import os
import argparse
import torch
import torch.nn.functional as F
from torch_geometric.data import Batch, Data
from RIdiffusion import HEGNN, RIdiffusion, seq_recovery

nt_types = ['A', 'U', 'G', 'C']

def prepare_graph(data: Data) -> Data:
    # Remove heavy optional attrs safely
    for attr in ['distances', 'edge_dist']:
        try:
            if attr in data:
                del data[attr]
        except Exception:
            if hasattr(data, attr):
                delattr(data, attr)
    mu_r_norm = getattr(data, 'mu_r_norm', None)
    extra_x_feature = data.x[:, 4:] if mu_r_norm is None else torch.cat([data.x[:, 4:], mu_r_norm], dim=1)
    graph = Data(
        x=data.x[:, :4],
        extra_x=extra_x_feature,
        pos=data.pos,
        edge_index=data.edge_index,
        edge_attr=data.edge_attr,
        ss=data.ss[:data.x.shape[0], :] if hasattr(data, 'ss') else None,
        sasa=data.x[:, 4] if data.x.shape[1] > 4 else None
    )
    return graph


def load_prior(run_dir: str, variant: str):
    # Try variant prior first, fallback to baseline
    var_path = os.path.join(run_dir, 'saved_priors', variant, 'prior_top_testset.pt')
    base_path = os.path.join(run_dir, 'saved_priors', 'baseline', 'prior_top_testset.pt')
    saved = {}
    try:
        saved = torch.load(var_path, map_location='cpu')
        print(f"Loaded prior from {var_path}")
    except Exception as e:
        print(f"Warning: failed to load {var_path}: {e}; fallback to baseline.")
        try:
            saved = torch.load(base_path, map_location='cpu')
            print(f"Loaded baseline prior from {base_path}")
        except Exception as e2:
            print(f"Error: failed to load baseline prior {base_path}: {e2}")
            saved = {}
    prior_top = saved.get('prior_top', None)
    test_ids = saved.get('test_ids', None)
    if isinstance(test_ids, torch.Tensor):
        test_ids = test_ids.tolist()
    return prior_top, test_ids


def evaluate_variant_seq_recovery(diffusion: RIdiffusion, prior_top: torch.Tensor, test_ids, graph_root: str, w: float):
    if prior_top is None or test_ids is None or len(test_ids) == 0:
        print("Error: prior_top or test_ids not available; skip evaluation.")
        return float('nan')
    device = next(diffusion.parameters()).device
    total_seq_rec = 0.0
    count = 0
    with torch.no_grad():
        for i, tid in enumerate(test_ids):
            gpath = os.path.join(graph_root, tid)
            try:
                data = torch.load(gpath, map_location=device)
            except Exception as e:
                print(f"Error loading graph {gpath}: {e}")
                continue
            prepared = prepare_graph(data)
            input_graph = Batch.from_data_list([prepared]).to(device)
            if i >= prior_top.shape[0]:
                print(f"Warning: prior_top out of range at {i}; skip.")
                continue
            fixed_prior_bits_graph = prior_top[i].unsqueeze(0).to(device=device, dtype=torch.float32)
            loss, pred_logits = diffusion.forward(
                input_graph,
                logit=True,
                prior_bits_per_graph=fixed_prior_bits_graph,
                w=w
            )
            pred_seq = F.one_hot(pred_logits.argmax(dim=1), num_classes=4).float()
            recovery, ind = seq_recovery(input_graph, pred_seq)
            total_seq_rec += float(recovery)
            count += 1
    return (total_seq_rec / count) if count > 0 else float('nan')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ckpt_path', type=str, required=True, help='A checkpoint path within the run directory')
    parser.add_argument('--prior_weight', type=float, default=0.5, help='Fusion weight for prior bias w')
    parser.add_argument('--graph_root', type=str, default='graph_dataset/dataset_0.9', help='Root for test graphs; joined with test_ids')
    parser.add_argument('--variants', type=str, default='top,mid,tail,mixed', help='Comma-separated variants to evaluate')
    args = parser.parse_args()

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    # Load base checkpoint to build model/config
    try:
        ckpt = torch.load(args.ckpt_path, map_location=device)
    except Exception as e:
        print(f"Error loading checkpoint {args.ckpt_path}: {e}")
        return
    run_dir = os.path.dirname(args.ckpt_path)

    config = ckpt['config']
    config['noise_type'] = 'uniform'
    # Infer num_qubits from prior_embedding weight if possible
    expected_num_qubits = config.get('num_qubits', 9)
    model_sd = ckpt.get('model', None)
    ema_sd = ckpt.get('ema', None)
    try:
        if model_sd is not None and 'prior_embedding.weight' in model_sd:
            expected_num_qubits = int(model_sd['prior_embedding.weight'].shape[1])
        elif ema_sd is not None:
            if 'ema_model.prior_embedding.weight' in ema_sd:
                expected_num_qubits = int(ema_sd['ema_model.prior_embedding.weight'].shape[1])
            elif 'online_model.prior_embedding.weight' in ema_sd:
                expected_num_qubits = int(ema_sd['online_model.prior_embedding.weight'].shape[1])
    except Exception as e:
        print(f"Warning: failed to infer num_qubits: {e}")
    config['num_qubits'] = int(expected_num_qubits)
    print(f"Using num_qubits={config['num_qubits']}")

    gnn = HEGNN(config, input_feat_dim=config['input_feat_dim'], hidden_channels=config['hidden_dim'], edge_attr_dim=config['edge_attr_dim'], dropout=config['drop_out'], n_layers=config['depth'], update_edge=config['update_edge'], embedding=config['embedding'], embedding_dim=config['embedding_dim'], embed_ss=config['embed_ss'], norm_feat=config['norm_feat'])
    diffusion = RIdiffusion(model=gnn, config=config, timesteps=config['timesteps']).to(device)

    variants = [v.strip() for v in args.variants.split(',') if v.strip()]

    for variant in variants:
        # Load variant checkpoint if available
        variant_ckpt_path = os.path.join(run_dir, f"best_{variant}.pt")
        try:
            ckpt_v = torch.load(variant_ckpt_path, map_location=device)
            print(f"Loaded variant checkpoint: {variant_ckpt_path}")
        except Exception as e:
            print(f"Warning: failed to load {variant_ckpt_path}: {e}; fallback to base ckpt.")
            ckpt_v = ckpt
        try:
            diffusion.load_state_dict(ckpt_v['model'])
        except Exception as e:
            print(f"Warning: failed to load model state_dict for variant {variant}: {e}")
        diffusion.eval()
        torch.set_grad_enabled(False)

        # Load prior_top and aligned test_ids for this variant
        prior_top, test_ids = load_prior(run_dir, variant)
        if prior_top is not None:
            prior_top = prior_top.to(device=device, dtype=torch.float32)
            if prior_top.shape[1] != config['num_qubits']:
                print(f"Warning: prior_top dim {prior_top.shape[1]} != num_qubits {config['num_qubits']}")
        graph_root = args.graph_root or config.get('graph_root', 'graph_dataset/dataset_0.9')
        avg_seq_rec = evaluate_variant_seq_recovery(diffusion, prior_top, test_ids, graph_root, w=args.prior_weight)
        print(f"[Prior-Fusion] Variant={variant} | Seq-Level Recovery: {avg_seq_rec:.4f}")

if __name__ == '__main__':
    main()