import os
import argparse
import torch
import torch.nn.functional as F
from torch_geometric.loader import DataLoader

from ema_pytorch import EMA

from RIdiffusion import HEGNN, RIdiffusion, seq_recovery
from dataset_src.large_dataset import RNAsolo


def list_graph_ids(base_dir: str, split: str):
    """
    Returns relative file paths (under base_dir) for .pt graphs in the given split.
    Example: base_dir='graph_dataset/dataset_0.8/', split='test_0.8' ->
    ['test_0.8/XYZ.pt', ...]
    """
    split_dir = os.path.join(base_dir, split)
    if not os.path.isdir(split_dir):
        raise FileNotFoundError(f"Graph split directory not found: {split_dir}. Run generate_graph_ss.py first.")
    files = [f for f in os.listdir(split_dir) if f.endswith('.pt')]
    return [os.path.join(split, f) for f in files]


def build_config_from_ckpt(ckpt: dict):
    if 'config' not in ckpt:
        raise KeyError("Checkpoint missing 'config'. Please re-train with config saved.")
    return ckpt['config']


def build_models_from_ckpt(ckpt_path: str, device: torch.device):
    ckpt = torch.load(ckpt_path, map_location=device)
    config = build_config_from_ckpt(ckpt)
    # Adjust config for CPU-only environments to prevent CUDA-only tensors in hyp modules
    if not torch.cuda.is_available():
        try:
            config = dict(config)  # shallow copy
            config['cuda'] = -1
            config['device'] = 'cpu'
        except Exception:
            pass

    # Recreate HEGNN and Diffusion as in train.py
    gnn = HEGNN(config,
                input_feat_dim=config['input_feat_dim'],
                hidden_channels=config['hidden_dim'],
                edge_attr_dim=config['edge_attr_dim'],
                dropout=config['drop_out'],
                n_layers=config['depth'],
                update_edge=config['update_edge'],
                embedding=config['embedding'],
                embedding_dim=config['embedding_dim'],
                embed_ss=config['embed_ss'],
                norm_feat=config['norm_feat'])

    diffusion = RIdiffusion(model=gnn,
                            timesteps=config['timesteps'],
                            loss_type=config['loss_type'],
                            objective=config['objective'],
                            config=config).to(device)

    diffusion.load_state_dict(ckpt['model'])

    ema_model = None
    if 'ema' in ckpt:
        # instantiate EMA and load saved state
        ema = EMA(diffusion, beta=0.995)
        try:
            ema.load_state_dict(ckpt['ema'])
            ema_model = ema.ema_model
        except Exception as e:
            print(f"[inference] Warning: failed to load EMA state: {e}")

    return diffusion, ema_model, config


@torch.no_grad()
def evaluate(dataloader, model, device):
    total_correct = 0
    total_nodes = 0
    y_true_all = []
    y_pred_all = []
    per_seq_rec_list = []

    model.eval()
    for batch in dataloader:
        data = batch.to(device)
        # forward with logits for recovery and F1
        _, pred_logits = model.forward(data, logit=True)
        pred = pred_logits.argmax(dim=1)

        # recovery (node-level)
        pred_seq = F.one_hot(pred, num_classes=4).float()
        _, ind = seq_recovery(data, pred_seq)
        total_correct += ind.sum().item()
        total_nodes += ind.shape[0]

        # recovery (sequence-level): mean over each graph in current batch
        if hasattr(data, 'batch'):
            unique_graphs = torch.unique(data.batch)
            for gid in unique_graphs:
                mask = (data.batch == gid)
                if mask.any():
                    per_seq_rec_list.append(ind[mask].float().mean().item())

        # labels for macro F1 (node-level classification)
        y_true = data.x[:, :4].argmax(dim=1)
        y_true_all.append(y_true.cpu())
        y_pred_all.append(pred.cpu())

    rec_node = total_correct / max(total_nodes, 1)
    rec_seq = (sum(per_seq_rec_list) / len(per_seq_rec_list)) if len(per_seq_rec_list) > 0 else float('nan')

    y_true_all = torch.cat(y_true_all).numpy() if len(y_true_all) > 0 else []
    y_pred_all = torch.cat(y_pred_all).numpy() if len(y_pred_all) > 0 else []

    f1_macro = float('nan')
    if len(y_true_all) > 0:
        try:
            from sklearn.metrics import precision_recall_fscore_support
            _, _, f1_macro, _ = precision_recall_fscore_support(
                y_true_all, y_pred_all, labels=[0, 1, 2, 3], average='macro'
            )
        except Exception as e:
            print(f"[inference] Warning: sklearn not available or error computing F1: {e}")

    return rec_node, rec_seq, f1_macro


def main():
    parser = argparse.ArgumentParser(description='Inference script: recovery score and macro F1 on test set')
    parser.add_argument('--ckpt_path', type=str,
                        default='weight/weight_custom.pt',
                        help='Path to checkpoint saved by train.py')
    parser.add_argument('--graph_root', type=str, default='graph_dataset/dataset_0.8',
                        help='Root folder of preprocessed graphs (generated by generate_graph_ss.py)')
    parser.add_argument('--test_split', type=str, default='test_0.8',
                        help='Subfolder under graph_root for test graphs')
    parser.add_argument('--batch_size', type=int, default=1, help='PyG DataLoader batch size')
    parser.add_argument('--use_ema', type=int, default=1, choices=[0, 1], help='Use EMA model parameters if available')
    args = parser.parse_args()

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    # Prepare test dataset
    try:
        test_ids = list_graph_ids(args.graph_root, args.test_split)
    except FileNotFoundError as e:
        print(str(e))
        print("Hint: run 'python3 generate_graph_ss.py' to produce test graphs.")
        return

    if len(test_ids) == 0:
        print(f"[inference] No .pt graphs found under {os.path.join(args.graph_root, args.test_split)}.\n"
              f"Please ensure you ran generate_graph_ss.py and that dataset_src/{args.test_split} contains PDB files.")
        return

    test_ds = RNAsolo(test_ids, baseDIR=args.graph_root + '/')
    test_loader = DataLoader(test_ds, batch_size=args.batch_size, shuffle=False)

    # Build model from checkpoint
    diffusion, ema_model, config = build_models_from_ckpt(args.ckpt_path, device)
    model = ema_model if (args.use_ema == 1 and ema_model is not None) else diffusion

    # Evaluate
    rec_node, rec_seq, f1_macro = evaluate(test_loader, model, device)

    # Output
    print(f"Test recovery (node-level): {rec_node:.4f}")
    if rec_seq == rec_seq:  # not NaN
        print(f"Test recovery (sequence-level): {rec_seq:.4f}")
    else:
        print("Test recovery (sequence-level): NaN (no test samples)")
    if not (f1_macro != f1_macro):  # check for NaN
        print(f"Test macro F1: {f1_macro:.4f}")
    else:
        print("Test macro F1: NaN (sklearn missing or no test samples)")


if __name__ == '__main__':
    main()