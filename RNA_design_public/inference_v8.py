import os
import json
import csv
import argparse
from datetime import datetime

import torch
import torch.nn.functional as F
from torch_geometric.loader import DataLoader

from ema_pytorch import EMA

from RIdiffusion import HEGNN, RIdiffusion, seq_recovery
from dataset_src.large_dataset import RNAsolo
from layers.quantum_prior import QuantumPriorBias


def list_graph_ids(base_dir: str, split: str):
    split_dir = os.path.join(base_dir, split)
    if not os.path.isdir(split_dir):
        raise FileNotFoundError(f"Graph split directory not found: {split_dir}. Run generate_graph_ss.py first.")
    files = [f for f in os.listdir(split_dir) if f.endswith('.pt')]
    return [os.path.join(split, f) for f in files]


def build_config():
    return {
        'manifold': 'PoincareBall',
        'num_layers': 2,
        'c': 1.0,
        'dropout': 0.0,

        'act': 'relu',
        'dim': 128,
        'bias': True,
        'task': 'rec',
        'cuda': -1,
        'device': 'cuda',

        'hidden_dim': 128,
        'depth': 4,
        'drop_out': 0.1,
        'update_edge': True,
        'embedding': True,
        'embedding_dim': 128,
        'embed_ss': -1,
        'norm_feat': False,

        'input_feat_dim': 128,
        'edge_attr_dim': 128,

        'timesteps': 500,
        'noise_type': 'uniform',
        'loss_type': 'CE',
        'objective': 'pred_x0',
    }


def decode_base(idx: int) -> str:
    # Dataset uses order ['A','U','G','C']
    mapping = ['A', 'U', 'G', 'C']
    try:
        return mapping[idx]
    except Exception:
        return 'N'


@torch.no_grad()
def run_inference(args):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    # Prepare dataset
    graph_ids = list_graph_ids(args.graph_root, args.split)
    ds = RNAsolo(graph_ids, baseDIR=args.graph_root + '/')
    loader = DataLoader(ds, batch_size=args.batch_size, shuffle=False)

    # Build model
    config = build_config()
    gnn = HEGNN(
        config,
        input_feat_dim=config['input_feat_dim'],
        hidden_channels=config['hidden_dim'],
        edge_attr_dim=config['edge_attr_dim'],
        dropout=config['drop_out'],
        n_layers=config['depth'],
        embedding=config['embedding'],
        embedding_dim=config['embedding_dim'],
        embed_ss=config['embed_ss'],
        norm_feat=config['norm_feat'],
    )

    diffusion = RIdiffusion(
        model=gnn,
        timesteps=config['timesteps'],
        loss_type=config['loss_type'],
        objective=config['objective'],
        config=config,
    ).to(device)

    # Load checkpoint
    ckpt = torch.load(args.ckpt_path, map_location=device)
    diffusion.load_state_dict(ckpt['model'])

    # EMA (optional)
    model_for_eval = diffusion
    if args.use_ema and ('ema' in ckpt):
        try:
            ema = EMA(diffusion, beta=args.ema_decay)
            ema.load_state_dict(ckpt['ema'])
            model_for_eval = ema.ema_model
        except Exception:
            model_for_eval = diffusion

    # Prior bias module (optional – only if checkpoint contains it)
    prior_module = None
    if args.use_quantum_prior_bias and ('prior' in ckpt):
        sample = ds.get(0)
        extra_dim = int(sample.extra_x.shape[1])
        prior_module = QuantumPriorBias(
            extra_feat_dim=extra_dim,
            bottleneck_dim=args.prior_bottleneck,
            n_qubits=args.prior_n_qubits,
            n_layers=args.prior_n_layers,
        ).to(device)
        try:
            prior_module.load_state_dict(ckpt['prior'])
        except Exception:
            prior_module = None

    # Output directory
    os.makedirs(args.out_dir, exist_ok=True)
    csv_path = os.path.join(args.out_dir, 'inference.csv')
    fasta_path = os.path.join(args.out_dir, 'pred_sequences.fasta')
    metrics_path = os.path.join(args.out_dir, 'inference_metrics.json')

    # Prepare writers
    csv_file = open(csv_path, 'w', newline='')
    csv_writer = csv.writer(csv_file)
    csv_writer.writerow(['graph_id', 'node_idx', 'pred_idx', 'pred_char', 'prob_A', 'prob_U', 'prob_G', 'prob_C'])
    fasta_file = open(fasta_path, 'w') if args.save_fasta else None

    total_correct = 0
    total_nodes = 0

    # Inference loop (batch_size defaults to 1)
    for batch_i, batch in enumerate(loader):
        data = batch.to(device)

        # Get logits from model
        loss_dummy, pred_logits = model_for_eval.forward(data, logit=True)

        # Inject prior bias if available
        if prior_module is not None and hasattr(data, 'extra_x') and data.extra_x is not None:
            bias = prior_module(data.extra_x)
            pred_logits = pred_logits + args.prior_strength * bias

        probs = F.softmax(pred_logits, dim=1)
        pred_idx = probs.argmax(dim=1)

        # Write per-node outputs
        base_probs = ['A', 'U', 'G', 'C']
        for i in range(pred_idx.shape[0]):
            p = probs[i].detach().cpu().tolist()
            csv_writer.writerow([
                graph_ids[batch_i],
                i,
                int(pred_idx[i].item()),
                decode_base(int(pred_idx[i].item())),
                p[0], p[1], p[2], p[3],
            ])

        # FASTA per-graph (batch_size==1 assumed)
        if args.save_fasta:
            seq = ''.join(decode_base(int(idx)) for idx in pred_idx.detach().cpu().tolist())
            fasta_file.write(f">{graph_ids[batch_i]}\n{seq}\n")

        # Metric vs ground truth
        pred_seq_onehot = F.one_hot(pred_idx, num_classes=4).float()
        recovery, ind = seq_recovery(data, pred_seq_onehot)
        total_correct += ind.sum().item()
        total_nodes += ind.shape[0]

    # Close writers
    csv_file.close()
    if fasta_file is not None:
        fasta_file.close()

    metrics = {
        'num_graphs': len(graph_ids),
        'total_nodes': total_nodes,
        'recovery': (total_correct / max(total_nodes, 1)) if total_nodes > 0 else None,
        'ckpt_path': args.ckpt_path,
        'use_ema': bool(args.use_ema and ('ema' in ckpt)),
        'used_prior_bias': bool(prior_module is not None),
    }
    with open(metrics_path, 'w') as f:
        json.dump(metrics, f, indent=2)

    print(f"Inference completed. Outputs saved to: {args.out_dir}")
    print(f"- Predictions CSV: {csv_path}")
    if args.save_fasta:
        print(f"- FASTA: {fasta_path}")
    print(f"- Metrics: {metrics_path}")


def main():
    parser = argparse.ArgumentParser(description='v8 inference: load trained checkpoint and predict sequences')
    parser.add_argument('--graph_root', type=str, default='graph_dataset/dataset_0.8')
    parser.add_argument('--split', type=str, default='test_0.8')
    parser.add_argument('--batch_size', type=int, default=1)
    parser.add_argument('--ckpt_path', type=str, default='./weight/weight_v8.pt')
    parser.add_argument('--out_dir', type=str, default='./runs/inference_v8')
    # EMA
    parser.add_argument('--use_ema', type=int, default=1, choices=[0, 1])
    parser.add_argument('--ema_decay', type=float, default=0.995)
    # Prior bias
    parser.add_argument('--use_quantum_prior_bias', type=int, default=1, choices=[0, 1])
    parser.add_argument('--prior_bottleneck', type=int, default=8)
    parser.add_argument('--prior_n_qubits', type=int, default=4)
    parser.add_argument('--prior_n_layers', type=int, default=2)
    parser.add_argument('--prior_strength', type=float, default=0.2)
    # Output options
    parser.add_argument('--save_fasta', type=int, default=1, choices=[0, 1])

    args = parser.parse_args()

    # Stamp run dir with timestamp if default
    if args.out_dir.endswith('inference_v8'):
        ts = datetime.now().strftime('%Y%m%d-%H%M%S')
        args.out_dir = os.path.join(args.out_dir, ts)

    run_inference(args)


if __name__ == '__main__':
    main()