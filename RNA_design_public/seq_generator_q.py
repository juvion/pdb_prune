import torch
from torch_geometric.data import Batch, Data
import os
from ema_pytorch import EMA
import numpy as np
from RIdiffusion import HEGNN, RIdiffusion
from generate_graph_ss import prepare_graph, pdb2graph
import argparse
from hyp_utils import HNN
from qcbm3 import QCBM
parser = argparse.ArgumentParser()


# pdb_dir = "./input_pdb/demo"
# pdb_dir="./input_pdb/CASP15_structure"

parser.add_argument("--pdb_dir", type=str, default="./input_pdb/demo", help="Path to input PDB directory")
args = parser.parse_args()
pdb_dir= args.pdb_dir

print("pdb_dir:", pdb_dir)

nt_types = ['A', 'U', 'G', 'C']
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# 统一定义 ckpt_path 并据此推导 run_dir
ckpt_path = './runs/RIdiffusion_q_debug_dataset_0.9_20251122-203109/best_base.pt'
ckpt = torch.load(ckpt_path, map_location=device)
run_dir = os.path.dirname(ckpt_path)

config = ckpt['config']
config['noise_type'] = 'uniform'
# Auto-align num_qubits with checkpoint prior_embedding in_features to avoid size mismatch
ema_sd = ckpt.get('ema', None)
expected_num_qubits = config.get('num_qubits', 9)
if ema_sd is not None:
    if 'ema_model.prior_embedding.weight' in ema_sd:
        expected_num_qubits = int(ema_sd['ema_model.prior_embedding.weight'].shape[1])
    elif 'online_model.prior_embedding.weight' in ema_sd:
        expected_num_qubits = int(ema_sd['online_model.prior_embedding.weight'].shape[1])
config['num_qubits'] = int(expected_num_qubits)
print(f"Using num_qubits={config['num_qubits']} (inferred from checkpoint)")




gnn = HEGNN(config, input_feat_dim=config['input_feat_dim'], hidden_channels=config['hidden_dim'], edge_attr_dim=config['edge_attr_dim'], dropout=config['drop_out'], n_layers=config['depth'], update_edge=config['update_edge'], embedding=config['embedding'], embedding_dim=config['embedding_dim'], embed_ss=config['embed_ss'], norm_feat=config['norm_feat'])

diffusion = RIdiffusion(model=gnn, config=config, timesteps=config['timesteps']).to(device)
diffusion = EMA(diffusion)

diffusion.load_state_dict(ckpt['ema'])
diffusion = diffusion.to(device)
# 加载训练阶段保存的 prior_top 和 test_ids，替代 QCBM 采样
prior_path = os.path.join(run_dir, 'saved_priors','prior_top_testset.pt')
prior_top = None
test_ids = None
try:
    saved = torch.load(prior_path, map_location='cpu')
    prior_top = saved.get('prior_top', None)
    test_ids = saved.get('test_ids', None)
    if isinstance(test_ids, torch.Tensor):
        test_ids = test_ids.tolist()
    if prior_top is None or test_ids is None:
        print(f"Error: missing prior_top or test_ids in {prior_path}")
    else:
        prior_top = prior_top.to(device=device, dtype=torch.float32)
        if prior_top.shape[1] != config['num_qubits']:
            print(f"Warning: prior_top dim {prior_top.shape[1]} != num_qubits {config['num_qubits']}")
        if prior_top.shape[0] != len(test_ids):
            print(f"Warning: prior_top count {prior_top.shape[0]} != test_ids count {len(test_ids)}")
except Exception as e:
    print(f"Error loading prior_top from {prior_path}: {e}")


def prepare_graph(data):
    # Safeguard: handle None input and optional attributes
    if data is None:
        return None
    # Safely remove optional heavy attributes if present
    try:
        if 'distances' in data:
            del data['distances']
    except Exception:
        if hasattr(data, 'distances'):
            delattr(data, 'distances')
    try:
        if 'edge_dist' in data:
            del data['edge_dist']
    except Exception:
        if hasattr(data, 'edge_dist'):
            delattr(data, 'edge_dist')

    mu_r_norm = getattr(data, 'mu_r_norm', None)
    if mu_r_norm is not None:
        extra_x_feature = torch.cat([data.x[:, 4:], mu_r_norm], dim=1)
    else:
        # Fallback: if mu_r_norm missing, only use existing extra features
        extra_x_feature = data.x[:, 4:]

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






pdb_files = [os.path.join(pdb_dir, f) for f in os.listdir(pdb_dir) if f.endswith('.pdb')]
# 与 prior_top 对齐的数量，尽量一一对应
pair_count = len(pdb_files)
if prior_top is None:
    print("Error: prior_top not loaded; cannot perform prior-guided generation.")
    pair_count = 0
else:
    if prior_top.shape[0] != pair_count:
        print(f"Warning: {pair_count} PDB files vs {prior_top.shape[0]} prior samples; pairing first {min(pair_count, prior_top.shape[0])} by index.")
        pair_count = min(pair_count, prior_top.shape[0])

output_file = './pdb.fasta'
error_files = []
all_recoveries = []

with open(output_file, 'w') as file:
    for idx, pdb_file in enumerate(pdb_files[:pair_count]):
        # try:
        graph = pdb2graph(pdb_file, normalize_path='./mean_attr.pt')
        if graph is None:
            print(f"Error: pdb2graph returned None for {pdb_file}")
            error_files.append(pdb_file)
            continue
        prepared_graph = prepare_graph(graph)
        if prepared_graph is None:
            print(f"Error: prepare_graph returned None for {pdb_file}")
            error_files.append(pdb_file)
            continue
        input_graph = Batch.from_data_list([prepared_graph]).to(device)
        original_sequence = ''.join([nt_types[idx.item()] for idx in input_graph.x.argmax(dim=1)])

        pdb_filename = os.path.splitext(os.path.basename(pdb_file))[0]
        file.write(f'>Original_{pdb_filename}\n')
        file.write(f'{original_sequence}\n')
        
        print(f'>Original_{pdb_filename}')
        print(original_sequence)
        
        # 针对当前 PDB，固定其先验比特（来自 prior_top[idx]），在 100 次迭代内复用
        if prior_top is None or idx >= prior_top.shape[0]:
            print("Error: prior_top unavailable or out of range; skip this PDB.")
            continue
        fixed_prior_bits_graph = prior_top[idx].unsqueeze(0)  # [1, num_qubits]
        
        recoveries = []
        for i in range(10):
            # 使用固定的先验比特进行采样（不再进行 QCBM 采样）
            prob, sample_graph = diffusion.ema_model.ddim_sample(
                input_graph,
                prior_bits_per_graph=fixed_prior_bits_graph,
                w=config.get('prior_weight', 0.5)
            )
            sampled_sequence = ''.join([nt_types[idx.item()] for idx in sample_graph.argmax(dim=1)])
            recovery = (prob.argmax(dim=1) == input_graph.x.argmax(dim=1)).sum().item() / input_graph.x.shape[0]
            recoveries.append(recovery)
            all_recoveries.append(recovery)

            print(f'Iteration {i+1}: Sampled: {sampled_sequence}')
            print(f'Iteration {i+1}: Original: {original_sequence}')
            print(f'Iteration {i+1}: Recovery: {recovery}\n')

            file.write(f'>seq{i+1}_{pdb_filename}--{recovery}\n')
            file.write(f'{sampled_sequence}\n')
        
        avg_recovery = sum(recoveries) / len(recoveries) if recoveries else 0.0
        print(f'Average Recovery over {len(recoveries)} iterations: {avg_recovery}\n')
    # overall average across all pdb files and all iterations
    overall_avg_recovery = sum(all_recoveries) / len(all_recoveries) if all_recoveries else 0.0
    print(f'Overall Average Recovery across {len(all_recoveries)} iterations (all PDBs): {overall_avg_recovery}\n')
        
        
        # except Exception as e:
        #     print(f"Error processing {pdb_file}: {e}")
        #     error_files.append(pdb_file)

if error_files:
    error_log_file = './error_log.txt'
    with open(error_log_file, 'w') as ef:
        for error_file in error_files:
            ef.write(f"{error_file}\n")
    print(f"Errors occurred in {len(error_files)} files. Check {error_log_file} for details.")
