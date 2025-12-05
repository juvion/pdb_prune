import torch
from torch_geometric.data import Batch, Data
import os
from ema_pytorch import EMA
import numpy as np
from RIdiffusion import HEGNN, RIdiffusion, seq_recovery
from generate_graph_ss import prepare_graph, pdb2graph
import argparse
from hyp_utils import HNN
from qcbm3 import QCBM
import torch.nn.functional as F
parser = argparse.ArgumentParser()


# pdb_dir = "./input_pdb/demo"
# pdb_dir="./input_pdb/CASP15_structure"

parser.add_argument("--pdb_dir", type=str, default="./input_pdb/demo", help="Path to input PDB directory")
# 新增：对齐训练评估的variant与prior权重
parser.add_argument("--variant", type=str, default="top", choices=["top", "mid", "tail", "mixed", "baseline"], help="Checkpoint variant to evaluate against")
parser.add_argument("--prior_weight", type=float, default=0.5, help="Fusion weight for QCBM prior bias (align with training)")
parser.add_argument("--ckpt_path", type=str, default="./runs/RIdiffusion_q_debug_dataset_0.9_20251122-230319/best_mid.pt", help="Path to checkpoint to evaluate")
args = parser.parse_args()
pdb_dir= args.pdb_dir

print("pdb_dir:", pdb_dir)
print("variant:", args.variant)
print("prior_weight:", args.prior_weight)
print("ckpt_path:", args.ckpt_path)

nt_types = ['A', 'U', 'G', 'C']
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# 统一定义 ckpt_path 并据此推导 run_dir
ckpt_path = args.ckpt_path
try:
    ckpt = torch.load(ckpt_path, map_location=device)
except Exception as e:
    print(f"Error loading checkpoint {ckpt_path}: {e}")
    raise
run_dir = os.path.dirname(ckpt_path)
# 若存在所选variant的checkpoint，则优先加载该variant
variant_ckpt_path = os.path.join(run_dir, f"best_{args.variant}.pt")
if os.path.exists(variant_ckpt_path):
    try:
        ckpt = torch.load(variant_ckpt_path, map_location=device)
        ckpt_path = variant_ckpt_path
        print(f"Loaded variant checkpoint: {variant_ckpt_path}")
    except Exception as e:
        print(f"Warning: failed to load variant checkpoint {variant_ckpt_path}, fallback to default: {e}")

config = ckpt['config']
config['noise_type'] = 'uniform'
# 更稳健地推断 num_qubits：优先从 model 的 prior_embedding 权重维度读取
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
    print(f"Warning: failed to infer num_qubits from checkpoint weights: {e}")
config['num_qubits'] = int(expected_num_qubits)
print(f"Using num_qubits={config['num_qubits']} (inferred from checkpoint)")




gnn = HEGNN(config, input_feat_dim=config['input_feat_dim'], hidden_channels=config['hidden_dim'], edge_attr_dim=config['edge_attr_dim'], dropout=config['drop_out'], n_layers=config['depth'], update_edge=config['update_edge'], embedding=config['embedding'], embedding_dim=config['embedding_dim'], embed_ss=config['embed_ss'], norm_feat=config['norm_feat'])

# 使用与训练一致：直接加载model权重，不包EMA
diffusion = RIdiffusion(model=gnn, config=config, timesteps=config['timesteps']).to(device)
try:
    diffusion.load_state_dict(ckpt['model'])
    print("Loaded model state_dict from checkpoint['model']")
except Exception as e:
    print(f"Warning: failed to load model params from ckpt['model']: {e}")

diffusion = diffusion.to(device)
# 评估模式与禁用梯度（与训练评估一致）
diffusion.eval()
torch.set_grad_enabled(False)
# 加载训练阶段保存的 prior_top 和 test_ids，替代 QCBM 采样
prior_path = os.path.join(run_dir, 'saved_priors', args.variant,'prior_top_testset.pt')
prior_top = None
test_ids = None
saved = {}
try:
    saved = torch.load(prior_path, map_location='cpu')
except Exception as e_var:
    print(f"Warning: failed to load variant prior at {prior_path}, try baseline: {e_var}")
    prior_path_baseline = os.path.join(run_dir, 'saved_priors', 'baseline','prior_top_testset.pt')
    try:
        saved = torch.load(prior_path_baseline, map_location='cpu')
        print(f"Loaded baseline prior from {prior_path_baseline}")
    except Exception as e_base:
        print(f"Error loading baseline prior at {prior_path_baseline}: {e_base}")
        saved = {}
# 读取 prior_top/test_ids 并转换类型
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
# 基于保存的 test_ids 构建从样本基名到 prior 索引的映射
id_basename_to_idx = {}
if test_ids:
    try:
        id_basename_to_idx = {os.path.splitext(os.path.basename(tid))[0]: i for i, tid in enumerate(test_ids)}
    except Exception as e:
        print(f"Error building test_ids mapping: {e}")

# 使用与训练一致的测试集图文件进行评估，避免 PDB 转图的预处理差异
if prior_top is None or test_ids is None or len(test_ids) == 0:
    print("Error: prior_top or test_ids not available; cannot evaluate on test set graphs.")
else:
    graph_root = config.get('graph_root', 'graph_dataset/dataset_0.9')
    test_graph_paths = [os.path.join(graph_root, tid) for tid in test_ids]

    output_file = './pdb.fasta'
    error_files = []
    processed_count = 0
    # 节点级加权平均（与训练一致）
    total_correct = 0
    total_nodes = 0

    with open(output_file, 'w') as file:
        for i, gpath in enumerate(test_graph_paths):
            try:
                data = torch.load(gpath, map_location=device)
            except Exception as e:
                print(f"Error loading graph {gpath}: {e}")
                error_files.append(gpath)
                continue
            if data is None:
                print(f"Error: Loaded None for {gpath}")
                error_files.append(gpath)
                continue
            # 使用与训练一致的数据结构：裁剪到4类并保留额外特征
            prepared = prepare_graph(data)
            input_graph = Batch.from_data_list([prepared]).to(device)
            original_sequence = ''.join([nt_types[ix.item()] for ix in input_graph.x[:, :4].argmax(dim=1)])

            pdb_filename = os.path.splitext(os.path.basename(gpath))[0]
            file.write(f'>Original_{pdb_filename}\n')
            file.write(f'{original_sequence}\n')
            print(f'>Original_{pdb_filename}')
            print(original_sequence)

            # 通过索引与 prior_top 一一对应
            if i >= prior_top.shape[0]:
                print(f"Warning: prior_top out of range for index {i}; skip {gpath}.")
                error_files.append(gpath)
                continue
            fixed_prior_bits_graph = prior_top[i].unsqueeze(0)

            with torch.no_grad():
                loss, pred_logits = diffusion.forward(
                    input_graph,
                    logit=True,
                    prior_bits_per_graph=fixed_prior_bits_graph,
                    w=args.prior_weight
                )
            pred_seq = F.one_hot(pred_logits.argmax(dim=1), num_classes=4).float()
            recovery, ind = seq_recovery(input_graph, pred_seq)
            predicted_sequence = ''.join([nt_types[ix.item()] for ix in pred_logits.argmax(dim=1)])
            # 累计节点级统计
            total_correct += ind.sum().item()
            total_nodes += ind.shape[0]

            print(f'Predicted: {predicted_sequence}')
            print(f'Original: {original_sequence}')
            print(f'Recovery: {float(recovery):.6f}\n')

            file.write(f'>Pred_{pdb_filename}--{float(recovery):.6f}\n')
            file.write(f'{predicted_sequence}\n')
            processed_count += 1
        overall_avg_recovery = (total_correct / max(total_nodes, 1)) if total_nodes > 0 else 0.0
        print(f'Overall Node-Weighted Recovery across {total_nodes} nodes: {overall_avg_recovery:.6f}\n')
        print(f'汇总：成功处理的序列数 {processed_count} / {len(test_graph_paths)}，跳过 {len(error_files)}')

    if error_files:
        error_log_file = './error_log.txt'
        with open(error_log_file, 'w') as ef:
            for error_file in error_files:
                ef.write(f"{error_file}\n")
        print(f"Errors occurred in {len(error_files)} files. Check {error_log_file} for details.")

if error_files:
    error_log_file = './error_log.txt'
    with open(error_log_file, 'w') as ef:
        for error_file in error_files:
            ef.write(f"{error_file}\n")
    print(f"Errors occurred in {len(error_files)} files. Check {error_log_file} for details.")
