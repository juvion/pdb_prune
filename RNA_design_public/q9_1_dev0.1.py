# Q8: QCBM-RNA Design with Prior Sample Dataset and Sequence-level Recovery Rate
# 基于q7.py修改，主要变化：
# 1. 生成先验样本数据集，与训练集RNA序列一一对应
# 2. 修改RNAModelWithQCBM训练逻辑：先验样本embedding与RNA序列node embedding相加
# 3. 训练后按序列粒度计算recovery rate（参考public2.py逻辑）
# 4. 计算目标概率分布时候不再考虑没有采样的比特串样本，直接用softmax函数作用在唯一先验样本对应的分数向量上，
# 得到target_probs，因此target_probs向量的大小就是唯一先验样本的数量。电路预测的概率分布pred_probs则从qcbm_distribution
# 中提取唯一先验样本对应的概率数值，，得到最终的电路预测概率分布。
# 5. 计算电路预测概率分布与目标概率分布的交叉熵损失，作为训练QCBM的损失函数。

# 基本代码框架来自q8_6_fuxian.py
# 在edge embedding上添加先验样本


import torch
from dataclasses import dataclass, field
from torch.utils.data import Dataset
from torch.utils.data import DataLoader
import pandas as pd
import random
import os
import torch
import torch.nn as nn
import numpy as np
from torch_scatter import scatter_sum, scatter_softmax
import torch.nn.functional as F
from typing import List
from torch.optim import Adam
from scipy.optimize import minimize
from tqdm import tqdm
from Bio import SeqIO
import math
import pennylane as qml
from pennylane import numpy as pnp
import sys
sys.path.append('./package_score')
from sequence_similarity_score import RNASequenceScorerOptimized
from pdb_to_fasta import SequenceExtractor
import tempfile

# Define function to read FASTA files using Biopython
#def read_fasta_biopython(file_path):
#    sequences = {}
#    for record in SeqIO.parse(file_path, "fasta"):
#        sequences[record.id] = str(record.seq)
#    return sequences
def read_fasta_biopython(file_path):
    sequences = {}
    try:
        # Robust read: ignore non-UTF8 bytes and continue
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as handle:
            for record in SeqIO.parse(handle, "fasta"):
                sequences[record.id] = str(record.seq)
    except Exception as e:
        print(f"Warning: failed to read FASTA {file_path}: {e}")
    return sequences



# 对文件列表进行排序以确保可复现性
# 原始：train_file_list = sorted(os.listdir("./RNAdesignv1/train/seqs"))
# 过滤隐藏文件和非FASTA扩展名，避免解析到 AppleDouble 资源文件（比如 ._pdbXXXX.fasta）
allowed_exts = {".fasta", ".fa", ".fna", ".txt"}
train_seqs_dir = "./RNAdesignv1/train/seqs"
train_file_list = sorted([
    f for f in os.listdir(train_seqs_dir)
    if not f.startswith('.') and os.path.splitext(f)[1].lower() in allowed_exts
])
content_dict = {
    "pdb_id": [],
    "seq": []
}
for file in tqdm(train_file_list):
    file_path = os.path.join(train_seqs_dir, file)
    sequences = read_fasta_biopython(file_path)
    # 如果解析不到任何记录，则跳过该文件并给出警告
    if not sequences:
        print(f"Warning: no FASTA records parsed in {file}; skipping.")
        continue
    # 取第一条记录
    pdb_id, seq = next(iter(sequences.items()))
    content_dict["pdb_id"].append(pdb_id)
    content_dict["seq"].append(seq)

data = pd.DataFrame(content_dict)

# 数据划分函数，使用指定的随机种子
def split_data(data, split_seed):
    np.random.seed(split_seed)
    random.seed(split_seed)
    split = np.random.choice(['train', 'valid', 'test'], size=len(data), p=[0.7, 0.2, 0.1])
    data_copy = data.copy()
    data_copy['split'] = split
    train_data = data_copy[data_copy['split']=='train']
    valid_data = data_copy[data_copy['split']=='valid']
    test_data = data_copy[data_copy['split']=='test']
    return train_data, valid_data, test_data

@dataclass
class DataConfig:
    train_npy_data_dir: str = './RNAdesignv1/train/coords'
    train_data_path: str = 'public_train_data.csv'
    valid_npy_data_dir: str = './RNAdesignv1/train/coords'
    valid_data_path: str = 'public_valid_data.csv'
    test_npy_data_dir: str = './RNAdesignv1/train/coords'
    test_data_path: str = 'public_test_data.csv'

@dataclass
class ModelConfig:
    smoothing: float = 0.1
    hidden: int = 128
    vocab_size: int = 4
    k_neighbors: int = 30
    dropout: float = 0.1
    node_feat_types: List[str] = field(default_factory=lambda: ['angle', 'distance', 'direction'])
    edge_feat_types: List[str] = field(default_factory=lambda: ['orientation', 'distance', 'direction'])
    num_encoder_layers: int = 3
    num_decoder_layers: int = 3

@dataclass
class QCBMConfig:
    num_qubits: int = 11  # 量子比特数
    num_layers: int = 3  # 量子电路层数
    qcbm_epochs: int = 300  # QCBM训练轮数
    qcbm_lr: float = 0.05  # QCBM学习率

@dataclass
class TrainConfig:
    batch_size: int = 16
    epoch: int = 100
    lr: float = 0.001
    output_dir: str = 'ckpts/q7'
    ckpt_path: str = 'ckpts/q7/best_q7.pt'

@dataclass
class Config:
    pipeline: str = 'train'
    seed: int = 2025
    device: str = 'cuda:0' if torch.cuda.is_available() else 'cpu'
    data_config: DataConfig = DataConfig()
    model_config: ModelConfig = ModelConfig()
    qcbm_config: QCBMConfig = QCBMConfig()
    train_config: TrainConfig = TrainConfig()

# RNA数据集类
class RNADataset(Dataset):
    def __init__(self, data_path, npy_dir):
        super(RNADataset, self).__init__()
        self.data = pd.read_csv(data_path)
        self.npy_dir = npy_dir
        self.seq_list = self.data['seq'].to_list()
        self.name_list = self.data['pdb_id'].to_list()

    def __len__(self):
        return len(self.name_list)

    def __getitem__(self, idx):
        seq = self.seq_list[idx]
        pdb_id = self.name_list[idx]
        coords = np.load(os.path.join(self.npy_dir, pdb_id + '.npy'))

        feature = {
            "name": pdb_id,
            "seq": seq,
            "coords": {
                "P": coords[:, 0, :],
                "O5'": coords[:, 1, :],
                "C5'": coords[:, 2, :],
                "C4'": coords[:, 3, :],
                "C3'": coords[:, 4, :],
                "O3'": coords[:, 5, :],
            }
        }

        return feature
# 设置随机种子
def seeding(seed):
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

# 特征化函数
def featurize(batch):
    alphabet = 'AUCG'
    B = len(batch)
    # 原来的实现使用序列长度作为L_max，导致当坐标长度与序列长度不一致时出现广播错误
    # 这里改为对每个样本使用有效长度 = min(序列长度, 坐标长度)，并在批次中用这些有效长度的最大值作为L_max
    seq_lengths = np.array([len(b['seq']) for b in batch], dtype=np.int32)
    coord_lengths = np.array([len(b['coords']["P"]) for b in batch], dtype=np.int32)
    eff_lengths = np.minimum(seq_lengths, coord_lengths)
    L_max = int(eff_lengths.max())

    X = np.zeros([B, L_max, 6, 3], dtype=np.float32)
    S = np.zeros([B, L_max], dtype=np.int32)
    names = []

    # Build the batch
    for i, b in enumerate(batch):
        x = np.stack([np.nan_to_num(b['coords'][c], nan=0.0) for c in ["P", "O5'", "C5'", "C4'", "C3'", "O3'"]], 1)
        l_seq = len(b['seq'])
        l_coord = x.shape[0]
        l_eff = int(min(l_seq, l_coord))
        if l_seq != l_coord:
            print(f"Warning: length mismatch for {b['name']}: seq={l_seq}, coords={l_coord}. Using min={l_eff}.")
        # Trim to effective length
        x_trim = x[:l_eff, :, :]
        # Pad to batch L_max
        pad_len = L_max - l_eff
        if pad_len > 0:
            x_pad = np.pad(x_trim, [[0, pad_len], [0,0], [0,0]], 'constant', constant_values=(np.nan, ))
        else:
            x_pad = x_trim
        X[i,:,:,:] = x_pad
        # Sequence indices, trimmed to effective length
        seq_trim = b['seq'][:l_eff]
        try:
            indices = np.asarray([alphabet.index(a) for a in seq_trim], dtype=np.int32)
        except ValueError:
            # 处理非AUCG字符：回退到'A'（0）
            indices = np.asarray([alphabet.index(a) if a in alphabet else 0 for a in seq_trim], dtype=np.int32)
            print(f"Warning: non-AUCG character(s) found in {b['name']}; mapped to 'A'.")
        S[i, :l_eff] = indices
        names.append(b['name'])

    # Mask 和 NaN处理，与原始逻辑保持一致
    mask = np.isfinite(np.sum(X,(2,3))).astype(np.float32)
    numbers = np.sum(mask, axis=1).astype(np.int32)
    S_new = np.zeros_like(S)
    X_new = np.zeros_like(X)+np.nan
    for i, n in enumerate(numbers):
        X_new[i,:n,::] = X[i][mask[i]==1]
        S_new[i,:n] = S[i][mask[i]==1]

    X = X_new
    S = S_new
    isnan = np.isnan(X)
    mask = np.isfinite(np.sum(X,(2,3))).astype(np.float32)
    X[isnan] = 0.
    # Conversion
    S = torch.from_numpy(S).to(dtype=torch.long)
    X = torch.from_numpy(X).to(dtype=torch.float32)
    mask = torch.from_numpy(mask).to(dtype=torch.float32)
    return X, S, mask, eff_lengths, names


def gather_edges(edges, neighbor_idx):
    neighbors = neighbor_idx.unsqueeze(-1).expand(-1, -1, -1, edges.size(-1))
    return torch.gather(edges, 2, neighbors)

def gather_nodes(nodes, neighbor_idx):
    neighbors_flat = neighbor_idx.view((neighbor_idx.shape[0], -1))
    neighbors_flat = neighbors_flat.unsqueeze(-1).expand(-1, -1, nodes.size(2))
    neighbor_features = torch.gather(nodes, 1, neighbors_flat)
    neighbor_features = neighbor_features.view(list(neighbor_idx.shape)[:3] + [-1])
    return neighbor_features

def gather_nodes_t(nodes, neighbor_idx):
    idx_flat = neighbor_idx.unsqueeze(-1).expand(-1, -1, nodes.size(2))
    return torch.gather(nodes, 1, idx_flat)

def cat_neighbors_nodes(h_nodes, h_neighbors, E_idx):
    h_nodes = gather_nodes(h_nodes, E_idx)
    return torch.cat([h_neighbors, h_nodes], -1)

class MPNNLayer(nn.Module):
    def __init__(self, num_hidden, num_in, dropout=0.1, num_heads=None, scale=30):
        super(MPNNLayer, self).__init__()
        self.num_hidden = num_hidden
        self.num_in = num_in
        self.scale = scale
        self.dropout = nn.Dropout(dropout)
        self.norm1 = nn.LayerNorm(num_hidden)
        self.norm2 = nn.LayerNorm(num_hidden)

        self.W1 = nn.Linear(num_hidden + num_in, num_hidden, bias=True)
        self.W2 = nn.Linear(num_hidden, num_hidden, bias=True)
        self.W3 = nn.Linear(num_hidden, num_hidden, bias=True)
        self.act = nn.ReLU()

        self.dense = nn.Sequential(
            nn.Linear(num_hidden, num_hidden*4),
            nn.ReLU(),
            nn.Linear(num_hidden*4, num_hidden)
        )

    def forward(self, h_V, h_E, edge_idx, batch_id=None):
        src_idx, dst_idx = edge_idx[0], edge_idx[1]
        h_message = self.W3(self.act(self.W2(self.act(self.W1(h_E)))))
        dh = scatter_sum(h_message, src_idx, dim=0) / self.scale
        h_V = self.norm1(h_V + self.dropout(dh))
        dh = self.dense(h_V)
        h_V = self.norm2(h_V + self.dropout(dh))
        return h_V

class Normalize(nn.Module):
    def __init__(self, features, epsilon=1e-6):
        super(Normalize, self).__init__()
        self.gain = nn.Parameter(torch.ones(features))
        self.bias = nn.Parameter(torch.zeros(features))
        self.epsilon = epsilon

    def forward(self, x, dim=-1):
        mu = x.mean(dim, keepdim=True)
        sigma = torch.sqrt(x.var(dim, keepdim=True) + self.epsilon)
        gain = self.gain
        bias = self.bias
        if dim != -1:
            shape = [1] * len(mu.size())
            shape[dim] = self.gain.size()[0]
            gain = gain.view(shape)
            bias = bias.view(shape)
        return gain * (x - mu) / (sigma + self.epsilon) + bias

# 特征维度定义
feat_dims = {
    'node': {
        'angle': 12,
        'distance': 80,
        'direction': 9,
    },
    'edge': {
        'orientation': 4,
        'distance': 96,
        'direction': 15,
    }
}

def nan_to_num(tensor, nan=0.0):
    idx = torch.isnan(tensor)
    tensor[idx] = nan
    return tensor

def _normalize(tensor, dim=-1):
    return nan_to_num(
        torch.div(tensor, torch.norm(tensor, dim=dim, keepdim=True)))

class RNAFeatures(nn.Module):
    def __init__(self, edge_features, node_features, node_feat_types=[], edge_feat_types=[], num_rbf=16, top_k=30, augment_eps=0., dropout=0.1, args=None):
        super(RNAFeatures, self).__init__()
        """Extract RNA Features"""
        self.edge_features = edge_features
        self.node_features = node_features
        self.top_k = top_k
        self.augment_eps = augment_eps
        self.num_rbf = num_rbf
        self.dropout = nn.Dropout(dropout)
        self.node_feat_types = node_feat_types
        self.edge_feat_types = edge_feat_types

        node_in = sum([feat_dims['node'][feat] for feat in node_feat_types])
        edge_in = sum([feat_dims['edge'][feat] for feat in edge_feat_types])
        self.node_embedding = nn.Linear(node_in,  node_features, bias=True)
        self.edge_embedding = nn.Linear(edge_in, edge_features, bias=True)
        self.norm_nodes = Normalize(node_features)
        self.norm_edges = Normalize(edge_features)

    def _dist(self, X, mask, eps=1E-6):
        mask_2D = torch.unsqueeze(mask,1) * torch.unsqueeze(mask,2)
        dX = torch.unsqueeze(X,1) - torch.unsqueeze(X,2)
        D = (1. - mask_2D)*10000 + mask_2D* torch.sqrt(torch.sum(dX**2, 3) + eps)

        D_max, _ = torch.max(D, -1, keepdim=True)
        D_adjust = D + (1. - mask_2D) * (D_max+1)
        D_neighbors, E_idx = torch.topk(D_adjust, min(self.top_k, D_adjust.shape[-1]), dim=-1, largest=False)
        return D_neighbors, E_idx

    def _rbf(self, D):
        D_min, D_max, D_count = 0., 20., self.num_rbf
        D_mu = torch.linspace(D_min, D_max, D_count, device=D.device)
        D_mu = D_mu.view([1,1,1,-1])
        D_sigma = (D_max - D_min) / D_count
        D_expand = torch.unsqueeze(D, -1)
        return torch.exp(-((D_expand - D_mu) / D_sigma)**2)

    def _get_rbf(self, A, B, E_idx=None, num_rbf=16):
        if E_idx is not None:
            D_A_B = torch.sqrt(torch.sum((A[:,:,None,:] - B[:,None,:,:])**2,-1) + 1e-6)
            D_A_B_neighbors = gather_edges(D_A_B[:,:,:,None], E_idx)[:,:,:,0]
            RBF_A_B = self._rbf(D_A_B_neighbors)
        else:
            D_A_B = torch.sqrt(torch.sum((A[:,:,None,:] - B[:,:,None,:])**2,-1) + 1e-6)
            RBF_A_B = self._rbf(D_A_B)
        return RBF_A_B

    def _quaternions(self, R):
        diag = torch.diagonal(R, dim1=-2, dim2=-1)
        Rxx, Ryy, Rzz = diag.unbind(-1)
        magnitudes = 0.5 * torch.sqrt(torch.abs(1 + torch.stack([
              Rxx - Ryy - Rzz,
            - Rxx + Ryy - Rzz,
            - Rxx - Ryy + Rzz
        ], -1)))
        _R = lambda i,j: R[:,:,:,i,j]
        signs = torch.sign(torch.stack([
            _R(2,1) - _R(1,2),
            _R(0,2) - _R(2,0),
            _R(1,0) - _R(0,1)
        ], -1))
        xyz = signs * magnitudes
        w = torch.sqrt(F.relu(1 + diag.sum(-1, keepdim=True))) / 2.
        Q = torch.cat((xyz, w), -1)
        Q = F.normalize(Q, dim=-1)
        return Q

    def _orientations_coarse(self, X, E_idx, eps=1e-6):
        V = X.clone()
        X = X[:,:,:6,:].reshape(X.shape[0], 6*X.shape[1], 3)
        dX = X[:,1:,:] - X[:,:-1,:]
        U = _normalize(dX, dim=-1)
        u_0, u_1 = U[:,:-2,:], U[:,1:-1,:]
        n_0 = _normalize(torch.cross(u_0, u_1), dim=-1)
        b_1 = _normalize(u_0 - u_1, dim=-1)

        # select C3'
        n_0 = n_0[:,4::6,:]
        b_1 = b_1[:,4::6,:]
        X = X[:,4::6,:]

        Q = torch.stack((b_1, n_0, torch.cross(b_1, n_0)), 2)
        Q = Q.view(list(Q.shape[:2]) + [9])
        Q = F.pad(Q, (0,0,0,1), 'constant', 0) # [16, 464, 9]

        Q_neighbors = gather_nodes(Q, E_idx) # [16, 464, 30, 9]
        P_neighbors = gather_nodes(V[:,:,0,:], E_idx) # [16, 464, 30, 3]
        O5_neighbors = gather_nodes(V[:,:,1,:], E_idx)
        C5_neighbors = gather_nodes(V[:,:,2,:], E_idx)
        C4_neighbors = gather_nodes(V[:,:,3,:], E_idx)
        O3_neighbors = gather_nodes(V[:,:,5,:], E_idx)

        Q = Q.view(list(Q.shape[:2]) + [3,3]).unsqueeze(2) # [16, 464, 1, 3, 3]
        Q_neighbors = Q_neighbors.view(list(Q_neighbors.shape[:3]) + [3,3]) # [16, 464, 30, 3, 3]

        dX = torch.stack([P_neighbors,O5_neighbors,C5_neighbors,C4_neighbors,O3_neighbors], dim=3) - X[:,:,None,None,:] # [16, 464, 30, 3]
        dU = torch.matmul(Q[:,:,:,None,:,:], dX[...,None]).squeeze(-1) # [16, 464, 30, 3] 邻居的相对坐标
        B, N, K = dU.shape[:3]
        E_direct = _normalize(dU, dim=-1)
        E_direct = E_direct.reshape(B, N, K,-1)
        R = torch.matmul(Q.transpose(-1,-2), Q_neighbors)
        E_orient = self._quaternions(R)

        dX_inner = V[:,:,[0,2,3],:] - X.unsqueeze(-2)
        dU_inner = torch.matmul(Q, dX_inner.unsqueeze(-1)).squeeze(-1)
        dU_inner = _normalize(dU_inner, dim=-1)
        V_direct = dU_inner.reshape(B,N,-1)
        return V_direct, E_direct, E_orient

    def _dihedrals(self, X, eps=1e-7):
        # P, O5', C5', C4', C3', O3'
        X = X[:,:,:6,:].reshape(X.shape[0], 6*X.shape[1], 3)

        # Shifted slices of unit vectors
        # https://iupac.qmul.ac.uk/misc/pnuc2.html#220
        # https://x3dna.org/highlights/torsion-angles-of-nucleic-acid-structures
        # alpha:   O3'_{i-1} P_i O5'_i C5'_i
        # beta:    P_i O5'_i C5'_i C4'_i
        # gamma:   O5'_i C5'_i C4'_i C3'_i
        # delta:   C5'_i C4'_i C3'_i O3'_i
        # epsilon: C4'_i C3'_i O3'_i P_{i+1}
        # zeta:    C3'_i O3'_i P_{i+1} O5'_{i+1}
        # What's more:
        #   chi: C1' - N9
        #   chi is different for (C, T, U) and (A, G) https://x3dna.org/highlights/the-chi-x-torsion-angle-characterizes-base-sugar-relative-orientation

        dX = X[:, 5:, :] - X[:, :-5, :] # O3'-P, P-O5', O5'-C5', C5'-C4', ...
        U = F.normalize(dX, dim=-1)
        u_2 = U[:,:-2,:]  # O3'-P, P-O5', ...
        u_1 = U[:,1:-1,:] # P-O5', O5'-C5', ...
        u_0 = U[:,2:,:]   # O5'-C5', C5'-C4', ...
        # Backbone normals
        n_2 = F.normalize(torch.cross(u_2, u_1), dim=-1)
        n_1 = F.normalize(torch.cross(u_1, u_0), dim=-1)

        # Angle between normals
        cosD = (n_2 * n_1).sum(-1)
        cosD = torch.clamp(cosD, -1+eps, 1-eps)
        D = torch.sign((u_2 * n_1).sum(-1)) * torch.acos(cosD)

        D = F.pad(D, (3,4), 'constant', 0)
        D = D.view((D.size(0), D.size(1) //6, 6))
        return torch.cat((torch.cos(D), torch.sin(D)), 2) # return D_features

    def forward(self, X, S, mask):
        if self.training and self.augment_eps > 0:
            X = X + self.augment_eps * torch.randn_like(X)

        # Build k-Nearest Neighbors graph
        B, N, _,_ = X.shape
        # P, O5', C5', C4', C3', O3'
        atom_P = X[:, :, 0, :]
        atom_O5_ = X[:, :, 1, :]
        atom_C5_ = X[:, :, 2, :]
        atom_C4_ = X[:, :, 3, :]
        atom_C3_ = X[:, :, 4, :]
        atom_O3_ = X[:, :, 5, :]

        X_backbone = atom_P
        D_neighbors, E_idx = self._dist(X_backbone, mask)

        mask_bool = (mask==1)
        mask_attend = gather_nodes(mask.unsqueeze(-1), E_idx).squeeze(-1)
        mask_attend = (mask.unsqueeze(-1) * mask_attend) == 1
        edge_mask_select = lambda x: torch.masked_select(x, mask_attend.unsqueeze(-1)).reshape(-1,x.shape[-1])
        node_mask_select = lambda x: torch.masked_select(x, mask_bool.unsqueeze(-1)).reshape(-1, x.shape[-1])

        # node features
        h_V = []
        # angle
        V_angle = node_mask_select(self._dihedrals(X))
        # distance
        node_list = ['O5_-P', 'C5_-P', 'C4_-P', 'C3_-P', 'O3_-P']
        V_dist = []

        for pair in node_list:
            atom1, atom2 = pair.split('-')
            V_dist.append(node_mask_select(self._get_rbf(vars()['atom_' + atom1], vars()['atom_' + atom2], None, self.num_rbf).squeeze()))
        V_dist = torch.cat(tuple(V_dist), dim=-1).squeeze()
        # direction
        V_direct, E_direct, E_orient = self._orientations_coarse(X, E_idx)
        V_direct = node_mask_select(V_direct)
        E_direct, E_orient = list(map(lambda x: edge_mask_select(x), [E_direct, E_orient]))

        # edge features
        h_E = []
        # dist
        edge_list = ['P-P', 'O5_-P', 'C5_-P', 'C4_-P', 'C3_-P', 'O3_-P']
        E_dist = []
        for pair in edge_list:
            atom1, atom2 = pair.split('-')
            E_dist.append(edge_mask_select(self._get_rbf(vars()['atom_' + atom1], vars()['atom_' + atom2], E_idx, self.num_rbf)))
        E_dist = torch.cat(tuple(E_dist), dim=-1)

        if 'angle' in self.node_feat_types:
            h_V.append(V_angle)
        if 'distance' in self.node_feat_types:
            h_V.append(V_dist)
        if 'direction' in self.node_feat_types:
            h_V.append(V_direct)

        if 'orientation' in self.edge_feat_types:
            h_E.append(E_orient)
        if 'distance' in self.edge_feat_types:
            h_E.append(E_dist)
        if 'direction' in self.edge_feat_types:
            h_E.append(E_direct)

        # Embed the nodes
        h_V = self.norm_nodes(self.node_embedding(torch.cat(h_V, dim=-1)))
        h_E = self.norm_edges(self.edge_embedding(torch.cat(h_E, dim=-1)))

        # prepare the variables to return
        S = torch.masked_select(S, mask_bool)
        shift = mask.sum(dim=1).cumsum(dim=0) - mask.sum(dim=1)
        src = shift.view(B,1,1) + E_idx
        src = torch.masked_select(src, mask_attend).view(1,-1)
        dst = shift.view(B,1,1) + torch.arange(0, N, device=src.device).view(1,-1,1).expand_as(mask_attend)
        dst = torch.masked_select(dst, mask_attend).view(1,-1)
        E_idx = torch.cat((dst, src), dim=0).long()

        sparse_idx = mask.nonzero()
        X = X[sparse_idx[:,0], sparse_idx[:,1], :, :]
        batch_id = sparse_idx[:,0]
        return X, S, h_V, h_E, E_idx, batch_id

class SimpleQCBM:
    """量子电路玻尔兹曼机，基于PennyLane框架实现"""
    def __init__(self, num_qubits, num_layers, total_sequences, device='cpu'):
        self.num_qubits = num_qubits
        self.num_layers = num_layers
        self.total_sequences = total_sequences  # 训练集中RNA序列的总数
        self.device = device

        # 创建PennyLane量子设备，设置shots参数用于采样
        self.qml_device = qml.device('default.qubit', wires=num_qubits, shots=self.total_sequences )

        # 参数化量子门参数 - 使用numpy数组而不是torch参数
        self.rotation_params = pnp.random.randn(num_layers, num_qubits, 3) * 0.1
        self.entangling_params = pnp.random.randn(num_layers, max(1, num_qubits-1)) * 0.1

        # 创建量子电路用于采样
        self.qcircuit_sample = qml.QNode(self._circuit_sample, self.qml_device, interface='numpy')
        # 创建量子电路用于概率计算
        #self.qcircuit_probs = qml.QNode(self._circuit_probs, self.qml_device, interface='numpy')

    def _circuit_sample(self, params):
        """定义参数化量子电路用于采样"""
        # 参数化量子电路层
        for layer in range(self.num_layers):
            # 单量子比特旋转门
            for qubit in range(self.num_qubits):
                qml.RX(params[layer, qubit, 0], wires=qubit)
                qml.RY(params[layer, qubit, 1], wires=qubit)
                qml.RZ(params[layer, qubit, 2], wires=qubit)

            # 纠缠门（CNOT）
            if self.num_qubits > 1:
                for i in range(self.num_qubits - 1):
                    qml.CNOT(wires=[i, i+1])

        # 返回采样结果
        return qml.sample(wires=range(self.num_qubits))



    def get_sample_counts(self, num_samples=None):
        """获取采样结果的counts字典"""
        if num_samples is None:
            num_samples = self.total_sequences

        # 进行采样
        raw_samples = self.qcircuit_sample(self.rotation_params)

        # 统计每个比特串的出现次数
        counts = {}
        for sample in raw_samples:
            # 将比特串转换为字符串作为key
            bit_string = ''.join([str(int(bit)) for bit in sample])
            counts[bit_string] = counts.get(bit_string, 0) + 1

        return counts

    def sample(self, num_samples=None, device=None):
        """从QCBM中采样并返回比特串格式"""
        if device is None:
            device = self.device

        # 如果没有指定采样数，使用total_sequences
        if num_samples is None:
            num_samples = self.total_sequences

        # 进行原始采样

        raw_samples = self.qcircuit_sample(self.rotation_params)

        # 转换为torch张量格式
        bit_strings = []
        for sample in raw_samples:
            bit_tensor = torch.tensor(sample, dtype=torch.float32, device=device)
            bit_strings.append(bit_tensor)

        result = torch.stack(bit_strings).to(device)
        return result

    def get_parameters(self):
        """获取所有参数的扁平化向量"""
        return pnp.concatenate([self.rotation_params.flatten(), self.entangling_params.flatten()])

    def set_parameters(self, params_vector):
        """从扁平化向量设置参数"""
        # 计算参数分割点
        rotation_size = self.num_layers * self.num_qubits * 3

        # 重新整形参数
        self.rotation_params = params_vector[:rotation_size].reshape(self.num_layers, self.num_qubits, 3)
        self.entangling_params = params_vector[rotation_size:].reshape(self.num_layers, max(1, self.num_qubits-1))

    def parameters(self):
        """返回参数列表，兼容原有接口"""
        return [self.rotation_params, self.entangling_params]

    def state_dict(self):
        """返回模型状态字典，用于保存模型"""
        return {
            'rotation_params': self.rotation_params,
            'entangling_params': self.entangling_params,
            'num_qubits': self.num_qubits,
            'num_layers': self.num_layers,
            'total_sequences': self.total_sequences
        }

    def load_state_dict(self, state_dict):
        """从状态字典加载模型参数"""
        self.rotation_params = state_dict['rotation_params']
        self.entangling_params = state_dict['entangling_params']
        self.num_qubits = state_dict['num_qubits']
        self.num_layers = state_dict['num_layers']
        self.total_sequences = state_dict['total_sequences']


# 带QCBM的RNA模型
class RNAModelWithQCBM(nn.Module):
    def __init__(self, model_config, qcbm_config, device='cpu'):
        super(RNAModelWithQCBM, self).__init__()
        self.device = device
        self.model_config = model_config
        self.qcbm_config = qcbm_config
        self.hidden_dim = model_config.hidden

        # RNA特征提取
        self.features = RNAFeatures(
            edge_features=model_config.hidden,
            node_features=model_config.hidden,
            node_feat_types=model_config.node_feat_types,
            edge_feat_types=model_config.edge_feat_types,
            top_k=model_config.k_neighbors,
            dropout=model_config.dropout
        )

        # 编码器层
        # Encoder layers
        self.encoder_layers = nn.ModuleList([
            MPNNLayer(self.hidden_dim, self.hidden_dim*2, dropout=model_config.dropout)
            for _ in range(model_config.num_encoder_layers)
        ])

        # Decoder layers
        self.decoder_layers = nn.ModuleList([
            MPNNLayer(self.hidden_dim, self.hidden_dim*2, dropout=model_config.dropout)
            for _ in range(model_config.num_decoder_layers)
        ])

        # 输出层
        self.readout = nn.Linear(model_config.hidden, model_config.vocab_size)

        # 先验样本嵌入层
        self.prior_embedding = nn.Linear(qcbm_config.num_qubits, model_config.hidden)

    def forward(self, X, S, mask, prior_samples=None, sequence_indices=None):
        # 特征提取
        X_feat, S_feat, h_V, h_E, E_idx, batch_id = self.features(X, S, mask)


        # 如果提供了先验样本，将其嵌入并添加到边特征中
        if prior_samples is not None and sequence_indices is not None:
            # 根据sequence_indices为每个碱基分配对应的先验样本
            batch_prior_samples = prior_samples[sequence_indices]

            # 确保在正确的设备上
            batch_prior_samples = batch_prior_samples.to(self.device)
            #print(batch_prior_samples)

            # 处理先验样本：[batch_size, num_qubits] -> [batch_size, hidden_dim]
            prior_embeddings = self.prior_embedding(batch_prior_samples)  # [batch_size, hidden_dim]

            # 为每条边分配对应序列的先验嵌入
            # E_idx[0]包含边的源节点索引，通过batch_id可以找到对应的序列
            edge_batch_ids = batch_id[E_idx[0]]  # 获取每条边对应的序列ID
            print('batch_id',batch_id,batch_id.shape)
            print('E_idx[0]',E_idx[0],E_idx[0].shape)
            print("edge_batch_ids",edge_batch_ids,edge_batch_ids.shape)
            edge_prior_embeddings = prior_embeddings[edge_batch_ids]  # [num_edges, hidden_dim]
            
            # 将先验嵌入添加到对应的边嵌入
            h_E = h_E + edge_prior_embeddings  # 逐元素相加

        # Encoder层
        for enc_layer in self.encoder_layers:
            h_EV = torch.cat([h_E, h_V[E_idx[0]], h_V[E_idx[1]]], dim=-1)
            h_V = enc_layer(h_V, h_EV, E_idx, batch_id)

        # Decoder层
        for dec_layer in self.decoder_layers:
            h_EV = torch.cat([h_E, h_V[E_idx[0]], h_V[E_idx[1]]], dim=-1)
            h_V = dec_layer(h_V, h_EV, E_idx, batch_id)

        # 输出
        logits = self.readout(h_V)
        log_probs = F.log_softmax(logits, dim=-1)

        return log_probs, S_feat, batch_id

    def compute_sequence_recovery_scores(self, X, S, mask, lengths, prior_samples, sequence_indices):
        """计算序列级别的recovery rate分数"""
        with torch.no_grad():
            # 获取预测结果
            log_probs, S_feat, batch_id = self.forward(X, S, mask, prior_samples, sequence_indices)

            # 获取预测的碱基
            probs = F.softmax(log_probs, dim=-1)
            predicted = probs.argmax(dim=-1)

            # 按序列计算recovery rate
            recovery_scores = []
            start_idx = 0
            for length in lengths:
                end_idx = start_idx + length.item()
                pred_seq = predicted[start_idx:end_idx]
                true_seq = S_feat[start_idx:end_idx]

                # 计算该序列的recovery rate
                correct = (pred_seq == true_seq).sum().float()
                total = length
                recovery = correct / total
                recovery_scores.append(recovery.cpu().numpy())

                start_idx = end_idx

            return recovery_scores

    def compute_package_scores(self, X, S, mask, lengths, prior_samples, sequence_indices, pdb_ids):
        """使用package_score目录下的打分机制计算先验样本分数"""
        with torch.no_grad():
            # 获取预测结果
            log_probs, S_feat, batch_id = self.forward(X, S, mask, prior_samples, sequence_indices)

            # 获取预测的碱基
            probs = F.softmax(log_probs, dim=-1)
            predicted = probs.argmax(dim=-1)

            # 碱基映射：0->A, 1->U, 2->G, 3->C
            base_map = {0: 'A', 1: 'U', 2: 'G', 3: 'C'}

            # 按序列计算package_score分数
            package_scores = []
            start_idx = 0

            for i, length in enumerate(lengths):
                end_idx = start_idx + length.item()
                pred_seq = predicted[start_idx:end_idx]
                true_seq = S_feat[start_idx:end_idx]

                # 转换为碱基序列字符串
                pred_sequence = ''.join([base_map[int(base)] for base in pred_seq.cpu().numpy()])
                true_sequence = ''.join([base_map[int(base)] for base in true_seq.cpu().numpy()])

                try:
                    # 创建临时文件保存预测序列
                    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as pred_file:
                        pred_file.write(f">predicted_seq\n{pred_sequence}\n")
                        pred_file_path = pred_file.name

                    # 创建临时文件保存参考序列
                    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as ref_file:
                        ref_file.write(f">reference_seq\n{true_sequence}\n")
                        ref_file_path = ref_file.name

                    # 使用RNASequenceScorer计算分数
                    scorer = RNASequenceScorerOptimized()

                    # 加载序列
                    pred_seq_loaded = scorer.load_fasta_sequence(pred_file_path)
                    ref_seq_loaded = scorer.load_fasta_sequence(ref_file_path)

                    # 创建简单的碱基对矩阵（对角线为1，表示自配对）
                    seq_len = len(true_sequence)
                    ref_matrix = np.eye(seq_len)

                    # 计算综合分数
                    combined_score_result = scorer.calculate_combined_score(
                        pred_seq_loaded, ref_seq_loaded, ref_matrix, lambda_param=0.5
                    )

                    # 从字典中提取combined_score值
                    package_scores.append(combined_score_result['combined_score'])

                    # 清理临时文件
                    os.unlink(pred_file_path)
                    os.unlink(ref_file_path)

                except Exception as e:
                    print(f"Error computing package score for sequence {i}: {e}")
                    # 如果出错，使用简单的序列相似度作为备选
                    similarity = (pred_seq == true_seq).sum().float() / length
                    package_scores.append(similarity.cpu().numpy())

                start_idx = end_idx

            return package_scores

# QCBM训练器
class QCBMTrainer:
    def __init__(self, qcbm_model, rna_model_with_qcbm, qcbm_config, total_sequences, device='cpu'):
        self.qcbm_model = qcbm_model
        self.rna_model_with_qcbm = rna_model_with_qcbm
        self.qcbm_config = qcbm_config
        self.total_sequences = total_sequences  # 训练集中RNA序列的总数
        self.device = device

        # 生成先验样本数据集（与训练集RNA序列一一对应）
        self.prior_samples = self._generate_prior_dataset()
        print(self.prior_samples)

        # 存储目标概率分布
        self.target_distribution = None

    def _generate_prior_dataset(self):
        """生成先验样本数据集"""
        print(f"Generating {self.total_sequences} prior samples...")
        prior_samples = self.qcbm_model.sample(self.total_sequences, device=self.device)
        print(f"Generated prior samples shape: {prior_samples.shape}")
        return prior_samples

    def generate_prior_samples(self, num_samples):
        """根据指定数量生成先验样本"""
        print(f"Generating {num_samples} prior samples with current QCBM parameters...")
        prior_samples = self.qcbm_model.sample(num_samples, device=self.device)
        print(f"Generated prior samples shape: {prior_samples.shape}")
        return prior_samples

    def train_qcbm_epoch(self, train_loader, optimizer):
        """训练一个epoch"""
        self.rna_model_with_qcbm.train()

        # 1. 使用先验数据集和训练集训练RNAModelWithQCBM
        total_loss = 0
        num_batches = 0

        # 创建序列索引映射
        sequence_index = 0

        for batch in tqdm(train_loader, desc="Training RNAModelWithQCBM"):
            # 清零梯度
            optimizer.zero_grad()

            X, S, mask, lengths, names = batch
            X = X.to(self.device)
            S = S.to(self.device)
            mask = mask.to(self.device)

            # 为当前batch的每个碱基分配对应的先验样本索引
            batch_size = sum(lengths)  # 当前batch中的总碱基数
            sequence_indices = []

            for length in lengths:
                # 为当前序列的每个碱基分配相同的先验样本索引
                seq_indices = [sequence_index % self.total_sequences] * length.item()
                sequence_indices.extend(seq_indices)
                sequence_index += 1

            sequence_indices = torch.tensor(sequence_indices, device=self.device)

            # 前向传播
            log_probs, S_feat, batch_id = self.rna_model_with_qcbm(
                X, S, mask, self.prior_samples, sequence_indices
            )

            # 计算损失
            criterion = nn.CrossEntropyLoss()
            loss = criterion(log_probs, S_feat)

            # 反向传播
            loss.backward()

            # 更新参数
            optimizer.step()

            total_loss += loss.item()
            num_batches += 1

        avg_loss = total_loss / num_batches
        print(f"RNAModelWithQCBM training loss: {avg_loss:.4f}")

        # 2. 计算先验样本的分数（序列级别的recovery rate）
        self.rna_model_with_qcbm.eval()
        prior_scores = self._compute_prior_scores(train_loader)

        # 3. 构建目标概率分布并优化QCBM参数
        self._optimize_qcbm_parameters(prior_scores)

       # 4. 重新生成先验样本集合（使用更新后的QCBM参数）
        print("Regenerating prior samples with updated QCBM parameters...")
        self.prior_samples = self._generate_prior_dataset()

        return avg_loss

    def _compute_prior_scores(self, train_loader):
        """计算先验样本的分数（序列级别recovery rate）"""
        print("Computing prior sample scores...")

        # 存储所有分数
        all_scores = []
        sequence_index = 0

        with torch.no_grad():
            for batch in tqdm(train_loader, desc="Computing scores"):
                X, S, mask, lengths, names = batch
                X = X.to(self.device)
                S = S.to(self.device)
                mask = mask.to(self.device)

                # 为当前batch的每个碱基分配对应的先验样本索引
                sequence_indices = []

                for length in lengths:
                    seq_idx = sequence_index % self.total_sequences
                    # 为当前序列的每个碱基分配相同的先验样本索引
                    seq_indices = [seq_idx] * length.item()
                    sequence_indices.extend(seq_indices)
                    sequence_index += 1

                sequence_indices = torch.tensor(sequence_indices, device=self.device)

                # 获取当前批次的pdb_ids
                batch_pdb_ids = names  # names包含了当前batch的pdb_ids

                # 计算package_score分数
                package_scores = self.rna_model_with_qcbm.compute_package_scores(
                    X, S, mask, lengths, self.prior_samples, sequence_indices, batch_pdb_ids
                )

                # 将分数添加到总列表中
                all_scores.extend(package_scores)

        # 合并所有分数，确保转换为float类型
        all_scores = torch.tensor([float(score) for score in all_scores], device=self.device)

        # 处理重复采样：对相同的先验样本取平均分数
        unique_samples = {}
        for i, sample in enumerate(self.prior_samples):
            # 将样本转换为tuple作为key
            sample_key = tuple(sample.cpu().numpy())
            if sample_key not in unique_samples:
                unique_samples[sample_key] = []
            unique_samples[sample_key].append(all_scores[i])

        # 计算唯一先验样本的平均分数并按十进制排序
        unique_sample_keys = list(unique_samples.keys())

        # 按十进制值排序
        unique_sample_keys.sort(key=lambda x: int(''.join([str(int(bit)) for bit in x]), 2))
        print(unique_sample_keys)

        unique_scores = []

        for sample_key in unique_sample_keys:
            scores_list = unique_samples[sample_key]
            avg_score = torch.stack(scores_list).mean()
            unique_scores.append(avg_score)

        # 将分数转换为张量
        unique_scores = torch.stack(unique_scores)

        # 直接对唯一先验样本的分数向量应用softmax得到目标概率分布
        target_probs = F.softmax(unique_scores, dim=0)

        print(f"Number of unique prior samples: {len(unique_sample_keys)}")
        print(f"Target probabilities shape: {target_probs.shape}")

        # 存储唯一样本的键值，用于后续提取QCBM概率
        self.unique_sample_keys = unique_sample_keys

        return target_probs

    def _optimize_qcbm_parameters(self, prior_scores):
        """使用基于非梯度的scipy.optimize优化QCBM参数以匹配先验分布"""
        print("Optimizing QCBM parameters with gradient-free optimization...")

        # prior_scores已经是归一化的目标概率分布
        if isinstance(prior_scores, torch.Tensor):
            self.target_distribution = prior_scores.detach().cpu().numpy()
        else:
            self.target_distribution = prior_scores

        # 定义目标函数
        def objective_function(params_vector):
            # 设置QCBM参数
            self.qcbm_model.set_parameters(params_vector)

            # 1. 让QCBM采样，计算counts字典
            counts = self.qcbm_model.get_sample_counts(num_samples=self.total_sequences)

            # 2. 从counts字典中提取目标概率分布中样本对应的采样次数
            pred_probs = []
            for sample_key in self.unique_sample_keys:
                # 将sample_key转换为字符串格式
                sample_str = ''.join([str(int(bit)) for bit in sample_key])
                # 使用counts.get方法，没有出现的样本采样次数为0
                count = counts.get(sample_str, 0)
                pred_probs.append(count)

            # 3. 将采样次数转换为概率分布
            pred_probs = pnp.array(pred_probs, dtype=float)
            total_samples = pnp.sum(pred_probs)
            if total_samples > 0:
                pred_probs = pred_probs / total_samples
            else:
                # 如果所有目标样本都没有被采样到，使用均匀分布
                pred_probs = pnp.ones(len(self.unique_sample_keys)) / len(self.unique_sample_keys)

            # 4. 计算交叉熵损失
            eps = 1e-8
            pred_probs = pred_probs + eps  # 避免log(0)

            if isinstance(self.target_distribution, torch.Tensor):
                target_np = self.target_distribution.cpu().numpy()
            else:
                target_np = self.target_distribution

            cross_entropy_loss = -pnp.sum(target_np * pnp.log(pred_probs))

            return cross_entropy_loss

        # 获取初始参数
        initial_params = self.qcbm_model.get_parameters()

        # 使用scipy.optimize进行优化
        print(f"Starting optimization with {len(initial_params)} parameters...")

        # 定义callback函数来显示优化过程
        iteration_count = [0]  # 使用列表来在闭包中修改值

        def optimization_callback(xk):
            iteration_count[0] += 1
            current_loss = objective_function(xk)
            print(f"Iteration {iteration_count[0]:3d}: Loss = {current_loss:.6f}")

            # 每10次迭代显示参数范围
            if iteration_count[0] % 10 == 0:
                param_min, param_max = np.min(xk), np.max(xk)
                param_mean, param_std = np.mean(xk), np.std(xk)
                print(f"  Parameters - Min: {param_min:.4f}, Max: {param_max:.4f}, Mean: {param_mean:.4f}, Std: {param_std:.4f}")

        # 使用COBYLA方法进行优化（无梯度优化）
        result = minimize(
            objective_function,
            initial_params,
            method='COBYLA',
            callback=optimization_callback,
            options={
                'maxiter': 200,  # 减少优化迭代次数
                'disp': True,
                'ftol': 1e-6
            }
        )

        # 设置最优参数
        self.qcbm_model.set_parameters(result.x)

        print(f"QCBM optimization completed. Final loss: {result.fun:.6f}")
        print(f"Optimization success: {result.success}")
        # COBYLA优化器的结果对象可能没有nit属性
        if hasattr(result, 'nit'):
            print(f"Number of iterations: {result.nit}")
        else:
            print("Number of iterations: Not available for this optimizer")

    def train_qcbm(self, train_loader, optimizer, num_epochs=30):
        """训练QCBM多个epoch"""
        for epoch in range(num_epochs):
            print(f"\n=== QCBM Epoch {epoch + 1}/{num_epochs} ===")
            avg_loss = self.train_qcbm_epoch(train_loader, optimizer)
            print(f"Epoch {epoch + 1} completed. Average loss: {avg_loss:.4f}")

def run_single_experiment(split_seed, model_seed=2025):
    """运行单次实验"""
    # 设置随机种子
    seeding(model_seed)

    # 配置
    config = Config()
    config.seed = model_seed

    # 数据划分
    train_data, valid_data, test_data = split_data(data, split_seed)
    train_data.to_csv(config.data_config.train_data_path, index=False)
    valid_data.to_csv(config.data_config.valid_data_path, index=False)
    test_data.to_csv(config.data_config.test_data_path, index=False)

    # 创建数据集
    train_dataset = RNADataset(
        data_path=config.data_config.train_data_path,
        npy_dir=config.data_config.train_npy_data_dir,
    )
    valid_dataset = RNADataset(
        data_path=config.data_config.valid_data_path,
        npy_dir=config.data_config.valid_npy_data_dir,
    )
    test_dataset = RNADataset(
        data_path=config.data_config.test_data_path,
        npy_dir=config.data_config.test_npy_data_dir,
    )

    # 创建数据加载器
    g = torch.Generator()
    g.manual_seed(config.seed)

    train_loader = DataLoader(train_dataset,
            batch_size=config.train_config.batch_size,
            shuffle=True,
            num_workers=0,
            generator=g,
            collate_fn=featurize)

    valid_loader = DataLoader(valid_dataset,
            batch_size=config.train_config.batch_size,
            shuffle=False,
            num_workers=0,
            collate_fn=featurize)

    test_loader = DataLoader(test_dataset,
            batch_size=config.train_config.batch_size,
            shuffle=False,
            num_workers=0,
            collate_fn=featurize)

    # 计算训练集中RNA序列的总数
    total_sequences = len(train_dataset)
    print(f"Total training sequences: {total_sequences}")

    # 创建模型
    qcbm_model = SimpleQCBM(
        num_qubits=config.qcbm_config.num_qubits,
        num_layers=config.qcbm_config.num_layers,
        total_sequences=total_sequences,
        device=config.device
    )

    rna_model_with_qcbm = RNAModelWithQCBM(
        model_config=config.model_config,
        qcbm_config=config.qcbm_config,
        device=config.device
    ).to(config.device)

    # 创建优化器
    optimizer = Adam(rna_model_with_qcbm.parameters(), config.train_config.lr)

    # 创建QCBM训练器
    qcbm_trainer = QCBMTrainer(
        qcbm_model=qcbm_model,
        rna_model_with_qcbm=rna_model_with_qcbm,
        qcbm_config=config.qcbm_config,
        total_sequences=total_sequences,
        device=config.device
    )

    # 创建输出目录
    if not os.path.exists(config.train_config.output_dir):
        os.makedirs(config.train_config.output_dir)

    # 训练循环
    best_valid_recovery = 0

    # 初始化记录列表
    training_losses = []
    validation_recoveries = []

    for epoch in range(config.train_config.epoch):
        print(f"\n{'='*60}")
        print(f"Epoch {epoch + 1}/{config.train_config.epoch}")
        print(f"{'='*60}")

        # 训练QCBM一个epoch
        avg_loss = qcbm_trainer.train_qcbm_epoch(train_loader, optimizer)
        print(f'Epoch {epoch + 1}/{config.train_config.epoch}, Training Loss: {avg_loss:.4f}')

        # 记录训练损失
        training_losses.append(avg_loss)

        # 验证
        rna_model_with_qcbm.eval()

        # 使用当前epoch训练好的QCBM重新采样验证集大小的先验样本
        valid_dataset_size = len(valid_dataset)
        validation_prior_samples = qcbm_trainer.generate_prior_samples(valid_dataset_size)

        with torch.no_grad():
            recovery_list = []
            sequence_index = 0

            for batch in tqdm(valid_loader, desc="Validation"):
                X, S, mask, lengths, names = batch
                X = X.to(config.device)
                S = S.to(config.device)
                mask = mask.to(config.device)

                # 为当前batch分配先验样本索引
                sequence_indices = []
                for length in lengths:
                    seq_idx = sequence_index % valid_dataset_size
                    seq_indices = [seq_idx] * length.item()
                    sequence_indices.extend(seq_indices)
                    sequence_index += 1

                sequence_indices = torch.tensor(sequence_indices, device=config.device)

                # 计算recovery rate
                recovery_scores = rna_model_with_qcbm.compute_sequence_recovery_scores(
                    X, S, mask, lengths, validation_prior_samples, sequence_indices
                )
                recovery_list.extend(recovery_scores)

            valid_recovery = np.mean(recovery_list)
            print(f'Epoch {epoch + 1}/{config.train_config.epoch}, Validation Recovery: {valid_recovery:.4f}')

            # 记录验证恢复率
            validation_recoveries.append(valid_recovery)

            # 保存最佳模型
            if valid_recovery > best_valid_recovery:
                best_valid_recovery = valid_recovery
                torch.save({
                    'rna_model': rna_model_with_qcbm.state_dict(),
                    'qcbm_model': qcbm_model.state_dict(),
                    'epoch': epoch,
                    'valid_recovery': valid_recovery
                }, config.train_config.ckpt_path)
                print(f"New best model saved with validation recovery: {valid_recovery:.4f}")

    # 测试
    print("\n" + "="*60)
    print("Testing...")
    print("="*60)

    # 加载最佳模型
    checkpoint = torch.load(config.train_config.ckpt_path, map_location=config.device)
    rna_model_with_qcbm.load_state_dict(checkpoint['rna_model'])
    qcbm_model.load_state_dict(checkpoint['qcbm_model'])

    rna_model_with_qcbm.eval()

    # 使用最后一个epoch的QCBM重新采样测试集大小的先验样本
    test_dataset_size = len(test_dataset)
    test_prior_samples = qcbm_trainer.generate_prior_samples(test_dataset_size)

    with torch.no_grad():
        recovery_list = []
        sequence_index = 0

        for batch in tqdm(test_loader, desc="Testing"):
            X, S, mask, lengths, names = batch
            X = X.to(config.device)
            S = S.to(config.device)
            mask = mask.to(config.device)

            # 为当前batch分配先验样本索引
            sequence_indices = []
            for length in lengths:
                seq_idx = sequence_index % test_dataset_size
                seq_indices = [seq_idx] * length.item()
                sequence_indices.extend(seq_indices)
                sequence_index += 1

            sequence_indices = torch.tensor(sequence_indices, device=config.device)

            # 计算recovery rate
            recovery_scores = rna_model_with_qcbm.compute_sequence_recovery_scores(
                X, S, mask, lengths, test_prior_samples, sequence_indices
            )
            recovery_list.extend(recovery_scores)

        test_recovery = np.mean(recovery_list)
        print(f'Test Recovery: {test_recovery:.4f}')

    # 保存训练记录到文件
    training_history = {
        'epoch': list(range(1, len(training_losses) + 1)),
        'training_loss': training_losses,
        'validation_recovery': validation_recoveries
    }

    # 保存为CSV文件
    history_df = pd.DataFrame(training_history)
    history_csv_path = os.path.join(config.train_config.output_dir, f'training_history_seed_{split_seed}.csv')
    history_df.to_csv(history_csv_path, index=False)
    print(f'Training history saved to: {history_csv_path}')

    # 保存为numpy文件
    #history_npy_path = os.path.join(config.train_config.output_dir, f'training_history_seed_{split_seed}.npy')
    #np.save(history_npy_path, training_history)
    #print(f'Training history saved to: {history_npy_path}')

    return {
        'split_seed': split_seed,
        'model_seed': model_seed,
        'best_val_recovery': best_valid_recovery,
        'test_recovery': test_recovery
    }

if __name__ == "__main__":
    print("Starting Q9_1 experiments with prior sample dataset and sequence-level recovery rate...")
    num_experiments = 5
    split_seeds = [42,123, 456, 789, 1011]
    #split_seeds = [42]
    model_init_seed = 2025

    results = []
    # 初始化空的DataFrame用于追加结果
    results_df = pd.DataFrame(columns=['split_seed', 'best_val_recovery', 'test_recovery'])

    for i, split_seed in enumerate(split_seeds):
        print(f"\n{'='*60}")
        print(f"Running Experiment {i+1}/{num_experiments} with split_seed={split_seed}")
        print(f"{'='*60}")

        result = run_single_experiment(split_seed, model_init_seed)
        results.append(result)

        # 将当前实验结果追加到DataFrame中
        new_row = pd.DataFrame({
            'split_seed': [result['split_seed']],
            'best_val_recovery': [result['best_val_recovery']],
            'test_recovery': [result['test_recovery']]
        })
        results_df = pd.concat([results_df, new_row], ignore_index=True)

        # 保存当前的results_df到文件（追加模式）
        results_df.to_csv('qcbm_rna_results_q9_1.csv', index=False)
        print(f"Results updated and saved to qcbm_rna_results_q8_10.csv")

        print(f"\nExperiment {i+1} Results:")
        print(f"Best Validation Recovery: {result['best_val_recovery']:.4f}")
        print(f"Test Recovery: {result['test_recovery']:.4f}")

    # 计算统计结果
    val_recoveries = [r['best_val_recovery'] for r in results]
    test_recoveries = [r['test_recovery'] for r in results]

    print(f"\n{'='*60}")
    print("FINAL RESULTS SUMMARY")
    print(f"{'='*60}")
    print(f"Validation Recovery - Mean: {np.mean(val_recoveries):.4f}, Std: {np.std(val_recoveries):.4f}")
    print(f"Test Recovery - Mean: {np.mean(test_recoveries):.4f}, Std: {np.std(test_recoveries):.4f}")

    results_df.to_csv('qcbm_rna_results_q9_1.csv', index=False)
    print(f"\nResults saved to qcbm_rna_results_q9_1.csv")
    print(f"\nAll Q8 experiments completed!")
