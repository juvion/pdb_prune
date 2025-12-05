import os
from dataset_src.RNA_graph_generate_ss import RNA_imem, dataset_argument
from torch.optim import Adam
from torch_geometric.data import Batch, Data
from dataset_src.utils import NormalizeRNA, get_stat
from Bio.PDB import PDBParser
import torch.nn.functional as F
import torch
from tqdm import tqdm
import traceback
import argparse
import numpy as np
# Compatibility shim: NumPy >=1.20 removed np.int; Biopython SASA still uses it
if not hasattr(np, "int"):
    np.int = int
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import logging

nucleotides_type = ['A', 'U', 'G', 'C']

logger = logging.getLogger("generate_graph_ss")

def quick_scan_atom_names(pdb_path, max_lines=100000):
    """Fast text scan of key atom names to aid debugging when graph generation returns None.
    Returns a dict with counts of occurrences for important atom name variants.
    """
    keys = {"C4'": 0, "C4*": 0, "C1'": 0, "C1*": 0, "N9": 0, "N1": 0}
    try:
        with open(pdb_path, 'r') as f:
            for i, line in enumerate(f):
                if i >= max_lines:
                    break
                # Check typical record lines
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue
                for k in keys:
                    if k in line:
                        keys[k] += 1
        return keys
    except Exception as e:
        return {"error": str(e)}

def get_struc2ndRes(pdb_filename):
    integer_encoded = []
    struc_2nds_res_alphabet = ['S', 'M', 'I', 'B', 'H', 'K', 'E', 'X']
    char_to_int = dict((c, i) for i, c in enumerate(struc_2nds_res_alphabet))
    p = PDBParser()
    structure = p.get_structure('random_id', pdb_filename)
    model = structure[0]
    model_residues = [(chain.id, residue.id[1]) for chain in model for residue in chain if residue.id[0] == ' ']
    one_hot_list = torch.zeros(len(model_residues), len(struc_2nds_res_alphabet))
    # st_file_name = pdb_filename.replace(".pdb", ".st").replace(f'{key}', f"{key}_ss")
    st_file_name = pdb_filename.replace(".pdb", ".st")

    st_file = st_file_name
    if not os.path.exists(st_file):
        return one_hot_list 
    else:
        with open(st_file, 'r') as f:
            counter = 0
            for line in f:
                if counter != 5:
                    counter += 1
                    continue
                else:
                    counter += 1
                    ss = line.strip()
                    current_position = 0
                    for cha in ss:
                        integer_encoded.append(char_to_int[cha])
                        one_hot = F.one_hot(torch.tensor(integer_encoded[-1]), num_classes=8)
                        one_hot_list[current_position] = one_hot
                        current_position += 1
        return one_hot_list

def prepare_graph(data):
    del data['distances']
    del data['edge_dist']
    mu_r_norm = data.mu_r_norm

    extra_x_feature = torch.cat([data.x[:, 4:], mu_r_norm], dim=1)
    graph = Data(
        x=data.x[:, :4],
        extra_x=extra_x_feature,
        pos=data.pos,
        edge_index=data.edge_index,
        edge_attr=data.edge_attr,
        ss=data.ss[:data.x.shape[0], :],
        sasa=data.x[:, 4]
    )
    return graph

def pdb2graph(filename, normalize_path='dataset_src/mean_attr.pt', if_transform=False):
    dataset_arg = dataset_argument(n=666)
    dataset = RNA_imem(dataset_arg['root'], dataset_arg['name'], split='test',
                       divide_num=dataset_arg['divide_num'], divide_idx=dataset_arg['divide_idx'],
                       c_4prime_max_neighbors=dataset_arg['c_4prime_max_neighbors'],
                       set_length=dataset_arg['set_length'],
                       struc_2nds_res_path=dataset_arg['struc_2nds_res_path'],
                       random_sampling=True, diffusion=True)
    rec_info = dataset.get_receptor_inference(filename)
    if rec_info is None:
        logger.debug(f"get_receptor_inference returned None for {filename}")
        return None
    rec, rec_coords, c_4prime_coords, p_coords, n_coords= rec_info
    struc_2nd_res = get_struc2ndRes(filename)

    rec_graph = dataset.get_c4prime_graph(rec, c_4prime_coords, p_coords, n_coords, rec_coords, struc_2nd_res)

    if rec_graph:
        if if_transform:
            normalize_transform = NormalizeRNA(filename=normalize_path)
            graph = normalize_transform(rec_graph)
        else:
            graph = rec_graph
        return graph
    else:
        return None

def get_graph(
    filename,
    pdb_dir,
    save_dir,
    all_dir,
    exclude_ls=None,
    normalize_path='dataset_src/mean_attr.pt',
    use_transform=False,
):
    exclude_ls = exclude_ls or []
    os.makedirs(save_dir, exist_ok=True)
    os.makedirs(all_dir, exist_ok=True)

    out_path = os.path.join(save_dir, filename.replace('.pdb', '.pt'))
    if os.path.exists(out_path):
        return True
    try:
        if filename in exclude_ls:
            print(filename + "  excluded")
            return False
        graph = pdb2graph(
            os.path.join(pdb_dir, filename),
            normalize_path=normalize_path,
            if_transform=use_transform,
        )
        if graph:
            torch.save(graph, out_path)
            torch.save(graph, os.path.join(all_dir, filename.replace('.pdb', '.pt')))
            return True
        else:
            print("err")
            return False
    except (IndexError, KeyError, ValueError, torch.serialization.pickle.UnpicklingError) as e:
        print(f"err: {filename} - {e}")
        return False

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='将PDB转换为图数据（支持二级结构特征）')
    parser.add_argument('--pdb_dir', type=str, default=None,
                        help='输入PDB文件夹路径；若未提供，则使用默认dataset_src下的各split')
    parser.add_argument('--out_dir', type=str, default='graph_dataset',
                        help='输出根目录，默认graph_dataset')
    parser.add_argument('--dataset_id', type=str, default='dataset_0.8',
                        help='输出数据集子目录名，例如dataset_0.8')
    parser.add_argument('--split', type=str, default=None,
                        help='自定义split名称（仅在指定--pdb_dir时有效），例如demo')
    parser.add_argument('--use_transform', action='store_true',
                        help='对图特征执行NormalizeRNA变换（需要mean_attr.pt）')
    parser.add_argument('--normalize_path', type=str, default='dataset_src/mean_attr.pt',
                        help='NormalizeRNA所用均值文件路径')
    parser.add_argument('--verbose', action='store_true',
                        help='打印详细调试信息（逐文件状态、失败原因扫描）')
    args = parser.parse_args()

    # Configure logging
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='[%(asctime)s] %(levelname)s %(message)s'
    )

    logger.info(
        f"Args: pdb_dir={args.pdb_dir}, out_dir={args.out_dir}, dataset_id={args.dataset_id}, "
        f"split={args.split}, use_transform={args.use_transform}, normalize_path={args.normalize_path}"
    )
    if args.use_transform:
        logger.info(f"Normalize file exists: {os.path.exists(args.normalize_path)} -> {args.normalize_path}")

    out_dir = args.out_dir
    data_set_id = args.dataset_id
    all_dir = os.path.join(out_dir, data_set_id, 'all')
    os.makedirs(all_dir, exist_ok=True)

    def process_split(pdb_dir: str, split_name: str):
        save_dir = os.path.join(out_dir, data_set_id, split_name)
        os.makedirs(save_dir, exist_ok=True)
        filename_list = sorted([i for i in os.listdir(pdb_dir) if i.endswith('.pdb')])
        print(f"[generate_graph_ss] Split '{split_name}': 在 '{pdb_dir}' 找到 {len(filename_list)} 个PDB文件，输出目录 '{save_dir}'。")
        logger.info(f"Split '{split_name}': num_files={len(filename_list)}")
        if args.verbose:
            logger.debug(f"Sample files: {filename_list[:5]}")

        saved_count = 0
        for filename in tqdm(filename_list, desc=f"Processing {split_name}"):
            out_path = os.path.join(save_dir, filename.replace('.pdb', '.pt'))
            if os.path.exists(out_path):
                if args.verbose:
                    logger.debug(f"Skip existing: {filename}")
                continue
            try:
                graph = pdb2graph(
                    os.path.join(pdb_dir, filename),
                    normalize_path=args.normalize_path,
                    if_transform=args.use_transform
                )
                if graph:
                    torch.save(graph, out_path)
                    torch.save(graph, os.path.join(all_dir, filename.replace('.pdb', '.pt')))
                    saved_count += 1
                    if args.verbose:
                        logger.info(f"Saved: {filename} -> {out_path}")
                else:
                    # Provide quick diagnostics for failure
                    if args.verbose:
                        pdb_path = os.path.join(pdb_dir, filename)
                        scan = quick_scan_atom_names(pdb_path)
                        st_exists = os.path.exists(pdb_path.replace('.pdb', '.st'))
                        logger.warning(
                            f"Graph None for {filename} | atom_scan={scan} | st_exists={st_exists}"
                        )
            except Exception:
                print(f"[generate_graph_ss] 处理 {filename} 失败：")
                traceback.print_exc()
                logger.exception(f"Exception while processing {filename}")
        print(f"[generate_graph_ss] Split '{split_name}' 完成：成功保存 {saved_count} 个图文件，失败 {len(filename_list)-saved_count} 个。")
        logger.info(
            f"Split '{split_name}' finished: saved={saved_count}, failed={len(filename_list)-saved_count}"
        )

    # 执行流程：单目录或默认数据集结构
    if args.pdb_dir:
        if not os.path.isdir(args.pdb_dir):
            raise FileNotFoundError(f"--pdb_dir 不存在：{args.pdb_dir}")
        split_name = args.split or 'custom'
        process_split(args.pdb_dir, split_name)
    else:
        for key in ['test_0.8', 'validation_0.8', 'train_0.8']:
            pdb_dir = os.path.join('dataset_src', key)
            process_split(pdb_dir, key)
