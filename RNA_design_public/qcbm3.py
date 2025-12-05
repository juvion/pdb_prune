import numpy as np
import pennylane as qml
from pennylane import numpy as pnp
import torch
import torch.nn.functional as F
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from tqdm import tqdm
import random

# 设置随机种子
np.random.seed(42)
random.seed(42)
torch.manual_seed(42)

class QCBM:
    """量子电路玻尔兹曼机"""
    
    def __init__(self, num_qubits=9, num_layers=3, device_name='default.qubit', fixed_shots=4096):
        self.num_qubits = num_qubits
        self.num_layers = num_layers
        self.device_name = device_name
        self.fixed_shots = int(fixed_shots)
        # 预创建：无shots的概率设备 + 固定shots的采样设备
        self.device_probs = qml.device(device_name, wires=num_qubits, shots=None)
        self.device_sample = qml.device(device_name, wires=num_qubits, shots=self.fixed_shots)
        
        # 参数数量：每层每比特包含4个旋转参数 + 每层包含 C(num_qubits,2) 个 IsingXX 纠缠角度参数（全连接）
        self.num_params = self.num_layers * self.num_qubits * 4 + self.num_layers * (self.num_qubits * (self.num_qubits - 1) // 2)
        
        # 初始化参数
        self.params = pnp.random.uniform(0, 2*pnp.pi, self.num_params, requires_grad=True)
        
        # 定义概率分布电路（绑定无shots设备，返回qml.probs）
        @qml.qnode(self.device_probs, interface='autograd')
        def circuit(params):
            wires = list(range(self.num_qubits))
            rot_size = self.num_layers * self.num_qubits * 4
            rotation_weights = params[:rot_size].reshape(self.num_layers, self.num_qubits, 4)
            pair_count = self.num_qubits * (self.num_qubits - 1) // 2
            entangler_weights = params[rot_size:].reshape(self.num_layers, pair_count)
            for l in range(self.num_layers):
                for i in range(self.num_qubits):
                    qml.RX(rotation_weights[l, i, 0], wires=i)
                    qml.RZ(rotation_weights[l, i, 1], wires=i)
                k = 0
                for a in range(self.num_qubits):
                    for b in range(a + 1, self.num_qubits):
                        qml.IsingXX(entangler_weights[l, k], wires=[a, b])
                        k += 1
                for i in range(self.num_qubits):
                    qml.RZ(rotation_weights[l, i, 2], wires=i)
                    qml.RX(rotation_weights[l, i, 3], wires=i)
            return qml.probs(wires=wires)
        self.circuit = circuit
        
        # 采样电路（返回qml.sample），与固定shots设备绑定的QNode
        def _circuit_sample(params):
            wires = list(range(self.num_qubits))
            rot_size = self.num_layers * self.num_qubits * 4
            rotation_weights = params[:rot_size].reshape(self.num_layers, self.num_qubits, 4)
            pair_count = self.num_qubits * (self.num_qubits - 1) // 2
            entangler_weights = params[rot_size:].reshape(self.num_layers, pair_count)
            for l in range(self.num_layers):
                for i in range(self.num_qubits):
                    qml.RX(rotation_weights[l, i, 0], wires=i)
                    qml.RZ(rotation_weights[l, i, 1], wires=i)
                k = 0
                for a in range(self.num_qubits):
                    for b in range(a + 1, self.num_qubits):
                        qml.IsingXX(entangler_weights[l, k], wires=[a, b])
                        k += 1
                for i in range(self.num_qubits):
                    qml.RZ(rotation_weights[l, i, 2], wires=i)
                    qml.RX(rotation_weights[l, i, 3], wires=i)
            return qml.sample(wires=wires)
        self._circuit_sample = _circuit_sample
        self.qnode_sample = qml.QNode(self._circuit_sample, self.device_sample, interface='numpy')
    
    def get_probabilities(self, params=None):
        if params is None:
            params = self.params
        return self.circuit(params)
    
    def sample(self, num_samples=1000, params=None):
        if params is None:
            params = self.params
        if int(num_samples) == int(self.fixed_shots):
            raw_samples = self.qnode_sample(params)
        else:
            # 新建设备与QNode，按需设置shots=num_samples，使用采样电路
            temp_device = qml.device(self.device_name, wires=self.num_qubits, shots=int(num_samples))
            temp_qnode = qml.QNode(self._circuit_sample, temp_device, interface='numpy')
            raw_samples = temp_qnode(params)
        binary_samples = [[int(bit) for bit in sample] for sample in raw_samples]
        return np.array(binary_samples)
    
    def get_sample_counts(self, num_samples=1000, params=None):
        if params is None:
            params = self.params
        if int(num_samples) == int(self.fixed_shots):
            raw = self.qnode_sample(params)
        else:
            # 按需创建临时设备与QNode以匹配可变shots
            temp_device = qml.device(self.device_name, wires=self.num_qubits, shots=int(num_samples))
            temp_qnode = qml.QNode(self._circuit_sample, temp_device, interface='numpy')
            raw = temp_qnode(params)
        counts = {}
        for sample in raw:
            s = ''.join(str(int(b)) for b in sample)
            counts[s] = counts.get(s, 0) + 1
        return counts
    
    def set_parameters(self, params):
        self.params = pnp.array(params, requires_grad=True)

def kl_divergence(p, q, eps=1e-12):
    """计算KL散度"""
    # 添加小的epsilon避免log(0)
    q_safe = pnp.clip(q, eps, 1.0)
    q_safe = q_safe / pnp.sum(q_safe)
    p_safe = pnp.clip(p, eps, 1.0)
    return pnp.sum(p_safe * pnp.log(p_safe / q_safe))


def cross_entropy_loss(p, q, eps=1e-12):
    """计算交叉熵损失"""
    q_safe = pnp.clip(q, eps, 1.0)
    q_safe = q_safe / pnp.sum(q_safe)
    return -pnp.sum(p * pnp.log(q_safe))

def train_qcbm(target_probs, num_qubits=9, num_layers=3, max_iterations=1000, learning_rate=0.01):
    """训练QCBM以拟合目标概率分布"""
    
    print(f"开始训练QCBM: {num_qubits} qubits, {num_layers} layers")
    print(f"目标分布维度: {len(target_probs)}")
    print(f"预期维度: {2**num_qubits}")
    
    # 确保目标概率分布归一化
    target_probs = np.array(target_probs)
    target_probs = target_probs / np.sum(target_probs)
    target_probs_pnp = pnp.array(target_probs)
    
    # 初始化QCBM
    qcbm = QCBM(num_qubits=num_qubits, num_layers=num_layers)
    
    # 记录训练历史
    loss_history = []
    kl_history = []
    
    # 定义目标函数
    def objective_function(params):
        try:
            # 获取当前概率分布
            params_pnp = pnp.array(params, requires_grad=True)
            current_probs = qcbm.get_probabilities(params_pnp)
            
            # 计算交叉熵损失（使用pennylane.numpy以便解析梯度）
            loss = kl_divergence(target_probs_pnp, current_probs)
            
            return loss
        except Exception as e:
            print(f"目标函数计算错误: {e}")
            return 1e6
    
    # 解析梯度函数（autograd）
    grad_fn = qml.grad(objective_function)
    
    # SciPy 兼容的 Jacobian 封装：确保输入为 autograd 数组，输出为 numpy.ndarray
    def jac_wrapper(x):
        x_pnp = pnp.array(x, requires_grad=True)
        g = grad_fn(x_pnp)
        return np.array(g, dtype=float)
    
    # 使用scipy优化
    print("开始优化...")
    
    # 尝试多种优化方法
    #methods = ['L-BFGS-B', 'SLSQP', 'TNC']
    methods = ['SLSQP']
    best_result = None
    best_loss = float('inf')
    
    for method in methods:
        print(f"\n尝试优化方法: {method}")
        try:
            iter_counter = {'count': 0}
            def optimization_callback(xk):
                iter_counter['count'] += 1
                try:
                    current_probs = qcbm.get_probabilities(pnp.array(xk))
                    current_loss = cross_entropy_loss(target_probs_pnp, current_probs)
                    current_kl = kl_divergence(target_probs_pnp, current_probs)
                    print(f"[{method}] Iter {iter_counter['count']}: loss={float(current_loss):.6f}, KL={float(current_kl):.6f}")
                except Exception as e_cb:
                    print(f"[{method}] Iter {iter_counter['count']}: 回调计算错误: {e_cb}")
            result = minimize(
                objective_function,
                np.array(qcbm.params, dtype=float),
                method=method,
                jac=jac_wrapper,
                callback=optimization_callback,
                options={
                    'maxiter': max_iterations,
                    'disp': True,
                    'ftol': 1e-9,
                    'gtol': 1e-9
                }
            )
            
            if result.fun < best_loss:
                best_loss = result.fun
                best_result = result
                print(f"找到更好的结果，损失: {best_loss:.6f}")
                
        except Exception as e:
            print(f"优化方法 {method} 失败: {e}")
            continue
    
    if best_result is None:
        print("所有优化方法都失败了")
        return qcbm, [], []
    
    # 设置最佳参数
    qcbm.set_parameters(best_result.x)
    
    # 计算最终结果
    final_probs = qcbm.get_probabilities()
    final_loss = cross_entropy_loss(target_probs_pnp, final_probs)
    final_kl = kl_divergence(target_probs_pnp, final_probs)
    
    print(f"\n训练完成!")
    print(f"最终交叉熵损失: {float(final_loss):.6f}")
    print(f"最终KL散度: {float(final_kl):.6f}")
    print(f"优化成功: {best_result.success}")
    print(f"优化信息: {best_result.message}")
    
    return qcbm, [final_loss], [final_kl]

def plot_comparison(target_probs, learned_probs, save_path='qcbm_comparison.png'):
    """绘制目标分布和学习分布的比较图"""
    plt.figure(figsize=(15, 5))
    
    # 子图1：概率分布比较
    plt.subplot(1, 3, 1)
    x = np.arange(len(target_probs))
    plt.plot(x, target_probs, 'b-', label='目标分布', alpha=0.7)
    plt.plot(x, learned_probs, 'r--', label='学习分布', alpha=0.7)
    plt.xlabel('状态索引')
    plt.ylabel('概率')
    plt.title('概率分布比较')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 子图2：散点图
    plt.subplot(1, 3, 2)
    plt.scatter(target_probs, learned_probs, alpha=0.6)
    plt.plot([0, max(target_probs)], [0, max(target_probs)], 'r--', alpha=0.8)
    plt.xlabel('目标概率')
    plt.ylabel('学习概率')
    plt.title('目标 vs 学习概率')
    plt.grid(True, alpha=0.3)
    
    # 子图3：误差分布
    plt.subplot(1, 3, 3)
    errors = np.abs(target_probs - learned_probs)
    plt.hist(errors, bins=50, alpha=0.7, edgecolor='black')
    plt.xlabel('绝对误差')
    plt.ylabel('频次')
    plt.title('误差分布')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()
    
    # 计算统计指标
    mse = np.mean((target_probs - learned_probs)**2)
    mae = np.mean(np.abs(target_probs - learned_probs))
    max_error = np.max(np.abs(target_probs - learned_probs))
    
    print(f"\n统计指标:")
    print(f"均方误差 (MSE): {mse:.8f}")
    print(f"平均绝对误差 (MAE): {mae:.8f}")
    print(f"最大绝对误差: {max_error:.8f}")

def main():
    """主函数"""
    # 用户提供的目标概率分布
    target_distribution = [0.00237882, 0.00300125, 0.00158726, 0.00183445, 0.00162958, 0.00139402, 0.00121911, 0.00076471, 0.00330065, 0.00085499, 0.00104196, 0.00106246, 0.00374068, 0.00054028, 0.00094598, 0.00159915, 0.00119503, 0.00119279, 0.00437206, 0.00058078, 0.00251297, 0.00265011, 0.00242964, 0.00305757, 0.00350357, 0.00157878, 0.00084112, 0.00341341, 0.00020158, 0.00231721, 0.00280185, 0.00143969, 0.00044163, 0.00179225, 0.00039262, 0.00181754, 0.00126620, 0.00323503, 0.00169229, 0.00106724, 0.00283863, 0.00312701, 0.00143899, 0.00175013, 0.00299158, 0.00187167, 0.00057466, 0.00066291, 0.00120290, 0.00131333, 0.00123480, 0.00174579, 0.00245570, 0.00076471, 0.00417570, 0.00154486, 0.00186010, 0.00171818, 0.00146771, 0.00073691, 0.00239104, 0.00267900, 0.00061233, 0.00232188, 0.00497059, 0.00155804, 0.00184359, 0.00113340, 0.00168320, 0.00184043, 0.00175189, 0.00070357, 0.00164137, 0.00232384, 0.00148900, 0.00465889, 0.00036555, 0.00065285, 0.00206393, 0.00079252, 0.00125933, 0.00185366, 0.00190202, 0.00076471, 0.00167842, 0.00124053, 0.00265854, 0.00253917, 0.00122626, 0.00278107, 0.00135927, 0.00199678, 0.00065369, 0.00218481, 0.00298945, 0.00395913, 0.00183940, 0.00326111, 0.00143632, 0.00316329, 0.00245570, 0.00472600, 0.00117720, 0.00129182, 0.00099914, 0.00150753, 0.00148946, 0.00076471, 0.00188545, 0.00494191, 0.00122626, 0.00194806, 0.00037659, 0.00299846, 0.00260873, 0.00035812, 0.00114034, 0.00158120, 0.00229310, 0.00247848, 0.00216684, 0.00221112, 0.00036953, 0.00171264, 0.00304178, 0.00360572, 0.00227979, 0.00515615, 0.00140382, 0.00190094, 0.00091719, 0.00170122, 0.00378230, 0.00511984, 0.00232656, 0.00115654, 0.00128436, 0.00242847, 0.00142626, 0.00321433, 0.00225866, 0.00076471, 0.00507477, 0.00253808, 0.00229827, 0.00186010, 0.00195648, 0.00061233, 0.00129300, 0.00166020, 0.00087509, 0.00069279, 0.00076471, 0.00154344, 0.00076471, 0.00152717, 0.00176827, 0.00053221, 0.00289834, 0.00329856, 0.00144062, 0.00069678, 0.00150078, 0.00186010, 0.00223947, 0.00077874, 0.00106724, 0.00161889, 0.00194847, 0.00481832, 0.00144196, 0.00092333, 0.00497506, 0.00273520, 0.00073560, 0.00180975, 0.00323144, 0.00177238, 0.00223317, 0.00098163, 0.00326676, 0.00135415, 0.00106919, 0.00192295, 0.00221277, 0.00054794, 0.00310072, 0.00020158, 0.00479499, 0.00056217, 0.00061233, 0.00091060, 0.00196637, 0.00025577, 0.00162322, 0.00100097, 0.00106724, 0.00215538, 0.00033663, 0.00193197, 0.00058048, 0.00070357, 0.00088682, 0.00194125, 0.00331024, 0.00195163, 0.00192125, 0.00215765, 0.00234662, 0.00033934, 0.00044578, 0.00228676, 0.00298260, 0.00291334, 0.00498681, 0.00200261, 0.00091719, 0.00460877, 0.00106724, 0.00300668, 0.00128436, 0.00032452, 0.00022706, 0.00028132, 0.00271181, 0.00426508, 0.00186010, 0.00044392, 0.00273046, 0.00198869, 0.00129133, 0.00212390, 0.00417214, 0.00178607, 0.00181810, 0.00367841, 0.00142867, 0.00182078, 0.00282223, 0.00146291, 0.00092517, 0.00101960, 0.00136424, 0.00066926, 0.00067742, 0.00245083, 0.00100788, 0.00487488, 0.00193029, 0.00206701, 0.00086653, 0.00095915, 0.00227038, 0.00248738, 0.00173716, 0.00187495, 0.00259198, 0.00123318, 0.00379474, 0.00136852, 0.00087361, 0.00191280, 0.00138739, 0.00347513, 0.00265081, 0.00208283, 0.00244963, 0.00032452, 0.00052246, 0.00356276, 0.00499482, 0.00409274, 0.00475590, 0.00137907, 0.00020158, 0.00079099, 0.00084514, 0.00257206, 0.00175873, 0.00076471, 0.00080840, 0.00092143, 0.00085602, 0.00083947, 0.00245570, 0.00226063, 0.00227947, 0.00414991, 0.00447960, 0.00071706, 0.00074465, 0.00129443, 0.00052246, 0.00158919, 0.00261215, 0.00252584, 0.00354632, 0.00390155, 0.00331684, 0.00324199, 0.00438154, 0.00072801, 0.00148946, 0.00168025, 0.00028630, 0.00075679, 0.00213809, 0.00070357, 0.00183157, 0.00230322, 0.00100557, 0.00164327, 0.00248163, 0.00165057, 0.00228276, 0.00106724, 0.00399380, 0.00270148, 0.00106724, 0.00190570, 0.00312516, 0.00295952, 0.00335833, 0.00302379, 0.00245570, 0.00308232, 0.00222367, 0.00399051, 0.00076471, 0.00130684, 0.00201157, 0.00450401, 0.00129803, 0.00061233, 0.00046382, 0.00258606, 0.00377919, 0.00057125, 0.00258111, 0.00415070, 0.00042280, 0.00124184, 0.00119266, 0.00169805, 0.00454731, 0.00365491, 0.00391879, 0.00244158, 0.00132834, 0.00373565, 0.00300949, 0.00067742, 0.00061233, 0.00096089, 0.00036953, 0.00178720, 0.00506303, 0.00336862, 0.00046382, 0.00213725, 0.00310612, 0.00266911, 0.00261140, 0.00052246, 0.00109081, 0.00402253, 0.00301452, 0.00161889, 0.00233822, 0.00256156, 0.00429740, 0.00281320, 0.00467648, 0.00346582, 0.00344303, 0.00203260, 0.00053730, 0.00033663, 0.00334671, 0.00088682, 0.00181095, 0.00094840, 0.00236167, 0.00316354, 0.00498896, 0.00177071, 0.00052246, 0.00174648, 0.00221697, 0.00160402, 0.00396240, 0.00042280, 0.00243143, 0.00188932, 0.00256018, 0.00068602, 0.00186010, 0.00093882, 0.00389375, 0.00089955, 0.00140896, 0.00186010, 0.00486757, 0.00076471, 0.00061233, 0.00086192, 0.00241620, 0.00234775, 0.00259370, 0.00337524, 0.00325409, 0.00287559, 0.00052246, 0.00174642, 0.00064788, 0.00108467, 0.00066291, 0.00199387, 0.00301454, 0.00157748, 0.00134924, 0.00081044, 0.00216973, 0.00340753, 0.00120299, 0.00233593, 0.00163806, 0.00039262, 0.00106724, 0.00128419, 0.00194937, 0.00126894, 0.00107704, 0.00081031, 0.00124184, 0.00148946, 0.00087460, 0.00170261, 0.00152217, 0.00318764, 0.00223278, 0.00186010, 0.00442082, 0.00053730, 0.00195179, 0.00126782, 0.00136474, 0.00191390, 0.00290106, 0.00106724, 0.00168860, 0.00159858, 0.00054794, 0.00298842, 0.00147407, 0.00134221, 0.00098327, 0.00107949, 0.00514830, 0.00114453, 0.00156043, 0.00088682, 0.00186010, 0.00311494, 0.00475381, 0.00514277, 0.00521550, 0.00082030, 0.00226105, 0.00148946, 0.00046382, 0.00148946, 0.00049816, 0.00082030, 0.00091719, 0.00196477, 0.00060909, 0.00191691, 0.00199387, 0.00233822, 0.00258845, 0.00261893, 0.00061233, 0.00176745, 0.00186392, 0.00066291, 0.00428542, 0.00067767, 0.00129685, 0.00065189, 0.00067742, 0.00462305, 0.00109323, 0.00061233, 0.00173730, 0.00398145, 0.00209755, 0.00082403, 0.00163931, 0.00188612, 0.00197341, 0.00281030, 0.00185622, 0.00203859, 0.00276975, 0.00080840, 0.00166680, 0.00061233, 0.00240694, 0.00067742, 0.00112551, 0.00347162]
    
    print("QCBM训练脚本")
    print("=" * 50)
    
    # 验证目标分布
    print(f"目标分布长度: {len(target_distribution)}")
    print(f"目标分布总和: {sum(target_distribution):.6f}")
    print(f"目标分布最小值: {min(target_distribution):.8f}")
    print(f"目标分布最大值: {max(target_distribution):.8f}")
    
    # 训练QCBM
    qcbm, loss_history, kl_history = train_qcbm(
        target_distribution,
        num_qubits=9,
        num_layers=6,
        max_iterations=200,
        learning_rate=0.01
    )
    
    # 获取学习到的分布
    learned_probs = qcbm.get_probabilities()
    learned_probs = np.array(learned_probs)
    # 使用 QCBM 电路直接采样，shots=1024，并统计每个样本（状态索引）的采样次数
    shots = 10240
    samples = qcbm.sample(num_samples=shots)
    print(samples[0:10])
    indices = [int(''.join(str(int(bit)) for bit in row), 2) for row in samples]
    sample_counts = np.bincount(indices, minlength=2**qcbm.num_qubits)
    print(f"采样 {shots} 次基于量子电路，按索引统计次数：")
    print(sample_counts)
    print(f"采样形状: {np.array(samples).shape}")
    print(f"前5个样本:")
    num_qubits_disp = int(np.log2(len(learned_probs)))
    for i in range(min(5, len(samples))):
        s_idx = indices[i]
        binary_str = format(s_idx, f'0{num_qubits_disp}b')
        print(f"  样本 {i+1}: {binary_str} (十进制: {s_idx})")

if __name__ == "__main__":
    main()