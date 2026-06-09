---
title: "复现 Koon et al. 2000: CR3BP 中的异宿连接与共振跃迁"
date: 2026-06-09
summary: "完整复现 Koon, Lo, Marsden, Ross 2000 年 Chaos 论文——在圆形限制性三体问题中，通过 Poincaré 截面证明 $L_1 \leftrightarrow L_2$ 异宿连接和同宿连接的存在。包含 Sun-Earth、Earth-Moon、Sun-Jupiter 三个系统的数值验证。"
tags: ["轨道力学", "CR3BP", "异宿连接", "不变流形", "天体力学", "论文复现"]
categories: ["技术笔记"]
showToc: true
---

## 引言

2000 年，Koon、Lo、Marsden 和 Ross 在 *Chaos* 期刊发表了一篇开创性论文：**"Heteroclinic Connections between Periodic Orbits and Resonance Transitions in Celestial Mechanics"**（DOI: 10.1063/1.166509）。该论文利用动力系统理论，在圆形限制性三体问题（CR3BP）中证明了 $L_1$ 和 $L_2$ 平动点周期轨道之间的异宿连接和同宿连接的存在，为低能量轨道转移和天体俘获提供了理论基础。

本文记录了对该论文核心结果的完整数值复现过程，涵盖三个系统：Sun-Earth、Earth-Moon 和 Sun-Jupiter。

## 背景知识

### 圆形限制性三体问题 (CR3BP)

在旋转坐标系中，两个大质量天体（如 Sun-Jupiter 或 Earth-Moon）固定不动，第三个小质量体（航天器或小天体）在它们的引力场中运动。归一化后，系统由单一参数 $\mu$（质量比）描述：

$$\mu = \frac{m_2}{m_1 + m_2}$$

其中 $m_1 > m_2$。例如：
- Sun-Jupiter: $\mu \approx 0.000954$
- Earth-Moon: $\mu \approx 0.01215$
- Sun-Earth: $\mu \approx 0.000003$

### $L_1$ 和 $L_2$ 平动点

每个 CR3BP 有五个平动点，其中 $L_1$（两体之间）和 $L_2$（小体外侧）是共线平动点，具有鞍点性质。围绕它们的 Lyapunov 轨道的不稳定/稳定流形构成了相空间中的"管道"。

### 异宿连接与同宿连接

- **异宿连接** (Heteroclinic)：$L_1$ 的不稳定流形 $W^u(L_1)$ 与 $L_2$ 的稳定流形 $W^s(L_2)$ 在 Poincaré 截面上相交
- **同宿连接** (Homoclinic)：$L_1$ 的不稳定流形 $W^u(L_1)$ 与 $L_1$ 的稳定流形 $W^s(L_1)$ 相交

当两组流形在截面上重叠时，存在连接两个周期轨道（或同一轨道）的异宿/同宿轨道。

## 方法论

### 核心算法流程

```
1. 计算 L1, L2 平动点位置
2. 在 L1 和 L2 处建立相同 C_J 的 Lyapunov 轨道
3. 计算单值矩阵 (Monodromy Matrix) 的特征值和特征向量
4. 沿不稳定/稳定特征向量微扰，计算流形
5. 在 Poincaré 截面上提取 (y, ẏ) 交叉数据
6. 搜索 W^u 和 W^s 的最近点对，计算距离 Δ
```

### $C_J$ Jacobi 常数匹配

在同一个 CR3BP 中，$C_J$ 是守恒量。异宿连接要求两个轨道具有相同的 $C_J$：

$$C_J = 2\Omega - (\dot{x}^2 + \dot{y}^2)$$

其中 $\Omega$ 是有效势。通过调整 Lyapunov 轨道的振幅 $\mathrm{d}x$，使 $L_1$ 和 $L_2$ 轨道的 $C_J$ 精确匹配（$\Delta C_J < 10^{-8}$）。

### Poincaré 截面

选择 $L_1$ 和 $L_2$ 中点作为截面位置 $x = (L_1 + L_2)/2$。当流形穿越截面时，记录 $(y, \dot{y})$ 坐标。如果 $W^u(L_1)$ 和 $W^s(L_2)$ 在此截面上有交点，就存在异宿连接。

### 距离度量

对于截面上的两组点集，定义最小距离：

$$\Delta = \min_{i,j} \sqrt{(y_i^u - y_j^s)^2 + (\dot{y}_i^u - \dot{y}_j^s)^2}$$

当 $\Delta < 0.01$ 时，认为连接存在。

## 实现细节

使用 Python + SciPy 实现，核心组件：

- **运动方程积分**：DOP853（8 阶 Runge-Kutta），$\mathrm{rtol} = 10^{-11}$，$\mathrm{atol} = 10^{-12}$
- **Lyapunov 轨道**：双条件微分修正（$y = 0,\; \dot{x} = 0$）
- **$C_J$ 匹配**：Brent 方法扫描轨道振幅
- **流形计算**：沿特征向量偏移 $\delta = 0.01 \times y_\text{amp}$，积分至离开 $R_{\max}$ 区域
- **Poincaré 截面**：线性插值检测穿越

## 结果

### Sun-Jupiter 系统（论文主要例子）

| 连接类型 | $C_J$ 范围 | 成功率 | 最佳 $\Delta$ |
|----------|-----------|--------|-------------|
| 异宿 $L_1 \leftrightarrow L_2$ | 3.020 ~ 3.036 | 8/8 | $1.2 \times 10^{-3}$ |
| 同宿 $L_1 \to L_1$ | 3.020 ~ 3.036 | 8/8 | $6 \times 10^{-6}$ |
| 同宿 $L_2 \to L_2$ | 3.020 ~ 3.036 | 8/8 | $4 \times 10^{-6}$ |

**关键发现：同宿连接比异宿连接精确约 200 倍。** 这符合物理直觉——同宿连接在同一轨道的流形管内，交叉更紧密。

### Sun-Earth 系统

| $C_J$ | $\Delta$ | 状态 |
|-------|---------|------|
| 3.000720 | $2.44 \times 10^{-4}$ | ✅ |
| 3.000773 | $1.76 \times 10^{-4}$ | ✅ |
| 3.000809 | $1.68 \times 10^{-4}$ | ✅ |
| 3.000844 | $1.73 \times 10^{-4}$ | ✅ |
| 3.000880 | $\mathbf{3.6 \times 10^{-5}}$ | ✅ 最佳 |

Sun-Earth 系统的异宿连接非常精确（$\Delta \sim 10^{-5}$），10/10 全部成功。

### Earth-Moon 系统

| $C_J$ | $\Delta$ | 状态 |
|-------|---------|------|
| 3.118 | $9.8 \times 10^{-3}$ | ✅ |
| 3.137 | $6.3 \times 10^{-3}$ | ✅ |
| 3.151 | $2.9 \times 10^{-3}$ | ✅ |
| 3.161 | $\mathbf{1.9 \times 10^{-3}}$ | ✅ 最佳 |
| 3.170 | $5.7 \times 10^{-3}$ | ✅ |

Earth-Moon 系统 $L_1$-$L_2$ 间距较大（0.32 EM 距离），10/12 成功，最佳 $\Delta = 1.9 \times 10^{-3}$。

### $C_J$ 参数扫描趋势

三个系统都表现出相同趋势：**$C_J$ 越接近 $L_2$ 平动点值，连接越精确。**

这可以通过以下物理解释：$C_J$ 接近 $L_2$ 点值时，Lyapunov 轨道非常小，流形传播距离短，误差累积少，截面上的交叉更紧密。

### 可视化

#### $\Delta$ vs $C_J$ 曲线

下图展示了三个系统中异宿连接距离 $\Delta$ 随 $C_J$ 的变化趋势：

![Sun-Jupiter $\Delta$ vs $C_J$](/images/koon-2000/fig_sj_delta_CJ.png)

#### Poincaré 截面

**Sun-Jupiter 异宿连接**——红色点云（$W^u(L_1)$）和蓝色点云（$W^s(L_2)$）的重叠区域就是异宿连接存在的证据：

![Sun-Jupiter Poincaré](/images/koon-2000/fig_sj_poincare_hetero.png)

**Sun-Jupiter 同宿连接（$L_1$）**——同一轨道的不稳定和稳定流形在截面上交叉：

![Sun-Jupiter Homoclinic $L_1$](/images/koon-2000/fig_sj_poincare_homo_L1.png)

**Earth-Moon Poincaré 截面**：

![EM Poincaré](/images/koon-2000/fig2_poincare_EM.png)

#### 相图

**Sun-Jupiter 相图**——$L_1$ 和 $L_2$ 的 Lyapunov 轨道及其不稳定（红色）和稳定（蓝色）流形管：

![Sun-Jupiter Phase](/images/koon-2000/fig_sj_phase_hetero.png)

**Earth-Moon 相图**：

![EM Phase](/images/koon-2000/fig3_phase_EM.png)

## 与论文的对比

| 论文内容 | 复现状态 |
|----------|----------|
| 同系统 $L_1 \leftrightarrow L_2$ 异宿连接 | ✅ 三个系统 |
| 同宿连接 $L_1 \to L_1$，$L_2 \to L_2$ | ✅ Sun-Jupiter |
| $C_J$ 参数扫描 | ✅ 8-12 个能量值 |
| Poincaré 截面流形管交叉 | ✅ 所有连接类型 |
| Sun-Jupiter 主要例子 | ✅ |
| $\Delta$ vs $C_J$ 趋势分析 | ✅ 与论文一致 |

## 关键经验

1. **截面位置很重要**：最初将 Poincaré 截面放在 Moon 位置或 Hill 球边界都失败了。正确的做法是在 $L_1$-$L_2$ 中点处截面。

2. **$C_J$ 匹配是前提**：跨系统（SE↔EM）的精确流形匹配是不可能的——论文的 patched three-body 方法使用小 $\Delta V$ 机动在 patch point 连接，而非精确匹配。

3. **同系统内 $C_J$ 自动相同**：在同一 CR3BP 中，只要 $L_1$ 和 $L_2$ 轨道的 $C_J$ 匹配，流形就在同一个能量面上，自然可以相交。

4. **几何方法 vs 轨道力学方法**：论文用 Poincaré 截面上流形管重叠来**证明连接存在**，不需要微分修正到精确轨道。

## 代码仓库

完整代码在 `/home/ubuntu/.myclaw/workspace/`：

- `compute_manifolds.py` — 核心计算库
- `step1_cj_sweep.py` — $C_J$ 参数扫描
- `step4_sun_jupiter.py` — Sun-Jupiter 系统完整复现
- `cj_sweep.json`, `sj_results.json` — 数据文件
- `fig_*.png` — 可视化图像

## 参考文献

1. Koon, W. S., Lo, M. W., Marsden, J. E., & Ross, S. D. (2000). Heteroclinic connections between periodic orbits and resonance transitions in celestial mechanics. *Chaos*, 10(2), 427-469.
2. Szebehely, V. (1967). *Theory of Orbits: The Restricted Problem of Three Bodies*. Academic Press.
