---
title: "复现 Hargraves & Paris 1987：用非线性规划和配点法做轨迹优化"
date: 2026-07-06
summary: "复现经典论文 Direct Trajectory Optimization Using Nonlinear Programming and Collocation：实现三次 Hermite 状态插值、控制线性插值、中心点配点缺陷，并用 brachistochrone 解析解验证数值收敛。"
tags: ["轨迹优化", "最优控制", "直接配点法", "非线性规划", "论文复现", "Python"]
categories: ["技术笔记"]
showToc: true
---

## 为什么复现这篇论文

Hargraves 和 Paris 在 1987 年发表的 **Direct Trajectory Optimization Using Nonlinear Programming and Collocation** 是直接轨迹优化路线中的一篇经典短文。它的核心思想很工程化：

> 不先推导庞大的协态方程，而是直接把状态、控制和时间离散成有限维变量，再用非线性规划求解。

这条路线后来成为很多轨迹优化工具的基础思路。现代工程中常见的 direct collocation、direct transcription、Hermite-Simpson、pseudospectral method，本质上都沿着这条路继续发展。

我这次复现的目标不是复制 Boeing 当年的 NPDOT 程序，而是复现论文中最关键的数值转录方法：

- 状态变量用三次 Hermite 多项式表示；
- 控制变量在节点之间线性插值；
- 动力学方程在每段中心点强制满足；
- 轨迹优化问题转成 NLP；
- 用 SQP 类优化器求解。

代码和验证结果保存在本地：

```text
/home/ubuntu/repro/direct-collocation
```

## 论文方法：把最优控制问题变成 NLP

论文考虑一般形式的轨迹优化问题：

```text
minimize J = Φ[x(E), u(E), ω, E]
subject to x' = f(x, u, ω, t)
```

其中状态、控制、事件时间、设计参数一起构成 NLP 的决策变量。论文把所有变量收集为一个向量：

```text
P = [Z, E, u]
```

然后把动力学缺陷、边界条件和路径约束都写成非线性约束：

```text
C(P) = 0 或 lower <= C(P) <= upper
```

最终得到标准非线性规划问题。

### 三次 Hermite 状态插值

对每一个离散段，状态用三次多项式表示：

```text
x(s) = C0 + C1 s + C2 s^2 + C3 s^3,   s ∈ [0, 1]
```

端点状态为 `x_i, x_{i+1}`，端点导数由动力学给出：

```text
f_i     = f(x_i, u_i)
f_{i+1} = f(x_{i+1}, u_{i+1})
```

段长为 `h`。根据 Hermite 插值，中心点状态为：

```text
x_c = 0.5 * (x_i + x_{i+1}) + h/8 * (f_i - f_{i+1})
```

中心点处由多项式导出的斜率为：

```text
slope_c = 1.5 * (x_{i+1} - x_i) / h - 0.25 * (f_i + f_{i+1})
```

控制变量线性插值得到中心控制：

```text
u_c = 0.5 * (u_i + u_{i+1})
```

于是每一段的 collocation defect 是：

```text
defect = slope_c - f(x_c, u_c) = 0
```

这就是我复现代码中的核心。

## Python 实现

使用 Python + SciPy 实现 NLP。核心函数如下：

```python
def defects(z, p):
    x, y, v, theta, tf = unpack(z, p)
    h = tf / p.n_segments
    cons = []
    for i in range(p.n_segments):
        yi = np.array([x[i], y[i], v[i]])
        yj = np.array([x[i + 1], y[i + 1], v[i + 1]])
        fi = dynamics(yi, theta[i], p)
        fj = dynamics(yj, theta[i + 1], p)

        yc = 0.5 * (yi + yj) + h / 8.0 * (fi - fj)
        thetac = 0.5 * (theta[i] + theta[i + 1])
        fc = dynamics(yc, thetac, p)

        slope_c = 1.5 * (yj - yi) / h - 0.25 * (fi + fj)
        cons.extend(slope_c - fc)
    return np.asarray(cons)
```

优化器使用 `scipy.optimize.minimize(method="SLSQP")`。这和论文里使用 NPSOL 的路线一致：都是 SQP 类非线性规划求解器。

## 验证问题：brachistochrone

论文提到 NPDOT 验证过四类例子：

1. Van der Pol；
2. Brachistochrone；
3. Supersonic interceptor minimum-time climb；
4. Advanced booster trajectory。

但这篇论文只有 5 页，很多工程例子的气动、推进、燃耗表并没有直接给出。因此最适合独立验证的是 brachistochrone，因为无约束版本有解析 cycloid 解。

我设置的验证问题为：

```text
x(0) = 0, y(0) = 1, v(0) = 0
x(tf) = 1, y(tf) = 0
g = 1
minimize tf
```

动力学为：

```text
x_dot = v cos(theta)
y_dot = v sin(theta)
v_dot = -g sin(theta)
```

这里 `theta` 是控制角，`y` 正方向向上。为了下降，优化器会选择负的 `theta`。

## 解析解对照

无约束 brachistochrone 的解析解是摆线：

```text
x = a (phi - sin phi)
y_drop = a (1 - cos phi)
tf = phi * sqrt(a/g)
```

对于起点 `(0, 1)`、终点 `(1, 0)`、`g = 1`，求得：

```text
phi_f = 2.412011143914
a     = 0.572917037532
tf    = 1.825682189397
```

这给了我们一个很干净的 benchmark：数值配点结果应该随着段数增加逐渐收敛到 `1.825682189397`。

## 收敛结果

使用不同段数运行直接配点 NLP：

| segments | success | final time | analytic time | abs error | rel error | max defect residual |
|---:|---:|---:|---:|---:|---:|---:|
| 5 | True | 1.825708717 | 1.825682189 | 2.653e-05 | 1.453e-05 | 3.525e-11 |
| 10 | True | 1.825693832 | 1.825682189 | 1.164e-05 | 6.377e-06 | 1.129e-10 |
| 20 | True | 1.825687905 | 1.825682189 | 5.715e-06 | 3.131e-06 | 7.892e-11 |

可以看到两个结果：

1. 终端时间随着段数增加收敛到解析解；
2. 动力学配点残差保持在 `1e-10` 量级。

这说明 Hermite midpoint collocation 的转录本身是正确工作的。

下图是 `N=20` 时，数值轨迹与解析摆线的对比：

![Brachistochrone analytic vs collocation](/images/direct-collocation/unconstrained_vs_analytic.png)

单独看数值节点和 Hermite 轨迹：

![Unconstrained N20 trajectory](/images/direct-collocation/unconstrained_N20_trajectory.png)

速度和控制角变化：

![Unconstrained N20 state and control](/images/direct-collocation/unconstrained_N20_state_control.png)

## 路径约束验证

论文中特别强调直接法的一个优势：容易处理路径约束。为此我也做了一个 constrained brachistochrone-like 版本，在节点和中心点都检查路径约束：

```text
y >= smooth obstacle / floor
```

对应代码中，节点约束和中心点约束都会进入 NLP：

```python
vals.extend(y - floor_y(x, p))
...
for each segment:
    yc = Hermite midpoint state
    vals.append(yc[1] - floor_y(yc[0], p))
```

收敛结果：

| segments | success | final time | max defect residual |
|---:|---:|---:|---:|
| 5 | True | 1.825708717 | 8.431e-13 |
| 10 | True | 1.825693831 | 7.761e-11 |
| 20 | True | 1.825687905 | 1.026e-10 |

当前这个路径约束设置没有显著改变最优轨迹，因此数值和无约束情况接近。但这个实验验证了两点：

- 路径约束能以普通 NLP inequality 的形式加入；
- 约束可以同时在节点和段中心点检查，这和论文描述一致。

![Constrained N20 trajectory](/images/direct-collocation/constrained_N20_trajectory.png)

## 论文中工程例子的复现边界

复现过程中一个重要原则是：**论文没给的数据，不要自己编。**

### Supersonic interceptor minimum-time climb

论文中给出的可用信息包括：

```text
初始：sea level, Mach 0.38
终端：altitude 20 km, Mach 1.0
初始重量：42,000 lb
控制：pitch function
初始猜测：final range = 360,000 ft
初始猜测：final mass = 1204 slugs
初始控制：pitch = 0.18 rad
8 段 Chebyshev 分布结果：time to climb = 325.2 s
15 段非均匀分布结果：time to climb = 317.3 s
```

但它没有给出完整的：

- 气动系数表；
- 推力随 Mach/高度变化表；
- 燃耗数据；
- NACA-1962 atmosphere 的具体实现细节；
- Ref. 35 中使用的原始数据。

所以我没有声称复现了 `317.3 s`。准确复现这个例子，需要继续找到 Bryson、Desai、Hoffman 1969 年那篇 *Energy-State Approximation in Performance Optimization of Supersonic Aircraft* 的数据表。

### Advanced booster

论文给出的信息是：

```text
Mach 0-2：第一推进阶段
Mach 2-20：第二推进阶段
Mach >20：火箭阶段
三维球形非旋转地球动力学
路径约束：动压、高度、攻角
NPDOT 使用 22 段
最终比 baseline 重量提升 10%
```

但同样缺少车辆、气动、推进、目标轨道和 baseline 轨迹数据。因此只能记录论文报告值，不能精确数值复现。

## 这次复现得到的经验

### 1. 直接配点法的核心很短

真正核心的数学公式只有几行：

```text
x_c = 0.5 * (x_i + x_{i+1}) + h/8 * (f_i - f_{i+1})
slope_c = 1.5 * (x_{i+1} - x_i) / h - 0.25 * (f_i + f_{i+1})
defect = slope_c - f(x_c, u_c)
```

但这几行把连续时间最优控制问题变成了有限维 NLP。

### 2. 路径约束是直接法的强项

间接法处理路径约束时通常会牵涉复杂的切换结构、乘子条件和边界弧。直接法则可以直接写成：

```text
g(x_i, u_i) >= 0
g(x_c, u_c) >= 0
```

然后交给 NLP 求解器。

### 3. 论文复现要区分“方法复现”和“数据复现”

这篇论文的方法已经可以完整复现；但工程例子因为缺数据，不能完整复现论文中的数值表。把这两件事区分开很重要：

- 方法复现：实现 Hermite collocation + NLP，并用可验证问题证明正确；
- 数据复现：需要论文或引用文献中完整的气动、推进、环境模型和约束数据。

### 4. SLSQP 能跑通，但不是最终形态

SciPy SLSQP 足够验证方法，但对于更大的飞行器轨迹优化问题，建议迁移到：

- CasADi + IPOPT；
- JAX / CasADi 自动微分；
- 稀疏 Jacobian；
- mesh refinement。

这会更接近现代轨迹优化工具链。

## 下一步

后续可以沿两条路线继续：

1. **算法路线**：实现自动网格细化，比较 10、20、40、80 段下的误差和计算时间；
2. **工程路线**：寻找 Ref. 35 的超音速飞机数据表，尝试复现论文报告的 `317.3 s` minimum-time climb。

这篇论文的价值在于，它把最优控制问题从“解析推导协态方程”推进到了“构造有限维 NLP 并交给优化器”。这也是现在很多工程轨迹优化工具仍在使用的基本思想。

## 参考文献

1. Hargraves, C. R., & Paris, S. W. (1987). Direct Trajectory Optimization Using Nonlinear Programming and Collocation. *Journal of Guidance, Control, and Dynamics*, 10(4), 338–342.
2. Bryson, A. E., & Ho, Y. C. (1969). *Applied Optimal Control*. Blaisdell.
3. Bryson, A. E., Desai, M. N., & Hoffman, W. C. (1969). Energy-State Approximation in Performance Optimization of Supersonic Aircraft. *Journal of Aircraft*, 6(6), 481–488.
