# 按 Featherstone Table 5.1 逐行推导二连杆 6D RNEA

> **目标**：以《Rigid Body Dynamics Algorithms》（Roy Featherstone，2008）为依据，针对**固定基座、两个标准旋转关节（2R）、无外部接触力**的串联机器人，从 Table 5.1 第 1 行走到最后 1 行。所有运动向量、加速度、空间力始终使用 **6×1**，所有空间变换、空间惯性始终使用 **6×6**；明确给出 `XJ`、`XT`、`Xup`、`S1`、`S2`、`vJ1`、`vJ2`、`cJ1`、`cJ2` 等量的**每个元素**，先符号推导，再代入一组数值。
>
> **来源范围**：本文件依据原书相关公式搭建 2R 示例；**书中的 Table 5.1 是通用算法，下面的 2R 符号展开、具体参数以及数值计算是按该算法自行推导的示例，并非原书中直接给出的二连杆算例。**
>
> **页码说明**：文内的“正文页”指书本印刷页码，“PDF 页”指阅读器从 1 开始计数的页码；这份 279 页的上传 PDF，在此章节范围内满足 **PDF 页 = 正文页 + 7**。

## 目录

1. [原书公式与准确页码索引](#1-原书公式与准确页码索引)
2. [建立 2R 模型和六维符号约定](#2-建立-2r-模型和六维符号约定)
3. [Table 5.1：准备关节数据 `XJ`、`XT`、`S`、`vJ`、`cJ`](#3-table-51准备关节数据-xjxtsvjcj)
4. [Table 5.1：初始化 `v0`、`a0`](#4-table-51初始化-v0a0)
5. [Table 5.1：第一次循环，`i=1`，第一连杆](#5-table-51第一次循环-i1第一连杆)
6. [Table 5.1：第一次循环，`i=2`，第二连杆](#6-table-51第一次循环-i2第二连杆)
7. [空间惯性矩阵 `I1` 和 `I2` 为什么长这样](#7-空间惯性矩阵-i1-和-i2-为什么长这样)
8. [Table 5.1：本体空间力 `f1B`、`f2B`](#8-table-51本体空间力-f1bf2b)
9. [Table 5.1：第二次循环，`i=2→1`](#9-table-51第二次循环-i21)
10. [最终 `τ`、`M(q)`、`C(q,dq)`、`G(q)`](#10-最终-τmqcq-dqgq)
11. [代入数值，把每个中间变量全部算出来](#11-代入数值把每个中间变量全部算出来)
12. [一份可独立运行的 MATLAB 符号复现程序](#12-一份可独立运行的-matlab-符号复现程序)
13. [与 Table 5.1 的逐行核对表](#13-与-table-51-的逐行核对表)

---

## 1. 原书公式与准确页码索引

| 本文需要的知识 | 原书位置 | 正文页 | PDF 页 | 公式或表编号 |
|---|---|---:|---:|---|
| 三维反对称叉乘矩阵 $[r]_\times$ | §2.8 Coordinate Transforms | 21 | 28 | **(2.23)** |
| 六维空间运动坐标变换 | §2.8 Coordinate Transforms | 22 | 29 | **(2.24)** |
| 六维空间力坐标变换 | §2.8 Coordinate Transforms | 22 | 29 | **(2.25)** |
| 六维运动叉乘算子 | §2.9 Spatial Cross Products | 25 | 32 | **(2.31)** |
| 六维力叉乘算子 | §2.9 Spatial Cross Products | 25 | 32 | **(2.32)** |
| 空间动量关系 $\mathsf H=\mathcal Iv$ | §2.13 Inertia | 32 | 39 | **(2.61)** |
| 质心处空间惯性 | §2.13 Inertia | 32 | 39 | **(2.62)** |
| **任意坐标原点的 6×6 空间惯性矩阵** | **§2.13 Inertia** | **33** | **40** | **(2.63)** |
| 旋转关节的局部运动子空间 $S$、关节变换 | §4.4 Joint Models | 79 | 86 | **Table 4.1** |
| 六维刚体空间力表达式 | §5.3 RNEA | 93 | 100 | **(5.9)** |
| 关节广义力投影 | §5.3 RNEA | 94 | 101 | **(5.11)** |
| 连杆坐标系的速度递推 | §5.3 RNEA | 95 | 102 | **(5.14)** |
| 含 $\mathring S_i$ 的加速度递推 | §5.3 RNEA | 95 | 102 | **(5.15)** |
| $v_{J_i}=S_i\dot q_i$ | §5.3 RNEA | 95 | 102 | **(5.16)** |
| $c_{J_i}=\mathring S_i\dot q_i$ | §5.3 RNEA | 95 | 102 | **(5.17)** |
| 整理后的速度递推 | §5.3 RNEA | 95 | 102 | **(5.18)** |
| 整理后的加速度递推 | §5.3 RNEA | 95 | 102 | **(5.19)** |
| 反向力递推 | §5.3 RNEA | 95 | 102 | **(5.20)** |
| **前向循环 + 反向循环伪代码** | **§5.3 RNEA** | **96** | **103** | **Table 5.1** |

**使用建议：翻到 PDF 第 103 页（书中第 96 页），把 Table 5.1 放在一边，对照本文件第 3–9 节逐行计算。**

原书 Table 5.1 除了上述数学式，还通过 `jcalc` 求关节运动学参数：

$$
[X_{J_i},\ S_i,\ v_{J_i},\ c_{J_i}]=\operatorname{jcalc}(\mathrm{jtype}_i,q_i,\dot q_i),
$$

并计算

$$
{}^iX_{\lambda(i)}=X_{J_i}X_{T_i}.
$$

其中 $X_{T_i}$ 是**固定几何变换**，$X_{J_i}$ 是**随关节角变化的运动变换**。这两个变换不能遗漏。

---

## 2. 建立 2R 模型和六维符号约定

### 2.1 二连杆的几何假设

```text
             Y0 ↑
                |
                |                           C2 ● ——→（连杆2末端方向）
                |                         /
                |                    L2 /   关节2旋转轴：+Z2
                |                       /
                |          O2 ●────────
                |           /    ∠q2 是相对连杆1的关节角
                |       L1 /
                |         /
                | O1 ●──────→ X0
                +------------------------------------------→
                    关节1旋转轴：+Z1

  俯视/图示在 XY 平面，Z 垂直纸面指向外侧。
  实际 q1、q2 任意；图只表示连杆间的连接关系，不按比例。
```

三个连杆固连坐标系：

- $\{0\}$：固定基座，原点 $O_1$；
- $\{1\}$：原点 $O_1$，局部 $x_1$ 沿第一连杆；
- $\{2\}$：原点 $O_2$，局部 $x_2$ 沿第二连杆；
- $O_1O_2=l_1$，第二连杆长度 $l_2$；
- 第一、二连杆的质心都位于其**各自的局部 x 轴**，距关节原点分别为 $r_1,r_2$；
- 两关节均绕各自的 $+z$ 轴旋转，旋转角 $q_1,q_2$，角速度 $\dot q_1,\dot q_2$，角加速度 $\ddot q_1,\ddot q_2$；
- 世界重力加速度 $^0\boldsymbol g=[0,-g,0]^T$；
- 没有末端载荷或外部接触力，暂时没有关节摩擦。

使用简写：

$$
\begin{aligned}
c_1=\cos q_1,\quad&s_1=\sin q_1,\\
c_2=\cos q_2,\quad&s_2=\sin q_2,\\
c_{12}=\cos(q_1+q_2),\quad&s_{12}=\sin(q_1+q_2).
\end{aligned}
$$

> 注意：这里的 $c_1,s_1$ 是余弦/正弦简写，**不是**质心向量，也不是运动子空间 $S_1$。质心向量记为 $\boldsymbol r_{C_i}$。

### 2.2 为什么自由度是 2，向量却是 6 维？

关节变量：

$$
q=\begin{bmatrix}q_1\\q_2\end{bmatrix}\in\mathbb R^2.
$$

对**任一根刚体**，空间运动向量按“**角速度在前、线速度在后**”排列：

$$
\boxed{v_i=\begin{bmatrix}\omega_{ix}\\\omega_{iy}\\\omega_{iz}\\u_{ix}\\u_{iy}\\u_{iz}\end{bmatrix}\in\mathbb R^6}
$$

对应空间加速度：

$$
\boxed{a_i=\begin{bmatrix}\alpha_{ix}\\\alpha_{iy}\\\alpha_{iz}\\a_{ix}\\a_{iy}\\a_{iz}\end{bmatrix}\in\mathbb R^6}
$$

对应空间力：

$$
\boxed{f_i=\begin{bmatrix}n_{ix}\\n_{iy}\\n_{iz}\\F_{ix}\\F_{iy}\\F_{iz}\end{bmatrix}\in\mathbb R^6}
$$

这里 $a_i$ 是**Featherstone 的空间加速度**；其线加速度部分并非在所有情况下都等于连杆坐标原点的经典线加速度，尤其不能忽视虚拟基座重力加速度处理。

---

## 3. Table 5.1：准备关节数据 `XJ`、`XT`、`S`、`vJ`、`cJ`

**原书位置**：Table 4.1（正文 p.79 / PDF p.86）；空间变换式 (2.24)（正文 p.22 / PDF p.29）；关节速度 (5.16)、偏置加速度 (5.17)（正文 p.95 / PDF p.102）；Table 5.1（正文 p.96 / PDF p.103）。

### 3.1 两个关节的 $S_1,S_2$ 必须明确写出来

每个关节都是**单自由度旋转关节，轴在局部 +Z 方向**：

$$
\boxed{S_1=\begin{bmatrix}0\\0\\1\\0\\0\\0\end{bmatrix},\qquad S_2=\begin{bmatrix}0\\0\\1\\0\\0\\0\end{bmatrix}}
$$

尺寸：$S_1,S_2\in\mathbb R^{6\times1}$。

前三行是角速度方向；后三行为零，表示**该关节相对运动没有沿关节轴以外的纯平移自由度**。

### 3.2 两个关节的相对空间速度 $v_{J1},v_{J2}$

依据原书 **(5.16)**：

$$
v_{J_i}=S_i\dot q_i.
$$

逐个写：

$$
\boxed{v_{J1}=\begin{bmatrix}0\\0\\\dot q_1\\0\\0\\0\end{bmatrix},\qquad
v_{J2}=\begin{bmatrix}0\\0\\\dot q_2\\0\\0\\0\end{bmatrix}}
$$

**不要把 $v_{J2}$ 与连杆 2 的总速度 $v_2$ 混为一谈**：$v_{J2}$ 只有第二关节自身的相对转动；$v_2$ 还要加上第一连杆传递来的运动。

### 3.3 两个关节的偏置加速度 $c_{J1},c_{J2}$

依据原书 **(5.17)**：

$$
c_{J_i}=\mathring S_i\dot q_i.
$$

这里 $S_i$ 在各自的连杆固连坐标系中是常量，因此

$$
\mathring S_1=\mathring S_2=\begin{bmatrix}0\\0\\0\\0\\0\\0\end{bmatrix}.
$$

得到明确的数值：

$$
\boxed{c_{J1}=\begin{bmatrix}0\\0\\0\\0\\0\\0\end{bmatrix},\qquad
c_{J2}=\begin{bmatrix}0\\0\\0\\0\\0\\0\end{bmatrix}}
$$

$\mathring S_i=0$ **不代表** $v_i\times v_{J_i}=0$；后者必须独立计算。

### 3.4 根据原书的空间变换 (2.24) 构建 $X_{J_i}$ 和 $X_{T_i}$

先定义绕 Z 轴的主动旋转矩阵（将局部向量表示转入父坐标系）：

$$
R_z(q)=\begin{bmatrix}\cos q&-\sin q&0\\\sin q&\cos q&0\\0&0&1\end{bmatrix}.
$$

Featherstone 的**父系 → 子系**空间运动变换使用 $E_i=R_z(q_i)^T$：

$$
E_i=\begin{bmatrix}c_i&s_i&0\\-s_i&c_i&0\\0&0&1\end{bmatrix}.
$$

关节转动部分：

$$
X_{J_i}=\begin{bmatrix}E_i&0\\0&E_i\end{bmatrix}.
$$

对关节 1：关节原点 $O_1$ 与基座原点重合，故

$$
\boxed{X_{T1}=I_6,\qquad X_{J1}=\begin{bmatrix}
c_1&s_1&0&0&0&0\\
-s_1&c_1&0&0&0&0\\
0&0&1&0&0&0\\
0&0&0&c_1&s_1&0\\
0&0&0&-s_1&c_1&0\\
0&0&0&0&0&1
\end{bmatrix}}
$$

因此：

$$\boxed{X_1=X_{J1}X_{T1}=X_{J1}.}$$

对关节 2：在第一连杆坐标系中，$O_2$ 相对 $O_1$ 的**固定位置**为

$$
p_{12}=\begin{bmatrix}l_1\\0\\0\end{bmatrix}.
$$

使用三维叉乘矩阵

$$
[p_{12}]_\times=\begin{bmatrix}0&0&0\\0&0&-l_1\\0&l_1&0\end{bmatrix},
$$

固定平移部分为

$$
X_{T2}=\begin{bmatrix}I_3&0\\-[p_{12}]_\times&I_3\end{bmatrix}
=\boxed{\begin{bmatrix}
1&0&0&0&0&0\\
0&1&0&0&0&0\\
0&0&1&0&0&0\\
0&0&0&1&0&0\\
0&0&l_1&0&1&0\\
0&-l_1&0&0&0&1
\end{bmatrix}}.
$$

第二关节转动部分：

$$
X_{J2}=\begin{bmatrix}E_2&0\\0&E_2\end{bmatrix}
=\begin{bmatrix}
c_2&s_2&0&0&0&0\\
-s_2&c_2&0&0&0&0\\
0&0&1&0&0&0\\
0&0&0&c_2&s_2&0\\
0&0&0&-s_2&c_2&0\\
0&0&0&0&0&1
\end{bmatrix}.
$$

根据 Table 5.1：

$$
\boxed{X_2=X_{J2}X_{T2}}
$$

做 6×6 矩阵乘法，得到：

$$
\boxed{X_2=\begin{bmatrix}
c_2&s_2&0&0&0&0\\
-s_2&c_2&0&0&0&0\\
0&0&1&0&0&0\\
0&0&l_1s_2&c_2&s_2&0\\
0&0&l_1c_2&-s_2&c_2&0\\
0&-l_1&0&0&0&1
\end{bmatrix}.}
$$

其中 $X_i\in\mathbb R^{6\times6}$。$X_2$ 左下角的 $l_1s_2,l_1c_2$ 表示：**父连杆角速度会使下一关节原点产生线速度**。

### 3.5 对 Table 5.1 `jcalc` 输出的总结

| Table 5.1 变量 | 关节 1 | 关节 2 | 维度 |
|---|---|---|---|
| `jtype(i)` | `R`（局部 +Z） | `R`（局部 +Z） | — |
| $q_i$ | $q_1$ | $q_2$ | 标量 |
| $S_i$ | $[0,0,1,0,0,0]^T$ | $[0,0,1,0,0,0]^T$ | 6×1 |
| $v_{J_i}$ | $[0,0,\dot q_1,0,0,0]^T$ | $[0,0,\dot q_2,0,0,0]^T$ | 6×1 |
| $\mathring S_i$ | 零 6×1 向量 | 零 6×1 向量 | 6×1 |
| $c_{J_i}$ | 零 6×1 向量 | 零 6×1 向量 | 6×1 |
| $X_{J_i}$ | $\operatorname{diag}(E_1,E_1)$ | $\operatorname{diag}(E_2,E_2)$ | 6×6 |
| $X_{T_i}$ | $I_6$ | $\operatorname{xlt}([l_1,0,0]^T)$ | 6×6 |
| $X_i$ | $X_{J1}$ | $X_{J2}X_{T2}$ | 6×6 |

---

## 4. Table 5.1：初始化 `v0`、`a0`

**原书位置**：重力作为虚拟基座加速度的讨论在 §5.3（正文 p.94 / PDF p.101）；Table 5.1（正文 p.96 / PDF p.103）。

机器人基座固定：

$$
\boxed{v_0=\begin{bmatrix}0\\0\\0\\0\\0\\0\end{bmatrix}.}
$$

世界重力加速度：

$$
{}^0\boldsymbol g=\begin{bmatrix}0\\-g\\0\end{bmatrix}.
$$

由于原书采用虚拟基座加速度，设：

$$
\boxed{a_0=-\begin{bmatrix}0\\0\\0\\0\\-g\\0\end{bmatrix}
=\begin{bmatrix}0\\0\\0\\0\\g\\0\end{bmatrix}.}
$$

它是算法技巧，**不是**说固定基座真的向上加速。

---

## 5. Table 5.1：第一次循环，`i=1`，第一连杆

**原书位置**：(5.16)–(5.19)，正文 p.95 / PDF p.102；Table 5.1，正文 p.96 / PDF p.103。

### 5.1 计算 $v_1$（原书 5.18）

$$
v_1=X_1v_0+v_{J1}.
$$

代入 $v_0=0$：

$$
\boxed{v_1=\begin{bmatrix}0\\0\\\dot q_1\\0\\0\\0\end{bmatrix}.}
$$

这很合理：连杆 1 的坐标原点就在固定的第一关节处，因此原点线速度为零；绕 Z 轴的角速度是 $\dot q_1$。

### 5.2 计算 $a_1$（原书 5.19）

$$
a_1=X_1a_0+S_1\ddot q_1+c_{J1}+\operatorname{crm}(v_1)v_{J1}.
$$

四项逐个代入：

$$
X_1a_0=\begin{bmatrix}0\\0\\0\\gs_1\\gc_1\\0\end{bmatrix},\qquad
S_1\ddot q_1=\begin{bmatrix}0\\0\\\ddot q_1\\0\\0\\0\end{bmatrix},
$$

$$
c_{J1}=\boldsymbol0_6,\qquad
\operatorname{crm}(v_1)v_{J1}=\boldsymbol0_6.
$$

最后：

$$
\boxed{a_1=\begin{bmatrix}
0\\0\\\ddot q_1\\gs_1\\gc_1\\0
\end{bmatrix}.}
$$

**不要因为 `cJ1=0` 就删去所有关节的叉乘项**；只是在本例连杆 1 的当前状态下，空间速度叉乘碰巧为零。

---

## 6. Table 5.1：第一次循环，`i=2`，第二连杆

### 6.1 计算 $v_2$（原书 5.18）

$$
v_2=X_2v_1+v_{J2}.
$$

第一项：

$$
X_2v_1=\begin{bmatrix}
0\\0\\\dot q_1\\l_1s_2\dot q_1\\l_1c_2\dot q_1\\0
\end{bmatrix}.
$$

第二项：

$$
v_{J2}=\begin{bmatrix}0\\0\\\dot q_2\\0\\0\\0\end{bmatrix}.
$$

相加得到：

$$
\boxed{v_2=\begin{bmatrix}
0\\0\\\dot q_1+\dot q_2\\l_1s_2\dot q_1\\l_1c_2\dot q_1\\0
\end{bmatrix}.}
$$

解释：$v_2$ 的第三分量是两关节角速度之和；第四、五分量是关节 1 旋转给第二关节原点带来的线速度。

### 6.2 计算 $a_2$（原书 5.19）

$$
a_2=X_2a_1+S_2\ddot q_2+c_{J2}+\operatorname{crm}(v_2)v_{J2}.
$$

依次计算：

**第一项**（父连杆空间加速度传递）：

$$
\boxed{X_2a_1=\begin{bmatrix}
0\\0\\\ddot q_1\\l_1s_2\ddot q_1+gs_{12}\\l_1c_2\ddot q_1+gc_{12}\\0
\end{bmatrix}.}
$$

**第二项**（关节 2 自身角加速度）：

$$
\boxed{S_2\ddot q_2=\begin{bmatrix}0\\0\\\ddot q_2\\0\\0\\0\end{bmatrix}.}
$$

**第三项**（标准旋转关节局部运动子空间恒定）：

$$\boxed{c_{J2}=\begin{bmatrix}0\\0\\0\\0\\0\\0\end{bmatrix}.}$$

**第四项**（空间叉乘，原书 (2.31)，正文 p.25 / PDF p.32）：

定义

$$
w=\dot q_1+\dot q_2,\quad u_x=l_1s_2\dot q_1,\quad u_y=l_1c_2\dot q_1.
$$

$\operatorname{crm}(v_2)$ 的具体值是

$$
\operatorname{crm}(v_2)=\begin{bmatrix}
0&-w&0&0&0&0\\
w&0&0&0&0&0\\
0&0&0&0&0&0\\
0&0&u_y&0&-w&0\\
0&0&-u_x&w&0&0\\
-u_y&u_x&0&0&0&0
\end{bmatrix}.
$$

所以：

$$
\boxed{\operatorname{crm}(v_2)v_{J2}
=\begin{bmatrix}
0\\0\\0\\l_1c_2\dot q_1\dot q_2\\-l_1s_2\dot q_1\dot q_2\\0
\end{bmatrix}.}
$$

将四项相加：

$$
\boxed{a_2=\begin{bmatrix}
0\\
0\\
\ddot q_1+\ddot q_2\\
l_1s_2\ddot q_1+l_1c_2\dot q_1\dot q_2+gs_{12}\\
l_1c_2\ddot q_1-l_1s_2\dot q_1\dot q_2+gc_{12}\\
0
\end{bmatrix}.}
$$

---

## 7. 空间惯性矩阵 $I_1$ 和 $I_2$ 为什么长这样

**原书位置**：§2.13（正文 p.32–33 / PDF p.39–40）；空间动量 (2.61)、质心处惯性 (2.62)、任意原点惯性 **(2.63)**。

### 7.1 论文中的 $I_i$ 对应 $I_{xx}$、$I_{yy}$ 还是 $I_{zz}$？

截图对应论文 [*Three High Performance Global Tracking Composite Adaptive Controllers for Fully Actuated Euler–Lagrange Systems: Experimental Validation*](https://arxiv.org/html/2609.27010v1#S5)。论文 §5 的 Table 1 明确将 $I_1,I_2$ 标为关于质心的惯量；式 (64)–(67) 给出二连杆动力学和参数组合。

**在本文与 `RNEA_SYM.m` 的建模约定下，关节轴为局部 $+z$，因此论文的标量 $I_i$ 对应质心转动惯量张量的 $I_{zz,i}$，也就是代码中的 `I_zz(i)`。** 这个轴分量的对应关系是根据本文的坐标轴约定推导的；论文的平面示意图本身没有标出三维 $x,y,z$ 轴。

设质心处、沿连杆坐标轴方向表达的转动惯量张量为

$$
\bar I_{C_i}=\begin{bmatrix}
I_{xx,i}&I_{xy,i}&I_{xz,i}\\
I_{xy,i}&I_{yy,i}&I_{yz,i}\\
I_{xz,i}&I_{yz,i}&I_{zz,i}
\end{bmatrix}.
$$

因为 $\omega_i=[0,0,\omega_{iz}]^T$，旋转动能为

$$
T_{\mathrm{rot},i}=\frac12\omega_i^T\bar I_{C_i}\omega_i
=\frac12 I_{zz,i}\omega_{iz}^2.
$$

从而

$$
\boxed{I_i^{\mathrm{paper}}=e_z^T\bar I_{C_i}e_z=I_{zz,i},\qquad e_z=[0,0,1]^T.}
$$

若另一套坐标系把同一关节轴定义为 $x$ 或 $y$，则它分别对应 $I_{xx}$ 或 $I_{yy}$；对于任意单位轴向量 $e$，对应标量是 $e^T\bar I_Ce$。因此应依据关节轴选择分量，不能仅根据连杆在图中的横竖方向选择。

为避免“$I_1$”既表示论文标量又表示 MATLAB 的 6×6 矩阵，本节之后采用下列区分：

| 符号或变量 | 含义 | 尺寸 |
|---|---|---|
| $I_i^{\mathrm{paper}}=I_{zz,i}$ | 质心处、平行于关节轴的转动惯量 | 标量 |
| $\bar I_{C_i}$ / `mat_Ic(:,:,i)` | 质心转动惯量张量 | 3×3 |
| $J_i=I_{zz,i}+m_i l_{c_i}^2$ | 关于关节原点的 Z 轴转动惯量 | 标量 |
| $\mathcal I_i$ / `mat_I(:,:,i)` / `I_1`、`I_2` | 关于关节原点的空间惯量 | 6×6 |
| $M(q)$ | 整个二连杆的关节空间惯性 | 2×2 |

这里 $l_{c_i}=r_i$，与第 2–6 节的质心距离记号一致。**论文的 $I_i^{\mathrm{paper}}$ 不包含 $m_i l_{c_i}^2$；平行轴项在构造 $\mathcal I_i$ 或组合论文参数时加入一次。**

### 7.2 从质心惯量构造关节原点处的空间惯量

原书是：

$$
\mathsf H_i=\mathcal I_i v_i\qquad\text{(2.61)}
$$

其中 $\mathsf H_i$ 表示空间动量。若坐标原点设在质心 $C_i$，那么：

$$
\mathcal I_{C_i}=\begin{bmatrix}\bar I_{C_i}&0\\0&m_i\mathbf 1_3\end{bmatrix}\qquad\text{(2.62)}.
$$

但本例坐标系原点设在**关节**，与质心不重合，因此要用原书 **(2.63)**：

$$
\boxed{\mathcal I_i=\begin{bmatrix}
\bar I_{C_i}+m_i[r_{C_i}]_\times[r_{C_i}]_\times^T&m_i[r_{C_i}]_\times\\
m_i[r_{C_i}]_\times^T&m_i\mathbf 1_3
\end{bmatrix}.}
$$

这一步的物理原因是：关节坐标原点线速度为零，不代表远离关节的质心线速度为零。$\omega\times r_C$ 会产生质心运动，进而产生线动量和关于关节原点的角动量。

本例假设：

$$
r_{C_1}=\begin{bmatrix}r_1\\0\\0\end{bmatrix},\quad
r_{C_2}=\begin{bmatrix}r_2\\0\\0\end{bmatrix},\quad
h_1=m_1r_1,\quad h_2=m_2r_2.
$$

因为 $r_{C_i}=[r_i,0,0]^T$，有

$$
[r_{C_i}]_\times=\begin{bmatrix}0&0&0\\0&0&-r_i\\0&r_i&0\end{bmatrix},\qquad
m_i[r_{C_i}]_\times[r_{C_i}]_\times^T=\operatorname{diag}(0,m_ir_i^2,m_ir_i^2).
$$

**保留全部质心惯量分量时**，完整六维矩阵是

$$
\boxed{\mathcal I_i=\begin{bmatrix}
I_{xx,i}&I_{xy,i}&I_{xz,i}&0&0&0\\
I_{xy,i}&I_{yy,i}+m_ir_i^2&I_{yz,i}&0&0&-h_i\\
I_{xz,i}&I_{yz,i}&I_{zz,i}+m_ir_i^2&0&h_i&0\\
0&0&0&m_i&0&0\\
0&0&h_i&0&m_i&0\\
0&-h_i&0&0&0&m_i
\end{bmatrix},\qquad h_i=m_ir_i.}
$$

这与当前代码调用 `cross_3x1(mat_c(:,i))` 后构造的 `mat_I(:,:,i)` 一致。右上块为 $m_i[r_{C_i}]_\times$，左下块为它的转置，保证空间惯量对称；也可对照 [Featherstone 的 `mcI` 实现](https://royfeatherstone.org/spatial/v2/sourceText/mcI.txt)。

### 7.3 后续逐分量展开采用的简化

为保持第 8–11 节的空间力展开清晰，下面额外假设连杆局部轴为惯性主轴，即惯量积为零。此时：

$$
J_{O_i}=\operatorname{diag}(J_{xi},J_{yi},J_i),\qquad
J_{xi}=I_{xx,i},\quad J_{yi}=I_{yy,i}+m_ir_i^2,\quad
J_i=I_{zz,i}+m_ir_i^2=I_i^{\mathrm{paper}}+m_i l_{c_i}^2.
$$

最终得到完整六维矩阵：

$$
\boxed{\mathcal I_1=\begin{bmatrix}
J_{x1}&0&0&0&0&0\\
0&J_{y1}&0&0&0&-h_1\\
0&0&J_1&0&h_1&0\\
0&0&0&m_1&0&0\\
0&0&h_1&0&m_1&0\\
0&-h_1&0&0&0&m_1
\end{bmatrix},}
$$

$$
\boxed{\mathcal I_2=\begin{bmatrix}
J_{x2}&0&0&0&0&0\\
0&J_{y2}&0&0&0&-h_2\\
0&0&J_2&0&h_2&0\\
0&0&0&m_2&0&0\\
0&0&h_2&0&m_2&0\\
0&-h_2&0&0&0&m_2
\end{bmatrix}.}
$$

**务必注意**：上面稀疏结构额外使用了“质心在局部 X 轴、惯量积为零”的假设。论文只提供平面模型所需的 $I_{zz,i}$，不能据此把 $I_{xx,i}$、$I_{yy,i}$ 设为相同数值，也不能由该表断言惯量积为零。若保留一般惯量积，空间力的 X、Y 力矩分量可能非零；对于本文固定平面 2R 运动，最终 Z 轴关节力矩仍只依赖 $I_{zz,1}$、$I_{zz,2}$。第 12 节用一般对称质心惯量验证这一点。

---

## 8. Table 5.1：本体空间力 $f_{1}^{B}$、$f_{2}^{B}$

**原书位置**：§5.3 式 **(5.9)**，正文 p.93 / PDF p.100；Table 5.1，正文 p.96 / PDF p.103；力叉乘 (2.32)，正文 p.25 / PDF p.32。

以下逐分量表达式使用 §7.3 的惯量积为零假设；其中始终有 $J_i=I_{zz,i}+m_i r_i^2$，而非 $J_i=I_i^{\mathrm{paper}}$。

Table 5.1 每次前向循环计算：

$$
\boxed{f_i^B=\mathcal I_i a_i+\operatorname{crf}(v_i)\mathcal I_i v_i}
$$

其中

$$
\operatorname{crf}(v)=-\operatorname{crm}(v)^T
=\begin{bmatrix}[\omega]_\times&[u]_\times\\0&[\omega]_\times\end{bmatrix}.
$$

表中写的 $f_i$ 在循环初始化阶段等于 $f_i^B$（没有外力），反向递推时再累加子连杆传递载荷。不要把两个量混为一谈。

### 8.1 第一连杆 $f_1^B$：两个 6×6 乘法

第一部分：

$$
\mathcal I_1a_1=\begin{bmatrix}
0\\0\\J_1\ddot q_1+h_1gc_1\\m_1gs_1\\h_1\ddot q_1+m_1gc_1\\0
\end{bmatrix}.
$$

第二部分：

$$
\mathcal I_1v_1=\begin{bmatrix}0\\0\\J_1\dot q_1\\0\\h_1\dot q_1\\0\end{bmatrix}.
$$

再计算力叉乘：

$$
\operatorname{crf}(v_1)\mathcal I_1v_1
=\begin{bmatrix}0\\0\\0\\-h_1\dot q_1^2\\0\\0\end{bmatrix}.
$$

得到：

$$
\boxed{f_1^B=\begin{bmatrix}
0\\0\\J_1\ddot q_1+h_1gc_1\\m_1gs_1-h_1\dot q_1^2\\h_1\ddot q_1+m_1gc_1\\0
\end{bmatrix}.}
$$

### 8.2 第二连杆 $f_2^B$：输出完整的 6 个分量

把 $v_2,a_2,\mathcal I_2$ 代入 (5.9)，可得：

$$
\boxed{f_2^B=\begin{bmatrix}0\\0\\n_{2z}\\F_{2x}\\F_{2y}\\0\end{bmatrix},}
$$

其中：

$$
\boxed{\begin{aligned}
n_{2z}={}&J_2(\ddot q_1+\ddot q_2)+h_2l_1c_2\ddot q_1\\
&+h_2l_1s_2\dot q_1^2+h_2gc_{12},
\end{aligned}}
$$

$$
\boxed{\begin{aligned}
F_{2x}={}&m_2l_1s_2\ddot q_1-m_2l_1c_2\dot q_1^2\\
&+m_2gs_{12}-h_2(\dot q_1+\dot q_2)^2,
\end{aligned}}
$$

$$
\boxed{\begin{aligned}
F_{2y}={}&h_2(\ddot q_1+\ddot q_2)+m_2l_1c_2\ddot q_1\\
&+m_2l_1s_2\dot q_1^2+m_2gc_{12}.
\end{aligned}}
$$

空间力的第 3 行是绕 Z 轴的力矩；第 4、5 行是局部 X、Y 方向的力。这些分量都是 6D 计算后的结果，不是先把空间向量降成三维。

---

## 9. Table 5.1：第二次循环，`i=2→1`

**原书位置**：关节力矩式 **(5.11)**（正文 p.94 / PDF p.101），反向递推 **(5.20)**（正文 p.95 / PDF p.102），Table 5.1（正文 p.96 / PDF p.103）。

### 9.1 `i=2`：先计算第二关节力矩 $\tau_2$

第二连杆是末端，无子连杆，故

$$f_2=f_2^B.$$

Table 5.1 计算

$$
\tau_2=S_2^Tf_2
=\begin{bmatrix}0&0&1&0&0&0\end{bmatrix}
\begin{bmatrix}0\\0\\n_{2z}\\F_{2x}\\F_{2y}\\0\end{bmatrix}.
$$

因此

$$
\boxed{\begin{aligned}
\tau_2={}&(J_2+l_1h_2c_2)\ddot q_1+J_2\ddot q_2\\
&+l_1h_2s_2\dot q_1^2+h_2gc_{12}.
\end{aligned}}
$$

### 9.2 `i=2`：向父连杆传递空间力

书中 Table 5.1 的力变换用与速度变换对偶的坐标变换；在我们的列向量约定下，它具体等于 $X_2^T$：

$$
\boxed{f_1=f_1^B+X_2^Tf_2.}
$$

前面已知

$$
X_2^T=\begin{bmatrix}
c_2&-s_2&0&0&0&0\\
s_2&c_2&0&0&0&-l_1\\
0&0&1&l_1s_2&l_1c_2&0\\
0&0&0&c_2&-s_2&0\\
0&0&0&s_2&c_2&0\\
0&0&0&0&0&1
\end{bmatrix}.
$$

因此

$$
\boxed{X_2^Tf_2=\begin{bmatrix}
0\\0\\n_{2z}+l_1s_2F_{2x}+l_1c_2F_{2y}\\
c_2F_{2x}-s_2F_{2y}\\
s_2F_{2x}+c_2F_{2y}\\0
\end{bmatrix}.}
$$

为什么第 3 行不是只有 $n_{2z}$？因为第二连杆对父连杆除了传递力矩，还传递力；力通过 $O_1O_2=l_1$ 的力臂也会形成绕 Z 轴的力矩。

### 9.3 `i=1`：计算第一关节力矩 $\tau_1$

$$
\tau_1=S_1^Tf_1
=\begin{bmatrix}0&0&1&0&0&0\end{bmatrix}(f_1^B+X_2^Tf_2).
$$

因此

$$
\boxed{\tau_1=J_1\ddot q_1+h_1gc_1+n_{2z}+l_1s_2F_{2x}+l_1c_2F_{2y}.}
$$

将上一节的 $n_{2z},F_{2x},F_{2y}$ 全部代入，按 $\ddot q$、$\dot q$、$g$ 分组：

$$
\boxed{\begin{aligned}
\tau_1={}&(J_1+J_2+m_2l_1^2+2l_1h_2c_2)\ddot q_1\\
&+(J_2+l_1h_2c_2)\ddot q_2\\
&-l_1h_2s_2(2\dot q_1\dot q_2+\dot q_2^2)\\
&+g[(h_1+m_2l_1)c_1+h_2c_{12}].
\end{aligned}}
$$

这就是完整的六维空间向量 RNEA 对平面 2R 示例计算出的第一关节力矩。

---

## 10. 最终 $\tau$、$M(q)$、$C(q,\dot q)$、$G(q)$

### 10.1 按第 2–6 节水平零位约定整理 RNEA

把两关节结果组合：

$$
\boxed{\tau=\begin{bmatrix}\tau_1\\\tau_2\end{bmatrix}
=M(q)\ddot q+C(q,\dot q)\dot q+G(q).}
$$

其中：

$$
\boxed{M(q)=\begin{bmatrix}
J_1+J_2+m_2l_1^2+2l_1h_2c_2&J_2+l_1h_2c_2\\
J_2+l_1h_2c_2&J_2
\end{bmatrix}.}
$$

$$
\boxed{C(q,\dot q)=l_1h_2s_2
\begin{bmatrix}-\dot q_2&-(\dot q_1+\dot q_2)\\\dot q_1&0\end{bmatrix}.}
$$

$$
\boxed{G(q)=g\begin{bmatrix}
(h_1+m_2l_1)c_1+h_2c_{12}\\h_2c_{12}
\end{bmatrix}.}
$$

这部分矩阵分组是**对上述 RNEA 结果的代数整理**，并不是说 Table 5.1 必须先求 $M,C,G$ 才能计算力矩。恰恰相反：RNEA 直接计算 $\tau$，不需要显式构造这三个矩阵。

### 10.2 与论文式 (64) 的参数逐项对应

将 $r_i=l_{c_i}$、$h_i=m_i l_{c_i}$、$J_i=I_{zz,i}+m_i l_{c_i}^2$ 代入，得到论文式 (67) 的参数组合：

$$
\boxed{\begin{aligned}
\theta_1&=J_1+J_2+m_2l_1^2
=I_{zz,1}+m_1l_{c_1}^2+I_{zz,2}+m_2l_{c_2}^2+m_2l_1^2,\\
\theta_2&=l_1h_2=l_1m_2l_{c_2},\\
\theta_3&=J_2=I_{zz,2}+m_2l_{c_2}^2,\\
\theta_4&=h_1+m_2l_1=m_1l_{c_1}+m_2l_1,\\
\theta_5&=h_2=m_2l_{c_2},\\
\theta_6&=f_{v_1},\qquad \theta_7=f_{v_2}.
\end{aligned}}
$$

因此 RNEA 的惯性矩阵恰好变为

$$
\boxed{M=\begin{bmatrix}
\theta_1+2\theta_2\cos q_2&\theta_3+\theta_2\cos q_2\\
\theta_3+\theta_2\cos q_2&\theta_3
\end{bmatrix},}
$$

与[论文 §5 式 (64)、(67)](https://arxiv.org/html/2609.27010v1#S5)一致。这里再次可见：若把论文的 $I_2$ 当成关于关节原点的 $J_2$ 再加入 $m_2l_{c_2}^2$，就会重复计算平行轴项。

### 10.3 论文式 (65) 的重力项为什么是正弦？

第 2–6 节的 $q_1$ 是从世界水平 $+x_0$ 方向量起的角；论文图 1 与势能式 (65) 对应从竖直向下方向量起的第一关节角。为了保留前六节的定义，用 $Q_1,Q_2$ 标记论文中的角度：

$$
\boxed{q_1=Q_1-\frac\pi2,\qquad q_2=Q_2.}
$$

这是恒定零位偏移，所以 $\dot q_i=\dot Q_i$、$\ddot q_i=\ddot Q_i$，$M$、$C$ 的形式不变。本文水平零位下的势能（忽略加性常数）为

$$
U_{\mathrm{doc}}=g\theta_4\sin q_1+g\theta_5\sin(q_1+q_2).
$$

代入零位转换后

$$
\boxed{U_{\mathrm{paper}}=-g\theta_4\cos Q_1-g\theta_5\cos(Q_1+Q_2),}
$$

这就是截图中的式 (65)。对论文角度求偏导，得到

$$
\boxed{G_{\mathrm{paper}}=\frac{\partial U_{\mathrm{paper}}}{\partial Q}
=g\begin{bmatrix}
\theta_4\sin Q_1+\theta_5\sin(Q_1+Q_2)\\
\theta_5\sin(Q_1+Q_2)
\end{bmatrix}.}
$$

所以文档中的余弦重力项与论文的正弦重力项是同一物理模型在不同零位下的表达。若直接把论文角度输入当前 MDH 表，需要把第一行的转角写为 `Q1 - sym(pi)/2`，第二行仍为 `Q2`。

### 10.4 加入论文式 (66) 的黏性摩擦

前文 RNEA 求的是无关节摩擦的刚体力矩。论文还包含 Rayleigh 耗散函数

$$
F_R=\frac12\theta_6\dot Q_1^2+\frac12\theta_7\dot Q_2^2,
\qquad
\tau_{\mathrm{fric}}=\frac{\partial F_R}{\partial\dot Q}
=\begin{bmatrix}\theta_6\dot Q_1\\\theta_7\dot Q_2\end{bmatrix}.
$$

因此与论文完整模型对应的驱动力矩为

$$
\boxed{\begin{aligned}
\tau_1={}&(\theta_1+2\theta_2\cos Q_2)\ddot Q_1
+(\theta_3+\theta_2\cos Q_2)\ddot Q_2\\
&-\theta_2\sin Q_2(2\dot Q_1\dot Q_2+\dot Q_2^2)\\
&+g[\theta_4\sin Q_1+\theta_5\sin(Q_1+Q_2)]
+\theta_6\dot Q_1,\\
\tau_2={}&(\theta_3+\theta_2\cos Q_2)\ddot Q_1+\theta_3\ddot Q_2
+\theta_2\sin Q_2\dot Q_1^2\\
&+g\theta_5\sin(Q_1+Q_2)+\theta_7\dot Q_2.
\end{aligned}}
$$

摩擦项在关节力矩层加入；它与前向循环的本体空间力 $f_i^B$ 是不同的计算环节。

---

## 11. 代入数值，把每个中间变量全部算出来

§11.1–11.8 保留便于手算的自选细杆算例，**不是论文实物参数**；§11.9 使用论文 Table 1 的数据进行对照。

这一节的目的就是回答“**S、vJ、cJ 等到底是多少**”。除了给出通式，我们再选择一个运动状态，把所有六维变量具体算成数字。

### 11.1 选择一个非静止状态

为保证三角函数也能简单手算，取：

| 参数 | 选定值 | 单位或意义 |
|---|---:|---|
| $q_1$ | $0$ | rad |
| $q_2$ | $\pi/2$ | rad |
| $\dot q_1$ | $1$ | rad/s |
| $\dot q_2$ | $2$ | rad/s |
| $\ddot q_1$ | $1$ | rad/s² |
| $\ddot q_2$ | $-1$ | rad/s² |
| $l_1$ | $1$ | m |
| $l_2$ | $0.5$ | m |
| $m_1$ | $1$ | kg |
| $m_2$ | $2$ | kg |
| $r_1$ | $0.5$ | m |
| $r_2$ | $0.25$ | m |
| $g$ | $10$ | m/s² |
| $J_1$ | $1/3$ | kg·m²，关于第一关节 |
| $J_2$ | $1/6$ | kg·m²，关于第二关节 |

相应的质心 Z 轴惯量是 $I_{zz,1}=J_1-m_1r_1^2=1/12$、$I_{zz,2}=J_2-m_2r_2^2=1/24$，而非 $1/3$、$1/6$。这里文档角度 $q_1=0$ 对应论文角度 $Q_1=\pi/2$。

以理想均匀细杆为例，其他惯量可选 $J_{x1}=J_{x2}=0$，$J_{y1}=1/3$，$J_{y2}=1/6$；轴向转动惯量取零是细杆零半径的理想化近似，**不是对一般机器人零惯量的假设**。这里计算关节 Z 轴力矩时，$J_{x_i},J_{y_i}$ 不影响最终值。

$$
h_1=m_1r_1=\frac12,\qquad h_2=m_2r_2=\frac12,
$$

$$c_1=1,\quad s_1=0,\quad c_2=0,\quad s_2=1,\quad c_{12}=0,\quad s_{12}=1.$$

### 11.2 关节运动子空间和关节速度：显式数字

$$
\boxed{S_1=S_2=\begin{bmatrix}0\\0\\1\\0\\0\\0\end{bmatrix}}
$$

$$
\boxed{v_{J1}=\begin{bmatrix}0\\0\\1\\0\\0\\0\end{bmatrix}\ \mathrm{(SI)},\qquad
v_{J2}=\begin{bmatrix}0\\0\\2\\0\\0\\0\end{bmatrix}\ \mathrm{(SI)}.}
$$

$$
\boxed{c_{J1}=c_{J2}=\begin{bmatrix}0\\0\\0\\0\\0\\0\end{bmatrix}.}
$$

> 六维空间向量的前三项是角量，后三项是线量，严格说整个六维向量的各分量具有不同单位；因此上面的“SI”只说明各分量分别采用 SI 单位，而不是整体一个统一单位。

### 11.3 数值空间变换

$$
\boxed{X_1=I_6.}
$$

由 $q_2=\pi/2$、$l_1=1$ 得：

$$
\boxed{X_2=\begin{bmatrix}
0&1&0&0&0&0\\
-1&0&0&0&0&0\\
0&0&1&0&0&0\\
0&0&1&0&1&0\\
0&0&0&-1&0&0\\
0&-1&0&0&0&1
\end{bmatrix}.}
$$

### 11.3.1 数值 $X_{J2}$、$X_{T2}$ 与空间惯性矩阵

为了与 Table 5.1 中 `Xup = XJ*XT` 一一对应，还可以把第二关节两部分分别写成数字：

$$
\boxed{X_{J2}=\begin{bmatrix}
0&1&0&0&0&0\\
-1&0&0&0&0&0\\
0&0&1&0&0&0\\
0&0&0&0&1&0\\
0&0&0&-1&0&0\\
0&0&0&0&0&1
\end{bmatrix},\quad
X_{T2}=\begin{bmatrix}
1&0&0&0&0&0\\
0&1&0&0&0&0\\
0&0&1&0&0&0\\
0&0&0&1&0&0\\
0&0&1&0&1&0\\
0&-1&0&0&0&1
\end{bmatrix}.}
$$

对应的两根连杆空间惯性矩阵的**每个数值**是：

$$
\boxed{\mathcal I_1=\begin{bmatrix}
0&0&0&0&0&0\\
0&\frac13&0&0&0&-\frac12\\
0&0&\frac13&0&\frac12&0\\
0&0&0&1&0&0\\
0&0&\frac12&0&1&0\\
0&-\frac12&0&0&0&1
\end{bmatrix}.}
$$

$$
\boxed{\mathcal I_2=\begin{bmatrix}
0&0&0&0&0&0\\
0&\frac16&0&0&0&-\frac12\\
0&0&\frac16&0&\frac12&0\\
0&0&0&2&0&0\\
0&0&\frac12&0&2&0\\
0&-\frac12&0&0&0&2
\end{bmatrix}.}
$$

检查一下：**两根连杆质量不同，因此右下角的平动惯性块分别是 $1I_3$ 和 $2I_3$；它们没有被写成同一个质量。**

### 11.4 初始化后的 $v_0,a_0$

$$
\boxed{v_0=\begin{bmatrix}0\\0\\0\\0\\0\\0\end{bmatrix},\qquad
a_0=\begin{bmatrix}0\\0\\0\\0\\10\\0\end{bmatrix}.}
$$

### 11.5 第一次循环，$i=1$

$$
v_1=X_1v_0+v_{J1}
\quad\Rightarrow\quad
\boxed{v_1=\begin{bmatrix}0\\0\\1\\0\\0\\0\end{bmatrix}.}
$$

$$
a_1=X_1a_0+S_1\ddot q_1+c_{J1}+\operatorname{crm}(v_1)v_{J1}
\quad\Rightarrow\quad
\boxed{a_1=\begin{bmatrix}0\\0\\1\\0\\10\\0\end{bmatrix}.}
$$

本体空间力：

$$
\boxed{f_1^B=\begin{bmatrix}
0\\0\\\frac{16}{3}\\-\frac12\\\frac{21}{2}\\0
\end{bmatrix}.}
$$

### 11.6 第一次循环，$i=2$

$$
X_2v_1=\begin{bmatrix}0\\0\\1\\1\\0\\0\end{bmatrix},
\qquad v_{J2}=\begin{bmatrix}0\\0\\2\\0\\0\\0\end{bmatrix}.
$$

相加：

$$
\boxed{v_2=\begin{bmatrix}0\\0\\3\\1\\0\\0\end{bmatrix}.}
$$

此时：

$$
X_2a_1=\begin{bmatrix}0\\0\\1\\11\\0\\0\end{bmatrix},
\qquad S_2\ddot q_2=\begin{bmatrix}0\\0\\-1\\0\\0\\0\end{bmatrix},
$$

$$
\operatorname{crm}(v_2)v_{J2}=\begin{bmatrix}0\\0\\0\\0\\-2\\0\end{bmatrix}.
$$

三项相加（$c_{J2}=0$）：

$$
\boxed{a_2=\begin{bmatrix}0\\0\\0\\11\\-2\\0\end{bmatrix}.}
$$

本体空间力：

$$
\boxed{f_2^B=\begin{bmatrix}
0\\0\\\frac12\\\frac{35}{2}\\2\\0
\end{bmatrix}.}
$$

### 11.7 第二次循环，$i=2\to1$

关节 2：

$$
\boxed{\tau_2=S_2^Tf_2=\frac12\ \mathrm{N\,m}.}
$$

传给第一连杆的空间力：

$$
\boxed{X_2^Tf_2=\begin{bmatrix}
0\\0\\18\\-2\\\frac{35}{2}\\0
\end{bmatrix}.}
$$

累加第一连杆自身的本体空间力：

$$
\boxed{f_1=f_1^B+X_2^Tf_2=\begin{bmatrix}
0\\0\\\frac{70}{3}\\-\frac52\\28\\0
\end{bmatrix}.}
$$

最后：

$$
\boxed{\tau_1=S_1^Tf_1=\frac{70}{3}\ \mathrm{N\,m}\approx23.333\ \mathrm{N\,m}.}
$$

**数值结果**：

$$
\boxed{\tau=\begin{bmatrix}\frac{70}{3}\\\frac12\end{bmatrix}\ \mathrm{N\,m}.}
$$

### 11.8 把整张 Table 5.1 的数值过程放在一起

| 量 | 第一关节/连杆 $i=1$ | 第二关节/连杆 $i=2$ |
|---|---|---|
| $S_i$ | $[0,0,1,0,0,0]^T$ | $[0,0,1,0,0,0]^T$ |
| $v_{J_i}$ | $[0,0,1,0,0,0]^T$ | $[0,0,2,0,0,0]^T$ |
| $c_{J_i}$ | $[0,0,0,0,0,0]^T$ | $[0,0,0,0,0,0]^T$ |
| $v_i$ | $[0,0,1,0,0,0]^T$ | $[0,0,3,1,0,0]^T$ |
| $a_i$ | $[0,0,1,0,10,0]^T$ | $[0,0,0,11,-2,0]^T$ |
| $f_i^B$ | $[0,0,16/3,-1/2,21/2,0]^T$ | $[0,0,1/2,35/2,2,0]^T$ |
| $f_i$（反向递推后） | $[0,0,70/3,-5/2,28,0]^T$ | $[0,0,1/2,35/2,2,0]^T$ |
| $\tau_i$ | $70/3\ \mathrm{N\,m}$ | $1/2\ \mathrm{N\,m}$ |

---

### 11.9 使用论文 Table 1 的实物参数对照

以下数据来自[论文 §5 Table 1](https://arxiv.org/html/2609.27010v1#S5)，惯量在本建模约定下对应 $I_{zz}$：

| 参数 | 数值 | 参数 | 数值 |
|---|---:|---|---:|
| $l_1$（m） | 0.45 | $l_2$（m） | 0.45 |
| $l_{c_1}$（m） | 0.091 | $l_{c_2}$（m） | 0.048 |
| $m_1$（kg） | 23.902 | $m_2$（kg） | 3.88 |
| $I_{zz,1}$（kg·m²） | 1.266 | $I_{zz,2}$（kg·m²） | 0.093 |
| $f_{v_1}$（N·m·s/rad） | 2.288 | $f_{v_2}$（N·m·s/rad） | 0.175 |
| $g$（m/s²） | 9.81 | | |

由这些数值重新组合得到

$$
\begin{aligned}
J_1&=1.463932462,\quad J_2=0.10193952,\\
h_1&=2.175082,\quad h_2=0.18624,\\
\theta&=[2.351571982,\ 0.083808,\ 0.10193952,\ 3.921082,\ 0.18624,\ 2.288,\ 0.175]^T.
\end{aligned}
$$

论文 Table 2 以有限位数列出对应参数；这里保留由 Table 1 显示值直接计算的结果，不要求它与 Table 2 的有限位数值逐位一致。

取与前一算例相同的角速度、角加速度，并取论文角度 $Q_1=Q_2=\pi/2$，即文档角度 $q_1=0,q_2=\pi/2$。无摩擦刚体力矩为

$$
\boxed{\tau_{\mathrm{rigid}}=[40.044982882,\ 0.083808]^T\ \mathrm{N\,m}.}
$$

加上黏性摩擦 $[2.288\times1,\ 0.175\times2]^T$ 后为

$$
\boxed{\tau_{\mathrm{paper}}=[42.332982882,\ 0.433808]^T\ \mathrm{N\,m}.}
$$

论文没有在该表中提供 $I_{xx}$、$I_{yy}$ 或惯量积，因此无法仅凭这张表唯一恢复两个完整的 3×3 质心惯量张量；但这不妨碍计算平面 2R 的关节力矩。

---

## 12. 一份可独立运行的 MATLAB 符号复现程序

将下面代码保存为 `demo_2r_table51.m` 并运行，需要 MATLAB Symbolic Math Toolbox。程序从一般对称的 3×3 质心惯量张量构造 6×6 空间惯量，全程执行六维 RNEA；保留惯量积，避免把论文未提供的惯量分量擅自设为相同数值。

程序包含五类校验：水平零位的解析力矩、最终力矩对其他惯量分量的零依赖、论文零位下的 M/C/G/摩擦形式、独立拉格朗日推导，以及自选细杆算例和论文实物参数算例。数值状态的第一关节角始终以文档零位输入，`q1=0` 对应论文 `Q1=pi/2`。

```matlab
clear; clc;

%% q1 是文档的水平零位角；Q1 是论文的竖直向下零位角
syms q1 q2 Q1 Q2 dq1 dq2 ddq1 ddq2 l1 g real
syms m1 m2 lc1 lc2 b1 b2 real
syms Ixx1 Ixy1 Ixz1 Iyy1 Iyz1 Izz1 real
syms Ixx2 Ixy2 Ixz2 Iyy2 Iyz2 Izz2 real

skew = @(r) [0,-r(3),r(2); r(3),0,-r(1); -r(2),r(1),0];
Rz = @(q) [cos(q),-sin(q),0; sin(q),cos(q),0; 0,0,1];
crm = @(v) [skew(v(1:3)),sym(zeros(3));
            skew(v(4:6)),skew(v(1:3))];
crf = @(v) -crm(v).';
spatialI = @(m,c,Ic) [Ic+m*skew(c)*skew(c).',m*skew(c);
                     m*skew(c).',m*sym(eye(3))];

%% 质心张量的 Izz 对应论文标量惯量；I_1、I_2 为 6x6 矩阵
Ic1 = [Ixx1,Ixy1,Ixz1; Ixy1,Iyy1,Iyz1; Ixz1,Iyz1,Izz1];
Ic2 = [Ixx2,Ixy2,Ixz2; Ixy2,Iyy2,Iyz2; Ixz2,Iyz2,Izz2];
I_1 = spatialI(m1,[lc1;0;0],Ic1);
I_2 = spatialI(m2,[lc2;0;0],Ic2);
mat_I = {I_1,I_2};

%% 关节运动学和空间变换
S = sym([0;0;1;0;0;0]);
mat_vJ = [S*dq1,S*dq2];
mat_cJ = sym(zeros(6,2));
E1 = Rz(q1).'; E2 = Rz(q2).';
XJ1 = blkdiag(E1,E1); XJ2 = blkdiag(E2,E2);
XT1 = sym(eye(6));
XT2 = [sym(eye(3)),sym(zeros(3)); -skew([l1;0;0]),sym(eye(3))];
mat_X = {XJ1*XT1,XJ2*XT2};
v0 = sym(zeros(6,1));
a0 = sym([0;0;0;0;g;0]);
dq = [dq1;dq2]; ddq = [ddq1;ddq2];

%% 6D RNEA 前向循环
mat_v = cell(2,1); mat_a = cell(2,1); mat_fB = cell(2,1);
for i = 1:2
    if i == 1
        vp = v0; ap = a0;
    else
        vp = mat_v{i-1}; ap = mat_a{i-1};
    end
    mat_v{i} = simplify(mat_X{i}*vp+mat_vJ(:,i));
    mat_a{i} = simplify(mat_X{i}*ap+S*ddq(i)+mat_cJ(:,i) ...
                       +crm(mat_v{i})*mat_vJ(:,i));
    mat_fB{i} = simplify(mat_I{i}*mat_a{i} ...
                        +crf(mat_v{i})*(mat_I{i}*mat_v{i}));
end

%% 反向循环：本体空间力累加为子树空间力，再投影关节力矩
mat_f = mat_fB;
tauRigid = sym(zeros(2,1));
for i = 2:-1:1
    tauRigid(i) = S.'*mat_f{i};
    if i > 1
        mat_f{i-1} = simplify(mat_f{i-1}+mat_X{i}.'*mat_f{i});
    end
end
tauRigid = simplify(tauRigid);

%% 检验水平零位的解析式及所依赖的惯量分量
h1 = m1*lc1; h2 = m2*lc2;
J1 = Izz1+m1*lc1^2; J2 = Izz2+m2*lc2^2;
M = [J1+J2+m2*l1^2+2*l1*h2*cos(q2),J2+l1*h2*cos(q2);
     J2+l1*h2*cos(q2),J2];
C = l1*h2*sin(q2)*[-dq2,-(dq1+dq2);dq1,0];
Gdoc = g*[(h1+m2*l1)*cos(q1)+h2*cos(q1+q2);h2*cos(q1+q2)];
assert(isequal(simplify(tauRigid-M*ddq-C*dq-Gdoc),sym(zeros(2,1))));
assert(isequal(simplify(jacobian(tauRigid,ddq)-M),sym(zeros(2))));
otherInertias = [Ixx1,Ixy1,Ixz1,Iyy1,Iyz1,Ixx2,Ixy2,Ixz2,Iyy2,Iyz2];
assert(isequal(simplify(jacobian(tauRigid,otherInertias)),sym(zeros(2,10))));

%% 论文式 (64)-(67)：零位转换与黏性摩擦
theta = [J1+J2+m2*l1^2;l1*h2;J2;h1+m2*l1;h2;b1;b2];
Mpaper = [theta(1)+2*theta(2)*cos(Q2),theta(3)+theta(2)*cos(Q2);
          theta(3)+theta(2)*cos(Q2),theta(3)];
Cpaper = theta(2)*sin(Q2)*[-dq2,-(dq1+dq2);dq1,0];
Upaper = -g*theta(4)*cos(Q1)-g*theta(5)*cos(Q1+Q2);
Gpaper = jacobian(Upaper,[Q1,Q2]).';
friction = [b1*dq1;b2*dq2];
tauPaper = simplify(subs(tauRigid,[q1,q2],[Q1-sym(pi)/2,Q2])+friction);
assert(isequal(simplify(tauPaper-Mpaper*ddq-Cpaper*dq-Gpaper-friction), ...
               sym(zeros(2,1))));

%% 独立拉格朗日校验：质心速度和质心 Izz 构造动能
pC1 = lc1*[sin(Q1);-cos(Q1);0];
pC2 = l1*[sin(Q1);-cos(Q1);0]+lc2*[sin(Q1+Q2);-cos(Q1+Q2);0];
uC1 = jacobian(pC1,[Q1,Q2])*dq;
uC2 = jacobian(pC2,[Q1,Q2])*dq;
K = (m1*(uC1.'*uC1)+m2*(uC2.'*uC2) ...
     +Izz1*dq1^2+Izz2*(dq1+dq2)^2)/2;
U = g*(m1*pC1(2)+m2*pC2(2));
p = jacobian(K,dq).';
tauLagrange = jacobian(p,[Q1,Q2])*dq+jacobian(p,dq)*ddq ...
              -jacobian(K,[Q1,Q2]).'+jacobian(U,[Q1,Q2]).'+friction;
assert(isequal(simplify(U-Upaper),sym(0)));
assert(isequal(simplify(tauLagrange-tauPaper),sym(zeros(2,1))));

%% 自选细杆算例：精确有理数，惯量积为零
vars = [q1,q2,dq1,dq2,ddq1,ddq2,l1,g,m1,m2,lc1,lc2, ...
        Ixx1,Ixy1,Ixz1,Iyy1,Iyz1,Izz1,Ixx2,Ixy2,Ixz2,Iyy2,Iyz2,Izz2];
vals = [0,sym(pi)/2,1,2,1,-1,1,10,1,2,sym(1)/2,sym(1)/4, ...
        0,0,0,sym(1)/12,0,sym(1)/12,0,0,0,sym(1)/24,0,sym(1)/24];
expectedVB = sym([0;0;1;0;0;0;0;0;3;1;0;0]);
expectedAB = sym([0;0;1;0;10;0;0;0;0;11;-2;0]);
expectedFB = [0;0;sym(16)/3;-sym(1)/2;sym(21)/2;0; ...
              0;0;sym(1)/2;sym(35)/2;2;0];
assert(isequal(simplify(subs([mat_v{1};mat_v{2}],vars,vals)),expectedVB));
assert(isequal(simplify(subs([mat_a{1};mat_a{2}],vars,vals)),expectedAB));
assert(isequal(simplify(subs([mat_fB{1};mat_fB{2}],vars,vals)),expectedFB));
assert(isequal(simplify(subs(tauRigid,vars,vals)),[sym(70)/3;sym(1)/2]));

names = {'S1','S2','vJ1','vJ2','cJ1','cJ2','XJ1','XJ2','XT1','XT2', ...
         'X1','X2','I_1','I_2','v0','a0','v1','v2','a1','a2', ...
         'f1B','f2B','f1','f2','tauRigid'};
expressions = {S,S,mat_vJ(:,1),mat_vJ(:,2),mat_cJ(:,1),mat_cJ(:,2), ...
               XJ1,XJ2,XT1,XT2,mat_X{1},mat_X{2},I_1,I_2,v0,a0, ...
               mat_v{1},mat_v{2},mat_a{1},mat_a{2},mat_fB{1},mat_fB{2}, ...
               mat_f{1},mat_f{2},tauRigid};
for k = 1:numel(names)
    fprintf('%s =\n',names{k});
    disp(simplify(subs(expressions{k},vars,vals)));
end

%% 论文 Table 1 实物参数：不赋值未提供的 Ixx、Iyy、惯量积
paperVars = [q1,q2,dq1,dq2,ddq1,ddq2,l1,g,m1,m2,lc1,lc2,Izz1,Izz2,b1,b2];
paperVals = [0,sym(pi)/2,1,2,1,-1,sym('0.45'),sym('9.81'), ...
             sym('23.902'),sym('3.88'),sym('0.091'),sym('0.048'), ...
             sym('1.266'),sym('0.093'),sym('2.288'),sym('0.175')];
thetaNum = double(subs(theta,paperVars,paperVals));
rigidNum = double(subs(tauRigid,paperVars,paperVals));
totalNum = double(subs(tauRigid+friction,paperVars,paperVals));
assert(norm(rigidNum-[40.044982882;0.083808],inf)<1e-9);
assert(norm(totalNum-[42.332982882;0.433808],inf)<1e-9);
fprintf('Paper theta =\n'); disp(thetaNum);
fprintf('Paper rigid-body torque =\n'); disp(rigidNum);
fprintf('Paper total torque with friction =\n'); disp(totalNum);
disp('PASS: 6D RNEA, inertia components, paper model, Lagrange model, both examples.');
```

**数值处理**：手算示例使用符号 pi 和有理数；论文 Table 1 的小数通过 `sym('数值字符串')` 输入。表中有限精度数据不能代表实物参数的无限精度真值，结果按其显示值计算。未提供的 Ixx、Iyy 和惯量积保留为符号，最终关节力矩中它们会消去。

**验证记录**：本节代码已在当前项目的 MATLAB 环境中实际运行通过，确认六维 RNEA、论文参数化表达和独立拉格朗日推导一致，两组数值算例也均通过校验。

---

## 13. 与 Table 5.1 的逐行核对表

| 原书 Table 5.1 的行 | 本例实际内容 | 结果在哪里 |
|---|---|---|
| `v0=0` | 固定基座六维速度零向量 | §4 |
| `a0=-ag` | $[0,0,0,0,g,0]^T$ | §4 |
| `jcalc(...)` | 两个旋转关节的 $XJ,S,vJ,cJ$ | §3 |
| `Xup = XJ*XT` | $X_1=X_{J1}$；$X_2=X_{J2}X_{T2}$ | §3.4 |
| `v = Xup*v_parent + vJ` | 先求 $v_1$，再求 $v_2$ | §5.1、§6.1 |
| `a = Xup*a_parent + S*ddq + cJ + v×vJ` | 先求 $a_1$，再求 $a_2$，叉乘没有省略 | §5.2、§6.2 |
| `f = I*a + v×*I*v` | $f_1^B,f_2^B$ | §7–8 |
| `tau = S^T*f` | 先求 $\tau_2$，再求 $\tau_1$ | §9 |
| `f_parent += Xup^T*f` | $f_1=f_1^B+X_2^Tf_2$ | §9.2 |
| **输出 `tau`** | $[\tau_1,\tau_2]^T$ | §10、§11 |

### 与论文对应时需要区分的量

1. **$v_{J_i}$ 不是 $v_i$。** $v_{J_i}$ 只表示这个关节自身的相对运动；$v_i=X_iv_{i-1}+v_{J_i}$ 才是完整的连杆空间速度。
2. **$f_i^B$ 不是最终的 $f_i$。** 前者仅对应当前刚体；反向递推必须把所有下游连杆传来的空间力加进去。
3. **$\mathcal I_i$ 不等于机器人关节空间惯性矩阵 $M(q)$。** $\mathcal I_i$ 是单根刚体的 6×6 空间惯性，$M(q)$ 是整个 2R 系统的 2×2 关节空间惯性。
4. **论文的 $I_i^{\mathrm{paper}}=I_{zz,i}$ 是质心标量惯量。** 它不等于关节原点惯量 $J_i=I_{zz,i}+m_i l_{c_i}^2$，也不等于代码中的 6×6 `I_1`、`I_2`。
5. **文档与论文的第一关节零位不同。** $q_1=Q_1-\pi/2$，故水平零位的余弦重力项转换为论文的正弦重力项；惯性矩阵不受该常量偏移影响。
6. **论文还包含黏性摩擦。** 完整驱动力矩为 $\tau_{\mathrm{rigid}}+[\theta_6\dot Q_1,\theta_7\dot Q_2]^T$。

**文献来源**：Roy Featherstone, *Rigid Body Dynamics Algorithms*, Springer, 2008，重点参见 §2.8、§2.9、§2.13、§4.4、§5.3。原书页码索引沿用本文件前六节；空间惯量结构另对照 [Featherstone `mcI` 官方实现](https://royfeatherstone.org/spatial/v2/sourceText/mcI.txt)。论文对照依据 Luis Cervantes-Pérez 等，[*Three High Performance Global Tracking Composite Adaptive Controllers for Fully Actuated Euler–Lagrange Systems: Experimental Validation*，arXiv:2609.27010v1](https://arxiv.org/html/2609.27010v1#S5)，§5 图 1、Table 1–2、式 (64)–(67)。$I_i^{\mathrm{paper}}$ 与 $I_{zz,i}$ 的对应及两套角度的转换依据本文建模约定推导。
