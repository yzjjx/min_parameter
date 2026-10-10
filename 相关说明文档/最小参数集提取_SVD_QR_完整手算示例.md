# 从采样力矩中提取机器人动力学最小参数集：SVD + 两次 QR 完整手算

> 对应你给出的 MATLAB 代码：从 `tau_samples` 开始，依次构造 `theta_full`、`W_symbolic`，使用 SVD 判断秩，使用两次带列主元 QR 选择独立列和独立行，最后通过 `parameter_map` 得到 `pi_param`。
>
> **本文用一个 2 关节、2 次采样、3 个原始参数单项式、2 个独立参数组合的例子逐步推导。** 例子特意设计得易于手算，重点是演示矩阵算法，不代表某台具体机器人的完整 Newton–Euler 动力学模型。

---

## 0. 先给出本例最终要证明的结论

设三个原始动力学参数为 $h_1,h_2,h_3$（可理解为具有相同量纲的惯性类参数或参数单项式），机器人有两个关节，定义 $a=\ddot q_1$、$b=\ddot q_2$。构造一个简单的力矩方程：

$$
\boxed{\begin{aligned}
\tau_1&=a h_1+b h_2+(a+b)h_3,\\
\tau_2&=b h_1+a h_2+(a+b)h_3.
\end{aligned}}
$$

乍看有三个未知量 $h_1,h_2,h_3$，但是展开重组：

$$
\begin{aligned}
\tau_1&=a(h_1+h_3)+b(h_2+h_3),\\
\tau_2&=b(h_1+h_3)+a(h_2+h_3).
\end{aligned}
$$

因此只需要两个固定参数组合：

$$
\boxed{\pi_1=h_1+h_3,\qquad\pi_2=h_2+h_3}
$$

并可以写成：

$$
\boxed{
\begin{bmatrix}\tau_1\\\tau_2\end{bmatrix}
=\underbrace{\begin{bmatrix}a&b\\b&a\end{bmatrix}}_{Y_{\mathrm{base}}(a,b)}
\underbrace{\begin{bmatrix}\pi_1\\\pi_2\end{bmatrix}}_{\pi}
}
$$

我们接下来**故意假装事先不知道** $\pi_1,\pi_2$ 是什么，只按照你的 MATLAB 代码，从两次采样的力矩公式中自动找出它们。

---

## 1. 准备输入：`tau_samples` 从哪里来？

对应代码：

```matlab
for sample_idx = 1:n_samples
    for joint_idx = 1:joint_count
        row_idx = (sample_idx-1)*joint_count + joint_idx;
        ...
    end
end
```

取：

```matlab
joint_count = 2;   % 两个关节
n_samples   = 2;   % 两次采样
```

对 $a=\ddot q_1,b=\ddot q_2$ 采样：

| 采样编号 | $a$ | $b$ | $\tau_1$ | $\tau_2$ |
|---|---:|---:|---|---|
| 1 | 1 | 0 | $h_1+h_3$ | $h_2+h_3$ |
| 2 | 0 | 1 | $h_2+h_3$ | $h_1+h_3$ |

所以 `tau_samples` 是一个 **2 行 × 2 列**的符号矩阵：

$$
\boxed{\texttt{tau\_samples}=\begin{bmatrix}
 h_1+h_3 & h_2+h_3\\
 h_2+h_3 & h_1+h_3
\end{bmatrix}}
$$

- **行**：关节编号，第一行表示关节 1，第二行表示关节 2。
- **列**：采样编号，第一列表示第 1 次采样，第二列表示第 2 次采样。
- 采样改变的是 $a,b$，**没有**给未知参数 $h_1,h_2,h_3$ 赋数值。否则后面无法提取它们的符号系数。

此时你的 `row_idx` 对应如下：

| `sample_idx` | `joint_idx` | `row_idx=(sample_idx-1)*2+joint_idx` | 对应力矩 |
|---:|---:|---:|---|
| 1 | 1 | 1 | $\tau_1^{(1)}=h_1+h_3$ |
| 1 | 2 | 2 | $\tau_2^{(1)}=h_2+h_3$ |
| 2 | 1 | 3 | $\tau_1^{(2)}=h_2+h_3$ |
| 2 | 2 | 4 | $\tau_2^{(2)}=h_1+h_3$ |

这样力矩按列堆叠后的向量为：

$$
T=\begin{bmatrix}
\tau_1^{(1)}\\\tau_2^{(1)}\\\tau_1^{(2)}\\\tau_2^{(2)}
\end{bmatrix}
=\begin{bmatrix}
h_1+h_3\\h_2+h_3\\h_2+h_3\\h_1+h_3
\end{bmatrix}\in\mathbb R^{4\times1}.
$$

---

## 2. `coeffs`：如何从力矩公式提取参数单项式？

对应代码：

```matlab
theta_full = sym(zeros(0, 1));
W_symbolic = sym(zeros(joint_count*n_samples, 0));
term_names = {};

for sample_idx = 1:n_samples
    for joint_idx = 1:joint_count
        row_idx = (sample_idx-1)*joint_count + joint_idx;
        [term_coefficients, parameter_terms] = ...
            coeffs(expand(tau_samples(joint_idx, sample_idx)), parameter_symbols);
        ...
    end
end
```

设置：

```matlab
parameter_symbols = [h1, h2, h3];
```

对于第 1 次采样、第 1 个关节：

$$\tau_1^{(1)}=h_1+h_3=1\cdot h_1+1\cdot h_3.$$

使用 `coeffs` 后，可以识别出：

| 参数单项式 `parameter_terms` | 系数 `term_coefficients` |
|---|---:|
| $h_1$ | 1 |
| $h_3$ | 1 |

对于第 1 次采样、第 2 个关节：

$$\tau_2^{(1)}=h_2+h_3=1\cdot h_2+1\cdot h_3.$$

提取到：

| 参数单项式 | 系数 |
|---|---:|
| $h_2$ | 1 |
| $h_3$ | 1 |

**注意：** MATLAB 的 `coeffs` 不保证按这里展示的顺序返回项。为了便于手算，下面统一将原始参数列重排为 $[h_1,h_2,h_3]$。你的程序中实际列顺序由首次出现时的 `term_names` 决定；即使顺序不同，原理完全一样。

### 2.1 `ismember` 为什么重要？

你的代码：

```matlab
term_name = char(parameter_terms(term_idx));
[found, column_idx] = ismember(term_name, term_names);

if ~found
    column_idx = numel(term_names) + 1;
    term_names{column_idx} = term_name;
    theta_full(column_idx, 1) = parameter_terms(term_idx);
    W_symbolic(:, column_idx) = sym(zeros(joint_count*n_samples, 1));
end

W_symbolic(row_idx, column_idx) = ...
    W_symbolic(row_idx, column_idx) + term_coefficients(term_idx);
```

- 遇到新的参数单项式：给它新建一列。
- 再次遇到同一参数单项式：找到它原来的列，不重复建列。
- 每行的系数写入该参数对应的列；没有出现的参数其系数为 0。

本例提取完成后，原始参数向量为：

$$
\boxed{\theta_{\mathrm{full}}=\begin{bmatrix}h_1\\h_2\\h_3\end{bmatrix}}
$$

系数矩阵为：

$$
\boxed{
W_{\mathrm{symbolic}}=
\begin{array}{c|ccc}
 &h_1&h_2&h_3\\\hline
\tau_1^{(1)}&1&0&1\\
\tau_2^{(1)}&0&1&1\\
\tau_1^{(2)}&0&1&1\\
\tau_2^{(2)}&1&0&1
\end{array}}
$$

去掉行列标签，这就是：

$$
\boxed{W_{\mathrm{symbolic}}=
\begin{bmatrix}
1&0&1\\
0&1&1\\
0&1&1\\
1&0&1
\end{bmatrix}\in\mathbb R^{4\times3}}
$$

核对一次矩阵乘法：

$$
W_{\mathrm{symbolic}}\theta_{\mathrm{full}}=
\begin{bmatrix}
1&0&1\\0&1&1\\0&1&1\\1&0&1
\end{bmatrix}
\begin{bmatrix}h_1\\h_2\\h_3\end{bmatrix}
=\begin{bmatrix}h_1+h_3\\h_2+h_3\\h_2+h_3\\h_1+h_3\end{bmatrix}=T.
$$

**第一阶段结束：** 此时有 3 个参数单项式，还不知道几个是独立的。

---

## 3. 数值代入和列归一化：`W_scaled`

对应代码：

```matlab
W = double(subs(W_symbolic, known_symbols, known_values));
column_scale = vecnorm(W, 2, 1);
column_scale(column_scale == 0) = 1;
W_scaled = W ./ column_scale;
```

这个例子中 $W_{\mathrm{symbolic}}$ 已经全是常数，没有待代入的重力、几何量，因此直接有：

$$
W=\begin{bmatrix}1&0&1\\0&1&1\\0&1&1\\1&0&1\end{bmatrix}.
$$

将它按列写成：

$$
W=[w_1\ w_2\ w_3],\quad
w_1=\begin{bmatrix}1\\0\\0\\1\end{bmatrix},\quad
w_2=\begin{bmatrix}0\\1\\1\\0\end{bmatrix},\quad
w_3=\begin{bmatrix}1\\1\\1\\1\end{bmatrix}.
$$

先用肉眼发现：

$$\boxed{w_3=w_1+w_2}.$$

因此第三列不是新增的独立信息，整个矩阵最多只有 2 列独立。

但实际机器人的 $W$ 可能是 $56\times9$、$1000\times70$，不可能一直靠肉眼判断，因此才需要 SVD 和 QR。

### 3.1 每一列的二范数

$$
\|w_1\|_2=\sqrt{1^2+0^2+0^2+1^2}=\sqrt2,
$$

$$
\|w_2\|_2=\sqrt{0^2+1^2+1^2+0^2}=\sqrt2,
$$

$$
\|w_3\|_2=\sqrt{1^2+1^2+1^2+1^2}=2.
$$

所以：

$$\texttt{column\_scale}=[\sqrt2,\sqrt2,2].$$

归一化以后：

$$
\boxed{\widetilde W=W_{\mathrm{scaled}}=
\begin{bmatrix}
1/\sqrt2&0&1/2\\
0&1/\sqrt2&1/2\\
0&1/\sqrt2&1/2\\
1/\sqrt2&0&1/2
\end{bmatrix}}
$$

其中每列的二范数都是 1。

> 列归一化是为了改善数值秩判断的尺度问题，不会改变精确的列线性相关性。它**不是**简单把每列除以第一项，也不意味着所有数都变成整数。

---

## 4. SVD：为什么最终判断 `n_base = 2`？

对应代码：

```matlab
singular_values = svd(W_scaled);

if isempty(singular_values) || singular_values(1) == 0
    svd_threshold = 0;
    n_base = 0;
else
    svd_threshold = rank_tol * singular_values(1);
    n_base = sum(singular_values > svd_threshold);
end
```

SVD 的标准形式为：

$$\widetilde W=U\Sigma V^T.$$

非零奇异值的数量等于矩阵的秩。为了看到本例奇异值从哪里来，我们可以不用求完整的 $U,V$，只分析：

$$
\widetilde W^T\widetilde W=
\begin{bmatrix}
1&0&1/\sqrt2\\
0&1&1/\sqrt2\\
1/\sqrt2&1/\sqrt2&1
\end{bmatrix}.
$$

它的特征值为：

$$\lambda_1=2,\qquad\lambda_2=1,\qquad\lambda_3=0.$$

奇异值是这些非负特征值的平方根：

$$
\boxed{\sigma_1=\sqrt2,\qquad\sigma_2=1,\qquad\sigma_3=0.}
$$

例如设：

```matlab
rank_tol = 1e-10;
```

则：

$$
\texttt{svd\_threshold}=10^{-10}\sqrt2\approx1.4142\times10^{-10}.
$$

逐个比较：

| 奇异值 | 大约的数值 | 大于阈值？ | 是否计入秩 |
|---|---:|---|---|
| $\sigma_1$ | 1.4142 | 是 | 是 |
| $\sigma_2$ | 1.0000 | 是 | 是 |
| $\sigma_3$ | 0 | 否 | 否 |

因此：

$$\boxed{\texttt{n\_base}=2.}$$

这一步只给出一个**数量**：原始 3 个参数单项式对应的系数矩阵，只有 2 个线性独立方向。

**SVD 并没有自动告诉我们 $\pi_1=h_1+h_3$、$\pi_2=h_2+h_3$；它只告诉我们要找 2 个独立组合。**

---

## 5. 第一次带列主元 QR：从 3 列选出 2 列

对应代码：

```matlab
[~, R_sample, column_order] = qr(W_scaled, 'vector');
base_columns = sort(column_order(1:n_base));
```

普通 QR 写成 $W=QR$；带列置换的 QR 写成：

$$
\widetilde W(:,p)=QR,
$$

其中 $p$ 是列索引的排列，表示 MATLAB 选择列主元的顺序。

我们已知 $n_{\mathrm{base}}=2$，现在需要从 $w_1,w_2,w_3$ 中找到任意两列独立列。

由于：

$$w_3=w_1+w_2,$$

以下三种选法都可以：

- 第 1、2 列：$[w_1\ w_2]$；
- 第 1、3 列：$[w_1\ w_3]$；
- 第 2、3 列：$[w_2\ w_3]$。

它们都能张成同一个二维列空间。

### 5.1 为了手算清晰，先选第 1、2 列

**这是下面的手工演示所采用的一种合法 QR 选基结果，并非断言 MATLAB 必然返回这个顺序。** 本例归一化后各列长度相同，发生主元并列时，实际数值实现可能选出不同的基列；第 8 节会展示另一种合法结果。

假设：

```matlab
base_columns = [1, 2];
```

则：

$$
\boxed{W_B=W(:,[1,2])=
\begin{bmatrix}
1&0\\
0&1\\
0&1\\
1&0
\end{bmatrix}\in\mathbb R^{4\times2}.}
$$

它是原始 $W$ 的列空间的一组基，而第三列满足：

$$W(:,3)=W_B\begin{bmatrix}1\\1\end{bmatrix}.$$

这一步完成的是**选择哪些系数列作为独立列**，还没有求出最终的参数映射矩阵。


### 5.2 第一次 QR 的正交化过程，也可以手算

假设主元顺序为 $[1,2,3]$。针对**归一化矩阵**，有：

$$
\widetilde w_1=\frac{w_1}{\sqrt2},\qquad
\widetilde w_2=\frac{w_2}{\sqrt2},\qquad
\widetilde w_3=\frac{w_3}{2}.
$$

先取第一列作为正交基：

$$q_1=\widetilde w_1,\qquad r_{11}=\|\widetilde w_1\|_2=1.$$

因为 $w_1^Tw_2=0$，第二列对第一列的投影为 0：

$$r_{12}=q_1^T\widetilde w_2=0.$$

所以第二个正交基直接是：

$$q_2=\widetilde w_2,\qquad r_{22}=1.$$

第三列分别投影到 $q_1,q_2$：

$$r_{13}=q_1^T\widetilde w_3=\frac{1}{\sqrt2},\qquad
r_{23}=q_2^T\widetilde w_3=\frac{1}{\sqrt2}.$$

减去这两个投影：

$$
\widetilde w_3-r_{13}q_1-r_{23}q_2
=\widetilde w_3-\frac{\widetilde w_1+\widetilde w_2}{\sqrt2}
=0.
$$

**第三列的正交化剩余向量为零，这就是第三列线性相关的直接证据。** 对应一种合法的上三角矩阵：

$$
R=\begin{bmatrix}
1&0&1/\sqrt2\\
0&1&1/\sqrt2\\
0&0&0
\end{bmatrix}.
$$

实际 MATLAB 采用数值稳定的 QR 实现，并可能置换列，所以打印出的 `R_sample` 不一定正好是上面的形式；但它一定体现只有两个独立列方向（在计算精度允许的情况下）。


---

## 6. 第二次带列主元 QR：为什么要从 4 行中再挑 2 行？

对应代码：

```matlab
[~, ~, row_order] = qr(W_scaled(:, base_columns).', 'vector');
exact_rows = row_order(1:n_base);
A_exact = W_symbolic(exact_rows, :);
```

第一次选完列后的矩阵是 $4\times2$：

$$W_B=\begin{bmatrix}1&0\\0&1\\0&1\\1&0\end{bmatrix}.$$

我们想求某个矩阵 $P\in\mathbb R^{2\times3}$，满足：

$$\boxed{W=W_BP.}$$

因为 $W_B$ 不是方阵，不方便用一个方阵求解进行手工推导。因此，挑出 2 个独立行，形成可逆的 $2\times2$ 子矩阵。

### 6.1 转置后进行列主元 QR 的原因

$$
W_B^T\in\mathbb R^{2\times4}.
$$

**从 $W_B^T$ 挑独立列，等价于从 $W_B$ 挑独立行。**

我们的例子里，第 1 行 $[1,0]$ 和第 2 行 $[0,1]$ 显然独立；第 3 行和第 2 行重复，第 4 行和第 1 行重复。

为了方便手算，假设第二次 QR 选中：

```matlab
exact_rows = [1, 2];
```

那么：

$$
\boxed{A_{\mathrm{exact}}=W_{\mathrm{symbolic}}([1,2],:)=
\begin{bmatrix}
1&0&1\\
0&1&1
\end{bmatrix}\in\mathbb R^{2\times3}.}
$$

再取出第一次选中的两列：

$$
\boxed{A_B=A_{\mathrm{exact}}(:,[1,2])=
\begin{bmatrix}1&0\\0&1\end{bmatrix}\in\mathbb R^{2\times2}.}
$$

这里 $\det(A_B)=1\ne0$，因此可逆。

**第二次 QR 的作用，不是再次决定参数个数，而是找到能解出列间关系的独立行。**


如果想进一步看转置操作的具体形式，将归一化后的前两列取出并转置：

$$
\widetilde W_B^T=
\begin{bmatrix}
1/\sqrt2&0&0&1/\sqrt2\\
0&1/\sqrt2&1/\sqrt2&0
\end{bmatrix}.
$$

它的第 1、2 列相互正交；第 3 列与第 2 列相同，第 4 列与第 1 列相同。所以选其第 1、2 列，就相当于选原始 $W_B$ 的第 1、2 行。由此得到上面的 $A_B$ 可逆方阵。

> 实际 `row_order` 也可能选第 4 行代替第 1 行，或第 3 行代替第 2 行。只要选出的子矩阵可逆，计算出的列间关系就是一样的（行顺序不同也不影响正确结果）。

---

## 7. `parameter_map`：真正求出参数是如何组合的

你的代码最重要的一行：

```matlab
parameter_map = simplify(A_exact(:, base_columns) \ A_exact);
```

### 7.1 这行 MATLAB 究竟在求什么？

定义：

$$
P=\texttt{parameter\_map}.
$$

希望满足：

$$A_BP=A_{\mathrm{exact}}.$$

前面已经求出：

$$
A_B=\begin{bmatrix}1&0\\0&1\end{bmatrix},\qquad
A_{\mathrm{exact}}=\begin{bmatrix}1&0&1\\0&1&1\end{bmatrix}.
$$

所以：

$$
\begin{bmatrix}1&0\\0&1\end{bmatrix}P
=\begin{bmatrix}1&0&1\\0&1&1\end{bmatrix}.
$$

用数学形式写是：

$$P=A_B^{-1}A_{\mathrm{exact}}.$$

这里 $A_B$ 就是单位矩阵，因此：

$$
\boxed{P=\begin{bmatrix}1&0&1\\0&1&1\end{bmatrix}.}
$$

MATLAB 的反斜杠 `A_B \ A_exact` 表示**解线性方程** $A_BP=A_{\mathrm{exact}}$，而不是按元素相除；实践中比显式写 `inv(A_B)*A_exact` 更合适。

### 7.2 再验证：所有行都满足 $W=W_BP$

$$
W_BP=
\begin{bmatrix}
1&0\\0&1\\0&1\\1&0
\end{bmatrix}
\begin{bmatrix}1&0&1\\0&1&1\end{bmatrix}
=\begin{bmatrix}
1&0&1\\0&1&1\\0&1&1\\1&0&1
\end{bmatrix}=W.
$$

说明不仅选出的 2 行成立，全部 4 行都成立。

### 7.3 最后一行如何生成 `pi_param`？

对应你的代码：

```matlab
pi_param = simplify(parameter_map * theta_full);
```

代入：

$$
\pi=P\theta_{\mathrm{full}}
=\begin{bmatrix}1&0&1\\0&1&1\end{bmatrix}
\begin{bmatrix}h_1\\h_2\\h_3\end{bmatrix}
=\begin{bmatrix}h_1+h_3\\h_2+h_3\end{bmatrix}.
$$

因此：

$$
\boxed{\pi_1=h_1+h_3,\qquad\pi_2=h_2+h_3.}
$$

**这就是从原始 3 个参数单项式提取出的 2 个独立参数组合。**

最后将原来的力矩表达式变成：

$$
T=W\theta_{\mathrm{full}}=W_BP\theta_{\mathrm{full}}=W_B\pi.
$$

这说明采样力矩只与 $\pi_1,\pi_2$ 有关，没必要把 $h_1,h_2,h_3$ 都当作独立参数来辨识。

---

## 8. 重点：为什么你实际运行得到的 `pi_param` 可能与手算顺序不同？

由于有多个合法的独立列组合，带列主元 QR **不保证**选择 $[1,2]$。

例如假设 MATLAB 实际选择：

```matlab
base_columns = [2, 3];
```

那么基矩阵为：

$$
W_B=[w_2\ w_3]=
\begin{bmatrix}
0&1\\
1&1\\
1&1\\
0&1
\end{bmatrix}.
$$

由 $w_3=w_1+w_2$ 得：

$$w_1=-w_2+w_3.$$

于是所有列在新基下的系数为：

$$
\boxed{P=\begin{bmatrix}-1&1&0\\1&0&1\end{bmatrix}.}
$$

对应的新参数组合：

$$
\boxed{
\pi=\begin{bmatrix}-h_1+h_2\\h_1+h_3\end{bmatrix}}
$$

这与第 7 节的 $[h_1+h_3,\ h_2+h_3]^T$ 长得不一样，但仍然完全正确，因为：

$$
\tau_1=b(-h_1+h_2)+(a+b)(h_1+h_3)
=a h_1+b h_2+(a+b)h_3,
$$

$$
\tau_2=a(-h_1+h_2)+(a+b)(h_1+h_3)
=b h_1+a h_2+(a+b)h_3.
$$

**结论：最小参数集的维数（本例为 2）是关键；最小参数向量的具体表达形式不唯一。** 不同的 QR 基列选择会给出不同、但彼此等价的参数组合。

---

## 9. 为什么最小参数集不能直接辨识全部原始参数？

在本例中：

$$
\tau_1=a(h_1+h_3)+b(h_2+h_3),
\qquad
\tau_2=b(h_1+h_3)+a(h_2+h_3).
$$

假设同时做以下参数变化：

$$
h_1'=h_1-\delta,\quad
h_2'=h_2-\delta,\quad
h_3'=h_3+\delta.
$$

则：

$$
h_1'+h_3'=h_1+h_3,\qquad
h_2'+h_3'=h_2+h_3.
$$

无论 $a,b$ 取什么值，关节力矩始终相同！

因此单靠此力矩模型，$h_1,h_2,h_3$ 无法分别唯一辨识。真正能唯一对应输入输出关系的是两个组合 $\pi_1,\pi_2$（当然，想从**有限实验数据**可靠辨识这两个组合，还必须保证采样的回归矩阵具有足够的秩，并考虑噪声）。

这也是你研究**机器人自适应控制、复合学习、动力学在线辨识**时需要先提取最小参数集的原因：消除那些在动力学方程中无法独立区分的参数方向。

---

## 10. 可直接运行的 MATLAB 完整验证代码

下面这段脚本沿用你的循环、`coeffs`、SVD、两次 QR、参数映射的全部核心结构。为了便于独立运行，本例无需 `known_symbols`、`known_values`，因此直接把已是纯数值的 `W_symbolic` 转成 `double`。

```matlab
clear; clc;

%% 0) 定义一个简单的双关节符号动力学模型
syms h1 h2 h3 a b real
joint_count = 2;
n_samples   = 2;
rank_tol    = 1e-10;

% h1、h2、h3 都是待辨识的惯性类参数单项式
theta_expected = [h1; h2; h3];
parameter_symbols = [h1, h2, h3];

% 两个关节的原始符号力矩
tau_symbolic = [
    a*h1 + b*h2 + (a+b)*h3;
    b*h1 + a*h2 + (a+b)*h3
];

%% 1) 采样加速度：列代表不同采样
% 第1列：(a,b)=(1,0)
% 第2列：(a,b)=(0,1)
acc_samples = [1, 0;
               0, 1];

tau_samples = sym(zeros(joint_count, n_samples));
for k = 1:n_samples
    tau_samples(:, k) = subs(tau_symbolic, [a, b], acc_samples(:, k).');
end

disp('===== tau_samples =====');
disp(tau_samples);

%% 2) 使用 coeffs 构造 theta_full、W_symbolic
theta_full = sym(zeros(0, 1));
W_symbolic = sym(zeros(joint_count*n_samples, 0));
term_names = {};

for sample_idx = 1:n_samples
    for joint_idx = 1:joint_count
        row_idx = (sample_idx-1)*joint_count + joint_idx;

        [term_coefficients, parameter_terms] = coeffs( ...
            expand(tau_samples(joint_idx, sample_idx)), parameter_symbols);

        for term_idx = 1:numel(parameter_terms)
            term_name = char(parameter_terms(term_idx));
            [found, column_idx] = ismember(term_name, term_names);

            if ~found
                column_idx = numel(term_names) + 1;
                term_names{column_idx} = term_name;
                theta_full(column_idx, 1) = parameter_terms(term_idx);
                W_symbolic(:, column_idx) = ...
                    sym(zeros(joint_count*n_samples, 1));
            end

            W_symbolic(row_idx, column_idx) = ...
                W_symbolic(row_idx, column_idx) + term_coefficients(term_idx);
        end
    end
end

parameter_count = numel(theta_full);

disp('===== theta_full（实际列顺序以程序输出为准）=====');
disp(theta_full);
disp('===== W_symbolic =====');
disp(W_symbolic);

% 验证原始采样力矩的重构
tau_stacked = reshape(tau_samples, [], 1);
assert(isequal(simplify(W_symbolic*theta_full - tau_stacked), ...
               sym(zeros(joint_count*n_samples, 1))));

%% 3) 列归一化 + SVD 求数值秩
% 本例已没有待代入的 g、L 等符号，直接转为数值
W = double(W_symbolic);
column_scale = vecnorm(W, 2, 1);
column_scale(column_scale == 0) = 1;
W_scaled = W ./ column_scale;

singular_values = svd(W_scaled);
if isempty(singular_values) || singular_values(1) == 0
    svd_threshold = 0;
    n_base = 0;
else
    svd_threshold = rank_tol * singular_values(1);
    n_base = sum(singular_values > svd_threshold);
end

fprintf('\n===== SVD =====\n');
disp(singular_values);
fprintf('数值秩 n_base = %d\n', n_base);

%% 4) 第一次 QR：找独立列
[~, R_sample, column_order] = qr(W_scaled, 'vector');
base_columns = sort(column_order(1:n_base));

disp('===== 第一次 QR 的主元列选择顺序 =====');
disp(column_order);
disp('===== base_columns =====');
disp(base_columns);

%% 5) 第二次 QR：找独立行
if n_base == 0
    parameter_map = sym(zeros(0, parameter_count));
else
    [~, ~, row_order] = qr(W_scaled(:, base_columns).', 'vector');
    exact_rows = row_order(1:n_base);

    A_exact = W_symbolic(exact_rows, :);
    A_B = A_exact(:, base_columns);

    disp('===== exact_rows =====');
    disp(exact_rows);
    disp('===== A_exact =====');
    disp(A_exact);
    disp('===== A_B（可逆方阵） =====');
    disp(A_B);

    %% 6) 求参数映射矩阵 P
    parameter_map = simplify(A_B \ A_exact);
end

%% 7) 求最小参数组合
pi_param = simplify(parameter_map * theta_full);

disp('===== parameter_map =====');
disp(parameter_map);
disp('===== pi_param =====');
disp(pi_param);

%% 8) 关键验证：是否能完整重建原始系数矩阵
W_base_symbolic = W_symbolic(:, base_columns);
residual_symbolic = simplify(W_symbolic - W_base_symbolic*parameter_map);

disp('===== W - W_base * P（应为零矩阵） =====');
disp(residual_symbolic);

assert(isequal(residual_symbolic, sym(zeros(size(W_symbolic)))), ...
       '参数映射无法重构原始采样系数矩阵！');

% 再验证力矩（按列堆叠形式）
assert(isequal(simplify(tau_stacked-W_base_symbolic*pi_param), ...
               sym(zeros(joint_count*n_samples, 1))));

fprintf('验证通过：W = W_base * parameter_map\n');

%% 9) 额外：不依赖 QR 选列顺序的手工标准答案
% 用 [h1+h3; h2+h3] 直接重组原始任意状态力矩
pi_hand = [h1+h3; h2+h3];
Y_hand = [a, b;
          b, a];

assert(isequal(simplify(tau_symbolic - Y_hand*pi_hand), ...
               sym(zeros(joint_count,1))));

disp('===== 手工标准答案 pi_hand =====');
disp(pi_hand);
disp('===== 任意状态的符号回归矩阵 Y_hand =====');
disp(Y_hand);
```

### 运行时应该关注哪些输出？

- `n_base` 应当是 **2**。
- `singular_values` 约为 $[1.4142,\ 1.0000,\ 0]^T$，最后一个可能是接近机器精度的极小数。
- `parameter_map` 的具体形式取决于 `theta_full` 的顺序和 `base_columns` 的选择，**不要求必须等于第 7 节的矩阵**。
- `pi_param` 应有 **2 个表达式**，可能与手工答案形式不同，但应当等价。
- `W - W_base*P` 应为全零矩阵；这比仅检查 `n_base` 更重要。
- 手工标准答案中的 $Y_{\mathrm{hand}}$ 和 $\pi_{\mathrm{hand}}$ 对**任意** $a,b$ 都满足原始动力学方程。

---

## 11. 把例子一一对应回你原来的 2 关节、28 次采样代码

| 本文手算示例 | 你原来的机器人程序 | 作用 |
|---|---|---|
| `joint_count = 2` | `joint_count = 2` | 每次采样有 2 个关节力矩方程 |
| `n_samples = 2` | `n_samples = 28` | 采样次数 |
| 3 个参数单项式 | 你此前得到 9 个 | `theta_full` 的维数 |
| $W\in\mathbb R^{4\times3}$ | $W\in\mathbb R^{56\times9}$ | 堆叠后的采样系数矩阵 |
| `n_base = 2` | 你此前得到 `n_base = 5` | SVD 识别出的数值秩 |
| `base_columns` 选 2 列 | 选 5 列 | 第一次 QR 选择独立列 |
| `exact_rows` 选 2 行 | 选 5 行 | 第二次 QR 选择独立行 |
| $A_B\in\mathbb R^{2\times2}$ | $A_B\in\mathbb R^{5\times5}$ | 用于符号线性方程求解的方阵 |
| $P\in\mathbb R^{2\times3}$ | $P\in\mathbb R^{5\times9}$ | 把原始单项式映射到独立参数组合 |
| $\pi\in\mathbb R^{2}$ | $\pi\in\mathbb R^{5}$ | 最终用于动力学辨识的最小参数向量 |

### 需要特别区分的三件事

1. **SVD 确定有几个：** $n_{\mathrm{base}}=\operatorname{rank}(W)$，前提是采样足够丰富、阈值合理。
2. **QR 决定选哪几个：** 选择基列与独立采样行；它不直接产生固定参数组合。
3. **符号线性求解确定怎么组合：** 解 $A_BP=A_{\mathrm{exact}}$，再算 $\pi=P\theta_{\mathrm{full}}$。

最后要记住：这段代码从**有限采样**中推断参数结构。对于真实机器人，应该再用未参与选基的随机样本验证 $W=W_BP$，最好结合原始符号动力学表达式验证对一般 $q,\dot q,\ddot q$ 也成立。**随机样本求出的数值秩不能单独替代解析的最小参数集证明。**

### 一行公式概括全部代码

$$
\boxed{
\underbrace{T}_{\text{采样力矩}}
=\underbrace{W}_{\text{完整采样系数矩阵}}
\underbrace{\theta_{\mathrm{full}}}_{\text{全部参数单项式}}
=\underbrace{W_B}_{\text{独立系数列}}
\underbrace{P\theta_{\mathrm{full}}}_{\displaystyle\pi_{\mathrm{param}}}
}
$$

**目的不是把原始物理参数全部求出来，而是找到能完整描述关节力矩的最少独立参数组合。**
