# 机器人最小参数集构建
## 1 文件说明

## 2 执行步骤（动力学RNEA）
代码使用MATLAB进行机器人最小参数集的构建，构建流程如下：    
1、首先使用RNEA构建机器人tau方程  
未知量输入：质量m，质心c，连杆长度l，惯量矩阵I；  
输出为tau矩阵、M矩阵、C矩阵、G矩阵  
该流程参考书本：《Rigid Body Dynamics Algorithms》    
主要分为三个步骤：  
1、从基座到末端，求解每个刚体的速度和加速度；  
2、根据牛顿欧拉方程，求解刚体运动所需空间力；  
3、从末端到基座传递空间力，并且投影得到关节力矩    

## 3 参数集构建步骤

假如一个二关节机器人，在某个随机数下的力矩表达式为：

$$
\begin{aligned}
\tau_1 &= 3m_1 + 2m_2 l_{c2} + 5I_{zz1} \\
\tau_2 &= 4m_1 - m_2 l_{c2} + 2I_{zz1}
\end{aligned}
$$

可以发现，两个关节力矩共同包含三个参数单项式，分别为：

$$
\theta_{\text{full}} =
\begin{bmatrix}
m_1 \\
m_2 l_{c2} \\
I_{zz1}
\end{bmatrix}
$$

那么就可以得到矩阵：

$$
\underbrace{
\begin{bmatrix}
\tau_1 \\
\tau_2
\end{bmatrix}}_{\text{力矩}}
=
\underbrace{
\begin{bmatrix}
3 & 2 & 5 \\
4 & -1 & 2
\end{bmatrix}}_{W_{\text{symbolic}}}
\underbrace{
\begin{bmatrix}
m_1 \\
m_2 l_{c2} \\
I_{zz1}
\end{bmatrix}}_{\theta_{\text{full}}}

$$

假设进行28次随机数的输入，2关节机器人，那么每次采样都会产生两个关节的力矩表达式

$$
\tau^{(k)} =
\begin{bmatrix}
\tau_1^{(k)} \\
\tau_2^{(k)}
\end{bmatrix}
$$

经过28次输入之后，可以将所有参数表达为：

$$
T =
\begin{bmatrix}
\tau_1^{(1)} \\
\tau_2^{(1)} \\
\tau_1^{(2)} \\
\tau_2^{(2)} \\
\vdots \\
\tau_1^{(28)} \\
\tau_2^{(28)}
\end{bmatrix}
\in \mathbb{R}^{56}
$$

代码中的row_idx就是为了计算这些力矩放在矩阵的哪一行

之后，用SVD判断矩阵存在几个线性独立组合，即最后的 $\theta_{full}$ 大概是几维

然后第一次QR分解：

按理说，只要输入的tau包含q，dq，ddq，即使带有摩擦参数，也是可以得到pi矩阵和Y矩阵