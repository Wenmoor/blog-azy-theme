---
title: FFT 和 FWT 学习笔记
date: 2025.2.2
tags: FFT,FWT
categories: 学习笔记
---

尝试自己理解一下原理。

---

首先是 FFT，解决的问题是将多项式和点值表达式在 $O(n\log n)$ 时间内互相转化。

设要转化的多项式为 $F(x)$，次数为 $n$，那么我们知道代入 $n+1$ 个点值可以唯一确定这个多项式。

此处选用的点值为单位根，记 $\omega_n$ 为 $n$ 次单位根，即将平面直角坐标系上的单位圆均分成 $n$ 份，几何意义是逆时针旋转 $\dfrac{2\pi}{n} \mathrm{rad}$。

由此容易得到一些性质，$\omega_{2n}^{2k} = \omega_{n}^{k}$，$\omega_{n}^{k+\frac{n}{2}} = -\omega_n^k$ 等。

现在设多项式为 $F(x) = a_0+a_1x^1+a_2x^2+\cdots+a_{n-1}x^{n-1}$，要求的是 $F(\omega_n^0),F(\omega_n^1),\cdots,F(\omega_n^{n-1})$ 的值。

考虑分离奇偶项，设 $F_0(x) = a_0 + a_2x + a_4x^2+\cdots$，$F_1(x) = a_1+a_3x+a_5x^2+\cdots$。

那么有 $F(x) = F_0(x^2)+xF_1(x^2)$，尝试代入 $\omega_n^k$。

分类讨论一下，若 $k<\dfrac{n}{2}$，则有：

$$\begin{aligned}
    F(\omega_n^k)&=F_0(\omega_n^{2k})+\omega_n^kF_1(\omega_n^{2k})\\
    &=F_0(\omega_{\frac{n}{2}}^{k})+\omega_n^kF_1(\omega_{\frac{n}{2}}^{k})
\end{aligned}$$

容易发现变成了两个规模减半的子问题。

若 $k \ge \dfrac{n}{2}$，则设 $k' = k - \dfrac{n}{2}$，有：

$$\begin{aligned}
    F(\omega_n^k)&=F_0(\omega_{\frac{n}{2}}^{k})+\omega_n^{k'+\frac{n}{2}}F_1(\omega_{\frac{n}{2}}^{k})\\
    &=F_0(\omega_{\frac{n}{2}}^{k'})-\omega_n^{k'}F_1(\omega_{\frac{n}{2}}^{k'})
\end{aligned}$$

形式相同。

优化运算可以使用位运算置换，每次是将最低位奇偶分类，可以看作是二进制位逆序后的基数排序，所以一个数字的二进制位的逆序就是最终的位置，免去了递归计算的大常数。

根据定义式：$F(\omega_n^x)=\sum\limits_{i=0}^{n-1}a_i\omega_n^{xi}$ 可得这是一个线性变换，等价于乘范德蒙德矩阵。

由于这个矩阵的元素非常特殊，它的逆矩阵也有特殊的性质，就是每一项取倒数，再除以变换的长度 $n$，就能得到它的逆矩阵。

NTT 即将 $\omega_n$ 用 $g^{\frac{p - 1}{n}}$ 代替，$g$ 是模数 $p$ 的原根，只因它们有相似的性质。

---

FWT 部分。

仿照 FFT 直接设计线性运算，设 $c(i,j)$ 表示原多项式的第 $j$ 项对 FWT 后多项式的第 $i$ 项的贡献系数。

联立 $\mathrm{FWT}(A)_i\times\mathrm{FWT}(B)_i=\mathrm{FWT}(A\times B)_i$ 和 $\mathrm{FWT}(A)_i = \sum\limits_{j=0}^{n}c(i,j)A_j$ 可得 $c(i,j)c(i,k)=c(i,j\oplus k)$ 是等价条件，其中 $\oplus$ 是任意位运算，而位运算的一大优点在于其每位独立，故 $c(i,j)$ 相当于 $i,j$ 每一位的变换系数相乘。

计算 $\mathrm{FWT}(A)_i = \sum\limits_{j=0}^{n-1}c(i,j)A_j$ 仍然考虑分治，记 $i'$ 表示 $i$ 去掉最高位后的值，$i_0$ 表示 $i$ 的最高位，$j$ 同理。

$$\begin{aligned}
    \mathrm{FWT}(A)_i &= \sum\limits_{j=0}^{n-1}c(i,j)A_j\\
    &= \sum\limits_{j=0}^{\frac{n}{2}-1}c(i,j)A_j + \sum\limits_{j=\frac{n}{2}}^{n-1}c(i,j)A_j\\
    &= c(i_0,0)\sum\limits_{j=0}^{\frac{n}{2}-1}c(i',j')A_j + c(i_0,1)\sum\limits_{j=\frac{n}{2}}^{n-1}c(i',j')A_j\\
\end{aligned}$$

转化为两个规模减半的子问题，根据 $i_0$ 的值计算即可。

推演 $c([0,1],[0,1])$ 的方法：

以 and 卷积为例。

首先有 $c(i,0)c(i,j)=c(i,0)$，故有 $c(i,0)=0$ 或 $c(i,j)=1$。

由于这是一个线性运算，矩阵必须满秩，所以不能存在一行或一列均为 $0$。

由此得到合法的矩阵可能是 $\begin{bmatrix} 0 & 1\\ 1 & 1 \end{bmatrix}$ 或者 $\begin{bmatrix} 1 & 1\\ 0 & 1 \end{bmatrix}$，逆运算直接乘上矩阵的逆即可。

- 子集卷积

解决的问题是计算 $F_i = \sum\limits_{j \operatorname{and} k = 0 \wedge j \operatorname{or} k = i} a_jb_k$。

或卷积解决不了与为 $0$ 的条件，考虑将其转化为 $|j|+|k|=|i|$，那么只需要在 $F$ 内部维护一个多项式即可。

时间复杂度 $O(n^22^n)$。