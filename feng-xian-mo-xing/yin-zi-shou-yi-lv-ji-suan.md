# 因子收益率计算

在对因子暴露进行处理后，我们需要在每个交易日截面上，用t天的股票因子暴露对股票t+1天的收益进行回归，估计收益率参数。

结构化风险模型给出了任一股票收益率的线性分解形式：

$$
r_j=x_cf_c+\sum_iX_{ni}f_i+\sum_sX_{ns}f_s+u_n
$$

其中$$f_c$$ 为国家因子，而$$f_i$$为行业因子，$$f_s$$为风格因子,$$u_n$$为个股特质收益。

> 1. 是否加入国家因子需要进一步考虑，因为加入国家因子后，行业因子与国家因子存在多重共线性问题，需要限制行业的加权收益为0:
>
>    $$
>    \sum_iw_if_i=0
>    $$
>
> 2. 针对线性分解形式需要注意的是，线性分解并不包含常数项，因此仅针对风险因子进行回归，而没有截距项（原因？）

那么对于N只股票的组合而言，组合的收益率向量为：

$$
R=Xf+U
$$

对于风险模型而言，需要在已知股票收益率R和因子载荷矩阵X的情况下，对因子收益率向量$$f$$进行估计。

## 估计方法

### 最小二乘法\(OLS\)

我们需要找到收益率向量$$f$$ 使得残差平方和达到最小，即：

$$
\begin{align*}
Min\ Q&=\sum_{i=1}^{N}\varepsilon_i^2=\sum_{i=1}^{N}(r_i-\hat{r}_i)^2\\
&=（R-X\hat{f})'(R-X\hat{f})\\
&= R'R-R'X\hat{f}-\hat{f}'X'R+\hat{f}'X'X\hat{f}'\\
&= R'R-2\hat{f}'X'R+\hat{f}'X'X\hat{f}
\end{align*}
$$

令$$\frac{\partial Q}{\partial\hat{f}}=0$$，得到$$-X'R+X'X\hat{f}=0$$，因此$$\hat{f}=(X'X)^{-1}XR$$

OL S 最小二乘估计法仅当不同股票的残差序列$$\varepsilon_{it}$$方差相同时,$$\hat{f}$$才是最优估计。

> 通常情况下，金融时间序列的数据均存在较明显的异方差性，即每只股票的$$\varepsilon_{it}$$方差是不相同的。为了解决异方差性，通常采用广义最小二乘（GLS）估计方法。

### 广义最小二乘（GLS）

广义最小二乘法假设的$$\varepsilon_{it}$$方差不相同，即:

$$
Var(U)=\Sigma=\left| \begin{array}{}
   \sigma_1^2 & \sigma_{12} &\cdots& \sigma_{1n} \\
   \sigma_{21} & \sigma_2^2 &\cdots& \sigma_{2n} \\
   \cdots &\cdots&\cdots&\cdots\\
   \sigma_{n1} &\cdots&\cdots&  \sigma_n^2
  \end{array} \right|
$$

由于矩阵$$\Sigma$$是正定阵，因此可以写成$$\Sigma=KK'$$，其中$$K$$为非奇异矩阵

对于$$R=Xf+U$$,可得到：

$$
K^{-1}R=K^{-1}Xf+K^{-1}U
$$

令$$R^*=K^{-1}R,X^*=K^{-1}X,U^*=K^{-1}U,$$则可以得到：

$$
R^*=X^*f+U^*
$$

由于$$E(U^*)=0$$,并且

$$
\begin{align*}
Var(U^*)&=Var(K^{-1}U)\\
&=K^{-1}Var(U)K^{-1}\\
&=K^{-1}KK'(K')^{-1}\\
&=I

\end{align*}
$$

满足最小二乘估计（OLS）的同方差性条件。据OLS 中$$f$$ 的估计结果，可以得到

$$
\begin{align*}
\hat{f}_{GLS}&=((X^*)'X^*)^{-1}(X^*)'Y^*\\
&=((K^{-1}X)'K^{-1}X)^{-1}(K^{-1}X)'(K^{-1}R)\\
&=(X'(K^{-1})'K^{-1}X)^{-1}X'(K^{-1})'(K^{-1}R)\\

\end{align*}
$$

因为$$(K^{-1})'K^{-1}=(K')^{-1}K^{-1}=(KK')^{-1}=\Sigma^{-1}$$故

$$
\hat{f}_{GLS}=(X'\Sigma^{-1}X)^{-1}X'\Sigma^{-1}R
$$

广义最小二乘（GLS）估计法在已知残值波动率矩阵$$\Sigma$$ 的条件下，可得到$$f$$ 的无偏估计量$$\hat{f}_{GLS}$$。

而在结构化风险模型中，由于假设残差之间不存在相关性，因此可利用广义最小二乘GLS方法的特殊形式加权最小二乘法（WLS）处理。

### 加权最小二乘法（WLS）

加权最小二乘法WLS假设的$$\varepsilon_{it}$$方差不相同，但$$\varepsilon_{it}$$之间协方差为0，即：

$$
Var(U)=\Sigma=\left| \begin{array}{}
   \sigma_1^2 &  & & \\
    & \sigma_2^2 & &   \\
     & &\cdots& \\
    & & &  \sigma_n^2
  \end{array} \right|
$$

将上式中$$\Sigma$$写作

$$
\Sigma=\sigma_w^2\left| \begin{array}{}
   1/w_1 &  & & \\
    &1/w_2 & &   \\
     & &\cdots& \\
    & & &  1/w_n
  \end{array} \right|
$$

令$$W=diag(w_1,w_2,\cdots,w_n)$$,则$$\Sigma=\sigma_w^2W^{-1}$$,则$$\Sigma^{-1}=(1/\sigma_w^2)W$$。因此，根据GLS中$$\hat{f}_{GLS}$$的估计表达式可得：

$$
\begin{align*}
\hat{f}_{WLS}&=(X'\Sigma^{-1}X)^{-1}X'\Sigma^{-1}R\\
&=\sigma_w^2(X'WX)^{-1}X'(1/\sigma_w^2)WR\\
&=(X'WX)^{-1}X'WR\\
\end{align*}
$$

### 国家因子带来的线性限制

由于上述回归问题由于存在一个线性约束，不能直接回归求解，我们可以采用以下几种方法，

第一种是将约束条件带入回归方程，减少一个回归变量；

第二种方法是直接求一个带约束的优化问题，即最小化加权残差平方和。

我们采用第二种方法进行求解。

