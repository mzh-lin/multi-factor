# 因子风险估计

 对股票组合风险的预测主要包括两部分：因子收益协方差矩阵的预测与特异性收益协方差矩阵的风险建模。这一章主要讨论如何预测因子的协方差矩阵。



对于因子收益协方差矩阵$$F$$的预测，最简单的方法莫过于计算样本内因子收益率的历史协方差矩阵，这种方法简单假设因子收益率的波动是稳定的，这在实际应用中容易产生较大偏差。

Barra USE4主要的预测步骤包括：

1. 日度协方差估计（Newey-West 调整）`GEM2(2008)`
2. Eigenfactor 调整（特征值调整）`USE4(2011)`
3. Volatility Regime调整（波动率偏误调整）`USE4(2011)`

对于特异性收益方差矩阵的调整

1. 日度协方差估计（Newey-West 调整）`GEM2(2008)`
2. 结构化调整 `EUE3(2009)`
3. 贝叶斯压缩 `USE4(2011)`
4. Volatility Regime调整（波动率偏误调整）`USE4(2011)`

### 风险矩阵的估计：偏误估计量

$$
\begin{gather*}
b_{k,t}=\frac{r_{k,t~t+\delta}}{\sigma_{k,t}}\\
r_{k,t~t+\delta}:k资产在当前第t个截面到未来t+\delta个截面的收益率\\
\sigma_{k,t}：在截面t上，模型对资产k在截面t~t+\delta区间收益率波动率的预测值\\
B_k=\sqrt{\frac{1}{T-1}\sum_{t=1}^{T}(b_{k,t}-\bar{b}_k)^2}\\
\bar{b}_k:标准化收益b_{k,t}在T期内的均值
\end{gather*}
$$

### 因子收益率协方差矩阵：EWMA

利用因子日频收益率数据$f$，指数衰减加权移动平均计算因子收益率协方差矩阵$F^{Raw}$:
$$
F^{Raw}_{a,b}=cov(f_a,f_b)_t=\frac{\sum_{s=0}^h\lambda_{t-s}(f_{a,t-s}-\bar{f_a})(f_{b,t-s}-\bar{f_b})}{\sum_{s=0}^h\lambda_{t-s}}\\
\lambda_{t-s}=0.5^{s/\tau}\\
cov(f_a,f_b)_t:第t期，因子a与因子b之间的协方差\\
\lambda:指数衰减权重\\
f_{a,t-s},f_{b,t-s}:第t-s期，因子a和因子b之间收益率\\
\bar{f_{a}},\bar{f_{b}}:从截面t-h到截面t的区间内,因子a和因子b的收益率的加权均值\\
参数:rolling\ window(h)=252,半衰期(\tau)=90
$$

### NW调整

$F^{Raw}$存在自相关性，需要利用Newey-West调整:
$$
F^{NW}=21\cdot\hat{\Omega}=21\cdot[F^{Raw}+\sum_{d=1}^D(1-\frac{d}{D+1})(\hat{\Omega}_d+\hat{\Omega}_d')],\\
\hat{\Omega}_d=\frac{\sum_{t=1}^{T-d}\lambda^{T-d-t}f_tf_{t+d}'}{\sum_{t=1}^{T-d}\lambda^{T-d-t}}
\\
f_t:截面t上所有因子的收益率序列\\
\hat{\Omega}_d:滞后期为d时的自协方差矩阵
$$
参数：滞后期$D=2$

### Eigenfactor 调整（特征值调整）！！！

Shepard指出模型低估风险，即采样得到的风险小于真实风险。我们利用蒙特卡洛模拟，估计低估程度，假定低估程度=采样/真实=模拟/采样

1. 首先对协方差矩阵进行特征值分解:
   $$
   D_0=U_0'F^{NW}U_0\\
   D_0:特征值构成的对角矩阵\\
   U_0:正交矩阵,矩阵第k列表示D_0中第k个特征值对应的特征向量
   $$

2. 蒙特卡洛模拟

   ![image-20191121213433835](E:\GitHub\multi-factor\feng-xian-mo-xing\yin-zi-feng-xian-gu-ji.assets\image-20191121213433835.png)

3. 计算特征值偏误

   ![image-20191121213533483](E:\GitHub\multi-factor\feng-xian-mo-xing\yin-zi-feng-xian-gu-ji.assets\image-20191121213533483.png)

