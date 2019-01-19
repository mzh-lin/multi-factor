# 因子定义

##### 样本空间选择

与大多海外成熟市场不同，A 股市场中市值最大、流动性最好的股票池并不能完整的反应A 股市场的全部特征。因此，我们对风险因子检验的样本空间选取全部A 股标的池。

##### 因子选取

针对Barra CNE5 我们进行相应因子定义的描述，其中共同因子主要分为行业因子及风险因子

###### 国家因子

Barra USE4认为加入国家因子会造成多重共线性，因为对于所有的股票$$n$$,行业因子的暴露如下:

$$

n

$$

因此要限制行业的加权收益为0

###### 行业因子

行业因子是风险模型的重要部分，通过对A 股市场全部个股的行业划分，反应了个股所属行业的独有特点。我们在中信一级29个行业分类的基础上，将非银金融行业划分为证券、保险和多元金融3 个子行业。因此，我们的风险模型包含31个行业因子，具体为:

| **石油石化** | **煤炭** | **有色金属** | **电力及公用事业** | **钢铁** | **基础化工** |
| :---: | :---: | :---: | :---: | :---: | :---: |
| **建筑** | **建材** | **轻工制造** | **机械** | **电力设备** | **国防军工** |
| **汽车** | **商贸零售** | **餐饮旅游** | **家电** | **纺织服装** | **医药** |
| **食品饮料** | **农林牧渔** | **银行** | **证券** | **保险** | **多元金融** |
| **房地产** | **交通运输** | **电子元器件** | **通信** | **计算机** | **传媒** |
| **综合** |  |  |  |  |  |

###### 风格因子

风格因子是共同因子中另一重要部分，风格因子总共包含9 大类因子、20 个小因子，其中大类因子包含Beta、Momentum、Size、Earnings Yield、Volatility、Growth、Value、Leverage 和Liquidity，具体为：

![](/assets/import.png)
