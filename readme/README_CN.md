# pyCEA
[![](https://img.shields.io/badge/主页-pyCEA%2FREADME_CN-orange)](https://github.com/Koni2020/pyCEA/blob/master/README.md)
![](https://img.shields.io/badge/Python-3.10-blue)
![](https://img.shields.io/badge/状态-Building-green)
## 目录
- [pyCEA简介](#pycea是什么)
- [依赖项](#依赖项)
- [安装](#安装)
- [使用](#usage)
- [文档](DOCUMENTATION_CN.md)
- [README](../README.md)

## pyCEA是什么？

pyCEA 是 Python Compound Event Analysis 的缩写。pyCEA以事件分析为基础，识别变量间的事件及其复合，事件复合连级指时间上变量间同步发生。例如：厄尔尼诺—干旱—野火连级事件；热浪—干旱复合事件；干旱—洪涝复合事件。\
它与R包[`CoinCalc`](https://github.com/JonatanSiegmund/CoinCalc)在一定程度上相似，但存在不同。
- 单独事件上，pyCEA考虑更多参数，允许变量事件中断。它能够识别特定阈值下，每个变量事件出发生的时间，强度，峰值，持续时间。
- 在复合事件上，还能识别事件间的连级（复合事件），并统计其概率。

## 依赖项

若要运行pyCEA,以下的包需要安装到当前Python环境。

* [numpy](https://numpy.org/)
* [pandas]()
* [numba]()

## 安装

pyCEA的安装很简单，在当前的Python环境内使用pip即可安装\
```pip install pyCEA```

## 使用

pyCEA目前支持两种输入数据类型，（1）原始$S_{mn}$信号矩阵,m为采样点，n为时间\
以下是一个简单的例子使用，更多信息可以参考tutorial里的示例脚本\
1 输入一个mxn的表格或者矩阵的信号，行为采样点，列是变量。

```python
# 导入所需要的库
import numpy as np
from pyCEA import CEA
import pandas as pd

# 读入data/demo.csv的数据,该数据列表示不同变量，行代表采样时间。
ts = pd.read_csv("./data/demo.csv", index_col=0, header=0)

cea = CEA(ts, delta=3, threshold=[-np.inf, -0.5], tau=3) # 关注小于-0.5即干旱部分, 窗口为3的干旱连级
cea.run_cea() #运行复合事件分析

cea.summary() # 输出运行结果到终端
cea.event_trip_info.to_excel("变量单次事件信息.xlsx")
```

2 输入一个mxn的 bool 表格或者矩阵。

```python
from pyCEA import CEA
import numpy as np

ts = np.random.choice([True, False], [720, 3]) # 生成一个bool矩阵，bool代事件发生

cea = CEA(ts, delta=3, is_binary_array=True) # 如果输入已经是一个bool矩阵，\
# 那么需要设置“is_binary”为TRUE,且不需要提供阈值参数“threshold”
```
     
## 参考文献
