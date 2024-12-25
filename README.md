# What is the pyCEA?
pyCEA stands for Python Compound Event Analysis. pyCEA is capable not only of identifying events within a signal under specific thresholds but also of detecting cascading events (compound events). The signals can be of various types, such as: El Niño–drought–wildfire cascading events, heatwave–drought compound events, or drought–flood compound
# Dependencies
For the installation of pyLiang, the following packages are required:
* [numpy](https://numpy.org/)
* [pandas]()
* [numba]()
# Installation
pyLiang can be installed using pip\
```pip install pyCEA```
# Usage
A quick example of pyLiang usage is as follow. 
```python
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
# README.md
- zh_CN[简体中文](readme/README_CN.md)