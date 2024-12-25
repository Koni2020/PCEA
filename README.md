# pyCEA
[![](https://img.shields.io/badge/Host-pyCEA%2FREADME-orange)](https://github.com/Koni2020/pyCEA/blob/master/README.md)
![](https://img.shields.io/badge/Python-3.10-blue)
## Catalog
- [Introduction of pyECA](#what-is-the-pycea)
- [Installation](#installation)
- [Usage](#usage)
- [Documentation](#usage)
- [README_ZN (简体中文)](readme/README_CN.md)
## What is the pyCEA?
pyCEA stands for Python Compound Event Analysis. pyCEA is capable not only of identifying events within a signal under specific thresholds but also of detecting cascading events (compound events). The signals can be of various types, such as: El Niño–drought–wildfire cascading events, heatwave–drought compound events, or drought–flood compound
# Dependencies
For the installation of pyECA, the following packages are required:
* [numpy](https://numpy.org/)
* [pandas]()
* [numba]()
## Installation
pyECA can be installed using pip\
```pip install pyCEA```
## Usage
A quick example of pyLiang usage is as follow. 
```python
import numpy as np
from pyCEA import CEA
import pandas as pd

# Read the data in data/demo.csv, where the columns represent different 
# variables and the rows represent sampling times.
ts = pd.read_csv("./data/demo.csv", index_col=0, header=0)

cea = CEA(ts, delta=3, threshold=[-np.inf, -0.5], tau=3) # 关注小于-0.5即干旱部分, 窗口为3的干旱连级
cea.run_cea() # run compound event analysis

cea.summary() # print the running results to the terminal
cea.event_trip_info.to_excel("event_statistics.xlsx")
```
