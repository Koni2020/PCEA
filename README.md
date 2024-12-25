# pyCEA
[![](https://img.shields.io/badge/Host-pyCEA%2FREADME-orange)](https://github.com/Koni2020/pyCEA/blob/master/README.md)
![](https://img.shields.io/badge/Python-3.10-blue)
![](https://img.shields.io/badge/Status-Building-green)
## Catalog
- [Introduction of pyECA](#what-is-the-pycea)
- [Installation](#installation)
- [Dependencies](#dependencies)
- [Usage](#usage)
- [Documentation](readme/DOCUMENTATION_CN.md)
- [README_ZN [简体中文]](readme/README_CN.md)
## What is the pyCEA?

pyCEA is the abbreviation for Python Compound Event Analysis. 
pyCEA is based on event analysis. It can identify events within signals based on specific thresholds and calculate the occurrence time, intensity, peak, and duration of those events. 
It can also detect event chains (compound events) and calculate their probabilities. An event chain refers to two events that occur synchronously in time. 
Examples include: El Niño—drought—wildfire chain events; heatwave—drought compound events; drought—flood compound events.
## Dependencies
For the installation of pyECA, the following packages are required:
* [numpy](https://numpy.org/)
* [pandas]()
* [numba]()
## Installation
pyECA can be installed using pip\
```pip install pyCEA```
## Usage
Two quick examples of pyECA usage is as following. 
1. The data from data/demo.csv is read in, 
where the columns represent different variables and the rows represent sampling times.
```python
import numpy as np
from pyCEA import CEA
import pandas as pd

# Read the data in data/demo.csv, where the columns represent different 
# variables and the rows represent sampling times.
ts = pd.read_csv("./data/demo.csv", index_col=0, header=0)

cea = CEA(ts, delta=3, threshold=[-np.inf, -0.5], tau=3) 
cea.run_cea() # run compound event analysis

cea.summary() # print the running results to the terminal
cea.event_trip_info.to_excel("event_statistics.xlsx")
```
2. Input a boolean matrix of size m x n.

```python
from pyCEA import CEA
import numpy as np

ts = np.random.choice([True, False], [720, 3]) # Generate a boolean matrix, 
# where True represents the occurrence of an event.

cea = CEA(ts, delta=3, is_binary_array=True) # If the input is already a boolean matrix,
# then you need to set "is_binary" to TRUE, 
# and the threshold parameter "threshold" is not required.
```
More details can be seen in [`Demo.ipynb`](tutorial/compound_event_analysis.ipynb)