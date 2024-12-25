## pyCEA是什么？

pyCEA 是Python Compound Event Analysis 的缩写。pyCEA不仅能够识别特定阈值下信号内的事件还能识别事件的连级（复合事件）。信号可以为多种
类型，例如：厄尔尼诺—干旱—野火连级事件；热浪—干旱复合事件；干旱—洪涝复合事件。

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

## 函数说明

**函数**: `cea`

1. **`si`** (DataFrame | ndarray):  
   
   - 输入信号数据，必选。
   - 输入信号数据要求为$m_{time} \times n_{var}$的矩阵,行与列分别为采样点及变量。其可以是 Pandas 的 DataFrame 或 NumPy 数组。若
     输入原始信号矩阵`arr[float]`，则事件筛选的阈值必须提供。若输入为二进制矩阵`arr[bool]`，每个值代表是否超过阈值，
     则事件筛选的阈值不需要提供。

2. **`threshold`** (ndarray | list | tuple, 默认值: `[-np.inf, -0.5]`):  
   
   - 用于事件筛选的阈值，可选。  
   - 如果指定，函数会根据 `operator` 参数进行筛选。如果提供 ndarray, 那么必须为维度为2，shape为`(n_var, 2)`的矩阵或者维度为1
     shape为`(2, )`的矩阵。若维度

3. **`tau`** (int | list | tuple | ndarray, 默认值: `np.ones(n_var)`):  
   
   - 复合时间窗口大小，可选。  
   - 表示复合事件的时间窗口，即变量$B$的事件开始时间，是否在变量$A$结束时间$\pm \tau$窗口内。

4. **`delta`** (int | tuple | ndarray, 默认值: `np.ones(n_var)`):  
   
   - 。  
   - 用于避免连续事件被误识别为同一事件。

5. **`operator`** (list | tuple, 默认值: `('ge', 'le')`):  
   
   - 阈值条件的操作符。  
     - `'ge'`: 大于等于  
     - `'le'`: 小于等于  
   - 可组合用于多条件筛选。

6. **`max_gap_length`** (int | list | tuple | ndarray, 默认值: `None`):  
   
   - 允许单次事件中断的最大长度。  
   - 用于过滤掉间隔过长的事件序列。

7. **`max_gap`** (int | list | tuple | ndarray, 默认值: `None`):  
   
   - 限制信号间的最大差异。  
   - 当信号间的差异大于 `max_gap` 时，认为是独立的事件。

8. **`is_binary_array`** (bool, 默认值: `False`):  
   
   - 是否将输入信号视为二值化数组。  
   - 如果为 `True`，表示输入信号已经被处理为二值化形式。

9. **`get_compound_event`** (bool, 默认值: `True`):  
   
   - 是否获取复合事件。  
   
   - 若为 `True`，函数会返回复合事件的统计量。
     
     # 参考文献
