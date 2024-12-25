**函数**: `cea`

1. **`si`** (DataFrame | ndarray):  
   
   - 输入信号数据，必选。
   - 输入信号数据要求为 $m_{time} \times n_{var}$ 的矩阵,行与列分别为采样点及变量。其可以是 Pandas 的 DataFrame 或 NumPy 数组。若
     输入原始信号矩阵`arr<float>`，则事件筛选的阈值必须提供。若输入为二进制矩阵`arr<bool>`，每个值代表是否超过阈值，
     则事件筛选的阈值不需要提供。

2. **`threshold`** (ndarray | list | tuple, 默认值: `[-np.inf, -0.5]`):  
   
   - 用于事件筛选的阈值，可选。  
   - 如果指定，函数会根据 `operator` 参数进行筛选。如果提供 ndarray, 那么必须为维度为2，shape为`(n_var, 2)`的矩阵或者维度为1
     shape为`(2, )`的矩阵。若维度为1，则默认所有变量使用该阈值。若提供为list和tuple类型，则长度必须为2，且默认所有变量使用该阈值。

3. **`tau`** (int | list | tuple | ndarray, 默认值: `np.ones(n_var)`):  
   
   - 复合时间窗口大小，可选。  
   - 表示复合事件的时间窗口，即变量$B$的事件开始时间，是否在变量$A$结束时间$\pm \tau$窗口内。若提供为int类型，则默认所有变量$n_{var}-1$间
   使用该大小。若提供为`list<int>`或`tuple<int>`，那么长度应该为$n_{var}-1$；若长度为1，则默认所有变量使用该窗口大小。

4. **`delta`** (int | tuple | ndarray, 默认值: `np.ones(n_var)`):  
   
   - 单次事件最小持续时间，可选。
   - 表示单次事件最小持续事件，例如洪涝最小持续事件。若提供为int类型，则默认所有变量$n_{var}$间
   使用该大小。若提供为`list<int>`或`tuple<int>`，那么长度应该为$n_{var}$；若长度为1，则默认所有变量使用该窗口大小。

5. **`operator`** (list | tuple, 默认值: `('ge', 'le')`):  
   
   - 阈值条件的区间操作符。  
     - `'ge'`: 大于等于  
     - `'le'`: 小于等于  
     - `'g'`: 大于
     - `'g'`： 小于
   - 可组合用于多条件筛选。

6. **`max_gap_length`** (int | list | tuple | ndarray, 默认值: `None`):  
   
   - 允许单次事件中断的最大长度。  
   - 用于过滤掉间隔过长的事件序列。

7. **`max_gap`** (int | list | tuple | ndarray, 默认值: `None`):  
   
   - 限制信号间的最大差异。  
   - 当信号间的差异大于 `max_gap` 时，认为是独立的事件。

8. **`is_binary_array`** (bool, 默认值: `False`):  
   
   - 是否将输入信号视为二值化数组。  
    - `True`：表示输入信号已经被处理为二值化形式。
    - `False`：表示输入信号是原始值，需要提供阈值。

9. **`get_compound_event`** (bool, 默认值: `True`):  
   
   - 是否获取复合事件。  
   
   - `True`：函数会返回复合事件的统计量。
   - `False`：表示输入信号是原始值，需要提供阈值。