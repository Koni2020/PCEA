PCEA
=====

.. image:: https://img.shields.io/badge/Host-pyCEA%2FREADME-orange
   :target: https://github.com/Koni2020/pyCEA/blob/master/README.md

.. image:: https://img.shields.io/badge/Python-3.10-blue

.. image:: https://img.shields.io/badge/Status-Building-green

Catalog
-------

- `What is PCEA? <#what-is-pycea>`_
- `Installation <#installation>`_
- `Dependencies <#dependencies>`_
- `Usage <#usage>`_
- `Documentation <readme/DOCUMENTATION_CN.md>`_
- `README_ZN [Chinese] <readme/README_CN.md>`_

What is PCEA?
--------------


pyCEA (Python Compound Event Analysis) identifies events in signals using specific thresholds and computes their occurrence time, intensity, peak, and duration. It also detects event chains (compound events) and calculates their probabilities, such as:

- ENSO-drought-wildfire chain events
- Heatwave-drought compound events
- Drought-flood compound events

Compound events are the superposition of events within a specific time window, exhibiting six types of relationships.

Dependencies
------------

The following Python packages are required:

- `numpy <https://numpy.org/>`_
- `pandas <https://pandas.pydata.org/>`_

Installation
------------

Install pyCEA via pip::

   pip install PCEA

Usage
-----

Here are two quick examples of using pyCEA:

1. Reading data from a CSV file::

   import numpy as np
   from PCEA import CEA
   import pandas as pd

   # Read data from demo.csv, where columns represent variables and rows represent sampling times.
   ts = pd.read_csv("./data/demo.csv", index_col=0, header=0)
   ts = ts.iloc[:, [0, 1]]
   cea = CEA(ts, delta=6, threshold=[-np.inf, -0.5])
   cea.run_cea(save_path='./data/resuls.xlsx')  # Run compound event analysis

2. Input a boolean matrix::

   from PCEA import CEA
   import numpy as np

   ts = np.random.choice([True, False], [720, 3])  # Generate a boolean matrix

   cea = CEA(ts, delta=3, is_binary_array=True)  # Set "is_binary_array" to True for boolean input
   cea.run_cea(save_path='./data/resuls.xlsx')  # Run compound event analysis

More examples can be found in `Jin et al. (2025) <tutorial/compound_event_analysis.ipynb>`_.

References
----------

- `Donges J F, Schleussner C F, Siegmund J F, et al. Event coincidence analysis for quantifying statistical interrelationships between event time series: On the role of flood events as triggers of epidemic outbreaks[J]. The European Physical Journal Special Topics, 2016, 225: 471-487. <https://link.springer.com/article/10.1140/epjst/e2015-50233-y>`_
