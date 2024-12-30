#!/usr/bin/env python
# 
# Written by Satoshi Ono
#
#
import pandas as pd
import numpy as np
from pathlib import Path
import kkkit
import rangeinfo
# import warnings
# warnings.simplefilter(action='ignore', category=FutureWarning)

np.seterr(divide='ignore')

ri = rangeinfo.RangeInfo('temp_s')
ri.read_ttp_inp()
nvst = len(ri.vs_order)
vinterval = ri.vinterval

# Define column names
cols = ['pot']
for i in range(nvst):
    cols.append(str(i))
    
# Initialize empty DataFrame
df = pd.DataFrame(index=[], columns=cols)

# Read and process histogram files
p = Path(".")
# for s in sorted(p.glob("**/md.vmcmd.hst")):
#     tmp = pd.read_table(s, sep=r'\s+', header=None, names=cols)
#     if not tmp.empty:  # Skip if `tmp` is empty
#         df = pd.concat([df, tmp])
dataframes = []  # 
for s in sorted(p.glob("**/md.vmcmd.hst")):
    tmp = pd.read_table(s, sep=r'\s+', header=None, names=cols)
    if not tmp.empty:  # Skip if `tmp` is empty
        dataframes.append(tmp)

if dataframes:  # 
    df = pd.concat(dataframes)
else:
    df = pd.DataFrame(columns=cols)  # 

# Group by 'pot' and sum the values
df = df.groupby('pot').sum()

# Total count for normalization
total = df.sum().sum()

# Normalize and apply logarithm
for i in range(nvst):
    column = str(i)
    df[column] = df[column].astype(float) / total  # Ensure column is numeric
    df[column] = np.log(df[column].replace(0, np.nan))  # Replace 0 with NaN to avoid log(0)
    
# Replace -inf with NaN (already handled by replacing 0 with NaN above, but extra safety)
df = df.replace(-np.inf, np.nan)

# Save each column to a separate file
for i in range(nvst):
    column = str(i)
    df[column].dropna().to_csv(f'e{i+1}.pdd', sep=' ', float_format='%.8E',header=False)
