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
dataframes = []  # 空のリストを初期化
for s in sorted(p.glob("**/md.vmcmd.hst")):
    tmp = pd.read_table(s, sep=r'\s+', header=None, names=cols)
    if not tmp.empty:  # Skip if `tmp` is empty
        dataframes.append(tmp)

if dataframes:  # 空リストでない場合にのみ結合
    df = pd.concat(dataframes)
else:
    df = pd.DataFrame(columns=cols)  # 空の DataFrame を初期化

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


# #!/usr/bin/env python
# # 
# # Written by Satoshi Ono
# #
# #
# import pandas as pd
# import numpy as np
# from pathlib import Path
# import kkkit
# import rangeinfo

# np.seterr(divide='ignore')

# ri = rangeinfo.RangeInfo('temp_s')
# ri.read_ttp_inp()
# nvst = len(ri.vs_order)
# vinterval = ri.vinterval


# #cols=['pot','v0','v1','v2','v3','v4','v5','v6','v7']
# cols=['pot']
# for i in range(nvst):
#     cols.append(str(i))
    
# p=Path(".")

# df=pd.DataFrame(index=[],columns=cols)
# for s in sorted(p.glob("**/md.vmcmd.hst")):
# #    tmp=np.loadtxt(str(s),unpack=True,dtype=int).T
#     tmp=pd.read_table(s,sep='\s+',header=None,names=cols)
# #    tmp=pd.read_table(s,delim_whitespace=True,header=None,names=cols)
#     df=pd.concat([df,tmp])
# df=df.groupby('pot').sum()
# total=df.sum().sum()

# for i in range(nvst):
#     df[str(i)]=np.log(df[str(i)]/total)
# df=df.replace(-np.inf,np.nan)
# for i in range(nvst):
#     df[str(i)].dropna().to_csv('e{0:d}.pdd'.format(i+1),sep=' ',float_format='%.8E')
