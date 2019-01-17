import pandas as pd
df = pd.read_csv('log')
DSs = df.DS.unique()
for DS in DSs:
    tmpdf = df[df.DS == DS]
    for i in [500, 900, 1000]:
        print(DS, i, tmpdf[str(i)].mean(), tmpdf[str(i)].std())

