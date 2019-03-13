import pandas as pd
import sys

if __name__ == '__main__':
    df = pd.read_csv(sys.argv[1])
    DSs = df.DS.unique()
    for DS in DSs:
        tmpdf = df[df.DS == DS]
        for tpe in ['n', 's']:
            m, s = tmpdf[(tmpdf.type == tpe)]['result'].mean(), tmpdf[(tmpdf.type == tpe)]['result'].std() 
            print(DS, tpe, m, s)
