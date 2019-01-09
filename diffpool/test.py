import pandas as pd
from glob import glob
import numpy as np
import sys
"""
files = glob(f'{sys.argv[1]}/*/*')

for f in files:
    print(f)
    accs100 = []
    accs200 = []
    accs500 = []
    accs1000 = []
    cnt = 0
    for event in tf.train.summary_iterator(f):

        for value in event.summary.value:

            if value.tag == 'acc/val_acc':
                val_acc = value.simple_value
            else:
                continue
            

            cnt += 1
            if (cnt-100) % 1000 == 0:
                accs100.append(val_acc)

    print(cnt)

    print(np.mean(accs100), np.std(accs100))

     # print(value.tag)
     # if value.HasField('simple_value'):
  # print(value.simple_value)
"""

if __name__ == '__main__':
    # DS = sys.argv[1]

    df = pd.read_csv('log_')
    # df = df[df.DS == DS]
    gcs = df.gc.unique()
    types = df.method.unique()
    for gc in gcs:
        for tpe in types:
            tmpdf = df[(df.gc == gc) & (df.method == tpe)]
            for i in range(10, 110, 10):
                print(gc, tpe, i, tmpdf[str(i)+'-mean'].mean(), tmpdf[str(i) + '-mean'].std())
