import numpy as np
from glob import glob
import sys
files = glob(sys.argv[1])
num=288
scores = [[] for i in range(num)]
for f in files:
    with open(f, 'r') as f:
        for idx, line in enumerate(f):
            scores[idx%288].append(float(line.strip().split()[-1]))


res = [np.mean(scores[idx]) for idx in range(num)]
amax = np.argmax(res)
print(amax, np.mean(scores[amax]), np.std(scores[amax]), len(scores[amax]))

