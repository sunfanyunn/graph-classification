f = open('../data/MUTAG_nodelabels.txt', 'r')
print(f.readline())
node_labels = []
for i in range(188):
    num = int(f.readline()[:-1])
    for j in range(num):
        node_labels.append(int(f.readline()[:-1]))

for i in range(8):
    print(node_labels.count(i))

