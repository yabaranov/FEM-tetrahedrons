from calendar import c
import matplotlib.pyplot as plt

x = []
y = []
z = []

with open('../../res/output/nodes.txt', 'r') as f:
    for line in f.readlines():
        (xt,yt,zt) = line.split()
        x.append(float(xt))
        y.append(float(yt))
        z.append(float(zt))

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

ax.scatter(x, y, z, c=['black'])

plt.show()