from matplotlib.ft2font import VERTICAL
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
import numpy as np
from matplotlib.widgets import Slider, Button
from functools import partial

x = []
y = []
z = []
elements = []

def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)

def draw(index):
    ax.collections.pop()
    elem = elements[index]
    tmp = []
    for ind in elem:
        tmp.append([ x[ind], y[ind], z[ind] ])
    v = np.array(tmp)
    verts = []
    for i in range(4):
        verts.append([ v[i], v[(i + 1) % 4], v[(i + 2) % 4] ])

    cmap = get_cmap(len(elements))
    ax.add_collection3d(Poly3DCollection(verts, linewidths=1, facecolors=cmap(index), edgecolors='r', alpha=.25))
    plt.show()

def set_slider(s,val):
    s.val = round(val)
    s.poly.xy[2] = s.val,1
    s.poly.xy[3] = s.val,0
    s.valtext.set_text(s.valfmt % s.val)
    draw(s.val)

with open('../../res/output/nodes.txt', 'r') as f:
    for line in f.readlines():
        (xt,yt,zt) = line.split()
        x.append(float(xt))
        y.append(float(yt))
        z.append(float(zt))

with open('../../res/output/elements.txt', 'r') as f:
    for line in f.readlines():
        (p1, p2, p3, p4) = line.split();
        elements.append([int(p1), int(p2), int(p3), int(p4)])

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

#ax.scatter(x, y, z)

#axfreq = plt.axes([0.25, 0.1, 0.65, 0.03])
#elem_slider = Slider(
#    ax=axfreq,
#    label='Elem number',
#    valmin=0,
#    valmax = len(elements),
#    orientation='horizontal',
#    valinit=0.,
#    valfmt="%i"
#)
#
#elem_slider.on_changed(partial(set_slider,elem_slider))

for j in range(len(elements)):
    tmp = []
    for ind in elements[j]:
        tmp.append([ x[ind], y[ind], z[ind] ])
    v = np.array(tmp)
    verts = []
    for i in range(4):
        verts.append([ v[i], v[(i + 1) % 4], v[(i + 2) % 4] ])

    cmap = get_cmap(len(elements))
    ax.add_collection3d(Poly3DCollection(verts, linewidths=1, facecolors=cmap(j), edgecolors='black', alpha=0.25))

ax.scatter(x, y, z)
draw(0)
plt.show()