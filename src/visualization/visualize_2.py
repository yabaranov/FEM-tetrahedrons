import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.widgets import Button, Slider
import matplotlib.tri as mtri


def SameSide(v1, v2, v3, v4, p):
    normal = np.cross(v2-v1, v3-v1)
    return (np.dot(normal, v4-v1) * np.dot(normal, p-v1) >= 0)
 
  
def pointInside(verts, p):   
    return SameSide(verts[1], verts[0], verts[3], verts[5], p) and SameSide(verts[5], verts[4], verts[7], verts[1], p) and SameSide(verts[0], verts[1], verts[4], verts[2], p) and SameSide(verts[3], verts[2], verts[7], verts[1], p) and SameSide(verts[1], verts[3], verts[5], verts[0], p) and SameSide(verts[0], verts[4], verts[2], verts[1], p)
x = []
y = []
z = []
values = []
subdomains = []

with open('../../res/output/solution.txt', 'r') as f:
    for line in f.readlines():
        (xt, yt, zt, valuest) = line.split()
        x.append(float(xt))
        y.append(float(yt))
        z.append(float(zt))
        values.append(float(valuest))

with open('../../res/output/hexahedron subdomains.txt', 'r') as f:
    for line in f.readlines():
        subdomain = line.split()
        nodes = []
        for i in range(8):
            node = []
            for j in range(3):
                node.append(float(subdomain[3 * i + j]))
            nodes.append(node)
        subdomains.append(nodes)

unique_x = sorted(list(set(x)))
unique_y = sorted(list(set(y)))
unique_z = sorted(list(set(z)))

fig = plt.figure()

#--------------------------------------------------------------------------------------------------

def getMaskx(ymid, zmid, unx):    
    mask = []
    for i in range(len(ymid)):
        flag = True
        p = [unx, ymid[i], zmid[i]]   
        for j in range(len(subdomains)):
            if pointInside(np.array(subdomains[j]), np.array(p)):    
                flag = False
                break
        mask.append(flag)
    return mask
    
def getDatax(xVal):
    y_graph = []
    z_graph = []
    values_graph = []
    
    for i in range(len(x)):
        if(x[i]==xVal):
            z_graph.append(z[i])
            y_graph.append(y[i])
            values_graph.append(values[i])
            
    return y_graph, z_graph, values_graph


axx = fig.add_subplot(1, 3, 1)
axx.set_xlim([unique_y[0], unique_y[len(unique_y)-1]])
axx.set_ylim([unique_z[0], unique_z[len(unique_z)-1]])
axx.set_xlabel('y', fontsize=24)  
axx.set_ylabel('z', fontsize=24)

axx_slider = fig.add_axes([0.07, 0.02, 0.23, 0.03])
x_slider = Slider(
    ax=axx_slider,
    label='x-section',
    valmin=0,
    valmax=len(unique_x) - 1,
    valinit=0,
    valfmt="%i"
)

datax = getDatax(unique_x[int(len(unique_x) / 2)])

triangx = mtri.Triangulation(datax[0] , datax[1])

ymidx = np.array(datax[0])[triangx.triangles].mean(axis=1)
zmidx = np.array(datax[1])[triangx.triangles].mean(axis=1)
triangx.set_mask(getMaskx(ymidx, zmidx, unique_x[int(len(unique_x) / 2)]))

countourx = axx.tricontourf(triangx, datax[2])
axx.tricontour(countourx, colors='black')
colbarx = fig.colorbar(countourx, ax=axx)
axx.set_title('x = ' + str(unique_x[int(len(unique_x) / 2)]), fontsize=28)


def updatex(val):
    try:
        data = getDatay(unique_x[int(val)])
        
        triang = mtri.Triangulation(data[0] , data[1])
    
        ymid = np.array(data[0])[triang.triangles].mean(axis=1)
        zmid = np.array(data[1])[triang.triangles].mean(axis=1)
        triang.set_mask(getMaskx(ymid, zmid, unique_x[int(val)]))
    
        axx.collections.clear()    
        countour = axx.tricontourf(triang, data[2])
        axx.tricontour(countour, colors='black')
        global colbarx
        colbarx.remove()
        colbarx = fig.colorbar(countour, ax=axx)   
        axx.set_title('x = ' + str(unique_x[int(val)]), fontsize=28)
        plt.draw() 
    except:
        print("Error! x = " + str(unique_x[int(val)]))
        pass


x_slider.on_changed(updatex)

#--------------------------------------------------------------------------------------------------
def getMasky(xmid, zmid, uny):    
    mask = []
    for i in range(len(xmid)):
        flag = True
        p = [xmid[i], uny, zmid[i]]   
        for j in range(len(subdomains)):
            if pointInside(np.array(subdomains[j]), np.array(p)):      
                flag = False
                break
        mask.append(flag)
    return mask
    
def getDatay(yVal):
    x_graph = []
    z_graph = []
    values_graph = []
    
    for i in range(len(y)):
        if(y[i]==yVal):
            x_graph.append(x[i])
            z_graph.append(z[i])
            values_graph.append(values[i])
            
    return x_graph, z_graph, values_graph

axy = fig.add_subplot(1, 3, 2)
axy.set_xlim([unique_x[0], unique_x[len(unique_x)-1]])
axy.set_ylim([unique_z[0], unique_z[len(unique_z)-1]])
axy.set_xlabel('x', fontsize=24)  
axy.set_ylabel('z', fontsize=24)

axy_slider = fig.add_axes([0.38, 0.02, 0.23, 0.03])
y_slider = Slider(
    ax=axy_slider,
    label='y-section',
    valmin=0,
    valmax=len(unique_y) - 1,
    valinit=0,
    valfmt="%i"
)

datay = getDatay(unique_y[0])

triangy = mtri.Triangulation(datay[0] , datay[1])

xmidy = np.array(datay[0])[triangy.triangles].mean(axis=1)
ymidy = np.array(datay[1])[triangy.triangles].mean(axis=1)
triangy.set_mask(getMasky(xmidy, ymidy, unique_y[0]))

countoury = axy.tricontourf(triangy, datay[2])
axy.tricontour(countoury, colors='black')
colbary = fig.colorbar(countoury, ax=axy)
axy.set_title('y = ' + str(unique_y[0]), fontsize=28)


def updatey(val):
    data = getDatay(unique_y[int(val)])
     
    triang = mtri.Triangulation(data[0] , data[1])
    
    xmid = np.array(data[0])[triang.triangles].mean(axis=1)
    ymid = np.array(data[1])[triang.triangles].mean(axis=1)
    triang.set_mask(getMasky(xmid, ymid, unique_y[int(val)]))

    axy.collections.clear()    
    countour = axy.tricontourf(triang, data[2])
    axy.tricontour(countour, colors='black')
    global colbary
    colbary.remove()
    colbary = fig.colorbar(countour, ax=axy)   
    axy.set_title('y = ' + str(unique_y[int(val)]), fontsize=28)
    plt.draw() 


y_slider.on_changed(updatey)

#--------------------------------------------------------------------------------------------------
def getMaskz(xmid, ymid, unz):    
    mask = []
    for i in range(len(xmid)):
        flag = True
        p = [xmid[i], ymid[i], unz]   
        for j in range(len(subdomains)):
            if pointInside(np.array(subdomains[j]), np.array(p)):       
                flag = False
                break
        mask.append(flag)
    return mask
    
def getDataz(zVal):
    x_graph = []
    y_graph = []
    values_graph = []
    
    for i in range(len(z)):
        if(z[i]==zVal):
            x_graph.append(x[i])
            y_graph.append(y[i])
            values_graph.append(values[i])
            
    return x_graph, y_graph, values_graph

axz = fig.add_subplot(1, 3, 3)
axz.set_xlim([unique_x[0], unique_x[len(unique_x)-1]])
axz.set_ylim([unique_y[0], unique_y[len(unique_y)-1]])
axz.set_xlabel('x', fontsize=24)  
axz.set_ylabel('y', fontsize=24)

axz_slider = fig.add_axes([0.69, 0.02, 0.23, 0.03])
z_slider = Slider(
    ax=axz_slider,
    label='z-section',
    valmin=0,
    valmax=len(unique_z) - 1,
    valinit=0,
    valfmt="%i"
)


dataz = getDataz(unique_z[0])

triangz = mtri.Triangulation(dataz[0] , dataz[1])

xmidz = np.array(dataz[0])[triangz.triangles].mean(axis=1)
ymidz = np.array(dataz[1])[triangz.triangles].mean(axis=1)
triangz.set_mask(getMaskz(xmidz, ymidz, unique_z[0]))

countourz = axz.tricontourf(triangz, dataz[2])
axz.tricontour(countourz, colors='black')
colbarz = fig.colorbar(countourz, ax=axz)
axz.set_title('z = ' + str(unique_z[0]), fontsize=28)


def updatez(val):
    data = getDataz(unique_z[int(val)])
     
    triang = mtri.Triangulation(data[0] , data[1])
    
    xmid = np.array(data[0])[triang.triangles].mean(axis=1)
    ymid = np.array(data[1])[triang.triangles].mean(axis=1)
    triang.set_mask(getMaskz(xmid, ymid, unique_z[int(val)]))

    axz.collections.clear()    
    countour = axz.tricontourf(triang, data[2])
    axz.tricontour(countour, colors='black')
    global colbarz
    colbarz.remove()
    colbarz = fig.colorbar(countour, ax=axz)   
    axz.set_title('z = ' + str(unique_z[int(val)]), fontsize=28)
    plt.draw() 


z_slider.on_changed(updatez)

plt.show()
