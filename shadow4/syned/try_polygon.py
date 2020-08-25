from time import time
import numpy as np
import matplotlib.path as mpltPath

# regular polygon for testing
lenpoly = 15
polygon = [[np.sin(x)+0.5,np.cos(x)+0.5] for x in np.linspace(0,2*np.pi,lenpoly)[:-1]]
xp = np.zeros(len(polygon)+1)
yp = np.zeros_like(xp)

print(len(polygon),polygon)
for i in range(len(polygon)):
    print(i,polygon[i])
    xp[i] = polygon[i][0]
    yp[i] = polygon[i][1]
xp[-1] = xp[0]
yp[-1] = yp[0]

from srxraylib.plot.gol import plot, set_qt
set_qt()
plot(xp,yp)

def ray_tracing(x,y,poly):
    n = len(poly)
    inside = False
    p2x = 0.0
    p2y = 0.0
    xints = 0.0
    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside

def ray_tracing_numpy(x,y,poly):
    n = len(poly)
    inside = np.zeros(len(x),np.bool_)
    p2x = 0.0
    p2y = 0.0
    xints = 0.0
    p1x,p1y = poly[0]
    print(">>>>--",p1x,p1y,poly[0])
    for i in range(n+1):
        p2x,p2y = poly[i % n]

        idx = np.nonzero((y > min(p1y,p2y)) & (y <= max(p1y,p2y)) & (x <= max(p1x,p2x)))[0]

        print(">>>>", i, idx, p2x, p2y, poly[i % n])
        if len(idx > 0):
            idx = idx[0]
            if p1y != p2y:
                xints = (y[idx]-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
            if p1x == p2x:
                inside[idx] = ~inside[idx]
            else:
                idxx = idx[x[idx] <= xints]
                inside[idxx] = ~inside[idxx]

        p1x,p1y = p2x,p2y
    return inside

# def ray_tracing_mult(x,y,poly):
#     return [ray_tracing(xi, yi, poly[:-1,:]) for xi,yi in zip(x,y)]

print(ray_tracing(0.5,0.5,polygon))
print(ray_tracing(0,2,polygon))
print(ray_tracing_numpy([0.5],[0.5],polygon))