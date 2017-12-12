#!/bin/python3
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def grad_x(x, y, z):
	return 2 * (0.1 - pow(x, 2) + 2 * y*z - x)*(-2 * x - 1) - 6 * (-0.2 + pow(y, 2) - 3 * x*z - y)*z - 4 * (0.3 - pow(z, 2) - 2 * x*y - z)*y

def grad_y(x, y, z):
	return 4 * (0.1 - pow(x, 2) + 2 * y*z - x)*z + 2 * (-0.2 + pow(y, 2) - 3 * x*z - y)*(2 * y - 1) - 4 * (0.3 - pow(z, 2) - 2 * x*y - z)*x

def grad_z(x, y, z):
	return 4 * (0.1 - pow(x, 2) + 2 * y*z - x)*y - 6 * (-0.2 + pow(y, 2) - 3 * x*z - y)*x + 2 * (0.3 - pow(z, 2) - 2 * x*y - z)*(-2 * z - 1);

xk = 2
xk_old = 0
yk = 2
yk_old = 0
zk = 2
zk_old = 0

lambd = math.exp(-5)
sigm = 0.1
eps = 1e-4

xk_old = xk
yk_old = yk
zk_old = zk
xk = xk_old - lambd * grad_x(xk_old, yk_old, zk_old)
yk = yk_old - lambd * grad_y(xk_old, yk_old, zk_old)
zk = zk_old - lambd * grad_z(xk_old, yk_old, zk_old)


count = 0

xk_arr = []
yk_arr = []
zk_arr = []

while max(abs(xk - xk_old), abs(yk - yk_old), abs(zk - zk_old)) > eps:
	count += 1
	
	xk_old = xk
	yk_old = yk
	zk_old = zk

	xk = xk_old - lambd * grad_x(xk_old, yk_old, zk_old)
	yk = yk_old - lambd * grad_y(xk_old, yk_old, zk_old)
	zk = zk_old - lambd * grad_z(xk_old, yk_old, zk_old)

	xk_arr.append(xk)
	yk_arr.append(yk)
	zk_arr.append(zk)


print("xk: {}".format(xk))
print("yk: {}".format(yk))
print("zk: {}".format(zk))

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(xs=xk_arr, ys=yk_arr, zs=zk_arr, marker='o')
plt.show()
plt.plot(xk_arr, yk_arr, '-o')
plt.show()
plt.plot(xk_arr, zk_arr, '-o')
plt.show()
plt.plot(yk_arr, zk_arr, '-o')
plt.show()
