#!/bin/python3
import math

def grad_x(x, y):
	return 2*(10*(x**3) + 3*(x**2)*(y+2) - 2*x*(3*(y**2) - 3*y+7) + (y**3) - 8*y +2)

def grad_y(x, y):
	return 2*((x**3) +(x**2)*(3-6*y) + x*(3*(y**2)-8) + 10*(y**3) - 3*y + 16*y - 16)

xk = 0
xk_old = 0
yk = 0
yk_old = 0

# lambd = [e^-4; e^-9]
lambd = math.exp(-5)
sigm = 0.1
eps = 1e-4

xk_old = xk
yk_old = yk
xk = xk_old - lambd * grad_x(xk_old, yk_old)
yk = yk_old - lambd * grad_y(xk_old, yk_old)


count = 1

while max(abs(xk - xk_old), abs(yk - yk_old)) > eps and count < 5000:
	count += 1
	print("iter: {}".format(count))
	print("xk: {}".format(xk))
	print("yk: {}".format(yk))
	print("xk_old: {}".format(xk_old))
	print("yk_old: {}".format(yk_old))
	print("grad_x: {}".format(grad_x(xk_old, yk_old)))
	
	xk_old = xk
	yk_old = yk

	xk = xk_old - lambd * grad_x(xk_old, yk_old)
	yk = yk_old - lambd * grad_y(xk_old, yk_old)