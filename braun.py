import matplotlib.pyplot as plt

def f(x, y):
	return x**2 - 2*(y**2) - x*y + 2*x - y + 1

def g(x, y):
	return 2*(x**2) - y**2 + x*y + 3*y - 5

# proizvodnay f po x
def f_x(x, y):
	return 2*x -y + 2

def f_y(x, y):
	return -4*y - x - 1

def g_x(x, y):
	return 4*x + y

def g_y(x, y):
	return -2*y + x + 3

eps = 1e-8

xk = 0
yk = 0
x_ = 0
qk = 1
pk = 1

count = 0

xk_arr = []
yk_arr = []

while abs(max(qk, pk)) > eps and count < 10:
	count += 1

	x_ = xk - (f(xk, yk)/f_x(xk, yk))
	qk = g(x_, yk)*f_x(xk, yk) / ( f_x(xk, yk)*g_y(xk, yk) - f_y(xk, yk)*g_x(x_, yk) )
	pk = (f(xk, yk) - qk*f_y(xk, yk))/f_x(xk, yk)
	xk = xk - pk
	yk = yk - qk

	xk_arr.append(xk)
	yk_arr.append(yk)

print("count: {}".format(count))
print("xk: {}".format(xk))
print("yk: {}".format(yk))

plt.plot(xk_arr, yk_arr, '-o')
plt.show()
