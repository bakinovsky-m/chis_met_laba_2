#!/bin/python3

def f(x):
	return x**2 + 2*x - 3

def razd_razn(x_arr, n):
	if n == 0:
		return 1
	res = 0
	if n == 1:
		chis = f(x_arr[1]) - f(x_arr[0])
	else:
		chis = razd_razn(x_arr[1:], n - 1) - razd_razn(x_arr[:-1], n - 1)
	znam = x_arr[n] - x_arr[0]
	res = chis / znam
	return res

def newton(x, x_arr, n):
	res = f(x_arr[0])
	sk = 1
	for i in range(1, n):
		rr = razd_razn(x_arr, i)
		sk *= x - x_arr[i-1]
		res += rr * sk
	return res

eps = 1e-4

y_arr = []
x0 = -3

count = 0

y_arr.append(-4)
y_arr.append(-3.9)
y_arr.append(-3.5)
y_arr.append(-3.2)
y_arr.append(-3.1)
y_arr.append(-2.8)
y_arr.append(-2.7)
y_arr.append(-2.3)
y_arr.append(-2)
y_arr.append(-1.9)
y_arr.append(-1.8)
qwe = y_arr.copy()
qwe.reverse()

new_n = 1

done = False
while not done:
	print("newton: {}".format(newton(f(x0), y_arr, new_n)))
	if abs(newton(f(x0), y_arr, new_n) - x0) < eps:
		done = True

	new_n += 1
	count += 1
