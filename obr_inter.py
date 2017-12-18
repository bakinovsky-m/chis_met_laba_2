#!/bin/python3

def f(x):
	# return 1/(x-1) - 1
	return x**2 + 2*x - 3

# def razd_razn(x_arr, n):
# 	x0 = x_arr[0]
# 	x1 = x_arr[1]
# 	if n <= 1:
# 		return (f(x1) - f(x0)) / (x1 - x0)
# 	else:
# 		# return razd_razn((x1, x2), n - 1) - razd_razn((x0, x1), n - 1) / x2 - x0
# 		return (razd_razn(x_arr[1:], n - 1) - razd_razn(x_arr, n - 1)) / x_arr[n] - x0

# def razd_razn(x_arr, n):
# 	if n == 0:
# 		return 1
# 	res = 0
# 	for j in range(n):
# 		znam = 1
# 		for i in range(n):
# 			if i != j:
# 				znam *= x_arr[j] - x_arr[i]
# 		res += f(x_arr[j]) / znam
# 	return res
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



def razd_razn_1(x0, x1):
	return (f(x1) - f(x0)) / (x1 - x0)
def razd_razn_2(x0, x1, x2):
	return (razd_razn_1(x1, x2) - razd_razn_1(x0, x1)) / (x2 - x0)
def razd_razn_3(x0, x1,x2,x3):
	chis = razd_razn_2(x1,x2,x3) - razd_razn_2(x0,x1,x2)
	znam = x3 - x0
	return chis / znam
def razd_razn_4(x0, x1,x2,x3,x4):
	chis = razd_razn_3(x1,x2,x3,x4) - razd_razn_3(x0,x1,x2,x3)
	znam = x4 - x0
	return chis / znam
def nnewton(x, x0, x1, x2, x3, x4):
	p1 = f(x0)
	p2 = razd_razn_1(x0, x1) * (x - x0)
	p3 = razd_razn_2(x0,x1,x2)*(x-x0)*(x-x1)
	p4 = razd_razn_3(x0,x1,x2,x3)*(x-x0)*(x-x1)*(x-x2)
	p5 = razd_razn_4(x0,x1,x2,x3,x4)*(x-x0)*(x-x1)*(x-x2)*(x-x3)
	return p1 + p2 + p3 + p4 + p5

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

# steps = 20
steps = len(range(-5, 5))
# for i in range(-5, 5):
# 	# y_arr.append(f(i)/steps)
# 	# y_arr.append(i/steps)
# 	y_arr.append(i)
# 	# print(i/steps)
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

while count < 10:
	print("f(x): {}".format(f(x0)))
	print("newton: {}".format(newton(f(x0), y_arr, new_n)))
	print("nnewton: {}".format(nnewton(f(x0), y_arr[0],y_arr[1],y_arr[2], y_arr[3], y_arr[4])))

	new_n += 1
	count += 1
