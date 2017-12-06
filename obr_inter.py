#!/bin/python3

# f(2) = 0
def f(x):
	return 1/(x-1) - 1

eps = 1e-4

x0 = 1.5
x1 = 10

count = 0

# while abs(x0 - x1) > eps and count < 10:
while f(x0)*f(x1) < 0 and count < 10:
	print("count: {}".format(count))
	print("x0: {}".format(x0))
	print("x1: {}".format(x1))
	x2 = x1 - (f(x1) * (x1 - x0))/(f(x1)-f(x0))
	print("x2: {}".format(x2))
	print()
	# x0 = x1
	x1 = x2

	count += 1