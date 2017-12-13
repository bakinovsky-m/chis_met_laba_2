#!/bin/python3

import matplotlib.pyplot as plt

def f(x):
	return 1/(x-1) - 1

eps = 1e-4

x0 = 1.5
x1 = 10

count = 0

x_arr = []

while f(x0)*f(x1) < 0 and count < 10:
	x2 = x1 - (f(x1) * (x1 - x0))/(f(x1)-f(x0))
	x1 = x2

	count += 1

	x_arr.append(x1)

print("x: {}".format(x1))

plt.plot(x_arr, '-o')
plt.show()
