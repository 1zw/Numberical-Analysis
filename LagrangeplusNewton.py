# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 00:52:53

@author: wily_elite
"""
import matplotlib.pyplot as plt
import numpy as np
def Lagrange_interpolation():
    # 1. lagrange interpolation to verify Runge phenomenon  
    x = np.arange(-5,5,0.1)
    y = 1/(1+x**2)
    x_la = np.arange(-5,5,1)
    y_la_x = 1/(1+x_la**2)
    px = 0
    for k in range(10):
        lc = 1
        for j in range(10):
            if k!=j:
                lc *= (x-x_la[j])/(x_la[k]-x_la[j])
        px += lc*y_la_x[k]
    
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    plt.plot(x,y)
    plt.scatter(x_la,y_la_x,marker='.',color='red',s=100)
    plt.plot(x,px)
    plt.ylim(0,1)
    plt.legend(['原函数图像','拉格朗日插值图像','插值点'])
    plt.show()
def Newton_interpolation():
    # 2. Newton interpolation code
    x_newton = np.array([1,4,9])
    y_newton = np.array([1,2,3])
    # ask for f(5),obviously,the answer is sqrt(5)
    y_newton_div= dict()
    for i in range(1,len(x_newton)):
        hi = []
        for k in range(len(x_newton)):
            if i+k <= len(x_newton)-1:
                hi.append(x_newton[i+k] - x_newton[k])
        fi = []
        if y_newton_div == dict():       
            for k in range(len(y_newton)-1):
                fi.append(y_newton[k+1] - y_newton[k])
        else:
            for z in range(len(y_newton_div[i-1])-1):
                fi.append(y_newton_div[i-1][z+1]-y_newton_div[i-1][z])
        y_newton_div[i] = [yi/xi for xi,yi in zip(hi,fi)]
    y_newton_div[0] = list(y_newton)
    x = 5
    y_inter_newton = 0
    for i in range(len(x_newton)):
        wi = 1
        for j in range(i):
            wi *= (x-x_newton[j])
        y_inter_newton += wi*y_newton_div[i][0]
    print(y_inter_newton)
if __name__ == '__main__':
    Lagrange_interpolation()
    Newton_interpolation()