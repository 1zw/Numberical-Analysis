# -*- coding: utf-8 -*-
from cProfile import label
import numpy as np
import matplotlib.pyplot as plt
# Fit the data by cubic and fourth-degree polynomials
data = np.array([[0.0, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0],
                 [1.0, 0.41, 0.50, 0.61, 0.91, 2.02, 2.4]])
degree = [3,4]
x = np.linspace(np.min(data[0,:]),np.max(data[0,:]),100)
def fitted_y(params,deg,side):
        global x
        y = np.zeros(100)
        for i,param in enumerate(params):
            if side=="inside":
                y += param*x**(deg-i)
            else:
                y += param*x**i
        return y
def inside():
    for deg in degree:
        params = np.polyfit(data[0,:],data[1,:],deg = deg)
        # print(params)
        y = fitted_y(params,deg,side="inside")
        plt.plot(x,y,label=str(deg)+"_polynomials_fitted")
    plt.scatter(data[0,:],data[1,:],label="points")
    # plt.plot(data[0,:],data[1,:],label = 'origin')
    plt.legend()
    plt.savefig("./img/inside.png")
    plt.close()
def outside():
    
    for deg in degree:
        fai_origin = np.zeros((deg+1,data.shape[1]))
        g_orgin = data[1,:]
        for i in range(0,deg+1):
            fai_origin[i] = data[0,:]**i
            # g_orgin[i] = data[1,:]**i
        # print (fai_origin)
        G = np.zeros((deg,deg))
        g = np.zeros((deg,1))
        for i in range(G.shape[0]):
            for j in range(G.shape[1]):
                G[i][j] = np.dot(fai_origin[i],fai_origin[j])
            g[i][0] = np.dot(fai_origin[i],g_orgin)
        params = np.dot(np.linalg.inv(G),g)
        params.reshape(1,-1)[0]
        print(params)
        y_fitted = fitted_y(params.reshape(1,-1)[0],deg,side="outside")
        plt.plot(x,y_fitted,label=str(deg)+"_polynomials_fitted")
    plt.scatter(data[0,:],data[1,:],label="points")
    plt.legend()
    plt.savefig("./img/outside.png")
    plt.close()     
if __name__ == '__main__':
    inside() # python 内置多项式拟合
    outside()# outside手写拟合