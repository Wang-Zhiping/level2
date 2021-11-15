import numpy as np
import matplotlib.pyplot as plt

A=np.array([[   -2  ,   1  ,   0    ,   0    ,    0],
            [    1  ,  -2  ,   1    ,   0    ,    0],
            [    0  ,   1  ,  -2    ,   1    ,    0],
            [    0  ,   0  ,   1    ,  -2    ,    1],
            [    0  ,   0  ,   0    ,   1    ,   -2]])

def getr(ii,A):
    r = 0
    for i in range(len(A)):
        if (ii != i):
            r += abs(A[ii][i])
    return r

def p_circle(r,ii,fig):
    a = A[ii][ii]
    b = 0
    theta = np.arange(0, 2*np.pi, 0.01)
    x = a + r * np.cos(theta)
    y = b + r * np.sin(theta)
    plt.plot(x, y,alpha = 0.4,linestyle=':',label = '{}'.format(i),linewidth = '{}'.format(len(A)-i+1),marker=',') 

if __name__ == "__main__":

    Gershgorin = plt.figure()
    for i in range(len(A)):
        r = getr(i,A)
        p_circle(r,i,Gershgorin)
    plt.axis('equal')
    plt.legend(['G-circle of 1',"G-circle of 2","G-circle of 3","G-circle of 4","G-circle of 5"] ,loc='upper right' )
    plt.title('Gershgorin circle')
    plt.show()

