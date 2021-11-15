
    
import numpy as np
import matplotlib.pyplot as plt

col = ['c','darkmagenta','orangE','royalblue','teal']

A=np.eye(8)
for i in range(len(A)):
    A[i][i]=-2
    if (0<i):
        A[i][i-1]=1
    if (4>i):
        A[i][i+1]=1


def Jacobi(A,eps,Eigv,detail):
    for ii in range(1,100):
        simga = 0
        apq = 0
        tan2 = 0
        for i in range(len(A)):
            for j in range(len(A)):
                if (i!=j):
                    simga +=A[i][j]             #step1
                    if (abs(A[i][j])>abs(apq)):
                        apq = A[i][j]           #step2
                        (p,q) = (i,j)
        if ((A[q][q]-A[p][p])!=0):
            tan2 = 2*A[p][q]/(A[q][q]-A[p][p])
            theta = np.math.atan(tan2) / 2
        else:
            theta = np.math.pi/4
        cos = np.math.cos(theta)
        sin = np.math.sin(theta)
        tan = np.math.tan(theta)
        Q=np.eye(len(A))
        Q[p][p]=cos
        Q[q][q]=cos
        Q[p][q]=sin
        Q[q][p]=-sin
        A=Product(Product(Trans(Q),A),Q)
        # for i in range(len(A)):
        #     if ((i !=p)and(i != q)):
        #         A[i][p] = A[p][i]*cos - A[q][i]*sin 
        #         A[i][q] = A[p][i]*sin + A[q][i]*cos
        #         A[p][i] = A[i][p]
        #         A[q][i] = A[i][q]
        # app = A[p][p]
        # aqq = A[q][q]
        # A[p][p] = sin*sin*app + cos*cos*aqq + 2*sin*cos*app
        # A[p][p] = sin*sin*aqq + cos*cos*app - 2*sin*cos*app
        # A[q][q] = A[q][q] + tan*A[p][q]
        E = Product(Eigv,Q)
        for i in range(len(A)):
            for j in range(len(A)):
                Eigv[i][j] = E[i][j]
        #print(np.round(A,10))
        if (detail):
            print("times={}\n{}".format(ii,np.round(A,4)))
            print(np.round(Eigv,4))
        sum = 0
        for i in range(len(A)):
            for j in range(len(A)):
                if (i != j):
                    sum += abs(A[i][j])
        if (sum < eps):
            break
    for i in range(len(A)):
        value[0][i]=(A[i][i])
    for i in range(len(A)):
        for j in range(len(A)):
            B[i][j] = A[i][j]
    #print(np.round(Eigv,3))

    return (value)

def Product(A, B):
    C = np.zeros([len(A),len(B[0])])
    for i in range(len(A)):
        for j in range(len(B[0])):
            for k in range(len(B)):
                C[i][j] += A[i][k]*B[k][j]
    return C

def Trans(A):
    C=np.zeros((len(A[0]),len(A)))
    for i in range (len(A[0])):
        for k in range(len(A)):
            C[i][k]=A[k][i]
    return C

if  __name__ == "__main__":
    Eigv=np.eye(len(A))
    B=np.eye(len(A))
    G=np.eye(len(A))
    yy = np.zeros(5)
    xx = np.zeros(5)
    value = np.zeros([1,len(A)])
    Jacobi(A,1.e-10,Eigv,False)
    print("\t\t\t\t\t\t\tEigenvalue")
    print(np.round(value,10))
    print("\t\t\t\t\t\t\tEigenvector")
    print(np.round(Eigv,6))
    omega = np.sqrt(-value)[0].copy()

    print("\n\n\n")
    
    dotx0 = np.zeros([1,len(A)])
    x0 = np.zeros([1,len(A)])
    x0[0][0] = 1
    x = Product(Trans(Eigv),x0)
    x = (Trans(x))[0].copy()
    print("gamma1",x)
    gamma1 = x
    x = Product(Trans(Eigv),dotx0)
    x = (Trans(x))[0]
    print("wgamma2",x)
    wgamma2 = (x)
    # plt.figure(figsize=(15, 5),dpi=100)
    # x = np.linspace(0,25, 1000) 
    # for i in range(len(Eigv)):
    #     y = 0
    #     for j in range(len(Eigv)):
    #         y += Eigv[i][j]*(gamma1[j]*np.cos(omega[j]*x)+wgamma2[j]/omega[j]*np.sin(omega[j]*x))
    #     plt.plot(x, y,color = col[i%5],label = 'm_{}'.format(i))
    # plt.legend(['m_1','m_2','m_3','m_4','m_5'] ,loc='upper right' )
    # plt.ylim(-2,2)
    # plt.show()
