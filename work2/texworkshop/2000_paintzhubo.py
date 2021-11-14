
    
import numpy as np
import matplotlib.pyplot as plt
import time
import copy
N = 2000
col = ['c','darkmagenta','orangE','royalblue','teal','r']

A=np.eye(N)

for i in range(N):
    A[i][i]=-2
    if (0<i):
        A[i][i-1]=1
    if (N-1>i):
        A[i][i+1]=1
        
print(A)

def Jacobi(A,eps,Eigv,detail):
    for ii in range(1,1000*len(A)):
        simga = 0
        apq = 0
        tan2 = 0
        N=len(A)
        for i in range(N):
            for j in range(i+1,N):
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
        cos2 = np.math.cos(2*theta)
        sin2 = np.math.sin(2*theta)
        Q=np.eye(N)
        Q[p][p]=cos
        Q[q][q]=cos
        Q[p][q]=sin
        Q[q][p]=-sin
        # A=Product(Product(Trans(Q),A),Q)
        app = float(A[p][p])
        aqq = float(A[q][q])
        apq = float(A[p][q])
        for j in range(len(A)):
            if (j!=q)and(j!=p):
                apj = float(A[p][j])
                aqj = float(A[q][j])
                A[p][j] = cos*apj - sin*aqj
                A[q][j] = sin * apj + cos * aqj
                A[j][p] = A[p][j]
                A[j][q] = A[q][j]

        A[p][p] = aqq*sin*sin + app*cos*cos -apq *sin2
        A[q][q] = app*sin*sin + aqq*cos*cos +apq *sin2
        A[p][q] = float(apq)*cos2 + (app-aqq)/2*(sin2)
        A[q][p] = A[p][q]

        #E = copy.copy(Product(Eigv,Q))
        for i in range(N):
            aip = Eigv[i][p]
            aiq = Eigv[i][q]
            Eigv[i][p] = cos*aip-sin*aiq
            Eigv[i][q]= sin*aip+cos*aiq
        #E=Eigv@Q
        # for i in range(N):
        #     for j in range(N):
        #         Eigv[i][j] = E[i][j]
        #print(np.round(A,10))
        if (detail):
            print("times={}\n{}".format(ii,np.round(A,4)))
            print(np.round(Eigv,4))
        sum = 0
        for i in range(N):
            for j in range(N):
                if (i != j):
                    sum += abs(A[i][j])
        if (sum < eps):
            break
    for i in range(N):
        value[0][i]=(A[i][i])
    for i in range(N):
        for j in range(N):
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
    Eigv=np.eye(N)
    B=np.eye(len(A))
    G=np.eye(len(A))
    yy = np.zeros(N)
    xx = np.zeros(N)
    value = np.zeros([1,N])
    time1 = time.time()
    #Jacobi(A.copy(),1.e-10,Eigv,False)
    time2 = time.time()

    time3 = time.time()
    (value,Eigv)=np.linalg.eig(A)
    time4 = time.time()


    print("\t\t\t\t\t\t\tEigenvalue")
    print(np.round(value,10))
    print("\t\t\t\t\t\t\tEigenvector")
    print(np.round(Eigv,4))

    # value = np.sort(value[0])

    # print(a-value)


    print("\n\n\n")
    dotx0 = np.zeros([1,N])
    x0 = np.zeros([1,N])
    x0[0][0] = 1
    x = Product(Trans(Eigv),x0)
    x = (Trans(x))[0].copy()
    print("gamma1",x)
    gamma1 = x
    x = Product(Trans(Eigv),dotx0)
    x = (Trans(x))[0]
    print("wgamma2",x)
    wgamma2 = (x)

    
    dotx0 = np.zeros([1,N])
    x0 = np.zeros([1,N])
    x0[0][0] = 1
    x = ((Eigv.T)@x0[0])
    x = x
    print("gamma1",x)
    gamma1 = x
    x = ((Eigv.T)@dotx0[0])
    x = x
    print("wgamma2",x)
    wgamma2 = (x)
    omega = np.sqrt(-value).copy()
    print("omega=",omega)
    plt.figure(figsize=(15, 8),dpi=100)
    # x = np.linspace(0,25, 1000) 
    # for i in range(len(Eigv)):
    #         y = 0
    #         for j in range(len(Eigv)):
    #             y += Eigv[i][j]*(gamma1[j]*np.cos(omega[j]*x)+wgamma2[j]/omega[j]*np.sin(omega[j]*x))
    #         plt.plot(x, y,color = col[i%5],label = 'm_{}'.format(i))
    #     #plt.legend(['m_1','m_2','m_3','m_4','m_5'] ,loc='upper right' )
    # plt.ylim(-2,2)
    
    xx = np.linspace(0,N, N)
    yy = []
    x=10000
    y = 0
    yy.clear()
    for j in range(len(Eigv)):
        y += Eigv[i][j]*(gamma1[j]*np.cos(omega[j]*x)+wgamma2[j]/omega[j]*np.sin(omega[j]*x))
        yy.append(y)
        y = 0
    plt.plot(xx, yy,alpha=0.2)
    plt.scatter(xx, yy,color="",edgecolor=col[0])
        #plt.legend(['m_1','m_2','m_3','m_4','m_5'] ,loc='upper right' )
    plt.show()


