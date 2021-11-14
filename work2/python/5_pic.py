
    
import numpy as np
import matplotlib.pyplot as plt

col = ['olive','darkmagenta','orangE','royalblue','teal']

A=np.array([[   -2  ,   1  ,   0    ,   0    ,    0],
            [    1  ,  -2  ,   1    ,   0    ,    0],
            [    0  ,   1  ,  -2    ,   1    ,    0],
            [    0  ,   0  ,   1    ,  -2    ,    1],
            [    0  ,   0  ,   0    ,   1    ,   -2]])

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
    Eigv=np.eye(5)
    B=np.eye(len(A))
    G=np.eye(len(A))
    yy = np.zeros(5)
    xx = np.zeros(5)
    value = np.zeros([1,5])
    Jacobi(A.copy(),1.e-15,Eigv,False)
    print("\t\t\t\t\t\t\tEigenvalue")
    print(np.round(value,10))
    print("\t\t\t\t\t\t\tEigenvector")
    print(np.round(Eigv,10))
    omega = np.sqrt(-value)[0].copy()

    print("\n\n\n")
    
    dotx0 = np.array([[1,-1, 0, 1,-1]])
    x0 = np.array([[-1, 1, 0, -1, 1]])
    x = Product(x0,Eigv.copy())
    print("gamma1",x)
    gamma1 = x[0]
    x = Product(dotx0,Eigv.copy())
    print("wgamma2",x)
    wgamma2 = (x[0])
    plt.figure(figsize=(15, 15),dpi=100)
    plt.subplot(2,  1,  1)  
    x = np.linspace(0,25, 1000) 
    for i in range(len(Eigv)):
        y = 0

        for j in range(len(Eigv)):
            y += (Eigv[i][j]*gamma1[j]*np.cos(omega[j]*x)+
            Eigv[i][j]*wgamma2[j]/omega[j]*np.sin(omega[j]*x))

        plt.plot(x, y,color = col[i%5],label = 'm_{}'.format(i))
        # plt.scatter(x,y,color='r',marker='x')
    
    plt.legend(['m_1','m_2','m_3','m_4','m_5'] ,loc='upper right' )
    plt.ylim(-2,2)
    plt.subplot(2,  1,  2)
    xxx=[]
    yyy=[]
    col = ['c','darkmagenta','orangE','royalblue','teal']
    for j in range(len(A)):
        xxx.clear()
        yyy.clear()
        for i in range(len(A)):
            xxx.append(i+1)
            yyy.append(Eigv[i][j])
        plt.plot(xxx, yyy,color = col[j%5],alpha=0.5,linestyle='--',label = 'modle_{}'.format(j))
        plt.scatter(xxx, yyy,marker='o',color="",edgecolors = col[j%5])
        plt.ylim(-1,1)
        plt.legend(['modle_1','modle_2','modle_3','modle_4','modle_5'] ,loc='upper right' )
    plt.savefig('plot1.jpg')
    plt.show()
    print(omega)
