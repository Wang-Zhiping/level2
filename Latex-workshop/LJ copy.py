import numpy as np 
import matplotlib.pyplot as plt
epsilon=-0.99
N=100
gamma=21.7
solve=[]
kk=1
col=['orange','y','forestgreen','royalblue','m','r']
ee=np.math.e
def fun(x):           
    f = np.math.sqrt(round((epsilon)-(vx(x)),10)) 
    return f

def T_2n(a, b, n, T_n):       
    if n<1:                
        print('n should larger than 1')
    h = (b - a)/n             
    sum_f = 0.            
    for k in range(0, n):
        sum_f = sum_f + fun(a + (k + 0.5)*h)
    T_2n = T_n/2. + sum_f*h/2.
    return T_2n

#Lennard-Jones
def inout(epsilon):
    epsilon=epsilon/4
    y2=(1-np.math.sqrt(1+4*epsilon))/2
    y1=(1+np.math.sqrt(1+4*epsilon))/2
    x1,x2=1.1,1.3
    for i in range(1,1000):
        x1=x1-(x1**-6-y1)/(-6*x1**-7)
        x2=x2-(x2**-6-y2)/(-6*x2**-7)

    solve.append(x1)
    solve.append(x2)
    return solve[0],solve[1]


def Romberg(a, b, err_min):
    kmax = 10
    tm = np.zeros(kmax,dtype = float)      
    tm1 = np.zeros(kmax,dtype = float)     
    tm[0] = 0.5*(b-a)*(fun(a) + fun(b))  
    #print(tm)
    err = 1.
    k = 0
    np.set_printoptions(precision = 9)
    while(err>err_min and k <kmax-2):  
        n = 2**k                      
        m = 1
        tm1[0] = T_2n(a, b, n, tm[0]) 
        while(err>err_min and m <= (k+1)):  
            tm1[m] = tm1[m-1]+(tm1[m-1]-(tm[m-1]))/(4.**m-1)
            result = tm1[m]
            err1 = abs(tm1[m]-tm[m-1])
            err2 = abs(tm1[m]-tm1[m-1])
            err = min(err1,err2)
            m = m+1
        tm = np.copy(tm1)
        k = k+1
        #print(tm)     
    return result

#fx=0
def f(n):
    fx=-(n+0.5)*np.math.pi
    fx=fx+Romberg(solve[0], solve[1], 1.e-12)*gamma
    return fx

def vx(x):
    return (x**(-12)-x**(-6))*4

def paint(ii):
	
    x1=inout(result[ii])
    n1=[vx(x1[0]),vx(x1[1])]
    plt.plot(x1,n1,color=col[i%6])
    solve.clear()

def paint2(i):
    e=result[i]/4
    #Lennard −Jones 势
    y2=(1-np.math.sqrt(1+4*e))/2
    y1=(1+np.math.sqrt(1+4*e))/2
    x1,x2=1.1,1.1
    for o in range(1,1000):
        x1=x1-(x1**-6-y1)/(-6*x1**-7)
        x2=x2-(x2**-6-y2)/(-6*x2**-7)
    xx=np.linspace(x1,x2,10000)
    #p=np.sqrt(np.abs(result[i]-vx(xx)))
    p=(np.abs(result[i]-vx(xx))**0.5)
    plt.plot(xx,p,color=col[i%6])
    plt.plot(xx,-p,color=col[i%6])
    


if __name__=="__main__":
    epsl=[]
    fnnn=[]
    result=[]
    for n in range(0,N):
        epsl=[]
        fnnn=[]
        while (1):
            
            solve=[]
            
            inout(epsilon)
            if (abs(f(n))<1.e-2):
                print(f(n),n,epsilon,solve)
                fnnn.append(f(n))
                epsl.append(epsilon)
            epsilon=epsilon+-epsilon/1000
            if epsilon>-0.0001:
                break
        print(fnnn)
        print(n,epsl)
        if (len(epsl)==0):
            print("MAX=",n-1)
            break
        mi=fnnn[0]
        for i in range(0,len(fnnn)):
            if abs(mi)>abs(fnnn[i]):
                mi=fnnn[i]
        for i in range(0,len(epsl)):
            if mi==fnnn[i]:
                result.append(epsl[i])
        print("r",result)
        epsilon=result[-1]
    print(result)
    #势能图
    fig = plt.figure(num=1)
    solve.clear()
    x=np.linspace(1,3,1000)
    v=4*(x**(-12)-x**(-6))
    plt.plot(x,v)
    for i in range(0,len(result)):
        paint(i)
    plt.title("D.L.J Potential Energy-D.r Diagram") 
    plt.ylabel("D-L.J Potential Energy(v)") 
    plt.xlabel("D.r(x)") 
    
    #相图
    fig = plt.figure(num=2)
    for i in range(0,len(result)):
        paint2(i)
    plt.title("Momentum-D.r Diagram") 
    plt.ylabel("Momentum") 
    plt.xlabel("D.r(x)") 

    
    fig = plt.figure(num=3)
    plt.subplot(2,  1,  1)  
    # 绘制第一个图像 

    x=np.linspace(1,3,1000)
    v=4*(x**(-12)-x**(-6))
    plt.plot(x,v)
    for i in range(0,len(result)):
        paint(i)
    plt.title("D.L.J Potential Energy-D.r Diagram") 
    plt.ylabel("D-L.J Potential Energy(v)") 
    plt.xlabel("D.r(x)") 
    plt.xlim(0.9,3.1)
    plt.ylim(-1,0)
    # 将第二个 subplot 激活，并绘制第二个图像
    plt.subplot(2,  1,  2) 
    for i in range(0,len(result)):
        paint2(i)
    plt.title("Momentum-D.r Diagram") 
    plt.ylabel("Momentum") 
    plt.xlabel("D.r(x)") 
    plt.xlim(0.9,3.1)
    # 展示图像
    plt.show()


