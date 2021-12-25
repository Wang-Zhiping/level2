#include <stdio.h>
#include <math.h>

#define M (2000)
#define Mu (2)
#define Lambda (0.5)
#define Epsilon (0)
#define Omega (0)


double h(double t){
    return Epsilon*sin(Omega*t);
}

double f1(double t,double n,double m){
    return n*(Lambda-(n+(m/(Mu+h(t))))/M);
}
double f2(double t,double n,double m){
    return m*(Lambda-(m+n/(Mu-h(t)))/M);
}

void RK_4(double t,double n,double m,double tm,double h,FILE *fp){
    double k11,k12,k13,k14,k21,k22,k23,k24,fai = 0;
    while(t <= tm ){
        // printf("%f,%f,%f\n",t,n,m);
        fai = Mu*(n*n+m*m)/2/M/M-Lambda*Mu/M*(n+m)+n*m/M/M;
        fwrite(&t, sizeof(double) , 1, fp );
        fwrite(&n, sizeof(double) , 1, fp );
        fwrite(&m, sizeof(double) , 1, fp );
        fwrite(&fai, sizeof(double) , 1, fp );
        k11 = f1(t,n,m);
        k21 = f2(t,n,m);
        k12 = f1(t+h/2,n+h/2*k11,m+h/2*k21);
        k22 = f2(t+h/2,n+h/2*k11,m+h/2*k21);
        k13 = f1(t+h/2,n+h/2*k12,m+h/2*k22);
        k23 = f2(t+h/2,n+h/2*k12,m+h/2*k22);
        k14 = f1(t+h  ,n+h*k13  ,m+h*k23);
        k24 = f2(t+h  ,n+h*k13  ,m+h*k23);
        n += h/6*(k11+2*k12+2*k13+k14);
        m += h/6*(k21+2*k22+2*k23+k24);
        t += h;

    }
    
}


int main(){
    int N = 240;
    FILE *fp;
    fp = fopen( "fig1_rk4.bin" , "w" );
    int i,j;
    for (i = 0;i<1200;i=i+N){
        RK_4(1,50,i,200,0.02,fp);
        RK_4(1,1200,i,200,0.02,fp);
        RK_4(1,i,50,200,0.02,fp);
        RK_4(1,i,1200,200,0.02,fp);
    }
    RK_4(1,90,50,200,0.02,fp);
    RK_4(1,100,150,200,0.02,fp);
    RK_4(1,1200,1200,200,0.02,fp);
    printf("\n%d\n",2+(i/N)*(4));
    fclose(fp);

    return 0;
}