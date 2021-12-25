#include <stdio.h>
#include <math.h>

#define M (2000)
#define Lambda (0.5)
#define Epsilon (0.01)
#define L1 (0.6)
#define L2 (0.1)



double hh(double t,double Omega){
    return Epsilon*sin(Omega*t);
}

double f(double Omega,double t,double x,double mu){
    double hp = hh(t,Omega)*Lambda*Lambda/(mu+1)/(mu+1);
    return sqrt(2)*hp-(mu-1)*x*Lambda/(mu+1);
}

void figaa(double mu,double Omega,double t,double x,double tm,double h,FILE* fp){
    double k1,k2,k3= 0;
    double tau =  (mu+1)/Lambda/(mu-1);
    double max = 0;
    double num = 0;
    double chid =0;
    while(t <= tm ){
        k1 = f(Omega,t,x,mu);
        k2 = f(Omega,t+0.5*h , x + 0.5*h*k1,mu);
        k3 = f(Omega,t+h , x -h*k1 + h*2*k2,mu);
        x += h*(k1+4*k2+k3)/6;
        t += h;
        chid = M*x/Epsilon/sqrt(1+Omega*Omega*tau*tau);
        if ((t>tm*0.995)&&(f(Omega,t,x,mu)*f(Omega,t+h,x+ h * f(Omega,t,x,mu),mu)<0)) {
            max += fabs(chid);
            num++;
            }
        chid = max/num;
    }

    fwrite(&Omega, sizeof(double) , 1, fp );
    fwrite(&chid, sizeof(double) , 1, fp );
}

void figbb(double mu,double Omega,double t,double x,double tm,double h,FILE* fp){
    double k1,k2,k3= 0;
    double tau =  (mu+1)/Lambda/(mu-1);
    double max = 0;
    double num = 0;
    double chid =0;
    while(t <= tm ){
        k1 = f(Omega,t,x,mu);
        k2 = f(Omega,t+0.5*h , x + 0.5*h*k1,mu);
        k3 = f(Omega,t+h , x -h*k1 + h*2*k2,mu);
        x += h*(k1+4*k2+k3)/6;
        t += h;
        chid = M*x/Epsilon/sqrt(1+Omega*Omega*tau*tau);
        if ((t>tm*0.999)&&(f(Omega,t,x,mu)*f(Omega,t+h,x+ h * f(Omega,t,x,mu),mu)<0)) {
            max += fabs(chid);
            num++;
            }
        chid = max/num;
    }

    fwrite(&mu, sizeof(double) , 1, fp );
    fwrite(&chid, sizeof(double) , 1, fp );
}

void figa(double Mu ,FILE *fp){
    double tau =  (Mu+1)/(Lambda*(Mu-1));
    for (int i =1;i<200;i++){
        double omega = i/1000.;
        double xs = (sqrt(2)*Lambda*Lambda*M/(Mu+1)/(Mu+1)* tau)/(1+omega*omega*tau*tau);
        fwrite(&omega, sizeof(double) , 1, fp );
        fwrite(&xs, sizeof(double) , 1, fp );
    }
    fwrite(&Mu, sizeof(double) , 1, fp );
    fwrite(&Mu, sizeof(double) , 1, fp );

}

void figb(double omega ,FILE *fp){

    for (int i =100;i<299;i++){
        double Mu = i/100.;
        double tau =  (Mu+1)/(Lambda*(Mu-1));
        double xs = (sqrt(2)*Lambda*Lambda*M/(Mu+1)/(Mu+1)* tau)/(1+omega*omega*tau*tau);
        fwrite(&Mu, sizeof(double) , 1, fp );
        fwrite(&xs, sizeof(double) , 1, fp );
    }
    fwrite(&omega, sizeof(double) , 1, fp );
    fwrite(&omega, sizeof(double) , 1, fp );
}




int main(){
    int N = 240;
    FILE *fp;
    fp = fopen( "fig2.bin" , "w" );
    figa(1.75,fp);
    figa(2,fp);
    figa(2.25,fp);
    figb(0.05,fp);
    figb(0.1,fp);
    figb(0.2,fp);



    for (int j = 0;j<3;j++){
        double mu = 1.75+0.25*j;
    for (int i = 0;i<10;i++){
        double omega = 0.01+i*0.02;
        figaa(mu,omega,0,1,200000,0.1,fp);
    }}

    for (int j = 0;j<3;j++){
        double omega;
        if (j==0) omega = 0.05;
        if (j==1) omega = 0.1;
        if (j==2) omega = 0.2;
    for (int i = 0;i<10;i++){
        double mu = 1.05+i*0.2;
        figbb(mu,omega,0,1,200000,0.1,fp);
    }}


    fclose(fp);
    return 0;
}