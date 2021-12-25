#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define M (2000)
#define Mu (2)
#define Lambda (0.5)
#define Epsilon (0)
#define Omega (0)
#define L1 (0.6)
#define L2 (0.1)


double ht(double t){
    return Epsilon*sin(Omega*t);
}

double* f(double t,double A[],int nn){
    double* B=(double*)malloc(sizeof(double)*nn);
    double q1=A[0];
    double q2=A[1];
    double p1=A[2];
    double p2=A[3];
    B[0] = L1*(2*p1-1)*q1 - L2*q1 - 1*(2*p1-1)*q1*q1/M - 1*(2*p2-1)*q2*q1/Mu/M;
    B[1] = L1*(2*p2-1)*q2 - L2*q2 - 1*(2*p2-1)*q2*q2/M - 1*(2*p1-1)*q1*q2/M/Mu;
    B[2] = L1*(1-p1)*p1 + L2*(p1-1) + 2*(p1-1)*q1*p1/M + 1*(2*p1*p2-p1-p2)*q2/Mu/M;
    B[3] = L1*(1-p2)*p2 + L2*(p2-1) + 2*(p2-1)*q2*p2/M + 1*(2*p1*p2-p1-p2)*q1/Mu/M;
    return B;
}



double* make(double a,double b,double c,double d,int n){
    double* B=(double*)malloc(sizeof(double)*n);
    B[0] = a;
    B[1] = b;
    B[2] = c;
    B[3] = d;
    return B;
}

void improvedeuler(double* tem,double t,double q1,double q2,double p1 ,double p2,double tm,double h,FILE *fp){
    int ii = 0;
    while (t<tm){
        // printf("%f %f %f %f %f\n",t,q1,q2,p1,p2);

        fwrite(&t, sizeof(double) , 1, fp );
        fwrite(&q1, sizeof(double) , 1, fp );
        fwrite(&q2, sizeof(double) , 1, fp );
        fwrite(&p1, sizeof(double) , 1, fp );
        fwrite(&p2, sizeof(double) , 1, fp );

        tem[0*20001+20000-ii] = t;
        tem[1*20001+20000-ii] = q1;
        tem[2*20001+20000-ii] = q2;
        tem[3*20001+20000-ii] = p1;
        tem[4*20001+20000-ii] = p2;

        double q1_ = q1 + h*f(t,make(q1,q2,p1,p2,4),4)[0];
        double q2_ = q2 + h*f(t,make(q1,q2,p1,p2,4),4)[1];
        double p1_ = p1 + h*f(t,make(q1,q2,p1,p2,4),4)[2];
        double p2_ = p2 + h*f(t,make(q1,q2,p1,p2,4),4)[3];

        double q1__ = q1 + h/2*(f(t,make(q1,q2,p1,p2,4),4)[0]+f(t,make(q1_,q2_,p1_,p2_,4),4)[0]);
        double q2__ = q2 + h/2*(f(t,make(q1,q2,p1,p2,4),4)[1]+f(t,make(q1_,q2_,p1_,p2_,4),4)[1]);
        double p1__ = p1 + h/2*(f(t,make(q1,q2,p1,p2,4),4)[2]+f(t,make(q1_,q2_,p1_,p2_,4),4)[2]);
        double p2__ = p2 + h/2*(f(t,make(q1,q2,p1,p2,4),4)[3]+f(t,make(q1_,q2_,p1_,p2_,4),4)[3]);
        
        q1 = q1__;
        q2 = q2__;
        p1 = p1__;
        p2 = p2__;
        t += h;
        ii += 1;
    }
}


void improvedeuler2(double* tem,double t,double q1,double q2,double p1 ,double p2,double tm,double h,FILE *fp){
    int ii = 0;
    while (t<tm){
        // printf("--%f %f %f %f %f\n",t,q1,q2,p1,p2);
        fwrite(&t, sizeof(double) , 1, fp );
        fwrite(&q1, sizeof(double) , 1, fp );
        fwrite(&q2, sizeof(double) , 1, fp );
        fwrite(&p1, sizeof(double) , 1, fp );
        fwrite(&p2, sizeof(double) , 1, fp );
        
        tem[0*20001+ii] = t;
        tem[1*20001+ii] = q1;
        tem[2*20001+ii] = q2;
        p1 = tem[3*20001+ii];
        p2 = tem[4*20001+ii];

        double q1_ = q1 + h*f(t,make(q1,q2,p1,p2,4),4)[0];
        double q2_ = q2 + h*f(t,make(q1,q2,p1,p2,4),4)[1];
        double p1_ = p1 + h*f(t,make(q1,q2,p1,p2,4),4)[2];
        double p2_ = p2 + h*f(t,make(q1,q2,p1,p2,4),4)[3];

        double q1__ = q1 + h/2*(f(t,make(q1,q2,p1,p2,4),4)[0]+f(t,make(q1_,q2_,p1_,p2_,4),4)[0]);
        double q2__ = q2 + h/2*(f(t,make(q1,q2,p1,p2,4),4)[1]+f(t,make(q1_,q2_,p1_,p2_,4),4)[1]);
        double p1__ = p1 + h/2*(f(t,make(q1,q2,p1,p2,4),4)[2]+f(t,make(q1_,q2_,p1_,p2_,4),4)[2]);
        double p2__ = p2 + h/2*(f(t,make(q1,q2,p1,p2,4),4)[3]+f(t,make(q1_,q2_,p1_,p2_,4),4)[3]);
        
        q1 = q1__;
        q2 = q2__;
        p1 = p1__;
        p2 = p2__;
        t += h;
        ii += 1;
    }
}


void improvedeuler3(double* tem,double t,double q1,double q2,double p1 ,double p2,double tm,double h,FILE *fp){
    int ii = 0;
    while (t<tm){
        // printf("--%f %f %f %f %f\n",t,q1,q2,p1,p2);
        tem[3*20001+20000-ii] = p1;
        tem[4*20001+20000-ii] = p2;
        q1 = tem[1*20001+20000-ii];
        q2 = tem[2*20001+20000-ii];
        fwrite(&t, sizeof(double) , 1, fp );
        fwrite(&q1, sizeof(double) , 1, fp );
        fwrite(&q2, sizeof(double) , 1, fp );
        fwrite(&p1, sizeof(double) , 1, fp );
        fwrite(&p2, sizeof(double) , 1, fp );

        double q1_ = q1 - h*f(t,make(q1,q2,p1,p2,4),4)[0];
        double q2_ = q2 - h*f(t,make(q1,q2,p1,p2,4),4)[1];
        double p1_ = p1 - h*f(t,make(q1,q2,p1,p2,4),4)[2];
        double p2_ = p2 - h*f(t,make(q1,q2,p1,p2,4),4)[3];

        double q1__ = q1 - h/2*(f(t,make(q1,q2,p1,p2,4),4)[0]+f(t,make(q1_,q2_,p1_,p2_,4),4)[0]);
        double q2__ = q2 - h/2*(f(t,make(q1,q2,p1,p2,4),4)[1]+f(t,make(q1_,q2_,p1_,p2_,4),4)[1]);
        double p1__ = p1 - h/2*(f(t,make(q1,q2,p1,p2,4),4)[2]+f(t,make(q1_,q2_,p1_,p2_,4),4)[2]);
        double p2__ = p2 - h/2*(f(t,make(q1,q2,p1,p2,4),4)[3]+f(t,make(q1_,q2_,p1_,p2_,4),4)[3]);
        
        q1 = q1__;
        q2 = q2__;
        p1 = p1__;
        p2 = p2__;
        t += h;
        ii += 1;
    }
}


int main(){
    FILE *fp;
    double tm = 200;
// 1
    // 正向求Q
    fp = fopen( "fig3_2.bin" , "w" );
    double tem[5*20001] = {0};
    improvedeuler(tem,0,1.e-3,1000,1,1,tm,0.01,fp);
    fclose(fp);


    for (int i = 0; i<10;i++){
        // 反向求P
        fp = fopen( "fig3_1.bin" , "w" );
        improvedeuler3(tem,0,0,1000,7/12.0,1,tm,0.01,fp);
        fclose(fp);

        // 正向求Q

        fp = fopen( "fig3_2.bin" , "w" );
        improvedeuler2(tem,0,2000./3,2000./3,1,1,tm,0.01,fp);
        fclose(fp);
    }


// 2
    // 正向求Q
    fp = fopen( "fig3_4.bin" , "w" );
    improvedeuler(tem,0,1000,1.e-3,1,1,tm,0.01,fp);
    fclose(fp);


    for (int i = 0; i<10;i++){
        // 反向求P
        fp = fopen( "fig3_3.bin" , "w" );
        improvedeuler3(tem,0,1000,0,1,7/12.0,tm,0.01,fp);
        fclose(fp);

        // 正向求Q

        fp = fopen( "fig3_4.bin" , "w" );
        improvedeuler2(tem,0,2000./3,2000./3,1,1,tm,0.01,fp);
        fclose(fp);
    }

    return 0;
}