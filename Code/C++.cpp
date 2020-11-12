/************************************************************************************************/
/*** Topic: Wang-Buzsaki model with Runge-Kutta 4th Order Method for one neuron    Ali-Seif   ***/
/*** Version Release 17.12 rev 11256                                                          ***/
/*** Date: 11/10/2020                                                                         ***/
/*** Code implemented in CodeBlocks C++ compiler (v. 17.12),                                  ***/
/*** MSI: PX60 6QD/ DDR4                                                                      ***/
/*** Run under a Intel® Core™ i7-6700HQ CPU @ 2.60GHz × 64 based processor with 16 GB RAM     ***/
/************************************************************************************************/
#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;
//_________________________________Calculate alpha and betas_________________________________//

double alpha_n(double v) {
   return   -0.01 * (v + 34.0) / (exp(-0.1 * (v + 34.0)) - 1.0);
}

double beta_n(double v) {
   return   0.125 * exp(-(v + 44.0) / 80.0);
}

double alpha_m(double v) {
   return   -0.1 * (v + 35.0) / (exp(-0.1 * (v + 35.0)) - 1.0);
}

double beta_m(double v) {
   return   4.0 * exp(-(v + 60.0) / 18.0);
}

double alpha_h(double v) {
   return   0.07 * exp(-(v + 58.0) / 20.0);
}

double beta_h(double v) {
   return   1.0 / (exp(-0.1 * (v + 28.0)) + 1.0);
}
double F(double v) {
   int  teta=0;
   return   1.0 / (1.0+exp(-0.5 * (v - teta)));
}

//__________________________Calculate infinite activation variables__________________________//

double n_inf(double v) {
   return   alpha_n(v) / (alpha_n(v) + beta_n(v));
}

double h_inf(double v) {
   return   alpha_h(v) / (alpha_h(v) + beta_h(v));
}

double s_inf(double v) {
   int  alpha_s=12;
   double beta_s=0.1;
   return   alpha_s / (alpha_s + beta_s);
}

double m_inf(double v) {
   return   alpha_m(v) / (alpha_m(v) + beta_m(v));
}

//__________________________________Calculation of currents__________________________________//

double INa(double v,double h) {
    int   g_Na=         35;
    int   E_Na=         55;
   return    g_Na*h*pow(m_inf(v),3)*(v-E_Na);
}

double IK(double v,double n) {
    int   g_K=          9;
    int   E_K=          -90;
   return    g_K*pow(n,4)*(v-E_K);
}

double Il(double v) {
    int   E_l=          -65;
    float g_l=          0.1;
   return        g_l*(v-E_l);
}

double Isyn(double v,double s,double k) {
    int   E_s=          -75;
    float g_s=          0.1;
    float N=          1;
   return        (k/N)*g_s*s*(v-E_s);
}

//___________________________________Differential Equations__________________________________//

double dvdt(double t, double v,double n,double h,double s,double k) {
    float Iapp =        0.75;
    float C_m =         1.0;
   return   (1/C_m)*(Iapp -(INa(v,h)+IK(v,n)+Il(v)+Isyn(v,s,k)));
}

double dndt(double t,double n, double v) {
    int   phi=          5;
   return   phi*((alpha_n(v)*(1-n))-beta_n(v)*n);
}

double dhdt(double t, double h, double v) {
    int   phi=          5;
   return   phi*((alpha_h(v)*(1-h))-beta_h(v)*h);
}
double dsdt(double t, double s, double v) {
    int  alpha_s=12;
   double beta_s=0.1;

   return   ((alpha_s*F(v)*(1-s))-beta_s*s);
}


//__________________________________Runge-Kutta calculations_________________________________//

double rk4thOrder_v(double t0, double v, double dt,double n,double h,double s,double k) {
    double  k1, k2, k3, k4;
            k1=     dt*dvdt(t0, v,n,h,s,k);
            k2=     dt*dvdt((t0+dt/2), (v+k1/2),n,h,s,k);
            k3=     dt*dvdt((t0+dt/2), (v+k2/2),n,h,s,k);
            k4=     dt*dvdt((t0+dt), (v+k3),n,h,s,k);
            v=      v+double((1.0/6.0)*(k1+2*k2+2*k3+k4));
   return   v;
}

double rk4thOrder_n(double t0, double v, double dt, double n) {

    double  k1, k2, k3, k4;
            k1=     dt*dndt(t0, n,v);
            k2=     dt*dndt((t0+dt/2), (n+k1/2),v);
            k3=     dt*dndt((t0+dt/2), (n+k2/2),v);
            k4=     dt*dndt((t0+dt), (n+k3),v);
            n=      n+double((1.0/6.0)*(k1+2*k2+2*k3+k4));
   return   n;
}

double rk4thOrder_h(double t0, double v, double dt,double h) {

    double  k1, k2, k3, k4;
            k1=     dt*dhdt(t0, h,v);
            k2=     dt*dhdt((t0+dt/2), (h+k1/2),v);
            k3=     dt*dhdt((t0+dt/2), (h+k2/2),v);
            k4=     dt*dhdt((t0+dt), (h+k3),v);
            h=      h+double((1.0/6.0)*(k1+2*k2+2*k3+k4));
   return   h;
}
double rk4thOrder_s(double t0, double v, double dt,double s) {

    double  k1, k2, k3, k4;
            k1=     dt*dsdt(t0, s,v);
            k2=     dt*dsdt((t0+dt/2), (s+k1/2),v);
            k3=     dt*dsdt((t0+dt/2), (s+k2/2),v);
            k4=     dt*dsdt((t0+dt), (s+k3),v);
            s=      s+double((1.0/6.0)*(k1+2*k2+2*k3+k4));
   return   s;
}

//*******************************************************************************************//
//                                                                                           //
//________________________________The principle of the program_______________________________//
//                                                                                           //
//*******************************************************************************************//
int main() {

    double  t0=0,t_final=200, dt=0.01;                                          //Initial values
    int    number=3;                                                            //Number of neurons

    float v[number] =    {0.0};                                                 //Definition of potential
    float n[number] =    {0.0};
    float h[number] =    {0.0};
    float s[number] =    {0.0};

    for (int i=1 ; i<=number ;i++){
       v[i]=-63.0;                                                              //Initial values of potential
       n[i]=n_inf(v[i]);                                                        //Infinite initial value
       h[i]=h_inf(v[i]);                                                        //Infinite initial value
       s[i]=s_inf(v[i]);                                                        //Infinite initial values
    }

    //cout<<s[0]<<'\t'<<s[1]<<'\t'<<s[2]<<'\t'<<s[3];

    ofstream temp("temp.txt", ios::out | ios::trunc);

    for (t0 ; t0<=t_final ;t0=t0+dt){
        //temp << t0 << '\t';
        /*
        for (int i=1 ; i<=number ;i++){
            v[i]=rk4thOrder_v(t0, v[i], dt, n[i], h[i],s[i]);
            v[0]=v[3];
            n[i]=rk4thOrder_n(t0, v[i], dt ,n[i]);
            h[i]=rk4thOrder_h(t0, v[i], dt ,h[i]);
            s[i]=rk4thOrder_s(t0, v[i-1], dt ,s[i]);

            temp << v[i] << '\t';
        }
        */
        v[1]=rk4thOrder_v(t0, v[1], dt, n[1], h[1],s[1],0.4);
        n[1]=rk4thOrder_n(t0, v[1], dt ,n[1]);
        h[1]=rk4thOrder_h(t0, v[1], dt ,h[1]);
        s[1]=rk4thOrder_s(t0, v[3], dt ,s[1]);

        v[2]=rk4thOrder_v(t0, v[2], dt, n[2], h[2],s[2],0.5);
        n[2]=rk4thOrder_n(t0, v[2], dt ,n[2]);
        h[2]=rk4thOrder_h(t0, v[2], dt ,h[2]);
        s[2]=rk4thOrder_s(t0, v[1], dt ,s[2]);

        v[3]=rk4thOrder_v(t0, v[3], dt, n[3], h[3],s[3],0.6);
        n[3]=rk4thOrder_n(t0, v[3], dt ,n[3]);
        h[3]=rk4thOrder_h(t0, v[3], dt ,h[3]);
        s[3]=rk4thOrder_s(t0, v[2], dt ,s[3]);

        temp<<t0 <<'\t'<<v[1]<<'\t'<<v[2]<<'\t'<<v[3]<<endl;
        //temp  << endl;

    }

    temp.close();
    cout << "\nFinish" << endl;

    return 0;
}
