#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include <iomanip>
#include <omp.h>  
#include <algorithm>
#include "TDMA_solver.h"
#pragma once
using namespace std ; 
//會用到的參數
constexpr int Nx = 160 ;
constexpr int Ny = 160 ;
constexpr double dx = 1.0 / double(Nx);
constexpr double dy = 1.0 / double(Ny);
constexpr double rho = 1.0;
const int n = Nx * Ny;
double Coe_P[Nx][Ny] , Coe_W[Nx][Ny] , Coe_E[Nx][Ny] , Coe_S[Nx][Ny] , Coe_N[Nx][Ny] , Source[Nx][Ny] ; 
//infact we will define the Boundary Condition and the size od the vector 
extern double u_star[Nx+1][Ny+2] ; //x-velocity predicted value  J + 1 
extern double v_star[Nx+2][Ny+1] ; //y-velocity predcted valure  I + 1
extern double a_P[Nx+1][Ny] ; //xpart velocity pmain coefficient  
extern double b_P[Nx][Ny+1] ; //y-part velocity main coefficient
extern double p_prime[Nx][Ny] ; //pressure correction vlaue
double p_prime_old[Nx][Ny] , p_prime_old2[Nx][Ny] ; //pressure correcrion old value


//initialize
void initialze() {
    #pragma omp parallel
    {
        #pragma omp for collapse(2) schedule(static)
        for (int I = 1; I <= Nx; ++I) {
            for (int J = 1; J <= Ny; ++J) {
                Coe_P[I-1][J-1] = 0.0;
                Coe_W[I-1][J-1] = 0.0;
                Coe_E[I-1][J-1] = 0.0;
                Coe_S[I-1][J-1] = 0.0;
                Coe_N[I-1][J-1] = 0.0;
                Source[I-1][J-1] = 0.0;
                //First step : initialize the coefficient for solving the presure poisson equation
            }
        }
    }
}

//setting coefficient 
void coefficient(){
    // 先複製當前解到x_old用於計算誤差
    #pragma omp parallel default(shared)
    {
    #pragma omp  for collapse(2) schedule(static)
    for(int I = 2; I <= Nx-1; I++) {
        for(int J = 2; J <= Ny-1; J++) {
            Coe_P[I-1][J-1] = (((dy*dy)/(rho*(a_P[I][J-1])))+((dy*dy)/(rho*(a_P[I-1][J-1])))+((dx*dx)/(rho*(b_P[I-1][J])))+((dx*dx)/(rho*(b_P[I-1][J-1])))) ;
            Coe_E[I-1][J-1] = -((dy*dy)/(rho*(a_P[I][J-1]))) ;
            Coe_W[I-1][J-1] = -((dy*dy)/(rho*(a_P[I-1][J-1]))) ;
            Coe_N[I-1][J-1] = -((dx*dx)/(rho*(b_P[I-1][J]))) ;
            Coe_S[I-1][J-1] = -((dx*dx)/(rho*(b_P[I-1][J-1])));
            Source[I-1][J-1] = dy*u_star[I-1][J] - dy*u_star[I][J] + dx*v_star[I][J-1] - dx*v_star[I][J] ;
            
        }
    }
    //right boundary nodes
    #pragma omp  for schedule(static)
    for(int J = 2; J <= Ny-1; J++) {
        Coe_P[Nx-1][J-1] =(((dy*dy)/(rho*(a_P[Nx-1][J-1])))+((dx*dx)/(rho*(b_P[Nx-1][J])))+((dx*dx)/(rho*(b_P[Nx-1][J-1]))));
        Coe_W[Nx-1][J-1] = -((dy*dy)/(rho*(a_P[Nx-1][J-1]))) ; 
        Coe_N[Nx-1][J-1] = -((dx*dx)/(rho*(b_P[Nx-1][J]))) ;
        Coe_S[Nx-1][J-1] = -((dx*dx)/(rho*(b_P[Nx-1][J-1]))) ;
        Source[Nx-1][J-1] = dy*u_star[Nx-1][J] - dy*u_star[Nx][J] + dx*v_star[Nx][J-1] - dx*v_star[Nx][J] ;
    }
    //left boundary nodes
    #pragma omp  for schedule(static)
    for(int J = 2; J <= Ny-1; J++) {
        Coe_P[1-1][J-1] =(((dy*dy)/(rho*(a_P[1-1+1][J-1])))+((dx*dx)/(rho*(b_P[1-1][J-1+1])))+((dx*dx)/(rho*(b_P[1-1][J-1]))));
        Coe_E[1-1][J-1] = -((dy*dy)/(rho*(a_P[1][J-1]))) ; 
        Coe_N[1-1][J-1] = -((dx*dx)/(rho*(b_P[1-1][J]))) ;
        Coe_S[1-1][J-1] = -((dx*dx)/(rho*(b_P[1-1][J-1]))) ;
        Source[1-1][J-1] = dy*u_star[1-1][J] - dy*u_star[1][J] + dx*v_star[1][J-1] - dx*v_star[1][J]  ;
    }
    //top boundary nodes
    #pragma omp  for schedule(static)
    for(int I = 2; I <= Nx-1; I++) {
        Coe_P[I-1][Ny-1] =(((dy*dy)/(rho*(a_P[I][Ny-1])))+((dy*dy)/(rho*(a_P[I-1][Ny-1])))+((dx*dx)/(rho*(b_P[I-1][Ny-1])))) ;
        Coe_E[I-1][Ny-1] = -((dy*dy)/(rho*(a_P[I][Ny-1]))) ; 
        Coe_W[I-1][Ny-1] = -((dy*dy)/(rho*(a_P[I-1][Ny-1]))) ;
        Coe_S[I-1][Ny-1] = -((dx*dx)/(rho*(b_P[I-1][Ny-1]))) ;
        Source[I-1][Ny-1] = dy*u_star[I-1][Ny] - dy*u_star[I][Ny] + dx*v_star[I][Ny-1] - dx*v_star[I][Ny] ;
    }
    //bottom boundary nodes
    #pragma omp  for schedule(static)
    for(int I = 2; I <= Nx-1 ; I++) {
        Coe_P[I-1][1-1] = (((dy*dy)/(rho*(a_P[I][1-1])))+((dy*dy)/(rho*(a_P[I-1][1-1]))) +((dx*dx)/(rho*(b_P[I-1][1])))) ;
        Coe_E[I-1][1-1] =-((dy*dy)/(rho*(a_P[I][1-1])));
        Coe_W[I-1][1-1] =-((dy*dy)/(rho*(a_P[I-1][1-1]))) ;
        Coe_N[I-1][1-1] =-((dx*dx)/(rho*(b_P[I-1][1])));
        Source[I-1][1-1] = dy*u_star[I-1][1] - dy*u_star[I][1] + dx*v_star[I][1-1] - dx*v_star[I][1] ;
    }
    {
        Coe_P[1-1][1-1] = 1.0;
        Coe_E[1-1][1-1] = 0.0 ;
        Coe_N[1-1][1-1] = 0.0;
        Source[1-1][1-1] = 0.0 ;
    }
    {
        Coe_P[Nx-1][1-1] =(((dy*dy)/(rho*(a_P[Nx-1][1-1])))+((dx*dx)/(rho*(b_P[Nx-1][1])))) ;
        Coe_W[Nx-1][1-1] = -((dy*dy)/(rho*(a_P[Nx-1][1-1])));
        Coe_N[Nx-1][1-1] = -((dx*dx)/(rho*(b_P[Nx-1][1])));
        Source[Nx-1][1-1] = dy*u_star[Nx-1][1] - dy*u_star[Nx][1] + dx*v_star[Nx][1-1] - dx*v_star[Nx][1] ;
    }
    {
        Coe_P[1-1][Ny-1] = (((dy*dy)/(rho*(a_P[1][Ny-1])))+((dx*dx)/(rho*(b_P[1-1][Ny-1])))) ;
        Coe_E[1-1][Ny-1] =  -((dy*dy)/(rho*(a_P[1][Ny-1]))) ;
        Coe_S[1-1][Ny-1] =  -((dx*dx)/(rho*(b_P[1-1][Ny-1])));
        Source[1-1][Ny-1] =  dy*u_star[1-1][Ny] - dy*u_star[1][Ny] + dx*v_star[1][Ny-1] - dx*v_star[1][Ny] ;
    }
    {
        Coe_P[Nx-1][Ny-1] =(((dy*dy)/(rho*(a_P[Nx-1][Ny-1])))+((dx*dx)/(rho*(b_P[Nx-1][Ny-1])))); //coefficient for main elements
        Coe_W[Nx-1][Ny-1] = -((dy*dy)/(rho*(a_P[Nx-1][Ny-1]))) ;
        Coe_S[Nx-1][Ny-1] = -((dx*dx)/(rho*(b_P[Nx-1][Ny-1])));
        Source[Nx-1][Ny-1] = dy*u_star[Nx-1][Ny] - dy*u_star[Nx][Ny] + dx*v_star[Nx][Ny-1] - dx*v_star[Nx][Ny] ;
    }
    }
}
void ADII_1(){
    //Non-need to itteration and this is the directly solving process.
    #pragma omp parallel default(shared)
    {
        //first : given value to the old value
        #pragma omp for collapse(2) schedule(static) //average cutting workking for threads 
        for(int I = 1 ; I <= Nx ; I++){
            for(int J = 1 ; J <= Ny ; J++){
                p_prime_old[I-1][J-1] = p_prime[I-1][J-1]; 
            }
        }
        //find the coefficeints of the computing process D , beta , alpha , x , S , boundary condition 
        #pragma omp for schedule(static) //average cutting workking for threads 
        for(int J = 1 ; J <= Ny ; J++){//define the index from 0~Nx-2
            vector <double> beta (Nx) ; vector <double> D (Nx) ; vector <double> alpha (Nx) ;
            vector <double> x(Nx) ; vector <double> S (Nx) ; //for neumann condition , det the boundary condition is zero through the confficients setting to get neumann Condition 
            for(int I = 1 ; I <= Nx ; I++){
            beta[I-1] = ((I == 1) ? 0.0 : Coe_W[I-1][J-1]) ;
            D[I-1] =  ((I == 1 && J == 1) ? 1.0 :  Coe_P[I-1][J-1]) ; 
            alpha[I-1] =   (((I == 1 && J == 1) || I == Nx ) ? 0.0 :  Coe_E[I-1][J-1]) ; 
            x[I-1] = p_prime[I-1][J-1] ;
            if(I == 1 && J == 1){ Source[I-1][J-1] = 0.0 ; }
            else if(J == 1 || J == Ny){
                S[I-1] = Source[I-1][J-1] + ((J == 1) ?  - Coe_N[I-1][J-1]*p_prime_old[I-1][J]: - Coe_S[I-1][J-1]*p_prime_old[I-1][J-2]); 
            }else{S[I-1] = Source[I-1][J-1] - Coe_S[I-1][J-1]*p_prime_old[I-1][J-2] - Coe_N[I-1][J-1]*p_prime_old[I-1][J] ; }
            }
            TDMA_FVM(D , beta , alpha , x , S , Nx) ; //x[NX] = pressure correction value 
            for(int I = 1 ; I <= Nx ; I++){
                p_prime[I-1][J-1] = x[I-1] ;
            }
        }
        //second : givee the value to the process value 
        #pragma omp for collapse(2) schedule(static) //average cutting workking for threads 
        for(int I = 1 ; I <= Nx ; I++){
            for(int J = 1 ; J <= Ny ; J++){
                p_prime_old2[I-1][J-1] = p_prime[I-1][J-1]; 
            }
        }
        #pragma omp for schedule(static) //average cutting workking for threads 
        for(int I = 1 ; I <=Nx ; I++){//define the index from 1 ~ Ny
            vector <double> beta(Ny) ; vector <double> D(Ny) ; vector <double> alpha (Ny) ;
            vector <double> x (Ny) ; vector <double> S (Ny) ; //For Neumann Condition , We first set the B.C. is zero and through coefficient setting to schieve the Neumann Condition  
            for(int J = 1 ; J <= Ny ; J++){
            beta[J-1] = ((J == 1) ? 0.0 : Coe_S[I-1][J-1]) ; 
            D[J-1] = ((J == 1 && I == 1) ? 1.0 : Coe_P[I-1][J-1]) ;  
            alpha[J-1]= (((J == 1 && I == 1) || J == Ny ) ? 0.0 : Coe_N[I-1][J-1] ) ; 
            x[J-1] = p_prime[I-1][J-1] ; 
            if(I == 1 && J == 1){ Source[I-1][J-1] = 0.0 ; }
            else if(I == 1 || I == Nx){S[J-1] = Source[I-1][J-1]+( (I == Nx) ? - Coe_W[I-1][J-1]*p_prime_old2[I-2][J-1] : - Coe_E[I-1][J-1]*p_prime_old2[I][J-1] );}
            else{S[J-1] = Source[I-1][J-1] - Coe_W[I-1][J-1]*p_prime_old2[I-2][J-1] - Coe_E[I-1][J-1]*p_prime_old2[I][J-1] ; }
            }
            TDMA_FVM(D , beta , alpha , x , S , Ny) ; 
            for(int J = 1 ; J <= Ny ; J++){
                p_prime[I-1][J-1] =  x[J-1]  ;
            }
        }
    } 
}


void step3(){
    initialze() ; 
    coefficient() ; 
    ADII_1() ;
}