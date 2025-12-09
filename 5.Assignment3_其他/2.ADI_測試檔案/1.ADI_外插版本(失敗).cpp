#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <vector> 
#include <omp.h> 
#include "0.P_correction_eq_TDMA.h" 
#include "TDMA_solver.h"//adding the TDMA to solve the predicted velocity 
//this code is taiking about "deferred correction and adding the domain boundary for velocity predicted value".
using namespace std ;
int G ; 
double min_step7 ;
const int max_G = 3000000 ; //最大步數設定
double omega_u = 0.0002 ; 
double omega_v = 0.0002 ; 
double omega_P = 0.0005 ; 
double u_star[Nx+1][Ny+2] , u_star_old[Nx+1][Ny+2] , u_star_old2[Nx+1][Ny+2] , u_prime[Nx+1][Ny] , u[Nx+1][Ny] ; //x-elocity predicted value and correctin value and correct value 
double v_star[Nx+2][Ny+1] , v_star_old[Nx+2][Ny+1] , v_star_old2[Nx+2][Ny+1] , v_prime[Nx][Ny+1] , v[Nx][Ny+1] ; //y-elocity predicted value and correctin value and correct value 
double p_star[Nx][Ny] , p_prime[Nx][Ny] , p[Nx][Ny] ; //pressure predicted value and pressure correction value and  correct pressure
const double Re =100.0 ; 
//Boundary Condition 
const double u_left = 0.0 ;const double u_right = 0.0 ; const double u_bottom = 0.0 ; const double u_top = 1.0 ; 
const double v_left = 0.0 ;const double v_right = 0.0 ; const double v_bottom = 0.0 ; const double v_top = 0.0 ; 
//質量流率收斂判據
const float CE_tolerance = 1e-10 ; 
bool CE_physical ; 
//coefficient of the discrete diffuion term 
const double D_e = (1.0/Re)*(dy/dx) ; 
const double D_w = (1.0/Re)*(dy/dx) ; 
const double D_n = (1.0/Re)*(dx/dy) ; 
const double D_s = (1.0/Re)*(dx/dy) ;
//r_f is the variable of the flux limiter
//convection term about u velocity
//domain of the r_f_u is the i : 1 ~ Nx+1 , J : 1 ~ Ny 
double r_w_u[Nx+1][Ny] , r_e_u[Nx+1][Ny] , r_s_u[Nx+1][Ny] , r_n_u[Nx+1][Ny] ; 
double alpha_w_u[Nx+1][Ny] , alpha_e_u[Nx+1][Ny] , alpha_s_u[Nx+1][Ny] , alpha_n_u[Nx+1][Ny] ; 
double F_w_u[Nx+1][Ny] , F_e_u[Nx+1][Ny] , F_s_u[Nx+1][Ny] , F_n_u[Nx+1][Ny] ; 
//convection term about v velocity
//domain of the r_f_v is the I = 1 ~ Nx , j  = 1 ~ Ny+1
double r_w_v[Nx][Ny+1] , r_e_v[Nx][Ny+1] , r_s_v[Nx][Ny+1] , r_n_v[Nx][Ny+1] ; 
double alpha_w_v[Nx][Ny+1] , alpha_e_v[Nx][Ny+1], alpha_s_v[Nx][Ny+1] , alpha_n_v[Nx][Ny+1] ; 
double F_w_v[Nx][Ny+1] , F_e_v[Nx][Ny+1] , F_s_v[Nx][Ny+1] , F_n_v[Nx][Ny+1] ; 
//coefficient of the discrete x-momentum equation
double a_P[Nx+1][Ny] , a_W[Nx+1][Ny] , a_E[Nx+1][Ny] , a_S[Nx+1][Ny] , a_N[Nx+1][Ny] , S1[Nx+1][Ny] ; 
//coefficient of the discrete y-momentum equaiton 
double b_P[Nx][Ny+1] , b_W[Nx][Ny+1] , b_E[Nx][Ny+1] , b_S[Nx][Ny+1] , b_N[Nx][Ny+1] , S2[Nx][Ny+1] ; 
//TVD-MUSCL funtion : upstream and downstream coefficient funtion : 
double D(double x) {
    /*return 0.5 ; //Central Diffeence Scheme */
    return 3.0/8.0 ; //2.QUICK Scheme 
    /*if (x <= 0) return 0.0;
    else if (x > 0 && x < (1.0/3.0)) return 0.0;
    else if (x >= (1.0/3.0) && x < 3.0) return 0.25;
    else return 1.0; //3.TVD-MUSCL*/
}

double U(double x) {
    /*return 0.5 ; //Central Diffeence Scheme */
    return 6.0/8.0 ; //2.QUICK Scheme 
    /*if (x <= 0) return 1.0;
    else if (x > 0 && x < (1.0/3.0)) return 1.0;
    else if (x >= (1.0/3.0) && x < 3.0) return 0.0;
    else return 0.0;//3.TVD-MUSCL*/
}

double UU(double x) {
    /*return 0.0 ; //Central Diffeence Scheme */
    return -1.0/8.0 ; //2.QUICK Scheme 
     /*if (x <= 0) return 0.0;
    else if (x > 0 && x < (1.0/3.0)) return -1.0;
    else if (x >= (1.0/3.0) && x < 3.0) return -0.25;
    else return 0.0;//3.TVD-MUSCL*/
}
//Simple algorithm : 
//step1 : initialize what ? 
//: Given  the predicted value. 
//So the variables about the predicted form of the discrete x-momentum equation and the predicted form of the discrete y-momentum equaition need to be initialized .
// u_star,u_star_old,v_star,v_star_old,p_star,
//F,alpha,r_f......for w,e,s,n side and u and v vice computational node need to be initialized ...
//We have initialize 29 Varibles  
//step2 : compute velocity predixted value 
//using u_star_old[I-1][j] , u_star[I-1][j] , v_star_old[i][J-1] , v_star_old[i][J-1] //watch out : the x-velocity need add the bottom layer and top layer for expolation and the y-velocity need add the left and right layer for expolation 
//step3 : compute the pressure correction equaition using inner iteration 
//step4 : update  three things : correct pressure , velocity corretion and corrct velocity 
//step5 : converge creteria : mass flow rate or the the value comparred with U. Ghia. Benchmark 
//step6 : update the predixted pressure and update predicted velocity 
//and the algorithm loop is the SIMPLE algorithm loop.....
void step1(){
    #pragma omp parallel default(shared)
    {
        //1.u_star,u_star_old,v_star,v_star_old,p_star initialize 
        //I = [1,NX+1] , j = [0.Ny+1] 
        #pragma omp for collapse (2) schedule(static)//(1)
        for(int i = 1 ; i <= Nx+1 ; i++){
            for(int J = 1 ; J <= Ny ; J++){
                //initialize the u_star 
                if(i == 1 || i == Nx+1){
                    u_star[i-1][J] = ((i == 1) ? u_left : u_right) ; 
                }else {
                    u_star[i-1][J]=0.0 ; //second we decide the inner node
                }
            }
        }
        //given value for the extrapolation
        #pragma omp for schedule(static)
        for(int i = 1 ; i <= Nx+1 ; i++){
            u_star[i-1][0] = 2*u_bottom - u_star[i-1][1] ; 
            u_star[i-1][Ny+1] = 2*u_top - u_star[i-1][Ny] ;
        }
        #pragma omp for collapse (2) schedule(static)//(1)
        for(int i = 1 ; i <= Nx+1 ; i++){
            for(int J = 0 ; J <= Ny+1 ; J++){
                u_star_old[i-1][J] = u_star[i-1][J] ; 
                u_star_old2[i-1][J] = u_star[i-1][J] ; 
            }
        }
        #pragma omp for collapse (2) schedule(static)//(2)
        //i = [0,Nx+1] , j = [1,Ny+1]
        for(int I = 1 ; I <= Nx ; I++){
            for(int j = 1 ; j <= Ny+1 ; j++){
                //initialize the v_ star 
                if(j == 1 || j == Ny+1){
                    v_star[I][j-1] = ((j == 1) ? v_bottom : v_top) ;
                }else {
                    v_star[I][j-1] = 0.0 ; //second we decide the inner node
                }
            }
        }
        //given value for the extrapolation
        #pragma omp for schedule(static)
        for(int j = 1 ; j <= Ny+1 ; j++){
            v_star[0][j-1] = 2*v_left - v_star[1][j-1] ; 
            v_star[Nx+1][j-1] = 2*v_right - v_star[Nx][j-1] ;
        }
        #pragma omp for collapse (2) schedule(static)//(2)
        //i = [0,Nx+1] , j = [1,Ny+1]
        for(int I = 0 ; I <= Nx+1 ; I++){
            for(int j = 1 ; j <= Ny+1 ; j++){
                v_star_old[I][j-1] = v_star[I][j-1] ;
                v_star_old2[I][j-1] = v_star[I][j-1] ;
            }
        }
        #pragma omp for collapse (2) schedule(static)//(3)
        //I = [1,Nx] , J = [1,Ny]//initilaize the pressure predicted value 
        for(int I = 1 ; I <= Nx ; I++){
            for(int J = 1 ; J <= Ny ; J++){
                p_star[I-1][J-1] = 0.0 ; 
            }
        }
        #pragma omp for collapse (2) schedule(static)//(4)
        for(int i = 1 ; i <= Nx+1 ; i++){
            for(int J = 1 ; J <= Ny ; J++){
                F_w_u[i-1][J-1] = ((i == 1 || i == Nx+1) ? 0.0 : dy*0.5*(u_star_old[i][J]+u_star_old[i-1][J])) ;
                alpha_w_u[i-1][J-1] = ((F_w_u[i-1][J-1] > 0) ? 1.0 : 0.0 ) ;
                r_w_u[i-1][J-1] = ((i == 1 || i == Nx+1 || i == 2) ? 0.0 : 0.0) ; //the finial term in step 1 give 0.0 
                F_e_u[i-1][J-1] = ((i == 1 || i == Nx+1) ? 0.0 : dy*0.5*(u_star_old[i-2][J]+u_star_old[i-1][J])) ;
                alpha_e_u[i-1][J-1] = ((F_e_u[i-1][J-1] > 0) ? 1.0 : 0.0 ) ;
                r_e_u[i-1][J-1] = ((i == 1 || i == Nx+1 || i == Nx) ? 0.0 : 0.0) ; //the finial term in step 1 give 0.0 
                F_n_u[i-1][J-1] = ((i == 1 || i == Nx+1) ? 0.0 : dx*0.5*(v_star_old[i-1][J]+v_star_old[i][J])) ;
                alpha_n_u[i-1][J-1] = ((F_n_u[i-1][J-1] > 0) ? 1.0 : 0.0 ) ;
                r_n_u[i-1][J-1] = ((i == 1 || i == Nx+1 || J == 1) ? 0.0 : 0.0 ) ;  
                F_s_u[i-1][J-1] = ((i == 1 || i == Nx+1) ? 0.0 : dx*0.5*(v_star_old[i-1][J-1]+v_star_old[i][J-1])) ;
                alpha_s_u[i-1][J-1] = ((F_s_u[i-1][J-1] > 0) ? 1.0 : 0.0 ) ;
                r_s_u[i-1][J-1] = ((i == 1 || i == Nx+1 || J == Ny) ? 0.0 : 0.0 ) ; 
            }
        }
        #pragma omp for collapse (2) schedule(static)//(5)
        for(int I = 1 ; I <= Nx ; I++){
            for(int j = 1 ; j <= Ny+1 ; j++){
                F_w_v[I-1][j-1] = ((j == 1 || j == Ny+1) ? 0.0 : dy*0.5*(u_star_old[I-1][j]+u_star_old[I-1][j-1])) ;
                alpha_w_v[I-1][j-1] = ((F_w_v[I-1][j-1] > 0) ? 1.0 : 0.0 ) ; 
                r_w_v[I-1][j-1] = ((j == 1 || j == Ny+1 || I == 1) ? 0.0 : 0.0) ; 
                F_e_v[I-1][j-1] = ((j == 1 || j == Ny+1) ? 0.0 : dy*0.5*(u_star_old[I][j]+u_star_old[I][j-1])) ; 
                alpha_e_v[I-1][j-1] = ((F_e_v[I-1][j-1] > 0) ? 1.0 : 0.0 ) ; 
                r_e_v[I-1][j-1] = ((j == 1 || j == Ny+1 || I == Nx) ? 0.0 : 0.0) ;
                F_s_v[I-1][j-1] = ((j == 1 || j == Ny+1) ? 0.0 : dx*0.5*(v_star_old[I][j-1]+v_star_old[I][j-2])) ; 
                alpha_s_v[I-1][j-1] = ((F_s_v[I-1][j-1] > 0) ? 1.0 : 0.0 ) ;
                r_s_v[I-1][j-1] = ((j ==1 || j == Ny+1 || j == 2) ? 0.0 : 0.0) ;  
                F_n_v[I-1][j-1] = ((j == 1 || j == Ny+1) ? 0.0 : dx*0.5*(v_star_old[I][j-1]+v_star_old[I][j])) ; 
                alpha_n_v[I-1][j-1] = ((F_n_v[I-1][j-1] > 0) ? 1.0 : 0.0 ) ; 
                r_n_v[I-1][j-1] = ((j ==1 || j == Ny+1 || j == Ny) ? 0.0 : 0.0) ; 
            }
        }

    

    }
}

void step2(){//solving x part .... 
    #pragma omp parallel for collapse(2) //Openmp 執行係數部分加速運算
    //compute prrdictednform of the x-momentum equation 
    for(int i = 1 ; i <= Nx+1 ; i++){
        for(int J = 1 ; J <=  Ny ; J++){
            //you can see that we set the top and bottom layer to set the unexisted point..
            a_P[i-1][J-1] = 0.0 ; a_W[i-1][J-1] = 0.0 ; a_E[i-1][J-1] = 0.0 ; a_S[i-1][J-1] = 0.0 ; a_N[i-1][J-1] = 0.0 ; S1[i-1][J-1] = 0.0 ; //set the coefficient of the discrete z-momentum eqaution = 0.0 ..
            if(i == 1){
                u_star[i-1][J] = u_left ; 
            }else if(i == Nx+1){
                u_star[i-1][J] = u_right ; 
            }else{
                //special case is for north side and south side J = 1 about south side and J = Ny about north side ...
                //L.H.S
                //(1)west side ..despite of the i == 2 , for L.H.S shoud  use upwind scheme for deferredd correction 
                //so...
                a_W[i-1][J-1] += -D_w + (-F_w_u[i-1][J-1])*alpha_w_u[i-1][J-1];
                a_P[i-1][J-1] += +D_w + (-F_w_u[i-1][J-1])*(1-alpha_w_u[i-1][J-1]);
                //(2)east side
                a_E[i-1][J-1] += -D_e + (F_e_u[i-1][J-1])*(1-alpha_e_u[i-1][J-1]);
                a_P[i-1][J-1] += +D_e + (F_e_u[i-1][J-1])*(alpha_e_u[i-1][J-1]);
                //(3)south side
                a_S[i-1][J-1] += ( (J == 1) ? 0 : -D_s + (-F_s_u[i-1][J-1])*alpha_s_u[i-1][J-1]) ; 
                a_P[i-1][J-1] += ((J == 1) ? 2*D_s : +D_s + (-F_s_u[i-1][J-1])*(1-alpha_s_u[i-1][J-1])) ;
                //(4)north side 
                a_N[i-1][J-1] += ((J == Ny) ? 0 : -D_n + (F_n_u[i-1][J-1])*(1-alpha_n_u[i-1][J-1])) ; 
                a_P[i-1][J-1] += ((J == Ny) ? 2*D_n : +D_n + (F_n_u[i-1][J-1])*alpha_n_u[i-1][J-1]) ;
                //R.H.S(source term)----deferred correction
                //(1)west side
                double A_P1, A_WW1, A_W1, A_E1 ; //tracient variables//because the deferred correction , it will note that upwind stream coefficient nee to be minused 1                
                A_P1 = ((1-alpha_w_u[i-1][J-1])*(U(r_w_u[i-1][J-1])-1)+alpha_w_u[i-1][J-1]*(D(r_w_u[i-1][J-1])))*(-F_w_u[i-1][J-1]); 
                A_WW1 = alpha_w_u[i-1][J-1]*UU(r_w_u[i-1][J-1])*(-F_w_u[i-1][J-1]); 
                A_W1 = ((1-alpha_w_u[i-1][J-1])*(D(r_w_u[i-1][J-1]))+alpha_w_u[i-1][J-1]*(U(r_w_u[i-1][J-1])-1))*(-F_w_u[i-1][J-1]); 
                A_E1 = (1-alpha_w_u[i-1][J-1])*UU(r_w_u[i-1][J-1])*(-F_w_u[i-1][J-1]); 
                S1[i-1][J-1] += -(-p_star[i-2][J-1])*dy/rho  + ((i == 2) ? 0 : 
                -u_star_old[i-1][J]*A_P1
                -u_star_old[i-3][J]*A_WW1
                -u_star_old[i-2][J]*A_W1
                -u_star_old[i][J]*A_E1 ); 
                //(2)east side
                double A_P2, A_W2, A_E2, A_EE2 ; //trancient variables
                A_P2 = F_e_u[i-1][J-1]*((1-alpha_e_u[i-1][J-1])*D(r_e_u[i-1][J-1])+alpha_e_u[i-1][J-1]*(U(r_e_u[i-1][J-1])-1)); 
                A_W2 = F_e_u[i-1][J-1]*alpha_e_u[i-1][J-1]*UU(r_e_u[i-1][J-1]); 
                A_E2 = F_e_u[i-1][J-1]*((1-alpha_e_u[i-1][J-1])*(U(r_e_u[i-1][J-1])-1)+alpha_e_u[i-1][J-1]*(D(r_e_u[i-1][J-1]))); 
                A_EE2 = F_e_u[i-1][J-1]*(1-alpha_e_u[i-1][J-1])*UU(r_e_u[i-1][J-1]);
                S1[i-1][J-1] += -(p_star[i-1][J-1])*dy/rho + ((i == Nx) ? 0 : 
                -u_star_old[i-1][J]*A_P2 
                -u_star_old[i-2][J]*A_W2
                -u_star_old[i][J]*A_E2
                -u_star_old[i+1][J]*A_EE2 ); 
                //(3)south side 
                double A_P3, A_SS3, A_S3, A_N3 ;
                A_P3 = ((1-alpha_s_u[i-1][J-1])*(U(r_s_u[i-1][J-1])-1)+alpha_s_u[i-1][J-1]*(D(r_s_u[i-1][J-1])))*(-F_s_u[i-1][J-1]); 
                A_SS3 = alpha_s_u[i-1][J-1]*UU(r_s_u[i-1][J-1])*(-F_s_u[i-1][J-1]); 
                A_S3 = ((1-alpha_s_u[i-1][J-1])*(D(r_s_u[i-1][J-1]))+alpha_s_u[i-1][J-1]*(U(r_s_u[i-1][J-1])-1))*(-F_s_u[i-1][J-1]); 
                A_N3 = (1-alpha_s_u[i-1][J-1])*UU(r_s_u[i-1][J-1])*(-F_s_u[i-1][J-1]); 
                S1[i-1][J-1] += ((J == 1) ? -u_bottom*(-F_s_u[i-1][J-1]-2*D_s) : 
                -u_star_old[i-1][J]*A_P3 
                -u_star_old[i-1][J-2]*A_SS3 //have been considered the extrapolation
                -u_star_old[i-1][J-1]*A_S3
                -u_star_old[i-1][J+1]*A_N3 );
                //(4)norht side
                double A_P4, A_S4, A_N4, A_NN4 ; //trancient variables
                A_P4 = F_n_u[i-1][J-1]*((1-alpha_n_u[i-1][J-1])*D(r_n_u[i-1][J-1])+alpha_n_u[i-1][J-1]*(U(r_n_u[i-1][J-1])-1)); 
                A_S4 = F_n_u[i-1][J-1]*alpha_n_u[i-1][J-1]*UU(r_n_u[i-1][J-1]); 
                A_N4 = F_n_u[i-1][J-1]*((1-alpha_n_u[i-1][J-1])*(U(r_n_u[i-1][J-1])-1)+alpha_n_u[i-1][J-1]*(D(r_n_u[i-1][J-1]))); 
                A_NN4 = F_n_u[i-1][J-1]*(1-alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1]);
                S1[i-1][J-1] += ((J == Ny) ? -u_top*(F_n_u[i-1][J-1]-2*D_n) :  
                -u_star_old[i-1][J]*A_P4
                -u_star_old[i-1][J-1]*A_S4
                -u_star_old[i-1][J+1]*A_N4
                -u_star_old[i-1][J+2]*A_NN4) ;   
            }
        }
    }
    //Using ADI iteration to solve (i = 2 : Nx ; J = 1 : Ny)
    //ADI-1 scaning direction from left to right and from bottom to top 
    //parameters setting : beta, D, alpha x, S, U_left, U_right , n
    #pragma omp parallel default(shared)
    {
        #pragma omp for schedule(static) //average cutting workking for threads 
        for(int J = 1 ; J <= Ny ; J++){//define the index from 0~Nx-2
        vector <double> beta (Nx-1) ;//left hand side elemets 
        vector <double> D (Nx-1) ; //primary elements
        vector <double> alpha (Nx-1) ; //right and side elements
        vector <double> x (Nx-1) ; //solutions 
        vector <double> S (Nx-1) ; //Source Vector
        double U_left = u_left ; 
        double U_right = u_right ; 
        for(int k = 2 ; k <= Nx ; k++){
            beta[k-2] = a_W[k-1][J-1] ; 
            D[k-2] = a_P[k-1][J-1] ; 
            alpha[k-2] = a_E[k-1][J-1] ; 
            x[k-2] = u_star[k-1][J] ; 
            S[k-2] = S1[k-1][J-1]-u_star_old[k-1][J-1]*a_S[k-1][J-1]-u_star_old[k-1][J+1]*a_N[k-1][J-1];
        }
        TDMA(D , beta , alpha , x , S , Nx-1 , U_left , U_right) ; 
            for(int i = 2 ; i <= Nx ; i++){
                u_star[i-1][J] =   x[i-2] ;
            }
        }
        #pragma omp for collapse(2) schedule(static) //average cutting workking for threads 
        for(int i = 1 ; i <= Nx+1 ; i++){
            for(int J = 0 ; J <= Ny+1 ; J++){
                u_star_old2[i-1][J] = u_star[i-1][J] ; 
            }
        }//Global elements update to the old value2 
        //ADI-2 scaning direction from bottom to top and from left to right
        //parameter setting : beta , D , alpha x , S , U_left , U_right , n 
        #pragma omp for schedule(static) //average cutting workking for threads 
        for(int i = 2 ; i <=Nx ; i++){//define the index from 1 ~ Ny
            vector <double> beta (Ny) ; vector <double> D (Ny) ; vector <double> alpha (Ny) ;
            vector <double> x (Ny) ; vector <double> S (Ny) ; //we would not worry about the index > Ny because we haved set the extrapolation 
            for(int k = 1 ; k <= Ny ; k++){
            beta[k-1] = a_S[i-1][k-1] ; 
            D[k-1] = a_P[i-1][k-1] ; 
            alpha[k-1] = a_N[i-1][k-1] ; 
            x[k-1] = u_star[i-1][k] ; 
            S[k-1] = S1[i-1][k-1]-u_star_old2[i-2][k]*a_W[i-1][k-1]-u_star_old2[i][k]*a_E[i-1][k-1];
            }
            TDMA_FVM(D , beta , alpha , x , S , Ny) ; 
            for(int J = 1 ; J <= Ny ; J++){
                u_star[i-1][J] =  x[J-1] ;
                u_star[i-1][J] = omega_u*u_star[i-1][J] + (1-omega_u)*u_star_old2[i-1][J] ; //SOR relaxation
            }
        }
    }//this part is about ADIiteration to solve the x-momentum equation 


    //solving y part .....
    #pragma omp parallel for collapse(2) 
    for(int j = 1 ; j <= Ny+1 ; j++){//writting the range of the coefficient not about the v_star 
        for(int I = 1 ; I <= Nx ; I++){
            b_P[I-1][j-1] = 0.0 ; b_W[I-1][j-1] = 0.0 ; b_E[I-1][j-1] = 0.0 ; b_S[I-1][j-1] = 0.0 ; b_N[I-1][j-1] = 0.0 ; S2[I-1][j-1] = 0.0 ; 
            if(j == 1){
                v_star[I][j-1] = v_bottom ; 
            }else if(j == Ny+1){
                v_star[I][j-1] = v_top ; 
            }else{
                //L.H.S : Using Upwind schemne for non-special case
                //(1)south side 
                b_S[I-1][j-1] += -D_s + (-F_s_v[I-1][j-1])*1*(alpha_s_v[I-1][j-1]);
                b_P[I-1][j-1] += +D_s + (-F_s_v[I-1][j-1])*1*(1-alpha_s_v[I-1][j-1]);
                //(2)north side
                b_N[I-1][j-1] += -D_n + (F_n_v[I-1][j-1])*(1-alpha_n_v[I-1][j-1]) ; 
                b_P[I-1][j-1] += +D_n + (F_n_v[I-1][j-1])*(alpha_n_v[I-1][j-1]) ; 
                //(3)west side
                b_W[I-1][j-1] += ((I == 1) ? 0.0 : -D_w + (-F_w_v[I-1][j-1])*1*(alpha_w_v[I-1][j-1])) ;
                b_P[I-1][j-1] += ((I == 1) ? 2*D_w : D_w + (-F_w_v[I-1][j-1])*1*(1-alpha_w_v[I-1][j-1])) ;
                //(4)east side
                b_E[I-1][j-1] += ((I == Nx) ? 0.0 : -D_e + (F_e_v[I-1][j-1])*1*(1-alpha_e_v[I-1][j-1])) ;
                b_P[I-1][j-1] += ((I == Nx) ? 2*D_e : D_e + (F_e_v[I-1][j-1])*1*(alpha_e_v[I-1][j-1])) ;
                //R.H.S(source term) -------deferred correction 
                //(1)south side  
                double B_P1 , B_SS1 , B_S1 , B_N1 ; 
                B_P1 = ((1-alpha_s_v[I-1][j-1])*(U(r_s_v[I-1][j-1])-1)+alpha_s_v[I-1][j-1]*D(r_s_v[I-1][j-1]))*(-F_s_v[I-1][j-1]) ;
                B_SS1 = (alpha_s_v[I-1][j-1]*UU(r_s_v[I-1][j-1]))*(-F_s_v[I-1][j-1]) ; 
                B_S1 = ((1-alpha_s_v[I-1][j-1])*D(r_s_v[I-1][j-1])+alpha_s_v[I-1][j-1]*(U(r_s_v[I-1][j-1])-1))*(-F_s_v[I-1][j-1]) ;
                B_N1 = ((1-alpha_s_v[I-1][j-1])*UU(r_s_v[I-1][j-1]))*(-F_s_v[I-1][j-1]) ; 
                S2[I-1][j-1] += -(-p_star[I-1][j-2])*dx/rho + ((j == 2) ? 0.0 : 
                -v_star_old[I][j-1]*B_P1
                -v_star_old[I][j-3]*B_SS1
                -v_star_old[I][j-2]*B_S1
                -v_star_old[I][j]*B_N1                 
                );
                //(2)north side 
                double B_P2 , B_S2 , B_N2 , B_NN2 ; 
                B_P2 = ((1-alpha_n_v[I-1][j-1])*D(r_n_v[I-1][j-1])+alpha_n_v[I-1][j-1]*(U(r_n_v[I-1][j-1])-1))*F_n_v[I-1][j-1] ; 
                B_S2 = alpha_n_v[I-1][j-1]*UU(r_n_v[I-1][j-1])*F_n_v[I-1][j-1];  
                B_N2 = ((1-alpha_n_v[I-1][j-1])*(U(r_n_v[I-1][j-1])-1)+alpha_n_v[I-1][j-1]*(D(r_n_v[I-1][j-1])))*F_n_v[I-1][j-1] ; 
                B_NN2 = (1-alpha_n_v[I-1][j-1])*UU(r_n_v[I-1][j-1])*F_n_v[I-1][j-1];  
                S2[I-1][j-1] += -(p_star[I-1][j-1])*dx/rho + ((j == Ny) ? 0.0 : 
                -v_star_old[I][j-1]*B_P2
                -v_star_old[I][j-2]*B_S2
                -v_star_old[I][j]*B_N2
                -v_star_old[I][j+1]*B_NN2
                );
                //(3)west side 
                double B_P3 , B_WW3 , B_W3 , B_E3 ; 
                B_P3 = ((1-alpha_w_v[I-1][j-1])*(U(r_w_v[I-1][j-1])-1)+alpha_w_v[I-1][j-1]*D(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1]);
                B_WW3 = (alpha_w_v[I-1][j-1]*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1]); 
                B_W3 = ((1-alpha_w_v[I-1][j-1])*D(r_w_v[I-1][j-1])+alpha_w_v[I-1][j-1]*(U(r_w_v[I-1][j-1])-1))*(-F_w_v[I-1][j-1]); 
                B_E3 = ((1-alpha_w_v[I-1][j-1])*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1]);  
                S2[I-1][j-1] += ((I == 1) ?  -v_left*(-2*D_w+(-F_w_v[I-1][j-1])) : 
                -v_star_old[I][j-1]*B_P3 
                -v_star_old[I-2][j-1]*B_WW3
                -v_star_old[I-1][j-1]*B_W3
                -v_star_old[I+1][j-1]*B_E3
               );
               //(4)east side 
               double B_P4 , B_W4 , B_E4 , B_EE4 ; 
               B_P4 = ((1-alpha_e_v[I-1][j-1])*D(r_e_v[I-1][j-1])+alpha_e_v[I-1][j-1]*(U(r_e_v[I-1][j-1])-1))*F_e_v[I-1][j-1]; 
               B_W4 = alpha_e_v[I-1][j-1]*UU(r_e_v[I-1][j-1])*F_e_v[I-1][j-1]; 
               B_E4 = ((1-alpha_e_v[I-1][j-1])*(U(r_e_v[I-1][j-1])-1)+alpha_e_v[I-1][j-1]*(D(r_e_v[I-1][j-1])))*F_e_v[I-1][j-1]; 
               B_EE4 = (1-alpha_e_v[I-1][j-1])*UU(r_e_v[I-1][j-1])*F_e_v[I-1][j-1]; 
               S2[I-1][j-1] += ((I == Nx) ? -v_right*(-2*D_e+F_e_v[I-1][j-1]) : 
               -v_star_old[I][j-1]*B_P4 
               -v_star_old[I-1][j-1]*B_W4
               -v_star_old[I+1][j-1]*B_E4
               -v_star_old[I+2][j-1]*B_EE4
               ) ; 
            }
        }
    }
    //using ADI + TDMA algorithm to compute the y-Velocity predicted value 
    #pragma omp parallel default(shared) // create the thream team and share the data memory for different threads
    {   //computing range : I = [1,Nx] , j = [2,Ny]
        //ADI-1 scanning direction is from left to right and from bottom to top
        #pragma omp for schedule(static) //average cutting the working spaxc for my threads
        for(int j = 2 ; j <= Ny ; j++){
            vector <double> beta(Nx) ; vector<double> D(Nx) ; vector<double> alpha(Nx) ; 
            vector <double> x(Nx) ; vector<double> S(Nx) ; //Solving Range from 1 to Nx because it is stiraged at the center of the Grid
            for(int k = 1 ; k <= Nx ; k++){//For ADI-1 west and east using the new value and north and south using old value
            beta[k-1] = b_W[k-1][j-1] ; 
            D[k-1] = b_P[k-1][j-1] ; 
            alpha[k-1] = b_E[k-1][j-1] ; 
            x[k-1] = v_star[k][j-1] ; 
            S[k-1] = S2[k-1][j-1] -v_star_old[k][j-2]*b_S[k-1][j-1]-v_star_old[k][j]*b_N[k-1][j-1] ;   
            }
            TDMA_FVM(D , beta , alpha , x , S , Nx ) ;
            for(int I = 1 ; I <= Nx ; I++){
                v_star[I][j-1] = x[I-1]  ; 
                //predicted value doing relaxation factor
            }
        }
        #pragma omp for collapse(2) schedule(static) //average cutting workking for threads 
        for(int I = 0 ; I <= Nx+1 ; I++){
            for(int j = 1 ; j <= Ny+1 ; j++){
                v_star_old2[I][j-1] = v_star[I][j-1] ; 
            }
        } 
        #pragma omp for schedule(static)
        //ADI-2 from bottom to top and from left to right scaning 
        for(int I = 1 ; I <= Nx ; I++){
            vector <double> beta(Ny-1) ; 
            vector <double> D(Ny-1) ; 
            vector <double> alpha(Ny-1) ; 
            vector <double> x(Ny-1) ; 
            vector <double> S(Ny-1) ; 
            double U_left = v_bottom  ; //Boundary Condition 
            double U_right = v_top ;    //BOundary Condition
            for(int k = 2 ; k <= Ny ; k++){
                beta[k-2] = b_S[I-1][k-1] ; 
                D[k-2] = b_P[I-1][k-1] ; 
                alpha[k-2] = b_N[I-1][k-1] ; 
                x[k-2] = v_star[I][k-1] ; 
                S[k-2] = S2[I-1][k-1] - v_star_old2[I-1][k-1]*b_W[I-1][k-1] - v_star_old2[I+1][k-1]*b_E[I-1][k-1] ; 
            }
            TDMA(D , beta , alpha , x , S , Ny-1 , U_left , U_right) ;
            for(int j = 2 ; j <= Ny ; j++){
                v_star[I][j-1] = x[j-2]  ; 
                v_star[I][j-1] = omega_v * v_star[I][j-1] + (1-omega_v)*v_star_old[I][j-1] ;
            }
        }
    }
}
//step4 : updating three things : velocity correction , correct pressure , and correct velocity ; /....
void step4 (){
    //(1)updating velocity correction : 
    #pragma omp parallel default(shared)
    {
        #pragma omp for collapse (2) schedule(static)
        //computing the y-part
        for(int J = 1 ; J <= Ny ; J++){
            for(int i = 1 ; i <= Nx+1 ; i++){
                double delta_Px = ((i == 1 || i == Nx+1) ? 0.0 : -(p_prime[i-1][J-1] - p_prime[i-2][J-1])) ; 
                u_prime[i-1][J-1] = ((i == 1 || i == Nx+1) ? 0.0 : delta_Px*(dy/(rho*a_P[i-1][J-1]))) ;
                //updating correct velocity ......
                u[i-1][J-1] = ((i == 1 || i == Nx+1) ? 0.0 : u_star[i-1][J] + u_prime[i-1][J-1]) ;  
            }
        }//Using SIMPLE omission...
        #pragma omp for collapse (2) schedule(static)
        //computing the x-part
        for(int I = 1 ; I <= Nx ; I++){
            for(int j = 1 ; j <= Ny+1 ; j++){
                double delta_Py = ((j == 1 || j == Ny+1) ? 0.0 : -(p_prime[I-1][j-1] - p_prime[I-1][j-2])) ;
                v_prime[I-1][j-1] = ((j == 1 || j == Ny+1) ? 0.0 : delta_Py*(dx/(rho*b_P[I-1][j-1]))) ; 
                //updating the correct velocity 
                v[I-1][j-1] = ((j == 1 || j == Ny+1) ? 0.0 : v_star[I][j-1] + v_prime[I-1][j-1]) ; 
            }
        }//Using SIMPLE omission...
        #pragma omp for collapse (2) schedule(static)
        //updating the correct pressure
        for(int I = 1 ; I <= Nx ; I++){
            for(int J = 1 ; J <= Ny ; J++){
                p[I-1][J-1] = p_star[I-1][J-1] + omega_P * p_prime[I-1][J-1] ; 
            }
        }//pressure predicted and corrected method
    }
}

double step5(int m){
    double sum = 0.0 ;
    sum += fabs(u[80][0]-0) ; 
    sum += fabs(u[80][8]-(-0.41165)) ;
    sum += fabs(u[80][9]-(-0.42901)) ; 
    sum += fabs(u[80][10]-(-0.43643)) ;  
    sum += fabs(u[80][15]-(-0.40435)) ;  
    sum += fabs(u[80][27]-(-0.3305)) ;  
    sum += fabs(u[80][44]-(-0.22855)) ;
    sum += fabs(u[80][71]-(-0.07404)) ;
    sum += fabs(u[80][79]-(-0.03039)) ;
    sum += fabs(u[80][98]-(0.08183)) ; 
    sum += fabs(u[80][117]-(0.20087)) ;  
    sum += fabs(u[80][135]-(0.33556)) ;
    sum += fabs(u[80][151]-(0.46036)) ;
    sum += fabs(u[80][153]-(0.45992)) ;
    sum += fabs(u[80][154]-(0.4612)) ;
    sum += fabs(u[80][155]-(0.48223)) ;
    sum += fabs(u[80][159]-(1)) ;
    sum += fabs(u[80][0]-0) ; 
    sum += fabs(v[0][80]-0) ;
    sum += fabs(v[9][80]-(0.42447)) ; 
    sum += fabs(v[10][80]-(0.43329)) ;
    sum += fabs(v[12][80]-(0.43648)) ; 
    sum += fabs(v[14][80]-(0.42951)) ;
    sum += fabs(v[24][80]-(0.35368)) ; 
    sum += fabs(v[35][80]-(0.28066)) ;
    sum += fabs(v[37][80]-(0.2728)) ;
    sum += fabs(v[80][80]-(0.00945)) ;
    sum += fabs(v[129][80]-(-0.30018)) ; 
    sum += fabs(v[137][80]-(-0.36214)) ;
    sum += fabs(v[145][80]-(-0.41442)) ; 
    sum += fabs(v[151][80]-(-0.52876)) ;
    sum += fabs(v[152][80]-(-0.55408)) ; 
    sum += fabs(v[154][80]-(-0.55069)) ;
    sum += fabs(v[155][80]-(-0.49774)) ;
    sum += fabs(v[159][80]-(0.0)) ;
    /*double sum = 0.0 ;
    sum += fabs(u[40][0]-0) ; 
    sum += fabs(u[40][3]-(-0.41165)) ;
    sum += fabs(u[40][4]-(-0.42901)) ; 
    sum += fabs(u[40][5]-(-0.43643)) ;  
    sum += fabs(u[40][7]-(-0.40435)) ;  
    sum += fabs(u[40][13]-(-0.3305)) ;  
    sum += fabs(u[40][22]-(-0.22855)) ;
    sum += fabs(u[40][35]-(-0.07404)) ;
    sum += fabs(u[40][39]-(-0.03039)) ;
    sum += fabs(u[40][48]-(0.08183)) ; 
    sum += fabs(u[40][58]-(0.20087)) ;  
    sum += fabs(u[40][67]-(0.33556)) ;
    sum += fabs(u[40][75]-(0.46036)) ;
    sum += fabs(u[40][76]-(0.45992)) ;
    sum += fabs(u[40][77]-(0.4612)) ;
    sum += fabs(u[40][77]-(0.48223)) ;
    sum += fabs(u[40][79]-(1)) ;
    sum += fabs(u[40][0]-0) ; 
    sum += fabs(v[0][40]-0) ;
    sum += fabs(v[4][40]-(0.42447)) ; 
    sum += fabs(v[5][40]-(0.43329)) ;
    sum += fabs(v[5][40]-(0.43648)) ; 
    sum += fabs(v[7][40]-(0.42951)) ;
    sum += fabs(v[12][40]-(0.35368)) ; 
    sum += fabs(v[17][40]-(0.28066)) ;
    sum += fabs(v[18][40]-(0.2728)) ;
    sum += fabs(v[39][40]-(0.00945)) ;
    sum += fabs(v[64][40]-(-0.30018)) ; 
    sum += fabs(v[69][40]-(-0.36214)) ;
    sum += fabs(v[72][40]-(-0.41442)) ; 
    sum += fabs(v[76][40]-(-0.52876)) ;
    sum += fabs(v[76][40]-(-0.55408)) ; 
    sum += fabs(v[77][40]-(-0.55069)) ;
    sum += fabs(v[77][40]-(-0.49774)) ;
    sum += fabs(v[79][40]-(0.0)) ;*/
    return sum ;
}

double d(int I, int J , int m){
    double u_e = u_star[I][J] ;
    double u_w = u_star[I-1][J] ;
    double v_n = v_star[I][J] ;
    double v_s = v_star[I][J-1] ;
    return rho*(u_e - u_w)*dy + rho*(v_n - v_s)*dx ;
} 

double step5star(int m){
    double sum = 0.0 ;
    #pragma omp parallel for collapse(2) reduction(+:sum)  // *** 加這行 ***
    for(int I = 1 ; I <= Nx ; I++){
        for(int J = 1 ; J <= Ny ; J++){
            sum = sum + fabs(d(I,J,m)) ;
        }
    }
    return sum ;
}


void step6(){
    #pragma omp parallel default(shared)
    {
        //1.u_star,u_star_old,v_star,v_star_old,p_star updating
        //I = [1,NX+1] , j = [0.Ny+1] 
        #pragma omp for collapse (2) schedule(static)//(1)
        for(int i = 1 ; i <= Nx+1 ; i++){
            for(int J = 1 ; J <= Ny ; J++){
                //initialize the u_star 
                if(i == 1 || i == Nx+1){
                    u_star[i-1][J] = ((i == 1) ? u_left : u_right) ; 
                }else {
                    u_star[i-1][J] = u[i-1][J-1] ; //second we decide the inner node
                }
            }
        }
        //given value for the extrapolation
        #pragma omp for schedule(static)
        for(int i = 1 ; i <= Nx+1 ; i++){
            u_star[i-1][0] = 2*u_bottom - u_star[i-1][1] ; 
            u_star[i-1][Ny+1] = 2*u_top - u_star[i-1][Ny] ;
        }
        #pragma omp for collapse (2) schedule(static)//(1)
        for(int i = 1 ; i <= Nx+1 ; i++){
            for(int J = 0 ; J <= Ny+1 ; J++){
                u_star_old[i-1][J] = u_star[i-1][J] ; 
            }
        }
        #pragma omp for collapse (2) schedule(static)//(2)
        //i = [0,Nx+1] , j = [1,Ny+1]
        for(int I = 1 ; I <= Nx ; I++){
            for(int j = 1 ; j <= Ny+1 ; j++){
                //initialize the v_ star 
                if(j == 1 || j == Ny+1){
                    v_star[I][j-1] = ((j == 1) ? v_bottom : v_top) ;
                }else {
                    v_star[I][j-1] = v[I-1][j-1] ; //second we decide the inner node
                }
            }
        }
        //given value for the extrapolation
        #pragma omp for schedule(static)
        for(int j = 1 ; j <= Ny+1 ; j++){
            v_star[0][j-1] = 2*v_left - v_star[1][j-1] ; 
            v_star[Nx+1][j-1] = 2*v_right - v_star[Nx][j-1] ;
        }
        #pragma omp for collapse (2) schedule(static)//(2)
        //i = [0,Nx+1] , j = [1,Ny+1]
        for(int I = 0 ; I <= Nx+1 ; I++){
            for(int j = 1 ; j <= Ny+1 ; j++){
                v_star_old[I][j-1] = v_star[I][j-1] ;
            }
        }
        #pragma omp for collapse (2) schedule(static)//(3)
        //I = [1,Nx] , J = [1,Ny]//updating the pressure predicted value 
        for(int I = 1 ; I <= Nx ; I++){
            for(int J = 1 ; J <= Ny ; J++){
                p_star[I-1][J-1] = p[I-1][J-1] ; 
            }
        }
        #pragma omp for collapse (2) schedule(static)//(4)updating the coefficient of convection term and alpha and side anaylized factor ----u 
        for(int i = 1 ; i <= Nx+1 ; i++){
            for(int J = 1 ; J <= Ny ; J++){
                r_w_u[i-1][J-1] = 0.0 ; r_e_u[i-1][J-1] = 0.0 ;  r_s_u[i-1][J-1] = 0.0 ; r_n_u[i-1][J-1] = 0.0 ;    
                F_w_u[i-1][J-1] = ((i == 1 || i == Nx+1) ? 0.0 : dy*.5*(u_star_old[i][J]+u_star_old[i-1][J])) ;
                alpha_w_u[i-1][J-1] = ((F_w_u[i-1][J-1] > 0) ? 1.0 : 0.0 ) ;
                double A = alpha_w_u[i-1][J-1] ;  
                if(i >= 3 && i <= Nx){
                    double P = u_star_old[i-1][J] ; double W = u_star_old[i-2][J] ; double WW = u_star_old[i-3][J] ; double E = u_star_old[i][J] ;
                    double D1 = P*A + (1-A)*W ; 
                    double U1 = W*A + (1-A)*P ; 
                    double UU1 = WW*A + (1-A)*E ; 
                    r_w_u[i-1][J-1] = (U1-UU1)/(D1-U1+1e-30) ; 
                }else{r_w_u[i-1][J-1] = 0.0 ;}
                F_e_u[i-1][J-1] = ((i == 1 || i == Nx+1) ? 0.0 : dy*0.5*(u_star_old[i-2][J]+u_star_old[i-1][J])) ;
                alpha_e_u[i-1][J-1] = ((F_e_u[i-1][J-1] > 0) ? 1.0 : 0.0 ) ; 
                if(i >= 2 && i <= Nx-1 ){
                    double P = u_star_old[i-1][J] ; double W = u_star_old[i-2][J] ; double E = u_star_old[i][J] ; double EE = u_star_old[i+1][J] ;
                    double B = alpha_e_u[i-1][J-1] ;
                    double D2 = B*E + (1-B)*P ; 
                    double U2 = B*P + (1-B)*E ; 
                    double UU2 = B*W + (1-B)*EE ; 
                    r_e_u[i-1][J-1] = (U2-UU2)/(D2-U2+1e-30) ;  
                }else{r_e_u[i-1][J-1]= 0.0 ;}
                F_n_u[i-1][J-1] = ((i == 1 || i == Nx+1) ? 0.0 : dx*0.5*(v_star_old[i-1][J]+v_star_old[i][J])) ;
                alpha_n_u[i-1][J-1] = ((F_n_u[i-1][J-1] > 0) ? 1.0 : 0.0 ) ;
                if(i>= 2 && i <= Nx && J <= Ny-1){
                    double P = u_star_old[i-1][J] ; double S = u_star_old[i-1][J-1] ; double N = u_star_old[i-1][J+1] ; double NN = u_star_old[i-1][J+2] ;
                    double C = alpha_n_u[i-1][J-1] ; 
                    double D3 = C*N + (1-C)*P ; 
                    double U3 = C*P + (1-C)*N ; 
                    double UU3 = C*S + (1-C)*NN ;
                    r_n_u[i-1][J-1] =(U3-UU3)/(D3-U3+1e-30) ;  
                }else{r_n_u[i-1][J-1] =0.0 ;}
                F_s_u[i-1][J-1] = ((i == 1 || i == Nx+1) ? 0.0 : dx*0.5*(v_star_old[i-1][J-1]+v_star_old[i][J-1])) ;
                alpha_s_u[i-1][J-1] = ((F_s_u[i-1][J-1] > 0) ? 1.0 : 0.0 ) ;
                if(i>= 2 && i <= Nx && J >= 2){
                    double P = u_star_old[i-1][J] ; double S = u_star_old[i-1][J-1] ; double SS = u_star_old[i-1][J-2] ; double N = u_star_old[i-1][J+1]  ;
                    double D = alpha_s_u[i-1][J-1] ; 
                    double D4 = D*P + (1-D)*S ;
                    double U4 = D*S + (1-D)*P ; 
                    double UU4 = D*SS + (1-D)*N ;
                    r_s_u[i-1][J-1] = (U4-UU4)/(D4-U4+1e-30)  ; 
                }else{r_s_u[i-1][J-1] = 0.0 ;}
            }
        }
        #pragma omp for collapse (2) schedule(static)//(5)
        for(int I = 1 ; I <= Nx ; I++){
            for(int j = 1 ; j <= Ny+1 ; j++){
                r_w_v[I-1][j-1] = 0.0 ; r_e_v[I-1][j-1] = 0.0 ;  r_s_v[I-1][j-1] = 0.0 ; r_n_v[I-1][j-1] = 0.0 ;
                F_w_v[I-1][j-1] = ((j == 1 || j == Ny+1) ? 0.0 : dy*0.5*(u_star_old[I-1][j]+u_star_old[I-1][j-1])) ;
                alpha_w_v[I-1][j-1] = ((F_w_v[I-1][j-1] > 0) ? 1.0 : 0.0 ) ; 
                if(j >= 2 && j <= Ny && I >= 2){
                    double P = v_star_old[I][j-1] ; double W = v_star_old[I-1][j-1] ; double WW = v_star_old[I-2][j-1] ; double E = v_star_old[I+1][j-1] ;
                    double A = alpha_w_v[I-1][j-1] ; double D1 = P*A + (1-A)*W ; double U1 = W*A + (1-A)*P ; double UU1 = WW*A + (1-A)*E ;
                    r_w_v[I-1][j-1] = (U1-UU1)/(D1-U1+1e-30) ; 
                }else{r_w_v[I-1][j-1] = 0.0 ;}
                F_e_v[I-1][j-1] = ((j == 1 || j == Ny+1) ? 0.0 : dy*0.5*(u_star_old[I][j]+u_star_old[I][j-1])) ; 
                alpha_e_v[I-1][j-1] = ((F_e_v[I-1][j-1] > 0) ? 1.0 : 0.0 ) ; 
                if(j >= 2 && j <= Ny && I <= Nx-1){
                    double P = v_star_old[I][j-1] ; double W = v_star_old[I-1][j-1] ; double E = v_star_old[I+1][j-1] ; double EE = v_star_old[I+2][j-1] ;
                    double B = alpha_e_v[I-1][j-1] ; double D2 = B*E + (1-B)*P ; double U2 = B*P + (1-B)*E ; double UU2 = B*W + (1-B)*EE ;
                    r_e_v[I-1][j-1] = (U2-UU2)/(D2-U2+1e-30) ;
                }else{r_e_v[I-1][j-1]=0.0;}
                F_s_v[I-1][j-1] = ((j == 1 || j == Ny+1) ? 0.0 : dx*0.5*(v_star_old[I][j-1]+v_star_old[I][j-2])) ; 
                alpha_s_v[I-1][j-1] = ((F_s_v[I-1][j-1] > 0) ? 1.0 : 0.0 ) ;
                if(j >= 3 && j <= Ny){
                    double P = v_star_old[I][j-1] ; double S = v_star_old[I][j-2] ; double SS = v_star_old[I][j-3] ; double N = v_star_old[I][j] ;
                    double C = alpha_s_v[I-1][j-1] ; double D3 = C*P + (1-C)*S ; double U3 = C*S + (1-C)*P ; double UU3 = C*SS + (1-C)*N ; 
                    r_s_v[I-1][j-1] = (U3-UU3)/(D3-U3+1e-30);
                }else{r_s_v[I-1][j-1] = 0.0 ; }
                F_n_v[I-1][j-1] = ((j == 1 || j == Ny+1) ? 0.0 : dx*0.5*(v_star_old[I][j-1]+v_star_old[I][j])) ; 
                alpha_n_v[I-1][j-1] = ((F_n_v[I-1][j-1] > 0) ? 1.0 : 0.0 ) ; 
                if(j >= 2 && j <= Ny-1){
                    double P = v_star_old[I][j-1] ; double S = v_star_old[I][j-2]; double N = v_star_old[I][j] ; double NN = v_star_old[I][j+1] ;
                    double D = alpha_n_v[I-1][j-1] ; double D4 = D*N + (1-D)*P ; double U4 = D*P + (1-D)*N ; double UU4 = D*S + (1-D)*NN ;
                    r_n_v[I-1][j-1] = (U4-UU4)/(D4-U4+1e-30) ; 
                }else{r_n_v[I-1][j-1] = 0.0 ;}
            }
        }
    }
}


void output_point_fields(int m) {
    int nodes_i = Nx + 1;
    int nodes_j = Ny + 1;

    std::ostringstream name;
    name << "3.CDs_Re" << Re << "_Grid" << Nx+1 << "x" << Ny+1
         << "_DC_step" << std::setfill('0') << std::setw(6) << m 
         << "_error" << min_step7 <<
          ".vtk";

    std::ofstream out(name.str());
    
    out << "# vtk DataFile Version 3.0\n";
    out << "3.CDs_Re = 5000_Grid"<<Nx+1<<"x"<<Nx+1<<"\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_GRID\n";
    out << "DIMENSIONS " << nodes_i << " " << nodes_j << " 1\n";

    // ------------------------------
    // Point coordinates
    // ------------------------------
    out << "POINTS " << nodes_i * nodes_j << " double\n";
    for (int j = 0; j < nodes_j; j++) {
        for (int i = 0; i < nodes_i; i++) {
            double x = i * dx;
            double y = j * dy;
            out << x << " " << y << " 0.0\n";
        }
    }

    // ---------------------------------
    // (A) Pressure → point data (interpolation)
    // ---------------------------------
    out << "POINT_DATA " << nodes_i * nodes_j << "\n";
    out << "SCALARS Pressure double 1\n";
    out << "LOOKUP_TABLE default\n";

    for (int j = 0; j < nodes_j; j++) {
        for (int i = 0; i < nodes_i; i++) {

            double p_val = 0.0;

            // Bilinear interpolation over surrounding 4 cells
            int I = i;
            int J = j;

            if (I >= 1 && I < Nx && J >= 1 && J < Ny) {
                p_val = 0.25 * (
                      p_star[I-1][J-1]
                    + p_star[I][J-1]
                    + p_star[I-1][J]
                    + p_star[I][J]
                );
            }
            else {
                // boundary：nearest valid cell
                int ii = std::min(std::max(I,1), Nx) - 1;
                int jj = std::min(std::max(J,1), Ny) - 1;
                p_val = p_star[ii][jj];
            }

            out << p_val << "\n";
        }
    }

    // ---------------------------------
    // (B) Velocity → point data (interpolation)
    // ---------------------------------
    out << "VECTORS Velocity double\n";

    for (int j = 0; j < nodes_j; j++) {
        for (int i = 0; i < nodes_i; i++) {

            double u_val = 0.0, v_val = 0.0;

            // interpolate u (vertical faces)
            if (i > 0 && i < Nx && j > 0 && j < Ny) {
                u_val = 0.25 * (
                      u_star[i-1][j]
                    + u_star[i][j]
                    + u_star[i-1][j+1]
                    + u_star[i][j+1]
                );
            }

            // interpolate v (horizontal faces)
            if (i > 0 && i < Nx && j > 0 && j < Ny) {
                v_val = 0.25 * (
                      v_star[i][j-1]
                    + v_star[i][j]
                    + v_star[i+1][j-1]
                    + v_star[i+1][j]
                );
            }

            out << u_val << " " << v_val << " 0.0\n";
        }
    }
    out.close();
    cout << "VTK POINT_DATA has been output: " << name.str() << endl;
}
int main(){
    //using OpenMP to Go Fast !!
    int num_threads = omp_get_max_threads();
    omp_set_num_threads(num_threads);
    cout << "CDs Using " << num_threads << " threads" << endl;
    cout << "Grid Size:" << Nx << "x" << Ny << endl;
    cout << "Grid spacing : " << dx << "x " << dy << endl;
    cout << "Boundary Conditions: " ;
    cout << "u: left = " << u_left << ", right = " << u_right << ", bottom = " << u_bottom << ", top = " << u_top << endl;
    cout << "v: left = " << v_left << ", right = " << v_right << ", bottom = " << v_bottom << ", top = " << v_top << endl;
    cout << "Re = " << Re << ",rho = " << rho << endl ;
    CE_physical = false;
    G = 0 ;
    step1() ;
    min_step7 = 50.0 ;
    double min_G = 0.0  ; 
    for(G = 0 ; G <= max_G ; G++){
        if( G % 100 == 0){output_point_fields(G) ;}
        step2() ; step3() ; step4() ; //...
        //step5:判斷收斂步驟
        if(min_step7 > step5(G)){ //每一步都判斷
                min_step7 = step5(G) ;
                min_G = G ; 
                if( min_step7 < 2.0 ){output_point_fields(min_G) ;}
            }
        if(G%10== 0){
            cout << "Itteration step = " << G << endl ;
            cout << "(G = "<< G << " )step: (1)Ghia = " << step5(G) << "_(2)Mass = " << step5star (G)<< endl ;
            cout << "currnt minerror(G =" << min_G << ")= "<< min_step7 << endl ;
        }
        if(G>0 && step5(G) < 0.1){
            CE_physical = true ;
            cout << "Haved converged at step G = " << G << endl ;
            cout << "(G = "<< G << " )step: (1)Ghia = " << step5(G) << "_(2)Mass = " << step5star (G)<< endl ;
            cout << "currnt minerror(G =" << G << ")= "<< min_step7 << endl ;
            output_point_fields(G) ;
            break ;
        }
        step6() ; 
    }
    if(!CE_physical){
        cout << "Not converged after max_G = " << max_G << " steps." << endl ;
        cout << "(G = "<< G << " )step: (1)Ghia = " << step5(G) << "_(2)Mass = " << step5star (G)<< endl ;
        output_point_fields(max_G) ;
    }
    cout << "Grid Size : " << Nx << "x" << Ny << "Completed." << endl ;
    return 0 ;
}