// Solve the Navier-Stokes equations using the Finite Volume Method with the SIMPLE algorithm
// Add deferred correction and use alternating iteration
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <vector> 
#include <omp.h>  // Add after #include <vector>
#include "0.P_correction_eq_TDMA.h" 
using namespace std ;
// Parameter settings
// Parameters to be used
// Pressure field has Nx cells in the x-direction
// x-velocity field has Nx+1 cells in the x-direction
// y-velocity field has Ny+1 cells in the y-direction
int G ;  // G is the current calculation step starting from 0, different from iteration step
const int  max_G = 500000 ; // Maximum calculation steps
double min_step7 ; 
// x-velocity field: i in [0,Nx], j in [0,Ny-1]
// y-velocity field: i in [0,Nx-1], j in [0,Ny]
// Pressure field: i in [0,Nx-1], j in [0,Ny-1]
double u_star[Nx+1][Ny] ; // Predicted x-velocity field
double u_star_old[Nx+1][Ny] ;
double u_star_old2[Nx+1][Ny] ;
double v_star[Nx][Ny+1] ; // Predicted y-velocity field
double v_star_old[Nx][Ny+1] ;
double v_star_old2[Nx][Ny+1] ;
double p_star[Nx][Ny] ; // Predicted pressure field
double u_prime[Nx+1][Ny] ; // Corrected x-velocity field
double v_prime[Nx][Ny+1] ; // Corrected y-velocity field
double p_prime[Nx][Ny] ; // Corrected pressure field
double u[Nx+1][Ny] ; // Correct x-velocity field
// Create matrix for old correct x-velocity field
double v[Nx][Ny+1] ; // Correct y-velocity field
// Create matrix for old correct y-velocity field
double p[Nx][Ny] ; // Correct pressure field

const double  Re = 5000.0 ;
// Boundary conditions for x-velocity field
const double u_left  = 0.0 ; // No-penetration wall condition
const double u_right = 0.0 ; // No-penetration wall condition
const double u_bottom = 0.0 ; // No-slip condition
const double u_top = 1.0 ; // No-slip condition
// Boundary conditions for y-velocity field
const double v_left  = 0.0 ; // No-slip condition
const double v_right = 0.0 ; // No-slip condition
const double v_bottom = 0.0 ; // No-penetration wall condition
const double v_top = 0.0 ; // No-penetration wall condition
const double omega_P = 0.09 ; // According to Patankar's approach, the relaxation factor for pressure: alpha_P = 0.8
double omega_u  = 0.01 ; 
double omega_v  = 0.01 ;
double max_error ;
const float CE_tolerance = 1e-10; // Convergence criterion: continuity equation for predicted values
bool CE_physical ;
//
// Diffusion coefficients of the discretized equations (shared by x and y components)
const double D_e = (1.0/Re)*(dy/dx);
const double D_w = (1.0/Re)*(dy/dx);
const double D_n = (1.0/Re)*(dy/dx);
const double D_s = (1.0/Re)*(dy/dx);
double r_w_u[Nx+1][Ny] ; // Interface decision coefficient for x-velocity field
double r_e_u[Nx+1][Ny] ;
double r_s_u[Nx+1][Ny] ;
double r_n_u[Nx+1][Ny] ;
double r_w_v[Nx][Ny+1] ; // Interface decision coefficient for y-velocity field
double r_e_v[Nx][Ny+1] ;
double r_s_v[Nx][Ny+1] ;
double r_n_v[Nx][Ny+1] ;
// Convection coefficients of the general discretized form of the x-momentum integral equation // The factor alpha is determined by F
double F_w_u[Nx+1][Ny] ;
double alpha_w_u[Nx+1][Ny] ;
double F_e_u[Nx+1][Ny] ;
double alpha_e_u[Nx+1][Ny] ;
double F_s_u[Nx+1][Ny] ;
double alpha_s_u[Nx+1][Ny] ;
double F_n_u[Nx+1][Ny] ;
double alpha_n_u[Nx+1][Ny] ;
// Linear combination coefficients of the general discretized form of the x-momentum integral equation
double a_P[Nx+1][Ny] ;
double a_WW[Nx+1][Ny] ;
double a_W[Nx+1][Ny] ;
double a_E[Nx+1][Ny] ;
double a_EE[Nx+1][Ny] ;
double a_SS[Nx+1][Ny] ;
double a_S[Nx+1][Ny] ;
double a_N[Nx+1][Ny] ;
double a_NN[Nx+1][Ny] ; 
double S1[Nx+1][Ny] ; // Inhomogeneous term of the x-momentum equation
// Convection coefficients of the general discretized form of the y-momentum integral equation // The factor alpha is determined by F
double F_w_v[Nx][Ny+1] ;
double alpha_w_v[Nx][Ny+1] ;
double F_e_v[Nx][Ny+1] ;
double alpha_e_v[Nx][Ny+1] ;
double F_s_v[Nx][Ny+1] ;
double alpha_s_v[Nx][Ny+1] ;
double F_n_v[Nx][Ny+1] ;
double alpha_n_v[Nx][Ny+1] ;
// Linear combination coefficients of the general discretized form of the y-momentum integral equation
double b_P[Nx][Ny+1] ;
double b_WW[Nx][Ny+1] ;
double b_W[Nx][Ny+1] ;
double b_E[Nx][Ny+1] ;
double b_EE[Nx][Ny+1] ;
double b_SS[Nx][Ny+1] ;
double b_S[Nx][Ny+1] ;
double b_N[Nx][Ny+1] ;
double b_NN[Nx][Ny+1] ;
double S2[Nx][Ny+1] ; // Inhomogeneous term of the y-momentum equation
// Remember to update these matrices immediately after the predicted velocity field is updated


//(1)
// CDs upstream and downstream iteration coefficients
double D(double x){
	return 0.5  ;
}
double U(double x){
    return 0.5 ;
}
double UU(double x){
	return 0.0 ;

}

/*//(2)
// QUICK upstream and downstream iteration coefficients
double D(double x){
	return 3.0/8.0     ;
}
double U(double x){
    return 6.0/8.0 ;
}
double UU(double x){
	return -1.0/8.0 ;

}*/


/*//(3)TVD-MUSCL upstream and downstream iteration coefficients
double D(double x) {
    if (x <= 0) return 0.0;
    else if (x > 0 && x < (1.0/3.0)) return 0.0;
    else if (x >= (1.0/3.0) && x < 3.0) return 0.25;
    else return 1.0;
}

double U(double x) {
    if (x <= 0) return 1.0;
    else if (x > 0 && x <= (1.0/3.0)) return 2.0;
    else if (x > (1.0/3.0) && x <= 3.0) return 1.0;
    else return 0.0;
}

double UU(double x) {
    if (x <= 0) return 0.0;
    else if (x > 0 && x < (1.0/3.0)) return -1.0;
    else if (x >= (1.0/3.0) && x < 3.0) return -0.25;
    else return 0.0;
}*/


// SIMPLE algorithm:
// step1: Initialize predicted x-velocity field, predicted y-velocity field, and predicted pressure field
void initialize(){
    #pragma omp parallel  // *** Need to add this line! ***
    {
    // Initialize predicted x-velocity field
     #pragma omp parallel for collapse(2) 
    for(int i = 1 ; i <= Nx+1 ; i++){
        // Domain of predicted x-velocity field: i in [0,Nx], j in [0,Ny-1]
        for(int J= 1 ; J <= Ny ; J++){
            u_star[i-1][J-1] = 0.0 ;
            if(i == 1) u_star[i-1][J-1] = u_left ; // No-penetration wall condition
            if(i == Nx + 1) u_star[i-1][J-1] = u_right ; // No-penetration wall condition
        }
    }
    // Initialize old correct x-velocity field -> assigned from initial predicted values, later updated from true correct values
     #pragma omp parallel for collapse(2) 
    for(int i = 1 ; i <= Nx+1 ; i++){
        for(int J = 1 ; J <= Ny ; J++){
            u_star_old[i-1][J-1] = u_star[i-1][J-1] ;
            u_star_old2[i-1][J-1] = u_star_old[i-1][J-1];
        }
    }
    // Initialize predicted y-velocity field
     #pragma omp parallel for collapse(2) 
    for(int I = 1 ; I <= Nx ; I++){
        // Domain of predicted y-velocity field: I in [0,Nx-1], j in [0,Ny]
        for(int j = 1 ; j <= Ny+1 ; j++){
            v_star[I-1][j-1] = 0.0 ;
            if(j == 1)v_star[I-1][j-1] = v_bottom ; // No-penetration wall condition
            if(j == Ny + 1)v_star[I-1][j-1] = v_top ; // No-penetration wall condition
        }
    }
    // Initialize old correct y-velocity field -> assigned from initial predicted values, later updated from true correct values
     #pragma omp parallel for collapse(2) 
    for(int I = 1 ; I <= Nx ; I++){
        for(int j = 1 ; j <= Ny+1 ; j++){
            v_star_old[I-1][j-1] = v_star[I-1][j-1] ;
            v_star_old2[I-1][j-1] =v_star_old[I-1][j-1] ;  
        }
    }
    // Initialize predicted pressure field
     #pragma omp parallel for collapse(2) 
    for(int I = 1 ; I <= Nx ; I++){
        for(int J = 1 ; J<= Ny ; J++){
            p_star[I-1][J-1] =  0.0  ; // Boundary values follow a Neumann condition enforced by the iterative scheme
            p_prime[I-1][J-1] = 0.0 ;
        }
    }
    // Matrices required for the predictor equations: interface ratios -> convective fluxes -> linear-combination coefficients
    // Initialization (and subsequent updates) must follow that order: interface ratios, convective terms, then linear coefficients
    // Initialize interface ratio matrices
    // x-velocity interface ratios (north and south faces require special handling)
     #pragma omp parallel for collapse(2) 
    for(int i = 1 ; i <= Nx+1 ; i++){
        for(int J = 1 ; J <= Ny ; J++){
            // Cells that fall outside the upwind stencil simply set ratios to zero
            if(i == 1 || i == Nx+1){
                r_w_u[i-1][J-1] = 0.0 ;
                r_e_u[i-1][J-1] = 0.0 ;
                r_n_u[i-1][J-1] = 0.0 ;
                r_s_u[i-1][J-1] = 0.0 ;
            }else {
                // Handle each boundary separately (west, east, south, north)
                // West boundary
                if(i == 2){
                    r_w_u[i-1][J-1] = -1.0 ; // West boundary uses upwind interpolation
                }else{
                    r_w_u[i-1][J-1] = 0.0 ; // Default initialization
                }
                // East boundary
                if(i == Nx){
                    r_e_u[i-1][J-1] = -1.0 ; // East boundary uses upwind interpolation
                }else{
                    r_e_u[i-1][J-1] = 0.0 ;// Default initialization
                }
                // South boundary
                if(J == 1){
                    r_s_u[i-1][J-1] = 0.0 ; // Bottom-most cells do not use a south-face ratio
                }else if(J == 2){
                    r_s_u[i-1][J-1] = -1.0 ;// Upwind scheme with safeguards against division by zero
                }else if(J == Ny){
                    r_s_u[i-1][J-1] = -1.0 ;// Upwind scheme with safeguards against division by zero
                }else{
                    r_s_u[i-1][J-1] = 0.0 ;// Default initialization
                }
                // North boundary
                if(J == 1){
                    r_n_u[i-1][J-1] = -1.0 ; // Upwind scheme with safeguards against division by zero
                }else if (J == Ny-1){
                    r_n_u[i-1][J-1] = -1.0 ;// Upwind scheme with safeguards against division by zero
                }else if(J == Ny){
                    r_n_u[i-1][J-1] = 0.0 ; // Top-most cells do not use a north-face ratio
                }else{
                    r_n_u[i-1][J-1] = 0.0 ;// Default initialization
                }
            }
        }    
    }// Remember to refresh interface ratios immediately after updating predictors
    // y-velocity interface ratios (east/west faces require special handling)
     #pragma omp parallel for collapse(2) 
    for(int I = 1 ; I <= Nx ; I++){
        for(int j = 1 ; j <= Ny+1 ; j++){
            // Cells that fall outside the upwind stencil simply set ratios to zero
            if(j == 1 || j == Ny+1){
                r_w_v[I-1][j-1] = 0.0 ;
                r_e_v[I-1][j-1] = 0.0 ;
                r_n_v[I-1][j-1] = 0.0 ;
                r_s_v[I-1][j-1] = 0.0 ;
            }else {
                // Handle each boundary in turn (south, north, west, east)
                // South boundary
                if(j == 2){
                    r_s_v[I-1][j-1] = -1.0 ; // Bottom boundary uses upwind interpolation
                }else{
                    r_s_v[I-1][j-1] = 0.0 ; // Default initialization
                }
                // North boundary
                if(j == Ny){
                    r_n_v[I-1][j-1] = -1.0 ; // Top boundary uses upwind interpolation
                }else{
                    r_n_v[I-1][j-1] = 0.0 ; // Default initialization
                }
                // West boundary
                if(I == 1){
                    r_w_v[I-1][j-1] = 0.0 ; // West-most cells do not require a west-face ratio
                }else if (I == 2){
                    r_w_v[I-1][j-1] = -1.0 ;// Upwind scheme with safeguards against division by zero
                }else if(I == Nx){
                    r_w_v[I-1][j-1] = -1.0 ;// Upwind scheme with safeguards against division by zero
                }else{
                    r_w_v[I-1][j-1] = 0.0 ;// Default initialization
                }
                // East boundary
                if(I == 1){
                    r_e_v[I-1][j-1] = -1.0 ; // Upwind scheme with safeguards against division by zero
                }else if (I == Nx-1){
                    r_e_v[I-1][j-1] = -1.0 ; // Upwind scheme with safeguards against division by zero
                }else if(I == Nx){
                    r_e_v[I-1][j-1] = 0.0 ; // East-most cells do not require an east-face ratio
                }else {
                    r_e_v[I-1][j-1] = 0.0 ; // Default initialization
                }
            }
        }
    }
    // Initialize the convective-flux coefficient matrices
    // x-momentum convective coefficients in the discrete form
     #pragma omp parallel for collapse(2) 
    for(int i = 1 ; i <= Nx+1 ; i++){
        for(int J = 1 ; J <= Ny ; J++){
            if(i == 1 || i == Nx+1) {
                F_w_u[i-1][J-1] = 0.0 ;
                F_e_u[i-1][J-1] = 0.0 ;
                F_n_u[i-1][J-1] = 0.0 ;
                F_s_u[i-1][J-1] = 0.0 ;
            }else{
                F_e_u[i-1][J-1] = 0.5*dy*(u_star[i-1][J-1]+u_star[i][J-1]); // Control-volume centers align with storage locations
                F_w_u[i-1][J-1] = 0.5*dy*(u_star[i-1][J-1]+u_star[i-2][J-1]);
                F_n_u[i-1][J-1] = 0.5*dx*(v_star[i-1][J]+v_star[i-2][J]);
                F_s_u[i-1][J-1] = 0.5*dx*(v_star[i-1][J-1]+v_star[i-2][J-1]);
            }
        }
    }// alpha_f_u flags are stored at the centers of control volumes
     #pragma omp parallel for collapse(2) 
    for(int i = 1 ; i <= Nx+1 ; i++){
        for(int J = 1 ; J <= Ny ; J++){
            if(F_w_u[i-1][J-1] > 0) alpha_w_u[i-1][J-1] = 1.0 ;
            else alpha_w_u[i-1][J-1] = 0.0 ;
            if(F_e_u[i-1][J-1] > 0) alpha_e_u[i-1][J-1] = 1.0 ;
            else alpha_e_u[i-1][J-1] = 0.0 ;
            if(F_n_u[i-1][J-1] > 0) alpha_n_u[i-1][J-1] = 1.0 ;
            else alpha_n_u[i-1][J-1] = 0.0 ;
            if(F_s_u[i-1][J-1] > 0) alpha_s_u[i-1][J-1] = 1.0 ;
            else alpha_s_u[i-1][J-1] = 0.0 ;
        }
    }
    // y-momentum convective coefficients in the discrete form
     #pragma omp parallel for collapse(2) 
    for(int I = 1 ; I <= Nx ; I++){
        for(int j = 1 ; j <= Ny+1 ; j++){
            if(j == 1 || j == Ny+1) {
                F_w_v[I-1][j-1] = 0.0 ;
                F_e_v[I-1][j-1] = 0.0 ;
                F_n_v[I-1][j-1] = 0.0 ;
                F_s_v[I-1][j-1] = 0.0 ;
            }else{
                F_w_v[I-1][j-1] = 0.5*dy*(u_star[I-1][j-1]+u_star[I-1][j-2]); 
                F_e_v[I-1][j-1] = 0.5*dy*(u_star[I][j-1]+u_star[I][j-2]);
                F_n_v[I-1][j-1] = 0.5*dx*(v_star[I-1][j-1]+v_star[I-1][j]);// Control-volume centers align with storage locations
                F_s_v[I-1][j-1] = 0.5*dx*(v_star[I-1][j-1]+v_star[I-1][j-2]);
            }
        }
    }// alpha_f_v flags are stored at the centers of control volumes
     #pragma omp parallel for collapse(2) 
    for(int I = 1 ; I <= Nx ; I++){
        for(int j = 1 ; j <= Ny+1 ; j++){
            if(F_w_v[I-1][j-1] > 0) alpha_w_v[I-1][j-1] = 1.0 ;
            else alpha_w_v[I-1][j-1] = 0.0 ;
            if(F_e_v[I-1][j-1] > 0) alpha_e_v[I-1][j-1] = 1.0 ;
            else alpha_e_v[I-1][j-1] = 0.0 ;
            if(F_n_v[I-1][j-1] > 0) alpha_n_v[I-1][j-1] = 1.0 ;
            else alpha_n_v[I-1][j-1] = 0.0 ;
            if(F_s_v[I-1][j-1] > 0) alpha_s_v[I-1][j-1] = 1.0 ;
            else alpha_s_v[I-1][j-1] = 0.0 ;
        }
    }
}
}

//step2: update the predicted velocity fields (u and v components)
//step2: update the predicted velocity fields (u and v components)
void step2(int m){
    //update the linear combination coefficients for the predicted u-momentum equation
    #pragma omp parallel for collapse(2)
    for(int i = 1 ; i <= Nx+1 ; i++){
        for(int J = 1 ; J <= Ny ; J++){
            //reset all neighbor coefficients to zero before accumulation
            a_P[i-1][J-1] = 0.0 ;
            a_WW[i-1][J-1] = 0.0 ;
            a_W[i-1][J-1] = 0.0 ;
            a_E[i-1][J-1] = 0.0 ;
            a_EE[i-1][J-1] = 0.0 ;
            a_SS[i-1][J-1] = 0.0 ;
            a_S[i-1][J-1] = 0.0 ;
            a_N[i-1][J-1] = 0.0 ;
            a_NN[i-1][J-1] = 0.0 ;
            S1[i-1][J-1] = 0.0 ;
            if(i == 1){
                a_P[i-1][J-1] += 1.0 ;
                S1[i-1][J-1] += u_left ;//no-slip wall; west boundary has zero pressure gradient
                u_star[i-1][J-1] = u_left ;
            }else if(i == Nx+1){
                a_P[i-1][J-1] += 1.0 ;
                S1[i-1][J-1] += u_right ;//no-slip wall; east boundary has zero pressure gradient
                u_star[i-1][J-1] = u_right ;
            }// Special cases: west-column second row, east-column second row, and first/second rows near bottom/top boundaries
            // Second column adjacent to the west boundary
            else{ // For all other interior cells, handle west/east/south/north faces individually
                // West boundary
                if(i == 2){// Second column next to the west boundary uses an upwind scheme on the west face
                    a_W[i-1][J-1] += alpha_w_u[i-1][J-1]*(-F_w_u[i-1][J-1])-D_w;// Subtract west-face diffusion when applying boundary treatment
                    a_P[i-1][J-1] += (1-alpha_w_u[i-1][J-1])*(-F_w_u[i-1][J-1])+D_w;
                    S1[i-1][J-1] += -(-p_star[i-2][J-1])*dy/rho ;
                }else{// General west-face upwind stencil involves W and P
                    //a_WW[i-1][J-1] += (alpha_w_u[i-1][J-1]*UU(r_w_u[i-1][J-1]))*(-F_w_u[i-1][J-1]) ;
                    a_W[i-1][J-1] += (alpha_w_u[i-1][J-1]*1)*(-F_w_u[i-1][J-1])-D_w ; // Subtract west-face diffusion when applying boundary treatment
                    //a_E[i-1][J-1] += ((1-alpha_w_u[i-1][J-1])*UU(r_w_u[i-1][J-1]))*(-F_w_u[i-1][J-1]);
                    a_P[i-1][J-1] += ((1-alpha_w_u[i-1][J-1])*1)*(-F_w_u[i-1][J-1])+D_w ;
                    S1[i-1][J-1] += -(-p_star[i-2][J-1])*dy/rho
                    -u_star_old[i-3][J-1]*(alpha_w_u[i-1][J-1]*UU(r_w_u[i-1][J-1]))*(-F_w_u[i-1][J-1])
                    -u_star_old[i-2][J-1]*((1-alpha_w_u[i-1][J-1])*D(r_w_u[i-1][J-1])+alpha_w_u[i-1][J-1]*(U(r_w_u[i-1][J-1])-1))*(-F_w_u[i-1][J-1])
                    -u_star_old[i][J-1]*((1-alpha_w_u[i-1][J-1])*UU(r_w_u[i-1][J-1]))*(-F_w_u[i-1][J-1])
                    -u_star_old[i-1][J-1]*((1-alpha_w_u[i-1][J-1])*(U(r_w_u[i-1][J-1])-1)+alpha_w_u[i-1][J-1]*D(r_w_u[i-1][J-1]))*(-F_w_u[i-1][J-1]) 
                    ; 
                }// Deferred correction complete
                // East boundary
                if(i == Nx){// Second column next to the east boundary uses an upwind scheme on the east face
                    a_E[i-1][J-1] += (1-alpha_e_u[i-1][J-1])*(F_e_u[i-1][J-1])-D_e ;// Subtract east-face diffusion when applying boundary treatment
                    a_P[i-1][J-1] += alpha_e_u[i-1][J-1]*(F_e_u[i-1][J-1])+D_e; 
                    S1[i-1][J-1] += -(p_star[i-1][J-1])*dy/rho ;
                }else{// General east-face upwind stencil involves E and P
                    //a_EE[i-1][J-1] += (1-alpha_e_u[i-1][J-1])*UU(r_e_u[i-1][J-1])*(F_e_u[i-1][J-1]) ;
                    a_E[i-1][J-1] += ((1-alpha_e_u[i-1][J-1])*1)*(F_e_u[i-1][J-1])-D_e ;// Subtract east-face diffusion when applying boundary treatment
                    //a_W[i-1][J-1] += (alpha_e_u[i-1][J-1])*UU(r_e_u[i-1][J-1])*(F_e_u[i-1][J-1]) ;
                    a_P[i-1][J-1] += (alpha_e_u[i-1][J-1]*1)*(F_e_u[i-1][J-1])+D_e;   
                    S1[i-1][J-1] += -(p_star[i-1][J-1])*dy/rho 
                    -u_star_old[i+1][J-1]*(1-alpha_e_u[i-1][J-1])*UU(r_e_u[i-1][J-1])*(F_e_u[i-1][J-1])//EE
                    -u_star_old[i][J-1]*(alpha_e_u[i-1][J-1]*D(r_e_u[i-1][J-1])+(1-alpha_e_u[i-1][J-1])*(U(r_e_u[i-1][J-1])-1))*(F_e_u[i-1][J-1])//E
                    -u_star_old[i-2][J-1]*(alpha_e_u[i-1][J-1])*UU(r_e_u[i-1][J-1])*(F_e_u[i-1][J-1])//W
                    -u_star_old[i-1][J-1]*((1-alpha_e_u[i-1][J-1])*D(r_e_u[i-1][J-1])+alpha_e_u[i-1][J-1]*(U(r_e_u[i-1][J-1])-1))*(F_e_u[i-1][J-1])//p
                    ;   
                }// Deferred correction complete
                // South boundary: treat three different rows
                if(J == 1){// First row adjacent to the bottom boundary
                    a_P[i-1][J-1] += 2*D_s ;
                    S1[i-1][J-1] += -(u_bottom)*(-F_s_u[i-1][J-1])+2*D_s*(u_bottom) ;// No-slip condition; south boundary pressure equals the cell value via the Neumann condition
                }else if(J == 2){// Second row just above the bottom boundary uses extrapolation instead of the SS node; south-face upwind couples S and P
                    //a_SS[i-1][J-1] += (alpha_s_u[i-1][J-1]*UU(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1]) ;
                    a_S[i-1][J-1] += ((alpha_s_u[i-1][J-1])*1)*(-F_s_u[i-1][J-1])-D_s ;// Subtract south-face diffusion when applying boundary treatment
                    //a_N[i-1][J-1] += ((1-alpha_s_u[i-1][J-1])*UU(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1]) ;
                    a_P[i-1][J-1] += ((1-alpha_s_u[i-1][J-1])*1)*(-F_s_u[i-1][J-1])+D_s ; 
                    S1[i-1][J-1] += 
                    -(2*u_bottom - u_star_old[i-1][J-2])*(alpha_s_u[i-1][J-1]*UU(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1]) //SS
                    -u_star_old[i-1][J-2]*((1-alpha_s_u[i-1][J-1])*D(r_s_u[i-1][J-1])+(alpha_s_u[i-1][J-1])*(U(r_s_u[i-1][J-1])-1))*(-F_s_u[i-1][J-1])//S
                    -u_star_old[i-1][J]*((1-alpha_s_u[i-1][J-1])*UU(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1])//N
                    -u_star_old[i-1][J-1]*((1-alpha_s_u[i-1][J-1])*(U(r_s_u[i-1][J-1])-1)+alpha_s_u[i-1][J-1]*D(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1])//P
                    ;// Deferred correction complete
                }else if(J == Ny){
                    // Extrapolate the north node back to the current cell
                     //a_SS[i-1][J-1] += (alpha_s_u[i-1][J-1]*UU(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1]) ;
                    a_S[i-1][J-1] += ((alpha_s_u[i-1][J-1])*1)*(-F_s_u[i-1][J-1])-D_s ;// Subtract south-face diffusion when applying boundary treatment
                    //a_N[i-1][J-1] += ((1-alpha_s_u[i-1][J-1])*UU(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1]) ;
                    a_P[i-1][J-1] += ((1-alpha_s_u[i-1][J-1])*1)*(-F_s_u[i-1][J-1])+D_s ; 
                    S1[i-1][J-1] += 
                    -u_star_old[i-1][J-3]*(alpha_s_u[i-1][J-1]*UU(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1]) //SS
                    -u_star_old[i-1][J-2]*((1-alpha_s_u[i-1][J-1])*D(r_s_u[i-1][J-1])+(alpha_s_u[i-1][J-1])*(U(r_s_u[i-1][J-1])-1))*(-F_s_u[i-1][J-1])//S
                    -(2*u_top -u_star_old[i-1][J-1])*((1-alpha_s_u[i-1][J-1])*UU(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1])//N
                    -u_star_old[i-1][J-1]*((1-alpha_s_u[i-1][J-1])*(U(r_s_u[i-1][J-1])-1)+alpha_s_u[i-1][J-1]*D(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1])//P
                    ;// Deferred correction complete
                }else{// General south-face upwind stencil involves S and P
                    //a_SS[i-1][J-1] += (alpha_s_u[i-1][J-1]*UU(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1]) ;
                    a_S[i-1][J-1] += ((alpha_s_u[i-1][J-1])*1)*(-F_s_u[i-1][J-1])-D_s ;// Subtract south-face diffusion when applying boundary treatment
                    //a_N[i-1][J-1] += ((1-alpha_s_u[i-1][J-1])*UU(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1]) ;
                    a_P[i-1][J-1] += ((1-alpha_s_u[i-1][J-1])*1)*(-F_s_u[i-1][J-1])+D_s ; 
                    S1[i-1][J-1] += 
                    -u_star_old[i-1][J-3]*(alpha_s_u[i-1][J-1]*UU(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1]) //SS
                    -u_star_old[i-1][J-2]*((1-alpha_s_u[i-1][J-1])*D(r_s_u[i-1][J-1])+(alpha_s_u[i-1][J-1])*(U(r_s_u[i-1][J-1])-1))*(-F_s_u[i-1][J-1])//S
                    -u_star_old[i-1][J]*((1-alpha_s_u[i-1][J-1])*UU(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1])//N
                    -u_star_old[i-1][J-1]*((1-alpha_s_u[i-1][J-1])*(U(r_s_u[i-1][J-1])-1)+alpha_s_u[i-1][J-1]*D(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1])//P
                    ;// Deferred correction complete
                }
                // North boundary requires three special cases: J == 1, J == Ny, and J == Ny-1
                if(J == Ny){
                    a_P[i-1][J-1] += 2*D_n ;
                    S1[i-1][J-1] += -(u_top)*F_n_u[i-1][Ny-1] + 2*D_n*(u_top);// No-slip condition; north boundary pressure equals the cell value via the Neumann condition
                }else if(J == 1){
                     //a_NN[i-1][J-1] += (1-alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1]) ;
                    a_N[i-1][J-1] += ((1-alpha_n_u[i-1][J-1])*1)*(F_n_u[i-1][J-1]) - D_n ;// Subtract north-face diffusion when applying boundary treatment
                    //a_S[i-1][J-1] += (alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1]) ;
                    a_P[i-1][J-1] += (alpha_n_u[i-1][J-1]*1)*(F_n_u[i-1][J-1])+D_n ;
                    S1[i-1][J-1] += 
                    -u_star_old[i-1][J+1]*(1-alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1]) 
                    -u_star_old[i-1][J]*(alpha_n_u[i-1][J-1]*D(r_n_u[i-1][J-1])+(1-alpha_n_u[i-1][J-1])*(U(r_n_u[i-1][J-1])-1))*(F_n_u[i-1][J-1]) 
                    -(2*u_bottom - u_star_old[i-1][J-1])*(alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1])
                    -u_star_old[i-1][J-1]*((1-alpha_n_u[i-1][J-1])*D(r_n_u[i-1][J-1])+alpha_n_u[i-1][J-1]*(U(r_n_u[i-1][J-1])-1))*(F_n_u[i-1][J-1])
                    ;
                }else if(J == Ny-1){
                    //a_NN[i-1][J-1] += (1-alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1]) ;
                    a_N[i-1][J-1] += ((1-alpha_n_u[i-1][J-1])*1)*(F_n_u[i-1][J-1]) - D_n ;// Subtract north-face diffusion when applying boundary treatment
                    //a_S[i-1][J-1] += (alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1]) ;
                    a_P[i-1][J-1] += (alpha_n_u[i-1][J-1]*1)*(F_n_u[i-1][J-1])+D_n ;
                    S1[i-1][J-1] += 
                    -(2*u_top - u_star_old[i-1][J])*(1-alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1]) 
                    -u_star_old[i-1][J]*(alpha_n_u[i-1][J-1]*D(r_n_u[i-1][J-1])+(1-alpha_n_u[i-1][J-1])*(U(r_n_u[i-1][J-1])-1))*(F_n_u[i-1][J-1]) 
                    -u_star_old[i-1][J-2]*(alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1])
                    -u_star_old[i-1][J-1]*((1-alpha_n_u[i-1][J-1])*D(r_n_u[i-1][J-1])+alpha_n_u[i-1][J-1]*(U(r_n_u[i-1][J-1])-1))*(F_n_u[i-1][J-1])
                    ;
                }else{
                    //a_NN[i-1][J-1] += (1-alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1]) ;
                    a_N[i-1][J-1] += ((1-alpha_n_u[i-1][J-1])*1)*(F_n_u[i-1][J-1]) - D_n ;// Subtract north-face diffusion when applying boundary treatment
                    //a_S[i-1][J-1] += (alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1]) ;
                    a_P[i-1][J-1] += (alpha_n_u[i-1][J-1]*1)*(F_n_u[i-1][J-1])+D_n ;
                    S1[i-1][J-1] += 
                    -u_star_old[i-1][J+1]*(1-alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1]) 
                    -u_star_old[i-1][J]*(alpha_n_u[i-1][J-1]*D(r_n_u[i-1][J-1])+(1-alpha_n_u[i-1][J-1])*(U(r_n_u[i-1][J-1])-1))*(F_n_u[i-1][J-1]) 
                    -u_star_old[i-1][J-2]*(alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1])
                          -u_star_old[i-1][J-1]*((1-alpha_n_u[i-1][J-1])*D(r_n_u[i-1][J-1])+alpha_n_u[i-1][J-1]*(U(r_n_u[i-1][J-1])-1))*(F_n_u[i-1][J-1])
                    ;
                 }
            }
        }
    }
    //update the linear combination coefficients for the predicted v-momentum equation
    #pragma omp parallel default(shared)
    {
        #pragma omp for schedule(static) //average cutting workking for threads 
        for(int J = 1 ; J <= Ny ; J++){//define the index from 0~Nx-2
        vector <double> beta (Nx-1) ;//left hand side elemets 
        vector <double> D (Nx-1) ; //primary elements
        vector <double> alpha (Nx-1) ; //right and side elements
        vector <double> x (Nx-1) ; //solutions 
        vector <double> S (Nx-1) ; //Source Vector 
        for(int k = 2 ; k <= Nx ; k++){
            beta[k-2] = a_W[k-1][J-1] ; 
            D[k-2] = a_P[k-1][J-1] ; 
            alpha[k-2] = a_E[k-1][J-1] ; 
            x[k-2] = u_star[k-1][J-1] ; 
            if(J == 1){S[k-2] = S1[k-1][J-1] -u_star_old[k-1][J]*a_N[k-1][J-1];}
            else if(J == Ny){ S[k-2] = S1[k-1][J-1]-u_star_old[k-1][J-2]*a_S[k-1][J-1] ; }
            else {S[k-2] = S1[k-1][J-1]-u_star_old[k-1][J-2]*a_S[k-1][J-1]-u_star_old[k-1][J]*a_N[k-1][J-1]; ; }
        }
        TDMA(D , beta , alpha , x , S , Nx-1 , u_star[1-1][J-1] , u_star[Nx+1-1][J-1]) ; 
            for(int i = 2 ; i <= Nx ; i++){
                u_star[i-1][J-1] = x[i-2] ;
            }
        }
        #pragma omp for collapse(2) schedule(static) //average cutting workking for threads 
        for(int i = 1 ; i <= Nx+1 ; i++){
            for(int J = 1 ; J <= Ny ; J++){
                u_star_old2[i-1][J-1] = u_star[i-1][J-1] ; 
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
            x[k-1] = u_star[i-1][k-1] ; 
            S[k-1] = S1[i-1][k-1]-u_star_old2[i-2][k-1]*a_W[i-1][k-1]-u_star_old2[i][k-1]*a_E[i-1][k-1];
            }
            TDMA_FVM(D , beta , alpha , x , S , Ny) ; 
            for(int J = 1 ; J <= Ny ; J++){
                u_star[i-1][J-1] = x[J-1] ;
            }
        }
        #pragma omp for collapse(2) schedule(static)
        for(int i = 2 ; i <= Nx ; i++){
            for(int J = 1 ; J <= Ny ; J++){
                u_star[i-1][J-1] = omega_u*u_star[i-1][J-1] + (1-omega_u)*u_star_old[i-1][J-1] ;
            }
        }
    }
    #pragma omp parallel for collapse(2)
    for(int I = 1 ; I <= Nx ; I++){
        for(int j = 1 ; j <= Ny+1 ; j++){
            //reset all neighbor coefficients to zero before accumulation
            b_P[I-1][j-1] = 0.0 ;
            b_WW[I-1][j-1] = 0.0 ;
            b_W[I-1][j-1] = 0.0 ;
            b_E[I-1][j-1] = 0.0 ;
            b_EE[I-1][j-1] = 0.0 ;
            b_SS[I-1][j-1] = 0.0 ;
            b_S[I-1][j-1] = 0.0 ;
            b_N[I-1][j-1] = 0.0 ;
            b_NN[I-1][j-1] = 0.0 ;
            S2[I-1][j-1] = 0.0 ;
            if(j == 1){
                b_P[I-1][j-1] += 1.0 ;
                S2[I-1][j-1] += v_bottom ;//no-slip wall; bottom boundary has zero pressure gradient
                v_star[I-1][j-1] = v_bottom ;
            }else if(j == Ny+1){
                b_P[I-1][j-1] += 1.0 ;
                S2[I-1][j-1] += v_top ;//no-slip wall; top boundary has zero pressure gradient
                v_star[I-1][j-1] = v_top ;
            }//special handling needed for first/second columns near west and east walls as well as second rows near top and bottom
            //second row next to the bottom boundary
            else{ //all remaining cells are processed by checking west, east, south, and north faces individually
                //south boundary
                if(j == 2){//second row above the bottom boundary uses upwind on the south face
                    b_S[I-1][j-1] += alpha_s_v[I-1][j-1]*(-F_s_v[I-1][j-1])-D_s;//subtract south-face diffusion when enforcing the boundary
                    b_P[I-1][j-1] += (1-alpha_s_v[I-1][j-1])*(-F_s_v[I-1][j-1])+D_s;
                    S2[I-1][j-1] += -(-p_star[I-1][j-2])*dx/rho ;
                }else{//interior south faces use the TVD-MUSCL stencil that involves the south node and the current cell
                    //b_SS[I-1][j-1] += (alpha_s_v[I-1][j-1]*UU(r_s_v[I-1][j-1]))*(-F_s_v[I-1][j-1]) ;
                    b_S[I-1][j-1] += (alpha_s_v[I-1][j-1]*1)*(-F_s_v[I-1][j-1])-D_s ; //subtract south-face diffusion when enforcing the boundary
                    //b_N[I-1][j-1] += ((1-alpha_s_v[I-1][j-1])*UU(r_s_v[I-1][j-1]))*(-F_s_v[I-1][j-1]) ;
                    b_P[I-1][j-1] += ((1-alpha_s_v[I-1][j-1])*1)*(-F_s_v[I-1][j-1])+D_s ;
                    S2[I-1][j-1] += -(-p_star[I-1][j-2])*dx/rho 
                    -v_star_old[I-1][j-3]*(alpha_s_v[I-1][j-1]*UU(r_s_v[I-1][j-1]))*(-F_s_v[I-1][j-1]) //SS
                    -v_star_old[I-1][j-2]*((1-alpha_s_v[I-1][j-1])*D(r_s_v[I-1][j-1])+alpha_s_v[I-1][j-1]*(U(r_s_v[I-1][j-1])-1))*(-F_s_v[I-1][j-1]) //S
                    -v_star_old[I-1][j]*((1-alpha_s_v[I-1][j-1])*UU(r_s_v[I-1][j-1]))*(-F_s_v[I-1][j-1]) //N
                    -v_star_old[I-1][j-1]*((1-alpha_s_v[I-1][j-1])*(U(r_s_v[I-1][j-1])-1)+alpha_s_v[I-1][j-1]*D(r_s_v[I-1][j-1]))*(-F_s_v[I-1][j-1]) //P
                    ;//deferred correction completed
                }
                //north boundary
                if(j == Ny){//second row below the top wall uses upwind on the north face
                    b_N[I-1][j-1] += (1-alpha_n_v[I-1][j-1])*(F_n_v[I-1][j-1])-D_n ;//subtract north-face diffusion when enforcing the boundary
                    b_P[I-1][j-1] += alpha_n_v[I-1][j-1]*(F_n_v[I-1][j-1])+D_n; 
                    S2[I-1][j-1] += -(p_star[I-1][j-1])*dx/rho ;
                }else{//other north faces employ the TVD-MUSCL stencil coupling the current and north cells
                    //b_NN[I-1][j-1] += (1-alpha_n_v[I-1][j-1])*UU(r_n_v[I-1][j-1])*(F_n_v[I-1][j-1]) ;
                    b_N[I-1][j-1] += ((1-alpha_n_v[I-1][j-1])*1)*(F_n_v[I-1][j-1]) - D_n ;//subtract north-face diffusion when enforcing the boundary
                    //b_S[I-1][j-1] += (alpha_n_v[I-1][j-1])*UU(r_n_v[I-1][j-1])*(F_n_v[I-1][j-1]) ;
                    b_P[I-1][j-1] += (alpha_n_v[I-1][j-1]*1)*(F_n_v[I-1][j-1]) + D_n ;
                    S2[I-1][j-1] += -(p_star[I-1][j-1])*dx/rho
                    -v_star_old[I-1][j+1]*(1-alpha_n_v[I-1][j-1])*UU(r_n_v[I-1][j-1])*(F_n_v[I-1][j-1]) //NN 
                    -v_star_old[I-1][j]*(alpha_n_v[I-1][j-1]*D(r_n_v[I-1][j-1])+(1-alpha_n_v[I-1][j-1])*(U(r_n_v[I-1][j-1])-1))*(F_n_v[I-1][j-1]) //N
                    -v_star_old[I-1][j-2]*(alpha_n_v[I-1][j-1])*UU(r_n_v[I-1][j-1])*(F_n_v[I-1][j-1]) //S
                    -v_star_old[I-1][j-1]*((1-alpha_n_v[I-1][j-1])*D(r_n_v[I-1][j-1])+alpha_n_v[I-1][j-1]*(U(r_n_v[I-1][j-1])-1))*(F_n_v[I-1][j-1]) //P
                    ; //deferred correction completed
                }
                //west boundary (three cases: I = 1, 2, Nx)
                if(I == 2){//second column near the west wall uses extrapolation to replace the WW value
                    //b_WW[I-1][j-1] += (alpha_w_v[I-1][j-1]*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1]) ;
                    b_W[I-1][j-1] += ((alpha_w_v[I-1][j-1])*1)*(-F_w_v[I-1][j-1])-D_w ; //subtract west-face diffusion when enforcing the boundary
                    //b_E[I-1][j-1] += ((1-alpha_w_v[I-1][j-1])*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1]);
                    b_P[I-1][j-1] += ((1-alpha_w_v[I-1][j-1])*1)*(-F_w_v[I-1][j-1])+D_w ;
                    S2[I-1][j-1] += 
                    -(2*v_left - v_star_old[I-2][j-1])*(alpha_w_v[I-1][j-1]*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1]) //WW
                    -v_star_old[I-2][j-1]*((1-alpha_w_v[I-1][j-1])*D(r_w_v[I-1][j-1])+(alpha_w_v[I-1][j-1])*(U(r_w_v[I-1][j-1])-1))*(-F_w_v[I-1][j-1]) //W
                    -v_star_old[I][j-1]*((1-alpha_w_v[I-1][j-1])*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1]) //E
                    -v_star_old[I-1][j-1]*((1-alpha_w_v[I-1][j-1])*(U(r_w_v[I-1][j-1])-1)+(alpha_w_v[I-1][j-1])*D(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1]) //P
                    ;//deferred correction completed
                }else if(I == 1){//first column next to the west wall uses upwind on the west face
                    b_P[I-1][j-1] +=  2*D_w;//subtract west-face diffusion when enforcing the boundary
                    S2[I-1][j-1] += -(-F_w_v[I-1][j-1])*v_left+2*D_w*v_left;
                 }else if(I == Nx){//E->P
                     //b_WW[I-1][j-1] += (alpha_w_v[I-1][j-1]*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1]) ;
                    b_W[I-1][j-1] += ((alpha_w_v[I-1][j-1])*1)*(-F_w_v[I-1][j-1])-D_w ; //subtract west-face diffusion when enforcing the boundary
                    //b_E[I-1][j-1] += ((1-alpha_w_v[I-1][j-1])*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1]);
                    b_P[I-1][j-1] += ((1-alpha_w_v[I-1][j-1])*1)*(-F_w_v[I-1][j-1])+D_w ;
                    S2[I-1][j-1] += 
                    -v_star_old[I-3][j-1]*(alpha_w_v[I-1][j-1]*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1]) //WW
                    -v_star_old[I-2][j-1]*((1-alpha_w_v[I-1][j-1])*D(r_w_v[I-1][j-1])+(alpha_w_v[I-1][j-1])*(U(r_w_v[I-1][j-1])-1))*(-F_w_v[I-1][j-1]) //W
                    -(2*v_right - v_star_old[I-1][j-1])*((1-alpha_w_v[I-1][j-1])*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1]) //E
                    -v_star_old[I-1][j-1]*((1-alpha_w_v[I-1][j-1])*(U(r_w_v[I-1][j-1])-1)+(alpha_w_v[I-1][j-1])*D(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1]) //P
                    ;//deferred correction completed
                }else{//remaining west faces use the TVD-MUSCL stencil
                    //b_WW[I-1][j-1] += (alpha_w_v[I-1][j-1]*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1]) ;
                    b_W[I-1][j-1] += ((alpha_w_v[I-1][j-1])*1)*(-F_w_v[I-1][j-1])-D_w ; //subtract west-face diffusion when enforcing the boundary
                    //b_E[I-1][j-1] += ((1-alpha_w_v[I-1][j-1])*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1]);
                    b_P[I-1][j-1] += ((1-alpha_w_v[I-1][j-1])*1)*(-F_w_v[I-1][j-1])+D_w ;
                    S2[I-1][j-1] += 
                    -v_star_old[I-3][j-1]*(alpha_w_v[I-1][j-1]*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1]) //WW
                    -v_star_old[I-2][j-1]*((1-alpha_w_v[I-1][j-1])*D(r_w_v[I-1][j-1])+(alpha_w_v[I-1][j-1])*(U(r_w_v[I-1][j-1])-1))*(-F_w_v[I-1][j-1]) //W
                    -v_star_old[I][j-1]*((1-alpha_w_v[I-1][j-1])*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1]) //E
                    -v_star_old[I-1][j-1]*((1-alpha_w_v[I-1][j-1])*(U(r_w_v[I-1][j-1])-1)+(alpha_w_v[I-1][j-1])*D(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1]) //P
                    ;//deferred correction completed
                }
                //east boundary (three cases: I = 1, Nx-1, Nx)
                if(I == 1){//second column near the east wall replaces the EE value by extrapolation (W->P)
                    //b_EE[I-1][j-1] += (1-alpha_e_v[I-1][j-1])*UU(r_e_v[I-1][j-1])*(F_e_v[I-1][j-1]) ;
                    b_E[I-1][j-1] += ((1-alpha_e_v[I-1][j-1])*1)*(F_e_v[I-1][j-1]) - D_e ;//subtract east-face diffusion when enforcing the boundary
                    //b_W[I-1][j-1] += (alpha_e_v[I-1][j-1])*UU(r_e_v[I-1][j-1])*(F_e_v[I-1][j-1]) ;
                    b_P[I-1][j-1] += (alpha_e_v[I-1][j-1]*1)*(F_e_v[I-1][j-1]) + D_e ;
                    S2[I-1][j-1] += 
                    -v_star_old[I+1][j-1]*(1-alpha_e_v[I-1][j-1])*UU(r_e_v[I-1][j-1])*(F_e_v[I-1][j-1])//EE
                    -v_star_old[I][j-1]*(alpha_e_v[I-1][j-1]*D(r_e_v[I-1][j-1])+(1-alpha_e_v[I-1][j-1])*(U(r_e_v[I-1][j-1])-1))*(F_e_v[I-1][j-1])//E
                    -(2*v_left - v_star_old[I-1][j-1])*(alpha_e_v[I-1][j-1])*UU(r_e_v[I-1][j-1])*(F_e_v[I-1][j-1]) //W
                    -v_star_old[I-1][j-1]*((1-alpha_e_v[I-1][j-1])*D(r_e_v[I-1][j-1])+alpha_e_v[I-1][j-1]*(U(r_e_v[I-1][j-1])-1))*(F_e_v[I-1][j-1])//P
                    ;//deferred correction completed
                }else if(I == Nx){//second column near the east wall uses upwind on the east face
                    b_P[I-1][j-1] += 2*D_e ;//subtract east-face diffusion when enforcing the boundary
                    S2[I-1][j-1] += -(F_e_v[I-1][j-1])*v_right +2*D_e*v_right;
                }else if(I == Nx-1){//EE node does not exist; extrapolate to the E location
                    //b_EE[I-1][j-1] += (1-alpha_e_v[I-1][j-1])*UU(r_e_v[I-1][j-1])*(F_e_v[I-1][j-1]) ;
                    b_E[I-1][j-1] += ((1-alpha_e_v[I-1][j-1])*1)*(F_e_v[I-1][j-1]) - D_e ;//subtract east-face diffusion when enforcing the boundary
                    //b_W[I-1][j-1] += (alpha_e_v[I-1][j-1])*UU(r_e_v[I-1][j-1])*(F_e_v[I-1][j-1]) ;
                    b_P[I-1][j-1] += (alpha_e_v[I-1][j-1]*1)*(F_e_v[I-1][j-1]) + D_e ;
                    S2[I-1][j-1] += 
                    -(2*v_right-v_star_old[I][j-1])*(1-alpha_e_v[I-1][j-1])*UU(r_e_v[I-1][j-1])*(F_e_v[I-1][j-1])//EE
                    -v_star_old[I][j-1]*(alpha_e_v[I-1][j-1]*D(r_e_v[I-1][j-1])+(1-alpha_e_v[I-1][j-1])*(U(r_e_v[I-1][j-1])-1))*(F_e_v[I-1][j-1])//E
                    -v_star_old[I-2][j-1]*(alpha_e_v[I-1][j-1])*UU(r_e_v[I-1][j-1])*(F_e_v[I-1][j-1]) //W
                    -v_star_old[I-1][j-1]*((1-alpha_e_v[I-1][j-1])*D(r_e_v[I-1][j-1])+alpha_e_v[I-1][j-1]*(U(r_e_v[I-1][j-1])-1))*(F_e_v[I-1][j-1])//P
                    ;//deferred correction completed
                }else{//remaining east faces use the TVD-MUSCL stencil
                    //b_EE[I-1][j-1] += (1-alpha_e_v[I-1][j-1])*UU(r_e_v[I-1][j-1])*(F_e_v[I-1][j-1]) ;
                    b_E[I-1][j-1] += ((1-alpha_e_v[I-1][j-1])*1)*(F_e_v[I-1][j-1]) - D_e ;//subtract east-face diffusion when enforcing the boundary
                    //b_W[I-1][j-1] += (alpha_e_v[I-1][j-1])*UU(r_e_v[I-1][j-1])*(F_e_v[I-1][j-1]) ;
                    b_P[I-1][j-1] += (alpha_e_v[I-1][j-1]*1)*(F_e_v[I-1][j-1]) + D_e ;
                    S2[I-1][j-1] += 
                    -v_star_old[I+1][j-1]*(1-alpha_e_v[I-1][j-1])*UU(r_e_v[I-1][j-1])*(F_e_v[I-1][j-1])//EE
                    -v_star_old[I][j-1]*(alpha_e_v[I-1][j-1]*D(r_e_v[I-1][j-1])+(1-alpha_e_v[I-1][j-1])*(U(r_e_v[I-1][j-1])-1))*(F_e_v[I-1][j-1])//E
                    -v_star_old[I-2][j-1]*(alpha_e_v[I-1][j-1])*UU(r_e_v[I-1][j-1])*(F_e_v[I-1][j-1]) //W
                    -v_star_old[I-1][j-1]*((1-alpha_e_v[I-1][j-1])*D(r_e_v[I-1][j-1])+alpha_e_v[I-1][j-1]*(U(r_e_v[I-1][j-1])-1))*(F_e_v[I-1][j-1])//P
                    ;//deferred correction completed
                }
            }
        }
    }
    //compute the y-velocity predicted value 
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
            x[k-1] = v_star[k-1][j-1] ; 
            S[k-1] = S2[k-1][j-1] -v_star_old[k-1][j-2]*b_S[k-1][j-1]-v_star_old[k-1][j]*b_N[k-1][j-1] ;   
            } 
            TDMA_FVM(D , beta , alpha , x , S , Nx ) ;
            for(int I = 1 ; I <= Nx ; I++){
                v_star[I-1][j-1] = x[I-1] ; 
            }
        }
        #pragma omp for collapse(2) schedule(static) //average cutting workking for threads 
        for(int I = 1 ; I <= Nx ; I++){
            for(int j = 1 ; j <= Ny+1 ; j++){
                v_star_old2[I-1][j-1] = v_star[I-1][j-1] ; 
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
            for(int k = 2 ; k <= Ny ; k++){
                beta[k-2] = b_S[I-1][k-1] ; 
                D[k-2] = b_P[I-1][k-1] ; 
                alpha[k-2] = b_N[I-1][k-1] ; 
                x[k-2] = v_star[I-1][k-1] ; 
                if(I == 1){ S[k-2] = S2[I-1][k-1] - v_star_old2[I][k-1]*b_E[I-1][k-1] ; }
                else if(I == Nx){ S[k-2] = S2[I-1][k-1] - v_star_old2[I-2][k-1]*b_W[I-1][k-1] ; }
                else{ S[k-2] = S2[I-1][k-1] - v_star_old2[I-2][k-1]*b_W[I-1][k-1] - v_star_old2[I][k-1]*b_E[I-1][k-1] ;}
                 
            }
            TDMA(D , beta , alpha , x , S , Ny-1 , 0.0 , 0.0 ) ;
            for(int j = 2 ; j <= Ny ; j++){
                v_star[I-1][j-1] = x[j-2]  ; 
            }
        }
        #pragma omp for collapse(2) schedule(static) //average cutting workking for threads 
        for(int I = 1 ; I <= Nx ; I++){
            for(int j = 1 ; j <= Ny+1 ; j++){
                v_star[I-1][j-1] = omega_v*v_star[I-1][j-1]+(1-omega_v)*v_star_old[I-1][j-1] ; 
            }
        } 
    }

}
//u_star and v_star are now available at step N+1/2
//Function: mass flux immediately after the first velocity predictor update
double d(int I, int J , int m){
    //I,J reference the cell-center coordinates; m is the pressure iteration counter
    double u_e = u_star[I+1-1][J-1] ;//velocity predictor is evaluated at N+1/2
    double u_w = u_star[I-1][J-1] ;//velocity predictor is evaluated at N+1/2
    double v_n = v_star[I-1][J+1-1] ;//velocity predictor is evaluated at N+1/2
    double v_s = v_star[I-1][J-1] ;//velocity predictor is evaluated at N+1/2
    return rho*(u_e - u_w)*dy + rho*(v_n - v_s)*dx ;
} //used after the first update of the velocity predictor


//step4: update the corrected pressure field
void step4(int m){
    #pragma omp parallel for collapse(2) 
    for(int I = 1 ; I <= Nx ; I++){ //pressure values (corrected, predicted, and corrected) reside at cell centers
        for(int J = 1 ; J <= Ny ; J++){
           p[I-1][J-1] = p_star[I-1][J-1] + omega_P*p_prime[I-1][J-1] ;
        }//p_star is evaluated at step N; p_prime at step N+1
    }//p_prime is evaluated at step N+1
}

//step5: update velocity corrections; boundary velocities are fixed by the boundary conditions and should not iterate
void step5(int m){
    #pragma omp parallel for collapse(2) 
    //update correction for the x-velocity component
    for(int i = 1 ; i <= Nx+1 ; i++){
        for(int J = 1 ; J <= Ny ; J++){
            if(i == 1 || i == Nx+1){
                u_prime[i-1][J-1] = 0.0 ;//boundary velocities already satisfy the physical constraints
            }else{//pressure correction: current cell minus the west neighbor
                u_prime[i-1][J-1] = -(p_prime[i-1][J-1] - p_prime[i-2][J-1])*dy/(rho*a_P[i-1][J-1]) ; ;
            }

        }
    }
    //update correction for the y-velocity component
    #pragma omp parallel for collapse(2) 
    for(int I = 1 ; I <= Nx ; I++){
        for(int j = 1 ; j <= Ny+1 ; j++){
            if(j == 1 || j == Ny+1){
                v_prime[I-1][j-1] = 0.0 ;//boundary velocities already satisfy the physical constraints
            }else{//pressure correction: current cell minus the south neighbor
                v_prime[I-1][j-1] = -(p_prime[I-1][j-1] - p_prime[I-1][j-2])*dx/(rho*b_P[I-1][j-1]) ;
            }
        }
    }
}

double step5star(int m){
    double sum = 0.0 ;
    #pragma omp parallel for collapse(2) reduction(+:sum)  // accumulate the residuals in parallel
    for(int I = 1 ; I <= Nx ; I++){
        for(int J = 1 ; J <= Ny ; J++){
            sum = sum + fabs(d(I,J,m)) ;
        }
    }
    return sum ;
}

//step6: update the corrected velocity fields and refresh the stored old values
void step6(int m){//boundary velocities remain fixed by the boundary conditions
    //update corrected x-velocity
    #pragma omp parallel for collapse(2) 
    for(int i = 1 ; i <= Nx+1 ; i++){
        for(int J = 1 ; J <= Ny ; J++){
            if(i == 1 ){
                u[i-1][J-1] = 0.0 ;//boundary velocity imposed by the wall
            }else if (i == Nx+1){
                u[i-1][J-1] = 0.0 ;//boundary velocity imposed by the wall
            }else{
                u[i-1][J-1] = u_prime[i-1][J-1]+u_star[i-1][J-1];//u_star at N+1/2, u_prime at N+1, u_old at N, final u at N+1
            }
        }
    }
    //update corrected y-velocity
    #pragma omp parallel for collapse(2) 
    for(int I = 1 ; I <= Nx ; I++){
        for(int j = 1 ; j <= Ny+1 ; j++){
            if(j == 1 ){
                v[I-1][j-1] = 0.0 ;//impermeable wall
            }else if (j == Ny+1){
                v[I-1][j-1] = 0.0 ;//impermeable wall
            }else{
                v[I-1][j-1] = v_prime[I-1][j-1]+v_star[I-1][j-1];//v_star at N+1/2, v_prime at N+1, v_old at N, final v at N+1
            }
        }
    }
}
double step7(int m){
    //81

    double sum = 0.0 ;
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
    sum += fabs(v[79][40]-(0.0)) ;
    /*
    //113
    double sum = 0.0 ;
    sum += fabs(u[56][0]-0) ; 
    sum += fabs(u[56][6]-(-0.41165)) ;
    sum += fabs(u[56][6]-(-0.42901)) ; 
    sum += fabs(u[56][7]-(-0.43643)) ;  
    sum += fabs(u[56][11]-(-0.40435)) ;  
    sum += fabs(u[56][19]-(-0.3305)) ;  
    sum += fabs(u[56][31]-(-0.22855)) ;
    sum += fabs(u[56][50]-(-0.07404)) ;
    sum += fabs(u[56][55]-(-0.03039)) ;
    sum += fabs(u[56][69]-(0.08183)) ; 
    sum += fabs(u[56][82]-(0.20087)) ;  
    sum += fabs(u[56][95]-(0.33556)) ;
    sum += fabs(u[56][106]-(0.46036)) ;
    sum += fabs(u[56][107]-(0.45992)) ;
    sum += fabs(u[56][108]-(0.4612)) ;
    sum += fabs(u[56][109]-(0.48223)) ;
    sum += fabs(u[56][111]-(1)) ;
    
    sum += fabs(v[0][56]-0) ;
    sum += fabs(v[6][56]-(0.42447)) ; 
    sum += fabs(v[7][56]-(0.43329)) ;
    sum += fabs(v[8][56]-(0.43648)) ; 
    sum += fabs(v[10][56]-(0.42951)) ;
    sum += fabs(v[17][56]-(0.35368)) ; 
    sum += fabs(v[25][56]-(0.28066)) ;
    sum += fabs(v[26][56]-(0.2728)) ;
    sum += fabs(v[55][56]-(0.00945)) ;
    sum += fabs(v[90][56]-(-0.30018)) ; 
    sum += fabs(v[96][56]-(-0.36214)) ;
    sum += fabs(v[101][56]-(-0.41442)) ; 
    sum += fabs(v[105][56]-(-0.52876)) ;
    sum += fabs(v[106][56]-(-0.55408)) ; 
    sum += fabs(v[107][56]-(-0.55069)) ;
    sum += fabs(v[108][56]-(-0.49774)) ;
    sum += fabs(v[111][56]-(0.0)) ;*/
    /*
    //161
    double sum = 0.0 ;
        sum += fabs(u[80][0]-0) ; 
    sum += fabs(u[80][8]-(-0.18109)) ;
    sum += fabs(u[80][9]-(-0.20196)) ; 
    sum += fabs(u[80][10]-(-0.22220)) ;  
    sum += fabs(u[80][15]-(-0.29730)) ;  
    sum += fabs(u[80][27]-(-0.38289)) ;  
    sum += fabs(u[80][44]-(-0.27805)) ;
    sum += fabs(u[80][71]-(-0.10648)) ;
    sum += fabs(u[80][79]-(-0.06080)) ;
    sum += fabs(u[80][98]-(0.05702)) ; 
    sum += fabs(u[80][117]-(0.18719)) ;  
    sum += fabs(u[80][135]-(0.33304)) ;
    sum += fabs(u[80][151]-(0.46604)) ;
    sum += fabs(u[80][153]-(0.51117)) ;
    sum += fabs(u[80][154]-(0.57492)) ;
    sum += fabs(u[80][155]-(0.65928)) ;
    sum += fabs(u[80][159]-(1)) ;
    
    sum += fabs(v[0][80]-0) ;
    sum += fabs(v[9][80]-(0.27485)) ; 
    sum += fabs(v[10][80]-(0.29012)) ;
    sum += fabs(v[12][80]-(0.30353)) ; 
    sum += fabs(v[14][80]-(0.32627)) ;
    sum += fabs(v[24][80]-(0.37095)) ; 
    sum += fabs(v[35][80]-(0.33075)) ;
    sum += fabs(v[37][80]-(0.32235)) ;
    sum += fabs(v[80][80]-(0.02526)) ;
    sum += fabs(v[129][80]-(-0.31966)) ; 
    sum += fabs(v[137][80]-(-0.42655)) ;
    sum += fabs(v[145][80]-(-0.51550)) ; 
    sum += fabs(v[151][80]-(-0.39188)) ;
    sum += fabs(v[152][80]-(-0.33714)) ; 
    sum += fabs(v[154][80]-(-0.27669)) ;
    sum += fabs(v[155][80]-(-0.21388)) ;
    sum += fabs(v[159][80]-(0.0)) ;*/
    /*//227
    double sum = 0.0 ;
    sum += fabs(u[113][0]-0) ; 
    sum += fabs(u[113][12]-(-0.41165)) ;
    sum += fabs(u[113][14]-(-0.42901)) ; 
    sum += fabs(u[113][15]-(-0.43643)) ;  
    sum += fabs(u[113][22]-(-0.40435)) ;  
    sum += fabs(u[113][38]-(-0.3305)) ;  
    sum += fabs(u[113][63]-(-0.22855)) ;
    sum += fabs(u[113][102]-(-0.07404)) ;
    sum += fabs(u[113][113]-(-0.03039)) ;
    sum += fabs(u[113][139]-(0.08183)) ; 
    sum += fabs(u[113][165]-(0.20087)) ;  
    sum += fabs(u[113][192]-(0.33556)) ;
    sum += fabs(u[113][215]-(0.46036)) ;
    sum += fabs(u[113][217]-(0.45992)) ;
    sum += fabs(u[113][218]-(0.4612)) ;
    sum += fabs(u[113][220]-(0.48223)) ;
    sum += fabs(u[113][225]-(1)) ;
    
    sum += fabs(v[0][113]-0) ;
    sum += fabs(v[14][113]-(0.42447)) ; 
    sum += fabs(v[15][113]-(0.43329)) ;
    sum += fabs(v[17][113]-(0.43648)) ; 
    sum += fabs(v[21][113]-(0.42951)) ;
    sum += fabs(v[35][113]-(0.35368)) ; 
    sum += fabs(v[51][113]-(0.28066)) ;
    sum += fabs(v[52][113]-(0.2728)) ;
    sum += fabs(v[113][113]-(0.00945)) ;
    sum += fabs(v[181][113]-(-0.30018)) ; 
    sum += fabs(v[194][113]-(-0.36214)) ;
    sum += fabs(v[204][113]-(-0.41442)) ; 
    sum += fabs(v[214][113]-(-0.52876)) ;
    sum += fabs(v[215][113]-(-0.55408)) ; 
    sum += fabs(v[217][113]-(-0.55069)) ;
    sum += fabs(v[218][113]-(-0.49774)) ;
    sum += fabs(v[225][113]-(0.0)) ;*/
    return sum ;
}
    
//step8: refresh the pressure predictor at cell centers
void step8(int m){
    #pragma omp parallel for collapse(2)  // use OpenMP to refresh predictions in parallel
    for(int I = 1 ; I <= Nx ; I++){
        for(int J = 1 ; J <= Ny ; J++){
            p_star[I-1][J-1] = p[I-1][J-1] ; //both p and p_star are now evaluated at step N+1
        }
    }
}    
    
//step9: update velocity predictors, interface ratios, and convective coefficients
void step9(int m){
    //update x-velocity predictor; boundary values come directly from the imposed conditions
    #pragma omp parallel for collapse(2)
    for(int i = 1 ; i <= Nx+1 ; i++){
        for(int J = 1 ; J <= Ny ; J++){
            if(i == 1){
                u_star[i-1][J-1] = u_left;
            }else if (i == Nx+1){
                u_star[i-1][J-1] = u_right ;
            }else{
                //internal cells, including top and bottom faces, inherit the corrected velocity
                u_star[i-1][J-1] = u[i-1][J-1]; //u evaluated at N+1; new u_star also at N+1
                u_star_old[i-1][J-1] = u_star[i-1][J-1]; //store the predictor immediately
            }
        }
    }
    //update y-velocity predictor; boundary values come directly from the imposed conditions
    #pragma omp parallel for collapse(2)
    for(int I = 1 ; I <= Nx ; I++){
        for(int j = 1 ; j <= Ny+1 ; j++){
            if(j == 1){
                v_star[I-1][j-1] = v_bottom ;
            }else if (j == Ny+1){
                v_star[I-1][j-1] = v_top ;
            }else{
                //internal cells, including left and right faces, inherit the corrected velocity
                v_star[I-1][j-1] = v[I-1][j-1]; //v evaluated at N+1; new v_star also at N+1
                v_star_old[I-1][j-1] = v_star[I-1][j-1]; //store the predictor immediately
            }
        }
    }
    //update convective coefficients
    //x-momentum convection terms in the general finite-volume discretization
    #pragma omp parallel for collapse(2)
    for(int i = 1 ; i <= Nx+1 ; i++){
        for(int J = 1 ; J <= Ny ; J++){
            if(i == 1 || i == Nx+1) {
                F_w_u[i-1][J-1] = 0.0 ;
                F_e_u[i-1][J-1] = 0.0 ;
                F_n_u[i-1][J-1] = 0.0 ;
                F_s_u[i-1][J-1] = 0.0 ;
            }else{//recompute convective fluxes
                F_e_u[i-1][J-1] = 0.5*dy*(u_star[i-1][J-1]+u_star[i][J-1]); //cell centers align with the storage location
                F_w_u[i-1][J-1] = 0.5*dy*(u_star[i-1][J-1]+u_star[i-2][J-1]);
                F_n_u[i-1][J-1] = 0.5*dx*(v_star[i-1][J]+v_star[i-2][J]);
                F_s_u[i-1][J-1] = 0.5*dx*(v_star[i-1][J-1]+v_star[i-2][J-1]);
            }
        }
    }//alpha_f_u switch; indices denote control-volume centers
    #pragma omp parallel for collapse(2)
    for(int i = 1 ; i <= Nx+1 ; i++){
        for(int J = 1 ; J <= Ny ; J++){
            if(F_w_u[i-1][J-1] > 0) alpha_w_u[i-1][J-1] = 1.0 ;
            else alpha_w_u[i-1][J-1] = 0.0 ;
            if(F_e_u[i-1][J-1] > 0) alpha_e_u[i-1][J-1] = 1.0 ;
            else alpha_e_u[i-1][J-1] = 0.0 ;
            if(F_n_u[i-1][J-1] > 0) alpha_n_u[i-1][J-1] = 1.0 ;
            else alpha_n_u[i-1][J-1] = 0.0 ;
            if(F_s_u[i-1][J-1] > 0) alpha_s_u[i-1][J-1] = 1.0 ;
            else alpha_s_u[i-1][J-1] = 0.0 ;
        }
    }
    //y-momentum convection terms in the general finite-volume discretization
    #pragma omp parallel for collapse(2)
    for(int I = 1 ; I <= Nx ; I++){
        for(int j = 1 ; j <= Ny+1 ; j++){
            if(j == 1 || j == Ny+1) {
                F_w_v[I-1][j-1] = 0.0 ;
                F_e_v[I-1][j-1] = 0.0 ;
                F_n_v[I-1][j-1] = 0.0 ;
                F_s_v[I-1][j-1] = 0.0 ;
            }else{
                F_w_v[I-1][j-1] = 0.5*dy*(u_star[I-1][j-1]+u_star[I-1][j-2]); 
                F_e_v[I-1][j-1] = 0.5*dy*(u_star[I][j-1]+u_star[I][j-2]);
                F_n_v[I-1][j-1] = 0.5*dx*(v_star[I-1][j-1]+v_star[I-1][j]);//cell centers align with the storage location
                F_s_v[I-1][j-1] = 0.5*dx*(v_star[I-1][j-1]+v_star[I-1][j-2]);
            }
        }
    }//alpha_f_v switch; indices denote control-volume centers
    #pragma omp parallel for collapse(2)
    for(int I = 1 ; I <= Nx ; I++){
        for(int j = 1 ; j <= Ny+1 ; j++){
            if(F_w_v[I-1][j-1] > 0) alpha_w_v[I-1][j-1] = 1.0 ;
            else alpha_w_v[I-1][j-1] = 0.0 ;
            if(F_e_v[I-1][j-1] > 0) alpha_e_v[I-1][j-1] = 1.0 ;
            else alpha_e_v[I-1][j-1] = 0.0 ;
            if(F_n_v[I-1][j-1] > 0) alpha_n_v[I-1][j-1] = 1.0 ;
            else alpha_n_v[I-1][j-1] = 0.0 ;
            if(F_s_v[I-1][j-1] > 0) alpha_s_v[I-1][j-1] = 1.0 ;
            else alpha_s_v[I-1][j-1] = 0.0 ;
        }
    }
    //update interface ratio indicators
    #pragma omp parallel for collapse(2)
    for(int i = 1 ; i <= Nx+1 ; i++){
        for(int J = 1 ; J <= Ny ; J++){
            //cells outside the stencil set the ratio to zero
            if(i == 1 || i == Nx+1){
                r_w_u[i-1][J-1] = 0.0 ;
                r_e_u[i-1][J-1] = 0.0 ;
                r_n_u[i-1][J-1] = 0.0 ;
                r_s_u[i-1][J-1] = 0.0 ;
            }else {
                //treat each boundary case-by-case
                //west boundary
                if(i == 2){
                    r_w_u[i-1][J-1] = -1.0 ; //west boundary uses upwind weighting
                }else{
                    double u_W = u_star[i-2][J-1]; 
                    double u_P = u_star[i-1][J-1];
                    double u_WW = u_star[i-3][J-1];
                    double u_E = u_star[i][J-1] ;
                    double D = u_W + alpha_w_u[i-1][J-1]*(u_P - u_W) ;
                    double U = u_P + alpha_w_u[i-1][J-1]*(u_W - u_P) ;
                    double UU = u_E + alpha_w_u[i-1][J-1]*(u_WW - u_E) ;
                    r_w_u[i-1][J-1] = (U-UU)/(D-U+1e-30) ; //initial MUSCL ratio
                }
                //east boundary
                if(i == Nx){
                    r_e_u[i-1][J-1] = -1.0 ; //east boundary uses upwind weighting
                }else{
                    double u_P = u_star[i-1][J-1] ;
                    double u_E = u_star[i][J-1] ;
                    double u_W = u_star[i-2][J-1] ;
                    double u_EE = u_star[i+1][J-1] ;
                    double D = u_P + alpha_e_u[i-1][J-1]*(u_E - u_P) ;
                    double U = u_E + alpha_e_u[i-1][J-1]*(u_P - u_E) ;
                    double UU = u_EE + alpha_e_u[i-1][J-1]*(u_W - u_EE) ;
                    r_e_u[i-1][J-1] = (U-UU)/(D-U+1e-30) ; //initial MUSCL ratio
                }
                //south boundary
                if(J == 1){
                    r_s_u[i-1][J-1] = 0.0 ; //bottom-most control volumes do not need the south ratio
                }else if(J == 2){
                    double u_S = u_star[i-1][J-2] ;
                    double u_P = u_star[i-1][J-1] ;
                    double u_SS = 2*u_bottom-u_S ; //mirror south neighbor around the boundary
                    double u_N = u_star[i-1][J] ;
                    double D = u_S + alpha_s_u[i-1][J-1]*(u_P - u_S) ;
                    double U = u_P + alpha_s_u[i-1][J-1]*(u_S - u_P) ;
                    double UU = u_N + alpha_s_u[i-1][J-1]*(u_SS - u_N) ;
                    r_s_u[i-1][J-1] =  (U-UU)/(D-U+1e-30) ; //initial MUSCL ratio
                }else if(J == Ny){
                    double u_S = u_star[i-1][J-2] ;
                    double u_P = u_star[i-1][J-1] ;
                    double u_SS = u_star[i-1][J-3] ;
                    double u_N = 2*u_top - u_P ; //mirror the point across the top boundary
                    double D = u_S + alpha_s_u[i-1][J-1]*(u_P - u_S) ;
                    double U = u_P + alpha_s_u[i-1][J-1]*(u_S - u_P) ;
                    double UU = u_N + alpha_s_u[i-1][J-1]*(u_SS - u_N) ;
                    r_s_u[i-1][J-1] =  (U-UU)/(D-U+1e-30) ; //initial MUSCL ratio
                }else{
                    double u_S = u_star[i-1][J-2] ;
                    double u_P = u_star[i-1][J-1] ;
                    double u_SS = u_star[i-1][J-3] ;
                    double u_N = u_star[i-1][J] ;
                    double D = u_S + alpha_s_u[i-1][J-1]*(u_P - u_S) ;
                    double U = u_P + alpha_s_u[i-1][J-1]*(u_S - u_P) ;
                    double UU = u_N + alpha_s_u[i-1][J-1]*(u_SS - u_N) ;
                    r_s_u[i-1][J-1] =  (U-UU)/(D-U+1e-30) ; //initial MUSCL ratio
                }
                //north boundary
                if(J == 1){
                    double u_P = u_star[i-1][J-1] ;
                    double u_N = u_star[i-1][J] ;
                    double u_S = 2*u_bottom-u_P ; //mirror the point across the bottom boundary
                    double u_NN = u_star[i-1][J+1] ;
                    double D = u_P + alpha_n_u[i-1][J-1]*(u_N - u_P) ;
                    double U = u_N + alpha_n_u[i-1][J-1]*(u_P - u_N) ;
                    double UU = u_NN + alpha_n_u[i-1][J-1]*(u_S - u_NN) ;
                    r_n_u[i-1][J-1] =  (U-UU)/(D-U+1e-30) ; //initial MUSCL ratio
                }else if (J == Ny-1){
                    double u_P = u_star[i-1][J-1] ;
                    double u_N = u_star[i-1][J] ;
                    double u_S = u_star[i-1][J-2] ;
                    double u_NN = 2*u_top - u_N ; //mirror the north neighbor across the boundary
                    double D = u_P + alpha_n_u[i-1][J-1]*(u_N - u_P) ;
                    double U = u_N + alpha_n_u[i-1][J-1]*(u_P - u_N) ;
                    double UU = u_NN + alpha_n_u[i-1][J-1]*(u_S - u_NN) ;
                    r_n_u[i-1][J-1] =  (U-UU)/(D-U+1e-30) ; //initial MUSCL ratio
                }else if(J == Ny){
                    r_n_u[i-1][J-1] = 0.0 ; //top-most control volumes do not need the north ratio
                }else{
                    double u_P = u_star[i-1][J-1] ;
                    double u_N = u_star[i-1][J] ;
                    double u_S = u_star[i-1][J-2] ;
                    double u_NN = u_star[i-1][J+1] ;
                    double D = u_P + alpha_n_u[i-1][J-1]*(u_N - u_P) ;
                    double U = u_N + alpha_n_u[i-1][J-1]*(u_P - u_N) ;
                    double UU = u_NN + alpha_n_u[i-1][J-1]*(u_S - u_NN) ;
                    r_n_u[i-1][J-1] =  (U-UU)/(D-U+1e-30) ; //initial MUSCL ratio
                }
            }
        }    
    }//after updating predictors, refresh the interface ratios immediately
    //interface ratios for the y-velocity (east and west faces require special handling)
    #pragma omp parallel for collapse(2)
    for(int I = 1 ; I <= Nx ; I++){
        for(int j = 1 ; j <= Ny+1 ; j++){
            //cells outside the stencil set the ratio to zero
            if(j == 1 || j == Ny+1){
                r_w_v[I-1][j-1] = 0.0 ;
                r_e_v[I-1][j-1] = 0.0 ;
                r_n_v[I-1][j-1] = 0.0 ;
                r_s_v[I-1][j-1] = 0.0 ;
            }else {
                //process each boundary configuration explicitly
                //south boundary
                if(j == 2){
                    r_s_v[I-1][j-1] = -1.0 ; //south boundary uses upwind weighting
                }else{
                    double v_S = v_star[I-1][j-2] ;
                    double v_P = v_star[I-1][j-1] ;
                    double v_SS = v_star[I-1][j-3] ;
                    double v_N = v_star[I-1][j] ;
                    double D = v_S + alpha_s_v[I-1][j-1]*(v_P - v_S) ;
                    double U = v_P + alpha_s_v[I-1][j-1]*(v_S - v_P) ;
                    double UU = v_N + alpha_s_v[I-1][j-1]*(v_SS - v_N) ;
                    r_s_v[I-1][j-1] = (U-UU)/(D-U+1e-30) ; //avoid division by zero
                }
                //north boundary
                if(j == Ny){
                    r_n_v[I-1][j-1] = -1.0 ; //north boundary uses upwind weighting
                }else{
                    double v_P = v_star[I-1][j-1] ;
                    double v_N = v_star[I-1][j] ;
                    double v_S = v_star[I-1][j-2] ;
                    double v_NN = v_star[I-1][j+1] ;
                    double D = v_P + alpha_n_v[I-1][j-1]*(v_N - v_P) ;
                    double U = v_N + alpha_n_v[I-1][j-1]*(v_P - v_N) ;
                    double UU = v_NN + alpha_n_v[I-1][j-1]*(v_S - v_NN) ;
                    r_n_v[I-1][j-1] = (U-UU)/(D-U+1e-30) ;
                }
                //west boundary
                if(I == 1){
                    r_w_v[I-1][j-1] = 0.0 ; //left-most control volumes do not need the west ratio
                }else if (I == 2){
                    double v_W = v_star[I-2][j-1] ;
                    double v_P = v_star[I-1][j-1] ;
                    double v_WW = 2*v_left-v_W ; //mirror the west neighbor across the boundary
                    double v_E = v_star[I][j-1] ;
                    double D = v_W + alpha_w_v[I-1][j-1]*(v_P - v_W) ;
                    double U = v_P + alpha_w_v[I-1][j-1]*(v_W - v_P) ;
                    double UU = v_E + alpha_w_v[I-1][j-1]*(v_WW - v_E) ;
                    r_w_v[I-1][j-1] = (U-UU)/(D-U+1e-30) ;
                }else if(I == Nx){
                    double v_W = v_star[I-2][j-1] ;
                    double v_P = v_star[I-1][j-1] ;
                    double v_WW = v_star[I-3][j-1] ;
                    double v_E = 2*v_right - v_P ; //mirror the point across the east boundary
                    double D = v_W + alpha_w_v[I-1][j-1]*(v_P - v_W) ;
                    double U = v_P + alpha_w_v[I-1][j-1]*(v_W - v_P) ;
                    double UU = v_E + alpha_w_v[I-1][j-1]*(v_WW - v_E) ;
                    r_w_v[I-1][j-1] = (U-UU)/(D-U+1e-30) ;
                }else{
                    double v_W = v_star[I-2][j-1] ;
                    double v_P = v_star[I-1][j-1] ;
                    double v_WW = v_star[I-3][j-1] ;
                    double v_E = v_star[I][j-1] ;
                    double D = v_W + alpha_w_v[I-1][j-1]*(v_P - v_W) ;
                    double U = v_P + alpha_w_v[I-1][j-1]*(v_W - v_P) ;
                    double UU = v_E + alpha_w_v[I-1][j-1]*(v_WW - v_E) ;
                    r_w_v[I-1][j-1] = (U-UU)/(D-U+1e-30) ;
                }
                //east boundary
                if(I == 1){
                    double v_P = v_star[I-1][j-1] ;
                    double v_E = v_star[I][j-1] ;
                    double v_W = 2*v_left - v_P ; //mirror the point across the west boundary
                    double v_EE = v_star[I+1][j-1] ;
                    double D = v_P + alpha_e_v[I-1][j-1]*(v_E - v_P) ;
                    double U = v_E + alpha_e_v[I-1][j-1]*(v_W - v_E) ;
                    double UU = v_EE + alpha_e_v[I-1][j-1]*(v_W - v_EE) ;
                    r_e_v[I-1][j-1] = (U-UU)/(D-U+1e-30) ;
                }else if (I == Nx-1){
                    double v_P = v_star[I-1][j-1] ;
                    double v_E = v_star[I][j-1] ;
                    double v_W = v_star[I-2][j-1] ;
                    double v_EE = 2*v_right - v_E ; //mirror the east neighbor across the boundary
                    double D = v_P + alpha_e_v[I-1][j-1]*(v_E - v_P) ;
                    double U = v_E + alpha_e_v[I-1][j-1]*(v_W - v_E) ;
                    double UU = v_EE + alpha_e_v[I-1][j-1]*(v_W - v_EE) ;
                    r_e_v[I-1][j-1] = (U-UU)/(D-U+1e-30) ;
                }else if(I == Nx){
                    r_e_v[I-1][j-1] = 0.0 ; //right-most control volumes do not need the east ratio
                }else {
                    double v_P = v_star[I-1][j-1] ;
                    double v_E = v_star[I][j-1] ;
                    double v_W = v_star[I-2][j-1] ;
                    double v_EE = v_star[I+1][j-1] ;
                    double D = v_P + alpha_e_v[I-1][j-1]*(v_E - v_P) ;
                    double U = v_E + alpha_e_v[I-1][j-1]*(v_W - v_E) ;
                    double UU = v_EE + alpha_e_v[I-1][j-1]*(v_W - v_EE) ;
                    r_e_v[I-1][j-1] = (U-UU)/(D-U+1e-30) ;
                }
            }
        }
    }
}

//Export finite-volume cell data to VTK
void output_combined_fields(int m) {
    int nodes_i = Nx + 1;  // number of nodes along x
    int nodes_j = Ny + 1;  // number of nodes along y
    int cells_i = Nx;      // number of cells along x
    int cells_j = Ny;      // number of cells along y
    
    std::ostringstream name;
    name << "3.TVD_Re" << Re << "_Grid" << Nx+1 << "x" << Ny+1
         << "_DC_step" << std::setfill('0') << std::setw(6) << m 
         << "_error" << min_step7 <<
          ".vtk";

    std::ofstream out(name.str());
    
    out << "# vtk DataFile Version 3.0\n";
    out << "3.TVD_Re = 5000_Grid "<< Nx+1 << "x" << Ny+1 << "\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_GRID\n";
    out << "DIMENSIONS " << nodes_i << " " << nodes_j << " 1\n";
    
    // write node coordinates
    out << "POINTS " << nodes_i * nodes_j << " double\n";
    for(int j = 0; j < nodes_j; j++) {
        for(int i = 0; i < nodes_i; i++) {
            double x_coord = i * dx;
            double y_coord = j * dy;
            out << scientific << setprecision(6) 
                << x_coord << " " << y_coord << " 0.0\n";
        }
    }
    
    out << "CELL_DATA " << cells_i * cells_j << "\n";
    
    // 1. pressure field
    out << "SCALARS Pressure double 1\n";
    out << "LOOKUP_TABLE default\n";
    for(int J = 1; J <= cells_j; J++) {
        for(int I = 1; I <= cells_i; I++) {
            out << scientific << setprecision(6) 
                << p_star[I-1][J-1] << "\n";
        }
    }
    
    // 2. velocity field (CELL_DATA already declared)
    out << "VECTORS Velocity double\n";
    for(int J = 1; J <= cells_j; J++) {
        for(int I = 1; I <= cells_i; I++) {
            double u_val = 0.5 * (u_star[I-1][J-1] + u_star[I][J-1]);
            double v_val = 0.5 * (v_star[I-1][J-1] + v_star[I-1][J]);
            
            out << scientific << setprecision(6)
                << u_val << " " << v_val << " 0.0\n";
        }
    }
    
    out.close();
    cout << "VTK document has been output: " << name.str() << endl;
}

int main(){
    // initialize OpenMP with the available resources
    int num_threads = omp_get_max_threads();
    omp_set_num_threads(num_threads);
    cout << "Using " << num_threads << " threads" << endl;
    cout << "Using FVD with SIMPLE algorithm to solve 2D incompressible Navier-Stokes equations" << endl;
    cout << "Grid Size:" << Nx+1 << "x" << Ny+1 << endl;
    cout << "Grid spacing : " << dx << "x " << dy << endl;
    cout << "Boundary Conditions: " ;
    cout << "u: left = " << u_left << ", right = " << u_right << ", bottom = " << u_bottom << ", top = " << u_top << endl;
    cout << "v: left = " << v_left << ", right = " << v_right << ", bottom = " << v_bottom << ", top = " << v_top << endl;
    cout << "Re = " << Re << ",rho = " << rho << endl ;
    CE_physical = false;
    G = 0 ;
    initialize() ;
    min_step7 = 50.0 ;
    double min_G = 0.0  ; 
    for(G = 0 ; G <= max_G ; G++){
        
        step2(G) ;//update the velocity predictor (first update)
        step3() ;//update the pressure correction
        step4(G) ;//update the corrected pressure field
        step5(G) ;//update the velocity correction
        step6(G) ;//update the corrected velocity field
        //step7: convergence check
        if(min_step7 > step7(G)){ //evaluate each iteration
                min_step7 = step7(G) ;
                min_G = G ; 
                //if(G > 1000 && min_step7 < 2.55   ){output_combined_fields(min_G) ;}
               
            }
        if(G%1000 == 0){
            cout << "Itteration step = " << G << endl ;
            cout << "Non-equilibrium mass flow rate for (G = "<< G << " )step:" << step5star(G) << endl ;
            cout << "currnt minerror(G =" << min_G << ")= "<< min_step7 << endl ; 
        }
        if(G>0 && step5star(G) < 1e-10){
            CE_physical = true ;
            cout << "Haved converged at step G = " << G << endl ;
            cout << "Non-equilibrium mass flow rate for (G = "<< G << " )step:" << step5star(G) <<endl ;
            cout << "currnt minerror(G =" << G << ")= "<< min_step7 << endl ;
            output_combined_fields(G) ;
            break ;
        }
        step8(G) ;//refresh the pressure predictor
        step9(G) ;//refresh velocity predictors, convection coefficients, and interface ratios
    }
    if(!CE_physical){
        cout << "Haved converged at step G = " << G << endl ;
        cout << "Non-equilibrium mass flow rate for (G = "<< G << " )step:" << step5star(G)  <<endl ;
        cout << "currnt minerror(G =" << G << ")= "<< min_step7 << endl ;
        output_combined_fields(min_G) ;
    }
    cout << "Grid Size : " << Nx << "x" << Ny << "Completed." << endl ;
    return 0 ;
}



