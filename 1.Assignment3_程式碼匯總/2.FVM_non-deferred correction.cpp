// Finite volume SIMPLE solver for the Navier-Stokes equations.
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <vector> 
#include <omp.h>  
#include "0.P_correction_eq.h" 
using namespace std ;
// Parameter definitions.
// Grid configuration used by the solver.
// Pressure grid has Nx cells along the x direction.
// u-velocity grid has Nx+1 cells along the x direction.
// v-velocity grid has Ny+1 cells along the y direction.
int G ;  // Current time step (starts from 0), different from inner iterations.
double min_step7 ; 
const int  max_G = 10000000;// Maximum allowed time steps.
// u field indices: i in [0, Nx], j in [0, Ny-1].
// v field indices: i in [0, Nx-1], j in [0, Ny].
// Pressure indices: i in [0, Nx-1], j in [0, Ny-1].
double u_star[Nx+1][Ny] ; // Predicted u-velocity field.
double u_star_old[Nx+1][Ny] ;
double v_star[Nx][Ny+1] ; // Predicted v-velocity field.
double v_star_old[Nx][Ny+1] ;
double p_star[Nx][Ny] ; // Predicted pressure field.
double u_prime[Nx+1][Ny] ; // u-velocity correction.
double v_prime[Nx][Ny+1] ; // v-velocity correction.
double p_prime[Nx][Ny] ; // Pressure correction.
double u[Nx+1][Ny] ; // Corrected u-velocity field.
// Previous corrected u-velocity field.
double u_old[Nx+1][Ny] ;// Needs initialization before iterations.
double v[Nx][Ny+1] ; // Corrected v-velocity field.
// Previous corrected v-velocity field.
double v_old[Nx][Ny+1] ;// Needs initialization before iterations.
double p[Nx][Ny] ; // Corrected pressure field.

const double  Re = 5000.0 ;
// Boundary conditions for the u velocity.
const double u_left  = 0.0 ;// No penetration at the wall.
const double u_right = 0.0 ;// No penetration at the wall.
const double u_bottom = 0.0 ;// No-slip condition.
const double u_top = 1.0 ;// No-slip lid velocity.
// Boundary conditions for the v velocity.
const double v_left  = 0.0 ;// No-slip condition.
const double v_right = 0.0 ;// No-slip condition.
const double v_bottom = 0.0 ;// No penetration at the wall.
const double v_top = 0.0 ;// No penetration at the wall.
const double omega_P = 0.5 ; // Relaxation factor following Patankar recommendations.
const double omega_v  = 0.3 ;
const double omega_u  = 0.3 ;
double max_error ;
const float CE_tolerance = 1e-10; // Convergence check: predicted continuity residual.
bool CE_physical ;

// Convection coefficients shared by u and v discretizations.
const double D_e = (1.0/Re)*(dy/dx);
const double D_w = (1.0/Re)*(dy/dx);
const double D_n = (1.0/Re)*(dy/dx);
const double D_s = (1.0/Re)*(dy/dx);
double r_w_u[Nx+1][Ny] ; // Face selection coefficients for u.
double r_e_u[Nx+1][Ny] ;
double r_s_u[Nx+1][Ny] ;
double r_n_u[Nx+1][Ny] ;
double r_w_v[Nx][Ny+1] ; // Face selection coefficients for v.
double r_e_v[Nx][Ny+1] ;
double r_s_v[Nx][Ny+1] ;
double r_n_v[Nx][Ny+1] ;
// u-momentum convection coefficients; alpha factors determined by F.
double F_w_u[Nx+1][Ny] ;
double alpha_w_u[Nx+1][Ny] ;
double F_e_u[Nx+1][Ny] ;
double alpha_e_u[Nx+1][Ny] ;
double F_s_u[Nx+1][Ny] ;
double alpha_s_u[Nx+1][Ny] ;
double F_n_u[Nx+1][Ny] ;
double alpha_n_u[Nx+1][Ny] ;
// u-momentum linear combination coefficients.
double a_P[Nx+1][Ny] ;
double a_WW[Nx+1][Ny] ;
double a_W[Nx+1][Ny] ;
double a_E[Nx+1][Ny] ;
double a_EE[Nx+1][Ny] ;
double a_SS[Nx+1][Ny] ;
double a_S[Nx+1][Ny] ;
double a_N[Nx+1][Ny] ;
double a_NN[Nx+1][Ny] ; 
double S1[Nx+1][Ny] ; // Source term for u-momentum equation.
// v-momentum convection coefficients; alpha factors determined by F.
double F_w_v[Nx][Ny+1] ;
double alpha_w_v[Nx][Ny+1] ;
double F_e_v[Nx][Ny+1] ;
double alpha_e_v[Nx][Ny+1] ;
double F_s_v[Nx][Ny+1] ;
double alpha_s_v[Nx][Ny+1] ;
double F_n_v[Nx][Ny+1] ;
double alpha_n_v[Nx][Ny+1] ;
// v-momentum linear combination coefficients.
double b_P[Nx][Ny+1] ;
double b_WW[Nx][Ny+1] ;
double b_W[Nx][Ny+1] ;
double b_E[Nx][Ny+1] ;
double b_EE[Nx][Ny+1] ;
double b_SS[Nx][Ny+1] ;
double b_S[Nx][Ny+1] ;
double b_N[Nx][Ny+1] ;
double b_NN[Nx][Ny+1] ;
double S2[Nx][Ny+1] ;// Source term for v-momentum equation.
// Update these matrices immediately after each predictor step.


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


// SIMPLE algorithm workflow.
// Step 1: initialize u, v, and pressure predictor fields.
void initialize(){
    #pragma omp parallel  // Parallel initialization block.
    {
    // Initialize predicted u-velocity field.
     #pragma omp parallel for collapse(2) 
    for(int i = 1 ; i <= Nx+1 ; i++){
        // u-star defined with i in [0, Nx], j in [0, Ny-1].
        for(int J= 1 ; J <= Ny ; J++){
            u_star[i-1][J-1] = 0.0 ;
            if(i == 1) u_star[i-1][J-1] = u_left ; // No-penetration wall.
            if(i == Nx + 1) u_star[i-1][J-1] = u_right ; // No-penetration wall.
        }
    }
    // Initialize previous corrected u field using the predictor.
     #pragma omp parallel for collapse(2) 
    for(int i = 1 ; i <= Nx+1 ; i++){
        for(int J = 1 ; J <= Ny ; J++){
            u_old[i-1][J-1] = u_star[i-1][J-1] ;
        }
    }
    // Initialize predicted v-velocity field.
     #pragma omp parallel for collapse(2) 
    for(int I = 1 ; I <= Nx ; I++){
        // v-star defined with i in [0, Nx-1], j in [0, Ny].
        for(int j = 1 ; j <= Ny+1 ; j++){
            v_star[I-1][j-1] = 0.0 ;
            if(j == 1)v_star[I-1][j-1] = v_bottom ;// No-penetration wall.
            if(j == Ny + 1)v_star[I-1][j-1] = v_top ;// No-penetration wall.
        }
    }
    // Initialize previous corrected v field using the predictor.
     #pragma omp parallel for collapse(2) 
    for(int I = 1 ; I <= Nx ; I++){
        for(int j = 1 ; j <= Ny+1 ; j++){
            v_old[I-1][j-1] = v_star[I-1][j-1] ;
        }
    }
    // Initialize predicted pressure and correction.
     #pragma omp parallel for collapse(2) 
    for(int I = 1 ; I <= Nx ; I++){
        for(int J = 1 ; J<= Ny ; J++){
            p_star[I-1][J-1] =  0.0  ; // Neumann boundary imposed via solver.
            p_prime[I-1][J-1] = 0.0 ;
        }
    }
    // Predictor uses: face-selection factors -> convection coefficients -> linear combination coefficients.
    // Initialization order: face-selection, then convection, then linear coefficients.
    // Face-selection factors for the u field (special treatment near north/south boundaries).
     #pragma omp parallel for collapse(2) 
    for(int i = 1 ; i <= Nx+1 ; i++){
        for(int J = 1 ; J <= Ny ; J++){
            // Trivial cells get zero face factors.
            if(i == 1 || i == Nx+1){
                r_w_u[i-1][J-1] = 0.0 ;
                r_e_u[i-1][J-1] = 0.0 ;
                r_n_u[i-1][J-1] = 0.0 ;
                r_s_u[i-1][J-1] = 0.0 ;
            }else {
                // Handle each boundary separately (west/east/south/north cases).
                // West boundary.
                if(i == 2){
                    r_w_u[i-1][J-1] = -1.0 ; // Upwind at west face.
                }else{
                    r_w_u[i-1][J-1] = 1.0 ; // Interior initialization.
                }
                // East boundary.
                if(i == Nx){
                    r_e_u[i-1][J-1] = -1.0 ; // Upwind at east face.
                }else{
                    r_e_u[i-1][J-1] = 1.0 ;// Interior initialization.
                }
                // South boundary.
                if(J == 1){
                    r_s_u[i-1][J-1] = 0.0 ; // Bottom cells do not need south factor.
                }else if(J == 2){
                    r_s_u[i-1][J-1] = -1.0 ;// Upwind for stability.
                }else if(J == Ny){
                    r_s_u[i-1][J-1] = -1.0 ;// Upwind near top boundary.
                }else{
                    r_s_u[i-1][J-1] = 1.0 ;// Interior initialization.
                }
                // North boundary.
                if(J == 1){
                    r_n_u[i-1][J-1] = -1.0 ; // Upwind near bottom.
                }else if (J == Ny-1){
                    r_n_u[i-1][J-1] = -1.0 ;// Upwind near top.
                }else if(J == Ny){
                    r_n_u[i-1][J-1] = 0.0 ; // Top cells do not need north factor.
                }else{
                    r_n_u[i-1][J-1] = 1.0 ;// Interior initialization.
                }
            }
        }    
    }// After each predictor update, refresh these face-selection factors.
    // Face-selection factors for the v field (special handling east/west).
     #pragma omp parallel for collapse(2) 
    for(int I = 1 ; I <= Nx ; I++){
        for(int j = 1 ; j <= Ny+1 ; j++){
            // Trivial cells receive zero factors.
            if(j == 1 || j == Ny+1){
                r_w_v[I-1][j-1] = 0.0 ;
                r_e_v[I-1][j-1] = 0.0 ;
                r_n_v[I-1][j-1] = 0.0 ;
                r_s_v[I-1][j-1] = 0.0 ;
            }else {
                // Evaluate south/north/west/east boundaries separately.
                // South boundary.
                if(j == 2){
                    r_s_v[I-1][j-1] = -1.0 ; // Upwind at bottom face.
                }else{
                    r_s_v[I-1][j-1] = 1.0 ; // Interior initialization.
                }
                // North boundary.
                if(j == Ny){
                    r_n_v[I-1][j-1] = -1.0 ; // Upwind at top face.
                }else{
                    r_n_v[I-1][j-1] = 1.0 ; // Interior initialization.
                }
                // West boundary.
                if(I == 1){
                    r_w_v[I-1][j-1] = 0.0 ; // Left cells do not need west factor.
                }else if (I == 2){
                    r_w_v[I-1][j-1] = -1.0 ;// Upwind near inlet.
                }else if(I == Nx){
                    r_w_v[I-1][j-1] = -1.0 ;// Upwind near outlet.
                }else{
                    r_w_v[I-1][j-1] = 1.0 ;// Interior initialization.
                }
                // East boundary.
                if(I == 1){
                    r_e_v[I-1][j-1] = -1.0 ; // Upwind near inlet.
                }else if (I == Nx-1){
                    r_e_v[I-1][j-1] = -1.0 ; // Upwind near outlet.
                }else if(I == Nx){
                    r_e_v[I-1][j-1] = 0.0 ; // Right cells do not need east factor.
                }else {
                    r_e_v[I-1][j-1] = 1.0 ; // Interior initialization.
                }
            }
        }
    }
    // Initialize convection coefficient matrices.
    // u-momentum convection coefficients.
     #pragma omp parallel for collapse(2) 
    for(int i = 1 ; i <= Nx+1 ; i++){
        for(int J = 1 ; J <= Ny ; J++){
            if(i == 1 || i == Nx+1) {
                F_w_u[i-1][J-1] = 0.0 ;
                F_e_u[i-1][J-1] = 0.0 ;
                F_n_u[i-1][J-1] = 0.0 ;
                F_s_u[i-1][J-1] = 0.0 ;
            }else{
                F_e_u[i-1][J-1] = 0.5*dy*(u_star[i-1][J-1]+u_star[i][J-1]); // Cell center aligns with storage.
                F_w_u[i-1][J-1] = 0.5*dy*(u_star[i-1][J-1]+u_star[i-2][J-1]);
                F_n_u[i-1][J-1] = 0.5*dx*(v_star[i-1][J]+v_star[i-2][J]);
                F_s_u[i-1][J-1] = 0.5*dx*(v_star[i-1][J-1]+v_star[i-2][J-1]);
            }
        }
    }// Alpha factors for u, indexed by cell centers.
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
    // v-momentum convection coefficients.
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
                F_n_v[I-1][j-1] = 0.5*dx*(v_star[I-1][j-1]+v_star[I-1][j]);// Cell center aligns with storage.
                F_s_v[I-1][j-1] = 0.5*dx*(v_star[I-1][j-1]+v_star[I-1][j-2]);
            }
        }
    }// Alpha factors for v, indexed by cell centers.
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

// Step 2: update predicted u and v velocity fields.
void step2(int m){
    // Update the u-momentum predictor coefficients.
    #pragma omp parallel for collapse(2)
    for(int i = 1 ; i <= Nx+1 ; i++){
        for(int J = 1 ; J <= Ny ; J++){
            // Reset all neighbor coefficients before accumulation.
            a_P[i-1][J-1] = 0.0 ; // Reset pressure coefficients.
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
                S1[i-1][J-1] += u_left ;// No-penetration wall; no pressure gradient.
                u_star[i-1][J-1] = u_left ;
            }else if(i == Nx+1){
                a_P[i-1][J-1] += 1.0 ;
                S1[i-1][J-1] += u_right ;// No-penetration wall; no pressure gradient.
                u_star[i-1][J-1] = u_right ;
            }// Special treatments: near west/east boundaries and first/last rows.
            // West boundary second column.
            else{ // For remaining cells handle west/east/south/north faces separately.
                // West boundary.
                if(i == 2){// Second column uses upwind on the west face.
                    a_W[i-1][J-1] += alpha_w_u[i-1][J-1]*(-F_w_u[i-1][J-1])-D_w;// Remove west diffusive term explicitly.
                    a_P[i-1][J-1] += (1-alpha_w_u[i-1][J-1])*(-F_w_u[i-1][J-1])+D_w;
                    S1[i-1][J-1] += -(-p_star[i-2][J-1])*dy/rho ;
                }else{// Interior west faces use TVD-MUSCL.
                    a_WW[i-1][J-1] += (alpha_w_u[i-1][J-1]*UU(r_w_u[i-1][J-1]))*(-F_w_u[i-1][J-1]) ; // Upwind for west faces.
                    a_W[i-1][J-1] += ((1-alpha_w_u[i-1][J-1])*D(r_w_u[i-1][J-1])+alpha_w_u[i-1][J-1]*U(r_w_u[i-1][J-1]))*(-F_w_u[i-1][J-1]) - D_w ; // Subtract west diffusive term when processing this face.
                    a_E[i-1][J-1] += ((1-alpha_w_u[i-1][J-1])*UU(r_w_u[i-1][J-1]))*(-F_w_u[i-1][J-1]);
                    a_P[i-1][J-1] += ((1-alpha_w_u[i-1][J-1])*U(r_w_u[i-1][J-1])+alpha_w_u[i-1][J-1]*D(r_w_u[i-1][J-1]))*(-F_w_u[i-1][J-1])+D_w ;
                    S1[i-1][J-1] += -(-p_star[i-2][J-1])*dy/rho ;
                }
                // East boundary.
                if(i == Nx){// Second-to-last column uses upwind to the east.
                    a_E[i-1][J-1] += (1-alpha_e_u[i-1][J-1])*(F_e_u[i-1][J-1])-D_e ;// Subtract east diffusive term explicitly.
                    a_P[i-1][J-1] += alpha_e_u[i-1][J-1]*(F_e_u[i-1][J-1])+D_e; 
                    S1[i-1][J-1] += -(p_star[i-1][J-1])*dy/rho ;
                }else{// Interior east faces use TVD-MUSCL.
                    a_EE[i-1][J-1] += (1-alpha_e_u[i-1][J-1])*UU(r_e_u[i-1][J-1])*(F_e_u[i-1][J-1]) ; // Upwind for east faces.
                    a_E[i-1][J-1] += (alpha_e_u[i-1][J-1]*D(r_e_u[i-1][J-1])+(1-alpha_e_u[i-1][J-1])*U(r_e_u[i-1][J-1]))*(F_e_u[i-1][J-1]) - D_e ;// Remove east diffusive contribution handled separately.
                    a_W[i-1][J-1] += (alpha_e_u[i-1][J-1])*UU(r_e_u[i-1][J-1])*(F_e_u[i-1][J-1]) ;
                    a_P[i-1][J-1] += ((1-alpha_e_u[i-1][J-1])*D(r_e_u[i-1][J-1])+alpha_e_u[i-1][J-1]*U(r_e_u[i-1][J-1]))*(F_e_u[i-1][J-1])+D_e;   
                    S1[i-1][J-1] += -(p_star[i-1][J-1])*dy/rho ;
                }
                // South boundary (three cases).
                if(J == 1){// First row adjacent to bottom wall.
                    a_P[i-1][J-1] += 2*D_s ;
                    S1[i-1][J-1] += -(u_bottom)*(-F_s_u[i-1][J-1]) +2*D_s*(u_bottom) ;// No-slip; apply Neumann pressure extrapolation at the south boundary.
                }else if(J == 2){// Second row uses extrapolation instead of the SS node.
                    
                    a_S[i-1][J-1] += ((1-alpha_s_u[i-1][J-1])*D(r_s_u[i-1][J-1])+(alpha_s_u[i-1][J-1])*U(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1])-D_s;
                    a_N[i-1][J-1] += ((1-alpha_s_u[i-1][J-1])*UU(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1]) ;
                    a_P[i-1][J-1] += ((1-alpha_s_u[i-1][J-1])*U(r_s_u[i-1][J-1])+(alpha_s_u[i-1][J-1])*D(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1]) + D_s ;
                    S1[i-1][J-1] +=  -2*(alpha_s_u[i-1][J-1]*UU(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1])*u_bottom 
                    -u_star_old[i-1][J-2]*(-(alpha_s_u[i-1][J-1]*UU(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1])) 
                    ; 
                }else if(J == Ny){
                    // No neighboring cell; use extrapolation.
                    
                    a_SS[i-1][J-1]+= (alpha_s_u[i-1][J-1]*UU(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1]) ;
                    a_S[i-1][J-1] += ((1-alpha_s_u[i-1][J-1])*D(r_s_u[i-1][J-1])+(alpha_s_u[i-1][J-1])*U(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1])-D_s ;
                    a_P[i-1][J-1]+= ((1-alpha_s_u[i-1][J-1])*U(r_s_u[i-1][J-1])+alpha_s_u[i-1][J-1]*D(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1])+D_s;
                    S1[i-1][J-1] += -2*u_top*((1-alpha_s_u[i-1][J-1])*UU(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1])
                    -u_star_old[i-1][J-1]*(-(1-alpha_s_u[i-1][J-1])*UU(r_s_u[i-1][J-1])*(-F_s_u[i-1][J-1])) 
                    ;
                }else{
                    a_SS[i-1][J-1] += (alpha_s_u[i-1][J-1]*UU(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1]) ;
                    a_S[i-1][J-1] += ((1-alpha_s_u[i-1][J-1])*D(r_s_u[i-1][J-1])+(alpha_s_u[i-1][J-1])*U(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1])-D_s ;// Remove south diffusive term already accounted for.
                    a_N[i-1][J-1] += ((1-alpha_s_u[i-1][J-1])*UU(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1]) ;
                    a_P[i-1][J-1] += ((1-alpha_s_u[i-1][J-1])*U(r_s_u[i-1][J-1])+alpha_s_u[i-1][J-1]*D(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1])+D_s ; 
                }
                // North boundary special cases: J==1, J==Ny, J==Ny-1.
                if(J == Ny){
                    a_P[i-1][J-1] += 2*D_n ;
                    S1[i-1][J-1] += -(u_top)*F_n_u[i-1][Ny-1] + 2*D_n*(u_top);// No-slip; apply Neumann pressure extrapolation at the north boundary.
                }else if(J == 1){
                    
                    a_NN[i-1][J-1] += (1-alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1]) ;
                    a_N[i-1][J-1] += (alpha_n_u[i-1][J-1]*D(r_n_u[i-1][J-1])+(1-alpha_n_u[i-1][J-1])*U(r_n_u[i-1][J-1]))*(F_n_u[i-1][J-1]) - D_n ;// Remove north diffusive term when handling this face.
                    a_P[i-1][J-1] += ((1-alpha_n_u[i-1][J-1])*D(r_n_u[i-1][J-1])+(alpha_n_u[i-1][J-1])*U(r_n_u[i-1][J-1]))*(F_n_u[i-1][J-1])+ D_n ;   
                    S1[i-1][J-1] += -2*(alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1])*u_bottom 
                    -u_star_old[i-1][J-1]*(-(alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1]))
                    ;
                }else if(J == Ny-1){
                    // No north cell; apply extrapolation.
                    //a_NN[i-1][J-1] += (1-alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1]) ;
                    
                    a_N[i-1][J-1] += (alpha_n_u[i-1][J-1]*D(r_n_u[i-1][J-1])+(1-alpha_n_u[i-1][J-1])*U(r_n_u[i-1][J-1]))*(F_n_u[i-1][J-1]) - D_n  ;// Remove north diffusive contribution explicitly.
                    a_S[i-1][J-1] += (alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1]) ;
                    a_P[i-1][J-1] += ((1-alpha_n_u[i-1][J-1])*D(r_n_u[i-1][J-1])+alpha_n_u[i-1][J-1]*U(r_n_u[i-1][J-1]))*(F_n_u[i-1][J-1]) +D_n ;
                    S1[i-1][J-1] += -2*u_top*(1-alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1])
                    -u_star_old[i-1][J]*(-(1-alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1])) 
                    ;
                }else{
                    a_NN[i-1][J-1] += (1-alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1]) ;
                    a_N[i-1][J-1] += (alpha_n_u[i-1][J-1]*D(r_n_u[i-1][J-1])+(1-alpha_n_u[i-1][J-1])*U(r_n_u[i-1][J-1]))*(F_n_u[i-1][J-1]) - D_n ;// Remove north diffusive term when handling this face.
                    a_S[i-1][J-1] += (alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1]) ;
                    a_P[i-1][J-1] += ((1-alpha_n_u[i-1][J-1])*D(r_n_u[i-1][J-1])+alpha_n_u[i-1][J-1]*U(r_n_u[i-1][J-1]))*(F_n_u[i-1][J-1])+D_n ;
                }
                double WW = 0.0 , W = 0.0 , E = 0.0 , EE = 0.0 ;
                double SS = 0.0 , S = 0.0 , N = 0.0 , NN = 0.0 ;
                E = u_star[i][J-1] ;
                W = u_star[i-2][J-1] ;
                // West-west neighbor.
                if(i == 2){
                    WW = 0.0 ;
                }else{
                    WW = u_star[i-3][J-1] ;
                }
                // East-east neighbor.
                if(i == Nx){
                    EE = 0.0 ;
                }else{
                    EE = u_star[i+1][J-1] ;
                }
                // South-south neighbor.
                if(J == 1 || J == 2){
                    SS = 0.0 ;
                }else{
                    SS = u_star[i-1][J-3] ;
                }
                // South neighbor.
                if(J == 1){
                    S = 0.0 ;
                }else{
                    S = u_star[i-1][J-2] ;
                }
                // North neighbor.
                if(J == Ny){
                    N = 0.0 ;
                }else{
                    N = u_star[i-1][J] ;
                }
                // North-north neighbor.
                if(J == Ny || J == Ny-1){
                    NN = 0.0 ;
                }else {
                    NN = u_star[i-1][J+1] ;
                }
                u_star[i-1][J-1] = (-a_WW[i-1][J-1]*WW - a_W[i-1][J-1]*W - a_E[i-1][J-1]*E - a_EE[i-1][J-1]*EE
                    - a_SS[i-1][J-1]*SS - a_S[i-1][J-1]*S - a_N[i-1][J-1]*N - a_NN[i-1][J-1]*NN
                    + S1[i-1][J-1]) / a_P[i-1][J-1];
                // Relaxation relation for the predicted velocity.
                u_star[i-1][J-1] = (1-omega_u)*u_star_old[i-1][J-1] + omega_u*u_star[i-1][J-1] ;
            }
        }
    }
   
    // Update the v-momentum predictor coefficients.
    #pragma omp parallel for collapse(2)
    for(int I = 1 ; I <= Nx ; I++){
        for(int j = 1 ; j <= Ny+1 ; j++){
            // Reset stencil coefficients before accumulating contributions.
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
                S2[I-1][j-1] += v_bottom ;// Enforce no-penetration; no pressure gradient on the bottom boundary
                v_star[I-1][j-1] = v_bottom ;
            }else if(j == Ny+1){
                b_P[I-1][j-1] += 1.0 ;
                S2[I-1][j-1] += v_top ;// Enforce no-penetration; no pressure gradient on the top boundary
                v_star[I-1][j-1] = v_top ;
            }// Special treatment for near-boundary rows (left/right first and second rows, second rows adjacent to top/bottom)
            // Bottom boundary second row
            else{ // Remaining cells are categorized by west, east, south, and north boundary cases
                // South boundary
                if(j == 2){// Second row above the bottom uses upwind on the south face
                    b_S[I-1][j-1] += alpha_s_v[I-1][j-1]*(-F_s_v[I-1][j-1])-D_s;// Remove south-face diffusive term when handling the boundary
                    b_P[I-1][j-1] += (1-alpha_s_v[I-1][j-1])*(-F_s_v[I-1][j-1])+D_s;
                    S2[I-1][j-1] += -(-p_star[I-1][j-2])*dx/rho ;
                }else{// Other interior points use TVD-MUSCL on the south face
                    b_SS[I-1][j-1] += (alpha_s_v[I-1][j-1]*UU(r_s_v[I-1][j-1]))*(-F_s_v[I-1][j-1]) ;
                    b_S[I-1][j-1] += ((1-alpha_s_v[I-1][j-1])*D(r_s_v[I-1][j-1])+alpha_s_v[I-1][j-1]*U(r_s_v[I-1][j-1]))*(-F_s_v[I-1][j-1]) - D_s ; // Remove south-face diffusive term when handling the boundary
                    b_N[I-1][j-1] += ((1-alpha_s_v[I-1][j-1])*UU(r_s_v[I-1][j-1]))*(-F_s_v[I-1][j-1]) ;
                    b_P[I-1][j-1] += ((1-alpha_s_v[I-1][j-1])*U(r_s_v[I-1][j-1])+alpha_s_v[I-1][j-1]*D(r_s_v[I-1][j-1]))*(-F_s_v[I-1][j-1])+D_s ;
                    S2[I-1][j-1] += -(-p_star[I-1][j-2])*dx/rho ;
                }
                // North boundary
                if(j == Ny){// Second row below the top uses upwind on the north face
                    b_N[I-1][j-1] += (1-alpha_n_v[I-1][j-1])*(F_n_v[I-1][j-1])-D_n ;// Remove north-face diffusive term when handling the boundary
                    b_P[I-1][j-1] += alpha_n_v[I-1][j-1]*(F_n_v[I-1][j-1])+D_n; 
                    S2[I-1][j-1] += -(p_star[I-1][j-1])*dx/rho ;
                }else{// Other interior points use TVD-MUSCL on the north face
                    b_NN[I-1][j-1] += (1-alpha_n_v[I-1][j-1])*UU(r_n_v[I-1][j-1])*(F_n_v[I-1][j-1]) ;
                    b_N[I-1][j-1] += (alpha_n_v[I-1][j-1]*D(r_n_v[I-1][j-1])+(1-alpha_n_v[I-1][j-1])*U(r_n_v[I-1][j-1]))*(F_n_v[I-1][j-1]) - D_n ;// Remove north-face diffusive term when handling the boundary
                    b_S[I-1][j-1] += (alpha_n_v[I-1][j-1])*UU(r_n_v[I-1][j-1])*(F_n_v[I-1][j-1]) ;
                    b_P[I-1][j-1] += ((1-alpha_n_v[I-1][j-1])*D(r_n_v[I-1][j-1])+alpha_n_v[I-1][j-1]*U(r_n_v[I-1][j-1]))*(F_n_v[I-1][j-1]) + D_n ;
                    S2[I-1][j-1] += -(p_star[I-1][j-1])*dx/rho ;
                }
                // West boundary (consider I = 1, 2, Nx separately)
                if(I == 2){// Second column from the west uses extrapolation instead of the WW node
                    
                    b_W[I-1][j-1] += ((1-alpha_w_v[I-1][j-1])*D(r_w_v[I-1][j-1])+(alpha_w_v[I-1][j-1])*U(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1])-D_w;
                    b_E[I-1][j-1] += ((1-alpha_w_v[I-1][j-1])*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1]) ;
                    b_P[I-1][j-1] += ((1-alpha_w_v[I-1][j-1])*U(r_w_v[I-1][j-1])+(alpha_w_v[I-1][j-1])*D(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1])+D_w;
                    S2[I-1][j-1] +=  -2*(alpha_w_v[I-1][j-1]*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1])*v_left 
                    -v_star_old[I-2][j-1]*(-(alpha_w_v[I-1][j-1]*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1]))
                    ;
                }else if(I == 1){// Westernmost column adopts upwind on the west face
                    b_P[I-1][j-1] +=  2*D_w;// Remove west-face diffusive term when handling the boundary
                    S2[I-1][j-1] += -(-F_w_v[I-1][j-1])*v_left+2*D_w*v_left;
                }else if(I == Nx){// E -> P
                    // When I = Nx the east control volume does not exist; extrapolate from interior values
                    
                    b_WW[I-1][j-1] += (alpha_w_v[I-1][j-1]*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1]);
                    b_W[I-1][j-1] += ((1-alpha_w_v[I-1][j-1])*D(r_w_v[I-1][j-1])+alpha_w_v[I-1][j-1]*U(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1]) -D_w ;
                    b_P[I-1][j-1] += ((1-alpha_w_v[I-1][j-1])*U(r_w_v[I-1][j-1])+alpha_w_v[I-1][j-1]*D(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1]) +D_w;
                    S2[I-1][j-1] += -2*v_right*(((1-alpha_w_v[I-1][j-1])*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1])) 
                    -v_star_old[I-1][j-1]*(-(1-alpha_w_v[I-1][j-1])*UU(r_w_v[I-1][j-1])*(-F_w_v[I-1][j-1]))
                    ;
                }else{// Other interior points use TVD-MUSCL on the west face
                    b_WW[I-1][j-1] += (alpha_w_v[I-1][j-1]*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1]) ;
                    b_W[I-1][j-1] += ((1-alpha_w_v[I-1][j-1])*D(r_w_v[I-1][j-1])+(alpha_w_v[I-1][j-1])*U(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1]) - D_w ; // Remove west-face diffusive term when handling the boundary
                    b_E[I-1][j-1] += ((1-alpha_w_v[I-1][j-1])*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1]);
                    b_P[I-1][j-1] += ((1-alpha_w_v[I-1][j-1])*U(r_w_v[I-1][j-1])+(alpha_w_v[I-1][j-1])*D(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1])+D_w ;
                }
                // East boundary (consider I = 1, Nx-1, Nx separately)
                if(I == 1){// Second column from the east uses extrapolation instead of the EE node
                    
                    b_EE[I-1][j-1] += ((1-alpha_e_v[I-1][j-1])*UU(r_e_v[I-1][j-1]))*F_e_v[I-1][j-1] ;
                    b_E[I-1][j-1] += ((alpha_e_v[I-1][j-1])*D(r_e_v[I-1][j-1])+(1-alpha_e_v[I-1][j-1])*U(r_e_v[I-1][j-1]))*(F_e_v[I-1][j-1])-D_e ;// Remove east-face diffusive term when handling the boundary
                    b_P[I-1][j-1] += ((alpha_e_v[I-1][j-1])*U(r_e_v[I-1][j-1])+(1-alpha_e_v[I-1][j-1])*D(r_e_v[I-1][j-1]))*(F_e_v[I-1][j-1])  +D_e ;
                    S2[I-1][j-1] +=  -2*v_left*(alpha_e_v[I-1][j-1]*UU(r_e_v[I-1][j-1]))*(F_e_v[I-1][j-1]) 
                    -v_star_old[I-1][j-1]*(-(alpha_e_v[I-1][j-1]*UU(r_e_v[I-1][j-1]))*(F_e_v[I-1][j-1]))
                    ;
                }else if(I == Nx){// Easternmost column adopts upwind on the east face
                    b_P[I-1][j-1] += 2*D_e ;// Remove east-face diffusive term when handling the boundary
                    S2[I-1][j-1] += -(F_e_v[I-1][j-1])*v_right +2*D_e*v_right;
                }else if(I == Nx-1){// EE node does not exist; extrapolate to E
                    
                    b_E[I-1][j-1] += ((1-alpha_e_v[I-1][j-1])*U(r_e_v[I-1][j-1])+alpha_e_v[I-1][j-1]*D(r_e_v[I-1][j-1]))*F_e_v[I-1][j-1] 
                    - D_e ;
                    b_W[I-1][j-1] += (alpha_e_v[I-1][j-1]*UU(r_e_v[I-1][j-1]))*F_e_v[I-1][j-1] ;
                    b_P[I-1][j-1] += ((1-alpha_e_v[I-1][j-1])*D(r_e_v[I-1][j-1])+alpha_e_v[I-1][j-1]*U(r_e_v[I-1][j-1]))*F_e_v[I-1][j-1] + D_e ;
                    S2[I-1][j-1] += -2*v_right*((1-alpha_e_v[I-1][j-1])*UU(r_e_v[I-1][j-1]))*F_e_v[I-1][j-1]
                    -v_star_old[I][j-1]*(-(1-alpha_e_v[I-1][j-1])*UU(r_e_v[I-1][j-1])*F_e_v[I-1][j-1])
                    ; 
                }else{// Other interior points use TVD-MUSCL on the east face
                    b_EE[I-1][j-1] += (1-alpha_e_v[I-1][j-1])*UU(r_e_v[I-1][j-1])*(F_e_v[I-1][j-1]) ;
                    b_E[I-1][j-1] += (alpha_e_v[I-1][j-1]*D(r_e_v[I-1][j-1])+(1-alpha_e_v[I-1][j-1])*U(r_e_v[I-1][j-1]))*(F_e_v[I-1][j-1]) - D_e ;// Remove east-face diffusive term when handling the boundary
                    b_W[I-1][j-1] += (alpha_e_v[I-1][j-1])*UU(r_e_v[I-1][j-1])*(F_e_v[I-1][j-1]) ;
                    b_P[I-1][j-1] += ((1-alpha_e_v[I-1][j-1])*D(r_e_v[I-1][j-1])+alpha_e_v[I-1][j-1]*U(r_e_v[I-1][j-1]))*(F_e_v[I-1][j-1])+D_e;
                }
                double WW=0.0 , W=0.0 , E=0.0 , EE=0.0 ;
                double SS=0.0 , S=0.0 , N=0.0 , NN=0.0 ;
                S = v_star[I-1][j-2] ;
                N = v_star[I-1][j] ;
                //SS:
                if(j == 2){
                    SS = 0.0 ;
                }else{
                    SS = v_star[I-1][j-3] ;
                }
                //NN:
                if(j == Ny){
                    NN = 0.0 ;
                }else{
                    NN = v_star[I-1][j+1] ;
                }
                //WW:
                if(I == 1 || I == 2){
                    WW = 0.0 ;
                }else{
                    WW = v_star[I-3][j-1] ;
                }
                //W:
                if(I == 1){
                    W = 0.0 ;
                }else{
                    W = v_star[I-2][j-1] ;
                }
                //E:
                if(I == Nx){
                    E = 0.0 ;
                }else{
                    E = v_star[I][j-1] ;
                }
                //EE:
                if(I == Nx || I == Nx-1){
                    EE = 0.0 ;
                }else{
                    EE = v_star[I+1][j-1] ;
                }
                v_star[I-1][j-1] = (-b_WW[I-1][j-1]*WW - b_W[I-1][j-1]*W - b_E[I-1][j-1]*E - b_EE[I-1][j-1]*EE
                            - b_SS[I-1][j-1]*SS - b_S[I-1][j-1]*S - b_N[I-1][j-1]*N - b_NN[I-1][j-1]*NN
                            + S2[I-1][j-1]) / b_P[I-1][j-1] ;
                // Apply under-relaxation to the predicted velocity.
                v_star[I-1][j-1] = (1-omega_u)*v_star_old[I-1][j-1] + omega_u*v_star[I-1][j-1] ;
            }
        }
    }
}
// Obtained u_star and v_star at the N + 1/2 stage.
// Function: mass flux after the first update of the velocity predictor.
double d(int I, int J , int m){
    // I,J correspond to cell centers; m is the SIMPLE iteration counter, distinct from the velocity time level.
    double u_e = u_star[I+1-1][J-1] ;// Velocity predictor evaluated at N+1/2
    double u_w = u_star[I-1][J-1] ;// Velocity predictor evaluated at N+1/2
    double v_n = v_star[I-1][J+1-1] ;// Velocity predictor evaluated at N+1/2
    double v_s = v_star[I-1][J-1] ;// Velocity predictor evaluated at N+1/2
    return rho*(u_e - u_w)*dy + rho*(v_n - v_s)*dx ;
} // Used to judge imbalance after the first velocity predictor update

// Step 4: update the corrected pressure field.
void step4(int m){
    #pragma omp parallel for collapse(2) 
    for(int I = 1 ; I <= Nx ; I++){ // Pressure field (true, predicted, corrected) is stored at cell centers
        for(int J = 1 ; J <= Ny ; J++){
           p[I-1][J-1] = p_star[I-1][J-1] + omega_P*p_prime[I-1][J-1] ;
        }// p_star[I-1][J-1] evaluated at step N; p_prime[I-1][J-1] at N+1
    }// p_prime[I-1][J-1] evaluated at step N+1
}

// Step 5: update the velocity corrections. Boundary velocities are fixed by conditions, so corrections vanish on physical boundaries.
void step5(int m){
    #pragma omp parallel for collapse(2) 
    // Update x-component velocity correction
    for(int i = 1 ; i <= Nx+1 ; i++){
        for(int J = 1 ; J <= Ny ; J++){
            if(i == 1 || i == Nx+1){
                u_prime[i-1][J-1] = 0.0 ;// Boundary velocities already satisfy conditions and are not iterated
            }else{// Pressure correction gradient: current minus western neighbor
                u_prime[i-1][J-1] = -(p_prime[i-1][J-1] - p_prime[i-2][J-1])*dy/(rho*a_P[i-1][J-1]) ;
            }

        }
    }
    // Update y-component velocity correction
    #pragma omp parallel for collapse(2) 
    for(int I = 1 ; I <= Nx ; I++){
        for(int j = 1 ; j <= Ny+1 ; j++){
            if(j == 1 || j == Ny+1){
                v_prime[I-1][j-1] = 0.0 ;// Boundary velocities already satisfy conditions and are not iterated
            }else{// Pressure correction gradient: current minus southern neighbor
                v_prime[I-1][j-1] = -(p_prime[I-1][j-1] - p_prime[I-1][j-2])*dx/(rho*b_P[I-1][j-1]) ;
            }
        }
    }
}

// Used to compute the non-equilibrium mass flux.
double step5star(int m){
    double sum = 0.0 ;
    #pragma omp parallel for collapse(2) reduction(+:sum)  // *** Added for parallel reduction ***
    for(int I = 1 ; I <= Nx ; I++){
        for(int J = 1 ; J <= Ny ; J++){
            sum = sum + fabs(d(I,J,m)) ;
        }
    }
    return sum ;
}


// Step 6: update the corrected velocity field and refresh the previous field.
void step6(int m){// As noted earlier, physical boundary velocities are prescribed and not iterated
    // Update the corrected velocity for the x component
    #pragma omp parallel for collapse(2) 
    for(int i = 1 ; i <= Nx+1 ; i++){
        for(int J = 1 ; J <= Ny ; J++){
            if(i == 1 ){
                u[i-1][J-1] = 0.0 ;// Boundary velocity satisfied by conditions; no iteration
            }else if (i == Nx+1){
                u[i-1][J-1] = 0.0 ;// Boundary velocity satisfied by conditions; no iteration
            }else{
                u[i-1][J-1] = u_prime[i-1][J-1]+u_star[i-1][J-1];// u_star at N+1/2; u_prime at N+1; u_old at N; u at N+1
            }
        }
    }
    // Refresh previous x-component velocity
    #pragma omp parallel for collapse(2) 
    for(int i = 1 ; i <= Nx+1 ; i++){
        for(int J = 1 ; J <= Ny ; J++){
            u_old[i-1][J-1] = u[i-1][J-1] ; // u at N+1; u_old now stores the current correct velocity
        }
    }
    // Update the corrected velocity for the y component
    #pragma omp parallel for collapse(2) 
    for(int I = 1 ; I <= Nx ; I++){
        for(int j = 1 ; j <= Ny+1 ; j++){
            if(j == 1 ){
                v[I-1][j-1] = 0.0 ;// Boundary velocity satisfied by conditions; no iteration (no-penetration)
            }else if (j == Ny+1){
                v[I-1][j-1] = 0.0 ;// Boundary velocity satisfied by conditions; no iteration (no-penetration)
            }else{
                v[I-1][j-1] = v_prime[I-1][j-1]+v_star[I-1][j-1];// v_star at N+1/2; v_prime at N+1; v_old at N; v at N+1
            }
        }
    }
    // Refresh previous y-component velocity
    #pragma omp parallel for collapse(2) 
    for(int I = 1 ; I <= Nx ; I++){
        for(int j = 1 ; j <= Ny+1 ; j++){
            v_old[I-1][j-1] = v[I-1][j-1] ; // v at N+1; v_old now stores the current correct velocity
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
    
// Step 8: update the pressure predictor at cell centers.
void step8(int m){
    #pragma omp parallel for collapse(2)  // *** Added for parallel execution ***
    for(int I = 1 ; I <= Nx ; I++){
        for(int J = 1 ; J <= Ny ; J++){
            p_star[I-1][J-1] = p[I-1][J-1] ; // p at N+1; p_star at N+1
        }
    }
}    
    
// Step 9: refresh velocity predictors, interface limiters, and convective coefficients.
void step9(int m){
    // Update u-component predictor; boundary velocities are imposed directly.
    #pragma omp parallel for collapse(2)
    for(int i = 1 ; i <= Nx+1 ; i++){
        for(int J = 1 ; J <= Ny ; J++){
            if(i == 1){
                u_star[i-1][J-1] = u_left;
            }else if (i == Nx+1){
                u_star[i-1][J-1] = u_right ;
            }else{
                // Interior points, including top/bottom faces, take the corrected value.
                u_star[i-1][J-1] = u[i-1][J-1]; // u at N+1; u_star at N+1
                u_star_old[i-1][J-1] = u_star[i-1][J-1]; // Cache the predictor immediately
            }
        }
    }
    // Update v-component predictor; boundary velocities are imposed directly.
    #pragma omp parallel for collapse(2)
    for(int I = 1 ; I <= Nx ; I++){
        for(int j = 1 ; j <= Ny+1 ; j++){
            if(j == 1){
                v_star[I-1][j-1] = v_bottom ;
            }else if (j == Ny+1){
                v_star[I-1][j-1] = v_top ;
            }else{
                // Interior points, including left/right faces, take the corrected value.
                v_star[I-1][j-1] = v[I-1][j-1]; // v at N+1; v_star at N+1
                v_star_old[I-1][j-1] = v_star[I-1][j-1]; // Cache the predictor immediately
            }
        }
    }
    // Update convective coefficients
    // Convective coefficients for the u-momentum discrete equation
    #pragma omp parallel for collapse(2)
    for(int i = 1 ; i <= Nx+1 ; i++){
        for(int J = 1 ; J <= Ny ; J++){
            if(i == 1 || i == Nx+1) {
                F_w_u[i-1][J-1] = 0.0 ;
                F_e_u[i-1][J-1] = 0.0 ;
                F_n_u[i-1][J-1] = 0.0 ;
                F_s_u[i-1][J-1] = 0.0 ;
            }else{// Recompute convective flux coefficients
                F_e_u[i-1][J-1] = 0.5*dy*(u_star[i-1][J-1]+u_star[i][J-1]); // Cell centers coincide with storage points
                F_w_u[i-1][J-1] = 0.5*dy*(u_star[i-1][J-1]+u_star[i-2][J-1]);
                F_n_u[i-1][J-1] = 0.5*dx*(v_star[i-1][J]+v_star[i-2][J]);
                F_s_u[i-1][J-1] = 0.5*dx*(v_star[i-1][J-1]+v_star[i-2][J-1]);
            }
        }
    }// alpha_f_u flags use control-volume center indices
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
    // Convective coefficients for the v-momentum discrete equation
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
                F_n_v[I-1][j-1] = 0.5*dx*(v_star[I-1][j-1]+v_star[I-1][j]);// Cell centers coincide with storage points
                F_s_v[I-1][j-1] = 0.5*dx*(v_star[I-1][j-1]+v_star[I-1][j-2]);
            }
        }
    }// alpha_f_v flags use control-volume center indices
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
    // Update interface limiters
    #pragma omp parallel for collapse(2)
    for(int i = 1 ; i <= Nx+1 ; i++){
        for(int J = 1 ; J <= Ny ; J++){
            // Cells that do not require limiter evaluation are set to zero
            if(i == 1 || i == Nx+1){
                r_w_u[i-1][J-1] = 0.0 ;
                r_e_u[i-1][J-1] = 0.0 ;
                r_n_u[i-1][J-1] = 0.0 ;
                r_s_u[i-1][J-1] = 0.0 ;
            }else {
                // Handle each boundary case separately (west, east, south, north)
                // West boundary
                if(i == 2){
                    r_w_u[i-1][J-1] = -1.0 ; // West boundary uses upwind scheme
                }else{
                    double u_W = u_star[i-2][J-1]; 
                    double u_P = u_star[i-1][J-1];
                    double u_WW = u_star[i-3][J-1];
                    double u_E = u_star[i][J-1] ;
                    double D = u_W + alpha_w_u[i-1][J-1]*(u_P - u_W) ;
                    double U = u_P + alpha_w_u[i-1][J-1]*(u_W - u_P) ;
                    double UU = u_E + alpha_w_u[i-1][J-1]*(u_WW - u_E) ;
                    r_w_u[i-1][J-1] = (U-UU)/(D-U+1e-30) ; // Initial formulation
                }
                // East boundary
                if(i == Nx){
                    r_e_u[i-1][J-1] = -1.0 ; // East boundary uses upwind scheme
                }else{
                    double u_P = u_star[i-1][J-1] ;
                    double u_E = u_star[i][J-1] ;
                    double u_W = u_star[i-2][J-1] ;
                    double u_EE = u_star[i+1][J-1] ;
                    double D = u_P + alpha_e_u[i-1][J-1]*(u_E - u_P) ;
                    double U = u_E + alpha_e_u[i-1][J-1]*(u_P - u_E) ;
                    double UU = u_EE + alpha_e_u[i-1][J-1]*(u_W - u_EE) ;
                    r_e_u[i-1][J-1] = (U-UU)/(D-U+1e-30) ; // Initial formulation
                }
                // South boundary
                if(J == 1){
                    r_s_u[i-1][J-1] = 0.0 ; // Bottom boundary cells do not need the south-face limiter
                }else if(J == 2){
                    double u_S = u_star[i-1][J-2] ;
                    double u_P = u_star[i-1][J-1] ;
                    double u_SS = 2*u_bottom-u_S ; // 2*boundary value minus the southern node
                    double u_N = u_star[i-1][J] ;
                    double D = u_S + alpha_s_u[i-1][J-1]*(u_P - u_S) ;
                    double U = u_P + alpha_s_u[i-1][J-1]*(u_S - u_P) ;
                    double UU = u_N + alpha_s_u[i-1][J-1]*(u_SS - u_N) ;
                    r_s_u[i-1][J-1] =  (U-UU)/(D-U+1e-30) ; // Initial formulation
                    r_s_u[i-1][J-1] =  -1.0  ; // Use upwind scheme
                }else if(J == Ny){
                    double u_S = u_star[i-1][J-2] ;
                    double u_P = u_star[i-1][J-1] ;
                    double u_SS = u_star[i-1][J-3] ;
                    double u_N = 2*u_top - u_P ; // 2*boundary value minus the current node
                    double D = u_S + alpha_s_u[i-1][J-1]*(u_P - u_S) ;
                    double U = u_P + alpha_s_u[i-1][J-1]*(u_S - u_P) ;
                    double UU = u_N + alpha_s_u[i-1][J-1]*(u_SS - u_N) ;
                    r_s_u[i-1][J-1] =  (U-UU)/(D-U+1e-30) ; // Initial formulation
                    r_s_u[i-1][J-1] =  -1.0  ; // Use upwind scheme
                }else{
                    double u_S = u_star[i-1][J-2] ;
                    double u_P = u_star[i-1][J-1] ;
                    double u_SS = u_star[i-1][J-3] ;
                    double u_N = u_star[i-1][J] ;
                    double D = u_S + alpha_s_u[i-1][J-1]*(u_P - u_S) ;
                    double U = u_P + alpha_s_u[i-1][J-1]*(u_S - u_P) ;
                    double UU = u_N + alpha_s_u[i-1][J-1]*(u_SS - u_N) ;
                    r_s_u[i-1][J-1] =  (U-UU)/(D-U+1e-30) ; // Initial formulation
                }
                // North boundary
                if(J == 1){
                    double u_P = u_star[i-1][J-1] ;
                    double u_N = u_star[i-1][J] ;
                    double u_S = 2*u_bottom-u_P ; // 2*boundary value minus the current node
                    double u_NN = u_star[i-1][J+1] ;
                    double D = u_P + alpha_n_u[i-1][J-1]*(u_N - u_P) ;
                    double U = u_N + alpha_n_u[i-1][J-1]*(u_P - u_N) ;
                    double UU = u_NN + alpha_n_u[i-1][J-1]*(u_S - u_NN) ;
                    r_n_u[i-1][J-1] =  (U-UU)/(D-U+1e-30) ; // Initial formulation
                    r_n_u[i-1][J-1] =  -1.0  ; // Use upwind scheme
                }else if (J == Ny-1){
                    double u_P = u_star[i-1][J-1] ;
                    double u_N = u_star[i-1][J] ;
                    double u_S = u_star[i-1][J-2] ;
                    double u_NN = 2*u_top - u_N ; // 2*boundary value minus the north node
                    double D = u_P + alpha_n_u[i-1][J-1]*(u_N - u_P) ;
                    double U = u_N + alpha_n_u[i-1][J-1]*(u_P - u_N) ;
                    double UU = u_NN + alpha_n_u[i-1][J-1]*(u_S - u_NN) ;
                    r_n_u[i-1][J-1] =  (U-UU)/(D-U+1e-30) ; // Initial formulation
                    r_n_u[i-1][J-1] =  -1.0  ; // Use upwind scheme
                }else if(J == Ny){
                    r_n_u[i-1][J-1] = 0.0 ; // Top boundary cells do not need the north-face limiter
                }else{
                    double u_P = u_star[i-1][J-1] ;
                    double u_N = u_star[i-1][J] ;
                    double u_S = u_star[i-1][J-2] ;
                    double u_NN = u_star[i-1][J+1] ;
                    double D = u_P + alpha_n_u[i-1][J-1]*(u_N - u_P) ;
                    double U = u_N + alpha_n_u[i-1][J-1]*(u_P - u_N) ;
                    double UU = u_NN + alpha_n_u[i-1][J-1]*(u_S - u_NN) ;
                    r_n_u[i-1][J-1] =  (U-UU)/(D-U+1e-30) ; // Initial formulation
                }
            }
        }    
    }// Remember to refresh limiters immediately after updating predictors
    // Interface limiters for the v-component velocity field (special handling on east/west boundaries)
    #pragma omp parallel for collapse(2)
    for(int I = 1 ; I <= Nx ; I++){
        for(int j = 1 ; j <= Ny+1 ; j++){
            // Cells that do not require limiter evaluation are set to zero
            if(j == 1 || j == Ny+1){
                r_w_v[I-1][j-1] = 0.0 ;
                r_e_v[I-1][j-1] = 0.0 ;
                r_n_v[I-1][j-1] = 0.0 ;
                r_s_v[I-1][j-1] = 0.0 ;
            }else {
                // Evaluate boundary cases sequentially (south, north, east, west)
                // South boundary
                if(j == 2){
                    r_s_v[I-1][j-1] = -1.0 ; // Bottom boundary uses an upwind scheme on the south face
                }else{
                    double v_S = v_star[I-1][j-2] ;
                    double v_P = v_star[I-1][j-1] ;
                    double v_SS = v_star[I-1][j-3] ;
                    double v_N = v_star[I-1][j] ;
                    double D = v_S + alpha_s_v[I-1][j-1]*(v_P - v_S) ;
                    double U = v_P + alpha_s_v[I-1][j-1]*(v_S - v_P) ;
                    double UU = v_N + alpha_s_v[I-1][j-1]*(v_SS - v_N) ;
                    r_s_v[I-1][j-1] = (U-UU)/(D-U+1e-30) ; // Avoid zero denominators
                }
                // North boundary
                if(j == Ny){
                    r_n_v[I-1][j-1] = -1.0 ; // Top boundary uses an upwind scheme on the north face
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
                // West boundary
                if(I == 1){
                    r_w_v[I-1][j-1] = 0.0 ; // Westernmost cells do not need the west-face limiter
                }else if (I == 2){
                    double v_W = v_star[I-2][j-1] ;
                    double v_P = v_star[I-1][j-1] ;
                    double v_WW = 2*v_left-v_W ; // 2*boundary value minus the west node
                    double v_E = v_star[I][j-1] ;
                    double D = v_W + alpha_w_v[I-1][j-1]*(v_P - v_W) ;
                    double U = v_P + alpha_w_v[I-1][j-1]*(v_W - v_P) ;
                    double UU = v_E + alpha_w_v[I-1][j-1]*(v_WW - v_E) ;
                    r_w_v[I-1][j-1] = (U-UU)/(D-U+1e-30) ;
                    r_w_v[I-1][j-1] = -1.0 ; // Use upwind scheme
                }else if(I == Nx){
                    double v_W = v_star[I-2][j-1] ;
                    double v_P = v_star[I-1][j-1] ;
                    double v_WW = v_star[I-3][j-1] ;
                    double v_E = 2*v_right - v_P ; // 2*boundary value minus the current node
                    double D = v_W + alpha_w_v[I-1][j-1]*(v_P - v_W) ;
                    double U = v_P + alpha_w_v[I-1][j-1]*(v_W - v_P) ;
                    double UU = v_E + alpha_w_v[I-1][j-1]*(v_WW - v_E) ;
                    r_w_v[I-1][j-1] = (U-UU)/(D-U+1e-30) ;
                    r_w_v[I-1][j-1] = -1.0 ; // Use upwind scheme
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
                // East boundary
                if(I == 1){
                    double v_P = v_star[I-1][j-1] ;
                    double v_E = v_star[I][j-1] ;
                    double v_W = 2*v_left - v_P ; // 2*boundary value minus the current node
                    double v_EE = v_star[I+1][j-1] ;
                    double D = v_P + alpha_e_v[I-1][j-1]*(v_E - v_P) ;
                    double U = v_E + alpha_e_v[I-1][j-1]*(v_W - v_E) ;
                    double UU = v_EE + alpha_e_v[I-1][j-1]*(v_W - v_EE) ;
                    r_e_v[I-1][j-1] = (U-UU)/(D-U+1e-30) ;
                    r_e_v[I-1][j-1] = -1.0 ; // Use upwind scheme
                }else if (I == Nx-1){
                    double v_P = v_star[I-1][j-1] ;
                    double v_E = v_star[I][j-1] ;
                    double v_W = v_star[I-2][j-1] ;
                    double v_EE = 2*v_right - v_E ; // 2*boundary value minus the east node
                    double D = v_P + alpha_e_v[I-1][j-1]*(v_E - v_P) ;
                    double U = v_E + alpha_e_v[I-1][j-1]*(v_W - v_E) ;
                    double UU = v_EE + alpha_e_v[I-1][j-1]*(v_W - v_EE) ;
                    r_e_v[I-1][j-1] = (U-UU)/(D-U+1e-30) ;
                    r_e_v[I-1][j-1] = -1.0 ; // Use upwind scheme
                }else if(I == Nx){
                    r_e_v[I-1][j-1] = 0.0 ; // Easternmost cells do not need the east-face limiter
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

// Output cell data for the finite-volume method.
void output_combined_fields(int m) {
    int nodes_i = Nx + 1;  // Number of nodes along x
    int nodes_j = Ny + 1;  // Number of nodes along y
    int cells_i = Nx;      // Number of cells along x
    int cells_j = Ny;      // Number of cells along y
    std::ostringstream name;
    name << "3.CDs_Re" << Re << "_Grid" << Nx+1 << "x" << Ny+1
         << "_nonDC_step" << std::setfill('0') << std::setw(6) << m 
         << "_error" << min_step7 <<
          ".vtk";

    std::ofstream out(name.str());
    
    out << "# vtk DataFile Version 3.0\n";
    out << "3.CDs_Re = 5000_Grid" << Nx+1 << "x" << Ny+1 << "\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_GRID\n";
    out << "DIMENSIONS " << nodes_i << " " << nodes_j << " 1\n";
    
    // Output node coordinates
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
    
    // 1. Pressure field
    out << "SCALARS Pressure double 1\n";
    out << "LOOKUP_TABLE default\n";
    for(int J = 1; J <= cells_j; J++) {
        for(int I = 1; I <= cells_i; I++) {
            out << scientific << setprecision(6) 
                << p_star[I-1][J-1] << "\n";
        }
    }
    
    // 2. Velocity field (no additional CELL_DATA declaration required)
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
/* //THis is for Tecplot output 
void output_combined_fields(int m) {
    const int cells_i = Nx; // cells in x
    const int cells_j = Ny; // cells in y

    std::ostringstream name;
    name << "3.QUICK_Re" << Re << "_Grid" << Nx << "x" << Ny
                 << "_nonDC_step" << std::setfill('0') << std::setw(6) << m 
                 << "_error" << min_step7 <<
                    ".dat";

    std::ofstream out(name.str());
    
    if (!out) {
        std::cerr << "Error: Cannot open file: " << name.str() << std::endl;
        return;
    }

    // Tecplot header
    out << "TITLE = \"SIMPLE-QUICK CFD Simulation\"\n";
    out << "VARIABLES = \"X\", \"Y\", \"P\", \"U\", \"V\", \"Vel_Mag\", \"Vorticity\"\n";
    out << "ZONE T=\"Time_Step_" << m << "\", "
        << "I=" << cells_i << ", J=" << cells_j << ", "
        << "DATAPACKING=POINT\n"; // Remove trailing comma

    // AUXDATA (quoting strings is more robust)
    out << "AUXDATA Re=\""   << Re        << "\"\n";
    out << "AUXDATA Grid=\"" << Nx << "x" << Ny << "\"\n";
    out << "AUXDATA dx=\""   << dx        << "\"\n";
    out << "AUXDATA dy=\""   << dy        << "\"\n";

    out << std::scientific << std::setprecision(8);

    for (int j = 0; j < cells_j; ++j) {
        for (int i = 0; i < cells_i; ++i) {
            // Cell center
            const double xC = (i + 0.5) * dx;
            const double yC = (j + 0.5) * dy;

            // p at center
            const double pC = p[i][j]; // Assume p has size [Nx][Ny]

            // Interpolate u,v from faces to the cell center (MAC: u(i+0.5,j), v(i,j+0.5))
            // Requires u sized [Nx+1][Ny] and v sized [Nx][Ny+1]
            const double uC = 0.5 * (u[i][j] + u[i+1][j]);
            const double vC = 0.5 * (v[i][j] + v[i][j+1]);

            const double vel_mag = std::sqrt(uC*uC + vC*vC);

            // Vorticity at the center using first differences consistent with the MAC grid
            double vort = 0.0;
            // Need i+1 <= Nx-1 => i <= Nx-2, j-1 >= 0 => j >= 1
            if (i <= cells_i - 2 && j >= 1) {
                const double dvdx = (v[i+1][j] - v[i][j]) / dx; // v: Nx × (Ny+1), use the same j index
                const double dudy = (u[i][j]   - u[i][j-1]) / dy; // u: (Nx+1) × Ny, use the same i index
                vort = dvdx - dudy;
            }
            // One-sided differences could be applied near boundaries to avoid zeros

            out << xC << " " << yC << " "
                << pC << " " << uC << " " << vC << " "
                << vel_mag << " " << vort << "\n";
        }
    }

    out.close();
    std::cout << "Tecplot file successfully written: " << name.str() << std::endl;
}
*/

int main(){
    // *** Added initialization for thread count ***
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
        
        step2(G) ;// Update velocity predictors (first update)
        step3() ;// Update pressure correction
        step4(G) ;// Update corrected pressure
        step5(G) ;// Update velocity corrections
        step6(G) ;// Update corrected velocity field
        // Step 7: convergence assessment
        if(min_step7 > step7(G)){ // Evaluate every iteration
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
        step8(G) ;// Update pressure predictor
        step9(G) ;// Update velocity predictor, convective coefficients, and limiters
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

