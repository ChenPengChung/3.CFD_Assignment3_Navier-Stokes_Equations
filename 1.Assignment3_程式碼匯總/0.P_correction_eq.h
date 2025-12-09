// This code implements the inner iteration routine for the pressure correction equation.
// All four boundaries use Neumann conditions.
// The solver uses Jacobi iterations combined with under-relaxation.
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include <iomanip>
#include <omp.h> 
using namespace std ;
// Parameters used in the solver.
constexpr int Nx = 80;
constexpr int Ny = 80;
constexpr double dx = 1.0 / double(Nx);
constexpr double dy = 1.0 / double(Ny);
constexpr double rho = 1.0;
int n = Nx * Ny;
vector<double> a; // 1D array for left-hand side coefficients.
vector<double> b; // 1D array for right-hand side source terms.
vector<double> x;
vector<double> x_old;
vector<double> x_old2;
vector<double> x_nuse ;
extern double u_star[Nx+1][Ny] ; // Predicted x-velocity component.
extern double v_star[Nx][Ny+1] ; // Predicted y-velocity component.
extern double a_P[Nx+1][Ny] ; // Cell-centered coefficients for x-velocity.
extern double b_P[Nx][Ny+1] ; // Cell-centered coefficients for y-velocity.
extern double p_prime[Nx][Ny] ; // Pressure correction field.
extern int G;
const double lamda = 0.1 ;
int F, max_F = 4000;
double maxerror; 
const float tolerance = 1e-4;
bool steadystate;



// Initialize matrices.
void initial(vector<double>& a, vector<double>& b , vector<double>& x) {
    // Initialize all entries to zero.
    #pragma omp parallel for
    for(int i = 1; i <= n; i++) {
        a[i] = 0.0;
        b[i] = 0.0;
        x[i] = 0.0;
        x_old[i] = 0.0;
    }
}

void GS(vector<double>& a, vector<double>& b , vector<double>& x , int n) {
    // Copy current solution into x_old for error evaluation.
     #pragma omp parallel for
    for(int k = 1; k <= n; k++) {
        x_old[k] = x[k];
    }
    
    // Note: Gauss-Seidel has data dependencies, so interior updates cannot be parallelized directly.
    // Use Red-Black ordering to enable parallelization.
    
    // Red nodes (I + J even).
    #pragma omp parallel for collapse(2)    // Key section.
    for(int I = 2; I <= Nx-1; I++) {
        for(int J = 2; J <= Ny-1; J++) {
            if((I+J) % 2 == 0) {  // Red node
                int idx = (J-1)*Nx + I ;
                double sum2 = 0.0 ;
                 a[idx] = (((dy*dy)/(rho*(a_P[I][J-1])))
                        +((dy*dy)/(rho*(a_P[I-1][J-1])))
                        +((dx*dx)/(rho*(b_P[I-1][J])))
                        +((dx*dx)/(rho*(b_P[I-1][J-1])))) ;// Diagonal coefficient
                 a[idx+1] = -((dy*dy)/(rho*(a_P[I][J-1]))) ; // East coefficient
                 a[idx-1] = -((dy*dy)/(rho*(a_P[I-1][J-1]))) ; // West coefficient
                 a[idx+Nx] =-((dx*dx)/(rho*(b_P[I-1][J]))) ;// North coefficient
                 a[idx-Nx] =-((dx*dx)/(rho*(b_P[I-1][J-1])));// South coefficient
                b[idx] = dy*u_star[I-1][J-1] - dy*u_star[I][J-1] + dx*v_star[I-1][J-1] - dx*v_star[I-1][J] ;
                sum2 = a[idx+1]*x[idx+1] + a[idx-1]*x[idx-1] 
                     + a[idx+Nx]*x[idx+Nx] + a[idx-Nx]*x[idx-Nx];
                x[idx] = (b[idx] - sum2) / a[idx];
                x[idx] = x[idx]*lamda + x_old[idx]*(1-lamda);
            }
        }
    }
    
    // Black nodes (I + J odd).
    #pragma omp parallel for collapse(2)
    for(int I = 2; I <= Nx-1; I++) {
        for(int J = 2; J <= Ny-1; J++) {
            if((I+J) % 2 == 1) {  // Black node
                int idx = (J-1)*Nx + I ;
                double sum2 = 0.0 ;
                a[idx] = (((dy*dy)/(rho*(a_P[I][J-1])))
                         +((dy*dy)/(rho*(a_P[I-1][J-1])))
                         +((dx*dx)/(rho*(b_P[I-1][J])))
                         +((dx*dx)/(rho*(b_P[I-1][J-1])))) ;// Diagonal coefficient
                a[idx+1] = -((dy*dy)/(rho*(a_P[I][J-1]))) ; // East coefficient
                a[idx-1] = -((dy*dy)/(rho*(a_P[I-1][J-1]))) ; // West coefficient
                a[idx+Nx] =-((dx*dx)/(rho*(b_P[I-1][J]))) ;// North coefficient
                a[idx-Nx] =-((dx*dx)/(rho*(b_P[I-1][J-1])));// South coefficient
                b[idx] = dy*u_star[I-1][J-1] - dy*u_star[I][J-1] + dx*v_star[I-1][J-1] - dx*v_star[I-1][J] ;
                sum2 += a[idx+1] *x[idx+1] ;
                sum2 += a[idx-1] *x[idx-1] ;
                sum2 += a[idx+Nx] *x[idx+Nx] ;
                sum2 += a[idx-Nx] *x[idx-Nx] ;
                x[idx] = (b[idx] - sum2) / a[idx];
                x[idx] = x[idx]*lamda + x_old[idx]*(1-lamda);
            }
        }
    }
    
    // Right boundary (excluding corners).
    #pragma omp parallel for
    for(int J = 2; J <= Ny-1; J++) {
        int idx = Nx*(J-1)+Nx; 
        a[idx] = ( 
                         +((dy*dy)/(rho*(a_P[Nx-1][J-1])))
                         +((dx*dx)/(rho*(b_P[Nx-1][J])))
                        +((dx*dx)/(rho*(b_P[Nx-1][J-1])))) ;// Diagonal coefficient
        a[idx-1] = -((dy*dy)/(rho*(a_P[Nx-1][J-1]))) ; // West coefficient
        a[idx+Nx] = -((dx*dx)/(rho*(b_P[Nx-1][J]))) ;// North coefficient
        a[idx-Nx] = -((dx*dx)/(rho*(b_P[Nx-1][J-1]))) ;// South coefficient
        b[idx] = dy*u_star[Nx-1][J-1] - dy*u_star[Nx][J-1] + dx*v_star[Nx-1][J-1] - dx*v_star[Nx-1][J] ;
        double sum2 = 0.0 ;
        sum2 += a[idx-1] * x[idx-1] ;
        sum2 += a[idx+Nx] * x[idx+Nx] ;
        sum2 += a[idx-Nx] * x[idx-Nx] ;
        x[idx] = (b[idx] - sum2) / a[idx];
        x[idx] = x[idx]*lamda + x_old[idx]*(1-lamda);
    }
    
    // Left boundary (excluding corners).
    #pragma omp parallel for
    for(int J = 2; J <= Ny-1; J++) {
        int idx = Nx*(J-1)+1;
        a[idx] = (((dy*dy)/(rho*(a_P[1-1+1][J-1])))
                         +((dx*dx)/(rho*(b_P[1-1][J-1+1])))
                         +((dx*dx)/(rho*(b_P[1-1][J-1])))) ;// Diagonal coefficient
        a[idx+1] = -((dy*dy)/(rho*(a_P[1][J-1]))) ; // East coefficient
        a[idx+Nx] = -((dx*dx)/(rho*(b_P[1-1][J]))) ;// North coefficient
        a[idx-Nx] = -((dx*dx)/(rho*(b_P[1-1][J-1]))) ;// South coefficient
        b[idx] = dy*u_star[1-1][J-1] - dy*u_star[1][J-1] + dx*v_star[1-1][J-1] - dx*v_star[1-1][J]  ;
        double sum2 = 0.0 ;
        sum2 += a[idx+1] * x[idx+1] ;
        sum2 += a[idx+Nx] * x[idx+Nx] ;
        sum2 += a[idx-Nx] * x[idx-Nx] ;
        x[idx] = (b[idx] - sum2) / a[idx];
        x[idx] = x[idx]*lamda + x_old[idx]*(1-lamda);
    }
    
    // Top boundary (excluding corners).
    #pragma omp parallel for
    for(int I = 2; I <= Nx-1; I++) {
        int idx = n - Nx + I ;
        a[idx] = (((dy*dy)/(rho*(a_P[I][Ny-1])))
                    +((dy*dy)/(rho*(a_P[I-1][Ny-1])))
                    +((dx*dx)/(rho*(b_P[I-1][Ny-1])))) ;// Diagonal coefficient
        a[idx+1] = -((dy*dy)/(rho*(a_P[I][Ny-1]))) ; 
        a[idx-1] = -((dy*dy)/(rho*(a_P[I-1][Ny-1]))) ;
        a[idx-Nx] = -((dx*dx)/(rho*(b_P[I-1][Ny-1]))) ;
        b[idx] =dy*u_star[I-1][Ny-1] - dy*u_star[I][Ny-1] + dx*v_star[I-1][Ny-1] - dx*v_star[I-1][Ny] ;
        double sum2 = 0.0 ;
        sum2 += a[idx+1] * x[idx+1] ;
        sum2 += a[idx-1] * x[idx-1] ;
        sum2 += a[idx-Nx] * x[idx-Nx] ;
        x[idx] = (b[idx] - sum2) / a[idx];
        x[idx] = x[idx]*lamda + x_old[idx]*(1-lamda) ; 
    }
    
    // Bottom boundary (excluding corners).
    #pragma omp parallel for
    for(int I = 2; I <= Nx-1 ; I++) {
        a[I] = (((dy*dy)/(rho*(a_P[I][1-1])))
                    +((dy*dy)/(rho*(a_P[I-1][1-1])))
                    +((dx*dx)/(rho*(b_P[I-1][1])))) ;// Diagonal coefficient
        a[I+1] =-((dy*dy)/(rho*(a_P[I][1-1])));
        a[I-1] =-((dy*dy)/(rho*(a_P[I-1][1-1]))) ;
        a[I+Nx] =-((dx*dx)/(rho*(b_P[I-1][1])));
        b[I] = dy*u_star[I-1][1-1] - dy*u_star[I][1-1] + dx*v_star[I-1][1-1] - dx*v_star[I-1][1] ;
        double sum2 = 0.0 ;
        sum2 += a[I+1] * x[I+1] ;
        sum2 += a[I-1] * x[I-1] ;
        sum2 += a[I+Nx] * x[I+Nx] ;
        x[I] = (b[I] - sum2) / a[I];
        x[I] = x[I]*lamda + x_old[I]*(1-lamda);
    }
    
    // Four corner points handled sequentially.
    {
        a[1] = ((dy*dy)/(rho*(a_P[1][1-1])) 
                         //+((dy*dy)/(rho*SAFE(a_P[I-1][J-1])))
                        + ((dx*dx)/(rho*(b_P[1-1][1]))));
                         //+((dx*dx)/(rho*SAFE(b_P[I-1][J-1])))) // Diagonal coefficient
        a[1+1] = -((dy*dy)/(rho*(a_P[1][1-1]))) ; // East coefficient
        //a[1][1-1] = -((dy*dy)/(rho*SAFE(a_P[I-1][J-1]))) ; // West coefficient
        a[1+Nx] =-((dx*dx)/(rho*(b_P[1-1][1]))) ;// North coefficient
        //a[1][1-Nx] =-((dx*dx)/(rho*SAFE(b_P[I-1][J-1])));// South coefficient
        b[1] = dy*u_star[1-1][1-1] - dy*u_star[1][1-1] + dx*v_star[1-1][1-1] - dx*v_star[1-1][1] ;
        double sum2 = 0.0 ;
        sum2 += a[1+1] * x[1+1] ;
        sum2 += a[1+Nx] * x[1+Nx] ;
        x[1] = (b[1] - sum2) / a[1];
        x[1] = x[1]*lamda + x_old[1]*(1-lamda);
    }
    {
        a[Nx] = (
                +((dy*dy)/(rho*(a_P[Nx-1][1-1])))
            +((dx*dx)/(rho*(b_P[Nx-1][1])))) ;// Diagonal coefficient
        a[Nx-1] = -((dy*dy)/(rho*(a_P[Nx-1][1-1])));
        a[2*Nx] = -((dx*dx)/(rho*(b_P[Nx-1][1])));
        b[Nx] = dy*u_star[Nx-1][1-1] - dy*u_star[Nx][1-1] + dx*v_star[Nx-1][1-1] - dx*v_star[Nx-1][1] ;
        double sum2 = 0.0 ;
        sum2 += a[Nx-1] * x[Nx-1] ;
        sum2 += a[2*Nx] * x[2*Nx] ;
        x[Nx] = (b[Nx] - sum2) / a[Nx];
        x[Nx] = x[Nx]*lamda + x_old[Nx]*(1-lamda);
    }
    {
        a[n-Nx+1] = (((dy*dy)/(rho*(a_P[1][Ny-1])))
            +((dx*dx)/(rho*(b_P[1-1][Ny-1])))) ;// Diagonal coefficient
        a[n-Nx+2] =  -((dy*dy)/(rho*(a_P[1][Ny-1]))) ;
        a[n-2*Nx+1] = -((dx*dx)/(rho*(b_P[1-1][Ny-1])));
        b[n-Nx+1] =  dy*u_star[1-1][Ny-1] - dy*u_star[1][Ny-1] + dx*v_star[1-1][Ny-1] - dx*v_star[1-1][Ny] ;
        double sum2 = 0.0 ;
        sum2 += a[n-Nx+2] * x[n-Nx+2] ;
        sum2 += a[n-2*Nx+1] * x[n-2*Nx+1] ;
        x[n-Nx+1] = (b[n-Nx+1] - sum2) / a[n-Nx+1];
        x[n-Nx+1] = x[n-Nx+1]*lamda + x_old[n-Nx+1]*(1-lamda);
    }
    {
        a[n] =( +((dy*dy)/(rho*(a_P[Nx-1][Ny-1])))+((dx*dx)/(rho*(b_P[Nx-1][Ny-1])))) ;// Diagonal coefficient
        a[n-1] = -((dy*dy)/(rho*(a_P[Nx-1][Ny-1]))) ;
        a[n-Nx] =-((dx*dx)/(rho*(b_P[Nx-1][Ny-1])));
        b[n] = dy*u_star[Nx-1][Ny-1] - dy*u_star[Nx][Ny-1] + dx*v_star[Nx-1][Ny-1] - dx*v_star[Nx-1][Ny] ;
        double sum2 = 0.0 ;
        sum2 += a[n-1] * x[n-1] ;
        sum2 += a[n-Nx] * x[n-Nx] ;
        x[n] = (b[n] - sum2) / a[n];
        x[n] = x[n]*lamda + x_old[n]*(1-lamda);
    }
    
    // Compute maximum error via reduction.
    maxerror = 0;
    #pragma omp parallel for reduction(max:maxerror)
    for(int k = 1; k <= n; k++) {
        double error = fabs(x[k] - x_old[k]);
        if(maxerror < error) {
            maxerror = error;
        }
    }
}

void setreference(int m) {
    // Compute the domain average.
    double mean = 0.0;
    #pragma omp parallel for reduction(+:mean)
    for(int k=1; k<=n; ++k) {
        mean += x[k];
    }
    mean /= n;
    
    // Subtract the average.
    #pragma omp parallel for
    for(int k=1; k<=n; ++k) {
        x[k] -= mean;
    }
}

void step3() {
    // Optional: set OpenMP thread count.
    // omp_set_num_threads(8);  // Adjust to match available CPU cores.
    // Resize working arrays.
    a.assign(n+2,  0.0);
    b.assign(n+2, 0.0);
    x.assign(n+2, 0.0);
    x_old.assign(n+2, 0.0);
    x_old2.assign(n+2, 0.0);
    x_nuse.assign(n+2, 0.0);
    
    
    steadystate = false;
    initial(a, b, x);
    
    #pragma omp parallel for
    for (int k=1; k<=n; ++k) {
        x_old[k] = x[k];
    }
    
    for(F = 0; F < max_F ; F++) {
        if(F > 0) {
            GS(a, b, x, n);
            setreference(F) ;
        } 
        if(F > 100 && maxerror < tolerance) {
            steadystate = true;
            if(G%100 == 0){cout << "G = " << G << "_F = " << F << "_Maxerror = " << maxerror << "achieve the inner itteration converged !!! "<< endl;} 
            break;
        }
    }
    if(!steadystate && G%100 == 0) {
        cout << "\nWarning: Reached maximum iteration count, but not converged to steady state!" << endl;
    }
    // Map the 1D solution back to the 2D pressure field.
    #pragma omp parallel for collapse(2)
    for(int I = 1; I <= Nx; I++) {
        for(int J = 1; J <= Ny; J++) {
            p_prime[I-1][J-1] = x[(J-1)*Nx+ I];
        }
    }
}