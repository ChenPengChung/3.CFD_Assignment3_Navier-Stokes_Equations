// Finite difference SIMPLE solver for the 2D lid-driven cavity using a central difference scheme.
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <vector>
#include <omp.h>
using namespace std ;

// Parameters used in the solver.
const int Nx = 161;
const int Ny = 161;
const int n = Nx * Ny ;
const int n1 = (Nx-2)*(Ny-2) ; 
const double dx = 1.0/double(Nx-1) ;
const double dy = 1.0/double(Ny-1) ;
const double dt = 0.001 ;

int t , T ,G , K ;
const int  max_t = 10000000 ;
const int  max_G = 10000000 ;
double pu_star[Nx-1][Ny]; // Predicted x-momentum field.
double pv_star[Nx][Ny-1]; // Predicted y-momentum field.
double p_star[Nx][Ny]; // Predicted pressure field.
double p[Nx][Ny]; // Corrected pressure field.
double *x  ;
double *x_old ;
double p_prime[Nx][Ny]; // Pressure correction field.
double **a ;// Coefficient matrix.
double *b ;// Right-hand side vector.

const double  rho = 1.0 ;
const double  Re = 100 ;
const double  mu = (1.0/Re) ;

const double pu_bottom = 0.0 ;
const double pu_up = rho ;
const double pv_left = 0.0 ;
const double pv_right = 0.0 ;


double a1 = 2.0/(dx*dx) + 2.0/(dy*dy);
double b1 = -1.0/(dx*dx) ;
double c1 = -1.0/(dy*dy) ;

const double omega = 0.8 ; 
double max_error ;
const float tolerance = 1e-6;
const float CE_tolerance = 1e-10;
bool CE_physical ;
bool steadystate ;

void allocate_memory(){
    x = new double[n1+2];
    x_old = new double[n1+2];
    b = new double[n1+2];
    
    a = new double*[n1+2];
    for(int i = 0; i < n1+2; i++){
        a[i] = new double[n1+2];
    }
    
    cout << "Memory allocated successfully for matrix size " << n1+2 << " x " << n1+2 << endl;
}

void deallocate_memory(){
    delete[] x;
    delete[] x_old;
    delete[] b;
    
    for(int i = 0; i < n1+2; i++){
        delete[] a[i];
    }
    delete[] a;
    
    cout << "Memory deallocated successfully" << endl;
}

double pressure(double x_coord , double y_coord){
    int I = x_coord/dx ;
    int J = y_coord/dy ;
    return p[I][J] ;
}

double x_Velocity(double x_coord , double y_coord){
    int i = (x_coord-(dx/2.0))/dx ;
    int J = y_coord/dy ;
    return pu_star[i][J]/rho;
}

double y_Velocity(double x_coord , double y_coord){
    int I = x_coord/dx ;
    int j = (y_coord-(dy/2.0))/dy ;
    return pv_star[I][j] /rho;
}

double d(int I , int J , int m){
    if(I < 2 || I > Nx-1 || J < 2 || J > Ny-1) return 0.0;
    
    double dudx = (pu_star[I-1][J-1] - pu_star[I-2][J-1])/dx;
    double dvdy = (pv_star[I-1][J-1] - pv_star[I-1][J-2])/dy;
    
    return dudx + dvdy;
}

double get_max_residual(int k , int m) {
    double max_residual = 0.0;
    
    // Parallel residual computation.
    #pragma omp parallel for collapse(2) reduction(max:max_residual)
    for(int I = 2 ; I <= Nx-1; I++ ) {
        for(int J = 2 ; J <= Ny-1; J++) {
            double residual = abs(d(I,J,m));
            if(residual > max_residual) {
                max_residual = residual;
            }
        }
    }
    return max_residual;
}

void time_init(int k){
    // Parallel initialization.
    #pragma omp parallel
    {
        // Initialize x-momentum field.
        #pragma omp for collapse(2)
        for(int i = 0; i < Nx-1; i++){ 
            for(int J = 0; J < Ny; J++){
                pu_star[i][J] = 0.0;
                if(J == 0){
                    pu_star[i][J] = pu_bottom ;
                }
                if(J == Ny-1){
                    pu_star[i][J] = pu_up ;
                }
            }
        }
        
        // Initialize y-momentum field.
        #pragma omp for collapse(2)
        for(int I = 0; I < Nx; I++){
            for(int j = 0; j < Ny-1; j++){
                pv_star[I][j] = 0.0;
                if(I == 0){
                    pv_star[I][j] = pv_left ;
                }
                if(I == Nx-1){
                    pv_star[I][j] = pv_right ;
                }
            }
        }
        
        // Initialize pressure fields.
        #pragma omp for collapse(2)
        for(int I = 0; I < Nx; I++){
            for(int J = 0; J < Ny; J++){
                p[I][J] = 0.0;
                p_star[I][J] = 0.0;
                p_prime[I][J] = 0.0;
            }
        }
    }
}

void comp_itter(double **a, double *b, int n1) {
    // Parallel matrix initialization.
    #pragma omp parallel for
    for(int i = 1 ; i <= n1 ; i++){
        for(int j = 1 ; j <= n1 ; j++){
            a[i][j] = 0.0 ;
        }
        b[i] = 0.0 ;
        x[i] = 0.0 ;
        x_old[i] = 0.0 ;
    }
    a[n1+1][n1+1] = 1.0 ;
    a[0][0] = 0.0 ;
}

void CD(double **a, double *b , int n , int m ){
    // Set interior coefficients in parallel.
    #pragma omp parallel for collapse(2)
    for(int I = 3 ; I <= Nx-2 ; I++){
        for(int J = 3 ; J <= Ny-2 ; J++){
            int idx = (J-2)*(Nx-2) + (I-1) ;
            b[idx] = -rho * d(I,J,m) ;
            a[idx][idx] = a1 ;
            a[idx][idx-1] = b1 ;
            a[idx][idx+1] = b1 ;
            a[idx][idx-(Nx-2)] = c1 ;
            a[idx][idx+(Nx-2)] = c1 ;
        }
    }
    
    // Boundary handling.
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            for(int J = 3 ; J <= Ny-2 ; J++){
                int idx = (J-2)*(Nx-2) + (2-1) ;
                b[idx] = -rho * d(2,J,m) ;
                a[idx][idx] = a1 ;
                a[idx][idx+1] = b1 ;
                a[idx][idx-(Nx-2)] = c1 ;
                a[idx][idx+(Nx-2)] = c1 ;
            }
        }
        
        #pragma omp section
        {
            for(int J = 3 ; J <= Ny-2 ; J++){
                int idx = (J-2)*(Nx-2) + ((Nx-1)-1) ;
                b[idx] = -rho * d((Nx-1),J,m) ;
                a[idx][idx] = a1 ;
                a[idx][idx-1] = b1 ;
                a[idx][idx-(Nx-2)] = c1 ;
                a[idx][idx+(Nx-2)] = c1 ;
            }
        }
        
        #pragma omp section
        {
            for(int I = 3 ; I <= Nx-2 ; I++){
                int idx = (2-2)*(Nx-2) + (I-1) ;
                b[idx] = -rho * d(I,2,m) ;
                a[idx][idx] = a1 ;
                a[idx][idx-1] = b1 ;
                a[idx][idx+1] = b1 ;
                a[idx][idx+(Nx-2)] = c1 ;
            }
        }
        
        #pragma omp section
        {
            for(int I = 3 ; I <= Nx-2 ; I++){
                int idx = ((Ny-1)-2)*(Nx-2) + (I-1) ;
                b[idx] = -rho * d(I,(Ny-1),m) ;
                a[idx][idx] = a1 ;
                a[idx][idx-1] = b1 ;
                a[idx][idx+1] = b1 ;
                a[idx][idx-(Nx-2)] = c1 ;
            }
        }
    }
    
    // Corner handling.
    int p1 = (2-2)*(Nx-2) + (2-1) ; 
    b[p1] = -rho * d(2,2,m) ;
    a[p1][p1] = a1 ;
    a[p1][p1+1] = b1 ;
    a[p1][p1+(Nx-2)] = c1 ;
    
    int p2 = (2-2)*(Nx-2) + ((Nx-1)-1) ;
    b[p2] = -rho * d(Nx-1,2,m) ;
    a[p2][p2] = a1 ;
    a[p2][p2-1] = b1 ;
    a[p2][p2+(Nx-2)] = c1 ;
    
    int p3 = ((Ny-1)-2)*(Nx-2) + (2-1) ;
    b[p3] = -rho * d(2,Ny-1,m) ;
    a[p3][p3] = a1 ;
    a[p3][p3+1] = b1 ;
    a[p3][p3-(Nx-2)] = c1 ;
    
    int p4 = ((Ny-1)-2)*(Nx-2) + ((Nx-1)-1) ;
    b[p4] = -rho * d(Nx-1,Ny-1,m) ;
    a[p4][p4] =  a1 ;
    a[p4][p4-1] = b1 ;
    a[p4][p4-(Nx-2)] = c1 ;
}

void momentum_time(int k , int m){ 
    // Parallel update of the x-momentum field.
    #pragma omp parallel for collapse(2)
    for(int i = 1 ; i < Nx-2 ; i++){
        for(int J = 1 ; J < Ny-1 ; J++){
            double u_old = pu_star[i][J]/rho;
            double u_e = pu_star[i+1][J]/rho;
            double u_w = pu_star[i-1][J]/rho;
            double u_n = (J < Ny-1) ? pu_star[i][J+1]/rho : u_old;
            double u_s = (J > 0) ? pu_star[i][J-1]/rho : u_old;
            
            double conv_x = u_old * (u_e - u_w)/(2*dx);
            
            double v_avg = 0.0;
            if(i > 0 && i < Nx-1 && J > 0 && J < Ny) {
                v_avg = 0.25*(pv_star[i][J] + pv_star[i+1][J] + 
                             pv_star[i][J-1] + pv_star[i+1][J-1])/rho;
            }
            double conv_y = v_avg * (u_n - u_s)/(2*dy);
            
            double diff_x = mu/rho * (u_e - 2*u_old + u_w)/(dx*dx);
            double diff_y = mu/rho * (u_n - 2*u_old + u_s)/(dy*dy);
            
            double dp_dx = (p_star[i+1][J] - p_star[i][J])/dx;
            
            pu_star[i][J] = pu_star[i][J] + dt*(-conv_x - conv_y + diff_x + diff_y - dp_dx/rho) * rho;
        }
    }
    
    // x-momentum boundary conditions.
    #pragma omp parallel for
    for(int i = 0; i < Nx-1; i++){
        pu_star[i][0] = pu_bottom;
        pu_star[i][Ny-1] = pu_up;
    }
    
    // Parallel update of the y-momentum field.
    #pragma omp parallel for collapse(2)
    for(int I = 1 ; I < Nx-1 ; I++){
        for(int j = 1 ; j < Ny-2 ; j++){
            double v_old = pv_star[I][j]/rho;
            double v_n = pv_star[I][j+1]/rho;
            double v_s = pv_star[I][j-1]/rho;
            double v_e = (I < Nx-1) ? pv_star[I+1][j]/rho : v_old;
            double v_w = (I > 0) ? pv_star[I-1][j]/rho : v_old;
            
            double conv_y = v_old * (v_n - v_s)/(2*dy);
            
            double u_avg = 0.0;
            if(I > 0 && I < Nx && j > 0 && j < Ny-1) {
                u_avg = 0.25*(pu_star[I][j] + pu_star[I][j+1] + 
                             pu_star[I-1][j] + pu_star[I-1][j+1])/rho;
            }
            double conv_x = u_avg * (v_e - v_w)/(2*dx);
            
            double diff_x = mu/rho * (v_e - 2*v_old + v_w)/(dx*dx);
            double diff_y = mu/rho * (v_n - 2*v_old + v_s)/(dy*dy);
            
            double dp_dy = (p_star[I][j+1] - p_star[I][j])/dy;
            
            pv_star[I][j] = pv_star[I][j] + dt*(-conv_x - conv_y + diff_x + diff_y - dp_dy/rho) * rho;
        }
    }
    
    // y-momentum boundary conditions.
    #pragma omp parallel for
    for(int j = 0; j < Ny-1; j++){
        pv_star[0][j] = pv_left;
        pv_star[Nx-1][j] = pv_right;
    }
}

void Jacobi(double **a, double *b, double *x, int n1 , int m) {
    // Copy previous solution.
    #pragma omp parallel for
    for(int k = 1; k <= n1; k++) {
        x_old[k] = x[k];
    }
    
    // Set interior coefficients in parallel.
   // For uneven workloads consider dynamic scheduling.
   #pragma omp parallel for collapse(2) schedule(static)
    for(int I = 3 ; I <= Nx-2 ; I++){
        for(int J = 3 ; J <= Ny-2 ; J++){
            int idx = (J-2)*(Nx-2) + (I-1) ;
            double sum = 0.0 ;
            // Five-point stencil.
            sum += a[idx][idx+1] * x_old[idx+1] ;
            sum += a[idx][idx-1] * x_old[idx-1] ;
            sum += a[idx][idx+(Nx-2)] * x_old[idx+(Nx-2)] ;
            sum += a[idx][idx-(Nx-2)] * x_old[idx-(Nx-2)] ;
            x[idx] = (1-omega)*x_old[idx] + (omega/a[idx][idx]) * (b[idx] - sum) ;
        }
    }
    // Handle all four boundaries in a single loop.
#pragma omp parallel
{
    // Left and right boundaries.
    #pragma omp for nowait
    for(int J = 3; J <= Ny-2; J++) {
        // Left boundary.
        int idx = (J-2)*(Nx-2) + (2-1);
        double sum = a[idx][idx+1] * x_old[idx+1] + 
                     a[idx][idx+(Nx-2)] * x_old[idx+(Nx-2)] +
                     a[idx][idx-(Nx-2)] * x_old[idx-(Nx-2)];
        x[idx] = (1-omega)*x_old[idx] + (omega/a[idx][idx]) * (b[idx] - sum);
        
        // Right boundary.
        idx = (J-2)*(Nx-2) + ((Nx-1)-1);
        sum = a[idx][idx-1] * x_old[idx-1] + 
              a[idx][idx+(Nx-2)] * x_old[idx+(Nx-2)] +
              a[idx][idx-(Nx-2)] * x_old[idx-(Nx-2)];
        x[idx] = (1-omega)*x_old[idx] + (omega/a[idx][idx]) * (b[idx] - sum);
    }
    
    // Top and bottom boundaries.
    #pragma omp for nowait
    for(int I = 3; I <= Nx-2; I++) {
        // Bottom boundary.
        int idx = (2-2)*(Nx-2) + (I-1);
        double sum = a[idx][idx+1] * x_old[idx+1] +
                     a[idx][idx-1] * x_old[idx-1] +
                     a[idx][idx+(Nx-2)] * x_old[idx+(Nx-2)];
        x[idx] = (1-omega)*x_old[idx] + (omega/a[idx][idx]) * (b[idx] - sum);
        
        // Top boundary.
        idx = ((Ny-1)-2)*(Nx-2) + (I-1);
        sum = a[idx][idx+1] * x_old[idx+1] +
              a[idx][idx-1] * x_old[idx-1] +
              a[idx][idx-(Nx-2)] * x_old[idx-(Nx-2)];
        x[idx] = (1-omega)*x_old[idx] + (omega/a[idx][idx]) * (b[idx] - sum);
    }
}
#pragma omp parallel sections
{
    #pragma omp section
    {
        int p1 = (2-2)*(Nx-2) + (2-1);
        double sum = a[p1][p1+1] * x_old[p1+1] + a[p1][p1+(Nx-2)] * x_old[p1+(Nx-2)];
        x[p1] = (1-omega)*x_old[p1] + (omega/a[p1][p1]) * (b[p1] - sum);
    }
    
    #pragma omp section
    {
        int p2 = (2-2)*(Nx-2) + ((Nx-1)-1);
        double sum = a[p2][p2-1] * x_old[p2-1] + a[p2][p2+(Nx-2)] * x_old[p2+(Nx-2)];
        x[p2] = (1-omega)*x_old[p2] + (omega/a[p2][p2]) * (b[p2] - sum);
    }
    
    #pragma omp section
    {
        int p3 = ((Ny-1)-2)*(Nx-2) + (2-1);
        double sum = a[p3][p3+1] * x_old[p3+1] + a[p3][p3+(Nx-2)] * x_old[p3+(Nx-2)];
        x[p3] = (1-omega)*x_old[p3] + (omega/a[p3][p3]) * (b[p3] - sum);
    }
    
    #pragma omp section
    {
        int p4 = ((Ny-1)-2)*(Nx-2) + ((Nx-1)-1);
        double sum = a[p4][p4-1] * x_old[p4-1] + a[p4][p4+(Nx-2)] * x_old[p4+(Nx-2)];
        x[p4] = (1-omega)*x_old[p4] + (omega/a[p4][p4]) * (b[p4] - sum);
    }
}

    // Compute maximum error (using reduction).
    max_error = 0.0;
    #pragma omp parallel for reduction(max:max_error)
    for(int i = 1 ; i <= n1 ; i++){ 
        double error = fabs(x[i] - x_old[i]);
        if(error > max_error) max_error = error ;
    }
    
    // Update pressure correction field.
    #pragma omp parallel for collapse(2)
    for(int I = 2 ; I <= Nx-1 ; I++){
        for(int J = 2 ; J <= Ny-1 ; J++){
            p_prime[I-1][J-1] = x[(J-2)*(Nx-2) + (I-1)] ;
        }
    } 
    
    // Apply boundary conditions.
    #pragma omp parallel
    {
        #pragma omp for
        for(int I = 2 ; I <= Nx-1 ; I++){ 
            p_prime[I-1][0] = p_prime[I-1][1] ;
            p_prime[I-1][Ny-1] = p_prime[I-1][Ny-2] ;
        }
        
        #pragma omp for
        for(int J = 2 ; J <= Ny-1 ; J++){
            p_prime[0][J-1] = p_prime[1][J-1] ;
            p_prime[Nx-1][J-1] = p_prime[Nx-2][J-1] ;
        }
    }
    
    p_prime[0][0] = p_prime[1][1];
    p_prime[Nx-1][0] = p_prime[Nx-2][1];
    p_prime[0][Ny-1] = p_prime[1][Ny-2];
    p_prime[Nx-1][Ny-1] = p_prime[Nx-2][Ny-2];
}

void pressure_time(int k , int m){
    // Parallel pressure update.
    #pragma omp parallel
    {
        #pragma omp for collapse(2)
        for(int I = 2 ; I <= Nx-1 ; I++){
            for(int J = 2 ; J <= Ny-1 ; J++){
                p[I-1][J-1] = p_prime[I-1][J-1] + p_star[I-1][J-1] ;
            }
        }
        
        // Pressure boundary conditions.
        #pragma omp for
        for(int I = 2 ; I <= Nx-1 ; I++){ 
            p[I-1][0] = p[I-1][1] ;
            p[I-1][Ny-1] = p[I-1][Ny-2] ;
        }
        
        #pragma omp for
        for(int J = 2 ; J <= Ny-1 ; J++){
            p[0][J-1] = p[1][J-1] ;
            p[Nx-1][J-1] = p[Nx-2][J-1] ;
        }
    }
    
    p[0][0] = p[1][1];
    p[Nx-1][0] = p[Nx-2][1];
    p[0][Ny-1] = p[1][Ny-2];
    p[Nx-1][Ny-1] = p[Nx-2][Ny-2];
    
    // Parallel velocity correction.
    #pragma omp parallel
    {
        // Correct x-velocity component.
        #pragma omp for collapse(2)
        for(int i = 1; i < Nx-2; i++){
            for(int J = 1; J < Ny-1; J++){
                double dp_dx = (p_prime[i+1][J] - p_prime[i][J])/dx;
                pu_star[i][J] -= dt * dp_dx;
            }
        }
        
        // Correct y-velocity component.
        #pragma omp for collapse(2)
        for(int I = 1; I < Nx-1; I++){
            for(int j = 1; j < Ny-2; j++){
                double dp_dy = (p_prime[I][j+1] - p_prime[I][j])/dy;
                pv_star[I][j] -= dt * dp_dy;
            }
        }
        
        // Reapply boundary conditions.
        #pragma omp for
        for(int i = 0; i < Nx-1; i++){
            pu_star[i][0] = pu_bottom;
            pu_star[i][Ny-1] = pu_up;
        }
        
        #pragma omp for
        for(int j = 0; j < Ny-1; j++){
            pv_star[0][j] = pv_left;
            pv_star[Nx-1][j] = pv_right;
        }
    }
}

void refresh_pressure(int k , int m ){
    #pragma omp parallel for collapse(2)
    for(int I = 1 ; I <= Nx ; I++){
        for(int J = 1 ; J <= Ny ; J++){
            p_star[I-1][J-1] = p[I-1][J-1] ;
        }
    }
}

void output_combined_fields(int m) {
    ostringstream name;
    name << "1.Central Difference_Re = "<< Re << "_Grid"<< Nx << "x" << Ny << setfill('0') << setw(6) << m << ".vtk";
    ofstream out(name.str().c_str());
    
    out << "# vtk DataFile Version 3.0\n";
    out << "SIMPLE_CD_Combined_Fields\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_POINTS\n";
    out << "DIMENSIONS " << Nx << " " << Ny << " 1\n";
    out << "ORIGIN 0 0 0\n";
    out << "SPACING " << dx << " " << dy << " 1\n";
    
    out << "POINT_DATA " << Nx * Ny << "\n";
    out << "SCALARS Pressure double 1\n";
    out << "LOOKUP_TABLE default\n";
    for(int j = 0; j < Ny; j++) {
        for(int i = 0; i < Nx; i++) {
            out << scientific << p[i][j] << "\n";
        }
    }
    
    out << "VECTORS Velocity double\n";
    for(int j = 0; j < Ny; j++) {
        for(int i = 0; i < Nx; i++) {
            double u = 0.0, v = 0.0;
            
            if(i > 0 && i < Nx-1) {
                u = 0.5 * (pu_star[i-1][j] + pu_star[i][j])/rho;
            }
            if(j > 0 && j < Ny-1) {
                v = 0.5 * (pv_star[i][j-1] + pv_star[i][j])/rho;
            }
            
            out << scientific << u << " " << v << " 0.0\n";
        }
    }
    
    out.close();
    cout << "VTK file has been output: " << name.str() << endl;
}

int main() {
    // Configure OpenMP thread count.
    int num_threads = omp_get_max_threads();
    omp_set_num_threads(num_threads);
    
    cout << "========================================" << endl;
    cout << "SIMPLE_CD OpenMP Parallel Version" << endl;
    cout << "Using " << num_threads << " threads" << endl;
    cout << "========================================" << endl;
    cout << "Mesh size : " << Nx << " x " << Ny << endl;
    cout << "Grid space: dx=" << dx << ", dy=" << dy << endl;
    cout << "Boundary conditions (u,v) : " << endl;
    cout << "Left boundary = (0,0) ; Right boundary = (0,0) ; Bottom boundary = (0,0) ; Top boundary = (1,0) " << endl ;
    cout << "Re = " << Re << ", viscosity = " << mu << endl ;
    
    cout << "Allocating memory for pressure correction matrix..." << endl;
    allocate_memory();
    
    CE_physical = false;
    t = 0 ;
    time_init(n) ;
    
    for(int t = 0 ; t <= max_t ; t++){
        if(t%10000 == 0){
            cout << "Time step :" << t << endl;
            cout << "the next time step , mass imbalance = " << get_max_residual(n , t) << endl ;
            output_combined_fields(t) ;
        }
        
        momentum_time(n , t) ;
        
        steadystate = false ;
        G = 0 ;
        comp_itter(a,b,n1) ;
        CD(a,b,n1,t);
        
        for(G = 0 ; G <= max_G ; G++){
            if(t%10000 == 0 && G%5000 == 0){
                cout <<"       " << "Iteration:" << G << endl ;
                cout <<"       " << "Maximum iteration error = " << max_error << endl ;
            }
            Jacobi(a,b,x,n1,t) ;
            if(G>0 && max_error < tolerance){
                steadystate = true ;
                K = G ;
                if(t%10000 == 0){
                    cout << "       " << "Pressure correction convergence achieved!" << endl;
                    cout << "       " << "Iteration:" << K << endl ;
                    cout << "       " << "Maximum iteration error = " << max_error << endl ;
                }
                break ; 
            }
        }
        
        if(!steadystate){
            K = G-1 ;
            cout << "       " << "Maximum number of iterations reached, but pressure correction convergence not achieved!" << endl;
            cout << "       " << "Iteration:" << K << endl ;
            cout << "       " << "Maximum iteration error = " << max_error << endl ;
            T = t ;
            break ;
        }else{
            pressure_time(n , t) ;
            refresh_pressure(n , t) ;   
            
            if(t>0 && get_max_residual(n , t) < CE_tolerance){
                CE_physical = true ;
                T = t ;
                cout << "Continuous equation convergence achieved!" << endl ;
                cout << "Time step:" << T << endl ;
                cout << "the next time step , mass imbalance = " << get_max_residual(n , t) << endl ;
                output_combined_fields(T) ;
                break ;
            }
        }
    }
    
    if(!CE_physical){
        if(!steadystate){
            T = T-1 ;
            cout << "Maximum number of iterations reached, but continuous equation convergence not achieved!" << endl ;
            cout << "Time step:" << T << endl ;
            cout << "the next time step , mass imbalance = " << get_max_residual(n , T) << endl ;
        }else{
            T = t-1 ;
            cout << "Maximum time step reached, but continuous equation convergence not achieved!" << endl ;
            cout << "Time step:" << T << endl ;
            cout << "the next time step , mass imbalance = " << get_max_residual(n , T) << endl ;
        }
    }
    
    cout << "Grid size " << Nx << "x" << Ny << " computation completed" << endl;
    
    deallocate_memory();
    
    return 0 ;
}