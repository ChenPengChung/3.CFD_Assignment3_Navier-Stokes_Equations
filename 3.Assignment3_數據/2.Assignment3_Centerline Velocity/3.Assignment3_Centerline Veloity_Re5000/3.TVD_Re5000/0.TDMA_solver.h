#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <vector> 
#pragma once
using namespace std ; 
//writting TDMA solver
inline void TDMA(
    const vector<double>& D,       //primary delements 
    const vector<double>& beta2,   //left hand side elemets 
    const vector<double>& alpha2,  //right hand side elememts 
    vector<double>& x,             // solutions 
    const vector<double>& C,       // sourceterm 
    int n,
    double U_left, double U_right) // left boundary and right boundary condition 
{
    if (n <= 0) return;
    if ((int)D.size()!=n || (int)beta2.size()!=n || (int)alpha2.size()!=n ||
        (int)C.size()!=n || (int)x.size()!=n) {
        throw std::invalid_argument("All input vector sizes must equal n.");
    }

    vector<double> A(n, 0.0);   
    vector<double> Cp(n, 0.0);  

    // forward elimination 
    double denom = D[0];
    if (denom == 0.0) throw std::runtime_error("Singular pivot at row 0.");
    A[0]  = (n > 1) ? (alpha2[0] / denom) : 0.0;
    Cp[0] = (C[0] - beta2[0]*U_left) / denom;

    
    for (int i = 1; i < n-1; ++i) {
        denom = D[i] - beta2[i]*A[i-1];
        if (denom == 0.0) throw std::runtime_error("Singular pivot at row " + std::to_string(i));
        A[i]  = alpha2[i] / denom;
        Cp[i] = (C[i] - beta2[i]*Cp[i-1]) / denom;          
    }

    //backward substituition 
    if (n >= 2) {
        int i = n-1;
        denom = D[i] - beta2[i]*A[i-1];
        if (denom == 0.0) throw std::runtime_error("Singular pivot at last row.");
        Cp[i] = (C[i] - beta2[i]*Cp[i-1] - alpha2[i]*U_right) / denom; 
    }

    
    x[n-1] = Cp[n-1];
    for (int i = n-2; i >= 0; --i) {
        x[i] = Cp[i] - A[i]*x[i+1];
    }
}


inline void TDMA_FVM(
    const vector<double>& D,
    const vector<double>& beta2,
    const vector<double>& alpha2,
    vector<double>& x,
    const vector<double>& C,
    int n)
{
    if (n <= 0) return;
    if ((int)D.size()!=n || (int)beta2.size()!=n || (int)alpha2.size()!=n ||
        (int)C.size()!=n || (int)x.size()!=n) {
        throw std::invalid_argument("All input vector sizes must equal n.");
    }

    vector<double> A(n, 0.0), Cp(n, 0.0);

    // forward i=0
    double denom = D[0];
    if (denom == 0.0) throw std::runtime_error("Singular pivot at row 0.");
    A[0]  = (n>1) ? (alpha2[0] / denom) : 0.0;
    Cp[0] = C[0] / denom;

    // forward i=1..n-1
    for (int i = 1; i < n; ++i) {
        denom = D[i] - beta2[i]*A[i-1];
        if (denom == 0.0) throw std::runtime_error("Singular pivot at row " + std::to_string(i));
        A[i]  = (i < n-1) ? (alpha2[i] / denom) : 0.0;     
        Cp[i] = (C[i] - beta2[i]*Cp[i-1]) / denom;        
    }

    // backward
    x[n-1] = Cp[n-1];
    for (int i = n-2; i >= 0; --i) {
        x[i] = Cp[i] - A[i]*x[i+1];
    }
}

//註記：TDMA(D,beta,alpha,x,C,n,U_left,Uright)  , TDMA_FVM(D,beta,alpha,x,C,n)
//D  ; primary elements 
//beta ; left hand side elements 
//alpha ; right hand side elements 
//x ; solutions 
//C ; source term 
//n ; number of the elements of the vector
