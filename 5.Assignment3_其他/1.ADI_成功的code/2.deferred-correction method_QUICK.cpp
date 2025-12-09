//利用有限體積法配合SIMPLE算法求解Navier-Stokes方程組
//加入deferred correction 並使用交替迭代
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <vector> 
#include <omp.h>  // 加在 #include <vector> 之後 
#include "0.P_correction_eq_TDMA.h" 
using namespace std ;
//參數設定
//會用到的參數
//壓力場在ｘ方向有Nx格
//x分量速度場在x方向有Nx+1格
//y分量速度場在y方向有Ny+1格
int G ;  //G為當前計算步數從0開始，與迭代步數意義不同
const int  max_G = 700000;//最大計算步數
double min_step7 ; 
//x分量速度場 : i \int [0,Nx] , j \int [0,Ny-1]
//y分量速度場 : i \int [0,Nx-1] , j \int [0,Ny]
//壓力場 : i \int [0,Nx-1] , j \int [0,Ny-1]
double u_star[Nx+1][Ny] ; //x分量速度場預估值
double u_star_old[Nx+1][Ny] ;
double v_star[Nx][Ny+1] ; //y分量速度場預估值
double v_star_old[Nx][Ny+1] ;
double p_star[Nx][Ny] ; //壓力場預估值
double u_prime[Nx+1][Ny] ; //x分量速度場修正值
double v_prime[Nx][Ny+1] ; //y分量速度場修正值
double p_prime[Nx][Ny] ; //壓力場修正值
double u[Nx+1][Ny] ; //x分量速度場正確值
//創建舊x分量速度場正確值 矩陣
double u_old[Nx+1][Ny] ;//舊x分量速度場正確值需要進行初始化
double v[Nx][Ny+1] ; //y分量速度場正確值
//創建舊y分量速度場正確值 矩陣
double v_old[Nx][Ny+1] ;//舊y分量速度場正確值需要進行初始化
double p[Nx][Ny] ; //壓力場正確值

const double  Re = 5000.0 ;
//x分量速度場邊界條件
const double u_left  = 0.0 ;//壁面不可穿透條件
const double u_right = 0.0 ;//壁面不可穿透條件
const double u_bottom = 0.0 ;//無滑移條件
const double u_top = 1.0 ;//無滑移條件
//y分量速度場邊界條件
const double v_left  = 0.0 ;//無滑移條件
const double v_right = 0.0 ;//無滑移條件
const double v_bottom = 0.0 ;//壁面不可穿透條件
const double v_top = 0.0 ;//壁面不可穿透條件
const double omega_P = 0.6 ; //根據Pateanca的做法，鬆弛因子我應該設pressure 的relation of predicted-correction method :　　ａｌｐｈａ＿Ｐ　＝　０．８
double omega_u  = 0.2 ; 
double omega_v  = 0.2 ;
double max_error ;
const float CE_tolerance = 1e-10; //收斂判據：預估值的連續方程式
bool CE_physical ;

//離散方程的對流部分係數(x,y兩分量可共用)
const double D_e = (1.0/Re);
const double D_w = (1.0/Re);
const double D_n = (1.0/Re);
const double D_s = (1.0/Re);
double r_w_u[Nx+1][Ny] ; //x分量速度場介面判斷係數
double r_e_u[Nx+1][Ny] ;
double r_s_u[Nx+1][Ny] ;
double r_n_u[Nx+1][Ny] ;
double r_w_v[Nx][Ny+1] ; //y分量速度場介面判斷係數
double r_e_v[Nx][Ny+1] ;
double r_s_v[Nx][Ny+1] ;
double r_n_v[Nx][Ny+1] ;
//x分量動量積分方程的一般離散格式的對流部分係數//相應的因子\alpha需要由Ｆ判斷
double F_w_u[Nx+1][Ny] ;
double alpha_w_u[Nx+1][Ny] ;
double F_e_u[Nx+1][Ny] ;
double alpha_e_u[Nx+1][Ny] ;
double F_s_u[Nx+1][Ny] ;
double alpha_s_u[Nx+1][Ny] ;
double F_n_u[Nx+1][Ny] ;
double alpha_n_u[Nx+1][Ny] ;
//x分量動量積分方程的一般離散格式的線性組合係數
double a_P[Nx+1][Ny] ;
double a_WW[Nx+1][Ny] ;
double a_W[Nx+1][Ny] ;
double a_E[Nx+1][Ny] ;
double a_EE[Nx+1][Ny] ;
double a_SS[Nx+1][Ny] ;
double a_S[Nx+1][Ny] ;
double a_N[Nx+1][Ny] ;
double a_NN[Nx+1][Ny] ; 
double S1[Nx+1][Ny] ; //x分量動量方程的非齊性項
//y分量動量積分方程的一般離散格式的對流部分係數//相應的因子\alpha需要由Ｆ判斷
double F_w_v[Nx][Ny+1] ;
double alpha_w_v[Nx][Ny+1] ;
double F_e_v[Nx][Ny+1] ;
double alpha_e_v[Nx][Ny+1] ;
double F_s_v[Nx][Ny+1] ;
double alpha_s_v[Nx][Ny+1] ;
double F_n_v[Nx][Ny+1] ;
double alpha_n_v[Nx][Ny+1] ;
//y分量動量積分方程的一般離散格式的線性組合係數
double b_P[Nx][Ny+1] ;
double b_WW[Nx][Ny+1] ;
double b_W[Nx][Ny+1] ;
double b_E[Nx][Ny+1] ;
double b_EE[Nx][Ny+1] ;
double b_SS[Nx][Ny+1] ;
double b_S[Nx][Ny+1] ;
double b_N[Nx][Ny+1] ;
double b_NN[Nx][Ny+1] ;
double S2[Nx][Ny+1] ;//y分量動量方程的非齊性項
//以上矩陣記得更新，一但速度場預估值更新完，就要馬上更新

/*//TVD-MUIDSCL 上游下游迭代係數
double D(double x) {
    if (x <= 0) return 0.0;
    else if (x > 0 && x < (1.0/3.0)) return 0.0;
    else if (x >= (1.0/3.0) && x < 3.0) return 0.25;
    else return 1.0;
}

double U(double x) {
    if (x <= 0) return 1.0;
    else if (x > 0 && x < (1.0/3.0)) return 1.0;
    else if (x >= (1.0/3.0) && x < 3.0) return 0.0;
    else return 0.0;
}

double UU(double x) {
    if (x <= 0) return 0.0;
    else if (x > 0 && x < (1.0/3.0)) return -1.0;
    else if (x >= (1.0/3.0) && x < 3.0) return -0.25;
    else return 0.0;
}*/



//QUICK 上游下游迭代係數
double D(double x){
	return 3.0/8.0 ;
}
double U(double x){
    return 6.0/8.0 ;
}
double UU(double x){
	return -1.0/8.0 ;
}

//SIMPLE算法：
//step1：x分量速度場預估值初始化、y分量速度場預估值初始化、壓力場預估值初始化
void initialize(){
    #pragma omp parallel  // *** 需要加這行! ***
    {
    //x分量速度場預估值初始化
     #pragma omp parallel for collapse(2) 
    for(int i = 1 ; i <= Nx+1 ; i++){
        //x分量速度場預估值的定義域為i \int [0,Nx] , j \int [0,Ny-1]
        for(int J= 1 ; J <= Ny ; J++){
            u_star[i-1][J-1] = 0.0 ;
            if(i == 1) u_star[i-1][J-1] = u_left ; //壁面不可穿透條件
            if(i == Nx + 1) u_star[i-1][J-1] = u_right ; //壁面不可穿透條件
        }
    }
    //舊x分量速度場正確值初始化->由一開始的預估值賦值，後續更新由真正的正確值賦值
     #pragma omp parallel for collapse(2) 
    for(int i = 1 ; i <= Nx+1 ; i++){
        for(int J = 1 ; J <= Ny ; J++){
            u_old[i-1][J-1] = u_star[i-1][J-1] ;
        }
    }
    //y分量速度場預估值初始化
     #pragma omp parallel for collapse(2) 
    for(int I = 1 ; I <= Nx ; I++){
        //y分量速度場預估值的定義譽為I \int [0,Nx-1] , j \int [0,Ny]
        for(int j = 1 ; j <= Ny+1 ; j++){
            v_star[I-1][j-1] = 0.0 ;
            if(j == 1)v_star[I-1][j-1] = v_bottom ;//壁面不可穿透條件
            if(j == Ny + 1)v_star[I-1][j-1] = v_top ;//壁面不可穿透條件
        }
    }
    //舊y分量速度場正確值初始化->由一開始的預估值賦值，後續更新由真正的正確值賦值
     #pragma omp parallel for collapse(2) 
    for(int I = 1 ; I <= Nx ; I++){
        for(int j = 1 ; j <= Ny+1 ; j++){
            v_old[I-1][j-1] = v_star[I-1][j-1] ;
        }
    }
    //壓力場預估值初始化
     #pragma omp parallel for collapse(2) 
    for(int I = 1 ; I <= Nx ; I++){
        for(int J = 1 ; J<= Ny ; J++){
            p_star[I-1][J-1] =  0.0  ; //邊界條件為Neumann條件，由迭代方程實現
            p_prime[I-1][J-1] = 0.0 ;
        }
    }
    //預估值方程會用到的矩陣，分別為介面判斷係數矩陣->對流部分係數矩陣->線性組合係數矩陣
    //更新，或者初始化順序如下，介面判斷係數矩陣->對流部分係數矩陣->線性組合係數矩陣
    //初始化介面判斷係數
    //x分量速度場介面判斷係數(北邊界判斷係數，南邊界判斷係數需要特殊討論)
     #pragma omp parallel for collapse(2) 
    for(int i = 1 ; i <= Nx+1 ; i++){
        for(int J = 1 ; J <= Ny ; J++){
            //不需要討論的cell->介面判斷係數直接設為0
            if(i == 1 || i == Nx+1){
                r_w_u[i-1][J-1] = 0.0 ;
                r_e_u[i-1][J-1] = 0.0 ;
                r_n_u[i-1][J-1] = 0.0 ;
                r_s_u[i-1][J-1] = 0.0 ;
            }else {
                //按照各個邊界討論（西邊界特殊討論＊1東邊界特輸討論*南邊界特殊討論*北邊界特殊討論*3)
                //西邊界
                if(i == 2){
                    r_w_u[i-1][J-1] = -1.0 ; //西邊界採用迎風格式
                }else{
                    r_w_u[i-1][J-1] = 1.0 ; //初始化設計
                }
                //東邊界
                if(i == Nx){
                    r_e_u[i-1][J-1] = -1.0 ; //東邊界採用迎風格式
                }else{
                    r_e_u[i-1][J-1] = 1.0 ;//初始化設計
                }
                //南邊界
                if(J == 1){
                    r_s_u[i-1][J-1] = 0.0 ; //下邊界計算點用不到南邊界判斷因子
                }else if(J == 2){
                    r_s_u[i-1][J-1] = -1.0 ;//迎風格式//需要特殊討論//避免除以零風險
                }else if(J == Ny){
                    r_s_u[i-1][J-1] = -1.0 ;//迎風格式//需要特殊討論//避免除以零風險
                }else{
                    r_s_u[i-1][J-1] = 1.0 ;//初始化設計
                }
                //北邊界
                if(J == 1){
                    r_n_u[i-1][J-1] = -1.0 ; //迎風格式//需要特殊討論//避免除以零風險
                }else if (J == Ny-1){
                    r_n_u[i-1][J-1] = -1.0 ;//迎風格式//需要特殊討論//避免除以零風險
                }else if(J == Ny){
                    r_n_u[i-1][J-1] = 0.0 ; //上邊界計算點用不到北邊界判斷因子
                }else{
                    r_n_u[i-1][J-1] = 1.0 ;//初始化設計
                }
            }
        }    
    }//記得預估值更新後，立即更新介面判斷係數
    //y分量速度場介面判斷係數(東邊界判斷細數，西邊界判斷係數需要特殊討論)
     #pragma omp parallel for collapse(2) 
    for(int I = 1 ; I <= Nx ; I++){
        for(int j = 1 ; j <= Ny+1 ; j++){
            //不需要討論的cell->介面判斷係數直接設為0
            if(j == 1 || j == Ny+1){
                r_w_v[I-1][j-1] = 0.0 ;
                r_e_v[I-1][j-1] = 0.0 ;
                r_n_v[I-1][j-1] = 0.0 ;
                r_s_v[I-1][j-1] = 0.0 ;
            }else {
                //按照各個邊界依序討論（南邊界特殊討論*1北邊界特殊討論*1東邊界特殊討論*3西邊界特殊討論*3)
                //南邊界
                if(j == 2){
                    r_s_v[I-1][j-1] = -1.0 ; //下邊界計算點的南邊界採用迎風格式
                }else{
                    r_s_v[I-1][j-1] = 1.0 ; //初始化設計
                }
                //北邊界
                if(j == Ny){
                    r_n_v[I-1][j-1] = -1.0 ; //上邊界計算點的北邊界採用迎風格式
                }else{
                    r_n_v[I-1][j-1] = 1.0 ; //初始化設計
                }
                //西邊界
                if(I == 1){
                    r_w_v[I-1][j-1] = 0.0 ; //西邊界計算點用不到西邊界判斷因子
                }else if (I == 2){
                    r_w_v[I-1][j-1] = -1.0 ;//迎風格式//需要特殊討論//避免除以零風險
                }else if(I == Nx){
                    r_w_v[I-1][j-1] = -1.0 ;//迎風格式//需要特殊討論//避免除以零風險
                }else{
                    r_w_v[I-1][j-1] = 1.0 ;//初始化設計
                }
                //東邊界
                if(I == 1){
                    r_e_v[I-1][j-1] = -1.0 ; //迎風格式//需要特殊討論//避免除以零風險
                }else if (I == Nx-1){
                    r_e_v[I-1][j-1] = -1.0 ; //迎風格式//需要特殊討論//避免除以零風險
                }else if(I == Nx){
                    r_e_v[I-1][j-1] = 0.0 ; //東邊界計算點用不到東邊界判斷因子
                }else {
                    r_e_v[I-1][j-1] = 1.0 ; //初始化設計
                }
            }
        }
    }
    //初始化對流部分係數矩陣
    //x分量動量積分方程的一般離散格式的對流部分係數
     #pragma omp parallel for collapse(2) 
    for(int i = 1 ; i <= Nx+1 ; i++){
        for(int J = 1 ; J <= Ny ; J++){
            if(i == 1 || i == Nx+1) {
                F_w_u[i-1][J-1] = 0.0 ;
                F_e_u[i-1][J-1] = 0.0 ;
                F_n_u[i-1][J-1] = 0.0 ;
                F_s_u[i-1][J-1] = 0.0 ;
            }else{
                F_e_u[i-1][J-1] = 0.5*dy*(u_star[i-1][J-1]+u_star[i][J-1]); //網格中心點剛好對應到變量儲存點
                F_w_u[i-1][J-1] = 0.5*dy*(u_star[i-1][J-1]+u_star[i-2][J-1]);
                F_n_u[i-1][J-1] = 0.5*dx*(v_star[i-1][J]+v_star[i-2][J]);
                F_s_u[i-1][J-1] = 0.5*dx*(v_star[i-1][J-1]+v_star[i-2][J-1]);
            }
        }
    }//alpha_f_u判斷因子，下標代表控制體積的中心點座標
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
    //y分量動量積分方程的一般離散格式的對流部分係數
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
                F_n_v[I-1][j-1] = 0.5*dx*(v_star[I-1][j-1]+v_star[I-1][j]);//網格中心點剛好對應到變量儲存點
                F_s_v[I-1][j-1] = 0.5*dx*(v_star[I-1][j-1]+v_star[I-1][j-2]);
            }
        }
    }//alpha_f_v判斷因子，其下標代表控制體積中心點座標
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

//step2：x分量速度場預估值更新、y分量速度場預估值更新
void step2(int m){
    //更新x分量動量積分方程的一般離散格式的預估值形式的線性組合係數
    #pragma omp parallel for collapse(2)
    for(int i = 1 ; i <= Nx+1 ; i++){
        for(int J = 1 ; J <= Ny ; J++){
            //各個鄰計算點的線組係數先賦值為0，再來做疊加
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
                S1[i-1][J-1] += u_left ;//物面不可穿透條件//左邊界計算點沒有壓力梯度
                u_star[i-1][J-1] = u_left ;
            }else if(i == Nx+1){
                a_P[i-1][J-1] += 1.0 ;
                S1[i-1][J-1] += u_right ;//物面不可穿透條件//左邊界計算點沒有壓力梯度
                u_star[i-1][J-1] = u_right ;
            }//誰需要特殊討論，左邊界第二排，右邊界第二排，下邊界第一排＆第二排，上邊界第一排＆第二排
            //左邊第二排
            else{ //排除上二種特殊情形後，下面的cell分西邊界，東邊界，南邊界，北邊界，來討論
                //西邊界
                if(i == 2){//左邊界第二排計算點的西側邊界採用迎風格式
                    a_W[i-1][J-1] += alpha_w_u[i-1][J-1]*(-F_w_u[i-1][J-1])-D_w;//討論西逼界時，再減去西邊界擴散部分係數
                    a_P[i-1][J-1] += (1-alpha_w_u[i-1][J-1])*(-F_w_u[i-1][J-1])+D_w;
                    S1[i-1][J-1] += -(-p_star[i-2][J-1])*dy/rho ;
                }else{//西邊界迎風項為w,p 
                    //a_WW[i-1][J-1] += (alpha_w_u[i-1][J-1]*UU(r_w_u[i-1][J-1]))*(-F_w_u[i-1][J-1]) ;
                    a_W[i-1][J-1] += (alpha_w_u[i-1][J-1]*1)*(-F_w_u[i-1][J-1]) - D_w ; //討論西邊界時，再減去西側擃散部分係數
                    //a_E[i-1][J-1] += ((1-alpha_w_u[i-1][J-1])*UU(r_w_u[i-1][J-1]))*(-F_w_u[i-1][J-1]);
                    a_P[i-1][J-1] += ((1-alpha_w_u[i-1][J-1])*1)*(-F_w_u[i-1][J-1])+D_w ;
                    S1[i-1][J-1] += -(-p_star[i-2][J-1])*dy/rho
                    -u_star_old[i-3][J-1]*(alpha_w_u[i-1][J-1]*UU(r_w_u[i-1][J-1]))*(-F_w_u[i-1][J-1])
                    -u_star_old[i-2][J-1]*((1-alpha_w_u[i-1][J-1])*D(r_w_u[i-1][J-1])+alpha_w_u[i-1][J-1]*(U(r_w_u[i-1][J-1])-1))*(-F_w_u[i-1][J-1])
                    -u_star_old[i][J-1]*((1-alpha_w_u[i-1][J-1])*UU(r_w_u[i-1][J-1]))*(-F_w_u[i-1][J-1])
                    -u_star_old[i-1][J-1]*((1-alpha_w_u[i-1][J-1])*(U(r_w_u[i-1][J-1])-1)+alpha_w_u[i-1][J-1]*D(r_w_u[i-1][J-1]))*(-F_w_u[i-1][J-1]) 
                    ; 
                }//deferred correction 完成!!
                //東邊界
                if(i == Nx){//右邊界第二排計算點的東側邊界採用迎風格式
                    a_E[i-1][J-1] += (1-alpha_e_u[i-1][J-1])*(F_e_u[i-1][J-1])-D_e ;//討論東邊界時，再扣掉東邊界擴散部分係數
                    a_P[i-1][J-1] += alpha_e_u[i-1][J-1]*(F_e_u[i-1][J-1])+D_e; 
                    S1[i-1][J-1] += -(p_star[i-1][J-1])*dy/rho ;
                }else{//東邊界迎風項牽涉東點及本點
                    //a_EE[i-1][J-1] += (1-alpha_e_u[i-1][J-1])*UU(r_e_u[i-1][J-1])*(F_e_u[i-1][J-1]) ;
                    a_E[i-1][J-1] += ((1-alpha_e_u[i-1][J-1])*1)*(F_e_u[i-1][J-1]) - D_e ;//討論東邊界時，再扣掉東邊界擴散部分係數
                    //a_W[i-1][J-1] += (alpha_e_u[i-1][J-1])*UU(r_e_u[i-1][J-1])*(F_e_u[i-1][J-1]) ;
                    a_P[i-1][J-1] += (alpha_e_u[i-1][J-1]*1)*(F_e_u[i-1][J-1])+D_e;   
                    S1[i-1][J-1] += -(p_star[i-1][J-1])*dy/rho 
                    -u_star_old[i+1][J-1]*(1-alpha_e_u[i-1][J-1])*UU(r_e_u[i-1][J-1])*(F_e_u[i-1][J-1])
                    -u_star_old[i][J-1]*(alpha_e_u[i-1][J-1]*D(r_e_u[i-1][J-1])+(1-alpha_e_u[i-1][J-1])*(U(r_e_u[i-1][J-1])-1))*(F_e_u[i-1][J-1])
                    -u_star_old[i-2][J-1]*(alpha_e_u[i-1][J-1])*UU(r_e_u[i-1][J-1])*(F_e_u[i-1][J-1])
                    -u_star_old[i-1][J-1]*((1-alpha_e_u[i-1][J-1])*D(r_e_u[i-1][J-1])+alpha_e_u[i-1][J-1]*(U(r_e_u[i-1][J-1])-1))*(F_e_u[i-1][J-1])
                    ;   
                }//deferred corrction 完成!!
                //南邊界(要討論三個)
                if(J == 1){//下邊界第一件計算點
                    a_P[i-1][J-1] += 2*D_s ;
                    S1[i-1][J-1] += -(u_bottom)*(-F_s_u[i-1][J-1]) +2*D_s*(u_bottom) ;//無滑移條件//南邊界中心點上的壓力取值為本點壓力，根據南邊界上的條件Neumann條件啟示，將本點上的壓力插值到南邊界上的壓力
                }else if(J == 2){//下邊界第二排計算點的南側邊界使用外插法取代南南點上的取//南邊界迎風項牽涉南點以及本點
                    a_S[i-1][J-1] += ((alpha_s_u[i-1][J-1])*1)*(-F_s_u[i-1][J-1])-D_s;
                   //a_N[i-1][J-1] += ((1-alpha_s_u[i-1][J-1])*UU(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1]) ;
                    a_P[i-1][J-1] += ((1-alpha_s_u[i-1][J-1])*1)*(-F_s_u[i-1][J-1]) + D_s ;
                    S1[i-1][J-1] +=  -2*(alpha_s_u[i-1][J-1]*UU(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1])*u_bottom
                    +(alpha_s_u[i-1][J-1]*UU(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1])*u_star_old[i-1][J-2]
                    -u_star_old[i-1][J-2]*((1-alpha_s_u[i-1][J-1])*D(r_s_u[i-1][J-1])+(alpha_s_u[i-1][J-1])*(U(r_s_u[i-1][J-1])-1))*(-F_s_u[i-1][J-1])
                    -u_star_old[i-1][J]*((1-alpha_s_u[i-1][J-1])*UU(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1])
                    -u_star_old[i-1][J-1]*((1-alpha_s_u[i-1][J-1])*(U(r_s_u[i-1][J-1])-1)+(alpha_s_u[i-1][J-1])*D(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1])
                    ;//deferred correction 完成!!
                }else if(J == Ny){//南邊界迎風項牽涉南點以及本點
                    //北邊界計算點不存在，使用外插法，插值過去
                    //a_SS[i-1][J-1]+= (alpha_s_u[i-1][J-1]*UU(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1]) ;
                    a_S[i-1][J-1] += ((alpha_s_u[i-1][J-1])*1)*(-F_s_u[i-1][J-1])-D_s ;
                    a_P[i-1][J-1]+= ((1-alpha_s_u[i-1][J-1])*1)*(-F_s_u[i-1][J-1])+D_s;
                    S1[i-1][J-1] += -2*u_top*((1-alpha_s_u[i-1][J-1])*UU(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1])
                    +u_star_old[i-1][J-1]*((1-alpha_s_u[i-1][J-1])*UU(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1])
                    -u_star_old[i-1][J-3]*(alpha_s_u[i-1][J-1]*UU(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1])
                    -u_star_old[i-1][J-2]*((1-alpha_s_u[i-1][J-1])*D(r_s_u[i-1][J-1])+(alpha_s_u[i-1][J-1])*(U(r_s_u[i-1][J-1])-1))*(-F_s_u[i-1][J-1])
                    -u_star_old[i-1][J-1]*((1-alpha_s_u[i-1][J-1])*(U(r_s_u[i-1][J-1])-1)+alpha_s_u[i-1][J-1]*D(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1]) ; //deferred correction 完成 !!
                }else{//s南邊界迎風項牽涉南點以及本]點
                    //a_SS[i-1][J-1] += (alpha_s_u[i-1][J-1]*UU(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1]) ;
                    a_S[i-1][J-1] += ((alpha_s_u[i-1][J-1])*1)*(-F_s_u[i-1][J-1])-D_s ;//討論南邊界時，再把南邊界的擴散部分係數扣掉
                    //a_N[i-1][J-1] += ((1-alpha_s_u[i-1][J-1])*UU(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1]) ;
                    a_P[i-1][J-1] += ((1-alpha_s_u[i-1][J-1])*1)*(-F_s_u[i-1][J-1])+D_s ; 
                    S1[i-1][J-1] += 
                    -u_star[i-1][J-3]*(alpha_s_u[i-1][J-1]*UU(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1]) 
                    -u_star[i-1][J-2]*((1-alpha_s_u[i-1][J-1])*D(r_s_u[i-1][J-1])+(alpha_s_u[i-1][J-1])*(U(r_s_u[i-1][J-1])-1))*(-F_s_u[i-1][J-1])
                    -u_star_old[i-1][J]*((1-alpha_s_u[i-1][J-1])*UU(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1])
                    -u_star_old[i-1][J-1]*((1-alpha_s_u[i-1][J-1])*(U(r_s_u[i-1][J-1])-1)+alpha_s_u[i-1][J-1]*D(r_s_u[i-1][J-1]))*(-F_s_u[i-1][J-1])
                    ;//deferred correction 完成
                }
                //北邊界分成J == 1 && J == Ny && J == Ny-1 三種特殊討論
                if(J == Ny){
                    a_P[i-1][J-1] += 2*D_n ;
                    S1[i-1][J-1] += -(u_top)*F_n_u[i-1][Ny-1] + 2*D_n*(u_top);//無滑移條件//北邊界中心點上的壓力取值為本點壓力，根據北邊界上的條件Neumann條件啟示，將本點上的壓力插值到北邊界上的壓力
                }else if(J == 1){
                    //a_NN[i-1][J-1] += (1-alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1]) ;
                    a_N[i-1][J-1] += ((1-alpha_n_u[i-1][J-1])*1)*(F_n_u[i-1][J-1]) - D_n ;//討論北邊界時，再把北邊界的擴散部分係數扣掉
                    a_P[i-1][J-1] += ((alpha_n_u[i-1][J-1])*1)*(F_n_u[i-1][J-1]) + D_n ;   
                    S1[i-1][J-1] += -2*(alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1])*u_bottom
                    +u_star_old[i-1][J-1]*(alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1]) 
                    -u_star_old[i-1][J+1]*(1-alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1])
                    -u_star_old[i-1][J]*(alpha_n_u[i-1][J-1]*D(r_n_u[i-1][J-1])+(1-alpha_n_u[i-1][J-1])*(U(r_n_u[i-1][J-1])-1))*(F_n_u[i-1][J-1]) 
                    -u_star_old[i-1][J-1]*((1-alpha_n_u[i-1][J-1])*D(r_n_u[i-1][J-1])+(alpha_n_u[i-1][J-1])*(U(r_n_u[i-1][J-1])-1))*(F_n_u[i-1][J-1])   
                    ;//deferred correction 完成
                }else if(J == Ny-1){
                    //北邊界計算點不存在，使用外插法，插值過去
                    //a_NN[i-1][J-1] += (1-alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1]) ;
                    a_N[i-1][J-1] += (1-alpha_n_u[i-1][J-1])*1*(F_n_u[i-1][J-1]) - D_n  ;//討論北邊界時，再把北邊界的擴散部分係數扣掉
                    //a_S[i-1][J-1] += (alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1]) ;
                    a_P[i-1][J-1] += (alpha_n_u[i-1][J-1]*1)*(F_n_u[i-1][J-1]) +D_n ;
                    S1[i-1][J-1] += -2*u_top*(1-alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1]) 
                    +u_star_old[i-1][J]*(1-alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1])
                    -u_star_old[i-1][J]*(alpha_n_u[i-1][J-1]*D(r_n_u[i-1][J-1])+(1-alpha_n_u[i-1][J-1])*(U(r_n_u[i-1][J-1])-1))*(F_n_u[i-1][J-1]) 
                    -u_star_old[i-1][J-2]*(alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1]) 
                    -u_star_old[i-1][J-1]*((1-alpha_n_u[i-1][J-1])*D(r_n_u[i-1][J-1])+alpha_n_u[i-1][J-1]*(U(r_n_u[i-1][J-1])-1))*(F_n_u[i-1][J-1]) 
                    ;//deferred correction 完成 
                    }else{
                    //a_NN[i-1][J-1] += (1-alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1]) ;
                    a_N[i-1][J-1] += ((1-alpha_n_u[i-1][J-1])*1)*(F_n_u[i-1][J-1]) - D_n ;//討論北邊界時，再把北邊界的擴散部分係數扣掉
                    //a_S[i-1][J-1] += (alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1]) ;
                    a_P[i-1][J-1] += (alpha_n_u[i-1][J-1]*1)*(F_n_u[i-1][J-1])+D_n ;
                    S1[i-1][J-1] += 
                    -u_star_old[i-1][J+1]*(1-alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1]) 
                    -u_star_old[i-1][J]*(alpha_n_u[i-1][J-1]*D(r_n_u[i-1][J-1])+(1-alpha_n_u[i-1][J-1])*(U(r_n_u[i-1][J-1])-1))*(F_n_u[i-1][J-1]) 
                    -u_star_old[i-1][J-2]*(alpha_n_u[i-1][J-1])*UU(r_n_u[i-1][J-1])*(F_n_u[i-1][J-1])
                    -u_star_old[i-1][J-1]*((1-alpha_n_u[i-1][J-1])*D(r_n_u[i-1][J-1])+alpha_n_u[i-1][J-1]*(U(r_n_u[i-1][J-1])-1))*(F_n_u[i-1][J-1])
                    ;
                 }
                double WW = 0.0 , W = 0.0 , E = 0.0 , EE = 0.0 ;
                double SS = 0.0 , S = 0.0 , N = 0.0 , NN = 0.0 ;
                E = u_star[i][J-1] ;
                W = u_star[i-2][J-1] ;
                //WW:
                if(i == 2){
                    WW = 0.0 ;
                }else{
                    WW = u_star[i-3][J-1] ;
                }
                //EE:
                if(i == Nx){
                    EE = 0.0 ;
                }else{
                    EE = u_star[i+1][J-1] ;
                }
                //SS:
                if(J == 1 || J == 2){
                    SS = 0.0 ;
                }else{
                    SS = u_star[i-1][J-3] ;
                }
                //S :
                if(J == 1){
                    S = 0.0 ;
                }else{
                    S = u_star[i-1][J-2] ;
                }
                //N :
                if(J == Ny){
                    N = 0.0 ;
                }else{
                    N = u_star[i-1][J] ;
                }
                //NN:
                if(J == Ny || J == Ny-1){
                    NN = 0.0 ;
                }else {
                    NN = u_star[i-1][J+1] ;
                }
                u_star[i-1][J-1] = (-a_WW[i-1][J-1]*WW - a_W[i-1][J-1]*W - a_E[i-1][J-1]*E - a_EE[i-1][J-1]*EE
                    - a_SS[i-1][J-1]*SS - a_S[i-1][J-1]*S - a_N[i-1][J-1]*N - a_NN[i-1][J-1]*NN
                    + S1[i-1][J-1]) / a_P[i-1][J-1];
                //速度場預估值的鬆弛迭代關係式
                u_star[i-1][J-1] = (1-omega_u)*u_star_old[i-1][J-1] + omega_u*u_star[i-1][J-1] ;
            }
        }
    }
   
    //更新y分量動量積分方程的一般離散格式的預估值形式的線性組合係數
    #pragma omp parallel for collapse(2)
    for(int I = 1 ; I <= Nx ; I++){
        for(int j = 1 ; j <= Ny+1 ; j++){
            //各個鄰計算點的線組係數先賦值為0，再來做疊加
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
                S2[I-1][j-1] += v_bottom ;//物面不可穿透條件//下邊界計算點沒有壓力梯度
                v_star[I-1][j-1] = v_bottom ;
            }else if(j == Ny+1){
                b_P[I-1][j-1] += 1.0 ;
                S2[I-1][j-1] += v_top ;//物面不可穿透條件//上邊界計算點沒有壓力梯度
                v_star[I-1][j-1] = v_top ;
            }//誰需要特殊討論，左邊界第一排＆第二排，右邊界第一排＆第二排，上邊界第二排，下邊界第二排
            //下邊界第二排
            else{ //排除上二種特殊情形後，下面的cell分西邊界，東邊界，南邊界，北邊界，來討論
                //南邊界
                if(j == 2){//下邊界第二排計算點的南側邊界採用迎風格式
                    b_S[I-1][j-1] += alpha_s_v[I-1][j-1]*(-F_s_v[I-1][j-1])-D_s;//討論南邊界時，再減去南邊界擴散部分係數
                    b_P[I-1][j-1] += (1-alpha_s_v[I-1][j-1])*(-F_s_v[I-1][j-1])+D_s;
                    S2[I-1][j-1] += -(-p_star[I-1][j-2])*dx/rho ;
                }else{//其他點南側邊界採用TVD-MUSCL格式//營風格是牽涉南點以及本點
                    //b_SS[I-1][j-1] += (alpha_s_v[I-1][j-1]*UU(r_s_v[I-1][j-1]))*(-F_s_v[I-1][j-1]) ;
                    b_S[I-1][j-1] += (alpha_s_v[I-1][j-1]*1)*(-F_s_v[I-1][j-1]) - D_s ; //討論南邊界時，再減去南邊界擴散部分係數
                    //b_N[I-1][j-1] += ((1-alpha_s_v[I-1][j-1])*UU(r_s_v[I-1][j-1]))*(-F_s_v[I-1][j-1]) ;
                    b_P[I-1][j-1] += ((1-alpha_s_v[I-1][j-1])*1)*(-F_s_v[I-1][j-1])+D_s ;
                    S2[I-1][j-1] += -(-p_star[I-1][j-2])*dx/rho 
                    -v_star_old[I-1][j-3]*(alpha_s_v[I-1][j-1]*UU(r_s_v[I-1][j-1]))*(-F_s_v[I-1][j-1]) 
                    -v_star_old[I-1][j-2]*((1-alpha_s_v[I-1][j-1])*D(r_s_v[I-1][j-1])+alpha_s_v[I-1][j-1]*(U(r_s_v[I-1][j-1])-1))*(-F_s_v[I-1][j-1])
                    -v_star_old[I-1][j]*((1-alpha_s_v[I-1][j-1])*UU(r_s_v[I-1][j-1]))*(-F_s_v[I-1][j-1]) 
                    -v_star_old[I-1][j-1]*((1-alpha_s_v[I-1][j-1])*(U(r_s_v[I-1][j-1])-1)+alpha_s_v[I-1][j-1]*D(r_s_v[I-1][j-1]))*(-F_s_v[I-1][j-1])
                    ;//deferred correction 完成 
                }
                //北邊界
                if(j == Ny){//上邊界第二排計算點的北側邊界採用迎風格式
                    b_N[I-1][j-1] += (1-alpha_n_v[I-1][j-1])*(F_n_v[I-1][j-1])-D_n ;//討論北邊界時，再扣掉北邊界擴散部分係數
                    b_P[I-1][j-1] += alpha_n_v[I-1][j-1]*(F_n_v[I-1][j-1])+D_n; 
                    S2[I-1][j-1] += -(p_star[I-1][j-1])*dx/rho ;
                }else{//其他點北側邊界採用TVD-MUSCL格式//北邊界千牽涉本點以及北點
                    //b_NN[I-1][j-1] += (1-alpha_n_v[I-1][j-1])*UU(r_n_v[I-1][j-1])*(F_n_v[I-1][j-1]) ;
                    b_N[I-1][j-1] += ((1-alpha_n_v[I-1][j-1])*1)*(F_n_v[I-1][j-1]) - D_n ;//討論北邊界時，再扣掉北邊界擴散部分係數
                    //b_S[I-1][j-1] += (alpha_n_v[I-1][j-1])*UU(r_n_v[I-1][j-1])*(F_n_v[I-1][j-1]) ;
                    b_P[I-1][j-1] += (alpha_n_v[I-1][j-1]*1)*(F_n_v[I-1][j-1]) + D_n ;
                    S2[I-1][j-1] += -(p_star[I-1][j-1])*dx/rho
                    -v_star_old[I-1][j+1]*(1-alpha_n_v[I-1][j-1])*UU(r_n_v[I-1][j-1])*(F_n_v[I-1][j-1]) 
                    -v_star_old[I-1][j]*(alpha_n_v[I-1][j-1]*D(r_n_v[I-1][j-1])+(1-alpha_n_v[I-1][j-1])*(U(r_n_v[I-1][j-1])-1))*(F_n_v[I-1][j-1]) 
                    -v_star_old[I-1][j-1]*(alpha_n_v[I-1][j-1])*UU(r_n_v[I-1][j-1])*(F_n_v[I-1][j-1]) 
                    -v_star_old[I-1][j-1]*((1-alpha_n_v[I-1][j-1])*D(r_n_v[I-1][j-1])+alpha_n_v[I-1][j-1]*(U(r_n_v[I-1][j-1])-1))*(F_n_v[I-1][j-1]) 
                    ; //deferred correction 完成 !!
                }
                //西邊界(分三種情況考慮：I = 1,2,Nx)
                if(I == 2){//左邊界第二排計算點的西側邊界使用外插法取代西西點上的取值WW->W
                    b_W[I-1][j-1] += ((alpha_w_v[I-1][j-1])*1)*(-F_w_v[I-1][j-1])-D_w;
                    //b_E[I-1][j-1] += ((1-alpha_w_v[I-1][j-1])*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1]) ;
                    b_P[I-1][j-1] += ((1-alpha_w_v[I-1][j-1])*1)*(-F_w_v[I-1][j-1])+D_w;
                    S2[I-1][j-1] +=  -2*(alpha_w_v[I-1][j-1]*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1])*v_left 
                    +v_star_old[I-2][j-1]*(alpha_w_v[I-1][j-1]*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1])
                    -v_star_old[I-2][j-1]*((1-alpha_w_v[I-1][j-1])*D(r_w_v[I-1][j-1])+(alpha_w_v[I-1][j-1])*(U(r_w_v[I-1][j-1])-1))*(-F_w_v[I-1][j-1])
                    -v_star_old[I][j-1]*((1-alpha_w_v[I-1][j-1])*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1])
                    -v_star_old[I-1][j-1]*((1-alpha_w_v[I-1][j-1])*(U(r_w_v[I-1][j-1])-1)+(alpha_w_v[I-1][j-1])*D(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1])
                    ; //deferred correction完成 !!
                }else if(I == 1){//左邊界第一排計算點的西側邊界採用迎風格式
                    b_P[I-1][j-1] +=  2*D_w;//討論西邊界時，再減去西邊界擴散部分係數
                    S2[I-1][j-1] += -(-F_w_v[I-1][j-1])*v_left+2*D_w*v_left;
                }else if(I == Nx){//E->P
                    //由於當I = Nx，時，東邊界計算點不存在，使用外插法，插值過去，插到本點過去
                    //b_WW[I-1][j-1] += (alpha_w_v[I-1][j-1]*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1]);
                    b_W[I-1][j-1] += (alpha_w_v[I-1][j-1]*1)*(-F_w_v[I-1][j-1]) -D_w ;
                    b_P[I-1][j-1] += ((1-alpha_w_v[I-1][j-1])*1)*(-F_w_v[I-1][j-1]) +D_w;
                    S2[I-1][j-1] += -2*v_right*(((1-alpha_w_v[I-1][j-1])*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1])) 
                    +v_star_old[I-1][j-1]*((1-alpha_w_v[I-1][j-1])*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1])
                    -v_star_old[I-3][j-1]*(alpha_w_v[I-1][j-1]*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1])
                    -v_star_old[I-2][j-1]*((1-alpha_w_v[I-1][j-1])*D(r_w_v[I-1][j-1])+alpha_w_v[I-1][j-1]*(U(r_w_v[I-1][j-1])-1))*(-F_w_v[I-1][j-1])
                    -v_star_old[I-1][j-1]*((1-alpha_w_v[I-1][j-1])*(U(r_w_v[I-1][j-1])-1)+alpha_w_v[I-1][j-1]*D(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1])
                    ;//deferred correction完成 !!
                }else{//其他點西側邊界採用TVD-MUSCL格式
                    //b_WW[I-1][j-1] += (alpha_w_v[I-1][j-1]*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1]) ;
                    b_W[I-1][j-1] += ((alpha_w_v[I-1][j-1])*1)*(-F_w_v[I-1][j-1])-D_w ; //討論西邊界時，再減去西側擃散部分係數
                    //b_E[I-1][j-1] += ((1-alpha_w_v[I-1][j-1])*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1]);
                    b_P[I-1][j-1] += ((1-alpha_w_v[I-1][j-1])*1)*(-F_w_v[I-1][j-1])+D_w ;
                    S2[I-1][j-1] += 
                    -v_star_old[I-3][j-1]*(alpha_w_v[I-1][j-1]*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1]) 
                    -v_star_old[I-2][j-1]*((1-alpha_w_v[I-1][j-1])*D(r_w_v[I-1][j-1])+(alpha_w_v[I-1][j-1])*(U(r_w_v[I-1][j-1])-1))*(-F_w_v[I-1][j-1])
                    -v_star_old[I][j-1]*((1-alpha_w_v[I-1][j-1])*UU(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1])
                    -v_star_old[I-1][j-1]*((1-alpha_w_v[I-1][j-1])*(U(r_w_v[I-1][j-1])-1)+(alpha_w_v[I-1][j-1])*D(r_w_v[I-1][j-1]))*(-F_w_v[I-1][j-1])
                    ;//deferred correction完成 !!
                }
                //東邊界(分三種情況考慮：I = 1,Nx-1,Nx)
                if(I == 1){//右邊界第二排計算點的東側邊界使用外插法取代東東點上的取值 W->P : 西點上的直插到本點上//東邊界deferred correction 牽涉東點以及本點
                    //b_EE[I-1][j-1] += ((1-alpha_e_v[I-1][j-1])*UU(r_e_v[I-1][j-1]))*F_e_v[I-1][j-1] ;
                    b_E[I-1][j-1] += ((1-alpha_e_v[I-1][j-1])*1)*(F_e_v[I-1][j-1])-D_e ;//討論東邊界時，再扣掉東邊界擴散部分係數
                    b_P[I-1][j-1] += ((alpha_e_v[I-1][j-1])*1)*(F_e_v[I-1][j-1])  +D_e ;
                    S2[I-1][j-1] +=  -2*v_left*(alpha_e_v[I-1][j-1]*UU(r_e_v[I-1][j-1]))*(F_e_v[I-1][j-1]) 
                    +v_star_old[I-1][j-1]*(alpha_e_v[I-1][j-1]*UU(r_e_v[I-1][j-1]))*(F_e_v[I-1][j-1])
                    -v_star_old[I+1][j-1]*((1-alpha_e_v[I-1][j-1])*UU(r_e_v[I-1][j-1]))*F_e_v[I-1][j-1]
                    -v_star_old[I][j-1]*((alpha_e_v[I-1][j-1])*D(r_e_v[I-1][j-1])+(1-alpha_e_v[I-1][j-1])*(U(r_e_v[I-1][j-1])-1))*(F_e_v[I-1][j-1])
                    -v_star_old[I-1][j-1]*((alpha_e_v[I-1][j-1])*(U(r_e_v[I-1][j-1])-1)+(1-alpha_e_v[I-1][j-1])*D(r_e_v[I-1][j-1]))*(F_e_v[I-1][j-1])
                    ;//deferredm correction 完成 !!
                }else if(I == Nx){//右邊界第二排計算點的東側邊界採用迎風格式
                    b_P[I-1][j-1] += 2*D_e ;//討論東邊界時，再扣掉東邊界擴散部分係數
                    S2[I-1][j-1] += -(F_e_v[I-1][j-1])*v_right +2*D_e*v_right;
                }else if(I == Nx-1){//EE點不存在，差值到E點
                    b_E[I-1][j-1] += ((1-alpha_e_v[I-1][j-1])*1)*F_e_v[I-1][j-1] 
                    - D_e ;
                    //b_W[I-1][j-1] += (alpha_e_v[I-1][j-1]*UU(r_e_v[I-1][j-1]))*F_e_v[I-1][j-1] ;
                    b_P[I-1][j-1] += (alpha_e_v[I-1][j-1]*1)*F_e_v[I-1][j-1] + D_e ;
                    S2[I-1][j-1] += -2*v_right*((1-alpha_e_v[I-1][j-1])*UU(r_e_v[I-1][j-1]))*F_e_v[I-1][j-1] 
                    +v_star_old[I][j-1]*((1-alpha_e_v[I-1][j-1])*UU(r_e_v[I-1][j-1]))*F_e_v[I-1][j-1]
                    -v_star_old[I][j-1]*((1-alpha_e_v[I-1][j-1])*(U(r_e_v[I-1][j-1])-1)+alpha_e_v[I-1][j-1]*D(r_e_v[I-1][j-1]))*F_e_v[I-1][j-1]
                    -v_star_old[I-2][j-1]*(alpha_e_v[I-1][j-1]*UU(r_e_v[I-1][j-1]))*F_e_v[I-1][j-1] 
                    -v_star_old[I-1][j-1]*((1-alpha_e_v[I-1][j-1])*D(r_e_v[I-1][j-1])+alpha_e_v[I-1][j-1]*(U(r_e_v[I-1][j-1])-1))*F_e_v[I-1][j-1] 
                    ;
                }else{//其他點東側邊界採用TVD-MUSCL格式
                    //b_EE[I-1][j-1] += (1-alpha_e_v[I-1][j-1])*UU(r_e_v[I-1][j-1])*(F_e_v[I-1][j-1]) ;
                    b_E[I-1][j-1] += ((1-alpha_e_v[I-1][j-1])*1)*(F_e_v[I-1][j-1]) - D_e ;//討論東邊界時，再扣掉東邊界擴散部分係數
                    //b_W[I-1][j-1] += (alpha_e_v[I-1][j-1])*UU(r_e_v[I-1][j-1])*(F_e_v[I-1][j-1]) ;
                    b_P[I-1][j-1] += (alpha_e_v[I-1][j-1]*1)*(F_e_v[I-1][j-1]) + D_e ;
                    S2[I-1][j-1] += 
                    -v_star_old[I+1][j-1]*(1-alpha_e_v[I-1][j-1])*UU(r_e_v[I-1][j-1])*(F_e_v[I-1][j-1]) 
                    -v_star_old[I][j-1]*(alpha_e_v[I-1][j-1]*D(r_e_v[I-1][j-1])+(1-alpha_e_v[I-1][j-1])*(U(r_e_v[I-1][j-1])-1))*(F_e_v[I-1][j-1]) 
                    -v_star_old[I-2][j-1]*(alpha_e_v[I-1][j-1])*UU(r_e_v[I-1][j-1])*(F_e_v[I-1][j-1]) 
                    -v_star_old[I-1][j-1]*((1-alpha_e_v[I-1][j-1])*D(r_e_v[I-1][j-1])+alpha_e_v[I-1][j-1]*(U(r_e_v[I-1][j-1])-1))*(F_e_v[I-1][j-1])
                    ;//deferred correction 完成
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
                //速度場預估值的鬆弛迭代關係式
                v_star[I-1][j-1] = (1-omega_u)*v_star_old[I-1][j-1] + omega_u*v_star[I-1][j-1] ;
            }
        }
    }
}
//得到了N+\frac{1}{2}計算步的u_star、v_star
//函數：第一次更新速度場預估值後的mass flux
double d(int I, int J , int m){
    //I,J實際上為主計算點座標，m為方程式的迭代步數，與當前速度場的計算步數不同
    double u_e = u_star[I+1-1][J-1] ;//此時速度場的計算步數為N+1/2
    double u_w = u_star[I-1][J-1] ;//此時速度場的計算步數為N+1/2
    double v_n = v_star[I-1][J+1-1] ;//此時速度場的計算步數為N+1/2
    double v_s = v_star[I-1][J-1] ;//此時速度場的計算步數為N+1/2
    return rho*(u_e - u_w)*dy + rho*(v_n - v_s)*dx ;
} //上式之判斷用在速度場預估值第一次更新後

//step4：更新壓力場正確值
void step4(int m){
    #pragma omp parallel for collapse(2) 
    for(int I = 1 ; I <= Nx ; I++){ //壓力場不管正確值、預估值、修正值，都取在主計算點
        for(int J = 1 ; J <= Ny ; J++){
           p[I-1][J-1] = p_star[I-1][J-1] + omega_P*p_prime[I-1][J-1] ;
        }//p_star[I-1][J-1]的計算步數為N; p_prime[I-1][J-1]為N+1;
    }//p_prime[I-1][J-1]的計算步數為N+1 ;
}

//step5：更新速度場修正值，需要注意，因為物理邊界上的速度場已經由邊界條件更新，不做迭代，因此。物理邊界上的速度場(即i = 1 && i = Nx+1)的速度場修正值為0.0
void step5(int m){
    #pragma omp parallel for collapse(2) 
    //更新x分量速度場修正值
    for(int i = 1 ; i <= Nx+1 ; i++){
        for(int J = 1 ; J <= Ny ; J++){
            if(i == 1 || i == Nx+1){
                u_prime[i-1][J-1] = 0.0 ;//物理邊界上的速度場已經由邊界條件更新，不做迭代
            }else{//壓力場修正值：本點-西鄰計算點
                u_prime[i-1][J-1] = -(p_prime[i-1][J-1] - p_prime[i-2][J-1])*dy/(rho*a_P[i-1][J-1]) ; ;
            }

        }
    }
    //更新y量速度場修正值
    #pragma omp parallel for collapse(2) 
    for(int I = 1 ; I <= Nx ; I++){
        for(int j = 1 ; j <= Ny+1 ; j++){
            if(j == 1 || j == Ny+1){
                v_prime[I-1][j-1] = 0.0 ;//物理邊界上的速度場已經由邊界條件更新，不做迭代
            }else{//壓力場修正值：本點-南鄰計算點
                v_prime[I-1][j-1] = -(p_prime[I-1][j-1] - p_prime[I-1][j-2])*dx/(rho*b_P[I-1][j-1]) ;
            }
        }
    }
}

//step6:更新速度場正確值&&重新賦值舊速度場正確值
void step6(int m){//如前所述，物理邊界上的速度場不應該迭代，因此，物理邊界上的速度場必然由滿足邊界條件的速度場預估值給定
    //更新x分量素速度場正確值
    #pragma omp parallel for collapse(2) 
    for(int i = 1 ; i <= Nx+1 ; i++){
        for(int J = 1 ; J <= Ny ; J++){
            if(i == 1 ){
                u[i-1][J-1] = 0.0 ;//物理邊界上的速度場已經由邊界條件更新，不做迭代
            }else if (i == Nx+1){
                u[i-1][J-1] = 0.0 ;//物理邊界上的速度場已經由邊界條件更新，不做迭代
            }else{
                u[i-1][J-1] = u_prime[i-1][J-1]+u_star[i-1][J-1];//u_star[i-1][J-1]的計算步數為N+1/2; u_prime[i-1][J-1]的計算步數為N+1; u_old[i-1][J-1]的計算步數為N ; u[i-1][J-1]的計算步數為N+1 ; 
            }
        }
    }
    //重新賦值舊x分量速度場正確值
    #pragma omp parallel for collapse(2) 
    for(int i = 1 ; i <= Nx+1 ; i++){
        for(int J = 1 ; J <= Ny ; J++){
            u_old[i-1][J-1] = u[i-1][J-1] ; //u[i-1][J-1]的計算步數為N+1 ; u_old[i-1][J-1]的計算步數為N+1(當前x分量速度場正確值);
        }
    }
    //更新y分量素速度場正確值
    #pragma omp parallel for collapse(2) 
    for(int I = 1 ; I <= Nx ; I++){
        for(int j = 1 ; j <= Ny+1 ; j++){
            if(j == 1 ){
                v[I-1][j-1] = 0.0 ;//物理邊界上的速度場已經由邊界條件更新，不做迭代//物面不可穿透條件
            }else if (j == Ny+1){
                v[I-1][j-1] = 0.0 ;//物理邊界上的速度場已經由邊界條件更新，不做迭代//物面不可穿透條件
            }else{
                v[I-1][j-1] = v_prime[I-1][j-1]+v_star[I-1][j-1];//v_star[I-1][j-1]的計算步數為N+1/2; v_prime[I-1][j-1]的計算步數為N+1; v_old[I-1][j-1]的計算步數為N ; v[I-1][j-1]的計算步數為N+1 ; 
            }
        }
    }
    //重新賦值舊y分量速度場正確值
    #pragma omp parallel for collapse(2) 
    for(int I = 1 ; I <= Nx ; I++){
        for(int j = 1 ; j <= Ny+1 ; j++){
            v_old[I-1][j-1] = v[I-1][j-1] ; //v[I-1][j-1]的計算步數為N+1 ; v_old[I-1][j-1]的計算步數為N+1(當前y分量速度場正確值);
        }
    }
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
//step7 : 計算控制表面的質量通量（取非平衡質量流動的總和）
double step7(int m){
    /*//81
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
    return sum ;*/


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
    sum += fabs(v[111][56]-(0.0)) ;
    //161
    /*double sum = 0.0 ;
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
    sum += fabs(v[159][80]-(0.0)) ;*/
    //227
    /*double sum = 0.0 ;
    sum += fabs(u[80][0]-0) ; 
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
}
    
//step8 : 更新主計算點的壓力場預估值
void step8(int m){
    #pragma omp parallel for collapse(2)  // *** 加這行 ***
    for(int I = 1 ; I <= Nx ; I++){
        for(int J = 1 ; J <= Ny ; J++){
            p_star[I-1][J-1] = p[I-1][J-1] ; //p[I-1][J-1]的計算步數為N+1 ; p_star[I-1][J-1]的計算步數為N+1 ;
        }
    }
}    
    
//step9 : 更新速度場預估值，更新介面判斷係數，更新對流部分係數
void step9(int m){
    //更新x分量速度場預估值：需要注意的是，物理邊界上的速度場由邊界條件直接給定不需迭代賦值
    #pragma omp parallel for collapse(2)
    for(int i = 1 ; i <= Nx+1 ; i++){
        for(int J = 1 ; J <= Ny ; J++){
            if(i == 1){
                u_star[i-1][J-1] = u_left;
            }else if (i == Nx+1){
                u_star[i-1][J-1] = u_right ;
            }else{
                //其他內部點，包含上下邊界點由正確值賦予取值
                u_star[i-1][J-1] = u[i-1][J-1]; //u[i-1][J-1]的計算步數為N+1 ; u_star[i-1][J-1]的計算步數為N+1 ;
                u_star_old[i-1][J-1] = u_star[i-1][J-1]; //馬上紀錄速度場預估值
            }
        }
    }
    //更新y分量速度場預估值：需要注意的是，物理邊界上的速度場由邊界條件直接給定不需迭代賦值
    #pragma omp parallel for collapse(2)
    for(int I = 1 ; I <= Nx ; I++){
        for(int j = 1 ; j <= Ny+1 ; j++){
            if(j == 1){
                v_star[I-1][j-1] = v_bottom ;
            }else if (j == Ny+1){
                v_star[I-1][j-1] = v_top ;
            }else{
                //其他內部點，包含左右邊界點由正確值賦予取值
                v_star[I-1][j-1] = v[I-1][j-1]; //v[I-1][j-1]的計算步數為N+1 ; v_star[I-1][j-1]的計算步數為N+1 ;
                v_star_old[I-1][j-1] = v_star[I-1][j-1]; //馬上紀錄速度場預估值
            }
        }
    }
    //更新對流項係數
    //x分量動量積分方程的一般離散格式的對流部分係數
    #pragma omp parallel for collapse(2)
    for(int i = 1 ; i <= Nx+1 ; i++){
        for(int J = 1 ; J <= Ny ; J++){
            if(i == 1 || i == Nx+1) {
                F_w_u[i-1][J-1] = 0.0 ;
                F_e_u[i-1][J-1] = 0.0 ;
                F_n_u[i-1][J-1] = 0.0 ;
                F_s_u[i-1][J-1] = 0.0 ;
            }else{//重新計算，對流項係數
                F_e_u[i-1][J-1] = 0.5*dy*(u_star[i-1][J-1]+u_star[i][J-1]); //網格中心點剛好對應到變量儲存點
                F_w_u[i-1][J-1] = 0.5*dy*(u_star[i-1][J-1]+u_star[i-2][J-1]);
                F_n_u[i-1][J-1] = 0.5*dx*(v_star[i-1][J]+v_star[i-2][J]);
                F_s_u[i-1][J-1] = 0.5*dx*(v_star[i-1][J-1]+v_star[i-2][J-1]);
            }
        }
    }//alpha_f_u判斷因子，下標代表控制體積的中心點座標
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
    //y分量動量積分方程的一般離散格式的對流部分係數
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
                F_n_v[I-1][j-1] = 0.5*dx*(v_star[I-1][j-1]+v_star[I-1][j]);//網格中心點剛好對應到變量儲存點
                F_s_v[I-1][j-1] = 0.5*dx*(v_star[I-1][j-1]+v_star[I-1][j-2]);
            }
        }
    }//alpha_f_v判斷因子，其下標代表控制體積中心點座標
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
    //更新介面判斷因子
    #pragma omp parallel for collapse(2)
    for(int i = 1 ; i <= Nx+1 ; i++){
        for(int J = 1 ; J <= Ny ; J++){
            //不需要討論的cell->介面判斷係數直接設為0
            if(i == 1 || i == Nx+1){
                r_w_u[i-1][J-1] = 0.0 ;
                r_e_u[i-1][J-1] = 0.0 ;
                r_n_u[i-1][J-1] = 0.0 ;
                r_s_u[i-1][J-1] = 0.0 ;
            }else {
                //按照各個邊界討論（西邊界特殊討論＊1東邊界特輸討論*南邊界特殊討論*北邊界特殊討論*3)
                //西邊界
                if(i == 2){
                    r_w_u[i-1][J-1] = -1.0 ; //西邊界採用迎風格式
                }else{
                    double u_W = u_star[i-2][J-1]; 
                    double u_P = u_star[i-1][J-1];
                    double u_WW = u_star[i-3][J-1];
                    double u_E = u_star[i][J-1] ;
                    double D = u_W + alpha_w_u[i-1][J-1]*(u_P - u_W) ;
                    double U = u_P + alpha_w_u[i-1][J-1]*(u_W - u_P) ;
                    double UU = u_E + alpha_w_u[i-1][J-1]*(u_WW - u_E) ;
                    r_w_u[i-1][J-1] = (U-UU)/(D-U+1e-30) ; //初始化設計
                }
                //東邊界
                if(i == Nx){
                    r_e_u[i-1][J-1] = -1.0 ; //東邊界採用迎風格式
                }else{
                    double u_P = u_star[i-1][J-1] ;
                    double u_E = u_star[i][J-1] ;
                    double u_W = u_star[i-2][J-1] ;
                    double u_EE = u_star[i+1][J-1] ;
                    double D = u_P + alpha_e_u[i-1][J-1]*(u_E - u_P) ;
                    double U = u_E + alpha_e_u[i-1][J-1]*(u_P - u_E) ;
                    double UU = u_EE + alpha_e_u[i-1][J-1]*(u_W - u_EE) ;
                    r_e_u[i-1][J-1] = (U-UU)/(D-U+1e-30) ; //初始化設計
                }
                //南邊界
                if(J == 1){
                    r_s_u[i-1][J-1] = 0.0 ; //下邊界計算點用不到南邊界判斷因子
                }else if(J == 2){
                    double u_S = u_star[i-1][J-2] ;
                    double u_P = u_star[i-1][J-1] ;
                    double u_SS = 2*u_bottom-u_S ; //2*邊界條件 - 南點
                    double u_N = u_star[i-1][J] ;
                    double D = u_S + alpha_s_u[i-1][J-1]*(u_P - u_S) ;
                    double U = u_P + alpha_s_u[i-1][J-1]*(u_S - u_P) ;
                    double UU = u_N + alpha_s_u[i-1][J-1]*(u_SS - u_N) ;
                    r_s_u[i-1][J-1] =  (U-UU)/(D-U+1e-30) ; //初始化設計
                }else if(J == Ny){
                    double u_S = u_star[i-1][J-2] ;
                    double u_P = u_star[i-1][J-1] ;
                    double u_SS = u_star[i-1][J-3] ;
                    double u_N = 2*u_top - u_P ; //2*邊界條件 - 本點
                    double D = u_S + alpha_s_u[i-1][J-1]*(u_P - u_S) ;
                    double U = u_P + alpha_s_u[i-1][J-1]*(u_S - u_P) ;
                    double UU = u_N + alpha_s_u[i-1][J-1]*(u_SS - u_N) ;
                    r_s_u[i-1][J-1] =  (U-UU)/(D-U+1e-30) ; //初始化設計
                }else{
                    double u_S = u_star[i-1][J-2] ;
                    double u_P = u_star[i-1][J-1] ;
                    double u_SS = u_star[i-1][J-3] ;
                    double u_N = u_star[i-1][J] ;
                    double D = u_S + alpha_s_u[i-1][J-1]*(u_P - u_S) ;
                    double U = u_P + alpha_s_u[i-1][J-1]*(u_S - u_P) ;
                    double UU = u_N + alpha_s_u[i-1][J-1]*(u_SS - u_N) ;
                    r_s_u[i-1][J-1] =  (U-UU)/(D-U+1e-30) ; //初始化設計
                }
                //北邊界
                if(J == 1){
                    double u_P = u_star[i-1][J-1] ;
                    double u_N = u_star[i-1][J] ;
                    double u_S = 2*u_bottom-u_P ; //2*邊界條件 - 本點
                    double u_NN = u_star[i-1][J+1] ;
                    double D = u_P + alpha_n_u[i-1][J-1]*(u_N - u_P) ;
                    double U = u_N + alpha_n_u[i-1][J-1]*(u_P - u_N) ;
                    double UU = u_NN + alpha_n_u[i-1][J-1]*(u_S - u_NN) ;
                    r_n_u[i-1][J-1] =  (U-UU)/(D-U+1e-30) ; //初始化設計
                }else if (J == Ny-1){
                    double u_P = u_star[i-1][J-1] ;
                    double u_N = u_star[i-1][J] ;
                    double u_S = u_star[i-1][J-2] ;
                    double u_NN = 2*u_top - u_N ; //2*邊界條件 - 北點
                    double D = u_P + alpha_n_u[i-1][J-1]*(u_N - u_P) ;
                    double U = u_N + alpha_n_u[i-1][J-1]*(u_P - u_N) ;
                    double UU = u_NN + alpha_n_u[i-1][J-1]*(u_S - u_NN) ;
                    r_n_u[i-1][J-1] =  (U-UU)/(D-U+1e-30) ; //初始化設計
                }else if(J == Ny){
                    r_n_u[i-1][J-1] = 0.0 ; //上邊界計算點用不到北邊界判斷因子
                }else{
                    double u_P = u_star[i-1][J-1] ;
                    double u_N = u_star[i-1][J] ;
                    double u_S = u_star[i-1][J-2] ;
                    double u_NN = u_star[i-1][J+1] ;
                    double D = u_P + alpha_n_u[i-1][J-1]*(u_N - u_P) ;
                    double U = u_N + alpha_n_u[i-1][J-1]*(u_P - u_N) ;
                    double UU = u_NN + alpha_n_u[i-1][J-1]*(u_S - u_NN) ;
                    r_n_u[i-1][J-1] =  (U-UU)/(D-U+1e-30) ; //初始化設計
                }
            }
        }    
    }//記得預估值更新後，立即更新介面判斷係數
    //y分量速度場介面判斷係數(東邊界判斷細數，西邊界判斷係數需要特殊討論)
    #pragma omp parallel for collapse(2)
    for(int I = 1 ; I <= Nx ; I++){
        for(int j = 1 ; j <= Ny+1 ; j++){
            //不需要討論的cell->介面判斷係數直接設為0
            if(j == 1 || j == Ny+1){
                r_w_v[I-1][j-1] = 0.0 ;
                r_e_v[I-1][j-1] = 0.0 ;
                r_n_v[I-1][j-1] = 0.0 ;
                r_s_v[I-1][j-1] = 0.0 ;
            }else {
                //按照各個邊界依序討論（南邊界特殊討論*1北邊界特殊討論*1東邊界特殊討論*3西邊界特殊討論*3)
                //南邊界
                if(j == 2){
                    r_s_v[I-1][j-1] = -1.0 ; //下邊界計算點的南邊界採用迎風格式
                }else{
                    double v_S = v_star[I-1][j-2] ;
                    double v_P = v_star[I-1][j-1] ;
                    double v_SS = v_star[I-1][j-3] ;
                    double v_N = v_star[I-1][j] ;
                    double D = v_S + alpha_s_v[I-1][j-1]*(v_P - v_S) ;
                    double U = v_P + alpha_s_v[I-1][j-1]*(v_S - v_P) ;
                    double UU = v_N + alpha_s_v[I-1][j-1]*(v_SS - v_N) ;
                    r_s_v[I-1][j-1] = (U-UU)/(D-U+1e-30) ; //避免分母為0
                }
                //北邊界
                if(j == Ny){
                    r_n_v[I-1][j-1] = -1.0 ; //上邊界計算點的北邊界採用迎風格式
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
                //西邊界
                if(I == 1){
                    r_w_v[I-1][j-1] = 0.0 ; //西邊界計算點用不到西邊界判斷因子
                }else if (I == 2){
                    double v_W = v_star[I-2][j-1] ;
                    double v_P = v_star[I-1][j-1] ;
                    double v_WW = 2*v_left-v_W ; //2*邊界條件 - 西點
                    double v_E = v_star[I][j-1] ;
                    double D = v_W + alpha_w_v[I-1][j-1]*(v_P - v_W) ;
                    double U = v_P + alpha_w_v[I-1][j-1]*(v_W - v_P) ;
                    double UU = v_E + alpha_w_v[I-1][j-1]*(v_WW - v_E) ;
                    r_w_v[I-1][j-1] = (U-UU)/(D-U+1e-30) ;
                }else if(I == Nx){
                    double v_W = v_star[I-2][j-1] ;
                    double v_P = v_star[I-1][j-1] ;
                    double v_WW = v_star[I-3][j-1] ;
                    double v_E = 2*v_right - v_P ; //2*邊界條件 - 本點
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
                //東邊界
                if(I == 1){
                    double v_P = v_star[I-1][j-1] ;
                    double v_E = v_star[I][j-1] ;
                    double v_W = 2*v_left - v_P ; //2*邊界條件 - 本點
                    double v_EE = v_star[I+1][j-1] ;
                    double D = v_P + alpha_e_v[I-1][j-1]*(v_E - v_P) ;
                    double U = v_E + alpha_e_v[I-1][j-1]*(v_W - v_E) ;
                    double UU = v_EE + alpha_e_v[I-1][j-1]*(v_W - v_EE) ;
                    r_e_v[I-1][j-1] = (U-UU)/(D-U+1e-30) ;
                }else if (I == Nx-1){
                    double v_P = v_star[I-1][j-1] ;
                    double v_E = v_star[I][j-1] ;
                    double v_W = v_star[I-2][j-1] ;
                    double v_EE = 2*v_right - v_E ; //2*邊界條件 - 東點
                    double D = v_P + alpha_e_v[I-1][j-1]*(v_E - v_P) ;
                    double U = v_E + alpha_e_v[I-1][j-1]*(v_W - v_E) ;
                    double UU = v_EE + alpha_e_v[I-1][j-1]*(v_W - v_EE) ;
                    r_e_v[I-1][j-1] = (U-UU)/(D-U+1e-30) ;
                }else if(I == Nx){
                    r_e_v[I-1][j-1] = 0.0 ; //東邊界計算點用不到東邊界判斷因子
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
/*void output_combined_fields(int m) {
    const int cells_i = Nx; // cells in x
    const int cells_j = Ny; // cells in y

    std::ostringstream name;
    name << "3.QUICK_Re" << Re << "_Grid" << Nx << "x" << Ny
         << "_DC_step" << std::setfill('0') << std::setw(6) << m 
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
        << "DATAPACKING=POINT\n"; // 移除尾逗號

    // AUXDATA（字串加引號較穩）
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
            const double pC = p[i][j]; // 假設 p 尺寸為 [Nx][Ny]

            // u, v 由面值內插到 cell center（MAC: u (i+0.5,j), v(i,j+0.5)）
            // 需要 u 尺寸為 [Nx+1][Ny]，v 尺寸為 [Nx][Ny+1]
            const double uC = 0.5 * (u[i][j] + u[i+1][j]);
            const double vC = 0.5 * (v[i][j] + v[i][j+1]);

            const double vel_mag = std::sqrt(uC*uC + vC*vC);

            // vorticity at center, 與 MAC 網格一致的差分（一次差分）
            double vort = 0.0;
            // 需要 i+1 <= Nx-1 => i <= Nx-2，j-1 >= 0 => j >= 1
            if (i <= cells_i - 2 && j >= 1) {
                const double dvdx = (v[i+1][j] - v[i][j]) / dx; // v: Nx × (Ny+1)，用同一 j
                const double dudy = (u[i][j]   - u[i][j-1]) / dy; // u: (Nx+1) × Ny，用同一 i
                vort = dvdx - dudy;
            }
            // 也可對邊界用單邊差分以避免 0

            out << xC << " " << yC << " "
                << pC << " " << uC << " " << vC << " "
                << vel_mag << " " << vort << "\n";
        }
    }

    out.close();
    std::cout << "Tecplot file successfully written: " << name.str() << std::endl;
}
*/

//有限體積法輸出cell data
void output_combined_fields(int m) {
    int nodes_i = Nx + 1;  // x方向節點數
    int nodes_j = Ny + 1;  // y方向節點數
    int cells_i = Nx;      // x方向網格數
    int cells_j = Ny;      // y方向網格數
    
    std::ostringstream name;
    name << "3.QUICK_Re" << Re << "_Grid" << Nx << "x" << Ny
         << "_DC_step" << std::setfill('0') << std::setw(6) << m 
         << "_error" << min_step7 <<
          ".vtk";

    std::ofstream out(name.str());
    
    out << "# vtk DataFile Version 3.0\n";
    out << "3.QUICK_Re = 5000_Grid 161x161\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_GRID\n";
    out << "DIMENSIONS " << nodes_i << " " << nodes_j << " 1\n";
    
    // 輸出節點座標
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
    
    // 1. 輸出壓力場
    out << "SCALARS Pressure double 1\n";
    out << "LOOKUP_TABLE default\n";
    for(int J = 1; J <= cells_j; J++) {
        for(int I = 1; I <= cells_i; I++) {
            out << scientific << setprecision(6) 
                << p_star[I-1][J-1] << "\n";
        }
    }
    
    // 2. 輸出速度場（不需要再宣告 CELL_DATA）
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
    // *** 新增這3行 ***
    int num_threads = omp_get_max_threads();
    omp_set_num_threads(num_threads);
    cout << "Using " << num_threads << " threads" << endl;
    cout << "Using FVD with SIMPLE algorithm to solve 2D incompressible Navier-Stokes equations" << endl;
    cout << "Grid Size:" << Nx << "x" << Ny << endl;
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
        
        step2(G) ;//速度場預估值更新！(第一次更新！)
        step3() ;//壓力場修正值更新！
        step4(G) ;//壓力場正確值更新！
        step5(G) ;//速度場修正值更新！
        step6(G) ;//速度場正確值更新
        //step7:判斷收斂步驟
        if(min_step7 > step7(G)){ //每一步都判斷
                min_step7 = step7(G) ;
                min_G = G ; 
                //if(G > 1000 && min_step7 < 2.55   ){output_combined_fields(min_G) ;}
               
            }
        if(G%100 == 0){
            cout << "Itteration step = " << G << endl ;
            cout << "Non-equilibrium mass flow rate for (G = "<< G << " )step:" << step5star(G) << endl ;
            cout << "currnt minerror(G =" << min_G << ")= "<< min_step7 << endl ; 
        }
        if(G>0 && step7(G) < 0.1){
            CE_physical = true ;
            cout << "Haved converged at step G = " << G << endl ;
            cout << "Non-equilibrium mass flow rate for (G = "<< G << " )step:" << step5star(G) <<endl ;
            cout << "currnt minerror(G =" << G << ")= "<< min_step7 << endl ;
            output_combined_fields(G) ;
            break ;
        }
        step8(G) ;//壓力場預估值更新
        step9(G) ;//速度場預估值更新，對流部分係數更新，介面判斷係數更新
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






