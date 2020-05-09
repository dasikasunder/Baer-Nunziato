/*
 * hype1d.cc
 *      Author: sunder
 */

//----------------------------------------------------------------------------
// Commonly used C/C++ header files
//----------------------------------------------------------------------------

#include <iostream>
#include <string>
#include <cmath>
#include <iomanip>
#include <cassert>
#include <ctime>
#include <vector>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <array>
#include <map>
#include <utility>
#define BOOST_DISABLE_ASSERTS
#include "boost/multi_array.hpp"

using namespace boost; 

//----------------------------------------------------------------------------
// Various constants used throught the code
//----------------------------------------------------------------------------

const int nVar = 7;             // Number of components in the PDE system //
const int nLin = 3;             // Number of linear degenerate fields in the PDE system //

const double GAMMA_S       = 3.0;     // Specific heat ratio of the solid phase //
const double GAMMA_G       = 1.4;     // Specific heat ratio of the gas phase //
const double PI_S          = 100.0;   // Stiffness constant of the solid phase //
const double PI_G          = 0.0;     // Stiffness constant of the gas phase //
const double prs_floor     = 1.0e-12; // Pressure floor value //
const double rho_floor     = 1.0e-14; // Density floor value //
const double small_num     = 1.0e-12; // Effective small number in the code //
const int  dofs_per_cell   = 2;       // Number of degrees of freedom polynomial expansion in a cell //

typedef multi_array<double, 2> Matrix;
typedef multi_array<double, 1> Vector;

//----------------------------------------------------------------------------
// Various types of boundary conditions
//----------------------------------------------------------------------------

enum bndry_type{inflow, periodic, reflective, transmissive};

//----------------------------------------------------------------------------
// Structure defining various critical parameters controlling the simulation
//----------------------------------------------------------------------------

struct  AppCtx {
    double x_min;                      /* x-coordinate of the domain begining */
    double x_max;                      /* x-coordinate of the domain ending */
    double N_cells;                    /* No. of cells in the domain */
    double CFL;                        /* CFL condition, should be less than 1.0 */
    double InitialTime = 0.0;          /* Initial time of the simulation */
    double FinalTime;                  /* Final time of the simulation */
    bool reconstruct_primitive = true; /* Reconstruct primitive variables */
    enum bndry_type left_boundary;     /* Boundary condition on the left face */
    enum bndry_type right_boundary;    /* Boundary condition on the right face */
};

/* -------------------------------------------------------- Linear-Algebra Functions -------------------------------------------------------- */

//----------------------------------------------------------------------------
// Find the minimum value in a vector V of size N 
//----------------------------------------------------------------------------

double MinVal(const Vector& V, const int& N) {

    double min = V[0]; 
    
    for (int i = 1; i < N; ++i) {
        if (V[i] < min)
            min = V[i]; 
    }
    
    return min; 
}

//----------------------------------------------------------------------------
// Find the maximum value in a vector V of size N 
//---------------------------------------------------------------------------

double MaxVal(const Vector& V, const int& N) {

    double max = V[0]; 
    
    for (int i = 1; i < N; ++i) {
        if (V[i] > max)
            max = V[i]; 
    }
    
    return max; 
}

//----------------------------------------------------------------------------
// Given two vectors V1 and V2 of size N, find the L2 residue between them  
//---------------------------------------------------------------------------

double residue(const Vector& V1, const Vector& V2, const int& N) {
    
    double sum = 0.0;
    double diff; 
    
    for (int i = 0; i < N; ++i) {
        diff = V1[i] - V2[i];
        sum += diff*diff; 
    }
    
    return std::sqrt(sum); 
}

//----------------------------------------------------------------------------
// Matrix-Vector Multiplication
// Does the operation: y := A*x
// A -> (m x n) matrix
// x -> (n) vector
// y -> (m) vector
//----------------------------------------------------------------------------

void MatVecMult(const Matrix& A, const Vector& x, Vector& y, int m, int n) {
    for (int i = 0; i < m; i++ ) {
        y[i]= 0.0;
        for (int j = 0; j < n; j++ )
            y[i]+= A[i][j]*x[j];
    }
}

//----------------------------------------------------------------------------
// Matrix-Matrix Multiplication
// Does the operation: C := A*B
// A -> (n x m) matrix
// B -> (m x p) matrix
// C -> (n x p) matrix
//----------------------------------------------------------------------------

void MatMatMult(const Matrix& A, const Matrix& B, Matrix& C, int n, int m, int p) {
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < p; ++j) {
            C[i][j] = 0.0; 
            for (int k = 0; k < m; ++k) {
                C[i][j] += A[i][k]*B[k][j]; 
            }
        }
    }
}

/* -------------------------------------------------------- PDE Functions -------------------------------------------------------- */

//----------------------------------------------------------------------------
// Convert a conserved variable to primitive variable
//----------------------------------------------------------------------------

void Cons2Prim(const Vector& Q, Vector& V) {

    double phi[2]; double rho[2]; double vx[2]; double E[2]; double p[2];

    phi[0] = Q[6];
    phi[1] = 1.0 - Q[6];

    rho[0] = Q[0]/phi[0];
    rho[1] = Q[3]/phi[1];

    vx[0] = Q[1]/Q[0];
    vx[1] = Q[4]/Q[3];

    E[0] = Q[2]/phi[0];
    E[1] = Q[5]/phi[1];

    p[0] = (GAMMA_S - 1.0)*( E[0] - 0.5*rho[0]*vx[0]*vx[0] ) - GAMMA_S*PI_S;
    p[1] = (GAMMA_G - 1.0)*( E[1] - 0.5*rho[1]*vx[1]*vx[1] ) - GAMMA_G*PI_G;

    // Now assign them to vector V

    V[0] = rho[0]; // Density of solid phase
    V[1] = vx[0];  // velocity of solid phase
    V[2] = p[0];   // pressure of solid phase
    V[3] = rho[1]; // Density of gas phase
    V[4] = vx[1];  // velocity of gas phase
    V[5] = p[1];   // pressure of gas phase
    V[6] = phi[0]; // Volume fraction of solid phase
}

//----------------------------------------------------------------------------
// Convert a primitive variable to conserved variable
//----------------------------------------------------------------------------

void Prim2Cons(const Vector& V, Vector& Q) {

    Q[0] = V[6]*V[0]; // phi_s*rho_s
    Q[1] = Q[0]*V[1]; // phi_s*rho_s*vx_s
    Q[2] = Q[0]*( 0.5*V[1]*V[1] + ( V[2] + GAMMA_S*PI_S)/(V[0]*(GAMMA_S-1.0)) ); // phi_s*E_s

    Q[3] = (1.0-V[6])*V[3]; // phi_g*rho_g
    Q[4] = Q[3]*V[4];       // phi_g*rho_g*vx_g
    Q[5] = Q[3]*( 0.5*(V[4]*V[4]) + (V[5] + GAMMA_G*PI_G)/(V[3]*(GAMMA_G-1.0)) ); // phi_g*E_g

    Q[6] = V[6]; // phi_s
}

//----------------------------------------------------------------------------
// Find the conservative flux components F in the given normal direction
//----------------------------------------------------------------------------

double PDEFlux(const Vector& Q, const double& x, Vector& F) {

    double phi[2]; double rho[2]; double vx[2]; double E[2]; double p[2];

    phi[0] = Q[6];
    phi[1] = 1.0 - Q[6];

    rho[0] = Q[0]/phi[0];
    rho[1] = Q[3]/phi[1];

    vx[0] = Q[1]/Q[0];
    vx[1] = Q[4]/Q[3];

    E[0] = Q[2]/phi[0];
    E[1] = Q[5]/phi[1];

    p[0] = (GAMMA_S - 1.0)*( E[0] - 0.5*rho[0]*vx[0]*vx[0] ) - GAMMA_S*PI_S;
    p[1] = (GAMMA_G - 1.0)*( E[1] - 0.5*rho[1]*vx[1]*vx[1] ) - GAMMA_G*PI_G;

    // Check if the input state is physically admissible

    if (rho[0] < rho_floor) {
        std::cerr << "Negative density (solid) = " << rho[0] << std::endl;
        std::cerr << "At x = " << x << std::endl;
        std::exit(EXIT_FAILURE);
    }

    if (rho[1] < rho_floor) {
        std::cerr << "Negative density (gas) = " << rho[1] << std::endl;
        std::cerr << "At x = " << x << std::endl;
        std::exit(EXIT_FAILURE);
    }

    if ((p[0] + PI_S)  < prs_floor) {
        std::cerr << "Negative pressure (solid) = " << p[0] + PI_S << std::endl;
        std::cerr << "At x = " << x << std::endl;
        std::exit(EXIT_FAILURE);
    }

    if ((p[1] + PI_G)  < prs_floor) {
        std::cerr << "Negative pressure (gas) = " << p[1] + PI_G << std::endl;
        std::cerr << "At x = " << x << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Now find the fluxes

    F[0] = phi[0]*rho[0]*vx[0];
    F[1] = phi[0]*(rho[0]*vx[0]*vx[0] + p[0]);
    F[2] = phi[0]*vx[0]*(E[0] + p[0]);
    F[3] = phi[1]*rho[1]*vx[1];
    F[4] = phi[1]*(rho[1]*vx[1]*vx[1] + p[1]);
    F[5] = phi[1]*vx[1]*(E[1] + p[1]);
    F[6] = 0.0;

    // Return maximum eigen value
    
    double s_s = std::abs(vx[0]) + std::sqrt(GAMMA_S*(p[0] + PI_S)/rho[0]);
    double s_g = std::abs(vx[1]) + std::sqrt(GAMMA_G*(p[1] + PI_G)/rho[1]);

    return std::max(s_s, s_g);
}

//----------------------------------------------------------------------------
// Find the non-conservative flux components nF
//----------------------------------------------------------------------------

void PDENonConsFlux(const Vector& Q, const Vector& grad_Q_x, Vector& nF) {
    // First extract the primitive variables

    // u_I = u_s; p_I = p_g (standard approximation)

    double phi[2]; double rho[2]; double vx[2]; double E[2];

    phi[0] = Q[6];
    phi[1] = 1.0 - Q[6];

    rho[1] = Q[3]/phi[1];

    vx[0] = Q[1]/Q[0];
    vx[1] = Q[4]/Q[3];

    E[1] = Q[5]/phi[1];

    double p_I = (GAMMA_G - 1.0)*( E[1] - 0.5*rho[1]*(vx[1]*vx[1]) ) - GAMMA_G*PI_G;
    double u_I = vx[0];
    double phi_x = grad_Q_x[6];

    // Now find the fluxes

    nF[0] = 0.0;
    nF[1] = -p_I*phi_x;
    nF[2] = -p_I*u_I*phi_x;
    nF[3] = 0.0;
    nF[4] = p_I*phi_x;
    nF[5] = p_I*u_I*phi_x;
    nF[6] = u_I*phi_x;
}

//----------------------------------------------------------------------------
// For a given state Q, find all the eigenvalues L
//----------------------------------------------------------------------------

void PDEEigenvalues(const Vector& Q, Vector& L) {
    
    double phi[2]; double rho[2]; double vx[2]; double E[2]; double p[2];

    phi[0] = Q[6];
    phi[1] = 1.0 - Q[6];

    rho[0] = Q[0]/phi[0];
    rho[1] = Q[3]/phi[1];

    vx[0] = Q[1]/Q[0];
    vx[1] = Q[4]/Q[3];

    E[0] = Q[2]/phi[0];
    E[1] = Q[5]/phi[1];

    p[0] = (GAMMA_S - 1.0)*( E[0] - 0.5*rho[0]*vx[0]*vx[0] ) - GAMMA_S*PI_S;
    p[1] = (GAMMA_G - 1.0)*( E[1] - 0.5*rho[1]*vx[1]*vx[1] ) - GAMMA_G*PI_G;
    
    double a_s = std::sqrt(GAMMA_S*(p[0] + PI_S)/rho[0]);
    double a_g = std::sqrt(GAMMA_G*(p[1] + PI_G)/rho[1]);
    
    L[0] = vx[1];       // u_g
    L[1] = vx[0];       // u_s        
    L[2] = vx[0];       // u_s 
    L[3] = vx[1] - a_g; // u_g - a_g
    L[4] = vx[1] + a_g; // u_g + a_g
    L[5] = vx[0] - a_s; // u_s - a_s 
    L[6] = vx[0] + a_s; // u_s + a_s
}

//----------------------------------------------------------------------------
// Form the non-conservative matrix B
//----------------------------------------------------------------------------

void PDEMatrixB(const Vector& Q, Matrix& B) {

    // First extract the primitive variables

    // u_I = u_s; p_I = p_g (standard approximation)


    double phi[2]; double rho[2]; double vx[2]; double E[2];

    phi[0] = Q[6];
    phi[1] = 1.0 - Q[6];

    rho[1] = Q[3]/phi[1];

    vx[0] = Q[1]/Q[0];
    vx[1] = Q[4]/Q[3];

    E[1] = Q[5]/phi[1];

    double p_I = (GAMMA_G - 1.0)*( E[1] - 0.5*rho[1]*(vx[1]*vx[1]) ) - GAMMA_G*PI_G;
    double u_I = vx[0];

    for (int i = 0; i < nVar; ++i)
        for (int j = 0; j < nVar; ++j)
            B[i][j] = 0.0;

    B[1][6] = -p_I;
    B[2][6] = -p_I*u_I;
    B[4][6] =  p_I;
    B[5][6] =  p_I*u_I;
    B[6][6] =  u_I;
}

//----------------------------------------------------------------------------
// The Roe matrix B̃(Ua,Ub) between two generic states Qa and Qb is numerically
// via Gaussian quadrature, by the function RoeMatrix
//----------------------------------------------------------------------------

void RoeMatrix(const Vector& Qa, const Vector& Qb, Matrix& BRoe) {

    // 3 - point quadrature points in [0,1] for path integrals

    const int N_gps = 3;
    const double s_gp[] = {0.1127016653792583, 0.5, 0.8872983346207417};
    const double w_gp[] = {0.2777777777777778, 0.4444444444444444, 0.2777777777777778};
    
    Matrix B(extents[nVar][nVar]);
    Vector Q(extents[nVar]);

    // First make the BRoe matrix zero

    for (int i = 0; i < nVar; ++i)
        for (int j = 0; j < nVar; ++j)
            BRoe[i][j] = 0.0;

    for (int q = 0; q < N_gps; ++q) {

        for (int c = 0; c < nVar; ++c)
            Q[c] = Qa[c] + s_gp[q]*(Qb[c] - Qa[c]);

        PDEMatrixB(Q, B);

        for (int i = 0; i < nVar; ++i) {
            for (int j = 0; j < nVar; ++j) {
                BRoe[i][j] += w_gp[q]*B[i][j];
            }
        }
    }
}

//----------------------------------------------------------------------------
// For a given state Q, find the intermediate states 
//----------------------------------------------------------------------------

void PDEIntermediateFields(const Vector& Q, Matrix& Lambda, Matrix& RS, Matrix& LS) {
    
    int i, j; 
    double phi_s, phi_g, rho_s, rho_g, u_s, u_g, E_s, E_g, p_s, p_g; // Primitive variables
    Matrix DQ_DV(extents[nVar][nVar]); Matrix DV_DQ(extents[nVar][nVar]);      // Jacobian matrices
    Matrix RS_p(extents[nVar][nLin]), LS_p(extents[nLin][nVar]);               // Right and left eigenvectors in primitive variables 
    double tempa, tempb, tempc; 
    
    // First extract states 
    
    phi_s = Q[6];
    phi_g = 1.0 - Q[6];

    rho_s = Q[0]/phi_s;
    rho_g = Q[3]/phi_g;

    u_s = Q[1]/Q[0];
    u_g = Q[4]/Q[3];

    E_s = Q[2]/phi_s;
    E_g = Q[5]/phi_g;

    p_s = (GAMMA_S - 1.0)*( E_s - 0.5*rho_s*u_s*u_s ) - GAMMA_S*PI_S;
    p_g = (GAMMA_G - 1.0)*( E_g - 0.5*rho_g*u_g*u_g ) - GAMMA_G*PI_G;
    
    // First fill the diagonal intermediate characteristic fields matrix Lambda
    
    for (i = 0; i < nLin; ++i) 
        for (j = 0; j < nLin; ++j)
            Lambda[i][j] = 0.0; 
    
    Lambda[0][0] = u_g; Lambda[1][1] = u_s; Lambda[2][2] = u_s;
    
    // Fill the Jacobian matrices 
    
    for (i = 0; i < nVar; ++i) { 
        for (j = 0; j < nVar; ++j) {
            DQ_DV[i][j] = 0.0; DV_DQ[i][j] = 0.0; 
        }
    }
    
    //-------------------------------------------------------------------------------------------------------
    
    DQ_DV[0][0] = phi_s; 
    DQ_DV[0][6] = rho_s; // Row 0 
    
    DQ_DV[1][0] = phi_s*u_s; 
    DQ_DV[1][1] = phi_s*rho_s; 
    DQ_DV[1][6] = rho_s*u_s;  // Row 1
    
    DQ_DV[2][0] = 0.5*phi_s*u_s*u_s; 
    DQ_DV[2][1] = phi_s*rho_s*u_s; 
    DQ_DV[2][2] = phi_s/(GAMMA_S - 1.0);
    DQ_DV[2][6] = (GAMMA_S*PI_S + p_s + 0.5*rho_s*u_s*u_s*(GAMMA_S - 1.0))/(GAMMA_S - 1.0); // Row 2
    
    DQ_DV[3][3] = phi_g; 
    DQ_DV[3][6] = -rho_g; // Row 3 
    
    DQ_DV[4][3] = u_g*phi_g;
    DQ_DV[4][4] = rho_g*phi_g; 
    DQ_DV[4][6] = -rho_g*u_g; // Row 4 
    
    DQ_DV[5][3] = 0.5*u_g*u_g*phi_g;
    DQ_DV[5][4] = rho_g*u_g*phi_g;
    DQ_DV[5][5] = phi_g/(GAMMA_G - 1.0); 
    DQ_DV[5][6] = -(GAMMA_G*PI_G + p_g + 0.5*rho_g*u_g*u_g*(GAMMA_G - 1.0))/(GAMMA_G - 1.0); // Row 5
    
    DQ_DV[6][6] = 1.0; // Row 6 
     
    //-------------------------------------------------------------------------------------------------------

    DV_DQ[0][0] = 1.0/phi_s;
    DV_DQ[0][6] = -rho_s/phi_s; // Row 0
    
    DV_DQ[1][0] = -u_s/(phi_s*rho_s);
    DV_DQ[1][1] = 1.0/(phi_s*rho_s); // Row 1
    
    DV_DQ[2][0] = u_s*u_s*(GAMMA_S - 1.0)/(2.0*phi_s);
    DV_DQ[2][1] = -u_s*(GAMMA_S - 1.0)/phi_s;
    DV_DQ[2][2] = (GAMMA_S - 1.0)/phi_s;
    DV_DQ[2][6] = -(GAMMA_S*PI_S + p_s)/phi_s; // Row 2
    
    DV_DQ[3][3] =  1.0/phi_g;
    DV_DQ[3][6] =  rho_g/phi_g; // Row 3 
    
    DV_DQ[4][3] = -u_g/(rho_g*phi_g);
    DV_DQ[4][4] = 1.0/(rho_g*phi_g); 
    
    DV_DQ[5][3] = u_g*u_g*(GAMMA_G - 1.0)/(2.0*phi_g);
    DV_DQ[5][4] = -u_g*(GAMMA_G - 1.0)/(phi_g);
    DV_DQ[5][5] = (GAMMA_G - 1.0)/(phi_g); 
    DV_DQ[5][6] = (GAMMA_G*PI_G + p_g)/(phi_g);
    
    DV_DQ[6][6] = 1.0; 
    
    //-------------------------------------------------------------------------------------------------------
    
    // Form the right eigenvectors in primitive variables 
    
    for (i = 0; i < nVar; ++i) 
        for (j = 0; j < nLin; ++j) 
            RS_p[i][j] = 0.0;  
    
    // First eigenvector corresponding to u_g (contact discontinuity)
        
    RS_p[3][0] = 1.0; 
    
    // Second eigenvector corresponding to u_s (contact discontinuity)
    
    RS_p[0][1] = 1.0; 
    
    // Third eigenvector corresponding to u_s (material interface)
    
    tempa = p_g + PI_G;
    tempb = u_g - u_s; 
    tempc = GAMMA_G*tempa - rho_g*tempb*tempb; 
    
    RS_p[2][2] = (p_g - p_s)/phi_s;
    RS_p[3][2] = -rho_g*(GAMMA_G*tempa - tempc)/(phi_g*tempc);
    RS_p[4][2] =  GAMMA_G*(tempa*tempb)/(phi_g*tempc);
    RS_p[5][2] = -GAMMA_G*rho_g*tempa*tempb*tempb/(phi_g*tempc);
    RS_p[6][2] = 1.0; 
    
    //-------------------------------------------------------------------------------------------------------
    
    // Form the left eigenvectors in primitive variables
    
    for (i = 0; i < nLin; ++i) 
        for (j = 0; j < nVar; ++j) 
            LS_p[i][j] = 0.0;
    
    // First eigenvector corresponding to u_g (contact discontinuity)
        
    LS_p[0][3] = 1.0; 
    LS_p[0][5] = -rho_g/(GAMMA_G*tempa);
    
    // Second eigenvector corresponding to u_s (contact discontinuity)
    
    LS_p[1][0] = 1.0; 
    LS_p[1][2] = -rho_s/(GAMMA_S*(p_s + PI_S));
    LS_p[1][6] = rho_s*(p_g - p_s)/(GAMMA_S*phi_s*(p_s + PI_S));
    
    // Third eigenvector corresponding to u_s (material interface)
    
    LS_p[2][6] = 1.0; 
    
    //-------------------------------------------------------------------------------------------------------
    
    // Now obtain the eigenvector matrices in conserved variables 
    
    MatMatMult(LS_p, DV_DQ, LS, nLin, nVar, nVar);
    MatMatMult(DQ_DV, RS_p, RS, nVar, nVar, nLin);
}

//----------------------------------------------------------------------------
// HLLEM Non-Conservative Flux 
//----------------------------------------------------------------------------

double HLLEMNC(const Vector& QL, const Vector& QR, const double& x, Vector& Dm, Vector& Dp) {
        
    // Declare variables 
    
    int i, j, iHLL;
    const int MaxIter = 25; 
    Vector QM(extents[nVar]), FL(extents[nVar]), FR(extents[nVar]), FM(extents[nVar]); 
    Vector LL(extents[nVar]), LR(extents[nVar]), LM(extents[nVar]); 
    Vector PathInt(extents[nVar]), QHLL(extents[nVar]), QOld(extents[nVar]);
    Vector AD(extents[nVar]), Q_jump(extents[nVar]), Work1(extents[nVar]), Work2(extents[nVar]), Work3(extents[nLin]), Work4(extents[nLin]); 
    Matrix BRoe1(extents[nVar][nVar]), BRoe2(extents[nVar][nVar]);
    Matrix RS(extents[nVar][nLin]), LS(extents[nLin][nVar]), Lambda(extents[nLin][nLin]), Id(extents[nLin][nLin]);
    Matrix Delta(extents[nLin][nLin]), Lap(extents[nLin][nLin]), Lam(extents[nLin][nLin]);
    double sL, sR, smaxL, smaxR, res;
    const double tol = 1.0e-12;
    const double flattener = 1.0; 
    
    // Compute the intermediate state
    
    for (i = 0; i < nVar; ++i)
        QM[i] = 0.5*(QL[i] + QR[i]);
    
    // Compute the fluxes FL and FR 
    
    smaxL = PDEFlux(QL, x, FL);
    smaxR = PDEFlux(QR, x, FR);
    
    // Compute the eigenvalues of QM, QL and QR 

    PDEEigenvalues(QM, LM);
    PDEEigenvalues(QL, LL);
    PDEEigenvalues(QR, LR);
    
    // Compute the left and right  signal speed 
    
    sL = std::min(0.0, std::min(MinVal(LM, nVar), MinVal(LL, nVar)) ); 
    sR = std::max(0.0, std::max(MaxVal(LM, nVar), MaxVal(LR, nVar)) );
    
    // Compute the Roe matrix 
    
    RoeMatrix(QL, QR, BRoe1);

    // Compute the path integral from QL to QR

    for (i = 0; i < nVar; ++i)
        Q_jump[i] = QR[i] - QL[i]; 
    
    MatVecMult(BRoe1, Q_jump, Work1, nVar, nVar);
    
    for (i = 0; i < nVar; ++i) 
        PathInt[i] = Work1[i] + FR[i] - FL[i];
    
    // Initial guess 
        
    for (i = 0; i < nVar; ++i)
        QHLL[i] = ( QR[i]*sR - QL[i]*sL - PathInt[i] ) / (sR-sL); 
    
    // Start non-linear Picard iterations 
    
    for (iHLL = 0; iHLL < MaxIter; ++iHLL) {
    
        // Save old state 
        
        for (i = 0; i < nVar; ++i)
            QOld[i] = QHLL[i];
        
        RoeMatrix(QL, QHLL, BRoe1); // Compute the Roe matrix between QL and Qstar
        RoeMatrix(QHLL, QR, BRoe2); // Compute the Roe matrix between Qstar and QR
        
        for (i = 0; i < nVar; ++i)
            Q_jump[i] = QHLL[i] - QL[i]; 
    
        MatVecMult(BRoe1, Q_jump, Work1, nVar, nVar);
        
        for (i = 0; i < nVar; ++i)
            Q_jump[i] = QR[i] - QHLL[i]; 
    
        MatVecMult(BRoe2, Q_jump, Work2, nVar, nVar);
        
        for (i = 0; i < nVar; ++i) 
            PathInt[i] = Work1[i] + Work2[i] + FR[i] - FL[i];
            
        // Update the HLL state
        
        for (i = 0; i < nVar; ++i)
            QHLL[i] = ( QR[i]*sR - QL[i]*sL - PathInt[i] ) / (sR-sL); 
            
        res = residue(QHLL, QOld, nVar); 
        
        if (res < tol)
            break; 
    }
    
    if (res > tol) {
        std::cerr << "HLLEMNC solver: Non-linear iterations did not converge. Residue = " << res << std::endl; 
        std::exit(EXIT_FAILURE); 
    }
    
    for (i = 0; i < nVar; ++i) {
        
        Dm[i] = -sL/(sR-sL)*PathInt[i] + sL*sR/(sR-sL)*(QR[i]-QL[i]);
        Dp[i] = +sR/(sR-sL)*PathInt[i] - sL*sR/(sR-sL)*(QR[i]-QL[i]);
    
    }
    
    //-------------------------------------------------------------------------------------------------------
    // Add the anti-diffusive term of HLLEM-RS
    
    // Form identity matrix 
    
    for (i = 0; i < nLin; ++i) {
        for (j = 0; j < nLin; ++j) {
            if (i == j)
                Id[i][j] = 1.0; 
            else
                Id[i][j] = 0.0;
        }
    }
     
    // Compute intermediate eigenvalues and eigenvectors 

    PDEIntermediateFields(QM,Lambda,RS,LS);
    
    // Compute Lambda+ and Lambda- 
    
    for (i = 0; i < nLin; ++i) {
        for (j = 0; j < nLin; ++j) {
            Lap[i][j] = 0.5*(Lambda[i][j] + std::abs(Lambda[i][j]));
            Lam[i][j] = 0.5*(Lambda[i][j] - std::abs(Lambda[i][j]));
        }
    }
    
    // Compute Delta*
    
    for (i = 0; i < nLin; ++i) {
        for (j = 0; j < nLin; ++j) {
            Delta[i][j] = Id[i][j] - Lam[i][j]/(sL - 1.0e-14) - Lap[i][j]/(sR + 1.0e-14);
        }
    }
    
    for (i = 0; i < nVar; ++i)
        Q_jump[i] = QR[i] - QL[i]; 
    
    MatVecMult(LS, Q_jump, Work3, nLin, nVar);
    MatVecMult(Delta, Work3, Work4, nLin, nLin);
    MatVecMult(RS, Work4, AD, nVar, nLin);
    
    for (i = 0; i < nVar; ++i)
        AD[i] = (sR*sL)/(sR - sL)*AD[i]; 

    for (i = 0; i < nVar; ++i) {
        Dm[i] = Dm[i] - flattener*AD[i]; 
        Dp[i] = Dp[i] + flattener*AD[i];
    }
    
    return std::max(smaxL, smaxR); 
}

//----------------------------------------------------------------------------
// LLF (Rusanov) Non-Conservative Flux 
//----------------------------------------------------------------------------

double LLFNC(const Vector& QL, const Vector& QR, const double& x, Vector& F, Vector& D) {

    Vector FL(extents[nVar]), FR(extents[nVar]), Q_jump(extents[nVar]);
    Matrix B(extents[nVar][nVar]);

    double s_max_l = PDEFlux(QL, x, FL);
    double s_max_r = PDEFlux(QR, x, FR);

    double s_max = std::max(s_max_l, s_max_r);

    for (int c = 0; c < nVar; ++c) {
        Q_jump[c] = QR[c] - QL[c];
        F[c] = 0.5*(FR[c] + FL[c] - s_max*(Q_jump[c]));
    }

    RoeMatrix(QL, QR, B);
    
    // Multiply B with the jump term

    MatVecMult(B, Q_jump, D, nVar, nVar);

    return s_max;
}

//----------------------------------------------------------------------------
// Minmod limiter function 
//----------------------------------------------------------------------------

double sgn(double a) {
    if (a > 0.0) return 1.0;
    else if ( a == 0.0) return 0.0;
    else return -1.0;
}

double minmod_limiter(double a, double b) {
    return 0.5*(sgn(a)+ sgn(b))*std::min(std::abs(a), std::abs(b));
}

//----------------------------------------------------------------------------
// Initial Condition Function
//----------------------------------------------------------------------------

Vector initial_condition(double x) {

    Vector V(extents[nVar]);

    if (x <=  0.0) {
        V[0] = 800.0;
        V[1] = 1.0;
        V[2] = 1.0;
        V[3] = 2.0;
        V[4] = 1.0;
        V[5] = 1.0;
        V[6] = 0.99;
    }

    else {

        V[0] = 1000.0;
        V[1] = 1.0;
        V[2] = 1.0;
        V[3] = 1.0;
        V[4] = 1.0;
        V[5] = 1.0;
        V[6] = 0.01;
    }

    return V; 
}

//----------------------------------------------------------------------------
// Path Conservative Finite Volume Method Class
//----------------------------------------------------------------------------

class HyPE_1D {

    // Typedefs

    typedef boost::multi_array<Vector, 1> array_type;
    typedef array_type::index index;

    boost::multi_array<std::vector<double>, 2> U;  // Conserved variables at cells 
    boost::multi_array<Vector, 1> Dp;              // Non-Conservative fluxes at cell faces
    boost::multi_array<Vector, 1> Dm;              // Non-Conservative fluxes at cell faces
    boost::multi_array<Vector, 1> RHS;             // RHS term for each face 
    boost::multi_array<Vector, 1> U_L;             // Value of conserved variable at cell left face 
    boost::multi_array<Vector, 1> U_R;             // Value of conserved variable at cell right face 
    boost::multi_array<Vector, 1> F_L;             // Value of conserved variable at cell left face 
    boost::multi_array<Vector, 1> F_R;             // Value of conserved variable at cell right face 
    boost::multi_array<double, 1> x;               // Cell centers

    AppCtx Params;
    int N_ph;
    double dx;
    double dt;
    double time;
    int time_step;
    int rk_stage;

    void initialize();
    void apply_boundary_conditions();
    void limit_solution();
    void compute_rhs(double);
    void solve();
    void plot() const;

public:
    HyPE_1D(AppCtx);
    void run();
};

//----------------------------------------------------------------------------
// Constructor
//----------------------------------------------------------------------------

HyPE_1D::HyPE_1D(AppCtx params) :
    U(boost::extents[params.N_cells + 6][nVar]),
    Dp(boost::extents[params.N_cells + 1]),
    Dm(boost::extents[params.N_cells + 1]),
    RHS(boost::extents[params.N_cells]),
    U_L(boost::extents[params.N_cells + 2]),
    U_R(boost::extents[params.N_cells + 2]),
    F_L(boost::extents[params.N_cells + 2]),
    F_R(boost::extents[params.N_cells + 2]),
    x(boost::extents[params.N_cells + 6]),
    Params(params),
    N_ph(3),
    dt(0.0),
    time(0.0),
    time_step(0),
    rk_stage(0)
    {

    // Initialize the grid

    dx = (Params.x_max - Params.x_min)/static_cast<double>(Params.N_cells);
    boost::array<array_type::index, 1> bases_1d = {{-N_ph}};
    x.reindex(bases_1d);

    for (int i = -N_ph; i < Params.N_cells + N_ph; i++)
        x[i] = Params.x_min + ((i+1)-0.5)*dx;

    boost::array<array_type::index, 2> bases_2d = {{-N_ph, 0}};

    U.reindex(bases_2d);


    for (int i = -3; i < Params.N_cells + 3; ++i)
        for (int c = 0; c < nVar; ++c)
                U[i][c].resize(dofs_per_cell, 0.0);

    boost::array<array_type::index, 1> bases_1d_face = {{-1}};

    U_L.reindex(bases_1d_face);
    U_R.reindex(bases_1d_face);
    F_L.reindex(bases_1d_face);
    F_R.reindex(bases_1d_face);
    
    for (int i = 0; i < Params.N_cells; ++i) {
        
        RHS[i].resize(extents[nVar]);
    }
    
    for (int i = -1; i < Params.N_cells+1; ++i) {
        
        U_L[i].resize(extents[nVar]);
        U_R[i].resize(extents[nVar]);
        F_L[i].resize(extents[nVar]);
        F_R[i].resize(extents[nVar]);
        
    }
    
    for (int i = 0; i < Params.N_cells + 1; ++i) {
            
        Dp[i].resize(extents[nVar]);
        Dm[i].resize(extents[nVar]);
        
    }
}

//----------------------------------------------------------------------------
// Initialize the solution using midpoint rule 
//----------------------------------------------------------------------------

void HyPE_1D::initialize() {
    
    std::cout << "Initializing the solution" << std::endl;

    Vector Q(extents[nVar]), V(extents[nVar]);

    // Loop through all the cells

    for (int i = 0; i < Params.N_cells; i++) {

        V = initial_condition(x[i]);

        Prim2Cons(V, Q);
        
        for (int c = 0; c < nVar; ++c)
            U[i][c][0] = Q[c];
    }

    std::cout << "Done!" << std::endl;
}

//----------------------------------------------------------------------------
// Apply boundary conditions
//----------------------------------------------------------------------------

void HyPE_1D::apply_boundary_conditions() {

    int oned_begin, oned_end, ilhs, irhs;

    // ---------------------- Left boundary ----------------------

    oned_begin = 0; oned_end = Params.N_cells-1;

    for (int i = 0; i < N_ph; ++i) {

        // Outflow/Transmissive boundary

        if (Params.left_boundary == transmissive) {

            ilhs = oned_begin - i - 1;
            irhs = oned_begin;

            for (int c = 0; c < nVar; ++c) 
                U[ilhs][c][0] = U[irhs][c][0];
            
        }

        if (Params.left_boundary == reflective) {

            ilhs = oned_begin - i - 1;
            irhs = oned_begin + i;

            for (int c = 0; c < nVar; ++c) 
                U[ilhs][c][0] = U[irhs][c][0];

            U[ilhs][1][0] = -U[ilhs][1][0];
            U[ilhs][4][0] = -U[ilhs][4][0];
        }

        if (Params.left_boundary == periodic) {
            ilhs = oned_begin - i - 1;
            irhs = oned_end - i;

            for (int c = 0; c < nVar; ++c)
                U[ilhs][c][0] = U[irhs][c][0];
            
        }
    }

    // ---------------------- Right boundary ----------------------

    for (int i = 0; i < N_ph; ++i) {

        // Outflow/Transmissive boundary

        if (Params.left_boundary == transmissive) {

            ilhs = oned_end + i + 1;
            irhs = oned_end;

            for (int c = 0; c < nVar; ++c) 
                U[ilhs][c][0] = U[irhs][c][0];
        }

        if (Params.left_boundary == reflective) {

            ilhs = oned_end + i + 1;
            irhs = oned_end - i;

            for (int c = 0; c < nVar; ++c) 
                U[ilhs][c][0] = U[irhs][c][0];

            U[ilhs][1][0] = -U[ilhs][1][0];
            U[ilhs][4][0] = -U[ilhs][4][0];
        }

        if (Params.left_boundary == periodic) {

            ilhs = oned_end + i + 1;
            irhs = oned_begin + i;

            for (int c = 0; c < nVar; ++c) 
                U[ilhs][c][0] = U[irhs][c][0];
        }
    }
}

//----------------------------------------------------------------------------
// Find the RHS in each cell
//----------------------------------------------------------------------------

void HyPE_1D::compute_rhs(double t) {

    double s, left_slope, right_slope, s_max = 0.0;
    const double r1_dx = 1./dx;
    Vector Q(extents[nVar]);
    Vector grad_Q(extents[nVar]);
    Vector nF(extents[nVar]);
    Vector Q_L(extents[nVar]); Vector Q_R(extents[nVar]);
    
    apply_boundary_conditions();

    for (int i = -1; i < Params.N_cells+1; ++i) {

        for (int c = 0; c < nVar; ++c) {

            left_slope  = U[i][c][0] - U[i-1][c][0];
            right_slope = U[i+1][c][0] - U[i][c][0];
            
            U[i][c][1] = minmod_limiter(left_slope, right_slope);
            
            U_L[i][c] = U[i][c][0] - 0.5*U[i][c][1];
            U_R[i][c] = U[i][c][0] + 0.5*U[i][c][1];
        }
        
        PDEFlux(U_L[i], x[i]-0.5*dx, F_L[i]); 
        PDEFlux(U_R[i], x[i]+0.5*dx, F_R[i]); 
    }

    // Find upwind flux

    for (int i = 0; i < Params.N_cells + 1; ++i) {

        for (int c = 0; c < nVar; ++c) {
            Q_L[c] = U_R[i-1][c];
            Q_R[c] = U_L[i][c];
        }
    
        s = HLLEMNC(Q_L, Q_R, x[i]-0.5*dx, Dm[i], Dp[i]); 

        if (s > s_max)
            s_max = s;
    }

    // Find RHS

    for (int i = 0; i < Params.N_cells; ++i) {

        for (int c = 0; c < nVar; ++c) {
            RHS[i][c] = -r1_dx*(F_R[i][c] - F_L[i][c] +
                                     Dm[i+1][c] + Dp[i][c]);
        }
    }

    // Add the smooth part of the non-conservative term

    for (int i = 0; i < Params.N_cells; ++i) {

        for (int c = 0; c < nVar; ++c) {
            Q[c] = U[i][c][0];
            grad_Q[c] = r1_dx*U[i][c][1]; 
        }

        PDENonConsFlux(Q, grad_Q, nF);

        for (int c = 0; c < nVar; ++c)
            RHS[i][c] += -1.0*nF[c];
    }

    // Find the time step size (only if rk_stage = 1)

    if (rk_stage == 1) {
        dt = Params.CFL*dx/s_max;

        // If time step exceeds the final time, reduce it accordingly

        if((time + dt)>Params.FinalTime)
            dt = Params.FinalTime - time;
    }
}

//----------------------------------------------------------------------------
// Update Solution using SSPRK22 method
//----------------------------------------------------------------------------

void HyPE_1D::solve() {

    std::cout << "Solving using SSPRK (2,2) method" << std::endl;

    boost::multi_array<double, 2> U_old(extents[Params.N_cells][nVar]);

    while (time < Params.FinalTime) {

        printf ("time = %4.3e, dt = %4.3e, final time = %4.3e\n", time, dt, Params.FinalTime);

        //  Stage 1

        rk_stage = 1;
        compute_rhs(time);

        for (int i = 0; i < Params.N_cells; ++i){

            for (int c = 0; c < nVar; ++c) {

                U_old[i][c] =   U[i][c][0];

                U[i][c][0] +=   dt*RHS[i][c];

            }
        }
        
        //  Stage 2

        rk_stage = 2;
        compute_rhs(time + dt);


        for (int i = 0; i < Params.N_cells; ++i) {

            for (int c = 0; c < nVar; ++c)
                U[i][c][0] =   0.5*(U_old[i][c] +   U[i][c][0] + dt*RHS[i][c]);

        }

        time += dt;
        time_step++;
    }

    printf ("time = %4.3e, dt = %4.3e, final time = %4.3e\n", time, dt, Params.FinalTime);
}

//----------------------------------------------------------------------------
// Plot solution as csv file
//----------------------------------------------------------------------------

void HyPE_1D::plot() const {

    Vector V(extents[nVar]), Q(extents[nVar]);

    std::ofstream out_data;
    const std::string filename = "sol.csv";
    out_data.open (filename);
    out_data.flags( std::ios::dec | std::ios::scientific );
    out_data.precision(6);

    out_data << "x, RHO_s, V_s, P_s, RHO_g, V_g, P_g, PHI_s" << std::endl;

    for (int i = 0; i < Params.N_cells; ++i) {

        for (int c = 0; c < nVar; ++c) {
            Q[c] = U[i][c][0];
        }

        Cons2Prim(Q, V);

        out_data << x[i] << ",";

        for (int c = 0; c < nVar; ++c) {

            if (c == nVar -1)
                out_data << V[c];

            else
                out_data << V[c] << ",";
        }

        out_data << std::endl;

    }

    out_data.close();
}

//----------------------------------------------------------------------------
// Put everything together and run the problem
//----------------------------------------------------------------------------

void HyPE_1D::run() {
    
    auto start = std::chrono::system_clock::now();

    //-------------------------------------
    
    initialize();
    solve();
    plot();
    
    //-------------------------------------
    
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    std::cout << "No. of time steps = " << time_step << std::endl; 
    std::cout << "Time taken = " << elapsed_seconds.count() << std::endl;
}

//----------------------------------------------------------------------------
// Main function of the code 
//----------------------------------------------------------------------------

int main() {

    AppCtx Params;

    Params.x_min = -0.5;
    Params.x_max =  0.5;
    Params.CFL   = 0.7;
    Params.InitialTime = 0.0;
    Params.FinalTime = 0.25;
    Params.N_cells = 400;
    Params.left_boundary  = transmissive;
    Params.right_boundary = transmissive;

    HyPE_1D Baer_Nunziato(Params);

    Baer_Nunziato.run();
    
    return 0;
}
