#ifndef _DG1D
#define _DG1D 1

#include <stdbool.h>

#define LINCON 1
#define BURGER 2
#define EULER 3

// Available test cases
#define SOD 10
#define MSOD 11
#define SSOD 12
#define BLAST 20
#define LAX 30
#define LOWD 40
#define SHUOSHER 50
#define PULSE 60
#define POLY 70
#define CONTACT 80
#define NOH 90
#define RELBLAST 100
// Type of bc
#define FREE 0
#define FIXED 1
#define PERIODIC 2

REAL mass0, mass1[2][2], mass2[3][3], mass3[4][4], mass4[5][5], mass5[6][6];
REAL cfl, cfl_f, cfl_u, dt, finaltime;
REAL dx_init;
REAL XS; /* Shock position */
REAL xmin0, xmax0;
REAL xmin, xmax;
REAL d_left, u_left, u1_left, u2_left, p_left;
REAL d_right, u_right, u1_right, u2_right, p_right;
REAL Mfact;
UINT ALE;
REAL dxmin, dxmax;
UINT test_case;
UINT bc_left, bc_right;
UINT pos_lim;
UINT tvd_lim;
UINT predictor_method;
UINT error_order;
UINT adaptation;
REAL hmin;
REAL hmax;
UINT siter;
REAL chi;
REAL eta;
UINT n_extra_cells;
UINT free_cell_index;
UINT free_face_index;

// We list here the type of predictors used.
enum { taylor_1p = 0, cerk_nodal_2p = 12, cerk_nodal_3p = 13 };

struct CELL {
    int do_adapt;
    int level;
    REAL x, xl, xr, wl, wr, h, *xg;
    UINT p, ng, ngll;
    REAL **U, **Re;
    struct FACE *lface, *rface;
    struct CELL *lcell, *rcell;
    bool active;
};
typedef struct CELL CELL;

struct FACE {
    REAL x; // location of face
    REAL w; // velocity of face
    CELL *lcell, *rcell;
    bool active;
};
typedef struct FACE FACE;

CELL *Init();
FACE *InitFaces(CELL *);
void TimeStep(CELL *);
void SaveSol(CELL *);
void Flux(CELL *cell, FACE *face, const int predictor);
void Update(CELL *, FACE *);
void ApplyLimiter(CELL *);
void Result(CELL *, double time, int iter, int final);
void MeshVel(FACE *);

void GaussInit();
void GaussPoints(CELL *);
void GetGaussPoints(double xl, double xr, int ng, double xgauss[ng]);
REAL ShapeFun(REAL, CELL *, UINT);
REAL ShapeFunDeriv(REAL, CELL *, UINT);

void UatGLL(CELL *cell, REAL **U);
void Uvect(CELL *cell, REAL x, REAL *U);
void EulerFlux(REAL *U, REAL w, REAL *flux);
void Jacobian(REAL *U, REAL A[][NVAR]);
void RoeFlux(REAL *Ul, REAL *Ur, REAL w, REAL *flux);
void ALERoeFlux(REAL *Ul, REAL *Ur, REAL w, REAL *flux);
void LFFlux(REAL *Ul, REAL *Ur, REAL w, REAL *flux);
void ECUSPFlux(REAL *Ul, REAL *Ur, REAL *flux);
void HLLCFlux(double *Ul, double *Ur, double w, double *flux);
void HLLCPSFlux(double *Ul, double *Ur, double w, double *flux);
void HLLFlux(double *Ul, double *Ur, double w, double *flux);
void AUSMDVFlux(REAL *Ul, REAL *Ur, REAL *flux);
void LFCFlux(REAL *Ul, REAL *Ur, REAL *flux);
double sound_speed(REAL *U);
void initialize_cell_points(
    CELL *cell, double xl, double xr, double wl, double wr);

void EigMat(REAL *, REAL[][NVAR], REAL[][NVAR]);
void Multi(REAL[][NVAR], REAL *);
void MultiplyMatrix(REAL R1[][NVAR], REAL R2[][NVAR], REAL A[NVAR][NVAR]);

void MoveGrid(CELL *, FACE *);
REAL GetMeshVel(CELL *cell, REAL x);

// Initial Conditions
void InitCondShocktube(REAL x, REAL *U);
void InitCondBlast(REAL x, REAL *U);
void InitCondShuOsher(REAL x, REAL *U);
void InitCondPulse(REAL x, REAL *U);
void InitCondPoly(REAL x, REAL *V);

// Error Analysis
void l2error(CELL *cell, double time, double error[NVAR], int order);

// Predictors Declarations
void get_predictor(int predictor,
                   double tstep,
                   int n_time_points,
                   double t[n_time_points],
                   CELL *cell_p,
                   int nodes,
                   double x_nodes[nodes],
                   double u_p[nodes][n_time_points][NVAR]);

void taylor_1p_predictor(double tstep,
                         int n_time_points,
                         double t[n_time_points],
                         CELL *cell_p,
                         int nodes,
                         double x_nodes[nodes],
                         double u_p[nodes][n_time_points][NVAR]);

void cerk_nodal_2p_predictor(double tstep,
                             int n_time_points,
                             double t[n_time_points],
                             CELL *cell_p,
                             int nodes,
                             double x_nodes[nodes],
                             double u_p[nodes][n_time_points][NVAR]);

void cerk_nodal_3p_predictor(double tstep,
                             int n_time_points,
                             double t[n_time_points],
                             CELL *cell_p,
                             int nodes,
                             double x_nodes[nodes],
                             double u_p[nodes][n_time_points][NVAR]);

void nodal_to_modal(int modes,
                    double nodal[modes][NVAR],
                    double modal[modes][NVAR]);

// Cell Adaptation Function
void adapt(int nc, CELL cell[nc], int nf, FACE face[nf]);

// For converting Conserevd variables to Primitives using Newton-Raphson
void func(double conserved[], double f[2], double p);
double InitialGuess(double conserved[]);
double newton(double conserved[]);
void con2prim(double conserved[], double primitive[]);
void prim2con(double primitive[], double conserved[]);
#endif
