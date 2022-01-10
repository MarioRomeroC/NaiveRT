#ifndef DEFINE_H_INCLUDED
#define DEFINE_H_INCLUDED

//SUBGRID
//The values of R and z MUST be equispaced to this code to work.
//#define N_R 100
//#define N_z 100

#define THETA_0 5 //in degrees, starting theta (I avoided 0 because phi is ambiguous) 
#define THETA_STEP 10 //in degrees. Theta loop will iterate as theta = THETA_0 + i*THETA_STEP, as long as theta < 180º
//N_THETA = floor( (180º - THETA_0)/THETA_STEP )

#define PHI_0 0 //in degrees, starting phi
#define PHI_STEP 10 //in degress. Phi loop will iterate as phi = PHI_0 + j*PHI_STEP, as long as phi < 180º
//N_PHI = floor( (180º - PHI_0)/PHI_STEP )

//PRECISION
#define r_number double //Change 'double' by float if you want less accuracy (but less issues with memory)

#define METHOD 2 //0 = Euler, 1 = Midpoint, 2 = Dormand-Prince

#define FIXED_STEP 1 //1 // 0 to adaptive stepsize.
//#define N_FIXED_STEPS 50 //Only used if the stepsize is fixed

#define DP_TOLERANCE 0.05 //Relative tolerance of the Dormand-Prince integrator (Stepsize is chosen to have LESS relative error than this value)
#define APPROX_TOLERANCE 0.01 //Relative error between two numbers in approx function. If |A-B|/|B| < this number, then the code assumes that A is approximately equal to B.

//ERROR HANDLING
#define SOFT_ZERO 1e-8 //Needed for s_to_R method

//CONSTANTS
#define PI 3.14159265359
#define LIGHT_SPEED 1.0//2.99792458E10 //cm/s. Speed of light

#endif
