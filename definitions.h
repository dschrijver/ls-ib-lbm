#ifndef DEFINITIONS_H
#define DEFINITIONS_H

// ---------------------------
//     Boundary conditions        
// ---------------------------
#define XPERIODIC
// #define YPERIODIC
#define ZPERIODIC

// --- Half-way Bounce-Back ---
// #define LEFT_BOUNCEBACK_NOSLIP
// #define RIGHT_BOUNCEBACK_NOSLIP
// #define BOTTOM_BOUNCEBACK_NOSLIP
// #define TOP_BOUNCEBACK_NOSLIP
// #define BACK_BOUNCEBACK_NOSLIP
// #define FRONT_BOUNCEBACK_NOSLIP

// --- Non-Equilibrium Bounce-Back ---
// #define LEFT_NEBB_VELOCITY
// #define LEFT_U_VELOCITY 0.0
// #define LEFT_V_VELOCITY 0.0
// #define LEFT_W_VELOCITY 0.0

// #define RIGHT_NEBB_VELOCITY
// #define RIGHT_U_VELOCITY 0.0
// #define RIGHT_V_VELOCITY 0.0
// #define RIGHT_W_VELOCITY 0.0

#define BOTTOM_NEBB_VELOCITY
#define BOTTOM_U_VELOCITY 0.0
#define BOTTOM_V_VELOCITY 0.0
#define BOTTOM_W_VELOCITY 0.0

#define TOP_NEBB_VELOCITY
#define TOP_U_VELOCITY 0.0
#define TOP_V_VELOCITY 0.0
#define TOP_W_VELOCITY 0.0

// #define BACK_NEBB_VELOCITY
// #define BACK_U_VELOCITY 0.0
// #define BACK_V_VELOCITY 0.0
// #define BACK_W_VELOCITY 0.0

// #define FRONT_NEBB_VELOCITY
// #define FRONT_U_VELOCITY 0.0
// #define FRONT_V_VELOCITY 0.0
// #define FRONT_W_VELOCITY 0.0


// --------------------------
//     Initial conditions    
// --------------------------
#define INI_POISEUILLE
#define WIDTH 3

// --------------------
//     Logic checks        
// --------------------
#if defined(LEFT_NEBB_VELOCITY) || defined(RIGHT_NEBB_VELOCITY) || defined(BOTTOM_NEBB_VELOCITY) || defined(TOP_NEBB_VELOCITY) || defined(BACK_NEBB_VELOCITY) || defined(FRONT_NEBB_VELOCITY)
#define WETNODE
#endif

#endif