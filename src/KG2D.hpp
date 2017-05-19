/*******************************************************
 * Copyright (c) 2017, J.D. Bodyfelt
 * All rights reserved.
 *
 * This file is distributed under 3-clause BSD license.
 * The complete license agreement can be obtained at:
 * <add URL>
 ********************************************************/
#ifndef __KG2D_HPP__
#define __KG2D_HPP__
/*************** Library Imports ************************/
#include "symplktk.hpp"
#include <random>

/*************** Macro Definitions **********************/
// NOTE: Reset w compile flag "-D<VAR>=<VAL>"
#ifndef NX
#define NX 256 // Site Count on x-dimension
#endif
#ifndef NY
#define NY 256 // Site Count on y-dimension
#endif
#ifndef NZ
#define NZ 128 // Statistical Population
#endif

#ifndef PBC
#define PBC true // Periodic Boundaries
#endif

#define k_loop for(size_t k=0; k<N; k++)
#define xyz_loop                                                               \
  for (size_t z = 0; z < NZ; z++)                                              \
    for (size_t y = 0; y < NY; y++)                                            \
      for (size_t x = 0; x < NX; x++)

/******************** Statistics *************************/
size_t inline IDX(size_t x, size_t y, size_t z);
FType mean(FType *A, size_t nA);
FType stdev(FType *A, size_t nA);
/*************** Class Definitions **********************/
class KG2D : public SI_ALG, public TICTOC {
public:
  size_t N = NX * NY * NZ;
  FType *u, *p, *eps;
  FType *zeta, *stencil, *E;
  FType W = 4.0, sigma = 2.0, E0;

  KG2D(size_t RECL_, FType dt_, FType tmx_, bool log_, FType E0_, size_t LX,
       size_t LY);
  ~KG2D(void);
  void INFO(void);
  // Virtual Inherits from SI_ALG
  void STENCIL_JACOBI5(FType *A);
  void STEP_A(FType tau); // U Update
  void STEP_B(FType tau); // P Update
  void STEP_C(FType tau); // Corrector Update
  void NRG_DENSITY(FType *dE);
  void NRG(void);
  void PRINT_NRG(void);
  void PRINT_NRG_ERR(void);
  void Dump2File(std::string fid, bool wipe);
  void TextDump(FType *A, std::string fid, size_t idx);
};
/*******************************************************/

#endif