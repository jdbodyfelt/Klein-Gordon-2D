#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>

typedef std::chrono::high_resolution_clock hi_res_clock;
/*************** Macro Definitions **********************/
#ifndef DP
#define DP true  // Default double precision.
#endif
#if DP
#define FPE DBL_MIN  // Need to maybe template this...
typedef double FType;
#define FP f64
#else
#define FPE FLT_MIN
typedef float FType;
#define FP f32
#endif

#ifndef SI_ALG
#define SI_ALG SBAB2C
#endif
/*************** Class Definitions **********************/
class TICTOC {
 protected:
  hi_res_clock::time_point TIC;
  size_t ss_k;

 public:
  bool logf;
  size_t REC, RECL;
  FType t, dt, tmax, ss_dt, ss_t;

  TICTOC(size_t RECL_, FType dt_, FType tmx_, bool log_);
  ~TICTOC(void);
  FType TOC(bool);
  bool Update(void);
};
/********************************************************/
class SBAB2C {
 protected:
  const static size_t Nc = 3;
  FType C[Nc] = {1.0 / 6.0, 0.5, 2.0 / 3.0};
  FType G = 1.0 / 72.0;
  FType Geps = 1.0;

 public:
  SBAB2C(FType dt_);
  ~SBAB2C(void);
  virtual void STEP_A(FType tau) = 0;
  virtual void STEP_B(FType tau) = 0;
  virtual void STEP_C(FType tau) = 0;
  void STEP_FWD(void);
  void STEP_BCK(void);
  void STEP_FULL(void);
};
/********************************************************/
