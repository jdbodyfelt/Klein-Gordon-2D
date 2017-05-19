#include "symplktk.hpp"
/********************************************************/
TICTOC::TICTOC(size_t RECL_, FType dt_, FType tmx_, bool log_)
    : t(0.0), RECL(RECL_), dt(dt_), tmax(tmx_), logf(log_) {
  TIC = hi_res_clock::now();
  REC = logf ? 0 : 1;
  ss_dt = logf ? exp(log(tmax / dt) / (FType)(RECL - 1))
               : tmax / (dt * (FType)(RECL - 1));
  ss_t = logf ? pow(ss_dt, (FType)REC) * dt : ss_dt * dt * (FType)REC;
  printf("🢡 TICTOC Initialized.\n");
}
/********************************************************/
TICTOC::~TICTOC(void) {}
/********************************************************/
FType TICTOC::TOC(bool reset = false) {
  using namespace std::chrono;
  typedef duration<FType> dFType;
  if (reset) TICTOC(RECL, dt, tmax, logf);
  hi_res_clock::time_point NOW = hi_res_clock::now();
  dFType span = duration_cast<dFType>(NOW - TIC);
  return span.count();
}
/********************************************************/
bool TICTOC::Update(void) {
  t += dt;
  //printf("t=%.8e\n", t);
  if (t < ss_t) {
    return false;
  } else {
    ++REC;
    ss_t = logf ? pow(ss_dt, REC) * dt : ss_dt * dt * REC;
    return true;
  }
}
/********************************************************/
SBAB2C::SBAB2C(FType dt_) {
  for (size_t k = 0; k < Nc; k++) C[k] *= dt_;
  G *= pow(Geps, 2) * pow(dt_, 3);
  printf("🢡 SBAB2C Initialized.\n");
}
/********************************************************/
SBAB2C::~SBAB2C(void){};  // printf("🢡 SBAB2C Destroyed.\n"); }
/********************************************************/
void SBAB2C::STEP_FWD(void) { STEP_C(0.5 * G); }
/********************************************************/
void SBAB2C::STEP_BCK(void) { STEP_C(-0.5 * G); }
/********************************************************/
void SBAB2C::STEP_FULL(void) {
  STEP_B(C[0]);
  STEP_A(C[1]);
  STEP_B(C[2]);
  STEP_A(C[1]);
  STEP_B(C[0]);
  STEP_C(G);
}
/********************************************************/