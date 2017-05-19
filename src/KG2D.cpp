#include "KG2D.hpp"
/********************* Indexing **************************/
size_t inline IDX(size_t x, size_t y, size_t z) {
  return x + NX * (y + NY * z);
}
/******************** Averaging **************************/
FType mean(FType *A, size_t nA) {
  FType mu = 0.0;
  for (size_t n = 0; n < nA; n++) mu += A[n];
  return mu / ((FType)nA);
}  // Average
/**************** Standard Deviation *********************/
FType stdev(FType *A, size_t nA) {
  FType mu = mean(A, nA);
  FType sig = 0.0;
  for (size_t n = 0; n < nA; n++) sig += pow(A[n] - mu, 2.0);
  return sqrt(sig / nA);
}
/********************************************************/
KG2D::KG2D(size_t RECL_, FType dt_, FType tmx_, bool log_, FType E0_,
           size_t LX = NX, size_t LY = NY)
    : SI_ALG(dt_), TICTOC(RECL_, dt_, tmx_, log_), E0(E0_) {
  // Allocate
  size_t B = sizeof(FType), MB = B * (5 * N + NZ);
  printf("Est. Req. RAM: %g MB\n", 1.3 * MB / pow(1024, 2.0));
  u = (FType *)calloc(N, B);
  p = (FType *)calloc(N, B);
  eps = (FType *)calloc(N, B);
  zeta = (FType *)calloc(N, B);
  stencil = (FType *)calloc(N, B);
  E = (FType *)calloc(NZ, B);
  // Populate
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine genny(seed);
  std::uniform_real_distribution<FType> distro(0.5, 1.5);
  k_loop eps[k] = distro(genny);
  // Get packet dims.
  size_t LXh = LX / 2 + LX % 2, LYh = LY / 2 + LY % 2;
  size_t NXh = NX / 2 + NX % 2, NYh = NY / 2 + NY % 2;
  // Assign packet.
  FType dE = sqrt(2.0 * E0 / (FType)(4.0 * LXh * LYh));
  for (size_t z = 0; z < NZ; z++)
    for (size_t y = NYh - LYh; y < NYh + LYh; y++)
      for (size_t x = NXh - LXh; x < NXh + LXh; x++) p[IDX(x, y, z)] = dE;
}
/********************************************************/
KG2D::~KG2D(void) {
  free(u);
  free(p);
  free(eps);
  free(E);
}
/********************************************************/
void KG2D::INFO(void) {
  printf("•INTEGRATION:\t");
  printf("dt=%g  tmax=%g  RECL=%lu\n", dt, tmax, RECL);
  printf("•LATTICE:\t");
  printf("(NX,NY,NZ)=(%d,%d,%d)\t", NX, NY, NZ);
  printf("W=%g  σ=%g\n", W, sigma);
  printf("\n");
}
/********************************************************/
void KG2D::STENCIL_JACOBI5(FType *A) {
  xyz_loop {
    stencil[IDX(x, y, z)] = -4.0 * A[IDX(x, y, z)];
    // Left Edge
    if (x > 0)
      stencil[IDX(x, y, z)] += A[IDX(x - 1, y, z)];
    else if (PBC)
      stencil[IDX(x, y, z)] += A[IDX(NX - 1, y, z)];
    // Top Edge
    if (y > 0)
      stencil[IDX(x, y, z)] += A[IDX(x, y - 1, z)];
    else if (PBC)
      stencil[IDX(x, y, z)] += A[IDX(x, NY - 1, z)];
    // Right Edge
    if (x < NX - 1)
      stencil[IDX(x, y, z)] += A[IDX(x + 1, y, z)];
    else if (PBC)
      stencil[IDX(x, y, z)] += A[IDX(0, y, z)];
    // Bottom Edge
    if (y < NY - 1)
      stencil[IDX(x, y, z)] += A[IDX(x, y + 1, z)];
    else if (PBC)
      stencil[IDX(x, y, z)] += A[IDX(x, 0, z)];
  }
}
/********************************************************/
// Symplectic Steps
void KG2D::STEP_A(FType tau) { k_loop u[k] += tau * p[k]; }
/********************************************************/
void KG2D::STEP_B(FType tau) {
  k_loop zeta[k] = (eps[k] + pow(abs(u[k]), sigma)) * u[k];
  STENCIL_JACOBI5(u);
  k_loop p[k] -= tau * (zeta[k] - (stencil[k] / W));
}
/********************************************************/
void KG2D::STEP_C(FType tau) {
  STENCIL_JACOBI5(u);
  k_loop zeta[k] = (eps[k] + pow(abs(u[k]), sigma)) * u[k] - (stencil[k] / W);
  STENCIL_JACOBI5(zeta);
  k_loop zeta[k] *= (eps[k] + (sigma + 1.0) * pow(abs(u[k]), sigma));
  k_loop p[k] += 2.0 * tau * (zeta[k] - stencil[k] / W);
}
/********************************************************/
void KG2D::NRG_DENSITY(FType *dE) {
  FType s2 = sigma + 2.0, Wd = 4.0 * W;
  k_loop dE[k] =
      (pow(p[k], 2) + eps[k] * pow(u[k], 2.0)) / 2.0 + pow(abs(u[k]), s2) / s2;
  xyz_loop {
    if (x > 0)
      dE[IDX(x, y, z)] += pow(u[IDX(x, y, z)] - u[IDX(x - 1, y, z)], 2.0) / Wd;
    else if (PBC)
      dE[IDX(x, y, z)] += pow(u[IDX(x, y, z)] - u[IDX(NX - 1, y, z)], 2.0) / Wd;
    else
      dE[IDX(x, y, z)] += pow(u[IDX(x, y, z)], 2.0) / Wd;
    //
    if (y > 0)
      dE[IDX(x, y, z)] += pow(u[IDX(x, y, z)] - u[IDX(x, y - 1, z)], 2.0) / Wd;
    else if (PBC)
      dE[IDX(x, y, z)] += pow(u[IDX(x, y, z)] - u[IDX(x, NY - 1, z)], 2.0) / Wd;
    else
      dE[IDX(x, y, z)] += pow(u[IDX(x, y, z)], 2.0) / Wd;
    //
    if (x < NX - 1)
      dE[IDX(x, y, z)] += pow(u[IDX(x + 1, y, z)] - u[IDX(x, y, z)], 2.0) / Wd;
    else if (PBC)
      dE[IDX(x, y, z)] += pow(u[IDX(0, y, z)] - u[IDX(x, y, z)], 2.0) / Wd;
    else
      dE[IDX(x, y, z)] += pow(-u[IDX(x, y, z)], 2.0) / Wd;
    //
    if (y < NY - 1)
      dE[IDX(x, y, z)] += pow(u[IDX(x, y + 1, z)] - u[IDX(x, y, z)], 2.0) / Wd;
    else if (PBC)
      dE[IDX(x, y, z)] += pow(u[IDX(x, 0, z)] - u[IDX(x, y, z)], 2.0) / Wd;
    else
      dE[IDX(x, y, z)] += pow(-u[IDX(x, y, z)], 2.0) / Wd;
  }
}
/********************************************************/
void KG2D::NRG(void) {
  FType *dE = (FType *)calloc(N, sizeof(FType));
  NRG_DENSITY(dE);
  for (size_t z = 0; z < NZ; z++) E[z] = 0.0;
  xyz_loop E[z] += dE[IDX(x, y, z)];
}
/********************************************************/
void KG2D::PRINT_NRG(void) {
  NRG();
  printf("E0=%.8e\n", E0);
  for (size_t z = 0; z < NZ; z++) printf("E[%lu]=%.8e\n", z, E[z]);
}
/********************************************************/
void KG2D::PRINT_NRG_ERR(void) {
  NRG();
  FType dE[NZ];
  for (size_t z = 0; z < NZ; z++) dE[z] = log10(abs((E[z] - E0) / E0) + FPE);
  FType mu = mean(dE, NZ);
  FType sig = stdev(dE, NZ);
  printf("%.4e\t%.4e\t%.4e\t%.4e\n", t, TOC(false), mu, sig);
}
/********************************************************/
void KG2D::Dump2File(std::string fid, bool wipe) {
  // std::ostringstream ss;
  // ss << REC;
  // std::string SREC = ss.str();
  // std::string kU = "U" + SREC, kP = "P" + SREC;
  // if (wipe) {
  //   remove(fid.c_str());
  //   saveArray(kU.c_str(), u, fid.c_str(), false);
  // } else {
  //   saveArray(kU.c_str(), u, fid.c_str(), true);
  // }
  // saveArray(kP.c_str(), p, fid.c_str(), true);
}
/********************************************************/
void KG2D::TextDump(FType *A, std::string fid, size_t z) {
  FILE *pFile = fopen(fid.c_str(), "w");
  for (size_t x = 0; x < NX; x++) {
    for (size_t y = 0; y < NY; y++) fprintf(pFile, "%.8e\t", A[IDX(x, y, z)]);
    fprintf(pFile, "\n");
  }
  fclose(pFile);
}
/********************************************************/