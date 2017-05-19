#include "KG2D.hpp"

int main(void) {
  size_t L = 16;

  KG2D X(128, 5e-2, 1e2, false, L * L / 2.0, L, L);
  X.W = 4.0;
  X.sigma = 2.0;

  X.INFO();
  X.PRINT_NRG_ERR();
  X.STEP_FWD();

  while (X.REC < X.RECL) {
    X.STEP_FULL();
    if (X.Update()) {
      X.STEP_BCK();
      X.PRINT_NRG_ERR();
      X.STEP_FWD();
    }
  }
  return 0;
}