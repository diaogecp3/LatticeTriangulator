/*********************************************************
 * Utility functions.
 *********************************************************/


/*
 * Compute average value for range [start, end) of array a.
 */
float average(float[] a, int n, int start, int end) {
  assert start >= 0 && end <= n && end > start;
  float sum = 0.0;
  for (int i = start; i < end; ++i) {
    sum += a[i];
  }
  return sum / (end - start);
}


/*
 * Exception handler.
 */
void exceptionHandler() {
  P.savePts("data/pts");
}

boolean notAbsZero(float x) {
  return x <= -0.00001 || x >= 0.00001;
}

boolean isAbsZero(float x) {
  return x > -0.00001 && x < 0.00001;
}

boolean notZero(float x) {
  // assume that x >= 0
  return x >= 0.00001;
}

boolean isZero(float x) {
  // assume that x >= 0
  return x < 0.00001;
}