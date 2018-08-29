float average(float[] a, int n, int start, int end) {
  assert start >= 0 && end <= n && end > start;
  float sum = 0.0;
  for (int i = start; i < end; ++i) {
    sum += a[i];
  }
  return sum / (end - start);
}

void exceptionHandler() {
  P.savePts("data/pts");
}