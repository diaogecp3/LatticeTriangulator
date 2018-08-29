// geometry util functions

vec normalToTriangle(pt A, pt B, pt C) {
  vec N = N(A, B, C);
  N = U(N);
  return N;
}

void showNormalToTriangle(pt A, pt B, pt C, float d, float r) {
  vec N = normalToTriangle(A, B, C);
  pt D = P(A, B, C);
  arrow(D, V(d, N), r);
}