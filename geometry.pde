/*********************************************************
 * Geometry utility functions.
 *********************************************************/

class vec2 {
  float x, y;
  vec2() {
    x = y = 0.0;
  }
  vec2(float x_, float y_) {
    x = x_;
    y = y_;
  }
  vec2 set(vec2 v) {
    x = v.x;
    y = v.y;
    return this;
  }
  vec2 set(float x_, float y_) {
    x = x_;
    y = y_;
    return this;
  }
}

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



vec2 solveLinearEquationsInTwoVars(float a, float b, float c, float d, float e, float f) {
  float det = a * d - b * c;
  //println("det = ", det);
  if (abs(det) < 0.000001) return null;
  float x =(d * e - b * f) / det;
  float y = (a * f - c * e) / det;
  //System.out.format("x = %f, y = %f\n", x, y);
  return new vec2(x, y);
}


void constructAndSolveLE(vec A, vec B, vec C, vec D, vec2 yz0, vec2 yz1) {
  float dotAC = dot(A, C);
  float dotBD = dot(B, D);
  vec2 t0 = solveLinearEquationsInTwoVars(C.y, C.z, D.y, D.z, dotAC, dotBD);
  vec2 t1 = solveLinearEquationsInTwoVars(C.y, C.z, D.y, D.z, dotAC - C.x, dotBD - D.x);
  yz0.set(t0);
  yz1.set(t1);
  return;
}

void intersectionTwoPlanes(pt A, vec C, pt B, vec D, pt p0, pt p1) {
  // assume that two planes are not parallel
  vec a, b, c, d;
  if (notAbsZero(C.y * D.z - D.y * C.z)) {
    vec2 yz0 = new vec2(), yz1 = new vec2();
    a = new vec(A.x, A.y, A.z);
    b = new vec(B.x, B.y, B.z);
    c = new vec(C.x, C.y, C.z);
    d = new vec(D.x, D.y, D.z);
    constructAndSolveLE(a, b, c, d, yz0, yz1);
    p0.set(0.0, yz0.x, yz0.y);
    p1.set(1.0, yz1.x, yz1.y);
  } else if (notAbsZero(C.x * D.z - D.x * C.z)) {
    vec2 xz0 = new vec2(), xz1 = new vec2();
    a = new vec(A.y, A.x, A.z);
    b = new vec(B.y, B.x, B.z);
    c = new vec(C.y, C.x, C.z);
    d = new vec(D.y, D.x, D.z);
    constructAndSolveLE(a, b, c, d, xz0, xz1);
    p0.set(xz0.x, 0.0, xz0.y);
    p1.set(xz1.x, 1.0, xz1.y);
  } else if (notAbsZero(C.x * D.y - D.x * C.y)) {
    vec2 xy0 = new vec2(), xy1 = new vec2();
    a = new vec(A.z, A.x, A.y);
    b = new vec(B.z, B.x, B.y);
    c = new vec(C.z, C.x, C.y);
    d = new vec(D.z, D.x, D.y);
    constructAndSolveLE(a, b, c, d, xy0, xy1);
    p0.set(xy0.x, xy0.y, 0.0);
    p1.set(xy1.x, xy1.y, 1.0);
  } else {
    println("Can't find intersection of these two planes!");
    exit();
  }
  return;
}


float distanceToLine(pt P, pt A, pt B) {
  vec AP = V(A, P), AB = V(A, B);
  if (parallel(AP, AB)) return 0;
  vec UAB = U(AB);
  float d = n(A(AP, -dot(AP, UAB), UAB));
  //println("distance to line is ", d);
  return d;
}


boolean emptyIntersectionTwoDisks(pt ca, vec va, float ra, pt cb, vec vb, float rb) {
  if (parallel(va, vb)) {  // two planes are parallel
    if (abs(dot(V(ca, cb), va)) < 0.0001) {  // two disks on the same plane
      if (d(ca, cb) > ra + rb) return true;
      else return false;
    } else {
      return true;
    }
  } else {  // two planes are not parallel
    // find intersection line of two supporting planes
    pt p0 = new pt(), p1 = new pt();
    intersectionTwoPlanes(ca, va, cb, vb, p0, p1);
    if (distanceToLine(ca, p0, p1) > ra || distanceToLine(cb, p0, p1) > rb) {
      return true;
    } else {
      return false;
    }
  }
}

boolean emptyIntersectionTwoDisks(Disk d0, Disk d1) {
  return emptyIntersectionTwoDisks(d0.c, d0.d, d0.r, d1.c, d1.d, d1.r);
}