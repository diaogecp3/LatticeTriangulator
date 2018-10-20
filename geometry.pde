/******************************************************************************
 * Geometric utility functions.
 *
 * Intersection tests between geometric primitives.
 * Construction of basic geometric primitives.
 * And so on.
 ******************************************************************************/


vec normalToTriangle(pt A, pt B, pt C) {
  return U(N(A, B, C));
}

void showNormalToTriangle(pt A, pt B, pt C, float d, float r) {
  vec N = normalToTriangle(A, B, C);
  pt D = P(A, B, C);
  arrow(D, V(d, N), r);
}

void showPlane(pt p, vec n, float s) {
  vec vi = constructNormal(n);
  vec vj = N(n, vi);
  pt p0 = P(p, s, vi, s, vj);
  pt p1 = P(p, -s, vi, s, vj);
  pt p2 = P(p, s, vi, -s, vj);
  pt p3 = P(p, -s, vi, -s, vj);
  beginShape(TRIANGLE_STRIP);
  vertex(p0);
  vertex(p1);
  vertex(p2);
  vertex(p3);
  endShape();
}


/*
 * Construct a vector normal to given vector v. v is not necessarily a unit
 * vector. A unit vector will be returned.
 */
vec constructNormal(vec v) {
  if (isAbsZero(v.y) && isAbsZero(v.z)) {  // v is parallel to (1, 0, 0)
    return U(new vec(-v.z, 0.0, v.x));  // cross product of v and (0, 1, 0)
  } else {
    return U(new vec(0.0, v.z, -v.y));  // cross product of v and (1, 0, 0)
  }
}



/*
 * This is a helper function to find the two intersection points of two planes.
 */
void constructAndSolveLE(vec A,                                      // in
                         vec B,                                      // in
                         vec C,                                      // in
                         vec D,                                      // in
                         vec2 yz0,                                   // out
                         vec2 yz1) {                                 // out
  float dotAC = dot(A, C);
  float dotBD = dot(B, D);
  vec2 t0 = solveLinearEquationsInTwoVars(C.y, C.z, D.y, D.z, dotAC, dotBD);
  vec2 t1 = solveLinearEquationsInTwoVars(C.y, C.z, D.y, D.z, dotAC - C.x, dotBD - D.x);
  yz0.set(t0);
  yz1.set(t1);
  return;
}


/*
 * Compute the intersection line of two planes, given that they are not parallel.
 */
void intersectionTwoPlanes(pt A,                                     // in
                           vec C,                                    // in
                           pt B,                                     // in
                           vec D,                                    // in
                           pt p0,                                    // out
                           pt p1) {                                  // out
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

pt intersectionTwoLines(pt pa, pt pb, pt pc, pt pd) {
  vec v0 = U(pa, pb);  // unit vector
  vec v1 = U(pc, pd);  // unit vector
  vec vac = V(pa, pc);

  vec c0 = N(v0, v1);
  vec c1 = N(vac, v1);
  float x = c1.norm() / c0.norm();

  /* In case when cross(v0, v1) and cross(vac, v1) don't point in the same direction. */
  if (c0.x * c1.x < 0 || c0.y * c1.y < 0 || c0.z * c1.z < 0) x = -x;

  return P(pa, x, v0);
}


boolean emptyIntersectionLineDisk(pt a, pt b, pt c, float r, vec n) {
  vec ab = V(a, b);
  vec ac = V(a, c);
  float dnab = dot(n, ab);
  float dnac = dot(n, ac);
  if (notAbsZero(dnab)) {  // dot(N, AB) != 0
    float t = dnac / dnab;
    pt p = P(a, t, ab);
    if (d(p, c) <= r) return false;
    else return true;
  } else {  // dot(N, AB) == 0
    if (notAbsZero(dnac)) {  // dot(N, AC) != 0
      return true;  // line AB parallel to the supporting plane of circle C
    } else {  // line AB is on the same plane as circle C
      if (distanceToLine(c, a, b) <= r) return false;
      else return true;
    }
  }
}


// solve a*cos(theta) + b*sin(theta) = c
float[] solveLinearEquationInCosSin(float a, float b, float c) {
  assert a * a + b * b > 0.0000001;  // at least one of {a, b} is non-zero
  float[] thetas = new float[2];
  if (notAbsZero(a) && notAbsZero(b)) {  // a != 0 and b != 0
    //println("a != 0 and b != 0");
    println("a = ", a, "b = ", b);
    /* How to decide tmp0 here? According to the sign of b. */
    float tmp0 = asin(c / (sqrt(a * a + b * b)));  // [-pi/2, pi/2]
    //float tmp0 = asin(-c / (sqrt(a * a + b * b)));  // [-pi/2, pi/2]
    if (b < 0) tmp0 = -tmp0;  // [-pi/2, pi/2]
    float tmp1 = atan(a / b);  // (-pi/2, pi/2)
    thetas[0] = tmp0 - tmp1;  // (-pi, pi)
    thetas[1] = PI - tmp0 - tmp1;  // (0, 2pi)
    thetas[0] = thetas[0] < 0 ? thetas[0] + TWO_PI : thetas[0];  // [0, pi), (pi, 2pi)
  } else if (isAbsZero(a)) {  // a == 0 and b != 0
    println("a == 0");
    thetas[0] = asin(c/b);  // [-pi/2, pi/2]
    thetas[1] = PI - thetas[0];  // [pi/2, 3pi/2]
    if (thetas[0] < 0) thetas[0] += TWO_PI;  // [3pi/2, 0), [0, pi/2]
  } else {  // a != 0 and b == 0
    println("b == 0");
    thetas[0] = acos(c/a);  // [0, pi]
    thetas[1] = TWO_PI - thetas[0];  // [pi, 2pi]
  }
  return thetas;
}

pt[] pivotPlaneAroundLineHitCircle(pt c, float r, vec n, pt a, pt b, vec vi, vec vj) {
  boolean empty = emptyIntersectionLineDisk(a, b, c, r, n);
  println("line AB intersects circle C?", empty ? "no" : "yes");

  vec ca = V(c, a);
  vec ab = V(a, b);
  if (vi == null) vi = constructNormal(n);
  float dnab = dot(n, ab);
  float dnca = dot(n, ca);
  pt[] ps = new pt[2];
  if (isAbsZero(dnab) && isAbsZero(dnca)) {
    println("Line AB is on the supporting plane of Circle C");
    ps[0] = P(c, r, vi);
    ps[1] = P(c, -r, vi);
    return ps;
  }
  if (vj == null) vj = N(n, vi);  // cross product
  float[] thetas;  // valid theta should be in [-pi, pi]?
  if (notAbsZero(dnab)) {
    // println("dot(N, AB) != 0");
    vec v = A(ca, -dnca/dnab, ab);
    float fa = dot(vi, v);
    float fb = dot(vj, v);
    thetas = solveLinearEquationInCosSin(fa, fb, r);
  } else {
    float fa = r * dot(vi, ab);
    float fb = r * dot(vj, ab);
    thetas = solveLinearEquationInCosSin(fa, fb, 0);
  }

  assert thetas.length == 2;
  // println("theta 1 =", thetas[0], "theta 2 =", thetas[1]);
  ps[0] = P(c, r * cos(thetas[0]), vi, r * sin(thetas[0]), vj);
  ps[1] = P(c, r * cos(thetas[1]), vi, r * sin(thetas[1]), vj);
  return ps;
}

pt[] tangentPlaneThreeCircles(pt c0, float r0, vec n0, vec vi0, vec vj0,
                              pt c1, float r1, vec n1, vec vi1, vec vj1,
                              pt c2, float r2, vec n2, vec vi2, vec vj2) {
  if (vi0 == null) vi0 = constructNormal(n0);
  if (vj0 == null) vj0 = N(n0, vi0);
  if (vi1 == null) vi1 = constructNormal(n1);
  if (vj1 == null) vj1 = N(n1, vi1);
  if (vi2 == null) vi2 = constructNormal(n2);
  if (vj2 == null) vj2 = N(n2, vi2);

  /* Initialize a triangle with each vertex from each circle. */
  pt pa = P(c0, r0, vi0);
  pt pb = P(c1, r1, vi1);
  pt pc = P(c2, r2, vi2);
  vec ta = vj0;
  vec tb = vj1;
  vec tc = vj2;
  vec n = normalToTriangle(pa, pb, pc);

  int maxIter = 5;
  int iter = 0;
  while (iter < maxIter) {
    boolean update = false;
    // move A? pivot around BC
    if (notAbsZero(dot(n, ta)) || dot(n, V(pb, c0)) > 0) {
      pt[] candidates = pivotPlaneAroundLineHitCircle(c0, r0, n0, pb, pc, vi0, vj0);
      pa = candidates[0];
      n = N(pb, pc, pa);
      if (dot(n, V(pb, c0)) > 0) {
        pa = candidates[1];
        n = N(pb, pc, pa);
      }
      n.normalize();
      ta = N(n0, U(c0, pa));
      update = true;
    }

    // move B? pivot around CA
    if (notAbsZero(dot(n, tb)) || dot(n, V(pc, c1)) > 0) {
      pt[] candidates = pivotPlaneAroundLineHitCircle(c1, r1, n1, pc, pa, vi1, vj1);
      pb = candidates[0];
      n = N(pc, pa, pb);
      if (dot(n, V(pc, c1)) > 0) {
        pb = candidates[1];
        n = N(pc, pa, pb);
      }
      n.normalize();
      tb = N(n1, U(c1, pb));
      update = true;
    }

    // move C? pivot around AB
    if (notAbsZero(dot(n, tc)) || dot(n, V(pa, c2)) > 0) {
      pt[] candidates = pivotPlaneAroundLineHitCircle(c2, r2, n2, pa, pb, vi2, vj2);
      pc = candidates[0];
      n = N(pa, pb, pc);
      if (dot(n, V(pa, c2)) > 0) {
        pc = candidates[1];
        n = N(pa, pb, pc);
      }
      n.normalize();
      tc = N(n2, U(c2, pc));
      update = true;
    }

    if (!update) break;
    iter++;
  }  // end while
  println("iter =", iter);
  pt[] ps = new pt[3];
  ps[0] = pa;
  ps[1] = pb;
  ps[2] = pc;
  return ps;
}