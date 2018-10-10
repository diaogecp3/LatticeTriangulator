/******************************************************************************
 * Geometric utility functions.
 *
 * Intersection tests/computation between geometric primitives.
 * Construction of basic geometric primitives.
 ******************************************************************************/


vec normalToTriangle(pt A, pt B, pt C) {
  return U(N(A, B, C));
}

void showNormalToTriangle(pt A, pt B, pt C, float d, float r) {
  vec N = normalToTriangle(A, B, C);
  pt D = P(A, B, C);
  arrow(D, V(d, N), r);
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