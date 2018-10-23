/******************************************************************************
 * Geometric utility functions.
 *
 * Intersection tests between geometric primitives.
 * Construction of basic geometric primitives.
 * And so on.
 ******************************************************************************/


boolean show3T = true;
boolean show2T = true;
int methodEP = 0;  // 0: basic, 1: heuristic normal, 2: heuristic plane


/*
 * Construct the normal to triangle (A, B, C).
 */
vec normalToTriangle(pt A, pt B, pt C) {
  return U(N(A, B, C));
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
boolean intersectionTwoPlanes(pt A,                                     // in
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
    return false;
  }
  return true;
}

/*
 * Compute the distance from point P to line AB.
 */
float distanceToLine(pt P, pt A, pt B) {
  vec AP = V(A, P), AB = V(A, B);
  if (parallel(AP, AB)) return 0;
  vec UAB = U(AB);
  float d = n(A(AP, -dot(AP, UAB), UAB));
  return d;
}

/*
 * Return true if disk (ca, va, ra) and disk (cb, vb, rb) don't intersect. The
 * algorithm is: check if the two supporting planes intersect, and compute the
 * distances from the two centers to the intersection line, respectively, and
 * check if the two distances is smaller than or equal to ra/rb.
 */
boolean emptyIntersectionTwoDisks(pt ca, vec va, float ra, pt cb, vec vb, float rb) {
  if (parallel(va, vb)) {  // two planes are parallel
    if (isAbsZero(dot(V(ca, cb), va))) {  // two disks on the same plane
      if (d(ca, cb) > ra + rb) return true;
      else return false;
    } else {
      return true;
    }
  } else {  // two planes are not parallel
    // find intersection line of two supporting planes
    pt p0 = new pt(), p1 = new pt();
    if (!intersectionTwoPlanes(ca, va, cb, vb, p0, p1)) return true;
    if (distanceToLine(ca, p0, p1) > ra || distanceToLine(cb, p0, p1) > rb) {
      return true;
    } else {
      return false;
    }
  }
}

/*
 * Return true if disk d0 and disk d1 don't intersect.
 */
boolean emptyIntersectionTwoDisks(Disk d0, Disk d1) {
  return emptyIntersectionTwoDisks(d0.c, d0.d, d0.r, d1.c, d1.d, d1.r);
}

/*
 * Compute the intesection point of line (pa, pb) and line (pc, pd). Assume that
 * these two lines are on the same plane.
 */
pt intersectionTwoLines(pt pa, pt pb, pt pc, pt pd) {
  vec v0 = U(pa, pb);  // unit vector
  vec v1 = U(pc, pd);  // unit vector
  vec vac = V(pa, pc);

  vec c0 = N(v0, v1);
  vec c1 = N(vac, v1);
  // float x = c1.norm() / c0.norm();
  float x = sqrt(dot(c1, c1) / dot(c0, c0));

  /* In case when cross(v0, v1) and cross(vac, v1) don't point in the same direction. */
  if (c0.x * c1.x < 0 || c0.y * c1.y < 0 || c0.z * c1.z < 0) x = -x;

  return P(pa, x, v0);
}

/*
 * Return true if line (a, b) and disk (c, r, n) don't intersect.
 */
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

/*
 * Return true if circle (c, r, n) intersects plane (p, d). If true, (p0, p1)
 * will be the intersection line. p0 and p1 will be on the circle.
 */
boolean intersectionCirclePlane(pt c, float r, vec n, pt p, vec d, pt p0, pt p1) {
  pt pa = new pt();
  pt pb = new pt();
  if (!intersectionTwoPlanes(c, n, p, d, pa, pb)) return false;
  vec vab = V(pa, pb);
  vec vca = V(c, pa);
  float a = dot(vab, vab);
  float b = 2 * dot(vca, vab);
  float e = dot(vca, vca) - r * r;
  float[] xs = solveQuadraticEquation(a, b, e);
  if (xs == null) return false;
  if (p0 != null) p0.set(P(pa, xs[0], vab));
  if (p1 != null) p1.set(P(pa, xs[1], vab));
  return true;
}

/*
 * Return the points of contact between the plane that passes through line (a, b)
 * and the circle (c, r, n). (vi, vj), which is optional, is the local frame on
 * the circle.
 */
pt[] pivotPlaneAroundLineHitCircle(pt c, float r, vec n, pt a, pt b, vec vi, vec vj) {
  pt[] ps = new pt[2];
  boolean empty = emptyIntersectionLineDisk(a, b, c, r, n);
  if (!empty) {
    println("Line AB intersects disk C");
    ps[0] = P(c, r, vi);
    ps[1] = P(c, -r, vi);
    return ps;
  }

  vec ca = V(c, a);
  vec ab = V(a, b);
  if (vi == null) vi = constructNormal(n);
  float dnab = dot(n, ab);
  float dnca = dot(n, ca);
  if (isAbsZero(dnab) && isAbsZero(dnca)) {
    println("Line AB is on the supporting plane of Circle C");
    ps[0] = P(c, r, vi);
    ps[1] = P(c, -r, vi);
    return ps;
  }
  if (vj == null) vj = N(n, vi);  // cross product
  float[] thetas;
  if (notAbsZero(dnab)) {
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
  ps[0] = P(c, r * cos(thetas[0]), vi, r * sin(thetas[0]), vj);
  ps[1] = P(c, r * cos(thetas[1]), vi, r * sin(thetas[1]), vj);
  return ps;
}

/*
 * Return the points of contact that define the tangent plane of the three
 * circles (c0, r0, n0), (c1, r1, n1), and (c2, r2, n2). (vi0, vj0), (vi1, vj1),
 * (vi2, vj2), which are optional, are local frames on the three circles,
 * respectively.
 */
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
  pt p0, p1, p2;
  vec t0, t1, t2, n;
  // println("methodEP = ", methodEP);
  switch (methodEP) {
    case 1:
      vec tmpN1 = normalToTriangle(c0, c1, c2);
      p0 = P(c0, r0, U(A(tmpN1, -dot(tmpN1, n0), n0)));
      p1 = P(c1, r1, U(A(tmpN1, -dot(tmpN1, n1), n1)));
      p2 = P(c2, r2, U(A(tmpN1, -dot(tmpN1, n2), n2)));
      t0 = U(N(n0, V(c0, p0)));
      t1 = U(N(n1, V(c1, p1)));
      t2 = U(N(n2, V(c2, p2)));
      n = normalToTriangle(p0, p1, p2);
      break;
    case 2:
      vec tmpN2 = normalToTriangle(c0, c1, c2);
      p0 = new pt();
      intersectionCirclePlane(c0, r0, n0, c0, tmpN2, p0, null);
      p1 = new pt();
      intersectionCirclePlane(c1, r1, n1, c0, tmpN2, p1, null);
      p2 = new pt();
      intersectionCirclePlane(c2, r2, n2, c0, tmpN2, p2, null);
      t0 = U(N(n0, V(c0, p0)));
      t1 = U(N(n1, V(c1, p1)));
      t2 = U(N(n2, V(c2, p2)));
      n = normalToTriangle(p0, p1, p2);
      break;
    default:
      p0 = P(c0, r0, vi0);
      p1 = P(c1, r1, vi1);
      p2 = P(c2, r2, vi2);
      t0 = vj0;
      t1 = vj1;
      t2 = vj2;
      n = normalToTriangle(p0, p1, p2);
  }

  int maxIter = 10;
  int iter = 0;
  while (iter < maxIter) {
    boolean update = false;
    // move A, pivot around BC
    if (notAbsZero(dot(n, t0)) || dot(n, V(p1, c0)) > 0) {
      pt[] candidates = pivotPlaneAroundLineHitCircle(c0, r0, n0, p1, p2, vi0, vj0);
      p0 = candidates[0];
      n = N(p1, p2, p0);
      if (dot(n, V(p1, c0)) > 0) {
        p0 = candidates[1];
        n = N(p1, p2, p0);
      }
      n.normalize();
      t0 = N(n0, U(c0, p0));
      update = true;
    }

    // move B, pivot around CA
    if (notAbsZero(dot(n, t1)) || dot(n, V(p2, c1)) > 0) {
      pt[] candidates = pivotPlaneAroundLineHitCircle(c1, r1, n1, p2, p0, vi1, vj1);
      p1 = candidates[0];
      n = N(p2, p0, p1);
      if (dot(n, V(p2, c1)) > 0) {
        p1 = candidates[1];
        n = N(p2, p0, p1);
      }
      n.normalize();
      t1 = N(n1, U(c1, p1));
      update = true;
    }

    // move C, pivot around AB
    if (notAbsZero(dot(n, t2)) || dot(n, V(p0, c2)) > 0) {
      pt[] candidates = pivotPlaneAroundLineHitCircle(c2, r2, n2, p0, p1, vi2, vj2);
      p2 = candidates[0];
      n = N(p0, p1, p2);
      if (dot(n, V(p0, c2)) > 0) {
        p2 = candidates[1];
        n = N(p0, p1, p2);
      }
      n.normalize();
      t2 = N(n2, U(c2, p2));
      update = true;
    }

    if (!update) break;
    iter++;
  }  // end while
  // println("iter =", iter);
  pt[] ps = new pt[3];
  ps[0] = p0;
  ps[1] = p1;
  ps[2] = p2;
  return ps;
}

/*
 * Construct and show the exact convex hull defined by three circles (c0, r0, n0),
 * (c1, r1, n2), and (c2, r2, n2). (vi0, vj0), (vi1, vj1), (vi2, vj2), which are
 * optional, are local frames on the three circles, respectively.
 */
void exactCHThreeCircles(pt c0, float r0, vec n0, vec vi0, vec vj0,
                         pt c1, float r1, vec n1, vec vi1, vec vj1,
                         pt c2, float r2, vec n2, vec vi2, vec vj2) {
  if (vi0 == null) vi0 = constructNormal(n0);
  if (vj0 == null) vj0 = N(n0, vi0);
  if (vi1 == null) vi1 = constructNormal(n1);
  if (vj1 == null) vj1 = N(n1, vi1);
  if (vi2 == null) vi2 = constructNormal(n2);
  if (vj2 == null) vj2 = N(n2, vi2);

  ArrayList<pt> points = new ArrayList<pt>();
  pt[] ps;
  ps = tangentPlaneThreeCircles(c0, r0, n0, vi0, vj0, c1, r1, n1, vi1, vj1, c2, r2, n2, vi2, vj2);
  for (pt p : ps) points.add(p);
  ps = tangentPlaneThreeCircles(c0, r0, n0, vi0, vj0, c2, r2, n2, vi2, vj2, c1, r1, n1, vi1, vj1);
  for (pt p : ps) points.add(p);

  int[] oppo = {3, 5, 4};

  pt[] centers = {c0, c1, c2};
  float[] radii = {r0, r1, r2};
  vec[] normals = {n0, n1, n2};
  vec[] vis = {vi0, vi1, vi2};
  vec[] vjs = {vj0, vj1, vj2};

  int[][] correspondences = new int[3][4];
  int nSamples = 8;
  for (int i = 0; i < 3; ++i) {
    int j = (i + 1) % 3;
    pt pa0 = points.get(i);
    pt pa1 = points.get(oppo[i]);
    pt pb1 = points.get(oppo[j]);
    pt pb0 = points.get(j);

    pt ca = centers[i], cb = centers[j];
    float ra = radii[i], rb = radii[j];
    vec na = normals[i], nb = normals[j];

    float aa = acos(dot(V(ca, pa0), V(ca, pa1)) / (ra * ra));  // [0, PI]
    if (dot(na, N(ca, pa0, pa1)) < 0) aa = TWO_PI - aa;

    float ab = acos(dot(V(cb, pb0), V(cb, pb1)) / (rb * rb));
    if (dot(nb, N(cb, pb0, pb1)) < 0) ab = TWO_PI - ab;

    // TODO: may pick the arc with bigger angle to subdivide
    // Assume that arc b has a bigger angle than arc a
    correspondences[i][0] = j;
    correspondences[i][1] = i;
    correspondences[i][2] = oppo[j];
    correspondences[i][3] = oppo[i];

    vec via = vis[i], vja = vjs[i];
    vec vib = vis[j], vjb = vjs[j];
    float angle = acos(dot(V(cb, pb0), vib) / rb);
    if (dot(V(cb, pb0), vjb) < 0) angle = TWO_PI - angle;
    float da = ab / nSamples;
    // println("ab = ", ab);
    for (int k = 1; k < nSamples; ++k) {
      angle += da;
      float c = cos(angle), s = sin(angle);
      pt p0 = P(cb, rb * c, vib, rb * s, vjb);
      pt p1 = P(p0, -rb * s, vib, rb * c, vjb);
      pt[] candidates = pivotPlaneAroundLineHitCircle(ca, ra, na, p0, p1, via, vja);
      pt candidate = candidates[0];
      if (dot(V(p0, ca), N(p0, candidate, p1)) > 0) {
        candidate = candidates[1];
      }
      points.add(p0);
      points.add(candidate);
    }
  }

  // Show the mesh
  {
    fill(red, 100);
    disk(c0, n0, r0);
    disk(c1, n1, r1);
    disk(c2, n2, r2);
  }
  assert points.size() == 6 * nSamples;
  if (show3T) {
    fill(orange, 100);
    showTriangle(points.get(0), points.get(1), points.get(2));
    showTriangle(points.get(3), points.get(4), points.get(5));
    fill(cyan, 100);
    showNormalToTriangle(points.get(0), points.get(1), points.get(2), 20, 2);
    showNormalToTriangle(points.get(3), points.get(4), points.get(5), 20, 2);
  }
  if (show2T) {
    fill(purple, 100);
    stroke(0);
    strokeWeight(2);
    for (int i = 0; i < 3; ++i) {
      beginShape(TRIANGLE_STRIP);
      int ia = correspondences[i][0];
      int ib = correspondences[i][1];
      vertex(points.get(ia));
      vertex(points.get(ib));
      int start = 6 + i * 2 * (nSamples - 1);
      int size = 2 * (nSamples - 1);
      for (int j = 0; j < size; ++j) {
        vertex(points.get(start + j));
      }
      ia = correspondences[i][2];
      ib = correspondences[i][3];
      vertex(points.get(ia));
      vertex(points.get(ib));
      endShape();
    }
    noStroke();
    for (int i = 0; i < 3; ++i) {
      int start = 6 + i * 2 * (nSamples - 1);
      int size = 2 * (nSamples - 1);
      fill(green, 100);
      for (int j = 0; j < size; j += 2) {
        show(points.get(start + j), 2);
      }
      fill(blue, 100);
      for (int j = 1; j < size; j += 2) {
        show(points.get(start + j), 2);
      }
    }
  }
}