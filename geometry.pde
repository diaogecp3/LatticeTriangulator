/******************************************************************************
 * Geometric utility functions.
 *
 * Intersection tests between geometric primitives.
 * Construction of basic geometric primitives.
 * And so on.
 ******************************************************************************/


int methodIterSP = 0;  // 0: basic, 1: heuristic normal, 2: heuristic plane

class DebugIterSPInfo {  // debug info about iterative supporting plane of 3 circles
  int iter = 0;
  DebugIterSPInfo() {}
}

/*
 * Construct the normal of triangle (A, B, C).
 */
vec normalOfTriangle(pt A, pt B, pt C) {
  return U(N(A, B, C));
}

/*
 * Compute the circumradius of a triangle (pa, pb, pc).
 */
float circumradiusOfTriangle(pt pa, pt pb, pt pc) {
  float a = d(pa, pb);
  float b = d(pb, pc);
  float c = d(pc, pa);
  return (a*b*c)/(sqrt((a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c)));
}

/*
 * Compute the circumcenter of a triangle (pa, pb, pc).
 */
pt circumcenterOfTriangle(pt pa, pt pb, pt pc) {
  vec ab = V(pa, pb);
  vec ac = V(pa, pc);
  vec abac = N(ab, ac);
  vec v = V(n2(ac), N(abac, ab), n2(ab), N(ac, abac));
  v.div(2 * n2(abac));
  return P(pa, v);
}

/*
 * Construct a vector normal to given vector v. v is not necessarily a unit
 * vector. A unit vector will be returned.
 */
vec constructNormal(vec v) {
  if (isZero(v.y) && isZero(v.z)) {  // v is parallel to (1, 0, 0)
    return U(new vec(-v.z, 0.0, v.x));  // cross product of v and (0, 1, 0)
  } else {
    return U(new vec(0.0, v.z, -v.y));  // cross product of v and (1, 0, 0)
  }
}

/*
 * Compute the intersection line of two planes, defined by:
 * a0x + b0y + c0z = d0
 * a1x + b1y + c1z = d1
 */
boolean intersectionTwoPlanes(float a0, float b0, float c0, float d0,
                              float a1, float b1, float c1, float d1,
                              pt p0, pt p1) {
  if (notZero(b0 * c1 - c0 * b1)) {
    vec2 t0 = solveLinearEquationsInTwoVars(b0, c0, b1, c1, d0, d1);
    vec2 t1 = solveLinearEquationsInTwoVars(b0, c0, b1, c1, d0 - a0, d1 - a1);
    p0.set(0.0, t0.x, t0.y);
    p1.set(1.0, t1.x, t1.y);
  } else if (notZero(a0 * c1 - c0 * a1)) {
    vec2 t0 = solveLinearEquationsInTwoVars(a0, c0, a1, c1, d0, d1);
    vec2 t1 = solveLinearEquationsInTwoVars(a0, c0, a1, c1, d0 - b0, d1 - b1);
    p0.set(t0.x, 0.0, t0.y);
    p1.set(t1.x, 1.0, t1.y);
  } else if (notZero(a0 * b1 - b0 * a1)) {
    vec2 t0 = solveLinearEquationsInTwoVars(a0, b0, a1, b1, d0, d1);
    vec2 t1 = solveLinearEquationsInTwoVars(a0, b0, a1, b1, d0 - c0, d1 - c1);
    p0.set(t0.x, t0.y, 0.0);
    p1.set(t1.x, t1.y, 1.0);
  } else {
    println("No intersection between these two planes!");
    return false;
  }
  return true;
}

/*
 * Compute the intersection line of two planes.
 */
boolean intersectionTwoPlanes(pt pa,                                     // in
                              vec va,                                    // in
                              pt pb,                                     // in
                              vec vb,                                    // in
                              pt p0,                                    // out
                              pt p1) {                                  // out
  float a0 = va.x, b0 = va.y, c0 = va.z, d0 = dot(pa, va);
  float a1 = vb.x, b1 = vb.y, c1 = vb.z, d1 = dot(pb, vb);
  return intersectionTwoPlanes(a0, b0, c0, d0, a1, b1, c1, d1, p0, p1);
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
    if (isZero(dot(V(ca, cb), va))) {  // two disks on the same plane
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
  if (notZero(dnab)) {  // dot(N, AB) != 0
    float t = dnac / dnab;
    pt p = P(a, t, ab);
    if (d(p, c) <= r) return false;
    else return true;
  } else {  // dot(N, AB) == 0
    if (notZero(dnac)) {  // dot(N, AC) != 0
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
 * Compute the intersection points between a line (o, d) and a sphere (c, r).
 * d is not necessarily a unit vector.
 */
pt[] intersectionLineSphere(pt o, vec d, pt c, float r) {
  vec co = V(c, o);
  float c2 = dot(d, d);
  float c1 = 2 * dot(co, d);
  float c0 = dot(co, co) - r * r;
  float[] ts = solveQuadraticEquation(c2, c1, c0);
  if (ts == null) return null;
  return new pt[] {P(o, ts[0], d), P(o, ts[1], d)};
}

/*
 * Return the points of contact between the plane that passes through line (a, b)
 * and the circle (c, r, n). (vi, vj), which is optional, is the local frame on
 * the circle. May be deprecated.
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
  if (isZero(dnab) && isZero(dnca)) {
    println("Line AB is on the supporting plane of Circle C");
    ps[0] = P(c, r, vi);
    ps[1] = P(c, -r, vi);
    return ps;
  }
  if (vj == null) vj = N(n, vi);  // cross product
  float[] thetas;
  if (notZero(dnab)) {
    vec v = A(ca, -dnca/dnab, ab);
    float fa = dot(vi, v);
    float fb = dot(vj, v);
    thetas = solveLinearEquationInCosSin(fa, fb, r);
  } else {
    float fa = dot(vi, ab);
    float fb = dot(vj, ab);
    thetas = solveLinearEquationInCosSin(fa, fb, 0);
  }

  assert thetas.length == 2;
  ps[0] = P(c, r * cos(thetas[0]), vi, r * sin(thetas[0]), vj);
  ps[1] = P(c, r * cos(thetas[1]), vi, r * sin(thetas[1]), vj);
  return ps;
}

/*
 * Return the points of contact between the plane that passes through line (a, b)
 * and the circle (c, r, n). (vi, vj), which is optional, is the local frame on
 * the circle.
 */
pt[] contactsOnSupportingPlaneOfLineCircle(pt c, float r, vec n, pt a, pt b, vec vi, vec vj) {
  pt[] ps = new pt[2];
  boolean empty = emptyIntersectionLineDisk(a, b, c, r, n);
  if (!empty) {
    println("Line AB intersects disk C");
    ps[0] = P(c, r, vi);
    ps[1] = P(c, -r, vi);
    return ps;
  }

  vec ac = V(a, c);
  vec ab = V(a, b);
  if (vi == null) vi = constructNormal(n);
  float dnab = dot(n, ab);
  float dnac = dot(n, ac);
  if (isZero(dnab) && isZero(dnac)) {
    println("Line AB is on the supporting plane of Circle C");
    ps[0] = P(c, r, vi);
    ps[1] = P(c, -r, vi);
    return ps;
  }

  if (vj == null) vj = N(n, vi);  // cross product
  float t = dnac / dnab;
  pt tmp = P(a, t, ab);
  float cosa = r / d(c, tmp);  // 0 < cosa <= 1
  float alpha = acos(cosa);  // [0, PI/2)
  float sina = sin(alpha);
  vec v0 = U(c, tmp);
  vec v1 = R(v0, HALF_PI, vi, vj);
  float rc = r * cosa;
  float rs = r * sina;
  ps[0] = P(c, rc, v0, rs, v1);
  ps[1] = P(c, rc, v0, -rs, v1);
  return ps;
}

/*
 * Construct and show the exact convex hull defined by an edge (a, b) and a
 * circle (c, r, n). (vi, vj), which is optional, is local frame on the circle.
 */
void exactCHEdgeCircle(pt a, pt b, pt c, float r, vec n, vec vi, vec vj) {
  if (vi == null) vi = constructNormal(n);
  if (vj == null) vj = N(n, vi);
  vec ca = V(c, a);
  vec cb = V(c, b);
  pt[] ps = pivotPlaneAroundLineHitCircle(c, r, n, a, b, vi, vj);

  int nSamples = 16;
  float da = TWO_PI / nSamples;
  int nSamples0 = -1;
  int nSamples1 = -1;
  ArrayList<pt> points = new ArrayList<pt>();
  float dnca = dot(n, ca);
  float dncb = dot(n, cb);
  if (dnca * dncb > 0) {  // a and b on the same side of plane (c, n)
    // println("A and B on the same side of plane (C, N)");
    if (dnca > 0) {  // a and b on the positive side of plane (c, n)
      n.rev();
      vj.rev();
    }
    vec ab = V(a, b);
    if (dot(ab, n) < 0) {  // make sure dot(ab, n) >= 0
      pt tmp = a;
      a = b;
      b = tmp;
      ab.rev();
      // Do we need ca and cb any more? If yes, then we need to swap them.
    }
    if (dot(ab, N(b, ps[0], ps[1])) < 0) {  // make sure a is on the positive side of plane (b, ps[0], ps[1])
      pt tmp = ps[0];
      ps[0] = ps[1];
      ps[1] = tmp;
    }
    vec v0 = V(c, ps[0]);
    vec v1 = V(c, ps[1]);
    float angle0 = acos(dot(v0, vi) / r);
    if (dot(v0, vj) < 0) angle0 = TWO_PI - angle0;
    float angle1 = acos(dot(v1, vi) / r);
    if (dot(v1, vj) < 0) angle1 = TWO_PI - angle1;
    float angle01 = acos(dot(v0, v1) / (r * r));
    if (dot(N(v0, v1), n) < 0) angle01 = TWO_PI - angle01;
    float angle10 = TWO_PI - angle01;

    // connect arc ps[1]-ps[0] (w.r.t. -n) to a
    float angle = angle1;
    nSamples0 = int(angle01 / da);  // note that here we use angle01
    for (int i = 1; i < nSamples0; ++i) {  // insert nSamples0 - 1 points
      angle -= da;  // notice here we use minus!
      points.add(P(c, r * cos(angle), vi, r * sin(angle), vj));
    }

    // connect arc ps[0]-ps[1] (w.r.t. -n) to b
    angle = angle0;
    nSamples1 = int(angle10 / da);
    for (int i = 1; i < nSamples1; ++i) {  // insert nSamples1 - 1 points
      angle -= da;  // notice here we use minus!
      points.add(P(c, r * cos(angle), vi, r * sin(angle), vj));
    }

    {
      fill(magenta, 200);
      beginShape(TRIANGLE_FAN);
      vertex(a);
      vertex(ps[1]);
      for (int i = 0; i < nSamples0 - 1; ++i) vertex(points.get(i));
      vertex(ps[0]);
      endShape();
      show(ps[1], 3);

      fill(orange, 200);
      beginShape(TRIANGLE_FAN);
      vertex(b);
      vertex(ps[0]);
      for (int i = nSamples0 - 1; i < points.size(); ++i) vertex(points.get(i));
      vertex(ps[1]);
      endShape();
      show(ps[0], 3);

      fill(cyan, 200);
      showTriangle(a, b, ps[1]);
      showTriangle(a, ps[0], b);
    }
  } else if (dnca * dncb < 0) {  // a and b on the different side of plane (c, n)
    if (emptyIntersectionLineDisk(a, b, c, r, n)) {
      println("A and B on different side of plane (C, N) but AB doesn't intersect disk C");
      if (dnca > 0) {
        pt tmp = a;
        a = b;
        b = tmp;
      }
      vec ab = V(a, b);
      if (dot(ab, N(b, ps[0], ps[1])) < 0) {  // make sure a is on the positive side of plane (b, ps[0], ps[1])
        pt tmp = ps[0];
        ps[0] = ps[1];
        ps[1] = tmp;
      }
      vec v0 = V(c, ps[0]);
      vec v1 = V(c, ps[1]);
      float angle0 = acos(dot(v0, vi) / r);
      if (dot(v0, vj) < 0) angle0 = TWO_PI - angle0;
      float angle1 = acos(dot(v1, vi) / r);
      if (dot(v1, vj) < 0) angle1 = TWO_PI - angle1;
      float angle01 = acos(dot(v0, v1) / (r * r));
      if (dot(N(v0, v1), n) < 0) angle01 = TWO_PI - angle01;

      // connect arc ps[1]-ps[0] (w.r.t. -n) to a
      float angle = angle1;
      nSamples0 = int(angle01 / da);
      for (int i = 1; i < nSamples0; ++i) {
        angle -= da;
        points.add(P(c, r * cos(angle), vi, r * sin(angle), vj));
      }
      // connect arc ps[0]-ps[1] (w.r.t. n) to b
      // reuse the points generated for the previous arc

      {
        fill(magenta, 200);
        beginShape(TRIANGLE_FAN);
        vertex(a);
        vertex(ps[1]);
        for (int i = 0; i < nSamples0 - 1; ++i) vertex(points.get(i));
        vertex(ps[0]);
        endShape();
        show(ps[1], 3);

        fill(orange, 200);
        beginShape(TRIANGLE_FAN);
        vertex(b);
        vertex(ps[0]);
        for (int i = nSamples0 - 2; i >= 0; --i) vertex(points.get(i));
        vertex(ps[1]);
        endShape();
        show(ps[0], 3);

        fill(cyan, 200);
        showTriangle(a, b, ps[1]);
        showTriangle(a, ps[0], b);
      }
    } else {  // line (a, b) intersects disk (c, r, n)
      println("A and B on different side of plane (C, N) and AB intersects disk C");
      if (dnca > 0) {
        pt tmp = a;
        a = b;
        b = tmp;
      }
      float angle = 0;
      for (int i = 0; i < nSamples; ++i, angle += da) {
        points.add(P(c, r * cos(angle), vi, r * sin(angle), vj));
      }
      {
        fill(magenta, 200);
        beginShape(TRIANGLE_FAN);
        vertex(a);
        for (int i = nSamples - 1; i >= 0; --i) vertex(points.get(i));
        vertex(points.get(nSamples - 1));
        endShape();

        fill(orange, 200);
        beginShape(TRIANGLE_FAN);
        vertex(b);
        for (int i = 0; i < nSamples; ++i) vertex(points.get(i));
        vertex(points.get(0));
        endShape();
      }
    }
  }

  {
    fill(red, 200);
    disk(c, n, r);
    fill(green, 200);
    show(a, 3);
    fill(blue, 200);
    show(b, 3);
  }
}

/*
 * Return the 3 points of contact that define the oriented supporting plane of
 * the 3 circles (c0, r0, n0), (c1, r1, n1), and (c2, r2, n2). (vi0, vj0),
 * (vi1, vj1) and (vi2, vj2), which are optional, are local frames on the 3 circles,
 * respectively. This function uses an iterative method.
 */
pt[] supPlaneThreeCirclesIter(pt c0, float r0, vec n0, vec vi0, vec vj0,
                            pt c1, float r1, vec n1, vec vi1, vec vj1,
                            pt c2, float r2, vec n2, vec vi2, vec vj2,
                            DebugIterSPInfo dInfo) {
  if (vi0 == null) vi0 = constructNormal(n0);
  if (vj0 == null) vj0 = N(n0, vi0);
  if (vi1 == null) vi1 = constructNormal(n1);
  if (vj1 == null) vj1 = N(n1, vi1);
  if (vi2 == null) vi2 = constructNormal(n2);
  if (vj2 == null) vj2 = N(n2, vi2);

  /* Initialize a triangle with each vertex from each circle. */
  pt p0, p1, p2;
  vec t0, t1, t2, n;
  // println("methodIterSP = ", methodIterSP);
  switch (methodIterSP) {
    case 1:
      vec tmpN1 = normalOfTriangle(c0, c1, c2);
      p0 = P(c0, r0, U(A(tmpN1, -dot(tmpN1, n0), n0)));
      p1 = P(c1, r1, U(A(tmpN1, -dot(tmpN1, n1), n1)));
      p2 = P(c2, r2, U(A(tmpN1, -dot(tmpN1, n2), n2)));
      t0 = U(N(n0, V(c0, p0)));
      t1 = U(N(n1, V(c1, p1)));
      t2 = U(N(n2, V(c2, p2)));
      n = normalOfTriangle(p0, p1, p2);
      break;
    case 2:
      vec tmpN2 = normalOfTriangle(c0, c1, c2);
      p0 = new pt();
      intersectionCirclePlane(c0, r0, n0, c0, tmpN2, p0, null);
      p1 = new pt();
      intersectionCirclePlane(c1, r1, n1, c0, tmpN2, p1, null);
      p2 = new pt();
      intersectionCirclePlane(c2, r2, n2, c0, tmpN2, p2, null);
      t0 = U(N(n0, V(c0, p0)));
      t1 = U(N(n1, V(c1, p1)));
      t2 = U(N(n2, V(c2, p2)));
      n = normalOfTriangle(p0, p1, p2);
      break;
    default:
      p0 = P(c0, r0, vi0);
      p1 = P(c1, r1, vi1);
      p2 = P(c2, r2, vi2);
      t0 = vj0;
      t1 = vj1;
      t2 = vj2;
      n = normalOfTriangle(p0, p1, p2);
  }

  int maxIter = 10;
  int iter = 0;
  while (iter < maxIter) {
    boolean update = false;
    // move A, pivot around BC
    if (notZero(dot(n, t0)) || dot(n, V(p1, c0)) > 0) {
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
    if (notZero(dot(n, t1)) || dot(n, V(p2, c1)) > 0) {
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
    if (notZero(dot(n, t2)) || dot(n, V(p0, c2)) > 0) {
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
  if (dInfo != null) dInfo.iter = iter;

  pt[] ps = new pt[] {p0, p1, p2};
  // println("ps.length =", ps.length);
  return ps;
}

/*
 * Construct and show the exact convex hull defined by two circles (c0, r0, n0)
 * and (c1, r1, n2). (vi0, vj0) and (vi1, vj1), which are optional, are local
 * frames on the two circles, respectively.
 */
void exactCHTwoCircles(pt c0, float r0, vec n0, vec vi0, vec vj0,
                       pt c1, float r1, vec n1, vec vi1, vec vj1) {
  if (vi0 == null) vi0 = constructNormal(n0);
  if (vj0 == null) vj0 = N(n0, vi0);
  if (vi1 == null) vi1 = constructNormal(n1);
  if (vj1 == null) vj1 = N(n1, vi1);

  int nSamples = 16;
  ArrayList<pt> points = new ArrayList<pt>();
  float da = TWO_PI / nSamples;
  float a = 0;
  for (int i = 0; i < nSamples; ++i, a += da) {
    float s = sin(a);
    float c = cos(a);
    pt p0 = P(c0, r0 * c, vi0, r0 * s, vj0);
    pt p1 = P(p0, -r0 * s, vi0, r0 * c, vj0);
    pt[] candidates = pivotPlaneAroundLineHitCircle(c1, r1, n1, p0, p1, vi1, vj1);
    pt candidate = candidates[0];
    if (dot(V(p0, c1), N(p0, candidate, p1)) > 0) {
      candidate = candidates[1];
    }
    points.add(p0);
    points.add(candidate);
  }

  assert points.size() == 2 * nSamples;

  // show the mesh
  {
    fill(red, 100);
    disk(c0, n0, r0);
    fill(green, 100);
    disk(c1, n1, r1);
  }
  {
    fill(violet, 100);
    stroke(0);
    strokeWeight(2);
    beginShape(QUAD_STRIP);
    for (pt p : points) vertex(p);
    vertex(points.get(0));
    vertex(points.get(1));
    endShape();
  }
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
  ps = supPlaneThreeCirclesIter(c0, r0, n0, vi0, vj0, c1, r1, n1, vi1, vj1, c2, r2, n2, vi2, vj2, null);
  for (pt p : ps) points.add(p);
  ps = supPlaneThreeCirclesIter(c0, r0, n0, vi0, vj0, c2, r2, n2, vi2, vj2, c1, r1, n1, vi1, vj1, null);
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

    // TODO: may pick the arc with bigger angle to split
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
      pt p0 = P(cb, rb * c, vib, rb * s, vjb);  // sampled point
      pt p1 = P(p0, -rb * s, vib, rb * c, vjb);  // sampled point + tangent at that point
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
  if (show3RT) {
    fill(orange, 100);
    showTriangle(points.get(0), points.get(1), points.get(2));
    showTriangle(points.get(3), points.get(4), points.get(5));
    // fill(cyan, 100);
    // showNormalToTriangle(points.get(0), points.get(1), points.get(2), 20, 2);
    // showNormalToTriangle(points.get(3), points.get(4), points.get(5), 20, 2);
  }
  if (show2RT) {
    fill(violet, 100);
    stroke(0);
    strokeWeight(2);
    for (int i = 0; i < 3; ++i) {
      beginShape(QUAD_STRIP);
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

/*
 * Return true if the 4 points are coplanar.
 */
boolean coplanarFourPoints(pt pa, pt pb, pt pc, pt pd) {
  vec vb = U(pa, pb);
  vec vc = U(pa, pc);
  vec vd = U(pa, pd);
  float mix = m(vb, vc, vd);
  if (isZero(mix)) return true;
  else {
    println("mix product =", mix);
    return false;
  }
}

/*
 * Compute the distance from a point p to a round cone, which is defined by the
 * union of two balls and the tangential cone connecting them. The distance is
 * negative if the point is inside the round cone.
 * https://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm
 */
float roundConeDist(pt p, pt s1, float r1, pt s2, float r2) {
  vec s1s2 = V(s1, s2);
  vec s1p = V(s1, p);
  float h = s1s2.norm();
  vec2 q = new vec2(N(s1p, s1s2).norm() / h, dot(V(s1s2).div(h), s1p));

  float b = (r1 - r2) / h;
  float a = sqrt(1.0 - b * b);
  float k = dot(q, V(-b, a));

  if (k < 0.0) return q.norm() - r1;
  if (k > a * h) return V(q.x, q.y-h).norm() - r2;

  return dot(q, V(a,b)) - r1;
}
