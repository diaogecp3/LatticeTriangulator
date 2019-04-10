/******************************************************************************
 * Utility functions.
 *
 * Statistics, rendering, math, etc.
 ******************************************************************************/

float gEpsilon = 0.000001;
float gEpsilonBig = 0.0001;
boolean exitDraw = false;

/* Statistics functions below. */

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
 * Compute accuracy for range [start, end) of array a. If array valids is not
 * null, only the valid values (indicated by array valids) of array a will be
 * considered.
 */
float accuracy(boolean[] a, int n, int start, int end, boolean[] valids) {
  assert start >= 0 && end <= n && end > start;
  float sum = 0.0;
  if (valids == null) {
    for (int i = start; i < end; ++i) {
      sum += (a[i] ? 1 : 0);
    }
    return sum / n;
  } else {
    int k = 0;
    for (int i = start; i < end; ++i) {
      if (valids[i]) {
        sum += (a[i] ? 1 : 0);
        k++;
      }
    }
    return sum / k;
  }
}


/* Math functions below. */

/*
 * Return true iff |x| is not close to 0.
 */
boolean notZero(float x) {
  return x <= -gEpsilon || x >= gEpsilon;
}

/*
 * Return true iff |x| >= e.
 */
boolean notZero(float x, float e) {
  return x <= -e || x >= e;
}


/*
 * Return true iff |x| is close to 0.
 */
boolean isZero(float x) {
  return x > -gEpsilon && x < gEpsilon;
}

/*
 * Return true iff |x| < e.
 */
boolean isZero(float x, float e) {
  return x > -e && x < e;
}

/*
 * Return true iff x is non-positive.
 */
boolean isNonPositive(float x) {
  return x < gEpsilon;
}

/*
 * Return acos of clamped value of x.
 */
float acosClamp(float x) {
  return acos(constrain(x, -1.0, 1.0));
}

/*
 * Return asin of clamped value of x.
 */
float asinClamp(float x) {
  return asin(constrain(x, -1.0, 1.0));
}

boolean isNaNVec(vec v) {
  return Float.isNaN(v.x) || Float.isNaN(v.y) || Float.isNaN(v.z);
}

boolean isNaNPt(pt p) {
  return Float.isNaN(p.x) || Float.isNaN(p.y) || Float.isNaN(p.z);
}

boolean isZeroVec(vec v) {
  return isZero(v.x) && isZero(v.y) && isZero(v.z);
}

boolean isZeroVec(vec v, float e) {
  return isZero(v.x, e) && isZero(v.y, e) && isZero(v.z, e);
}

boolean isZeroPt(pt p) {
  return isZero(p.x) && isZero(p.y) && isZero(p.z);
}

boolean isZeroPt(pt p, float e) {
  return isZero(p.x, e) && isZero(p.y, e) && isZero(p.z, e);
}

boolean samePt(pt p, pt q) {
  return isZero(p.x - q.x) && isZero(p.y - q.y) && isZero(p.z - q.z);
}

boolean samePt(pt p, pt q, float e) {
  return isZero(p.x - q.x, e) && isZero(p.y - q.y, e) && isZero(p.z - q.z, e);
}

boolean isApproximately(float a, float b, float epsilon) {
  return a + epsilon > b && a - epsilon < b;
}

boolean isApproximately(double a, double b, double epsilon) {
  return a + epsilon > b && a - epsilon < b;
}

boolean isApproximately(vec a, vec b, float epsilon) {
  return isApproximately(a.x, b.x, epsilon) && isApproximately(a.y, b.y, epsilon) && isApproximately(a.z, b.z, epsilon);
}

boolean isApproximately(MinMaxF a, MinMaxF b, float epsilon) {
  return isApproximately(a.min, b.min, epsilon) && isApproximately(a.max, b.max, epsilon);
}

float logb(float base, float val) { return log(val) / log(base); }  // log base b

double logb(double base, double val) { return Math.log(val) / Math.log(base); }


/*
 * Return the index whose value is minimum in an array list. If there are more
 * than one minimums, return the lowest index.
 */
int argmin(ArrayList<Integer> a) {
  int ret = 0, min = a.get(0);
  for (int i = 1; i < a.size(); ++i) {
    if (a.get(i) < min) {
      min = a.get(i);
      ret = i;
    }
  }
  return ret;
}

int argAbsMax(vec v) {
  float a = abs(v.x);
  float b = abs(v.y);
  float c = abs(v.z);
  if (a >= b && a >= c) return 0;
  if (b >= a && b >= c) return 1;
  return 2;
}

float min3(float a, float b, float c) { return min(a, min(b, c)); }
float max3(float a, float b, float c) { return max(a, max(b, c)); }

MinMaxI intersect(MinMaxI a, MinMaxI b) {
  MinMaxI interval = new MinMaxI( a.min>b.min?a.min:b.min, a.max<b.max?a.max:b.max );
  if (interval.min > interval.max) interval.max = interval.min;
  return interval;
}

MinMaxI bound(MinMaxI a, MinMaxI b) {
  MinMaxI interval = new MinMaxI( a.min<b.min?a.min:b.min, a.max>b.max?a.max:b.max );
  return interval;
}

MinMaxF intersect(MinMaxF a, MinMaxF b) {
  MinMaxF interval = new MinMaxF(a.min>b.min?a.min:b.min, a.max<b.max?a.max:b.max);
  if (interval.min > interval.max) interval.max = interval.min;
  return interval;
}

MinMaxF bound(MinMaxF a, MinMaxF b) {
  MinMaxF interval = new MinMaxF(a.min<b.min?a.min:b.min, a.max>b.max?a.max:b.max);
  return interval;
}

/*
 * Return true if a and b point the same direction. Assume that a and b are colinear.
 */
boolean sameDirection(vec a, vec b) {
  int i = argAbsMax(a);  // pick the dimension with largest abs value for numerical stability
  return a.get(i) * b.get(i) > 0;
}

/*
 * Generate a random vector, with each element sampled from [-n, n).
 */
vec random3(float n) {
  return new vec(random(-n, n), random(-n, n), random(-n, n));
}

/*
 * Solve the following linear equations in two variables.
 *
 * ax + by = e
 * cx + dy = f
 */
vec2 solveLinearEquationsInTwoVars(float a, float b, float c, float d, float e, float f) {
  float det = a * d - b * c;
  if (isZero(det)) return null;
  float x =(d * e - b * f) / det;
  float y = (a * f - c * e) / det;
  return new vec2(x, y);
}

/*
 * Solve ax^2 + bx + c = 0, return two solutions if there exists.
 */
float[] solveQuadraticEquation(float a, float b, float c) {
  float det = b * b - 4 * a * c;
  if (det < 0) return null;
  float[] xs = new float[2];
  float sdet = sqrt(det);
  float aa = 2 * a;
  xs[0] = (-b + sdet) / aa;
  xs[1] = (-b - sdet) / aa;
  return xs;
}

/*
 * Solve a*cos(theta) + b*sin(theta) = c, return two thetas. The two thetas will
 * be in [0, 2pi].
 */
float[] solveLinearEquationInCosSin(float a, float b, float c) {
  if (a * a + b * b < 0.0000001) {
    println("In a*cos(theta) + b*sin(theta) = c, a and b might be too small!");
  }
  float[] thetas = new float[2];
  if (notZero(a) && notZero(b)) {  // a != 0 and b != 0
    // println("a = ", a, "b = ", b);
    /* How to decide (the sign of) tmp0 here? According to the sign of b. */
    float tmp0 = asinClamp(c / (sqrt(a * a + b * b)));  // [-pi/2, pi/2]
    if (b < 0) tmp0 = -tmp0;  // [-pi/2, pi/2]
    float tmp1 = atan(a / b);  // (-pi/2, pi/2)
    thetas[0] = tmp0 - tmp1;  // (-pi, pi)
    thetas[1] = PI - tmp0 - tmp1;  // (0, 2pi)
    thetas[0] = thetas[0] < 0 ? thetas[0] + TWO_PI : thetas[0];  // [0, pi), (pi, 2pi)
  } else if (isZero(a)) {  // a == 0 and b != 0
    println("In solving a*cos(theta) + b*sin(theta) = c: a == 0");
    thetas[0] = asinClamp(c/b);  // [-pi/2, pi/2]
    thetas[1] = PI - thetas[0];  // [pi/2, 3pi/2]
    if (thetas[0] < 0) thetas[0] += TWO_PI;  // [3pi/2, 0), [0, pi/2]
  } else {  // a != 0 and b == 0
    println("In solving a*cos(theta) + b*sin(theta) = c: b == 0");
    thetas[0] = acosClamp(c/a);  // [0, pi]
    thetas[1] = TWO_PI - thetas[0];  // [pi, 2pi]
  }
  return thetas;
}

/*
 * Solve a*cos(theta) + b*sin(theta) = 0, return two thetas. The two thetas will
 * be in [0, 2pi].
 */
float[] solveLinearEquationInCosSin(float a, float b) {
  if (isZero(a)) return new float[] {0.0, PI};
  if (isZero(b)) return new float[] {HALF_PI, HALF_PI + PI};
  float t = - a / b;
  if (Float.isNaN(t)) {
    println("t is NaN!");
    return null;
  }
  float theta = atan(t);
  if (theta >= 0) return new float[] {theta, theta + PI};
  else return new float[] {theta + PI, theta + TWO_PI};
}


/* Rendering functions below. */

/*
 * Show a triangle.
 */
void showTriangle(pt a, pt b, pt c) {
  beginShape(TRIANGLES);
  vertex(a);
  vertex(b);
  vertex(c);
  endShape();
}

/*
 * Show triangles.
 */
void showTriangles(ArrayList<Triangle> triangles, pt[] points) {
  int n = triangles.size();
  beginShape(TRIANGLES);
  for (int i = 0; i < n; ++i) {
    if (triangles.get(i) == null) continue;
    vertex(points[triangles.get(i).a]);
    vertex(points[triangles.get(i).b]);
    vertex(points[triangles.get(i).c]);
  }
  endShape();
}

/*
 * Show the normal to triangle (a, b, c) as an arrow. d controls the length of
 * the arrow and r controls the size of the arrow tip.
 */
void showTriangleNormal(pt a, pt b, pt c, float d, float r) {
  vec n = normalOfTriangle(a, b, c);
  pt m = P(a, b, c);
  arrow(m, V(d, n), r);
}

/*
 * Show normals of all triangles.
 */
void showTriangleNormals(ArrayList<Triangle> triangles, pt[] points) {
  for (Triangle tri : triangles) {
    if (tri == null) continue;
    showTriangleNormal(points[tri.a], points[tri.b], points[tri.c], 10, 1);
  }
}

/*
 * Show the ball defined by (c, r).
 */
void showBall(pt c, float r) {
  pushMatrix();
  translate(c.x, c.y, c.z);
  sphere(r);
  popMatrix();
}

/*
 * Show the line defined by (o, d). d is not necessarily a unit vector.
 */
void showLine(pt o, vec d) {
  pt p0 = P(o, -100.0, d);
  pt p1 = P(o, 100.0, d);
  line(p0.x, p0.y, p0.z, p1.x, p1.y, p1.z);
}

/*
 * Show the line defined by (a, b).
 */
void showLine(pt a, pt b) {
  vec d = V(a, b);
  showLine(a, d);
}

/*
 * Show the line segment defined by (a, b).
 */
void showSegment(pt a, pt b) {
  line(a.x, a.y, a.z, b.x, b.y, b.z);
}

/*
 * Show the poly line defined by a list of points.
 */
void showPolyLine(ArrayList<pt> ps) {
  noFill();
  beginShape();
  for (pt p : ps) vertex(p);
  endShape();
  fill(black);
}

/*
 * Show the poly-loop defined by a list of points.
 */
void showLoop(ArrayList<pt> ps) {
  noFill();
  beginShape();
  for (pt p : ps) vertex(p);
  endShape(CLOSE);
  fill(black);
}

/*
 * Show the oriented poly-loop defined by a list of points.
 */
void showOrientedLoop(ArrayList<pt> ps) {
  pt prev = ps.get(ps.size()-1);  // last point
  pt center = new pt();
  for (int i = 0; i < ps.size(); ++i) {
    pt cur = ps.get(i);
    arrow(prev, V(prev, cur), 1);
    prev = cur;
    center.add(cur);
  }
  center.div(ps.size());
  vec normal = U(N(ps.get(0), ps.get(1), ps.get(2)));
  arrow(center, V(8, normal), 1);
}

/*
 * Show the plane defined by (p, n). The plane is a finite square centered at p,
 * with side length being 2s.
 */
void showPlane(pt p, vec n, float s) {
  vec vi = constructNormal(n);
  vec vj = N(n, vi);
  beginShape(QUAD);
  vertex(P(p, s, vi, s, vj));
  vertex(P(p, s, vi, -s, vj));
  vertex(P(p, -s, vi, -s, vj));
  vertex(P(p, -s, vi, s, vj));
  endShape();
}

/*
 * Show the plane defined by (a, b, c). The plane is a finite square centered at
 * the center of the triangle (a, b, c), with side length being 2s.
 */
void showPlane(pt a, pt b, pt c, float s) {
  vec n = normalOfTriangle(a, b, c);
  pt p = P(a, b, c);
  showPlane(p, n, s);
}

/*
 * Show the circle defined by (c, n, r) where c is the center, n is the normal,
 * r is the radius. Remember to set stroke color before calling this function.
 */
void showCircle(pt c, vec n, float r) {
  noFill();
  float a = 0;
  float da = TWO_PI / 36;
  vec vi = constructNormal(n);
  vec vj = N(n, vi);
  beginShape();
  for (int i = 0; i < 36; ++i, a += da) {
    vertex(P(c, r * cos(a), vi, r * sin(a), vj));
  }
  endShape(CLOSE);
  fill(black);
}

/*
 * Show the circumcircle defined by triangle (pa, pb, pc).
 */
void showCircumcircleOfTriangle(pt pa, pt pb, pt pc, pt center, vec normal, Float radius) {
  if (center == null) center = circumcenterOfTriangle(pa, pb, pc);
  if (normal == null) normal = normalOfTriangle(pa, pb, pc);
  if (radius == null) radius = d(center, pa);
  showCircle(center, normal, radius);
}

/*
 * Show the oblique cone with apex p, and a circular base (c, n, r).
 */
void showObliqueCone(pt p, pt c, vec n, float r) {
  int nSamples = 120;
  float a = 0;
  float da = TWO_PI / nSamples;
  vec vi = constructNormal(n);
  vec vj = N(n, vi);
  beginShape(TRIANGLE_FAN);
  vertex(p);
  for (int i = 0; i <= nSamples; ++i, a += da) {
    vertex(P(c, r * cos(a), vi, r * sin(a), vj));
  }
  endShape();
}

/*
 * Spherical linear interpolation.
 */
vec slerp(vec U, float t, vec W) {
  float a = angle(U, W);
  float su = sin(a * (1.0 - t));
  float sw = sin(a * t);
  float s = sin(a);
  return V(su / s, U, sw / s, W);
}

/*
 * Sample a point on arc, centered at C, between A and B. t is the parameter
 * of the arc.
 */
pt onArc(pt A, pt C, pt B, float t) {
  vec U = V(C, A);
  vec W = V(C, B);
  vec V = slerp(U, t, W);
  return P(C, V);
}

/*
 * Draw an arc on sphere, centered at C, from A to B using weignt w.
 */
void showArc(pt A, pt C, pt B, float w) {
  pt PP = P(A);
  for(float t = 0.05; t <= 1.01; t += 0.05) {
    pt QQ = onArc(A, C, B, t);
    stub(PP, QQ, w, w);
    PP.set(QQ);
  }
}

/*
 * Show some arcs on sphere (c, r). Each arc is (ps[i], ps[i+1]) with weight w.
 */
void showArcs(pt[] ps, int nv, pt c, float r, float w0, float w1, float w2) {
  assert nv % 2 == 0;
  fill(red);
  for (int i = 0; i < nv; i += 2) {
    show(ps[i], w0);
  }

  fill(green);
  for (int i = 0; i < nv; i += 2) {
    showArc(ps[i], c, ps[i+1], w1);
  }

  fill(blue);
  for (int i = 1; i < nv; i += 2) {
    show(ps[i], w2);
  }
}

/*
 * Show the (minor) arc from A to B on sphere (C, r) approximated by a polyline.
 */
void showPolyArc(pt a, pt b, pt c, float r) {
  vec u = V(c, a), w = V(c, b);
  float angle = acosClamp(dot(u, w) / (r * r));
  float s = sin(angle);

  ArrayList<pt> points = new ArrayList<pt>();
  points.add(a);
  for(float t = 0.1; t < 1.01; t += 0.1) {  // 10 iterations
    float su = sin(angle * (1 - t));
    float sw = sin(angle * t);
    vec v = V(su/s, u, sw/s, w);  // SLERP
    points.add(P(c, v));
  }

  stroke(0);
  strokeWeight(1);
  noFill();
  beginShape();
  for (pt p : points) vertex(p);
  endShape();
  fill(black);
  noStroke();
}


/* Misc functions below. */

/*
 * Exception handler.
 */
void exceptionHandler() {
  gPoints.savePts("data/pts_unnamed");
  gRingSet.save("data/rs_unnamed");
  exitDraw = true;
}

/*
 * Keep only one point if multiple points are almost the same. Assume that points
 * are in cyclic order. pids are the IDs of points.
 */
void removeDuplicates(ArrayList<pt> points, ArrayList<Integer> pids) {
  if (points == null && pids == null) return;
  assert points.size() == pids.size();
  int i = 1;
  while (i < points.size()) {
    if (samePt(points.get(i), points.get(i-1))) {
      points.remove(i);
      pids.remove(i);
    } else i++;
  }
  if (samePt(points.get(i-1), points.get(0))) {
    points.remove(i-1);
    pids.remove(i-1);
  }
}