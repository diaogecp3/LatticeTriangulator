/******************************************************************************
 * Utility functions.
 *
 * Statistics, rendering, math, etc.
 ******************************************************************************/


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
  return x <= -0.00001 || x >= 0.00001;
}

/*
 * Return true iff |x| is close to 0.
 */
boolean isZero(float x) {
  return x > -0.00001 && x < 0.00001;
}

/*
 * Return true iff x is non-positive.
 */
boolean isNonPositive(float x) {
  return x < 0.00001;
}

/*
 * Return clamped value of x w.r.t. [a, b].
 */
float clamp(float x, float a, float b) {
  return min(max(x, a), b);
}

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
  // assert a * a + b * b > 0.0000001;  // at least one of {a, b} is non-zero
  if (a * a + b * b < 0.0000001) {
    println("In a*cos(theta) + b*sin(theta) = c, a and b might be too small!");
  }
  float[] thetas = new float[2];
  if (notZero(a) && notZero(b)) {  // a != 0 and b != 0
    // println("a = ", a, "b = ", b);
    /* How to decide (the sign of) tmp0 here? According to the sign of b. */
    float tmp0 = asin(c / (sqrt(a * a + b * b)));  // [-pi/2, pi/2]
    if (b < 0) tmp0 = -tmp0;  // [-pi/2, pi/2]
    float tmp1 = atan(a / b);  // (-pi/2, pi/2)
    thetas[0] = tmp0 - tmp1;  // (-pi, pi)
    thetas[1] = PI - tmp0 - tmp1;  // (0, 2pi)
    thetas[0] = thetas[0] < 0 ? thetas[0] + TWO_PI : thetas[0];  // [0, pi), (pi, 2pi)
  } else if (isZero(a)) {  // a == 0 and b != 0
    println("In solving a*cos(theta) + b*sin(theta) = c: a == 0");
    thetas[0] = asin(c/b);  // [-pi/2, pi/2]
    thetas[1] = PI - thetas[0];  // [pi/2, 3pi/2]
    if (thetas[0] < 0) thetas[0] += TWO_PI;  // [3pi/2, 0), [0, pi/2]
  } else {  // a != 0 and b == 0
    println("In solving a*cos(theta) + b*sin(theta) = c: b == 0");
    thetas[0] = acos(c/a);  // [0, pi]
    thetas[1] = TWO_PI - thetas[0];  // [pi, 2pi]
  }
  return thetas;
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
void showTriangles(ArrayList<Triangle> triangles, pt[] G) {
  int n = triangles.size();
  beginShape(TRIANGLES);
  for (int i = 0; i < n; ++i) {
    if (triangles.get(i) == null) continue;
    vertex(G[triangles.get(i).a]);
    vertex(G[triangles.get(i).b]);
    vertex(G[triangles.get(i).c]);
  }
  endShape();
}

/*
 * Show normals of triangles, one normal per triangle face starting at its center.
 */
void showTriangleNormals(ArrayList<Triangle> triangles, pt[] G) {
  int n = triangles.size();
  for (int i = 0; i < n; ++i) {
    if (triangles.get(i) == null) continue;
    pt A = G[triangles.get(i).a];
    pt B = G[triangles.get(i).b];
    pt C = G[triangles.get(i).c];
    showNormalToTriangle(A, B, C, 10, 1);
  }
}

/*
 * Show inner vertices.
 */
void showInnerVertices(ArrayList<Vertex> vertices, pt[] G) {
  int n = vertices.size();
  for (int i = 0; i < n; ++i) {
    if (vertices.get(i).isInner) {
      show(G[vertices.get(i).id], 2);
    }
  }
}

/*
 * Show manifold edges, i.e. those with exactly 2 adjacent faces.
 */
void showManifoldEdges(boolean[][] manifoldMask, pt[] G, int nv) {
  for (int i = 0; i < nv; ++i) {
    for (int j = i + 1; j < nv; ++j) {
      if (manifoldMask[i][j] && manifoldMask[j][i]) {
        collar(G[i], V(G[i], G[j]), 1, 1);
      }
    }
  }
}

/*
 * Show the normal to triangle (A, B, C) as an arrow. d controls the length of
 * the arrow and r controls the size of the arrow tip.
 */
void showNormalToTriangle(pt A, pt B, pt C, float d, float r) {
  vec N = normalToTriangle(A, B, C);
  pt D = P(A, B, C);
  arrow(D, V(d, N), r);
}

/*
 * Show the plane defined by (p, n). s controls the size of the plane shown.
 */
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
 * Show the circle defined by (c, n, r) where c is the center, n is the normal,
 * r is the radius.
 */
void showCircle(pt c, vec n, float r) {
  noFill();
  stroke(0);
  strokeWeight(3);
  float a = 0;
  float da = TWO_PI / 36;
  vec vi = constructNormal(n);
  vec vj = N(n, vi);
  beginShape();
  for (int i = 0; i < 36; ++i, a += da) {
    vertex(P(c, r * cos(a), vi, r * sin(a), vj));
  }
  endShape(CLOSE);
}

/*
 * Show the circumcircle defined by triangle (pa, pb, pc).
 */
void showCircumcircleOfTriangle(pt pa, pt pb, pt pc, pt center, vec normal, Float radius) {
  if (center == null) center = circumcenterOfTriangle(pa, pb, pc);
  if (normal == null) normal = normalToTriangle(pa, pb, pc);
  if (radius == null) radius = d(center, pa);
  showCircle(center, normal, radius);
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
    drawArc(ps[i], c, ps[i+1], r, w1);
  }

  fill(blue);
  for (int i = 1; i < nv; i += 2) {
    show(ps[i], w2);
  }
}

/* Misc functions below. */

/*
 * Exception handler.
 */
void exceptionHandler() {
  P.savePts("data/pts_unnamed");
  rs.save("data/rs_unnamed");
  exitDraw = true;
}

/*
 * Convert a 2D point array into a 1D point array. This is used for ring set
 * processing.
 */
pt[] convertTo1DArray(pt[][] points, int nr, int nc) {
  pt[] G = new pt[nr * nc];
  int k = 0;
  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < nc; ++j) {
      G[k++] = points[i][j];
    }
  }
  return G;
}