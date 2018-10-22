/******************************************************************************
 * Utility functions.
 ******************************************************************************/


boolean exitDraw = false;

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

/*
 * Exception handler.
 */
void exceptionHandler() {
  P.savePts("data/pts_unnamed");
  rs.saveRings("data/rs_unnamed");
  exitDraw = true;
}

/*
 * Return true iff |x| is not close to 0.
 */
boolean notAbsZero(float x) {
  return x <= -0.00001 || x >= 0.00001;
}

/*
 * Return true iff |x| is close to 0.
 */
boolean isAbsZero(float x) {
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

/*
 * Solve the following linear equations in two variables.
 *
 * ax + by = e
 * cx + dy = f
 */
vec2 solveLinearEquationsInTwoVars(float a, float b, float c, float d, float e, float f) {
  float det = a * d - b * c;
  if (isAbsZero(det)) return null;
  float x =(d * e - b * f) / det;
  float y = (a * f - c * e) / det;
  return new vec2(x, y);
}

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