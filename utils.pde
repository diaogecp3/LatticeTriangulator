/*********************************************************
 * Utility functions.
 *********************************************************/

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


/*
 * Exception handler.
 */
void exceptionHandler() {
  P.savePts("data/pts_unnamed");
  exitDraw = true;
}

boolean notAbsZero(float x) {
  return x <= -0.00001 || x >= 0.00001;
}

boolean isAbsZero(float x) {
  return x > -0.00001 && x < 0.00001;
}

boolean notZero(float x) {
  // assume that x >= 0
  return x >= 0.00001;
}

boolean isZero(float x) {
  // assume that x >= 0
  return x < 0.00001;
}

vec random3(float n) {
  return new vec(random(-n, n), random(-n, n), random(-n, n));
}

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
 * Show triangles.
 */
void showTriangles(ArrayList<Triangle> triangles, pt[] G) {
  int n = triangles.size();
  beginShape(TRIANGLES);
  for (int i = 0; i < n; ++i) {
    pt A = G[triangles.get(i).a];
    pt B = G[triangles.get(i).b];
    pt C = G[triangles.get(i).c];
    vertex(A);
    vertex(B);
    vertex(C);
  }
  endShape();
}

/*
 * Show normals of triangles, one normal per triangle face.
 */
void showTriangleNormals(ArrayList<Triangle> triangles, pt[] G) {
  int n = triangles.size();
  for (int i = 0; i < n; ++i) {
    pt A = G[triangles.get(i).a];
    pt B = G[triangles.get(i).b];
    pt C = G[triangles.get(i).c];
    showNormalToTriangle(A, B, C, 10, 1);
  }
}

/*
 * Show inner vertices
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