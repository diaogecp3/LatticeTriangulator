/*********************************************************
 * Test functions.
 *********************************************************/


int numTests = 1000;
int numPointsPerTest = 64;
boolean showResults = false;

/*
 * Test the performance of convex hull generation for points on a sphere.
 * The number of tests done is n. For each test, input size is m. If showResult
 * is true, 3D results will be shown on screen.
 */
void testCH(int n, int m, boolean showResults) {
  float[] times = new float[n];  // in millisecons
  int p = n / 10;
  for (int i = 0; i < n; ++i) {
    generatePointsOnSphere(P, centerOfSphere, radiusOfSphere, m);
    long startTime = System.nanoTime(); 
    ArrayList<Triangle> triangles = generateConvexHull(P.G, P.nv);
    long endTime = System.nanoTime();
    times[i] = (endTime - startTime) / 1000000.0;
    if (i % p == 0) System.out.format("duration of %d-th test = %f ms. \n", i, times[i]);
    if (showResults) {
      fill(red); showTriangles(triangles, P.G);
      fill(cyan); showTriangleNormals(triangles, P.G);
    }
  }
  
  float avg = average(times, n, 1, n);
  System.out.format("Generate a convex hull for %d point, %d tests in total," +
    " average time (ignore the first one) = %f ms. \n", m, n, avg);
}


void testThreeRingTriangle(int n, int np) {
  boolean[] successes = new boolean[n];
  float[] times = new float[n];
  for (int i = 0; i < n; ++i) {
    //int numPointsPerGroup = int(random(3, np+1));
    RingSet ringSet = new RingSet(centerOfSphere, radiusOfSphere, 3, np);
    ringSet.init();
    ringSet.generatePoints(1.0);
    ringSet.convexHullRef = new ArrayList<Triangle>();
    ringSet.convexHullRef.add(new Triangle(0, 1, 2));
    long st = System.nanoTime();
    ArrayList<Triangle> tris = ringSet.generateThreeRingTriangles();
    long ed = System.nanoTime();
    times[i] = (ed - st) / 1000000.0;
    pt[] points = ringSet.get1DPointArray();
    if (tris.get(0) == null) successes[i] = false;
    else {
      boolean success = passQualityTest(tris, points, points.length);
      successes[i] = success;
    }
    if (successes[i] == false) {
      ringSet.saveRings("data/rs_fail_" + i);
    }
  }
  float avgTime = average(times, n, 0, n);
  float succRate = accuracy(successes, n, 0, n);
  System.out.format("Three-ring triangle generation (np = %d, n = %d):" +
    " success rate = %f, average time = %f ms.", np, n, succRate, avgTime);
}




boolean passConvexityTest(ArrayList<Triangle> triangles, pt[] G, int nv) {
  int n = triangles.size();
  for (int i = 0; i < n; ++i) {
    if (triangles.get(i) == null) continue;
    int a = triangles.get(i).a;
    int b = triangles.get(i).b;
    int c = triangles.get(i).c;
    pt A = G[a];
    pt B = G[b];
    pt C = G[c];
    vec normal = N(A, B, C);  // non-unit normal
    for (int j = 0; j < nv; ++j) {
      if (j == a || j == b || j == c) continue;
      pt D = G[j];
      vec AD = V(A, D);
      if (notZero(dot(normal, AD))) {
        return false;
      }
    }
  }
  return true;
}


boolean passQualityTest(ArrayList<Triangle> triangles, pt[] G, int nv) {
  if (!passConvexityTest(triangles, G, nv)) return false;
  return true;
}


void testIntersectionTwoDisks() {
  pt ca = new pt(0.5, 0, 0);
  vec va = new vec(0, 0, 1);
  float ra = 1.0;
  
  pt cb = new pt(0.2, 0, 1.2);
  vec vb = new vec(1, 0, 0);
  float rb = 1.0;
  
  pt p0 = new pt(), p1 = new pt();
  intersectionTwoPlanes(ca, va, cb, vb, p0, p1);
  System.out.format("p0 = (%f, %f, %f), p1 = (%f, %f, %f) \n", p0.x, p0.y, p0.z, p1.x, p1.y, p1.z);
  
  boolean empty = emptyIntersectionTwoDisks(ca, va, ra, cb, vb, rb);
  if (empty) println("empty intersection");
  else println("non-empty intersection");
  return;
}


void oneConvexHullTest() {
  if (showPointSet) {
    fill(green);
    P.drawBalls(2);
  }
  if (generateCH) {
    long startTime = System.nanoTime();
    ArrayList<Triangle> triangles = generateConvexHull(P.G, P.nv);
    long endTime = System.nanoTime();
    timeCH = (endTime - startTime) / 1000000.0;
    fill(red);
    stroke(0);
    numTriangles = triangles.size();
    showTriangles(triangles, P.G);
  } else {
    numTriangles = -1;
  }
}

void oneConvexHullWithHolesTest() {
  rs.generatePoints(attenuation);
  if (showRingSet) rs.showRings();
  if (generateCH) {
    pt[] pointArray = rs.get1DPointArray();
    if (regenerateCH) {  // regenerate a convex hull
      long startTime = System.nanoTime();
      rs.generateTriangleMesh();
      long endTime = System.nanoTime();
      timeCH = (endTime - startTime) / 1000000.0; 
    }
    fill(red);
    stroke(0);
    showTriangles(rs.triangles, pointArray);
    numTriangles = rs.triangles.size();
  } else {
    numTriangles = -1;
  }
}

void oneFastConvexHullWithHolesTest() {
  rs.generatePoints(attenuation);
  if (showRingSet) rs.showRings();
  if (generateCH) {
    pt[] pointArray = rs.get1DPointArray();
    long startTime = System.nanoTime();
    ArrayList<Triangle> triangles = rs.generateThreeRingTriangles();
    long endTime = System.nanoTime();
    timeCH = (endTime - startTime) / 1000000.0;
    fill(red);
    stroke(0);
    showTriangles(triangles, pointArray);
    fill(green);
    noStroke();
    showTriangleNormals(triangles, pointArray);
    numTriangles = triangles.size();
    if (debugFastCH) {
      rs.showDebug3RTriInfo();
    } else {
      boolean success = passQualityTest(triangles, pointArray, pointArray.length);
      if (!success) {
        println("not pass quality (convexity) test!");
      }
    }
  } else {
    numTriangles = -1;
  }
}

void oneSubdivisionTest() {
  debugCH = false;
  rs.generatePoints(attenuation);
  if (showRingSet) rs.showRings();
  if (generateCH) {
    ArrayList<pt> positions = rs.get1DPointArrayList();
    long startTime = System.nanoTime();
    ArrayList<Triangle> triangles = generateConvexHull(rs.get2DPointArray(),
                                                       rs.getNumRings(),
                                                       rs.getNumPointsPerRing());
    long endTime = System.nanoTime();
    timeCH = (endTime - startTime) / 1000000.0;

    startTime = System.nanoTime();
    TriangleMesh triMesh = new TriangleMesh(positions, triangles);
    triMesh.subdivide(subdivisionTimes, centerOfSphere, radiusOfSphere);
    endTime = System.nanoTime();
    timeSB = (endTime - startTime) / 1000000.0;

    triMesh.showTriangleMesh(red, true);
    // triMesh.showCornerPairs(blue, 3);
    // triMesh.showVertices(green, 1);
    numTriangles = triMesh.getNumTriangles();
  } else {
    numTriangles = -1;
  }
}

void oneHubTest() {
  hub.showHub(red, 150);
  hub.showBoundingBall(blue, 100);
}