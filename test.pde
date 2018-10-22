/*********************************************************
 * Test functions.
 *********************************************************/


int numTests = 10000;
int numBackup3RT = 0;
int numPointsPerTest = 64;
boolean showResults = false;
boolean saveFailures = false;

class DebugInfoConvexity {
  pt[] points;
  int a, b, c, d;
  DebugInfoConvexity(pt[] points, int a, int b, int c, int d) {
    this.points = points;
    this.a = a;
    this.b = b;
    this.c = c;
    this.d = d;
  }
  void display(color c0, color c1) {
    fill(c0);
    beginShape(TRIANGLES);
    vertex(points[a]);
    vertex(points[b]);
    vertex(points[c]);
    endShape();
    fill(c1);
    show(points[d], 3);
  }
}


/*
 * Test the performance of convex hull generation for points on a sphere.
 * The number of tests done is n. For each test, input size is m. If showResult
 * is true, 3D results will be shown on screen.
 */
void testConvexHull(int n, int m, boolean showResults) {
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


void testThreeRingTriangle(int n, int np, float attenuation) {
  boolean[] successes = new boolean[n];
  boolean[] valids = new boolean[n];
  float[] times = new float[n];
  for (int i = 0; i < n; ++i) {
    RingSet ringSet = new RingSet(centerOfSphere, radiusOfSphere, 3, np);
    ringSet.init();
    ringSet.generatePoints(attenuation);
    ringSet.refConvexHull = new TriangleMesh();
    ringSet.refConvexHull.triangles = new ArrayList<Triangle>();
    ringSet.refConvexHull.triangles.add(new Triangle(0, 1, 2));
    ringSet.refConvexHull.nt = 1;
    long st = System.nanoTime();
    ringSet.generateThreeRingTriangles();
    long ed = System.nanoTime();
    times[i] = (ed - st) / 1000000.0;
    pt[] points = ringSet.get1DPointArray();
    if (ringSet.threeRingTriangles.get(0) == null) {
      successes[i] = false;
      System.out.format("%d fails\n", i);
    } else {
      valids[i] = true;
      boolean success = passQualityTest(ringSet.threeRingTriangles, points, points.length);
      successes[i] = success;
    }
    if (successes[i] == false && saveFailures) {
      if (valids[i]) ringSet.saveRings("data/tmp/rs_3rt_low_quality_" + i);
      else ringSet.saveRings("data/tmp/rs_3rt_null_" + i);
    }
  }
  float avgTime = average(times, n, 0, n);
  float succRate = accuracy(successes, n, 0, n, null);
  float succRateValid = accuracy(successes, n, 0, n, valids);
  switch (method3RT) {
    case 0:
      System.out.format("Naive method:\n");
      break;
    case 1:
      System.out.format("Heuristic Search method:\n");
      break;
    case 2:
      System.out.format("Breadth First Search method:\n");
      break;
    case 3:
      System.out.format("Breadth First Search method with heuristics:\n");
      break;
    case 4:
      System.out.format("Approximated Extreme Plane method:\n");
      break;
    default:
      System.out.format("Naive method: ");
      break;
  }
  System.out.format("Three-ring triangle generation " +
                    "(n = %d, np = %d, attenuation = %f): " +
                    "success rate = %f, success rate (valid) = %f, " +
                    "average time = %f ms, times to use backup method = %d.\n",
                    n, np, attenuation, succRate, succRateValid, avgTime, numBackup3RT);
}


void testExtremePlaneThreeCircles(int n, int np, float attenuation) {
  float[] times = new float[n];
  for (int i = 0; i < n; ++i) {
    RingSet ringSet = new RingSet(centerOfSphere, radiusOfSphere, 3, np);
    ringSet.init();
    ringSet.generatePoints(attenuation);
    long st = System.nanoTime();
    ringSet.generateExtremePlaneThreeRings(0, 1, 2);
    long ed = System.nanoTime();
    times[i] = (ed - st) / 1000000.0;
  }
  float avgTime = average(times, n, 0, n);
  switch (methodEP) {
    case 1:
      System.out.format("Heuristic initialization using normal defined by 3 centers\n");
      break;
    case 2:
      System.out.format("Heuristic initialization using plane defined by 3 centers\n");
      break;
    default:
      System.out.format("Basic initialization using 3 x axes\n");
  }
  System.out.format("Extreme plane generation (n = %d, np = %d (doesn't matter), attenuation = %f): " +
                    "average time = %f ms.\n", n, np, attenuation, avgTime);
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
      if (dot(normal, AD) > 0) {
        System.out.format("dot(normal, AD) = %f, dot(U(normal), U(AD)) = %f " +
                          "norm(normal) = %f, norm(AD) = %f\n",
                          dot(normal, AD), dot(U(normal), U(AD)),
                          normal.norm(), AD.norm());
        DebugInfoConvexity dInfo = new DebugInfoConvexity(G, a, b, c, j);
        dInfo.display(green, orange);
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
      rs.generateTriangleMesh(0);
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
    if (regenerateCH) {
      long startTime = System.nanoTime();
      rs.generateTriangleMesh(methodTM);
      long endTime = System.nanoTime();
      timeCH = (endTime - startTime) / 1000000.0;
    }
    numTriangles = 0;

    if (methodTM == 1) {
      if (rs.threeRingTriangles != null) {
        fill(red);
        stroke(0);
        strokeWeight(2);
        showTriangles(rs.threeRingTriangles, pointArray);
        fill(#0AED62, 100);  // light green
        noStroke();
        showTriangleNormals(rs.threeRingTriangles, pointArray);
        numTriangles += rs.threeRingTriangles.size();
      }
      if (rs.twoRingTriangles != null) {
        fill(blue);
        stroke(0);
        showTriangles(rs.twoRingTriangles, pointArray);
        fill(#0AED62, 100);
        noStroke();
        showTriangleNormals(rs.twoRingTriangles, pointArray);
        numTriangles += rs.twoRingTriangles.size();
      }
    } else {
      fill(green);
      stroke(0);
      showTriangles(rs.triangles, pointArray);
      fill(blue, 100);
      noStroke();
      showTriangleNormals(rs.triangles, pointArray);
      numTriangles += rs.triangles.size();
    }

    if (debug3RT && method3RT == 1) {
      rs.showDebug3RTInfo();
    } else if (debug2RT) {
      rs.showDebug2RTInfo();
    } else {
      // boolean success = passQualityTest(rs.threeRingTriangles, pointArray, pointArray.length);
      // if (!success) {
      //   println("not pass quality (convexity) test!");
      // }
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
    numTriangles = triMesh.nt;
  } else {
    numTriangles = -1;
  }
}

void oneHubTest() {
  hub.showHub(red, 150);
  hub.showBoundingBall(blue, 100);
}


void onePivotPlaneAroundLineHitCircleTest() {
  assert rs.nRings >= 3;
  rs.generatePoints(attenuation);
  assert rs.points != null && rs.centers != null && rs.radii != null;
  pt c = rs.centers[0];
  float r = rs.radii[0];
  vec n = U(rs.c, c);
  pt a = rs.points[1][0];
  pt b = rs.points[2][0];
  pt[] contacts = pivotPlaneAroundLineHitCircle(c, r, n, a, b, null, null);
  assert contacts != null && contacts.length == 2;

  vec t0 = U(N(n, V(c, contacts[0])));
  vec t1 = U(N(n, V(c, contacts[1])));
  vec n0 = U(N(a, b, contacts[0]));
  vec n1 = U(N(a, b, contacts[1]));
  println("dot(n0, t0) =", dot(n0, t0), "(shoule be close to 0)");
  println("dot(n1, t1) =", dot(n1, t1), "(shoule be close to 0)");

  fill(orange, 100);
  showTriangle(a, b, contacts[0]);
  fill(purple, 100);
  showTriangle(a, b, contacts[1]);
  fill(red);
  disk(c, n, r);
  fill(blue);
  show(a, 3);
  show(b, 3);

  // show the two planes using point-normal method
  {
    // fill(green, 100);
    // vec ab = V(a, b);
    // vec cp0 = V(c, contacts[0]);
    // float alpha0 = -dot(cp0, ab) / dot(n, ab);
    // vec normal0 = A(cp0, alpha0, n);
    // showPlane(contacts[0], normal0, 20);
    // vec cp1 = V(c, contacts[1]);
    // float alpha1 = -dot(cp1, ab) / dot(n, ab);
    // vec normal1 = A(cp1, alpha1, n);
    // showPlane(contacts[1], normal1, 20);
  }
}

void tangentPlaneThreeCirclesTest() {
  assert rs.nRings >= 3;
  rs.generatePoints(attenuation);
  assert rs.points != null && rs.centers != null && rs.radii != null;

  pt c0 = rs.centers[0];
  pt c1 = rs.centers[1];
  pt c2 = rs.centers[2];
  float r0 = rs.radii[0];
  float r1 = rs.radii[1];
  float r2 = rs.radii[2];
  vec n0 = U(rs.c, c0);
  vec n1 = U(rs.c, c1);
  vec n2 = U(rs.c, c2);
  vec vi0 = rs.xAxes[0];
  vec vi1 = rs.xAxes[1];
  vec vi2 = rs.xAxes[2];
  pt[] ps = tangentPlaneThreeCircles(c0, r0, n0, vi0, null,
                                     c1, r1, n1, vi1, null,
                                     c2, r2, n2, vi2, null);
  {
    fill(red);
    //show(ps[0], 3);
    disk(c0, n0, r0);
    fill(green);
    //show(ps[1], 3);
    disk(c1, n1, r1);
    fill(blue);
    //show(ps[2], 3);
    disk(c2, n2, r2);
  }
  {
    fill(orange, 180);
    showTriangle(ps[0], ps[1], ps[2]);
    fill(cyan, 100);
    showNormalToTriangle(ps[0], ps[1], ps[2], 20, 4);
  }

  ps = tangentPlaneThreeCircles(c0, r0, n0, vi0, null,
                                c2, r2, n2, vi2, null,
                                c1, r1, n1, vi1, null);
  {
    fill(purple, 180);
    showTriangle(ps[0], ps[1], ps[2]);
    fill(magenta, 100);
    showNormalToTriangle(ps[0], ps[1], ps[2], 20, 4);
  }

  if (showRingSet) rs.showRings();
}

void exactCHThreeCirclesTest() {
  assert rs.nRings >= 3;
  rs.generatePoints(attenuation);
  assert rs.points != null && rs.centers != null && rs.radii != null;

  pt c0 = rs.centers[0];
  pt c1 = rs.centers[1];
  pt c2 = rs.centers[2];
  float r0 = rs.radii[0];
  float r1 = rs.radii[1];
  float r2 = rs.radii[2];
  vec n0 = U(rs.c, c0);
  vec n1 = U(rs.c, c1);
  vec n2 = U(rs.c, c2);
  vec vi0 = rs.xAxes[0];
  vec vi1 = rs.xAxes[1];
  vec vi2 = rs.xAxes[2];
  exactCHThreeCircles(c0, r0, n0, vi0, null,
                      c1, r1, n1, vi1, null,
                      c2, r2, n2, vi2, null);

}

void testCirclePlaneIntersection() {
  assert rs.nRings >= 3;
  pt p0 = new pt();
  pt p1 = new pt();
  rs.generatePoints(attenuation);

  pt c0 = rs.centers[0];
  float r0 = rs.radii[0];
  vec n0 = rs.normals[0];
  pt c1 = rs.centers[1];
  pt c2 = rs.centers[2];
  vec d = normalToTriangle(c0, c1, c2);

  fill(orange, 100);
  // showTriangle(c0, c1, c2);
  showPlane(c0, d, r0);
  fill(red, 100);
  disk(c0, n0, r0);

  if (intersectionCirclePlane(c0, r0, n0, c0, d, p0, p1)) {
    fill(green, 100);
    show(p0, 3);
    fill(blue, 100);
    show(p1, 3);
  }
}