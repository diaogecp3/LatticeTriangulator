/******************************************************************************
 * Test functions.
 ******************************************************************************/


/*
 * TODO:
 * 1). For multiple-tests functions, it may be better to reuse the random inputs
 *     on different methods.
 */


int numTests = 10;
int numBackup3RT = 0;
int numPointsPerTest = 64;
boolean saveFailures = false;

boolean showFirstCone = false;
boolean showSecondCone = false;
boolean showCoarseCorridor = false;

boolean showSpheres = false;

boolean showStereoProjection = false;

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
  void show(color c0, color c1) {
    fill(c0);
    beginShape(TRIANGLES);
    vertex(points[a]);
    vertex(points[b]);
    vertex(points[c]);
    endShape();
    fill(c1);
    showBall(points[d], 3);
  }
}

/*
 * Return true if the given triangle mesh passes convexity test, i.e. if the
 * given triangle mesh is a convex hull.
 */
boolean passConvexityTest(ArrayList<Triangle> triangles, pt[] points, int nv) {
  for (Triangle t : triangles) {
    if (t == null) continue;
    pt pa = points[t.a];
    pt pb = points[t.b];
    pt pc = points[t.c];
    vec normal = U(N(pa, pb, pc));
    for (int j = 0; j < nv; ++j) {
      if (j == t.a || j == t.b || j == t.c) continue;
      pt pd = points[j];
      vec AD = U(pa, pd);
      if (dot(normal, AD) > gEpsilon) {
        // DebugInfoConvexity dInfo = new DebugInfoConvexity(points, t.a, t.b, t.c, j);
        // dInfo.show(green, orange);
        return false;
      }
    }
  }
  return true;
}

boolean passConvexityTest(ArrayList<Triangle> triangles, ArrayList<pt> points) {
  int nv = points.size();
  for (Triangle t : triangles) {
    if (t == null) continue;
    pt pa = points.get(t.a);
    pt pb = points.get(t.b);
    pt pc = points.get(t.c);
    vec normal = U(N(pa, pb, pc));
    for (int j = 0; j < nv; ++j) {
      if (j == t.a || j == t.b || j == t.c) continue;
      pt pd = points.get(j);
      vec AD = U(pa, pd);
      if (dot(normal, AD) > gEpsilon) return false;
    }
  }
  return true;
}

/*
 * Return true if the given triangle mesh passes quality test.
 *
 * For the time being, the quality test only contains convexity test.
 */
boolean passQualityTest(ArrayList<Triangle> triangles, pt[] G, int nv) {
  if (!passConvexityTest(triangles, G, nv)) return false;
  return true;
}


/* Multiple-tests functions below. */


/*
 * Test the performance of convex hull generation for points on a sphere. The
 * number of tests is n. For each test, input size is m.
 */
void convexHullTests(int n, int m) {
  float[] times = new float[n];  // in millisecons
  int p = n / 10;
  for (int i = 0; i < n; ++i) {
    generatePointsOnSphere(gPoints, gSphereCenter, gSphereRadius, m);
    long startTime = System.nanoTime();
    ArrayList<Triangle> triangles = generateConvexHull(gPoints.G, gPoints.nv);
    long endTime = System.nanoTime();
    times[i] = (endTime - startTime) / 1000000.0;
    if (i % p == 0) System.out.format("duration of %d-th test = %f ms. \n", i, times[i]);
  }

  float avg = average(times, n, 1, n);
  System.out.format("Generate a convex hull for %d points (%d tests)," +
                    " average time (ignore the first one) = %f ms.\n", m, n, avg);
}


/*
 * Test the performance of ring set triangulation. n is the number of tests. For
 * each test, there are nr rings with np points on each ring. attenuation controls
 * the size of each ring.
 */
void ringSetTriangulationTests(int n, int nr, int np, float attenuation) {
  assert nr >= 3;
  float[] times = new float[n];
  int method = 0;
  for (int i = 0; i < n; ++i) {
    RingSet ringSet = new RingSet(gSphereCenter, gSphereRadius, nr, np);
    ringSet.init();
    ringSet.generatePoints(attenuation);
    long st = System.nanoTime();
    ringSet.generateTriangleMesh(method);
    long ed = System.nanoTime();
    times[i] = (ed - st) / 1000000.0;
  }
  float avgTime = average(times, n, 0, n);
  switch (method) {
    case 1:
      System.out.format("fast triangle mesh generation using 3RT and 2RT trick\n");
      break;
    default:
      System.out.format("convex hull generation\n");
  }
  System.out.format("Triangle mesh generation for ring set " +
                    "(n = %d, nr = %d, np = %d, attenuation = %f): " +
                    "average time = %f ms.\n", n, nr, np, attenuation, avgTime);
}


/*
 * Test the performance of three-ring-triangle generation. n is the number of
 * tests. np is the number of points sampled on each circle. attenuation controls
 * the size of each circle.
 */
void threeRingTriangleTests(int n, int np, float attenuation) {
  boolean[] successes = new boolean[n];
  boolean[] valids = new boolean[n];
  float[] times = new float[n];
  for (int i = 0; i < n; ++i) {
    RingSet ringSet = new RingSet(gSphereCenter, gSphereRadius, 3, np);
    ringSet.init();
    ringSet.generatePoints(attenuation);
    ringSet.refConvexHull = new TriangleMesh(ringSet.centers);
    ringSet.refConvexHull.addTriangle(new Triangle(0, 1, 2));
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
      if (valids[i]) ringSet.save("data/tmp/rs_3rt_low_quality_" + i);
      else ringSet.save("data/tmp/rs_3rt_null_" + i);
    }
  }
  float avgTime = average(times, n, 0, n);
  float succRate = accuracy(successes, n, 0, n, null);
  float succRateValid = accuracy(successes, n, 0, n, valids);
  switch (method3RT) {
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
      System.out.format("Approximated Supporting Plane method:\n");
      break;
    default:
      System.out.format("Naive method: \n");
  }
  System.out.format("Three-ring triangle generation " +
                    "(n = %d, np = %d, attenuation = %f): " +
                    "success rate = %f, success rate (valid) = %f, " +
                    "average time = %f ms, times to use backup method = %d.\n",
                    n, np, attenuation, succRate, succRateValid, avgTime, numBackup3RT);
}

/*
 * Test the performance of supporting-plane-of-three-circles. n is the number of
 * tests. attenuation controls the size of each circle.
 */
void supPlaneThreeCirclesTests(int n, float attenuation) {
  float[] times = new float[n];
  for (int i = 0; i < n; ++i) {
    RingSet ringSet = new RingSet(gSphereCenter, gSphereRadius, 3, 3);
    ringSet.init();
    ringSet.generatePoints(attenuation);
    long st = System.nanoTime();
    ringSet.oneSupPlaneThreeCircles(0, 1, 2);  // closed-form solution
    long ed = System.nanoTime();
    times[i] = (ed - st) / 1000000.0;
  }
  float avgTime = average(times, n, 0, n);
  System.out.format("Extreme plane generation (n = %d, attenuation = %f): " +
                    "average time = %f ms.\n", n, attenuation, avgTime);
}

void exactCHAllCirclesTests(int n, int nRings) {
  float[] times = new float[n];
  int f = 2 * nRings - 4;
  int m = 0;  // number of examples that don't satisfy F = 2V - 4
  for (int i = 0; i < n; ++i) {
    RingSet ringSet = new RingSet(gSphereCenter, gSphereRadius, nRings, 3);
    ringSet.init();
    ringSet.generatePoints(1.0);
    long st = System.nanoTime();
    ringSet.generateExTrisNaive();
    long ed = System.nanoTime();
    times[i] = (ed - st) / 1000000.0;
    if (int(ringSet.exTriPoints.size() / 3) != f) {
      m++;
      ringSet.save("data/tmp/rs_wrong_number_ex_tris_"+str(i));
    }
  }
  float avgTime = average(times, n, 0, n);
  System.out.format("number of tests = %d, number of circles = %d, average time = %f, " +
                    "number of examples that don't satisfy Euler's formula = %d\n",
                    n, nRings, avgTime, m);
}



/*
 * n = number of circles
 * m = m tests
 */
void incrementalConvexHullTests(int n, int m) {
  String[] lines = new String[m + 1];
  String file = "results/stats/IncCH_" + str(n);
  lines[0] = str(m);
  float[] times = new float[m];
  float avg = 0.0;
  int k = 0;
  for (int i = 0; i < m; ++i) {
    if (i % 10 == 0) println("i = ", i);
    RingSet ringSet = new RingSet(gSphereCenter, gSphereRadius * 10, n, 3);
    ringSet.init();
    ringSet.generatePoints(1.0);
    boolean[] success = new boolean[1];
    success[0] = true;
    long st = System.nanoTime();
    ringSet.generateExactCHIncremental(success);
    long ed = System.nanoTime();
    if (success[0]) {
      times[i] = (ed - st) / 1000000.0;
      avg += times[i];
      k += 1;
      if (k == 100) break;
    } else {
      println("fail!");
      times[i] = -1.0;
      ringSet.save("data/tmp/rs_Inc_CH_wrong_" + str(n) + "_" + str(i));
    }
    lines[i + 1] = str(times[i]);
  }

  avg /= k;
  lines[0] = str(k);
  System.out.format("Test Incremental Convex Hull: n = %d, m = %d, k = %d, average time = %f ms.\n",
                     n, m, k, avg);
  saveStrings(file, lines);
}


void incrementalConvexHullTests() {
  int m = 120;
  for (int n = 3; n < 100; n *= 2) {
    println("n = ", n);
    incrementalConvexHullTests(n, m);
  }
}




/* Single-test functions below. */


void convexHullTest() {
  if (showPointSet) {
    fill(green);
    gPoints.drawBalls(2);
  }
  if (generateCH) {
    long startTime = System.nanoTime();
    ArrayList<Triangle> triangles = generateConvexHull(gPoints.G, gPoints.nv);
    long endTime = System.nanoTime();
    timeTM = (endTime - startTime) / 1000000.0;
    fill(red);
    stroke(0);
    gNumTriangles = triangles.size();
    showTriangles(triangles, gPoints.G);
  } else {
    gNumTriangles = -1;
  }
}

void convexHullWithHolesTest() {
  gRingSet.generatePoints(gAttenuation);
  assert gRingSet.points != null;
  if (showRingSet) gRingSet.show();
  if (generateCH) {
    pt[] pointArray = gRingSet.get1DPointArray();
    if (regenerateCH) {  // regenerate a convex hull
      long startTime = System.nanoTime();
      gRingSet.generateTriangleMesh(0);
      long endTime = System.nanoTime();
      timeTM = (endTime - startTime) / 1000000.0;
    }
    fill(red);
    stroke(0);
    showTriangles(gRingSet.triangles, pointArray);
    gNumTriangles = gRingSet.triangles.size();
  } else {
    gNumTriangles = -1;
  }
}

void ringSetTriangulationTest() {
  gRingSet.generatePoints(gAttenuation);
  if (showRingSet) gRingSet.show();
  if (generateCH) {
    pt[] pointArray = gRingSet.get1DPointArray();
    if (regenerateCH) {
      long startTime = System.nanoTime();
      gRingSet.generateTriangleMesh(methodTM);
      long endTime = System.nanoTime();
      timeTM = (endTime - startTime) / 1000000.0;
    }
    gNumTriangles = 0;

    if (methodTM == 1) {
      if (show3RT) {
        assert gRingSet.threeRingTriangles != null;
        fill(blue);
        stroke(0);
        showTriangles(gRingSet.threeRingTriangles, pointArray);
        noStroke();
        fill(cyan, 100);
        showTriangleNormals(gRingSet.threeRingTriangles, pointArray);
        gNumTriangles += gRingSet.threeRingTriangles.size();
      }
      if (show2RT) {
        assert gRingSet.twoRingTriangles != null;
        fill(green);
        stroke(0);
        showTriangles(gRingSet.twoRingTriangles, pointArray);
        noStroke();
        fill(lime, 100);
        showTriangleNormals(gRingSet.twoRingTriangles, pointArray);
        gNumTriangles += gRingSet.twoRingTriangles.size();
      }
    } else {
      fill(magenta);
      stroke(0);
      showTriangles(gRingSet.triangles, pointArray);
      noStroke();
      fill(cyan, 100);
      showTriangleNormals(gRingSet.triangles, pointArray);
      gNumTriangles += gRingSet.triangles.size();
    }

    if (debug3RT && method3RT == 1) {
      gRingSet.showDebug3RTInfo();
    } else if (debug2RT) {
      gRingSet.showDebug2RTInfo();
    }
  } else {
    gNumTriangles = -1;
  }
}

void subdivisionTest() {
  debugCH = false;
  gRingSet.generatePoints(gAttenuation);
  if (showRingSet) gRingSet.show();
  if (generateCH) {
    ArrayList<pt> positions = gRingSet.get1DPointArrayList();
    long startTime = System.nanoTime();
    ArrayList<Triangle> triangles = generateConvexHull(gRingSet.get2DPointArray(),
                                                       gRingSet.getNumRings(),
                                                       gRingSet.getNumPointsPerRing());
    long endTime = System.nanoTime();
    timeTM = (endTime - startTime) / 1000000.0;

    startTime = System.nanoTime();
    TriangleMesh triMesh = new TriangleMesh(positions, triangles);
    triMesh.subdivide(subdivisionTimes);

    if (projectOnSphere) {
      triMesh.projectOnSphere(gSphereCenter, gSphereRadius);
    }

    endTime = System.nanoTime();
    timeSD = (endTime - startTime) / 1000000.0;

    triMesh.show(red, true);
    // triMesh.showCornerPairs(blue, 3);
    // triMesh.showVertices(green, 1);
    gNumTriangles = triMesh.nt;
  } else {
    gNumTriangles = -1;
  }
}

void pivotPlaneAroundLineHitCircleTest() {
  assert gRingSet.nRings >= 3;
  gRingSet.generatePoints(gAttenuation);
  assert gRingSet.points != null && gRingSet.centers != null && gRingSet.normals != null && gRingSet.radii != null;
  pt c = gRingSet.centers[0];
  float r = gRingSet.radii[0];
  vec n = gRingSet.normals[0];
  pt a = gRingSet.points[1][0];
  pt b = gRingSet.points[2][0];
  pt[] contacts = pivotPlaneAroundLineHitCircle(c, r, n, a, b, gRingSet.xAxes[0], gRingSet.yAxes[0]);
  assert contacts != null && contacts.length == 2;

  fill(violet, 200);
  showTriangle(a, b, contacts[0]);
  // fill(violet, 100);
  // showTriangle(a, b, contacts[1]);
  fill(red, 200);
  disk(c, n, r);
  fill(green, 200);
  show(a, 3);
  fill(blue, 200);
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

  // debug/visualize construction
  {
    vec t0 = U(N(n, V(c, contacts[0])));
    vec t1 = U(N(n, V(c, contacts[1])));
    vec n0 = U(N(a, b, contacts[0]));
    vec n1 = U(N(a, b, contacts[1]));
    // println("dot(n0, t0) =", dot(n0, t0), "(shoule be close to 0)");
    // println("dot(n1, t1) =", dot(n1, t1), "(shoule be close to 0)");

    // show the two planes that pass through (contacts[0], contacts[0] + t0)
    fill(pink, 200);
    showPlane(contacts[0], U(c, contacts[0]), r+10);
    fill(yellow, 200);
    showPlane(contacts[0], n, r+10);

    fill(cyan, 200);
    arrow(c, V(20, n), 2);
    fill(magenta, 200);
    arrow(c, V(c, contacts[0]), 2);
  }
}

void exactCHEdgeCircleTest() {
  // assert rs.nRings >= 3;
  // rs.generatePoints(attenuation);
  // assert rs.points != null && rs.centers != null && rs.normals != null &&
  //        rs.radii != null && rs.sameRadius == false;
  // exactCHEdgeCircle(rs.points[1][0], rs.points[2][0], rs.centers[0], rs.radii[0],
  //                   rs.normals[0], rs.xAxes[0], rs.yAxes[0]);

  {
    assert gEdgeCircle != null;
    exactCHEdgeCircle(gEdgeCircle.a, gEdgeCircle.b, gEdgeCircle.c, gEdgeCircle.r, gEdgeCircle.n, gEdgeCircle.vi, gEdgeCircle.vj);
  }
}

void exactCHTwoCirclesTest() {
  assert gRingSet.nRings >= 2;
  gRingSet.generatePoints(gAttenuation);
  assert gRingSet.points != null && gRingSet.centers != null && gRingSet.normals != null &&
         gRingSet.radii != null && gRingSet.sameRadius == false;
  exactCHTwoCircles(gRingSet.centers[0], gRingSet.radii[0], gRingSet.normals[0], gRingSet.xAxes[0], gRingSet.yAxes[0],
                    gRingSet.centers[1], gRingSet.radii[1], gRingSet.normals[1], gRingSet.xAxes[1], gRingSet.yAxes[1]);
}

void interactiveExactCHTwoCirclesTest() {
  if (gPoints.nv < 4) {
    println("Should use at least 4 points.");
    return;
  }
  gRingSet = new RingSet(gSphereCenter, gSphereRadius, gPoints.G, 4);
  exactCHTwoCircles(gRingSet.centers[0], gRingSet.radii[0], gRingSet.normals[0], gRingSet.xAxes[0], gRingSet.yAxes[0],
                    gRingSet.centers[1], gRingSet.radii[1], gRingSet.normals[1], gRingSet.xAxes[1], gRingSet.yAxes[1]);

}

void supPlaneThreeCirclesIterTest() {
  assert gRingSet.nRings >= 3;
  gRingSet.generatePoints(gAttenuation);
  assert gRingSet.points != null && gRingSet.centers != null && gRingSet.normals != null &&
         gRingSet.radii != null && gRingSet.sameRadius == false;

  {
    fill(red);
    disk(gRingSet.centers[0], gRingSet.normals[0], gRingSet.radii[0]);
    fill(green);
    disk(gRingSet.centers[1], gRingSet.normals[1], gRingSet.radii[1]);
    fill(blue);
    disk(gRingSet.centers[2], gRingSet.normals[2], gRingSet.radii[2]);
  }

  pt[] ps = supPlaneThreeCirclesIter(gRingSet.centers[0], gRingSet.radii[0], gRingSet.normals[0], gRingSet.xAxes[0], gRingSet.yAxes[0],
                                     gRingSet.centers[1], gRingSet.radii[1], gRingSet.normals[1], gRingSet.xAxes[1], gRingSet.yAxes[1],
                                     gRingSet.centers[2], gRingSet.radii[2], gRingSet.normals[2], gRingSet.xAxes[2], gRingSet.yAxes[2],
                                     null);
  {
    fill(orange, 180);
    showTriangle(ps[0], ps[1], ps[2]);
    fill(cyan, 100);
    showTriangleNormal(ps[0], ps[1], ps[2], 20, 4);
  }
  ps = supPlaneThreeCirclesIter(gRingSet.centers[0], gRingSet.radii[0], gRingSet.normals[0], gRingSet.xAxes[0], gRingSet.yAxes[0],
                                gRingSet.centers[2], gRingSet.radii[2], gRingSet.normals[2], gRingSet.xAxes[2], gRingSet.yAxes[2],
                                gRingSet.centers[1], gRingSet.radii[1], gRingSet.normals[1], gRingSet.xAxes[1], gRingSet.yAxes[1],
                                null);
  {
    fill(violet, 180);
    showTriangle(ps[0], ps[1], ps[2]);
    fill(magenta, 100);
    showTriangleNormal(ps[0], ps[1], ps[2], 20, 4);
  }

  if (showRingSet) gRingSet.show();
}

void exactCHThreeCirclesTest() {
  assert gRingSet.nRings >= 3;
  gRingSet.generatePoints(gAttenuation);
  assert gRingSet.points != null && gRingSet.centers != null && gRingSet.normals != null &&
         gRingSet.radii != null && gRingSet.sameRadius == false;
  exactCHThreeCircles(gRingSet.centers[0], gRingSet.radii[0], gRingSet.normals[0], gRingSet.xAxes[0], gRingSet.yAxes[0],
                      gRingSet.centers[1], gRingSet.radii[1], gRingSet.normals[1], gRingSet.xAxes[1], gRingSet.yAxes[1],
                      gRingSet.centers[2], gRingSet.radii[2], gRingSet.normals[2], gRingSet.xAxes[2], gRingSet.yAxes[2]);
}

void threeRingTriangleTest() {
  gRingSet.generatePoints(gAttenuation);
  if (showRingSet) gRingSet.show();
  if (generateCH) {
    pt[] pointArray = gRingSet.get1DPointArray();
    if (regenerateCH) {
      long startTime = System.nanoTime();
      gRingSet.generateThreeRingTriangles();
      long endTime = System.nanoTime();
      timeTM = (endTime - startTime) / 1000000.0;
    }
    gNumTriangles = 0;
    assert gRingSet.threeRingTriangles != null;
    fill(red);
    stroke(0);
    strokeWeight(2);
    showTriangles(gRingSet.threeRingTriangles, pointArray);
    fill(#0AED62, 100);  // light green
    noStroke();
    showTriangleNormals(gRingSet.threeRingTriangles, pointArray);
    gNumTriangles += gRingSet.threeRingTriangles.size();
  } else {
    gNumTriangles = -1;
  }
}

void interactiveNaiveCHTest() {
  if (gPoints.nv < 6) {
    println("Should use at least 6 points.");
    return;
  }

  int nv = gPoints.nv - gPoints.nv % 2;
  gRingSet = new RingSet(gSphereCenter, gSphereRadius, gPoints.G, nv);

  if (!gRingSet.isValid()) {
    validRS = false;
    return;
  } else {
    validRS = true;
  }

  if (showCircleSet) {
    gRingSet.showCircles(null);
  }

  if (showDiskSet) {
    gRingSet.showDisks(null);
  }

  gRingSet.generateExTrisNaive();
  if (show3RT) {
    fill(blue, 200);
    gRingSet.showExTris();
  }

  if (debugST) {  // show circumcircles, normals of supporting triangles
    ArrayList<pt> centers = gRingSet.debugSTInfo.circumcenters;
    ArrayList<Float> radii = gRingSet.debugSTInfo.circumradii;
    ArrayList<vec> normals = gRingSet.debugSTInfo.normals;
    for (int i = 0; i < centers.size(); ++i) {
      if (i % 2 == 0) fill(green, 100);
      else fill(blue, 100);
      disk(centers.get(i), normals.get(i), radii.get(i));
    }
    fill(cyan, 100);
    for (int i = 0; i < centers.size(); ++i) {
      arrow(centers.get(i), V(30, normals.get(i)), 4);
    }

    {  // show the plane defined by centers[0, 1, 2]
      pt c = P(gRingSet.centers[0], gRingSet.centers[1], gRingSet.centers[2]);
      vec n = U(N(gRingSet.centers[0], gRingSet.centers[1], gRingSet.centers[2]));
      fill(magenta, 100);
      showPlane(c, n, 70);
    }
  }

  float delta = TWO_PI / 30;
  gRingSet.generateExEdges(delta);
  if (show2RT) {
    fill(green, 200);
    stroke(0);
    gRingSet.showExEdges();
  }
}

void interactiveIncCHTest() {
  if (gPoints.nv < 4) {
    println("Should use at least 4 points.");
    return;
  }

  int nv = gPoints.nv - gPoints.nv % 2;
  gRingSet = new RingSet(gSphereCenter, gSphereRadius, gPoints.G, nv);

  if (!gRingSet.isValid()) {
    validRS = false;
    return;
  } else {
    validRS = true;
  }

  if (showCircleSet) {
    if (debugIncCH && nv > 4) gRingSet.showCircles(debugIncCHIter);
    else gRingSet.showCircles(null);
  }

  if (showDiskSet) {
    fill(red);
    if (debugIncCH && nv > 4) gRingSet.showDisks(debugIncCHIter);
    else gRingSet.showDisks(null);
  }

  gRingSet.generateExactCHIncremental(null);

  if (showTriangleFaces) {
    fill(blue);
    gRingSet.showIncTriangles();
  }

  if (showCorridorFaces) {
    fill(green);
    gRingSet.showIncCorridors();
  }

  if (showPolygons) {
    gRingSet.generateConvexTriMesh();
    fill(red);
    gRingSet.showPolygons();
  }

  if (debugIncCH) {
    gRingSet.showDebugIncCHInfo();
  }

  if (showApolloniusDiagram) {
    gRingSet.showApolloniusDiagram();
  }

  // if (debugApolloniusDiagram) {
  //   gRingSet.showADDebugInfo();
  // }

  // if (true) {
    // boolean b0 = gRingSet.incCorridors.get(0).isUniformlySampled();
    // gRingSet.incCorridors.get(0).testCorrespondence();
  // }

  /* Test stereo projection. */
  if (showStereoProjection && gRingSet.nRings >= 3) {
    // pt northPole = P(gRingSet.sphereCenter, gRingSet.sphereRadius, V(0, 0, 1));
    // StereoProjector sp = new StereoProjector(gRingSet.sphereCenter, gRingSet.sphereRadius, StereoType.NORTH_POLE_SOUTH_PLANE);
    StereoProjector sp = new StereoProjector(gRingSet.sphereCenter, gRingSet.sphereRadius, StereoType.CUSTOM_NORTH_POLE_SOUTH_PLANE);
    int maxCircleID = gRingSet.getMaxCircleID();
    pt northPole = gRingSet.contacts[maxCircleID];
    if (sp.sType == StereoType.CUSTOM_NORTH_POLE_SOUTH_PLANE) {
      sp.setNorthPole(northPole);
    }

    /*
    Circle cir0 = gRingSet.getCircle(0);
    Circle cir1 = gRingSet.getCircle(1);
    Circle cir2 = gRingSet.getCircle(2);

    Circle pCir0 = sp.project(cir0);
    Circle pCir1 = sp.project(cir1);
    Circle pCir2 = sp.project(cir2);

    RingSet.IncTriangle t0 = gRingSet.incTriangles.get(0);
    RingSet.IncTriangle t1 = gRingSet.incTriangles.get(1);
    Circle circumcir0 = t0.getCircumcircle();
    Circle circumcir1 = t1.getCircumcircle();

    Circle pCircumcir0 = sp.project(circumcir0);
    Circle pCircumcir1 = sp.project(circumcir1);

    {
      strokeWeight(3);
      stroke(red);
      cir0.show();
      pCir0.show();
      showSegment(northPole, pCir0.c);

      stroke(green);
      cir1.show();
      pCir1.show();
      showSegment(northPole, pCir1.c);

      stroke(blue);
      cir2.show();
      pCir2.show();
      showSegment(northPole, pCir2.c);

      stroke(magenta);
      circumcir0.show();
      pCircumcir0.show();
      showSegment(northPole, pCircumcir0.c);

      stroke(cyan);
      circumcir1.show();
      pCircumcir1.show();
      showSegment(northPole, pCircumcir1.c);
      noStroke();
      strokeWeight(1);
    }
    */

    Circle[] pCirs = new Circle[gRingSet.nRings];
    for (int i = 0; i < gRingSet.nRings; ++i) {
      Circle cir = gRingSet.getCircle(i);
      pCirs[i] = sp.project(cir);
    }
    Circle[] pCircumcirs = new Circle[gRingSet.incTriangles.size()];
    for (int i = 0; i < gRingSet.incTriangles.size(); ++i) {
      RingSet.IncTriangle t = gRingSet.incTriangles.get(i);
      Circle cir = t.getCircumcircle();
      pCircumcirs[i] = sp.project(cir);
    }

    ArrayList<ArrayList<pt>> apolloniusEdges = new ArrayList();
    for (RingSet.IncCorridor cor : gRingSet.incCorridors) {
      apolloniusEdges.add(cor.allPointsAD());
    }

    ArrayList<ArrayList<pt>> pApolloniusEdges = new ArrayList();
    for (ArrayList<pt> edge : apolloniusEdges) {
      ArrayList<pt> pEdge = new ArrayList<pt>();
      for (pt p : edge) {
        pEdge.add(sp.project(p));
      }
      pApolloniusEdges.add(pEdge);
    }

    {
      hint(DISABLE_DEPTH_TEST);
      strokeWeight(3);
      stroke(red);
      for (Circle cir : pCirs) cir.show();
      stroke(blue);
      for (Circle cir : pCircumcirs) cir.show();

      strokeWeight(2);
      stroke(green);
      for (ArrayList<pt> pEdge : pApolloniusEdges) {
        showPolyLine(pEdge);
      }
      noStroke();
      strokeWeight(1);

      hint(ENABLE_DEPTH_TEST);
    }

    {  // show the south plane
      fill(orange, 100);
      // pt southPole = P(sp.center, V(0, 0, -sp.radius));
      // showPlane(southPole, V(0, 0, 1), 50 * sp.radius);
      pt southPole = P(sp.center, -1.02, V(sp.center, sp.northPole));
      showPlane(southPole, U(sp.center, sp.northPole), 50 * sp.radius);
    }
  }
}

void corridorTest() {
  if (gPoints.nv < 6) {
    println("Should use at least 6 points.");
    return;
  }

  int nv = gPoints.nv - gPoints.nv % 2;
  RingSet rs = new RingSet(gSphereCenter, gSphereRadius, gPoints.G, nv);

  if (!rs.isValid()) {
    validRS = false;
    return;
  } else {
    validRS = true;
  }

  if (showCircleSet) {
    rs.showCircles(null);
  }

  if (showDiskSet) {
    rs.showDisks(null);
  }

  rs.generateExactCHIncremental(null);

  if (showTriangleFaces) {
    fill(blue);
    rs.showIncTriangles();
  }

  if (showCorridorFaces) {
    fill(green);
    rs.showIncCorridors();
  }

  rs.showCorridor(idxIncCor);
}

void meshFromExactCHTest() {
  if (gPoints.nv < 6) {
    println("Should use at least 6 points.");
    return;
  }

  int nv = gPoints.nv - gPoints.nv % 2;
  gRingSet = new RingSet(gSphereCenter, gSphereRadius, gPoints.G, nv);

  if (!gRingSet.isValid()) {
    validRS = false;
    return;
  } else {
    validRS = true;
  }

  gRingSet.setNumPointsPerRing(gNumPointsPerRing);
  gRingSet.generatePoints(gAttenuation);
  if (showRingSet) {
    gRingSet.show();
  }

  if (showDiskSet) {
    gRingSet.showDisks(null);
  }

  gRingSet.generateExactCHIncremental(null);
  gTriangleMesh = gRingSet.generateMeshFromExactCH(0);

  if (subdivisionTimes > 0) {
    gTriangleMesh.subdivide(subdivisionTimes);
  }

  if (showTriMesh) {
    if (gTriangleMesh != null) {
      gTriangleMesh.show(khaki, true);
    }
  } else {
    if (show3RT) {
      fill(navy, 230);
      gRingSet.showThreeRingTriangles();
    }
    if (show2RT) {
      fill(springGreen, 230);
      gRingSet.showTwoRingTriangles();
    }
  }

  if (showBeams) {
    fill(orange, 230);
    gRingSet.showBeams(40.0);
  }
}

void interactiveHubTest() {
  if (gPoints.nv < 4) {
    println("Should use at least 4 points.");
    return;
  }

  int nv = gPoints.nv - gPoints.nv % 2;
  gHub = new Hub(gSphereCenter, gRadiusInnerBall, gPoints.G, nv);
  if (!gHub.valid) {
    println("Hub is not valid!");
    return;
  }

  gHub.generateIntersectionCircles();
  gHub.liftCones(gGapWidth);
  if (showIntersectionCircles) {
    gHub.showIntersectionCircles();
  }

  gRingSet = gHub.circlesToRingSet(gNumPointsPerRing);
  gHub.generateBeamSamples(gNumPointsPerRing, gRingSet.xAxes);

  gBeamMesh = gHub.generateBeamMesh();
  if (showLiftedCones) {
    // gBeamMesh.show(cyan, showTriangleStrokes);
    gHub.showLiftedCones(cyan, 255);
  }

  if (showDiskSet) {
    gRingSet.showDisks(null);
  }

  gRingSet.generateExactCHIncremental(null);

  if (showTriangleFaces) {
    fill(blue);
    // fill(purple);
    gRingSet.showIncTriangles();
  }

  if (showCorridorFaces) {
    fill(green);
    // fill(purple);
    gRingSet.showIncCorridors();
  }

  gTriangleMesh = gRingSet.generateConvexTriMesh();

  /* Generate gap mesh. */
  // gHub.initGaps(gRingSet.borders);
  // gGapMesh = gHub.generateGapMesh();
  // if (showGapMesh) {
  //   gGapMesh.show(navy, showTriangleStrokes);
  // }
  // gTriangleMesh.augment(gGapMesh);  // merge convex hull mesh and gap mesh

  // gTriangleMesh.subdivide(subdivisionTimes);
  // if (projectOnHub) {
  //   gTriangleMesh.projectOnHub(gHub, projectMethod);
  // }

  if (showTriMesh) {
    gTriangleMesh.show(purple, showTriangleStrokes);
  }

  if (showHub) {
    gHub.show(lightSalmon, 130);
  }

  if (showBoundingSphere) {
    gHub.showBoundingSphere(khaki, 100);
  }
}

/* Test the case when mix(N1, N2, N3) = 0. */
void supPlaneThreeCirclesSpecialTest() {
  if (gPoints.nv < 6) return;
  pt[] points = gPoints.G;

  /* If we fix the 3 centers: */
  // vec v1 = U(V(1, 1, 0));
  // vec v2 = U(V(0, 1, 1));
  // vec v3 = U(V(-0.5, v1, -0.8, v2));
  // points[0].set(P(centerOfSphere, radiusOfSphere, v1));
  // points[2].set(P(centerOfSphere, radiusOfSphere, v2));
  // points[4].set(P(centerOfSphere, radiusOfSphere, v3));

  /* If we don't fix the 3 centers: */
  vec v1 = U(gSphereCenter, points[0]);
  vec v2 = U(gSphereCenter, points[2]);
  vec v3 = U(V(-0.5, v1, 0.8, v2));
  points[4].set(P(gSphereCenter, gSphereRadius, v3));

  gRingSet = new RingSet(gSphereCenter, gSphereRadius, points, 6);

  if (!gRingSet.isValid()) {
    validRS = false;
    return;
  } else {
    validRS = true;
  }

  if (showDiskSet) {
    fill(red);
    gRingSet.showDisks(null);
  }

  gRingSet.generateExactCHIncremental(null);

  if (showTriangleFaces) {
    fill(blue);
    gRingSet.showIncTriangles();
  }

  if (showCorridorFaces) {
    fill(green);
    gRingSet.showIncCorridors();
  }

  return;
}

/*
 * Make some nice pictures about the supporting plane of 3 circles.
 */
void supPlaneThreeCirclesTest() {
  if (gPoints.nv < 6) return;
  gRingSet = new RingSet(gSphereCenter, gSphereRadius, gPoints.G, 6);

  if (!gRingSet.isValid()) {
    validRS = false;
    return;
  } else {
    validRS = true;
  }

  if (showDiskSet) {
    fill(red);
    gRingSet.showDisks(null);
  }
  gRingSet.generateExactCHIncremental(null);

  if (showTriangleFaces) {
    fill(blue);
    gRingSet.showIncTriangles();
  }

  if (showCones) {
    fill(purple);
    float extDist = 0;
    gRingSet.showCones(extDist);
  }

  if (showFirstCone) {
    RingSet.IncTriangle t0 = gRingSet.incTriangles.get(0);
    stroke(blue);
    strokeWeight(3);
    t0.showCircumcircle();
    strokeWeight(1);
    noStroke();
    fill(cyan, 100);
    t0.showCone(gSphereCenter, 120);
    fill(navy, 100);
    t0.showFace();
  }

  if (showSecondCone) {
    RingSet.IncTriangle t1 = gRingSet.incTriangles.get(1);
    stroke(green);
    strokeWeight(3);
    t1.showCircumcircle();
    strokeWeight(1);
    noStroke();
    fill(lime, 100);
    t1.showCone(gSphereCenter, M(t1.normal), 120);
    fill(springGreen, 100);
    t1.showFace();
  }



  if (showCircleSet) {
    hint(DISABLE_DEPTH_TEST);
    stroke(black);
    gRingSet.showCircles(null);
    noStroke();
    hint(ENABLE_DEPTH_TEST);
  }
}

/*
 * Given two circles, tangent to each other, on a sphere, test whether the two
 * circle centers, the sphere center, and the point of tangency are coplanar.
 */
void geodesicDistanceTest() {
  if (gPoints.nv < 4) return;
  pt[] points = gPoints.G;
  vec n0 = U(gSphereCenter, points[0]);  // n1 is normal to n0
  vec n1 = U(N(gSphereCenter, points[0], points[2]));
  vec n2 = N(n1, n0);

  float rr = gSphereRadius + gSphereRadius;
  float r2 = gSphereRadius * gSphereRadius;
  float d = d(points[0], points[1]);
  float r0 = d * sqrt((rr-d)*(rr+d)) / rr;
  pt c0 = null;
  if (dot(V(gSphereCenter, points[1]), n0) > 0) {
    c0 = P(gSphereCenter, sqrt(r2 - r0 * r0), n0);
  } else {
    c0 = P(gSphereCenter, -sqrt(r2 - r0 * r0), n0);
  }
  points[3].set(P(c0, r0, n2));

  gRingSet = new RingSet(gSphereCenter, gSphereRadius, points, 4);

  if (showAuxPlane) {
    fill(orange, 100);
    showPlane(gSphereCenter, points[0], points[2], 100);
  }

  hint(DISABLE_DEPTH_TEST);
  pen(cyan, 3);
  show(gSphereCenter, points[0]);
  show(gSphereCenter, points[2]);
  strokeWeight(3);
  showCircle(gSphereCenter, n1, gSphereRadius);
  strokeWeight(1);
  noStroke();
  hint(ENABLE_DEPTH_TEST);

  if (showCircleSet) {
    gRingSet.showCircles(null);
  }
}


void ellipticConeTest() {
  if (gPoints.nv < 4) return;

  gRingSet = new RingSet(gSphereCenter, gSphereRadius, gPoints.G, 4);

  if (!gRingSet.isValid()) {
    validRS = false;
    return;
  } else {
    validRS = true;
  }

  if (showCircleSet) {
    gRingSet.showCircles(null);
  }

  if (showDiskSet) {
    gRingSet.showDisks(null);
  }

  gRingSet.generateExactCHIncremental(null);

  if (showCorridorFaces) {
    fill(green);
    gRingSet.showIncCorridors();
  }

  if (debugIncCH) {
    gRingSet.showDebugIncCHInfo();
  }

  if (showCoarseCorridor) {
    float del = TWO_PI / 10;
    fill(green);
    gRingSet.incCorridors.get(0).show(del);
    gRingSet.incCorridors.get(1).show(del);
    gRingSet.incCorridors.get(0).checkCoplanarity(del, false);
    gRingSet.incCorridors.get(1).checkCoplanarity(del, false);
  }
}


/*
 * Test about converting a hub into a triangle mesh.
 */
void interactiveHubToMeshTest() {
  if (gPoints.nv < 4) {
    println("Should use at least 4 points.");
    return;
  }

  int nv = gPoints.nv - gPoints.nv % 2;
  gHub = new Hub(gSphereCenter, gRadiusInnerBall, gPoints.G, nv);
  if (!gHub.valid) {
    println("Hub is not valid!");
    return;
  }

  vec[] xAxes = new vec[gHub.nNeighbors];
  for (int i = 0; i < xAxes.length; ++i) {
    xAxes[i] = constructNormal(gHub.tCones[i].normal);
    vec y = N(gHub.tCones[i].normal, xAxes[i]);
    xAxes[i] = R(xAxes[i], PI, xAxes[i], y);
  }

  // gTriangleMesh = gHub.generateTriMesh(gNumPointsPerRing, gGapWidth, xAxes, true);

  {
    gHub.setNumEdgesRegularPolygon(gNumPointsPerRing);
    gHub.setGapDistance(gGapWidth);
    gHub.setXDirestions(xAxes);
    gHub.setIsHalved(true);
    gTriangleMesh = gHub.generateTriMesh();
  }

  if (showTriMesh) {
    gTriangleMesh.show(cyan, true);
  }

  if (showHub) {
    gHub.show(lightSalmon, 130);
  }
}

void constructEllipticConeTest() {
  if (gPoints.nv < 4) {
    println("Should use at least 4 points.");
    return;
  }

  gRingSet = new RingSet(gSphereCenter, gSphereRadius, gPoints.G, 4);

  gRingSet.showEllipticCone(0, 1);

  if (showCircleSet) gRingSet.showCircles(2);
  if (showDiskSet) gRingSet.showDisks(2);
  if (showAuxPlane) {
    fill(orange, 180);
    showPlane(gSphereCenter, gPoints.G[0], gPoints.G[2], gSphereRadius);
  }

}


/*
 * Read a hub from an augmented file and convert this hub to a mesh.
 */
void hubToMeshTest() {
  gHub = new Hub();
  // gHub.loadAugFile("data/hub_aug_unnamed");
  gHub.loadAugFile("data/outrings.cir");
  gHub.setGapDistance(gGapWidth);
  gHub.setIsHalved(true);
  gTriangleMesh = gHub.generateTriMesh();

  if (showTriMesh) gTriangleMesh.show(cyan, true);
  if (showHub) gHub.show(lightSalmon, 130);
}

void latticeTest() {
  if (showLattice) gLattice.show();

  long startTime = System.nanoTime();
  ArrayList<Integer>[] adjLists = gLattice.adjacencyLists();

  {  // debug
    // fill(orange);
    // gLattice.balls[21].show();

    // gHub = gLattice.generateHub(21, adjLists[21]);
    // gHub.generateIntersectionCircles();
    // gRingSet = gHub.circlesToRingSet(gNumPointsPerRing);
    // gRingSet.generateExactCHIncremental(null);
    // if (showDiskSet) {
    //   fill(red);
    //   gRingSet.showDisks(null);
    // }
    // if (showTriangleFaces) {
    //   fill(blue);
    //   gRingSet.showIncTriangles();
    // }
    // if (showCorridorFaces) {
    //   fill(green);
    //   gRingSet.showIncCorridors();
    //   RingSet.IncCorridor cor = gRingSet.incCorridors.get(2);
    //   cor.showFace();
    //   println(cor.samples.size());
    //   fill(yellow);
    //   for (int i = 0; i < cor.samples.size(); ++i) {
    //     showBall(cor.samples.get(i), 1);
    //     println("sample position =", cor.samples.get(i));
    //   }
    //   fill(cyan);
    //   showBall(cor.vertices[0].position, 1);
    //   showBall(cor.vertices[1].position, 1);
    //   showBall(cor.vertices[2].position, 1);
    //   showBall(cor.vertices[3].position, 1);
    // }

    // BorderedTriangleMesh btm = gHub.generateConvexHullMesh(gNumPointsPerRing);
    // btm.show(cyan, magenta, true);
  }

  gTriangleMesh = gLattice.triangulate(adjLists);
  long endTime = System.nanoTime();
  timeTM = (endTime - startTime) / 1000000.0;

  if (showTriMesh) gTriangleMesh.show(cyan, true);
}

void convexGapTest() {
  assert gGap != null && gGap.points0 != null && gGap.points1 != null;

  gFocus = gGap.center();

  {
    fill(magenta);
    gGap.show();
  }

  gTriangleMesh = gGap.toTriMesh();
  if (gTriangleMesh == null) return;

  if (gTriangleMesh.isConvex() == false) {
    println("Not convex!");
  }
  gTriangleMesh.show(cyan, true);
}

void incCHTest() {
  assert gRingSet != null;

  gFocus = gRingSet.sphereCenter;
  {
    // int i = 1, j = 2, k = 4;
    // fill(red);
    // gRingSet.showDisk(i);
    // gRingSet.showDisk(j);
    // gRingSet.showDisk(k);

    // fill(cyan);
    // pt[] sp1 = gRingSet.supPlaneThreeCirclesIterative(i, j, k);
    // showPlane(sp1[0], sp1[1], sp1[2], 40);
    // fill(magenta);
    // pt[] sp2 = gRingSet.supPlaneThreeCirclesIterative(i, k, j);
    // showPlane(sp2[0], sp2[1], sp2[2], 40);

    // fill(lime);
    // pt[] sps = gRingSet.twoSupPlanesThreeCircles(i, j, k, null, null, false);
    // showPlane(sps[0], sps[1], sps[2], 40);
    // showPlane(sps[3], sps[4], sps[5], 40);
  }

  if (showDiskSet) {
    fill(red);
    gRingSet.showDisks(null);
  }

  if (showCircleSet) {
    gRingSet.showCircles(null);
  }

  // gRingSet.generateExTrisNaive();
  // if (show3RT) {
  //   fill(blue, 200);
  //   gRingSet.showExTris();
  // }

  gRingSet.generateExactCHIncremental(null);

  if (showTriangleFaces) {
    fill(blue);
    gRingSet.showIncTriangles();
  }

  if (showCorridorFaces) {
    fill(green);
    gRingSet.showIncCorridors();
  }

  fill(yellow, 100);
  gRingSet.showSphere();
}

void hubTest() {
  assert gHub != null && gHub.nNeighbors > 0;

  gFocus = gHub.ball.c;

  {
    // gHub.intersectionDistance(0, 1);
    // println("intersection distance =", gHub.maxIntersectDist);
  }

  gHub.generateIntersectionCircles();
  if (showIntersectionCircles) {
    gHub.showIntersectionCircles();
  }

  {
    // gRingSet = gHub.circlesToRingSet(gNumPointsPerRing);
    // gRingSet.generateExactCHIncremental(null);
    // if (showTriangleFaces) {
    //   fill(blue, 200);
    //   gRingSet.showIncTriangles();
    // }
    // if (showCorridorFaces) {
    //   fill(green, 200);
    //   gRingSet.showIncCorridors();
    // }
  }

  {
    // BorderedTriangleMesh btm = gHub.generateConvexHullMesh(gNumPointsPerRing);  // gNumPointsPerRing controls the resolution of each corridor

    // if (!btm.validBorders()) {
    //   println("borders are not valid (convex).");
    // }

    // gTriangleMesh = btm.triangleMesh;
    // if (showTriMesh) {
    //   gTriangleMesh.show(cyan, true);
    // }
  }

  if (showHub) {
    gHub.show(lightSalmon, 130);
  }

  if (showBoundingSphere) {
    gHub.showBoundingSphere(khaki, 200);
  }
}

/* Some other tests for small features. */

void circlePlaneIntersectionTest() {
  assert gRingSet.nRings >= 3;
  gRingSet.generatePoints(gAttenuation);

  pt c0 = gRingSet.centers[0];
  float r0 = gRingSet.radii[0];
  vec n0 = gRingSet.normals[0];
  pt c1 = gRingSet.centers[1];
  pt c2 = gRingSet.centers[2];
  vec d = normalOfTriangle(c0, c1, c2);

  fill(orange, 100);
  showPlane(c0, d, r0);
  fill(red, 100);
  disk(c0, n0, r0);

  pt[] ps = intersectionCirclePlane(c0, r0, n0, c0, d);
  if (ps != null) {
    fill(green, 100);
    show(ps[0], 3);
    fill(blue, 100);
    show(ps[1], 3);
  }
}

void hubLineIntersectionTest() {
  if (gPoints.nv < 4) {
    println("Should use at least 4 points.");
    return;
  }

  int nv = gPoints.nv - 2 - gPoints.nv % 2;
  gHub = new Hub(gSphereCenter, gRadiusInnerBall, gPoints.G, nv);
  if (!gHub.valid) {
    println("Hub is not valid!");
    return;
  }

  if (showHub) {
    gHub.show(lightSalmon, 100);
  }

  if (gPoints.nv % 2 == 1) return;

  pt o = gPoints.G[nv];
  vec d = U(o, gPoints.G[nv + 1]);

  stroke(cyan);
  strokeWeight(3);
  showLine(o, d);
  noStroke();

  gHub.closestIntersectionWithLine(o, d);

  return;
}

void roundConeDistTest() {
  if (gPoints.nv < 3) {
    println("Should use at least 3 points.");
    return;
  }

  int nv = 2;
  gHub = new Hub(gSphereCenter, gRadiusInnerBall, gPoints.G, nv);
  if (!gHub.valid) {
    println("Hub is not valid!");
    return;
  }

  gHub.show(lightSalmon, 100);

  pt p = gPoints.G[2];
  float d = roundConeDist(p, gHub.ball.c, gHub.ball.r, gHub.neighbors[0].c, gHub.neighbors[0].r);
  println("distance =", d);

  fill(cyan, 100);
  show(p, d);
}

void intersectionTwoPlanesTest() {
  if (gPoints.nv < 6) return;

  pt[] points = gPoints.G;
  pt pa = points[0];
  vec va = normalOfTriangle(points[0], points[1], points[2]);

  pt pb = points[3];
  vec vb = normalOfTriangle(points[3], points[4], points[5]);

  fill(red);
  showPlane(pa, va, 200);
  fill(blue);
  showPlane(pb, vb, 200);

  pt[] ps = intersectionTwoPlanes(pa, va, pb, vb);
  if (ps != null) {
    pen(green, 3);
    showLine(ps[0], U(ps[0], ps[1]));
    strokeWeight(1);
    noStroke();
  }
}

void intersectionTwoSpheresTest() {
  if (gPoints.nv < 4) return;
  pt[] points = gPoints.G;
  pt c1 = points[0], c2 = points[2];
  float r1 = d(c1, points[1]), r2 = d(c2, points[3]);

  if (showSpheres) {  // show the two spheres
    fill(red, 100);
    show(c1, r1);
    fill(blue, 100);
    show(c2, r2);
  }

  // Circle cir = intersectionSphereSphere(c1, r1, c2, r2);

  Circle cir = intersectionTwoInflatedSpheres(c1, r1, c2, r2);

  if (cir == null) return;

  stroke(cyan);
  strokeWeight(5);
  cir.show();
  strokeWeight(1);
  noStroke();
}


void stereoProjectionTest() {
  if (gPoints.nv < 4) return;

  gRingSet = new RingSet(gSphereCenter, gSphereRadius, gPoints.G, 4);

  if (showStereoProjection) {
    // pt northPoleWorld = P(gRingSet.sphereCenter, gRingSet.sphereRadius, V(0, 0, 1));
    // StereoProjector sp = new StereoProjector(gRingSet.sphereCenter, gRingSet.sphereRadius, StereoType.NORTH_POLE_SOUTH_PLANE);
    pt northPole = gRingSet.contacts[0];
    StereoProjector sp = new StereoProjector(gRingSet.sphereCenter, gRingSet.sphereRadius, StereoType.CUSTOM_NORTH_POLE_SOUTH_PLANE);
    if (sp.sType == StereoType.CUSTOM_NORTH_POLE_SOUTH_PLANE) {
      sp.setNorthPole(northPole);
    }

    pt pNorthPole = sp.project(northPole);  // should be south pole
    Circle cir0 = gRingSet.getCircle(0);
    Circle cir1 = gRingSet.getCircle(1);
    Circle pCir0 = sp.project(cir0);
    Circle pCir1 = sp.project(cir1);

    strokeWeight(3);
    stroke(red);
    cir0.show();
    pCir0.show();
    stroke(blue);
    cir1.show();
    pCir1.show();
    noStroke();
    strokeWeight(1);

    fill(orange);
    // pt p = P(gRingSet.sphereCenter, gRingSet.sphereRadius, V(0, 0, -1));
    // showPlane(p, V(0, 0, 1), 10 * gRingSet.sphereRadius);
    vec southPlaneNormal = U(sp.center, sp.northPole);
    pt southPole = P(sp.center, -1.0, V(sp.center, sp.northPole));
    showPlane(southPole, southPlaneNormal, 10 * gRingSet.sphereRadius);

    fill(cyan);
    showBall(northPole, 3);
    showBall(pCir0.c, 3);
    showBall(pCir1.c, 3);

    hint(DISABLE_DEPTH_TEST);
    stroke(blue);
    showSegment(northPole, pCir0.c);
    showSegment(northPole, pCir1.c);
    noStroke();
    hint(ENABLE_DEPTH_TEST);
  }
}

void steadyLatticeTest() {
  assert gSteadyLattice != null;
  if (showSteadyLattice) gSteadyLattice.show();

  gRanges = gSteadyLattice.getCubeRange(gCubeCenter, gCubeHalfLength);

  ArrayList<Ball> balls = new ArrayList<Ball>();
  ArrayList<Edge> beams = new ArrayList<Edge>();
  gSteadyLattice.traverse(gRanges, balls, beams);

  gLattice = new Lattice(balls, beams);
  ArrayList<Integer>[] adjLists = gLattice.adjacencyLists();
  gTriangleMesh = gLattice.triangulate(adjLists);

  if (showTriMesh) gTriangleMesh.show(cyan, false);
}
