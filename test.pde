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
 * Return true if the given triangle mesh passes convexity test, i.e. if the
 * given triangle mesh is a convex hull.
 */
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
    generatePointsOnSphere(gPoints, centerOfSphere, radiusOfSphere, m);
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
    RingSet ringSet = new RingSet(centerOfSphere, radiusOfSphere, nr, np);
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
    RingSet ringSet = new RingSet(centerOfSphere, radiusOfSphere, 3, np);
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
      System.out.format("Approximated Extreme Plane method:\n");
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
    RingSet ringSet = new RingSet(centerOfSphere, radiusOfSphere, 3, 3);
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
    RingSet ringSet = new RingSet(centerOfSphere, radiusOfSphere, nRings, 3);
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
    numTriangles = triangles.size();
    showTriangles(triangles, gPoints.G);
  } else {
    numTriangles = -1;
  }
}

void convexHullWithHolesTest() {
  gRingSet.generatePoints(attenuation);
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
    numTriangles = gRingSet.triangles.size();
  } else {
    numTriangles = -1;
  }
}

void ringSetTriangulationTest() {
  gRingSet.generatePoints(attenuation);
  if (showRingSet) gRingSet.show();
  if (generateCH) {
    pt[] pointArray = gRingSet.get1DPointArray();
    if (regenerateCH) {
      long startTime = System.nanoTime();
      gRingSet.generateTriangleMesh(methodTM);
      long endTime = System.nanoTime();
      timeTM = (endTime - startTime) / 1000000.0;
    }
    numTriangles = 0;

    if (methodTM == 1) {
      if (show3RT) {
        assert gRingSet.threeRingTriangles != null;
        fill(red);
        stroke(0);
        strokeWeight(2);
        showTriangles(gRingSet.threeRingTriangles, pointArray);
        fill(#0AED62, 100);  // light green
        noStroke();
        showTriangleNormals(gRingSet.threeRingTriangles, pointArray);
        numTriangles += gRingSet.threeRingTriangles.size();
      }
      if (show2RT) {
        assert gRingSet.twoRingTriangles != null;
        fill(blue);
        stroke(0);
        showTriangles(gRingSet.twoRingTriangles, pointArray);
        fill(#0AED62, 100);
        noStroke();
        showTriangleNormals(gRingSet.twoRingTriangles, pointArray);
        numTriangles += gRingSet.twoRingTriangles.size();
      }
    } else {
      fill(green);
      stroke(0);
      showTriangles(gRingSet.triangles, pointArray);
      fill(blue, 100);
      noStroke();
      showTriangleNormals(gRingSet.triangles, pointArray);
      numTriangles += gRingSet.triangles.size();
    }

    if (debug3RT && method3RT == 1) {
      gRingSet.showDebug3RTInfo();
    } else if (debug2RT) {
      gRingSet.showDebug2RTInfo();
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

void subdivisionTest() {
  debugCH = false;
  gRingSet.generatePoints(attenuation);
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
      triMesh.projectOnSphere(centerOfSphere, radiusOfSphere);
    }

    endTime = System.nanoTime();
    timeSD = (endTime - startTime) / 1000000.0;

    triMesh.show(red, true);
    // triMesh.showCornerPairs(blue, 3);
    // triMesh.showVertices(green, 1);
    numTriangles = triMesh.nt;
  } else {
    numTriangles = -1;
  }
}

void pivotPlaneAroundLineHitCircleTest() {
  assert gRingSet.nRings >= 3;
  gRingSet.generatePoints(attenuation);
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
  gRingSet.generatePoints(attenuation);
  assert gRingSet.points != null && gRingSet.centers != null && gRingSet.normals != null &&
         gRingSet.radii != null && gRingSet.sameRadius == false;
  exactCHTwoCircles(gRingSet.centers[0], gRingSet.radii[0], gRingSet.normals[0], gRingSet.xAxes[0], gRingSet.yAxes[0],
                    gRingSet.centers[1], gRingSet.radii[1], gRingSet.normals[1], gRingSet.xAxes[1], gRingSet.yAxes[1]);
}

void supPlaneThreeCirclesIterTest() {
  assert gRingSet.nRings >= 3;
  gRingSet.generatePoints(attenuation);
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
    showNormalToTriangle(ps[0], ps[1], ps[2], 20, 4);
  }
  ps = supPlaneThreeCirclesIter(gRingSet.centers[0], gRingSet.radii[0], gRingSet.normals[0], gRingSet.xAxes[0], gRingSet.yAxes[0],
                                gRingSet.centers[2], gRingSet.radii[2], gRingSet.normals[2], gRingSet.xAxes[2], gRingSet.yAxes[2],
                                gRingSet.centers[1], gRingSet.radii[1], gRingSet.normals[1], gRingSet.xAxes[1], gRingSet.yAxes[1],
                                null);
  {
    fill(violet, 180);
    showTriangle(ps[0], ps[1], ps[2]);
    fill(magenta, 100);
    showNormalToTriangle(ps[0], ps[1], ps[2], 20, 4);
  }

  if (showRingSet) gRingSet.show();
}

void exactCHThreeCirclesTest() {
  assert gRingSet.nRings >= 3;
  gRingSet.generatePoints(attenuation);
  assert gRingSet.points != null && gRingSet.centers != null && gRingSet.normals != null &&
         gRingSet.radii != null && gRingSet.sameRadius == false;
  exactCHThreeCircles(gRingSet.centers[0], gRingSet.radii[0], gRingSet.normals[0], gRingSet.xAxes[0], gRingSet.yAxes[0],
                      gRingSet.centers[1], gRingSet.radii[1], gRingSet.normals[1], gRingSet.xAxes[1], gRingSet.yAxes[1],
                      gRingSet.centers[2], gRingSet.radii[2], gRingSet.normals[2], gRingSet.xAxes[2], gRingSet.yAxes[2]);
}

void threeRingTriangleTest() {
  gRingSet.generatePoints(attenuation);
  if (showRingSet) gRingSet.show();
  if (generateCH) {
    pt[] pointArray = gRingSet.get1DPointArray();
    if (regenerateCH) {
      long startTime = System.nanoTime();
      gRingSet.generateThreeRingTriangles();
      long endTime = System.nanoTime();
      timeTM = (endTime - startTime) / 1000000.0;
    }
    numTriangles = 0;
    assert gRingSet.threeRingTriangles != null;
    fill(red);
    stroke(0);
    strokeWeight(2);
    showTriangles(gRingSet.threeRingTriangles, pointArray);
    fill(#0AED62, 100);  // light green
    noStroke();
    showTriangleNormals(gRingSet.threeRingTriangles, pointArray);
    numTriangles += gRingSet.threeRingTriangles.size();
  } else {
    numTriangles = -1;
  }
}

void exactCHNaiveTest() {
  if (gPoints.nv < 6) {
    println("Should use at least 6 points.");
    return;
  }

  int nv = gPoints.nv - gPoints.nv % 2;
  RingSet rs = new RingSet(centerOfSphere, radiusOfSphere, gPoints.G, nv);

  if (!rs.isValid()) {
    validRS = false;
    return;
  } else {
    validRS = true;
  }

  if (showCircleSet) {
    rs.showCircles();
  }

  if (showDiskSet) {
    rs.showDisks();
  }

  rs.generateExTrisNaive();
  if (show3RT) {
    rs.showExTris();
    numTriangles = rs.exTriPoints.size() / 3;
    numRings = rs.nRings;

    fill(cyan);
    for (int i = 0; i < rs.exTriPoints.size(); i += 3) {
      showNormalToTriangle(rs.exTriPoints.get(i), rs.exTriPoints.get(i + 1), rs.exTriPoints.get(i + 2), 20, 3);
    }
  }

  if (debugST) {  // show circumcircles, normals of supporting triangles
    ArrayList<pt> centers = rs.debugSTInfo.circumcenters;
    ArrayList<Float> radii = rs.debugSTInfo.circumradii;
    ArrayList<vec> normals = rs.debugSTInfo.normals;
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
      pt c = P(rs.centers[0], rs.centers[1], rs.centers[2]);
      vec n = U(N(rs.centers[0], rs.centers[1], rs.centers[2]));
      fill(magenta, 100);
      showPlane(c, n, 70);
    }
  }

  rs.generateExEdges();
  if (show2RT) {
    rs.showExEdges();
  }

  rs.showApproxCorridors(approxMethodCorridor);
}

void exactCHIncrementalTest() {
  if (gPoints.nv < 4) {
    println("Should use at least 4 points.");
    return;
  }

  int nv = gPoints.nv - gPoints.nv % 2;
  gRingSet = new RingSet(centerOfSphere, radiusOfSphere, gPoints.G, nv);

  if (!gRingSet.isValid()) {
    validRS = false;
    return;
  } else {
    validRS = true;
  }

  if (showCircleSet) {
    gRingSet.showCircles();
  }

  if (showDiskSet) {
    gRingSet.showDisks();
  }

  gRingSet.generateExactCHIncremental();

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

  // if (showApolloniusDiagram) {
  //   gRingSet.showApolloniusDiagram();
  // }

  // if (debugApolloniusDiagram) {
  //   gRingSet.showADDebugInfo();
  // }
}

void corridorTest() {
  if (gPoints.nv < 6) {
    println("Should use at least 6 points.");
    return;
  }

  int nv = gPoints.nv - gPoints.nv % 2;
  RingSet rs = new RingSet(centerOfSphere, radiusOfSphere, gPoints.G, nv);

  if (!rs.isValid()) {
    validRS = false;
    return;
  } else {
    validRS = true;
  }

  if (showCircleSet) {
    rs.showCircles();
  }

  if (showDiskSet) {
    rs.showDisks();
  }

  rs.generateExactCHIncremental();

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
  gRingSet = new RingSet(centerOfSphere, radiusOfSphere, gPoints.G, nv);

  if (!gRingSet.isValid()) {
    validRS = false;
    return;
  } else {
    validRS = true;
  }

  gRingSet.setNumPointsPerRing(numPointsPerRing);
  gRingSet.generatePoints(attenuation);
  if (showRingSet) {
    gRingSet.show();
  }

  if (showDiskSet) {
    gRingSet.showDisks();
  }

  gRingSet.generateExactCHIncremental();
  gRingSet.generateMeshFromExactCH(0);

  gTriangleMesh = gRingSet.triMesh;  // make it global, so that we can save it by pressing 'w'

  if (subdivisionTimes > 0) {
    gRingSet.triMesh.subdivide(subdivisionTimes);
  }

  if (showTriMesh) {
    if (gRingSet.triMesh != null) {
      gRingSet.triMesh.show(khaki, true);
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
  if (gPoints.nv < 2) {
    println("Should use at least 2 points.");
    return;
  }

  int nv = gPoints.nv - gPoints.nv % 2;
  gHub = new Hub(centerOfSphere, rInnerBall, gPoints.G, nv);
  if (!gHub.valid) {
    println("Hub is not valid!");
    return;
  }

  gHub.generateIntersectionCircles();
  if (showIntersectionCircles) {
    gHub.showIntersectionCircles();
  }

  gBeamMesh = gHub.generateBeamMesh(numPointsPerRing);
  if (showLiftedCones) {
    gBeamMesh.show(cyan, showTriangleStrokes);
    // gHub.showLiftedCones(cyan, 255);
  }

  gRingSet = gHub.circlesToRingSet(numPointsPerRing);

  if (showDiskSet) {
    gRingSet.showDisks();
  }

  gRingSet.generateExactCHIncremental();

  if (showTriangleFaces) {
    fill(blue);
    gRingSet.showIncTriangles();
  }

  if (showCorridorFaces) {
    fill(green);
    gRingSet.showIncCorridors();
  }

  gTriangleMesh = gRingSet.generateConvexTriMesh();

  /* Generate gap mesh. */
  gHub.initGaps(gRingSet.borders);
  gGapMesh = gHub.generateGapMesh();
  if (showGapMesh) {
    gGapMesh.show(navy, showTriangleStrokes);
  }

  gTriangleMesh.subdivide(subdivisionTimes);
  if (projectOnHub) {
    gTriangleMesh.projectOnHub(gHub, projectMethod);
  }

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
  vec v1 = U(centerOfSphere, points[0]);
  vec v2 = U(centerOfSphere, points[2]);
  vec v3 = U(V(-0.5, v1, 0.8, v2));
  points[4].set(P(centerOfSphere, radiusOfSphere, v3));

  gRingSet = new RingSet(centerOfSphere, radiusOfSphere, points, 6);

  if (!gRingSet.isValid()) {
    validRS = false;
    return;
  } else {
    validRS = true;
  }

  if (showDiskSet) gRingSet.showDisks();

  gRingSet.generateExactCHIncremental();

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
  gRingSet = new RingSet(centerOfSphere, radiusOfSphere, gPoints.G, 6);

  if (!gRingSet.isValid()) {
    validRS = false;
    return;
  } else {
    validRS = true;
  }

  if (showDiskSet) {
    fill(red);
    gRingSet.showDisks();
  }
  gRingSet.generateExactCHIncremental();

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
    t0.showCone(centerOfSphere, 120);
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
    t1.showCone(centerOfSphere, M(t1.normal), 120);
    fill(springGreen, 100);
    t1.showFace();
  }



  if (showCircleSet) {
    hint(DISABLE_DEPTH_TEST);
    stroke(black);
    gRingSet.showCircles();
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
  vec n0 = U(centerOfSphere, points[0]);  // n1 is normal to n0
  vec n1 = U(N(centerOfSphere, points[0], points[2]));
  vec n2 = N(n1, n0);

  float rr = radiusOfSphere + radiusOfSphere;
  float r2 = radiusOfSphere * radiusOfSphere;
  float d = d(points[0], points[1]);
  float r0 = d * sqrt((rr-d)*(rr+d)) / rr;
  pt c0 = null;
  if (dot(V(centerOfSphere, points[1]), n0) > 0) {
    c0 = P(centerOfSphere, sqrt(r2 - r0 * r0), n0);
  } else {
    c0 = P(centerOfSphere, -sqrt(r2 - r0 * r0), n0);
  }
  points[3].set(P(c0, r0, n2));

  gRingSet = new RingSet(centerOfSphere, radiusOfSphere, points, 4);

  if (showAuxPlane) {
    fill(orange, 100);
    showPlane(centerOfSphere, points[0], points[2], 100);
  }

  hint(DISABLE_DEPTH_TEST);
  pen(cyan, 10);
  show(centerOfSphere, points[0]);
  show(centerOfSphere, points[2]);
  strokeWeight(3);
  showCircle(centerOfSphere, n1, radiusOfSphere);
  strokeWeight(1);
  noStroke();
  hint(ENABLE_DEPTH_TEST);

  if (showCircleSet) {
    gRingSet.showCircles();
  }
}


void ellipticConeTest() {
  if (gPoints.nv < 4) return;

  gRingSet = new RingSet(centerOfSphere, radiusOfSphere, gPoints.G, 4);

  if (!gRingSet.isValid()) {
    validRS = false;
    return;
  } else {
    validRS = true;
  }

  if (showCircleSet) {
    gRingSet.showCircles();
  }

  if (showDiskSet) {
    gRingSet.showDisks();
  }

  gRingSet.generateExactCHIncremental();

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
void hubToMeshTest() {
  if (gPoints.nv < 4) {
    println("Should use at least 4 points.");
    return;
  }

  int nv = gPoints.nv - gPoints.nv % 2;
  gHub = new Hub(centerOfSphere, rInnerBall, gPoints.G, nv);
  if (!gHub.valid) {
    println("Hub is not valid!");
    return;
  }

  gTriangleMesh = gHub.generateTriMesh(numPointsPerRing);

  if (showTriMesh) {
    gTriangleMesh.show(cyan, true);
  }

  if (showHub) {
    gHub.show(lightSalmon, 130);
  }
}


/* Some other tests for small features. */

void circlePlaneIntersectionTest() {
  assert gRingSet.nRings >= 3;
  pt p0 = new pt();
  pt p1 = new pt();
  gRingSet.generatePoints(attenuation);

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

  if (intersectionCirclePlane(c0, r0, n0, c0, d, p0, p1)) {
    fill(green, 100);
    show(p0, 3);
    fill(blue, 100);
    show(p1, 3);
  }
}

void hubLineIntersectionTest() {
  if (gPoints.nv < 4) {
    println("Should use at least 4 points.");
    return;
  }

  int nv = gPoints.nv - 2 - gPoints.nv % 2;
  gHub = new Hub(centerOfSphere, rInnerBall, gPoints.G, nv);
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
  gHub = new Hub(centerOfSphere, rInnerBall, gPoints.G, nv);
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

  pt p0 = P();
  pt p1 = P();
  boolean intersect = intersectionTwoPlanes(pa, va, pb, vb, p0, p1);
  if (intersect) {
    pen(green, 3);
    showLine(p0, U(p0, p1));
    strokeWeight(1);
    noStroke();
  }
}
