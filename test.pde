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
    generatePointsOnSphere(P, centerOfSphere, radiusOfSphere, m);
    long startTime = System.nanoTime();
    ArrayList<Triangle> triangles = generateConvexHull(P.G, P.nv);
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
 * Test the performance of extreme-plane generation. n is the number of tests.
 * attenuation controls the size of each circle.
 */
void extremePlaneTests(int n, float attenuation) {
  float[] times = new float[n];
  int[] iters = new int[n];
  for (int i = 0; i < n; ++i) {
    RingSet ringSet = new RingSet(centerOfSphere, radiusOfSphere, 3, 3);
    ringSet.init();
    ringSet.generatePoints(attenuation);
    DebugEPInfo dInfo = new DebugEPInfo();
    long st = System.nanoTime();
    ringSet.generateExTriThreeRings(0, 1, 2, dInfo);
    long ed = System.nanoTime();
    times[i] = (ed - st) / 1000000.0;
    iters[i] = dInfo.iter;
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
      System.out.format("Basic initialization using 3 x-axes\n");
  }
  System.out.format("Extreme plane generation (n = %d, attenuation = %f): " +
                    "average time = %f ms.\n", n, attenuation, avgTime);

  int max = 0;
  for (int i = 0; i < n; ++i) {
    //System.out.format("%d-th: %d\n", i, iters[i]);
    if (iters[i] > max) max = iters[i];
  }
  System.out.format("max number of iterations = %d\n", max);
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
    P.drawBalls(2);
  }
  if (generateCH) {
    long startTime = System.nanoTime();
    ArrayList<Triangle> triangles = generateConvexHull(P.G, P.nv);
    long endTime = System.nanoTime();
    timeTM = (endTime - startTime) / 1000000.0;
    fill(red);
    stroke(0);
    numTriangles = triangles.size();
    showTriangles(triangles, P.G);
  } else {
    numTriangles = -1;
  }
}

void convexHullWithHolesTest() {
  rs.generatePoints(attenuation);
  if (showRingSet) rs.showRings();
  if (generateCH) {
    pt[] pointArray = rs.get1DPointArray();
    if (regenerateCH) {  // regenerate a convex hull
      long startTime = System.nanoTime();
      rs.generateTriangleMesh(0);
      long endTime = System.nanoTime();
      timeTM = (endTime - startTime) / 1000000.0;
    }
    fill(red);
    stroke(0);
    showTriangles(rs.triangles, pointArray);
    numTriangles = rs.triangles.size();
  } else {
    numTriangles = -1;
  }
}

void ringSetTriangulationTest() {
  rs.generatePoints(attenuation);
  if (showRingSet) rs.showRings();
  if (generateCH) {
    pt[] pointArray = rs.get1DPointArray();
    if (regenerateCH) {
      long startTime = System.nanoTime();
      rs.generateTriangleMesh(methodTM);
      long endTime = System.nanoTime();
      timeTM = (endTime - startTime) / 1000000.0;
    }
    numTriangles = 0;

    if (methodTM == 1) {
      if (show3RT) {
        assert rs.threeRingTriangles != null;
        fill(red);
        stroke(0);
        strokeWeight(2);
        showTriangles(rs.threeRingTriangles, pointArray);
        fill(#0AED62, 100);  // light green
        noStroke();
        showTriangleNormals(rs.threeRingTriangles, pointArray);
        numTriangles += rs.threeRingTriangles.size();
      }
      if (show2RT) {
        assert rs.twoRingTriangles != null;
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

void subdivisionTest() {
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
    timeTM = (endTime - startTime) / 1000000.0;

    startTime = System.nanoTime();
    TriangleMesh triMesh = new TriangleMesh(positions, triangles);
    triMesh.subdivide(subdivisionTimes, centerOfSphere, radiusOfSphere);
    endTime = System.nanoTime();
    timeSD = (endTime - startTime) / 1000000.0;

    triMesh.showTriangleMesh(red, true);
    // triMesh.showCornerPairs(blue, 3);
    // triMesh.showVertices(green, 1);
    numTriangles = triMesh.nt;
  } else {
    numTriangles = -1;
  }
}

void hubTest() {
  hub.showHub(red, 150);
  hub.showBoundingBall(blue, 100);
}

void pivotPlaneAroundLineHitCircleTest() {
  assert rs.nRings >= 3;
  rs.generatePoints(attenuation);
  assert rs.points != null && rs.centers != null && rs.normals != null && rs.radii != null;
  pt c = rs.centers[0];
  float r = rs.radii[0];
  vec n = rs.normals[0];
  pt a = rs.points[1][0];
  pt b = rs.points[2][0];
  pt[] contacts = pivotPlaneAroundLineHitCircle(c, r, n, a, b, rs.xAxes[0], rs.yAxes[0]);
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
    assert ec != null;
    exactCHEdgeCircle(ec.a, ec.b, ec.c, ec.r, ec.n, ec.vi, ec.vj);
  }
}

void exactCHTwoCirclesTest() {
  assert rs.nRings >= 2;
  rs.generatePoints(attenuation);
  assert rs.points != null && rs.centers != null && rs.normals != null &&
         rs.radii != null && rs.sameRadius == false;
  exactCHTwoCircles(rs.centers[0], rs.radii[0], rs.normals[0], rs.xAxes[0], rs.yAxes[0],
                    rs.centers[1], rs.radii[1], rs.normals[1], rs.xAxes[1], rs.yAxes[1]);
}

void tangentPlaneThreeCirclesIterTest() {
  assert rs.nRings >= 3;
  rs.generatePoints(attenuation);
  assert rs.points != null && rs.centers != null && rs.normals != null &&
         rs.radii != null && rs.sameRadius == false;

  {
    fill(red);
    disk(rs.centers[0], rs.normals[0], rs.radii[0]);
    fill(green);
    disk(rs.centers[1], rs.normals[1], rs.radii[1]);
    fill(blue);
    disk(rs.centers[2], rs.normals[2], rs.radii[2]);
  }

  pt[] ps = tangentPlaneThreeCirclesIter(rs.centers[0], rs.radii[0], rs.normals[0], rs.xAxes[0], rs.yAxes[0],
                                     rs.centers[1], rs.radii[1], rs.normals[1], rs.xAxes[1], rs.yAxes[1],
                                     rs.centers[2], rs.radii[2], rs.normals[2], rs.xAxes[2], rs.yAxes[2],
                                     null);
  {
    fill(orange, 180);
    showTriangle(ps[0], ps[1], ps[2]);
    fill(cyan, 100);
    showNormalToTriangle(ps[0], ps[1], ps[2], 20, 4);
  }
  ps = tangentPlaneThreeCirclesIter(rs.centers[0], rs.radii[0], rs.normals[0], rs.xAxes[0], rs.yAxes[0],
                                rs.centers[2], rs.radii[2], rs.normals[2], rs.xAxes[2], rs.yAxes[2],
                                rs.centers[1], rs.radii[1], rs.normals[1], rs.xAxes[1], rs.yAxes[1],
                                null);
  {
    fill(violet, 180);
    showTriangle(ps[0], ps[1], ps[2]);
    fill(magenta, 100);
    showNormalToTriangle(ps[0], ps[1], ps[2], 20, 4);
  }

  if (showRingSet) rs.showRings();
}

void exactCHThreeCirclesTest() {
  assert rs.nRings >= 3;
  rs.generatePoints(attenuation);
  assert rs.points != null && rs.centers != null && rs.normals != null &&
         rs.radii != null && rs.sameRadius == false;
  exactCHThreeCircles(rs.centers[0], rs.radii[0], rs.normals[0], rs.xAxes[0], rs.yAxes[0],
                      rs.centers[1], rs.radii[1], rs.normals[1], rs.xAxes[1], rs.yAxes[1],
                      rs.centers[2], rs.radii[2], rs.normals[2], rs.xAxes[2], rs.yAxes[2]);
}

void threeRingTriangleTest() {
  rs.generatePoints(attenuation);
  if (showRingSet) rs.showRings();
  if (generateCH) {
    pt[] pointArray = rs.get1DPointArray();
    if (regenerateCH) {
      long startTime = System.nanoTime();
      rs.generateThreeRingTriangles();
      long endTime = System.nanoTime();
      timeTM = (endTime - startTime) / 1000000.0;
    }
    numTriangles = 0;
    assert rs.threeRingTriangles != null;
    fill(red);
    stroke(0);
    strokeWeight(2);
    showTriangles(rs.threeRingTriangles, pointArray);
    fill(#0AED62, 100);  // light green
    noStroke();
    showTriangleNormals(rs.threeRingTriangles, pointArray);
    numTriangles += rs.threeRingTriangles.size();
  } else {
    numTriangles = -1;
  }
}

void exactCHNaiveTest() {
  if (P.nv < 6) {
    println("Should use at least 6 points.");
    return;
  }

  int nv = P.nv - P.nv % 2;
  RingSet rs = new RingSet(centerOfSphere, radiusOfSphere, P.G, nv);

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

  if (showAuxPlane) {
    fill(gray, 100);
    showPlane(centerOfSphere, new vec(0, 0, 1), radiusOfSphere + 5);
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
  if (P.nv < 6) {
    println("Should use at least 6 points.");
    return;
  }

  int nv = P.nv - P.nv % 2;
  RingSet rs = new RingSet(centerOfSphere, radiusOfSphere, P.G, nv);

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

  if (debugIncCH) {
    rs.showDebugIncCHInfo();
  }
}

void corridorTest() {
  if (P.nv < 6) {
    println("Should use at least 6 points.");
    return;
  }

  int nv = P.nv - P.nv % 2;
  RingSet rs = new RingSet(centerOfSphere, radiusOfSphere, P.G, nv);

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

  rs.showCorridor(idxIncCor);
}

void meshFromExactCHTest() {
  if (P.nv < 6) {
    println("Should use at least 6 points.");
    return;
  }

  int nv = P.nv - P.nv % 2;
  RingSet rs = new RingSet(centerOfSphere, radiusOfSphere, P.G, nv);

  if (!rs.isValid()) {
    validRS = false;
    return;
  } else {
    validRS = true;
  }

  rs.setNumPointsPerRing(numPointsPerRing);
  rs.generatePoints(attenuation);
  if (showRingSet) {
    rs.showRings();
  }

  if (showDiskSet) {
    rs.showDisks();
  }

  rs.generateExactCHIncremental();
  rs.generateMeshFromExactCH(0);

  if (show3RT) {
    fill(navy, 230);
    rs.showThreeRingTriangles();
  }

  if (show2RT) {
    fill(springGreen, 230);
    rs.showTwoRingTriangles();
  }

  if (showBeams) {
    fill(orange, 230);
    rs.showBeams(40.0);
  }
}