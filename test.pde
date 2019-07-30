/******************************************************************************
 * Test.
 ******************************************************************************/


int gNumTests = 100000;
int numBackup3RT = 0;
int numPointsPerTest = 64;
boolean saveFailures = false;

boolean showFirstCone = false;
boolean showSecondCone = false;
boolean showCoarseCorridor = false;

boolean gShowTwoSpheres = false;

boolean showStereoProjection = false;
boolean showUpArrow = true;

boolean gCreateGap = false;
boolean gUseTriQuadMesh = true;
boolean gProjectOnCircleAfterSub = true;

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
      if (dot(normal, AD) > gEpsilon) {
        {  // debug
          warningMsg += "dot(normal, AD) = " + dot(normal, AD) + "; ";
          fill(red);
          showTriangle(pa, pb, pc);
          fill(yellow);
          showBall(pa, 3);
          fill(green);
          showBall(pb, 3);
          fill(blue);
          showBall(pc, 3);
          fill(magenta);
          showBall(pd, 3);
          println("normal =", normal, " AD =", AD);
          // println("index of the wrong triangle =", triangles.indexOf(t));
        }
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



/******************************************************************************
 * Multiple-tests functions below.
 ******************************************************************************/



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
    // ringSet.init();
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
    // ringSet.init();
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
    // ringSet.init();
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
    // ringSet.init();
    // ringSet.generatePoints(1.0);
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
    // ringSet.init();
    // ringSet.generatePoints(1.0);
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

private float errorPointOnCone(pt p, pt apex, vec axis, float alpha) {
  vec ap = V(apex, p);
  float z = dot(ap, axis);
  float r = N(ap, axis).norm();
  return abs(r - z * tan(alpha));
}

/*
 * (u, a) defines a cone, and (v, b) defines another cone. u and v are unit
 * vectors. a and b are half-angles of cones. a or b can be bigger than PI/2.
 */
private float errorConeTangentToCone(vec u, float a, vec v, float b) {
  return abs(dot(u, v) - cos(a + b));
}


/*
 * Compare the two approaches to find the supporting plane of three given circles.
 * The first approach: directly works in 3D space, Ashish's method.
 * The second approach: stereographic projection, 2D Apollonius problem, inverse
 * of stereographic projection.
 * n is the number of tests. In each test, a valid set of circles will be randomly
 * generate and serves as input to the two approaches.
 * Average running time and accuracy of each approach are reported.
 */
void compareSupPlaneApproachTests(int n, String filename) {
  float avgTime0 = 0.0, avgTime1 = 0.0, diffAvgTime = 0.0;
  float avgError0 = 0.0, avgError1 = 0.0, diffAvgError = 0.0;
  String[] lines = new String[2 * n + 4];
  int validCounts = n;
  for (int i = 0; i < n; ++i) {
    RingSet rs = new RingSet(gSphereCenter, gSphereRadius, 3);
    int[] ringIDs = new int[] {0, 1, 2};
    float[] alphas = new float[3];
    alphas[0] = asinClamp(rs.radii[ringIDs[0]] / rs.sphereRadius);
    alphas[1] = asinClamp(rs.radii[ringIDs[1]] / rs.sphereRadius);
    alphas[2] = asinClamp(rs.radii[ringIDs[2]] / rs.sphereRadius);
    vec[] normals = new vec[] {rs.normals[ringIDs[0]], rs.normals[ringIDs[1]],
                               rs.normals[ringIDs[2]]};

    {
      // println("alphas =", alphas[0], ", ", alphas[1], ", ", alphas[2]);
      // println("normals =", normals[0], ", ", normals[1], ", ", normals[2]);
    }

    float time0 = 0.0, time1 = 0.0, timeDiff = 0.0;
    float error0 = 0.0, error1 = 0.0, errorDiff = 0.0;
    boolean valid0 = true, valid1 = true;
    {  // first method (3D)
      // compute time
      long startTime = System.nanoTime();
      pt[] ps = rs.twoSupPlanesThreeCircles(ringIDs[0], ringIDs[1], ringIDs[2], null, null, true);
      long endTime = System.nanoTime();
      time0 = (endTime - startTime) / 1000000.0;

      if (ps == null) {
        valid0 = false;
      } else {  // compute error
        for (int j = 0; j < 3; ++j) error0 += errorPointOnCone(ps[j], rs.sphereCenter, normals[j], alphas[j]);
        Circle ap = circumcircleOfTriangle(ps[0], ps[1], ps[2]);
        float alpha = asinClamp(ap.r / rs.sphereRadius);
        if (dot(ap.n, V(rs.sphereCenter, ap.c)) < 0) alpha = PI - alpha;
        for (int j = 0; j < 3; ++j) error0 += errorConeTangentToCone(ap.n, alpha, normals[j], alphas[j]);
      }
    }
    {  // second method (stereographic projection + 2D)
      long startTime = System.nanoTime();
      int f = rs.getMaxCircleID();  // the ID of the infinity circle
      pt northPole = rs.contacts[f];
      StereoProjector sp = new StereoProjector(rs.sphereCenter, rs.sphereRadius, StereoType.CUSTOM_NORTH_POLE_SOUTH_PLANE);
      sp.setNorthPole(northPole);
      pt[] ps = rs.oneSupPlaneThreeCirclesStereoApollonius(ringIDs[0], ringIDs[1], ringIDs[2], sp, f);
      long endTime = System.nanoTime();
      time1 = (endTime - startTime) / 1000000.0;

      if (ps == null) {
        valid1 = false;
      } else {  // compute error
        for (int j = 0; j < 3; ++j) error1 += errorPointOnCone(ps[j], rs.sphereCenter, normals[j], alphas[j]);
        Circle ap = circumcircleOfTriangle(ps[0], ps[1], ps[2]);
        float alpha = asinClamp(ap.r / rs.sphereRadius);
        if (dot(ap.n, V(rs.sphereCenter, ap.c)) < 0) alpha = PI - alpha;
        for (int j = 0; j < 3; ++j) error1 += errorConeTangentToCone(ap.n, alpha, normals[j], alphas[j]);
      }
    }
    {  // compare two errors
      if (valid0 && valid1 && error1 > error0 * 10000) {  // if the second approach produces very high error w.r.t. the first approach
        rs.save("data/tmp/rs_high_error_" + str(i));
      }
    }
    String timeLine, errorLine;
    if (valid0 && valid1) {
      timeDiff = time0 - time1;
      timeLine = str(time0) + ", " + str(time1) + ", " + str(timeDiff);
      errorDiff = error0 - error1;
      errorLine = str(error0) + ", " + str(error1) + ", " + str(errorDiff);
      // println("time line:", timeLine);
      // println("error line:", errorLine);
      avgTime0 += time0;
      avgTime1 += time1;
      avgError0 += error0;
      avgError1 += error1;
    } else {
      if (valid0) {
        timeLine = str(time0) + ", INVALID, INVALID";
        errorLine = str(error0) + ", INVALID, INVALID";
      } else if (valid1) {
        timeLine = "INVALID, " + str(time1) + ", INVALID";
        errorLine = "INVALID, " + str(error1) + ", INVALID";
      } else {
        timeLine = "INVALID, INVALID, INVALID";
        errorLine = "INVALID, INVALID, INVALID";
      }
      validCounts -= 1;  // decrease valid counts
    }
    lines[2 * i] = timeLine;
    lines[2 * i + 1] = errorLine;
  }
  avgTime0 /= validCounts;
  avgTime1 /= validCounts;
  avgError0 /= validCounts;
  avgError1 /= validCounts;
  diffAvgTime = avgTime0 - avgTime1;
  diffAvgError = avgError0 - avgError1;
  String timeStatLine = str(avgTime0) + ", " + str(avgTime1) + ", " + str(diffAvgTime);
  println("avg time and diff:", timeStatLine);
  String errorStatLine = str(avgError0) + ", " + str(avgError1) + ", " + str(diffAvgError);
  println("avg error and diff:", errorStatLine);
  lines[2 * n] = timeStatLine;
  lines[2 * n + 1] = errorStatLine;
  lines[2 * n + 2] = "valid counts: " + str(validCounts);
  lines[2 * n + 3] = "total counts: " + str(n);
  if (filename != null) saveStrings(filename + "_" + str(n), lines);
}


/******************************************************************************
 * Single-test functions below.
 ******************************************************************************/



void convexHullTest() {
  if (gShowPointSet) {
    fill(green);
    gPoints.drawBalls(2);
  }
  if (gGenerateCH) {
    long startTime = System.nanoTime();
    ArrayList<Triangle> triangles = generateConvexHull(gPoints.G, gPoints.nv);
    long endTime = System.nanoTime();
    gTimeMeshing = (endTime - startTime) / 1000000.0;
    fill(red);
    stroke(0);
    gNumTriangles = triangles.size();
    showTriangles(triangles, gPoints.G);
    noStroke();
  } else {
    gNumTriangles = -1;
  }
}

void convexHullWithHolesTest() {
  gRingSet.generatePoints(gAttenuation);
  assert gRingSet.points != null;
  if (gShowRingSet) gRingSet.show();
  if (gGenerateCH) {
    pt[] pointArray = gRingSet.get1DPointArray();
    if (gRegenerateCH) {  // regenerate a convex hull
      long startTime = System.nanoTime();
      gRingSet.generateTriangleMesh(0);
      long endTime = System.nanoTime();
      gTimeMeshing = (endTime - startTime) / 1000000.0;
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
  if (gShowRingSet) gRingSet.show();
  if (gGenerateCH) {
    pt[] pointArray = gRingSet.get1DPointArray();
    if (gRegenerateCH) {
      long startTime = System.nanoTime();
      gRingSet.generateTriangleMesh(methodTM);
      long endTime = System.nanoTime();
      gTimeMeshing = (endTime - startTime) / 1000000.0;
    }
    gNumTriangles = 0;

    if (methodTM == 1) {
      if (gShow3RT) {
        assert gRingSet.threeRingTriangles != null;
        fill(blue);
        stroke(0);
        showTriangles(gRingSet.threeRingTriangles, pointArray);
        noStroke();
        fill(cyan, 100);
        showTriangleNormals(gRingSet.threeRingTriangles, pointArray);
        gNumTriangles += gRingSet.threeRingTriangles.size();
      }
      if (gShow2RT) {
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
  if (gShowRingSet) gRingSet.show();
  if (gGenerateCH) {
    ArrayList<pt> positions = gRingSet.get1DPointArrayList();
    long startTime = System.nanoTime();
    ArrayList<Triangle> triangles = generateConvexHull(gRingSet.get2DPointArray(),
                                                       gRingSet.getNumRings(),
                                                       gRingSet.getNumPointsPerRing());
    long endTime = System.nanoTime();
    gTimeMeshing = (endTime - startTime) / 1000000.0;

    startTime = System.nanoTime();
    TriangleMesh triMesh = new TriangleMesh(positions, triangles);
    triMesh.subdivide(gSubdivisonTimes);

    if (gProjectOnSphere) {
      triMesh.projectOnSphere(gSphereCenter, gSphereRadius);
    }

    endTime = System.nanoTime();
    gTimeSubdivision = (endTime - startTime) / 1000000.0;

    triMesh.show(red, true);
    // triMesh.showCornerPairs(blue, 3);
    // triMesh.showVertices(green, 1);
    gNumTriangles = triMesh.nt;
  } else {
    gNumTriangles = -1;
  }
}

void supPlaneLineCircleTest() {
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
  assert gEdgeCircle != null;
  exactCHEdgeCircle(gEdgeCircle.a, gEdgeCircle.b, gEdgeCircle.c, gEdgeCircle.r,
                    gEdgeCircle.n, gEdgeCircle.vi, gEdgeCircle.vj);
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

  if (gShowRingSet) gRingSet.show();
}

void exactCHThreeCirclesIterTest() {
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
  if (gShowRingSet) gRingSet.show();
  if (gGenerateCH) {
    pt[] pointArray = gRingSet.get1DPointArray();
    if (gRegenerateCH) {
      long startTime = System.nanoTime();
      gRingSet.generateThreeRingTriangles();
      long endTime = System.nanoTime();
      gTimeMeshing = (endTime - startTime) / 1000000.0;
    }
    gNumTriangles = 0;
    assert gRingSet.threeRingTriangles != null;
    fill(red);
    stroke(0);
    strokeWeight(2);
    showTriangles(gRingSet.threeRingTriangles, pointArray);
    fill(lightGreen, 100);
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

  if (!gRingSet.isValid()) return;


  if (gShowCircleSet) {
    gRingSet.showCircles(null);
  }

  if (gShowDiskSet) {
    gRingSet.showDisks(null);
  }

  gRingSet.generateExTrisNaive();
  if (gShow3RT) {
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
  if (gShow2RT) {
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

  if (!gRingSet.isValid()) return;

  if (gShowCircleSet) {
    if (debugIncCH && nv > 4) gRingSet.showCircles(debugIncCHIter);
    else gRingSet.showCircles(null);
  }

  if (gShowDiskSet) {
    fill(red);
    if (debugIncCH && nv > 4) gRingSet.showDisks(debugIncCHIter);
    else gRingSet.showDisks(null);
  }

  gRingSet.generateExactCHIncremental(null);

  if (gShowTriangleFaces) {
    fill(blue);
    gRingSet.showIncTriangles();
  }

  if (gShowCorridorFaces) {
    fill(green);
    gRingSet.showIncCorridors();
  }

  if (gShowPolygons) {
    gRingSet.generateConvexTriMesh();
    fill(red);
    gRingSet.showPolygons();
  }

  if (debugIncCH) {
    gRingSet.showDebugIncCHInfo();
  }

  // if (gShowApolloniusDiagram) {
  //   gRingSet.showApolloniusDiagram();
  // }

  // if (debugApolloniusDiagram) {
  //   gRingSet.showADDebugInfo();
  // }

  // {
  //   boolean b0 = gRingSet.incCorridors.get(0).isUniformlySampled();
  //   gRingSet.incCorridors.get(0).testCorrespondence();
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

    ArrayList<ArrayList<Circle>> apolloniusEdgeCircles = new ArrayList();
    for (RingSet.IncCorridor cor : gRingSet.incCorridors) {
      apolloniusEdgeCircles.add(cor.apolloniusCircles());
    }
    ArrayList<ArrayList<pt>> pApolloniusEdges = new ArrayList();
    for (ArrayList<Circle> edgeCircles : apolloniusEdgeCircles) {
      ArrayList<pt> edge = new ArrayList<pt>();
      for (Circle cir : edgeCircles) {
        Circle pCir = sp.project(cir);
        edge.add(pCir.c);
      }
      pApolloniusEdges.add(edge);
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
      pt southPole = P(sp.center, -1.3, V(sp.center, sp.northPole));
      showPlane(southPole, U(sp.center, sp.northPole), 50 * sp.radius);
    }
  }
}

/* Pick a specific corridor and analyse its property. */
void corridorTest() {
  if (gPoints.nv < 6) {
    println("Should use at least 6 points.");
    return;
  }

  int nv = gPoints.nv - gPoints.nv % 2;
  gRingSet = new RingSet(gSphereCenter, gSphereRadius, gPoints.G, nv);

  if (!gRingSet.isValid()) return;


  if (gShowCircleSet) {
    gRingSet.showCircles(null);
  }

  if (gShowDiskSet) {
    gRingSet.showDisks(null);
  }

  gRingSet.generateExactCHIncremental(null);

  if (gShowTriangleFaces) {
    fill(blue);
    gRingSet.showIncTriangles();
  }

  if (gShowCorridorFaces) {
    fill(green);
    gRingSet.showIncCorridors();
  }

  gRingSet.showCorridor(gIdxIncCorridor);
}

void meshFromExactCHTest() {
  if (gPoints.nv < 6) {
    println("Should use at least 6 points.");
    return;
  }

  int nv = gPoints.nv - gPoints.nv % 2;
  gRingSet = new RingSet(gSphereCenter, gSphereRadius, gPoints.G, nv);

  if (!gRingSet.isValid()) return;


  gRingSet.setNumPointsPerRing(gNumPointsPerRing);
  gRingSet.generatePoints(gAttenuation);
  if (gShowRingSet) gRingSet.show();


  if (gShowDiskSet) {
    gRingSet.showDisks(null);
  }

  gRingSet.generateExactCHIncremental(null);
  gTriangleMesh = gRingSet.generateMeshFromExactCH(0);

  if (gSubdivisonTimes > 0) {
    gTriangleMesh.subdivide(gSubdivisonTimes);
  }

  if (gShowTriMesh) {
    if (gTriangleMesh != null) {
      gTriangleMesh.show(khaki, true);
    }
  } else {
    if (gShow3RT) {
      fill(navy, 230);
      gRingSet.showThreeRingTriangles();
    }
    if (gShow2RT) {
      fill(springGreen, 230);
      gRingSet.showTwoRingTriangles();
    }
  }

  if (gShowBeams) {
    fill(orange, 230);
    gRingSet.showBeams(40.0);
  }
}

private void hubTest() {
  assert gHub != null;
  gHub.generateIntersectionCircles();
  if (gShowIntersectionCircles) {
    gHub.showIntersectionCircles();
  }

  gRingSet = gHub.circlesToRingSet(gNumPointsPerRing);
  gRingSet.generateExactCHIncremental(null);
  boolean containHubCenter = gRingSet.pointInCrudestConvexHull(gHub.ball.c);

  if (gShowDiskSet) {
    fill(red);
    gRingSet.showDisks(null);
  }
  if (gShowTriangleFaces) {
    fill(blue);
    gRingSet.showIncTriangles();
  }
  if (gShowCorridorFaces) {
    fill(green);
    gRingSet.showIncCorridors();
  }

  /* Generate a triangle-quad mesh. */
  gTriQuadMesh = gRingSet.generateConvexTriQuadMesh();
  gTriQuadMesh.setColors(cyan, lime);  // set triangle and corridor colors

  /* Subdivide the triangle-quad mesh. */
  for (int i = 0; i < gSubdivisonTimes; ++i) {
    gTriQuadMesh = gTriQuadMesh.subdivide(SubdivideTypeTriangle.LOOP, SubdivideTypeQuad.DIAMOND, gProjectOnCircleAfterSub);
  }

  /* Project the mesh, face-by-face, on the hub. */
  ProjectType projType = null;
  if (gMethodProjection == 1) projType = ProjectType.RAY;
  else if (gMethodProjection == 2) projType = ProjectType.SPHERE_TRACING;

  if (gSubdivisonTimes > 0 && projType != null) {
    // gTriQuadMesh.projectOnHub(gHub, projType);
    gTriangleMesh = gTriQuadMesh.toTriangleMesh();
    gTriangleMesh.containHubCenter = containHubCenter;

    if (debugProjection) dProjectionInfo.reset();
    gTriangleMesh.projectOnHub(gHub, projType);
    if (debugProjection) {
      dProjectionInfo.printDists();
      dProjectionInfo.printTs();
      dProjectionInfo.show();
    }
  }



  if (gShowTriMesh) {
    if (gSubdivisonTimes > 0 && projType != null) {
      gTriangleMesh.show(cyan, gShowTriangleStrokes);
    } else {
      gTriQuadMesh.show(gShowTriangleStrokes);
    }
  }

  if (gShowHub) {
    gHub.show(lightSalmon, 130);
  }

  if (gShowBoundingSphere) {
    gHub.showBoundingSphere(khaki, 200);
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

  if (gShowIntersectionCircles) {
    gHub.showIntersectionCircles();
  }

  gRingSet = gHub.circlesToRingSet(gNumPointsPerRing);
  gHub.generateBeamSamples(gNumPointsPerRing, gRingSet.xAxes);

  /* Generate convex hull. */
  gRingSet.generateExactCHIncremental(null);
  boolean containHubCenter = gRingSet.pointInCrudestConvexHull(gHub.ball.c);

  if (gInsertInfCircle) {
    if (containHubCenter == false) {
      // println("insert a new circle!");
      gRingSet = gRingSet.insertInfinitesimalCircle();
      gRingSet.generateExactCHIncremental(null);

    }
  }

  if (gShowDiskSet) {
    fill(red);
    gRingSet.showDisks(null);
  }
  if (gShowTriangleFaces) {
    fill(blue);
    gRingSet.showIncTriangles();
  }
  if (gShowCorridorFaces) {
    fill(green);
    gRingSet.showIncCorridors();
  }

  if (!gUseTriQuadMesh) {  // use a triangle mesh for the convex hull
    gTriangleMesh = gRingSet.generateConvexTriMesh();

    if (gCreateGap) {
      /* Generate beam mesh. */
      gBeamMesh = gHub.generateBeamMesh();
      if (gShowLiftedCones) {
        // gBeamMesh.show(cyan, gShowTriangleStrokes);  // each beam is an n-sided polygon with small n
        gHub.showLiftedCones(cyan, 255);  // each beam is an n-sided cylinder with very large n
      }

      /* Generate gap mesh. */
      gHub.initGaps(gRingSet.borders);
      gGapMesh = gHub.generateGapMesh();
      if (gShowGapMesh) gGapMesh.show(orange, gShowTriangleStrokes);
      // gTriangleMesh.augment(gGapMesh);  // merge convex hull mesh and gap mesh
    } else {  // don'create a gap between each beam and the convex hull
      int nNeighbors = gHub.nNeighbors;
      gGapMesh = new TriangleMesh();
      for (int i = 0; i < nNeighbors; ++i) {
        ArrayList<pt> innerLoop = gRingSet.borders[i];
        ArrayList<pt> outerLoop = new ArrayList<pt>();
        pt[] samples = gHub.liftedCones[i].samples;
        assert samples != null;
        int n = samples.length / 2;
        for (int j = n; j < samples.length; ++j) outerLoop.add(samples[j]);
        ConvexGap gap = new ConvexGap(innerLoop, outerLoop);  // a beam is actually a big gap
        TriangleMesh tm = gap.toTriMesh();
        gGapMesh.augmentWithShift(tm.positions, tm.triangles);
        if (gShowGapMesh) gGapMesh.show(lightSalmon, gShowTriangleStrokes);
      }
    }

    gTriangleMesh.subdivide(gSubdivisonTimes);

    ProjectType projType = null;
    if (gMethodProjection == 1) projType = ProjectType.RAY;
    else if (gMethodProjection == 2) projType = ProjectType.SPHERE_TRACING;

    gTriangleMesh.containHubCenter = containHubCenter;
    gTriangleMesh.projectOnHub(gHub, projType);

    if (gShowTriMesh) {
      gTriangleMesh.show(purple, gShowTriangleStrokes);
    }
  } else {  // use a triangle-quad mesh for the convex hull
    /* Generate a triangle-quad mesh. */
    gTriQuadMesh = gRingSet.generateConvexTriQuadMesh();
    gTriQuadMesh.setColors(blue, green);  // set triangle and corridor colors

    /* Subdivide the triangle-quad mesh. */
    for (int i = 0; i < gSubdivisonTimes; ++i) {
      gTriQuadMesh = gTriQuadMesh.subdivide(SubdivideTypeTriangle.LOOP, SubdivideTypeQuad.DIAMOND, gProjectOnCircleAfterSub);
    }

    /* Project the mesh on the hub. */
    ProjectType projType = null;
    if (gMethodProjection == 1) projType = ProjectType.RAY;
    else if (gMethodProjection == 2) projType = ProjectType.SPHERE_TRACING;

    if (gSubdivisonTimes > 0 && projType != null) {
      // gTriQuadMesh.projectOnHub(gHub, projType);
      if (debugProjection) {
        dProjectionInfo.show();
        dProjectionInfo.printDists();
        dProjectionInfo.printPointsInside();
      }

      gTriangleMesh = gTriQuadMesh.toTriangleMesh();
      gTriangleMesh.containHubCenter = containHubCenter;
      gTriangleMesh.projectOnHub(gHub, projType);
    }

    if (gShowTriMesh) {
      if (gSubdivisonTimes > 0 && projType != null) {
        gTriangleMesh.show(cyan, gShowTriangleStrokes);
      } else {
        if (gSubdivisonTimes == 0) gTriQuadMesh.setColors(blue, green);
        else gTriQuadMesh.setColors(cyan, green);
        gTriQuadMesh.show(gShowTriangleStrokes);
      }
    }

    /* Construct beams. */
    int nNeighbors = gHub.nNeighbors;
    gGapMesh = new TriangleMesh();
    for (int i = 0; i < nNeighbors; ++i) {
      ArrayList<pt> innerLoop = gTriQuadMesh.borders[i];
      ArrayList<pt> outerLoop = new ArrayList<pt>();
      pt[] samples = gHub.liftedCones[i].samples;
      assert samples != null;
      int n = samples.length / 2;
      for (int j = n; j < samples.length; ++j) outerLoop.add(samples[j]);
      ConvexGap gap = new ConvexGap(innerLoop, outerLoop);  // a beam is actually a big gap
      TriangleMesh tm = gap.toTriMesh();

      if (gShowExplodedView) {
        vec dir = V(gDeltaExplodedView, gHub.tCones[i].normal);
        tm.translate(dir);
      }

      gGapMesh.augmentWithShift(tm.positions, tm.triangles);
      if (gShowGapMesh) gGapMesh.show(lightSalmon, gShowTriangleStrokes);
    }
  }

  if (gShowHub) {
    if (!gShowExplodedView) {
      gHub.show(lightSalmon, 255);  // alpha: 130, 255
    } else {
      gHub.showBeamsExplodedView(gDeltaExplodedView, 40, false, lightSalmon, 255);
    }
  }

  if (gShowBoundingSphere) {
    gHub.showBoundingSphere(khaki, 255);  // alpha: 100, 255
  }
}

void staticHubTest() {
  gFocus = gHub.ball.c;
  {
    // fill(black);
    // showBall(gFocus, 1);
  }
  hubTest();
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

  if (!gRingSet.isValid()) return;


  if (gShowDiskSet) {
    fill(red);
    gRingSet.showDisks(null);
  }

  gRingSet.generateExactCHIncremental(null);

  if (gShowTriangleFaces) {
    fill(blue);
    gRingSet.showIncTriangles();
  }

  if (gShowCorridorFaces) {
    fill(green);
    gRingSet.showIncCorridors();
  }

  return;
}

/*
 * Show the two cones corresponding to the two supporting planes of three given
 * circles on a sphere.
 */
void supPlaneThreeCirclesTest() {
  if (gPoints.nv < 6) return;
  gRingSet = new RingSet(gSphereCenter, gSphereRadius, gPoints.G, 6);

  if (!gRingSet.isValid()) return;


  if (gShowDiskSet) {
    fill(red);
    gRingSet.showDisks(null);
  }
  gRingSet.generateExactCHIncremental(null);

  if (gShowTriangleFaces) {
    fill(blue);
    gRingSet.showIncTriangles();
  }

  if (gShowCones) {
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

  if (gShowCircleSet) {
    hint(DISABLE_DEPTH_TEST);
    stroke(black);
    gRingSet.showCircles(null);
    noStroke();
    hint(ENABLE_DEPTH_TEST);
  }
}

/*
 * Given two circles, tangent to each other, on a sphere, test whether the two
 * cap centers (on sphere), the sphere center, and the point of tangency are
 * coplanar.
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

  if (gShowAuxPlane) {
    fill(orange, 100);
    showPlane(gSphereCenter, points[0], points[2], 100);
  }

  hint(DISABLE_DEPTH_TEST);
  pen(cyan, 3);
  showSegment(gSphereCenter, points[0]);
  showSegment(gSphereCenter, points[2]);
  strokeWeight(3);
  showCircle(gSphereCenter, n1, gSphereRadius);
  strokeWeight(1);
  noStroke();
  hint(ENABLE_DEPTH_TEST);

  if (gShowCircleSet) {
    gRingSet.showCircles(null);
  }
}


void ellipticConeTest() {
  if (gPoints.nv < 4) return;

  gRingSet = new RingSet(gSphereCenter, gSphereRadius, gPoints.G, 4);

  if (!gRingSet.isValid()) return;


  if (gShowCircleSet) {
    gRingSet.showCircles(null);
  }

  if (gShowDiskSet) {
    gRingSet.showDisks(null);
  }

  gRingSet.generateExactCHIncremental(null);

  if (gShowCorridorFaces) {
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

  gHub.setNumEdgesRegularPolygon(gNumPointsPerRing);
  gHub.setGapDistance(gGapWidth);
  gHub.setXDirestions(xAxes);
  gHub.setIsHalved(true);
  gTriangleMesh = gHub.generateTriMesh();

  if (gShowTriMesh) {
    gTriangleMesh.show(cyan, true);
  }

  if (gShowHub) {
    gHub.show(lightSalmon, 130);
  }
}

void twoEllipticConesTest() {
  if (gPoints.nv < 4) {
    println("Should use at least 4 points.");
    return;
  }

  gRingSet = new RingSet(gSphereCenter, gSphereRadius, gPoints.G, 4);

  if (gShowEllipticCone1) gRingSet.showEllipticCone(0, 1, 0, red, color(green));  // alpha: 200
  if (gShowEllipticCone2) gRingSet.showEllipticCone(0, 1, 1, red, color(blue));  // alpha: 200

  if (gShowDiskSet) {
    fill(red, 200);
    gRingSet.showDisks(2);
  }
  if (gShowCircleSet) gRingSet.showCircles(2);
  if (gShowAuxPlane) {
    fill(orange, 180);
    showPlane(gSphereCenter, gPoints.G[0], gPoints.G[2], gSphereRadius);
  }
}

/*
 * Read a hub from an augmented hub file and convert this hub to a mesh.
 */
void staticAugHubToMeshTest() {
  gHub = new Hub();
  gHub.loadAugFile("data/hub_aug_unnamed");
  gHub.setGapDistance(gGapWidth);
  gHub.setIsHalved(true);
  gTriangleMesh = gHub.generateTriMesh();

  if (gShowTriMesh) gTriangleMesh.show(cyan, true);
  if (gShowHub) gHub.show(lightSalmon, 130);
}

void convexGapTest() {
  assert gGap != null && gGap.points0 != null && gGap.points1 != null;
  warningMsg = "";
  gFocus = gGap.center();

  fill(magenta);
  gGap.show();

  gTriangleMesh = gGap.toTriMesh();
  if (gTriangleMesh == null) return;
  if (gTriangleMesh.isConvex() == false) warningMsg += "The gap mesh isn't convex; ";
  gTriangleMesh.show(cyan, true);

  {  // debug
    // ArrayList<pt> ps = gTriangleMesh.positions;
    // ArrayList<Triangle> ts = gTriangleMesh.triangles;

    // println("Look at triangles");
    // for (Triangle t : ts) println("triangle: ", t);

    // show first triangle
    // Triangle t0 = ts.get(0);
    // fill(red);
    // showTriangle(ps.get(t0.a), ps.get(t0.b), ps.get(t0.c));
    // println("t0 =", t0);

    // show bad triangle
    // Triangle t = ts.get(3);
    // fill(blue);
    // showTriangle(ps.get(t.a), ps.get(t.b), ps.get(t.c));
    // println("t =", t);

    // fill(yellow);
    // showBall(ps.get(t.a), 3);
    // fill(magenta);
    // showBall(ps.get(t.b), 3);
    // fill(green);
    // showBall(ps.get(t.c), 3);
  }
}

void staticIncCHTest() {
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

  if (gShowDiskSet) {
    fill(red);
    gRingSet.showDisks(null);
  }

  if (gShowCircleSet) {
    gRingSet.showCircles(null);
  }

  // gRingSet.generateExTrisNaive();
  // if (gShow3RT) {
  //   fill(blue, 200);
  //   gRingSet.showExTris();
  // }

  gRingSet.generateExactCHIncremental(null);

  if (gShowTriangleFaces) {
    fill(blue);
    gRingSet.showIncTriangles();
  }

  if (gShowCorridorFaces) {
    fill(green);
    gRingSet.showIncCorridors();
  }

  gTriangleMesh = gRingSet.generateConvexTriMesh();
  if (gShowTriMesh) gTriangleMesh.show(cyan, true);

  if (gShowPolygons) {
    fill(red);
    gRingSet.showPolygons();
  }

  //fill(yellow, 100);
  //gRingSet.showSphere();
}

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

  if (gShowHub) {
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

  if (gShowTwoSpheres) {  // show the two spheres
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
  gRingSet = new RingSet(gSphereCenter, gSphereRadius, gPoints.G, gPoints.nv - gPoints.nv % 2);

  if (gShowDiskSet) {
    fill(red);
    gRingSet.showDisks(null);
  }

  // if (showUpArrow) {
  //   fill(green);
  //   arrow(gRingSet.sphereCenter, V(0, 0, gRingSet.sphereRadius), 5);
  // }

  // compute cones of circles of the ring set
  int[] ringIDs = new int[] {0, 1, 2};
  float[] alphas = new float[3];
  alphas[0] = asinClamp(gRingSet.radii[ringIDs[0]] / gRingSet.sphereRadius);
  alphas[1] = asinClamp(gRingSet.radii[ringIDs[1]] / gRingSet.sphereRadius);
  alphas[2] = asinClamp(gRingSet.radii[ringIDs[2]] / gRingSet.sphereRadius);
  vec[] normals = new vec[] {gRingSet.normals[ringIDs[0]],
                             gRingSet.normals[ringIDs[1]],
                             gRingSet.normals[ringIDs[2]]};
  // println("alphas =", alphas[0], ", ", alphas[1], ", ", alphas[2]);
  // println("normals =", normals[0], normals[1], normals[2]);

  if (showStereoProjection) {
    int f = gRingSet.getMaxCircleID();
    pt northPole = gRingSet.contacts[f];
    // pt northPole = gRingSet.contacts[0];
    // pt northPole = P(gRingSet.sphereCenter, V(0, 0, gRingSet.sphereRadius));
    StereoProjector sp = new StereoProjector(gRingSet.sphereCenter, gRingSet.sphereRadius, StereoType.CUSTOM_NORTH_POLE_SOUTH_PLANE);
    if (sp.sType == StereoType.CUSTOM_NORTH_POLE_SOUTH_PLANE) {
      sp.setNorthPole(northPole);
    }

    Circle cir0 = gRingSet.getCircle(ringIDs[0]);
    Circle cir1 = gRingSet.getCircle(ringIDs[1]);
    Circle cir2 = gRingSet.getCircle(ringIDs[2]);
    Circle pCir0 = sp.project(cir0);
    Circle pCir1 = sp.project(cir1);
    Circle pCir2 = sp.project(cir2);

    {
      // check projection
      // strokeWeight(3);
      // stroke(red);
      // cir0.show();
      // pCir0.show();
      // stroke(blue);
      // cir1.show();
      // pCir1.show();
      // noStroke();
      // strokeWeight(1);

      // fill(tomato);
      // showBall(northPole, 3);

      // check inversion of a point
      // pt c1 = gRingSet.contacts[1];  // the center of the second spherical cap
      // pt pc1 = sp.project(c1);
      // pt d1 = sp.inverse(pc1);
      // if (!samePt(c1, d1)) {
      //   fill(cyan);
      //   showBall(c1, 3);
      //   fill(magenta);
      //   showBall(d1, 3);
      //   println("c1 =", c1, "d1 =", d1);
      //   println("diff of c1 and d1 =", V(c1, d1));
      // }

      // check inversion of a circle
      // Circle invCir1 = sp.inverse(pCir1);
      // strokeWeight(3);
      // stroke(pink);
      // invCir1.show();
      // noStroke();
      // strokeWeight(1);
      // println("cir1.c =", cir1.c, "invCir1.c =", invCir1.c);
    }

    Circle2 qCir0 = sp.to2D(pCir0);
    Circle2 qCir1 = sp.to2D(pCir1);
    Circle2 qCir2 = sp.to2D(pCir2);
    Circle2[] qCirx = constructApolloniuCircles(qCir0, qCir1, qCir2, -1, 1, 1);
    Circle pCirx0 = sp.to3D(qCirx[0]);
    Circle cirx0 = sp.inverse(qCirx[0]);
    Circle pCirx1 = sp.to3D(qCirx[1]);
    Circle cirx1 = sp.inverse(qCirx[1]);

    {  // show input circles and Apollonius circles (in 2D and on sphere)
      // strokeWeight(3);
      // stroke(red);
      // cir0.show();
      // pCir0.show();
      // stroke(green);
      // cir1.show();
      // pCir1.show();
      // stroke(snow);
      // cir2.show();
      // pCir2.show();

      // stroke(blue);
      // cirx0.show();
      // pCirx0.show();
      // cirx1.show();
      // pCirx1.show();
      // noStroke();
      // strokeWeight(1);
    }

    {  // first method to compute the three contact points of the supporting plane
      pt[] ps = gRingSet.oneSupPlaneThreeCirclesStereoApollonius(ringIDs[0], ringIDs[1], ringIDs[2], sp, f);
      {  // show points of tangency
        fill(lime);
        showBall(ps[0], 3);
        showBall(ps[1], 3);
        showBall(ps[2], 3);
      }
      Circle cirx = circumcircleOfTriangle(ps[0], ps[1], ps[2]);
      {  // show the Apollonius circle and the triangle
        stroke(blue);
        strokeWeight(3);
        cirx.show();
        strokeWeight(1);
        noStroke();
        // fill(navy);
        // showTriangle(ps[0], ps[1], ps[2]);
        // fill(cyan);
        // showTriangleNormal(ps[0], ps[1], ps[2], 50, 1);
      }

      // compute error
      float errorOnCone = 0, errorTangency = 0;
      for (int j = 0; j < 3; ++j) errorOnCone += errorPointOnCone(ps[j], gRingSet.sphereCenter, normals[j], alphas[j]);
      println("2D Apollonius + stereographic projection: error on cone =", errorOnCone);
      float alpha = asinClamp(cirx.r / gRingSet.sphereRadius);
      if (dot(cirx.n, V(gRingSet.sphereCenter, cirx.c)) < 0) alpha = PI - alpha;
      for (int j = 0; j < 3; ++j) errorTangency += errorConeTangentToCone(cirx.n, alpha, normals[j], alphas[j]);
      println("2D Apollonius + stereographic projection: error tangency =", errorTangency);
    }

    {  // second method to compute the three contact points of the supporting plane
      pt[] ps = gRingSet.twoSupPlanesThreeCircles(ringIDs[0], ringIDs[1], ringIDs[2], null, null, true);
      Circle cirx = circumcircleOfTriangle(ps[0], ps[1], ps[2]);
      float errorOnCone = 0, errorTangency = 0;
      for (int j = 0; j < 3; ++j) errorOnCone += errorPointOnCone(ps[j], gRingSet.sphereCenter, normals[j], alphas[j]);
      println("3D Apollonius: error on cone =", errorOnCone);
      float alpha = asinClamp(cirx.r / gRingSet.sphereRadius);
      if (dot(cirx.n, V(gRingSet.sphereCenter, cirx.c)) < 0) alpha = PI - alpha;
      for (int j = 0; j < 3; ++j) errorTangency += errorConeTangentToCone(cirx.n, alpha, normals[j], alphas[j]);
      println("3D Apollonius: error tangency =", errorTangency);
    }

    fill(orange);
    vec southPlaneNormal = U(sp.center, sp.northPole);
    pt southPole = P(sp.center, -1.001, V(sp.center, sp.northPole));
    showPlane(southPole, southPlaneNormal, 10 * gRingSet.sphereRadius);
  }
}

void latticeTest() {
  assert gLattice != null;
  if (gShowLattice) {
    fill(red);
    gLattice.show();
  }

  if (gShowBoundingSphere) {
    fill(yellow, 100);
    gLattice.showInflatingSpheres();
  }

  if (gShowTriMesh) {
    long startTime = System.nanoTime();
    gTriangleMesh = gLattice.triangulate();
    long endTime = System.nanoTime();
    gTimeMeshing = (endTime - startTime) / 1000000.0;

    if (gShowCHoCCs) {
      gLattice.showCHoCCs(blue, green, lightSalmon);
    } else {
      gTriangleMesh.show(cyan, true);
    }

    if (debugLattice) {
      gAvgVertexCount = average(gLattice.junctionVertexCounts);
    }
  } else if (gUseTriQuadMesh) {
    ProjectType projType = null;
    if (gMethodProjection == 1) projType = ProjectType.RAY;
    else if (gMethodProjection == 2) projType = ProjectType.SPHERE_TRACING;

    if (debugProjection) dProjectionInfo.reset();

    long startTime = System.nanoTime();
    gLattice.tessellate(gSubdivisonTimes, projType);
    long endTime = System.nanoTime();
    gTimeMeshing = (endTime - startTime) / 1000000.0;

    gLattice.showJunctions(cyan, true);
    if (gShowBeams) gLattice.showBeams(lightSalmon, true);

    if (debugProjection) dProjectionInfo.show();

  }
}

void steadyLatticeTest() {
  warningMsg = "";
  assert gSteadyLattice != null;
  if (gShowSteadyLattice) gSteadyLattice.show(null);

  if (gNavigateSteadyLattice) {  // use user-specified cube-range
    gRanges = gSteadyLattice.getCubeRange(gCubeCenter, gCubeHalfLength);
  } else {  // use full cube-range
    gRanges = gSteadyLattice.getFullRange();
  }

  ArrayList<Ball> balls = new ArrayList<Ball>();
  ArrayList<Edge> beams = new ArrayList<Edge>();
  gSteadyLattice.traverse(gRanges, balls, beams);

  if (balls.size() == 0 || beams.size() == 0) {
    warningMsg += "The selected subset of lattice is empty; ";
    return;
  }

  gLattice = new Lattice(balls, beams);
  latticeTest();
}