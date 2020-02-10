/*******************************************************************************
 * Helper functions to generate input (samples on a circle, points in 3D, etc.).
 * Initializations for the scene.
 ******************************************************************************/


/* Input methods for differnt objects. */
int inputMethodPointSet = 0;  // 0: read from file, 1: generate randomly, 2: from ring set
int inputMethodRingSet = 1;  // 0: read from file, 1: generate randomly
int inputMethodHub = 1;  // 0: read from file, 1: generate randomly
int inputMethodEdgeCircle = 1;  // 0: read from file, 1: generate randomly
int inputMethodTriangleMesh = 1;  // 0: read from file, 1: nothing happens
int inputMethodLattice = 0;  // 0: read from file, 1: generate manually
int inputMethodGap = 1;  // 0: read from file, 1: nothing happens
int inputMethodSteadyLattice = 1;  // 0: nothing happens, 1: generate manually

/* File paths for different objects. */
// String gPointSetPath = "data/point_set/ps_convex_hull_circles";
String gPointSetPath = "data/point_set/ps_hub_tessellation";
String gRingSetPath = "data/ring_set/rs_convex_hull_samples_on_circles";
String gHubPath = "data/hub/hub_9";
String gEdgeCirclePath = "data/edge_circle/ec_0";
String gTriangleMeshPath = "data/triangle_mesh/tm_0";
String gLatticePath = "data/lattice/lattice";
String gGapPath = "data/gap/gap_10";
String gSteadyLatticePath = "data/lattice/lattice";
String gCameraPath = "data/camera/cam_hub_tessellation";

/* File paths for results. */
String gSupPlaneTestsStatFile = "results/stats_supporting_plane_three_circles/stat_2";


/* Global variables related to gPoints. */
pts gPoints;  // the global points

/* Global variables related to gRingSet. */
float gAttenuationMin = 0.05;
float gAttenuationDelta = 0.05;
float gAttenuation = 1.0;
int gNumRings = 5;
int gNumPointsPerRing = 4;
RingSet gRingSet;  // the global ring set

/* Global variables related to gHub. */
float gRadiusInnerBall = 20;  // radius of the inner ball
int gNumNeighbors = 4;  // number of outer balls
Hub gHub;  // the global hub

/* Global variables related to gTriangleMesh. */
TriangleMesh gTriangleMesh;  // the global triangle mesh
TriangleMesh gBeamMesh;  // the global triangle mesh for lifted beams
TriangleMesh gGapMesh;  // the global triangle mesh for gaps

/* Global variables related to gTriQuadMesh. */
BorderedTriQuadMesh gTriQuadMesh;

/* Global variables related to gEdgeCircle. */
EdgeCircle gEdgeCircle;  // the global edge circle

/* Global variables related to gLattice. */
Lattice gLattice;

/* Global variables related to gGap. */
ConvexGap gGap;

/* Global variables related to gSteadyLattice. */
SteadyLattice gSteadyLattice;
MinMaxI[] gRanges = new MinMaxI[3];  // [iMin, iMax], [jMin, jMax], [kMin, kMax]
Idx3 gCubeCenter = new Idx3(0, 0, 0);
int gCubeHalfLength = 1;

/*
 * Choose which steady lattice model to use in Test 25.
 * 0: a regular lattice
 * 1: a rotated regular lattice
 * 2: a swirl regular lattice
 * 3: a swirl-rotation-translation lattice
 * ...
 * 10: the model in the convex hull paper (SPM 2019 submission 20)
 */
int gSteadyOption = 3;



/*
 * Randomly generate a point on the sphere with center c and radius r.
 */
pt generateOnePointOnSphere(pt c, float r) {
  float alpha = random(-HALF_PI, HALF_PI);
  float beta = random(0, TWO_PI);
  float dx = r * cos(alpha) * cos(beta);
  float dy = r * cos(alpha) * sin(beta);
  float dz = r * sin(alpha);
  return new pt(c.x + dx, c.y + dy, c.z + dz);
}

/*
 * Randomly generate a point inside the sphere with center c and radius r.
 */
pt generateOnePointInsideSphere(pt c, float r) {
  pt p = new pt();
  float r2 = r * r;
  int k = 100, i = 0;
  while (i < k) {
    p.x = random(c.x - r, c.x + r);
    p.y = random(c.y - r, c.y + r);
    p.z = random(c.z - r, c.z + r);
    if (n2(V(c, p)) <= r2) return p;
    i++;
  }
  return p;
}

/*
 * Randomly generate n points on the sphere with center c and radius r.
 */
pt[] generatePointsOnSphere(pt c, float r, int n) {
  pt[] points = new pt[n];
  for (int i = 0; i < n; ++i) {
    while (true) {
      pt p = generateOnePointOnSphere(c, r);
      boolean bad = false;
      for (int j = 0; j < i; ++j) {
        if (isNonPositive(d(points[j], p))) {
          bad = true;
          break;
        }
      }
      if (!bad) {
        points[i] = p;
        break;
      }
    }
  }
  return points;
}

/*
 * Randomly generate n points on the sphere with center c and radius r. Clear ps
 * and store the n points in it.
 */
void generatePointsOnSphere(pts ps, pt c, float r, int n) {
  ps.empty();
  pt[] points = generatePointsOnSphere(c, r, n);
  for (int i = 0; i < n; ++i) {
    ps.addPt(points[i]);
  }
}

/*
 * Sample points uniformly on a circle.
 */
private pt[] generatePointsForOneCircle(pt capCenter, float circleRadius,
                                        pt sphereCenter, float sphereRadius,
                                        vec xAxis, int n, pt circleCenter) {
  pt[] points = new pt[n];
  vec normal = U(sphereCenter, capCenter);
  vec yAxis = N(normal, xAxis);  // the order matters

  /* Only consider the circle that corresponds to a cone with its half-angle in [0, PI/2]. */
  float d = sqrt(sphereRadius * sphereRadius - circleRadius * circleRadius);
  circleCenter.set(P(sphereCenter, d, normal));

  float da = TWO_PI / n;
  float a = 0;
  for (int i = 0; i < n; ++i) {
    points[i] = P(circleCenter, circleRadius, R(xAxis, a, xAxis, yAxis));
    a += da;
  }
  return points;
}

/*
 * Sample np points uniformly on each of the nc circles which have the same
 * radius.
 * contact = the intersection between the sphere and the outward axis of a disk
 * center = the center of a disk (under the sphere)
 */
pt[][] generatePointsForCircles(pt[] capCenters, float circleRadius,
                                pt sphereCenter, float sphereRadius,
                                vec[] xAxes, int nc, int np, pt[] circleCenters) {
  pt[][] points = new pt[nc][np];
  for (int i = 0; i < nc; ++i) {
    circleCenters[i] = new pt();
    points[i] = generatePointsForOneCircle(capCenters[i], circleRadius,
                                           sphereCenter, sphereRadius, xAxes[i],
                                           np, circleCenters[i]);
  }
  return points;
}

/*
 * Sample np points uniformly on each of the nc circles which may have different
 * radii.
 */
pt[][] generatePointsForCircles(pt[] capCenters, float[] circleRadii,
                                pt sphereCenter, float sphereRadius,
                                vec[] xAxes, int nc, int np, pt[] circleCenters) {
  pt[][] points = new pt[nc][np];
  for (int i = 0; i < nc; ++i) {
    circleCenters[i] = new pt();
    points[i] = generatePointsForOneCircle(capCenters[i], circleRadii[i],
                                           sphereCenter, sphereRadius,
                                           xAxes[i], np, circleCenters[i]);
  }
  return points;
}

/*
 * Randomly generate n cap centers and n radii for corresponding disks such that
 * the n caps/disks are disjoint.
 */
pt[] generateCapCentersAndCircleRadii(pt sphereCenter, float sphereRadius, int n,
                                      float[] circleRadii, vec[] normals) {
  pt[] capCenters = new pt[n];
  float[] alphas = new float[n];
  for (int i = 0; i < n; ++i) {
    while (true) {
      pt p = generateOnePointOnSphere(sphereCenter, sphereRadius);
      vec normal = U(sphereCenter, p);
      float r = random(1, sphereRadius * 0.9);
      float alpha = asinClamp(r / sphereRadius);  // [0, PI/2]
      boolean isValid = true;
      for (int j = 0; j < i; ++j) {
        float theta = acosClamp(dot(normal, normals[j]));  // [0, PI]
        if (theta <= alpha + alphas[j]) {
          isValid = false;
          break;
        }
      }
      if (isValid) {
        capCenters[i] = p;
        normals[i] = normal;
        alphas[i] = alpha;
        circleRadii[i] = r;
        break;
      }
    }
  }
  return capCenters;
}

/*
 * Randomly generate n cap centers on a sphere such that the n corresponding
 * disks, each with the same radius rMax, are disjoint.
 */
pt[] generateCapCenters(pt sphereCenter, float sphereRadius, int n, float rMax,
                        vec[] normals) {
  assert rMax > 0 && rMax < sphereRadius;
  pt[] capCenters = new pt[n];
  float alpha = 2 * asinClamp(rMax/sphereRadius);
  for (int i = 0; i < n; ++i) {
    while (true) {
      pt p = generateOnePointOnSphere(sphereCenter, sphereRadius);
      vec normal = U(sphereCenter, p);
      boolean isValid = true;
      for (int j = 0; j < i; ++j) {
        float theta = acosClamp(dot(normal, normals[j]));  // [0, PI]
        if (theta <= alpha) {
          isValid = false;
          break;
        }
      }
      if (isValid) {
        capCenters[i] = p;
        normals[i] = normal;
        break;
      }
    }
  }
  return capCenters;
}

/*
 * Randomly generate n balls with their centers on a sphere.
 */
Ball[] generateBallsOnSphere(pt c, float r, float rMax, int n) {
  Ball[] balls = new Ball[n];
  assert rMax >= 4.0;
  for (int i = 0; i < n; ++i) {
    while (true) {
      pt px = generateOnePointOnSphere(c, r);
      float rx = random(4, rMax);
      boolean isValid = true;
      for (int j = 0; j < i; ++j) {
        if (d(px, balls[j].c) < rx + balls[j].r) {
          isValid = false;
          break;
        }
      }
      if (isValid) {
        Ball ball = new Ball(px, rx);
        balls[i] = ball;
        break;
      }
    }
  }
  return balls;
}

/*
 * Randomly generate a hub.
 */
Hub generateHub(pt c, float r0, float r1, int n) {
  Ball ball = new Ball(c, r0);
  Ball[] neighbors = generateBallsOnSphere(c, r1, r1 - 1.3 * r0, n);
  return new Hub(ball, neighbors, n);
}

/*
 * Manully design a steady lattice.
 */
SteadyLattice designSteadyLattice() {
  SteadyLattice lattice = new SteadyLattice();

  switch (gSteadyOption) {
    case 0:   // a regular lattice
      {
        int size = 5;  // u, v, w
        float fsize = (float)size;
        lattice.setG(MakeScaling(50.0 / fsize, P(0, 0, 0)));
        lattice.setU(size, MakeTranslation(V(1, 0, 0)));
        lattice.setV(size, MakeTranslation(V(0, 1, 0)));
        lattice.setW(size, MakeTranslation(V(0, 0, 1)));
        int a = lattice.addJoint(P(0,0,0), 0.2);
        lattice.addBeam(a, a, I(1, 0, 0));
        lattice.addBeam(a, a, I(0, 1, 0));
        lattice.addBeam(a, a, I(0, 0, 1));
      }
      break;
    case 1:  // a rotated regular lattice
      {
        int size = 5;
        float fsize = (float)size;
        lattice.setG(MakeScaling(50.0 / fsize, P(0, 0, 0)));
        lattice.setU(size, MakeTranslation(V(1, 0, 0)));
        lattice.setV(size, MakeRotation(HALF_PI / fsize, V(0, 0, 1), P(fsize*3, 0, fsize/2)));
        lattice.setW(size, MakeTranslation(V(0, 0, 2)));
        int a = lattice.addJoint(P(0, 0, 0), 0.2);
        lattice.addBeam(a, a, I(1, 0, 0));
        lattice.addBeam(a, a, I(0, 1, 0));
        lattice.addBeam(a, a, I(0, 0, 1));
      }
      break;
    case 2:  // a swirl regular lattice
      {
        int size = 6;  // default: 9
        float fsize = (float)size;
        lattice.setG(MakeScaling(50.0 / fsize, P(0, 0, 0)));
        lattice.setU(size, MakeTranslation(V(1, 0, 0)));
        lattice.setV(size, MakeSwirl(HALF_PI / fsize, V(0, 0, 1), P(fsize*3, 0, fsize/2), 0.95));
        lattice.setW(size, MakeTranslation(V(0, 0, 2)));
        int a = lattice.addJoint(P(0, 0, 0), 0.2);
        lattice.addBeam(a, a, I(1, 0, 0));
        lattice.addBeam(a, a, I(0, 1, 0));
        lattice.addBeam(a, a, I(0, 0, 1));
      }
      break;
    case 3:  // a swirl-rotation-translation lattice
      {
        int size = 7;  // default: 7
        float fsize = (float)size;
        lattice.setG(MakeScaling(50.0 / fsize, P(0, 0, 0)));
        lattice.setU(size, MakeTranslation(V(1, 0, 0)));
        lattice.setV(size, MakeRotation(PI/6/size, V(0, 1, 0), P(fsize*3, 0, fsize/2)));
        lattice.setW(size, MakeSwirl(PI/6/size, V(0, 1, 0), P(fsize/2, fsize*2, fsize/2), pow(.994, 99.0/fsize)));
        int a = lattice.addJoint(P(0, 0, 0), 0.2);
        lattice.addBeam(a, a, I(1, 0, 0));
        lattice.addBeam(a, a, I(0, 1, 0));
        lattice.addBeam(a, a, I(0, 0, 1));
      }
      break;

    case 10:  // the model in the convex hull paper (SPM 2019 submission 20)
      {
        lattice.setG(MakeScaling(0.25, P(0, 0, 0)));
        lattice.setU(1, new SwirlTransform());
        lattice.setV(15, MakeRotation(TWO_PI/16, V(1, 0, 0), P(0, 0, 0)));
        lattice.setW(8, MakeSwirl(-TWO_PI/32, V(1, 0, 0), P(90, 0, 0), 0.9));
        // lattice.restrictInBounds(true, false, true);  // not used in my code
        int a = lattice.addJoint(P(0, 20, 0), 1.005);
        lattice.addBeam(a, a, I(0, 1, 0));
        lattice.addBeam(a, a, I(0, 0, 1));
      }
      break;
  }
  return lattice;
}

/*
 * Manually design a lattice.
 */
Lattice designLattice() {
  Ball[] balls = new Ball[7];
  int[][] beams = new int[6][2];

  balls[0] = new Ball(P(0, 0, 0), 25);
  balls[1] = new Ball(P(-70, 50, 0), 20);
  balls[2] = new Ball(P(50, 50, 0), 30);
  balls[3] = new Ball(P(0, -90, 0), 20);
  balls[4] = new Ball(P(-120, 170, 0), 10);
  balls[5] = new Ball(P(140, 190, 0), 10);
  balls[6] = new Ball(P(20, -210, 0), 35);

  beams[0][0] = 0;
  beams[0][1] = 1;
  beams[1][0] = 0;
  beams[1][1] = 2;
  beams[2][0] = 0;
  beams[2][1] = 3;

  beams[3][0] = 1;
  beams[3][1] = 4;
  beams[4][0] = 2;
  beams[4][1] = 5;
  beams[5][0] = 3;
  beams[5][1] = 6;

  return new Lattice(balls, beams);
}

/*
 * Initialize all objects in the scene for all test cases.
 */
void initScene() {
  switch (inputMethodRingSet) {
    case 0:  // read from file
      gRingSet = new RingSet();
      gRingSet.load(gRingSetPath);
      break;
    case 1:  // generate randomly
      gRingSet = new RingSet(gSphereCenter, gSphereRadius, gNumRings, gNumPointsPerRing);
      break;
    default:
      println("Please use a valid input method for ring set");
  }
  if (gRingSet == null || !gRingSet.isValid()) {
    println("No valid ring set!");
  }

  gPoints = new pts();
  gPoints.declare();
  switch (inputMethodPointSet) {
    case 0:  // read from file
      gPoints.load(gPointSetPath);
      break;
    case 1:  // generate randomly
      generatePointsOnSphere(gPoints, gSphereCenter, gSphereRadius, 10);
      break;
    case 2:  // from ring set
      if (gRingSet != null) gPoints = gRingSet.toPointSet();
      break;
    default:
      println("Please use a valid input method for point set");
  }

  switch (inputMethodHub) {
    case 0:  // read from file
      gHub = new Hub();
      gHub.load(gHubPath);
      break;
    case 1:  // generate randomly
      gHub = generateHub(gSphereCenter, gRadiusInnerBall, gSphereRadius, gNumNeighbors);
      break;
    default:
      println("Please use a valid input method for hub");
  }

  gEdgeCircle = new EdgeCircle();
  switch (inputMethodEdgeCircle) {
    case 0:
      gEdgeCircle.load(gEdgeCirclePath);
      break;
    case 1:
      gEdgeCircle.init();
      break;
    default:
      println("Please use a valid input method for edge-circle");
  }

  gTriangleMesh = new TriangleMesh();
  switch (inputMethodTriangleMesh) {
    case 0:
      gTriangleMesh.load(gTriangleMeshPath);
      break;
    default:
      println("Please use a valid input method for triangle mesh");
  }

  gLattice = new Lattice();
  switch (inputMethodLattice) {
    case 0:
      gLattice.load(gLatticePath);
      break;
    case 1:
      gLattice = designLattice();
      break;
    default:
      println("Please use a valid input method for lattice");
  }

  gGap = new ConvexGap();
  switch (inputMethodGap) {
    case 0:
      gGap.load(gGapPath);
      // gGap.translate(V(gGap.center()).rev());
      // gGap.scale(10.0);
      break;
    default:
      println("Please use a valid input method for gap");
  }

  gSteadyLattice = new SteadyLattice();
  switch (inputMethodSteadyLattice) {
    case 1:
      gSteadyLattice = designSteadyLattice();
      break;
    default:
      println("Please use a valid input method for steady lattice");
  }

  Idx3 uvw = gSteadyLattice.repetitionCounts();
  assert uvw != null;
  gRanges[0] = new MinMaxI(0, uvw.i);
  gRanges[1] = new MinMaxI(0, uvw.j);
  gRanges[2] = new MinMaxI(0, uvw.k);
}
