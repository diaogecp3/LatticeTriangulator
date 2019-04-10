/*******************************************************************************
 * Helper functions to generate input (samples on a circle, points in 3D, etc.).
 * Initializations for the scene.
 ******************************************************************************/


/* Input methods for differnt objects. */
int inputMethodPointSet = 0;  // 0: read from file, 1: generate randomly, 2: from ring set
int inputMethodRingSet = 1;  // 0: read from file, 1: generate randomly
int inputMethodHub = 0;  // 0: read from file, 1: generate randomly
int inputMethodEdgeCircle = 1;  // 0: read from file, 1: generate randomly
int inputMethodTriangleMesh = 1;  // 0: read from file, 1: nothing happens
int inputMethodLattice = 1;  // 0: read from file, 1: generate manually
int inputMethodGap = 0;  // 0: read from file
int inputMethodSteadyLattice = 1;  // 0: read from file, 1: generate manually

/* File paths for different objects. */
String gPointSetPath = "data/point_set/ps_stereo_3";
String gRingSetPath = "data/ring_set/rs_0";
String gHubPath = "data/hub/hub_5";
String gEdgeCirclePath = "data/edge_circle/ec_0";
String gTriangleMeshPath = "data/triangle_mesh/tm_0";
String gLatticePath = "data/lattice/lattice_4";
String gGapPath = "data/gap/gap_5";
String gSteadyLatticePath = "data/lattice/steady_lattice_0";

/* Global variables related to gPoints. */
pts gPoints;  // the global points

/* Global variables related to gRingSet. */
float gAttenuationMin = 0.05;
float gAttenuationDelta = 0.05;
float gAttenuation = 1.0;
int gNumRings = 5;
int gNumPointsPerRing = 10;
RingSet gRingSet;  // the global ring set

/* Global variables related to gHub. */
float gRadiusInnerBall = 20;  // radius of the inner ball
int gNumNeighbors = 4;  // number of outer balls
Hub gHub;  // the global hub

/* Global variables related to gTriangleMesh. */
TriangleMesh gTriangleMesh;  // the global triangle mesh
TriangleMesh gBeamMesh;  // the global triangle mesh for lifted beams
TriangleMesh gGapMesh;  // the global triangle mesh for gaps

/* Global variables related to gEdgeCircle. */
EdgeCircle gEdgeCircle;  // the global edge circle

/* Global variables related to gLattice. */
Lattice gLattice;

/* Global variables related to gGap. */
ConvexGap gGap;

/* Global variables related to gSteadyLattice. */
SteadyLattice gSteadyLattice;
MinMaxI[] gRanges = new MinMaxI[3];  // [iMin, iMax], [jMin, jMax], [kMin, kMax]
Idx3 gCubeCenter = new Idx3(2, 2, 2);
int gCubeHalfLength = 1;


pt generateOnePointOnSphere(pt c, float r) {
  float alpha = random(-HALF_PI, HALF_PI);
  float beta = random(0, TWO_PI);
  float dx = r * cos(alpha) * cos(beta);
  float dy = r * cos(alpha) * sin(beta);
  float dz = r * sin(alpha);
  return new pt(c.x + dx, c.y + dy, c.z + dz);
}

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

void generatePointsOnSphere(pts P, pt c, float r, int n) {
  P.empty();
  pt[] points = generatePointsOnSphere(c, r, n);
  for (int i = 0; i < n; ++i) {
    P.addPt(points[i]);
  }
}


pt[] generateContactsOnSphere(pt C, float R, int nc, float rMax) {
  assert rMax < R;
  pt[] points = new pt[nc];
  Disk[] disks = new Disk[nc];
  float distance = sqrt(R * R - rMax * rMax);
  for (int i = 0; i < nc; ++i) {
    while (true) {
      pt p = generateOnePointOnSphere(C, R);
      vec normal = U(C, p);
      pt c = P(C, distance, normal);  // push p towards C
      Disk disk = new Disk(c, normal, rMax);
      boolean bad = false;
      for (int j = 0; j < i; ++j) {
        if (!emptyIntersectionTwoDisks(disk, disks[j])) {  // too expensive
          bad = true;
          break;
        }
      }
      if (!bad) {
        points[i] = p;
        disks[i] = disk;
        break;
      }
    }  // end while
  }
  return points;
}


// TODO: need to revisit
pt[] generatePointsForOneCircle(pt p,                                // in
                                float r,                             // in
                                pt C,                                // in
                                float R,                             // in
                                vec initDir,                         // in
                                int np,                              // in
                                pt center) {                         // out
  pt[] points = new pt[np];
  vec normal = U(C, p);
  vec v = I = initDir;
  vec J = N(normal, I);  // the order matters! make sure cross(I, J) same as normal
  center.set(P(C, sqrt(R * R - r * r), normal));  // correct?
  float da = TWO_PI / np;
  float a = 0;
  for (int i = 0; i < np; ++i) {
    points[i] = P(center, r, R(v, a, I, J));
    a += da;
  }
  return points;
}


// a contact = the intersection point between the medial axis of a tube and the sphere
// a center = the center of a circle
pt[][] generatePointsForCircles(pt[] contacts,                       // in
                                float r,                             // in
                                pt C,                                // in
                                float R,                             // in
                                vec[] initDirs,                      // in
                                int nc,                              // in
                                int np,                              // in
                                pt[] centers) {                      // out
  pt[][] points = new pt[nc][np];
  for (int i = 0; i < nc; ++i) {
    centers[i] = new pt();
    points[i] = generatePointsForOneCircle(contacts[i], r, C, R, initDirs[i],
                                           np, centers[i]);
  }
  return points;
}

pt[][] generatePointsForCircles(pt[] contacts,                       // in
                                float[] radii,                       // in
                                pt C,                                // in
                                float R,                             // in
                                vec[] initDirs,                      // in
                                int nc,                              // in
                                int np,                              // in
                                pt[] centers) {                      // out
  pt[][] points = new pt[nc][np];
  for (int i = 0; i < nc; ++i) {
    centers[i] = new pt();
    points[i] = generatePointsForOneCircle(contacts[i], radii[i], C, R,
                                           initDirs[i], np, centers[i]);
  }
  return points;
}


pt[] generateContactsAndRadii(pt center,                             // in
                              float radius,                          // in
                              int n,                                 // in
                              float[] radii,                         // out
                              vec[] normals) {                       // out
  pt[] contacts = new pt[n];
  float[] alphas = new float[n];
  for (int i = 0; i < n; ++i) {
    while (true) {
      pt p = generateOnePointOnSphere(center, radius);
      vec normal = U(center, p);
      float r = random(1, radius * 0.9);
      float alpha = asinClamp(r/radius);  // [0, PI/2]
      boolean isValid = true;
      for (int j = 0; j < i; ++j) {
        float theta = acosClamp(dot(normal, normals[j]));  // [0, PI]
        if (theta <= alpha + alphas[j]) {
          isValid = false;
          break;
        }
      }
      if (isValid) {
        contacts[i] = p;
        normals[i] = normal;
        alphas[i] = alpha;
        radii[i] = r;
        break;
      }
    }
  }
  return contacts;
}


pt[] generateContacts(pt center,                                     // in
                      float radius,                                  // in
                      int n,                                         // in
                      float rMax,                                    // in
                      vec[] normals) {                               // out
  assert rMax > 0 && rMax < radius;
  pt[] contacts = new pt[n];
  float alpha = 2 * asinClamp(rMax/radius);
  for (int i = 0; i < n; ++i) {
    while (true) {
      pt p = generateOnePointOnSphere(center, radius);
      vec normal = U(center, p);
      boolean isValid = true;
      for (int j = 0; j < i; ++j) {
        float theta = acosClamp(dot(normal, normals[j]));  // [0, PI]
        if (theta <= alpha) {
          isValid = false;
          break;
        }
      }
      if (isValid) {
        contacts[i] = p;
        normals[i] = normal;
        break;
      }
    }
  }
  return contacts;
}

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

Hub generateHub(pt c, float r0, float r1, int n) {
  Ball ball = new Ball(c, r0);
  Ball[] neighbors = generateBallsOnSphere(c, r1, r1 - 1.3 * r0, n);
  return new Hub(ball, neighbors, n);
}

/* Manully design a steady lattice. */
SteadyLattice generateSteadyLattice() {
  int size = 9;
  float fsize = float(size);

  SteadyLattice lattice = new SteadyLattice();

  lattice.setG(MakeScaling(50.0/size, P(0,0,0)));

  lattice.setU(size, MakeTranslation(V( 1, 0, 0 )));
  lattice.setV(size, MakeTranslation(V( 0, 1, 0 )));
  lattice.setW(size, MakeTranslation(V( 0, 0, 1 )));

  // lattice.setU(size, MakeTranslation(V( 1, 0, 0 )));
  // lattice.setV(size, MakeRotation(PI/6/size, V( 0, 1, 0 ), P(fsize*3,0,fsize/2.0)));
  // lattice.setW(size, MakeSwirl(PI/6/size, V( 0, 1, 0 ), P(fsize/2.0,fsize*2,fsize/2.0), pow(.994, 99.0/size)));

  int a = lattice.addJoint(P(0,0,0), .2);

  lattice.addBeam(a, a, I(1,0,0));
  lattice.addBeam(a, a, I(0,1,0));
  lattice.addBeam(a, a, I(0,0,1));

  return lattice;
}

Lattice generateTetRobot() {
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

void initScene() {
  switch (inputMethodRingSet) {
    case 0:  // read from file
      gRingSet = new RingSet();
      gRingSet.load(gRingSetPath);
      {
        // gRingSet.translate(V(gRingSet.sphereCenter).rev());
      }
      break;
    case 1:  // generate randomly
      gRingSet = new RingSet(gSphereCenter, gSphereRadius, gNumRings, gNumPointsPerRing);
      gRingSet.init();
      break;
    default:
      println("Please use a valid input method for ring set");
      exit();
  }
  if (!gRingSet.isValid()) {  // check validity of ring set
    println("ring set not valid!");
    exit();
  }

  gPoints = new pts();
  gPoints.declare();  // some points on a sphere
  switch (inputMethodPointSet) {
    case 0:  // read from file
      gPoints.loadPts(gPointSetPath);
      break;
    case 1:  // generate randomly
      generatePointsOnSphere(gPoints, gSphereCenter, gSphereRadius, 10);
      break;
    case 2:  // from ring set
      gPoints = gRingSet.toPointSet();
      break;
    default:
      println("Please use a valid input method for point set");
      exit();
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
      exit();
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
      exit();
  }

  gTriangleMesh = new TriangleMesh();
  switch (inputMethodTriangleMesh) {
    case 0:
      gTriangleMesh.load(gTriangleMeshPath);
      break;
    case 1:
      break;
    default:
      println("Please use a valid input method for triangle mesh");
      exit();
  }

  gLattice = new Lattice();
  switch (inputMethodLattice) {
    case 0:
      gLattice.load(gLatticePath);
      break;
    case 1:
      gLattice = generateTetRobot();
      break;
    default:
      println("Please use a valid input method for lattice");
      exit();
  }

  gGap = new ConvexGap();
  switch (inputMethodGap) {
    case 0:
      gGap.load(gGapPath);
      {
        // pt c = gGap.center();
        // gGap.translate(V(c).rev());
        // gGap.scale(10.0);
        // println("before, gGap.points0.size() =", gGap.points0.size(), "gGap.points1.size() =", gGap.points1.size());
        // gGap.removeDuplicatePoints(0);
        // gGap.removeDuplicatePoints(1);
        // println("after, gGap.points0.size() =", gGap.points0.size(), "gGap.points1.size() =", gGap.points1.size());
      }
      break;
    default:
      // println("Please use a valid input method for gap");
      // exit();
  }

  gSteadyLattice = new SteadyLattice();
  switch (inputMethodSteadyLattice) {
    case 0:
      break;
    case 1:
      gSteadyLattice = generateSteadyLattice();
      break;
    default:
  }
  Idx3 uvw = gSteadyLattice.repetitionCounts();
  assert uvw != null;
  gRanges[0] = new MinMaxI(0, uvw.i);
  gRanges[1] = new MinMaxI(0, uvw.j);
  gRanges[2] = new MinMaxI(0, uvw.k);
}

