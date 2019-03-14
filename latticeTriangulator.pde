import processing.pdf.*;


/*
 * 0: one convex-hull test
 * 1: one convex-hull-with-holes test
 * 2: one ring-set-triangulation test (2 methods)
 * 3: one subdivision test
 * 4: one exact-convex-hull-of-circles test (non-interactive)
 * 5: one pivot-plane-around-line-until-hit-circle test
 * 6: one exact-convex-hull-of-edge-circle test
 * 7: one exact-convex-hull-of-two-circles test
 * 8: one supporting-plane-of-three-circles-iter test (iterative method with 3 different initializations)
 * 9: one exact-convex-hull-of-three-circles test
 * 10: one three-ring-triangle test
 * 11: one construct-elliptic-cone test
 * 12: one interactive-naive-exact-convex-hull test
 * 13: one interactive-incremental-exact-convex-hull test
 * 14: one corridor test
 * 15: one mesh-from-exact-convex-hull test
 * 16: one interactive-hub test
 * 17: one supporting-plane-of-three-circles-special-case test
 * 18: one geodesic-distance test
 * 19: one supporting-plane-of-three-circles test
 * 20: one elliptic-cone test
 * 21: one interactive-hub-to-mesh test
 * 22: one hub-to-mesh test
 * 23: one lattice-to-mesh test
 * 24: one convex-gap test
 * 25: one hub test (non-interactive)
 * ...
 * 100: many convex-hull tests
 * 101: many ring-set-triangulation tests
 * 102: many three-ring-triangle tests
 * 103: many supporting-plane-of-three-circles tests
 * 104: many exact-convex-hull-all-circles tests
 * 105: many incremental-convex-hull tests
 * ...
 * 200: one circle-plane-intersection test
 * 201: one hub-line-intersection test
 * 202: one round-cone-distance test
 * 203: one intersection-between-two-planes test
 * 204: one intersection-between-two-spheres test
 */
int test = 23;

/* Input methods for differnt objects. */
int inputMethodPointSet = 0;  // 0: read from file, 1: generate randomly, 2: from ring set
int inputMethodRingSet = 0;  // 0: read from file, 1: generate randomly
int inputMethodHub = 0;  // 0: read from file, 1: generate randomly
int inputMethodEdgeCircle = 1;  // 0: read from file, 1: generate randomly
int inputMethodTriangleMesh = 1;  // 0: read from file, 1: nothing happens
int inputMethodLattice = 0;  // 0: read from file
int inputMethodGap = 0;  // 0: read from file

/* File paths for different objects. */
String gPointSetPath = "data/point_set/ps_arcs_16";
String gRingSetPath = "data/ring_set/rs_3";
String gHubPath = "data/hub/hub_3";
String gEdgeCirclePath = "data/edge_circle/ec_0";
String gTriangleMeshPath = "data/triangle_mesh/tm_0";
String gLatticePath = "data/lattice/lattice_4";
String gGapPath = "data/gap/gap_5";

boolean showSphere = true;
boolean showCenterOfSphere = true;
boolean showArcSet = true;
boolean showAuxPlane = false;

boolean generateCH = true;
boolean regenerateCH = true;  // for ring set, test shrink/grow
boolean showPointSet = true;

boolean animating = true;
boolean tracking = false;
boolean center = true;
boolean showFrame = false;
boolean showFocus = false;
boolean snappingPDF = false;

pt gFocus = new pt(0.0, 0.0, 0.0);  // focus point: the camera is looking at it (moved when 'f or 'F' are pressed)
pt Pick = new pt(0.0, 0.0, 0.0);

float gSphereRadius = 100;  // default: 100
pt gSphereCenter = new pt(0.0, 0.0, 0.0);

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

/* Global variables related to gEdgeCircle. */
EdgeCircle gEdgeCircle;  // the global edge circle

/* Global variables related to gLattice. */
Lattice gLattice;

ConvexGap gGap;

Camera gCamera = new Camera(500, -0.06 * TWO_PI, -0.04 * TWO_PI);

int gNumTriangles = -1;
float timeTM = 0.0;
float timeSD = 0.0;

void setup() {
  face0 = loadImage("data/Yaohong.jpg");  // load Yaohong's image
  face1 = loadImage("data/Jarek.jpg");  // load Jarek's image
  textureMode(NORMAL);
  size(900, 900, P3D); // P3D means that we will do 3D graphics
  noSmooth();  // LEAVE HERE FOR 3D PICK TO WORK!!!

  if (test >= 100 && test < 200) {
    noLoop();
    debugCH = false;
    debug3RT = false;
    debug2RT = false;
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

  if (regenerateCH == false) {
    gAttenuation = gAttenuationMin;
    gRingSet.generatePoints(gAttenuation);  // shrink all rings
    if (test == 1) {
      debugCH = false;
      gRingSet.generateTriangleMesh(0);  // generate a triangle mesh and store it
    }
    if (test == 2) {
      debug3RT = false;
      debug2RT = false;
      gRingSet.generateTriangleMesh(1);  // generate a triangle mesh and store it
    }
  }
}

void draw() {
  background(255);

  pushMatrix();  // to ensure that we can restore the standard view before writing on the canvas
  if (snappingPDF) beginRecord(PDF, "PDFimages/P" + nf(pictureCounter++, 3) + ".pdf");

  /* Set perspective. */
  float cameraZ = (height / 2.0) / tan(gCamera.fov / 2.0);
  camera(0, 0, cameraZ, 0, 0, 0, 0, 1, 0);  // sets a standard perspective
  perspective(gCamera.fov, float(width)/float(height), gCamera.zNear, gCamera.zFar);

  /* Set view. */
  translate(gCamera.dx, gCamera.dy, gCamera.dz);
  lights();  // turns on view-dependent lighting
  rotateX(gCamera.rx);  // rotates the scene around the focus
  rotateY(gCamera.ry);  // rotates the scene around the focus
  rotateX(HALF_PI);  // rotates frame around X to make X and Y basis vectors parallel to the floor
  if (center) translate(-gFocus.x, -gFocus.y, -gFocus.z);

  computeProjectedVectors(); // computes screen projections I, J, K of basis vectors (see bottom of pv3D): used for dragging in viewer's frame
  if (showFrame) showFrame(150); // X-red, Y-green, Z-blue arrows
  if (showFocus) {
    fill(red);
    showBall(gFocus, 10);
  }

  if (showCenterOfSphere) {
    noStroke();
    fill(black);
    showBall(gSphereCenter, 3); // show center of the sphere
  }

  switch (test) {
    case 0:
      convexHullTest();
      break;
    case 1:
      convexHullWithHolesTest();
      break;
    case 2:
      ringSetTriangulationTest();
      break;
    case 3:
      subdivisionTest();
      break;
    case 4:
      incCHTest();
      break;
    case 5:
      pivotPlaneAroundLineHitCircleTest();
      break;
    case 6:
      exactCHEdgeCircleTest();
      break;
    case 7:
      exactCHTwoCirclesTest();
      break;
    case 8:
      supPlaneThreeCirclesIterTest();
      break;
    case 9:
      exactCHThreeCirclesTest();
      break;
    case 10:
      threeRingTriangleTest();
      break;
    case 11:
      constructEllipticConeTest();
      break;
    case 12:
      interactiveNaiveCHTest();
      break;
    case 13:
      interactiveIncCHTest();
      break;
    case 14:
      corridorTest();
      break;
    case 15:
      meshFromExactCHTest();
      break;
    case 16:
      interactiveHubTest();
      break;
    case 17:
      supPlaneThreeCirclesSpecialTest();
      break;
    case 18:
      geodesicDistanceTest();
      break;
    case 19:
      supPlaneThreeCirclesTest();
      break;
    case 20:
      ellipticConeTest();
      break;
    case 21:
      interactiveHubToMeshTest();
      break;
    case 22:
      hubToMeshTest();
      break;
    case 23:
      latticeTest();
      break;
    case 24:
      convexGapTest();
      break;
    case 25:
      hubTest();
      break;

    case 100:  // many convex-hull tests
      convexHullTests(numTests, numPointsPerTest);
      break;
    case 101:  // many ring-set-triangulation tests
      ringSetTriangulationTests(numTests, gNumRings, gNumPointsPerRing, gAttenuation);
      break;
    case 102:  // many three-ring-triangle tests
      threeRingTriangleTests(numTests, gNumPointsPerRing, gAttenuation);
      break;
    case 103:  // many supporting-plane-of-three-circles tests
      supPlaneThreeCirclesTests(numTests, gAttenuation);
      break;
    case 104:
      exactCHAllCirclesTests(numTests, gNumRings);
      break;
    case 105:
      incrementalConvexHullTests();
      break;

    case 200:
      circlePlaneIntersectionTest();
      break;
    case 201:
      hubLineIntersectionTest();
      break;
    case 202:
      roundConeDistTest();
      break;
    case 203:
      intersectionTwoPlanesTest();
      break;
    case 204:
      intersectionTwoSpheresTest();
      break;
    default:
      println("Please enter a correct test number");
      exit();
  }

  // Show the big yellow sphere
  if (showSphere) {
    fill(yellow, 100);
    noStroke();
    show(gSphereCenter, gSphereRadius);  // this should be before pick
    Pick = pick(mouseX, mouseY);
  }

  if (picking) {
    gPoints.setPickToIndexOfVertexClosestTo(Pick); // id of vertex of P with closest screen projection to mouse (us in keyPressed 'x'...
    picking = false;
  }

  if (showSphere && showArcSet) {
    int nv = gPoints.nv - gPoints.nv % 2;
    showArcs(gPoints.G, nv, gSphereCenter, gSphereRadius, 4, 3, 4);
    fill(red, 100);
    gPoints.showPicked(6); // show currently picked vertex of P
    fill(orange, 100);
    show(pick(mouseX, mouseY), 5);  // show mouse position on the sphere
  }

  if (exitDraw) noLoop();

  if (snappingPDF) {
    endRecord();
    snappingPDF = false;
  }
  popMatrix(); // done with 3D drawing, restore front view for writing text

  fill(black);
  if (scribeText) {
    displayHeader();
    displayDebugText();
  }

  // show menu at bottom, only if not filming
  if (scribeText && !filming) displayFooter();
  if (filming && (animating || change)) {
    saveFrame("FRAMES/F" + nf(frameCounter++, 4) + ".tif");  // save next frame
  }
  change = false; // to avoid capturing frames when nothing happens (change is set uppn action)
}
