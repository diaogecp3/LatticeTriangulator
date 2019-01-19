import processing.pdf.*;


/*
 * 0: one convex-hull test
 * 1: one convex-hull-with-holes test
 * 2: one ring-set-triangulation test (2 methods)
 * 3: one subdivision test
 * 4: one hub test
 * 5: one pivot-plane-around-line-until-hit-circle test
 * 6: one exact-convex-hull-of-edge-circle test
 * 7: one exact-convex-hull-of-two-circles test
 * 8: one supporting-plane-of-three-circles-iter test (iterative method with 3 different initializations)
 * 9: one exact-convex-hull-of-three-circles test
 * 10: one three-ring-triangle test
 * 11: one triangle-mesh test
 * 12: one interactive-naive-exact-convex-hull test
 * 13: one interactive-incremental-exact-convex-hull test
 * 14: one corridor test
 * 15: one mesh-from-exact-convex-hull test
 * 16: one interactive-hub test
 * ...
 * 100: many convex-hull tests
 * 101: many ring-set-triangulation tests
 * 102: many three-ring-triangle tests
 * 103: many extreme-plane tests
 * ...
 * 200: one circle-plane-intersection test
 * 201: one hub-line-intersection test
 */
int test = 16;

float tan0 = 0, tan1 = 0;  // for debugging supporting triangle of 3 circles
boolean validRS = false;

int inputMethodPointSet = 0;  // 0: read from file, 1: generate randomly
int inputMethodRingSet = 0;  // 0: read from file, 1: generate randomly
int inputMethodHub = 1;  // 0: read from file, 1: generate randomly
int inputMethodEdgeCircle = 1;  // 0: read from file, 1: generate randomly
int inputMethodTriangleMesh = 1;  // 0: read from file, 1: nothing happens

boolean showSphere = true;
boolean showCenterOfSphere = true;
boolean showArcSet = true;
boolean showAuxPlane = false;
int approxMethodCorridor = 0;
boolean generateCH = false;
boolean regenerateCH = true;  // for ring set, test shrink/grow
boolean showPointSet = true;

float dz = 500;  // distance to camera. Manipulated with mouse wheel
float rx = -0.06 * TWO_PI, ry = -0.04 * TWO_PI;  // view angles manipulated when space pressed but not mouse
boolean animating = true;
boolean tracking = false;
boolean center = true;
boolean showFrame = false;
boolean snappingPDF = false;
boolean viewpoint = false;
float t = 0, s = 0;
pt Viewer = P();
pt F = P(0,0,0);  // focus point: the camera is looking at it (moved when 'f or 'F' are pressed)
pt Of = P(100,100,0), Ob = P(110,110,0);  // red point controlled by the user via mouseDrag: used for inserting vertices ...
pt Vf = P(0,0,0), Vb = P(0,0,0);
pt Pick=P();

float radiusOfSphere = 100;
pt centerOfSphere = new pt(0.0, 0.0, 0.0);

/* Global variables related to gPoints. */
pts gPoints;  // the global points

/* Global variables related to gRingSet. */
float rMax = 50;
float attenuationMin = 0.05;
float attenuationDelta = 0.05;
float attenuation = 1.0;
int numRings = 5;
int numPointsPerRing = 6;
RingSet gRingSet;  // the global ring set

/* Global variables related to gHub. */
float rInnerBall = 20;  // radius of the inner ball
float rSphereOfOuterBalls = radiusOfSphere;  // radius of the sphere where the centers of the outer balls lie
int nNeighbors = 4;  // number of outer balls
Hub gHub;  // the global hub

/* Global variables related to gTriangleMesh. */
TriangleMesh gTriangleMesh;  // the global triangle mesh

/* Global variables related to gEdgeCircle. */
EdgeCircle gEdgeCircle;  // the global edge circle

int numTriangles = -1;
float timeTM = 0.0;
float timeSD = 0.0;

void setup() {
  face0 = loadImage("data/Yaohong.jpg");  // load Yaohong's image
  face1 = loadImage("data/Jarek.jpg");  // load Jarek's image
  textureMode(NORMAL);
  size(900, 900, P3D); // P3D means that we will do 3D graphics
  noSmooth();  // LEAVE HERE FOR 3D PICK TO WORK!!!

  gPoints = new pts();
  gPoints.declare();  // some points on a sphere

  if (test >= 100 && test < 200) {
    noLoop();
    debugCH = false;
    debug3RT = false;
    debug2RT = false;
  }

  switch (inputMethodPointSet) {
    case 0:  // read from file
      gPoints.loadPts("data/point_set/ps_arcs_14");
      break;
    case 1:  // generate randomly
      generatePointsOnSphere(gPoints, centerOfSphere, radiusOfSphere, 10);
      break;
    default:
      println("Please use a valid input method for point set");
      exit();
  }

  switch (inputMethodRingSet) {
    case 0:  // read from file
      gRingSet = new RingSet(centerOfSphere, radiusOfSphere);
      gRingSet.load("data/rs_unnamed");
      break;
    case 1:  // generate randomly
      gRingSet = new RingSet(centerOfSphere, radiusOfSphere,
                       numRings, numPointsPerRing);
      gRingSet.init();
      break;
    default:
      println("Please use a valid input method for ring set");
      exit();
  }

  switch (inputMethodHub) {
    case 0:  // read from file
      gHub = new Hub();
      gHub.load("data/hub/hub_easy_3");
      break;
    case 1:  // generate randomly
      gHub = generateHub(centerOfSphere, rInnerBall, rSphereOfOuterBalls, nNeighbors);
      break;
    default:
      println("Please use a valid input method for hub");
      exit();
  }

  switch (inputMethodEdgeCircle) {
    case 0:
      gEdgeCircle = new EdgeCircle();
      gEdgeCircle.load("data/edge_circle/ec_0");
      break;
    case 1:
      gEdgeCircle = new EdgeCircle();
      gEdgeCircle.init();
      break;
    default:
      println("Please use a valid input method for hub");
      exit();
  }

  switch (inputMethodTriangleMesh) {
    case 0:
      gTriangleMesh = new TriangleMesh();
      gTriangleMesh.load("data/triangle_mesh/tm_0");
      break;
    default:
      gTriangleMesh = null;
  }

  // check validity of ring set
  if (!gRingSet.isValid()) {
    println("ring set not valid!");
    exit();
  }

  if (regenerateCH == false) {
    attenuation = attenuationMin;
    gRingSet.generatePoints(attenuation);  // shrink all rings
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

  // SET PERSPECTIVE
  float fov = PI / 3.0;
  float cameraZ = (height / 2.0) / tan(fov / 2.0);
  camera(0, 0, cameraZ, 0, 0, 0, 0, 1, 0);  // sets a standard perspective
  perspective(fov, 1.0, 1.0, 10000);

  // SET VIEW
  translate(0, 0, dz);  // puts origin of model at screen center and moves forward/away by dz
  lights();  // turns on view-dependent lighting
  rotateX(rx); rotateY(ry);  // rotates the model around the new origin (center of screen)
  rotateX(HALF_PI);  // rotates frame around X to make X and Y basis vectors parallel to the floor
  if (center) translate(-F.x, -F.y, -F.z);
  if (viewpoint) {Viewer = viewPoint(); viewpoint = false;} // sets Viewer to the current viewpoint when ',' is pressed
  computeProjectedVectors(); // computes screen projections I, J, K of basis vectors (see bottom of pv3D): used for dragging in viewer's frame
  if (showFrame) showFrame(150); // X-red, Y-green, Z-blue arrows

  if (showCenterOfSphere) {
    noStroke();
    fill(black);
    show(centerOfSphere, 3); // show center of the sphere
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
      hubTest();
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
      tangentPlaneThreeCirclesIterTest();
      break;
    case 9:
      exactCHThreeCirclesTest();
      break;
    case 10:
      threeRingTriangleTest();
      break;
    case 11:
      triangleMeshTest();
      break;
    case 12:
      exactCHNaiveTest();
      break;
    case 13:
      exactCHIncrementalTest();
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

    case 100:  // many convex-hull tests
      convexHullTests(numTests, numPointsPerTest);
      break;
    case 101:  // many ring-set-triangulation tests
      ringSetTriangulationTests(numTests, numRings, numPointsPerRing, attenuation);
      break;
    case 102:  // many three-ring-triangle tests
      threeRingTriangleTests(numTests, numPointsPerRing, attenuation);
      break;
    case 103:  // many extreme-plane(-of-three-circles) tests
      extremePlaneTests(numTests, attenuation);
      break;
    case 104:
      exactCHAllCirclesTests(numTests, numRings);
      break;

    case 200:
      circlePlaneIntersectionTest();
      break;
    case 201:
      hubLineIntersectionTest();
      break;
    default:
      println("Please enter a correct test number");
      exit();
  }

  // Show the big yellow sphere
  if (showSphere) {
    fill(yellow, 100);
    noStroke();
    show(centerOfSphere, radiusOfSphere);  // this should be before pick
    Pick = pick(mouseX, mouseY);
  }

  if (picking) {
    gPoints.setPickToIndexOfVertexClosestTo(Pick); // id of vertex of P with closest screen projection to mouse (us in keyPressed 'x'...
    picking = false;
  }

  /* Tests that should be done after picking. */
  switch (test) {
    case 16:
      break;
  }

  if (showSphere && showArcSet) {
    noStroke();
    int nv = gPoints.nv - gPoints.nv % 2;
    showArcs(gPoints.G, nv, centerOfSphere, radiusOfSphere, 4, 3, 4);
    fill(red, 100);
    gPoints.showPicked(6); // show currently picked vertex of P
    fill(orange, 100);
    show(pick(mouseX,mouseY), 5);  // show mouse position on the sphere
  }

  if (exitDraw) noLoop();

  if (snappingPDF) {
    endRecord();
    snappingPDF = false;
  }
  popMatrix(); // done with 3D drawing, restore front view for writing text

  // display header, including my face
  if (scribeText) {
    fill(black);
    displayHeader();
  }
  // display anything related to my project
  //scribeHeader("debug = " + str(debugCH), 2);
  //scribeHeader("number of faces I enter = " + numFaces, 3);
  if (scribeText && test >=12 && test <= 15) {
    scribeHeader("valid ring set? " + (validRS ? "yes" : "no"), 3);
  }
  if (scribeText && numTriangles != -1) {
    if (test < 12) {
      scribeHeader("#triangles =" + numTriangles, 4);
    } else {
      // scribeHeader("#triangles =" + numTriangles + "(expect " + str(2 * numRings - 4) + "), #rings =" + numRings, 4);
      // scribeHeader("tan0 = " + tan0 + ", tan1 = " + tan1, 5);
      // scribeHeader("arctan0 = " + atan(tan0) + ", arctan1 = " + atan(tan1), 6);
    }
  }

  //scribeHeader("time for triangle mesh generation = " + timeTM + "ms", 5);
  if (scribeText && subdivisionTimes > 0) {
    if (test == 3) {
      scribeHeader("time for subdivision = " + timeSD + "ms", 6);
    }
    if (test == 15) {
      scribeHeader("subdivision times = " + subdivisionTimes, 6);
    }
  }
  if (scribeText && debug3RT) {
    if (test == 2) {
      scribeHeader("number of steps for a three-ring triangle = " + numSteps3RT, 7);
    }
  }

  //scribeHeader("regenerate = " + str(regenerateCH), 8);
  //scribeHeader("fix penetration among 3-ring triangles = " + str(fix3RT), 9);

  if (scribeText && gRingSet.exTriPoints != null) {
    if (test == 11) {
      scribeHeader("#triangles =" + int(gRingSet.exTriPoints.size() / 3) + " #vertices =" + gRingSet.nRings, 10);
    }
  }

  // show menu at bottom, only if not filming
  if (scribeText && !filming) displayFooter();
  if (animating) {  // periodic change of time
    t += PI/360;
    if (t >= TWO_PI) t = 0;
    s = (cos(t) + 1.0) / 2;
  }
  if (filming && (animating || change)) {
    saveFrame("FRAMES/F" + nf(frameCounter++, 4) + ".tif");  // save next frame
  }
  change = false; // to avoid capturing frames when nothing happens (change is set uppn action)
}