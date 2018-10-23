import processing.pdf.*;


/*
 * 0: one convex-hull test
 * 1: one convex-hull-with-holes test
 * 2: one ring-set-triangulation test (2 methods)
 * 3: one subdivision test
 * 4: one hub test
 * 5: one pivot-plane-around-line-until-hit-circle test
 * 6: one tangent-plane-of-three-circles test (3 initialization methods)
 * 7: one exact-convex-hull-for-three-circles test
 * 8: one three-ring-triangle test
 * ...
 * 10: many convex-hull tests
 * 11: many ring-set-triangulation tests
 * 12: many three-ring-triangle tests
 * 13: many extreme-plane tests
 * ...
 * 20: one circle-plane-intersection test
 */
int test = 13;

int inputMethodPointSet = 0;  // 0: read from file, 1: generate randomly
int inputMethodRingSet = 0;  // 0: read from file, 1: generate randomly
int inputMethodHub = 0;  // 0: read from file, 1: generate randomly

boolean showYellowSphere = false;
boolean generateCH = false;
boolean regenerateCH = true;  // for ring set, test shrink/grow
boolean showRingSet = true;
boolean showPointSet = true;
int subdivisionTimes = 0;

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
pt centerOfSphere = new pt();

float rMax = 50;
float attenuationMin = 0.05;
float attenuationDelta = 0.05;
float attenuation = 1.0;
int numRings = 8;
int numPointsPerRing = 8;
RingSet rs;

float r0 = 20;
float r1 = 80;
int nNeighbors = 2;
Hub hub;

int numTriangles = -1;
float timeTM = 0.0;
float timeSD = 0.0;

void setup() {
  face0 = loadImage("data/Yaohong.jpg");  // load Yaohong's image
  face1 = loadImage("data/Jarek.jpg");  // load Jarek's image
  textureMode(NORMAL);
  size(900, 900, P3D); // P3D means that we will do 3D graphics
  noSmooth();  // LEAVE HERE FOR 3D PICK TO WORK!!!

  P.declare();

  if (test >= 10 && test < 20) {
    noLoop();
    debugCH = false;
    debug3RT = false;
    debug2RT = false;
  }

  switch (inputMethodPointSet) {
    case 0:  // read from file
      P.loadPts("data/point_set/ps_easy_0");
      break;
    case 1:  // generate randomly
      generatePointsOnSphere(P, centerOfSphere, radiusOfSphere, 10);
      break;
    default:
      println("Please use a valid input method for point set");
      exit();
  }

  switch (inputMethodRingSet) {
    case 0:  // read from file
      rs = new RingSet(centerOfSphere, radiusOfSphere);
      //rs.load("data/rs_unnamed");
      rs.load("data/ring_set/rs_medium_0");
      //rs.loadRings("data/ring_set/rs_3rt_bfs_null_1");
      //rs.loadRings("data/ring_set/rs_3rt_penetration_1");
      //rs.loadRings("data/ring_set/rs_3rt_wrong_fix_2");
      //rs.loadRings("data/ring_set/rs_2rt_fail");
      //rs.loadRings("data/ring_set/rs_plane_line_circle_0");
      //rs.loadRings("data/ring_set/rs_exact_CH_0");
      break;
    case 1:  // generate randomly
      rs = new RingSet(centerOfSphere, radiusOfSphere,
                       numRings, numPointsPerRing);
      rs.init();
      break;
    default:
      println("Please use a valid input method for ring set");
      exit();
  }

  switch (inputMethodHub) {
    case 0:  // read from file
      hub = new Hub();
      hub.load("data/hub/hub_easy_1");
      break;
    case 1:  // generate randomly
      hub = generateHub(centerOfSphere, r0, r1, nNeighbors);
      break;
    default:
      println("Please use a valid input method for hub");
      exit();
  }

  // check validity of ring set
  if (!rs.isValid()) {
    println("ring set not valid!");
    exit();
  }

  if (regenerateCH == false) {
    attenuation = attenuationMin;
    rs.generatePoints(attenuation);  // shrink all rings
    if (test == 1) {
      debugCH = false;
      rs.generateTriangleMesh(0);  // generate a triangle mesh and store it
    }
    if (test == 2) {
      debug3RT = false;
      debug2RT = false;
      rs.generateTriangleMesh(1);  // generate a triangle mesh and store it
    }
  }
}

void draw() {
  background(255);

  pushMatrix();  // to ensure that we can restore the standard view before writing on the canvas
  if (snappingPDF) beginRecord(PDF, "PDFimages/P" + nf(pictureCounter++,3) + ".pdf");

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

  noStroke();
  fill(magenta); show(centerOfSphere, 4); // show center of the sphere
  Pick = pick(mouseX,mouseY);

  if (picking) {
    P.setPickToIndexOfVertexClosestTo(Pick); // id of vertex of P with closest screen projection to mouse (us in keyPressed 'x'...
    picking = false;
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
      tangentPlaneThreeCirclesTest();
      break;
    case 7:
      exactCHThreeCirclesTest();
      break;
    case 8:
      threeRingTriangleTest();
      break;
    case 10:  // many convex-hull tests
      convexHullTests(numTests, numPointsPerTest);
      break;
    case 11:  // many ring-set-triangulation tests
      ringSetTriangulationTests(numTests, numRings, numPointsPerRing, attenuation);
      break;
    case 12:  // many three-ring-triangle tests
      threeRingTriangleTests(numTests, numPointsPerRing, attenuation);
      break;
    case 13:  // many extreme-plane(-of-three-circles) tests
      extremePlaneTests(numTests, attenuation);
      break;

    case 20:
      testIntersectionCirclePlane();
      break;
    default:
      println("Please enter a correct test number");
      exit();
  }

  // Show the big yellow sphere
  if (showYellowSphere) {
    fill(yellow, 100);
    noStroke();
    show(centerOfSphere, radiusOfSphere);
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
  scribeHeader("debug = " + str(debugCH), 2);
  scribeHeader("number of faces I enter = " + numFaces, 3);
  if (numTriangles != -1) {
    scribeHeader("number of triangles = " + numTriangles, 4);
    scribeHeader("time for triangle mesh generation = " + timeTM + "ms", 5);
  }
  if (test == 3 && subdivisionTimes > 0) {
    scribeHeader("time for subdivision = " + timeSD + "ms", 6);
  }
  if (test == 2 && debug3RT) {
    scribeHeader("number of steps for a three-ring triangle = " + numSteps3RT, 7);
  }
  scribeHeader("regenerate = " + str(regenerateCH), 8);
  scribeHeader("fix penetration among 3-ring triangles = " + str(fix3RT), 9);
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

void keyPressed() {
  if (key == '`') picking = true;
  if (key == '?') scribeText = !scribeText;
  if (key == '!') snapPicture();
  if (key == '@') snappingPDF = true;
  if (key == '~') filming = !filming;
  if (key == 'q') Q.copyFrom(P);
  if (key == 'p') P.projectOnSphere(100);
  if (key == 'e') { PtQ.copyFrom(Q); Q.copyFrom(P); P.copyFrom(PtQ); }
  if (key == '=') { bu = fu; bv = fv; }
  // if (key == '.') F=P.Picked(); // snaps focus F to the selected vertex of P (easier to rotate and zoom while keeping it in center)
  if (key == 'c') center = !center; // snaps focus F to the selected vertex of P (easier to rotate and zoom while keeping it in center)
  if (key == 't') tracking = !tracking; // snaps focus F to the selected vertex of P (easier to rotate and zoom while keeping it in center)
  // if (key == 'x' || key == 'z' || key == 'd' || key == 'a') P.setPickedIndexTo(pp); // picks the vertex of P that has closest projeciton to mouse
  if (key == 'x' || key == 'z' || key == 'd' || key == 'a') P.setPickToIndexOfVertexClosestTo(Pick);  // picks the vertex of P that has closest projeciton to mouse
  if (key == 'd') P.deletePicked();
  if (key == 'i') P.insertClosestProjection(Pick);  // Inserts new vertex in P that is the closeset projection of O
  if (key == 'W') { P.savePts("data/pts"); Q.savePts("data/pts2"); }  // save vertices to pts2
  if (key == 'L') { P.loadPts("data/pts"); Q.loadPts("data/pts2"); }  // loads saved model
  if (key == 'w') {  // save data
    P.savePts("data/pts_unnamed");
    rs.save("data/rs_unnamed");
    hub.save("data/hub_unnamed");
  }
  if (key == 'l') {  // load data
    P.loadPts("data/pts_unnamed");
    rs.load("data/rs_unnamed");
    hub.load("data/hub.unnamed");
  }
  // if (key == 'a') animating = !animating; // toggle animation
  if (key == ',') viewpoint = true;
  if (key == '>') showFrame = !showFrame;
  if (key == '#') exit();

  /* Following are Yaohong's keys. */
  if (key == '0') {
    debugCH = !debugCH;
    debug3RT = !debug3RT;
    debug2RT = !debug2RT;
  }
  if (key == 'o') showYellowSphere = !showYellowSphere;
  if (key == 'h') generateCH = !generateCH;
  if (key == 'r') {
    regenerateCH = !regenerateCH;
    if (regenerateCH == false) {
      rs.generatePoints(attenuationMin);  // shrink all rings
      if (test == 1) {
        debugCH = false;
        rs.generateTriangleMesh(0);  // generate a triangle mesh and store it
      }
      if (test == 2) {
        debug3RT = false;
        debug2RT = false;
        rs.generateTriangleMesh(1);  // generate a triangle mesh and store it
      }
    }
  }
  if (key == '+') {
    if (numTriangles >= 0) {
      numFaces = numTriangles + 1;
      numFaces3RT = numTriangles + 1;
    }
    subdivisionTimes++;
    rs.debug2RTInfo.numGlobalStep = min(rs.debug2RTInfo.numGlobalStep + 1, rs.nRings);
    rs.debug2RTInfo.numLocalStep = 1;
  }
  if (key == '-') {
    if (numFaces > 0) {
      numFaces--;
      numFaces3RT--;
    }
    subdivisionTimes = max(0, subdivisionTimes - 1);
    rs.debug2RTInfo.numGlobalStep = max(1, rs.debug2RTInfo.numGlobalStep - 1);
    rs.debug2RTInfo.numLocalStep = 1;
  }
  if (key == '*') {
    numSteps3RT++;
    rs.debug2RTInfo.numLocalStep = min(rs.debug2RTInfo.numLocalStep + 1, rs.nPointsPerRing);
  }
  if (key == '/') {
    numSteps3RT = max(1, numSteps3RT - 1);
    rs.debug2RTInfo.numLocalStep = max(1, rs.debug2RTInfo.numLocalStep - 1);
  }
  if (key == '1') {
    numFaces = 1;
    numSteps3RT = 1;
  }
  if (key == '2') {
    show2RT = !show2RT;
    fix3RT = show2RT;
  }
  if (key == '3') {
    show3RT = !show3RT;
  }

  if (key == '[') attenuation = min(1.0, attenuation + attenuationDelta);
  if (key == ']') attenuation = max(attenuationMin, attenuation - attenuationDelta);
  if (key == 'g') {
    showRingSet = !showRingSet;
    showPointSet = !showPointSet;
  }
  if (key == 'f' && test == 8) {
    fix3RT = !fix3RT;
  }
  if (key == 'm') {
    if (methodTM == 0) methodTM = 1;
    else {
      methodTM = 0;
      if (test == 2) {
        rs.threeRingTriangles = null;
        rs.twoRingTriangles = null;
      }
    }
  }

  change = true;
}

void mouseWheel(MouseEvent event) {
  dz -= event.getAmount();
  change = true;
}
void mouseReleased() { picking = false; }
void mousePressed() {
  if (!keyPressed || key == 'a') picking = true;
}
void mouseMoved() {
  if (keyPressed && key == ' ') {
    rx -= PI * (mouseY - pmouseY) / height;
    ry += PI * (mouseX - pmouseX) / width;
  }
  if (keyPressed && key == 's') dz += (float)(mouseY-pmouseY); // approach view (same as wheel)
}
void mouseDragged() {
  if (!keyPressed) {
    Of.add(ToIJ(V((float)(mouseX - pmouseX), (float)(mouseY - pmouseY), 0)));
  }
  if (keyPressed && key == CODED && keyCode == SHIFT) {
    Of.add(ToK(V((float)(mouseX - pmouseX), (float)(mouseY - pmouseY), 0)));
  }
  if (keyPressed && key == 'i') P.setPickedTo(Pick);
  if (keyPressed && key == 'x') P.movePickedTo(Pick);
  if (keyPressed && key == 'z') P.movePicked(ToK(V((float)(mouseX - pmouseX), (float)(mouseY - pmouseY), 0)));
  if (keyPressed && key == 'X') P.moveAll(ToIJ(V((float)(mouseX - pmouseX), (float)(mouseY - pmouseY), 0)));
  if (keyPressed && key == 'Z') P.moveAll(ToK(V((float)(mouseX - pmouseX), (float)(mouseY - pmouseY), 0)));
  if (keyPressed && key == 'f') { // move focus point on plane
    if (center) F.sub(ToIJ(V((float)(mouseX - pmouseX), (float)(mouseY - pmouseY), 0)));
    else F.add(ToIJ(V((float)(mouseX - pmouseX), (float)(mouseY - pmouseY), 0)));
  }
  if (keyPressed && key == 'F') { // move focus point vertically
    if (center) F.sub(ToK(V((float)(mouseX - pmouseX), (float)(mouseY - pmouseY), 0)));
    else F.add(ToK(V((float)(mouseX - pmouseX), (float)(mouseY - pmouseY), 0)));
  }
}

// **** Header, footer, help text on canvas
void displayHeader() { // Displays title and authors face on screen
  scribeHeader(title, 0);
  scribeHeaderRight(name);
  fill(white);
  image(face0, width - (face1.width + face0.width) / 2 - 10, 25,
        face0.width / 2, face0.height / 2 + 5);
  image(face1, width - face1.width / 2, 25,
        face1.width / 2, face1.height / 2 - 5);
}

void displayFooter() { // Displays help text at the bottom
  scribeFooter(guide, 1);
  scribeFooter(menu, 0);
}

String title = "Lattice Triangulator",
       name = "Yaohong Wu, Jarek Rossignac",
       menu = "?:help, !:picture, ~:(start/stop)capture, " +
              "space:rotate, s/wheel:closer, >:frame, #:quit",
       guide = "x/z:select&edit, e:swap, q/p:copy, l/L: load, w/W:write to file";