/******************************************************************************
 * Graphical user interface.
 ******************************************************************************/


void keyPressed() {
  if (key == '`') picking = true;
  if (key == '?') scribeText = !scribeText;
  if (key == '!') snapPicture();
  if (key == '@') snappingPDF = true;
  if (key == '~') filming = !filming;
  if (key == 'c') center = !center; // snaps focus F to the selected vertex of P (easier to rotate and zoom while keeping it in center)
  if (key == 't') tracking = !tracking; // snaps focus F to the selected vertex of P (easier to rotate and zoom while keeping it in center)
  if (key == 'x' || key == 'z' || key == 'd' || key == 'a') gPoints.setPickToIndexOfVertexClosestTo(gPick);  // picks the vertex of P that has closest projeciton to mouse
  if (key == 'd') {
    // P.deletePicked();
    gPoints.deletePickedPair();
    if (test == 14) debugIncCHIter = max(3, min(debugIncCHIter, gPoints.nv / 2));
    if (test == 25) gCubeHalfLength = max(gCubeHalfLength - 1, 0);
  }
  if (key == 'i') {
    if (test != 25) gPoints.addPt(gPick);  // append the new vertex gPick in P
    if (test == 25) gCubeCenter.i = max(gCubeCenter.i - 1, gCubeHalfLength);
  }
  if (key == 'j') {
    if (test == 25) gCubeCenter.j = max(gCubeCenter.j - 1, gCubeHalfLength);
  }
  if (key == 'k') {
    if (test == 25) gCubeCenter.k = max(gCubeCenter.k - 1, gCubeHalfLength);
  }
  if (key == 'w') {  // save data
    if (gPoints != null) gPoints.savePts("data/pts_unnamed");
    if (gRingSet != null) gRingSet.save("data/rs_unnamed");
    // if (gHub != null) gHub.save("data/hub_unnamed");
    if (gHub != null) gHub.saveAugFile("data/hub_aug_unnamed");
    if (gEdgeCircle != null) gEdgeCircle.save("data/ec_unnamed");
    if (gTriangleMesh != null) gTriangleMesh.save("data/tm_unnamed");
    if (gCamera != null) gCamera.save("data/cam_unnamed");
    if (gGap != null) gGap.save("data/gap_unnamed");
  }
  if (key == 'l') {
    if (gCamera != null) gCamera.load("data/point_set/cam_37");
  }
  if (key == '>') showFrame = !showFrame;

  /* Following are Yaohong's keys. */
  /* Keys: numbers. */
  if (key == '0') {
    debugCH = !debugCH;
    debug3RT = !debug3RT;
    debug2RT = !debug2RT;
    debugST = !debugST;
    if (test == 14) {
      debugIncCH = !debugIncCH;
      debugApolloniusDiagram = !debugApolloniusDiagram;
    }
    if (test == 24) debugLattice = !debugLattice;
  }
  if (key == '1') {
    gNumFaces = 1;
    numSteps3RT = 1;
    showDiskSet = !showDiskSet;
  }
  if (key == '2') {
    show2RT = !show2RT;
    fix3RT = show2RT;
    if (test == 11) showFirstCone = !showFirstCone;
  }
  if (key == '3') {
    show3RT = !show3RT;
    if (test == 11) showSecondCone = !showSecondCone;
  }
  if (key == '4') {
    showCorridorFaces = !showCorridorFaces;
  }
  if (key == '5') {
    showTriangleFaces = !showTriangleFaces;
  }
  if (key == '6') {
    showPolygons = !showPolygons;
  }
  if (key == '7') {
    simpleCorridor = !simpleCorridor;
  }
  if (key == '8') {
    showArcSet = !showArcSet;
  }
  if (key == '9') {
    if (test == 26 || test == 27) showAuxPlane = !showAuxPlane;
    if (test == 24 || test == 25) showFocus = !showFocus;
  }

  /* Keys: increase/decrease operators. */
  if (key == '+') {
    if (test == 14) {
      debugIncCHIter = min(debugIncCHIter + 1, int(gPoints.nv / 2) - 1);
    }
    if (test >= 14 && test <= 25)  {
      gNumPointsPerRing++;
    }
    if (test == 1 && gNumTriangles >= 0) {
      gNumFaces = gNumTriangles + 1;
    }
  }
  if (key == '-') {
    if (test == 14) {
      debugIncCHIter = max(debugIncCHIter - 1, 3);
    }
    if (test >= 15 && test <= 25) {
      gNumPointsPerRing = max(gNumPointsPerRing - 1, 3);
    }
    if (test == 1 && gNumFaces > 0) {
      gNumFaces--;
    }
  }
  if (key == '/') {
    if (test == 16) idxIncCor++;
    if (test == 20) projectMethod = (projectMethod + 1) % numProjectMethod;
  }
  if (key == '*') {
    if (test == 16) idxIncCor = max(0, idxIncCor - 1);
    if (test == 20) projectMethod = (projectMethod + numProjectMethod - 1) % numProjectMethod;
  }
  if (key == '[') {
    if (test == 1 || test == 3) gAttenuation = min(1.0, gAttenuation + gAttenuationDelta);
    if (test == 14) idxIncTri++;
    if (test == 4 || test == 19 || test == 20) subdivisionTimes++;
  }
  if (key == ']') {
    if (test == 1 || test == 3) gAttenuation = max(gAttenuationMin, gAttenuation - gAttenuationDelta);
    if (test == 14) idxIncTri = max(0, idxIncTri - 1);
    if (test == 4 || test == 19 || test == 20) subdivisionTimes = max(0, subdivisionTimes - 1);
  }

  /* Keys: lowercase letters. */
  if (key == 'o') {
    showSphere = !showSphere;
  }
  if (key == 'h') {
    generateCH = !generateCH;
  }
  if (key == 'r') {
    if (test == 14) debugIncCHCor = !debugIncCHCor;
    if (test == 1 || test == 3) {
      regenerateCH = !regenerateCH;
      if (regenerateCH == false) {
        gRingSet.generatePoints(gAttenuationMin);  // shrink all rings
        if (test == 1) {
          debugCH = false;
          gRingSet.generateTriangleMesh(0);  // generate a triangle mesh and store it
        }
        if (test == 3) {
          debug3RT = false;
          debug2RT = false;
          gRingSet.generateTriangleMesh(1);  // generate a triangle mesh and store it
        }
      }
    }
    if (test == 17) showCoarseCorridor = !showCoarseCorridor;
  }
  if (key == 'n') {
    if (test == 14) debugIncCHNewView = !debugIncCHNewView;
  }
  if (key == 'g') {
    showRingSet = !showRingSet;
    showPointSet = !showPointSet;
    showCircleSet = !showCircleSet;
    if (test == 24) showLattice = !showLattice;
  }
  if (key == 'b') {
    showBeams = !showBeams;
  }
  if (key == 'f') {
    fix3RT = !fix3RT;
  }
  if (key == 'm') {
    if (methodTM == 0) methodTM = 1;
    else {
      methodTM = 0;
      if (test == 3) {
        gRingSet.threeRingTriangles = null;
        gRingSet.twoRingTriangles = null;
      }
    }
  }

  /* Keys: uppercase letters. */
  if (key == 'C') {
    showCenterOfSphere = !showCenterOfSphere;
  }
  if (key == 'A') {
    showApolloniusDiagram = !showApolloniusDiagram;
  }
  if (key == 'T') {
    showTriMesh = !showTriMesh;
  }
  if (key == 'B') {
    showBoundingSphere = !showBoundingSphere;
  }
  if (key == 'O') {
    showIntersectionCircles = !showIntersectionCircles;
  }
  if (key == 'H') {
    showHub = !showHub;
  }
  if (key == 'P') {
    if (test == 2) projectOnSphere = !projectOnSphere;
    if (test == 14 || test == 28) showStereoProjection = !showStereoProjection;
    if (test == 20) projectOnHub = !projectOnHub;
  }
  if (key == 'L') {
    if (test == 20) showLiftedCones = !showLiftedCones;
  }
  if (key == 'G') {
    showGapMesh = !showGapMesh;
  }
  if (key == 'S') {
    showTriangleStrokes = !showTriangleStrokes;
    showCorridorStrokes = !showCorridorStrokes;
    if (test == 204) showSpheres = !showSpheres;
    if (test == 25) showSteadyLattice = !showSteadyLattice;
  }
  if (key == 'K') {
    showCones = !showCones;
  }
  if (key == 'I') {
    if (test == 25) gCubeCenter.i = min(gCubeCenter.i + 1, gSteadyLattice.repetitionCountU() - 1 - gCubeHalfLength);
  }
  if (key == 'J') {
    if (test == 25) gCubeCenter.j = min(gCubeCenter.j + 1, gSteadyLattice.repetitionCountV() - 1 - gCubeHalfLength);
  }
  if (key == 'K') {
    if (test == 25) gCubeCenter.k = min(gCubeCenter.k + 1, gSteadyLattice.repetitionCountW() - 1 - gCubeHalfLength);
  }
  if (key == 'D') {
    if (test == 25) gCubeHalfLength = min(gCubeHalfLength + 1, int((gSteadyLattice.minRepetitionCount() - 1) / 2));
  }
  if (key == 'e') {
    showEllipticCone1 = !showEllipticCone1;
  }
  if (key == 'E') {
    showEllipticCone2 = !showEllipticCone2;
  }

  change = true;
}

void mouseWheel(MouseEvent event) {
  gCamera.dz -= event.getAmount();
  change = true;
}

void mouseReleased() {
  picking = false;
}

void mousePressed() {
  if (!keyPressed || key == 'a') picking = true;
}

void mouseMoved() {
  if (keyPressed && key == ' ') {
    gCamera.rx -= PI * (mouseY - pmouseY) / height;
    gCamera.ry += PI * (mouseX - pmouseX) / width;
  }
  if (keyPressed && key == 's') gCamera.dz += (float)(mouseY-pmouseY); // approach view (same as wheel)
}

void mouseDragged() {
  if (!keyPressed) {
    gCamera.dx += gCamera.sxPan * (float)(mouseX - pmouseX);
    gCamera.dy += gCamera.syPan * (float)(mouseY - pmouseY);
  }

  if (keyPressed && key == 'i') gPoints.setPickedTo(gPick);
  if (keyPressed && key == 'x') gPoints.movePickedTo(gPick);
  if (keyPressed && key == 'z') gPoints.movePicked(ToK(V((float)(mouseX - pmouseX), (float)(mouseY - pmouseY), 0)));
  if (keyPressed && key == 'X') gPoints.moveAll(ToIJ(V((float)(mouseX - pmouseX), (float)(mouseY - pmouseY), 0)));
  if (keyPressed && key == 'Z') gPoints.moveAll(ToK(V((float)(mouseX - pmouseX), (float)(mouseY - pmouseY), 0)));
  if (keyPressed && key == 'f') {  // move focus point on plane
    if (center) gFocus.sub(ToIJ(V((float)(mouseX - pmouseX), (float)(mouseY - pmouseY), 0)));
    else gFocus.add(ToIJ(V((float)(mouseX - pmouseX), (float)(mouseY - pmouseY), 0)));
  }
  if (keyPressed && key == 'F') {  // move focus point vertically
    if (center) gFocus.sub(ToK(V((float)(mouseX - pmouseX), (float)(mouseY - pmouseY), 0)));
    else gFocus.add(ToK(V((float)(mouseX - pmouseX), (float)(mouseY - pmouseY), 0)));
  }
}

void displayDebugText() {
  if (test == 2) {
    scribeHeader("time for subdivision = " + timeSD + "ms", 6);
  }

  if (test == 14) {
    scribeHeader("valid ring set? " + (validRS ? "yes" : "no"), 3);
  }

  if (test == 19 || test == 20) {
    scribeHeader("subdivision times = " + subdivisionTimes, 6);
  }

  if (test == 20) {
    switch (projectMethod) {
      case 1:
        scribeHeader("Project method: shooting lines", 7);
        break;
      case 2:
        scribeHeader("Project method: sphere tracing to blended hub", 7);
        break;
      default:
        scribeHeader("Project method: shooting rays", 7);
    }
  }

  if (test == 24) {
    scribeHeader("number of balls = " + gLattice.nBalls, 2);
    scribeHeader("number of beams = " + gLattice.nBeams, 3);
    scribeHeader("corridor resolution = " + gNumPointsPerRing, 4);
    scribeHeader("time for triangulation = " + timeTM + "ms", 5);
  }
}

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

String title = "Lattice Triangulator";
String name = "Yaohong Wu, Jarek Rossignac";
String menu = "?:help, !:picture, ~:(start/stop)capture, space:rotate, s/wheel:closer";
String guide = "x:select, w:save data";
