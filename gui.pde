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
  if (key == 'x' || key == 'z' || key == 'd' || key == 'a') gPoints.setPickToIndexOfVertexClosestTo(Pick);  // picks the vertex of P that has closest projeciton to mouse
  if (key == 'd') {
    // P.deletePicked();
    gPoints.deletePickedPair();
    if (test >= 12 && test <= 15) {
      debugIncCHIter = max(3, min(debugIncCHIter, gPoints.nv / 2));
    }
  }
  if (key == 'i') {
    gPoints.addPt(Pick);  // append the new vertex Pick in P
  }
  if (key == 'w') {  // save data
    if (gPoints != null) gPoints.savePts("data/pts_unnamed");
    if (gRingSet != null) gRingSet.save("data/rs_unnamed");
    // if (gHub != null) gHub.save("data/hub_unnamed");
    if (gHub != null) gHub.saveAugFile("data/hub_aug_unnamed");
    if (gEdgeCircle != null) gEdgeCircle.save("data/ec_unnamed");
    if (gTriangleMesh != null) gTriangleMesh.save("data/tm_unnamed");
    if (gCamera != null) gCamera.save("data/cam_unnamed");
  }
  if (key == 'l') {
    // if (gCamera != null) gCamera.load("data/cam_unnamed");
    if (gCamera != null) gCamera.load("data/point_set/cam_25");
  }
  if (key == ',') viewpoint = true;
  if (key == '>') showFrame = !showFrame;

  /* Following are Yaohong's keys. */
  /* Keys: numbers. */
  if (key == '0') {
    debugCH = !debugCH;
    debug3RT = !debug3RT;
    debug2RT = !debug2RT;
    debugST = !debugST;
    if (test == 13 || test == 20) {
      debugIncCH = !debugIncCH;
      debugApolloniusDiagram = !debugApolloniusDiagram;
    }
  }
  if (key == '1') {
    numFaces = 1;
    numSteps3RT = 1;
    showDiskSet = !showDiskSet;
  }
  if (key == '2') {
    show2RT = !show2RT;
    fix3RT = show2RT;
    if (test == 19) showFirstCone = !showFirstCone;
  }
  if (key == '3') {
    show3RT = !show3RT;
    if (test == 19) showSecondCone = !showSecondCone;
  }
  if (key == '4') {
    if (test >= 12 && test <= 20) showCorridorFaces = !showCorridorFaces;
  }
  if (key == '5') {
    if (test >= 12 && test <= 20) showTriangleFaces = !showTriangleFaces;
  }
  if (key == '6') {
    showPolygons = !showPolygons;
  }

  if (key == '8') {
    if (test >= 12 && test <= 20) showArcSet = !showArcSet;
  }
  if (key == '9') {
    if ((test >= 12 && test <= 15) || test == 18) showAuxPlane = !showAuxPlane;
  }

  /* Keys: increase/decrease operators. */
  if (key == '+') {
    if (test == 13) {
      // debugIncCHIter = min(debugIncCHIter + 1, int(gPoints.nv / 2) - 1);
    }
    if (test == 13 || test == 15 || test == 16 || test == 21) {
      gNumPointsPerRing++;
    }
    if (gNumTriangles >= 0) {
      numFaces = gNumTriangles + 1;
      numFaces3RT = gNumTriangles + 1;
    }
    // gRingSet.debug2RTInfo.numGlobalStep = min(gRingSet.debug2RTInfo.numGlobalStep + 1, gRingSet.nRings);
    // gRingSet.debug2RTInfo.numLocalStep = 1;
  }
  if (key == '-') {
    if (test == 13) {
      // debugIncCHIter = max(debugIncCHIter - 1, 3);
    }
    if (test == 13 || test == 15 || test == 16 || test == 21) {
      gNumPointsPerRing = max(gNumPointsPerRing - 1, 3);
    }
    if (numFaces > 0) {
      numFaces--;
      numFaces3RT--;
    }
    // gRingSet.debug2RTInfo.numGlobalStep = max(1, gRingSet.debug2RTInfo.numGlobalStep - 1);
    // gRingSet.debug2RTInfo.numLocalStep = 1;
  }
  if (key == '/') {
    if (test == 13 || test == 14) idxIncCor++;
    if (test == 16) projectMethod = (projectMethod + 1) % numProjectMethod;
    // numSteps3RT = max(1, numSteps3RT - 1);
    // gRingSet.debug2RTInfo.numLocalStep = max(1, gRingSet.debug2RTInfo.numLocalStep - 1);
  }
  if (key == '*') {
    if (test == 13 || test == 14) idxIncCor = max(0, idxIncCor - 1);
    if (test == 16) projectMethod = (projectMethod + numProjectMethod - 1) % numProjectMethod;
    // numSteps3RT++;
    // gRingSet.debug2RTInfo.numLocalStep = min(gRingSet.debug2RTInfo.numLocalStep + 1, gRingSet.nPointsPerRing);
  }
  if (key == '[') {
    if (test == 1 || test == 2) gAttenuation = min(1.0, gAttenuation + gAttenuationDelta);
    if (test == 13) idxIncTri++;
    if (test == 3 || test == 15 || test == 16) subdivisionTimes++;
  }
  if (key == ']') {
    if (test == 1 || test == 2) gAttenuation = max(gAttenuationMin, gAttenuation - gAttenuationDelta);
    if (test == 13) idxIncTri = max(0, idxIncTri - 1);
    if (test == 3 || test == 15 || test == 16) subdivisionTimes = max(0, subdivisionTimes - 1);
  }

  /* Keys: lowercase letters. */
  if (key == 'o') {
    showSphere = !showSphere;
  }
  if (key == 'h') {
    generateCH = !generateCH;
  }
  if (key == 'r') {
    if (test == 13) debugIncCHCor = !debugIncCHCor;
    regenerateCH = !regenerateCH;
    if (regenerateCH == false) {
      gRingSet.generatePoints(gAttenuationMin);  // shrink all rings
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
    if (test == 20) {
      showCoarseCorridor = !showCoarseCorridor;
    }
  }
  if (key == 'n') {
    if (test == 13) debugIncCHNewView = !debugIncCHNewView;
  }
  if (key == 'g') {
    showRingSet = !showRingSet;
    showPointSet = !showPointSet;
    showCircleSet = !showCircleSet;
  }
  if (key == 'b') {
    showBeams = !showBeams;
  }
  if (key == 'f') {
    if (test == 10) fix3RT = !fix3RT;
  }
  if (key == 'm') {
    if (methodTM == 0) methodTM = 1;
    else {
      methodTM = 0;
      if (test == 2) {
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
    if (test == 3) projectOnSphere = !projectOnSphere;
    if (test == 16) projectOnHub = !projectOnHub;
  }
  if (key == 'L') {
    if (test == 16) showLiftedCones = !showLiftedCones;
  }
  if (key == 'G') {
    showGapMesh = !showGapMesh;
  }
  if (key == 'S') {
    showTriangleStrokes = !showTriangleStrokes;
  }
  if (key == 'K') {
    if (test == 19) showCones = !showCones;
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
    Of.add(ToIJ(V((float)(mouseX - pmouseX), (float)(mouseY - pmouseY), 0)));
  }
  if (keyPressed && key == CODED && keyCode == SHIFT) {
    Of.add(ToK(V((float)(mouseX - pmouseX), (float)(mouseY - pmouseY), 0)));
  }
  if (keyPressed && key == 'i') gPoints.setPickedTo(Pick);
  if (keyPressed && key == 'x') gPoints.movePickedTo(Pick);
  if (keyPressed && key == 'z') gPoints.movePicked(ToK(V((float)(mouseX - pmouseX), (float)(mouseY - pmouseY), 0)));
  if (keyPressed && key == 'X') gPoints.moveAll(ToIJ(V((float)(mouseX - pmouseX), (float)(mouseY - pmouseY), 0)));
  if (keyPressed && key == 'Z') gPoints.moveAll(ToK(V((float)(mouseX - pmouseX), (float)(mouseY - pmouseY), 0)));
  if (keyPressed && key == 'f') { // move focus point on plane
    if (center) F.sub(ToIJ(V((float)(mouseX - pmouseX), (float)(mouseY - pmouseY), 0)));
    else F.add(ToIJ(V((float)(mouseX - pmouseX), (float)(mouseY - pmouseY), 0)));
  }
  if (keyPressed && key == 'F') { // move focus point vertically
    if (center) F.sub(ToK(V((float)(mouseX - pmouseX), (float)(mouseY - pmouseY), 0)));
    else F.add(ToK(V((float)(mouseX - pmouseX), (float)(mouseY - pmouseY), 0)));
  }
}

void displayDebugText() {
  if (test == 13) {
    scribeHeader("valid ring set? " + (validRS ? "yes" : "no"), 3);
  }

  if (subdivisionTimes >= 0) {
    if (test == 3) {
      scribeHeader("time for subdivision = " + timeSD + "ms", 6);
    }
    if (test == 15 || test == 16) {
      scribeHeader("subdivision times = " + subdivisionTimes, 6);
    }
  }

  if (test == 16) {
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

  if (gRingSet.exTriPoints != null) {
    if (test == 11) {
      scribeHeader("#triangles =" + int(gRingSet.exTriPoints.size() / 3) + " #vertices =" + gRingSet.nRings, 10);
    }
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
