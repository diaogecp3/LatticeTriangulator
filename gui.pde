/******************************************************************************
 * Graphical user interface.
 ******************************************************************************/

String warningMsg = "";

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
    if (test == 25) gCubeHalfLength = max(gCubeHalfLength - 1, 1);
  }
  if (key == 'i') {
    if (test != 25) gPoints.addPt(gPick);  // append the new vertex gPick in P
    if (test == 25) gCubeCenter.i = max(gCubeCenter.i - 1, 0);
  }
  if (key == 'j') {
    if (test == 25) gCubeCenter.j = max(gCubeCenter.j - 1, 0);
  }
  if (key == 'k') {
    if (test == 25) gCubeCenter.k = max(gCubeCenter.k - 1, 0);
  }
  if (key == 'w') {  // save data
    if (gPoints != null) gPoints.save("data/pts_unnamed");
    if (gRingSet != null) gRingSet.save("data/rs_unnamed");
    if (gHub != null) gHub.save("data/hub_unnamed");
    // if (gHub != null) gHub.saveAugFile("data/hub_aug_unnamed");
    if (gLattice != null) gLattice.save("data/lattice_unnamed");
    if (gEdgeCircle != null) gEdgeCircle.save("data/ec_unnamed");
    if (gTriangleMesh != null) gTriangleMesh.save("data/tm_unnamed");
    if (gCamera != null) gCamera.save("data/cam_unnamed");
    if (gGap != null) gGap.save("data/gap_unnamed");
  }
  if (key == 'l') {
    if (gCamera != null) {
      gCamera.load("data/camera/cam_hub_tessellation");
    }
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
    if (test == 24 || test == 25) debugLattice = !debugLattice;
  }
  if (key == '1') {
    gNumFaces = 1;
    gNumSteps3RT = 1;
    gShowDiskSet = !gShowDiskSet;
  }
  if (key == '2') {
    gShow2RT = !gShow2RT;
    gFix3RT = gShow2RT;
    if (test == 11) showFirstCone = !showFirstCone;
  }
  if (key == '3') {
    gShow3RT = !gShow3RT;
    if (test == 11) showSecondCone = !showSecondCone;
  }
  if (key == '4') {
    gShowCorridorFaces = !gShowCorridorFaces;
  }
  if (key == '5') {
    gShowTriangleFaces = !gShowTriangleFaces;
  }
  if (key == '6') {
    gShowPolygons = !gShowPolygons;
  }
  if (key == '7') {
    gUseSimpleCorridor = !gUseSimpleCorridor;
  }
  if (key == '8') {
    gShowArcSet = !gShowArcSet;
  }
  if (key == '9') {
    if (test == 26 || test == 27) gShowAuxPlane = !gShowAuxPlane;
    if (test == 23 || test == 24 || test == 25) showFocus = !showFocus;
    if (test == 20) gCreateGap = !gCreateGap;
  }


  /* Keys: increase/decrease operators. */
  if (key == '+') {
    if (test == 14 && debugIncCH) {
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
    if (test == 14 && debugIncCH) {
      debugIncCHIter = max(debugIncCHIter - 1, 3);
    }
    if (test >= 14 && test <= 25) {
      gNumPointsPerRing = max(gNumPointsPerRing - 1, 2);
    }
    if (test == 1 && gNumFaces > 0) {
      gNumFaces--;
    }
  }
  if (key == '/') {
    if (test == 16) gIdxIncCorridor++;
    if (test == 20 || test == 23 || test == 24 || test == 25) gMethodProjection = (gMethodProjection + 1) % gNumProjectMethods;
  }
  if (key == '*') {
    if (test == 16) gIdxIncCorridor = max(0, gIdxIncCorridor - 1);
    if (test == 20 || test == 23 || test == 24 || test == 25) gMethodProjection = (gMethodProjection + gNumProjectMethods - 1) % gNumProjectMethods;
  }
  if (key == '[') {
    if (test == 1 || test == 3) gAttenuation = min(1.0, gAttenuation + gAttenuationDelta);
    if (test == 14) gIdxIncTriangle++;
    if (test == 4 || test == 19 || test == 20 || test == 23 || test == 24 || test == 25) gSubdivisonTimes++;
  }
  if (key == ']') {
    if (test == 1 || test == 3) gAttenuation = max(gAttenuationMin, gAttenuation - gAttenuationDelta);
    if (test == 14) gIdxIncTriangle = max(0, gIdxIncTriangle - 1);
    if (test == 4 || test == 19 || test == 20 || test == 23 || test == 24 || test == 25) gSubdivisonTimes = max(0, gSubdivisonTimes - 1);
  }

  /* Keys: lowercase letters. */
  if (key == 'p') {
    gProjectOnCircleAfterSub = !gProjectOnCircleAfterSub;
  }
  if (key == 'o') {
    gShowSphere = !gShowSphere;
  }
  if (key == 'h') {
    gGenerateCH = !gGenerateCH;
  }
  if (key == 'r') {
    if (test == 14) debugIncCHCor = !debugIncCHCor;
    if (test == 1 || test == 3) {
      gRegenerateCH = !gRegenerateCH;
      if (gRegenerateCH == false) {
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
    if (test == 20) gInsertInfCircle = !gInsertInfCircle;
  }
  if (key == 'g') {
    gShowRingSet = !gShowRingSet;
    gShowPointSet = !gShowPointSet;
    gShowCircleSet = !gShowCircleSet;
    if (test == 24 || test == 25) gShowLattice = !gShowLattice;
  }
  if (key == 'b') {
    if (test == 19 || test == 24 || test == 25) gShowBeams = !gShowBeams;
  }
  if (key == 'f') {
    gFix3RT = !gFix3RT;
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

    if (test == 18 || test == 20 || test == 24 || test == 25) {
      gMethodConvexGap = 1 - gMethodConvexGap;
    }
  }

  /* Keys: uppercase letters. */
  if (key == 'C') {
    gShowCenterOfSphere = !gShowCenterOfSphere;
  }
  if (key == 'A') {
    gShowApolloniusDiagram = !gShowApolloniusDiagram;
  }
  if (key == 'T') {
    gShowTriMesh = !gShowTriMesh;
  }
  if (key == 'Q') {
    if (test == 20 || test == 24 || test == 25) {
      gUseTriQuadMesh = !gUseTriQuadMesh;
    }
  }
  if (key == 'B') {
    gShowBoundingSphere = !gShowBoundingSphere;
  }
  if (key == 'O') {
    gShowIntersectionCircles = !gShowIntersectionCircles;
  }
  if (key == 'H') {
    gShowHub = !gShowHub;
    if (test == 24 || test == 25) gShowCHoCCs = !gShowCHoCCs;
  }
  if (key == 'P') {
    if (test == 2) gProjectOnSphere = !gProjectOnSphere;
    if (test == 14 || test == 28) showStereoProjection = !showStereoProjection;
  }
  if (key == 'L') {
    if (test == 20) gShowLiftedCones = !gShowLiftedCones;
  }
  if (key == 'G') {
    gShowGapMesh = !gShowGapMesh;
  }
  if (key == 'S') {
    gShowTriangleStrokes = !gShowTriangleStrokes;
    gShowCorridorStrokes = !gShowCorridorStrokes;
    if (test == 204) gShowTwoSpheres = !gShowTwoSpheres;
    if (test == 25) gShowSteadyLattice = !gShowSteadyLattice;
  }
  if (key == 'K') {
    gShowCones = !gShowCones;
  }
  if (key == 'I') {
    if (test == 25) gCubeCenter.i = min(gCubeCenter.i + 1, gSteadyLattice.repetitionCountU() - 1);
  }
  if (key == 'J') {
    if (test == 25) gCubeCenter.j = min(gCubeCenter.j + 1, gSteadyLattice.repetitionCountV() - 1);
  }
  if (key == 'K') {
    if (test == 25) gCubeCenter.k = min(gCubeCenter.k + 1, gSteadyLattice.repetitionCountW() - 1);
  }
  if (key == 'D') {
    if (test == 25) gCubeHalfLength = min(gCubeHalfLength + 1, gSteadyLattice.maxRepetitionCount() - 1);
  }
  if (key == 'e') {
    if (test == 27) gShowEllipticCone1 = !gShowEllipticCone1;
  }
  if (key == 'E') {
    if (test == 27) gShowEllipticCone2 = !gShowEllipticCone2;
    if (test == 20) gShowExplodedView = !gShowExplodedView;
  }
  if (key == 'N') {
    gNavigateSteadyLattice = !gNavigateSteadyLattice;
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

private void displayProjectMethod(int line) {
  switch (gMethodProjection) {
    case 1:
      scribeHeader("project to exact surface", line);
      break;
    case 2:
      scribeHeader("project to blended surface", line);
      break;
    default:
      scribeHeader("no projection", line);
  }
}

private void displayGapMethod(int line) {
  if (gMethodConvexGap == 0) {
    scribeHeader("method 1: wants two negative dot products", line);
  } else if (gMethodConvexGap == 1) {
    scribeHeader("method 2: wants smaller pivot angle", line);
  }
}

void displayDebugText() {
  /* Display warning message at the top. */
  if (warningMsg.length() > 0) scribeHeader("Warning: " + warningMsg, 1);

  int line = 2;
  if (test == 2) {
    scribeHeader("time for subdivision = " + gTimeSubdivision + "ms", line++);
  }

  if (test == 14) {
    scribeHeader("valid ring set? " + (gRingSet.valid ? "yes" : "no"), line++);
    scribeHeader("number of samples per circle = " + gNumPointsPerRing, line++);
  }

  if (test == 18) {
    scribeHeader("nv0 = " + gGap.points0.size(), line++);
    scribeHeader("nv1 = " + gGap.points1.size(), line++);
    scribeHeader("number of triangles = " + gTriangleMesh.nt, line++);
    if (gMethodConvexGap == 0) {
      scribeHeader("method 1: wants two negative dot products", line++);
    } else if (gMethodConvexGap == 1) {
      scribeHeader("method 2: wants smaller pivot angle", line++);
    }
  }

  if (test == 19 || test == 20) {
    scribeHeader("number of sides of each beam = " + gNumPointsPerRing, line++);
    scribeHeader("subdivision times = " + gSubdivisonTimes, line++);
    scribeHeader("insert an inf circle: " + (gInsertInfCircle ? "yes" : "no"), line++);
  }

  if (test == 20 || test == 23) {
    displayGapMethod(line++);
    scribeHeader("push vertices on circles: " + (gProjectOnCircleAfterSub ? "yes" : "no"), line++);
    displayProjectMethod(line++);
  }

  if (test == 24) {
    scribeHeader("debug? " + (debugLattice ? "yes" : "no"), line++);
    scribeHeader("valid lattice? " + (gLattice.valid ? "yes" : "no"), line++);
    scribeHeader("number of balls = " + gLattice.nBalls, line++);
    scribeHeader("number of beams = " + gLattice.nBeams, line++);
    scribeHeader("corridor resolution = " + gNumPointsPerRing, line++);
    scribeHeader("time for triangulation = " + gTimeMeshing + "ms", line++);
    scribeHeader("show beams = " + (gShowBeams ? "yes" : "no"), line++);
    displayGapMethod(line++);
    displayProjectMethod(line++);
  }

  if (test == 25) {
    scribeHeader("debug? " + (debugLattice ? "yes" : "no"), line++);
    scribeHeader("u, v, w = " + gSteadyLattice.repetitionCounts(), line++);
    scribeHeader("cube center = " + gCubeCenter, line++);
    scribeHeader("cube half length = " + gCubeHalfLength, line++);
    scribeHeader("number of balls = " + gLattice.nBalls, line++);
    scribeHeader("number of beams = " + gLattice.nBeams, line++);
    scribeHeader("corridor resolution = " + gNumPointsPerRing, line++);
    scribeHeader("time for triangulation = " + gTimeMeshing + "ms", line++);
    scribeHeader("subdivision times = " + gSubdivisonTimes, line++);
    scribeHeader("project on circles = " + (gProjectOnCircleAfterSub ? "yes" : "no"), line++);
    scribeHeader("show beams = " + (gShowBeams ? "yes" : "no"), line++);
    displayGapMethod(line++);
    displayProjectMethod(line++);
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
