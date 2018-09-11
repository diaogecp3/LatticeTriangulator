import processing.pdf.*;


boolean showYellowSphere = false;
boolean generateCH = false;
boolean showRingSet = true;
boolean showPointSet = true;
int inputMethodPointSet = 0;
int inputMethodRingSet = 0;
int test = 2;  // 0: convex hull 1: convex hull with holes 2: many convex-hull tests 3: many convex-hull-with-holes tests

float dz=500; // distance to camera. Manipulated with mouse wheel
float rx=-0.06*TWO_PI, ry=-0.04*TWO_PI;    // view angles manipulated when space pressed but not mouse
boolean twistFree=false, animating=true, tracking=false, center=true,
        gouraud=true, showNormals=false, showFrame=false, snappingPDF=false;
float t=0, s=0;
boolean viewpoint=false;
pt Viewer = P();
pt F = P(0,0,0);  // focus point:  the camera is looking at it (moved when 'f or 'F' are pressed
pt Of=P(100,100,0), Ob=P(110,110,0); // red point controlled by the user via mouseDrag : used for inserting vertices ...
pt Vf=P(0,0,0), Vb=P(0,0,0);
pt Pick=P();

float radiusOfSphere = 100;
pt centerOfSphere = new pt();

float rMax = 40;
float attenuation = 1.0;
int numGroups = 6;
int numPointsPerGroup = 8;

RingSet rs;

int numTriangles = -1;
float timeCH = 0.0;


void setup() {
  face0 = loadImage("data/Yaohong.jpg");  // load Yaohong's image
  face1 = loadImage("data/Jarek.jpg");  // load Jarek's image
  textureMode(NORMAL);
  size(900, 900, P3D); // P3D means that we will do 3D graphics
  noSmooth();  // LEAVE HERE FOR 3D PICK TO WORK!!!

  P.declare();
  
  if (test >= 2) {
    noLoop();
    debugCH = false;
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
      rs.loadPointGroups("data/ring_set/rs_easy_0");
      break;
    case 1:  // generate randomly
      rs = new RingSet(centerOfSphere, radiusOfSphere,
                       numGroups, numPointsPerGroup);
      rs.init();
      break;
    default:
      println("Please use a valid input method for ring set");
      exit();
  }
}

void draw() {
  background(255);
  
  {
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
    rotateX(PI / 2);  // rotates frame around X to make X and Y basis vectors parallel to the floor
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
      case 0:  // one convex-hull test
        oneConvexHullTest();
        break;
      case 1:  // one convex-hull-with-holes test
        oneConvexHullWithHolesTest();
        break;
      case 2:  // many convex-hull tests
        testCH(numTests, numPointsPerTest, showResults);
        break;
      case 3:  // many convex-hull-with-holes tests
        println("Many convex-hull-with-holes tests coming soon");
        exit();
        break;
      default:
        println("Please enter a correct test number");
        exit();
    }

    // Show the big yellow sphere
    if (showYellowSphere) {
      fill(yellow, 100);
      show(centerOfSphere, radiusOfSphere);
    }

    if (exitDraw) noLoop();

    if(snappingPDF) { endRecord(); snappingPDF = false; }
    popMatrix(); // done with 3D drawing. Restore front view for writing text on canvas
  }
  if(scribeText) {fill(black); displayHeader();} // dispalys header on canvas, including my face
  if(scribeText && !filming) displayFooter(); // shows menu at bottom, only if not filming
  if (animating) { t+=PI/180/2; if(t>=TWO_PI) t=0; s=(cos(t)+1.)/2; } // periodic change of time 
  if(filming && (animating || change)) saveFrame("FRAMES/F"+nf(frameCounter++,4)+".tif");  // save next frame to make a movie
  change=false; // to avoid capturing frames when nothing happens (change is set uppn action)
  scribeHeader("debug convex hull = " + String.valueOf(debugCH), 2);
  scribeHeader("number of steps I enter = " + numFacesShown, 3);
  if (numTriangles != -1) {
    scribeHeader("number of triangles = " + numTriangles, 4);
    scribeHeader("time to generate convex hull = " + timeCH + "ms", 5);
  }
}

void keyPressed() { 
  if(key=='`') picking=true; 
  if(key=='?') scribeText=!scribeText;
  if(key=='!') snapPicture();
  if(key=='@') snappingPDF=true;
  if(key=='~') filming=!filming;
  if(key=='|') showNormals=!showNormals;
  if(key=='G') gouraud=!gouraud;
  if(key=='q') Q.copyFrom(P);
  if(key=='p') P.projectOnSphere(100);
  if(key=='e') {PtQ.copyFrom(Q);Q.copyFrom(P);P.copyFrom(PtQ);}
  if(key=='=') {bu=fu; bv=fv;}
  // if(key=='.') F=P.Picked(); // snaps focus F to the selected vertex of P (easier to rotate and zoom while keeping it in center)
  if(key=='c') center=!center; // snaps focus F to the selected vertex of P (easier to rotate and zoom while keeping it in center)
  if(key=='t') tracking=!tracking; // snaps focus F to the selected vertex of P (easier to rotate and zoom while keeping it in center)
//  if(key=='x' || key=='z' || key=='d' || key=='a') P.setPickedIndexTo(pp); // picks the vertex of P that has closest projeciton to mouse
  if(key=='x' || key=='z' || key=='d' || key=='a') P.setPickToIndexOfVertexClosestTo(Pick); // picks the vertex of P that has closest projeciton to mouse
  if(key=='d') P.deletePicked();
  if(key=='i') P.insertClosestProjection(Pick); // Inserts new vertex in P that is the closeset projection of O
  if(key=='W') {P.savePts("data/pts"); Q.savePts("data/pts2");}  // save vertices to pts2
  if(key=='L') {P.loadPts("data/pts"); Q.loadPts("data/pts2");}   // loads saved model
  if(key=='w') {P.savePts("data/pts"); rs.savePointGroups("data/rs_unnamed");}   // save vertices to pts
  if(key=='l') {P.loadPts("data/pts"); rs.loadPointGroups("data/rs_unnamed");} 
//  if(key=='a') animating=!animating; // toggle animation
  if(key==',') viewpoint=true;
  if(key=='>') showFrame=!showFrame;
  if(key=='#') exit();
  change=true;
  
  if (key == '0') debugCH = !debugCH;
  if (key == 'o') showYellowSphere = !showYellowSphere;
  if (key == 'h') generateCH = !generateCH;
  if (key == '+') if (numTriangles >= 0) numFacesShown = numTriangles + 1;
  if (key == '-') if (numFacesShown > 0) numFacesShown--;
  if (key == '1') numFacesShown = 1;
  if (key == '[') attenuation = min(1.0, attenuation + 0.1);
  if (key == ']') attenuation = max(0.1, attenuation - 0.1);
  if (key == 'g') { showRingSet = !showRingSet; showPointSet = !showPointSet; }
}

void mouseWheel(MouseEvent event) {dz -= event.getAmount(); change=true;}
void mouseReleased() {picking=false;}
void mousePressed() {
   if (!keyPressed || key=='a') picking=true;
  }
  
void mouseMoved() {
  if (keyPressed && key==' ') {rx-=PI*(mouseY-pmouseY)/height; ry+=PI*(mouseX-pmouseX)/width;};
  if (keyPressed && key=='s') dz+=(float)(mouseY-pmouseY); // approach view (same as wheel)
  }
  
void mouseDragged() {
  if (!keyPressed) {Of.add(ToIJ(V((float)(mouseX-pmouseX),(float)(mouseY-pmouseY),0))); }
  if (keyPressed && key==CODED && keyCode==SHIFT) {Of.add(ToK(V((float)(mouseX-pmouseX),(float)(mouseY-pmouseY),0)));};
  if (keyPressed && key=='i') P.setPickedTo(Pick); 
  if (keyPressed && key=='x') P.movePickedTo(Pick);
  if (keyPressed && key=='z') P.movePicked(ToK(V((float)(mouseX-pmouseX),(float)(mouseY-pmouseY),0))); 
  if (keyPressed && key=='X') P.moveAll(ToIJ(V((float)(mouseX-pmouseX),(float)(mouseY-pmouseY),0))); 
  if (keyPressed && key=='Z') P.moveAll(ToK(V((float)(mouseX-pmouseX),(float)(mouseY-pmouseY),0))); 
  if (keyPressed && key=='f') { // move focus point on plane
    if(center) F.sub(ToIJ(V((float)(mouseX-pmouseX),(float)(mouseY-pmouseY),0))); 
    else F.add(ToIJ(V((float)(mouseX-pmouseX),(float)(mouseY-pmouseY),0))); 
    }
  if (keyPressed && key=='F') { // move focus point vertically
    if(center) F.sub(ToK(V((float)(mouseX-pmouseX),(float)(mouseY-pmouseY),0))); 
    else F.add(ToK(V((float)(mouseX-pmouseX),(float)(mouseY-pmouseY),0))); 
    }
  }  

// **** Header, footer, help text on canvas
void displayHeader() { // Displays title and authors face on screen
    scribeHeader(title,0); scribeHeaderRight(name); 
    fill(white);
    image(face0, width-face1.width/2-face0.width/2-10, 25, face0.width/2, face0.height/2+5);
    image(face1, width-face1.width/2, 25, face1.width/2, face1.height/2-5);
    }
void displayFooter() { // Displays help text at the bottom
    scribeFooter(guide,1); 
    scribeFooter(menu,0); 
    }

String title ="Lattice Triangulator", name ="Yaohong Wu, Jarek Rossignac",
       menu="?:help, !:picture, ~:(start/stop)capture, space:rotate, s/wheel:closer, >:frame, #:quit",
       guide="x/z:select&edit, e:swap, q/p:copy, l/L: load, w/W:write to file"; // user's guide