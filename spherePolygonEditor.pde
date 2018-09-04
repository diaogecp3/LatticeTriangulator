import processing.pdf.*;


// global variables that control what to show
boolean showYellowSphere = false;
boolean generateCH = false;
boolean showPolygon = false;
boolean generateInput = false;
boolean doTests = false;
boolean showCircles = true;


float dz=500; // distance to camera. Manipulated with wheel or when 
float rx=-0.06*TWO_PI, ry=-0.04*TWO_PI;    // view angles manipulated when space pressed but not mouse
//float rx=0, ry=0;    // view angles manipulated when space pressed but not mouse
boolean twistFree=false, animating=true, tracking=false, center=true,
  gouraud=true, showControlPolygon=true, showNormals=false, showFrame=false,
  snappingPDF=false, twoD=false;
float t=0, s=0;
boolean viewpoint=false;
pt Viewer = P();
pt F = P(0,0,0);  // focus point:  the camera is looking at it (moved when 'f or 'F' are pressed
pt Of=P(100,100,0), Ob=P(110,110,0); // red point controlled by the user via mouseDrag : used for inserting vertices ...
pt Vf=P(0,0,0), Vb=P(0,0,0);
pt Pick=P();
float radius = 100;
pt centerOfSphere = P(0, 0, 0);
float rMax = 60;
int nc = 4;
int np = 6;
pt[] centers;
vec[] initDirs;

void setup() {
  myFace = loadImage("data/pic.jpg");  // load image from file pic.jpg in folder data *** replace that file with your pic of your own face
  textureMode(NORMAL);
  size(900, 900, P3D); // p3D means that we will do 3D graphics
  P.declare(); Q.declare(); // P is a polyloop in 3D: declared in pts
  // PtQ.declare(); 
  // P.resetOnCircle(12,100); // used to get started if no model exists on file 
  P.loadPts("data/pts");  // loads saved model from file
  Q.loadPts("data/pts2");  // loads saved model from file
  noSmooth();  // LEAVE HERE FOR 3D PICK TO WORK!!!
  
  P2.declare();
  P2.resetOnCircle(6,100);
  
  if (generateInput) {
    int n = 64;
    generatePointsOnSphere(P, centerOfSphere, radius, n);
  }
  
  if (doTests) {
    noLoop();
    debugCH = false;
  }
  
  //testIntersectionTwoDisks();
  centers = generateTubeCentersInSphere(centerOfSphere, radius, rMax, nc);
  initDirs = generateInitDirs(centerOfSphere, centers, nc);
}

void draw() {
  background(255);
  if(twoD) {
    fill(green); P2.drawArcs(100,4); P2.drawBalls(5);
    }
  else {
    pushMatrix();   // to ensure that we can restore the standard view before writing on the canvas
    if(snappingPDF) beginRecord(PDF, "PDFimages/P"+nf(pictureCounter++,3)+".pdf"); 

    // SET PERSPECTIVE
    float fov = PI/3.0;
    float cameraZ = (height/2.0) / tan(fov/2.0);
    camera(0,0,cameraZ,0,0,0,0,1,0  );       // sets a standard perspective
    perspective(fov, 1.0, 1.0, 10000);

    // SET VIEW
    translate(0,0,dz); // puts origin of model at screen center and moves forward/away by dz
    lights();  // turns on view-dependent lighting
    rotateX(rx); rotateY(ry); // rotates the model around the new origin (center of screen)
    rotateX(PI/2); // rotates frame around X to make X and Y basis vectors parallel to the floor
    if(center) translate(-F.x,-F.y,-F.z);
    if(viewpoint) {Viewer = viewPoint(); viewpoint=false;} // sets Viewer to the current viewpoint when ',' is pressed
    computeProjectedVectors(); // computes screen projections I, J, K of basis vectors (see bottom of pv3D): used for dragging in viewer's frame    
    if(showFrame) showFrame(150); // X-red, Y-green, Z-blue arrows

    noStroke();
    //   pp=P.idOfVertexWithClosestScreenProjectionTo(Mouse()); // id of vertex of P with closest screen projection to mouse (us in keyPressed 'x'...

    
    fill(magenta); show(P(), 3); // show center of the sphere
    Pick=pick(mouseX,mouseY);

    if(picking) { 
      P.setPickToIndexOfVertexClosestTo(Pick); // id of vertex of P with closest screen projection to mouse (us in keyPressed 'x'...
      picking=false;
      }
      
    // fill(red,100); show(Pick,5); // SHOWING THE PICK POINT ON SURFACE UNDER THE MOUSE
    // fill(green,100); show(P(),5); // SHOW PICKED VERTEX

    // Show polygon (not used in Yaohong's project)
    if (showPolygon) {
      fill(green);
      P.drawArcs(100,4);
      P.drawBalls(2); // draw curve P as cones with ball ends
      fill(red,100); P.showPicked(3); // shows currently picked vertex in red   
    } else {
      fill(green);
      //P.drawBalls(1);
      fill(red, 100);
      //P.showPicked(2);
    }

    // 3D mouse demo
    fill(red); noStroke(); show(pick(mouseX,mouseY),5);

    if (doTests) {
      testCH(numTests, numPointsPerTest, showResults);
    } else if (generateCH) {
      ArrayList<Triangle> triangles = generateConvexHull(P.G, P.nv);
      fill(red); showTriangles(triangles, P.G);
    }
    
    float r = rMax * 1.0;
    if (showCircles) {
      pt[][] points = generatePointsForCircles(centers, r, centerOfSphere, initDirs, nc, np);
      showCircles(centers, points, nc, np);
    }

    // Show the big yellow sphere
    if (showYellowSphere) {
      fill(yellow, 100); show(P(),radius);  // center of sphere = (0, 0, 0)
    }

    if(snappingPDF) {endRecord(); snappingPDF=false;}
    popMatrix(); // done with 3D drawing. Restore front view for writing text on canvas
    }
  if(scribeText) {fill(black); displayHeader();} // dispalys header on canvas, including my face
  if(scribeText && !filming) displayFooter(); // shows menu at bottom, only if not filming
  if (animating) { t+=PI/180/2; if(t>=TWO_PI) t=0; s=(cos(t)+1.)/2; } // periodic change of time 
  if(filming && (animating || change)) saveFrame("FRAMES/F"+nf(frameCounter++,4)+".tif");  // save next frame to make a movie
  change=false; // to avoid capturing frames when nothing happens (change is set uppn action)
  //scribeHeader("pp="+pp,2);
  scribeHeader("number of faces = " + numFacesShown, 3);
  }
  
void keyPressed() {
  if(key=='2') twoD=!twoD; 
  if(key=='`') picking=true; 
  if(key=='?') scribeText=!scribeText;
  if(key=='!') snapPicture();
  if(key=='@') snappingPDF=true;
  if(key=='~') filming=!filming;
  if(key==']') showControlPolygon=!showControlPolygon;
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
  if(key=='w') P.savePts("data/pts");   // save vertices to pts
  if(key=='l') P.loadPts("data/pts"); 
//  if(key=='a') animating=!animating; // toggle animation
  if(key==',') viewpoint=true;
  if(key=='>') showFrame=!showFrame;
  if(key=='#') exit();
  change=true;
  
  if (key == '0') debugCH = !debugCH;
  if (key == 'o') showYellowSphere = !showYellowSphere;
  if (key == 'h') generateCH = !generateCH;
  if (key == '+') numFacesShown++;
  if (key == '-') numFacesShown--;
  if (key == '1') numFacesShown = 1;
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
//  if (keyPressed && key=='x') P.movePicked(ToIJ(V((float)(mouseX-pmouseX),(float)(mouseY-pmouseY),0))); 
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
    fill(white); image(myFace, width-myFace.width/2,25,myFace.width/2,myFace.height/2);
    }
void displayFooter() { // Displays help text at the bottom
    scribeFooter(guide,1); 
    scribeFooter(menu,0); 
    }

String title ="FAN: Polyloop Editor on a sphere", name ="Yaohong Wu, Jarek Rossignac",
       menu="?:help, !:picture, ~:(start/stop)capture, space:rotate, s/wheel:closer, >:frame, #:quit",
       guide="x/z:select&edit, e:swap, q/p:copy, l/L: load, w/W:write to file"; // user's guide