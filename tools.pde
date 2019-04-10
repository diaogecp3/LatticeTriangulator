/******************************************************************************
 * Tools.
 *
 * Adapted from the version provided by Prof. Jarek Rossignac at Georgia Tech.
 ******************************************************************************/


// ************************************ IMAGES & VIDEO
PImage face0; // picture of first author's face, should be under data/ in sketch folder
PImage face1; // picture of second author's face, should be under data/ in sketch folder
int pictureCounter = 0, frameCounter = 0;
boolean filming = false, change = false;

void snapPicture() {
  println("saving picture", pictureCounter);
  saveFrame("pictures/P" + nf(pictureCounter++,3) + ".png");
}

class Camera {
  float dz = 500;  // distance to camera, manipulated with mouse wheel
  float rx = -0.06 * TWO_PI;  // view angle, manipulated when space pressed
  float ry = -0.04 * TWO_PI;  // view angle, manipulated when space pressed

  float dx = 0;
  float dy = 0;

  final float sxPan = 0.2;  // scaling factor for panning in x direction
  final float syPan = 0.2;  // scaling factor for panning in y direction

  final float fov = PI / 3.0;  // field of view
  final float zNear = 1.0;
  final float zFar = 10000.0;

  Camera(float dz, float rx, float ry) {
    this.dz = dz;
    this.rx = rx;
    this.ry = ry;
  }

  void save(String file) {
    println("saving camera:", file);
    String[] lines = new String[3];
    lines[0] = str(dz);
    lines[1] = str(rx);
    lines[2] = str(ry);
    saveStrings(file, lines);
  }

  void load(String file) {
    println("loading camera:", file);
    String[] lines = loadStrings(file);
    dz = float(lines[0]);
    rx = float(lines[1]);
    ry = float(lines[2]);
  }
}

// ************************************ COLORS
// For more color and color names, see https://htmlcolorcodes.com/color-names/
color red = #FF0000, darkRed = #8B0000, firebrick = #B22222;
color pink = #FFC0CB, deepPink = #FF1493, hotPink = #FF69B4;
color orange = #FFA500, tomato = #FF6347, lightSalmon = #FFA07A;
color yellow = #FFFF00, gold = #FFD700, khaki = #F0E68C;
color violet = #EE82EE, magenta = #FF00FF, purple = #800080;
color green = #008000, lime = #00FF00, springGreen = #00FF7F, lightGreen = #90EE90, darkGreen = #006400;
color blue = #0000FF, cyan = #00FFFF, navy = #000080;
color brown = #A52A2A, chocolate = #D2691E, sandyBrown = #F4A460;
color white = #FFFFFF, snow = #FFFAFA, ivory = #FFFFF0;
color black = #000000, gray = #808080, silver = #C0C0C0;

void pen(color c, float w) {
  stroke(c);
  strokeWeight(w);
}

// ************************************ TEXT, TITLE, and USER's GUIDE
boolean scribeText = true;  // toggle for displaying of help text

void scribe(String S, float x, float y) {  // write on screen at (x,y) with current fill color
  fill(0);
  text(S, x, y);
  noFill();
}

void scribeHeader(String S, int i) {  // write black at line i
  fill(0);
  text(S, 10, 20 + i * 20);
  noFill();
}

void scribeHeaderRight(String S) {  // write black on screen top, right-aligned
  fill(0);
  text(S, width - 7.5 * S.length(), 20);
  noFill();
}

void scribeFooter(String S, int i) {  // write black on screen at line i from bottom
  fill(0);
  text(S, 10, height - 10 - i*20);
  noFill();
}

void scribeAtMouse(String S) {  // write on screen near mouse
  fill(0);
  text(S, mouseX, mouseY);
  noFill();
}

void scribeMouseCoordinates() {  // write mouse coordinates near mouse
  fill(black);
  text("(" + mouseX + "," + mouseY + ")", mouseX + 7, mouseY + 25);
  noFill();
}