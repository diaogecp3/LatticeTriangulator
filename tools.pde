// ************************************ IMAGES & VIDEO
PImage face0; // picture of first author's face, should be under data/ in sketch folder
PImage face1; // picture of second author's face, should be under data/ in sketch folder
int pictureCounter = 0, frameCounter = 0;
boolean filming = false, change = false;

void snapPicture() {
  saveFrame("pictures/P"+nf(pictureCounter++,3)+".png");
}

// ************************************ COLORS
// For more color and color names, see https://htmlcolorcodes.com/color-names/
color red = #FF0000, darkRed = #8B0000, firebrick = #B22222;
color pink = #FFC0CB, deepPink = #FF1493, hotPink = #FF69B4;
color orange = #FFA500, tomato = #FF6347, lightSalmon = #FFA07A;
color yellow = #FFFF00, gold = #FFD700, khaki = #F0E68C;
color violet = #EE82EE, magenta = #FF00FF, purple = #800080;
color green = #008000, lime = #00FF00, springGreen = #00FF7F;
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