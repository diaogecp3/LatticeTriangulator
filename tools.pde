// ************************************ IMAGES & VIDEO
int pictureCounter=0, frameCounter=0;
Boolean filming=false, change=false;
PImage face0; // picture of first author's face, should be under data/ in sketch folder
PImage face1; // picture of second author's face, should be under data/ in sketch folder
void snapPicture() {saveFrame("PICTURES/P"+nf(pictureCounter++,3)+".jpg"); }

// ******************************************COLORS
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

void pen(color c, float w) {stroke(c); strokeWeight(w);}

// ******************************** TEXT , TITLE, and USER's GUIDE
Boolean scribeText=true; // toggle for displaying of help text
void scribe(String S, float x, float y) {fill(0); text(S,x,y); noFill();} // writes on screen at (x,y) with current fill color
void scribeHeader(String S, int i) {fill(0); text(S,10,20+i*20); noFill();} // writes black at line i
void scribeHeaderRight(String S) {fill(0); text(S,width-7.5*S.length(),20); noFill();} // writes black on screen top, right-aligned
void scribeFooter(String S, int i) {fill(0); text(S,10,height-10-i*20); noFill();} // writes black on screen at line i from bottom
void scribeAtMouse(String S) {fill(0); text(S,mouseX,mouseY); noFill();} // writes on screen near mouse
void scribeMouseCoordinates() {fill(black); text("("+mouseX+","+mouseY+")",mouseX+7,mouseY+25); noFill();}

// **************************** FILE SELECTION FOR SAVING AND LOADING MODELS
//String fileName="data/points";

//String path="data/pts";
//void saveToFile(File selection) {
//  if (selection == null) println("Window was closed or the user hit cancel.");
//  else path=selection.getAbsolutePath();
//  println("    save path = "+path);
//  }
//
//void readFromFile(File selection) {
//  if (selection == null) println("Window was closed or the user hit cancel or file not found.");
//  else path=selection.getAbsolutePath();
//  println("    read path = "+path);
//  }
//
//
//void fileSelected(File selection) {
//  if (selection == null) println("Window was closed or the user hit cancel.");
//  else {
//    fileName = selection.getAbsolutePath();
//    println("User selected " + fileName);
//    }
//  }
//