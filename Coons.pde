/******** Editor of an Animated Coons Patch

Implementation steps:
**<01 Manual control of (u,v) parameters. 
**<02 Draw 4 boundary curves CT(u), CB(u), SL(v), CR(v) using proportional Neville
**<03 Compute and show Coons point C(u,v)
**<04 Display quads filed one-by-one for the animated Coons patch
**<05 Compute and show normal at C(u,v) and a ball ON the patch

*/
//**<01: mouseMoved; 'v', draw: uvShow()
float u=0, v=0; 
float bu=0.5, bv=0.5;  // back foot
float fu=0.55, fv=0.55; // front foot

void uvShow() { 
  fill(red);
  if(keyPressed && key=='v')  text("u="+u+", v="+v,10,30);
  noStroke(); fill(blue); ellipse(u*width,v*height,5,5); 
  }
/*
0 1 2 3 
11    4
10    5
9 8 7 6
*/
pt coons(pt[] P, float s, float t) {
  pt Lst = L( L(P[0],s,P[3]), t, L(P[9],s,P[6]) ) ;
  pt Lt = L( N( 0,P[0], 1./3,P[1],  2./3,P[2],  1,P[3], s) ,t, N(0,P[9], 1./3,P[8], 2./3,P[7], 1,P[6], s) ) ;
  pt Ls = L( N( 0,P[0], 1./3,P[11], 2./3,P[10], 1,P[9], t) ,s, N(0,P[3], 1./3,P[4], 2./3,P[5] ,1,P[6], t) ) ;
  return P(Ls,V(Lst,Lt));
  }
pt B(pt A, pt B, pt C, float s) {return L(L(A,s,B),s,L(B,s,C)); } 
pt B(pt A, pt B, pt C, pt D, float s) {return L(B(A,B,C,s),s,B(B,C,D,s)); } 
pt B(pt A, pt B, pt C, pt D, pt E, float s) {return L(B(A,B,C,D,s),s,B(B,C,D,E,s)); } 
pt N(float a, pt A, float b, pt B, float t) {return L(A,(t-a)/(b-a),B);}
pt N(float a, pt A, float b, pt B, float c, pt C, float t) {return N(a,N(a,A,b,B,t),c,N(b,B,c,C,t),t);}
pt N(float a, pt A, float b, pt B, float c, pt C, float d, pt D, float t) {return N(a,N(a,A,b,B,c,C,t),d,N(b,B,c,C,d,D,t),t);}
pt N(float a, pt A, float b, pt B, float c, pt C, float d, pt D, float e, pt E, float t) {return N(a,N(a,A,b,B,c,C,d,D,t),e,N(b,B,c,C,d,D,e,E,t),t);}

void drawBorders(pt[] P){
  float e=0.01;
  beginShape(); for(float t=0; t<1.001; t+=e) v(coons(P,0,t)); endShape();
  beginShape(); for(float t=0; t<1.001; t+=e) v(coons(P,1,t)); endShape();
  beginShape(); for(float t=0; t<1.001; t+=e) v(coons(P,t,0)); endShape();
  beginShape(); for(float t=0; t<1.001; t+=e) v(coons(P,t,1)); endShape();
  }
  
vec CoonsNormal(pt[] P, float u, float v, float e) { 
  vec Tu = V(coons(P,u-e,v),coons(P,u+e,v));
  vec Tv = V(coons(P,u,v-e),coons(P,u,v+e));
  return U(N(Tu,Tv));
  }

void shadeSurface(pt[] P, float e){ 
  for(float s=0; s<1.001-e; s+=e) for(float t=0; t<1.001-e; t+=e) 
  { beginShape(); v(coons(P,s,t)); v(coons(P,s+e,t)); v(coons(P,s+e,t+e)); v(coons(P,s,t+e)); endShape(CLOSE);}
  }

void shadeSurfaceTextured(pt[] P, float e){ 
  fill(white);
  for(float s=0; s<1.001-e; s+=e) for(float t=0; t<1.001-e; t+=e) 
  {beginShape(); texture(myFace); vTextured(coons(P,s,t),s,t); vTextured(coons(P,s+e,t),s+e,t); vTextured(coons(P,s+e,t+e),s+e,t+e); vTextured(coons(P,s,t+e),s,t+e); endShape(CLOSE);}
  }

void shadeSurfaceGouraud(pt[] P, float e, float ee){ 
  Boolean col=true;
  for(float s=0; s<1.001-e; s+=e) {
    col=!col;
    for(float t=0; t<1.001-e; t+=e) {
      if(col) fill(cyan); else fill(magenta); col=!col;
      beginShape(); 
      nv(CoonsNormal(P,s,t,ee)); v(coons(P,s,t)); 
      nv(CoonsNormal(P,s+e,t,ee)); v(coons(P,s+e,t)); 
      nv(CoonsNormal(P,s+e,t+e,ee)); v(coons(P,s+e,t+e)); 
      nv(CoonsNormal(P,s,t+e,ee)); v(coons(P,s,t+e)); 
      endShape(CLOSE);
      }
    }
  }

void showNormals(pt[] P, float e, float ee){ 
  for(float s=0; s<1.001-e; s+=e) for(float t=0; t<1.001-e; t+=e) show(coons(P,s,t),50,CoonsNormal(P,s,t,ee));
  }


void slide(pt[] P, float e) {
  float nfu=fu, nfv=fv;
  float z=coons(P,fu,fv).z;
  for(float a=0; a<=TWO_PI; a+=PI/20) {
    float nu=fu+e*cos(a), nv=fv+e*sin(a);
    float nz=coons(P,nu,nv).z;
    if(nz<z) {z=nz; nfu=nu; nfv=nv;}
    }
  fu=nfu; fv=nfv;
  }
  
void attractFront(pt[] P, float e) {
  float nfu=fu, nfv=fv;
  float z=d(coons(P,fu,fv),Of);
  for(float a=0; a<=TWO_PI; a+=PI/20) {
    float nu=fu+e*cos(a), nv=fv+e*sin(a);
    float nz=d(coons(P,nu,nv),Of);
    if(nz<z) {z=nz; nfu=nu; nfv=nv;}
    }
  fu=nfu; fv=nfv;
  }  

void attractBack(pt[] P, float r, float e) {
  float nbu=bu, nbv=bv;
  pt O = coons(PtQ.G,fu,fv);
  float z = abs(d(coons(PtQ.G,bu,bv),O)-r);
  for(float a=0; a<=TWO_PI; a+=PI/20) {
    float nu=bu+e*cos(a), nv=bv+e*sin(a);
    float nz=abs(d(coons(P,nu,nv),O)-r);
    if(nz<z) {z=nz; nbu=nu; nbv=nv;}
    }
  bu=nbu; bv=nbv;
  }  
  
void showMan(pt[] P, float h) {
  pt Of = coons(PtQ.G,fu,fv); vec Nf = CoonsNormal(PtQ.G,fu,fv,0.01); pt Ff = P(Of,5,Nf); pt Kf = P(Of,h,Nf); 
  pt Ob = coons(PtQ.G,bu,bv); vec Nb = CoonsNormal(PtQ.G,bu,bv,0.01); pt Fb = P(Ob,5,Nb); pt Kb = P(Ob,h,Nb); 
  float d=d(Kf,Kb)/2,b=sqrt(sq(h)-sq(d));
  vec V = V(Nf,Nb); pt B = P(P(Kf,Kb),b,V); 
  pt T = P(B,h,V);
  fill(red); show(Ff,5); collar(Of,V(h,Nf),1,5); show(Kf,5); collar(Kf,V(Kf,B),5,10);
  fill(blue); show(Fb,5); collar(Ob,V(h,Nb),1,5); show(Kb,5); collar(Kb,V(Kb,B),5,10);
  fill(orange); show(B,10); show(T,15); collar(B,V(B,T),10,15);
  }  
  
  
