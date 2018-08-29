void drawArc(pt A, pt C, pt B, float sr, float w, boolean small) { // draws arc on sphere of radius sr and center C from A to B using weignt w
  if(small) drawArc(A,C,B,sr,w);
  else {
    vec M = M(Slerp(V(C,A),0.5,V(C,B)));
    drawArc(A,C,P(C,M),sr,w);
    drawArc(P(C,M),C,B,sr,w);
    }
  }
  
void drawArc(pt A, pt C, pt B, float sr, float w) { // draws arc on sphere of radius sr and center C from A to B using weignt w
  float a = angle(U(C,A),U(C,B));
  pt PP = P(A);
  for(float t=0.05; t<=1.01; t+=0.05) {pt QQ = onArc(A,C,B,t); show(QQ,w); stub(PP,QQ,w,w); PP.setTo(QQ);}
  }
  
pt onArc(pt A, pt C, pt B, float t) { // draws point on arc of center C between A and B 
//  vec U = U(V(C,A)), W = U(V(C,B)), V=Slerp(U,t,W);
//  float ca=d(C,A), cb=d(C,B), m=cb/ca;
//  return P(C,ca*pow(m,t),V);
  vec U = V(C,A), W = V(C,B);
  vec V=Slerp(U,t,W); 
  return P(C,V);
  }
  
vec Slerp(vec U, float t, vec W) {
  float a = angle(U,W);
  float su=sin(a*(1.-t)), sw=sin(a*t), s = sin(a);
  return V(su/s,U,sw/s,W);
  } 
