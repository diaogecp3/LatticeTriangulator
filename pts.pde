/******************************************************************************
 * Points.
 *
 * Adapted from the version provided by Prof. Jarek Rossignac at Georgia Tech.
 ******************************************************************************/

int pp = 1;  // index of picked vertex

class pts {  // class for manipulaitng and sisplaying polyloops
  boolean loop = true;
  int pv = 0;  // picked vertex index,
  int iv = 0;  // insertion vertex index
  int nv = 0;  // number of vertices currently used in P
  int maxnv = 16000;  //  max number of vertices
  pt[] G = new pt[maxnv];  // geometry table (vertices)
  pts() {}
  pts declare() {for (int i=0; i<maxnv; i++) G[i]=P(); return this;}     // init all point objects
  pts empty() {nv=0; pv=0; return this;} // resets P so that we can start adding points
  pts addPt(pt P) { G[nv].set(P); pv=nv; nv++;  return this;} // adds a point at the end
  pts addPt(float x,float y) { G[nv].x=x; G[nv].y=y; pv=nv; nv++; return this;}
  pts copyFrom(pts Q) {empty(); nv=Q.nv; for (int v=0; v<nv; v++) G[v]=P(Q.G[v]); return this;}
  pts setToL(pts P, float t, pts Q) { // lerp (linear interpolation betwen P and Q
    empty();
    nv=min(P.nv,Q.nv);
    for (int v=0; v<nv; v++) G[v]=L(P.G[v],t,Q.G[v]);
    return this;}
  pts resetOnCircle(int k, float r) { // makes new polyloo[p with k  points on a circle around origin
    empty(); // resert P
    pt C = P(); // center of circle
    for (int i=0; i<k; i++) addPt(R(P(C,V(0,-r,0)),2.*PI*i/k,C)); // points on z=0 plane
    pv=0; // picked vertex ID is set to 0
    return this;
    }
  pts projectOnSphere(float r) { // makes new polyloo[p with k  points on a circle around origin
    for (int i=0; i<nv; i++) G[i].snapToSphere(r); // points on z=0 plane
    pv=0; // picked vertex ID is set to 0
    return this;
    }

  int idOfVertexWithClosestScreenProjectionTo(pt M) { // for picking a vertex with the mouse
    pp=0;
    for (int i=1; i<nv; i++) if (d(M,ToScreen(G[i]))<=d(M,ToScreen(G[pp]))) pp=i;
    return pp;
    }

  pt closestProjectionOf(pt M) {   // for picking inserting O. Returns projection but also CHANGES iv !!!!
    pt C = P(G[0]); float d=d(M,C);
    for (int i=1; i<nv; i++) if (d(M,G[i])<=d) {iv=i; C=P(G[i]); d=d(M,C); }
    for (int i=nv-1, j=0; j<nv; i=j++) {
       pt A = G[i], B = G[j];
       if (projectsBetween(M,A,B) && disToLine(M,A,B)<d) {d=disToLine(M,A,B); iv=i; C=projectionOnLine(M,A,B);}
       }
    return C;
    }
  void setPickToIndexOfVertexClosestTo(pt M) {   // for picking inserting O. Returns projection but also CHANGES iv !!!!
    pv=0;
    for (int i=1; i<nv; i++) if (d(M,G[i])<=d(M,G[pv])) pv=i;
    }

  pts insertPt(pt P) { // inserts new vertex after vertex with ID iv
    for(int v=nv-1; v>iv; v--) G[v+1].set(G[v]);
     iv++;
     G[iv].set(P);
     nv++; // increments vertex count
     return this;
     }

  pts insertClosestProjection(pt M) {
    pt P = closestProjectionOf(M); // also sets iv
    insertPt(M);
    return this;
    }

  pts deletePicked() {for(int i=pv; i<nv; i++) G[i].set(G[i+1]); pv=max(0,pv-1); nv--;  return this;}
  pts deletePickedPair() {
    if (nv == 0 || nv % 2 == 1) return this;
    int j = pv - pv % 2;
    for (int i = j; i < nv - 2; ++i) {
      G[i].set(G[i+2]);
    }
    pv = max(0, j - 1);
    nv -= 2;
    return this;
  }
  pts setPt(pt P, int i) { G[i].set(P); return this;}
  pts setPickedTo(pt P) { G[pv].set(P); return this;}
  pts showPicked() {show(G[pv],13); return this;}
  pts drawBalls(float r) {for (int v=0; v<nv; v++) show(G[v],r); return this;}
  pts showPicked(float r) {show(G[pv],r); return this;}
  pts drawClosedCurve(float r) {for (int v=0; v<nv-1; v++) stub(G[v],V(G[v],G[v+1]),r,r/2);  stub(G[nv-1],V(G[nv-1],G[0]),r,r/2); return this;}
  pts setPickedIndexTo(int pp) {pv=pp; return this;}
  pts movePicked(vec V) { G[pv].add(V); return this;}      // moves selected point (index p) by amount mouse moved recently
  pts movePickedTo(pt P) { G[pv].set(P); return this;}      // moves selected point (index p) by amount mouse moved recently
  pts moveAll(vec V) {for (int i=0; i<nv; i++) G[i].add(V); return this;};
  pt Picked() {return G[pv];}

  void savePts(String fn) {
    println("saving point set:", fn);
    String [] inppts = new String [nv+1];
    int s=0;
    inppts[s++]=str(nv);
    for (int i=0; i<nv; i++) {inppts[s++]=str(G[i].x)+","+str(G[i].y)+","+str(G[i].z);}
    saveStrings(fn,inppts);
  };

  void loadPts(String fn) {
    print("loading: " + fn + " ");
    String [] ss = loadStrings(fn);
    int s = 0;
    nv = int(ss[s++]);
    println("nv =", nv);
    for (int k = 0; k < nv; k++) {
      int i = k + s;
      float[] xyz = float(split(ss[i], ","));
      G[k].set(xyz[0],xyz[1],xyz[2]);
    }
    pv = 0;
  };

} // end of pts class
