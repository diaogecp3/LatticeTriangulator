/******************************************************************************
 * Frames, similarities, and range finder for steady lattices.
 *
 * Adapted from the version provided by Kelsey Kurzeja at Georgia Tech.
 ******************************************************************************/


class Frame3 {
  public pt p;
  public vec u, v, w;

  public Frame3() { this( P(0,0,0), V(1,0,0), V(0,1,0), V(0,0,1) ); }
  public Frame3(pt pp, vec pu, vec pv, vec pw) { p=pp; u=pu; v=pv; w=pw; }

  public Frame3 deepCopy() { return new Frame3(p.c(), u.c(), v.c(), w.c()); }

  public Frame3 set(Frame3 f) { p.set(f.p); u.set(f.u); v.set(f.v); w.set(f.w); return this; }
  public Frame3 set(pt pp, vec pu, vec pv, vec pw) { p.set(pp); u.set(pu); v.set(pv); w.set(pw); return this; }

  public pt uTip() { return p.c().add(u); }
  public pt vTip() { return p.c().add(v); }
  public pt wTip() { return p.c().add(w); }

  public boolean scalingIsUniform() {
    return isApproximately(u.norm(), v.norm(), .001) &&
           isApproximately(v.norm(), w.norm(), .001);  // Need to keep epsilon large, or else we will have numerical accuraccy issues.
  }
  public boolean hasSkew() {
    return !isApproximately(dot(u, v), 0, .0001) ||
           !isApproximately(dot(u, w), 0, .0001) ||
           !isApproximately(dot(v, w), 0, .0001);
  }
  public boolean hasReflection() { return dot(u, cross(v, w)) < 0; }
  public boolean isSimilarity() { return scalingIsUniform() && !hasSkew() && !hasReflection(); }  // Is this true?

  public pt toGlobalPoint(float x, float y, float z) { return p.c().add(u.c().mul(x)).add(v.c().mul(y)).add(w.c().mul(z)); }
  public pt toGlobalPoint(pt pp) { return toGlobalPoint(pp.x, pp.y, pp.z); }

  public vec toGlobalVector(float x, float y, float z) { return u.c().mul(x).add(v.c().mul(y)).add(w.c().mul(z)); }
  public vec toGlobalVector(vec vv) { return toGlobalVector(vv.x, vv.y, vv.z); }

  public Frame3 toGlobalFrame(pt pp, vec pu, vec pv, vec pw) { return new Frame3( toGlobalPoint(pp), toGlobalVector(pu), toGlobalVector(pv), toGlobalVector(pw) ); }
  public Frame3 toGlobalFrame(Frame3 f) { return toGlobalFrame(f.p, f.u, f.v, f.w); }

  public pt toLocalPoint(pt g) {
    float V = abs(signedTetVolume(u, v, w));
    float Vu = signedTetVolume( V(p,g), v, w );
    float Vv = signedTetVolume( V(p,g), w, u );
    float Vw = signedTetVolume( V(p,g), u, v );
    float x = Vu/V, y = Vv/V, z = Vw/V;
    return new pt(x, y, z);
  }
  public pt toLocalPoint(float gx, float gy, float gz) { return toLocalPoint(P(gx,gy,gz)); }

  public vec toLocalVector(vec g) { return toLocalPoint(p.c().add(g)).toVec(); }
  public vec toLocalVector(float gx, float gy, float gz) { return toLocalVector(V(gx,gy,gz)); }

  public Frame3 toLocalFrame(pt gp, vec gu, vec gv, vec gw) { return new Frame3(toLocalPoint(gp), toLocalVector(gu), toLocalVector(gv), toLocalVector(gw)); }
  public Frame3 toLocalFrame(Frame3 g) { return toLocalFrame(g.p, g.u, g.v, g.w); }

  public Frame3 transform(Frame3 transform, Frame3 reference) {
    set(reference.toGlobalFrame(transform.toGlobalFrame(reference.toLocalFrame(this))));
    return this;
  }
  public Frame3 transformLocal(Frame3 transform) {
    set(toGlobalFrame(transform));
    return this;
  }
  public Frame3 transformGlobal(Frame3 transform) {
    return transform(transform, new Frame3());
  }

  public Frame3 translateGlobal(float x, float y, float z) { p.add(x,y,z); return this; }
  public Frame3 translateGlobal(vec v) { p.add(v); return this; }
  public Frame3 scaleGlobal(float s) { p.mul(s); u.mul(s); v.mul(s); w.mul(s); return this; }
  public Frame3 scaleGlobal(float s, pt F) { p.sub(F); p.mul(s); p.add(F); u.mul(s); v.mul(s); w.mul(s); return this; }
  public Frame3 scaleGlobal(float sx, float sy, float sz) { p.x*=sx; p.y*=sy; p.z*=sz; u.x*=sx; u.y*=sy; u.z*=sz; v.x*=sx; v.y*=sy; v.z*=sz; w.x*=sx; w.y*=sy; w.z*=sz; return this; }
  public Frame3 rotateXGlobal(float rad) { p.rotX(rad); u.rotX(rad); v.rotX(rad); w.rotX(rad); return this; }  //**TODO: Can make faster by reusing the trig from each rotation
  public Frame3 rotateYGlobal(float rad) { p.rotY(rad); u.rotY(rad); v.rotY(rad); w.rotY(rad); return this; }
  public Frame3 rotateZGlobal(float rad) { p.rotZ(rad); u.rotZ(rad); v.rotZ(rad); w.rotZ(rad); return this; }
  public Frame3 rotateGlobal(float rad, vec axisDir) { p.rot(rad,axisDir); u.rot(rad,axisDir); v.rot(rad,axisDir); w.rot(rad,axisDir); return this; }
  public Frame3 rotateGlobal(float rad, vec axisDir, pt axisPoint) { p.rot(rad,axisDir,axisPoint); u.rot(rad,axisDir); v.rot(rad,axisDir); w.rot(rad,axisDir); return this; }

  public Frame3 translateLocal(float x, float y, float z) { p.add(u.c().mul(x)).add(v.c().mul(y)).add(w.c().mul(z)); return this; }
  public Frame3 translateLocal(vec v) { return translateLocal(v.x, v.y, v.z); }
  public Frame3 scaleLocal(float s) { u.mul(s); v.mul(s); w.mul(s); return this; }
  public Frame3 scaleLocal(float sx, float sy, float sz) { u.mul(sx); v.mul(sy); w.mul(sz); return this; }
  public Frame3 rotateLocal(float rad, vec axisDir) { u.rot(rad,axisDir); v.rot(rad,axisDir); w.rot(rad,axisDir); return this; }
  public Frame3 rotateXLocal(float rad) { rotateLocal(rad, u); return this; }
  public Frame3 rotateYLocal(float rad) { rotateLocal(rad, v); return this; }
  public Frame3 rotateZLocal(float rad) { rotateLocal(rad, w); return this; }

  public float[] matrixAsArray() { return new float[] { u.x,v.x,w.x,p.x, u.y,v.y,w.y,p.y, u.z,v.z,w.z,p.z }; }

  public String toString() { return p + "," + u + ", " + v + "," + w; }
}

enum SwirlType { INVALID, IDENTITY, TRANSLATION, ROTATION, SCALING, SCREW, SWIRL }

class SwirlTransform {
  public SwirlType type;
  public vec translation;
  public pt fixedPoint;  // In the case of rotation and screw, this is just a point on the axis of rotation.
  public vec axisDirection;
  public float rotationAngle;
  public float scaleFactor;

  public SwirlTransform() { type = SwirlType.IDENTITY; }

  public SwirlTransform deepCopy() {
    SwirlTransform n = new SwirlTransform();
    n.type = type;
    n.translation = (translation == null)? null : translation.c();
    n.fixedPoint = (fixedPoint == null)? null : fixedPoint.c();
    n.axisDirection = (axisDirection == null)? null : axisDirection.c();
    n.rotationAngle = rotationAngle;
    n.scaleFactor = scaleFactor;
    return n;
  }

  public pt transformPointNoCopy(pt p, float t) {
    if (type == SwirlType.IDENTITY) return p;
    else if (type == SwirlType.TRANSLATION) return p.add(translation.x*t, translation.y*t, translation.z*t);
    else if (type == SwirlType.ROTATION) return p.rot(rotationAngle*t, axisDirection, fixedPoint);
    else if (type == SwirlType.SCALING) return p.mul(pow(scaleFactor,t), fixedPoint);
    else if (type == SwirlType.SCREW) return p.rot(rotationAngle*t, axisDirection, fixedPoint).add(translation.x*t, translation.y*t, translation.z*t);
    else if (type == SwirlType.SWIRL) return p.rot(rotationAngle*t, axisDirection, fixedPoint).mul(pow(scaleFactor,t), fixedPoint);
    else return null;
  }

  public vec transformVectorNoCopy(vec v, float t) {
    if (type == SwirlType.IDENTITY) return v;
    else if (type == SwirlType.TRANSLATION) return v;
    else if (type == SwirlType.ROTATION) return v.rot(rotationAngle*t, axisDirection);
    else if (type == SwirlType.SCALING) return v.mul(pow(scaleFactor,t));
    else if (type == SwirlType.SCREW) return v.rot(rotationAngle*t, axisDirection);
    else if (type == SwirlType.SWIRL) return v.rot(rotationAngle*t, axisDirection).mul(pow(scaleFactor,t));
    else return null;
  }

  public Frame3 transformFrameNoCopy(Frame3 f, float t) {
    if (type == SwirlType.IDENTITY) return f;
    else if (type == SwirlType.TRANSLATION) return f.translateGlobal(translation.x*t, translation.y*t, translation.z*t);
    else if (type == SwirlType.ROTATION) return f.rotateGlobal(rotationAngle*t, axisDirection, fixedPoint);
    else if (type == SwirlType.SCALING) return f.scaleGlobal(pow(scaleFactor,t), fixedPoint);
    else if (type == SwirlType.SCREW) return f.rotateGlobal(rotationAngle*t, axisDirection, fixedPoint).translateGlobal(translation.x*t, translation.y*t, translation.z*t);
    else if (type == SwirlType.SWIRL) return f.rotateGlobal(rotationAngle*t, axisDirection, fixedPoint).scaleGlobal(pow(scaleFactor,t), fixedPoint);
    else return null;
  }

  public Ball transformBallNoCopy(Ball b, float t) {
    if (type == SwirlType.IDENTITY) return b;
    else if (type == SwirlType.TRANSLATION) return b.translate(translation.x*t, translation.y*t, translation.z*t);
    else if (type == SwirlType.ROTATION) return b.rotate(rotationAngle*t, axisDirection, fixedPoint);
    else if (type == SwirlType.SCALING) return b.scale(pow(scaleFactor,t), fixedPoint);
    else if (type == SwirlType.SCREW) return b.rotate(rotationAngle*t, axisDirection, fixedPoint).translate(translation.x*t, translation.y*t, translation.z*t);
    else if (type == SwirlType.SWIRL) return b.rotate(rotationAngle*t, axisDirection, fixedPoint).scale(pow(scaleFactor,t), fixedPoint);
    else return null;
  }

  public pt transformPoint(pt p, float t) { return transformPointNoCopy(p.c(), t); }
  public vec transformVector(vec v, float t) { return transformVectorNoCopy(v.c(), t); }
  public Frame3 transformFrame(Frame3 f, float t) { return transformFrameNoCopy(f.deepCopy(), t); }
  public Ball transformBall(Ball b, float t) { return transformBallNoCopy(b.deepCopy(), t); }

  public pt inversePointNoCopy(pt q, float t) { return transformPointNoCopy(q, -t); }
  public vec inverseVectorNoCopy(vec u, float t) { return transformVectorNoCopy(u, -t); }
  public Frame3 inverseFrameNoCopy(Frame3 f, float t) { return transformFrameNoCopy(f, -t); }
  public Ball inverseBallNoCopy(Ball b, float t) { return transformBallNoCopy(b, -t); }

  public pt inversePoint(pt q, float t) { return transformPointNoCopy(q.c(), -t); }
  public vec inverseVector(vec u, float t) { return transformVectorNoCopy(u.c(), -t); }
  public Frame3 inverseFrame(Frame3 f, float t) { return transformFrameNoCopy(f.deepCopy(), -t); }
  public Ball inverseBall(Ball b, float t) { return transformBallNoCopy(b.deepCopy(), -t); }
}

SwirlTransform MakeTranslation(vec translation) {
  SwirlTransform swirl = new SwirlTransform();
  swirl.type = SwirlType.TRANSLATION;
  swirl.translation = translation;
  return swirl;
}

SwirlTransform MakeRotation(float rad, vec axisDirection, pt pointOnAxis) {
  SwirlTransform swirl = new SwirlTransform();
  swirl.type = SwirlType.ROTATION;
  swirl.rotationAngle = rad;
  swirl.axisDirection = axisDirection;
  swirl.fixedPoint = pointOnAxis;
  return swirl;
}

SwirlTransform MakeScaling(float scaleFactor, pt fixedPoint) {
  SwirlTransform swirl = new SwirlTransform();
  swirl.type = SwirlType.SCALING;
  swirl.scaleFactor = scaleFactor;
  swirl.fixedPoint = fixedPoint;
  return swirl;
}

SwirlTransform MakeScrew(float rad, vec axisDirection, pt pointOnAxis, vec translation) {
  SwirlTransform swirl = new SwirlTransform();
  swirl.type = SwirlType.SCREW;
  swirl.rotationAngle = rad;
  swirl.axisDirection = axisDirection;
  swirl.fixedPoint = pointOnAxis;
  swirl.translation = translation;
  return swirl;
}

SwirlTransform MakeSwirl(float rad, vec axisDirection, pt fixedPoint, float scaleFactor) {
  SwirlTransform swirl = new SwirlTransform();
  swirl.type = SwirlType.SWIRL;
  swirl.rotationAngle = rad;
  swirl.axisDirection = axisDirection;
  swirl.fixedPoint = fixedPoint;
  swirl.scaleFactor = scaleFactor;
  return swirl;
}


//------------------------------------------------------------------------------
// RangeFinder
//------------------------------------------------------------------------------

// For now, assume no rotations > 360


//------------------------------------------------------------
// Normalized mappings (M)
//------------------------------------------------------------

float translationMap(pt Q, vec V) { return dot(Q,V) / V.magSq(); }
float scalingMap(pt Q, pt F, float s) { return logb(s, dist(F, Q)); }
float rotationMap(pt Q, vec O, pt F, vec V, float a, float b) {
  return angleAroundDirection(O, disp(F,Q), V)/a + TWO_PI*b/abs(a);
}


//------------------------------------------------------------
// Extents
//------------------------------------------------------------

class Extents {
  public MinMaxF translation, rotation, dilation;
  public Extents() { translation = rotation = dilation = null; }
}

boolean extentsAreSame(Extents A, Extents B) {
  float epsilon = .0001;
  if ( ((A.translation == null) != (B.translation == null)) || (A.translation != null && !isApproximately(A.translation, B.translation, epsilon)) ) return false;
  if ( ((A.rotation == null) != (B.rotation == null)) || (A.rotation != null && !isApproximately(A.rotation, B.rotation, epsilon)) ) return false;
  if ( ((A.dilation == null) != (B.dilation == null)) || (A.dilation != null && !isApproximately(A.dilation, B.dilation, epsilon)) ) return false;
  return true;
}


//------------------------------------------------------------
// Ball extents
//------------------------------------------------------------

MinMaxF translationExtents(pt Q, float q, vec V) {
  float c = translationMap(Q, V);
  float r = q / V.mag();
  return new MinMaxF(c-r, c+r);
}
MinMaxF translationExtents(Ball Q, vec V) { return translationExtents(Q.c, Q.r, V); }

MinMaxF scalingExtents(pt Q, float q, pt F, float s) {
  float d = dist(F, Q);
  float m1 = logb(s,d-q), m2 = logb(s,d+q);
  return new MinMaxF(min(m1,m2), max(m1,m2));
}
MinMaxF scalingExtents(Ball Q, pt F, float s) { return scalingExtents(Q.c, Q.r, F, s); }

MinMaxF rotationExtents(pt Q, float q, vec O, pt F, vec V, float a, float b) {
  float c = rotationMap(Q, O, F, V, a, b);
  float theta = asin(q/dist(Q, projectOnAxis(Q, F, V))) / abs(a);
  return new MinMaxF(c-theta, c+theta);
}
MinMaxF rotationExtents(Ball Q, vec O, pt F, vec V, float a, float b) { return rotationExtents(Q.c, Q.r, O, F, V, a, b); }


//------------------------------------------------------------
// Beam extents
//------------------------------------------------------------

MinMaxF translationExtents(Ball A, Ball B, vec V) {
  return bound(translationExtents(A, V), translationExtents(B, V));
}

MinMaxF scalingExtents(Ball A, Ball B, pt F, float s) {
  MinMaxF a = scalingExtents(A, F, s);
  MinMaxF b = scalingExtents(B, F, s);
  float d = logb(s, max(0, signedDistToBallCappedCone(F, A, B)));
  return new MinMaxF(min3(a.min, b.min, d), max3(a.max, b.max, d));
}

MinMaxF rotationExtents(Ball A, Ball B, vec O, pt F, vec V, float a, float b) {
  MinMaxF bound1 = rotationExtents(A, O, F, V, a, b);
  float c1 = (bound1.max + bound1.min) * 0.5;
  float delta = angleAroundDirection(disp(F, A.c), disp(F, B.c), V) / a;
  float c2 = c1 + delta;
  float theta2 = asin(B.r/dist(B.c,projectOnAxis(B.c, F, V))) / abs(a);
  MinMaxF bound2 = new MinMaxF(c2-theta2, c2+theta2);
  return bound(bound1, bound2);
}

Extents beamExtentsForTransform(SwirlTransform X, Ball A, Ball B) {
  Extents extents = new Extents();
  if (X.type == SwirlType.TRANSLATION || X.type == SwirlType.SCREW) {
    extents.translation = translationExtents(A, B, X.translation);
  }
  if (X.type == SwirlType.SCALING || X.type == SwirlType.SWIRL) {
    extents.dilation = scalingExtents(A, B, X.fixedPoint, X.scaleFactor);
  }
  if (X.type == SwirlType.ROTATION || X.type == SwirlType.SCREW || X.type == SwirlType.SWIRL) {
    extents.rotation = rotationExtents(A, B, perp(X.axisDirection), X.fixedPoint, X.axisDirection, X.rotationAngle, 0);
  }
  return extents;
}

boolean transformsAreSeparableForBeam(SwirlTransform X, SwirlTransform Y, Ball A, Ball B) {
  Ball XA = X.transformBall(A, 1), XB = X.transformBall(B, 1);
  Ball YA = Y.transformBall(A, 1), YB = Y.transformBall(B, 1);
  if (!extentsAreSame(beamExtentsForTransform(X, A, B), beamExtentsForTransform(X, YA, YB))) return false;  // Extents of beam(A,B) wrt X == extents of beam(YA,YB) wrt X
  if (!extentsAreSame(beamExtentsForTransform(Y, A, B), beamExtentsForTransform(Y, XA, XB))) return false;  // Extents of beam(A,B) wrt Y == extents of beam(XA,XB) wrt Y
  return true;
}

boolean transformsAreSeparableForBeam(SwirlTransform X, SwirlTransform Y, SwirlTransform Z, Ball A, Ball B) {
  return transformsAreSeparableForBeam(X, Y, A, B) &&
         transformsAreSeparableForBeam(X, Z, A, B) &&
         transformsAreSeparableForBeam(Y, Z, A, B);
}


//------------------------------------------------------------
// Primitive RangeFinders
//------------------------------------------------------------

MinMaxI translationRangeFinder(pt Q, float q, int n, MinMaxF baseExtents, vec V) {
  MinMaxF qExtents = translationExtents(Q, q, V);
  return new MinMaxI(max(0, ceil(qExtents.min-baseExtents.max)), min(n, floor(qExtents.max-baseExtents.min)));
}
MinMaxI translationRangeFinder(Ball Q, int n, MinMaxF baseExtents, vec V) {
  return translationRangeFinder(Q.c, Q.r, n, baseExtents, V);
}

MinMaxI scalingRangeFinder(pt Q, float q, int n, MinMaxF baseExtents, pt F, float s) {
  MinMaxF qExtents = scalingExtents(Q, q, F, s);
  return new MinMaxI(max(0, ceil(qExtents.min-baseExtents.max)), min(n, floor(qExtents.max-baseExtents.min)));
}
MinMaxI scalingRangeFinder(Ball Q, int n, MinMaxF baseExtents, pt F, float s) {
  return scalingRangeFinder(Q.c, Q.r, n, baseExtents, F, s);
}

MinMaxI[] rotationRangeFinder(pt Q, float q, MinMaxF baseExtents, vec O, pt F, vec V, float a, MinMaxI trimRange) {
  // baseExtents may be mapped into one or two branches: {0}, {-1,0}, {0,1}
  MinMaxF qExtents = rotationExtents(Q, q, O, F, V, a, 0);  // qExtents may be mapped into one or two branches: {0}, {-1,0}, {0,1}

  float branchSize = TWO_PI / a;
  int startBranch = (int)floor( (baseExtents.min+trimRange.min) / branchSize );
  int endBranch = (int)ceil( (baseExtents.max+trimRange.max) / branchSize );

  qExtents.min += startBranch*branchSize;
  qExtents.max += startBranch*branchSize;

  int numRanges = endBranch-startBranch+1;
  if (numRanges < 0) numRanges = 0;
  MinMaxI[] ranges = new MinMaxI[numRanges];

  for (int i = startBranch; i <= endBranch; i++) {
    ranges[i-startBranch] = new MinMaxI(max(trimRange.min, ceil(qExtents.min-baseExtents.max)), min(trimRange.max, floor(qExtents.max-baseExtents.min)));
    qExtents.min += branchSize;
    qExtents.max += branchSize;
  }
  return ranges;
}
MinMaxI[] rotationRangeFinder(Ball Q, MinMaxF baseExtents, vec O, pt F, vec V, float a, MinMaxI trimRange) {
  return rotationRangeFinder(Q.c, Q.r, baseExtents, O, F, V, a, trimRange);
}


//------------------------------------------------------------
// RangeFinder invocations
//------------------------------------------------------------

// Balls in range of ball
MinMaxI[] rangeFinder(pt Q, float q, int repetitions, SwirlTransform transform, pt C, float r) {
  SwirlTransform T = transform;
  if (T.type == SwirlType.TRANSLATION) {
    MinMaxF baseExtents = translationExtents(C,r,T.translation);
    return new MinMaxI[]{ translationRangeFinder(Q,q,repetitions,baseExtents,T.translation) };
  }
  else if (T.type == SwirlType.SCALING) {
    MinMaxF baseExtents = scalingExtents(C,r,T.fixedPoint,T.scaleFactor);
    return new MinMaxI[]{ scalingRangeFinder(Q,q,repetitions,baseExtents,T.fixedPoint,T.scaleFactor) };
  }
  else if (T.type == SwirlType.ROTATION) {
    vec O = perp(T.axisDirection);
    MinMaxF rotationBaseExtents = rotationExtents(C, r, O, T.fixedPoint, T.axisDirection, T.rotationAngle, 0);
    return rotationRangeFinder(Q, q, rotationBaseExtents, O, T.fixedPoint, T.axisDirection, T.rotationAngle, new MinMaxI(0,repetitions));
  }
  else if (T.type == SwirlType.SCREW) {
    MinMaxF translationBaseExtents = translationExtents(C,r,T.translation);
    MinMaxI translationRange = translationRangeFinder(Q,q,repetitions,translationBaseExtents,T.translation);

    vec O = perp(T.axisDirection);
    MinMaxF rotationBaseExtents = rotationExtents(C, r, O, T.fixedPoint, T.axisDirection, T.rotationAngle, 0);
    return rotationRangeFinder(Q, q, rotationBaseExtents, O, T.fixedPoint, T.axisDirection, T.rotationAngle, translationRange);
  }
  else if (T.type == SwirlType.SWIRL) {
    MinMaxF scalingBaseExtents = scalingExtents(C,r,T.fixedPoint,T.scaleFactor);
    MinMaxI scalingRange = scalingRangeFinder(Q,q,repetitions,scalingBaseExtents,T.fixedPoint,T.scaleFactor);

    vec O = perp(T.axisDirection);
    MinMaxF rotationBaseExtents = rotationExtents(C, r, O, T.fixedPoint, T.axisDirection, T.rotationAngle, 0);
    return rotationRangeFinder(Q, q, rotationBaseExtents, O, T.fixedPoint, T.axisDirection, T.rotationAngle, scalingRange);
  }
  return null;
}
MinMaxI[] rangeFinder(Ball Q, int repetitions, SwirlTransform transform, Ball A) { return rangeFinder(Q.c, Q.r, repetitions, transform, A.c, A.r); }
MinMaxI[] rangeFinder(pt Q, float q, int repetitions, SwirlTransform transform, Ball A) { return rangeFinder(Q, q, repetitions, transform, A.c, A.r); }
MinMaxI[] rangeFinder(Ball Q, int repetitions, SwirlTransform transform, pt C, float r) { return rangeFinder(Q.c, Q.r, repetitions, transform, C, r); }


// Beams in range of ball
MinMaxI[] rangeFinder(pt Q, float q, int repetitions, SwirlTransform transform, Ball A, Ball B) {
  SwirlTransform T = transform;
  if (T.type == SwirlType.TRANSLATION) {
    MinMaxF baseExtents = translationExtents(A, B, T.translation);
    return new MinMaxI[]{ translationRangeFinder(Q,q,repetitions,baseExtents,T.translation) };
  }
  else if (T.type == SwirlType.SCALING) {
    MinMaxF baseExtents = scalingExtents(A,B,T.fixedPoint,T.scaleFactor);
    return new MinMaxI[]{ scalingRangeFinder(Q,q,repetitions,baseExtents,T.fixedPoint,T.scaleFactor) };
  }
  else if (T.type == SwirlType.ROTATION) {
    vec O = perp(T.axisDirection);
    MinMaxF rotationBaseExtents = rotationExtents(A, B, O, T.fixedPoint, T.axisDirection, T.rotationAngle, 0);
    return rotationRangeFinder(Q, q, rotationBaseExtents, O, T.fixedPoint, T.axisDirection, T.rotationAngle, new MinMaxI(0,repetitions));
  }
  else if (T.type == SwirlType.SCREW) {
    MinMaxF translationBaseExtents = translationExtents(A,B,T.translation);
    MinMaxI translationRange = translationRangeFinder(Q,q,repetitions,translationBaseExtents,T.translation);

    vec O = perp(T.axisDirection);
    MinMaxF rotationBaseExtents = rotationExtents(A, B, O, T.fixedPoint, T.axisDirection, T.rotationAngle, 0);
    return rotationRangeFinder(Q, q, rotationBaseExtents, O, T.fixedPoint, T.axisDirection, T.rotationAngle, translationRange);
  }
  else if (T.type == SwirlType.SWIRL) {
    MinMaxF scalingBaseExtents = scalingExtents(A,B,T.fixedPoint,T.scaleFactor);
    MinMaxI scalingRange = scalingRangeFinder(Q,q,repetitions,scalingBaseExtents,T.fixedPoint,T.scaleFactor);

    vec O = perp(T.axisDirection);
    MinMaxF rotationBaseExtents = rotationExtents(A, B, O, T.fixedPoint, T.axisDirection, T.rotationAngle, 0);
    return rotationRangeFinder(Q, q, rotationBaseExtents, O, T.fixedPoint, T.axisDirection, T.rotationAngle, scalingRange);
  }
  return new MinMaxI[0];
}
MinMaxI[] rangeFinder(Ball Q, int repetitions, SwirlTransform transform, Ball A, Ball B) { return rangeFinder(Q.c, Q.r, repetitions, transform, A, B); }