/******************************************************************************
 * Projection of a vertex onto a hub or blended hub.
 ******************************************************************************/

private final int gMaxIterSphereTracing = 32;
private final float gUpperBound = 10000000;

/* Projection type. */
enum ProjectType {RAY, SPHERE_TRACING}

boolean debugProjection = false;

class DebugProjectionInfo {
  FloatList distsInside = new FloatList();
  FloatList tsInside = new FloatList();
  ArrayList<pt> pointsInside = new ArrayList<pt>();
  ArrayList<pt> targetsInside = new ArrayList<pt>();
  void reset() {
    distsInside.clear();
    tsInside.clear();
    pointsInside.clear();
    targetsInside.clear();
  }
  void printDists() {
    if (distsInside.size() == 0) return;
    print("Distances: ");
    for (float d : distsInside) print(d + ", ");
    println();
    int count = 0;
    for (int i = 0; i < distsInside.size(); ++i) {
      float d = distsInside.get(i);
      int j = 0;
      for (; j < i; ++j) {
        float f = distsInside.get(j);
        if (isZero(d - f)) break;
      }
      if (j == i) count++;
    }
    println("number of distinct distance =", count);
  }
  void printTs() {
    print("Parameters: ");
    for (float t : tsInside) {
      print(t);
      if (Float.isInfinite(t)) print("(Inf)");
      print(", ");
    }
    println();
  }
  void printPointsInside() {
    if (pointsInside.size() == 0) return;
    for (pt p : pointsInside) print(System.identityHashCode(p) + ", ");
    println();
    for (pt p : pointsInside) print(p + "; ");
    println();
  }
  void show() {
    for (int i = 0; i < pointsInside.size(); ++i) {
      pt p = pointsInside.get(i);
      pt q = targetsInside.get(i);
      fill(red);
      showBall(p, 1);
      fill(magenta);
      arrow(p, V(p, q), 1);
    }
  }
}

DebugProjectionInfo dProjectionInfo = new DebugProjectionInfo();

/*
 * Project a point on a hub by shooting a ray. p will be updated as its
 * projection. Only points outside the hub will be projected.
 */
void projectPointOnExactHub(Hub hub, pt p) {
  float dist = hub.distanceFrom(p);
  if (dist < 0.0001) return;  // skip if p is inside or on the hub
  vec d = U(p, hub.ball.c);
  Float t = hub.closestIntersectionWithRay(p, d);
  if (t != null && t < gUpperBound) p.add(t, d);
}

/*
 * Project a point on a hub by shooting a ray. p will be updated as its
 * projection. Points inside and outside the hub will be projected. Please use
 * this function only when the hub center is contained in the crudest convex
 * hull of the hub.
 */
void projectPointOnExactHubBiDir(Hub hub, pt p) {
  float dist = hub.distanceFrom(p);
  if (isZero(dist, 0.0001)) return;  // skip if p is on the hub
  vec d = U(p, hub.ball.c);
  if (dist < -0.0001) d.rev();  // flip the ray direction if p is inside the hub
  Float t = hub.closestIntersectionWithRay(p, d);
  if (t != null && t < gUpperBound) p.add(t, d);

  {
    // pt q = P(p);  // make a copy
    // if (t != null && t < gUpperBound) q.set(P(p, t, d));
    // if (debugProjection && dist < -0.0001) {
    //   dProjectionInfo.distsInside.append(dist);
    //   dProjectionInfo.tsInside.append(t);
    //   dProjectionInfo.pointsInside.add(p);
    //   dProjectionInfo.targetsInside.add(q);
    // }
    // p.set(q);
  }
}

/*
 * Project a point on the blended surface of a hub by shooting a ray. Because it
 * is difficult to compute the intersection between a point and the blended
 * surface, we use sphere tracing method to compute an approximation of the
 * intersection between a ray and the blended surface. Note that only points
 * outside the blended hub will be projected.
 */
void projectPointOnBlendedHub(Hub hub, pt p) {
  if (hub.blendedDistanceFrom(p) < 0.0001) return;  // skip if p is inside or on the blended hub
  vec d = U(p, hub.ball.c);
  pt q = P(p);  // make a copy
  for (int i = 0; i < gMaxIterSphereTracing; ++i) {
    float t = hub.blendedDistanceFrom(q);
    if (Float.isInfinite(t) || t > gUpperBound || t < 0.0001) break;
    q.add(t, d);
  }
  p.set(q);
}

/*
 * Project a point on the blended surface of a hub by shooting a ray. We use
 * sphere tracing method to compute an approximation of the intersection between
 * a ray and the blended surface. Note that both points inside and outside the
 * blended hub will be projected. Please make sure the hub center is contained
 * in the crudest convex hull of the hub. Warning: this function may not work
 * well.
 */
void projectPointOnBlendedHubBiDir(Hub hub, pt p) {
  float dist = hub.blendedDistanceFrom(p);
  if (isZero(dist, 0.0001)) return;  // skip if p is on the blended hub
  vec d = U(p, hub.ball.c);
  if (dist < -0.0001) d.rev();  // flip the ray direction if p is inside the blended hub
  pt q = P(p);
  for (int i = 0; i < gMaxIterSphereTracing; ++i) {
    float t = hub.blendedDistanceFrom(q);
    if (Float.isInfinite(t) || t > gUpperBound || t < 0.0001) break;
    q.add(t, d);
  }
  p.set(q);
}