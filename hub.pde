/******************************************************************************
 * Hub (a ball and tangential cones) processing.
 ******************************************************************************/

boolean showHub = true;
boolean showBoundingSphere = false;  // the bounding sphere of a hub
boolean showIntersectionCircles = false;  // the intersection circles between the bounding sphere and cones
boolean showLiftedCones = false;
boolean showGapMesh = false;


class Cone {
  pt apex;  // apex
  vec axis;  // in the direction of incresing radius
  float alpha;  // half angle

  Cone (pt apex, vec axis, float alpha) {
    this.apex = apex;
    this.axis = axis;
    this.alpha = alpha;
  }

  /*
   * Compute the intersections between a line and a two-cone. Return the two
   * parameters corresponding to the two intersections, if any.
   */
  float[] intersectLine(pt o, vec d) {
    float dv = dot(d, axis);
    vec co = V(apex, o);
    float cov = dot(co, axis);
    float cos = cos(alpha);
    float cos2 = cos * cos;

    float a = dv * dv - cos2;
    float b = 2 * (dv * cov - dot(d, co) * cos2);
    float c = cov * cov - dot(co, co) * cos2;

    float[] sols = solveQuadraticEquation(a, b, c);
    return sols;
  }
}

class TruncatedCone {
  pt c0;  // center of the lower base
  float r0;  // radius of the lower base
  pt c1;  // center of the upper base
  float r1;  // radius of the upper base

  vec normal;  // unit vector from c0 to c1
  Cone cone = null;  // the original cone
  float lowHeight = -1.0;
  float highHeight = -1.0;

  pt[] samples = null;

  TruncatedCone(pt c0, float r0, pt c1, float r1) {
    this.c0 = c0;
    this.r0 = r0;
    this.c1 = c1;
    this.r1 = r1;
    normal = U(c0, c1);
  }

  TruncatedCone(pt c0, float r0, pt c1, float r1, vec normal, Cone cone) {
    this.c0 = c0;
    this.r0 = r0;
    this.c1 = c1;
    this.r1 = r1;
    this.normal = normal;
    this.cone = cone;
  }

  /*
   * Compute the original cone. This function is called when needed.
   */
  void initCone() {
    if (cone != null) return;
    float d = d(c0, c1);
    float tan = abs(r1 - r0) / d;
    float alpha = atan(tan);  // [0, PI/2)
    vec axis = normal;
    if (r1 < r0) axis = M(normal);  // reverse the direction
    float x = r0 / tan;
    pt apex = P(c0, -x, axis);
    cone = new Cone(apex, axis, alpha);

    if (r0 < r1) {
      lowHeight = x;
      highHeight = x + d;
    } else {
      lowHeight = x - d;
      highHeight = x;
    }
  }

  void generateSamples() {
    int n = numPointsPerRing;  // global variable
    samples = new pt[2 * n];
    vec vi = constructNormal(normal);
    vec vj = N(normal, vi);
    float da = TWO_PI / n;
    float a = 0;
    for (int i = 0; i < n; ++i) {
      float cos = cos(a);
      float sin = sin(a);
      samples[i] = P(c0, r0 * cos, vi, r0 * sin, vj);
      samples[i + n] = P(c1, r1 * cos, vi, r1 * sin, vj);
      a += da;
    }
  }

  /*
   * Find the closest valid intersection with line (o, d). The intersection
   * point must be on the positive cone and has height in [lowHeight, highHeight].
   * If there are two valid, return the closer one.
   */
  Float closestIntersectionWithLine(pt o, vec d) {
    assert cone != null;
    float[] ts = cone.intersectLine(o, d);
    if (ts == null) return null;

    pt p0 = P(o, ts[0], d);
    pt p1 = P(o, ts[1], d);

    float d0 = dot(V(cone.apex, p0), cone.axis);
    float d1 = dot(V(cone.apex, p1), cone.axis);
    boolean valid0 = (d0 >= lowHeight && d0 <= highHeight);
    boolean valid1 = (d1 >= lowHeight && d1 <= highHeight);

    if (valid0) {
      if (valid1 && abs(ts[1]) < abs(ts[0])) return new Float(ts[1]);
      else return new Float(ts[0]);
    } else {
      if (valid1) return new Float(ts[1]);
      else return null;
    }
  }

  /*
   * Find the closest valid intersection with ray (o, d). The intersection
   * point must be on the positive cone and has height in [lowHeight, highHeight].
   * If there are two valid, return the closer one.
   */
  Float closestIntersectionWithRay(pt o, vec d) {
    assert cone != null;
    float[] ts = cone.intersectLine(o, d);
    if (ts == null) return null;

    pt p0 = ts[0] > 0 ? P(o, ts[0], d) : null;
    pt p1 = ts[1] > 0 ? P(o, ts[1], d) : null;

    boolean valid0 = false;
    if (p0 != null) {
      float d0 = dot(V(cone.apex, p0), cone.axis);
      valid0 = (d0 >= lowHeight && d0 <= highHeight);
    }

    boolean valid1 = false;
    if (p1 != null) {
      float d1 = dot(V(cone.apex, p1), cone.axis);
      valid1 = (d1 >= lowHeight && d1 <= highHeight);
    }

    if (valid0) {
      if (valid1 && ts[1] < ts[0]) return new Float(ts[1]);
      else return new Float(ts[0]);
    } else {
      if (valid1) return new Float(ts[1]);
      else return null;
    }
  }

  void show() {
    if (samples == null || samples.length == 0) generateSamples();
    int n = samples.length / 2;
    stroke(0);
    beginShape(QUAD_STRIP);
    for (int i = 0; i < n; ++i) {
      vertex(samples[n + i]);
      vertex(samples[i]);
    }
    vertex(samples[n]);
    vertex(samples[0]);
    endShape();
    noStroke();
  }
}

class ConicGap {
  ArrayList<pt> points0;  // points on inner irregular loop
  ArrayList<pt> points1;  // points on outer regular loop

  ConicGap(ArrayList<pt> points0, ArrayList<pt> points1) {
    this.points0 = points0;
    this.points1 = points1;
  }

  ArrayList<pt> positionList() {
    ArrayList<pt> posList = new ArrayList<pt>();
    posList.addAll(points0);
    posList.addAll(points1);
    return posList;
  }

  ArrayList<Triangle> gapHull() {
    int nv0 = points0.size();
    int nv1 = points1.size();
    assert nv0 >= 3 && nv1 >= 3;

    ArrayList<Triangle> triangles = new ArrayList<Triangle>();

    /* Find the first triangle. */
    pt pa = points0.get(0);
    pt pb = points0.get(1);
    int j = 0;
    int jLeft = nv1 - 1, jRight;
    for (; j < nv1; ++j) {
      jRight = (j + 1) % nv1;
      vec vLeft = V(pa, points1.get(jLeft));
      vec vRight = V(pa, points1.get(jRight));
      vec n = N(pa, pb, points1.get(j));  // no need to normalize
      if (dot(n, vLeft) < 0 && dot(n, vRight) < 0) break;
      jLeft = j;
    }
    assert j < nv1;
    triangles.add(new Triangle(0, 1, j + nv0));

    /* Traverse the loops. */
    int stop0 = 0;
    int stop1 = j;
    int i = 1;
    int iNext = (i + 1) % nv0;
    int jNext = (j + 1) % nv1;
    int jCount = 0;
    while (i != stop0 && jCount < nv1) {
      pa = points1.get(j);
      pb = points0.get(i);
      pt pc = points0.get(iNext);

      /* Check convexity. */
      boolean valid = false;
      {
        vec n = N(pa, pb, pc);
        vec vLeft = V(pa, points1.get((j + nv1 - 1) % nv1));
        vec vRight = V(pa, points1.get(jNext));
        if (dot(n, vLeft) < 0 && dot(n, vRight) < 0) valid = true;
      }

      if (valid) {
        triangles.add(new Triangle(j + nv0, i, iNext));
        i = iNext;
        iNext = (i + 1) % nv0;
      } else {
        triangles.add(new Triangle(j + nv0, i, jNext + nv0));
        j = jNext;
        jNext = (j + 1) % nv1;
        jCount++;
      }
    }

    while (i != stop0) {
      triangles.add(new Triangle(j + nv0, i, iNext));
      i = iNext;
      iNext = (i + 1) % nv0;
    }

    while (j != stop1) {
      triangles.add(new Triangle(j + nv0, i, jNext + nv0));
      j = jNext;
      jNext = (j + 1) % nv1;
    }
    return triangles;
  }

  TriangleMesh toTriMesh() {
    ArrayList<pt> posList = positionList();
    ArrayList<Triangle> triList = gapHull();
    return new TriangleMesh(posList, triList);
  }
}


/*
 * Hub class.
 *
 * A hub is the union of a ball and a set of tangential cones
 * defined by it and each of its neighboring balls.
 */
class Hub {
  Ball ball = null;  // inner ball
  Ball[] neighbors = null;  // outer balls
  int nNeighbors = 0;

  boolean valid = false;  // whether the hub is valid or not
  float maxIntersectDist = -1.0;  // the radius of the bounding sphere

  TruncatedCone[] tCones = null;  // beams
  Circle[] circles = null;  // intersection between cones and bounding sphere

  float gapWidth = 10.0;  // the width of the gap between the convex hull and beams
  TruncatedCone[] liftedCones = null;

  ConicGap[] gaps = null;

  Hub() {}

  Hub(Ball ball, Ball[] neighbors, int nNeighbors) {
    this.ball = ball;
    this.neighbors = neighbors;
    this.nNeighbors = nNeighbors;

    valid = isValid();
    if (valid) this.init();
  }

  Hub(pt c, float r, pt[] points, int nv) {
    ball = new Ball(c, r);
    nNeighbors = nv / 2;
    neighbors = new Ball[nNeighbors];
    for (int i = 0; i < nv; i += 2) {
      neighbors[i / 2] =  new Ball(points[i], d(points[i], points[i+1]));
    }

    valid = isValid();
    if (valid) this.init();
  }

  private boolean isValid() {
    for (int i = 0; i < nNeighbors - 1; ++i) {
      for (int j = i + 1; j < nNeighbors; ++j) {
        if (neighbors[i].intersectBall(neighbors[j])) return false;
      }
      if (ball.intersectBall(neighbors[i])) return false;
    }
    return true;
  }

  private void generateTruncatedCones() {
    tCones = new TruncatedCone[nNeighbors];
    for (int i = 0; i < nNeighbors; ++i) {
      vec v = V(ball.c, neighbors[i].c);
      float d = v.norm();
      v.div(d);  // normalize
      float cos = (ball.r - neighbors[i].r) / d;
      float sin = sin(acos(cos));
      pt c0 = P(ball.c, ball.r * cos, v);  // the center corresponding to the inner ball
      pt c1 = P(neighbors[i].c, neighbors[i].r * cos, v);  // the center corrsponding to the i-th outer ball
      float r0 = ball.r * sin;
      float r1 = neighbors[i].r * sin;

      tCones[i] = new TruncatedCone(c0, r0, c1, r1);
      tCones[i].initCone();
    }
  }

  private float intersectionDistance(int i, int j) {
    assert tCones[i] != null && tCones[j] != null;
    pt pa = ball.c;
    pt pb = neighbors[i].c;
    pt pc = neighbors[j].c;

    vec vba = U(pb, pa);  // x-axis
    vec vca = U(pc, pa);
    vec normal = U(N(vba, vca));
    vec vj = N(normal, vba);  // y-axis = rotated vba
    vec vcaRotated = R(vca, HALF_PI, vba, vj);

    pt pba10 = P(tCones[i].c1, -tCones[i].r1, vj);
    pt pba11 = P(tCones[i].c0, -tCones[i].r0, vj);
    pt pca00 = P(tCones[j].c1, tCones[j].r1, vcaRotated);
    pt pca01 = P(tCones[j].c0, tCones[j].r0, vcaRotated);

    pt px = intersectionTwoLines(pba10, pba11, pca00, pca01);

    {
      // fill(green);
      // show(px, 4);
    }

    return d(pa, px);
  }

  private float maximumIntersectionDistance() {
    float t = -1.0;
    for (int i = 0; i < nNeighbors - 1; ++i) {
      for (int j = i + 1; j < nNeighbors; ++j) {
        float d = intersectionDistance(i, j);
        t = max(t, d);
      }
    }
    return t;
  }

  private void init() {
    generateTruncatedCones();

    maxIntersectDist = maximumIntersectionDistance();
    maxIntersectDist += 10;  // slightly bigger such that circles are disjoint
  }

  private TruncatedCone liftCone(int i, pt c, float r) {
    pt c0 = P(c, gapWidth, tCones[i].normal);
    float r0 = r;
    if (tCones[i].r1 < tCones[i].r0) {
      r0 -= gapWidth * tan(tCones[i].cone.alpha);
    } else {
      r0 += gapWidth * tan(tCones[i].cone.alpha);
    }
    return new TruncatedCone(c0, r0, tCones[i].c1, tCones[i].r1, tCones[i].normal, tCones[i].cone);
  }

  private Circle intersectionCircleWithBoundingSphere(int i) {
    float r0 = ball.r * ball.r / tCones[i].r0;  // the radius of the base of ball.c
    float r1 = neighbors[i].r * neighbors[i].r / tCones[i].r1;  // the radius of the base of neighbors[i].c
    vec v = V(ball.c, neighbors[i].c);

    float dd = n2(v);
    float r0r0 = r0 * r0;
    float r1r1 = r1 * r1;
    float r0r1 = 2 * r0 * r1;

    float a = dd + r0r0 + r1r1 - r0r1;
    float b = -2 * r0r0 + r0r1;
    float c = r0r0 - maxIntersectDist * maxIntersectDist;

    float[] sols = solveQuadraticEquation(a, b, c);
    float t = -1.0;  // t should be in [0, 1]
    if (sols != null) {
      t = max(sols[0], sols[1]);
    }

    if (t >= 0 && t <= 1) {
      float r = (1 - t) * r0 + t * r1;
      pt center = P(ball.c, t * sqrt(dd), tCones[i].normal);
      return new Circle(center, tCones[i].normal, r);
    }

    return null;
  }

  void intersectionCircles() {
    if (maxIntersectDist <= 0.0) return;

    circles = new Circle[nNeighbors];
    liftedCones = new TruncatedCone[nNeighbors];
    for (int i = 0; i < nNeighbors; ++i) {
      circles[i] = intersectionCircleWithBoundingSphere(i);
      liftedCones[i] = liftCone(i, circles[i].c, circles[i].r);
    }
    return;
  }

  RingSet circlesToRingSet() {
    return new RingSet(ball.c, maxIntersectDist, circles, nNeighbors);
  }

  /*
   * Return the distance from point p to the hub.
   */
  float distanceFrom(pt p) {
    float d = Float.MAX_VALUE;
    for (int i = 0; i < nNeighbors; ++i) {
      d = min(d, roundConeDist(p, ball.c, ball.r, neighbors[i].c, neighbors[i].r));
    }
    return d;
  }

  /*
   * Return the distance from point p to the blended hub. The blending is
   * independent on the order of the neighbors.
   */
  float blendedDistanceFrom(pt p) {
    float d = 0;
    float k = 8;  // shouldn't be too large
    for (int i = 0; i < nNeighbors; ++i) {
      d += exp(-k * roundConeDist(p, ball.c, ball.r, neighbors[i].c, neighbors[i].r));
    }
    return -log(d) / k;
  }

  float closestIntersectionWithLine(pt o, vec d) {
    /* Check intersections between each cone and the line. */
    float t = Float.MAX_VALUE;
    for (int i = 0; i < nNeighbors; ++i) {
      Float tmp = tCones[i].closestIntersectionWithLine(o, d);
      if (tmp != null && abs(tmp) < abs(t)) t = tmp;
    }

    {  // check intersections between the inner ball and the line
      float[] ts = ball.intersectLine(o, d);
      if (ts != null) {
        if (abs(ts[0]) < abs(t)) t = ts[0];
        if (abs(ts[1]) < abs(t)) t = ts[1];
      }
    }

    return t;
  }

  float closestIntersectionWithRay(pt o, vec d) {
    /* Check intersections between each cone and the ray. */
    float t = Float.MAX_VALUE;
    for (int i = 0; i < nNeighbors; ++i) {
      Float tmp = tCones[i].closestIntersectionWithRay(o, d);
      if (tmp != null && tmp < t) t = tmp;
    }

    {  // check intersections between the inner ball and the ray
      float[] ts = ball.intersectLine(o, d);
      if (ts != null) {
        if (ts[0] > 0 && ts[0] < t) t = ts[0];
        if (ts[1] > 0 && ts[1] < t) t = ts[1];
      }
    }

    return t;
  }

  /*
   * Triangulate the i-th beam. The ID of the first point of this beam is k.
   * Augment the triangle mesh tm with new positions and new triangles. Return
   * the number of new positions.
   */
  private int triangulateBeam(int b, int k, TriangleMesh tm) {
    TruncatedCone beam = liftedCones[b];
    assert beam != null;
    if (beam.samples == null) beam.generateSamples();
    pt[] samples = beam.samples;

    ArrayList<pt> positions = new ArrayList<pt>();
    for (int i = 0; i < samples.length; ++i) {
      positions.add(samples[i]);
    }
    ArrayList<Triangle> triangles = new ArrayList<Triangle>();
    int n = samples.length / 2;
    for (int i = 0; i < n - 1; ++i) {
      int ia = i + k;
      int ib = ia + 1;
      int ic = ib + n;
      int id = ia + n;
      triangles.add(new Triangle(ia, ib, ic));
      triangles.add(new Triangle(ic, id, ia));
    }
    {
      int ia = n - 1 + k;
      int ib = k;
      int ic = ib + n;
      int id = ia + n;
      triangles.add(new Triangle(ia, ib, ic));
      triangles.add(new Triangle(ic, id, ia));
    }
    tm.augmentWithoutShift(positions, triangles);

    return samples.length;
  }

  /*
   * Triangulate the beams.
   */
  TriangleMesh triangulateBeams() {
    TriangleMesh tm = new TriangleMesh();
    int k = 0;
    for (int i = 0; i < nNeighbors; ++i) {
      k += triangulateBeam(i, k, tm);
    }
    return tm;
  }

  void initGaps(ArrayList<pt>[] innerLoops) {
    assert innerLoops.length == nNeighbors;
    assert liftedCones != null;
    gaps = new ConicGap[nNeighbors];

    for (int i = 0; i < nNeighbors; ++i) {
      ArrayList<pt> innerLoop = innerLoops[i];
      if (liftedCones[i].samples == null) liftedCones[i].generateSamples();
      ArrayList<pt> outerLoop = new ArrayList<pt>();
      pt[] samples = liftedCones[i].samples;
      int n = samples.length / 2;  // or numPointsPerRings
      for (int j = 0; j < n; ++j) outerLoop.add(samples[j]);
      gaps[i] = new ConicGap(innerLoop, outerLoop);
    }
  }

  TriangleMesh gapMesh() {
    if (gaps == null || gaps.length == 0) return null;

    TriangleMesh gapMesh = new TriangleMesh();
    for (int i = 0; i < nNeighbors; ++i) {
      TriangleMesh tm = gaps[i].toTriMesh();
      gapMesh.augmentWithShift(tm.positions, tm.triangles);
    }
    return gapMesh;
  }

  void showIntersectionCircles() {
    stroke(cyan);
    for (int i = 0; i < nNeighbors; ++i) {
      if (circles[i] != null) circles[i].show();
    }
  }

  void showBoundingSphere(color c, int a) {
    if (maxIntersectDist > 0) {
      fill(c, a);
      noStroke();
      show(ball.c, maxIntersectDist);
    }
  }

  void showLiftedCones() {
    for (int i = 0; i < nNeighbors; ++i) {
      if (liftedCones[i] != null) liftedCones[i].show();
    }
  }

  void showHub(color c, int a) {
    fill(c, a);
    noStroke();
    ball.showBall();
    for (int i = 0; i < nNeighbors; ++i) {
      tCones[i].show();
      // neighbors[i].showBall();
    }
    return;
  }

  void save(String file) {
    String[] lines = new String[2 + nNeighbors];
    int i = 0;
    lines[i++] = str(nNeighbors + 1);
    lines[i++] = str(ball.c.x) + "," + str(ball.c.y) + "," + str(ball.c.z) +
                 "," + str(ball.r);
    for (int j = 0; j < nNeighbors; ++j) {
      lines[i++] = str(neighbors[j].c.x) + "," + str(neighbors[j].c.y) + "," +
                   str(neighbors[j].c.z) + "," + str(neighbors[j].r);
    }
    saveStrings(file, lines);
    return;
  }

  void load(String file) {
    String[] lines = loadStrings(file);
    int i = 0;
    nNeighbors = int(lines[i++]) - 1;
    println("loading:", file, "number of neighbors =", nNeighbors);
    float[] tmp = float(split(lines[i++], ","));
    ball = new Ball(tmp[0], tmp[1], tmp[2], tmp[3]);
    neighbors = new Ball[nNeighbors];
    for (int j = 0; j < nNeighbors; ++j) {
      tmp = float(split(lines[i++], ","));
      neighbors[j] = new Ball(tmp[0], tmp[1], tmp[2], tmp[3]);
    }
    return;
  }
}
