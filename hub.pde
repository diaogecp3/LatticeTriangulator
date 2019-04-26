/******************************************************************************
 * Hub (a ball and tangential cones).
 ******************************************************************************/

boolean showHub = true;
boolean showBoundingSphere = false;  // the bounding sphere of a hub
boolean showIntersectionCircles = false;  // the intersection circles between the bounding sphere and cones
boolean showLiftedCones = false;
boolean showGapMesh = true;

float gGapWidth = 10;

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

  float h0 = 0.0;
  float h1 = 0.0;

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

  void halve() {
    c1 = P(c0, c1);
    r1 = (r0 + r1) / 2;
    h1 = (h0 + h1) / 2;
  }

  /*
   * Compute the original cone. This function is called when needed.
   */
  void initCone() {
    assert cone == null;
    float d = d(c0, c1);
    float tan = abs(r1 - r0) / d;
    float alpha = atan(tan);  // [0, PI/2)
    vec axis = normal;
    if (r1 < r0) axis = M(normal);  // reverse the direction
    float x = r0 / tan;
    pt apex = P(c0, -x, axis);
    cone = new Cone(apex, axis, alpha);

    h0 = x;
    h1 = (r0 < r1) ? x + d : x - d;

    assert h0 > 0.0 && h1 > 0.0;
  }

  void generateSamples(int numSamples, vec xAxis) {
    int n = numSamples;  // the number of samples on a base
    samples = new pt[2 * n];
    vec vi = xAxis == null ? constructNormal(normal) : xAxis;
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
    float lowHeight = min(h0, h1);
    float highHeight = max(h0, h1);
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
    float lowHeight = min(h0, h1);
    float highHeight = max(h0, h1);
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

  Float furthestIntersectionWithRay(pt o, vec d) {
    assert cone != null;
    float[] ts = cone.intersectLine(o, d);
    if (ts == null) return null;

    pt p0 = ts[0] > 0 ? P(o, ts[0], d) : null;
    pt p1 = ts[1] > 0 ? P(o, ts[1], d) : null;

    boolean valid0 = false;
    float lowHeight = min(h0, h1);
    float highHeight = max(h0, h1);
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
      if (valid1 && ts[1] > ts[0]) return new Float(ts[1]);
      else return new Float(ts[0]);
    } else {
      if (valid1) return new Float(ts[1]);
      else return null;
    }
  }

  void show(int numSamples, boolean showStroke) {
    generateSamples(numSamples, null);
    int n = samples.length / 2;
    if (showStroke) stroke(0);
    beginShape(QUAD_STRIP);
    for (int i = 0; i < n; ++i) {
      vertex(samples[n + i]);
      vertex(samples[i]);
    }
    vertex(samples[n]);
    vertex(samples[0]);
    endShape();
    if (showStroke) noStroke();
  }

  void showDebugInfo() {
    fill(snow, 255);
    showBall(c0, 5);
    fill(black, 255);
    showBall(c1, 5);
  }
}

TruncatedCone truncatedConeOfTwoBalls(Ball b0, Ball b1) {
  vec v = V(b0.c, b1.c);
  float d = v.norm();
  v.div(d);  // normalize
  float cos = (b0.r - b1.r) / d;
  float sin = sin(acosClamp(cos));
  pt c0 = P(b0.c, b0.r * cos, v);  // the center corresponding to ball b0
  pt c1 = P(b1.c, b1.r * cos, v);  // the center corresponding to ball b1
  float r0 = b0.r * sin;
  float r1 = b1.r * sin;
  return new TruncatedCone(c0, r0, c1, r1);
}

class ConvexGap {
  ArrayList<pt> points0;  // points on the inner loop
  ArrayList<pt> points1;  // points on the outer loop, orientation same as inner loop

  boolean normalized = false;  // true if this object is a normalized version of another one

  ConvexGap() {}

  ConvexGap(ArrayList<pt> points0, ArrayList<pt> points1) {
    this.points0 = points0;
    this.points1 = points1;
  }

  private ArrayList<pt> positionList() {
    ArrayList<pt> posList = new ArrayList<pt>();
    posList.addAll(points0);
    posList.addAll(points1);
    return posList;
  }

  void removeDuplicatePoints(int idx) {
    ArrayList<pt> points = idx == 0 ? points0 : points1;
    if (points == null) return;
    int i = 1;
    while (i < points.size()) {
      if (samePt(points.get(i), points.get(i-1))) {
        points.remove(i);
      } else i++;
    }
    if (samePt(points.get(i-1), points.get(0))) {
      points.remove(i-1);
    }
  }

  /*
   * Make the convex gap center at (0, 0, 0) and have appropriate scale.
   * c is the center of this object, s is the "target" scale.
   * If rawS is true, then the object is scaled by s.
   * If rawS is false, then the object is scaled to s.
   */
  ConvexGap normalize(pt c, float s, boolean rawS) {
    if (c == null) c = this.center();
    vec t = V(c).rev();
    ArrayList<pt> ps0 = new ArrayList<pt>();
    ArrayList<pt> ps1 = new ArrayList<pt>();
    for (pt p : points0) ps0.add(P(p, t));
    for (pt p : points1) ps1.add(P(p, t));

    if (rawS == false) {
      float r = 0.0;
      for (pt p : ps0) r += sqrt(dot(p, p));
      for (pt p : ps1) r += sqrt(dot(p, p));
      r /= (ps0.size() + ps1.size());
      s /= r;
    }

    for (pt p : ps0) p.mul(s);
    for (pt p : ps1) p.mul(s);
    return new ConvexGap(ps0, ps1);
  }

  void setNormalized(boolean n) {
    normalized = n;
  }

  private ArrayList<Triangle> gapHull() {
    int nv0 = points0.size();
    int nv1 = points1.size();
    // println(nv0, nv1);
    assert nv0 >= 3 && nv1 >= 3;

    ArrayList<Triangle> triangles = new ArrayList<Triangle>();

    /* Find the first triangle. */
    pt pa = points0.get(0);
    pt pb = points0.get(1);
    int j = 0;
    int jLeft = nv1 - 1, jRight;
    for (; j < nv1; ++j) {
      jRight = (j + 1) % nv1;
      vec vLeft = U(pa, points1.get(jLeft));
      vec vRight = U(pa, points1.get(jRight));
      vec n = U(N(pa, pb, points1.get(j)));
      {
        // println("n = ", n, "vLeft = ", vLeft, "vRight = ", vRight, "dot(n, vLeft) = ", dot(n, vLeft), "dot(n, vRight) = ", dot(n, vRight));
      }
      if (dot(n, vLeft) <= 0 && dot(n, vRight) <= 0) break;
      jLeft = j;
    }
    if (j == nv1) {
      if (normalized) return null;

      // println("Cannot find first triangle! Let's normalize the gap!");
      pt c = this.center();
      for (float k = 10.0; k < 1001.0; k *= 10.0) {
        ConvexGap gap = this.normalize(c, k, true);
        gap.setNormalized(true);
        ArrayList<Triangle> tris = gap.gapHull();
        if (tris != null) return tris;
      }

      println("Still can't find first triangle after normalization!");
      return null;
    }
    assert j < nv1;
    triangles.add(new Triangle(0, 1, j + nv0));

    {
      // fill(red);
      // showBall(points0.get(0), 3);
      // fill(green);
      // showBall(points0.get(1), 3);
      // fill(blue);
      // showBall(points1.get(j), 3);
    }

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
        // if (i == 1) {
          // fill(navy);
          // arrow(pa, vLeft, 5);
          // fill(chocolate);
          // arrow(pa, vRight, 5);
          // for (pt p : points1) println("p =", p);
          // println("pa =", pa, "points1.get(jNext) =", points1.get(jNext));
          // println("j =", j, "jNext =", jNext);
          // println("n =", n, "vLeft =", vLeft, "vRight =", vRight);
          // println("dot(n, vLeft) =", dot(n, vLeft), "dot(n, vRight) =", dot(n, vRight));
        // }
        if (dot(n, vLeft) <= 0 && dot(n, vRight) <= 0) valid = true;
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

  private int toGlobalPID(int i, ArrayList<Integer> pIDs0, ArrayList<Integer> pIDs1) {
    if (i < pIDs0.size()) return pIDs0.get(i);
    else return pIDs1.get(i - pIDs0.size());
  }

  /* Map each local vertex ID to its global vertex ID. */
  ArrayList<Triangle> gapHullGlobal(ArrayList<Integer> pIDs0, ArrayList<Integer> pIDs1) {
    ArrayList<Triangle> tris = gapHull();
    if (tris == null) return null;

    {  // debug
      // ArrayList<pt> ps = new ArrayList<pt>(points0);
      // ps.addAll(points1);
      // boolean pass = passConvexityTest(tris, ps);
      // if (!pass) {
      //   println("Gap hull fails convexity test!");
      //   return null;
      // }
    }

    for (Triangle t : tris) {
      t.set(toGlobalPID(t.a, pIDs0, pIDs1), toGlobalPID(t.b, pIDs0, pIDs1), toGlobalPID(t.c, pIDs0, pIDs1));
    }
    return tris;
  }

  TriangleMesh toTriMesh() {
    ArrayList<Triangle> triList = gapHull();
    if (triList == null) return null;
    ArrayList<pt> posList = positionList();
    return new TriangleMesh(posList, triList);
  }

  void show() {
    showOrientedLoop(points0);
    showOrientedLoop(points1);
  }

  pt center() {
    pt c = new pt();
    for (pt p : points0) c.add(p);
    for (pt p : points1) c.add(p);
    c.div(points0.size() + points1.size());
    return c;
  }

  float averageDistanceTo(pt c) {
    float r = 0;
    for (pt p : points0) r += d(p, c);
    for (pt p : points1) r += d(p, c);
    return r / (points0.size() + points1.size());
  }

  void translate(vec v) {
    for (pt p : points0) p.set(P(p, v));
    for (pt p : points1) p.set(P(p, v));
  }

  void scale(float s) {
    for (pt p : points0) p.set(P(s, p));
    for (pt p : points1) p.set(P(s, p));
  }

  void save(String file) {
    if (points0 == null || points1 == null) return;
    println("saving convex gap:", file);
    String[] lines = new String[points0.size() + points1.size() + 2];
    int i = 0;
    lines[i++] = str(points0.size());
    for (pt p : points0) lines[i++] = str(p.x) + "," + str(p.y) + "," + str(p.z);
    lines[i++] = str(points1.size());
    for (pt p : points1) lines[i++] = str(p.x) + "," + str(p.y) + "," + str(p.z);
    saveStrings(file, lines);
  }

  void load(String file) {
    println("loading convex gap:", file);
    points0 = new ArrayList<pt>();
    points1 = new ArrayList<pt>();
    String[] lines = loadStrings(file);
    int i = 0;
    int nv0 = int(lines[i++]);
    for (int j = 0; j < nv0; ++j) {
      float[] tmp = float(split(lines[i++], ","));
      points0.add(new pt(tmp[0], tmp[1], tmp[2]));
    }
    int nv1 = int(lines[i++]);
    for (int j = 0; j < nv1; ++j) {
      float[] tmp = float(split(lines[i++], ","));
      points1.add(new pt(tmp[0], tmp[1], tmp[2]));
    }
  }
}


/*
 * Hub class.
 *
 * A hub is the union of a ball and a set of tangential cones defined by it and
 * each of its neighboring balls.
 */
class Hub {
  Ball ball = null;  // inner ball
  Ball[] neighbors = null;  // outer balls
  int nNeighbors = 0;  // number of outer balls

  boolean valid = false;  // whether the hub is valid or not

  /* Intermediate variables in hub triangulation. */
  TruncatedCone[] tCones = null;  // each truncated cone connects the inner ball and an outer ball
  float maxIntersectDist = -1.0;  // the radius of the bounding sphere
  Circle[] circles = null;  // intersection between tCones and the bounding sphere
  TruncatedCone[] liftedCones = null;  // truncated tCones
  ConvexGap[] gaps = null;  // a gap is created between a lifted cone and its corresponding hole of the convex hull

  /* Parameters for hub triangulation. */
  boolean isHalved = false;  // whether each tCone is halved or not
  float inflationRatio = 1.03;  // the ratio of the radius of the bounding sphere over the max intersection distance, default: 1.05
  int numEdgesRegularPolygon = -1;  // the number of edges of a beam cross section
  vec[] xDirections = null;  // each beam has an I vector to generate samples
  float gapDistance = 0.0;  // how much each tCone should be lifted/chopped

  RingSet ringSet;

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
      float sin = sin(acosClamp(cos));
      pt c0 = P(ball.c, ball.r * cos, v);  // the center corresponding to the inner ball
      pt c1 = P(neighbors[i].c, neighbors[i].r * cos, v);  // the center corrsponding to the i-th outer ball
      float r0 = ball.r * sin;
      float r1 = neighbors[i].r * sin;

      tCones[i] = new TruncatedCone(c0, r0, c1, r1);
      tCones[i].initCone();
    }
  }

  float intersectionDistance(int i, int j) {
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
    if (px == null) {  // parallel lines
      // println("parallel lines");
      // println("difference =", V(pba11, pca01));
      // assert isZeroVec(V(pba11, pca01));
      px = P(pba11, pca01);
    }

    {
      // println("intersection point =", px);
      // fill(red);
      // showBall(px, 5);
      // if (i == 0 && j == 1) {
      //   println("pba10 =", pba10, "pba11 =", pba11);
      //   println("pca00 =", pca00, "pca01 =", pca01);
      //   strokeWeight(5);
      //   stroke(green);
      //   showSegment(pba10, pba11);
      //   stroke(blue);
      //   showSegment(pca00, pca01);
      //   strokeWeight(1);
      //   noStroke();
      // }
    }

    return d(pa, px);
  }

  private float maximumIntersectionDistance() {
    float t = -1.0;
    for (int i = 0; i < nNeighbors - 1; ++i) {
      for (int j = i + 1; j < nNeighbors; ++j) {
        float d = intersectionDistance(i, j);
        // println("i =", i, "j =", j, "d =", d);
        t = max(t, d);
      }
    }
    t = max(t, ball.r);
    // println("t =", t);
    return t;
  }

  private void init() {
    generateTruncatedCones();  // a truncated cone is a beam
    maxIntersectDist = maximumIntersectionDistance();
    maxIntersectDist *= inflationRatio;
  }

  void setInflationRatio(float ratio) {
    assert ratio > 1.0;
    inflationRatio = ratio;
  }

  void setNumEdgesRegularPolygon(int numEdges) {
    assert numEdges >= 3;
    numEdgesRegularPolygon = numEdges;
  }

  void setXDirestions(vec[] xAxes) {
    assert xAxes.length == nNeighbors;
    xDirections = xAxes;
  }

  void setGapDistance(float dist) {
    assert dist > 0.0;
    gapDistance = dist;
  }

  void setIsHalved(boolean halve) {
    isHalved = halve;
  }

  void halveBeams() {
    for (int i = 0; i < nNeighbors; ++i) {
      tCones[i].halve();
    }
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

  void generateIntersectionCircles() {
    if (maxIntersectDist <= 0.0) return;
    circles = new Circle[nNeighbors];
    for (int i = 0; i < nNeighbors; ++i) {
      circles[i] = intersectionCircleWithBoundingSphere(i);
    }

    for (Circle cir : circles) {
      if (cir == null) {
        println("circle is null!");
        save("data/hub_unnamed");
        break;
      }
    }
    return;
  }

  private TruncatedCone liftCone(int i, pt c, float r, float gapWidth) {
    pt c0 = P(c, gapWidth, tCones[i].normal);
    float r0 = r;
    if (tCones[i].r1 < tCones[i].r0) {
      r0 -= gapWidth * tan(tCones[i].cone.alpha);
    } else {
      r0 += gapWidth * tan(tCones[i].cone.alpha);
    }

    if (dot(V(c0, tCones[i].c1), tCones[i].normal) < 0) {
      // println("The inner base is lifted too much." +
      //         "The axis direction will be reversed." +
      //         "To fix this, the outer base will also be lifted.");
      pt c1 = P(tCones[i].c1, 2 * gapWidth, tCones[i].normal);
      float r1 = tCones[i].r1;
      if (tCones[i].r1 < tCones[i].r0) {
        r1 -= 2 * gapWidth * tan(tCones[i].cone.alpha);
      } else {
        r1 += 2 * gapWidth * tan(tCones[i].cone.alpha);
      }
      return new TruncatedCone(c0, r0, c1, r1, tCones[i].normal, tCones[i].cone);
    }
    return new TruncatedCone(c0, r0, tCones[i].c1, tCones[i].r1, tCones[i].normal, tCones[i].cone);
  }

  void liftCones(float gapWidth) {
    liftedCones = new TruncatedCone[nNeighbors];
    for (int i = 0; i < nNeighbors; ++i) {
      if (circles[i] != null) {
        liftedCones[i] = liftCone(i, circles[i].c, circles[i].r, gapWidth);
      }
    }
  }

  RingSet circlesToRingSet(int m) {
    return new RingSet(ball.c, maxIntersectDist, circles, nNeighbors, m);
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
    float k = 0.185;  // shouldn't be too large, good values: 0.2,
    for (int i = 0; i < nNeighbors; ++i) {
      float t = exp(-k * roundConeDist(p, ball.c, ball.r, neighbors[i].c, neighbors[i].r));
      assert t > 0.0;
      d += t;
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

  float furthestIntersectionWithRay(pt o, vec d) {
    float t = -1;
    for (int i = 0; i < nNeighbors; ++i) {
      Float tmp = tCones[i].furthestIntersectionWithRay(o, d);
      if (tmp != null && tmp > t) t = tmp;
    }

    {
      float[] ts = ball.intersectLine(o, d);
      if (ts != null) {
        if (ts[0] > 0 && ts[0] > t) t = ts[0];
        if (ts[1] > 0 && ts[1] > t) t = ts[1];
      }
    }

    return t;
  }

  /*
   * Generate vertices on lifted cones (i.e. beams). m is the number of edges of
   * a beam cross section. xAxis is the I direction of a beam.
   */
  void generateBeamSamples(int m, vec[] xAxes) {
    for (int i = 0; i < nNeighbors; ++i) {
      liftedCones[i].generateSamples(m, xAxes == null ? null : xAxes[i]);
    }
  }

  /*
   * Triangulate the i-th beam. k is the ID of the first point of this beam. The
   * triangle mesh tm will be augmented with new points and new triangles.
   * Return the number of new points.
   */
  private int triangulateBeam(int b, int k, TriangleMesh tm) {
    TruncatedCone beam = liftedCones[b];
    assert beam != null;

    pt[] samples = beam.samples;
    assert samples != null && samples.length > 0 && samples.length % 2 == 0;

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
   * Triangulate all the incident beams. Return a single mesh.
   */
  TriangleMesh generateBeamMesh() {
    TriangleMesh tm = new TriangleMesh();
    for (int i = 0, k = 0; i < nNeighbors; ++i) {
      k += triangulateBeam(i, k, tm);
    }
    return tm;
  }

  /*
   * Initialize gaps, each formed by a hole of the convex hull and the base of
   * its corresponding beam. innerLoops[i] is the i-th hole of the convex hull.
   */
  void initGaps(ArrayList<pt>[] innerLoops) {
    assert innerLoops != null && innerLoops.length == nNeighbors;
    assert liftedCones != null && liftedCones.length == nNeighbors;
    gaps = new ConvexGap[nNeighbors];

    for (int i = 0; i < nNeighbors; ++i) {
      ArrayList<pt> innerLoop = innerLoops[i];
      ArrayList<pt> outerLoop = new ArrayList<pt>();
      pt[] samples = liftedCones[i].samples;
      assert samples != null;
      int n = samples.length / 2;
      for (int j = 0; j < n; ++j) outerLoop.add(samples[j]);
      gaps[i] = new ConvexGap(innerLoop, outerLoop);
    }
  }

  /*
   * Triangulate all the gaps. Return a single mesh.
   */
  TriangleMesh generateGapMesh() {
    if (gaps == null || gaps.length == 0) return null;

    TriangleMesh gapMesh = new TriangleMesh();
    for (int i = 0; i < nNeighbors; ++i) {
      TriangleMesh tm = gaps[i].toTriMesh();
      gapMesh.augmentWithShift(tm.positions, tm.triangles);
    }
    return gapMesh;
  }

  BorderedTriangleMesh generateConvexHullMesh(int numEdges) {
    generateIntersectionCircles();
    RingSet rs = circlesToRingSet(numEdges);
    ringSet = rs;
    if (rs.nRings >= 3) rs.generateExactCHIncremental(null);
    TriangleMesh tm = rs.generateConvexTriMesh();
    ArrayList<Integer>[] borders = rs.getSwingLists();  // each border is a loop of vertex IDs
    return new BorderedTriangleMesh(tm, borders);
  }

  BorderedTriQuadMesh generateConvexHullTQMesh() {
    assert ringSet != null;  // please also make sure that ring set has already generated the convex hull
    return ringSet.generateConvexTriQuadMesh();
  }

  /*
   * Generate a triangle mesh that approximates the hub. This mesh is made of
   * convex hull mesh, gap mesh and beam mesh.
   *
   * Parameters:
   * numEdges: an integer, the number of edges on a beam cross section
   * gapWidth: a float, the width of each gap
   * xAxes: an array of vectors with size nNeighbors, each x-axis is the I direction of a beam
   * halve: whether or not to halve each beam
   */
  TriangleMesh generateTriMesh(int numEdges, float gapWidth, vec[] xAxes, boolean halve) {
    if (halve) halveBeams();

    generateIntersectionCircles();
    liftCones(gapWidth);
    RingSet rs = circlesToRingSet(numEdges);

    rs.generateExactCHIncremental(null);
    TriangleMesh convexHullMesh = rs.generateConvexTriMesh();

    generateBeamSamples(numEdges, xAxes);  // make sure to call this function before generating gap mesh and beam mesh
    initGaps(rs.borders);
    TriangleMesh gapMesh = generateGapMesh();
    TriangleMesh beamMesh = generateBeamMesh();

    /* Merge the 3 meshes into 1 mesh. */
    TriangleMesh triMesh = new TriangleMesh(convexHullMesh);  // copy
    triMesh.augmentWithShift(beamMesh);

    HashMap<pt, Integer> vIDs = new HashMap<pt, Integer>();
    {  // store the global vertex ID for each vertex of the gap mesh
      for (int i = 0; i < gapMesh.nv; ++i) {
        pt p = gapMesh.positions.get(i);
        ArrayList<pt> ps = triMesh.positions;
        for (int j = 0; j < triMesh.nv; ++j) {
          if (p == ps.get(j)) {
            vIDs.put(p, j);
            break;
          }
        }
      }
      assert vIDs != null && vIDs.size() > 0;
    }
    {  // merge the current triangle mesh with the gap mesh
      ArrayList<pt> ps = gapMesh.positions;
      for (int i = 0; i < gapMesh.nt; ++i) {
        Triangle t = gapMesh.triangles.get(i);
        int a = vIDs.get(ps.get(t.a));
        int b = vIDs.get(ps.get(t.b));
        int c = vIDs.get(ps.get(t.c));
        triMesh.addTriangle(new Triangle(a, b, c));
      }
    }

    return triMesh;
  }

  TriangleMesh generateTriMesh() {
    assert numEdgesRegularPolygon >= 3 && gapDistance > 0.0;
    return generateTriMesh(numEdgesRegularPolygon, gapDistance, xDirections, isHalved);
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
      showBall(ball.c, maxIntersectDist);
    }
  }

  void showLiftedCones(color c, int a) {
    fill(c, a);
    int n = gNumPointsPerRing;
    for (int i = 0; i < nNeighbors; ++i) {
      if (liftedCones[i] != null) {
        liftedCones[i].show(n, true);
        // liftedCones[i].showDebugInfo();
      }
    }
  }

  void show(color c, int a) {
    fill(c, a);
    noStroke();
    ball.show();
    for (int i = 0; i < nNeighbors; ++i) {
      tCones[i].show(40, false);
      // neighbors[i].show();
    }
    return;
  }

  void save(String file) {
    println("saving hub:", file);
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
  }

  void load(String file) {
    println("loading hub:", file);
    String[] lines = loadStrings(file);
    int i = 0;
    nNeighbors = int(lines[i++]) - 1;
    float[] tmp = float(split(lines[i++], ","));
    ball = new Ball(tmp[0], tmp[1], tmp[2], tmp[3]);
    neighbors = new Ball[nNeighbors];
    for (int j = 0; j < nNeighbors; ++j) {
      tmp = float(split(lines[i++], ","));
      neighbors[j] = new Ball(tmp[0], tmp[1], tmp[2], tmp[3]);
    }

    valid = isValid();
    if (valid) this.init();
    else println("The hub is not valid! At lease two balls intersect.");
  }

  void saveAugFile(String file) {
    if (xDirections == null) return;
    println("saving augmented hub:", file);
    String[] lines = new String[3 * nNeighbors + 3];
    int i = 0;
    lines[i++] = str(numEdgesRegularPolygon);
    lines[i++] = str(ball.r);
    lines[i++] = str(nNeighbors);
    for (int j = 0; j < nNeighbors; ++j) {
      lines[i++] = str(neighbors[j].c.x) + "," + str(neighbors[j].c.y) + "," +
                   str(neighbors[j].c.z);
      lines[i++] = str(neighbors[j].r);
      lines[i++] = str(xDirections[j].x) + "," + str(xDirections[j].y) + "," +
                   str(xDirections[j].z);
    }
    saveStrings(file, lines);
  }

  /*
   * Construct a hub given an augmented file. This augmented file contains
   * information that defines the hub, and information about triangulation (e.g.
   * number of edges of a beam base).
   */
  void loadAugFile(String file) {
    println("loading augmented hub:", file);
    String[] lines = loadStrings(file);
    int i = 0;
    // number of edges on each regular polygon
    numEdgesRegularPolygon = int(lines[i++]);
    // the center ball centered at (0, 0, 0)
    ball = new Ball(new pt(0.0, 0.0, 0.0), float(lines[i++]));
    // the number of beam balls
    nNeighbors = int(lines[i++]);
    neighbors = new Ball[nNeighbors];
    xDirections = new vec[nNeighbors];
    for (int j = 0; j < nNeighbors; ++j) {
      float[] tmp = float(split(lines[i++], ","));
      float r = float(lines[i++]);
      neighbors[j] = new Ball(tmp[0], tmp[1], tmp[2], r);
      tmp = float(split(lines[i++], ","));
      xDirections[j] = new vec(tmp[0], tmp[1], tmp[2]);
    }

    valid = isValid();
    if (valid) this.init();
    else println("The hub is not valid! At lease two balls intersect.");
  }
}



/*
 * Instructions on converting a hub into a triangle mesh.
 *
 * gHub = new Hub();
 * gHub.loadAugFile(file);  // the data file of a hub
 * float gapWidth = 3.0;  // can put this into the data file later
 * TriangleMesh tm = generateTriMesh(gHub.numEdgesRegularPolygon, gapWidth, gHub.xDirections, true);  // true means each beam will be halved
 */
