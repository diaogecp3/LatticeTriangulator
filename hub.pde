/******************************************************************************
 * Hub (a ball and tangential cones) processing.
 ******************************************************************************/

boolean showHub = true;
boolean showBoundingSphere = false;  // the bounding sphere of a hub
boolean showIntersectionCircles = false;  // the intersection circles between the bounding sphere and cones


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
   * Compute the intersection between a line and a two-cone. Return the two
   * parameters, if any.
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
  pt c1;  // center of the upper base
  float r0;  // radius of the lower base
  float r1;  // radius of the upper base
  vec normal;  // unit vector from c0 to c1

  Cone cone;  // the original cone
  float lowHeight;
  float highHeight;

  TruncatedCone(pt c0, pt c1, float r0, float r1) {
    this.c0 = c0;
    this.c1 = c1;
    this.r0 = r0;
    this.r1 = r1;
    normal = U(c0, c1);
    initCone();
  }

  void initCone() {
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

  void show() {
    int n = 30;
    pt[] points0 = new pt[n];
    pt[] points1 = new pt[n];
    vec vi = constructNormal(normal);
    vec vj = N(normal, vi);
    float da = TWO_PI / n;
    float a = 0;
    for (int i = 0; i < n; ++i) {
      float cos = cos(a);
      float sin = sin(a);
      points0[i] = P(c0, r0 * cos, vi, r0 * sin, vj);
      points1[i] = P(c1, r1 * cos, vi, r1 * sin, vj);
      a += da;
    }
    beginShape(QUAD_STRIP);
    for (int i = 0; i < n; ++i) {
      vertex(points1[i]);
      vertex(points0[i]);
    }
    vertex(points1[0]);
    vertex(points0[0]);
    endShape();
  }
}


/*
 * Hub class.
 *
 * A hub is the union of a ball and a set of tangential cones
 * defined by it and each of its neighboring balls.
 */
class Hub {
  Ball ball = null;
  Ball[] neighbors = null;
  int nNeighbors = 0;

  boolean valid = false;  // whether the hub is valid or not
  float maxIntersectDist = -1.0;  // the radius of the bounding sphere

  TruncatedCone[] tCones = null;
  Circle[] circles = null;

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

      tCones[i] = new TruncatedCone(c0, c1, r0, r1);
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
    maxIntersectDist += 3;  // slightly bigger such that circles are disjoint
  }

  private Circle intersectionCircleWithBoundingSphere(int i) {
    float r0 = ball.r * ball.r / tCones[i].r0;
    float r1 = neighbors[i].r * neighbors[i].r / tCones[i].r1;
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
    for (int i = 0; i < nNeighbors; ++i) {
      circles[i] = intersectionCircleWithBoundingSphere(i);
    }
    return;
  }

  RingSet circlesToRingSet() {
    return new RingSet(ball.c, maxIntersectDist, circles, nNeighbors);
  }

  float closestIntersectionWithLine(pt o, vec d) {
    float t = 1e20;
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

    {
      // fill(tomato, 255);
      // show(P(o, t, d), 4);
    }

    return t;
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
