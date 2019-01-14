/******************************************************************************
 * Hub (a ball and tangential cones) processing.
 ******************************************************************************/


/*
 * A tangential cone touches ball ba and bb, and another tangential cone touches
 * ball ba and bc. Find the minimum radius such that, if ball ba has a radius
 * bigger than that, the two cones with be disjoint outside of ball ba.
 */
float intersectionDistance(Ball ba, Ball bb, Ball bc) {
  pt pa = ba.c;
  pt pb = bb.c;
  pt pc = bc.c;
  float ra = ba.r;
  float rb = bb.r;
  float rc = bc.r;
  vec vba = V(pb, pa);  // used as axis-I
  vec vca = V(pc, pa);
  float dba = vba.norm();  // length of vba
  float dca = vca.norm();
  vba.div(dba);  // normalize vba
  vca.div(dca);  // normalize vca
  vec normal = N(vba, vca);  // cross(BA, BC) = cross(AB, AC)
  vec axisJ = N(normal, vba);
  vec vbaRotated = axisJ;
  vec vcaRotated = R(vca, HALF_PI, vba, axisJ);
  float cosb = (rb - ra) / dba;
  float cosc = (rc - ra) / dca;
  float sinb = sqrt(1 - cosb * cosb);  // sinb is always positive
  float sinc = sqrt(1 - cosc * cosc);  // sinc is always positive


  // TODO: need to revisit
  vec v0 = V(cosb, vba, -sinb, vbaRotated);
  vec v1 = V(cosc, vca, sinc, vcaRotated);

  // vec v0 = V(sinb, vba, -cosb, vbaRotated);
  // vec v1 = V(sinc, vca, cosc, vcaRotated);

  pt pba10 = P(pb, rb, v0);
  pt pba11 = P(pa, ra, v0);
  pt pca00 = P(pc, rc, v1);
  pt pca01 = P(pa, ra, v1);
  pt px = intersectionTwoLines(pba10, pba11, pca00, pca01);

  {  // show the intersection point
    fill(green, 100);
    show(px, 4);
  }

  return d(pa, px);
}

/*
 * Show the tangential cone touching ball ba and bb.
 */
void showTangentialCone(Ball ba, Ball bb) {
  int n = 36;
  pt[] points0 = new pt[n];
  pt[] points1 = new pt[n];
  pt ca = ba.c, cb = bb.c;
  float ra = ba.r, rb = bb.r;
  vec vab = V(ca, cb);
  float d = vab.norm();
  vab.div(d);  // normalize vab
  vec axisI = constructNormal(vab);
  vec axisJ = N(vab, axisI);
  /* Compute first contact point on ba (and on bb). */
  float cosa = (ra - rb) / d;
  float sina = sqrt(1 - cosa * cosa);
  vec vabRotated = R(vab, HALF_PI, axisJ, vab);
  vec v = V(cosa, vab, sina, vabRotated);
  points0[0] = P(ca, ra, v);
  points1[0] = P(cb, rb, v);
  pt c0 = P(ca, ra * cosa, vab);
  pt c1 = P(cb, rb * cosa, vab);
  float da = TWO_PI / n;
  float a = da;
  /* Compute the rest contacts. */
  for (int i = 1; i < n; ++i) {
    points0[i] = R(points0[0], a, axisI, axisJ, c0);
    points1[i] = R(points1[0], a, axisI, axisJ, c1);
    a += da;
  }
  beginShape(QUAD_STRIP);
  for (int i = 0; i < n; ++i) {
    v(points1[i]);
    v(points0[i]);
  }
  v(points1[0]);
  v(points0[0]);
  endShape();
  return;
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
  private float maxIntersectDist = -1.0;  // the radius of the bounding sphere

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

  private float maximumIntersectionDistance() {
    float t = -1.0;
    for (int i = 0; i < nNeighbors; ++i) {
      for (int j = i + 1; j < nNeighbors; ++j) {
        float d = intersectionDistance(ball, neighbors[i], neighbors[j]);
        t = max(t, d);
      }
    }
    return t;
  }

  private void init() {
    maxIntersectDist = maximumIntersectionDistance();
    maxIntersectDist += 5;  // slightly bigger
  }

  private Circle intersectionCircleWithBoundingSphere(int i) {
    float r0 = ball.r;
    float r1 = neighbors[i].r;
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
      // println("sols[0] = ", sols[0], "sols[1] = ", sols[1]);
      t = max(sols[0], sols[1]);
    }

    if (t >= 0 && t <= 1) {
      float r = (1 - t) * r0 + t * r1;
      vec n = U(v);
      pt center = P(ball.c, t * sqrt(dd), n);
      return new Circle(center, n, r);
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
      showTangentialCone(ball, neighbors[i]);
      //neighbors[i].showBall();
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
