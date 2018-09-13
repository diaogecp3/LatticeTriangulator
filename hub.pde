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
  float dba = vba.norm();
  float dca = vca.norm();
  vba = V(1.0/dba, vba);  // normalize vba
  vca = V(1.0/dca, vca);  // normalize vca
  vec normal = N(vba, vca);  // cross(BA, BC) = cross(AB, AC)
  vec axisJ = N(normal, vba);
  vec vbaRotated = R(vba, HALF_PI, vba, axisJ);
  vec vcaRotated = R(vca, HALF_PI, vca, axisJ);
  float cosb = (rb - ra) / dba;
  float cosc = (rc - ra) / dca;
  float sinb = sqrt(1 - cosb * cosb);  // sinb is always positive
  float sinc = sqrt(1 - cosc * cosc);  // sinc is always positive
  vec v0 = V(cosb, vba, -sinb, vbaRotated);
  vec v1 = V(cosc, vca, sinc, vcaRotated);
  pt pba10 = P(pb, rb, v0);
  pt pba11 = P(pa, ra, v0);
  pt pca00 = P(pc, rc, v1);
  pt pca01 = P(pa, ra, v1);
  pt px = intersectionTwoLines(pba10, pba11, pca00, pca01);
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
  vab = V(1.0/d, vab);  // normalize vab
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
  Ball ball;
  Ball[] neighbors;
  int nNeighbors;
  Hub(Ball ball, Ball[] neighbors, int nNeighbors) {
    this.ball = ball;
    this.neighbors = neighbors;
    this.nNeighbors = nNeighbors;
  }

  float maximumIntersectionDistance() {
    float t = -1.0;
    for (int i = 0; i < nNeighbors; ++i) {
      for (int j = i + 1; j < nNeighbors; ++j) {
        float d = intersectionDistance(ball, neighbors[i], neighbors[j]);
        t = max(t, d);
      }
    }
    return t;
  }

  void showBoundingBall(color c, int a) {
    float t = maximumIntersectionDistance();
    fill(c, a);
    noStroke();
    show(ball.c, t);
  }

  void showHub(color c, int a) {
    fill(c, a);
    noStroke();
    ball.showBall();
    for (int i = 0; i < nNeighbors; ++i) {
      showTangentialCone(ball, neighbors[i]);
      neighbors[i].showBall();
    }
    return;
  }
}