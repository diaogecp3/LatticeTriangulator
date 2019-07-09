/******************************************************************************
 * Gap (the gap between two convex polygonal loops).
 ******************************************************************************/

/*
 * 0: wants two negative dot products
 * 1: wants smaller pivot angle (slower but more robust?)
 */
int gMethodConvexGap = 1;

/*
 * Convex gap class.
 *
 * The convex hull used to fill the gap between two convex polygonal loops.
 */
class ConvexGap {
  ArrayList<pt> points0;  // points on the inner loop
  ArrayList<pt> points1;  // points on the outer loop, orientation same as inner loop

  boolean normalized = false;  // true if it has been normalized

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

  /*
   * Remove duplicate points in the idx-th point array.
   */
  void removeDuplicatePoints(int idx) {
    removeDuplicates(idx == 0 ? points0 : points1);
  }

  /*
   * Make the gap center at (0, 0, 0) and have appropriate scale.
   *
   * Parameters:
   * c: the center of the current gap
   * s: the target scale or scale factor
   * rawS: true: the gap is scaled by s; false: the gap is scaled to s
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

  private void setNormalized(boolean n) {
    normalized = n;
  }

  /*
   * Compute the convex hull of two point-pairs (i.e., four points in total).
   */
  private ArrayList<Triangle> gapHullTwoTwo() {
    pt a = points0.get(0);
    pt b = points0.get(1);
    pt c = points1.get(0);
    pt d = points1.get(1);
    int ic = 2, id = 3;
    vec u = U(N(a, c, b));
    vec v = U(c, d);
    if (dot(u, v) > 0) {  // swap indices of c and d
      ic = 3;
      id = 2;
    }
    ArrayList<Triangle> tris = new ArrayList<Triangle>();
    tris.add(new Triangle(1, 0, ic));
    tris.add(new Triangle(ic, id, 1));
    tris.add(new Triangle(0, 1, id));
    tris.add(new Triangle(id, ic, 0));
    return tris;
  }

  /*
   * Compute the convex hull of the two convex polygonal loops.
   */
  private ArrayList<Triangle> gapHull() {
    int nv0 = points0.size();
    int nv1 = points1.size();
    assert nv0 >= 2 && nv1 >= 2;
    if (nv0 == 2 && nv1 == 2) return gapHullTwoTwo();

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
        // println("dot(n, vLeft) =", dot(n, vLeft), "dot(n, vRight) =", dot(n, vRight));
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

    // println("j =" + j + ", nv1 =" + nv1);
    assert j < nv1;
    triangles.add(new Triangle(0, 1, j + nv0));

    /* Traverse the loops. */
    int stop0 = 0;
    int stop1 = j;
    int i = 1;
    int iNext = (i + 1) % nv0;
    int jNext = (j + 1) % nv1;
    int jCount = 0;

    vec nPrev = null;
    if (gMethodConvexGap != 0) nPrev = U(N(points0.get(0), points0.get(1), points1.get(j)));  // the normal of previous face

    while (i != stop0 && jCount < nv1) {
      pa = points1.get(j);
      pb = points0.get(i);
      pt pc = points0.get(iNext);  // the third vertex of the first candidate
      pt pd = points1.get(jNext);  // the third vertex of the second candidate
      vec n = U(N(pa, pb, pc));  // the normal of the first candidate
      vec m = U(N(pa, pb, pd));  // the normal of the second candidate

      /* Pick one candidate. */
      boolean valid = false;  // true if pick the first candidate (j, i, i + 1)
      // println("first choice, j =", j, "i =", i, "iNext =", iNext);
      // println("second choice, j =", j, "i =", i, "jNext =", jNext);

      if (gMethodConvexGap == 0) {  // method 1 wants two negative dot products
        vec v0 = U(pa, points1.get((j + nv1 - 1) % nv1));
        vec v1 = U(pa, points1.get(jNext));
        float dFirstChoice = dot(n, v0) + dot(n, v1);
        // println("dot(n, v0) =", dot(n, v0), "dot(n, v1) =", dot(n, v1));

        vec u0 = U(pb, points0.get((i + nv0 - 1) % nv0));
        vec u1 = U(pb, points0.get(iNext));
        float dSecondChoice = dot(m, u0) + dot(m, u1);
        // println("dot(m, u0) =", dot(m, u0), "dot(m, u1) =", dot(m, u1));
        valid = dFirstChoice < dSecondChoice;
      } else {  // method 2 wants smaller pivot angle (slower but more robust?)
        vec ba = V(pb, pa);  // axis
        float angle0 = acosClamp(dot(nPrev, n));
        float angle1 = acosClamp(dot(nPrev, m));
        // println("before adjust, angle0 =", angle0, ", angle1 =", angle1);

        /* Only adjust angles when they are large enough. */
        if (angle0 > gEpsilonLarge && dot(N(nPrev, n), ba) < 0) angle0 = TWO_PI - angle0;
        if (angle1 > gEpsilonLarge && dot(N(nPrev, m), ba) < 0) angle1 = TWO_PI - angle1;

        // println("after adjust, angle0 =", angle0, ", angle1 =", angle1);
        valid = angle0 < angle1;
        nPrev = valid ? n : m;
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

    for (Triangle t : tris) {
      t.set(toGlobalPID(t.a, pIDs0, pIDs1), toGlobalPID(t.b, pIDs0, pIDs1), toGlobalPID(t.c, pIDs0, pIDs1));
    }
    return tris;
  }

  /*
   * Construct a convex hull of the two poly loops and use a triangle mesh to
   * represent this convex hull.
   */
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