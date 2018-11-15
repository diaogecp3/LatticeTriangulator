/******************************************************************************
 * Ring set processing.
 ******************************************************************************/

/*
 * TODO:
 * 1. Clean up code.
 */


import java.util.Queue;
import java.util.LinkedHashMap;

boolean debugST = false;
boolean debug3RT = false;
boolean debug2RT = false;
boolean fix3RT = false;
boolean show2RT = false;
boolean show3RT = false;
int numFaces3RT = 1;
int numSteps3RT = 1;
int maxIterHSGlobal = 100;
int maxIterHSLocal = 100;


/*
 * Method for three-ring triangle generation.
 * 0: Naive method. Time: O(n^3). Space: O(1)
 * 1: Heuristic search. Time: O(n). Space: O(1)
 * 2: Breadth first search. Time: roughly O(min(3^d, n^3)) where d is the length
      from initial state to optimal state. Space: O(min(3^d, n^3)). Not converge?
 * 3: Breadth first search with heuristics. Time and space should be less than 2.
 * 4: Approximated extreme plane. Time: roughly O(1), Space: roughly O(1).
 */
int method3RT = 3;

/*
 * Method for triangle mesh generation for ring set.
 * 0: Generate a convex hull. Time: O(n^2)
 * 1: Generate a mesh using referenced convex hull. Time: O(r^2 + n), where r is
      the number of rings.
 */
int methodTM = 1;


/*
 * RingSet class.
 *
 * A ring set is a set of rings (i.e. discritized cicles). We assume that these
 * rings lie on a sphere.
 */
class RingSet {
  /*
   * Debug info about supporting triangle/plane of 3 given disjoint circles on a sphere.
   */
  class DebugSTInfo {
    ArrayList<pt> circumcenters = new ArrayList<pt>();
    ArrayList<Float> circumradii = new ArrayList<Float>();
    ArrayList<vec> normals = new ArrayList<vec>();
    DebugSTInfo() {}
  }

  /*
   * Debug info about three-ring-triangle generation, used in heuristic search.
   */
  class Debug3RTInfo {
    int idr0, idr1, idr2;
    int idp0, idp1, idp2;
    int numSteps;
    int numFaces;
    Debug3RTInfo() {
      idr0 = idr1 = idr2 = -1;
      idp0 = idp1 = idp2 = -1;
      numSteps = 0;
      numFaces = 0;
    }
  }

  /*
   * Debug info about two-ring-triangle generation.
   */
  class Debug2RTInfo {
    pt pa0, pa1, pb0, pb1;
    int numGlobalStep = 1;
    int numLocalStep = 1;
  }

  /*
   * Triangle node used in BFS for three-ring-triangle generation.
   */
  class TriangleNode {
    Triangle tri;
    ArrayList<TriangleNode> children = null;
    TriangleNode(Triangle tri) {
      this.tri = tri;
    }

    boolean computeChildren(pt[][] rings, HashSet<Triangle> exploredSet) {
      assert rings.length == 3;
      int np = rings[0].length;
      pt[] ps = new pt[3];
      ps[0] = rings[0][tri.get(0)];
      ps[1] = rings[1][tri.get(1)];
      ps[2] = rings[2][tri.get(2)];
      vec normal = N(ps[0], ps[1], ps[2]);  // not necessarily unit vector
      boolean stable = true;
      pt p;
      vec v;
      children = new ArrayList<TriangleNode>();
      for (int r = 0; r < 3; ++r) {
        int[] idxs = new int[3];
        idxs[1] = tri.get((r + 1) % 3);
        idxs[2] = tri.get((r + 2) % 3);

        idxs[0] = (tri.get(r) - 1 + np) % np;
        p = rings[r][idxs[0]];
        v = V(ps[r], p);
        if (dot(normal, v) > 0) {
          if (stable == true) stable = false;
          Triangle tri = new Triangle(idxs[(3 - r) % 3], idxs[(4 - r) % 3], idxs[(5 - r) % 3]);
          if (!exploredSet.contains(tri)) {
            children.add(new TriangleNode(tri));
            exploredSet.add(tri);
          }
        }

        idxs[0] = (tri.get(r) + 1) % np;
        p = rings[r][idxs[0]];
        v = V(ps[r], p);
        if (dot(normal, v) > 0) {
          if (stable == true) stable = false;
          Triangle tri = new Triangle(idxs[(3 - r) % 3], idxs[(4 - r) % 3], idxs[(5 - r) % 3]);
          if (!exploredSet.contains(tri)) {
            children.add(new TriangleNode(tri));
            exploredSet.add(tri);
          }
        }
      }

      return stable;
    }
  }

  /*
   * A corridor has four vertices with ID a, b, c, d, respectively.
   */
  class Corridor {
    int a, b, c, d;
    Corridor(int a, int b, int c, int d) {
      this.a = a;
      this.b = b;
      this.c = c;
      this.d = d;
    }
  }

  pt c;  // the center of the sphere where the ring set lies
  float r;  // the radius of the sphere where the ring set lies
  int nRings;  // number of rings/circles
  int nPointsPerRing;  // number of points on each ring
  boolean sameRadius;  // whether all rings have the same radius
  pt[] contacts;  // the intersections between the outward normals of rings and the sphere
  float[] radii;  // the radii of rings
  vec[] xAxes;  // the x axes, one for each ring, used to generate the first point on ring
  vec[] yAxes;  // the y axes, one for each ring
  vec[] normals;  // the normals, one for each ring
  pt[][] points;  // the generated points on rings
  pt[] centers;  // the centers of rings, one for each ring

  ArrayList<Triangle> triangles = null;  // triangle mesh with ring vertices
  TriangleMesh refConvexHull = null; // convex hull generated by contacts
  ArrayList<Triangle> threeRingTriangles = null;
  ArrayList<Triangle> twoRingTriangles = null;
  ArrayList<ArrayList<Integer>> splitLists = null;

  int nSamples = 30;
  float dAngle = 0;
  ArrayList<pt> exTriPoints = null;  // points for extreme triangles, every 3 consecutive points form a trianle
  ArrayList<Integer> exTriRIDs = null;
  ArrayList<Integer>[] swingLists = null;  // each ring has a swing list to sort the corners touching it
  ArrayList<Float>[] angleLists = null;
  LinkedHashMap<Corridor, ArrayList<pt>> exEdges;  // key: (pid a, b, c, d), value: a list of points sampled from corresponding arcs

  Debug3RTInfo debug3RTInfo = new Debug3RTInfo();  // for debug
  Debug2RTInfo debug2RTInfo = new Debug2RTInfo();  // for debug
  DebugSTInfo debugSTInfo = new DebugSTInfo();

  RingSet(pt c, float r) {
    this.c = c;
    this.r = r;
    sameRadius = false;
  }

  RingSet(pt c, float r, int nc, int np) {
    this.c = c;
    this.r = r;
    this.nRings = nc;
    this.nPointsPerRing = np;
    sameRadius = false;
  }

  RingSet(pt c, float r, int nc, int np, float rMax) {
    this.c = c;
    this.r = r;
    this.nRings = nc;
    this.nPointsPerRing = np;
    sameRadius = true;
    radii = new float[1];
    radii[0] = rMax;
  }

  RingSet(pt c, float r, pt[] ps, int nv) {  // the ring set may not be valid
    this.c = c;
    this.r = r;
    this.nRings = int(nv / 2);
    this.nPointsPerRing = 3;  // default
    sameRadius = false;
    float rr = r + r;
    float r2 = r * r;
    contacts = new pt[nRings];
    centers = new pt[nRings];
    radii = new float[nRings];
    normals = new vec[nRings];
    xAxes = new vec[nRings];
    yAxes = new vec[nRings];
    for (int i = 0; i < nv; i += 2) {
      int j = i / 2;
      contacts[j] = ps[i];
      float d = d(ps[i], ps[i+1]);
      radii[j] = d * sqrt((rr-d)*(rr+d)) / rr;
      normals[j] = U(c, ps[i]);
      xAxes[j] = constructNormal(normals[j]);
      yAxes[j] = N(normals[j], xAxes[j]);
      centers[j] = P(c, sqrt(r2 - radii[j] * radii[j]), normals[j]);
    }
  }

  private void generateXYAxes() {
    xAxes = new vec[nRings];
    yAxes = new vec[nRings];
    for (int i = 0; i < nRings; ++i) {
      xAxes[i] = constructNormal(normals[i]);
      yAxes[i] = N(normals[i], xAxes[i]);
    }
  }

  void init() {
    normals = new vec[nRings];
    if (!sameRadius) {
      radii = new float[nRings];
      contacts = generateContactsAndRadii(c, r, nRings, radii, normals);
    } else {
      contacts = generateContacts(c, r, nRings, radii[0], normals);
    }
    generateXYAxes();
  }

  int getNumRings() {
    return nRings;
  }

  int getNumPointsPerRing() {
    return nPointsPerRing;
  }

  pt[] get1DPointArray() {
    if (points == null) return null;
    pt[] pointArray = new pt[nRings * nPointsPerRing];
    int k = 0;
    for (int i = 0; i < nRings; ++i) {
      for (int j = 0; j < nPointsPerRing; ++j) {
        pointArray[k++] = points[i][j];
      }
    }
    return pointArray;
  }

  ArrayList<pt> get1DPointArrayList() {
    if (points == null) return null;
    ArrayList<pt> positions = new ArrayList<pt>();
    for (int i = 0; i < nRings; ++i) {
      for (int j = 0; j < nPointsPerRing; ++j) {
        positions.add(points[i][j]);
      }
    }
    return positions;
  }

  pt[][] get2DPointArray() {
    return points;
  }

  private pt getPoint(int id) {
    return points[id / nPointsPerRing][id % nPointsPerRing];
  }

  void generatePoints(float attenuation) {
    centers = new pt[nRings];
    if (!sameRadius) {
      float[] curRadii = new float[nRings];
      for (int i = 0; i < nRings; ++i) curRadii[i] = radii[i] * attenuation;
      points = generatePointsForCircles(contacts, curRadii, c, r, xAxes, nRings,
                                        nPointsPerRing, centers);
    } else {
      float curRadius = radii[0] * attenuation;
      points = generatePointsForCircles(contacts, curRadius, c, r, xAxes, nRings,
                                        nPointsPerRing, centers);
    }
  }

  /*
   * Generate a convex hull using contacts. This convex hull is used as a
   * reference in three-ring-triangle generation.
   */
  private void generateRefConvexHull() {
    assert contacts.length >= 3;
    ArrayList<Triangle> trianglesRefCH = generateConvexHull(contacts, nRings);
    refConvexHull = new TriangleMesh(contacts, trianglesRefCH);
    refConvexHull.setupSwingLists();  // for fixing penetration
  }

  private Triangle generateThreeRingTriangleNaive(pt[] points0,
                                                  pt[] points1,
                                                  pt[] points2) {
    int i = 0, j = 0, k = 0;
    vec[] vs = new vec[6];
    for (i = 0; i < nPointsPerRing; ++i) {
      pt pa = points0[i];
      vs[0] = V(pa, points0[(i-1+nPointsPerRing)%nPointsPerRing]);
      vs[1] = V(pa, points0[(i+1)%nPointsPerRing]);
      for (j = 0; j < nPointsPerRing; ++j) {
        pt pb = points1[j];
        vs[2] = V(pb, points1[(j-1+nPointsPerRing)%nPointsPerRing]);
        vs[3] = V(pb, points1[(j+1)%nPointsPerRing]);
        for (k = 0; k < nPointsPerRing; ++k) {
          pt pc = points2[k];
          vs[4] = V(pc, points2[(k-1+nPointsPerRing)%nPointsPerRing]);
          vs[5] = V(pc, points2[(k+1)%nPointsPerRing]);
          vec normal = N(pa, pb, pc);
          boolean isValid = true;
          for (vec v : vs) {
            if (dot(normal, v) > 0) {
              isValid = false;
              break;
            }
          }
          if (isValid) return new Triangle(i, j, k);
        }
      }
    }
    println("return null even using naive method");
    return null;
  }

  private int findStablePoint(pt[] points, int i, pt a, pt b, vec normal) {
    vec vn = normal;
    int d = dot(vn, V(points[i], points[(i + 1) % nPointsPerRing])) > 0 ? 1 : -1;
    int inext = (i + d + nPointsPerRing) % nPointsPerRing;
    int iprev = (i - d + nPointsPerRing) % nPointsPerRing;
    int steps1 = 0;
    while ((dot(vn, V(points[i], points[inext])) > 0 ||
            dot(vn, V(points[i], points[iprev])) > 0) &&
            (steps1 < maxIterHSLocal)) {
      if (debug3RT &&
        debug3RTInfo.numFaces == numFaces3RT - 1 &&
        debug3RTInfo.numSteps >= numSteps3RT) break;
      iprev = i;
      i = inext;
      inext = (inext + d + nPointsPerRing) % nPointsPerRing;
      vn = N(points[i], a, b);
      if (debug3RT) {
        debug3RTInfo.numSteps++;
      }
      steps1++;
    }
    normal.setTo(vn);  // normal will be updated
    return i;
  }

  private boolean findStableTriangleHS(pt[] points0, pt[] points1, pt[] points2,
                                       Triangle tri, vec normal) {
    int i = tri.a, j = tri.b, k = tri.c;
    boolean update = false;
    /* Find ring A/B/C's stable point respectively. */
    int inew = findStablePoint(points0, i, points1[j], points2[k], normal);
    if (inew != i) {
      i = inew;
      if (update == false) update = true;
    }
    int jnew = findStablePoint(points1, j, points2[k], points0[i], normal);
    if (jnew != j) {
      j = jnew;
      if (update == false) update = true;
    }
    int knew = findStablePoint(points2, k, points0[i], points1[j], normal);
    if (knew != k) {
      k = knew;
      if (update == false) update = true;
    }
    tri.set(inew, jnew, knew);
    return update;
  }

  /*
   * Generate a three-ring triangle given 3 rings using heuristic search. It may
   * not converge when rings are large, rings have just a few points.
   */
  private Triangle generateThreeRingTriangleHS(pt[] points0,
                                               pt[] points1,
                                               pt[] points2) {
    int i = 0, j = 0, k = 0;
    Triangle tri = new Triangle(i, j, k);
    vec normal = N(points0[i], points1[j], points2[k]);
    if (debug3RT) {
      debug3RTInfo.idp0 = i;
      debug3RTInfo.idp1 = j;
      debug3RTInfo.idp2 = k;
      debug3RTInfo.numSteps = 1;
    }
    int iter = 0;
    while (iter < maxIterHSGlobal) {
      if (debug3RT &&
          debug3RTInfo.numFaces == numFaces3RT - 1 &&
          debug3RTInfo.numSteps >= numSteps3RT) break;
      boolean update = findStableTriangleHS(points0, points1, points2, tri, normal);
      if (update == false) break;  // break when stable
      iter++;
    }  // end while

    if (debug3RT) {
      debug3RTInfo.idp0 = i;
      debug3RTInfo.idp1 = j;
      debug3RTInfo.idp2 = k;
    }

    if (iter < maxIterHSGlobal) {
      return new Triangle(tri.a, tri.b, tri.c);
    } else {  // use backup method
      println("cannot find stable three-ring triangle using HS");
      numBackup3RT++;
      //return generateThreeRingTriangleNaive(points0, points1, points2);
      return null;
    }
  }

  /*
   * Generate a three-ring triangle given 3 rings using BFS method with a hint,
   * i.e. an initial triangle.
   */
  private Triangle generateThreeRingTriangleBFSWithHint(pt[] points0,
                                                        pt[] points1,
                                                        pt[] points2,
                                                        Triangle tri) {
    pt[][] rings = new pt[3][];
    rings[0] = points0;
    rings[1] = points1;
    rings[2] = points2;
    TriangleNode triNode = new TriangleNode(tri);
    Queue<TriangleNode> queue = new LinkedList<TriangleNode>();
    queue.add(triNode);
    HashSet<Triangle> set = new HashSet<Triangle>();
    set.add(tri);
    boolean stable = triNode.computeChildren(rings, set);
    if (stable) return tri;
    while (queue.size() > 0) {
      TriangleNode node = queue.poll();
      ArrayList<TriangleNode> children = node.children;
      int nc = children.size();
      for (int i = 0; i < nc; ++i) {
        TriangleNode child = children.get(i);
        stable = child.computeChildren(rings, set);
        if (stable) {
          // if (set.size() > nPointsPerRing * nPointsPerRing * nPointsPerRing) {
          //   println("size of explored set (bigger than np^3) = ", set.size());
          // }
          return child.tri;
        }
        queue.add(child);
      }
    }
    println("cannot find stable three-ring triangle using BFS");  // this is possible!
    numBackup3RT++;
    //return generateThreeRingTriangleNaive(points0, points1, points2);
    return null;
  }

  /*
   * Generate a three-ring triangle using BFS method. This method can be slow.
   * It is slower than naive method when the number of points on each ring is
   * small. However, with a good hint, this method can be fast.
   */
  private Triangle generateThreeRingTriangleBFS(pt[] points0,
                                                pt[] points1,
                                                pt[] points2,
                                                int option) {
    Triangle tri = null;
    switch (option) {
      case 0:  // randomly generate an initial triangle
        tri = new Triangle(int(random(points0.length)),
                           int(random(points1.length)),
                           int(random(points2.length)));
        break;
      case 1:  // heuristics
        vec normal = N(points0[0], points1[0], points2[0]);
        tri = new Triangle(0, 0, 0);
        findStableTriangleHS(points0, points1, points2, tri, normal);
        break;
      default:
        tri = new Triangle(0, 0, 0);
        break;
    }
    return generateThreeRingTriangleBFSWithHint(points0, points1, points2, tri);
  }

  /*
   * Generate a three-ring triangle using approximated extreme plane method.
   * This method first finds the exact extreme plane/triangle using the three
   * circles, and then trys to find a triangle defined by sampled points that
   * can approximate the exact extreme triangle.
   */
  private Triangle generateThreeRingTriangleApprox(int rid0, int rid1, int rid2) {
    pt[] cs = {centers[rid0], centers[rid1], centers[rid2]};
    float[] rs = {radii[rid0], radii[rid1], radii[rid2]};
    vec[] ns = {normals[rid0], normals[rid1], normals[rid2]};
    vec[] vis = {xAxes[rid0], xAxes[rid1], xAxes[rid2]};
    vec[] vjs = {yAxes[rid0], yAxes[rid1], yAxes[rid2]};
    pt[] ps = tangentPlaneThreeCirclesIter(cs[0], rs[0], ns[0], vis[0], vjs[0],
                                       cs[1], rs[1], ns[1], vis[1], vjs[1],
                                       cs[2], rs[2], ns[2], vis[2], vjs[2],
                                       null);

    float da = TWO_PI / nPointsPerRing;
    int[][] candidates = new int[3][2];
    for (int i = 0; i < 3; ++i) {
      vec cp = V(cs[i], ps[i]);
      float angle = acos(dot(cp, vis[i]) / rs[i]);
      if (dot(cp, vjs[i]) < 0) angle = TWO_PI - angle;
      float tmp = angle / da;
      int prev = int(tmp);
      int next = prev + 1;
      if (tmp - prev > next - tmp) {
        candidates[i][0] = next % nPointsPerRing;
        candidates[i][1] = prev % nPointsPerRing;
      } else {
        candidates[i][0] = prev % nPointsPerRing;
        candidates[i][1] = next % nPointsPerRing;
      }
    }

    pt[] ps0 = points[rid0];
    pt[] ps1 = points[rid1];
    pt[] ps2 = points[rid2];
    vec[] vs = new vec[6];
    for (int i = 0; i < 2; ++i) {
      int ia = candidates[0][i];
      pt pa = ps0[ia];
      vs[0] = V(pa, ps0[(ia - 1 + nPointsPerRing) % nPointsPerRing]);
      vs[1] = V(pa, ps0[(ia + 1) % nPointsPerRing]);
      for (int j = 0; j < 2; ++j) {
        int ib = candidates[1][j];
        pt pb = ps1[ib];
        vs[2] = V(pb, ps1[(ib - 1 + nPointsPerRing) % nPointsPerRing]);
        vs[3] = V(pb, ps1[(ib + 1) % nPointsPerRing]);
        for (int k = 0; k < 2; ++k) {
          int ic = candidates[2][k];
          pt pc = ps2[ic];
          vs[4] = V(pc, ps2[(ic - 1 + nPointsPerRing) % nPointsPerRing]);
          vs[5] = V(pc, ps2[(ic + 1) % nPointsPerRing]);
          vec n = N(pa, pb, pc);
          boolean isValid = true;
          for (vec v : vs) {
            if (dot(n, v) > 0) {
              isValid = false;
              break;
            }
          }
          if (isValid) return new Triangle(ia, ib, ic);
        }
      }
    }

    println("cannot find stable three-ring triangle using Approximated Extreme Plane method");
    numBackup3RT++;
    //return generateThreeRingTriangleNaive(ps0, ps1, ps2);
    return null;
  }

  private void setupSplitLists() {
    splitLists = new ArrayList<ArrayList<Integer>>();
    for (int i = 0; i < nRings; ++i) {
      ArrayList<Integer> swingList = refConvexHull.swingLists.get(i);
      int nCorners = swingList.size();  // number of adjacent corners/triangles
      splitLists.add(new ArrayList<Integer>());
      for (int j = 0; j < nCorners; ++j) {
        int cid = swingList.get(j);
        int tid = cid / 3;  // triangle ID w.r.t. ref convex hull/set of three-ring triangles
        int pid = threeRingTriangles.get(tid).get(cid % 3);  // point ID (global)
        splitLists.get(i).add(pid);
      }
    }
  }

  private boolean cornerIntersectsCorner(pt a0, pt a1, pt a2, pt b0, pt b1, pt b2) {
    return edgeIntersectsTriangle(a0, a1, b0, b1, b2) ||
           edgeIntersectsTriangle(a0, a2, b0, b1, b2) ||
           edgeIntersectsTriangle(b0, b1, a0, a1, a2) ||
           edgeIntersectsTriangle(b0, b2, a0, a1, a2);
  }

  // TODO: current implementation may be optimized, especially triangle-edge intersection
  private void fixPenetration3RT() {
    assert threeRingTriangles != null;
    setupSplitLists();
    for (int i = 0; i < nRings; ++i) {
      ArrayList<Integer> splits = splitLists.get(i);
      ArrayList<Integer> swings = refConvexHull.swingLists.get(i);
      int ns = splits.size();

      /* Pick the index with minimum location to start. */
      int first = argmin(splits);  // index, [0, ns)
      int firstSplit = splits.get(first);  // point ID (global)
      HashMap<Integer, Integer> splitCount = new HashMap<Integer, Integer>();
      splitCount.put(firstSplit, 1);
      ArrayList<pt> cornerTriangles = new ArrayList<pt>();

      /* Fix any issue around this index. */
      int cid0 = swings.get(first);
      int tid0 = cid0 / 3;
      Triangle tri0 = threeRingTriangles.get(tid0);
      int ia0 = tri0.get((cid0 % 3)), ib0 = tri0.get((cid0 + 1) % 3), ic0 = tri0.get((cid0 + 2) % 3);
      pt pa0 = getPoint(ia0), pb0 = getPoint(ib0), pc0 = getPoint(ic0);
      cornerTriangles.add(pa0);
      cornerTriangles.add(pb0);
      cornerTriangles.add(pc0);
      int start = (first + 1) % ns;
      int end = (first + ns - 1) % ns;
      int count = 1;

      // if (debugFixPenetration) {
      //   pt tmp = getPoint(tri0.get(cid0 % 3));
      //   fill(red, 100);
      //   show(tmp, 5);
      // }

      while (count < ns) {  // search right, i.e. increment start pointer if possible
        int cid = swings.get(start);
        Triangle tri = threeRingTriangles.get(cid / 3);
        int ia = tri.get(cid % 3), ib = tri.get((cid + 1) % 3), ic = tri.get((cid + 2) % 3);
        pt pa = getPoint(ia), pb = getPoint(ib), pc = getPoint(ic);
        int key = splits.get(start);
        if (splitCount.containsKey(key)) {
          splitCount.put(key, splitCount.get(key) + 1);
          cornerTriangles.add(pa);
          cornerTriangles.add(pb);
          cornerTriangles.add(pc);
          start = (start + 1) % ns;
          count++;
          continue;
        }

        int curSize = cornerTriangles.size();
        int j = curSize - 3;
        for (; j >= 0; j -= 3) {  // the order may help improve performance
          pt pd = cornerTriangles.get(j);
          pt pe = cornerTriangles.get(j+1);
          pt pf = cornerTriangles.get(j+2);
          if (cornerIntersectsCorner(pa, pb, pc, pd, pe, pf)) break;
        }
        if (j >= 0) {  // intersection found!
          splitCount.put(key, splitCount.getOrDefault(key, 0) + 1);
          cornerTriangles.add(pa);
          cornerTriangles.add(pb);
          cornerTriangles.add(pc);
          start = (start + 1) % ns;
          count++;
        } else {  // intersection not found!
          break;
        }
      }  // end searching right

      // if (debugFixPenetration) {
      //   int cid = swings.get(start);
      //   pt tmp = getPoint(threeRingTriangles.get(cid / 3).get(cid % 3));
      //   fill(green, 100);
      //   show(tmp, 5);
      // }

      while (count < ns) {  // search left, i.e. decrement end pointer if possible
        int cid = swings.get(end);
        Triangle tri = threeRingTriangles.get(cid / 3);
        int ia = tri.get(cid % 3), ib = tri.get((cid + 1) % 3), ic = tri.get((cid + 2) % 3);
        pt pa = getPoint(ia), pb = getPoint(ib), pc = getPoint(ic);
        int key = splits.get(end);
        if (splitCount.containsKey(key)) {
          splitCount.put(key, splitCount.get(key) + 1);
          cornerTriangles.add(pa);
          cornerTriangles.add(pb);
          cornerTriangles.add(pc);
          end = (end + ns - 1) % ns;
          count++;
          continue;
        }
        int curSize = cornerTriangles.size();
        int j = curSize - 3;
        for (; j >= 0; j -= 3) {  // the order may help improve performance
          pt pd = cornerTriangles.get(j);
          pt pe = cornerTriangles.get(j+1);
          pt pf = cornerTriangles.get(j+2);
          if (cornerIntersectsCorner(pa, pb, pc, pd, pe, pf)) break;
        }
        if (j >= 0) {  // intersection found!
          splitCount.put(key, splitCount.getOrDefault(key, 0) + 1);
          cornerTriangles.add(pa);
          cornerTriangles.add(pb);
          cornerTriangles.add(pc);
          end = (end + ns - 1) % ns;
          count++;
        } else {  // intersection not found!
          break;
        }
      }  // end searching left

      // if (debugFixPenetration) {
      //   int cid = swings.get(end);
      //   pt tmp = getPoint(threeRingTriangles.get(cid / 3).get(cid % 3));
      //   fill(blue, 100);
      //   show(tmp, 5);
      // }

      /* Fix the first group. */
      int newSplit = firstSplit;  // make all corners/triangles touch newSplit
      int maxCount = splitCount.get(firstSplit);
      for (Integer k : splitCount.keySet()) {
        if (splitCount.get(k) > maxCount) {
          maxCount = splitCount.get(k);
          newSplit = k;
        }
      }
      // for (int j = 0; j < splitIdxs.size(); ++j) {
      for (int idx = (end + 1) % ns; idx != start; idx = (idx + 1) % ns) {
        int cid = swings.get(idx);
        int tid = cid / 3;
        int ccid = cid % 3;
        if (threeRingTriangles.get(tid).get(ccid) != newSplit) {
          splits.set(idx, newSplit);
          threeRingTriangles.get(tid).set(ccid, newSplit);
        }
      }

      /* Fix the rest. */
      int prevSplit = splits.get(start);
      start = (start + 1) % ns;
      int curSplit = splits.get(start);
      while (count < ns - 1) {
        if (curSplit < prevSplit) {
          splits.set(start, prevSplit);  // increase cur split to match prev split
          int cid = swings.get(start);
          threeRingTriangles.get(cid / 3).set(cid % 3, prevSplit);
        } else prevSplit = curSplit;
        start = (start + 1) % ns;
        curSplit = splits.get(start);
        count++;
      }
    }
  }

  /*
   * Generate three-ring triangles. A three-ring triangle is an "extreme"
   * triangle that touches three rings and have them on the same side.
   */
  void generateThreeRingTriangles() {
    if (refConvexHull == null) generateRefConvexHull();
    int nt = refConvexHull.nt;  // number of triangles of convex hull
    if (threeRingTriangles != null) {
      threeRingTriangles.clear();
    } else {
      threeRingTriangles = new ArrayList<Triangle>();
    }
    if (debug3RT && method3RT == 1) debug3RTInfo.numFaces = 0;
    for (int i = 0; i < nt; ++i) {
      if (debug3RT && method3RT == 1 && i >= numFaces3RT) break;
      Triangle face = refConvexHull.triangles.get(i);
      Triangle t = null;
      switch (method3RT) {
        case 1:
          t = generateThreeRingTriangleHS(points[face.a], points[face.b], points[face.c]);
          break;
        case 2:
          t = generateThreeRingTriangleBFS(points[face.a], points[face.b], points[face.c], 0);
          break;
        case 3:
          t = generateThreeRingTriangleBFS(points[face.a], points[face.b], points[face.c], 1);
          break;
        case 4:
          t = generateThreeRingTriangleApprox(face.a, face.b, face.c);
          break;
        default:
          t = generateThreeRingTriangleNaive(points[face.a], points[face.b], points[face.c]);
      }
      Triangle triangle = (t == null) ?
                          null :
                          new Triangle(face.a * nPointsPerRing + t.a, face.b * nPointsPerRing + t.b, face.c * nPointsPerRing + t.c);
      threeRingTriangles.add(triangle);
      if (debug3RT && method3RT == 1) {
        debug3RTInfo.idr0 = face.a;
        debug3RTInfo.idr1 = face.b;
        debug3RTInfo.idr2 = face.c;
        debug3RTInfo.numFaces++;
      }
    }  // end for

    if (fix3RT) {
      fixPenetration3RT();
    }
  }

  // TODO: need to revisit, bacause current method is quite simple and may create
  // overlapping triangles (with flipped normals)
  private void fillCorridor(int loA, int hiA, int loB, int hiB, int rA, int rB) {
    int i = loA, j = loB;
    pt pa = getPoint(loA), pb = getPoint(loB);
    int iBase = rA * nPointsPerRing, jBase = rB * nPointsPerRing;
    while(i != hiA || j != hiB) {
      pt pc = null, pd = null;
      int iNext = -1, jNext = -1;
      if (i != hiA) {  // increase i
        iNext = iBase + (i + 1) % nPointsPerRing;
        pc = getPoint(iNext);
      }
      if (j != hiB) {  // decrease j
        jNext = jBase + (j + nPointsPerRing - 1) % nPointsPerRing;
        pd = getPoint(jNext);
      }
      boolean moveA = true;
      if (pc != null && pd != null) {
        /* TODO: the criterion used to select a good triangle may be changed. */
        if (d(pa, pd) < d(pb, pc)) {
          twoRingTriangles.add(new Triangle(i, j, jNext));
          moveA = false;
        } else {
          twoRingTriangles.add(new Triangle(i, j, iNext));
        }
      } else if (pc != null) {
        twoRingTriangles.add(new Triangle(i, j, iNext));
      } else {  // pd != null
        twoRingTriangles.add(new Triangle(i, j, jNext));
        moveA = false;
      }
      if (moveA) i = iNext;
      else j = jNext;
    }
  }

  /*
   * Generate two-ring triangles. A two-ring triangle is a triangle that has two
   * vertices on a ring and another vertex on another ring.
   */
  void generateTwoRingTriangles() {
    twoRingTriangles = new ArrayList<Triangle>();
    boolean[][] visited = new boolean[nRings][nRings];
    for (int i = 0; i < nRings; ++i) {  // i = ID of ring A
      if (debug2RT && i >= debug2RTInfo.numGlobalStep) break;
      ArrayList<Integer> swings = refConvexHull.swingLists.get(i);
      int ns = swings.size();
      ArrayList<Integer> splits = splitLists.get(i);
      for (int j = 0, curC = swings.get(0); j < ns; ++j) {
        if (debug2RT && i == debug2RTInfo.numGlobalStep - 1 && j >= debug2RTInfo.numLocalStep) break;
        int swingC = swings.get((j+1) % ns);
        int prevC = prevCorner(curC);
        int nextSwingC = nextCorner(swingC);
        int k = refConvexHull.triangles.get(prevC / 3).get(prevC % 3);  // k = ID of ring B
        if (visited[i][k]) {
          curC = swingC;
          continue;
        }
        int loA = splits.get(j);
        int hiA = splits.get((j+1) % ns);
        int loB = threeRingTriangles.get(prevC / 3).get(prevC % 3);
        int hiB = threeRingTriangles.get(nextSwingC / 3).get(nextSwingC % 3);
        fillCorridor(loA, hiA, loB, hiB, i, k);
        curC = swingC;
        visited[i][k] = visited[k][i] = true;

        if (debug2RT && i == debug2RTInfo.numGlobalStep - 1 && j == debug2RTInfo.numLocalStep - 1) {
          debug2RTInfo.pa0 = getPoint(loA);
          debug2RTInfo.pa1 = getPoint(hiA);
          debug2RTInfo.pb0 = getPoint(loB);
          debug2RTInfo.pb1 = getPoint(hiB);
        }
      }
    }
  }

  private void generateTriangleMeshCH() {
    triangles = generateConvexHull(points, nRings, nPointsPerRing);
  }

  private void generateTriangleMeshFast() {
    if (!debug2RT || show2RT) fix3RT = true;
    generateThreeRingTriangles();
    if (!debug2RT || show2RT) {
      generateTwoRingTriangles();
    } else {
      twoRingTriangles = null;
    }
    triangles = new ArrayList<Triangle>();
    triangles.addAll(threeRingTriangles);
    if (twoRingTriangles != null) triangles.addAll(twoRingTriangles);
  }

  /*
   * Generate a triangle mesh based on the given option. When option is 1, use
   * 3RT and 2RT trick to generate the mesh. Otherwise, a convex hull will be
   * generated.
   */
  void generateTriangleMesh(int option) {
    switch (option) {
      case 1:
        generateTriangleMeshFast();
        break;
      default:
        generateTriangleMeshCH();
    }
  }

  /*
   * Generate an extreme triangle touching three rings given the three ring IDs.
   * Two such extreme triangles will be generated, i.e. 6 points will be returned.
   */
  private pt[] generateExTriThreeRingsTwo(int i, int j, int k) {
    float s1 = radii[i] / r, s2 = radii[j] / r, s3 = radii[k] / r;
    float a1 = asin(s1), a2 = asin(s2), a3 = asin(s3);  // TODO: may store these half-angles?
    float c1 = cos(a1), c2 = cos(a2), c3 = cos(a3);

    vec n23 = N(normals[j], normals[k]);
    vec n31 = N(normals[k], normals[i]);
    vec n12 = N(normals[i], normals[j]);
    float mix = d(normals[i], n23);
    if (isAbsZero(mix)) {
      println("mix product is 0");
      return null;
    }
    float invMix = 1.0 / mix;
    n23 = V(invMix, n23);
    n31 = V(invMix, n31);
    n12 = V(invMix, n12);
    vec va = V(c1, n23, c2, n31, c3, n12);
    vec vb = V(s1, n23, s2, n31, s3, n12);
    float bb = dot(vb, vb);
    float bb1 = bb - 1;
    if (isAbsZero(bb1)) {
      println("dot(b, b) = 1");
      return null;
    }
    float aa = dot(va, va);
    float ab = dot(va, vb);
    float det = ab * ab - (aa - 1) * bb1;
    if (det < 0) {
      println("determinant < 0");
      return null;
    }
    float sqrtDet = sqrt(det);
    float[] ts = new float[2];
    ts[0] = (ab + sqrtDet) / bb1;
    ts[1] = (ab - sqrtDet) / bb1;
    vec[] ns = new vec[2];  // the 2 normals
    pt[] ps = new pt[2];  // the 2 centers

    /* Compute the 2 centers and 2 normals. */
    for (int ii = 0; ii < 2; ++ii) {
      float cosa = 1.0 / sqrt(1 + ts[ii] * ts[ii]);  // cosa > 0
      ns[ii] = V(cosa, A(va, -ts[ii], vb));
      ps[ii] = P(c, r * cosa, ns[ii]);
    }

    /*
     * For each computed circle, check if all input circles on its negative side.
     * Since all input circles are already on the same side, we just need to
     * check if one input circle is on its negative side.
     */
    for (int ii = 0; ii < 2; ++ii) {
      vec v = V(ps[ii], centers[i]);
      if (dot(v, ns[ii]) > 0) ns[ii].rev();
    }

    pt[] pointsTangency = new pt[6];
    // first triangle (i, j, k)
    vec n1 = U(N(normals[i], ns[0]));  // the tangent at the contact point
    vec d1 = N(n1, normals[i]);
    pointsTangency[0] = P(centers[i], radii[i], d1);
    vec n2 = U(N(normals[j], ns[0]));
    vec d2 = N(n2, normals[j]);
    pointsTangency[1] = P(centers[j], radii[j], d2);
    vec n3 = U(N(normals[k], ns[0]));
    vec d3 = N(n3, normals[k]);
    pointsTangency[2] = P(centers[k], radii[k], d3);
    // second triangle (i, j, k)
    vec n4 = U(N(normals[i], ns[1]));
    vec d4 = N(n4, normals[i]);
    pointsTangency[3] = P(centers[i], radii[i], d4);
    vec n5 = U(N(normals[j], ns[1]));
    vec d5 = N(n5, normals[j]);
    pointsTangency[4] = P(centers[j], radii[j], d5);
    vec n6 = U(N(normals[k], ns[1]));
    vec d6 = N(n6, normals[k]);
    pointsTangency[5] = P(centers[k], radii[k], d6);

    /* Make sure the first center/normal corresponds to orientation (i, j, k). */
    if (dot(N(pointsTangency[0], pointsTangency[1], pointsTangency[2]), ns[0]) > 0) {
      pt p = pointsTangency[4];
      pointsTangency[4] = pointsTangency[5];
      pointsTangency[5] = p;
    } else {  // flip the first triplet, swap the two triplets
      println("swap the two triplets");
      pt[] tmpPointsTangency = new pt[6];
      tmpPointsTangency[0] = pointsTangency[3];
      tmpPointsTangency[1] = pointsTangency[4];
      tmpPointsTangency[2] = pointsTangency[5];
      tmpPointsTangency[3] = pointsTangency[0];
      tmpPointsTangency[4] = pointsTangency[2];
      tmpPointsTangency[5] = pointsTangency[1];
      pointsTangency = tmpPointsTangency;
    }

    if (debugST) {
      tan0 = ts[0];
      tan1 = ts[1];
      for (int ii = 0; ii < 2; ++ii) {
        debugSTInfo.circumcenters.add(ps[ii]);
        float a = atan(ts[ii]);  // [-PI/2, PI/2]
        if (a < 0) {
          a += PI;  // [PI/2, PI]
          System.out.format("%d: half angle > PI/2\n", ii);
          fill(pink);
          show(ps[ii], 5);
        }
        float radius = r * sin(a);
        assert radius >= 0;
        debugSTInfo.circumradii.add(radius);
        debugSTInfo.normals.add(ns[ii]);
      }
    }

    return pointsTangency;
  }

  /*
   * Generate an extreme triangle touching three rings given the three ring IDs.
   * TODO: this function is deprecated.
   */
  pt[] generateExTriThreeRings(int i, int j, int k, DebugEPInfo dInfo) {
    float s1 = radii[i] / r, s2 = radii[j] / r, s3 = radii[k] / r;
    float a1 = asin(s1), a2 = asin(s2), a3 = asin(s3);  // TODO: may store these half-angles?
    float c1 = cos(a1), c2 = cos(a2), c3 = cos(a3);
    float x1 = normals[i].x, y1 = normals[i].y, z1 = normals[i].z;
    float x2 = normals[j].x, y2 = normals[j].y, z2 = normals[j].z;
    float x3 = normals[k].x, y3 = normals[k].y, z3 = normals[k].z;

    float DD = -x3*y2*z1 + x2*y3*z1 + x3*y1*z2 - x1*y3*z2 - x2*y1*z3 + x1*y2*z3;

    if (isAbsZero(DD)) {
      println("3 normals are coplanar.");
      return null;
    }

    float A1 = (-c3*y2*z1 + c2*y3*z1 + c3*y1*z2 - c1*y3*z2 - c2*y1*z3 + c1*y2*z3) / DD;
    float B1 = (s3*y2*z1 - s2*y3*z1 - s3*y1*z2 + s1*y3*z2 + s2*y1*z3 - s1*y2*z3) / DD;
    float A2 = (c3*x2*z1 - c2*x3*z1 - c3*x1*z2 + c1*x3*z2 + c2*x1*z3 - c1*x2*z3) / DD;
    float B2 = (-s3*x2*z1 + s2*x3*z1 + s3*x1*z2 - s1*x3*z2 - s2*x1*z3 + s1*x2*z3) / DD;
    float A3 = (-c3*x2*y1 + c2*x3*y1 + c3*x1*y2 - c1*x3*y2 - c2*x1*y3 + c1*x2*y3) / DD;
    float B3 = (s3*x2*y1 - s2*x3*y1 - s3*x1*y2 + s1*x3*y2 + s2*x1*y3 - s1*x2*y3) / DD;

    float AB = 2 * (A1*B1 + A2*B2 + A3*B3);
    float AA = A1*A1 + A2*A2 + A3*A3 - 1.0;
    float BB = B1*B1 + B2*B2 + B3*B3 - 1.0;

    float[] sols = solveQuadraticEquation(BB, AB, AA);
    if (sols == null) {
      println("No solutions found!");
      return null;
    }
    float t = sols[0];
    float t2 = t * t;
    float cosa = 1 / sqrt(1 + t2), sina = t / sqrt(1 + t2);
    vec n = V(cosa * A1 + sina * B1, cosa * A2 + sina * B2, cosa * A3 + sina * B3);
    pt p = P(c, r * cosa, n);  // center of the tangent circle
    vec nijk = N(centers[i], centers[j], centers[k]);
    if (dot(nijk, V(centers[i], p)) < 0) {
      t = sols[1];
      t2 = t * t;
      cosa = 1 / sqrt(1 + t2);
      sina = t / sqrt(1 + t2);
      n.set(cosa * A1 + sina * B1, cosa * A2 + sina * B2, cosa * A3 + sina * B3);
      p = P(c, r * cosa, n);
    }
    // println("t = ", t);
    if (dot(nijk, n) < 0) n.rev();
    // System.out.format("n = (%f, %f, %f), p = (%f, %f, %f)\n", n.x, n.y, n.z, p.x, p.y, p.z);

    // Find contact points on the 3 circles
    pt[] ps = new pt[3];
    vec n1 = U(N(normals[i], n));
    vec d1 = N(n1, normals[i]);
    ps[0] = P(centers[i], radii[i], d1);
    vec n2 = U(N(normals[j], n));
    vec d2 = N(n2, normals[j]);
    ps[1] = P(centers[j], radii[j], d2);
    vec n3 = U(N(normals[k], n));
    vec d3 = N(n3, normals[k]);
    ps[2] = P(centers[k], radii[k], d3);

    {
      // float radius = r * sina;
      // fill(cyan, 100);
      // disk(p, n, radius);
      // fill(dgreen, 100);
      // arrow(p, V(50, n), 5);

      // fill(magenta, 100);
      // arrow(centers[i], V(radii[i], d1), 5);
      // arrow(centers[j], V(radii[j], d2), 5);
      // arrow(centers[k], V(radii[k], d3), 5);
    }
    return ps;
  }

  private boolean isValidExTri(pt pa, pt pb, pt pc, int ia, int ib, int ic) {
    vec normal = N(pa, pb, pc);

    /* Make sure all circle centers on the negative side of the current plane. */
    for (int i = 0; i < nRings; ++i) {
      if (i == ia || i == ib || i == ic) continue;
      if (dot(V(pa, centers[i]), normal) > 0) return false;
    }

    /* Make sure the cone defined by the circumcircle of the current triangle is empty. */
    vec va = V(pc, pa);
    vec vb = V(pc, pb);
    float radius = sqrt(n2(va) * n2(vb) * n2(M(va, vb)) / n2(normal)) / 2;  // radius of circumcircle
    float alpha = asin(radius / r);  // half angle
    vec unitNormal = U(normal);

    for (int i = 0; i < nRings; ++i) {
      if (i == ia || i == ib || i == ic) continue;
      float beta = asin(radii[i] / r);
      if (cos(alpha + beta) < dot(unitNormal, normals[i])) return false;
    }

    return true;
  }

  private void updateSwingList(int i, pt p, int pid) {
    assert swingLists[i] != null && angleLists[i] != null;

    vec vp = V(centers[i], p);
    float angle = acos(dot(vp, xAxes[i]) / radii[i]);
    if (dot(vp, yAxes[i]) < 0) angle = TWO_PI - angle;

    int n = angleLists[i].size();

    int lo = 0, hi = n - 1, mid;
    while (lo <= hi) {  // binary search
      mid = (lo + hi) / 2;
      if (angleLists[i].get(mid) > angle) hi = mid - 1;
      else lo = mid + 1;
    }
    angleLists[i].add(lo, angle);
    swingLists[i].add(lo, pid);
  }

  private void insertExTriPoint(int rid, pt p, int pid) {
    exTriPoints.add(p);
    exTriRIDs.add(rid);
    updateSwingList(rid, p, pid);
  }

  void generateExTris() {
    exTriPoints = new ArrayList<pt>();
    exTriRIDs = new ArrayList<Integer>();
    swingLists = new ArrayList[nRings];
    angleLists = new ArrayList[nRings];
    for (int i = 0; i < nRings; ++i) {
      swingLists[i] = new ArrayList<Integer>();
      angleLists[i] = new ArrayList<Float>();
    }
    if (debugST) {
      debugSTInfo.circumcenters.clear();
      debugSTInfo.circumradii.clear();
      debugSTInfo.normals.clear();
    }
    for (int i = 0; i < nRings; ++i) {
      for (int j = i + 1; j < nRings; ++j) {
        for (int k = j + 1; k < nRings; ++k) {
          // pt[] tmp = generateExTriThreeRings(i, j, k, null);
          // if (isValidExTri(tmp[0], tmp[1], tmp[2], i, j, k)) {
          //   int pid = exTriPoints.size();
          //   exTriPoints.add(tmp[0]);
          //   exTriRIDs.add(i);
          //   updateSwingList(i, tmp[0], pid);

          //   exTriPoints.add(tmp[1]);
          //   exTriRIDs.add(j);
          //   updateSwingList(j, tmp[1], pid + 1);

          //   exTriPoints.add(tmp[2]);
          //   exTriRIDs.add(k);
          //   updateSwingList(k, tmp[2], pid + 2);
          // }
          // tmp = generateExTriThreeRings(i, k, j, null);
          // if (isValidExTri(tmp[0], tmp[1], tmp[2], i, k, j)) {
          //   int pid = exTriPoints.size();
          //   exTriPoints.add(tmp[0]);
          //   exTriRIDs.add(i);
          //   updateSwingList(i, tmp[0], pid);

          //   exTriPoints.add(tmp[1]);
          //   exTriRIDs.add(k);
          //   updateSwingList(k, tmp[1], pid + 1);

          //   exTriPoints.add(tmp[2]);
          //   exTriRIDs.add(j);
          //   updateSwingList(j, tmp[2], pid + 2);
          // }

          pt[] tmp = generateExTriThreeRingsTwo(i, j, k);
          if (tmp == null || tmp.length != 6) {
            System.out.format("No ex tris among (%d, %d, %d)\n", i, j, k);
          }
          if (isValidExTri(tmp[0], tmp[1], tmp[2], i, j, k)) {
            int pid = exTriPoints.size();
            insertExTriPoint(i, tmp[0], pid);
            insertExTriPoint(j, tmp[1], pid + 1);
            insertExTriPoint(k, tmp[2], pid + 2);
          }
          if (isValidExTri(tmp[3], tmp[4], tmp[5], i, k, j)) {
            int pid = exTriPoints.size();
            insertExTriPoint(i, tmp[3], pid);
            insertExTriPoint(k, tmp[4], pid + 1);
            insertExTriPoint(j, tmp[5], pid + 2);
          }
        }
      }
    }
  }

  private ArrayList<pt> fillExEdgeWithPoints(float angle, float da, int n, int u, int v) {
    ArrayList<pt> edgePoints = new ArrayList<pt>();
    pt cu = centers[u], cv = centers[v];
    float ru = radii[u], rv = radii[v];
    vec nu = normals[u], nv = normals[v];
    vec viu = xAxes[u], viv = xAxes[v];
    vec vju = yAxes[u], vjv = yAxes[v];
    for (int k = 1; k < n; ++k) {
      angle += da;
      float c = cos(angle), s = sin(angle);
      pt p0 = P(cu, ru * c, viu, ru * s, vju);  // sampled point
      pt p1 = P(p0, -ru * s, viu, ru * c, vju);  // sampled point + tangent at that point
      pt[] candidates = pivotPlaneAroundLineHitCircle(cv, rv, nv, p0, p1, viv, vjv);
      pt candidate = candidates[0];
      if (dot(V(p0, cv), N(p0, candidate, p1)) > 0) {
        candidate = candidates[1];
      }
      edgePoints.add(p0);
      edgePoints.add(candidate);
    }
    return edgePoints;
  }

  private void fillExEdge(int u, int v, int ia, int ib, int ic, int id) {
    pt pa = exTriPoints.get(ia);
    pt pb = exTriPoints.get(ib);
    pt pc = exTriPoints.get(ic);
    pt pd = exTriPoints.get(id);

    pt cu = centers[u], cv = centers[v];
    float ru = radii[u], rv = radii[v];
    vec nu = normals[u], nv = normals[v];

    // compute angle for arc DA and arc BC, w.r.t. corresponding ring orientations
    vec va = V(cu, pa), vd = V(cu, pd);
    vec vb = V(cv, pb), vc = V(cv, pc);
    float au = acos(dot(va, vd) / (ru * ru));  // [0, PI]
    if (dot(nu, N(vd, va)) < 0) au = TWO_PI - au;
    float av = acos(dot(vb, vc) / (rv * rv));
    if (dot(nv, N(vb, vc)) < 0) av = TWO_PI - av;

    // Sample points from the bigger arc
    if (au > av) {
      float angle = acos(dot(vd, xAxes[u]) / ru);  // starting angle for D w.r.t. ring u
      if (dot(vd, yAxes[u]) < 0) angle = TWO_PI - angle;
      int n = int(au / dAngle);  // split arc into n pieces by inserting n-1 points
      ArrayList<pt> edgePoints = fillExEdgeWithPoints(angle, dAngle, n, u, v);
      exEdges.put(new Corridor(id, ic, ia, ib), edgePoints);
    } else {
      float angle = acos(dot(vb, xAxes[v]) / rv);
      if (dot(vb, yAxes[v]) < 0) angle = TWO_PI - angle;
      int n = int(av / dAngle);
      ArrayList<pt> edgePoints = fillExEdgeWithPoints(angle, dAngle, n, v, u);
      exEdges.put(new Corridor(ib, ia, ic, id), edgePoints);
    }
  }

  void generateExEdges() {
    assert nSamples >= 1;
    dAngle = TWO_PI / nSamples;
    int nExTriPoints = exTriPoints.size();
    boolean[][] visited = new boolean[nExTriPoints][nExTriPoints];
    exEdges = new LinkedHashMap<Corridor, ArrayList<pt>>();
    for (int i = 0; i < nRings; ++i) {
      ArrayList<Integer> swingList = swingLists[i];
      int n = swingList.size();
      for (int j = 0; j < n; ++j) {
        int pidD = swingList.get(j);
        int pidC = int(pidD / 3) * 3 + (pidD + 2) % 3;
        if (visited[pidD][pidC]) continue;
        int pidA = swingList.get((j + 1) % n);
        int pidB = int(pidA / 3) * 3 + (pidA + 1) % 3;
        int u = exTriRIDs.get(pidD);
        int v = exTriRIDs.get(pidC);
        fillExEdge(u, v, pidA, pidB, pidC, pidD);
        visited[pidD][pidC] = visited[pidC][pidD] = true;
        visited[pidA][pidB] = visited[pidB][pidA] = true;
      }
    }
  }

  void generateExactCH() {
    generateExTris();
    generateExEdges();
  }

  void showRings() {
    noStroke();
    fill(orange);
    for (int i = 0; i < nRings; ++i) {
      show(centers[i], 1);
    }
    fill(green);
    for (int i = 0; i < nRings; ++i) {
      arrow(centers[i], V(centers[i], points[i][0]), 2);
    }
    fill(blue);
    for (int i = 0; i < nRings; ++i) {
      for (int j = 0; j < nPointsPerRing; ++j) {
        show(points[i][j], 2);
      }
    }
    fill(cyan);
    for (int i = 0; i < nRings; ++i) {
      for (int j = 0; j < nPointsPerRing; ++j) {
        collar(points[i][j], V(points[i][j], points[i][(j + 1) % nPointsPerRing]), 1, 1);
      }
    }
    fill(purple);
    for (int i = 0; i < nRings; ++i) {
      arrow(centers[i], V(20, normals[i]), 4);
    }
    return;
  }

  void showCircles() {
    stroke(0);
    strokeWeight(3);
    noFill();
    float da = TWO_PI / 36;
    for (int i = 0; i < nRings; ++i) {
      float a = 0;
      beginShape();
      for (int j = 0; j < 36; ++j, a += da) {
        vertex(P(centers[i], radii[i] * cos(a), xAxes[i], radii[i] * sin(a), yAxes[i]));
      }
      endShape(CLOSE);
    }
  }

  void showDisks() {
    noStroke();
    fill(red, 200);
    for (int i = 0; i < nRings; ++i) {
      disk(centers[i], xAxes[i], yAxes[i], radii[i]);
    }
  }

  void showExEdges() {
    if (exTriPoints == null || exEdges == null) return;
    stroke(0);
    strokeWeight(2);
    fill(green, 200);
    for (Corridor edge : exEdges.keySet()) {
      // System.out.format("(%d, %d, %d, %d)\n", edge.a, edge.b, edge.c, edge.d);
      pt pa = exTriPoints.get(edge.a);
      pt pb = exTriPoints.get(edge.b);
      pt pc = exTriPoints.get(edge.c);
      pt pd = exTriPoints.get(edge.d);
      ArrayList<pt> edgePoints = exEdges.get(edge);
      beginShape(QUAD_STRIP);
      vertex(pa);
      vertex(pb);
      for (pt p : edgePoints) vertex(p);
      vertex(pc);
      vertex(pd);
      endShape(QUAD_STRIP);
    }
  }

  void showExTris() {
    if (exTriPoints == null) return;
    noStroke();
    fill(blue, 200);
    beginShape(TRIANGLES);
    for (int i = 0; i < exTriPoints.size(); i += 3) {
      vertex(exTriPoints.get(i));
      vertex(exTriPoints.get(i+1));
      vertex(exTriPoints.get(i+2));
    }
    endShape();
  }

  void showDebug3RTInfo() {
    noStroke();
    fill(red, 150);
    show(contacts[debug3RTInfo.idr0], 3);
    fill(green, 150);
    show(contacts[debug3RTInfo.idr1], 3);
    fill(blue, 150);
    show(contacts[debug3RTInfo.idr2], 3);
    fill(#BF6868, 200);  // light red
    show(points[debug3RTInfo.idr0][debug3RTInfo.idp0], 5);
    fill(#40935D, 200);  // light green
    show(points[debug3RTInfo.idr1][debug3RTInfo.idp1], 5);
    fill(#517EC9, 200);  // light blue
    show(points[debug3RTInfo.idr2][debug3RTInfo.idp2], 5);
  }

  void showDebug2RTInfo() {
    fill(red, 200); show(debug2RTInfo.pa0, 5);
    fill(yellow, 200); show(debug2RTInfo.pa1, 5);
    fill(green, 100); show(debug2RTInfo.pb0, 5);
    fill(blue, 100); show(debug2RTInfo.pb1, 5);
    fill(#8B7373, 200); show(centers[debug2RTInfo.numGlobalStep - 1], 5);  // center of current ring
  }

  void save(String file) {
    String[] lines = new String[2 + 3 * nRings];
    int i = 0;
    lines[i++] = str(nRings);
    lines[i++] = str(nPointsPerRing);
    for (int j = 0; j < nRings; ++j) {
      lines[i++] = str(contacts[j].x) + "," + str(contacts[j].y) + "," +
                   str(contacts[j].z);
      lines[i++] = str(radii[j]);
      lines[i++] = str(xAxes[j].x) + "," + str(xAxes[j].y) + "," +
                   str(xAxes[j].z);
    }
    saveStrings(file, lines);
    return;
  }

  void load(String file) {
    String[] lines = loadStrings(file);
    int i = 0;
    nRings = int(lines[i++]);
    nPointsPerRing = int(lines[i++]);
    println("loading:", file, "nc =", nRings, "np =", nPointsPerRing);
    contacts = new pt[nRings];
    radii = new float[nRings];
    xAxes = new vec[nRings];
    normals = new vec[nRings];
    yAxes = new vec[nRings];
    for (int j = 0; j < nRings; ++j) {
      float[] contact = float(split(lines[i++], ","));
      contacts[j] = new pt(contact[0], contact[1], contact[2]);
      normals[j] = U(c, contacts[j]);
      radii[j] = float(lines[i++]);
      float[] initDir = float(split(lines[i++], ","));
      xAxes[j] = new vec(initDir[0], initDir[1], initDir[2]);
      yAxes[j] = N(normals[j], xAxes[j]);
    }
    return;
  }

  boolean isValid() {
    for (int i = 0; i < nRings; ++i) {
      float alpha0 = asin(clamp(radii[i]/r, -1.0, 1.0));
      for (int j = i + 1; j < nRings; ++j) {
        float alpha1 = asin(clamp(radii[j]/r, -1.0, 1.0));
        float theta = acos(clamp(dot(normals[i], normals[j]), -1.0, 1.0));
        if (theta <= alpha0 + alpha1) return false;
      }
    }
    return true;
  }
}