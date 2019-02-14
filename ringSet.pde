/******************************************************************************
 * Ring set processing.
 ******************************************************************************/


import java.util.Queue;
import java.util.LinkedHashMap;
import java.util.Collections;

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

int idxIncCor = 0;  // the ID of the interested corridor (used in debug)
int idxIncTri = 0;  // the ID of the interested triangle (used in debug)

boolean showRingSet = true;
boolean showCircleSet = false;
boolean showDiskSet = false;
boolean showCones = false;
boolean showPolygons = false;

boolean debugIncCH = false;
int debugIncCHIter = 3;
int debugIncCHBoundarySize = 2;
boolean debugIncCHNewView = false;
boolean debugIncCHCor = false;
boolean showCorridorFaces = false;
boolean showTriangleFaces = false;
boolean showCorridorStrokes = false;
boolean showTriangleStrokes = true;
boolean showBeams = false;
boolean showApolloniusDiagram = false;
boolean debugApolloniusDiagram = false;
boolean showTriMesh = false;  // the mesh generated from exact convex hull

boolean simpleCorridor = false;

boolean validRS = false;

/*
 * Method for three-ring triangle generation.
 * 0: Naive method. Time: O(n^3). Space: O(1)
 * 1: Heuristic search. Time: O(n). Space: O(1)
 * 2: Breadth first search. Time: roughly O(min(3^d, n^3)) where d is the length
      from initial state to optimal state. Space: O(min(3^d, n^3)). Not converge?
 * 3: Breadth first search with heuristics. Time and space should be less than 2.
 * 4: Approximated supporting plane. Time: roughly O(1), Space: roughly O(1).
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
  class NaiveCorridor {
    int a, b, c, d;  // point IDs
    NaiveCorridor(int a, int b, int c, int d) {
      this.a = a;
      this.b = b;
      this.c = c;
      this.d = d;
    }
  }

  class IncVertex {
    pt position;
    int rid;
    IncVertex(pt position, int rid) {
      this.position = position;
      this.rid = rid;
    }
  }

  class IncFace {
    IncVertex[] vertices = null;
    /*
     * -1: unvisited, 0: unvisible, 1: visible, 2: unvisible but should be removed
     * For example, a corridor may not be visible but if its adjacent triangle is
     * removed, it should also be removed.
     */
    int visible = -1;
    IncFace() {}
    boolean isVisibleFromCircle(vec v, float f) {
      return false;
    }
    IncEdge[] getEdges() {
      return null;
    }
    void showFace() {}

    /*
     * Find edges in the current face that are complemental to an edge e from an
     * adjacent face.
     */
    IncEdge[] getComplementalEdges(IncEdge e) {
      return null;
    }
  }

  class IncEdge {
    IncVertex va;
    IncVertex vb;
    IncFace adjFace = null;
    IncEdge() {}
    IncEdge(IncVertex va, IncVertex vb) {
      this.va = va;
      this.vb = vb;
    }
    void setEndPoints(IncVertex va, IncVertex vb) {
      this.va = va;
      this.vb = vb;
    }
    IncFace getAdjFace() {
      return adjFace;
    }
    void setAdjFace(IncFace adjFace) {
      this.adjFace = adjFace;
    }

    void showEdge() {
      noStroke();
      arrow(va.position, V(va.position, vb.position), 3);
    }
  }

  class IncTriangle extends IncFace {
    IncEdge[] edges;
    vec normal = null;
    float angle = -1.0;
    pt pAD = null;  // the position of the corresponding Apollonius diagram vertex

    IncTriangle(IncVertex v0, IncVertex v1, IncVertex v2) {
      vertices = new IncVertex[3];
      edges = new IncEdge[3];
      for (int i = 0; i < 3; ++i) edges[i] = new IncEdge();
      vertices[0] = v0;
      vertices[1] = v1;
      vertices[2] = v2;
      edges[0].setEndPoints(v0, v1);
      edges[1].setEndPoints(v1, v2);
      edges[2].setEndPoints(v2, v0);
    }

    void setAdjFace(IncFace adjFace, int i) {
      assert i >= 0 && i <= 2;
      edges[i].setAdjFace(adjFace);
    }

    void setNormalAndAngle(vec v, float a) {
      normal = v;
      angle = a;

      // Compute pAD
      pAD = P(c, r, normal);
    }

    @Override
    boolean isVisibleFromCircle(vec v, float f) {
      if (acos(dot(v, normal)) < angle + f) return true;
      return false;
    }

    @Override
    IncEdge[] getEdges() {
      return edges;
    }

    @Override
    IncEdge[] getComplementalEdges(IncEdge e) {
      IncVertex va = e.va;
      int idx = 0;
      for(; idx < 3; ++idx) {
        if (vertices[idx] == va) break;
      }
      IncEdge[] compEdges = new IncEdge[2];
      compEdges[0] = edges[idx];
      compEdges[1] = edges[(idx + 1) % 3];
      return compEdges;
    }

    Triangle toTriangle(HashMap<pt, Integer> pids) {
      int a = pids.get(vertices[0].position);
      int b = pids.get(vertices[1].position);
      int c = pids.get(vertices[2].position);
      return new Triangle(a, b, c);
    }

    @Override
    void showFace() {
      assert normal != null;
      if (showTriangleStrokes) stroke(0);
      else noStroke();
      showTriangle(vertices[0].position, vertices[1].position, vertices[2].position);
      if (showTriangleStrokes) noStroke();  // restore state
    }

    void showCircumcircle() {
      showCircumcircleOfTriangle(vertices[0].position, vertices[1].position,
                                 vertices[2].position, null, normal, null);
    }

    void showCone(pt apex, float height) {
      pt p = P(apex, height, normal);
      float r = height * tan(angle);
      fan(p, V(p, apex), r);
    }

    void showCone(pt apex, vec n, float height) {
      pt p = P(apex, height, n);
      float r = abs(height * tan(angle));
      fan(p, V(p, apex), r);
    }

    void showADCircle() {
      showCircumcircle();
    }
  }


  /*
   * A corridor is a structure as shown below.
   *                 triangle 0
   *           D -----(edge 1)----- A
   *  circle L )                    ( circle R
   *           C -----(edge 0)----- B
   *                 triangle 1
   *
   * Are vertices A, B, C, D coplanar? I think so!
   */
  class IncCorridor extends IncFace {
    IncEdge[] edges;
    vec[] coneNormals;  // 2 normals, one for cone ABC, one for cone CDA
    float[] coneAngles;  // 2 half-angles, one for cone ABC, one for cone CDA

    float delta = TWO_PI / 20;  // the sampling density on an arc, default = TWO_PI / 20
    ArrayList<pt> samples;  // samples.size = 2 * (numSegments - 1), where numSegments = int(angle of bigger arc / dAngle)

    ArrayList<pt> psAD;  // points on the corresponding Apollonius diagram edge

    IncCorridor() {}
    IncCorridor(IncVertex v0, IncVertex v1, IncVertex v2, IncVertex v3) {
      vertices = new IncVertex[4];
      edges = new IncEdge[2];
      for (int i = 0; i < 2; ++i) edges[i] = new IncEdge();
      vertices[0] = v0;
      vertices[1] = v1;
      vertices[2] = v2;
      vertices[3] = v3;
      edges[0].setEndPoints(v1, v2);  // BC
      edges[1].setEndPoints(v3, v0);  // DA
      setupNormalsAndAngles();

      if (gNumPointsPerRing >= 3) {
         delta = TWO_PI / (gNumPointsPerRing + 1);
        //delta = 2 * (TWO_PI / gNumPointsPerRing);
      }
      else delta = TWO_PI / 20;
      // if (simpleCorridor) {
      //   delta = TWO_PI / 3;
      // } else {
      //   delta = TWO_PI / gNumPointsPerRing;
      // }

      // generateSamples();
      // if (showApolloniusDiagram) generatePointsAD();
    }

    /*
     * Compute the normal and circumradius of triangle ABC/CDA.
     */
    private void setupNormalsAndAngles() {
      coneNormals = new vec[2];
      coneAngles = new float[2];
      coneNormals[0] = normalOfTriangle(vertices[0].position, vertices[1].position,
                                       vertices[2].position);
      coneNormals[1] = normalOfTriangle(vertices[2].position, vertices[3].position,
                                       vertices[0].position);
      float r0 = circumradiusOfTriangle(vertices[0].position, vertices[1].position,
                                        vertices[2].position);
      vec v = V(vertices[0].position, c);  // the vector from A to the center of the sphere
      coneAngles[0] = asin(r0 / r);  // r, i.e. the radius of the sphere, can be accessed!
      if (dot(v, coneNormals[0]) > 0) {  // if the center is on the positive side
        coneAngles[0] = PI - coneAngles[0];  // the angle should be > PI/2
      }
      float r1 = circumradiusOfTriangle(vertices[2].position, vertices[3].position,
                                        vertices[0].position);
      coneAngles[1] = asin(r1 / r);
      if (dot(v, coneNormals[1]) > 0) {  // if the center is on the positive side
        coneAngles[1] = PI - coneAngles[1];  // the angle should be > PI/2
      }
      assert isZero(coneNormals[0].norm() - 1) && isZero(coneNormals[1].norm() - 1);
    }

    void setDelta(float delta) {
      this.delta = delta;
    }

    void setAdjFace(IncFace adjFace, int i) {
      assert i == 0 || i == 1;
      edges[i].setAdjFace(adjFace);
    }

    private void generateSamples(float delta) {
      int left = vertices[3].rid;  // left ring
      int right = vertices[0].rid;  // right ring
      pt c0 = centers[left], c1 = centers[right];
      float r0 = radii[left], r1 = radii[right];
      vec n0 = normals[left], n1 = normals[right];

      /* Compute angles for arc DC and arc BA w.r.t. corresponding ring orientations. */
      vec vd = V(c0, vertices[3].position);
      vec vc = V(c0, vertices[2].position);
      vec vb = V(c1, vertices[1].position);
      vec va = V(c1, vertices[0].position);

      float a0 = acos(clamp(dot(vd, vc) / (r0 * r0), -1.0, 1.0));  // [0, PI]
      if (dot(n0, N(vd, vc)) < 0) a0 = TWO_PI - a0;
      {
        // fill(yellow);
        // show(vertices[3].position, 5);
        // fill(cyan);
        // show(vertices[2].position, 5);
        // println("dot(vd, vc) = ", dot(vd, vc));
        // println("r0 * r0 = ", r0 * r0);
        // println("a0 = ", a0);
      }
      float a1 = acos(clamp(dot(vb, va) / (r1 * r1), -1.0, 1.0));
      if (dot(n1, N(vb, va)) < 0) a1 = TWO_PI - a1;
      {
        // fill(springGreen);
        // show(vertices[1].position, 5);
        // fill(firebrick);
        // show(vertices[0].position, 5);
        // println("dot(vb, va) = ", dot(vb, va));
        // println("r1 * r1 = ", r1 * r1);
        // println("a1 = ", a1);
      }

      if (a0 > a1) {
        float stAng = acos(clamp(dot(vd, xAxes[left]) / r0, -1.0, 1.0));  // starting theta for D w.r.t. the left ring
        if (dot(vd, yAxes[left]) < 0) stAng = TWO_PI - stAng;
        int n = int(a0 / delta);
        samples = fillExEdgeWithPoints(stAng, delta, n, left, right);
      } else {
        float stAng = acos(clamp(dot(vb, xAxes[right]) / r1, -1.0, 1.0));
        if (dot(vb, yAxes[right]) < 0) stAng = TWO_PI - stAng;
        int n = int(a1 / delta);
        samples = fillExEdgeWithPoints(stAng, delta, n, right, left);
        Collections.reverse(samples);  // reverse the list
      }
    }

    /*
     * If a circle can see an adjacent trangle of the corridor, then it can see
     * this corridor since an adjacent triangle resides the tangent plane of an
     * straight edge of the corridor. Intuitively, if triangle 0 is visible,
     * then edge DA needs to move down. If triangle 1 is visible, then edge BC
     * needs to move up. A corridor is called marginally visible if it is
     * adjacent to a visible triangle.
     *
     * The following visibility test is to determine if a circle is in the
     * "middle" of the corridor, providing that it can't see the two adjacent
     * triangles. A corridor is called centrally visible if it passes the
     * following test. Note that a corridor can be both marginally visible and
     * centrally visible.
     *
     * If a corridor is only marginally visible, then it doesn't affect the
     * boundary.
     */
    @Override
    boolean isVisibleFromCircle(vec v, float f) {
      if ((acos(dot(v, coneNormals[0])) < coneAngles[0] + f) || (acos(dot(v, coneNormals[1])) < coneAngles[1] + f)) {
        // println("corridor is centrally visible");
        return true;
      }
      return false;
    }

    @Override
    IncEdge[] getEdges() {
      return edges;
    }

    @Override
    IncEdge[] getComplementalEdges(IncEdge e) {
      IncVertex va = e.va;
      IncEdge[] compEdges = new IncEdge[1];
      if (va == vertices[0]) {  // return edge BC
        compEdges[0] = edges[0];
      } else if (va == vertices[2]) {
        compEdges[0] = edges[1];
      } else {
        this.showFace();
        fill(darkRed);
        e.showEdge();
        println("Wrong input edge in getComplementalEdges() for corridor face.");
      }
      return compEdges;
    }

    @Override
    void showFace() {
      show(delta);
      // show(TWO_PI / 40);
    }

    void show(float delta) {
      generateSamples(delta);
      if (showCorridorStrokes) stroke(0);
      else noStroke();

      // show correspondences
      beginShape(QUAD_STRIP);
      vertex(vertices[3].position);  // D
      vertex(vertices[0].position);  // A
      for (pt p : samples) vertex(p);  // samples
      vertex(vertices[2].position);  // C
      vertex(vertices[1].position);  // B
      endShape();

      if (showCorridorStrokes) noStroke();  // restore state
    }

    void showCircumcircles() {
      showCircumcircleOfTriangle(vertices[0].position, vertices[1].position,
                                 vertices[2].position, null, coneNormals[0], null);
      showCircumcircleOfTriangle(vertices[2].position, vertices[3].position,
                                 vertices[0].position, null, coneNormals[1], null);
    }

    void showADCircle() {
      int k = samples.size() / 2;
      if (k % 2 == 0) k -= 2;
      else k -= 1;
      pt pa = null, pb = null;
      if (k < 0) {
        pa = vertices[3].position;
        pb = vertices[0].position;
      } else {
        pa = samples.get(k);
        pb = samples.get(k+1);
      }

      /* Compute the normal of the supporting plane of AB. */
      int ridLeft = vertices[3].rid;
      pt cLeft = centers[ridLeft];
      float rLeft = radii[ridLeft];
      vec xLeft = xAxes[ridLeft];
      vec yLeft = yAxes[ridLeft];
      vec ca = V(cLeft, pa);
      float theta = acos(dot(ca, xLeft) / rLeft);
      if (dot(ca, yLeft) < 0) theta = TWO_PI - theta;
      vec t = V(-sin(theta), xLeft, cos(theta), yLeft);
      vec ab = U(pa, pb);
      vec nor = U(N(ab, t));

      /* Compute the center and radius of the circle defined by the supporting plane of AB. */
      float d = dot(V(c, pa), nor);  // d can be negative
      pt center = P(c, d, nor);
      float radius = sin(acos(d / r)) * r;

      showCircle(center, nor, radius);

      hint(DISABLE_DEPTH_TEST);
      beginShape(LINES);
      vertex(pa);
      vertex(pb);
      endShape();
      hint(ENABLE_DEPTH_TEST);
    }

    /* Show things related to debug, e.g. circumcircles of two adjacent triangles. */
    void showDebugInfo() {
      /* Show the circumcircles of the two adjacent triangles. */
      stroke(magenta);
      strokeWeight(6);
      for (int i = 0; i < 2; ++i) {
        IncFace f = edges[i].getAdjFace();
        if (f != null) {
          assert f instanceof IncTriangle;
          IncTriangle t = (IncTriangle)f;
          t.showCircumcircle();
        }
      }

      /* Show the circumcircles of triangle ABC and CDA. These two circles are actually the same. */
      stroke(cyan);
      strokeWeight(6);
      showCircumcircles();
      strokeWeight(1);

      // checkCoplanarity(null, null);
    }

    private void generateExtendedPoints(pt pa, pt pb, float d, ArrayList<pt> extPoints) {
      vec v = U(pa, pb);
      extPoints.add(P(pa, -d, v));
      extPoints.add(P(pb, d, v));
    }

    void checkCoplanarity(Float del, Boolean disableDepthTest) {
      if (del == null) {
        generateSamples(TWO_PI / 40);
      } else {
        generateSamples(del);
      }
      coplanarFourPoints(vertices[3].position, vertices[0].position, vertices[2].position, vertices[1].position);
      if (samples.size() > 0) {  // >= 2
        int n = samples.size();
        assert n >= 2;
        coplanarFourPoints(vertices[3].position, vertices[0].position, samples.get(0), samples.get(1));
        for (int i = 0; i < n - 2; i += 2) {
          coplanarFourPoints(samples.get(i), samples.get(i+1), samples.get(i+2), samples.get(i+3));
        }
        coplanarFourPoints(samples.get(n-2), samples.get(n-1), vertices[2].position, vertices[1].position);
      }

      /* Extend the line segment defined by each pair of corresponding points. */
      ArrayList<pt> extPoints = new ArrayList<pt>();
      float d = 1000;
      generateExtendedPoints(vertices[3].position, vertices[0].position, d, extPoints);
      for (int i = 0; i < samples.size(); i += 2) {
        generateExtendedPoints(samples.get(i), samples.get(i+1), d, extPoints);
      }
      generateExtendedPoints(vertices[2].position, vertices[1].position, d, extPoints);

      if (disableDepthTest == null || disableDepthTest == true) {
        hint(DISABLE_DEPTH_TEST);
      }
      stroke(0);
      strokeWeight(1);
      beginShape(LINES);
      for (int i = 0; i < extPoints.size(); i += 2) {
        vertex(extPoints.get(i));
        vertex(extPoints.get(i+1));
      }
      endShape();
      noStroke();
      if (disableDepthTest == null || disableDepthTest == true) {
        hint(ENABLE_DEPTH_TEST);
      }
    }

    private pt generateADPoint(pt pa, pt pb, pt cLeft, float rLeft, vec xLeft, vec yLeft) {
      /* Compute the tangent at pa w.r.t. the left circle. */
      vec ca = V(cLeft, pa);
      float theta = acos(dot(ca, xLeft) / rLeft);
      if (dot(ca, yLeft) < 0) theta = TWO_PI - theta;
      vec t = V(-sin(theta), xLeft, cos(theta), yLeft);
      /* Compute the normal of the supporting plane touching pa and pb. */
      vec ab = U(pa, pb);
      vec nor = U(N(ab, t));

      {  // visualize the normal at the middle of (pa, pb)
        // fill(chocolate);
        // pt mid = P(pa, pb);
        // arrow(mid, V(20, nor), 1);
      }

      return P(c, r, nor);
    }

    void generatePointsAD() {
      if (psAD == null) psAD = new ArrayList<pt>();
      else psAD.clear();

      int n = samples.size();
      int ridLeft = vertices[3].rid;
      pt cLeft = centers[ridLeft];
      float rLeft = radii[ridLeft];
      vec xLeft = xAxes[ridLeft];
      vec yLeft = yAxes[ridLeft];

      /* Generate a point on DA if there is no adjacent triangle of DA. */
      if (edges[1].getAdjFace() == null) {
        pt p = generateADPoint(vertices[3].position, vertices[0].position, cLeft,
                               rLeft, xLeft, yLeft);
        psAD.add(p);
      }

      for (int i = 0; i < n; i += 2) {
        pt pa = samples.get(i);
        pt pb = samples.get(i+1);
        pt p = generateADPoint(pa, pb, cLeft, rLeft, xLeft, yLeft);
        psAD.add(p);
      }

      /* Generate a point on CB if there is no adjacent triangle of CB. */
      if (edges[0].getAdjFace() == null) {
        pt p = generateADPoint(vertices[2].position, vertices[1].position, cLeft,
                               rLeft, xLeft, yLeft);
        psAD.add(p);
      }
    }

    void showPointsAD() {
      IncTriangle t0 = (IncTriangle)edges[1].getAdjFace();
      IncTriangle t1 = (IncTriangle)edges[0].getAdjFace();

      pt p0 = null, p1 = null;
      if (t0 != null) {
        p0 = t0.pAD;
        fill(cyan);
        showBall(p0, 3);
      }
      for (pt p : psAD) {
        if (p0 != null) {
          showPolyArc(p0, p, c, r);
        }
        fill(lime);
        showBall(p, 2);
        p0 = p;
      }
      if (t1 != null) {
        p1 = t1.pAD;
        showPolyArc(p0, p1, c, r);
        fill(cyan);
        showBall(p1, 3);
      }
    }

    /*
     * Generate a list of triangles which approximate the corridor. The list of
     * positions will be updated as sampled points will be inserted.
     */
    ArrayList<Triangle> toTriangles(ArrayList<pt> positions, HashMap<pt, Integer> pids) {
      generateSamples(delta);
      ArrayList<Triangle> triangles = new ArrayList<Triangle>();
      if (samples == null || samples.size() == 0) {
        int a = pids.get(vertices[0].position);
        int b = pids.get(vertices[1].position);
        int c = pids.get(vertices[2].position);
        int d = pids.get(vertices[3].position);
        triangles.add(new Triangle(a, b, c));
        triangles.add(new Triangle(c, d, a));
        return triangles;
      }

      int td = pids.get(vertices[3].position);  // ID of tmp vertex D
      int ta = pids.get(vertices[0].position);  // ID of tmp vertex A

      int k = positions.size();  // k is the next valid vertex ID
      for (int i = 0; i < samples.size(); i += 2) {
        pt tpc = samples.get(i);
        positions.add(tpc);  // position of tmp vertex C
        int tc = k;  // ID of tmp vertex C
        pids.put(tpc, tc);

        pt tpb = samples.get(i + 1);
        positions.add(tpb);  // position of tmp vertex B
        int tb = k + 1;  // ID of tmp vertex B
        pids.put(tpb, tb);

        triangles.add(new Triangle(ta, tb, tc));
        triangles.add(new Triangle(tc, td, ta));

        td = tc;
        ta = tb;
        k += 2;
      }

      {  // C, B, and the last two samples
        int tc = pids.get(vertices[2].position);
        int tb = pids.get(vertices[1].position);
        triangles.add(new Triangle(ta, tb, tc));
        triangles.add(new Triangle(tc, td, ta));
      }
      return triangles;
    }
  }

  class DebugIncCHInfo {
    boolean newView = false;

    // things in old view
    ArrayList<IncFace> visibleFaces = new ArrayList<IncFace>();
    ArrayList<IncEdge> boundary = null;
    ArrayList<IncEdge> tmpBoundary = null;
    ArrayList<IncCorridor> removedUnvisibleCorridors = new ArrayList<IncCorridor>();
    ArrayList<IncCorridor> oldCorridors = new ArrayList<IncCorridor>();

    // things in new view
    ArrayList<IncFace> remainingFaces = null;
    ArrayList<IncFace> newFaces = null;
  }

  /* Basic data members. */
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

  /* For triangulation guided by convex hull of contacts. */
  TriangleMesh refConvexHull = null;  // convex hull generated by contacts
  ArrayList<Triangle> triangles = null;  // triangles, each formed by 3 ring vertices
  ArrayList<Triangle> threeRingTriangles = null;  // a three-ring triangle is formed by vertices from 3 rings
  ArrayList<Triangle> twoRingTriangles = null;  // a two-ring triangle is formed by vertices from 2 rings
  ArrayList<Integer>[] splitLists = null;  // each ring has a a sorted list of vertices of adjacent 3RTs (vertex IDs are stored)

  /* For naive convex hull. */
  ArrayList<pt> exTriPoints = null;  // a supporting triangle is formed by 3 consecutive points
  ArrayList<Integer> exTriRIDs = null;  // ring IDs of points in exTriPoints
  private ArrayList<Integer>[] swingLists = null;  // each ring has a sorted list of vertices of adjacent supporting triangles (vertex IDs are stored)
  private ArrayList<Float>[] angleLists = null;  // each vertex on a circle has a corresponding angle value
  private LinkedHashMap<NaiveCorridor, ArrayList<pt>> exEdges;  // key: (pid a, b, c, d), value: a list of points sampled from corresponding arcs

  /* For incremental convex hull. */
  ArrayList<IncFace> faces;
  ArrayList<IncCorridor> incCorridors;
  ArrayList<IncTriangle> incTriangles;

  /* For snapping vertices to ring samples. */
  HashMap<IncVertex, Integer> vertexToSample = new HashMap<IncVertex, Integer>();

  /* For gap filling in hub triangulation. swingLists and angleLists are also used. */
  ArrayList<pt>[] borders = null;

  /* For debug. */
  Debug3RTInfo debug3RTInfo = new Debug3RTInfo();
  Debug2RTInfo debug2RTInfo = new Debug2RTInfo();
  DebugSTInfo debugSTInfo = new DebugSTInfo();
  DebugIncCHInfo debugIncCHInfo = new DebugIncCHInfo();

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
    sameRadius = false;  // default
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
      if (dot(V(c, ps[i+1]), normals[j]) > 0) {
        centers[j] = P(c, sqrt(r2 - radii[j] * radii[j]), normals[j]);
      } else {
        centers[j] = P(c, -sqrt(r2 - radii[j] * radii[j]), normals[j]);
      }
    }
  }

  RingSet(pt c, float r, Circle[] circles, int nc, int np) {
    this.c = c;
    this.r = r;
    this.nRings = nc;
    this.nPointsPerRing = np;
    sameRadius = false;  // default
    centers = new pt[nc];
    normals = new vec[nc];
    radii = new float[nc];
    contacts = new pt[nc];
    xAxes = new vec[nRings];
    yAxes = new vec[nRings];
    for (int i = 0; i < nc; ++i) {
      centers[i] = circles[i].c;
      normals[i] = circles[i].n;
      radii[i] = circles[i].r;
      contacts[i] = P(this.c, this.r, normals[i]);
      xAxes[i] = constructNormal(normals[i]);
      yAxes[i] = N(normals[i], xAxes[i]);
    }
  }

  void setNumPointsPerRing(int np) {
    nPointsPerRing = np;
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

  private pt getPointFromGlobalID(int id) {
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
   * Generate a three-ring triangle using approximated supporting plane method.
   * This method first finds the exact supporting plane/triangle using the three
   * circles, and then trys to find a triangle defined by sampled points that
   * can approximate the exact supporting plane/triangle.
   */
  private Triangle generateThreeRingTriangleApprox(int rid0, int rid1, int rid2) {
    pt[] cs = {centers[rid0], centers[rid1], centers[rid2]};
    float[] rs = {radii[rid0], radii[rid1], radii[rid2]};
    vec[] ns = {normals[rid0], normals[rid1], normals[rid2]};
    vec[] vis = {xAxes[rid0], xAxes[rid1], xAxes[rid2]};
    vec[] vjs = {yAxes[rid0], yAxes[rid1], yAxes[rid2]};
    pt[] ps = supPlaneThreeCirclesIter(cs[0], rs[0], ns[0], vis[0], vjs[0],
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

    println("cannot find stable three-ring triangle using Approximated Supporting Plane method");
    numBackup3RT++;
    //return generateThreeRingTriangleNaive(ps0, ps1, ps2);
    return null;
  }

  private void setupSplitLists() {
    splitLists = new ArrayList[nRings];
    for (int i = 0; i < nRings; ++i) {
      ArrayList<Integer> swingList = refConvexHull.swingLists.get(i);
      int nCorners = swingList.size();  // number of adjacent corners/triangles
      splitLists[i] = new ArrayList<Integer>();
      for (int j = 0; j < nCorners; ++j) {
        int cid = swingList.get(j);
        int tid = cid / 3;  // triangle ID w.r.t. ref convex hull/set of three-ring triangles
        int pid = threeRingTriangles.get(tid).get(cid % 3);  // point ID (global)
        splitLists[i].add(pid);
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
      ArrayList<Integer> splits = splitLists[i];
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
      pt pa0 = getPointFromGlobalID(ia0), pb0 = getPointFromGlobalID(ib0), pc0 = getPointFromGlobalID(ic0);
      cornerTriangles.add(pa0);
      cornerTriangles.add(pb0);
      cornerTriangles.add(pc0);
      int start = (first + 1) % ns;
      int end = (first + ns - 1) % ns;
      int count = 1;

      // if (debugFixPenetration) {
      //   pt tmp = getPointFromGlobalID(tri0.get(cid0 % 3));
      //   fill(red, 100);
      //   show(tmp, 5);
      // }

      while (count < ns) {  // search right, i.e. increment start pointer if possible
        int cid = swings.get(start);
        Triangle tri = threeRingTriangles.get(cid / 3);
        int ia = tri.get(cid % 3), ib = tri.get((cid + 1) % 3), ic = tri.get((cid + 2) % 3);
        pt pa = getPointFromGlobalID(ia), pb = getPointFromGlobalID(ib), pc = getPointFromGlobalID(ic);
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
      //   pt tmp = getPointFromGlobalID(threeRingTriangles.get(cid / 3).get(cid % 3));
      //   fill(green, 100);
      //   show(tmp, 5);
      // }

      while (count < ns) {  // search left, i.e. decrement end pointer if possible
        int cid = swings.get(end);
        Triangle tri = threeRingTriangles.get(cid / 3);
        int ia = tri.get(cid % 3), ib = tri.get((cid + 1) % 3), ic = tri.get((cid + 2) % 3);
        pt pa = getPointFromGlobalID(ia), pb = getPointFromGlobalID(ib), pc = getPointFromGlobalID(ic);
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
      //   pt tmp = getPointFromGlobalID(threeRingTriangles.get(cid / 3).get(cid % 3));
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

  /*
   * Fill the discretized corridor shown below.
   *        loA ---------------- loB
   *  ring rA )                  ( ring rB
   *        hiA ---------------- hiB
   *
   * TODO:
   * I may need to revisit this function. Bacause this function may create
   * overlapping triangles (with non-consistent normals).
   */
  private void fillCorridor(int loA, int hiA, int loB, int hiB, int rA, int rB, vec normal) {
    int i = loA, j = loB;
    pt pa = getPointFromGlobalID(loA), pb = getPointFromGlobalID(loB);  // pa and pb should be updated in the while loop
    int iBase = rA * nPointsPerRing, jBase = rB * nPointsPerRing;
    vec nPrev = normal;  // may be null
    while(i != hiA || j != hiB) {
      pt pc = null, pd = null;
      int iNext = -1, jNext = -1;
      if (i != hiA) {  // increase i
        iNext = iBase + (i + 1) % nPointsPerRing;
        pc = getPointFromGlobalID(iNext);
      }
      if (j != hiB) {  // decrease j
        jNext = jBase + (j + nPointsPerRing - 1) % nPointsPerRing;
        pd = getPointFromGlobalID(jNext);
      }
      boolean moveA = true;
      if (pc != null && pd != null) {
        /* TODO: the criterion used to select a good triangle may be changed. */
        vec n1 = U(N(pa, pb, pc));
        vec n2 = U(N(pa, pb, pd));
        if ((nPrev != null && dot(n2, nPrev) < dot(n1, nPrev)) ||
            (nPrev == null) && d(pa, pd) < d(pb, pc)) {
          twoRingTriangles.add(new Triangle(i, j, jNext));
          nPrev = n2;
          moveA = false;
        } else {
          twoRingTriangles.add(new Triangle(i, j, iNext));
          nPrev = n1;
        }
      } else if (pc != null) {
        twoRingTriangles.add(new Triangle(i, j, iNext));
      } else {  // pd != null
        twoRingTriangles.add(new Triangle(i, j, jNext));
        moveA = false;
      }
      if (moveA) {
        i = iNext;
        pa = pc;
      } else {
        j = jNext;
        pb = pd;
      }
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
      ArrayList<Integer> splits = splitLists[i];
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
        fillCorridor(loA, hiA, loB, hiB, i, k, null);
        curC = swingC;
        visited[i][k] = visited[k][i] = true;

        if (debug2RT && i == debug2RTInfo.numGlobalStep - 1 && j == debug2RTInfo.numLocalStep - 1) {
          debug2RTInfo.pa0 = getPointFromGlobalID(loA);
          debug2RTInfo.pa1 = getPointFromGlobalID(hiA);
          debug2RTInfo.pb0 = getPointFromGlobalID(loB);
          debug2RTInfo.pb1 = getPointFromGlobalID(hiB);
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
        generateTriangleMeshFast();  // deprecated
        break;
      default:
        generateTriangleMeshCH();
    }
  }



  /*
   * Generate two supporting planes, with each touching three rings (i, j, k).
   * An array of 6 points will be returned. The first 3 points correspond to
   * orientation (i, j, k) while the last 3 points correspond to orientation
   * (i, k, j). If oneSolution is true, then the first 3 points are correct while
   * the last 3 points may be wrong.
   * normalsEx and halfAnglesEx are output. They define the two supporting
   * triangles.
   */
  pt[] twoSupPlanesThreeCircles(int i, int j, int k, vec[] normalsEx,
                                float[] halfAnglesEx, boolean oneSolution) {
    float s1 = radii[i] / r, s2 = radii[j] / r, s3 = radii[k] / r;
    float a1 = asin(s1), a2 = asin(s2), a3 = asin(s3);  // TODO: may store these half-angles?
    float c1 = cos(a1), c2 = cos(a2), c3 = cos(a3);

    vec n23 = N(normals[j], normals[k]);
    vec n31 = N(normals[k], normals[i]);
    vec n12 = N(normals[i], normals[j]);
    float mix = d(normals[i], n23);
    float aa = 0, bb = 0, ab = 0, det = 0;
    float[] ts = new float[2];  // the 2 values of tan(alpha)
    vec[] ns = new vec[2];  // the 2 normals
    pt[] ps = new pt[2];  // the 2 centers
    if (notZero(mix)) {
      float invMix = 1.0 / mix;
      n23 = V(invMix, n23);
      n31 = V(invMix, n31);
      n12 = V(invMix, n12);
      vec va = V(c1, n23, c2, n31, c3, n12);
      vec vb = V(s1, n23, s2, n31, s3, n12);
      bb = dot(vb, vb);
      float bb1 = bb - 1;
      if (isZero(bb1)) {
        println("dot(b, b) = 1");
        return null;
      }
      aa = dot(va, va);
      ab = dot(va, vb);
      det = ab * ab - (aa - 1) * bb1;
      if (det < 0) {
        println("determinant < 0");
        return null;
      }
      float sqrtDet = sqrt(det);
      ts[0] = (ab + sqrtDet) / bb1;
      ts[1] = (ab - sqrtDet) / bb1;
      /* Compute the 2 centers and 2 normals. */
      for (int ii = 0; ii < 2; ++ii) {
        float cosa = 1.0 / sqrt(1 + ts[ii] * ts[ii]);  // cosa > 0
        ns[ii] = V(cosa, A(va, -ts[ii], vb));
        ps[ii] = P(c, r * cosa, ns[ii]);
      }
    } else {
      println("mix product is 0. But I fixed it.");
      vec va = V(c1, n23, c2, n31, c3, n12);
      vec vb = V(s1, n23, s2, n31, s3, n12);

      if (notZero(va.x)) ts[0] = va.x / vb.x;
      else if (notZero(va.y)) ts[0] = va.y / vb.y;
      else ts[0] = va.z / vb.z;
      ts[1] = ts[0];

      /* Compute the 2 centers and 2 normals. */
      float alpha = atan(ts[0]);
      float cosa = cos(alpha);
      float sina = sin(alpha);

      pt p1 = new pt();
      pt p2 = new pt();

      float d1 = cosa * c1 - sina * s1;
      float d2 = cosa * c2 - sina * s2;
      float d3 = cosa * c3 - sina * s3;
      if (intersectionTwoPlanes(normals[i].x, normals[i].y, normals[i].z, d1,
                                normals[j].x, normals[j].y, normals[j].z, d2,
                                p1, p2)) {
      } else if (intersectionTwoPlanes(normals[j].x, normals[j].y, normals[j].z, d2,
                                       normals[k].x, normals[k].y, normals[k].z, d3,
                                       p1, p2)) {
      } else if (intersectionTwoPlanes(normals[k].x, normals[k].y, normals[k].z, d3,
                                       normals[i].x, normals[i].y, normals[i].z, d1,
                                       p1, p2)) {
      } else {
        println("No intersection between any 2 planes of the 3 planes!");
        return null;
      }

      pt[] tmp = intersectionLineSphere(p1, V(p1, p2), P(), 1);
      if (tmp == null) {
        println("No intersection between line and unit sphere!");
        return null;
      }
      ns[0] = V(tmp[0].x, tmp[0].y, tmp[0].z);
      ns[1] = V(tmp[1].x, tmp[1].y, tmp[1].z);
      ps[0] = P(c, r * cosa, ns[0]);
      ps[1] = P(c, r * cosa, ns[1]);
    }

    /*
     * For each computed circle, make sure all input circles on its negative side.
     * Since all input circles are already on the same side, we just need to
     * check if one input circle is on its negative side.
     */
    for (int ii = 0; ii < 2; ++ii) {
      vec v = V(ps[ii], centers[i]);
      if (dot(v, ns[ii]) > 0) ns[ii].rev();
    }

    pt[] pointsTangency = new pt[6];
    /* Compute the first triangle (i, j, k). */
    vec n1 = U(N(normals[i], ns[0]));  // the tangent at the contact point
    vec d1 = N(n1, normals[i]);
    pointsTangency[0] = P(centers[i], radii[i], d1);
    vec n2 = U(N(normals[j], ns[0]));
    vec d2 = N(n2, normals[j]);
    pointsTangency[1] = P(centers[j], radii[j], d2);
    vec n3 = U(N(normals[k], ns[0]));
    vec d3 = N(n3, normals[k]);
    pointsTangency[2] = P(centers[k], radii[k], d3);
    /* Compute the second triangle (i, j, k). */
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
      if (oneSolution == false) {  // skip some computation when only one solution is needed
        pt p = pointsTangency[4];
        pointsTangency[4] = pointsTangency[5];
        pointsTangency[5] = p;
      }
    } else {  // flip the first triplet, swap the two triplets
      pt[] tmpPointsTangency = new pt[6];
      tmpPointsTangency[0] = pointsTangency[3];
      tmpPointsTangency[1] = pointsTangency[4];
      tmpPointsTangency[2] = pointsTangency[5];
      if (oneSolution == false) {
        tmpPointsTangency[3] = pointsTangency[0];
        tmpPointsTangency[4] = pointsTangency[2];
        tmpPointsTangency[5] = pointsTangency[1];
      }
      pointsTangency = tmpPointsTangency;
      {  // swap normals and tan(alpha)'s
        vec v = ns[0];
        ns[0] = ns[1];
        ns[1] = v;
        float f = ts[0];
        ts[0] = ts[1];
        ts[1] = f;
      }
    }

    if (normalsEx != null && halfAnglesEx != null) {
      assert normalsEx.length <= 2 && halfAnglesEx.length <= 2;
      normalsEx[0] = ns[0];
      halfAnglesEx[0] = atan(ts[0]);
      if (halfAnglesEx[0] < 0) halfAnglesEx[0] += PI;
      if (oneSolution == false) {
        normalsEx[1] = ns[1];
        halfAnglesEx[1] = atan(ts[1]);
        if (halfAnglesEx[1] < 0) halfAnglesEx[1] += PI;
      }
    }

    return pointsTangency;
  }

  /*
   * Generate an supporting plane touching 3 circles (i, j, k). 3 points that
   * define the plane will be returned. These 3 points are on the 3 circles,
   * respectively.
   * TODO: this function is deprecated.
   */
  pt[] oneSupPlaneThreeCircles(int i, int j, int k) {
    float s1 = radii[i] / r, s2 = radii[j] / r, s3 = radii[k] / r;
    float a1 = asin(s1), a2 = asin(s2), a3 = asin(s3);  // TODO: may store these half-angles?
    float c1 = cos(a1), c2 = cos(a2), c3 = cos(a3);
    float x1 = normals[i].x, y1 = normals[i].y, z1 = normals[i].z;
    float x2 = normals[j].x, y2 = normals[j].y, z2 = normals[j].z;
    float x3 = normals[k].x, y3 = normals[k].y, z3 = normals[k].z;

    float DD = -x3*y2*z1 + x2*y3*z1 + x3*y1*z2 - x1*y3*z2 - x2*y1*z3 + x1*y2*z3;

    if (isZero(DD)) {
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

    return ps;
  }

  private boolean isValidExTri(pt pa, pt pb, pt pc, int ia, int ib, int ic) {
    vec normal = N(pa, pb, pc);  // not unit normal

    /* Make sure all circle centers on the negative side of the current plane. */
    for (int i = 0; i < nRings; ++i) {
      if (i == ia || i == ib || i == ic) continue;
      if (dot(V(pa, centers[i]), normal) > 0) return false;
    }

    /* Make sure the cone defined by the circumcircle of the current triangle is empty. */
    vec va = V(pc, pa);
    vec vb = V(pc, pb);
    float radius = sqrt(n2(va) * n2(vb) * n2(M(va, vb)) / n2(normal)) / 2;  // radius of circumcircle
    float alpha = asin(radius / r);  // half angle in [0, PI/2] since x >= 0 and asin(x) in [-PI/2, PI/2]
    vec unitNormal = U(normal);

    /*
     * If the center of the sphere is on the positive side of the triangle,
     * then the half angle of the cone is in [PI/2, PI].
     */
    if (dot(V(pa, c), unitNormal) > 0) alpha = PI - alpha;

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

  void generateExTrisNaive() {
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
          pt[] tmp = twoSupPlanesThreeCircles(i, j, k, null, null, false);
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

  /*
   * The exact convex hull of two circles should be a corridor (and two disks,
   * of course). I split this corridor into two pieces.
   * TODO: why sometimes the two pieces are the same?
   */
  private void generateExactCHTwoCircles() {
    if (incCorridors == null) incCorridors = new ArrayList<IncCorridor>();
    else incCorridors.clear();
    if (faces == null) faces = new ArrayList<IncFace>();
    else faces.clear();

    int ridLeft = 0, ridRight = 1;
    // if (radii[0] < radii[1]) {
    //   ridLeft = 1;
    //   ridRight = 0;
    // }

    float rLeft = radii[ridLeft];
    float rRight = radii[ridRight];
    pt cLeft = centers[ridLeft];
    pt cRight = centers[ridRight];
    vec nLeft = normals[ridLeft];
    vec nRight = normals[ridRight];
    pt pd = P(cLeft, rLeft, xAxes[ridLeft]);

    pt pc = P(cLeft, -rLeft, xAxes[ridLeft]);
    // pt pc = P(cLeft, rLeft, yAxes[ridLeft]);

    /* Find the corresponding point of pd. */
    pt pa = null;
    {
      pt pdd = P(pd, rLeft, yAxes[ridLeft]);  // sampled point + tangent at that point
      pt[] candidates = contactsOnSupportingPlaneOfLineCircle(cRight, rRight, nRight, pd, pdd, xAxes[ridRight], yAxes[ridRight]);
      pa = candidates[0];
      if (dot(V(pd, cRight), N(pd, pa, pdd)) > 0) {
        pa = candidates[1];
      }
    }

    /* Find the corresponding point of pc. */
    pt pb = null;
    {
      pt pcc = P(pc, -rLeft, yAxes[ridLeft]);
      // pt pcc = P(pc, -rLeft, xAxes[ridLeft]);

      pt[] candidates = contactsOnSupportingPlaneOfLineCircle(cRight, rRight, nRight, pc, pcc, xAxes[ridRight], yAxes[ridRight]);
      pb = candidates[0];
      if (dot(V(pc, cRight), N(pc, pb, pcc)) > 0) {
        pb = candidates[1];
      }
    }

    assert pa != null && pb != null;
    IncVertex va = new IncVertex(pa, ridRight);
    IncVertex vb = new IncVertex(pb, ridRight);
    IncVertex vc = new IncVertex(pc, ridLeft);
    IncVertex vd = new IncVertex(pd, ridLeft);

    /* Create corridor (A, B, C, D). */
    IncCorridor cor0 = new IncCorridor(va, vb, vc, vd);

    /* Create corridor (D, C, B, A). */
    IncCorridor cor1 = new IncCorridor(vd, vc, vb, va);

    incCorridors.add(cor0);
    incCorridors.add(cor1);

    faces.add(cor0);
    faces.add(cor1);
  }

  private ArrayList<IncFace> generateInitFaces(boolean[] success) {
    ArrayList<IncFace> fs = new ArrayList<IncFace>();
    vec[] nors = new vec[2];
    float[] angs = new float[2];
    pt[] ps = twoSupPlanesThreeCircles(0, 1, 2, nors, angs, false);
    // assert ps != null && ps.length == 6;
    if (ps == null) {
      if (success != null) success[0] = false;
      return null;
    }

    IncVertex v0 = new IncVertex(ps[0], 0);
    IncVertex v1 = new IncVertex(ps[1], 1);
    IncVertex v2 = new IncVertex(ps[2], 2);
    IncVertex v3 = new IncVertex(ps[3], 0);
    IncVertex v4 = new IncVertex(ps[4], 2);
    IncVertex v5 = new IncVertex(ps[5], 1);
    IncTriangle tri0 = new IncTriangle(v0, v1, v2);
    IncTriangle tri1 = new IncTriangle(v3, v4, v5);
    tri0.setNormalAndAngle(nors[0], angs[0]);
    tri1.setNormalAndAngle(nors[1], angs[1]);
    /* corridor w.r.t. edge (0, 1) and edge (5, 3) */
    IncCorridor cor0 = new IncCorridor(v5, v1, v0, v3);
    /* corridor w.r.t. edge (1, 2) and edge (4, 5) */
    IncCorridor cor1 = new IncCorridor(v4, v2, v1, v5);
    /* corridor w.r.t. edge (2, 0) and edge (3, 4) */
    IncCorridor cor2 = new IncCorridor(v3, v0, v2, v4);

    tri0.setAdjFace(cor0, 0);
    tri0.setAdjFace(cor1, 1);
    tri0.setAdjFace(cor2, 2);
    tri1.setAdjFace(cor2, 0);
    tri1.setAdjFace(cor1, 1);
    tri1.setAdjFace(cor0, 2);

    cor0.setAdjFace(tri0, 0);
    cor0.setAdjFace(tri1, 1);
    cor1.setAdjFace(tri0, 0);
    cor1.setAdjFace(tri1, 1);
    cor2.setAdjFace(tri0, 0);
    cor2.setAdjFace(tri1, 1);

    fs.add(tri0);
    fs.add(tri1);
    fs.add(cor0);
    fs.add(cor1);
    fs.add(cor2);
    return fs;
  }

  private ArrayList<IncEdge> exploreVisibleRegion(IncEdge e, vec normal, float angle, boolean[] success) {
    IncFace face = e.getAdjFace();
    if (face == null) {
      if (success != null) success[0] = false;
      return null;
    }
    ArrayList<IncEdge> edges = new ArrayList<IncEdge>();
    if (face.visible != -1) {  // if face is visited
      // println("face is visited, face.visible is ", face.visible);
      // assert face.visible == 0;  // then it must be unvisible
      if (face.visible != 0) {
        if (success != null) {
          success[0] = false;
          return edges;
        }
      }

      edges.add(e);
      return edges;
    }

    /* Current face is unvisited, check if it is visible. */
    if (!face.isVisibleFromCircle(normal, angle)) {  // if face is unvisible
      face.visible = 0;
      edges.add(e);
      return edges;
    }

    /* Current face is visible, try to expand boundary. */
    face.visible = 1;
    IncEdge[] compEdges = face.getComplementalEdges(e);
    assert compEdges.length == 1 || compEdges.length == 2;

    for (int i = 0; i < compEdges.length; ++i) {
      edges.addAll(exploreVisibleRegion(compEdges[i], normal, angle, success));
    }
    return edges;
  }

  private IncTriangle generateIncTri(int r0, int r1, int r2, boolean[] success) {
    vec[] nors = new vec[1];
    float[] angs = new float[1];
    pt[] ps = twoSupPlanesThreeCircles(r0, r1, r2, nors, angs, true);
    if (ps == null) {
      if (success != null) success[0] = false;
      return null;
    }
    IncVertex v0 = new IncVertex(ps[0], r0);
    IncVertex v1 = new IncVertex(ps[1], r1);
    IncVertex v2 = new IncVertex(ps[2], r2);
    IncTriangle tri = new IncTriangle(v0, v1, v2);
    tri.setNormalAndAngle(nors[0], angs[0]);
    return tri;
  }

  private IncCorridor generateIncCor(IncTriangle tri0, IncTriangle tri1, int rLeft, int rRight, boolean[] success) {
    int idx0 = 0;
    for (; idx0 < 3; ++idx0) {
      if (tri0.vertices[idx0].rid == rRight) break;
    }

    if (idx0 == 3) {
      if (success != null) success[0] = false;
      return null;
    }

    IncVertex va = tri0.vertices[idx0];
    IncVertex vd = tri0.vertices[(idx0 + 1) % 3];

    int idx1= 0;
    for(; idx1 < 3; ++idx1) {
      if (tri1.vertices[idx1].rid == rLeft) break;
    }
    if (idx1 == 3) {
      if (success != null) success[0] = false;
      return null;
    }

    IncVertex vc = tri1.vertices[idx1];
    IncVertex vb = tri1.vertices[(idx1 + 1) % 3];

    IncCorridor cor = new IncCorridor(va, vb, vc, vd);

    tri0.setAdjFace(cor, idx0);
    tri1.setAdjFace(cor, idx1);

    cor.setAdjFace(tri1, 0);
    cor.setAdjFace(tri0, 1);
    return cor;
  }

  private ArrayList<IncFace> generateNewFaces(ArrayList<IncEdge> boundary, int k, boolean[] success) {
    assert boundary != null && boundary.size() >= 1 && k >= 0 && k < nRings;

    /*
     * If the first edge and the last edge share a corridor that should be removed,
     * left shift the edge list by 1.
     */
    {
      IncFace firstAdjFace = boundary.get(0).getAdjFace();
      IncFace lastAdjFace = boundary.get(boundary.size() - 1).getAdjFace();
      if (firstAdjFace.visible >= 1 && firstAdjFace instanceof IncCorridor && firstAdjFace == lastAdjFace) {
        ArrayList<IncEdge> tmpBoundary = new ArrayList<IncEdge>();
        for (int i = 1; i < boundary.size(); ++i) {
          tmpBoundary.add(boundary.get(i));
        }
        tmpBoundary.add(boundary.get(0));
        boundary = tmpBoundary;
        // println("In generateNewFaces() after left shift, boundary =", boundary);
        // System.out.format("In generateNewFaces() after left shift, first edge = (%d, %d)\n",
        //                   boundary.get(0).va.rid, boundary.get(0).vb.rid);
      }
    }

    ArrayList<IncFace> newFaces = new ArrayList<IncFace>();
    IncTriangle prevTri = null;
    boolean missOppoTriPrev = false;
    for (IncEdge e : boundary) {
      IncTriangle oppoTri = null;
      IncFace adjFace = e.getAdjFace();
      if (adjFace instanceof IncCorridor) {
        IncEdge[] compEdges = adjFace.getComplementalEdges(e);
        IncFace adjAdjFace = compEdges[0].getAdjFace();
        assert adjAdjFace instanceof IncTriangle;
        oppoTri = (IncTriangle)adjAdjFace;
      } else {
        assert adjFace instanceof IncTriangle;
        oppoTri = (IncTriangle)adjFace;
      }

      int r0 = e.va.rid;
      int r1 = e.vb.rid;
      IncTriangle tri = generateIncTri(r0, r1, k, success);
      if (tri == null) {
        if (success != null) success[0] = false;
        return null;
      }
      newFaces.add(tri);


      if (oppoTri.visible < 1) {  // if the opposite triangle will remain
        assert oppoTri != null;
        IncCorridor cor = generateIncCor(tri, oppoTri, r1, r0, success);  // create a corridor between tri and oppoTri
        newFaces.add(cor);
      } else {
        /*
         * If the opposite triangle will be removed, then one more corridor must
         * be created between two new triangles w.r.t. edge (r1, r0). This corridor
         * is generated when we see the second new triangle.
         */
        if (missOppoTriPrev) {
          assert prevTri != null;
          IncCorridor cor = generateIncCor(tri, prevTri, r1, r0, success);  // create a corridor between tri and prevTri
          newFaces.add(cor);
          missOppoTriPrev = false;
        } else {
          missOppoTriPrev = true;
        }
      }

      if (prevTri != null) {
        IncCorridor cor = generateIncCor(tri, prevTri, r0, k, success);
        newFaces.add(cor);
      }
      prevTri = tri;
    }

    /* Create a corridor between the first triangle and the last triangle. */
    {
      assert newFaces.get(0) instanceof IncTriangle;
      IncTriangle firstTri = (IncTriangle)(newFaces.get(0));
      assert prevTri != null;
      IncCorridor cor = generateIncCor(firstTri, prevTri, boundary.get(0).va.rid, k, success);
      newFaces.add(cor);
    }

    return newFaces;
  }

  /*
   * Generate exact convex hull of circles using an incremental algorithm.
   */
  void generateExactCHIncremental(boolean[] success) {
    if (nRings < 2) {
      println("Less than 2 circles!");
      return;
    }

    /* Generate the convex hull of 2 circles. */
    if (nRings == 2) {
      generateExactCHTwoCircles();
      return;
    }

    /* Generate the initial convex hull using the first 3 circles. */
    faces = generateInitFaces(success);
    if (success != null && success[0] == false) return;

    ArrayList<IncEdge> boundary = new ArrayList<IncEdge>();
    ArrayList<IncFace> newFaces = null;
    ArrayList<IncFace> remainingFaces = new ArrayList<IncFace>();

    for (int i = 3; i < nRings; ++i) {
      // System.out.format("Insert the %d-th (0-base) circle\n", i);
      if (debugIncCH && debugIncCHIter < i) {
        break;
      }

      /* Initialize boundary using the first visible face. */
      boundary.clear();
      vec normal = normals[i];
      float angle = asin(radii[i] / r);  // [0, PI / 2]
      assert angle >= 0 && angle <= HALF_PI;
      for (IncFace f : faces) {
        if (f == null) return;
        if (f.isVisibleFromCircle(normal, angle)) {
          IncEdge[] edges = f.getEdges();
          assert edges.length == 2 || edges.length == 3;
          for (int j = 0; j < edges.length; ++j) {
            boundary.add(edges[j]);
          }
          f.visible = 1;  // visible
          break;
        } else {
          f.visible = 0;  // unvisible
        }
      }

      /* Expand the boundary. */
      int bsize = boundary.size();  // current number of edges
      if (bsize != 2 && bsize != 3) {
        System.out.format("Invalid boundary size %d, store this example.\n", bsize);
        gPoints.savePts("data/pts_unnamed");
        if (success != null) success[0] = false;
        return;
      }
      assert bsize == 2 || bsize == 3;

      if (debugIncCH && debugIncCHIter == i) {
        if (debugIncCHInfo.tmpBoundary == null) debugIncCHInfo.tmpBoundary = new ArrayList<IncEdge>();
        else debugIncCHInfo.tmpBoundary.clear();
      }

      ArrayList<IncEdge> tmpBoundary = new ArrayList<IncEdge>();
      for (IncEdge e : boundary) {
        ArrayList<IncEdge> tmpEdges = exploreVisibleRegion(e, normal, angle, success);
        if (tmpEdges == null) {
          if (success != null) success[0] = false;
          return;
        }
        if (success != null && success[0] == false) { //<>//
          return;
        }
        tmpBoundary.addAll(tmpEdges); //<>//
      }
      boundary = tmpBoundary;
      // System.out.format("boundary size = %d\n", boundary.size());

      /*
       * Label the corridors which are marginally visible but not centrally
       * visible.
       */
      for (IncFace f : faces) {
        if (f instanceof IncTriangle && f.visible == 1) {
          IncTriangle tri = (IncTriangle)f;
          for (IncEdge e : tri.edges) {
            IncFace adjFace = e.getAdjFace();
            // println("adjFace.visible =", adjFace.visible);
            assert adjFace.visible != -1;  // this corridor should be visited since it is adjacent to a visible triangle
            if (adjFace.visible < 1) {  // unvisited (not possible) or unvisible
              adjFace.visible = 2;  // should be removed
            }
          }
        }
      }

      /*
       * Construct new faces using boundary. Note that the boundary may get left
       * shifted by 1 inside of generateNewFaces(), but it stays unchanged
       * outside of generateNewFaces().
       */
      newFaces = generateNewFaces(boundary, i, success);
      // System.out.format("%d new faces generated\n", newFaces.size());
      // println("after generateNewFaces, boundary =", boundary);
      // System.out.format("after generateNewFaces, first edge = (%d, %d)\n",
      //                   boundary.get(0).va.rid, boundary.get(0).vb.rid);

      /*
       * Collect remaning faces, i.e. unvisible triangles and unvisible (neither
       * marginally nor centrally) corridors. A corridor will be removed if it
       * is marginally visible or centrally visible.
       */
      remainingFaces.clear();
      for (IncFace f : faces) {
        if (f.visible < 1) {
          f.visible = -1;  // old face becomes unvisited for next iteration
          remainingFaces.add(f);
        }
      }

      /* Store information related to debugging. */
      if (debugIncCH && debugIncCHIter == i) {
        debugIncCHInfo.boundary = boundary;  // store boundary
        debugIncCHInfo.visibleFaces.clear();  // store visible faces
        debugIncCHInfo.removedUnvisibleCorridors.clear();  // store just marginally visible corridors
        debugIncCHInfo.oldCorridors.clear();  // store all old corridors
        for (IncFace f : faces) {
          if (f.visible == 1) {
            debugIncCHInfo.visibleFaces.add(f);
          } else if (f.visible == 2) {
            debugIncCHInfo.removedUnvisibleCorridors.add((IncCorridor)f);
          }
          if (f instanceof IncCorridor) {
            debugIncCHInfo.oldCorridors.add((IncCorridor)f);
          }
        }
      }

      /* Collect remaining faces and new faces. All faces are unvisited. */
      faces.clear();
      faces.addAll(remainingFaces);
      if (newFaces == null) {
        if (success != null) success[0] = false;
        return;
      }
      faces.addAll(newFaces);

      /* Store information related to debugging. */
      if (debugIncCH && debugIncCHIter == i) {
        debugIncCHInfo.remainingFaces = remainingFaces;
        debugIncCHInfo.newFaces = newFaces;
      }
    }

    {  // collect triangles and corridors separately
      // println("#faces =", faces.size());
      incTriangles = new ArrayList<IncTriangle>();
      incCorridors = new ArrayList<IncCorridor>();
      for (IncFace f : faces) {
        if (f instanceof IncTriangle) incTriangles.add((IncTriangle)f);
        else if (f instanceof IncCorridor) incCorridors.add((IncCorridor)f);
      }
    }

    {  // print something
      // for (IncFace f : faces) {
      //   if (f instanceof IncCorridor) {
      //     IncCorridor cor = (IncCorridor)f;
      //     if (notZero(cor.coneNormals[0].x - cor.coneNormals[1].x) ||
      //         notZero(cor.coneNormals[0].y - cor.coneNormals[1].y) ||
      //         notZero(cor.coneNormals[0].z - cor.coneNormals[1].z)) {
      //         println("The two normals of a corridor are different.");
      //         System.out.format("normal 1 = (%f, %f, %f), normal 2 = (%f, %f, %f)\n",
      //                           cor.coneNormals[0].x, cor.coneNormals[0].y, cor.coneNormals[0].z,
      //                           cor.coneNormals[1].x, cor.coneNormals[1].y, cor.coneNormals[1].z);
      //         noStroke();
      //         fill(ivory);
      //         show(cor.vertices[0].position, 3);
      //         show(cor.vertices[1].position, 3);
      //         show(cor.vertices[2].position, 3);
      //         show(cor.vertices[3].position, 3);
      //     }
      //   }
      // }
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
      pt[] candidates = contactsOnSupportingPlaneOfLineCircle(cv, rv, nv, p0, p1, viv, vjv);
      pt candidate = candidates[0];
      if (dot(V(p0, cv), N(p0, candidate, p1)) > 0) {
        candidate = candidates[1];
      }
      edgePoints.add(p0);
      edgePoints.add(candidate);
    }
    return edgePoints;
  }

  private void fillExEdge(int u, int v, int ia, int ib, int ic, int id, float delta) {
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
      int n = int(au / delta);  // split arc into n pieces by inserting n-1 points
      ArrayList<pt> edgePoints = fillExEdgeWithPoints(angle, delta, n, u, v);
      exEdges.put(new NaiveCorridor(id, ic, ia, ib), edgePoints);
    } else {
      float angle = acos(dot(vb, xAxes[v]) / rv);
      if (dot(vb, yAxes[v]) < 0) angle = TWO_PI - angle;
      int n = int(av / delta);
      ArrayList<pt> edgePoints = fillExEdgeWithPoints(angle, delta, n, v, u);
      exEdges.put(new NaiveCorridor(ib, ia, ic, id), edgePoints);
    }
  }

  void generateExEdges(float delta) {
    int nExTriPoints = exTriPoints.size();
    boolean[][] visited = new boolean[nExTriPoints][nExTriPoints];
    exEdges = new LinkedHashMap<NaiveCorridor, ArrayList<pt>>();
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
        fillExEdge(u, v, pidA, pidB, pidC, pidD, delta);
        visited[pidD][pidC] = visited[pidC][pidD] = true;
        visited[pidA][pidB] = visited[pidB][pidA] = true;
      }
    }
  }

  private int findNearestSample(pt p, int rid) {
    pt center = centers[rid];
    vec v = V(center, p);
    float angle = acos(dot(v, xAxes[rid]) / radii[rid]);  // [0, PI]
    if (dot(v, yAxes[rid]) < 0) angle = TWO_PI - angle;
    float delta = TWO_PI / nPointsPerRing;
    int i = int(angle / delta);
    i = min(i, nPointsPerRing - 1);
    int j = (i + 1) % nPointsPerRing;

    pt a = points[rid][i];
    pt b = points[rid][j];
    int offset = rid * nPointsPerRing;
    if (d(a, p) > d(b, p)) return j + offset;
    return i + offset;
  }

  private TriangleMesh convertToTriMesh() {
    assert threeRingTriangles != null && threeRingTriangles.size() > 0;
    assert twoRingTriangles != null && twoRingTriangles.size() > 0;

    if (triangles != null) triangles.clear();
    else triangles = new ArrayList<Triangle>();

    ArrayList<pt> positions = get1DPointArrayList();
    triangles.addAll(threeRingTriangles);
    triangles.addAll(twoRingTriangles);

    return new TriangleMesh(positions, triangles);
  }

  private TriangleMesh generateMeshSnapping() {
    if (threeRingTriangles != null) threeRingTriangles.clear();
    else threeRingTriangles = new ArrayList<Triangle>();
    if (twoRingTriangles != null) twoRingTriangles.clear();
    else twoRingTriangles = new ArrayList<Triangle>();
    if (vertexToSample != null) vertexToSample.clear();
    else vertexToSample = new HashMap<IncVertex, Integer>();

    assert incTriangles != null && incCorridors != null && points != null;

    /* Create 3-ring triangles. */
    for (IncTriangle tri : incTriangles) {
      Integer[] sampleIDs = new Integer[3];
      for (int i = 0; i < 3; ++i) {
        IncVertex v = tri.vertices[i];
        if (vertexToSample.containsKey(v)) sampleIDs[i] = vertexToSample.get(v);
        else {  // find the nearest sample of v, update hash map
          sampleIDs[i] = findNearestSample(v.position, v.rid);
          vertexToSample.put(v, sampleIDs[i]);
        }
      }
      threeRingTriangles.add(new Triangle(sampleIDs[0], sampleIDs[1], sampleIDs[2]));
    }

    /* Create 2-ring triangles. */
    for (IncCorridor cor : incCorridors) {
      IncVertex va = cor.vertices[0];
      IncVertex vb = cor.vertices[1];
      IncVertex vc = cor.vertices[2];
      IncVertex vd = cor.vertices[3];
      int rLeft = vd.rid;
      int rRight = va.rid;
      int sa = vertexToSample.get(va);
      int sb = vertexToSample.get(vb);
      int sc = vertexToSample.get(vc);
      int sd = vertexToSample.get(vd);
      IncTriangle tri0 = (IncTriangle)cor.edges[1].getAdjFace();
      fillCorridor(sd, sc, sa, sb, rLeft, rRight, tri0.normal);
    }

    return convertToTriMesh();
  }

  TriangleMesh generateMeshFromExactCH(int option) {
    assert points != null;
    TriangleMesh tm = null;
    switch (option) {
      case 1:
        break;
      default:
        tm = generateMeshSnapping();
    }
    return tm;
  }

  private void updateSwingListGivenCorridor(IncCorridor cor, HashMap<pt, Integer> pids) {
    int ridLeft = cor.vertices[2].rid;
    int ridRight = cor.vertices[0].rid;
    int pidC = pids.get(cor.vertices[2].position);
    int pidA = pids.get(cor.vertices[0].position);

    ArrayList<pt> pointsDC = new ArrayList<pt>();
    ArrayList<pt> pointsBA = new ArrayList<pt>();
    for (int i = 0; i < cor.samples.size(); i += 2) {
      pointsDC.add(cor.samples.get(i));
      pointsBA.add(cor.samples.get(i+1));
    }
    Collections.reverse(pointsBA);

    ArrayList<Integer> pidsDC = new ArrayList<Integer>();
    ArrayList<Integer> pidsBA = new ArrayList<Integer>();
    for (int i = 0; i < pointsDC.size(); ++i) {
      pidsDC.add(pids.get(pointsDC.get(i)));
      pidsBA.add(pids.get(pointsBA.get(i)));
    }

    ArrayList<Integer> slLeft = swingLists[ridLeft];
    ArrayList<Integer> slRight = swingLists[ridRight];
    int locLeft = 0;
    int locRight = 0;
    for (; locLeft < slLeft.size(); ++locLeft) {
      if (slLeft.get(locLeft) == pidC) break;
    }
    for (; locRight < slRight.size(); ++locRight) {
      if (slRight.get(locRight) == pidA) break;
    }
    assert locLeft < slLeft.size() && locRight < slRight.size();

    slLeft.addAll(locLeft, pidsDC);
    slRight.addAll(locRight, pidsBA);
  }

  TriangleMesh generateConvexTriMesh() {
    /* Create a list of positions. */
    HashMap<pt, Integer> pids = new HashMap<pt, Integer>();  // position -> ID
    ArrayList<pt> posList = new ArrayList<pt>();  // ID -> position

    /* Create a list of triangles, each being a tuple of IDs. */
    ArrayList<Triangle> triList = new ArrayList<Triangle>();

    /* Setup the swing list and angle list for each circle. */
    swingLists = new ArrayList[nRings];
    angleLists = new ArrayList[nRings];
    for (int i = 0; i < nRings; ++i) {
      swingLists[i] = new ArrayList<Integer>();
      angleLists[i] = new ArrayList<Float>();
    }

    /* Store the vertices of all supporting triangles. */
    Integer id = 0;
    assert faces != null;
    for (IncFace f : faces) {
      if (f instanceof IncTriangle) {
        for (IncVertex v : f.vertices) {
          updateSwingList(v.rid, v.position, id);
          pids.put(v.position, id++);
          posList.add(v.position);
        }
      }
    }

    {  // verify the swing lists
      // fill(cyan, 100);
      // for (int i = 0; i < nRings; ++i) {
      //   ArrayList<Integer> sl = swingLists[i];
      //   for (int j = 0; j < sl.size() - 1; ++j) {
      //     pt a = posList.get(sl.get(j));
      //     pt b = posList.get(sl.get(j + 1));
      //     arrow(a, V(a, b), 3);
      //   }
      // }
    }

    /* Convert each supporting triangle to a tuple of IDs. */
    if (incTriangles != null) {
      for (IncTriangle tri : incTriangles) {
        triList.add(tri.toTriangle(pids));
      }
    }

    // println("Before converting corridors, size of pids =", pids.size());
    /* Convert each corridor to a list of tuples of IDs. */
    if (incCorridors != null) {
      for (IncCorridor cor : incCorridors) {
        triList.addAll(cor.toTriangles(posList, pids));
        updateSwingListGivenCorridor(cor, pids);
      }
    }
    // println("After converting corridors, size of pids =", pids.size());

    {  // verify the swing lists
      // fill(cyan);
      // for (int i = 0; i < nRings; ++i) {
      //   ArrayList<Integer> sl = swingLists[i];
      //   for (int j = 0; j < sl.size() - 1; ++j) {
      //     pt a = posList.get(sl.get(j));
      //     pt b = posList.get(sl.get(j + 1));
      //     arrow(a, V(a, b), 3);
      //   }
      // }
    }

    {  // setup borders
      borders = new ArrayList[nRings];
      for (int i = 0; i < nRings; ++i) {
        borders[i] = new ArrayList<pt>();
        ArrayList<Integer> sl = swingLists[i];
        for (int j = 0; j < sl.size(); ++j) {
          borders[i].add(posList.get(sl.get(j)));
        }
      }
    }

    if (triList.size() == 0) return null;
    TriangleMesh tm = new TriangleMesh(posList, triList);
    return tm;
  }

  void show() {
    fill(orange);
    for (int i = 0; i < nRings; ++i) {
      showBall(centers[i], 1);
    }
    fill(green);
    for (int i = 0; i < nRings; ++i) {
      arrow(centers[i], V(centers[i], points[i][0]), 2);
    }
    fill(blue);
    for (int i = 0; i < nRings; ++i) {
      for (int j = 0; j < nPointsPerRing; ++j) {
        showBall(points[i][j], 2);
      }
    }
    fill(cyan);
    for (int i = 0; i < nRings; ++i) {
      for (int j = 0; j < nPointsPerRing; ++j) {
        collar(points[i][j], V(points[i][j], points[i][(j + 1) % nPointsPerRing]), 1, 1);
      }
    }
    return;
  }

  void showCircles(Integer numCircles) {
    if (numCircles == null) numCircles = nRings;
    stroke(0);
    strokeWeight(3);
    for (int i = 0; i < numCircles; ++i) {
      showCircle(centers[i], normals[i], radii[i]);
    }
    strokeWeight(1);  // restore state
    noStroke();  // restore state
  }

  void showDisks(Integer numDisks) {
    if (numDisks == null) numDisks = nRings;
    fill(red, 255);
    for (int i = 0; i < numDisks; ++i) {
      disk(centers[i], xAxes[i], yAxes[i], radii[i]);
    }
  }

  /*
   * A polygon corresponds to a hole of the exact convex hull. Typically, each
   * polygon is a irregular polygon.
   */
  void showPolygons() {
    for (int i = 0; i < nRings; ++i) {
      ArrayList<pt> ps = borders[i];
      beginShape(TRIANGLE_FAN);
      vertex(centers[i]);
      for (int j = 0; j < ps.size(); ++j) {
        vertex(ps.get(j));
      }
      vertex(ps.get(0));
      endShape();
    }
  }

  void showExEdges() {
    if (exTriPoints == null || exEdges == null) return;
    for (NaiveCorridor edge : exEdges.keySet()) {
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
    beginShape(TRIANGLES);
    for (int i = 0; i < exTriPoints.size(); i += 3) {
      vertex(exTriPoints.get(i));
      vertex(exTriPoints.get(i+1));
      vertex(exTriPoints.get(i+2));
    }
    endShape();
  }

  void showCorridor(int index) {
    if (incCorridors.size() <= 0) return;
    int idx = index % incCorridors.size();
    incCorridors.get(idx).showDebugInfo();
  }

  void showIncTriangles() {
    if (incTriangles == null || incTriangles.size() == 0) return;
    for (IncTriangle tri : incTriangles) tri.showFace();
  }

  void showIncCorridors() {
    if (incCorridors == null || incCorridors.size() == 0) return;
    for (IncCorridor cor : incCorridors) cor.showFace();
  }

  void showThreeRingTriangles() {
    if (threeRingTriangles == null) return;
    stroke(0);
    beginShape(TRIANGLES);
    for (Triangle t : threeRingTriangles) {
      for (int i = 0; i < 3; ++i) {
        int pid = t.get(i);
        vertex(points[pid / nPointsPerRing][pid % nPointsPerRing]);
      }
    }
    endShape();
  }

  void showTwoRingTriangles() {
    if (twoRingTriangles == null) return;
    stroke(0);
    beginShape(TRIANGLES);
    for (Triangle t : twoRingTriangles) {
      for (int i = 0; i < 3; ++i) {
        int pid = t.get(i);
        vertex(points[pid / nPointsPerRing][pid % nPointsPerRing]);
      }
    }
    endShape();
  }

  /*
   * A beam is a cylinder. The i-th beam is the cylinder defined by the i-th
   * circle and its normal. length is the length of each cylinder.
   */
  void showBeams(float length) {
    if (points == null) return;
    stroke(0);
    for (int i = 0; i < nRings; ++i) {
      beginShape(QUAD_STRIP);
      for (int j = 0; j < nPointsPerRing; ++j) {
        pt a = points[i][j];
        pt b = P(a, length, normals[i]);
        vertex(a);
        vertex(b);
      }
      vertex(points[i][0]);
      vertex(P(points[i][0], length, normals[i]));
      endShape();
    }
  }

  /*
   * Show the cones, each defined by a circle and the center of the sphere.
   */
  void showCones(float extendDist) {
    for (int i = 0; i < nRings; ++i) {
      pt p = P(centers[i], extendDist, normals[i]);
      float r_tmp = -1.0;
      if (dot(V(c, centers[i]), V(c, contacts[i])) > 0) {
        r_tmp = radii[i] + extendDist * tan(asin(radii[i] / r));
      } else {
        r_tmp = abs(radii[i] - extendDist * tan(asin(radii[i] / r)));
      }
      vec v = V(p, c);
      fan(p, v, r_tmp);
    }
  }

  /*
   * Show the elliptic cone defined by circle c0 and circle c1.
   */
  void showEllipticCone(int c0, int c1) {
    pt pc0 = centers[c0];
    pt pc1 = centers[c1];
    float r0 = radii[c0];
    float r1 = radii[c1];
    vec n0 = normals[c0];
    vec n1 = normals[c1];

    vec n = U(N(n0, n1));  // normal to plane (c, pc0, pc1)
    vec v0 = N(n, n0);  // unit vector

    pt a1 = P(pc0, r0, v0);
    pt b1 = P(pc0, -r0, v0);

    vec v1 = N(n, n1);  // unit vector
    pt a2 = P(pc1, r1, v1);
    pt b2 = P(pc1, -r1, v1);

    pt f = intersectionTwoLines(a1, b2, b1, a2);
    vec v = U(A(U(a1, b2), U(b1, a2)));

    {
      stroke(red);
      showLine(f, v);
      noStroke();
    }

    {
      fill(blue, 220);
      showObliqueCone(f, pc0, n0, r0);
      fill(green, 150);
      showObliqueCone(f, pc1, n1, r1);
    }

    {
      // fill(cyan);
      // showBall(a1, 3);
      // showBall(b1, 3);

      // fill(magenta);
      // showBall(a2, 3);
      // showBall(b2, 3);
    }
  }

  void showApolloniusDiagram() {
    assert incCorridors != null && incCorridors.size() != 0;

    for (IncCorridor cor : incCorridors) {
      cor.showPointsAD();
    }
  }

  /* Show debug information related to Apollonius diagram. */
  void showADDebugInfo() {
    int i = idxIncTri % incTriangles.size();
    int j = idxIncCor % incCorridors.size();
    strokeWeight(3);
    stroke(navy);
    incTriangles.get(i).showADCircle();
    stroke(springGreen);
    incCorridors.get(j).showADCircle();
    strokeWeight(1);
    noStroke();
  }

  void showDebug3RTInfo() {
    fill(red, 150);
    showBall(contacts[debug3RTInfo.idr0], 3);
    fill(green, 150);
    showBall(contacts[debug3RTInfo.idr1], 3);
    fill(blue, 150);
    showBall(contacts[debug3RTInfo.idr2], 3);
    fill(#BF6868, 200);  // light red
    showBall(points[debug3RTInfo.idr0][debug3RTInfo.idp0], 5);
    fill(#40935D, 200);  // light green
    showBall(points[debug3RTInfo.idr1][debug3RTInfo.idp1], 5);
    fill(#517EC9, 200);  // light blue
    showBall(points[debug3RTInfo.idr2][debug3RTInfo.idp2], 5);
  }

  void showDebug2RTInfo() {
    fill(red, 200);
    showBall(debug2RTInfo.pa0, 5);
    fill(yellow, 200);
    showBall(debug2RTInfo.pa1, 5);
    fill(green, 100);
    showBall(debug2RTInfo.pb0, 5);
    fill(blue, 100);
    showBall(debug2RTInfo.pb1, 5);
    fill(#8B7373, 200);
    showBall(centers[debug2RTInfo.numGlobalStep - 1], 5);  // center of current ring
  }

  void showDebugIncCHInfo() {
    if (nRings == 2) {

      incCorridors.get(0).showDebugInfo();
      incCorridors.get(1).showDebugInfo();

      {  // show the line connecting the two centers
        // vec v = U(centers[0], centers[1]);
        // float d = 1000;
        // pt p0 = P(centers[0], -d, v);
        // pt p1 = P(centers[1], d, v);
        // stroke(darkRed);
        // strokeWeight(2);
        // beginShape(LINES);
        // vertex(p0);
        // vertex(p1);
        // endShape();
      }
    }

    if (nRings < 4) return;

    // show remaining faces in both views
    for (IncFace f : debugIncCHInfo.remainingFaces) {
      if (f instanceof IncTriangle) fill(blue);
      else if (f instanceof IncCorridor) fill(green);
      f.showFace();
    }
    // println("# remaining faces =", debugIncCHInfo.remainingFaces.size());

    {  // the new circle
      fill(gold);
      noStroke();
      disk(centers[debugIncCHIter], xAxes[debugIncCHIter], yAxes[debugIncCHIter],
           radii[debugIncCHIter]);
      stroke(0);
      strokeWeight(3);
      showCircle(centers[debugIncCHIter], normals[debugIncCHIter], radii[debugIncCHIter]);
      strokeWeight(1);
      noStroke();
    }

    if (debugIncCHNewView) {
      // show new faces
      for (IncFace f : debugIncCHInfo.newFaces) {
        if (f instanceof IncTriangle) fill(cyan);  // navy
        else if (f instanceof IncCorridor) {
          if (f.vertices[0].rid == debugIncCHIter || f.vertices[2].rid == debugIncCHIter) {
            fill(lightGreen);  // spring green for corridor adjacent to new circle
          } else {
            fill(ivory);  // ivory for corridor not adjacent to new circle
          }
        }
        f.showFace();
      }
      // println("# new faces =", debugIncCHInfo.newFaces.size());
    } else {
      // show visible faces
      for (IncFace f : debugIncCHInfo.visibleFaces) {
        if (f instanceof IncTriangle) fill(gray);
        else if (f instanceof IncCorridor) {
          IncCorridor cor = (IncCorridor)f;
          IncTriangle t0 = (IncTriangle)cor.edges[0].getAdjFace();
          IncTriangle t1 = (IncTriangle)cor.edges[1].getAdjFace();
          if (t0.visible == 1 && t1.visible == 1) {
            fill(gray);  // fully visible
          } else {
            fill(ivory);  // partly visible
          }
        }
        f.showFace();
      }

      // println("# visible faces =", debugIncCHInfo.visibleFaces.size());

      // show boundary
      // fill(violet);
      // for (IncEdge e : debugIncCHInfo.boundary) e.showEdge();
      // println("boundary size =", debugIncCHInfo.boundary.size());

      // show marginally, but not centrally, visible corridors
      fill(ivory);
      for (IncCorridor cor : debugIncCHInfo.removedUnvisibleCorridors) {
        cor.showFace();
      }
      // println("# removed and unvisible corridors =",
      //         debugIncCHInfo.removedUnvisibleCorridors.size());

      // show a corridor
      if (debugIncCHCor) {
        int idx = idxIncCor % debugIncCHInfo.oldCorridors.size();
        debugIncCHInfo.oldCorridors.get(idx).showDebugInfo();
      }
    }
  }

  void save(String file) {
    println("saving ring set:", file);
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
    centers = new pt[nRings];
    for (int j = 0; j < nRings; ++j) {
      float[] contact = float(split(lines[i++], ","));
      contacts[j] = new pt(contact[0], contact[1], contact[2]);
      normals[j] = U(c, contacts[j]);
      radii[j] = float(lines[i++]);
      float[] initDir = float(split(lines[i++], ","));
      xAxes[j] = new vec(initDir[0], initDir[1], initDir[2]);
      yAxes[j] = N(normals[j], xAxes[j]);
      centers[j] = P(c, sqrt(r * r - radii[j] * radii[j]), normals[j]);  // assume each cone has half-angle <= PI/2
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

  pts toPointSet() {
    int nv = 2 * nRings;
    pt[] ps = new pt[nv];
    for (int i = 0; i < nRings; ++i) {
      ps[2 * i] = contacts[i];
      ps[2 * i + 1] = P(centers[i], radii[i], xAxes[i]);
    }
    pts pointSet = new pts();
    pointSet.declare();
    pointSet.nv = nv;
    for (int i = 0; i < nv; ++i) {
      pointSet.G[i].set(ps[i]);
    }
    return pointSet;
  }
}
