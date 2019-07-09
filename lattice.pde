/******************************************************************************
 * Lattice and steady lattice.
 *
 * The code for steady lattice is adapted from the version provided by Kelsey
 * Kurzeja at Georgia Tech.
 ******************************************************************************/

import java.util.List;

boolean debugLattice = false;

boolean gShowBeams = true;
boolean gShowJunctions = true;
boolean gShowCHoCCs = false;
boolean gShowLattice = true;
boolean gShowSteadyLattice = false;
boolean gNavigateSteadyLattice = false;

float gAvgVertexCount = 0.0;

/*
 * General lattice, which is basiclly a graph with each vertex being a ball.
 */
class Lattice {
  int nBalls;
  int nBeams;
  Ball[] balls;
  int[][] beams;  // each beam is a pair of indices

  boolean valid = false;

  /* Intermediate variables in lattice tessellation. */
  ArrayList<Integer>[] adjLists = null;
  Hub[] hubs = null;

  /* For debug. */
  class DebugBeamInfo {
    ArrayList<Integer> beamStarts = new ArrayList<Integer>();
    ArrayList<Integer> beamSizes = new ArrayList<Integer>();
    DebugBeamInfo() {}
    void reset() {
      beamStarts.clear();
      beamSizes.clear();
    }
  }
  DebugBeamInfo dBeamInfo = new DebugBeamInfo();
  IntList junctionVertexCounts = new IntList();  // junctionVertexCounts[i] is the vertex count of the i-th CHoCC
  ArrayList<Mesh> junctionMeshes = null;
  ArrayList<TriangleMesh> beamMeshes = null;

  Lattice() {}

  Lattice(Ball[] balls, int[][] beams) {
    this.balls = balls;
    this.beams = beams;
    nBalls = balls.length;
    nBeams = beams.length;
    init();
  }

  Lattice(ArrayList<Ball> ballList, ArrayList<Edge> beamList) {
    nBalls = ballList.size();
    nBeams = beamList.size();
    balls = new Ball[nBalls];
    beams = new int[nBeams][2];
    for (int i = 0; i < nBalls; ++i) balls[i] = ballList.get(i);
    for (int i = 0; i < nBeams; ++i) {
      beams[i][0] = beamList.get(i).a;
      beams[i][1] = beamList.get(i).b;
    }
    init();
  }

  ArrayList<Integer>[] generateAdjLists() {
    if (beams == null) return null;
    if (adjLists == null) adjLists = new ArrayList[nBalls];
    for (int i = 0; i < nBalls; ++i) adjLists[i] = new ArrayList<Integer>();
    for (int i = 0; i < nBeams; ++i) {
      int a = beams[i][0];
      int b = beams[i][1];
      adjLists[a].add(b);
      adjLists[b].add(a);
    }

    return adjLists;
  }

  Hub generateHub(int ballID, ArrayList<Integer> adjList) {
    Ball centerBall = balls[ballID];
    int nNeighbors = adjList.size();
    Ball[] outerBalls = new Ball[nNeighbors];
    for (int i = 0; i < nNeighbors; ++i) outerBalls[i] = balls[adjList.get(i)];
    return new Hub(centerBall, outerBalls, nNeighbors);
  }

  Hub[] generateHubs() {
    assert adjLists != null;

    if (hubs == null) hubs = new Hub[nBalls];
    for (int i = 0; i < nBalls; ++i) {
      hubs[i] = generateHub(i, adjLists[i]);
    }
    return hubs;
  }

  /*
   * Check validity of the lattice.
   */
  boolean isValid() {
    assert nBalls > 0 && nBeams > 0;
    assert adjLists != null && hubs != null;

    /* Check if each hub is valid. */
    for (int i = 0; i < nBalls; ++i) {
      if (hubs[i].valid == false) return false;
    }

    /* Check intersection of inflating spheres of two neighboring hubs. */
    boolean[] visited = new boolean[nBalls];  // all false
    for (int i = 0; i < nBalls; ++i) {
      Ball b0 = hubs[i].getBoundingSphere();
      ArrayList<Integer> adjList = adjLists[i];
      for (Integer neighbor : adjList) {
        int j = (int)neighbor;
        if (visited[j]) continue;
        Ball b1 = hubs[j].getBoundingSphere();
        if (b0.intersectBall(b1)) {
          warningMsg += "Two inflating spheres intersect; ";
          return false;
        }
      }
      visited[i] = true;
    }
    return true;
  }

  /*
   * Initialize adjacency lists, hubs, and check validity of the lattice.
   */
  void init() {
    generateAdjLists();
    generateHubs();
    valid = isValid();
  }

  private ArrayList<pt> extractVertices(ArrayList<pt> positions, ArrayList<Integer> pIDs) {
    ArrayList<pt> ps = new ArrayList<pt>();
    for (Integer pID : pIDs) {
      ps.add(positions.get(pID));
    }
    return ps;
  }

  /*
   * Compute a triangle mesh that approximates the surface of the lattice.
   */
  TriangleMesh triangulate() {
    assert adjLists != null && hubs != null;
    TriangleMesh triMesh = new TriangleMesh();

    /*
     * Need to store information about beams that are being explored.
     * key: the two ball IDs of a beam
     * value: the vertex IDs on the visited end of the beam
     */
    HashMap<Edge, ArrayList<Integer>> beamToLoop = new HashMap<Edge, ArrayList<Integer>>();

    Queue<Integer> queue = new ArrayDeque<Integer>();
    queue.add(0);
    int k = 0;  // for counting how many hubs are processed
    boolean[] visited = new boolean[nBalls];  // all false

    if (debugLattice) {
      dBeamInfo.reset();
      junctionVertexCounts.clear();
    }

    while (queue.size() > 0) {
      int ballID = queue.poll();
      if (visited[ballID]) continue;
      visited[ballID] = true;

      ArrayList<Integer> adjList = adjLists[ballID];
      Hub hub = hubs[ballID];
      BorderedTriangleMesh btm = hub.generateConvexHullMesh(gNumPointsPerRing);  // gNumPointsPerRing controls the resolution of each corridor

      if (debugLattice) {
        junctionVertexCounts.append(btm.triangleMesh.nv);
      }

      int offset = triMesh.nv;
      btm.shiftVertexIDs(offset);

      /* Augment current triangle mesh. */
      triMesh.augmentWithoutShift(btm.triangleMesh);

      ArrayList<Integer>[] borders = btm.borders;

      /* Push new balls into the queue and store the beams that need to be matched. */
      for (int i = 0; i < adjList.size(); ++i) {
        int bID = adjList.get(i);
        if (visited[bID] == false) {
          queue.add(bID);
        }

        int u = min(ballID, bID);
        int v = max(ballID, bID);
        Edge uv = new Edge(u, v);

        ArrayList<Integer> border = borders[i];
        if (beamToLoop.containsKey(uv)) {  // if the i-th beam is initialized, construct the whole beam and remove this beam from the hash map
          ArrayList<Integer> otherBorder = beamToLoop.get(uv);
          ArrayList<pt> ps0 = extractVertices(triMesh.positions, border);
          ArrayList<pt> ps1 = extractVertices(triMesh.positions, otherBorder);

          Collections.reverse(ps1);
          Collections.reverse(otherBorder);

          {  // remove duplicate vertices before feeding them to convex gap constructor
            removeDuplicates(ps0, border);
            removeDuplicates(ps1, otherBorder);
          }

          ConvexGap gap = new ConvexGap(ps0, ps1);
          ArrayList<Triangle> tris = gap.gapHullGlobal(border, otherBorder);

          if (tris == null) {
            println("Failure in filling a gap! Please check the saved gap and hub files.");
            stroke(magenta);
            strokeWeight(3);
            showLoop(ps0);
            showLoop(ps1);
            strokeWeight(1);
            noStroke();
            gGap = gap;
            gGap.save("data/gap_unnamed");
            hub.save("data/hub_unnamed");
            return triMesh;
          }
          if (debugLattice) {
            dBeamInfo.beamStarts.add(triMesh.nt);
            dBeamInfo.beamSizes.add(tris.size());
          }

          triMesh.augmentWithoutShift(tris);
          beamToLoop.remove(uv);
        } else {  // if the i-th beam is not initialized, initialize it
          beamToLoop.put(uv, border);
        }
      }
      k++;
    }

    return triMesh;
  }

  /*
   * Compute junction meshes and beam meshes that together approximate the
   * surface of the lattice.
   *
   * Parameters:
   * subdivisionTimes: subdivision levels
   * projType: projection type
   */
  void tessellate(int subdivisionTimes, ProjectType projType) {
    assert adjLists != null && hubs != null;
    if (junctionMeshes == null) junctionMeshes = new ArrayList<Mesh>();
    else junctionMeshes.clear();
    if (beamMeshes == null) beamMeshes = new ArrayList<TriangleMesh>();
    else beamMeshes.clear();

    /*
     * Need to store information about beams that are being explored.
     * key: the two ball IDs of a beam
     * value: the vertices on the visited end of the beam
     */
    HashMap<Edge, ArrayList<pt>> beamToLoop = new HashMap<Edge, ArrayList<pt>>();

    Queue<Integer> queue = new ArrayDeque<Integer>();
    queue.add(0);
    boolean[] visited = new boolean[nBalls];  // all false

    while (queue.size() > 0) {
      int ballID = queue.poll();
      if (visited[ballID]) continue;
      visited[ballID] = true;

      ArrayList<Integer> adjList = adjLists[ballID];
      Hub hub = hubs[ballID];
      if (hub.ringSet == null) hub.generateRingSet(gNumPointsPerRing);
      hub.ringSet.generateExactCHIncremental(null);

      /* Generate junction mesh. */
      BorderedTriQuadMesh tqm = hub.ringSet.generateConvexTriQuadMesh();

      /* Subdivide it and project some vertices onto corresponding circles. */
      for (int i = 0; i < subdivisionTimes; ++i) {
        tqm = tqm.subdivide(SubdivideTypeTriangle.LOOP, SubdivideTypeQuad.DIAMOND, gProjectOnCircleAfterSub);
      }

      /* Project onto the hub. */
      if (subdivisionTimes > 0 && projType != null) {
        // tqm.projectOnHub(hub, projType);
        TriangleMesh tm = tqm.toTriangleMesh();
        boolean containHubCenter = hub.ringSet.pointInCrudestConvexHull(hub.ball.c);
        tm.containHubCenter = containHubCenter;
        tm.projectOnHub(hub, projType);
        junctionMeshes.add(tm);
      } else {
        junctionMeshes.add(tqm);
      }

      /* Push new balls into the queue and store the beams that need to be matched. */
      ArrayList<pt>[] borders = tqm.borders;
      assert borders.length == adjList.size();
      for (int i = 0; i < adjList.size(); ++i) {
        int bID = adjList.get(i);
        if (visited[bID] == false) {
          queue.add(bID);
        }

        int u = min(ballID, bID);
        int v = max(ballID, bID);
        Edge uv = new Edge(u, v);

        ArrayList<pt> border = borders[i];
        removeDuplicates(border);

        if (beamToLoop.containsKey(uv)) {  // if the i-th beam is initialized, construct the whole beam and remove this beam from the hash map
          ArrayList<pt> otherBorder = beamToLoop.get(uv);
          Collections.reverse(otherBorder);

          ConvexGap gap = new ConvexGap(border, otherBorder);
          TriangleMesh tm = gap.toTriMesh();
          if (tm == null) {
            println("Failure in filling a gap! Please check the saved gap and hub files.");
            gap.save("data/gap_unnamed");
            hub.save("data/hub_unnamed");
            return;
          }
          // if (!tm.isConvex()) {
          //   println("The gap mesh isn't convex! Please check the saved gap and hub files.");
          //   gap.save("data/gap_unnamed");
          //   hub.save("data/hub_unnamed");
          //   return;
          // }

          /* Store the beam mesh. */
          beamMeshes.add(tm);

          beamToLoop.remove(uv);
        } else {  // if the i-th beam is not initialized, initialize it
          beamToLoop.put(uv, border);
        }
      }
    }
  }

  void showInflatingSpheres() {
    assert hubs != null;
    for (int i = 0; i < nBalls; ++i) {
      Ball b = hubs[i].getBoundingSphere();
      b.show();
    }
  }

  /*
   * Show all convex hulls of cospherical circles. Please make sure that all
   * sets of circles and their corresponding convex hulls have been computed.
   *
   * Parameters:
   * cTriangle: color for triangles
   * cCorridor: color for corridors
   * cBeam: color for beams
   */
  void showCHoCCs(color cTriangle, color cCorridor, color cBeam) {
    assert hubs != null;
    /* Show all junctions with each being a CHoCC. */
    for (int i = 0; i < nBalls; ++i) {
      Hub hub = hubs[i];
      RingSet rs = hub.ringSet;
      if (rs == null) continue;
      fill(cTriangle);
      rs.showIncTriangles();
      fill(cCorridor);
      rs.showIncCorridors();
    }

    /* Show all beams. */
    fill(cBeam);
    int numSamples = 12;
    for (int i = 0; i < nBalls; ++i) {
      ArrayList<Integer> adjList = adjLists[i];
      Circle[] circles = hubs[i].circles;
      for (int j = 0; j < adjList.size(); ++j) {
        int k = adjList.get(j);
        if (k < i) continue;
        Circle c0 = circles[j];
        int h = adjLists[k].indexOf((Integer)i);
        Circle c1 = (hubs[k].circles)[h];

        /* Construct a beam using circle c0 and c1. */
        TruncatedCone cone = new TruncatedCone(c0.c, c0.r, c1.c, c1.r);
        cone.show(numSamples, false);
      }
    }
  }

  void showJunctions(color cJunction, boolean showStroke) {
    if (junctionMeshes == null || junctionMeshes.size() == 0) return;
    for (Mesh junction : junctionMeshes) {
      junction.show(cJunction, showStroke);
    }
  }

  void showBeams(color cBeam, boolean showStroke) {
    if (beamMeshes == null || beamMeshes.size() == 0) return;
    for (TriangleMesh beam : beamMeshes) beam.show(cBeam, showStroke);
  }

  void show() {
    /* Show balls. */
    for (int i = 0; i < nBalls; ++i) balls[i].show();

    /* Show beams. */
    int numSamples = 20;
    for (int i = 0; i < nBeams; ++i) {
      int[] b = beams[i];
      TruncatedCone cone = truncatedConeOfTwoBalls(balls[b[0]], balls[b[1]]);
      cone.show(numSamples, false);
    }
  }

  void save(String file) {
    println("saving lattice:", file);
    String[] lines = new String[nBalls + nBeams + 2];
    int i = 0;

    lines[i++] = str(nBalls);
    for (int j = 0; j < nBalls; ++j) {
      lines[i++] = str(balls[j].c.x) + "," + str(balls[j].c.y) + "," +
                   str(balls[j].c.z) + "," + str(balls[j].r);
    }

    lines[i++] = str(nBeams);
    for (int j = 0; j < nBeams; ++j) {
      lines[i++] = str(beams[j][0]) + "," + str(beams[j][1]);
    }

    saveStrings(file, lines);
  }

  void load(String file) {
    println("loading lattice:", file);
    String[] lines = loadStrings(file);
    int i = 0;

    nBalls = int(lines[i++]);
    balls = new Ball[nBalls];
    for (int j = 0; j < nBalls; ++j) {
      float[] tmp = float(split(lines[i++], ","));
      balls[j] = new Ball(tmp[0], tmp[1], tmp[2], tmp[3]);
    }

    nBeams = int(lines[i++]);
    beams = new int[nBeams][2];
    for (int j = 0; j < nBeams; ++j) {
      int[] tmp = int(split(lines[i++], ","));
      beams[j][0] = tmp[0];
      beams[j][1] = tmp[1];
    }

    init();
  }
}

/*
 * Use two colors to show respectively junctions and beams of a lattice mesh.
 *
 * Parameters:
 * mesh: the triangle mesh that approximates a lattice, assuming that every
 *       junction (or beam) is approximated by a subarray of triangles
 * beamStarts: beamStarts[i] is the index of the first triangle of the i-th beam
 * beamSizes: beamSizes[i] is the number of triangles of the i-th beam
 * cJunction: color for each junction
 * cBeam: color for each beam
 * showJunctions: show junctions or not
 * showBeams: show beams or not
 * showStroke: show the stroke or not
 */
void showGraphMesh(TriangleMesh mesh, ArrayList<Integer> beamStarts,
                   ArrayList<Integer> beamSizes, color cJunction, color cBeam,
                   boolean showJunctions, boolean showBeams, boolean showStroke) {
  ArrayList<pt> positions = mesh.positions;
  ArrayList<Triangle> triangles = mesh.triangles;
  int i = 0;
  int j = 0;
  if (showStroke) stroke(0);
  beginShape(TRIANGLES);
  while (i < mesh.nt) {
    if (j < beamStarts.size() && i == beamStarts.get(j)) {
      if (showBeams) {  // show beam triangle mesh on demand
        fill(cBeam);
        int end = beamStarts.get(j) + beamSizes.get(j);
        while (i < end) {
          vertex(positions.get(triangles.get(i).a));
          vertex(positions.get(triangles.get(i).b));
          vertex(positions.get(triangles.get(i).c));
          i++;
        }
      } else {
        i = beamStarts.get(j) + beamSizes.get(j);
      }
      j++;
    } else {
      if (showJunctions) {  // show junction triangle mesh on demand
        fill(cJunction);
        vertex(positions.get(triangles.get(i).a));
        vertex(positions.get(triangles.get(i).b));
        vertex(positions.get(triangles.get(i).c));
      }
      i++;
    }
  }
  endShape();
  if (showStroke) noStroke();
}

/*
 * Steady lattice, which is defined by a base ball cluster, a set of beams, and similarities.
 */
class SteadyLattice {
  private SwirlTransform G;  // Base transform

  private int u, v, w;             // Repetition counts
  private SwirlTransform U, V, W;  // Steady transformations

  private List<Ball> unitJoints;  // a unit joint is a ball in the base cluster
  private List<List<JointId>> nextJoints;  // adjacency lists

  // private boolean restrictInBoundsU, restrictInBoundsV, restrictInBoundsW;

  private boolean hasLinearTimeBIQ = false;
  private boolean hasConstantTimeBIQ = false;

  public SteadyLattice() {
    G = new SwirlTransform();
    u = 0;  U = new SwirlTransform();
    v = 0;  V = new SwirlTransform();
    w = 0;  W = new SwirlTransform();
    unitJoints = new ArrayList<Ball>();
    nextJoints = new ArrayList<List<JointId>>();
    // restrictInBoundsU = true;
    // restrictInBoundsV = true;
    // restrictInBoundsW = true;
  }

  public SteadyLattice deepCopy() {
    SteadyLattice n = new SteadyLattice();
    n.G = G.deepCopy();
    n.u = u;  n.U = U.deepCopy();
    n.v = v;  n.V = V.deepCopy();
    n.w = w;  n.W = W.deepCopy();
    for (Ball joint : unitJoints) n.unitJoints.add(joint.deepCopy());
    for (List<JointId> nextList : nextJoints) {
      List<JointId> newNextList = new ArrayList<JointId>();
      n.nextJoints.add(newNextList);
      for (JointId id : nextList) newNextList.add(id.deepCopy());
    }
    return n;
  }

  public SteadyLattice setG(SwirlTransform pG) { G = pG; return this; }
  public SteadyLattice setU(int pu, SwirlTransform pU) { u=pu; U=pU; return this; }
  public SteadyLattice setV(int pv, SwirlTransform pV) { v=pv; V=pV; return this; }
  public SteadyLattice setW(int pw, SwirlTransform pW) { w=pw; W=pW; return this; }

  public SwirlTransform getU() { return U; }
  public SwirlTransform getV() { return V; }
  public SwirlTransform getW() { return W; }

  public int addJoint(Ball ball) {
    unitJoints.add(ball);
    nextJoints.add( new ArrayList<JointId>() );
    return unitJoints.size()-1;
  }
  public int addJoint(pt p, float r) { return addJoint(new Ball(p, r)); }

  public void addBeam(int a, int b, Idx3 offset) { nextJoints.get(a).add( new JointId(b, offset) ); }

  // public SteadyLattice restrictInBounds(boolean inU, boolean inV, boolean inW) {
  //   restrictInBoundsU = inU;
  //   restrictInBoundsV = inV;
  //   restrictInBoundsW = inW;
  //   return this;
  // }

  public Idx3 repetitionCounts() { return I(u,v,w); }
  public int repetitionCountU() { return u; }
  public int repetitionCountV() { return v; }
  public int repetitionCountW() { return w; }
  public Idx3 dimensions() { return I(u+1,v+1,w+1); }
  public int numUnitJoints() { return unitJoints.size(); }
  public Ball getUnitJointBall(int i) { return unitJoints.get(i); }
  public List<JointId> getNextJoints(int i) { return nextJoints.get(i); }

  // public boolean restrictInBoundsU() { return restrictInBoundsU; }
  // public boolean restrictInBoundsV() { return restrictInBoundsV; }
  // public boolean restrictInBoundsW() { return restrictInBoundsW; }

  // public boolean hasLinearTimeBIQ() { return hasLinearTimeBIQ; }
  // public boolean hasConstantTimeBIQ() { return hasConstantTimeBIQ; }

  public Ball getJointBall(int i, Idx3 index) { return ballAt(getUnitJointBall(i), index); }

  public pt pointAt(pt start, Idx3 index) {
    // Index is not restricted to be in lattice bounds.
    pt p = start.c();
    U.transformPointNoCopy(p, index.i);
    V.transformPointNoCopy(p, index.j);
    W.transformPointNoCopy(p, index.k);
    return G.transformPointNoCopy(p, 1);
  }

  public Ball ballAt(Ball start, Idx3 index) {
    // Index is not restricted to be in lattice bounds.
    Ball ball = start.deepCopy();
    U.transformBallNoCopy(ball, index.i);
    V.transformBallNoCopy(ball, index.j);
    W.transformBallNoCopy(ball, index.k);
    return G.transformBallNoCopy(ball, 1);
  }

  public Frame3 frameAt(Idx3 index) {
    // Index is not restricted to be in lattice bounds.
    Frame3 f = new Frame3();
    U.transformFrameNoCopy(f, index.i);
    V.transformFrameNoCopy(f, index.j);
    W.transformFrameNoCopy(f, index.k);
    return G.transformFrameNoCopy(f, 1);
  }

  public Ball inverseBallAt(Ball start, Idx3 index) {
    // Index is not restricted to be in lattice bounds.
    Ball ball = G.inverseBallNoCopy(start.deepCopy(), 1);
    W.inverseBallNoCopy(ball, index.k);
    V.inverseBallNoCopy(ball, index.j);
    U.inverseBallNoCopy(ball, index.i);
    return ball;
  }

  public pt pointAt(pt start, vec index) {
    // Index is not restricted to be in lattice bounds.
    pt p = start.c();
    U.transformPointNoCopy(p, index.x);
    V.transformPointNoCopy(p, index.y);
    W.transformPointNoCopy(p, index.z);
    return G.transformPointNoCopy(p, 1);
  }

  public Ball ballAt(Ball start, vec index) {
    // Index is not restricted to be in lattice bounds.
    Ball ball = start.deepCopy();
    U.transformBallNoCopy(ball, index.x);
    V.transformBallNoCopy(ball, index.y);
    W.transformBallNoCopy(ball, index.z);
    return G.transformBallNoCopy(ball, 1);
  }

  public Frame3 frameAt(vec index) {
    // Index is not restricted to be in lattice bounds.
    Frame3 f = new Frame3();
    U.transformFrameNoCopy(f, index.x);
    V.transformFrameNoCopy(f, index.y);
    W.transformFrameNoCopy(f, index.z);
    return G.transformFrameNoCopy(f, 1);
  }

  public Ball inverseBallAt(Ball start, vec index) {
    // Index is not restricted to be in lattice bounds.
    Ball ball = G.inverseBallNoCopy(start.deepCopy(), 1);
    W.inverseBallNoCopy(ball, index.z);
    V.inverseBallNoCopy(ball, index.y);
    U.inverseBallNoCopy(ball, index.x);
    return ball;
  }

  public pt unitJointCenterCentroid() {
    pt c = P(0,0,0);
    for (Ball joint : unitJoints)
      c.add(joint.c);
    return c.div(unitJoints.size());
  }

  public int minRepetitionCount() {
    return min(u, min(v, w));
  }

  public int maxRepetitionCount() {
    return max(u, max(v, w));
  }

  private int globalBallIdx(int b, int i, int j, int k) {
    return b + (i + j * u + k * u * v) * unitJoints.size();
  }

  public MinMaxI[] getCubeRange(Idx3 cubeCenter, int cubeHalfLength) {
    MinMaxI[] ranges = new MinMaxI[3];
    ranges[0] = new MinMaxI(max(cubeCenter.i - cubeHalfLength, 0), min(cubeCenter.i + cubeHalfLength + 1, u));
    ranges[1] = new MinMaxI(max(cubeCenter.j - cubeHalfLength, 0), min(cubeCenter.j + cubeHalfLength + 1, v));
    ranges[2] = new MinMaxI(max(cubeCenter.k - cubeHalfLength, 0), min(cubeCenter.k + cubeHalfLength + 1, w));
    return ranges;
  }

  public MinMaxI[] getFullRange() {
    MinMaxI[] ranges = new MinMaxI[3];
    ranges[0] = new MinMaxI(0, u);
    ranges[1] = new MinMaxI(0, v);
    ranges[2] = new MinMaxI(0, w);
    return ranges;
  }

  /*
   * Traverse the lattice within specific ranges. Store the balls and beams
   * within those ranges. Each range is of the form [a, b).
   */
  public void traverse(MinMaxI[] ranges, ArrayList<Ball> balls, ArrayList<Edge> beams) {
    int iMin = ranges[0].min, iMax = ranges[0].max;
    int jMin = ranges[1].min, jMax = ranges[1].max;
    int kMin = ranges[2].min, kMax = ranges[2].max;

    balls.clear();
    for (int k = kMin; k < kMax; ++k) {
      for (int j = jMin; j < jMax; ++j) {
        for (int i = iMin; i < iMax; ++i) {
          Idx3 offset = new Idx3(i, j, k);
          for (Ball b : unitJoints) {
            balls.add(ballAt(b, offset));
          }
        }
      }
    }

    beams.clear();
    int kRange = kMax - kMin;
    int jRange = jMax - jMin;
    int iRange = iMax - iMin;
    int nUnitJoints = unitJoints.size();

    // each beam should store local ball IDs
    for (int k = kMin; k < kMax; ++k) {
      int kLocal = k - kMin;
      for (int j = jMin; j < jMax; ++j) {
        int jLocal = j - jMin;
        for (int i = iMin; i < iMax; ++i) {
          int iLocal = i - iMin;
          for (int b = 0; b < nextJoints.size(); ++b) {
            // local ball ID for the first ball of the beam
            int bId0 = b + (iLocal + jLocal * iRange + kLocal * jRange * iRange) * nUnitJoints;
            List<JointId> adjList = nextJoints.get(b);
            for (int a = 0; a < adjList.size(); ++a) {
              JointId jId = adjList.get(a);
              Idx3 nextCloudId = jId.offset.c().add(i, j, k);
              if (nextCloudId.i >= iMax || nextCloudId.j >= jMax || nextCloudId.k >= kMax) continue;
              // local ball ID for the second ball of the beam
              int bId1 = jId.unitId + (nextCloudId.i - iMin + (nextCloudId.j - jMin) * iRange + (nextCloudId.k - kMin) * jRange * iRange) * unitJoints.size();
              beams.add(new Edge(bId0, bId1));
            }
          }
        }
      }
    }
  }

  public void show(MinMaxI[] ranges) {
    ArrayList<Ball> balls = new ArrayList<Ball>();
    ArrayList<Edge> beams = new ArrayList<Edge>();  // an edge is a pair of ball indices

    if (ranges == null) ranges = getFullRange();
    traverse(ranges, balls, beams);

    // fill(red);
    // for (Ball b : balls) b.show();
    stroke(green);
    for (Edge beam : beams) {
      showSegment(balls.get(beam.a).c, balls.get(beam.b).c);
    }
    noStroke();
  }
}

class JointId {
  private int unitId;  // ball ID
  private Idx3 offset;

  public JointId(int pUnitId) { unitId=pUnitId; offset=I(0,0,0); }
  public JointId(int pUnitId, Idx3 pOffset) { unitId=pUnitId; offset=pOffset; }

  public JointId deepCopy() { return new JointId(unitId, offset.c()); }

  public int unitId() { return unitId; }
  public Idx3 offset() { return offset.c(); }
}

/*
enum BeamType { LOCAL, GLOBAL, BOTH }
enum BeamDirection { OUT, IN, BOTH }

class UnitBeamIterator {
  private SteadyLattice lattice;

  private BeamType beamType;
  private BeamDirection beamDirection;

  private boolean restrictUnitId;
  private int restrictedUnitId;

  private int index, index2;
  private boolean searchingOut;

  public UnitBeamIterator(SteadyLattice pLattice, BeamType pBeamType, BeamDirection pBeamDirection) {
    lattice = pLattice;
    beamType = pBeamType;
    beamDirection = pBeamDirection;
    restrictUnitId = false;
    restrictedUnitId = -1;
    index = index2 = 0;
    searchingOut = false;
  }

  public UnitBeamIterator restrictUnitId(int id) {
    restrictUnitId = true;
    restrictedUnitId = id;
    return this;
  }

  private boolean updateIndices() {
    while ((!restrictUnitId || index == restrictedUnitId) && index < lattice.numUnitJoints()) {
      if (index2 >= lattice.getNextJoints(index).size()) {
        index++;
        index2 = 0;
      } else {
        if (beamType == BeamType.GLOBAL && lattice.getNextJoints(index).get(index2).offset().isZero()) index2++;      // Ignore local connections
        else if (beamType == BeamType.LOCAL && !lattice.getNextJoints(index).get(index2).offset().isZero()) index2++; // Ignore global connections
        else return true;
      }
    }
    return false;
  }

  public boolean init() {
    searchingOut = beamDirection == BeamDirection.OUT || beamDirection == BeamDirection.BOTH;
    index = (restrictUnitId)? restrictedUnitId : 0;
    index2 = 0;
    return updateIndices();
  }

  public boolean next() {
    if (beamDirection == BeamDirection.BOTH) {
      if (searchingOut) {
        searchingOut = false;
        return true;
      } else {
        searchingOut = true;
        index2++;
        return updateIndices();
      }
    } else {
      index2++;
      return updateIndices();
    }
  }

  public Idx3 currentOffset() {
    if (searchingOut) return lattice.getNextJoints(index).get(index2).offset().c();
    else return lattice.getNextJoints(index).get(index2).offset().c().mul(-1);
  }
  public int currentStartUnitId() { return index; }
  public int currentEndUnitId() { return lattice.getNextJoints(index).get(index2).unitId(); }
}

class BeamIterator {
  private SteadyLattice lattice;

  private Idx3 min, max, stride;

  private UnitBeamIterator unitBeamIterator;
  private Idx3Iterator idxIterator;

  private boolean restrictInBoundsU, restrictInBoundsV, restrictInBoundsW;  // If true, trim beams that have at least one invalid ball k-index.

  public BeamIterator(SteadyLattice pLattice, BeamType pBeamType, BeamDirection pBeamDirection) {
    lattice = pLattice;
    unitBeamIterator = new UnitBeamIterator(lattice, pBeamType, pBeamDirection);
    min = I(0,0,0);
    max = lattice.repetitionCounts();
    stride = I(1,1,1);
    idxIterator = new Idx3Iterator(min, max, stride);
    restrictInBoundsU = lattice.restrictInBoundsU();
    restrictInBoundsV = lattice.restrictInBoundsV();
    restrictInBoundsW = lattice.restrictInBoundsW();
  }

  public BeamIterator restrictUnitId(int id) {
    unitBeamIterator.restrictUnitId(id);
    return this;
  }

  public BeamIterator restrictBeamStartIndex(Idx3 pMin, Idx3 pMax) { return restrictBeamStartIndex(pMin, pMax, I(1,1,1)); }
  public BeamIterator restrictBeamStartIndex(Idx3 pMin, Idx3 pMax, Idx3 pStride) {
    min = pMin;
    max = pMax;
    stride = pStride;
    return this;
  }

  public BeamIterator restrictInBounds(boolean inU, boolean inV, boolean inW) {
    restrictInBoundsU = inU;
    restrictInBoundsV = inV;
    restrictInBoundsW = inW;
    return this;
  }

  private boolean resetIdxIterator() {
    Idx3 dims = lattice.repetitionCounts();
    Idx3 offset = unitBeamIterator.currentOffset();

    int minI = min.i, minJ = min.j, minK = min.k;
    int maxI = max.i, maxJ = max.j, maxK = max.k;
    if (restrictInBoundsU) {
      minI = max(minI, (offset.i >= 0)? 0 : -offset.i);
      maxI = min(maxI, (offset.i >= 0)? dims.i - offset.i : dims.i);
    }
    if (restrictInBoundsV) {
      minJ = max(minJ, (offset.j >= 0)? 0 : -offset.j);
      maxJ = min(maxJ, (offset.j >= 0)? dims.j - offset.j : dims.j);
    }
    if (restrictInBoundsW) {
      minK = max(minK, (offset.k >= 0)? 0 : -offset.k);
      maxK = min(maxK, (offset.k >= 0)? dims.k - offset.k : dims.k);
    }

    idxIterator.set(I(minI, minJ, minK), I(maxI, maxJ, maxK), stride);
    return idxIterator.init();
  }

  public boolean init() {
    if (unitBeamIterator.init()) do {
      if (resetIdxIterator())
        return true;
    } while (unitBeamIterator.next());
    return false;
  }

  public boolean next() {
    if (idxIterator.next()) return true;
    while (unitBeamIterator.next())
      if (resetIdxIterator())
        return true;
    return false;
  }

  public Idx3 currentOffset() { return unitBeamIterator.currentOffset(); }
  public int currentStartUnitId() { return unitBeamIterator.currentStartUnitId(); }
  public int currentEndUnitId() { return unitBeamIterator.currentEndUnitId(); }
  public Idx3 currentStartIndex() { return idxIterator.current(); }
  public Idx3 currentEndIndex() { return idxIterator.current().add(currentOffset()); }

  public Ball currentStartBall() { return lattice.getJointBall(currentStartUnitId(), currentStartIndex()); }
  public Ball currentEndBall() { return lattice.getJointBall(currentEndUnitId(), currentEndIndex()); }
  // public ConicalFrustum currentConicalFrustum() { return getConicalFrustumBetweenBalls(currentStartBall(), currentEndBall()); }

  public int currentTotalInstancesAlongW() {
    int dimsW = lattice.repetitionCounts().k, offsetW = unitBeamIterator.currentOffset().k;
    return min(max.k, (offsetW >= 0)? dimsW - offsetW : dimsW) -
           max(min.k, (offsetW >= 0)? 0 : -offsetW) +
           1;
  }
}

class JointIterator {
  private SteadyLattice lattice;

  private boolean restrictUnitId;
  private int restrictedUnitId;

  private int unitId;

  private Idx3Iterator idxIterator;

  public JointIterator(SteadyLattice pLattice) {
    lattice = pLattice;
    unitId = 0;
    restrictUnitId = false;
    restrictedUnitId = -1;
    idxIterator = new Idx3Iterator(I(0,0,0), lattice.repetitionCounts(), I(1,1,1));
  }

  public JointIterator restrictUnitId(int id) {
    restrictUnitId = true;
    restrictedUnitId = id;
    return this;
  }

  public JointIterator restrictIndex(Idx3 pMin, Idx3 pMax) { return restrictIndex(pMin, pMax, I(1,1,1)); }
  public JointIterator restrictIndex(Idx3 pMin, Idx3 pMax, Idx3 pStride) {
    idxIterator.set(pMin, pMax, pStride);
    return this;
  }

  public boolean init() {
    unitId = (restrictUnitId)? restrictedUnitId : 0;
    return lattice.numUnitJoints() > unitId && idxIterator.init();
  }

  public boolean next() {
    if (restrictUnitId) return idxIterator.next();
    unitId++;
    if (unitId < lattice.numUnitJoints()) return true;
    unitId = 0;
    return idxIterator.next();
  }

  public int currentUnitId() { return unitId; }
  public Idx3 currentIndex() { return idxIterator.current(); }
  public Ball currentBall() { return lattice.getJointBall(currentUnitId(), currentIndex()); }
}
*/