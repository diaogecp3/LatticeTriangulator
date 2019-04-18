/******************************************************************************
 * Lattice and steady lattice.
 *
 * The code for steady lattice is adapted from the version provided by Kelsey
 * Kurzeja at Georgia Tech.
 ******************************************************************************/

import java.util.List;

boolean debugLattice = false;
boolean showLattice = true;
boolean showSteadyLattice = false;


/*
 * General lattice, which is basiclly a graph with each vertex being a ball.
 */
class Lattice {
  int nBalls;
  int nBeams;
  Ball[] balls;
  int[][] beams;  // each beam is a pair of indices

  class DebugGapInfo {
    ArrayList<Integer> gapStarts = new ArrayList<Integer>();
    ArrayList<Integer> gapSizes = new ArrayList<Integer>();
    DebugGapInfo() {}
    void reset() {
      gapStarts.clear();
      gapSizes.clear();
    }
  }
  DebugGapInfo dGapInfo = new DebugGapInfo();

  Lattice() {}

  Lattice(Ball[] balls, int[][] beams) {
    this.balls = balls;
    this.beams = beams;
    nBalls = balls.length;
    nBeams = beams.length;
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
  }

  ArrayList<Integer>[] adjacencyLists() {
    if (beams == null) return null;
    ArrayList<Integer>[] adjLists = new ArrayList[nBalls];
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

  private ArrayList<pt> extractVertices(ArrayList<pt> positions, ArrayList<Integer> pIDs) {
    ArrayList<pt> ps = new ArrayList<pt>();
    for (Integer pID : pIDs) {
      ps.add(positions.get(pID));
    }
    return ps;
  }

  TriangleMesh triangulate(ArrayList<Integer>[] adjLists) {
    TriangleMesh triMesh = new TriangleMesh();
    HashMap<Edge, ArrayList<Integer>> beamToLoop = new HashMap<Edge, ArrayList<Integer>>();
    Queue<Integer> queue = new ArrayDeque<Integer>();
    queue.add(0);
    int k = 0;
    boolean[] visited = new boolean[nBalls];  // all false
    if (debugLattice) dGapInfo.reset();
    while (queue.size() > 0) {
      // if (k == 4) break;

      int ballID = queue.poll();
      // println("ballID =", ballID);
      if (visited[ballID]) continue;
      visited[ballID] = true;

      ArrayList<Integer> adjList = adjLists[ballID];
      // println("hub ID =", ballID);
      Hub hub = generateHub(ballID, adjList);
      BorderedTriangleMesh btm = hub.generateConvexHullMesh(gNumPointsPerRing);  // gNumPointsPerRing controls the resolution of each corridor

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
            // println("Fail!");
            stroke(magenta);
            strokeWeight(3);
            // fill(magenta, 100);
            // showOrientedLoop(ps0);
            // showOrientedLoop(ps1);
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
            dGapInfo.gapStarts.add(triMesh.nt);
            dGapInfo.gapSizes.add(tris.size());
          }
          if (showGapMesh) triMesh.augmentWithoutShift(tris);
          beamToLoop.remove(uv);
        } else {  // if the i-th beam is not initialized, initialize it
          beamToLoop.put(uv, border);
        }
      }
      k++;
    }

    return triMesh;
  }

  void show() {
    fill(red);
    for (int i = 0; i < nBalls; ++i) balls[i].show();
    stroke(green);
    strokeWeight(3);
    beginShape(LINES);
    for (int i = 0; i < nBeams; ++i) {
      int[] b = beams[i];
      vertex(balls[b[0]].c);
      vertex(balls[b[1]].c);
    }
    endShape();
    strokeWeight(1);
    noStroke();
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
  }
}

void showHubGapMesh(TriangleMesh mesh, ArrayList<Integer> gapStarts, ArrayList<Integer> gapSizes, color cHub, color cGap, boolean useStroke) {
  ArrayList<pt> positions = mesh.positions;
  ArrayList<Triangle> triangles = mesh.triangles;
  int i = 0;
  int j = 0;
  if (useStroke) stroke(0);
  beginShape(TRIANGLES);
  while (i < mesh.nt) {
    if (j < gapStarts.size() && i == gapStarts.get(j)) {
      fill(cGap);
      int end = gapStarts.get(j) + gapSizes.get(j);
      while (i < end) {
        vertex(positions.get(triangles.get(i).a));
        vertex(positions.get(triangles.get(i).b));
        vertex(positions.get(triangles.get(i).c));
        i++;
      }
      j++;
    } else {
      fill(cHub);
      vertex(positions.get(triangles.get(i).a));
      vertex(positions.get(triangles.get(i).b));
      vertex(positions.get(triangles.get(i).c));
      i++;
    }
  }
  endShape();
  if (useStroke) noStroke();
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

  private boolean restrictInBoundsU, restrictInBoundsV, restrictInBoundsW;

  private boolean hasLinearTimeBIQ = false;
  private boolean hasConstantTimeBIQ = false;

  public SteadyLattice() {
    G = new SwirlTransform();
    u = 0;  U = new SwirlTransform();
    v = 0;  V = new SwirlTransform();
    w = 0;  W = new SwirlTransform();
    unitJoints = new ArrayList<Ball>();
    nextJoints = new ArrayList<List<JointId>>();
    restrictInBoundsU = true;
    restrictInBoundsV = true;
    restrictInBoundsW = true;
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

  public int addJoint(Ball ball) {
    unitJoints.add(ball);
    nextJoints.add( new ArrayList<JointId>() );
    return unitJoints.size()-1;
  }
  public int addJoint(pt p, float r) { return addJoint(new Ball(p, r)); }

  public void addBeam(int a, int b, Idx3 offset) { nextJoints.get(a).add( new JointId(b, offset) ); }

  public SteadyLattice restrictInBounds(boolean inU, boolean inV, boolean inW) {
    restrictInBoundsU = inU;
    restrictInBoundsV = inV;
    restrictInBoundsW = inW;
    return this;
  }

  // public void updateBIQCostFlags() {
  //   hasLinearTimeBIQ = true;
  //   hasConstantTimeBIQ = true;
  //   BeamIterator beamIterator = new BeamIterator(this, BeamType.BOTH, BeamDirection.OUT);
  //   beamIterator.restrictBeamStartIndex(I(0, 0, 0), I(0, 0, 0));
  //   if (beamIterator.init()) do {
  //     if (hasConstantTimeBIQ && !transformsAreSeparableForBeam(W, V, U, beamIterator.currentStartBall(), beamIterator.currentEndBall()))
  //       hasConstantTimeBIQ = false;
  //     if (hasLinearTimeBIQ && !hasConstantTimeBIQ && !transformsAreSeparableForBeam(W, V, beamIterator.currentStartBall(), beamIterator.currentEndBall()))
  //       hasLinearTimeBIQ = false;
  //   } while (beamIterator.next());
  // }

  // public void inflate(int n) { inflate(n,n,n); }
  // public void inflate(int nx, int ny, int nz) {
  //   Frame3 G_frame = G.transformFrameNoCopy(new Frame3(), 1);
  //   Frame3 U_frame = U.transformFrameNoCopy(new Frame3(), -nx);
  //   Frame3 V_frame = V.transformFrameNoCopy(new Frame3(), -ny);
  //   Frame3 W_frame = W.transformFrameNoCopy(new Frame3(), -nz);

  //   Frame3 WinG = G_frame.toGlobalFrame(W_frame);
  //   G = decomposeSwirlTransform(new Frame3(), WinG, 1);
  //   G_frame = G.transformFrameNoCopy(new Frame3(), 1);

  //   Frame3 VinG = G_frame.toGlobalFrame(V_frame);
  //   G = decomposeSwirlTransform(new Frame3(), VinG, 1);
  //   G_frame = G.transformFrameNoCopy(new Frame3(), 1);

  //   Frame3 UinG = G_frame.toGlobalFrame(U_frame);
  //   G = decomposeSwirlTransform(new Frame3(), UinG, 1);

  //   u += 2*nx;
  //   v += 2*ny;
  //   w += 2*nz;
  // }

  public Idx3 repetitionCounts() { return I(u,v,w); }
  public int repetitionCountU() { return u; }
  public int repetitionCountV() { return v; }
  public int repetitionCountW() { return w; }
  public Idx3 dimensions() { return I(u+1,v+1,w+1); }
  public int numUnitJoints() { return unitJoints.size(); }
  public Ball getUnitJointBall(int i) { return unitJoints.get(i); }
  public List<JointId> getNextJoints(int i) { return nextJoints.get(i); }

  public boolean restrictInBoundsU() { return restrictInBoundsU; }
  public boolean restrictInBoundsV() { return restrictInBoundsV; }
  public boolean restrictInBoundsW() { return restrictInBoundsW; }

  public boolean hasLinearTimeBIQ() { return hasLinearTimeBIQ; }
  public boolean hasConstantTimeBIQ() { return hasConstantTimeBIQ; }

  public Ball getJointBall(int i, Idx3 index) { return ballAt(getUnitJointBall(i), index); }

  // public AABB3 getAABB() {
  //   AABB3 aabb = EmptyAABB3();
  //   JointIterator joints = new JointIterator(this);
  //   if (joints.init()) do {
  //     Ball ball = joints.currentBall();
  //     aabb.combine(ball.aabb());
  //   } while (joints.next());
  //   return aabb;
  // }

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

  private int globalBallIdx(int b, int i, int j, int k) {
    return b + (i + j * u + k * u * v) * unitJoints.size();
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

  public void show() {
    ArrayList<Ball> balls = new ArrayList<Ball>();
    ArrayList<Edge> beams = new ArrayList<Edge>();  // an edge is a pair of ball indices
    /*
    for (int k = 0; k < w; ++k) {
      for (int j = 0; j < v; ++j) {
        for (int i = 0; i < u; ++i) {
          Idx3 offset = new Idx3(i, j, k);
          for (Ball b : unitJoints) {
            balls.add(ballAt(b, offset));
          }
        }
      }
    }
    for (int k = 0; k < w; ++k) {
      for (int j = 0; j < v; ++j) {
        for (int i = 0; i < u; ++i) {  // current ball cloud ID is (i, j, k)
          int bIdOffset = (i + j * u + k * u * v) * unitJoints.size();
          for (int b = 0; b < nextJoints.size(); ++b) {  // look at each ball's adjacency list
            int bId0 = b + bIdOffset;
            List<JointId> adjList = nextJoints.get(b);
            for (int a = 0; a < adjList.size(); ++a) {  // look at each neighbor of b
              JointId jId = adjList.get(a);
              Idx3 nextCloudId = jId.offset.c().add(i, j, k);
              if (nextCloudId.i >= u || nextCloudId.j >= v || nextCloudId.k >= w) continue;
              int bId1 = jId.unitId + (nextCloudId.i + nextCloudId.j * u + nextCloudId.k * u * v) * unitJoints.size();
              beams.add(new Edge(bId0, bId1));
            }
          }
        }
      }
    }
    */
    MinMaxI[] ranges = new MinMaxI[3];
    // ranges[0] = new MinMaxI(max(gCubeCenter.i - gCubeHalfLength, 0), min(gCubeCenter.i + gCubeHalfLength + 1, u));
    // ranges[1] = new MinMaxI(max(gCubeCenter.j - gCubeHalfLength, 0), min(gCubeCenter.j + gCubeHalfLength + 1, v));
    // ranges[2] = new MinMaxI(max(gCubeCenter.k - gCubeHalfLength, 0), min(gCubeCenter.k + gCubeHalfLength + 1, w));
    ranges[0] = new MinMaxI(0, u);
    ranges[1] = new MinMaxI(0, v);
    ranges[2] = new MinMaxI(0, w);
    traverse(ranges, balls, beams);
    // println("# balls =", balls.size());
    // println("# beams =", beams.size());

    // fill(red);
    // for (Ball b : balls) b.show();

    stroke(green);
    for (Edge beam : beams) {
      showSegment(balls.get(beam.a).c, balls.get(beam.b).c);
    }
    noStroke();

    // ArrayList<TruncatedCone> cones = new ArrayList<TruncatedCone>();
    // for (Edge beam : beams) {
    //   Ball a = balls.get(beam.a);
    //   Ball b = balls.get(beam.b);
    //   cones.add(truncatedConeOfTwoBalls(a, b));
    // }
    // fill(blue);
    // for (TruncatedCone cone : cones) {
    //   cone.show(6, false);
    // }
  }

  public MinMaxI[] getCubeRange(Idx3 cubeCenter, int cubeHalfLength) {
    MinMaxI[] ranges = new MinMaxI[3];
    ranges[0] = new MinMaxI(max(cubeCenter.i - cubeHalfLength, 0), min(cubeCenter.i + cubeHalfLength + 1, u));
    ranges[1] = new MinMaxI(max(cubeCenter.j - cubeHalfLength, 0), min(cubeCenter.j + cubeHalfLength + 1, v));
    ranges[2] = new MinMaxI(max(cubeCenter.k - cubeHalfLength, 0), min(cubeCenter.k + cubeHalfLength + 1, w));
    return ranges;
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
