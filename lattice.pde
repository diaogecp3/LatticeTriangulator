/******************************************************************************
 * Lattice.
 ******************************************************************************/


boolean showLattice = false;

class Lattice {
  int nBalls;
  int nBeams;
  Ball[] balls;
  int[][] beams;

  Lattice() {}

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
    while (queue.size() > 0) {
      // if (k == 4) break;

      int ballID = queue.poll();
      if (visited[ballID]) continue;
      visited[ballID] = true;

      ArrayList<Integer> adjList = adjLists[ballID];
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
            // return triMesh;
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
