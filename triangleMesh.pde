/******************************************************************************
 * Triangle mesh processing.
 ******************************************************************************/

boolean projectOnSphere = true;
boolean projectOnHub = false;
int subdivisionTimes = 0;  // subdivision times

int projectMethod = 0;
int numProjectMethod = 3;

int nextCorner(int cid) {
    return cid - (cid % 3) + (cid + 1) % 3;
}

int prevCorner(int cid) {
    return cid - (cid % 3) + (cid + 2) % 3;
}

/*
 * TriangleMesh class.
 *
 * A triangle mesh is a set of triangles. Another data structures may used to
 * simplify processing, e.g. opposite corner table.
 */
class TriangleMesh {
  ArrayList<pt> positions;
  ArrayList<Triangle> triangles;
  int nv;
  int nt;

  ArrayList<Integer> oppositeTable;
  ArrayList<ArrayList<Integer>> swingLists;

  TriangleMesh() {
    positions = new ArrayList<pt>();
    triangles = new ArrayList<Triangle>();
    nv = nt = 0;
  }

  TriangleMesh(ArrayList<pt> positions) {
    this.positions = positions;
    nv = positions.size();
    triangles = new ArrayList<Triangle>();
    nt = 0;
  }

  TriangleMesh(pt[] positionArray) {
    nv = positionArray.length;
    positions = new ArrayList<pt>();
    for (int i = 0; i < nv; ++i) positions.add(positionArray[i]);
    triangles = new ArrayList<Triangle>();
    nt = 0;
  }

  TriangleMesh(ArrayList<pt> positions, ArrayList<Triangle> triangles) {
    this.positions = positions;
    this.triangles = triangles;
    nv = this.positions.size();
    nt = this.triangles.size();
  }

  TriangleMesh(pt[] positionArray, ArrayList<Triangle> triangles) {
    nv = positionArray.length;
    positions = new ArrayList<pt>();
    for (int i = 0; i < nv; ++i) positions.add(positionArray[i]);
    this.triangles = triangles;
    nt = triangles.size();
  }

  TriangleMesh(TriangleMesh triMesh) {
    nv = triMesh.nv;
    nt = triMesh.nt;

    positions = new ArrayList<pt>();
    triangles = new ArrayList<Triangle>();

    ArrayList<pt> ps = triMesh.positions;
    ArrayList<Triangle> ts = triMesh.triangles;
    for (int i = 0; i < nv; ++i) positions.add(ps.get(i));
    for (int i = 0; i < nt; ++i) triangles.add(ts.get(i));
  }

  void augmentWithoutShift(ArrayList<Triangle> newTriangles) {
    if (newTriangles == null) return;
    triangles.addAll(newTriangles);
    nt = triangles.size();
  }

  void augmentWithoutShift(ArrayList<pt> newPositions, ArrayList<Triangle> newTriangles) {
    positions.addAll(newPositions);
    triangles.addAll(newTriangles);
    nv = positions.size();
    nt = triangles.size();
  }

  void augmentWithoutShift(TriangleMesh triMesh) {
    augmentWithoutShift(triMesh.positions, triMesh.triangles);
  }

  /*
   * Augment the current triangle mesh with a new set of triangles, assuming
   * the new triangles are all formed by new vertices. "shift" means that the
   * vertex IDs stored in each new triangle will be shifted, i.e. from local
   * vertex ID to global vertex ID.
   */
  void augmentWithShift(ArrayList<pt> newPositions, ArrayList<Triangle> newTriangles) {
    positions.addAll(newPositions);
    for (Triangle tri : newTriangles) {
      triangles.add(new Triangle(tri.a + nv, tri.b + nv, tri.c + nv));
    }
    nv = positions.size();
    nt = triangles.size();
  }

  /*
   * Augment the current triangle mesh with a new triangle mesh, assuming these
   * two triangle meshes don't share common vertices. "shift" means that the
   * vertex IDs stored in each new triangle will be shifted, i.e. from local
   * vertex ID to global vertex ID.
   */
  void augmentWithShift(TriangleMesh triMesh) {
    augmentWithShift(triMesh.positions, triMesh.triangles);
  }

  /*
   * Merge two triangle meshes, which may have vertices in common.
   */
  void augment(TriangleMesh triMesh) {
    int k = nv;  // ID for new vertex
    HashMap<pt, Integer> vids = new HashMap<pt, Integer>();
    ArrayList<pt> otherVertices = triMesh.positions;
    ArrayList<Triangle> otherTriangles = triMesh.triangles;

    for (pt p : otherVertices) {
      boolean isShared = false;
      int j = 0;
      for (; j < nv; ++j) {
        pt q = positions.get(j);
        if (isZero(d(p, q))) {
          isShared = true;
          break;
        }
      }
      if (isShared) { // p is in current mesh
        vids.put(p, j);
      } else {
        vids.put(p, k++);
        positions.add(p);  // add a new vertex
      }
    }

    /* Add new triangles. */
    for (Triangle tri : otherTriangles) {
      int a = vids.get(otherVertices.get(tri.a));
      int b = vids.get(otherVertices.get(tri.b));
      int c = vids.get(otherVertices.get(tri.c));
      triangles.add(new Triangle(a, b, c));
    }

    nv = positions.size();
    nt = triangles.size();
  }

  void shiftVertexIDs(int offset) {
    for (Triangle t : triangles) {
      t.set(t.a + offset, t.b + offset, t.c + offset);
    }
  }

  void addTriangle(Triangle tri) {
    assert tri.a < nv && tri.b < nv && tri.c < nv;
    assert tri.a >= 0 && tri.b >= 0 && tri.c >= 0;
    triangles.add(tri);
    nt++;
  }

  private void setupOppositeTable() {
    if (oppositeTable == null) {
      oppositeTable = new ArrayList<Integer>();
    }
    assert nv > 0 && nt > 0;
    int nc = 3 * nt;

    /* Make sure opposite table is filled with -1 and has proper size. */
    for (int i = 0; i < oppositeTable.size(); ++i) oppositeTable.set(i, -1);
    for (int i = oppositeTable.size(); i < nc; ++i) oppositeTable.add(-1);

    ArrayList<HashMap<Integer, int[]>> maps = new ArrayList<HashMap<Integer, int[]>>();
    for (int i = 0; i < nv; ++i) {
      maps.add(new HashMap<Integer, int[]>());
    }
    for (int i = 0; i < nt; ++i) {  // for each triangle
      Triangle triangle = triangles.get(i);
      for (int j = 0; j < 3; ++j) {  // for each corner
        int ea = triangle.get((j + 1) % 3);
        int eb = triangle.get((j + 2) % 3);  // (j - 1) % 3 can produce negative value
        if (ea > eb) {  // make sure ea < eb
          int tmp = ea;
          ea = eb;
          eb = tmp;
        }
        if (maps.get(ea).get(eb) == null) {
          maps.get(ea).put(eb, new int[2]);
          maps.get(ea).get(eb)[0] = 3 * i + j;  // store first corner
        } else {
          maps.get(ea).get(eb)[1] = 3 * i + j;  // store second corner
          int c0 = maps.get(ea).get(eb)[0];
          int c1 = maps.get(ea).get(eb)[1];
          oppositeTable.set(c0, c1);
          oppositeTable.set(c1, c0);
        }
      }
    }
    return;
  }

  void setupSwingLists() {
    if (oppositeTable == null) setupOppositeTable();
    assert triangles != null && oppositeTable != null;
    swingLists = new ArrayList<ArrayList<Integer>>();
    for (int i = 0; i < nv; ++i) swingLists.add(new ArrayList<Integer>());
    int numCorners = nt * 3;
    boolean[] visited = new boolean[numCorners];  // all false
    for (int i = 0; i < numCorners; ++i) {  // i = corner ID
      if (visited[i]) continue;
      ArrayList<Integer> list = new ArrayList<Integer>();
      int c = i;
      do {
        //System.out.format("c = %d ", c);
        list.add(c);
        visited[c] = true;
        int opp = oppositeTable.get(nextCorner(c));
        //System.out.format("o(n(c)) = %d ", opp);
        if (opp == -1) break;  // c is the last corner in cw order
        c = nextCorner(opp);  // update c
        //System.out.format("\n");
      } while (c != i && !visited[c]);
      //System.out.format("\n");
      int vid = cornerIDToVertexID(i);
      list.addAll(swingLists.get(vid));
      swingLists.set(vid, list);
    }
  }

  void subdivide(int times) {
    setupOppositeTable();
    int triedTimes = 0;
    while (triedTimes < times) {
      int nc = 3 * nt;
      int vid = nv;
      int[] midpointIDs = new int[nc];
      /* Insert a midpoint on each edge excluding border edges. */
      for (int i = 0; i < nt; ++i) {
        int t = 3 * i;
        Triangle triangle = triangles.get(i);
        for (int j = 0; j < 3; ++j) {
          int c = t + j;  // corner ID
          int o = oppositeTable.get(c);  // ID of opposite corner of c
          if (o == -1) {  // no opposite corner
            midpointIDs[c] = -1;  // no need to subdivide the corresponding edge
            continue;
          }
          if (o < c) continue;  // skip visited corner
          int a = triangle.get((j + 1) % 3);
          int b = triangle.get((j + 2) % 3);
          pt pa = positions.get(a);
          pt pb = positions.get(b);
          pt mid = P(pa, pb);

          positions.add(mid);
          midpointIDs[c] = midpointIDs[o] = vid++;
        }
      }
      nv = vid;  // update the size of positions

      ArrayList<Triangle> newTriangles = new ArrayList<Triangle>();
      /* Construct new triangles. */
      for (int i = 0; i < nt; ++i) {
        int[] loop = new int[6];
        int k = -1;
        Triangle triangle = triangles.get(i);
        int t = 3 * i;
        for (int j = 0; j < 3; ++j) {
          int c = t + j;
          int jj = 2 * j;
          loop[jj] = triangle.get(j);
          loop[(jj + 3) % 6] = midpointIDs[c];
          if (midpointIDs[c] == -1) {  // this occurs at most once for a triangle
            k = jj;
          }
        }
        if (k == -1) {  // 4 triangles will be created
          newTriangles.add(new Triangle(loop[5], loop[0], loop[1]));
          newTriangles.add(new Triangle(loop[1], loop[2], loop[3]));
          newTriangles.add(new Triangle(loop[3], loop[4], loop[5]));
          newTriangles.add(new Triangle(loop[5], loop[1], loop[3]));
        } else {  // 3 triangles will be created
          newTriangles.add(new Triangle(loop[(k + 5) % 6], loop[k], loop[k + 1]));
          /* TODO: Pick the second and third triangles based on some criterion. */
          int a = loop[k + 1];
          int b = loop[(k + 2) % 6];
          int c = loop[(k + 4) % 6];
          int d = loop[(k + 5) % 6];
          pt pa = positions.get(a);
          pt pb = positions.get(b);
          pt pc = positions.get(c);
          pt pd = positions.get(d);
          if (n2(N(pa, pb, pc)) > n2(N(pb, pc, pd))) {
            newTriangles.add(new Triangle(a, b, c));
            newTriangles.add(new Triangle(c, d, a));
          } else {
            newTriangles.add(new Triangle(b, c, d));
            newTriangles.add(new Triangle(d, a, b));
          }
        }
      }  // end for loop
      triangles = newTriangles;  // update triangles
      nt = newTriangles.size();  // update number of triangles

      /* TODO: may update this table when creating new triangle mesh. */
      setupOppositeTable();  // update opposite table
      triedTimes++;
    }
  }

  void projectOnSphere(pt c, float r) {
    for (int i = 0; i < nv; ++i) {
      pt p = positions.get(i);
      if (isZero(d(p, c) - r)) continue;
      vec v = U(c, p);
      p.set(P(c, r, v));  // this may not be good if c isn't in the interior of the mesh
    }
  }

  private void projectOnHubRay(Hub hub) {
    for (int i = 0; i < nv; ++i) {
      pt o = positions.get(i);
      if (hub.distanceFrom(o) < 0.0001) continue;
      vec d = U(o, hub.ball.c);  // from current vertex to the center of the hub
      Float t = hub.closestIntersectionWithRay(o, d);
      if (t != null) {
        o.set(P(o, t, d));
      }
    }
  }

  private void projectOnHubLine(Hub hub) {
    for (int i = 0; i < nv; ++i) {
      pt o = positions.get(i);
      // if (hub.distanceFrom(o) < 0.0001) continue;
      vec d = U(hub.ball.c, o);
      Float t = hub.closestIntersectionWithLine(o, d);
      if (t != null) {
        o.set(P(o, t, d));
      }
    }
  }

  private void projectOnHubSphereTrace(Hub hub) {
    int maxIter = 32;
    for (int i = 0; i < nv; ++i) {
      pt o = positions.get(i);
      if (hub.distanceFrom(o) < 0.0001) continue;
      vec d = U(o, hub.ball.c);  // from current vertex to the center of the hub

      pt p = P(o);  // copy
      for (int j = 0; j < maxIter; ++j) {
        float dist = hub.blendedDistanceFrom(p);
        if (dist < 0.0001) break;
        if (Float.isInfinite(dist)) {
          println("dist =", dist, "p =", p, "j =", j, "i =", i);
          dist = 0.0;
          break;
        }
        p.add(dist, d);
      }

      o.set(p);
    }
  }

  void projectOnHub(Hub hub, int option) {
    switch (option) {
      case 1:
        projectOnHubLine(hub);
        break;
      case 2:
        projectOnHubSphereTrace(hub);
        break;
      default:
        projectOnHubRay(hub);
    }
  }

  /* Test if the triangle mesh is convex. */
  boolean isConvex() {
    return passConvexityTest(triangles, positions);
  }

  void show(color c, boolean useStroke) {
    fill(c, 255);
    if (useStroke) {
      stroke(0);
      strokeWeight(1);
    } else noStroke();
    beginShape(TRIANGLES);
    for (int i = 0; i < nt; ++i) {
      vertex(positions.get(triangles.get(i).a));
      vertex(positions.get(triangles.get(i).b));
      vertex(positions.get(triangles.get(i).c));
    }
    endShape();
    if (useStroke) noStroke();
  }

  void showVertices(color c, float r) {
    fill(c, 255);
    noStroke();
    for (int i = 0; i < nv; ++i) {
      showBall(positions.get(i), r);
    }
    return;
  }

  void showCornerPairs(color c, float w) {
    int nc = 3 * nt;
    stroke(c);
    strokeWeight(w);
    for (int i = 0; i < nc; ++i) {
      int j = oppositeTable.get(i);
      if (j == -1 || j < i) continue;
      int a = cornerIDToVertexID(i);
      int b = cornerIDToVertexID(j);
      pt pa = positions.get(a);
      pt pb = positions.get(b);
      line(pa.x, pa.y, pa.z, pb.x, pb.y, pb.z);
    }
    return;
  }

  private int cornerIDToVertexID(int cid) {
    return triangles.get(cid / 3).get(cid % 3);
  }

  void save(String file) {
    println("saving triangle mesh:", file);
    String[] lines = new String[2 + nv + nt];
    int i = 0;
    lines[i++] = str(nv);
    for (int j = 0; j < nv; ++j) {
      pt p = positions.get(j);
      lines[i++] = str(p.x) + "," + str(p.y) + "," + str(p.z);
    }
    lines[i++] = str(nt);
    for (int j = 0; j < nt; ++j) {
      Triangle tri = triangles.get(j);
      lines[i++] = str(tri.a) + "," + str(tri.b) + "," + str(tri.c);
    }
    saveStrings(file, lines);
    return;
  }

  void load(String file) {
    println("loading triangle mesh:", file);
    String[] lines = loadStrings(file);
    int i = 0;

    nv = int(lines[i++]);
    if (positions == null) positions = new ArrayList<pt>();
    else positions.clear();
    for (int j = 0; j < nv; ++j) {
      float[] pos = float(split(lines[i++], ","));
      positions.add(new pt(pos[0], pos[1], pos[2]));
    }

    nt = int(lines[i++]);
    if (triangles == null) triangles = new ArrayList<Triangle>();
    else triangles.clear();
    for (int j = 0; j < nt; ++j) {
      int[] vid = int(split(lines[i++], ","));
      triangles.add(new Triangle(vid[0], vid[1], vid[2]));
    }

    setupOppositeTable();
    return;
  }
}

/*
 * Triangle mesh with borders.
 */
class BorderedTriangleMesh {
  TriangleMesh triangleMesh;
  ArrayList<Integer>[] borders;

  BorderedTriangleMesh(TriangleMesh triangleMesh, ArrayList<Integer>[] borders) {
    this.triangleMesh = triangleMesh;
    this.borders = borders;
  }

  BorderedTriangleMesh(ArrayList<pt> positions, ArrayList<Triangle> triangles,
                       ArrayList<Integer>[] borders) {
    this.triangleMesh = new TriangleMesh(positions, triangles);
    this.borders = borders;
  }

  void shiftVertexIDs(int offset) {
    this.triangleMesh.shiftVertexIDs(offset);
    for (ArrayList<Integer> border : borders) {
      for (int i = 0; i < border.size(); ++i) {
        border.set(i, border.get(i) + offset);
      }
    }
  }

  boolean validBorders() {
    ArrayList<pt> positions = triangleMesh.positions;
    for (ArrayList<Integer> border : borders) {
      ArrayList<pt> points = new ArrayList<pt>();
      for (Integer i : border) {
        points.add(positions.get(i));
      }
      if (!isConvexLoop(points)) {
        {  // debug
          // fill(magenta, 100);
          // showOrientedLoop(points);
        }
        return false;
      }
    }
    return true;
  }

  void show(color cMesh, color cBorder, boolean useStroke) {
    triangleMesh.show(cMesh, useStroke);
    ArrayList<pt> positions = triangleMesh.positions;
    fill(cBorder);
    for (ArrayList<Integer> border : borders) {
      ArrayList<pt> loop = new ArrayList<pt>();
      for (int i = 0; i < border.size(); ++i) {
        loop.add(positions.get(border.get(i)));
      }
      showOrientedLoop(loop);
    }
  }
}
