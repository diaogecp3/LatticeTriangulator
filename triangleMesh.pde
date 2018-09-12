/******************************************************************************
 * Triangle mesh processing.
 ******************************************************************************/


class TriangleMesh {
  ArrayList<pt> positions;
  ArrayList<Triangle> triangles;
  ArrayList<Integer> oppositeTable;
  int nv, nt;

  TriangleMesh(ArrayList<pt> positions, ArrayList<Triangle> triangles) {
    this.positions = positions;
    this.triangles = triangles;
    oppositeTable = new ArrayList<Integer>();
    nv = this.positions.size();
    nt = this.triangles.size();
    setupOppositeTable();
  }

  TriangleMesh(ArrayList<pt> positions, ArrayList<Triangle> triangles,
               ArrayList<Integer> oppositeTable) {
    this.positions = positions;
    this.triangles = triangles;
    this.oppositeTable = oppositeTable;
    nv = this.positions.size();
    nt = this.triangles.size();
    setupOppositeTable();
  }

  int getNumTriangles() {
    return nt;
  }

  void setupOppositeTable() {
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

  private int cornerIDToVertexID(int cid) {
    return triangles.get(cid / 3).get(cid % 3);
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

  void showTriangleMesh(color c, boolean useStroke) {
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
  }

  void showVertices(color c, float r) {
    fill(c, 255);
    noStroke();
    for (int i = 0; i < nv; ++i) {
      show(positions.get(i), r);
    }
    return;
  }

  void subdivide(int times, pt center, float r) {
    int triedTimes = 0;
    while (triedTimes < times) {
      int nc = 3 * nt;
      int vid = nv;
      int[] midpointIDs = new int[nc];
      /* Insert midpoints. */
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
          mid = P(center, r, U(center, mid));  // lift to surface of sphere
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
}