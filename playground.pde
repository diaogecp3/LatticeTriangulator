/******************************************************************************
 * Playground.
 *
 * For simple tests.
 ******************************************************************************/


void testHashSetContains() {
  Triangle t0 = new Triangle(0, 1, 2);
  Triangle t1 = new Triangle(0, 1, 2);
  HashSet<Triangle> set = new HashSet<Triangle>();
  set.add(t0);
  if (set.contains(t0)) {
    println("set contains t0");
  }
  if (set.contains(t1)) {
    println("set contains t1");
  }
}

void testSwingLists(boolean randomInput, boolean saveInput) {
  pts points = new pts();
  points.declare();
  if (randomInput) {
    generatePointsOnSphere(points, gSphereCenter, gSphereRadius, 4);
  } else {
    points.loadPts("data/point_set/ps_test");
  }

  ArrayList<pt> positions = new ArrayList<pt>();
  for (int i = 0; i < points.nv; ++i) positions.add(points.G[i]);

  ArrayList<Triangle> triangles = generateConvexHull(points.G, points.nv);
  TriangleMesh triMesh = new TriangleMesh(positions, triangles);
  triMesh.setupSwingLists();

  for (int i = 0; i < triangles.size(); ++i) {
    System.out.format("triangle %d = (%d, %d, %d)\n", i, triangles.get(i).a,
                      triangles.get(i).b, triangles.get(i).c);
  }

  System.out.format("opposite table:\n");
  for (int i = 0; i < triMesh.oppositeTable.size(); ++i) {
    System.out.format("%d <-> %d ", i, triMesh.oppositeTable.get(i));
  }
  System.out.format("\n");

  for (int i = 0; i < triMesh.nv; ++i) {
    System.out.format("swing list of vertex %d: ", i);
    int m = triMesh.swingLists.get(i).size();
    for (int j = 0; j < m; ++j) {
      System.out.format("%d ", triMesh.swingLists.get(i).get(j));
    }
    System.out.format("\n");
  }

  if (saveInput) {
    points.savePts("data/point_set/ps_unnamed");
  }
}
