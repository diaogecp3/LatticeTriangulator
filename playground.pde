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
    generatePointsOnSphere(points, centerOfSphere, radiusOfSphere, 4);
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

void testIntersectionTwoDisks() {
  pt ca = new pt(0.5, 0, 0);
  vec va = new vec(0, 0, 1);
  float ra = 1.0;

  pt cb = new pt(0.2, 0, 1.2);
  vec vb = new vec(1, 0, 0);
  float rb = 1.0;

  pt p0 = new pt(), p1 = new pt();
  intersectionTwoPlanes(ca, va, cb, vb, p0, p1);
  System.out.format("p0 = (%f, %f, %f), p1 = (%f, %f, %f) \n", p0.x, p0.y, p0.z, p1.x, p1.y, p1.z);

  boolean empty = emptyIntersectionTwoDisks(ca, va, ra, cb, vb, rb);
  if (empty) println("empty intersection");
  else println("non-empty intersection");
  return;
}

void testIntersectionCirclePlane() {
  assert rs.nRings >= 3;
  pt p0 = new pt();
  pt p1 = new pt();
  rs.generatePoints(attenuation);

  pt c0 = rs.centers[0];
  float r0 = rs.radii[0];
  vec n0 = rs.normals[0];
  pt c1 = rs.centers[1];
  pt c2 = rs.centers[2];
  vec d = normalToTriangle(c0, c1, c2);

  fill(orange, 100);
  showPlane(c0, d, r0);
  fill(red, 100);
  disk(c0, n0, r0);

  if (intersectionCirclePlane(c0, r0, n0, c0, d, p0, p1)) {
    fill(green, 100);
    show(p0, 3);
    fill(blue, 100);
    show(p1, 3);
  }
}