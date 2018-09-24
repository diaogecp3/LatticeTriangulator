/******************************************************************************
 * Playground.
 *
 * Test basic Processing/Java features.
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