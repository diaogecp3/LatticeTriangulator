int numTests = 1000000;
int numPointsPerTest = 128;
boolean showResults = false;
float[] times = new float[numTests];  // in milliseconds

/*
 * Test the performance of convex hull generation for points on a sphere.
 * The number of tests done is n. For each test, input size is m. If showResult
 * is true, 3D results will be shown on screen.
 */
void testCH(int n, int m, boolean showResults) {
  for (int i = 0; i < n; ++i) {
    generatePointsOnSphere(P, centerOfSphere, radius, m);
    
    long startTime = System.nanoTime(); 
    ArrayList<Triangle> triangles = generateConvexHull(P.G, P.nv);
    long endTime = System.nanoTime();
    times[i] = (endTime - startTime) / 1000000.0;
    System.out.format("duration = %f ms. \n", times[i]);
    if (showResults) {
      fill(red); showTriangles(triangles, P.G);
      fill(cyan); showTriangleNormals(triangles, P.G);
    }
  }
  
  float avg = average(times, n, 1, n);
  System.out.format("Generate a convex hull for %d point, %d tests in total, average time = %f ms. \n", m, n, avg);
}