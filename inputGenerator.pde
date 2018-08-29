pt[] generatePointsOnSphere(pt c, float r, int n) {
  pt[] points = new pt[n];
  for (int i = 0; i < n; ++i) {
    boolean notGoodEnough = true;
    while (notGoodEnough) {
      float alpha = random(-HALF_PI, HALF_PI);
      float beta = random(0, TWO_PI);
      float dx = r * cos(alpha) * cos(beta);
      float dy = r * cos(alpha) * sin(beta);
      float dz = r * sin(alpha);
      pt p = new pt(c.x + dx, c.y + dy, c.z + dz);
      
      boolean bad = false;
      for (int j = 0; j < i; ++j) {
        if (d(points[j], p) < 0.001) {
          bad = true;
          break;
        }
      }
      if (!bad) {
        notGoodEnough = false;
        points[i] = p;
      }
    }
  }
  return points;
}

void generatePointsOnSphere(pts P, pt c, float r, int n) {
  P.empty();
  pt[] points = generatePointsOnSphere(c, r, n);
  for (int i = 0; i < n; ++i) {
    P.addPt(points[i]);
  }
}