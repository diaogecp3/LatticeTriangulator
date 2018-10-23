/*********************************************************
 * Input generation for convex hulls.
 *********************************************************/


pt generateOnePointOnSphere(pt c, float r) {
  float alpha = random(-HALF_PI, HALF_PI);
  float beta = random(0, TWO_PI);
  float dx = r * cos(alpha) * cos(beta);
  float dy = r * cos(alpha) * sin(beta);
  float dz = r * sin(alpha);
  return new pt(c.x + dx, c.y + dy, c.z + dz);
}

pt[] generatePointsOnSphere(pt c, float r, int n) {
  pt[] points = new pt[n];
  for (int i = 0; i < n; ++i) {
    while (true) {
      pt p = generateOnePointOnSphere(c, r);
      boolean bad = false;
      for (int j = 0; j < i; ++j) {
        if (isNonPositive(d(points[j], p))) {
          bad = true;
          break;
        }
      }
      if (!bad) {
        points[i] = p;
        break;
      }
    }  // end while
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


pt[] generateContactsOnSphere(pt C, float R, int nc, float rMax) {
  assert rMax < R;
  pt[] points = new pt[nc];
  Disk[] disks = new Disk[nc];
  float distance = sqrt(R * R - rMax * rMax);
  for (int i = 0; i < nc; ++i) {
    while (true) {
      pt p = generateOnePointOnSphere(C, R);
      vec normal = U(C, p);
      pt c = P(C, distance, normal);  // push p towards C
      Disk disk = new Disk(c, normal, rMax);
      boolean bad = false;
      for (int j = 0; j < i; ++j) {
        if (!emptyIntersectionTwoDisks(disk, disks[j])) {  // too expensive
          bad = true;
          break;
        }
      }
      if (!bad) {
        points[i] = p;
        disks[i] = disk;
        break;
      }
    }  // end while
  }
  return points;
}




pt[] generatePointsForOneCircle(pt p,                                // in
                                float r,                             // in
                                pt C,                                // in
                                float R,                             // in
                                vec initDir,                         // in
                                int np,                              // in
                                pt center) {                         // out
  pt[] points = new pt[np];
  vec normal = U(C, p);
  vec v = I = initDir;
  vec J = N(normal, I);  // the order matters! make sure cross(I, J) same as normal
  center.set(P(C, sqrt(R * R - r * r), normal));
  float da = TWO_PI / np;
  float a = 0;
  for (int i = 0; i < np; ++i) {
    points[i] = P(center, r, R(v, a, I, J));
    a += da;
  }
  return points;
}


// a contact = the intersection point between the medial axis of a tube and the sphere
// a center = the center of a circle
pt[][] generatePointsForCircles(pt[] contacts,                       // in
                                float r,                             // in
                                pt C,                                // in
                                float R,                             // in
                                vec[] initDirs,                      // in
                                int nc,                              // in
                                int np,                              // in
                                pt[] centers) {                      // out
  pt[][] points = new pt[nc][np];
  for (int i = 0; i < nc; ++i) {
    centers[i] = new pt();
    points[i] = generatePointsForOneCircle(contacts[i], r, C, R, initDirs[i],
                                           np, centers[i]);
  }
  return points;
}

pt[][] generatePointsForCircles(pt[] contacts,                       // in
                                float[] radii,                       // in
                                pt C,                                // in
                                float R,                             // in
                                vec[] initDirs,                      // in
                                int nc,                              // in
                                int np,                              // in
                                pt[] centers) {                      // out
  pt[][] points = new pt[nc][np];
  for (int i = 0; i < nc; ++i) {
    centers[i] = new pt();
    points[i] = generatePointsForOneCircle(contacts[i], radii[i], C, R,
                                           initDirs[i], np, centers[i]);
  }
  return points;
}


pt[] generateContactsAndRadii(pt center,                             // in
                              float radius,                          // in
                              int n,                                 // in
                              float[] radii,                         // out
                              vec[] normals) {                       // out
  pt[] contacts = new pt[n];
  float[] alphas = new float[n];
  for (int i = 0; i < n; ++i) {
    while (true) {
      pt p = generateOnePointOnSphere(center, radius);
      vec normal = U(center, p);
      float r = random(5, radius * 0.8);
      float alpha = asin(clamp(r/radius, -1.0, 1.0));  // [0, PI/2]
      boolean isValid = true;
      for (int j = 0; j < i; ++j) {
        float theta = acos(clamp(dot(normal, normals[j]), -1.0, 1.0));  // [0, PI]
        if (theta <= alpha + alphas[j]) {
          isValid = false;
          break;
        }
      }
      if (isValid) {
        contacts[i] = p;
        normals[i] = normal;
        alphas[i] = alpha;
        radii[i] = r;
        break;
      }
    }
  }
  return contacts;
}


pt[] generateContacts(pt center,                                     // in
                      float radius,                                  // in
                      int n,                                         // in
                      float rMax,                                    // in
                      vec[] normals) {                               // out
  assert rMax > 0 && rMax < radius;
  pt[] contacts = new pt[n];
  float alpha = 2 * asin(clamp(rMax/radius, -1.0, 1.0));
  for (int i = 0; i < n; ++i) {
    while (true) {
      pt p = generateOnePointOnSphere(center, radius);
      vec normal = U(center, p);
      boolean isValid = true;
      for (int j = 0; j < i; ++j) {
        float theta = acos(clamp(dot(normal, normals[j]), -1.0, 1.0));  // [0, PI]
        if (theta <= alpha) {
          isValid = false;
          break;
        }
      }
      if (isValid) {
        contacts[i] = p;
        normals[i] = normal;
        break;
      }
    }
  }
  return contacts;
}

Ball[] generateBallsOnSphere(pt c, float r, float rMax, int n) {
  Ball[] balls = new Ball[n];
  assert rMax >= 4.0;
  for (int i = 0; i < n; ++i) {
    while (true) {
      pt px = generateOnePointOnSphere(c, r);
      float rx = random(4, rMax);
      boolean isValid = true;
      for (int j = 0; j < i; ++j) {
        if (d(px, balls[j].c) < rx + balls[j].r) {
          isValid = false;
          break;
        }
      }
      if (isValid) {
        Ball ball = new Ball(px, rx);
        balls[i] = ball;
        break;
      }
    }
  }
  return balls;
}

Hub generateHub(pt c, float r0, float r1, int n) {
  Ball ball = new Ball(c, r0);
  Ball[] neighbors = generateBallsOnSphere(c, r1, r1 - 1.3 * r0, n);
  return new Hub(ball, neighbors, n);
}