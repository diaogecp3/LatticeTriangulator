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
        if (isZero(d(points[j], p))) {
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
        if (!emptyIntersectionTwoDisks(disk, disks[j])) {
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

vec constructNormal(vec v) {
  if (isAbsZero(v.y) && isAbsZero(v.z)) { // v is parallel to (1, 0, 0)
    return U(new vec(-v.z, 0.0, v.x));  // cross product of v and (0, 1, 0)
  } else {
    return U(new vec(0.0, v.z, -v.y));  // cross product of v and (1, 0, 0)
  }
}

vec[] generateInitDirs(pt C, pt[] contacts, int nc) {
  vec[] initDirs = new vec[nc];
  for (int i = 0; i < nc; ++i) {
    vec v = V(C, contacts[i]);
    initDirs[i] = constructNormal(v);
  }
  return initDirs;
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




// a center = the center of a disk
void showGroups(pt[] centers, pt[][] points, int nc, int np) {
  fill(orange);
  for (int i = 0; i < nc; ++i) {
    show(centers[i], 3);
  }
  fill(green);
  for (int i = 0; i < nc; ++i) {
    arrow(centers[i], V(centers[i], points[i][0]), 3);
  }
  fill(blue);
  for (int i = 0; i < nc; ++i) {
    for (int j = 0; j < np; ++j) {
      show(points[i][j], 3);
    }
  }
  fill(cyan);
  for (int i = 0; i < nc; ++i) {
    for (int j = 0; j < np; ++j) {
      collar(points[i][j], V(points[i][j], points[i][(j + 1) % np]), 1, 1);
    }
  }
  return;
}



pt[] generateContactsAndRadii(pt center,                             // in
                              float radius,                          // in
                              int n,                                 // in
                              float[] radii) {                       // out
  pt[] contacts = new pt[n];
  vec[] normals = new vec[n];
  float[] alphas = new float[n];
  for (int i = 0; i < n; ++i) {
    while (true) {
      pt p = generateOnePointOnSphere(center, radius);
      vec normal = U(center, p);
      float r = random(5, radius * 0.8);
      float alpha = asin(r/radius);  // [0, PI/2]
      boolean isValid = true;
      for (int j = 0; j < i; ++j) {
        float theta = acos(dot(normal, normals[j]));  // [0, PI]
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
                      float rMax) {                                  // in
  assert rMax > 0 && rMax < radius;
  pt[] contacts = new pt[n];
  vec[] normals = new vec[n];
  float alpha = 2 * asin(rMax/radius);
  for (int i = 0; i < n; ++i) {
    while (true) {
      pt p = generateOnePointOnSphere(center, radius);
      vec normal = U(center, p);
      boolean isValid = true;
      for (int j = 0; j < i; ++j) {
        float theta = acos(dot(normal, normals[j]));  // [0, PI]
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