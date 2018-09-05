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


pt[] generateContactsOnSphere(pt C, float R, float rMax, int nc) {
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

pt[] generatePointsForOneCircle(pt p, float r, pt C, float R, vec initDir, int np, pt center) {
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
pt[][] generatePointsForCircles(pt[] contacts, float r, pt C, float R, vec[] initDirs, int nc, int np, pt[] centers) {
  pt[][] points = new pt[nc][np];
  for (int i = 0; i < nc; ++i) {
    centers[i] = new pt();
    points[i] = generatePointsForOneCircle(contacts[i], r, C, R, initDirs[i], np, centers[i]);
  }
  return points;
}


// a center = the center of a disk
void showCircles(pt[] centers, pt[][] points, int nc, int np) {
  fill(red);
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