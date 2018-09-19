/******************************************************************************
 * Ring set processing.
 ******************************************************************************/


boolean debugFastCH = true;
int numFacesFastCH = 1;
int numStepsFastCH = 1;
int maxIterThreeRings = 100;
int maxIterOneRing = 100;

class Debug3RTriInfo {  // debug info about three-ring-triangle generation
  int idr0, idr1, idr2;
  int idp0, idp1, idp2;
  int numSteps;
  int numFacesShown;
  Debug3RTriInfo() {
    idr0 = idr1 = idr2 = -1;
    idp0 = idp1 = idp2 = -1;
    numSteps = 0;
    numFacesShown = 0;
  }
}


/*
 * RingSet class.
 *
 * A ring set is a set of rings (i.e. discritized cicles). We assume that these
 * rings lie on a sphere.
 */
class RingSet {
  pt c;  // the center of the sphere where the ring set lies
  float r;  // the radius of the sphere where the ring set lies
  int nc, np;  // number of rings/circles, and number of points on each ring
  boolean sameRadius;  // whether all rings have the same radius
  pt[] contacts;  // the intersections between the outward normals of rings and the sphere
  float[] radii;  // the radii of rings
  vec[] initDirs;  // the initial directions, one for each ring, used to generate the first point on ring

  pt[][] points;  // the generated points on rings
  pt[] centers;  // the centers of rings, one for each ring

  ArrayList<Triangle> triangles = null;  // triangle mesh with ring vertices
  ArrayList<Triangle> convexHullRef = null;  // convex hull generated by contacts

  Debug3RTriInfo debug3RTriInfo = new Debug3RTriInfo();  // for debug

  RingSet(pt c, float r) {
    this.c = c;
    this.r = r;
    sameRadius = false;
  }
  
  RingSet(pt c, float r, int nc, int np) {
    this.c = c;
    this.r = r;
    this.nc = nc;
    this.np = np;
    sameRadius = false;
  }
  
  RingSet(pt c, float r, int nc, int np, float rMax) {
    this.c = c;
    this.r = r;
    this.nc = nc;
    this.np = np;
    sameRadius = true;
    radii = new float[1];
    radii[0] = rMax;
  }
  
  void init() {
    if (!sameRadius) {
      radii = new float[nc];
      contacts = generateContactsAndRadii(c, r, nc, radii);
    } else {
      contacts = generateContacts(c, r, nc, radii[0]);
    }
    initDirs = generateInitDirs(c, contacts, nc);
  }
  
  void generatePoints(float attenuation) {
    centers = new pt[nc];
    if (!sameRadius) {
      float[] curRadii = new float[nc];
      for (int i = 0; i < nc; ++i) curRadii[i] = radii[i] * attenuation;
      points = generatePointsForCircles(contacts, curRadii, c, r, initDirs, nc,
                                        np, centers);
    } else {
      float curRadius = radii[0] * attenuation;
      points = generatePointsForCircles(contacts, curRadius, c, r, initDirs, nc,
                                        np, centers);
    }
  }

  void generateTriangleMesh() {
    triangles = generateConvexHull(points, nc, np);
  }

  void generateConvexHullRef() {
    assert contacts.length >= 3;
    convexHullRef = generateConvexHull(contacts, nc);
  }

  private int findStablePoint(pt[] points, int i, pt a, pt b, vec normal) {
    vec vn = normal;
    int d = dot(vn, V(points[i], points[(i + 1) % np])) > 0 ? 1 : -1;
    int inext = (i + d + np) % np;
    int iprev = (i - d + np) % np;
    int steps1 = 0;
    while ((dot(vn, V(points[i], points[inext])) > 0 ||
            dot(vn, V(points[i], points[iprev])) > 0) &&
            (steps1 < maxIterOneRing)) {
      if (debugFastCH &&
        debug3RTriInfo.numFacesShown == numFacesFastCH - 1 &&
        debug3RTriInfo.numSteps >= numStepsFastCH) break;
      iprev = i;
      i = inext;
      inext = (inext + d + np) % np;
      vn = N(points[i], a, b);
      if (debugFastCH) {
        debug3RTriInfo.numSteps++;
      }
      steps1++;
    }
    normal.setTo(vn);
    return i;
  }

  private Triangle generateThreeRingTriangle(pt[] points0,
                                             pt[] points1,
                                             pt[] points2) {
    int i = 1, j = 1, k = 1;
    vec normal = N(points0[i], points1[j], points2[k]);
    if (debugFastCH) {
      debug3RTriInfo.numSteps = 1;
    }
    int iter = 0;
    while (iter < maxIterThreeRings) {
      if (debugFastCH &&
          debug3RTriInfo.numFacesShown == numFacesFastCH - 1 &&
          debug3RTriInfo.numSteps >= numStepsFastCH) break;
      boolean noUpdate = true;
      /* Find ring A/B/C's stable point respectively. */
      int inew = findStablePoint(points0, i, points1[j], points2[k], normal);
      if (inew != i) {
        i = inew;
        if (noUpdate) noUpdate = false;
      }
      int jnew = findStablePoint(points1, j, points2[k], points0[i], normal);
      if (jnew != j) {
        j = jnew;
        if (noUpdate) noUpdate = false;
      }
      int knew = findStablePoint(points2, k, points0[i], points1[j], normal);
      if (knew != k) {
        k = knew;
        if (noUpdate) noUpdate = false;
      }
      if (noUpdate) break;
      iter++;
    }
    if (iter >= maxIterThreeRings) return null;
    return new Triangle(i, j, k);
  }

  ArrayList<Triangle> generateThreeRingTriangles() {
    int nt = convexHullRef.size();  // number of triangles of convex hull
    ArrayList<Triangle> threeRingTriangles = new ArrayList<Triangle>();
    debug3RTriInfo.numFacesShown = 0;
    for (int i = 0; i < nt; ++i) {
      if (debugFastCH && i >= numFacesFastCH) break;
      Triangle face = convexHullRef.get(i);
      Triangle t = generateThreeRingTriangle(points[face.a], points[face.b], points[face.c]);
      Triangle triangle = (t == null) ?
                          null :
                          new Triangle(face.a * np + t.a, face.b * np + t.b, face.c * np + t.c);
      threeRingTriangles.add(triangle);
      if (debugFastCH) {
        debug3RTriInfo.idr0 = face.a;
        debug3RTriInfo.idr1 = face.b;
        debug3RTriInfo.idr2 = face.c;
        debug3RTriInfo.numFacesShown++;
      }
    }
    return threeRingTriangles;
  }

  pt[][] get2DPointArray() {
    return points;
  }

  int getNumRings() {
    return nc;
  }

  int getNumPointsPerRing() {
    return np;
  }

  pt[] get1DPointArray() {
    if (points == null) return null;
    pt[] G = new pt[nc * np];
    int k = 0;
    for (int i = 0; i < nc; ++i) {
      for (int j = 0; j < np; ++j) {
        G[k++] = points[i][j];
      }
    }
    return G;
  }

  ArrayList<pt> get1DPointArrayList() {
    if (points == null) return null;
    ArrayList<pt> positions = new ArrayList<pt>();
    for (int i = 0; i < nc; ++i) {
      for (int j = 0; j < np; ++j) {
        positions.add(points[i][j]);
      }
    }
    return positions;
  }

  void showRings() {
    fill(orange);
    for (int i = 0; i < nc; ++i) {
      show(centers[i], 1);
    }
    fill(green);
    for (int i = 0; i < nc; ++i) {
      arrow(centers[i], V(centers[i], points[i][0]), 2);
    }
    fill(blue);
    for (int i = 0; i < nc; ++i) {
      for (int j = 0; j < np; ++j) {
        show(points[i][j], 2);
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
  
  void showDebug3RTriInfo() {
    fill(red, 150);
    show(contacts[debug3RTriInfo.idr0], 3);
    fill(green, 150);
    show(contacts[debug3RTriInfo.idr1], 3);
    fill(blue, 150);
    show(contacts[debug3RTriInfo.idr2], 3);
  }

  void saveRings(String file) {
    String[] lines = new String[2 + 3 * nc];
    int i = 0;
    lines[i++] = str(nc);
    lines[i++] = str(np);
    for (int j = 0; j < nc; ++j) {
      lines[i++] = str(contacts[j].x) + "," + str(contacts[j].y) + "," +
                   str(contacts[j].z);
      lines[i++] = str(radii[j]);
      lines[i++] = str(initDirs[j].x) + "," + str(initDirs[j].y) + "," +
                   str(initDirs[j].z);
    }
    saveStrings(file, lines);
    return;
  }

  void loadRings(String file) {
    String[] lines = loadStrings(file);
    int i = 0;
    nc = int(lines[i++]);
    np = int(lines[i++]);
    contacts = new pt[nc];
    radii = new float[nc];
    initDirs = new vec[nc];
    for (int j = 0; j < nc; ++j) {
      float[] contact = float(split(lines[i++], ","));
      contacts[j] = new pt(contact[0], contact[1], contact[2]);
      radii[j] = float(lines[i++]);
      float[] initDir = float(split(lines[i++], ","));
      initDirs[j] = new vec(initDir[0], initDir[1], initDir[2]);
    }
    return;
  }
}