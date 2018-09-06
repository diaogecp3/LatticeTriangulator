/*********************************************************
 * Convex hull generation in a special case where all
 * vertices are on a sphere.
 *********************************************************/


import java.util.LinkedList;

int numFacesShown = 0;
boolean debugCH = true;
float disturbance = 0.0001;



/*
 * Initialize a vertex list given a point (position) array. A vertex may
 * contain more info than a position.
 */
ArrayList<Vertex> convertToVertexList(pt[] G, int nv, int k) {
  ArrayList<Vertex> vertices = new ArrayList<Vertex>();
  if (k > 0) {
    for (int i = 0; i < nv; ++i) {
      vec epsilon = random3(disturbance);
      pt p = P(G[i], epsilon);
      vertices.add(new Vertex(i, p));
    }
  } else {
    for (int i = 0; i < nv; ++i) {
      vertices.add(new Vertex(i, G[i]));
    }
  }
  return vertices;
}


/*
 * Check if a vertex is inner given a list of front edges.
 */
boolean isInnerVertex(LinkedList<TaggedEdge> frontEdges, Vertex v) {
  for (TaggedEdge edge : frontEdges) {
    if (edge.tag == 0 && (edge.e.a == v.id || edge.e.b == v.id)) return false;
  }
  return true;
}

/*
 * Find the first valid triangle (i.e. on the surface of the convex hull)
 * from input vertices. Because we assume that all vertices will be on the
 * surface of the convex hull (actually on a sphere), the lowest 2 vertices
 * (w.r.t. z-axis) must be on a valid triangle. Check each other vertex and
 * pick the one that can help form a valid triangle.
 */
Triangle findFirstTriangle(ArrayList<Vertex> vertices,
                           LinkedList<TaggedEdge> frontEdges,
                           boolean[][] manifoldMask) {
  int a = 0, b = 1, n = vertices.size();
  if (vertices.get(a).position.z > vertices.get(b).position.z) {
    int tmp = a;
    a = b;
    b = tmp;
  }
  for (int i = 2; i < n; ++i) {
    if (vertices.get(i).position.z < vertices.get(a).position.z) {
      b = a;
      a = i;
    } else if (vertices.get(i).position.z < vertices.get(b).position.z) {
      b = i;
    }
  }

  Vertex va, vb, vc;
  va = vertices.get(a);
  vb = vertices.get(b);
  int c = -1;
  for (int i = 0; i < n; ++i) {
    if (i == a || i == b) continue;
    Vertex vi = vertices.get(i);
    vec N = N(va.position, vb.position, vi.position);
    
    int j = 0;
    for (; j < n; ++j) {
      if (j != a && j != b && j != i) break;
    }

    Vertex vj = vertices.get(j);
    vec aj = V(va.position, vj.position);
    boolean sign = dot(N, aj) > 0;
    boolean isValid = true;
    for (int k = 0; k < n; ++k) {
      if (k == a || k == b || k == i || k == j) continue;
      Vertex vk = vertices.get(k);
      vec ak = V(va.position, vk.position);
      boolean tmpSign = dot(N, ak) > 0;
      if (sign != tmpSign) {
        isValid = false;
        break;
      }
    }
    
    if (isValid == true) {
      c = i;
      if (sign) {
        int tmp = b;
        b = c;
        c = tmp;
      }
      break;
    }
  }
  //assert c >= 0;
  if (c < 0) {
    println("cannot find c in first triangle");
    exceptionHandler();
    if (debugCH) {
      fill(#DD08FA); show(vertices.get(a).position, 5);
      fill(#08FA31); show(vertices.get(b).position, 5);
      numFacesShown = 0;
      return null;
    } else exit();
  }
  vc = vertices.get(c);
  TaggedEdge e0, e1, e2;
  e0 = new TaggedEdge(a, b, c, 0);
  e1 = new TaggedEdge(b, c, a, 0);
  e2 = new TaggedEdge(c, a, b, 0);
  frontEdges.add(e0);
  frontEdges.add(e1);
  frontEdges.add(e2);
  va.outEdges.put(vb, e0);
  vb.outEdges.put(vc, e1);
  vc.outEdges.put(va, e2);
  manifoldMask[a][b] = manifoldMask[b][c] = manifoldMask[c][a] = true;
  return new Triangle(a, b, c);
}

/*
 * Pivot around edge e to find the next valid triangle face. We only need to
 * find the next vertex such that it and edge e form a triangle face of the
 * surface of the convex hull.
 */
int pivot(TaggedEdge e, ArrayList<Vertex> vertices, boolean[][] manifoldMask, Integer nVertices) {
  int a = e.e.a, b = e.e.b, c = e.c;
  pt A = vertices.get(a).position, B = vertices.get(b).position, C = vertices.get(c).position;
  vec AC = V(A, C), AB = V(A, B);
  vec normalizedAB = U(AB);
  vec dir0 = U(M(AC, V(dot(AC, normalizedAB), normalizedAB)));  // normalize(AC - (AC dot (normalize(AB)) times (normalize(AB)))

  float minCosTheta = 1.5f;
  int d = -1;
  int nv = nVertices == null ? vertices.size() : nVertices.intValue();
  int gid = -1;
  if (vertices.get(a).groupID != -1 && 
      vertices.get(a).groupID == vertices.get(b).groupID) {
    gid = vertices.get(a).groupID;
  }
  for (int i = 0; i < nv; ++i) {
    if (i == a || i == b || i == c) continue;
    if (vertices.get(i).isInner) continue;
    if (manifoldMask[a][i] || manifoldMask[i][b]) continue;
    if (gid >= 0 && vertices.get(i).groupID == gid) continue;
    pt D = vertices.get(i).position;
    vec AD = V(A, D);
    vec dir1 = U(M(AD, V(dot(AD, normalizedAB), normalizedAB)));
    float cosTheta = dot(dir0, dir1);
    if (cosTheta < minCosTheta) {  // this is OK, the angle should be [0, pi]
      minCosTheta = cosTheta;
      d = i;
    }
  }

  return d;
}


/*
 * Generate a convex hull from a list of vertices. Assume that all parameters
 * except the first one, are initialized by default values (empty/zeros, etc).
 * If no error or strange situation encountered, this function returns true,
 * otherwise returns false.
 */
boolean generateConvexHullFromVertices(ArrayList<Vertex> vertices,
                                       ArrayList<Triangle> triangles,
                                       LinkedList<TaggedEdge> frontEdges,
                                       boolean[][] manifoldMask,
                                       DebugInfo debugInfo) {
  // Find the first valid triangle
  triangles.add(findFirstTriangle(vertices, frontEdges, manifoldMask));
  debugInfo.numSteps = 1;
  
  // Expand the triangle mesh until the end
  while(frontEdges.size() > 0) {
    // Break after creating the first numFacesShown faces
    if (debugCH && debugInfo.numSteps >= numFacesShown) break;

    TaggedEdge e = frontEdges.poll();  // get and remove the first front edge
    if (e.tag != 0) continue;  // skip non-front edge
    int d = pivot(e, vertices, manifoldMask, null);

    if (d < 0) {
      println("No hit found, number of steps = " + debugInfo.numSteps);
      exceptionHandler();
      return false;
    }

    int a = e.e.a, b = e.e.b;
    Vertex va = vertices.get(a);
    Vertex vb = vertices.get(b);
    Vertex vd = vertices.get(d);

    triangles.add(new Triangle(a, d, b));  // create a new triangle facet
    e.tag = 2;  // the front edge becomes an inner edge  // TODO: seems that no need to tag edge? if use manifold mask
    va.outEdges.remove(vb);  // remove edge ab
    manifoldMask[a][d] = manifoldMask[d][b] = manifoldMask[b][a] = true;

    if (debugCH) {
      debugInfo.a = a;  // for debugging
      debugInfo.b = b;  // for debugging
      debugInfo.d = d;  // for debugging
    }

    // Check if there is a front edge connecting va and vd
    TaggedEdge e0 = null;  
    if (vd.outEdges.containsKey(va)) {  // can't just check next valid edge in front list
      e0 = vd.outEdges.get(va);
      vd.outEdges.remove(va);
    }
    if (va.outEdges.containsKey(vd)) { // Is this possible? not possible when no 4 or more points on a plane?
      println("edge AD already exists, number of steps = " + debugInfo.numSteps);
      exceptionHandler();
      return false;
    }

    // Check if there is a front edge connecting vd and vb
    TaggedEdge e1 = null;
    if (vb.outEdges.containsKey(vd)) {  // can't just check prev valid edge in front list
      e1 = vb.outEdges.get(vd);
      vb.outEdges.remove(vd);
    }
    if (vd.outEdges.containsKey(vb)) { // Is this possible?
      println("edge DB already exists, number of steps = " + debugInfo.numSteps);
      exceptionHandler();
      return false;
    }

    if (e0 == null) { // create a new e0 and push it to the front list
      e0 = new TaggedEdge(a, d, b, 0);
      va.outEdges.put(vd, e0);
      frontEdges.add(e0);
    } else {  // e0 already exists, remove it from front list
      e0.tag = 2;
    }

    if (e1 == null) { // create a new e1 and push it to the front list
      e1 = new TaggedEdge(d, b, a, 0);
      vd.outEdges.put(vb, e1);
      frontEdges.add(e1);
    } else {  // e1 already exists, remove it from front list
      e1.tag = 2;
    }

    /*
     * Check if va becomes an inner vertex when e0 is removed from front list.
     * Check if vb becomes an inner vertex when e1 is removed from front list.
     * Check if vd becomes an inner vertex when e0 and e1 are removed from front list.
     */
    if (e0.tag == 2) {
      if (isInnerVertex(frontEdges, va)) va.isInner = true;
    }
    if (e1.tag == 2) {
      if (isInnerVertex(frontEdges, vb)) vb.isInner = true;
    }
    if (e0.tag == 2 && e1.tag == 2) {
      if (isInnerVertex(frontEdges, vd)) vd.isInner = true;
    }

    if (debugCH) debugInfo.numSteps++;  // increase steps or number of faces to be shown
  }  // end while
  return true;
}

/*
 * Generate the convex hull for a given point array G whose size is nv. Assume
 * that all input points will be on the surface of the convex hull.
 */
ArrayList<Triangle> generateConvexHull(pt[] G, int nv) {
  assert nv >= 4;
  ArrayList<Vertex> vertices;
  ArrayList<Triangle> triangles;
  LinkedList<TaggedEdge> frontEdges; // maybe using an array would be better, since edge DA, BD would be prev/next to AB, and hence no need to use outEdges for each vertex on the front
  boolean[][] manifoldMask;
  DebugInfo debugInfo;
  int k = 0;
  while (true) {
    vertices = convertToVertexList(G, nv, k);  // vertices may differ every time
    triangles = new ArrayList<Triangle>();
    frontEdges = new LinkedList<TaggedEdge>();
    manifoldMask = new boolean[nv][nv];
    debugInfo = new DebugInfo();

    boolean success = generateConvexHullFromVertices(vertices, triangles,
      frontEdges, manifoldMask, debugInfo);
    
    boolean passQT = passQualityTest(triangles, G, nv);
    if (success && !passQT) {
      println("generate CH successfully but not pass quality test");
    }
    k++;
    if (success && passQT) break;
  }
  if (k > 1) System.out.format("number of times tried = %d\n", k);
  
  if (debugCH) {
    fill(blue); showFrontEdges(frontEdges, G);
    fill(black); showInnerVertices(vertices, G);
    fill(#9AFF05); showManifoldEdges(manifoldMask, G, nv);  // light green
    if (debugInfo.a >= 0) { fill(#970EED, 160); show(G[debugInfo.a], 3); }  // purple
    if (debugInfo.b >= 0) { fill(#21C2FA, 160); show(G[debugInfo.b], 3); }  // light blue
    if (debugInfo.d >= 0) { fill(#FFF705, 160); show(G[debugInfo.d], 3); }  // light yellow
  }
  
  return triangles;
}


pt[] convertTo1DArray(pt[][] points, int nc, int np) {
  pt[] G = new pt[nc * np];
  int k = 0;
  for (int i = 0; i < nc; ++i) {
    for (int j = 0; j < np; ++j) {
      G[k++] = points[i][j];
    }
  }
  return G;
}


ArrayList<Vertex> convertToVertexList(pt[][] points, pt[] contacts, int nc, int np, int k) {
  ArrayList<Vertex> vertices = new ArrayList<Vertex>();
  if (k > 0) {
    for (int i = 0; i < nc; ++i) {
      for (int j = 0; j < np; ++j) {
        vec epsilon = random3(disturbance);
        pt position = P(points[i][j], epsilon);
        Vertex vertex = new Vertex(i * np + j, position, i);
        vertices.add(vertex);
      }
    }
  } else {
    for (int i = 0; i < nc; ++i) {
      for (int j = 0; j < np; ++j) {
        Vertex vertex = new Vertex(i * np + j, points[i][j], i);
        vertices.add(vertex);
      }
    }
  }
  int nv = nc * np;
  for (int i = 0; i < nc; ++i) {
    Vertex vertex = new Vertex(nv + i, contacts[i], i);
    vertices.add(vertex);
  }
  return vertices;
}

void initFrontsAndManifoldMask(int nc, int np, ArrayList<Vertex> vertices, Front[] fronts, boolean[][] manifoldMask) {
  int nv = nc * np;
  for (int i = 0; i < nc; ++i) {
    LinkedList<TaggedEdge> edges = new LinkedList<TaggedEdge>();
    HashSet<Integer> groupIDs = new HashSet<Integer>();
    groupIDs.add(new Integer(i));
    int head = i * np;
    int contact = nv + i;
    for (int j = 0; j < np - 1; ++j) {
      int current = head + j;
      int next = current + 1;
      TaggedEdge e = new TaggedEdge(current, next, contact);
      manifoldMask[current][next] = true;
      edges.add(e);
      vertices.get(current).outEdges.put(vertices.get(next), e);
    }
    int tail = head + np - 1;
    TaggedEdge lastEdge = new TaggedEdge(tail, head, contact);
    edges.add(lastEdge);
    vertices.get(tail).outEdges.put(vertices.get(head), lastEdge);
    manifoldMask[tail][head] = true;
    fronts[i] = new Front(edges, groupIDs);
  }
  return;
}


boolean generateConvexHullFromVertices(ArrayList<Vertex> vertices,
                                       ArrayList<Triangle> triangles,
                                       Front[] fronts,
                                       boolean[][] manifoldMask,
                                       int nc, int np,
                                       DebugInfo debugInfo) {
  Front curFront = fronts[0];
  assert curFront.size() > 0;
  Integer nv = new Integer(nc * np);
  debugInfo.numSteps = 0;
  while (curFront.size() > 0) {
    if (debugCH && debugInfo.numSteps >= numFacesShown) break;
    TaggedEdge e = curFront.poll();
    if (e.tag != 0) continue;
    int d = pivot(e, vertices, manifoldMask, nv);
    
    if (d < 0) {
      println("No hit found, number of steps = " + debugInfo.numSteps);
      exceptionHandler();
      return false;
    }
    
    if (debugCH) {
      debugInfo.a = e.e.a;  // for debugging
      debugInfo.b = e.e.b;  // for debugging
      debugInfo.d = d;  // for debugging      
    }
    
    int a = e.e.a, b = e.e.b;
    Vertex va = vertices.get(a);
    Vertex vb = vertices.get(b);
    Vertex vd = vertices.get(d);

    triangles.add(new Triangle(a, d, b));  // create a new triangle facet
    e.tag = 2;  // the front edge becomes an inner edge  // TODO: seems that no need to tag edge? if use manifold mask
    va.outEdges.remove(vb);  // remove edge ab
    manifoldMask[a][d] = manifoldMask[d][b] = manifoldMask[b][a] = true;
    
    // Check if there is a front edge connecting va and vd
    TaggedEdge e0 = null;
    if (vd.outEdges.containsKey(va)) {
      e0 = vd.outEdges.get(va);
      vd.outEdges.remove(va);
    }
    
    if (va.outEdges.containsKey(vd)) { // Is this possible? not possible when no 4 or more points on a plane?
      println("edge AD already exists, number of steps = " + debugInfo.numSteps);
      exceptionHandler();
      return false;
    }
    
    // Check if there is a front edge connecting vd and vb
    TaggedEdge e1 = null;
    if (vb.outEdges.containsKey(vd)) {
      e1 = vb.outEdges.get(vd);
      vb.outEdges.remove(vd);
    }
    
    if (vd.outEdges.containsKey(vb)) { // Is this possible?
      println("edge DB already exists, number of steps = " + debugInfo.numSteps);
      exceptionHandler();
      return false;
    }
    
    
    if (e0 == null) { // create a new e0 and push it to the front list
      e0 = new TaggedEdge(a, d, b, 0);
      va.outEdges.put(vd, e0);
      curFront.add(e0);
    } else {  // e0 already exists, remove it from front list
      e0.tag = 2;
    }

    if (!curFront.containGroupID(vd.groupID)) {
      curFront.mergeWithFront(vd, fronts);
    }

    if (e1 == null) { // create a new e1 and push it to the front list
      e1 = new TaggedEdge(d, b, a, 0);
      vd.outEdges.put(vb, e1);
      curFront.add(e1);
    } else {  // e1 already exists, remove it from front list
      e1.tag = 2;
    }
    
    if (e0.tag == 2) {
      if (curFront.isInnerVertex(va)) va.isInner = true;
    }
    if (e1.tag == 2) {
      if (curFront.isInnerVertex(vb)) vb.isInner = true;
    }
    if (e0.tag == 2 && e1.tag == 2) {
      if (curFront.isInnerVertex(vd)) vd.isInner = true;
    }
    
    if (debugCH) debugInfo.numSteps++;
  }
  
  return true;
}

ArrayList<Triangle> generateConvexHull(pt[][] points, pt[] contacts, int nc, int np) {
  ArrayList<Vertex> vertices;
  ArrayList<Triangle> triangles;
  Front[] fronts;
  boolean[][] manifoldMask;
  DebugInfo debugInfo;
  int nv = nc * np;
  int k = 0;
  pt[] G = convertTo1DArray(points, nc, np);
  while (true) {
    vertices = convertToVertexList(points, contacts, nc, np, k);
    triangles = new ArrayList<Triangle>();
    fronts = new Front[nc];
    manifoldMask = new boolean[nv][nv];
    debugInfo = new DebugInfo();
    initFrontsAndManifoldMask(nc, np, vertices, fronts, manifoldMask);
    
    boolean success = generateConvexHullFromVertices(vertices, triangles, fronts, manifoldMask, nc, np, debugInfo);
    boolean passQT = passQualityTest(triangles, G, nv);
    
    if (success && !passQT) {
      println("generate CH (from a list of fronts) successfully but not pass quality test");
    }
    if (success && passQT) break;
    k++;
  }
  
  if (debugCH) {
    fill(blue); fronts[0].showEdges(G);
    fill(#A08888);  // light grey
    for (int i = 1; i < nc; ++i) {
      if (fronts[i].size() > 0) fronts[i].showEdges(G);
    }
    fill(black); showInnerVertices(vertices, G);
    fill(#9AFF05); showManifoldEdges(manifoldMask, G, nv);  // light green
    if (debugInfo.a >= 0) { fill(#970EED, 160); show(G[debugInfo.a], 3); }  // purple
    if (debugInfo.b >= 0) { fill(#21C2FA, 160); show(G[debugInfo.b], 3); }  // light blue
    if (debugInfo.d >= 0) { fill(#FFF705, 160); show(G[debugInfo.d], 3); }  // light yellow
  }
  
  return triangles;
}


/*
 * Show triangles.
 */
void showTriangles(ArrayList<Triangle> triangles, pt[] G) {
  int n = triangles.size();
  beginShape(TRIANGLES);
  for (int i = 0; i < n; ++i) {
    pt A, B, C;
    A = G[triangles.get(i).a];
    B = G[triangles.get(i).b];
    C = G[triangles.get(i).c];
    vertex(A);
    vertex(B);
    vertex(C);
  }
  endShape();
}

/*
 * Show front edges.
 */
void showFrontEdges(LinkedList<TaggedEdge> frontEdges, pt[] G) {
  for (TaggedEdge edge : frontEdges) {
    if (edge.tag != 0) continue;
    pt A, B;
    A = G[edge.e.a];
    B = G[edge.e.b];
    arrow(A, V(A, B), 2);
  }
}

/*
 * Show normals of triangles, one normal per triangle face.
 */
void showTriangleNormals(ArrayList<Triangle> triangles, pt[] G) {
  int n = triangles.size();
  for (int i = 0; i < n; ++i) {
    pt A, B, C;
    A = G[triangles.get(i).a];
    B = G[triangles.get(i).b];
    C = G[triangles.get(i).c];
    showNormalToTriangle(A, B, C, 10, 1);
  }
}

/*
 * Show inner vertices
 */
void showInnerVertices(ArrayList<Vertex> vertices, pt[] G) {
  int n = vertices.size();
  for (int i = 0; i < n; ++i) {
    if (vertices.get(i).isInner) {
      show(G[vertices.get(i).id], 2);
    }
  }
}

/*
 * Show manifold edges, i.e. those with exactly 2 adjacent faces.
 */
void showManifoldEdges(boolean[][] manifoldMask, pt[] G, int nv) {
  for (int i = 0; i < nv; ++i) {
    for (int j = i + 1; j < nv; ++j) {
      if (manifoldMask[i][j] && manifoldMask[j][i]) {
        collar(G[i], V(G[i], G[j]), 1, 1);
      }
    }
  }
}