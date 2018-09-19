/******************************************************************************
 * Convex hull generation in a special case where all vertices are on a sphere.
 ******************************************************************************/


/*
 * TODO:
 * 1). manifoldMask is a n by n sparse matrix. Use an array of hash map can
 *     decrease space complexity without losing too much speed.
 */


boolean debugCH = false;
int numFaces = 0;
float disturbance = 0.0001;

/*
 * Convert a point array to a vertex array list. Disturb positions if needed.
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
 * Convert a point array to a vertex array list. Associate a group ID
 * for each vertex. Disturb positions if needed.
 */
ArrayList<Vertex> convertToVertexList(pt[][] points, int nr, int nc, int k) {
  ArrayList<Vertex> vertices = new ArrayList<Vertex>();
  if (k > 0) {
    for (int i = 0; i < nr; ++i) {
      for (int j = 0; j < nc; ++j) {
        vec epsilon = random3(disturbance);
        pt position = P(points[i][j], epsilon);
        Vertex vertex = new Vertex(i * nc + j, position, i);
        vertices.add(vertex);
      }
    }
  } else {
    for (int i = 0; i < nr; ++i) {
      for (int j = 0; j < nc; ++j) {
        Vertex vertex = new Vertex(i * nc + j, points[i][j], i);
        vertices.add(vertex);
      }
    }
  }
  return vertices;
}

/*
 * Find the first valid triangle (i.e. on the surface of the convex hull)
 * from input vertices. Because we assume that all vertices will be on the
 * surface of the convex hull (actually on a sphere), the lowest 2 vertices
 * (w.r.t. z-axis) must be on a valid triangle. Check each other vertex and
 * pick the one that can help form a valid triangle.
 */
boolean findFirstTriangle(ArrayList<Vertex> vertices,                // in/out
                          ArrayList<Triangle> triangles,             // out
                          Front front,                               // out
                          boolean[][] manifoldMask,                  // out
                          DebugCHInfo debugInfo) {                     // out
  assert front.size() == 0 && triangles.size() == 0 && vertices.size() > 3;
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
  vb = vertices.get(b);  // b can be changed
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

  if (c < 0) {
    println("cannot find c to form the first triangle");
    exceptionHandler();
    exit();
  }
  vb = vertices.get(b);
  vc = vertices.get(c);
  vec normalABC = normalToTriangle(va.position, vb.position, vc.position);
  FrontEdge e0, e1, e2;
  e0 = new FrontEdge(a, b, c, normalABC);
  e1 = new FrontEdge(b, c, a, normalABC);
  e2 = new FrontEdge(c, a, b, normalABC);
  front.add(e0);
  front.add(e1);
  front.add(e2);
  va.outEdges.put(vb, e0);
  vb.outEdges.put(vc, e1);
  vc.outEdges.put(va, e2);
  manifoldMask[a][b] = manifoldMask[b][c] = manifoldMask[c][a] = true;
  triangles.add(new Triangle(a, b, c));
  debugInfo.numFaces = 1;
  return true;
}

/*
 * Initialize data structures (e.g. fronts, manifoldMask) related to convex
 * hull generation.
 */
void initFronts(int nGroups,                                         // in    
                int nPointsPerGroup,                                 // in
                ArrayList<Vertex> vertices,                          // in/out
                Front[] fronts,                                      // out
                boolean[][] manifoldMask,                            // out
                DebugCHInfo debugInfo) {                               // out
  for (int i = 0; i < nGroups; ++i) {
    LinkedList<FrontEdge> edges = new LinkedList<FrontEdge>();
    HashSet<Integer> groupIDs = new HashSet<Integer>();
    groupIDs.add(new Integer(i));
    int head = i * nPointsPerGroup;
    
    // Compute the normal of current face
    pt A = vertices.get(head).position;
    pt B = vertices.get(head + 1).position;
    pt C = vertices.get(head + 2).position;
    vec N = normalToTriangle(A, B, C);
    
    for (int j = 0; j < nPointsPerGroup - 1; ++j) {
      int current = head + j;
      int next = current + 1;
      FrontEdge e = new FrontEdge(current, next, -1, N);
      manifoldMask[current][next] = true;
      edges.add(e);
      vertices.get(current).outEdges.put(vertices.get(next), e);
    }
    int tail = head + nPointsPerGroup - 1;
    FrontEdge lastEdge = new FrontEdge(tail, head, -1, N);
    edges.add(lastEdge);
    vertices.get(tail).outEdges.put(vertices.get(head), lastEdge);
    manifoldMask[tail][head] = true;
    fronts[i] = new Front(edges, groupIDs);
  }
  debugInfo.numFaces = 0;
  return;
}

/*
 * Pivot around edge e to find the next valid triangle face.
 */
int pivot(FrontEdge e,                                               // in/out
          ArrayList<Vertex> vertices,                                // in/out
          boolean[][] manifoldMask,                                  // in/out
          vec N) {                                                   // out
  int a = e.a, b = e.b, c = e.c;
  Vertex va = vertices.get(a), vb = vertices.get(b);
  pt A = va.position, B = vb.position;
  vec normalABC = e.N;
  
  float maxCosTheta = -2.0f;
  int d = -1, gid = -1;
  int nv = vertices.size();
  if (va.groupID != -1 && va.groupID == vb.groupID) {
    gid = va.groupID;
  }
  for (int i = 0; i < nv; ++i) {
    if (i == a || i == b || i == c) continue;
    if (vertices.get(i).isInner) continue;
    if (manifoldMask[a][i] || manifoldMask[i][b]) continue;
    if (gid >= 0 && vertices.get(i).groupID == gid) continue;
    pt D = vertices.get(i).position;
    vec normalADB = normalToTriangle(A, D, B);
    float cosTheta = dot(normalABC, normalADB);
    if (cosTheta > maxCosTheta) {  // the angle should be [0, pi]
      maxCosTheta = cosTheta;
      d = i;
      N.setTo(normalADB);
    }
  }
  return d;
}

/*
 * Generate a convex hull from a list of vertices. If no error or strange
 * situation encountered, this function returns true, otherwise returns false.
 */
boolean generateConvexHullWithFronts(ArrayList<Vertex> vertices,     // in/out
                                     ArrayList<Triangle> triangles,  // out
                                     Front[] fronts,                 // out
                                     boolean[][] manifoldMask,       // out
                                     DebugCHInfo debugInfo) {          // out
  Front curFront = fronts[0];
  assert curFront.size() > 0;
  while (curFront.size() > 0) {
    if (debugCH && debugInfo.numFaces >= numFaces) break;
    FrontEdge e = curFront.poll();
    if (!e.isValid) continue;
    vec normalADB = new vec();
    int d = pivot(e, vertices, manifoldMask, normalADB);
    
    if (d < 0) {
      if (debugCH) println("No hit found, number of steps = " + debugInfo.numFaces);
      else println("No hit found");
      exceptionHandler();
      return false;
    }

    int a = e.a, b = e.b;
    Vertex va = vertices.get(a);
    Vertex vb = vertices.get(b);
    Vertex vd = vertices.get(d);

    if (debugCH) {
      debugInfo.a = a;  // for debugging
      debugInfo.b = b;  // for debugging
      debugInfo.d = d;  // for debugging      
    }

    triangles.add(new Triangle(a, d, b));  // create a new triangle face
    e.isValid = false;
    va.outEdges.remove(vb);  // remove edge ab
    manifoldMask[a][d] = manifoldMask[d][b] = manifoldMask[b][a] = true;
    
    // Check if there is a front edge connecting va and vd
    FrontEdge e0 = null;
    if (vd.outEdges.containsKey(va)) {
      e0 = vd.outEdges.get(va);
      vd.outEdges.remove(va);
      e0.isValid = false;
    } else {
      e0 = new FrontEdge(a, d, b, normalADB);
      va.outEdges.put(vd, e0);
      curFront.add(e0);
    }
    
    // Check if current front touches a new front
    if (vd.groupID != -1 && !curFront.containGroupID(vd.groupID)) {
      curFront.mergeWithFront(vd, fronts);
    }

    // Check if there is a front edge connecting vd and vb
    FrontEdge e1 = null;
    if (vb.outEdges.containsKey(vd)) {
      e1 = vb.outEdges.get(vd);
      vb.outEdges.remove(vd);
      e1.isValid = false;
    } else {
      e1 = new FrontEdge(d, b, a, normalADB);
      vd.outEdges.put(vb, e1);
      curFront.add(e1);
    }

    if (!e0.isValid) {
      if (curFront.isInnerVertex(va)) va.isInner = true;
    }
    if (!e1.isValid) {
      if (curFront.isInnerVertex(vb)) vb.isInner = true;
    }
    if (!e0.isValid && !e1.isValid) {
      if (curFront.isInnerVertex(vd)) vd.isInner = true;
    }
    
    if (debugCH) debugInfo.numFaces++;
  }
  
  return true;
}

/*
 * Generate the convex hull for a given point array G whose size is nv. Assume
 * that all input points will be on the surface of the convex hull.
 */
ArrayList<Triangle> generateConvexHull(pt[] G, int nv) {
  assert nv >= 3;
  if (nv == 3) {
    Triangle t0 = new Triangle(0, 1, 2);
    Triangle t1 = new Triangle(0, 2, 1);
    ArrayList<Triangle> ch = new ArrayList<Triangle>();
    ch.add(t0);
    ch.add(t1);
    return ch;
  }
  ArrayList<Vertex> vertices;
  ArrayList<Triangle> triangles;
  Front[] fronts;
  boolean[][] manifoldMask;
  DebugCHInfo debugInfo;
  int k = 0;
  while (true) {
    vertices = convertToVertexList(G, nv, k);  // vertices may differ every time
    triangles = new ArrayList<Triangle>();
    fronts = new Front[1];
    fronts[0] = new Front();
    manifoldMask = new boolean[nv][nv];
    debugInfo = new DebugCHInfo();
    boolean successInit = findFirstTriangle(vertices, triangles, fronts[0],
                                            manifoldMask, debugInfo);
    if (!successInit) { k++; continue; }
    boolean successCH = generateConvexHullWithFronts(vertices, triangles,
                                                     fronts, manifoldMask,
                                                     debugInfo);
    if (!successCH) { k++; continue; }
    boolean successQT = passQualityTest(triangles, G, nv);
    if (!successQT) {
      k++;
      println("fail quality test!");
      if (k > 15) {
        exceptionHandler();
        break;
      }
      continue;
    } else break;
  }
  if (k >= 1) System.out.format("number of times tried = %d\n", k + 1);
  if (debugCH) {
    fill(blue); fronts[0].showEdges(G);
    fill(black); showInnerVertices(vertices, G);
    fill(#9AFF05); showManifoldEdges(manifoldMask, G, nv);  // light green
    fill(cyan); showTriangleNormals(triangles, G);
    if (debugInfo.a >= 0) { fill(#970EED, 160); show(G[debugInfo.a], 3); }  // purple
    if (debugInfo.b >= 0) { fill(#21C2FA, 160); show(G[debugInfo.b], 3); }  // light blue
    if (debugInfo.d >= 0) { fill(#FFF705, 160); show(G[debugInfo.d], 3); }  // light yellow
  }
  
  return triangles;
}

/*
 * Generate the convex hull for a given point array points whose size
 * is nGroup by nPointsPerGroup. Assume that all input points will be
 * on the surface of the convex hull.
 */
ArrayList<Triangle> generateConvexHull(pt[][] points,                // in      
                                       int nGroup,                   // in
                                       int nPointsPerGroup) {        // in
  ArrayList<Vertex> vertices;
  ArrayList<Triangle> triangles;
  Front[] fronts;
  boolean[][] manifoldMask;
  DebugCHInfo debugInfo;
  int nv = nGroup * nPointsPerGroup;
  int k = 0;
  pt[] G = convertTo1DArray(points, nGroup, nPointsPerGroup);
  while (true) {
    vertices = convertToVertexList(points, nGroup, nPointsPerGroup, k);
    triangles = new ArrayList<Triangle>();
    fronts = new Front[nGroup];
    manifoldMask = new boolean[nv][nv];
    debugInfo = new DebugCHInfo();
    initFronts(nGroup, nPointsPerGroup, vertices, fronts, manifoldMask,
               debugInfo);    
    boolean successCH = generateConvexHullWithFronts(vertices, triangles,
                                                     fronts, manifoldMask,
                                                     debugInfo);
    if (!successCH) { k++; continue; }
    boolean successQT = passQualityTest(triangles, G, nv);
    if (!successQT) {
      k++;
      println("fail quality test!");
      if (k > 15) {
        exceptionHandler();
        break;
      }
      continue;
    } else break;
  }
  if (k >= 1) System.out.format("number of times tried = %d\n", k + 1);
  if (debugCH) {
    fill(blue); fronts[0].showEdges(G);
    fill(#A08888);  // light grey
    for (int i = 1; i < nGroup; ++i) {
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