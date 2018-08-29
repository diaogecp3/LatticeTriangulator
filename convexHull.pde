import java.util.Queue;
import java.util.LinkedList;

int numVerticesIgnored = 0;
int numFacesShown = 150;

boolean debugCH = true;
boolean showCH = true;

/*
 * Initialize a vertex list given a point (position) array. A vertex may
 * contain more info than a position.
 */
ArrayList<Vertex> convertToVertexList(pt[] G, int nv) {
  ArrayList<Vertex> vertices = new ArrayList<Vertex>();
  for (int i = 0; i < nv; ++i) {
    vertices.add(new Vertex(i, G[i]));
  }
  return vertices;
}


/*
 * Check if a vertex is inner given a list of front edges.
 */
boolean isInnerVertex(Queue<TaggedEdge> frontEdges, Vertex v) {
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
                           Queue<TaggedEdge> frontEdges,
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
int pivot(TaggedEdge e, ArrayList<Vertex> vertices, boolean[][] manifoldMask) {
  int a = e.e.a, b = e.e.b, c = e.c;
  pt A = vertices.get(a).position, B = vertices.get(b).position, C = vertices.get(c).position;
  vec AC = V(A, C), AB = V(A, B);
  vec normalizedAB = U(AB);
  vec dir0 = U(M(AC, V(dot(AC, normalizedAB), normalizedAB)));  // normalize(AC - (AC dot (normalize(AB)) times (normalize(AB)))

  float minCosTheta = 1.5f;
  int d = -1;
  int nv = vertices.size();
  for (int i = 0; i < nv; ++i) {
    if (i == a || i == b || i == c) continue;
    if (vertices.get(i).isInner) { // ignore inner vertices
      numVerticesIgnored++;
      continue;  
    }
    if (manifoldMask[a][i] || manifoldMask[i][b]) continue;
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
 * Generate the convex hull for a given point array G whose size is nv. Assume
 * that all input points will be on the surface of the convex hull.
 */
ArrayList<Triangle> generateConvexHull(pt[] G, int nv) {
  if (nv < 4) {
    print("Number of vertices less than 4!");
    return null;
  }
  ArrayList<Triangle> triangles = new ArrayList<Triangle>();
  Queue<TaggedEdge> frontEdges = new LinkedList<TaggedEdge>();
  ArrayList<Vertex> vertices = convertToVertexList(G, nv);
  boolean[][] manifoldMask = new boolean[nv][nv];  // initialize as false matrix

  // Find the first valid triangle
  triangles.add(findFirstTriangle(vertices, frontEdges, manifoldMask));

  int k = 1;
  int ia = -1, ib = -1, id = -1;
  // Expand the triangle mesh until the end
  while(frontEdges.size() > 0) {
    if (debugCH && k >= numFacesShown) break;  // stop when creating the first k faces 

    TaggedEdge e = frontEdges.poll();  // get and remove the first front edge
    if (e.tag != 0) continue;  // skip non-front edge
    int d = pivot(e, vertices, manifoldMask);

    if (debugCH) {
      ia = e.e.a;  // for debugging
      ib = e.e.b;  // for debugging
      id = d;  // for debugging
    }

    if (d < 0) {
      print("No hit found, number of steps = " + k);
      exceptionHandler();
      if (debugCH) break;
      else exit();
      numVerticesIgnored = 0;
      return triangles;
    }

    int a = e.e.a, b = e.e.b;
    Vertex va = vertices.get(a);
    Vertex vb = vertices.get(b);
    Vertex vd = vertices.get(d);

    triangles.add(new Triangle(a, d, b));  // create a new triangle facet
    e.tag = 2;  // the front edge becomes an inner edge  // TODO: seems that no need to tag edge or vertex? if use manifold mask
    manifoldMask[a][d] = manifoldMask[d][b] = manifoldMask[b][a] = true;

    va.outEdges.remove(vb);  // remove edge ab

    // Check if there is a front edge connecting va and vd
    TaggedEdge e0 = null;
    if (vd.outEdges.containsKey(va)) {
      e0 = vd.outEdges.get(va); 
      vd.outEdges.remove(va);
    }
    else if (va.outEdges.containsKey(vd)) { // is this possible? not possible when no 4 points on a plane?
      println("edge AD already exists, number of steps = " + k);
      exceptionHandler();
      if (debugCH) break;
      else exit();
      e0 = va.outEdges.get(vd);
      va.outEdges.remove(vd);  // remove edge ad and then add it back?
    }

    // Check if there is a front edge connecting vd and vb
    TaggedEdge e1 = null;
    if (vb.outEdges.containsKey(vd)) {
      e1 = vb.outEdges.get(vd);
      vb.outEdges.remove(vd);
    }
    else if(vd.outEdges.containsKey(vb)) { // is this possible?
      println("edge DB already exists, number of steps = " + k);
      exceptionHandler();
      if (debugCH) break;
      else exit();
      e1 = vd.outEdges.get(vb); 
      vd.outEdges.remove(vb);
    }

    if (e0 == null) { // mark e0 as front edge and push it to the front list
      e0 = new TaggedEdge(a, d, b, 0);  // create an edge from a to d
      va.outEdges.put(vd, e0);  // update the adjacent edges starting from a
      frontEdges.add(e0);
    } else {  // e0 already exists
      e0.tag = 2;  // remove this edge from front edges
    }

    if (e1 == null) { // mark e1 as front edge and push it to the front list
      e1 = new TaggedEdge(d, b, a, 0);  // create an edge from d to b
      vd.outEdges.put(vb, e1);  // update the adjacent edges starting from d
      frontEdges.add(e1);
    } else {  // e1 already exists
      e1.tag = 2;  // remove this edge from front edges
    }

    /*
     * Check if vd becomes an inner vertex when e0 and e1 are inner.
     * Check if va becomes an inner vertex when e0 is inner.
     * Check if vb becomes an inner vertex when e1 is inner.
     */
    if (e0.tag == 2 && e1.tag == 2) {
      if (isInnerVertex(frontEdges, vd)) vd.isInner = true;
    }
    if (e0.tag == 2) {
      if (isInnerVertex(frontEdges, va)) va.isInner = true;
    }
    if (e1.tag == 2) {
      if (isInnerVertex(frontEdges, vb)) vb.isInner = true;
    }
    
    k++;
  } // end while

  if (debugCH) {
    fill(blue); showFrontEdges(frontEdges, G);
    fill(black); showInnerVertices(vertices);
    fill(#9AFF05); showManifoldEdges(manifoldMask, G, nv);
    if (ia >= 0) { fill(#970EED, 100); show(G[ia], 3); }
    if (ib >= 0) { fill(#21C2FA, 100); show(G[ib], 3); }
    if (id >= 0) { fill(#FFF705, 100); show(G[id], 3); }
  }
  
  numVerticesIgnored = 0; // reset to 0
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
void showFrontEdges(Queue<TaggedEdge> frontEdges, pt[] G) {
  for (TaggedEdge edge : frontEdges) {
    if (edge.tag != 0) continue;
    pt A, B;
    A = G[edge.e.a];
    B = G[edge.e.b];
    arrow(A, V(A, B), 1);
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
void showInnerVertices(ArrayList<Vertex> vertices) {
  int n = vertices.size();
  for (int i = 0; i < n; ++i) {
    if (vertices.get(i).isInner) {
      show(vertices.get(i).position, 2);
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