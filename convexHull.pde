/*********************************************************
 * Convex hull generation in a special case where all
 * vertices are on a sphere.
 *********************************************************/


import java.util.Queue;
import java.util.LinkedList;

int numFacesShown = 300;

boolean debugCH = false;
boolean debugDisturbance = true;
float disturbance = 0.001;

class DebugInfo {
  int a, b, d;
  int numSteps;
  DebugInfo() {
    a = b = d = -1;
    numSteps = 1;
  }
}

/*
 * Initialize a vertex list given a point (position) array. A vertex may
 * contain more info than a position.
 */
ArrayList<Vertex> convertToVertexList(pt[] G, int nv) {
  ArrayList<Vertex> vertices = new ArrayList<Vertex>();
  if (debugDisturbance) {
    for (int i = 0; i < nv; ++i) {
      float x = G[i].x + random(-1.0, 1.0) * disturbance;
      float y = G[i].y + random(-1.0, 1.0) * disturbance;
      float z = G[i].z + random(-1.0, 1.0) * disturbance;
      pt tmp = new pt(x, y, z);
      vertices.add(new Vertex(i, tmp));
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
 * Generate a convex hull from a list of vertices. Assume that all parameters
 * except the first one, are initialized by default values (empty/zeros, etc).
 * If no error or strange situation encountered, this function returns true,
 * otherwise returns false.
 */
boolean generateConvexHullFromVertices(ArrayList<Vertex> vertices,
                                       ArrayList<Triangle> triangles,
                                       Queue<TaggedEdge> frontEdges,
                                       boolean[][] manifoldMask,
                                       DebugInfo debugInfo) {
  // Find the first valid triangle
  triangles.add(findFirstTriangle(vertices, frontEdges, manifoldMask));
  debugInfo.numSteps = 1;
  
  // Expand the triangle mesh until the end
  while(frontEdges.size() > 0) {
    // stop after creating the first numFacesShown faces
    if (debugCH && debugInfo.numSteps >= numFacesShown) break;

    TaggedEdge e = frontEdges.poll();  // get and remove the first front edge
    if (e.tag != 0) continue;  // skip non-front edge
    int d = pivot(e, vertices, manifoldMask);

    if (d < 0) {
      //println("No hit found, number of steps = " + debugInfo.numSteps);
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
    else if (va.outEdges.containsKey(vd)) { // Is this possible? not possible when no 4 or more points on a plane?
      println("edge AD already exists, number of steps = " + debugInfo.numSteps);
      exceptionHandler();
      return false;
      //e0 = va.outEdges.get(vd);
      //va.outEdges.remove(vd);  // remove edge AD and don't add it back
    }

    // Check if there is a front edge connecting vd and vb
    TaggedEdge e1 = null;
    if (vb.outEdges.containsKey(vd)) {
      e1 = vb.outEdges.get(vd);
      vb.outEdges.remove(vd);
    }
    else if(vd.outEdges.containsKey(vb)) { // Is this possible?
      println("edge DB already exists, number of steps = " + debugInfo.numSteps);
      exceptionHandler();
      return false;
      //e1 = vd.outEdges.get(vb); 
      //vd.outEdges.remove(vb); // remove edge DB and don't add it back
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
  if (nv < 4) {
    println("Number of vertices less than 4!");
    return null;
  }
  
  ArrayList<Triangle> triangles;
  ArrayList<Vertex> vertices;
  Queue<TaggedEdge> frontEdges; // maybe using an array would be better, since edge DA, BD would be prev/next to AB, and hence no need to use outEdges for each vertex on the front
  boolean[][] manifoldMask;
  DebugInfo debugInfo;
  int k = 1;
  while (true) {
    vertices = convertToVertexList(G, nv);  // vertices may differ every time
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
    if (success && passQT) break;
    //if (k % 10 == 0) println("number of times tried = ", k);
    k++;
  }
  //System.out.format("number of times tried = %d, number of triangles = %d\n", k, triangles.size());

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