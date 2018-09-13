/******************************************************************************
 * Self-defined classes used in this project.
 ******************************************************************************/


import java.util.LinkedList;
import java.util.HashMap;
import java.util.HashSet;


class DebugInfo {
  int a, b, d;
  int numSteps;
  DebugInfo() {
    a = b = d = -1;
    numSteps = 1;
  }
}

class Triangle {
  int a, b, c;  // index of vertex
  Triangle(int a, int b, int c) {
    this.a = a;
    this.b = b;
    this.c = c;
  }

  int get(int index) {
    assert index >= 0 && index < 3;
    switch (index) {
      case 0:
        return a;
      case 1:
        return b;
      case 2:
        return c;
    }
    return -1;
  }
}

class Edge {
  int a, b;  // index of vertex
  Edge(int a, int b) {
    this.a = a;
    this.b = b;
  }
}

class FrontEdge {
  int a, b;  // the id of two ends of the edge
  int c; // the id of opposite vertex w.r.t. front edge [a, b]
  boolean isValid;  // true:front, false: non-front 
  vec N;  // the outward normal defined by vertices A, B, C
  FrontEdge(int a, int b, int c, boolean isValid, vec N) {
    this.a = a;
    this.b = b;
    this.c = c;
    this.isValid = isValid;
    this.N = N;
  }
  FrontEdge(int a, int b, int c, vec N) {
    this.a = a;
    this.b = b;
    this.c = c;
    this.isValid = true;
    this.N = N;
  }
}

class Front {
  LinkedList<FrontEdge> edges;
  HashSet<Integer> groupIDs;
  Front() {
    edges = new LinkedList<FrontEdge>();
    groupIDs = new HashSet<Integer>();
  }
  
  Front(LinkedList<FrontEdge> edges) {
    this.edges = edges;
    groupIDs = new HashSet<Integer>();
  }
  
  Front(LinkedList<FrontEdge> edges, HashSet<Integer> groupIDs) {
    this.edges = edges;
    this.groupIDs = groupIDs;
  }
  
  LinkedList<FrontEdge> getEdges() {
    return edges;
  }
  
  HashSet<Integer> getGroupIDs() {
    return groupIDs;
  }
  
  int size() {
    return edges.size();
  }
  
  boolean empty() {
    return edges.size() == 0;
  }
  
  FrontEdge poll() {
    return edges.poll();
  }
  
  boolean add(FrontEdge e) {
    return edges.add(e);
  }
  
  void invalidate(FrontEdge e) {
    e.isValid = false;
  }
  
  FrontEdge get(int index) {
    return edges.get(index);
  }
  
  int indexOf(Vertex v) {
    int n = edges.size();
    int i = 0;
    for (; i < n; ++i) {
      if (edges.get(i).a == v.id) break;
    }
    return i;
  }
  
  void shiftK(int k) {
    while (k > 0) {
      FrontEdge e = edges.poll();
      edges.add(e);
      k--;
    }
    return;
  }
  
  boolean containGroupID(Integer gid) {
    return groupIDs.contains(gid);
  }
  
  void mergeWithFront(Front front) {
    edges.addAll(front.getEdges());
    front.getEdges().clear();
    for (Integer gid : front.getGroupIDs()) {
      groupIDs.add(gid);
    }
    front.getGroupIDs().clear();
  }
  
  void mergeWithFront(Vertex v, Front[] fronts) {
    int groupID = v.groupID;
    Front anotherFront = fronts[groupID];
    int k = anotherFront.indexOf(v);
    anotherFront.shiftK(k);
    this.mergeWithFront(anotherFront);
  }
  
  boolean isInnerVertex(Vertex v) {
    for (FrontEdge edge : edges) {
      if (edge.isValid && (edge.a == v.id || edge.b == v.id)) return false;
    }
    return true;
  }
  
  void showEdges(pt[] G) {
    for (FrontEdge edge : edges) {
      if (!edge.isValid) continue;
      pt A, B;
      A = G[edge.a];
      B = G[edge.b];
      arrow(A, V(A, B), 2);
    }
  }
}

class Vertex {
  int id;
  pt position;
  boolean isInner;
  int groupID;
  HashMap<Vertex, FrontEdge> outEdges;  // front edges going from this vertex
  Vertex(){}
  Vertex(int id, pt position) {
    this.id = id;
    this.position = position;
    isInner = false;
    groupID = -1;
    outEdges = new HashMap<Vertex, FrontEdge>();
  }
  Vertex(int id, pt position, int groupID) {
    this.id = id;
    this.position = position;
    isInner = false;
    this.groupID = groupID;
    outEdges = new HashMap<Vertex, FrontEdge>();
  }
}

class Disk {
  pt c;
  vec d;
  float r;
  Disk(pt c, vec d, float r) {
    this.c = c;
    this.d = d;
    this.r = r;
  }
}

class RingSet {
  pt c;
  float r;
  int nc, np;
  boolean sameRadius;
  
  pt[] contacts;
  float[] radii;
  vec[] initDirs;
  
  pt[][] points;
  pt[] centers;
  
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
  
  pt[][] get2DPointArray() {
    return points;
  }
  
  int getNumGroups() {
    return nc;
  }
  
  int getNumPointsPerGroup() {
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
  
  void showGroups() {
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
  
  void savePointGroups(String file) {
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

  void loadPointGroups(String file) {
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

/*
 * Ball class.
 */
class Ball {
  pt c;
  float r;
  Ball(pt c, float r) {
    this.c = c;
    this.r = r;
  }
  void showBall() {
    show(c, r);
  }
}

class vec2 {
  float x, y;
  vec2() {
    x = y = 0.0;
  }
  vec2(float x_, float y_) {
    x = x_;
    y = y_;
  }
  vec2 set(vec2 v) {
    x = v.x;
    y = v.y;
    return this;
  }
  vec2 set(float x_, float y_) {
    x = x_;
    y = y_;
    return this;
  }
}