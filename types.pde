/*********************************************************
 * Geometric data structures used in this project.
 *********************************************************/

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
  Triangle(int _a, int _b, int _c) { a = _a; b = _b; c = _c;}
}

class Edge {
  int a, b;  // index of vertex
  Edge(int _a, int _b) { a = _a; b = _b;}
}

class FrontEdge {
  int a, b;  // the id of two ends of the edge
  int c; // the id of opposite vertex w.r.t. edge e, it only makes senses when e is a front edge
  boolean isValid;  // true:front, false: non-front 
  vec N;  // the outward normal defined by vertices A, B, C
  FrontEdge(int _a, int _b, int _c, boolean _isValid, vec _N) {
    a = _a;
    b = _b;
    c = _c;
    isValid = _isValid;
    N = _N;
  }
  FrontEdge(int _a, int _b, int _c, vec _N) {
    a = _a;
    b = _b;
    c = _c;
    isValid = true;
    N = _N;
  }
}

class Front {
  LinkedList<FrontEdge> edges;
  HashSet<Integer> groupIDs;
  Front() {
    edges = new LinkedList<FrontEdge>();
    groupIDs = new HashSet<Integer>();
  }
  
  Front(LinkedList<FrontEdge> _edges) {
    edges = _edges;
    groupIDs = new HashSet<Integer>();
  }
  
  Front(LinkedList<FrontEdge> _edges, HashSet<Integer> _groupIDs) {
    edges = _edges;
    groupIDs = _groupIDs;
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
  HashMap<Vertex, FrontEdge> outEdges;  // front or border edges going from this vertex  // seems that no need to use this if usinsk
  Vertex(){}
  Vertex(int _id, pt _position) {
    id = _id;
    position = _position;
    isInner = false;
    groupID = -1;
    outEdges = new HashMap<Vertex, FrontEdge>();
  }
  Vertex(int _id, pt _position, int _groupID) {
    id = _id;
    position = _position;
    isInner = false;
    groupID = _groupID;
    outEdges = new HashMap<Vertex, FrontEdge>();
  }
}

class Disk {
  pt c;
  vec d;
  float r;
  Disk(pt c_, vec d_, float r_) {
    c = c_;
    d = d_;
    r = r_;
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
  
  RingSet(pt _c, float _r) {
    c = _c;
    r = _r;
    sameRadius = false;
  }
  
  RingSet(pt _c, float _r, int _nc, int _np) {
    c = _c;
    r = _r;
    nc = _nc;
    np = _np;
    sameRadius = false;
  }
  
  RingSet(pt _c, float _r, int _nc, int _np, float rMax) {
    c = _c;
    r = _r;
    nc = _nc;
    np = _np;
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
      points = generatePointsForCircles(contacts, curRadii, c, r, initDirs, nc, np, centers);
    } else {
      float curRadius = radii[0] * attenuation;
      points = generatePointsForCircles(contacts, curRadius, c, r, initDirs, nc, np, centers);
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
      lines[i++] = str(contacts[j].x) + "," + str(contacts[j].y) + "," + str(contacts[j].z);
      lines[i++] = str(radii[j]);
      lines[i++] = str(initDirs[j].x) + "," + str(initDirs[j].y) + "," + str(initDirs[j].z);
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