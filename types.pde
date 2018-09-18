/******************************************************************************
 * Simple self-defined classes used in this project.
 ******************************************************************************/


import java.util.LinkedList;
import java.util.HashMap;
import java.util.HashSet;


class DebugCHInfo {
  int a, b, d;
  int numSteps;
  DebugCHInfo() {
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