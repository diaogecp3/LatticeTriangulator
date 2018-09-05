/*********************************************************
 * Geometric data structures used in this project.
 *********************************************************/


import java.util.HashMap;
import java.util.HashSet;

class Triangle {
  int a, b, c;  // index of vertex
  Triangle(int _a, int _b, int _c) { a = _a; b = _b; c = _c;}
}

class Edge {
  int a, b;  // index of vertex
  Edge(int _a, int _b) { a = _a; b = _b;}
}

class TaggedEdge {
  Edge e;
  int c; // the id of opposite vertex w.r.t. edge e, it only makes senses when e is a front or border edge
  int tag;  // 0:front, 1:border, 2:inner
  TaggedEdge(int _a, int _b, int _c, int _tag) {
    e = new Edge(_a,_b);
    c = _c;
    tag = _tag;
  }
  TaggedEdge(int _a, int _b, int _c) {
    e = new Edge(_a, _b);
    c = _c;
    tag = 0;
  }
}

class Front {
  LinkedList<TaggedEdge> edges;
  HashSet<Integer> groupIDs;
  Front() {
    edges = new LinkedList<TaggedEdge>();
    groupIDs = new HashSet<Integer>();
  }
  Front(LinkedList<TaggedEdge> _edges, HashSet<Integer> _groupIDs) {
    edges = _edges;
    groupIDs = _groupIDs;
  }
  
  LinkedList<TaggedEdge> getEdges() {
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
  
  TaggedEdge poll() {
    return edges.poll();
  }
  
  boolean add(TaggedEdge e) {
    return edges.add(e);
  }
  
  void invalidate(TaggedEdge e) {
    e.tag = 2;
  }
  
  boolean remove(TaggedEdge e) {
    return edges.remove(e);
  }
  
  TaggedEdge get(int index) {
    return edges.get(index);
  }
  
  int indexOf(Vertex v) {
    int n = edges.size();
    int i = 0;
    for (; i < n; ++i) {
      if (edges.get(i).e.a == v.id) break;
    }
    return i;
  }
  
  void shiftK(int k) {
    while (k > 0) {
      TaggedEdge e = edges.poll();
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
    for (TaggedEdge edge : edges) {
      if (edge.tag == 0 && (edge.e.a == v.id || edge.e.b == v.id)) return false;
    }
    return true;
  }
  
  void showEdges(pt[] G) {
    for (TaggedEdge edge : edges) {
      if (edge.tag != 0) continue;
      pt A, B;
      A = G[edge.e.a];
      B = G[edge.e.b];
      arrow(A, V(A, B), 2);
    }
  }
}

class Vertex {
  int id;
  pt position;
  boolean isInner;
  int groupID;
  HashMap<Vertex, TaggedEdge> outEdges;  // front or border edges going from this vertex  // seems that no need to use this if usinsk
  Vertex(){}
  Vertex(int _id, pt _position) {
    id = _id;
    position = _position;
    isInner = false;
    groupID = -1;
    outEdges = new HashMap<Vertex, TaggedEdge>();
  }
  Vertex(int _id, pt _position, int _groupID) {
    id = _id;
    position = _position;
    isInner = false;
    groupID = _groupID;
    outEdges = new HashMap<Vertex, TaggedEdge>();
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