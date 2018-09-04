/*********************************************************
 * Geometric data structures used in this project.
 *********************************************************/


import java.util.HashMap;

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
}

class Vertex {
  int id;
  pt position;
  boolean isInner;
  HashMap<Vertex, TaggedEdge> outEdges;  // front or border edges going from this vertex  // seems that no need to use this if usinsk
  Vertex(){}
  Vertex(int _id, pt _position) {
    id = _id;
    position = _position;
    isInner = false;
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