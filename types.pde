/******************************************************************************
 * Some simple primitive types: disk, circle, ball, etc.
 ******************************************************************************/


import java.util.LinkedList;
import java.util.HashMap;
import java.util.HashSet;

class Triangle {
  int a, b, c;  // index of vertex
  Triangle(int a, int b, int c) {
    this.a = a;
    this.b = b;
    this.c = c;
  }

  void set(int a, int b, int c) {
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

  void set(int index, int value) {
    assert index >= 0 && index < 3;
    switch (index) {
      case 0:
        a = value;
        break;
      case 1:
        b = value;
        break;
      case 2:
        c = value;
        break;
    }
  }

  @Override
  public boolean equals(Object o) {
    if (o == this) return true;
    if (!(o instanceof Triangle)) return false;
    Triangle t = (Triangle)o;
    return a == t.a && b == t.b && c == t.c;
  }

  @Override
  public int hashCode() {
    return 31 * (31 * (31 + a) + b) + c;
  }
}

class Edge {
  int a, b;  // index of vertex
  Edge(int a, int b) {
    this.a = a;
    this.b = b;
  }

  @Override
  public boolean equals(Object o) {
    if (o == this) return true;
    if (!(o instanceof Edge)) return false;
    Edge e = (Edge)o;
    return a == e.a && b == e.b;
  }

  @Override
  public int hashCode() {
    return 31 * (31 + a) + b;
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

class Circle {
  pt c;
  vec n;
  float r;

  Circle(pt c, vec n, float r) {
    this.c = c;
    this.n = n;
    this.r = r;
  }

  void show() {
    showCircle(c, n, r);
  }
}

class Ball {
  pt c;
  float r;

  Ball(pt c, float r) {
    this.c = c;
    this.r = r;
  }

  Ball(float x, float y, float z, float r) {
    this.c = new pt(x, y, z);
    this.r = r;
  }

  public Ball copy() { return new Ball(c, r); }
  public Ball deepCopy() { return new Ball(c.c(), r); }

  boolean intersectBall(Ball other) {
    if (d(c, other.c) <= r + other.r) return true;
    return false;
  }

  float[] intersectLine(pt o, vec d) {
    vec co = V(c, o);
    float c2 = dot(d, d);
    float c1 = 2 * dot(co, d);
    float c0 = dot(co, co) - r * r;
    return solveQuadraticEquation(c2, c1, c0);
  }

  public Ball translate(vec v) { c.add(v); return this; }
  public Ball translate(float x, float y, float z) { c.add(x, y, z); return this; }
  public Ball rotateX(float rad) { c.rotX(rad); return this; }
  public Ball rotateY(float rad) { c.rotY(rad); return this; }
  public Ball rotateZ(float rad) { c.rotZ(rad); return this; }
  public Ball rotate(float rad, vec axis) { c.rot(rad, axis); return this; }
  public Ball rotate(float rad, vec axis, pt pointOnAxis) { c.rot(rad, axis, pointOnAxis); return this; }
  public Ball scale(float s) { c.mul(s); r *= s; return this; }
  public Ball scale(float s, pt f) { c.mul(s, f); r *= s; return this; }

  public Ball transform(Frame3 F) { return new Ball(F.toGlobalPoint(c), r * F.u.mag()); }  // Assume F is a similarity transform

  void show() {
    showBall(c, r);
  }
}

class vec2 {
  float x, y;
  vec2() {
    x = y = 0.0;
  }
  vec2(float x, float y) {
    this.x = x;
    this.y = y;
  }
  vec2 c() {
    return new vec2(x, y);
  }
  vec2 set(vec2 v) {
    x = v.x;
    y = v.y;
    return this;
  }
  vec2 set(float x, float y) {
    this.x = x;
    this.y = y;
    return this;
  }
  vec2 add(vec2 v) { x+=v.x; y+=v.y; return this; }
  vec2 add(float vx, float vy) { x+=vx; y+=vy; return this; }
  vec2 sub(vec2 v) { x-=v.x; y-=v.y; return this; }
  vec2 sub(float vx, float vy) { x-=vx; y-=vy; return this; }
  vec2 mul(float s) { x*=s; y*=s; return this; }
  vec2 mul(float sx, float sy) { x*=sx; y*=sy; return this; }
  float norm() {
    return sqrt(x * x + y * y);
  }

  String toString() {
    return "(" + x + ", " + y + ")";
  }
}

vec2 V(float a, float b) {
  return new vec2(a, b);
}

float dot(vec2 u, vec2 v) {
  return u.x * v.x + u.y * v.y;
}

vec2 disp(vec2 p, vec2 q) {
  return new vec2(q.x - p.x, q.y - p.y);
}

class EdgeCircle {
  pt a, b, c;
  float r;
  vec n, vi, vj;
  EdgeCircle() {}
  void init() {
    c = generateOnePointOnSphere(gSphereCenter, 0.5 * gSphereRadius);
    n = U(gSphereCenter, c);
    r = random(30, 50);
    a = generateOnePointInsideSphere(gSphereCenter, gSphereRadius);
    while (isZero(dot(V(c, a), n))) a = generateOnePointInsideSphere(gSphereCenter, gSphereRadius);
    b = generateOnePointInsideSphere(gSphereCenter, gSphereRadius);
    while (isZero(dot(V(c, b), n))) b = generateOnePointInsideSphere(gSphereCenter, gSphereRadius);
    vi = constructNormal(n);
    vj = N(n, vi);
  }

  void save(String file) {
    String[] lines = new String[7];
    int i = 0;
    lines[i++] = str(a.x) + "," + str(a.y) + "," + str(a.z);
    lines[i++] = str(b.x) + "," + str(b.y) + "," + str(b.z);
    lines[i++] = str(c.x) + "," + str(c.y) + "," + str(c.z);
    lines[i++] = str(r);
    lines[i++] = str(n.x) + "," + str(n.y) + "," + str(n.z);
    lines[i++] = str(vi.x) + "," + str(vi.y) + "," + str(vi.z);
    lines[i++] = str(vj.x) + "," + str(vj.y) + "," + str(vj.z);
    saveStrings(file, lines);
    return;
  }

  void load(String file) {
    String[] lines = loadStrings(file);
    int i = 0;
    float[] tmp = float(split(lines[i++], ","));
    a = new pt(tmp[0], tmp[1], tmp[2]);
    tmp = float(split(lines[i++], ","));
    b = new pt(tmp[0], tmp[1], tmp[2]);
    tmp = float(split(lines[i++], ","));
    c = new pt(tmp[0], tmp[1], tmp[2]);
    r = float(lines[i++]);
    tmp = float(split(lines[i++], ","));
    n = new vec(tmp[0], tmp[1], tmp[2]);
    tmp = float(split(lines[i++], ","));
    vi = new vec(tmp[0], tmp[1], tmp[2]);
    tmp = float(split(lines[i++], ","));
    vj = new vec(tmp[0], tmp[1], tmp[2]);
  }
}

class MinMaxI {
  public int min, max;  // inclusive, inclusive
  public MinMaxI(int pmin, int pmax) { min=pmin; max=pmax; }
  public MinMaxI copy() { return new MinMaxI(min, max); }
  public int size() { return max-min+1; }
  public String toString() { return min + ", " + max; }
}

class MinMaxF {
  public float min, max;
  public MinMaxF(float pmin, float pmax) { min = pmin; max = pmax; }
  public MinMaxF copy() { return new MinMaxF(min, max); }
  public float size() { return max-min; }
  public String toString() { return min + ", " + max; }
}

class Idx3 {
  public int i, j, k;
  public Idx3() { i = j = k = 0; }
  public Idx3(int pi, int pj, int pk) { i = pi; j = pj; k = pk; }
  public Idx3 c() { return new Idx3(i, j, k); }
  public boolean equals(Idx3 other) { return i==other.i && j==other.j && k==other.k; }
  public boolean isZero() { return i==0 && j==0 && k==0; }
  public Idx3 set(int pi, int pj, int pk) { i=pi; j=pj; k=pk; return this; }
  public Idx3 set(Idx3 v) { i=v.i; j=v.j; k=v.k; return this; }
  public Idx3 add(Idx3 v) { i+=v.i; j+=v.j; k+=v.k; return this; }
  public Idx3 add(int vi, int vj, int vk) { i+=vi; j+=vj; k+=vk; return this; }
  public Idx3 sub(Idx3 v) { i-=v.i; j-=v.j; k-=v.k; return this; }
  public Idx3 sub(int vi, int vj, int vk) { i-=vi; j-=vj; k-=vk; return this; }
  public Idx3 mul(int s) { i*=s; j*=s; k*=s; return this; }
  public Idx3 mul(int si, int sj, int sk) { i*=si; j*=sj; k*=sk; return this; }
  public Idx3 mul(Idx3 s) { i*=s.i; j*=s.j; k*=s.k; return this; }
  public Idx3 div(int s) { i/=s; j/=s; k/=s; return this; }
  public Idx3 div(int si, int sj, int sk) { i/=si; j/=sj; k/=sk; return this; }
  public Idx3 div(Idx3 s) { i/=s.i; j/=s.j; k/=s.k; return this; }
  public String toString() { return i + "," + j + "," + k; }
}

Idx3 I(int i, int j, int k) { return new Idx3(i, j, k); }

enum Idx3Order { IJK, IKJ, JIK, JKI, KIJ, KJI }

class Idx3Iterator {
  private Idx3 min, max, stride, current;

  public Idx3Iterator(Idx3 pMin, Idx3 pMax) { this(pMin, pMax, I(1,1,1)); }
  public Idx3Iterator(Idx3 pMin, Idx3 pMax, Idx3 pStride) {
    min = pMin.c();
    max = pMax.c();
    stride = pStride.c();
    current = min.c();
  }

  public Idx3Iterator set(Idx3 pMin, Idx3 pMax, Idx3 pStride) {
    min.set(pMin);
    max.set(pMax);
    stride.set(pStride);
    current.set(pMin);
    return this;
  }

  public boolean init() {
    current = min.c();
    return current.i <= max.i && current.j <= max.j && current.k <= max.k;
  }

  public boolean next() {
    return nextKJI();
  }

  public boolean nextIJK() {
    current.i += stride.i;
    if (current.i > max.i) {
      current.i = min.i;
      current.j += stride.j;
      if (current.j > max.j) {
        current.j = min.j;
        current.k += stride.k;
      }
    }
    return current.k <= max.k;
  }

  public boolean nextIKJ() {
    current.i += stride.i;
    if (current.i > max.i) {
      current.i = min.i;
      current.k += stride.k;
      if (current.k > max.k) {
        current.k = min.k;
        current.j += stride.j;
      }
    }
    return current.j <= max.j;
  }

  public boolean nextJIK() {
    current.j += stride.j;
    if (current.j > max.j) {
      current.j = min.j;
      current.i += stride.i;
      if (current.i > max.i) {
        current.i = min.i;
        current.k += stride.k;
      }
    }
    return current.k <= max.k;
  }

  public boolean nextJKI() {
    current.j += stride.j;
    if (current.j > max.j) {
      current.j = min.j;
      current.k += stride.k;
      if (current.k > max.k) {
        current.k = min.k;
        current.i += stride.i;
      }
    }
    return current.i <= max.i;
  }

  public boolean nextKIJ() {
    current.k += stride.k;
    if (current.k > max.k) {
      current.k = min.k;
      current.i += stride.i;
      if (current.i > max.i) {
        current.i = min.i;
        current.j += stride.j;
      }
    }
    return current.j <= max.j;
  }

  public boolean nextKJI() {
    current.k += stride.k;
    if (current.k > max.k) {
      current.k = min.k;
      current.j += stride.j;
      if (current.j > max.j) {
        current.j = min.j;
        current.i += stride.i;
      }
    }
    return current.i <= max.i;
  }

  public Idx3 current() { return current.c(); }
}