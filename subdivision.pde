/******************************************************************************
 * Subdivision of the polygonal mesh that approximates the convex hull of
 * co-spherical circles.
 * The polygonal mesh contains triangles and quads, and has polygonal border
 * loops for each circle.
 * The polygonal mesh becomes a triangle mesh with polygonal border loops after
 * subdivision.
 ******************************************************************************/

/* Subdivide type for a triangle. */
enum SubdivideTypeTriangle {LOOP}

/* Subdivide type for a quad. */
enum SubdivideTypeQuad {FAN, DIAMOND}

/* Projection type. */
enum ProjectType {RAY, SPHERE_TRACING}

/*
 * The vertex class used in subdivision.
 */
class VertexSub {
  pt position;
  int circleID;  // the ID of the circle where it lies, -1: no such circle
  vec normal;  // the face normal of the vertex if it doesn't belong to a circle, or the edge normal of the vertex if it belongs to a circle/disk
  VertexSub(pt p, int cID, vec n) {
    this.position = p;
    this.circleID = cID;
    this.normal = n;
  }
}


/*
 * Project point p, which has normal n, to circle (c, r). Point p is on the
 * plane of the circle.
 */
pt projectPointOnCircle(pt p, vec n, pt c, float r) {
  if (n == null) {
    return P(c, r, U(c, p));  // radial projection
  } else {
    vec cp = V(c, p);
    if (dot(cp, n) > 0) return P(c, r, cp.normalize());  // radial projection
    else {
      vec d = reflect(cp, n);
      return P(c, r, d.normalize());  // reflected radial projection
    }
  }
}


/*
 * The triangle class used in subdivision.
 */
class TriangleSub {
  VertexSub[] vertices;  // of length 3
  color col = 0;

  TriangleSub(VertexSub[] vertices) {
    this.vertices = vertices;
  }

  TriangleSub(VertexSub va, VertexSub vb, VertexSub vc, color c) {
    this.vertices = new VertexSub[] {va, vb, vc};
    this.col = c;
  }

  void setColor(color c) {
    col = c;
  }

  /* Subdivie a triangle into 4 triangles. */
  ArrayList<TriangleSub> subdivide(SubdivideTypeTriangle subTypeTri, Circle[] circles, ArrayList<pt>[] borders, boolean projectOnCircle) {
    ArrayList<TriangleSub> tris = new ArrayList<TriangleSub>();
    vec normal = U(N(vertices[0].position, vertices[1].position, vertices[2].position));
    if (subTypeTri == SubdivideTypeTriangle.LOOP) {
      VertexSub[] mids = new VertexSub[3];
      for (int i = 0; i < 3; ++i) {
        VertexSub va = vertices[i];
        VertexSub vb = vertices[(i + 1) % 3];
        if (va.circleID != -1 && va.circleID == vb.circleID) {  // on the same circle
          pt c = circles[va.circleID].c;
          float r = circles[va.circleID].r;
          vec n = circles[va.circleID].n;
          pt p = P(va.position, vb.position);  // p is the midpoint of edge (va, vb)

          vec edgeNormal = U(N(V(va.position, vb.position), n));
          if (projectOnCircle) {
            p = projectPointOnCircle(p, null, c, r);  // radial projection
            // p = projectPointOnCircle(p, edgeNormal, c, r);  // this doesn't work
          }

          mids[i] = new VertexSub(p, va.circleID, edgeNormal);
          borders[va.circleID].add(p);  // this border is not sorted after insertion
        } else {  // not on the same circle
          pt p = P(va.position, vb.position);
          mids[i] = new VertexSub(p, -1, normal);
        }
      }
      tris.add(new TriangleSub(vertices[0], mids[0], mids[2], this.col));
      tris.add(new TriangleSub(vertices[1], mids[1], mids[0], this.col));
      tris.add(new TriangleSub(vertices[2], mids[2], mids[1], this.col));
      tris.add(new TriangleSub(mids[0], mids[1], mids[2], this.col));
    } else {
      println("Please use a correct subdivision type for triangles.");
    }
    assert tris.size() == 4;
    return tris;
  }
}

/*
 * The quad class used in subdivision.
 */
class QuadSub {
  VertexSub[] vertices;  // of length 4
  color col = 0;

  QuadSub(VertexSub[] vertices) {
    this.vertices = vertices;
  }

  QuadSub(VertexSub va, VertexSub vb, VertexSub vc, VertexSub vd, color c) {
    this.vertices = new VertexSub[] {va, vb, vc, vd};
    this.col = c;
  }

  void setColor(color c) {
    col = c;
  }

  private VertexSub[] midVertices(Circle[] circles, ArrayList<pt>[] borders, boolean projectOnCircle) {
    VertexSub[] mids = new VertexSub[5];  // 4 edge mid + 1 face mid
    vec normal = U(N(vertices[0].position, vertices[1].position, vertices[2].position));
    for (int i = 0; i < 4; ++i) {
      VertexSub va = vertices[i];
      VertexSub vb = vertices[(i + 1) % 4];
      if (va.circleID != -1 && va.circleID == vb.circleID) {  // on the same circle
        pt c = circles[va.circleID].c;
        float r = circles[va.circleID].r;
        vec n = circles[va.circleID].n;
        pt p = P(va.position, vb.position);  // p is the midpoint of edge (va, vb)

        vec edgeNormal = U(N(V(va.position, vb.position), n));
        if (projectOnCircle) {
          p = projectPointOnCircle(p, null, c, r);  // radial projection
          // p = projectPointOnCircle(p, edgeNormal, c, r);  // this doesn't work
        }

        mids[i] = new VertexSub(p, va.circleID, edgeNormal);
        borders[va.circleID].add(p);  // this border is not sorted after insertion
      } else {  // not on the same circle
        pt p = P(va.position, vb.position);
        mids[i] = new VertexSub(p, -1, normal);
      }
    }

    {
      pt p = P(vertices[0].position, vertices[1].position, vertices[2].position, vertices[3].position);
      mids[4] = new VertexSub(p, -1, normal);
    }

    return mids;
  }

  /* Subdivide a quad into 8 triangles. */
  ArrayList<TriangleSub> subdivide(SubdivideTypeQuad subTypeQuad, Circle[] circles, ArrayList<pt>[] borders, boolean projectOnCircle) {
    ArrayList<TriangleSub> tris = new ArrayList<TriangleSub>();
    if (subTypeQuad == SubdivideTypeQuad.FAN) {
      VertexSub[] mids = midVertices(circles, borders, projectOnCircle);
      tris.add(new TriangleSub(vertices[0], mids[0], mids[4], col));
      tris.add(new TriangleSub(vertices[1], mids[1], mids[4], col));
      tris.add(new TriangleSub(vertices[2], mids[2], mids[4], col));
      tris.add(new TriangleSub(vertices[3], mids[3], mids[4], col));
      tris.add(new TriangleSub(mids[0], vertices[1], mids[4], col));
      tris.add(new TriangleSub(mids[1], vertices[2], mids[4], col));
      tris.add(new TriangleSub(mids[2], vertices[3], mids[4], col));
      tris.add(new TriangleSub(mids[3], vertices[0], mids[4], col));
    } else if (subTypeQuad == SubdivideTypeQuad.DIAMOND) {
      VertexSub[] mids = midVertices(circles, borders, projectOnCircle);
      tris.add(new TriangleSub(vertices[0], mids[0], mids[3], col));
      tris.add(new TriangleSub(vertices[1], mids[1], mids[0], col));
      tris.add(new TriangleSub(vertices[2], mids[2], mids[1], col));
      tris.add(new TriangleSub(vertices[3], mids[3], mids[2], col));
      tris.add(new TriangleSub(mids[0], mids[1], mids[4], col));
      tris.add(new TriangleSub(mids[1], mids[2], mids[4], col));
      tris.add(new TriangleSub(mids[2], mids[3], mids[4], col));
      tris.add(new TriangleSub(mids[3], mids[0], mids[4], col));
    } else {
      println("Please use a correct subdivision type for quads.");
    }
    assert tris.size() == 8;
    return tris;
  }
}

/*
 * Bordered triangle-quad mesh used to approximate the convex hull of circles.
 */
class BorderedTriQuadMesh {
  ArrayList<TriangleSub> triangles;
  ArrayList<QuadSub> quads;
  Circle[] circles;
  ArrayList<pt>[] borders;

  BorderedTriQuadMesh(ArrayList<TriangleSub> triangles, ArrayList<QuadSub> quads, Circle[] circles, ArrayList<pt>[] borders) {
    this.triangles = triangles;
    this.quads = quads;
    this.circles = circles;
    this.borders = borders;
  }

  void setColors(color ct, color cq) {
    for (TriangleSub tri : triangles) tri.setColor(ct);
    for (QuadSub quad : quads) quad.setColor(cq);
  }

  BorderedTriQuadMesh subdivide(SubdivideTypeTriangle subTypeTri, SubdivideTypeQuad subTypeQuad, boolean projectOnCircle) {
    ArrayList<TriangleSub> newTriangles = new ArrayList<TriangleSub>();
    for (TriangleSub tri : triangles) {
      ArrayList<TriangleSub> tris = tri.subdivide(subTypeTri, circles, borders, projectOnCircle);
      newTriangles.addAll(tris);
    }

    if (quads != null) {
      for (QuadSub quad : quads) {
        ArrayList<TriangleSub> tris = quad.subdivide(subTypeQuad, circles, borders, projectOnCircle);
        newTriangles.addAll(tris);
      }
    }

    /* Sort vertices on each border. */
    for (int i = 0; i < borders.length; ++i) {
      ArrayList<pt> newBorder = sortBorderLoop(circles[i], borders[i]);
      borders[i].clear();
      borders[i].addAll(newBorder);
    }

    return new BorderedTriQuadMesh(newTriangles, null, circles, borders);
  }

  /* The hub may be partially outside the mesh. */
  private void projectOnHubRay(Hub hub) {
    for (TriangleSub tri : triangles) {  // no need to consider quads since quads only contain vertices on circles
      for (VertexSub v : tri.vertices) {
        pt p = v.position;
        if (hub.distanceFrom(p) < 0.0001) continue;
        assert v.normal != null;
        vec d = U(p, hub.ball.c);
        if (dot(d, v.normal) < 0) {  // hub center cannot see vertex v
          Float t = hub.closestIntersectionWithRay(p, d);
          if (t != null) p.set(P(p, t, d));
        } else {  // hub center can see vertex v
          // println("hub center can see vertex v");
          /* Redirect the ray direction. */
          d.rev();
          vec r = reflect(d, v.normal);
          r.normalize();

          // Float t = hub.furthestIntersectionWithRay(p, v.normal);  // use vertex normal
          Float t = hub.furthestIntersectionWithRay(hub.ball.c, r);
          {  // debug
            // noStroke();
            // fill(magenta);
            // arrow(p, V(20, v.normal), 1);
          }
          // println("t =", t);
          // if (t != null) p.set(P(p, t, v.normal));
          if (t != null) p.set(P(hub.ball.c, t, r));
        }
      }
    }
  }

  /* Assume that the hub is completely inside the mesh. */
  private void projectOnHubSphereTrace(Hub hub) {
    int maxIter = 32;
    for (TriangleSub tri : triangles) {
      for (VertexSub v : tri.vertices) {
        pt o = v.position;
        if (hub.distanceFrom(o) < 0.0001) continue;
        vec d = U(o, hub.ball.c);
        pt p = P(o);
        for (int j = 0; j < maxIter; ++j) {
          float dist = hub.blendedDistanceFrom(p);
          if (dist < 0.0001) break;
          if (Float.isInfinite(dist)) {
            dist = 0.0;
            break;
          }
          p.add(dist, d);
        }
        o.set(p);
      }
    }
  }

  void projectOnHub(Hub hub, ProjectType option) {
    switch (option) {
      case SPHERE_TRACING:
        projectOnHubSphereTrace(hub);
        break;
      default:
        projectOnHubRay(hub);
    }
  }

  /* Only use triangle faces. */
  TriangleMesh toTriangleMesh() {
    assert quads == null || quads.size() == 0;

    ArrayList<pt> poss = new ArrayList<pt>();
    ArrayList<Triangle> tris = new ArrayList<Triangle>();

    for (TriangleSub tri : triangles) {
      int[] vids = new int[] {-1, -1, -1};
      for (int i = 0; i < 3; ++i) {
        pt p = tri.vertices[i].position;
        int j = 0;
        for (; j < poss.size(); ++j) {
          if (d(p, poss.get(j)) < 0.0001) {
            vids[i] = j;
            break;
          }
        }
        if (j == poss.size()) {  // find a new vertex
          vids[i] = j;
          poss.add(p);
        }
      }
      assert vids[0] != -1 && vids[1] != -1 && vids[2] != -1;
      tris.add(new Triangle(vids[0], vids[1], vids[2]));
    }

    return new TriangleMesh(poss, tris);
  }

  void show(color ct, color cq, color cc, color cb, boolean showStroke) {
    /* Show triangles. */
    if (showStroke) {
      stroke(0);
      strokeWeight(1);
    } else noStroke();
    fill(ct);
    beginShape(TRIANGLES);
    for (TriangleSub tri : triangles) {
      vertex(tri.vertices[0].position);
      vertex(tri.vertices[1].position);
      vertex(tri.vertices[2].position);
    }
    endShape();

    /* Show quads. */
    if (quads != null) {
      fill(cq);
      beginShape(QUADS);
      for (QuadSub quad : quads) {
        vertex(quad.vertices[0].position);
        vertex(quad.vertices[1].position);
        vertex(quad.vertices[2].position);
        vertex(quad.vertices[3].position);
      }
      endShape(QUADS);
    }
    if (showStroke) noStroke();

    /* Show circles. */
    stroke(cc);
    strokeWeight(3);
    // for (int i = 0; i < circles.length; ++i) circles[i].show();
    strokeWeight(1);
    noStroke();

    /* Show borders. */
    stroke(cb);
    strokeWeight(3);
    // for (int i = 0; i < borders.length; ++i) showLoop(borders[i]);
    strokeWeight(1);
    noStroke();
  }

  void show(boolean showStroke) {
    /* Show triangles. */
    if (showStroke) {
      stroke(0);
      strokeWeight(1);
    } else noStroke();
    beginShape(TRIANGLES);
    for (TriangleSub tri : triangles) {
      fill(tri.col);
      vertex(tri.vertices[0].position);
      vertex(tri.vertices[1].position);
      vertex(tri.vertices[2].position);
    }
    endShape();

    /* Show quads. */
    if (quads != null) {
      beginShape(QUADS);
      for (QuadSub quad : quads) {
        fill(quad.col);
        vertex(quad.vertices[0].position);
        vertex(quad.vertices[1].position);
        vertex(quad.vertices[2].position);
        vertex(quad.vertices[3].position);
      }
      endShape(QUADS);
    }
    if (showStroke) noStroke();
  }
}
