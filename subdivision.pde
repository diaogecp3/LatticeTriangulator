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

/*
 * The vertex class used in subdivision.
 */
class VertexSub {
  pt position;
  int circleID;  // the ID of the circle where it lies, -1: no such circle

  /*
   * The normal of a vertex is used in projection. Note that a vertex adjacent
   * to two faces has to consider two face normals.
   * Todo: merge two vertices if they are close to each other, and merge their
   * normals.
   */
  vec normal;

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
private pt projectPointOnCircle(pt p, vec n, pt c, float r) {
  if (n == null) {
    return P(c, r, U(c, p));  // radial projection, not always true
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
 * Project the midpoint of a chord (a, b) on circle with center c, radius r and
 * normal n. The projected point x should be on the positive side of (a, b),
 * i.e., triangle (a, b, x) should has the same normal as that of the circle.
 */
private pt projectChordMidpointOnCircle(pt a, pt b, pt c, float r, vec n) {
  if (!samePt(a, b)) {
    vec v = U(N(n, V(a, b)));
    return P(c, r, v);
  } else {  // a is close to b
    pt m = P(a, b);
    vec u = U(c, m);  // c is not close to m unless the circle is so small
    pt x = P(c, r, u);
    if (dot(N(a, b, x), n) < 0) x = P(c, -r, u);
    return x;
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

  /*
   * Subdivie a triangle into 4 triangles.
   */
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

          // vec edgeNormal = U(N(n, V(va.position, vb.position)));  // TODO: do we need to store edge normal?
          if (projectOnCircle) {
            p = projectChordMidpointOnCircle(va.position, vb.position, c, r, n);
            // p = projectPointOnCircle(p, null, c, r);  // radial projection
            // p = projectPointOnCircle(p, edgeNormal, c, r);  // this doesn't work
          }

          mids[i] = new VertexSub(p, va.circleID, null);
          borders[va.circleID].add(p);  // this border is not sorted after insertion
        } else {  // not on the same circle
          pt p = P(va.position, vb.position);
          mids[i] = new VertexSub(p, -1, normal);  // todo: this normal may not be the normal at that vertex
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

  /*
   * Create 4 edge midpoints and 1 face midpoints.
   */
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

        // vec edgeNormal = U(N(n, V(va.position, vb.position)));
        if (projectOnCircle) {
          p = projectChordMidpointOnCircle(va.position, vb.position, c, r, n);
          // p = projectPointOnCircle(p, null, c, r);  // radial projection
          // p = projectPointOnCircle(p, edgeNormal, c, r);  // this doesn't work
        }

        mids[i] = new VertexSub(p, va.circleID, null);
        borders[va.circleID].add(p);  // this border is not sorted after insertion
      } else {  // not on the same circle
        pt p = P(va.position, vb.position);
        mids[i] = new VertexSub(p, -1, normal);  // todo: this normal may not be the normal at that vertex
      }
    }

    {
      pt p = P(vertices[0].position, vertices[1].position, vertices[2].position, vertices[3].position);
      mids[4] = new VertexSub(p, -1, normal);
    }

    return mids;
  }

  /*
   * Subdivide a quad into 8 triangles.
   */
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
class BorderedTriQuadMesh extends Mesh {
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

  /*
   * Use ct as the color for triangles and cq as the color for quads.
   */
  void setColors(color ct, color cq) {
    for (TriangleSub tri : triangles) tri.setColor(ct);
    if (quads != null) {
      for (QuadSub quad : quads) quad.setColor(cq);
    }
  }

  /*
   * Subdivide each face of the mesh.
   *
   * Parameters:
   * subTypeTri: subdivision type for triangles
   * subTypeQuad: subdivision type for quads
   * projectOnCircle: project vertices on their corresponding circles or not
   */
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

  /*
   * Project vertices on a hub by shooting rays. Note that a vertex may be
   * outside or inside the hub.
   */
  private void projectOnHubRay(Hub hub) {
    if (debugProjection) dProjectionInfo.reset();

    for (TriangleSub tri : triangles) {
      for (VertexSub v : tri.vertices) {
        pt p = v.position;
        projectPointOnExactHub(hub, p);

        {
          // float dist = hub.distanceFrom(p);

          // // if (isZero(dist, 0.01)) continue;  // on the hub
          // if (dist < 0.0001) continue;  // inside or on

          // // if (dist < -0.0001) {  // inside the hub
          //   // fill(magenta, 100);
          //   // showBall(p, 3);
          // // }

          // vec d = U(p, hub.ball.c);
          // // if (dist < -0.01) d.rev();  // flip the ray direction if p is inside the hub
          // Float t = hub.closestIntersectionWithRay(p, d);
          // if (t != null) p.add(t, d);
        }

        {
          // /* The vertex normal is not correct. */
          // if (dot(d, v.normal) < 0) {  // hub center cannot see vertex v
          //   Float t = hub.closestIntersectionWithRay(p, d);
          //   if (t != null) p.add(t, d);
          // } else {  // hub center can see vertex v
          //   println("hub center can see a vertex");
          //   {  // debug
          //     // fill(yellow, 100);
          //     // showBall(p, 4);
          //   }

          //   /* Redirect the ray direction. */
          //   d.rev();
          //   vec r = reflect(d, v.normal);
          //   r.normalize();

          //   // Float t = hub.furthestIntersectionWithRay(p, v.normal);  // use vertex normal
          //   // if (t != null) p.add(t, v.normal);
          //   Float t = hub.furthestIntersectionWithRay(hub.ball.c, r);  // use reflected radial projection
          //   if (t != null) p.set(P(hub.ball.c, t, r));
          // }
        }

      }
    }
  }

  /*
   * Project vertices on a blended hub by tracing spheres. Assume that the hub
   * is completely inside the mesh.
   */
  private void projectOnHubSphereTrace(Hub hub) {
    for (TriangleSub tri : triangles) {
      for (VertexSub v : tri.vertices) {
        projectPointOnBlendedHub(hub, v.position);
      }
    }
  }

  /*
   * Project vertices on a hub. Projection is performed after subdivision. Thus,
   * only triangles need to be considered.
   */
  void projectOnHub(Hub hub, ProjectType option) {
    if (option == null) return;
    switch (option) {
      case RAY:  // project radially on the exact hub
        projectOnHubRay(hub);
        break;
      case SPHERE_TRACING:  // project radially on the blended hub
        projectOnHubSphereTrace(hub);
        break;
    }
  }

  private void tessellateAndMergeVerticesForPolygons(
    ArrayList<pt> poss,
    ArrayList<Triangle> tris) {
    int n = circles.length;
    for (int k = 0; k < n; ++k) {
      pt pc = circles[k].c;
      int ic = poss.size();  // vertex index of the center of the k-th circle
      poss.add(pc);

      ArrayList<Integer> vids = new ArrayList<Integer>();
      for (pt p : borders[k]) {
        int j = 0;
        for (; j < poss.size(); ++j) {
          if (samePt(p, poss.get(j))) {  // these two vertices are almost the same
            vids.add(j);
            break;
          }
        }
        if (j == poss.size()) {  // find a new vertex
          vids.add(j);
          poss.add(p);
        }
      }

      int m = borders[k].size();  // m points on the k-th border
      for (int i = 0; i <= m - 2; ++i) {
        tris.add(new Triangle(ic, vids.get(i), vids.get(i+1)));
      }
      tris.add(new Triangle(ic, vids.get(m-1), vids.get(0)));
    }
  }

  private void tessellateAndMergeVerticesForQuads(
    ArrayList<pt> poss,
    ArrayList<Triangle> tris) {
    if (quads == null || quads.size() == 0) return;
    for (QuadSub quad : quads) {
      int[] vids = new int[] {-1, -1, -1, -1};
      for (int i = 0; i < 4; ++i) {
        pt p = quad.vertices[i].position;
        int j = 0;
        for (; j < poss.size(); ++j) {
          if (samePt(p, poss.get(j))) {  // these two vertices are almost the same
            vids[i] = j;
            break;
          }
        }
        if (j == poss.size()) {  // find a new vertex
          vids[i] = j;
          poss.add(p);
        }
      }
      assert vids[0] != -1 && vids[1] != -1 && vids[2] != -1 && vids[3] != -1;
      tris.add(new Triangle(vids[0], vids[1], vids[2]));
      tris.add(new Triangle(vids[0], vids[2], vids[3]));
    }
  }

  private void mergeVerticesForTriangles(
    ArrayList<pt> poss,
    ArrayList<Triangle> tris) {
    for (TriangleSub tri : triangles) {
      int[] vids = new int[] {-1, -1, -1};
      for (int i = 0; i < 3; ++i) {
        pt p = tri.vertices[i].position;
        int j = 0;
        for (; j < poss.size(); ++j) {
          if (samePt(p, poss.get(j))) {  // these two vertices are almost the same
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
  }

  /*
   * Convert this triangle-quad mesh to a triangle mesh. Please make sure that
   * no quads exist in this triangle-quad mesh. Vertices are merged if they are
   * at the same location.
   */
  TriangleMesh toTriangleMesh() {
    assert quads == null || quads.size() == 0;

    ArrayList<pt> poss = new ArrayList<pt>();
    ArrayList<Triangle> tris = new ArrayList<Triangle>();

    // for (TriangleSub tri : triangles) {
    //   int[] vids = new int[] {-1, -1, -1};
    //   for (int i = 0; i < 3; ++i) {
    //     pt p = tri.vertices[i].position;
    //     int j = 0;
    //     for (; j < poss.size(); ++j) {
    //       if (samePt(p, poss.get(j))) {  // these two vertices are almost the same
    //         vids[i] = j;
    //         break;
    //       }
    //     }
    //     if (j == poss.size()) {  // find a new vertex
    //       vids[i] = j;
    //       poss.add(p);
    //     }
    //   }
    //   assert vids[0] != -1 && vids[1] != -1 && vids[2] != -1;
    //   tris.add(new Triangle(vids[0], vids[1], vids[2]));
    // }

    mergeVerticesForTriangles(poss, tris);
    return new TriangleMesh(poss, tris);
  }

  /*
   * Tessellate the polygon faces, each bounded by a polygonal border, and quads.
   * Collect all triangles (including those tessellated from polygons and quads)
   * into a TriangleMesh.
   */
  TriangleMesh tessellate() {
    ArrayList<pt> poss = new ArrayList<pt>();
    ArrayList<Triangle> tris = new ArrayList<Triangle>();

    tessellateAndMergeVerticesForPolygons(poss, tris);
    tessellateAndMergeVerticesForQuads(poss, tris);
    mergeVerticesForTriangles(poss, tris);

    return new TriangleMesh(poss, tris);
  }

  /*
   * Show the mesh and its borders.
   *
   * Parameters:
   * ct: color for triangles
   * cq: color for quads
   * cc: color for circles
   * cb: color for borders
   * showStroke: show the stroke or not
   */
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
    for (int i = 0; i < circles.length; ++i) circles[i].show();
    strokeWeight(1);
    noStroke();

    /* Show borders. */
    stroke(cb);
    strokeWeight(3);
    for (int i = 0; i < borders.length; ++i) showLoop(borders[i]);
    strokeWeight(1);
    noStroke();
  }

  /*
   * Show the mesh. Each face has its own color.
   */
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

  @Override
  void show(color colorMesh, boolean showStroke) {
    setColors(colorMesh, colorMesh);
    show(showStroke);
  }
}
