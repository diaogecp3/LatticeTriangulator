/******************************************************************************
 * Stereographic projection.
 ******************************************************************************/

enum StereoType {NORTH_POLE_SOUTH_PLANE, CUSTOM_NORTH_POLE_SOUTH_PLANE}

/* Steregraphic projection of p w.r.t. sphere (c, r). */
pt stereoProjectNorthPoleSouthPlane(pt c, float r, pt p) {
  pt q = new pt(0, 0, c.z - r);
  float t = 2 * r / (c.z + r - p.z);
  q.x = c.x + t * (p.x - c.x);
  q.y = c.y + t * (p.y - c.y);
  return q;
}

/* Stereographic projection of q w.r.t. sphere (c, r). */
pt inverseStereoProjectNorthPoleSouthPlane(pt c, float r, pt q) {
  // println("center =", c);
  float x = q.x - c.x;
  float y = q.y - c.y;

  float x2 = x * x;
  float y2 = y * y;
  float fr2 = 4 * r * r;
  float d = x2 + y2 + fr2;

  float px = fr2 * x / d;
  float py = fr2 * y / d;
  float pz = (x2 + y2 - fr2) * r / d;

  {
    // float k1 = px / q.x;
    // float k2 = py / q.y;
    // float k3 = (r - pz) / (2 * r);
    // println("k1, k2, k3 =", k1, k2, k3);
  }

  return new pt(px + c.x, py + c.y, pz + c.z);
}

class StereoProjector {
  pt center;
  float radius;
  StereoType sType;

  pt northPole = null;
  vec rotAxis = null;
  float rotAngle = 0;

  pt originSouthPlane = null;
  vec xAxisSouthPlane = null;
  vec yAxisSouthPlane = null;

  StereoProjector(pt c, float r, StereoType t) {
    center = c;
    radius = r;
    sType = t;

    if (sType == StereoType.NORTH_POLE_SOUTH_PLANE) {
      setNorthPole(P(center, V(0, 0, radius)));
    }
  }

  void setNorthPole(pt n) {
    northPole = n;
    vec cn = V(center, northPole);
    rotAxis = U(N(cn, V(0, 0, 1)));
    rotAngle = acosClamp(cn.z / radius);

    originSouthPlane = P(center, -1, cn);  // south pole
    xAxisSouthPlane = constructNormal(cn);
    yAxisSouthPlane = U(N(cn, xAxisSouthPlane));
  }

  pt project(pt p) {
    if (sType == StereoType.NORTH_POLE_SOUTH_PLANE) {
      return stereoProjectNorthPoleSouthPlane(center, radius, p);
    } else if (sType == StereoType.CUSTOM_NORTH_POLE_SOUTH_PLANE) {
      pt a = p.c().rot(rotAngle, rotAxis, center);  // same effect as rotating frame
      pt b = stereoProjectNorthPoleSouthPlane(center, radius, a);
      return b.rot(-rotAngle, rotAxis, center);  // same effect as rotating frame
    } else println("invalid stereo type!");
    return null;
  }

  pt inverse(pt q) {
    // vec d = V(northPole, q);
    // pt[] ps = intersectionLineSphere(northPole, d, center, radius);  // not efficient!
    // if (d(ps[0], northPole) < 0.001) return ps[1];
    // else return ps[0];

    if (sType == StereoType.NORTH_POLE_SOUTH_PLANE) {
      return inverseStereoProjectNorthPoleSouthPlane(center, radius, q);
    } else if (sType == StereoType.CUSTOM_NORTH_POLE_SOUTH_PLANE) {
      pt b = q.c().rot(rotAngle, rotAxis, center);
      pt a = inverseStereoProjectNorthPoleSouthPlane(center, radius, b);
      return a.rot(-rotAngle, rotAxis, center);
    } else println("invalid stereo type!");
    return null;
  }

  pt inverse(pt2 q) {
    return inverse(to3D(q));
  }

  pt to3D(pt2 p) {
    return P(originSouthPlane, p.x, xAxisSouthPlane, p.y, yAxisSouthPlane);
  }

  pt2 to2D(pt p) {
    vec d = V(originSouthPlane, p);
    return new pt2(dot(d, xAxisSouthPlane), dot(d, yAxisSouthPlane));
  }

  /* Represent a 2D circle, which lies on the projection plane, using 3D coordinates. */
  Circle to3D(Circle2 cir) {
    vec n = U(center, northPole);  // the up vector of the projection plane
    pt c = to3D(cir.center);
    return new Circle(c, n, cir.radius);
  }

  /* Represent a 3D circle, which lies on the projection plane, using 2D coordinates. */
  Circle2 to2D(Circle cir) {
    pt2 c = to2D(cir.c);
    return new Circle2(c, cir.r);
  }

  /* The circle cir must be on the sphere (center, radius). */
  Circle project(Circle cir) {
    if (sType == StereoType.NORTH_POLE_SOUTH_PLANE || sType == StereoType.CUSTOM_NORTH_POLE_SOUTH_PLANE) {
      /* Pick 3 points on circle, project them, construct the projected circle. */
      vec x = constructNormal(cir.n);
      vec y = N(cir.n, x);
      pt a = P(cir.c, cir.r, x);
      pt b = P(cir.c, cir.r, y);
      pt c = P(cir.c, -cir.r, x);
      pt pa = project(a);
      pt pb = project(b);
      pt pc = project(c);

      {  // debug
        // fill(magenta);
        // showBall(a, 5);
        // showBall(b, 5);
        // showBall(c, 5);

        // stroke(cyan);
        // strokeWeight(5);
        // showSegment(northPole, pa);
        // showSegment(northPole, pb);
        // showSegment(northPole, pc);
        // noStroke();
      }

      return circumcircleOfTriangle(pa, pb, pc);
    } else println("invalid stereo type!");
    return null;
  }

  /* Inverse a circle, which is represented in 2D and lies on the projection plane, back to the sphere. */
  Circle inverse(Circle2 cir) {
    if (sType == StereoType.NORTH_POLE_SOUTH_PLANE || sType == StereoType.CUSTOM_NORTH_POLE_SOUTH_PLANE) {
      pt2 c = cir.center;
      float r = cir.radius;
      pt2 p0 = new pt2(c.x + r, c.y);
      pt2 p1 = new pt2(c.x, c.y + r);
      pt2 p2 = new pt2(c.x - r, c.y);

      pt q0 = to3D(p0);
      pt q1 = to3D(p1);
      pt q2 = to3D(p2);

      pt pq0 = inverse(q0);
      pt pq1 = inverse(q1);
      pt pq2 = inverse(q2);

      return circumcircleOfTriangle(pq0, pq1, pq2);
    } else println("invalid stereo type!");
    return null;
  }

  /* Inverse a circle, which is represented in 3D and lies on the projection plane, back to the sphere. */
  Circle inverse(Circle cir) {
    Circle2 cir2 = to2D(cir);
    return inverse(cir2);
  }
}

/* Circle in 2D. */
class Circle2 {
  pt2 center = new pt2();
  float radius = 1;
  Circle2 () {};
  Circle2 (pt2 c, float r) {
    center = c;
    radius = r;
  }
  void setTo(Circle2 cir) {
    center = cir.center;
    radius = cir.radius;
  }
}

/*
 * Compute the two Apollonius circles of three given circles.
 * si are +/- 1 (s1 = +1 means C1 is outside of an Apollonius circle)
 * See https://rasmusfonseca.github.io/implementations/apollonius.html
 */
Circle2[] constructApolloniuCircles(Circle2 C1, Circle2 C2, Circle2 C3, int s1, int s2, int s3) {
  float x1 = C1.center.x;
  float y1 = C1.center.y;
  float r1 = C1.radius;
  float x2 = C2.center.x;
  float y2 = C2.center.y;
  float r2 = C2.radius;
  float x3 = C3.center.x;
  float y3 = C3.center.y;
  float r3 = C3.radius;

  float v11 = 2*x2 - 2*x1;
  float v12 = 2*y2 - 2*y1;
  float v13 = x1*x1 - x2*x2 + y1*y1 - y2*y2 - r1*r1 + r2*r2;
  float v14 = 2*s2*r2 - 2*s1*r1;

  float v21 = 2*x3 - 2*x2;
  float v22 = 2*y3 - 2*y2;
  float v23 = x2*x2 - x3*x3 + y2*y2 - y3*y3 - r2*r2 + r3*r3;
  float v24 = 2*s3*r3 - 2*s2*r2;

  float w12 = v12/v11;
  float w13 = v13/v11;
  float w14 = v14/v11;

  float w22 = v22/v21-w12;
  float w23 = v23/v21-w13;
  float w24 = v24/v21-w14;

  float P = -w23/w22;
  float Q = w24/w22;
  float M = -w12*P-w13;
  float N = w14 - w12*Q;

  float a = N*N + Q*Q - 1;
  float b = 2*M*N - 2*N*x1 + 2*P*Q - 2*Q*y1 + 2*s1*r1;
  float c = x1*x1 + M*M - 2*M*x1 + P*P + y1*y1 - 2*P*y1 - r1*r1;

  // Find roots of a quadratic equation
  float[] rs = solveQuadraticEquation(a, b, c);
  if (rs == null) return null;
  float[] xs = new float[] {M + N * rs[0], M + N * rs[1]};
  float[] ys = new float[] {P + Q * rs[0], P + Q * rs[1]};
  return new Circle2[] {new Circle2(new pt2(xs[0], ys[0]), abs(rs[0])),
                        new Circle2(new pt2(xs[1], ys[1]), abs(rs[1]))};
}
