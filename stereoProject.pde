/******************************************************************************
 * Stereographic projection.
 ******************************************************************************/

enum StereoType {NORTH_POLE_SOUTH_PLANE, CUSTOM_NORTH_POLE_SOUTH_PLANE}

pt stereoProjectNorthPoleSouthPlane(pt c, float r, pt p) {
  pt q = new pt(0, 0, -r);
  float t = 2 * r / (r - p.z);
  q.x = t * p.x;
  q.y = t * p.y;
  return q;
}

class StereoProjector {
  pt center;
  float radius;
  StereoType sType;

  pt northPole = null;
  vec rotAxis = null;
  float rotAngle = 0;

  StereoProjector(pt c, float r, StereoType t) {
    center = c;
    radius = r;
    sType = t;
  }

  void setNorthPole(pt n) {
    northPole = n;
    vec cn = V(center, northPole);
    rotAxis = U(N(cn, V(0, 0, 1)));
    rotAngle = acosClamp(cn.z / radius);
  }

  pt project(pt p) {
    if (sType == StereoType.NORTH_POLE_SOUTH_PLANE) {
      return stereoProjectNorthPoleSouthPlane(center, radius, p);
    } else if (sType == StereoType.CUSTOM_NORTH_POLE_SOUTH_PLANE) {
      pt a = p.c().rot(rotAngle, rotAxis);
      pt b = stereoProjectNorthPoleSouthPlane(center, radius, a);
      return b.rot(-rotAngle, rotAxis);
    } else println("invalid stereo type!");
    return null;
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
      return circumcircleOfTriangle(pa, pb, pc);
    } else println("invalid stereo type!");
    return null;
  }
}