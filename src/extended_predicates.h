#ifndef extended_predicates_h
#define extended_predicates_h

#include <stdio.h>
#include "implicit_point.h"

#include "delaunay.h"

#pragma intrinsic(fabs)

//------------------------
// Basic operations
//------------------------

//  Input: two points p and q, having coordinates p=(px,py,pz), q=(qx,qy,qz).
// Output: 1 -> if p and q have the same coordinates,
//         0 -> otherwise.
static inline uint32_t same_point(const double* p, const double* q) {
    return (p[0] == q[0] && p[1] == q[1] && p[2] == q[2]);
}

//-----------------------------------
// Signe-version of basic predicates
//-----------------------------------

static inline int signe_orient2d(const double* p, const double* q, const double* r) {
    return orient2d(p[0], p[1], q[0], q[1], r[0], r[1]);
}

static inline int signe_orient3d(const double* p, const double* q, const double* r, const double* s) {
    return orient3d(p[0], p[1], p[2], q[0], q[1], q[2], r[0], r[1], r[2], s[0], s[1], s[2]);
}

bool misAlignment(const double * p, const double * q, const double * r);

uint32_t same_half_plane(const double * p, const double * q, const double * v1, const double * v2);
uint32_t pointInInnerSegment(const double * p, const double * v1, const double * v2);
uint32_t pointInSegment(const double * p, const double * v1, const double * v2);
uint32_t innerSegmentsCross(const double* u1, const double* u2, const double* v1, const double* v2);

uint32_t pointInInnerTriangle(const double * p, const double * v1, const double * v2, const double * v3);
uint32_t pointInTriangle(const double * p, const double * v1, const double * v2, const double * v3);
uint32_t innerSegmentCrossesInnerTriangle(const double * u1, const double * u2, const double * v1, const double * v2, const double * v3);
uint32_t innerSegmentCrossesTriangle(const double * u1, const double * u2,
                          const double * v1, const double * v2, const double * v3);

int triangles_overlap(const double* v11, const double* v12, const double* v13, const double* v21, const double* v22, const double* v23);
int triangle_intersects_inner_tet(
    const double* v1, const double* v2, const double* v3,
    const double* t1, const double* t2, const double* t3, const double* t4);


#endif /* extended_predicates_h */
