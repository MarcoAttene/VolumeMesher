#include "extended_predicates.h"
#include <algorithm>


//----------------------
// Derivated predicates
//----------------------

// ----- Class1: the highest dimensional object is a POINT -------

//  Input: point p=(px,py,pz), q=(qx,qy,qz), r=(rx,ry,rz).
// Output: 1 -> if the 3 input points are NOT aligned,
//         0 -> otherwise
bool misAlignment(const double * p, const double * q, const double * r){

    // Projection on (x,y)-plane
    if( orient2d(p,q,r) )      return 1;

    // Projection on (y,z)-plane
    if( orient2d(p+1,q+1,r+1) )      return 1;

    // Projection on (x,z)-plane
    const double pxz[]={p[0],p[2]};
    const double qxz[]={q[0],q[2]};
    const double rxz[]={r[0],r[2]};
    return ( orient2d(pxz,qxz,rxz) );
}


// ----- Class2: the highest dimensional object is a SEGMENT -------


//  Input: point p=(px,py,pz); segment v1-v2 with v1=(v1x,v1y,v1z), v2=(v2x,v2y,v2z).
// Output: 1 -> if the point p belong to the interior of the segment v1-v2 (endpoints excluded),
//         0 -> otherwise.
uint32_t pointInInnerSegment(const double * p, const double * v1, const double * v2){
    return ( misAlignment(p, v1, v2)==0 && (
               ( v1[0]<v2[0] && v1[0]<p[0] && p[0]<v2[0] )
            || ( v1[0]>v2[0] && v1[0]>p[0] && p[0]>v2[0] )
            || ( v1[1]<v2[1] && v1[1]<p[1] && p[1]<v2[1] )
            || ( v1[1]>v2[1] && v1[1]>p[1] && p[1]>v2[1] )
            || ( v1[2]<v2[2] && v1[2]<p[2] && p[2]<v2[2] )
            || ( v1[2]>v2[2] && v1[2]>p[2] && p[2]>v2[2] )
    ));
}

//  Input: point p=(px,py,pz); segment v1-v2 with v1=(v1x,v1y,v1z), v2=(v2x,v2y,v2z).
// Output: 1 -> if the point p belong to the segment v1-v2 (endpoints included),
//         0 -> otherwise.
uint32_t pointInSegment(const double * p, const double * v1, const double * v2){
    return ( same_point(p, v1) || same_point(p, v2) || pointInInnerSegment(p,v1,v2) );
}

//  Input: point p, point q, and a segment v1-v2, througt their coordinates:
//         p=(px,py,pz), q=(qx,qy,qz), v1=(v1x,v1y,v1z), v2=(v2x,v2y,v2z).
// Output: 1 -> points p and q lies both on the same side of the straight line passing througt v1 and v2.
//         0 -> otherwise.
// Note. Points and segment must be complanar.
uint32_t same_half_plane(const double * p, const double * q, const double * v1, const double * v2){
    // Projection on (x,y)-plane
    if(signe_orient2d(p, v1, v2) != signe_orient2d(q, v1, v2))    return 0;

    // Projection on (y,z)-plane
    if(signe_orient2d(p + 1, v1 + 1, v2 + 1) != signe_orient2d(q + 1, v1 + 1, v2 + 1))    return 0;

    // Projection on (x,z)-plane
    const double pxz[]={p[0],p[2]};
    const double qxz[]={q[0],q[2]};
    const double v1xz[]={v1[0],v1[2]};
    const double v2xz[]={v2[0],v2[2]};
    return (signe_orient2d(pxz, v1xz, v2xz) == signe_orient2d(qxz, v1xz, v2xz));
}


//  Input: segments u1-u2 and v1-v2 through their coordinates:
//         u1=(u1x,u1y,u1z), u2=(u2x,u2y,u2z), v1=(v1x,v1y,v1z), v2=(v2x,v2y,v2z).
// Output: 1 -> segments properly intesects, i.e. intersection occours in both segments interior,
//         0 -> otherwise.
// Note. Collinear overlapping segments are not considered to be properly intersecting.
uint32_t innerSegmentsCross(const double * u1, const double * u2, const double * v1, const double * v2 ){

    // Segments have not to share any endpoints
    if(same_point(u1,v1) || same_point(u1,v2) ||
       same_point(u2,v1) || same_point(u2,v2)    )      return 0;

    // The 4 endpoints must be coplanar
    if( orient3d(u1,u2,v1,v2)!=0. )                     return 0;

    // Endpoints of one segment cannot stay either on the same side of the other one.
    if(same_half_plane(u1,u2,v1,v2) ||
       same_half_plane(v1,v2,u1,u2)   )                 return 0;

    // Each segment endpoint cannot be aligned with the other segment.
    if( !misAlignment(u1, v1, v2) )                     return 0;
    if( !misAlignment(u2, v1, v2) )                     return 0;
    if( !misAlignment(v1, u1, u2) )                     return 0;
    if( !misAlignment(v2, u1, u2) )                     return 0;

    // If the segment projected on one coordinate plane cross -> segmant cross.
    // Projection on (x,y)-plane
    if ( orient2d(u1,u2,v1) != 0. )                     return 1;
    if ( orient2d(v1,v2,u2) != 0. )                     return 1;
    // Projection on (y,z)-plane
    if ( orient2d(u1+1,u2+1,v1+1) != 0. )               return 1;
    if ( orient2d(v1+1,v2+1,u2+1) != 0. )               return 1;
    // Projection on (z,x)-plane
    const double u1xz[]={u1[0],u1[2]};
    const double u2xz[]={u2[0],u2[2]};
    const double v1xz[]={v1[0],v1[2]};
    const double v2xz[]={v2[0],v2[2]};
    if ( orient2d(u1xz,u2xz,v1xz) != 0. )               return 1;
    if ( orient2d(v1xz,v2xz,u2xz) != 0. )               return 1;

    return 0;
}


// ----- Class3: the highest dimensional object is a TRIANGLE -------


//  Input: point p and triangle <v1,v2,v3> through their coordinates:
//         p=(px,py,pz), v1=(v1x,v1y,v1z), v2=(v2x,v2y,v2z), v3=(v3x,v3y,v3z).
// Output: 1 -> point belong to the interior of the triangle,
//         0 -> otherwise.
uint32_t pointInInnerTriangle(const double * p, const double * v1, const double * v2, const double * v3){

    double o1, o2, oo2, oo4, oo6;

    // Projection on (x,y)-plane -> p VS v1
    o1 = signe_orient2d(p, v2, v3);
    o2 = signe_orient2d(v1, v2, v3);
    oo2 = o2;
    if ( o1 != o2 )         return 0;

    // Projection on (y,z)-plane -> p VS v1
    o1 = signe_orient2d(p+1, v2+1, v3+1);
    o2 = signe_orient2d(v1+1, v2+1, v3+1);
    oo4 = o2;
    if ( o1 != o2 )         return 0;

    // Projection on (x,z)-plane -> p VS v1
    const double pxz[]  = {p[0],p[2]};
    const double v1xz[] = {v1[0],v1[2]};
    const double v2xz[] = {v2[0],v2[2]};
    const double v3xz[] = {v3[0],v3[2]};
    o1 = signe_orient2d(pxz, v2xz, v3xz);
    o2 = signe_orient2d(v1xz, v2xz, v3xz);
    oo6 = o2;
    if ( o1 != o2 )         return 0;

    // Projection on (x,y)-plane -> p VS v2
    o1 = signe_orient2d(p, v3, v1);
    o2 = oo2;
    if ( o1 != o2 )         return 0;

    // Projection on (y,z)-plane -> p VS v2
    o1 = signe_orient2d(p+1, v3+1, v1+1);
    o2 = oo4;
    if ( o1 != o2 )         return 0;

    // Projection on (x,z)-plane -> p VS v2
    o1 = signe_orient2d(pxz, v3xz, v1xz);
    o2 = oo6;
    if ( o1 != o2 )         return 0;

    // Projection on (x,y)-plane -> p VS v3
    o1 = signe_orient2d(p, v1, v2);
    o2 = oo2;
    if ( o1 != o2 )         return 0;

    // Projection on (y,z)-plane -> p VS v3
    o1 = signe_orient2d(p+1, v1+1, v2+1);
    o2 = oo4;
    if ( o1 != o2 )         return 0;

    // Projection on (x,z)-plane -> p VS v3
    o1 = signe_orient2d(pxz, v1xz, v2xz);
    o2 = oo6;
    if ( o1 != o2 )         return 0;

    return 1;
}

//  Input: point p and triangle <v1,v2,v3> through their coordinates:
//         p=(px,py,pz), v1=(v1x,v1y,v1z), v2=(v2x,v2y,v2z), v3=(v3x,v3y,v3z).
// Output: 1 -> point belong to the  triangle (interior or boundary),
//         0 -> otherwise.
uint32_t pointInTriangle(const double * p, const double * v1, const double * v2, const double * v3){

    return ( pointInSegment(p,v1,v2) ||
             pointInSegment(p,v2,v3) ||
             pointInSegment(p,v3,v1) ||
             pointInInnerTriangle(p,v1,v2,v3) );

}

//  Input: segment u1-u2 and triangle <v1,v2,v3> through their coordinates:
//         u1=(u1x,u1y,u1z), u2=(u2x,u2y,u2z), v1=(v1x,v1y,v1z), v2=(v2x,v2y,v2z), v3=(v3x,v3y,v3z).
// Output: 1 -> segment and triangle properly intesects, i.e. intersection occours in both segment and triangle interior,
//        0 -> otherwise.
// Note. Collinear overlapping segments are not considered to be properly intersecting.
uint32_t innerSegmentCrossesInnerTriangle(const double * u1, const double * u2, const double * v1, const double * v2, const double * v3){

    // "out of the Box" check.
    double bound;
    bound = std::min(u1[0],u2[0]);                       // min(u1,u2) alogng x-axis
    if(v1[0]<=bound && v2[0]<=bound && v3[0]<=bound)         return 0;
    bound = std::max(u1[0],u2[0]);                       // max(u1,u2) alogng x-axis
    if(v1[0]>=bound && v2[0]>=bound && v3[0]>=bound)         return 0;
    bound = std::min(u1[1],u2[1]);                       // min(u1,u2) alogng y-axis
    if(v1[1]<=bound && v2[1]<=bound && v3[1]<=bound)         return 0;
    bound = std::max(u1[1],u2[1]);                       // max(u1,u2) alogng y-axis
    if(v1[1]>=bound && v2[1]>=bound && v3[1]>=bound)         return 0;
    bound = std::min(u1[2],u2[2]);                       // min(u1,u2) alogng z-axis
    if(v1[2]<=bound && v2[2]<=bound && v3[2]<=bound)         return 0;
    bound = std::max(u1[2],u2[2]);                       // max(u1,u2) alogng z-axis
    if(v1[2]>=bound && v2[2]>=bound && v3[2]>=bound)         return 0;

    const int orient_u1_tri = signe_orient3d(u1,v1,v2,v3);
    const int orient_u2_tri = signe_orient3d(u2,v1,v2,v3);

    // Check if triangle vertices and at least one of the segment endpoints are coplanar:
    // in this case there is no proper intersection.
    if( orient_u1_tri==0 || orient_u2_tri==0 )             return 0;

    // Endpoints of one segment cannot stay both in one of the same half-space defined by the triangle.
    if( orient_u1_tri == orient_u2_tri )                   return 0;


    // Since now, endpoints are one abouve and one below the triangle-plane.

    // Intersection between segment and triangle sides are not proper.
    // Check also if segment intersect the triangle-plane outside the triangle.
    const int orient_u_v1v2 = signe_orient3d(u1,u2,v1,v2);
    const int orient_u_v2v3 = signe_orient3d(u1,u2,v2,v3);

    if( orient_u_v1v2==0 || orient_u_v2v3==0 )              return 0;
    if( orient_u_v1v2 != orient_u_v2v3 )                    return 0;

    const int orient_u_v3v1 = signe_orient3d(u1,u2,v3,v1);

    if( orient_u_v3v1==0 )                                  return 0;
    if( orient_u_v3v1 != orient_u_v2v3 )                    return 0;


    // Finally, we have a proper intersection.
    return 1;
}

//  Input: segment u1-u2 and triangle <v1,v2,v3> through their coordinates:
//         u1=(u1x,u1y,u1z), u2=(u2x,u2y,u2z),
//         v1=(v1x,v1y,v1z), v2=(v2x,v2y,v2z), v3=(v3x,v3y,v3z).
// Output: 1 -> inner segment and triangle intesects,
//              i.e. intersection occours between segment interior and triangle,
//         0 -> otherwise.
uint32_t innerSegmentCrossesTriangle(const double * u1, const double * u2,
                          const double * v1, const double * v2, const double * v3){
    if( pointInInnerSegment(v1, u1, u2)                     ||
        pointInInnerSegment(v2, u1, u2)                     ||
        pointInInnerSegment(v3, u1, u2)                     ||
        innerSegmentsCross(v2, v3, u1, u2)                  ||
        innerSegmentsCross(v3, v1, u1, u2)                  ||
        innerSegmentsCross(v1, v2, u1, u2)                  ||
        innerSegmentCrossesInnerTriangle(u1, u2, v1, v2, v3)   )  return 1;
    return 0;
}




// SLOW FUNCTION - Use for checking purposes only !
//
// Returns TRUE is the area of the intersection of the two triangles is non zero

int triangles_overlap(const double* v11, const double* v12, const double* v13, const double* v21, const double* v22, const double* v23)
{
    // If not coplanar -> FALSE
    if (orient3d(v11, v12, v13, v21) || orient3d(v11, v12, v13, v22) || orient3d(v11, v12, v13, v23)) return 0;

    // If one has vertex into the other -> TRUE
    if (pointInInnerTriangle(v11, v21, v22, v23) || pointInInnerTriangle(v12, v21, v22, v23) || pointInInnerTriangle(v13, v21, v22, v23) ||
        pointInInnerTriangle(v21, v11, v12, v13) || pointInInnerTriangle(v22, v11, v12, v13) || pointInInnerTriangle(v23, v11, v12, v13)) return 1;

    // If edges corss -> TRUE
    if (innerSegmentsCross(v11, v12, v21, v22) || innerSegmentsCross(v12, v13, v21, v22) || innerSegmentsCross(v13, v11, v21, v22) ||
        innerSegmentsCross(v11, v12, v22, v23) || innerSegmentsCross(v12, v13, v22, v23) || innerSegmentsCross(v13, v11, v22, v23) ||
        innerSegmentsCross(v11, v12, v23, v21) || innerSegmentsCross(v12, v13, v23, v21) || innerSegmentsCross(v13, v11, v23, v21)) return 1;

    if (pointInTriangle(v11, v21, v22, v23) && pointInTriangle(v12, v21, v22, v23) && pointInTriangle(v13, v21, v22, v23)) return 1;
    if (pointInTriangle(v21, v11, v12, v13) && pointInTriangle(v22, v11, v12, v13) && pointInTriangle(v23, v11, v12, v13)) return 1;

    return 0;
}



static inline uint32_t pointInTriangleCC(const double*  p, const double*  v1, const double*  v2, const double*  v3) {
    return (orient3d(p, v1, v2, v3) == 0 && pointInTriangle(p, v1, v2, v3));
}

// SLOW FUNCTION - Use for checking purposes only !
//
// Let V be a triangle and T be a tetrahedron, and let's assume that none of the
// three vertices of V belongs to the interior int(T) of T.
// If I is the intersection of V and int(T), this function returns TRUE
// if the area of I is non zero.

int triangle_intersects_inner_tet(
    const double* v1, const double* v2, const double* v3,
    const double* t1, const double* t2, const double* t3, const double* t4)
{
    // Exclude face coplanarity : at least one positive and one negative orient3d
    const double o1 = orient3d(t1, v1, v2, v3);
    const double o2 = orient3d(t2, v1, v2, v3);
    const double o3 = orient3d(t3, v1, v2, v3);
    const double o4 = orient3d(t4, v1, v2, v3);
    if ((o1 <= 0 && o2 <= 0 && o3 <= 0 && o4 <= 0) || (o1 >= 0 && o2 >= 0 && o3 >= 0 && o4 >= 0)) return 0;

    // Inner tet edge intersects inner triangle
   if (innerSegmentCrossesInnerTriangle(t1, t2, v1, v2, v3) ||
       innerSegmentCrossesInnerTriangle(t2, t3, v1, v2, v3) ||
       innerSegmentCrossesInnerTriangle(t3, t4, v1, v2, v3) ||
       innerSegmentCrossesInnerTriangle(t1, t3, v1, v2, v3) ||
       innerSegmentCrossesInnerTriangle(t2, t4, v1, v2, v3) ||
       innerSegmentCrossesInnerTriangle(t4, t1, v1, v2, v3)) return 1;

    // Inner inner triangle edge intersects inner tet face
    if (innerSegmentCrossesInnerTriangle(v1, v2, t1, t2, t3) ||
        innerSegmentCrossesInnerTriangle(v1, v2, t2, t3, t4) ||
        innerSegmentCrossesInnerTriangle(v1, v2, t3, t4, t1) ||
        innerSegmentCrossesInnerTriangle(v1, v2, t1, t2, t4)) return 1;

    if (innerSegmentCrossesInnerTriangle(v2, v3, t1, t2, t3) ||
        innerSegmentCrossesInnerTriangle(v2, v3, t2, t3, t4) ||
        innerSegmentCrossesInnerTriangle(v2, v3, t3, t4, t1) ||
        innerSegmentCrossesInnerTriangle(v2, v3, t1, t2, t4)) return 1;

    if (innerSegmentCrossesInnerTriangle(v3, v1, t1, t2, t3) ||
        innerSegmentCrossesInnerTriangle(v3, v1, t2, t3, t4) ||
        innerSegmentCrossesInnerTriangle(v3, v1, t3, t4, t1) ||
        innerSegmentCrossesInnerTriangle(v3, v1, t1, t2, t4)) return 1;

    // All three triangle vertices are on tet faces
    if ((pointInTriangleCC(v1, t1, t2, t3) || pointInTriangleCC(v1, t2, t3, t4) || pointInTriangleCC(v1, t4, t1, t2) || pointInTriangleCC(v1, t1, t4, t3)) &&
        (pointInTriangleCC(v2, t1, t2, t3) || pointInTriangleCC(v2, t2, t3, t4) || pointInTriangleCC(v2, t4, t1, t2) || pointInTriangleCC(v2, t1, t4, t3)) &&
        (pointInTriangleCC(v3, t1, t2, t3) || pointInTriangleCC(v3, t2, t3, t4) || pointInTriangleCC(v3, t4, t1, t2) || pointInTriangleCC(v3, t1, t4, t3))) return 1;

    return 0;
}
