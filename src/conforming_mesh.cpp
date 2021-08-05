#include "conforming_mesh.h"
#include "extended_predicates.h"
#include "implicit_point.h"

#define INTERSECTION 1
#define IMPROPER_INTERSECTION 2
#define IMPROPER_INTERSECTION_COUNTED 3
#define PROPER_INTERSECTION_COUNTED 4
#define OVERLAP2D_F0 10
#define OVERLAP2D_F1 11
#define OVERLAP2D_F2 12
#define OVERLAP2D_F3 13
#define OVERLAP2D_F0_COUNTED 20
#define OVERLAP2D_F1_COUNTED 21
#define OVERLAP2D_F2_COUNTED 22
#define OVERLAP2D_F3_COUNTED 23

#define INSERT_TET_IN_LIST(tet,array,ind,marker) if(marker[tet]==0) array[ind++]=tet


//-------------------------------------
// Simplexes and Sub-Simplexes vertices
//-------------------------------------

static inline void extract_tetVrts(uint32_t* tet_vrts, uint64_t tet_ind,
                                   const TetMesh* mesh){
    std::memcpy(tet_vrts, mesh->tet_node + 4 * tet_ind, 4 * sizeof(uint32_t));
    //tet_vrts[0] = mesh->tet_node[4*tet_ind    ];
    //tet_vrts[1] = mesh->tet_node[4*tet_ind + 1];
    //tet_vrts[2] = mesh->tet_node[4*tet_ind + 2];
    //tet_vrts[3] = mesh->tet_node[4*tet_ind + 3];
}

// Given a tet = <v0,v1,v2,v3> the i-th face is opposite to vertex vi
// (ex. face_0 = <v1,v2,v3> is opposite to v0)
static inline void extract_tetFaceVrts(uint32_t* tetFace_vrts, uint64_t tet_ind,
                                       uint64_t ID, const TetMesh* mesh){
    tetFace_vrts[0] = mesh->tet_node[4*tet_ind + (ID+1)%4];
    tetFace_vrts[1] = mesh->tet_node[4*tet_ind + (ID+2)%4];
    tetFace_vrts[2] = mesh->tet_node[4*tet_ind + (ID+3)%4];
}

// Given a tet = <v0,v1,v2,v3> and the ID of two vertices returns the relative
// edge endpoints. (ex.  1,3 -> edge_1,3 = <v1,v3>)
static inline void extract_tetEdgeVrts(uint32_t* tetEdge_vrts, uint64_t tet_ind,
                                       uint64_t ID1, uint64_t ID2,
                                       const TetMesh* mesh){
    tetEdge_vrts[0] = mesh->tet_node[4*tet_ind + ID1];
    tetEdge_vrts[1] = mesh->tet_node[4*tet_ind + ID2];
}

//---------------------------------
// Geometric Predicates Interfaces
//---------------------------------
// Interfaces to use geometric predicates with vertex
// indices instead of coordinates

static inline uint32_t vrt_pointInInnerSegment(uint32_t pt_ind,
                                      uint32_t endpt0_ind, uint32_t endpt1_ind,
                                               const TetMesh* mesh){

    return (pt_ind!=endpt0_ind && pt_ind!=endpt1_ind &&
        pointInInnerSegment(mesh->vertices[pt_ind].coord, mesh->vertices[endpt0_ind].coord, mesh->vertices[endpt1_ind].coord));
}

static inline uint32_t vrt_pointInSegment(uint32_t pt_ind,
                                      uint32_t endpt0_ind, uint32_t endpt1_ind,
                                          const TetMesh* mesh){

    return (pt_ind==endpt0_ind || pt_ind==endpt1_ind ||
        pointInSegment(mesh->vertices[pt_ind].coord, mesh->vertices[endpt0_ind].coord, mesh->vertices[endpt1_ind].coord));
}

static inline uint32_t vrt_pointInInnerTriangle(uint32_t pt_ind,
            uint32_t tri_vrt0_ind, uint32_t tri_vrt1_ind, uint32_t tri_vrt2_ind,
                                                const TetMesh* mesh){

    return (pt_ind != tri_vrt0_ind &&
       pt_ind != tri_vrt1_ind &&
       pt_ind != tri_vrt2_ind &&
        pointInInnerTriangle(mesh->vertices[pt_ind].coord, 
            mesh->vertices[tri_vrt0_ind].coord, mesh->vertices[tri_vrt1_ind].coord, mesh->vertices[tri_vrt2_ind].coord));
}

static inline uint32_t vrt_pointInTriangle(uint32_t pt_ind,
            uint32_t tri_vrt0_ind, uint32_t tri_vrt1_ind, uint32_t tri_vrt2_ind,
                                           const TetMesh* mesh){

    return (pt_ind == tri_vrt0_ind ||
       pt_ind == tri_vrt1_ind ||
       pt_ind == tri_vrt2_ind ||
        pointInTriangle(mesh->vertices[pt_ind].coord,
            mesh->vertices[tri_vrt0_ind].coord, mesh->vertices[tri_vrt1_ind].coord, mesh->vertices[tri_vrt2_ind].coord));
}

static inline uint32_t vrt_innerSegmentsCross(uint32_t endpt0sgmA_ind,
                                              uint32_t endpt1sgmA_ind,
                                              uint32_t endpt0sgmB_ind,
                                              uint32_t endpt1sgmB_ind,
                                              const TetMesh* mesh){

    if(endpt0sgmA_ind == endpt0sgmB_ind || endpt0sgmA_ind == endpt1sgmB_ind ||
       endpt1sgmA_ind == endpt0sgmB_ind || endpt1sgmA_ind == endpt1sgmB_ind   )
         return 0;

    const double* e0A_coord = mesh->vertices[endpt0sgmA_ind].coord;
    const double* e1A_coord = mesh->vertices[endpt1sgmA_ind].coord;
    const double* e0B_coord = mesh->vertices[endpt0sgmB_ind].coord;
    const double* e1B_coord = mesh->vertices[endpt1sgmB_ind].coord;

    return innerSegmentsCross(e0A_coord, e1A_coord, e0B_coord, e1B_coord);
}

static inline uint32_t vrt_innerSegmentCrossesInnerTriangle(uint32_t endpt0_ind,
                                                            uint32_t endpt1_ind,
                                                          uint32_t tri_vrt0_ind,
                                                          uint32_t tri_vrt1_ind,
                                                          uint32_t tri_vrt2_ind,
                                                            const TetMesh* mesh){


    if(endpt0_ind==tri_vrt0_ind || endpt1_ind==tri_vrt0_ind ||
       endpt0_ind==tri_vrt1_ind || endpt1_ind==tri_vrt1_ind ||
       endpt0_ind==tri_vrt2_ind || endpt1_ind==tri_vrt2_ind    ) return 0;

    const double* e0_coord = mesh->vertices[endpt0_ind].coord;
    const double* e1_coord = mesh->vertices[endpt1_ind].coord;
    const double* t0_coord = mesh->vertices[tri_vrt0_ind].coord;
    const double* t1_coord = mesh->vertices[tri_vrt1_ind].coord;
    const double* t2_coord = mesh->vertices[tri_vrt2_ind].coord;
    return innerSegmentCrossesInnerTriangle(e0_coord, e1_coord,
                                            t0_coord, t1_coord, t2_coord);
}

static inline uint32_t vrt_innerSegmentCrossesTriangle(uint32_t endpt0_ind,
                                                       uint32_t endpt1_ind,
                                                       uint32_t tri_vrt0_ind,
                                                       uint32_t tri_vrt1_ind,
                                                       uint32_t tri_vrt2_ind,
                                                       const TetMesh* mesh){

    if(endpt0_ind==tri_vrt0_ind || endpt1_ind==tri_vrt0_ind ||
       endpt0_ind==tri_vrt1_ind || endpt1_ind==tri_vrt1_ind ||
       endpt0_ind==tri_vrt2_ind || endpt1_ind==tri_vrt2_ind    ) return 0;

    const double* e0_coord = mesh->vertices[endpt0_ind].coord;
    const double* e1_coord = mesh->vertices[endpt1_ind].coord;
    const double* t0_coord = mesh->vertices[tri_vrt0_ind].coord;
    const double* t1_coord = mesh->vertices[tri_vrt1_ind].coord;
    const double* t2_coord = mesh->vertices[tri_vrt2_ind].coord;
    return innerSegmentCrossesTriangle(e0_coord, e1_coord,
                                       t0_coord, t1_coord, t2_coord);
}

static inline uint32_t vrt_innerTrianglesCrosses(uint32_t triA_vrt0_ind,
                                                 uint32_t triA_vrt1_ind,
                                                 uint32_t triA_vrt2_ind,
                                                 uint32_t triB_vrt0_ind,
                                                 uint32_t triB_vrt1_ind,
                                                 uint32_t triB_vrt2_ind,
                                                 const TetMesh* mesh){
    if(vrt_innerSegmentCrossesInnerTriangle(triA_vrt0_ind,triA_vrt1_ind,
                            triB_vrt0_ind,triB_vrt1_ind,triB_vrt2_ind, mesh) ||
       vrt_innerSegmentCrossesInnerTriangle(triA_vrt1_ind,triA_vrt2_ind,
                            triB_vrt0_ind,triB_vrt1_ind,triB_vrt2_ind, mesh) ||
       vrt_innerSegmentCrossesInnerTriangle(triA_vrt2_ind,triA_vrt0_ind,
                            triB_vrt0_ind,triB_vrt1_ind,triB_vrt2_ind, mesh)   )
       return 1;

    return 0;
}

uint32_t vrt_signe_orient3d(uint32_t vrt1, uint32_t vrt2,
                            uint32_t vrt3, uint32_t vrt4, const TetMesh* mesh){

    if(vrt1==vrt2 || vrt1==vrt3 || vrt1==vrt4 ||
       vrt2==vrt3 || vrt2==vrt4 || vrt3==vrt4   ) return 0;

    const double* vrt1_coord = mesh->vertices[vrt1].coord;
    const double* vrt2_coord = mesh->vertices[vrt2].coord;
    const double* vrt3_coord = mesh->vertices[vrt3].coord;
    const double* vrt4_coord = mesh->vertices[vrt4].coord;
    return signe_orient3d(vrt1_coord, vrt2_coord, vrt3_coord, vrt4_coord);
}

// It is assumed that test verices and segment endpoints (that define a
// straight line) lie on the same plane.
static inline uint32_t vrt_same_half_plane(uint32_t test_vrt0, uint32_t test_vrt1,
                                uint32_t strLine4_vrt0, uint32_t strLine4_vrt1,
                                            const TetMesh* mesh){

  if(test_vrt0==strLine4_vrt0 || test_vrt0==strLine4_vrt1 ||
     test_vrt1==strLine4_vrt0 || test_vrt1==strLine4_vrt1   ) return 0;

  const double* test0_coord = mesh->vertices[test_vrt0].coord;
  const double* test1_coord = mesh->vertices[test_vrt1].coord;
  const double* sL0_coord = mesh->vertices[strLine4_vrt0].coord;
  const double* sL1_coord = mesh->vertices[strLine4_vrt1].coord;
  return (same_half_plane(test0_coord, test1_coord, sL0_coord, sL1_coord));
}

//----------------------------
// Ad-hok geometric Predicates
//----------------------------

static inline uint32_t tetvrtONconstraintplane(uint32_t vrt,
                                               const uint32_t* c_vrts,
                                               const TetMesh* mesh){
  return (vrt_signe_orient3d(c_vrts[0],c_vrts[1],c_vrts[2], vrt, mesh) == 0);
}

static inline uint32_t tetedgeONconstraintplane(const uint32_t* te_vrts,
                                                const uint32_t* c_vrts,
                                                const TetMesh* mesh){

  return (tetvrtONconstraintplane(te_vrts[0], c_vrts, mesh) &&
          tetvrtONconstraintplane(te_vrts[1], c_vrts, mesh)    );
}

static inline uint32_t tetfaceONconstraintplane(const uint32_t* tf_vrts,
                                                const uint32_t* c_vrts,
                                                const TetMesh* mesh){

  return (tetvrtONconstraintplane(tf_vrts[0], c_vrts, mesh) &&
          tetvrtONconstraintplane(tf_vrts[1], c_vrts, mesh) &&
          tetvrtONconstraintplane(tf_vrts[2], c_vrts, mesh)    );
}

//  Input: the indices of two vertices: test0, test1,
//         the indices of three verttices that define a plane: p0, p1, p2,
//         pointer to mesh: mesh.
// Output: return 1 if test0 and test1 have the same orient3d wrt the plane for
//         p0,p1,p2. (i.e. they belong to the same half-space).
// Note. if test0 or test1 belong to the plane for p0,p1,p2 the function returns
//       1 if both test vertices belong to the plane.
static inline uint32_t vrtsInSameHalfSpace(uint32_t test0, uint32_t test1,
                                          uint32_t p0, uint32_t p1, uint32_t p2,
                                          const TetMesh* mesh){
  return (vrt_signe_orient3d(p0,p1,p2, test0, mesh) ==
          vrt_signe_orient3d(p0,p1,p2, test1, mesh)    );
}

//-----------------
// Mesh Exploration
//-----------------

// Returns 1 if an unsigned int on 32 bit (vertex index type) IS AN ELEMENT OF
// an array of unsigned int on 32 bit,
// otherwise returns 0.
static inline bool arrayUINT32_contains_elem(const uint32_t* array,
                                                uint32_t lenght, uint32_t elem){
    for(uint32_t i=0; i<lenght; i++)
        if(array[i]==elem)  return true;

    return false;
}

//  Input: index of the vertex vrt: vrt,
//         ID of the tetrahedron tet: tet_ID (4 * tetrahedron_index),
//         pointer to the mesh.
// Output: returns 1 if vrt is a vertex of tet,
//         0 otherwise.
static inline bool is_vertex_ofTet(uint32_t vrt, uint64_t tet_ID,
                                      const TetMesh* mesh){
    return (mesh->tet_node[tet_ID  ] == vrt ||
            mesh->tet_node[tet_ID+1] == vrt ||
            mesh->tet_node[tet_ID+2] == vrt ||
            mesh->tet_node[tet_ID+3] == vrt   );
}

//  Input: endponit of a segment: v0,v1,
//         three vertices of a triangle: tri_vrts.
// Output: true if <v0,v1> is a side of the triangle, false otherwise.
bool segment_is_triSide(uint32_t v0, uint32_t v1, const uint32_t* tri_vrts){
  return (arrayUINT32_contains_elem(tri_vrts,3, v0) &&
          arrayUINT32_contains_elem(tri_vrts,3, v1)     );

}

//  Input: endponit of a segment: v0,v1,
//         three vertices of a triangle: tri_vrts.
// Output: return 1 if <v0,v1> is a side of the triangle or is included in a
//         side of the triangle, 0 otherwise.
//         if 1 is returned, by using opp_tri_vrt return the index of the vertex
//         of the triangle opposite to the side that contains <v0,v1>
bool segment_in_triSide(uint32_t v0, uint32_t v1, const uint32_t* tri_vrts,
                            uint32_t* opp_tri_vrt, const TetMesh* mesh){
  if(vrt_pointInSegment(v0,tri_vrts[0],tri_vrts[1],mesh)){
    if(vrt_pointInSegment(v1,tri_vrts[0],tri_vrts[1],mesh)){
      *opp_tri_vrt = tri_vrts[2];
      return true;
    }
    else return false;
  }

  if(vrt_pointInSegment(v0,tri_vrts[1],tri_vrts[2],mesh)){
    if(vrt_pointInSegment(v1,tri_vrts[1],tri_vrts[2],mesh)){
      *opp_tri_vrt = tri_vrts[0];
      return true;
    }
    else return false;
  }

  if(vrt_pointInSegment(v0,tri_vrts[2],tri_vrts[0],mesh)){
    if(vrt_pointInSegment(v1,tri_vrts[2],tri_vrts[0],mesh)){
      *opp_tri_vrt = tri_vrts[1];
      return true;
    }
    else return false;
  }

  return false;

}

//  Input: triangle <v0,v1,v2>,
//         ID of the tetrahedron tet: tet_ID (4 * tetrahedron_index),
//         pointer to the mesh.
// Output: returns 1 if <v0,v1,v2> is a face of tet,
//         0 otherwise.
static inline bool triangle_is_tetFace(uint32_t v0, uint32_t v1, uint32_t v2,
                                           uint64_t tet_ID, const TetMesh* mesh){
    return (is_vertex_ofTet(v0, tet_ID, mesh) &&
        is_vertex_ofTet(v1, tet_ID, mesh) &&
        is_vertex_ofTet(v2, tet_ID, mesh));
}

//  Input: trinagle <v0,v1,v2>,
//         ID of the tetrahedron tet: tet_ID (4 * tetrahedron_index),
//         pointer to the mesh.
// Output: returns the face_ID (0,1,2,3) of the face of vertices v0,v1,v2
//         i.e. returns the vertex_ID of the vertex different from v0,v1,v2.
static inline uint32_t tet_faceID(uint32_t v0, uint32_t v1, uint32_t v2,
                                  uint64_t tet_ID, const TetMesh* mesh){
    uint32_t* v = mesh->tet_node + tet_ID;
    if(v[0] != v0 && v[0] != v1 && v[0] != v2 ) return 0;
    if(v[1] != v0 && v[1] != v1 && v[1] != v2 ) return 1;
    if(v[2] != v0 && v[2] != v1 && v[2] != v2 ) return 2;
    if(v[3] != v0 && v[3] != v1 && v[3] != v2 ) return 3;

    ip_error("tet_faceID: FATAL ERROR no faces of tet correspond to vertices.\n");
}

//  Input: trinagle <v0,v1,v2>,
//         pointer to the mesh,
//         pointer of face_ID type.
// Output: returns 1 if <v0,v1,v2> is a face of one tetrahedron between those
//         incident in v0, 0 otherwise.
//         by using tet_face_ind returns the tetrahedron-face ID (0,1,2,3) that
//         of the face equal to the constraint.
// Note. Check if the triangle is a face of a tetrahedron incident in v0.
bool triangle_in_VT(uint32_t v0, uint32_t v1, uint32_t v2,
                                      TetMesh* mesh, uint64_t* tet_face_ind){

    bool found = false;
    uint64_t num_incTet_v0 = 0, tet_ID;
    uint64_t* incTet_v0 = NULL;
    incTet_v0 = mesh->incident_tetrahedra(v0, &num_incTet_v0);

    for(uint64_t i=0; i<num_incTet_v0; i++){
        tet_ID = 4*incTet_v0[i];
        if(triangle_is_tetFace(v0, v1, v2, tet_ID, mesh) ){
          *tet_face_ind = tet_ID + tet_faceID(v0,v1,v2, tet_ID, mesh);
          found = true;
          break;
         }
    }

    free(incTet_v0);
    return found;
}

//  Input: pointer to the mesh,
//         index of tetrahedron tet: tet_ind,
//         pointer to array of indices of the vertices of a
//         face of tet: face_vrts.
// Output: returns the vertex of tet opposite to the face defined by face_vrts
//         (i.e. the other vertex of tet).
uint32_t opposite_vertex_face(const TetMesh* mesh,
                                  uint64_t tet_ind, const uint32_t* face_vrts){

    const uint32_t *vp = mesh->tet_node + 4 * tet_ind;
    uint32_t v;
    v = *(vp++);
    if(v!=face_vrts[0] && v!=face_vrts[1] && v!=face_vrts[2]) return v;
    v = *(vp++);
    if(v!=face_vrts[0] && v!=face_vrts[1] && v!=face_vrts[2]) return v;
    v = *(vp++);
    if(v!=face_vrts[0] && v!=face_vrts[1] && v!=face_vrts[2]) return v;
    v = *(vp++);
    if(v!=face_vrts[0] && v!=face_vrts[1] && v!=face_vrts[2]) return v;
    return UINT32_MAX;
}

//  Input: pointer to the mesh,
//         index of tetrahedron tet: tet_ind,
//         index of the vertex v of tet: v_ind,
//         pointer to vertex index type: ptr.
// Output: by using ptr returns the 3 vertices of the face of tet opposite to v
//        (i.e. the other 3 vertices of tet).
void opposite_face_vertices(const TetMesh* mesh, uint64_t tet_ind,
                                          uint32_t v_ind, uint32_t* ptr){

    uint32_t j=0, curr_v_ind;
    curr_v_ind = mesh->tet_node[4*tet_ind];
    if(curr_v_ind!=v_ind)   ptr[j++] = curr_v_ind;
    curr_v_ind = mesh->tet_node[4*tet_ind +1];
    if(curr_v_ind!=v_ind)   ptr[j++] = curr_v_ind;
    curr_v_ind = mesh->tet_node[4*tet_ind +2];
    if(curr_v_ind!=v_ind)   ptr[j++] = curr_v_ind;
    curr_v_ind = mesh->tet_node[4*tet_ind +3];
    if(curr_v_ind!=v_ind)   ptr[j++] = curr_v_ind;
}

//  Input: pointer to the mesh,
//         index of tetrahedron tet: tet_ind,
//         index of vertex v of tet: v_ind.
// Output: returns the index of the tetrahedron that
//         share with tet the face opposite to v.
uint64_t adjTet_oppsiteTo_vertex(const TetMesh* mesh,
                                               uint64_t tet_ind,
                                               uint32_t v_ind){
    if(mesh->tet_node[4*tet_ind] == v_ind)
      return mesh->tet_neigh[4*tet_ind]>>2;
    if(mesh->tet_node[4*tet_ind + 1] == v_ind)
      return mesh->tet_neigh[4*tet_ind + 1]>>2;
    if(mesh->tet_node[4*tet_ind + 2] == v_ind)
      return mesh->tet_neigh[4*tet_ind + 2]>>2;
    if(mesh->tet_node[4*tet_ind + 3] == v_ind)
      return mesh->tet_neigh[4*tet_ind + 3]>>2;
    return UINT64_MAX;
}

//  Input: pointer to the mesh,
//         index of tetrahedron tet: tet_ind,
//         index of two vertices u_inds[0],u_inds[1] of tet: u_inds,
//         pointer to vertex index type: ptr.
// Output: by using ptr returns the 2 vertices of the side of tet
//         opposite to (u1,u2) (i.e. the other 2 vertices of tet).
void opposite_side_vertices(const TetMesh* mesh, uint64_t tet_ind,
                                         const uint32_t* u_inds, uint32_t* ptr){
    uint32_t j=0, v;
    v = mesh->tet_node[4*tet_ind];
    if(v!=u_inds[0] && v!=u_inds[1]) ptr[j++] = v;
    v = mesh->tet_node[4*tet_ind +1];
    if(v!=u_inds[0] && v!=u_inds[1])  ptr[j++] = v;
    v = mesh->tet_node[4*tet_ind +2];
    if(v!=u_inds[0] && v!=u_inds[1])  ptr[j++] = v;
    v = mesh->tet_node[4*tet_ind +3];
    if(v!=u_inds[0] && v!=u_inds[1])  ptr[j++] = v;
}

//  Input: tetrahedron (tet) index: tet_ind,
//         endpoints of a segment: edge=(e0,e1),
//         pointer to vertex index type: n_ret_vrts,
//         pointer to vertex index type: ptr,
//         pointer to mesh.
// Output: by using ptr returns the list of the vertices of tet
//         that belongs to the closure of edge (at most 2),
//         by using n_ret_vrts returns the number of the vertices of tet
//         that belongs to the closure of edge.
void TetVrts_InSegment(uint64_t tet_ind, const uint32_t* edge,
                                     uint32_t* n_ret_vrts, uint32_t* ptr,
                                     const TetMesh* mesh){

  uint32_t n_found=0, found[2];
  uint32_t tet_vrts[4];
  extract_tetVrts(tet_vrts, tet_ind, mesh);

  // Check if one of tet vertices is an endpoint of edge.
  if(arrayUINT32_contains_elem(tet_vrts,4, edge[0]) )  found[n_found++]=edge[0];
  if(arrayUINT32_contains_elem(tet_vrts,4, edge[1]) )  found[n_found++]=edge[1];

  if(n_found<2){
  // Check if one of tet vertices is inside the edge.
  uint32_t j_start=0;

    if(n_found==0)
      for(uint32_t i=0; i<4; i++)
        if(vrt_pointInInnerSegment(tet_vrts[i],edge[0],edge[1],mesh)){
          found[n_found++] = tet_vrts[i];
          j_start=i+1;
          break;
        }

   if(n_found==1)
      for(uint32_t j=j_start; j<4; j++)
        if(vrt_pointInInnerSegment(tet_vrts[j],edge[0],edge[1],mesh)){
          found[n_found++] = tet_vrts[j];
          break;
        }
  }

  // Returns:
  *n_ret_vrts = n_found;
  if(n_found==0)  return;
  // Here one vertex has been found for shure.
  ptr[0] = found[0];
  if(n_found==2)  ptr[1] = found[1];
}


//------------------
// ARRAY COMPOSITION
//------------------

//  Input: array of tetrahedra indices: arrayOfTetsA,
//         numbero of elements of arrayOofTetsA: lenghtA,
//         array of tetrahedra indices: arrayOfTetsB,
//         numbero of elements of arrayOofTetsB: lenghtB.
// Output: by using arrayOfTetsB returns the new array composed by the elements
//         of arrayOfTetsB and then those of arrayOfTetsA,
//         by using lenghtB returns the number of elemnets in the new array.
void enqueueTetsArray(
                                const uint64_t* arrayOfTetsA, uint64_t lenghtA,
                                uint64_t** arrayOfTetsB, uint64_t* lenghtB){

    if( *lenghtB == 0 ){
        *lenghtB = lenghtA;
        *arrayOfTetsB = (uint64_t*) malloc(sizeof(uint64_t) * (*lenghtB));
        for(uint64_t j=0; j<lenghtA; j++)
            (*arrayOfTetsB)[j] = arrayOfTetsA[j];
    }
    else{
        uint64_t tmp = *lenghtB;
        *lenghtB += lenghtA;
        *arrayOfTetsB = (uint64_t*) realloc(*arrayOfTetsB,
                                            sizeof(uint64_t) * (*lenghtB));
        for(uint64_t j=0; j<lenghtA; j++)
            (*arrayOfTetsB)[tmp + j] = arrayOfTetsA[j];
    }
}

/****************************************/
/** half-edges and virtual constraints **/
/****************************************/

//  Input: pointer to two half-edges: he1, he2.
// Output: returns 1 if he1 > he2 in lexicographic order (w.r.t. endpts),
//         returns -1 if he1 < he2 in lexicographic order (w.r.t. endpts),
//         returns 0 if he2 = he1 in lexicographic order (w.r.t. endpts).
int half_edges_compare(const void* void_he1, const void* void_he2){
  const half_edge_t* he1 = (half_edge_t*)void_he1;
  const half_edge_t* he2 = (half_edge_t*)void_he2;
  if( he1->endpts[0] > he2->endpts[0]) return 1;
  else if( he1->endpts[0] == he2->endpts[0]){
    if( he1->endpts[1] > he2->endpts[1]) return 1;
    else if( he1->endpts[1] == he2->endpts[1]) return 0;
    else return -1;
  }
  else return -1;
}

//  Input: pointer to constraints,
//         pointer to half-edges.
// Output: fills the array half_edges of struct half_edge_t, for each constraint
//         side one half edge is created by using side endpoints and the index
//         of the constraint.
void fill_half_edges(const constraints_t* constraints, half_edge_t* half_edges){
  for(uint32_t t=0; t<constraints->num_triangles; t++){
    uint32_t tri_id = 3*t; // t is the constraint index.
    uint32_t v0 = constraints->tri_vertices[tri_id  ];
    uint32_t v1 = constraints->tri_vertices[tri_id+1];
    uint32_t v2 = constraints->tri_vertices[tri_id+2];
    // for each half-edge is required that endpts[0] < endpts[1].
    if(v0 <= v1){
    half_edges[tri_id  ].endpts[0] = v0;
    half_edges[tri_id  ].endpts[1] = v1;
    }
    else{
      half_edges[tri_id  ].endpts[0] = v1;
      half_edges[tri_id  ].endpts[1] = v0;
    }
    half_edges[tri_id  ].tri_ind = t;

    if(v1 <= v2){
    half_edges[tri_id+1].endpts[0] = v1;
    half_edges[tri_id+1].endpts[1] = v2;
    }
    else{
      half_edges[tri_id+1].endpts[0] = v2;
      half_edges[tri_id+1].endpts[1] = v1;
    }
    half_edges[tri_id+1].tri_ind = t;

    if(v2 <= v0){
    half_edges[tri_id+2].endpts[0] = v2;
    half_edges[tri_id+2].endpts[1] = v0;
    }
    else{
      half_edges[tri_id+2].endpts[0] = v0;
      half_edges[tri_id+2].endpts[1] = v2;
    }
    half_edges[tri_id+2].tri_ind = t;
  }
}

//  Input: pointer to half-edges.
//         total number of half-edges.
// Output: by using half_edges return the array of half-edge sorted by
//         increasing lexicographic order.
void sort_half_edges(half_edge_t* half_edges, uint32_t num_half_edges){
  qsort(half_edges, num_half_edges, sizeof(half_edge_t), half_edges_compare);
}

//  Input: index of an half-edge w.r.t. the array half_edges: he,
//         index of the vector constraints->tri_vertices in which insert the
//         3 new vertices of the virtual constraint: pos,
//         pointer to constraints,
//         pointer to half-edges: half_edges,
//         pointer to mesh.
// Output: modifing constraints add the three vertices of a virtual constraint
//         in position pos, pos+1, pos+2 of the vector constraints->tri_vertices
void add_virtual_constraint(uint32_t he, uint32_t pos,
                            constraints_t* constraints,
                            const half_edge_t* half_edges, const TetMesh* mesh){
  const uint32_t tri_id = 3*half_edges[he].tri_ind;
  const uint32_t v0 = constraints->tri_vertices[tri_id  ];
  const uint32_t v1 = constraints->tri_vertices[tri_id+1];
  const uint32_t v2 = constraints->tri_vertices[tri_id+2];
  // Find a vertex (u) of a tetrahedron incident in one edge endpoints
  // half_edges[he].endpts[0] (or half_edges[he].endpts[1])
  // not aligned with the vertices of triangle half_edges[he].tri_id
  const uint64_t tet_ind = mesh->vertices[v0].inc_tet; // non-ghost tet
  uint32_t u = mesh->tet_node[4*tet_ind];
  uint32_t k=0;

  while(vrt_signe_orient3d(u, v0, v1, v2, mesh) == 0 ){
      u = mesh->tet_node[4*tet_ind + (++k)];
  }

  // Create the new (virtual) constraint:
  // <u, half_edges[he].endpts[0], half_edges[he].endpts[1]>
  // and insert its vertices in the constraints vertices vector.
  constraints->tri_vertices[pos  ] = half_edges[he].endpts[0];
  constraints->tri_vertices[pos+1] = half_edges[he].endpts[1];
  constraints->tri_vertices[pos+2] = u;
}

//  Input: first and last indices of two half-edges such that all the
//         consecutive half-edges from them have the sam endpoints: he0, he1,
//         pointer to half-edges,
//         pointer to constraints,
//         pointer to mesh.
// Output: return 1 if all constraints incident on the common edge are coplanar
//         and on the same side of the halfedge,
//         0 otherwise.
bool tri_onSameEdge_allCoPlanar(uint32_t he0, uint32_t he1,
                                    const half_edge_t* half_edges,
                                    const constraints_t* constraints,
                                    const TetMesh* mesh){
  // Find the vertex (u) of constraint half_edges[he].tri_ind different
  // from edge endpoints.
  uint32_t v0, v1, v2;
  uint32_t tri_id = 3*half_edges[he0].tri_ind;
  uint32_t u = constraints->tri_vertices[tri_id];
  uint32_t k=0;
  const uint32_t hp0 = half_edges[he0].endpts[0];
  const uint32_t hp1 = half_edges[he0].endpts[1];
  while(u==hp0 || u==hp1)
    u = constraints->tri_vertices[tri_id + (++k)];

  for(uint32_t he = he0+1; he<=he1; he++){
    tri_id = 3*half_edges[he].tri_ind;
    v0 = constraints->tri_vertices[tri_id  ];
    v1 = constraints->tri_vertices[tri_id+1];
    v2 = constraints->tri_vertices[tri_id+2];
    if(vrt_signe_orient3d(u, v0, v1, v2, mesh) != 0 ) return false;
  }

  // Assume that triangles are not degenerate
  // Calculate dominating normal
  // Return TRUE only if all the coplanar triangles are on the same side of the edge

  tri_id = 3 * half_edges[he0].tri_ind;
  v0 = constraints->tri_vertices[tri_id];
  v1 = constraints->tri_vertices[tri_id + 1];
  v2 = constraints->tri_vertices[tri_id + 2];
  const double* v0c = mesh->vertices[v0].coord;
  const double* v1c = mesh->vertices[v1].coord;
  const double* v2c = mesh->vertices[v2].coord;
  const int dom = genericPoint::maxComponentInTriangleNormal(v0c[0], v0c[1], v0c[2], v1c[0], v1c[1], v1c[2], v2c[0], v2c[1], v2c[2]);

  double e0c[2], e1c[2], oc[2];

  int i, j;
  for (i = j = 0; i < 3; i++) if (i != dom)
  {
      e0c[j] = mesh->vertices[hp0].coord[i];
      e1c[j] = mesh->vertices[hp1].coord[i];
      oc[j] = mesh->vertices[u].coord[i];
      j++;
  }

  const double base_or = orient2d(e0c, e1c, oc);
  const int base_or_sign = (base_or > 0) - (base_or < 0);

  for (uint32_t he = he0 + 1; he <= he1; he++) {
      tri_id = 3 * half_edges[he].tri_ind;
      v0 = constraints->tri_vertices[tri_id];
      v1 = constraints->tri_vertices[tri_id + 1];
      v2 = constraints->tri_vertices[tri_id + 2];
      u = v0;
      if (hp0 != v1 && hp1 != v1) u = v1;
      else if (hp0 != v2 && hp1 != v2) u = v2;
      for (i = j = 0; i < 3; i++) if (i != dom)
          oc[j++] = mesh->vertices[u].coord[i];
      const double new_or = orient2d(e0c, e1c, oc);
      const int new_or_sign = (new_or > 0) - (new_or < 0);
      if (new_or_sign != base_or_sign) return false;
  }

  return true;
}

//  Input: pointer to mesh,
//         pointer to constraints,
//         pointer to half edges.
// Output: by using constraints it adds needed vistual constraints to
//         constraints structure in order to ensure that during BSP subdivision
//         all constraints edges will be inserted.
uint32_t place_virtual_constraints(TetMesh* mesh, constraints_t* constraints,
                               half_edge_t* half_edges){
    // 1- Count needed virtual constraints.
    const uint32_t num_half_edges = 3*constraints->num_triangles;
    // Mark half-edges which requires a virtual constraint:
    // - those that have only one incident constraint,
    // - those that have all incident constraints coplanar.
    uint32_t* need_virtual_constraint = (uint32_t*) calloc(num_half_edges, sizeof(uint32_t));
    uint32_t num_virtual_constraints = 0;

    uint32_t he=0;    // half-edge index.
    while(he < num_half_edges -1){

      uint32_t onSameEdge = 1; // Counts how many constraints are
                               // incident on current half-edge.
      while(he + onSameEdge < num_half_edges &&
          half_edges_compare(&half_edges[he], &half_edges[he+onSameEdge]) == 0 )
             onSameEdge++;

      if(onSameEdge == 1 ||
         tri_onSameEdge_allCoPlanar(he, he+onSameEdge-1, half_edges, constraints, mesh) ){
        num_virtual_constraints++;
        need_virtual_constraint[he] = 1;
      }

      he += onSameEdge; // next half-edge with different endpoints.
    }

    // 2- Add virtual constraints.

    // Resize constraints_vrts
    uint32_t pos = 3*constraints->num_triangles;
    constraints->num_virtual_triangles = num_virtual_constraints;
    constraints->num_triangles += num_virtual_constraints;
    constraints->tri_vertices = (uint32_t*) realloc(constraints->tri_vertices,
                                3*constraints->num_triangles*sizeof(uint32_t));
    constraints->constr_group = (CONSTR_GROUP_T*) realloc(constraints->constr_group,
                                constraints->num_triangles*sizeof(CONSTR_GROUP_T));

    // Fill constraints->tri_vertices with virtual constraints vertices
    for(uint32_t he=0; he<num_half_edges; he++)
      if(need_virtual_constraint[he]==1){
        add_virtual_constraint(he, pos, constraints, half_edges, mesh);
        pos += 3;
      }

    free(need_virtual_constraint);

    return num_virtual_constraints;
}


/*************************************************************/
/** search intersections between constraints and tetrahedra **/
/*************************************************************/

//----------------------------------------------------------------------
// FUNCTIONS TO CLASSIFY GENERIC INTERSECTIONS INTO: PROPER AND IMPROPER
//----------------------------------------------------------------------

//  Input: the 3 vertices of the contraint-triangle: c[3],
//         the 3 vertices of the a tetrahedron face: t[3],
//         pointer to mesh.
// Output: returns 2 if the intersection is the whole tet-face (proper),
//                 1 if it is improper,
//                 0 if it is proper but not the whole face.
//         otherwise UINT32_MAX.
// Note. it is assumed that the vertices of tet-face lie on the constraint-plane.
uint32_t intersectionClass_tetFace_constr(const uint32_t* c, const uint32_t* t,
                                           const TetMesh* mesh){

  // Since the constraint-triangle is convex:
  // IF all the 3 vertices of tet-face belong to the constraint (with boundary)
  // THEN the whole tet-face is included in the constraint -> class. PROPER
  if(vrt_pointInTriangle(t[0], c[0],c[1],c[2], mesh) &&
     vrt_pointInTriangle(t[1], c[0],c[1],c[2], mesh) &&
     vrt_pointInTriangle(t[2], c[0],c[1],c[2], mesh)    )   return 2;

  // If there is a tet-edge contraint-side inner corissing -> class. IMPROPER
  if(vrt_innerSegmentsCross(t[0],t[1], c[0],c[1], mesh))   return 1;
  if(vrt_innerSegmentsCross(t[0],t[1], c[1],c[2], mesh))   return 1;
  if(vrt_innerSegmentsCross(t[0],t[1], c[2],c[0], mesh))   return 1;

  if(vrt_innerSegmentsCross(t[1],t[2], c[0],c[1], mesh))   return 1;
  if(vrt_innerSegmentsCross(t[1],t[2], c[1],c[2], mesh))   return 1;
  if(vrt_innerSegmentsCross(t[1],t[2], c[2],c[0], mesh))   return 1;

  if(vrt_innerSegmentsCross(t[2],t[0], c[0],c[1], mesh))   return 1;
  if(vrt_innerSegmentsCross(t[2],t[0], c[1],c[2], mesh))   return 1;
  if(vrt_innerSegmentsCross(t[2],t[0], c[2],c[0], mesh))   return 1;

  // Reaming tetrahedron are such that at least one vertex is outside the
  // constraint and no tet-edge inner cross any constraint-side.
  //
  // An intersection exists (with this triangle since the tet-vertex opoosite
  // to <tF_vrts[0],tF_vrts[1],tF_vrts[2]>) cannot lie on the constraint-plane.
  //
  // Any constraint vertex can belong to the tet-face (included boundary),
  // unless it concide with a tet-vertex.
  //
  // If a tet-vertex (a) is outside the constraint it must exists at least
  // another tet-vertex (b) on the constraint boundary or in the interior.
  // [ (b) cannot be in the interior, on the countrary <(a),(b)> will inner
  //   cross a constraint side, or a constraint vertex must lye on <(a),(b)>  ]
  // The same holds for the other tet-vert (different from (a) and (b)).

  // Otherwise...
  // the intersection is only a tet-edge or a tet-vertex -> class. PROPER
  return 0;
}

//  Input: 3 vertices of a contraint-triangle: c_vrts[3],
//         2 vertices of a tetrahedron edge: tE_vrts[2],
//         the other 2 vertices of the tetrahedron: opptE_vrts[2],
//         a flag (1 if opptE_vrts have same orient3d w.r.t. constraint-plane,
//         0 otherwise): opptE_same_or,
//         pointer to mesh.
// Output: returns 0 if the intersection is proper, 1 if it is improper,
//         otherwise UINT32_MAX.
// Note. It is assumed that the orient3d of opptE_vrts[0] and opptE_vrts[1] are
//       different from 0 and that the orient3d of tE_vrts[0], tE_vrts[1] is 0.
uint32_t intersectionClass_tetEdge_constr(const uint32_t* c_vrts,
                                           const uint32_t* tE_vrts,
                                           const uint32_t* opptE_vrts,
                                           uint32_t opptE_same_or,
                                           const TetMesh* mesh){

  const uint32_t tE0_in=vrt_pointInTriangle(tE_vrts[0], c_vrts[0],c_vrts[1],c_vrts[2], mesh);
  const uint32_t tE1_in=vrt_pointInTriangle(tE_vrts[1], c_vrts[0],c_vrts[1],c_vrts[2], mesh);

  if(opptE_same_or){
    // IF the 2 vertices of tet-edge vestices belong to
    // the constraint (boundary included) -> PROPER INTERSECTION.
    if(tE0_in && tE1_in)  return 0;
    // At least one vertices between tE_vrts[0] and tE_vrts[1] must stay in the
    // the constraint (otherwise there will be no intersection - IMPOSSIBLE)
    #ifdef DEBUG
    if(!tE0_in && !tE1_in)
      printf("[conforming_mesh.c]intersectionClass_tetEdge_constr: ERROR "
             "there is no intersection.\n");
    #endif

    // Note. a constraint vertex cannot stay in the interior of the edge
    //       <tE_vrts[0],tE_vrts[1]>.

    // IF <tE_vrts[0],tE_vrts[1]> inner crosses a constraint-side
    // THEN the intersection is IMPROPER.
    if(vrt_innerSegmentsCross(tE_vrts[0],tE_vrts[1],c_vrts[0],c_vrts[1],mesh) ||
       vrt_innerSegmentsCross(tE_vrts[0],tE_vrts[1],c_vrts[1],c_vrts[2],mesh) ||
       vrt_innerSegmentsCross(tE_vrts[0],tE_vrts[1],c_vrts[2],c_vrts[0],mesh)  )
        return 1;
    // The only possible disposition is that tE_vrts[0] or tE_vrts[1] lies on
    // the constraint boundary and the other one is outside the constraint:
    // PROPER INTERSECTION.
    return 0;
  }
  else{ // opptE_same_or == 0 ->
        // <opptE_vrts[0],opptE_vrts[1]> crosses constraint-plane

     // IF at least one between tE_vrts[0] and tE_vrts[1] lie in the constraint
     // interior (no boundary)
     // THEN the intersection is IMPROPER.
     if(vrt_pointInInnerTriangle(tE_vrts[0], c_vrts[0],c_vrts[1],c_vrts[2], mesh) ||
        vrt_pointInInnerTriangle(tE_vrts[1], c_vrts[0],c_vrts[1],c_vrts[2], mesh)   )
         return 1;
     // Note. neither tE_vrts[0] or tE_vrts[1] lies in the constraint interior.
     if(tE0_in && tE1_in){

       uint32_t cs01, cs12, cs20;
       cs01 = vrt_pointInSegment(tE_vrts[0], c_vrts[0],c_vrts[1], mesh) +
              vrt_pointInSegment(tE_vrts[1], c_vrts[0],c_vrts[1], mesh);
       cs12 = vrt_pointInSegment(tE_vrts[0], c_vrts[1],c_vrts[2], mesh) +
              vrt_pointInSegment(tE_vrts[1], c_vrts[1],c_vrts[2], mesh);
       cs20 = vrt_pointInSegment(tE_vrts[0], c_vrts[2],c_vrts[0], mesh) +
              vrt_pointInSegment(tE_vrts[1], c_vrts[2],c_vrts[0], mesh);

       // IF both tE_vrts[0] or tE_vrts[1] belong to the constraint boundary they
       // must lie on the same constraint side, otherwise IMPROPER intersection.
       if(cs01!=2 && cs12!=2 && cs20!=2) return 1;
       // Note. both tE_vrts[0] or tE_vrts[1] lie on the same constraint side.
       // IF the edge <opptE_vrts[0],opptE_vrts[1]> crosses the constraint
       // (boundary included)
       // THEN the intersection is IMPROPER.
       if(vrt_innerSegmentCrossesTriangle(opptE_vrts[0],opptE_vrts[1],c_vrts[0],c_vrts[1],c_vrts[2], mesh))
          return 1;
       // IF a constraint side (not aligned with  <tE_vrts[0],tE_vrts[1]>)
       // crosses the interior of
       // <opptE_vrts[0],opptE_vrts[1],tE_vrts[0]> or
       // <opptE_vrts[0],opptE_vrts[1],tE_vrts[1]>
       // THEN the intersection is IMPROPER.
       uint32_t comm_side0, comm_side1, not_comm;
       if(cs01==2) {comm_side0=c_vrts[0]; comm_side1=c_vrts[1]; not_comm=c_vrts[2];}
       else if(cs12==2) {comm_side0=c_vrts[1]; comm_side1=c_vrts[2]; not_comm=c_vrts[0];}
            else {comm_side0=c_vrts[2]; comm_side1=c_vrts[0]; not_comm=c_vrts[1];}

       if(vrt_innerSegmentCrossesInnerTriangle(comm_side0,not_comm, opptE_vrts[0],opptE_vrts[1],tE_vrts[0], mesh) ||
          vrt_innerSegmentCrossesInnerTriangle(comm_side1,not_comm, opptE_vrts[0],opptE_vrts[1],tE_vrts[0], mesh) ||
          vrt_innerSegmentCrossesInnerTriangle(comm_side0,not_comm, opptE_vrts[0],opptE_vrts[1],tE_vrts[1], mesh)||
          vrt_innerSegmentCrossesInnerTriangle(comm_side1,not_comm, opptE_vrts[0],opptE_vrts[1],tE_vrts[1], mesh)   )
            return 1;
       // Otherwise.. the intersection is PROPER
       return 0;
     }
     // Note. at least one between tE_vrts[0] and tE_vrts[1] belong to the
     // constraint boundary.
     if(tE0_in || tE1_in){
       // IF <tE_vrts[0],tE_vrts[1]> inner crosses a constraint-side
       // THEN the intersection is IMPROPER.
       if(vrt_innerSegmentsCross(tE_vrts[0],tE_vrts[1],c_vrts[0],c_vrts[1],mesh) ||
          vrt_innerSegmentsCross(tE_vrts[0],tE_vrts[1],c_vrts[1],c_vrts[2],mesh) ||
          vrt_innerSegmentsCross(tE_vrts[0],tE_vrts[1],c_vrts[2],c_vrts[0],mesh)  )
           return 1;
       // IF the edge <opptE_vrts[0],opptE_vrts[1]> crosses the constraint
       // (boundary included)
       // THEN the intersection is IMPROPER.
       if(vrt_innerSegmentCrossesTriangle(opptE_vrts[0],opptE_vrts[1],c_vrts[0],c_vrts[1],c_vrts[2], mesh))
          return 1;
       // IF a constraint side crosses the interior of
       // <opptE_vrts[0],opptE_vrts[1],tE_vrts[0]> or
       // <opptE_vrts[0],opptE_vrts[1],tE_vrts[1]>
       // THEN the intersection is IMPROPER.
       if(vrt_innerSegmentCrossesInnerTriangle(c_vrts[0],c_vrts[1], opptE_vrts[0],opptE_vrts[1],tE_vrts[0], mesh) ||
          vrt_innerSegmentCrossesInnerTriangle(c_vrts[1],c_vrts[2], opptE_vrts[0],opptE_vrts[1],tE_vrts[0], mesh) ||
          vrt_innerSegmentCrossesInnerTriangle(c_vrts[2],c_vrts[0], opptE_vrts[0],opptE_vrts[1],tE_vrts[0], mesh) ||
          vrt_innerSegmentCrossesInnerTriangle(c_vrts[0],c_vrts[1], opptE_vrts[0],opptE_vrts[1],tE_vrts[1], mesh) ||
          vrt_innerSegmentCrossesInnerTriangle(c_vrts[1],c_vrts[2], opptE_vrts[0],opptE_vrts[1],tE_vrts[1], mesh) ||
          vrt_innerSegmentCrossesInnerTriangle(c_vrts[2],c_vrts[0], opptE_vrts[0],opptE_vrts[1],tE_vrts[1], mesh)   )
           return 1;
       // Otherwise.. the intersection is PROPER
       return 0;
     }
     // OTHERWISE since an intersection exists it must be IMPROPER.
     return 1;
  }
}

//  Input: the 3 vertices of the contraint-triangle: c_vrts[3],
//         the vertex of the a tetrahedron : tV,
//         the 3 vertices of the tetrahedron face opposite to tV: opptF_vrts[3],
//         the 3 orient3d of the vertices of the tetrahedron face opposite to tV
//         w.r.t. the constraint_plane: or_opptF_vrts[3],
//         a flag (1 if opptF_vrts have same orient3d w.r.t. constraint-plane,
//         0 otherwise): opptF_same_or,
//         pointer to mesh.
// Output: returns 0 if the intersection is proper, 1 if it is improper,
//         otherwise UINT32_MAX.
// Note. It is assumed that the orient3d of opptF_vrts[0] opptF_vrts[1] and
//       opptE_vrts[2] are different from 0 and that the orient3d of tV is 0.
uint32_t intersectionClass_tetVrt_constr(const uint32_t* c_vrts, uint32_t tV,
                                         const uint32_t* opptF_vrts,
                                         const uint32_t* or_opptF_vrts,
                                          uint32_t opptF_same_or,

                                          const TetMesh* mesh){



  // IF <opptF_vrts[0],opptF_vrts[1],opptF_vrts[2]> do not croesses the
  // constraint-plane the only possible intersection is tV itself
  // THEN the intersection is PROPER.
  if(opptF_same_or) return 0;
  else{
    // IF tv is not in the constraint (boundary included)
    // THEN the intersection is IMPROPER
    if(!vrt_pointInTriangle(tV, c_vrts[0],c_vrts[1],c_vrts[2], mesh))
      return 1;
    // IF tv is in the constraint interior
    // THEN the intersection is IMPROPER
    if(vrt_pointInInnerTriangle(tV, c_vrts[0],c_vrts[1],c_vrts[2], mesh))
      return 1;

    // tV is on the constraint boundary

    // Find the disposition of the vertex of tetrahedron-face opposite to tV
    // w.r.t. the constraint.
    uint32_t tV_UNDER1, tV_UNDER2, tV_OVER;

    if(or_opptF_vrts[0] == or_opptF_vrts[1]){
        tV_UNDER1 = opptF_vrts[0];
        tV_UNDER2 = opptF_vrts[1];
        tV_OVER = opptF_vrts[2];
      }
      else if(or_opptF_vrts[1] == or_opptF_vrts[2]){
              tV_UNDER1 = opptF_vrts[1];
              tV_UNDER2 = opptF_vrts[2];
              tV_OVER = opptF_vrts[0];
            }
            else{
              tV_UNDER1 = opptF_vrts[2];
              tV_UNDER2 = opptF_vrts[0];
              tV_OVER = opptF_vrts[1];
            }

    // There is an IMPROPER improper in each of the following cases:
    //   (i) <tV_OVER, tV_UNDER1> crosses the constraint (boundary or interior)
    //  (ii) <tV_OVER, tV_UNDER2> crosses the constraint (boundary or interior)
    // (iii) the interior of <tV, tV_OVER, tV_UNDER1>
    //       is inner crossed by any constraint side,
    //  (iv) the interior of <tV, tV_OVER, tV_UNDER2>
    //       is inner crossed by any constraint side,
    //   (v) the interior of <tV_OVER, tV_UNDER1, tV_UNDER2>
    //       is inner crossed by any constraint side,


    // Condition (i)
    if(vrt_innerSegmentCrossesTriangle(tV_OVER,tV_UNDER1, c_vrts[0],c_vrts[1],c_vrts[2], mesh))
      return 1;

    // Condition (ii)
    if(vrt_innerSegmentCrossesTriangle(tV_OVER,tV_UNDER2, c_vrts[0],c_vrts[1],c_vrts[2], mesh))
      return 1;

    // Condition (iii)
    if(vrt_innerSegmentCrossesInnerTriangle(c_vrts[0], c_vrts[1], tV, tV_OVER, tV_UNDER1, mesh) ||
       vrt_innerSegmentCrossesInnerTriangle(c_vrts[1], c_vrts[2], tV, tV_OVER, tV_UNDER1, mesh) ||
       vrt_innerSegmentCrossesInnerTriangle(c_vrts[2], c_vrts[0], tV, tV_OVER, tV_UNDER1, mesh)   )
      return 1;

    // Condition (iv)
    if(vrt_innerSegmentCrossesInnerTriangle(c_vrts[0], c_vrts[1], tV, tV_OVER, tV_UNDER2, mesh) ||
       vrt_innerSegmentCrossesInnerTriangle(c_vrts[1], c_vrts[2], tV, tV_OVER, tV_UNDER2, mesh) ||
       vrt_innerSegmentCrossesInnerTriangle(c_vrts[2], c_vrts[0], tV, tV_OVER, tV_UNDER2, mesh)   )
      return 1;

    // Condition (v)
    if(vrt_innerSegmentCrossesInnerTriangle(c_vrts[0], c_vrts[1], tV_OVER, tV_UNDER1, tV_UNDER2, mesh) ||
       vrt_innerSegmentCrossesInnerTriangle(c_vrts[1], c_vrts[2], tV_OVER, tV_UNDER1, tV_UNDER2, mesh) ||
       vrt_innerSegmentCrossesInnerTriangle(c_vrts[2], c_vrts[0], tV_OVER, tV_UNDER1, tV_UNDER2, mesh)   )
      return 1;

    // Otherwise...  PROPER intersection.
    return 0;
  }

}

//  Input: 3 vertices of the contraint-triangle: constr_v[3],
//         array of indices of the tetrahedra that intersects
//         the constraint: tets,
//         number of tetrahedra that intersects the constraint: num_tets,
//         array of marker for the intersections
//         tetrahedra-constraint: mark_TetIntersection,
//         pointer to mesh.
// Output: by using mark_TetIntersection marks the improper
//         intersections tetrahedra-constraint.
// Note. It is assumed that each tetrahedron of tet_list intersects
//       (properly or improperly) the constraint.
void find_improperIntersection(const uint32_t* constr_v,
                               const uint64_t* tets,
                               uint64_t num_tets,
                               uint32_t* mark_TetIntersection,
                               const TetMesh* mesh){
  for(uint64_t tet_i=0; tet_i<num_tets; tet_i++){
      uint64_t tet_ind = tets[tet_i]; // index of tetrahedron tet
      uint64_t tet_ID = 4*tet_ind;

      // tet vertices
      uint32_t tet_v[4];
      extract_tetVrts(tet_v, tet_ind, mesh);

      // tet vertices orient w.r.t. constraint-plane.
      int or_tet_v[4];
      or_tet_v[0] = vrt_signe_orient3d(constr_v[0], constr_v[1], constr_v[2], tet_v[0], mesh);
      or_tet_v[1] = vrt_signe_orient3d(constr_v[0], constr_v[1], constr_v[2], tet_v[1], mesh);
      or_tet_v[2] = vrt_signe_orient3d(constr_v[0], constr_v[1], constr_v[2], tet_v[2], mesh);
      if(tet_v[3] == INFINITE_VERTEX) or_tet_v[3] = 3;    // Meaningless value
      else or_tet_v[3] = vrt_signe_orient3d(constr_v[0], constr_v[1], constr_v[2], tet_v[3], mesh);

      // Usefull strucutures for subsimplexs of tet.
      uint32_t tetFace_v[3];
      uint32_t tetEdge_v[2];
      uint32_t oppTetEdge_v[2];
      uint32_t tet_u;
      uint32_t oppTetFace_v[3];

      uint32_t iclass;

      // ghost-tet case
      if(tet_v[3] == INFINITE_VERTEX) continue;

      // (non-ghost) tet case

      // Coplanarity measure
      uint32_t num_coplan_v = 0;
      for(uint32_t i=0; i<4; i++)  if(or_tet_v[i] == 0) num_coplan_v++;

      // Case: a tet-face lies on the constraint plane.
      if(num_coplan_v == 3){ // For BSP cuts mapping these kind of intersections are usless,
                             // by the way they are used to split BSPfaces that partially
                             // (2D)overlap with a constraint. So, we map them with extra symbols.
        uint32_t k = 0;
        for(uint32_t i=0; i<4; i++)
          if(or_tet_v[i] == 0)    tetFace_v[k++] = tet_v[i];
          else tet_u = tet_v[i];

        iclass = intersectionClass_tetFace_constr(constr_v, tetFace_v, mesh);
        if(iclass == 1){ // IMPROPER INTERSECTION -> not interior, but face 2Doverlap.
          const uint32_t f = tet_faceID(tetFace_v[0], tetFace_v[1], tetFace_v[2], tet_ID, mesh);
          mark_TetIntersection[tet_ind] = 10 + f; // See MACRO at the beginning.
        }
        else if(iclass == 2){ // whole face PROPER INTERSECTION -> face 2Doverlap.
          const uint32_t f = tet_faceID(tetFace_v[0], tetFace_v[1], tetFace_v[2], tet_ID, mesh);
          mark_TetIntersection[tet_ind] = 10 + f; // See MACRO at the beginning.
        }

        continue; // No further intersection are possible.
      }

      // Assumption: no tet-face is aligned with the constraint.

      if(num_coplan_v==2){
        uint32_t k=0, j=0, same_or=0;
        int non_zero_or[2];
        for(uint32_t i=0; i<4; i++)
          if(or_tet_v[i] == 0) tetEdge_v[k++] = tet_v[i];
          else{
            oppTetEdge_v[j] = tet_v[i];
            non_zero_or[j++] = or_tet_v[i];
          }

        // MARCO: if tet is tangent to constraint no split will occur in any case
        if (non_zero_or[0] == non_zero_or[1]) continue; // same_or = 1; <------------------------

        iclass = intersectionClass_tetEdge_constr(constr_v, tetEdge_v, oppTetEdge_v, same_or, mesh);
        if(iclass == 1) mark_TetIntersection[tet_ind] = IMPROPER_INTERSECTION;

        continue; // No further intersection are possible.
      }

      // Assumption: no tet-face or tet-edge is aligned with the constraint.

      if(num_coplan_v == 1){
        uint32_t k=0, same_or=0;
        uint32_t or_oppTetFace_v[3];
        for(uint32_t i=0; i<4; i++)
          if(or_tet_v[i] == 0) tet_u = tet_v[i];
          else{
            oppTetFace_v[k] = tet_v[i];
            or_oppTetFace_v[k++] = or_tet_v[i];
          }

        // As above. If tet is tangent to constraint no split will occur in any case
        if (oppTetFace_v[0] == oppTetFace_v[1] &&
            oppTetFace_v[0] == oppTetFace_v[2]    ) continue; // same_or = 1; <-----------------------

        iclass = intersectionClass_tetVrt_constr(constr_v, tet_u, oppTetFace_v, or_oppTetFace_v, same_or, mesh);
        if(iclass == 1) mark_TetIntersection[tet_ind] = IMPROPER_INTERSECTION;

        continue; // No further intersection are possible.
      }

      // No tet-vertrices on the constraint-plane -> the intersection is improper.
      mark_TetIntersection[tet_ind] = IMPROPER_INTERSECTION;

  }
}

//---------------------------------------
// FUNCTIONS TO FIND PROPER INTERSECTIONS
//---------------------------------------

//  Input: the 3 vertices of the contraint-triangle:
//         c_vrts[0], c_vrts[1], c_vrts[2],
//         the 3 vertices of the a tetrahedron face:
//         tF_vrts[0], tF_vrts[1], tF_vrts[2],
//         pointer to mesh.
// Output: returns 1 if there is a proper inetrsection (of the whole face),
//         0 otherwise.
uint32_t properIntersection_tetFace_constr(const uint32_t* c_vrts,
                                           const uint32_t* tF_vrts,
                                           const TetMesh* mesh){

  // Necessary Cond. -> the 3 vertices of tet-face lies on the constraint plane.
  if(!tetfaceONconstraintplane(tF_vrts, c_vrts, mesh))  return 0;

  // Since the constraint-triangle is convex:
  // Sufficient Cond. -> IF all the 3 vertices of tet-face belong
  //                     to the constraint (boundary included)
  //                     THEN the whole tet-face is included in
  //                     the constraint                 -> proper intersection.
  if(vrt_pointInTriangle(tF_vrts[0], c_vrts[0],c_vrts[1],c_vrts[2], mesh) &&
     vrt_pointInTriangle(tF_vrts[1], c_vrts[0],c_vrts[1],c_vrts[2], mesh) &&
     vrt_pointInTriangle(tF_vrts[2], c_vrts[0],c_vrts[1],c_vrts[2], mesh)    )
     return 1;

  // Otherwise...
  //(No proper intersection between the whole face and the constraint)
  return 0;
}

//  Input: 3 vertices of a contraint-triangle: c_vrts[3],
//         2 vertices of a tetrahedron edge: tE_vrts[2],
//         the other 2 vertices of the tetrahedron: opptE_vrts[2],
//         pointer to mesh.
// Output: returns 1 if there is a proper inetrsection (of the whole edge),
//         0 otherwise.
// Note. It is assumed that there are no faces of the tetrahedron
//       that proper intersect the constraint plane.
// Note. In the case of a ghost-tet, it is assumed that
//       ghost vertex is opptE_vrts[1].
uint32_t properIntersection_tetEdge_constr(const uint32_t* c_vrts,
                                           const uint32_t* tE_vrts,
                                           const uint32_t* opptE_vrts,
                                           const TetMesh* mesh){

  // Necessary Cond. -> the 2 vertices of tet-edge lies on the constraint plane.
  if(!tetedgeONconstraintplane(tE_vrts, c_vrts, mesh))  return 0;

  // Necessary Cond. -> the 2 vertices of tet-edge vestices belong to
  //                    the constraint (boundary included).
  if(!vrt_pointInTriangle(tE_vrts[0], c_vrts[0],c_vrts[1],c_vrts[2], mesh) ||
     !vrt_pointInTriangle(tE_vrts[1], c_vrts[0],c_vrts[1],c_vrts[2], mesh)    )
     return 0;

  // ghost-tet: the constraint can not cross the face opposite to ghost-tet,
  //            since it is a face of the convex-hull of the mesh.
  if(opptE_vrts[1] == UINT32_MAX) return 1;
  // LEAK -> if the hull is locally flat it must be added a further case that
  // for current purpopsal of the algorithm is ignored.

  uint32_t or_opptE0 =
         vrt_signe_orient3d(c_vrts[0],c_vrts[1],c_vrts[2], opptE_vrts[0], mesh);
  uint32_t or_opptE1 =
      vrt_signe_orient3d(c_vrts[0],c_vrts[1],c_vrts[2], opptE_vrts[1], mesh);

  // Co-planar tet-face case:
  // IF one of opptE_vrts (cop_opp_vrt) lie on the constraint-plane
  // AND
  // IF both tE_vrts belong to one of the constraint sides
  // AND
  // IF cop_opp_vrt is in the same half-plane of the other constraint vrt w.r.t.
  //    the aligned edge-side.
  // THEN the constraint-tetrahedron instersection is the edge
  //      <tE_vrts[0],tE_vrts[1]>                        -> proper intersection.
  if(or_opptE0==0){
     uint32_t opp_c_vrt;
     if(segment_in_triSide(tE_vrts[0],tE_vrts[1], c_vrts, &opp_c_vrt, mesh))
       if(!vrt_same_half_plane(opptE_vrts[0],opp_c_vrt,
                                tE_vrts[0],tE_vrts[1], mesh)) return 1;
     // otherwise there is no proper intersection of dimension 1..
     return 0;
  }
  if(or_opptE1==0){
     uint32_t opp_c_vrt;
     if(segment_in_triSide(tE_vrts[0],tE_vrts[1], c_vrts, &opp_c_vrt, mesh))
        if(!vrt_same_half_plane(opptE_vrts[1],opp_c_vrt,
                                tE_vrts[0],tE_vrts[1], mesh)) return 1;
     // otherwise there is no proper intersection of dimension 1..
     return 0;
  }

  // Now we can assume that no tet-faces lie on the constraint-plane:
  // neither opptE_vrts[0] and opptE_vrts[1] can lie on the constraint-plane.

  // Sufficient Cond. -> IF both endpoints of opptE_vrts stay over (or under)
  //                     the constraint plane
  //                     THEN the only intersection in with tetEdge
  //                                                     -> proper intersection.
  if( or_opptE0 == or_opptE1 )    return 1;

  // Remaining tetrahedron-constraint disposition are such that:
  // - the edge tE_vrts belong to the constraint (boundary included),
  // - the edge opptE_vrts crosses the plane of the constraint.

  // Sufficient Cond. ->  IF both the conditions hold:
  //  (i) <tE_vrts[0],tE_vrts[1]> is aligned with a side of the constraint,
  // (ii) the constraint vertex not aligned with <tE_vrts[0], tE_vrts[1]>
  //      and opptE_vrts[1] do not stay in the half-space, defined by
  //      the plane for tE_vrts[0], tE_vrts[1], opptE_vrts[0].
  //                      THEN the tetrahedron intersects the constraint
  //                      only in <tE_vrts[0],tE_vrts[1]> -> proper intersection.

  // Note. since no constraint verices can belong to a tetrahedron edge,
  //       the alignement occours only when the edge is included in to one of
  //       the constraint sides.

  // Case: <tE_vrts[0],tE_vrts[1]> is aligned with <c_vrts[0], c_vrts[1]>
  if(vrt_pointInSegment(tE_vrts[0], c_vrts[0], c_vrts[1], mesh) &&
     vrt_pointInSegment(tE_vrts[1], c_vrts[0], c_vrts[1], mesh) &&
     (-1*or_opptE0 !=
      vrt_signe_orient3d(c_vrts[0],c_vrts[1],opptE_vrts[0], opptE_vrts[1], mesh)) )
     return 1;

  // Case: <tE_vrts[0],tE_vrts[1]> is aligned with <c_vrts[1], c_vrts[2]>
  if(vrt_pointInSegment(tE_vrts[0], c_vrts[1], c_vrts[2], mesh) &&
     vrt_pointInSegment(tE_vrts[1], c_vrts[1], c_vrts[2], mesh) &&
     (-1*or_opptE0 !=
      vrt_signe_orient3d(c_vrts[1],c_vrts[2],opptE_vrts[0], opptE_vrts[1], mesh)) )
     return 1;

  // Case: <tE_vrts[0],tE_vrts[1]> is aligned with <c_vrts[2], c_vrts[0]>
  if(vrt_pointInSegment(tE_vrts[0], c_vrts[2], c_vrts[0], mesh) &&
     vrt_pointInSegment(tE_vrts[1], c_vrts[2], c_vrts[0], mesh) &&
     (-1*or_opptE0 !=
      vrt_signe_orient3d(c_vrts[2],c_vrts[0],opptE_vrts[0], opptE_vrts[1], mesh)) )
     return 1;

  // Otherwise...
  // No proper intersection between the whole edge and the constraint)
  return 0;
}

//  Input: the 3 vertices of the contraint-triangle:
//         c_vrts[0], c_vrts[1], c_vrts[2],
//         the vertex of the a tetrahedron: tV,
//         the 3 vertices of the tetrahedron face opposite to tV:
//         opptF_vrts[0], opptF_vrts[1], opptF_vrts[2],
//         pointer to mesh.
// Output: returns 1 if there is a proper inetrsection (with the vertex tetVrt),
//         0 otherwise.
// Note. In the case of a ghost-tet, it is assumed that the ghost vertex
//       is always in opptF_vrts[2].
uint32_t properIntersection_tetVrt_constr(const uint32_t* c_vrts, uint32_t tV,
                                          const uint32_t* opptF_vrts,
                                          const TetMesh* mesh){

  // Necessary Cond. -> tV lies on the constraint plane.
  if(!tetvrtONconstraintplane(tV, c_vrts, mesh))   return 0;

  // Necessary Cond. -> tV belong to the constraint (boundary included).
  if(!vrt_pointInTriangle(tV, c_vrts[0],c_vrts[1],c_vrts[2], mesh))   return 0;

  // ghost-tet: the constraint can not cross the face opposite to ghost-tet,
  //            since it is a face of the convex-hull of the mesh.
  if(opptF_vrts[2] == UINT32_MAX)        return 1;

  uint32_t or_opptF0 =
        vrt_signe_orient3d(c_vrts[0],c_vrts[1],c_vrts[2], opptF_vrts[0], mesh);
  uint32_t or_opptF1 =
        vrt_signe_orient3d(c_vrts[0],c_vrts[1],c_vrts[2], opptF_vrts[1], mesh);
  uint32_t or_opptF2 =
        vrt_signe_orient3d(c_vrts[0],c_vrts[1],c_vrts[2], opptF_vrts[2], mesh);

  // Co-planar tet-faces:
  // Since constraint cannot have vertices inside the co-planar tet-face
  // IF there are no tet-edge constraint-side crossing
  // THEN the intersection is proper,
  // OTHERWISE is improper.
  int a;
  if(or_opptF0==0 && or_opptF1==0){a=0; }
  if(or_opptF1==0 && or_opptF2==0){a=1; }
  if(or_opptF2==0 && or_opptF0==0){a=2; }

  // Co-planar edges
  if(or_opptF0==0){a=3; }
  if(or_opptF1==0){a=4; }
  if(or_opptF2==0){a=5; }

  // Sufficient Cond. -> IF the 3 vertices of oppTetFace belong to the same
  //                     half-space defined by the constraint plane,
  //                     THEN the unique intersection between tet and
  //                     constraint is tV               -> proper intersection.
  if( or_opptF0 == or_opptF1 && or_opptF1 == or_opptF2 )  return 1;

  // Remaining tetrahedron-constraint disposition are such that:
  // - tV belong to the constraint (boundary or interior),
  // - two vertices (tV_UNDER1,tV_UNDER2) of the tetrahedron-face opposite to
  //   tV are in one half-plane defined by the constraint, while the other
  //   vertex (tV_OVER) is in the other half-plane.

  // Necessary Cond. -> tV belong to the constraint boundary.
  // (If tV is in the constraint interior, the intersection is not only tV)
  if(vrt_pointInInnerTriangle(tV,c_vrts[0],c_vrts[1],c_vrts[2], mesh)) return 0;

  // Remaining tetrahedron-constraint disposition are such that:
  // - tV is on the constraint boundary,
  // - two vertices (tV_UNDER1,tV_UNDER2) of the tetrahedron-face opposite to
  //   tV are in one half-plane defined by the constraint, while the other
  //   vertex (tV_OVER) is in the other half-plane.

  // Find the disposition of the vertex of tetrahedron-face opposite to tV
  // w.r.t. the constraint.
  uint32_t tV_UNDER1, tV_UNDER2, tV_OVER;

  if(or_opptF0 == or_opptF1){
      tV_UNDER1 = opptF_vrts[0];
      tV_UNDER2 = opptF_vrts[1];
      tV_OVER = opptF_vrts[2];
    }
    else if(or_opptF1 == or_opptF2){
            tV_UNDER1 = opptF_vrts[1];
            tV_UNDER2 = opptF_vrts[2];
            tV_OVER = opptF_vrts[0];
          }
          else{
            tV_UNDER1 = opptF_vrts[2];
            tV_UNDER2 = opptF_vrts[0];
            tV_OVER = opptF_vrts[1];
          }

  // Sufficient Cond. -> IF the following conditions hold:
  //   (i) tet-edge <tV_OVER, tV_UNDER1>
  //       do NOT cross the constraint (boundary or interior)
  //  (ii) tet-edge <tV_OVER, tV_UNDER2>
  //       do NOT cross the constraint (boundary or interior)
  // (iii) tet-face <tV, tV_OVER, tV_UNDER1>
  //       is NOT inner crossed by any constraint side,
  //  (iv) tet-face <tV, tV_OVER, tV_UNDER2>
  //       is NOT inner crossed by any constraint side,
  //   (v) tet-face <tV_OVER, tV_UNDER1, tV_UNDER2>
  //       is NOT inner crossed by any constraint side,
  //                      THEN the unique intersection between tet
  //                      and constraint is tV           -> proper intersection.

  // Condition (i)
  if(vrt_innerSegmentCrossesTriangle(tV_OVER, tV_UNDER1,
                                     c_vrts[0], c_vrts[1], c_vrts[2], mesh))
    return 0;


  // Condition (ii)
  if(vrt_innerSegmentCrossesTriangle(tV_OVER, tV_UNDER2,
                                     c_vrts[0], c_vrts[1], c_vrts[2], mesh))
    return 0;

  // Condition (iii)
  if(vrt_innerTrianglesCrosses(c_vrts[0], c_vrts[1], c_vrts[2],
                               tV, tV_OVER, tV_UNDER1, mesh)    )
    return 0;

  // Condition (iv)
  if(vrt_innerTrianglesCrosses(c_vrts[0],c_vrts[1],c_vrts[2],
                               tV,tV_OVER,tV_UNDER2, mesh)    )
    return 0;

  // Condition (v)
  if(vrt_innerTrianglesCrosses(c_vrts[0], c_vrts[1], c_vrts[2],
                               tV_OVER, tV_UNDER1, tV_UNDER2, mesh)    )
    return 0;

  // Otherwise... (The sufficient condition is verified -> proper intersection)
  return 1;
}

//-----------------------------------------
// FUNCTIONS TO FIND IMPROPER INTERSECTIONS
//-----------------------------------------

//  Input: 3 vertices of the contraint-triangle: constr_vrts[3],
//         array of indices of the tetrahedra that intersects
//         the constraint: tet_list,
//         number of tetrahedra that intersects the constraint: num_tet_list,
//         array of marker for the intersections
//         tetrahedra-constraint: mark_TetIntersection,
//         pointer to mesh.
// Output: by using mark_TetIntersection marks the improper
//         intersections tetrahedra-constraint.
// Note. It is assumed that each tetrahedron of tet_list intersects
//       (properly or improperly) the constraint.
void improperIntersection(const uint32_t* constr_vrts, const uint64_t* tet_list,
                          uint64_t num_tet_list, uint32_t* mark_TetIntersection,
                          const TetMesh* mesh){
    for(uint64_t list_ind=0; list_ind<num_tet_list; list_ind++){
        uint64_t tet_ind = tet_list[list_ind]; // tetrahedron tet

        // tet vertices
        uint32_t tetVrts[4];
        extract_tetVrts(tetVrts, tet_ind, mesh);

        // Usefull strucutures for subsimplex of tet.
        uint32_t tetFace_vrts[3];
        uint32_t tetEdge_vrts[2];
        uint32_t oppTetEdge_vrts[2];
        uint32_t tetVrt;
        uint32_t oppTetFace_vrts[3];

        // Property: proper intersections and improper intersections
        //           form a partition of intersections.

        // ghost-tet case
        if(tetVrts[3] == UINT32_MAX){

          // Proper intersection of dimension 2:
          extract_tetFaceVrts(tetFace_vrts, tet_ind, 3, mesh);
          if(properIntersection_tetFace_constr(constr_vrts, tetFace_vrts, mesh))
            continue;

          // Proper intersection of dimension 1:
          extract_tetEdgeVrts(tetEdge_vrts, tet_ind, 0, 1, mesh);
          extract_tetEdgeVrts(oppTetEdge_vrts, tet_ind, 2, 3, mesh);
          if(properIntersection_tetEdge_constr(constr_vrts, tetEdge_vrts,
                                               oppTetEdge_vrts, mesh)     )
            continue;

          extract_tetEdgeVrts(tetEdge_vrts, tet_ind, 1, 2, mesh);
          extract_tetEdgeVrts(oppTetEdge_vrts, tet_ind, 0, 3, mesh);
          if(properIntersection_tetEdge_constr(constr_vrts, tetEdge_vrts,
                                                oppTetEdge_vrts, mesh)    )
            continue;

          extract_tetEdgeVrts(tetEdge_vrts, tet_ind, 2, 0, mesh);
          extract_tetEdgeVrts(oppTetEdge_vrts, tet_ind, 1, 3, mesh);
          if(properIntersection_tetEdge_constr(constr_vrts, tetEdge_vrts,
                                                oppTetEdge_vrts, mesh)    )
            continue;

          // Proper intersection of dimension 0:
          tetVrt = tetVrts[0];
          extract_tetFaceVrts(oppTetFace_vrts, tet_ind, 0, mesh);
          if(properIntersection_tetVrt_constr(constr_vrts, tetVrt,
                                              oppTetFace_vrts, mesh) )
            continue;

          tetVrt = tetVrts[1];
          extract_tetFaceVrts(oppTetFace_vrts, tet_ind, 1, mesh);
          if(properIntersection_tetVrt_constr(constr_vrts, tetVrt,
                                              oppTetFace_vrts, mesh) )
            continue;

          tetVrt = tetVrts[2];
          extract_tetFaceVrts(oppTetFace_vrts, tet_ind, 2, mesh);
          if(properIntersection_tetVrt_constr(constr_vrts, tetVrt,
                                              oppTetFace_vrts, mesh) )
            continue;

          // If no proper intersection -> the intersection is improper.
          mark_TetIntersection[tet_ind] = IMPROPER_INTERSECTION;
          continue;
        }

        // normal tet case (i.e. non-ghost)

        // Proper intersection of dimension 2:
        extract_tetFaceVrts(tetFace_vrts, tet_ind, 3, mesh);
        if(properIntersection_tetFace_constr(constr_vrts, tetFace_vrts, mesh) )
          continue;

        extract_tetFaceVrts(tetFace_vrts, tet_ind, 0, mesh);
        if(properIntersection_tetFace_constr(constr_vrts, tetFace_vrts, mesh) )
          continue;

        extract_tetFaceVrts(tetFace_vrts, tet_ind, 1, mesh);
        if(properIntersection_tetFace_constr(constr_vrts, tetFace_vrts, mesh) )
          continue;

        extract_tetFaceVrts(tetFace_vrts, tet_ind, 2, mesh);
        if(properIntersection_tetFace_constr(constr_vrts, tetFace_vrts, mesh) )
          continue;

        // Proper intersection of dimension 1:
        extract_tetEdgeVrts(tetEdge_vrts, tet_ind, 0, 1, mesh);
        extract_tetEdgeVrts(oppTetEdge_vrts, tet_ind, 2, 3, mesh);
        if(properIntersection_tetEdge_constr(constr_vrts, tetEdge_vrts,
                                             oppTetEdge_vrts, mesh)     )
          continue;

        extract_tetEdgeVrts(tetEdge_vrts, tet_ind, 1, 2, mesh);
        extract_tetEdgeVrts(oppTetEdge_vrts, tet_ind, 3, 0, mesh);
        if(properIntersection_tetEdge_constr(constr_vrts, tetEdge_vrts,
                                             oppTetEdge_vrts, mesh)     )
          continue;

        extract_tetEdgeVrts(tetEdge_vrts, tet_ind, 2, 3, mesh);
        extract_tetEdgeVrts(oppTetEdge_vrts, tet_ind, 0, 1, mesh);
        if(properIntersection_tetEdge_constr(constr_vrts, tetEdge_vrts,
                                             oppTetEdge_vrts, mesh)     )
          continue;

        extract_tetEdgeVrts(tetEdge_vrts, tet_ind, 3, 0, mesh);
        extract_tetEdgeVrts(oppTetEdge_vrts, tet_ind, 1, 2, mesh);
        if(properIntersection_tetEdge_constr(constr_vrts, tetEdge_vrts,
                                             oppTetEdge_vrts, mesh)     )
          continue;

        extract_tetEdgeVrts(tetEdge_vrts, tet_ind, 1, 3, mesh);
        extract_tetEdgeVrts(oppTetEdge_vrts, tet_ind, 0, 2, mesh);
        if(properIntersection_tetEdge_constr(constr_vrts, tetEdge_vrts,
                                              oppTetEdge_vrts, mesh)     )
          continue;

        extract_tetEdgeVrts(tetEdge_vrts, tet_ind, 2, 0, mesh);
        extract_tetEdgeVrts(oppTetEdge_vrts, tet_ind, 3, 1, mesh);
        if(properIntersection_tetEdge_constr(constr_vrts, tetEdge_vrts,
                                             oppTetEdge_vrts, mesh)       )
          continue;

        // If there is a common edge between tet and constraint:
        // IMPROPER INTERSECTION
        if(segment_is_triSide(tetVrts[0], tetVrts[1], constr_vrts) ||
           segment_is_triSide(tetVrts[1], tetVrts[2], constr_vrts) ||
           segment_is_triSide(tetVrts[2], tetVrts[3], constr_vrts) ||
           segment_is_triSide(tetVrts[3], tetVrts[0], constr_vrts) ||
           segment_is_triSide(tetVrts[1], tetVrts[3], constr_vrts) ||
           segment_is_triSide(tetVrts[0], tetVrts[2], constr_vrts)   ){
          mark_TetIntersection[tet_ind] = IMPROPER_INTERSECTION;
          continue;
        }

        // Proper intersection of dimension 0:
        tetVrt = tetVrts[0];
        extract_tetFaceVrts(oppTetFace_vrts, tet_ind, 0, mesh);
        if(properIntersection_tetVrt_constr(constr_vrts, tetVrt,
                                            oppTetFace_vrts, mesh) )
          continue;

        tetVrt = tetVrts[1];
        extract_tetFaceVrts(oppTetFace_vrts, tet_ind, 1, mesh);
        if(properIntersection_tetVrt_constr(constr_vrts, tetVrt,
                                            oppTetFace_vrts, mesh) )
          continue;

        tetVrt = tetVrts[2];
        extract_tetFaceVrts(oppTetFace_vrts, tet_ind, 2, mesh);
        if(properIntersection_tetVrt_constr(constr_vrts, tetVrt,
                                            oppTetFace_vrts, mesh) )
          continue;

        tetVrt = tetVrts[3];
        extract_tetFaceVrts(oppTetFace_vrts, tet_ind, 3, mesh);
        if(properIntersection_tetVrt_constr(constr_vrts, tetVrt,
                                            oppTetFace_vrts, mesh) )
          continue;

        // If no proper intersection -> the intersection is improper.
        mark_TetIntersection[tet_ind] = IMPROPER_INTERSECTION;

    }
}


//------------------------------------------------------------------------
// FUNCTIONS TO FIND (GENERIC) INTERSECTIONS ALONG THE CONSTRAINT BOUNDARY
//------------------------------------------------------------------------

// elem0 can be only 1,2 or 3.
// Its value (i) means how many elements (connecting_vrts[i]) have to be assigned.
static inline void fill_connecting_vrts(uint32_t* connecting_vrts,
                                        uint32_t elem0, uint32_t elem1,
                                        uint32_t elem2, uint32_t elem3 ){
    connecting_vrts[0] = elem0;
    connecting_vrts[1] = elem1;
    if(elem0==1) return;
    connecting_vrts[2] = elem2;
    if(elem0==2) return;
    connecting_vrts[3] = elem3;
}

//  Input: pointer to the mesh,
//         index of a vertex that belong to the constraint in which we have to
//         find intersections: v_curr,
//         vertex index of the constraint side endponit we are
//         moving towards: v_stop,
//         index of the constraint vertex that do not belong to the side of the
//         constraint we are travelling: other_constr_vrt,
//         array of tetrahedra marker: mark_TetIntersection,
//         pointer to array of vertex index type: connecting_vrts[4],
//         pointer to tetrahedron index type: nextTet_ind,
//         pointer to a tetrahedron index type: num_intersecated_tet.
// Output: returns the indices of tetrahedra incident in v_curr and not
//         already visited,
//         by using num_intersecated_tet returns the number of tetrahedra
//         in the above array,
//         by using mark_TetIntersection marks the intersection as explained
//         at the beginning of the function insert_constraints,
//         by uising connecting_vrts returns the informations to continue the
//         travelling along the constraint side,
//         by using nextTet_ind returns the tetrahedron from which
//         continue the side travelling.
uint64_t* intersections_TetVrtOnConstraintSide(TetMesh* mesh,
                                               uint32_t v_curr, uint32_t v_stop,
                                               uint32_t other_constr_vrt,
                                               uint32_t* mark_TetIntersection,
                                               uint32_t* connecting_vrts,
                                               uint64_t* nextTet_ind,
                                               uint64_t* num_intersecated_tet){

    //   Consider the intersection between the constraint side (v_start,v_stop)
    //   and a tetrahedron incident in v_curr.
    //   There are 5 cases:
    //   (0) (v_curr,v_stop) intesect the tetrahedron only in v_curr.
    //   (1) (v_curr,v_stop) intersects exactly 1 tetrahedron in the
    //       interior of the face opposite to v_curr.
    //   (2) (v_curr,v_stop) intersects exactly 2 tetrahedra in the
    //       interior of theire common face.
    //   (3) (v_curr,v_stop) intersects n tetrahedra along a common edge.
    //   (3') like 3), but the other endpoint of the common edge is v_stop.

    uint64_t num_incTet;
    uint64_t* incTet = mesh->incident_tetrahedra(v_curr, &num_incTet);
    // ghost-tets are NOT returned.

    // Not already visited tet_
    uint64_t num_newTetIn_incTet = 0;
    uint64_t* newTetIn_incTet = (uint64_t*) malloc(sizeof(uint64_t)*num_incTet);
    for(uint64_t i=0; i<num_incTet; i++)
        INSERT_TET_IN_LIST(incTet[i], newTetIn_incTet,
                           num_newTetIn_incTet, mark_TetIntersection);
    // Mark all new tetrahedra: they have a generic intersection in v_curr.
    for(uint64_t i=0; i<num_newTetIn_incTet; i++)
        mark_TetIntersection[ newTetIn_incTet[i] ] = INTERSECTION;

    // Cycle over the tetrahedra incident in v_curr to fill connecting_vrts
    // and return information on how the constraint side exit from
    // the intersecated tet_

    uint64_t tet_ind;
    uint32_t v_oppf_inds[3]; // Indices of the vertices of
                             // the face opposite to v_curr.

    for(uint64_t i=0; i<num_incTet; i++){
        tet_ind = incTet[i];
        // Fill v_oppf_inds[0,1,2].
        opposite_face_vertices(mesh, tet_ind, v_curr, v_oppf_inds);

        // Check if the case (3'):
        // one of the vertices of the face opposite to v_curr is v_stop.
        if( arrayUINT32_contains_elem(v_oppf_inds, 3, v_stop) ){

          fill_connecting_vrts(connecting_vrts,1,v_stop,UINT32_MAX,UINT32_MAX);
          break;
        }

        // Check if the case (1):
        // (v_curr, v_stop) intersect the interior of the face
        // opposite to v_curr.
        if(vrt_innerSegmentCrossesInnerTriangle(v_curr, v_stop,
                                              v_oppf_inds[0], v_oppf_inds[1],
                                              v_oppf_inds[2],mesh)            ){

          fill_connecting_vrts(connecting_vrts,3,
                               v_oppf_inds[0],v_oppf_inds[1],v_oppf_inds[2]);
          break;
        }


        // Check if case (2):
        // (v_curr, v_stop) intersect the interior of an edge of the face
        // opposite to v_start.

        // Edge < v_oppf_inds[0] , v_oppf_inds[1] >
        if(vrt_innerSegmentsCross(v_curr,v_stop,
                                  v_oppf_inds[0],v_oppf_inds[1], mesh) ){

          fill_connecting_vrts(connecting_vrts,2,
                               v_oppf_inds[0],v_oppf_inds[1],UINT32_MAX);
          break;
        }

        // Edge < v_oppf_inds[1] , v_oppf_inds[2] >
        if(vrt_innerSegmentsCross(v_curr,v_stop, v_oppf_inds[1],v_oppf_inds[2],
                                  mesh)                                       ){

          fill_connecting_vrts(connecting_vrts,2,
                               v_oppf_inds[1],v_oppf_inds[2],UINT32_MAX);
          break;
        }

        // Edge < v_oppf_inds[2] , v_oppf_inds[0] >
        if(vrt_innerSegmentsCross(v_curr,v_stop, v_oppf_inds[2],v_oppf_inds[0],
                                  mesh)                                       ){

          fill_connecting_vrts(connecting_vrts,2,
                               v_oppf_inds[2],v_oppf_inds[0],UINT32_MAX);
          break;
        }


        // Check if case (3):
        // (v_curr, v_stop) intersect a vertex of the face opposite to v_curr.

        // Vertex v_oppf_inds[0]
        if( vrt_pointInInnerSegment(v_oppf_inds[0], v_curr, v_stop, mesh) ){

          fill_connecting_vrts(connecting_vrts,1,
                               v_oppf_inds[0],UINT32_MAX,UINT32_MAX);
          break;
        }

        // Vertex v_oppf_inds[1]
        if( vrt_pointInInnerSegment(v_oppf_inds[1], v_curr, v_stop, mesh) ){

          fill_connecting_vrts(connecting_vrts,1,
                               v_oppf_inds[1],UINT32_MAX,UINT32_MAX);
          break;
        }

        // Vertex v_oppf_inds[2]
        if(vrt_pointInInnerSegment(v_oppf_inds[2], v_curr, v_stop, mesh) ){

          fill_connecting_vrts(connecting_vrts,1,
                               v_oppf_inds[2],UINT32_MAX,UINT32_MAX);
          break;
        }

        // Case (0): intersects only v_curr -> visit the incTet[i+1]

    } // END Cycle over the tetrahedra incident in v_curr

    // Returns the number of new tetrahedra intersecated.
    *num_intersecated_tet = num_newTetIn_incTet;

    // Returns the tetrahedron from which the side travelling continue.
    // In case (3') return the tetrahedron itself.
    *nextTet_ind = tet_ind;

    free(incTet);
    incTet = NULL;
    return newTetIn_incTet;
}

//  Input: pointer to the mesh,
//         endpoints of the tetrahedron edge whose interior intersects
//         the constraint: tet_cutted_edge,
//         vertex index of the constraint side endponit we are
//         moving away: v_start,
//         vertex index of the constraint side endponit we are
//         moving towards: v_stop,
//         index of the constraint vertex that do not belong to the side of the
//         constraint we are travelling: other_constr_vrt,
//         array of tetrahedra marker: mark_TetIntersection,
//         pointer to a tetrahedron index type: num_intersecated_tet,
//         pointer to array of vertex index type: connecting_vrts[4],
//         pointer to tetrahedron index type: nextTet_ind, (contains the tet
//         from which start searching).
// Output: returns the indices of tetrahedra incident in the edge
//         tet_cutted_edge and not already visited,
//         by using num_intersecated_tet returns the number of tetrahedra in
//         the above array,
//         by using mark_TetIntersection marks the intersection as explained
//         at the beginning of the function insert_constraints,
//         by uising connecting_vrts returns the informations to continue
//         the travelling along the constraint side,
//         by using nextTet_ind returns the tetrahedron from which
//         continue the side travelling.
uint64_t* intersections_TetEdgeCrossConstraintSide(TetMesh* mesh,
                                                const uint32_t* tet_cutted_edge,
                                                  uint32_t v_start,
                                                  uint32_t v_stop,
                                                  uint32_t other_constr_vrt,
                                                uint32_t* mark_TetIntersection,
                                                  uint32_t* connecting_vrts,
                                                  uint64_t* nextTet_ind,
                                                uint64_t* num_intersecated_tet){

    //   Name p* the point in which (the interior of) tet_cutted_edge
    //   intersects the constraint side (v_start,v_stop).
    //   Consider the intersection between the constraint side (v_start,v_stop)
    //   and a tetrahedron incident in tet_cutted_edge.
    //   There are 6 cases:
    //   (0) (v_start,v_stop) intesect the tetrahedron only in p*.
    //   (1) (v_start,v_stop) intersects exactly 1 tetrahedron in the
    //       interior of one of the two faces opposite to an endpoints
    //       of tet_cutted_edge.
    //   (2) (v_start,v_stop) intersects exactly 1 tetrahedron in the interior
    //       of the edge whose endpoints are not those of tet_cutted_edge.
    //   (3) (v_start,v_stop) intersects exactly 2 tetrahedra along a common
    //       face and pass through the vertex, of that face,
    //       opposite to tet_cutted_edge.
    //  (3') like (3), but pass through v_stop.
    //   (4) (v_start,v_stop) intersects exactly 2 tetrahedra along a common
    //       face and cuts one the edge (different from tet_cutted_edge)
    //       of that face.

    uint64_t num_incTet;
    uint64_t* incTet = mesh->ETrelation(tet_cutted_edge,
                                  *nextTet_ind, &num_incTet);

    // Not already visited tet_
    uint64_t num_newTetIn_incTet = 0;
    uint64_t* newTetIn_incTet = (uint64_t*) malloc(sizeof(uint64_t)*num_incTet);
    for(uint64_t i=0; i<num_incTet; i++)
        INSERT_TET_IN_LIST(incTet[i], newTetIn_incTet, num_newTetIn_incTet,
                           mark_TetIntersection);
    // Mark all new tetrahedra:
    // they have an intersection in p* (improper or not).
    for(uint64_t i=0; i<num_newTetIn_incTet; i++)
        mark_TetIntersection[ newTetIn_incTet[i] ] = INTERSECTION;

    // Cycle over the tetrahedra incident in tet_cutted_edge to fill
    // connecting_vrts and return information on how the constraint
    // side exit from the intersecated tet_

    #ifdef DEBUG
    uint32_t all_tets_intersects_only_pstar = 1;
    #endif

    uint64_t tet_ind;
    uint32_t v_oppe_inds[2];    // Indices of the endpoints of the edge
                                // opposite to tet_cutted_edge.
    for(uint64_t ind=0; ind<num_incTet; ind++){
        tet_ind = incTet[ind];

        // Skip ghost-tet
        if(mesh->tet_node[4*tet_ind+3]==UINT32_MAX) continue;

        // Fill v_oppe_inds[0,1].
        opposite_side_vertices(mesh, tet_ind, tet_cutted_edge, v_oppe_inds);

        // Check if the case (3'):
        // one of the vertices of the side opposite to tet_cutted_edge is v_stop.
        if( arrayUINT32_contains_elem(v_oppe_inds, 2, v_stop) ){

          fill_connecting_vrts(connecting_vrts,1,v_stop,UINT32_MAX,UINT32_MAX);
          break;
        }


        // Check if the case (1):
        // (v_start, v_stop) intersect the interior of a face opposite to an
        // endpoint of tet_cutted_edge.

        // Face <v_oppe_inds[0], v_oppe_inds[1], tet_cutted_edge[0]>
        // Check also if it is the right diretction: against v_start.
        if(vrt_innerSegmentCrossesInnerTriangle(v_start,v_stop,
              v_oppe_inds[0],v_oppe_inds[1],tet_cutted_edge[0], mesh)    &&
            vrtsInSameHalfSpace(tet_cutted_edge[1],v_start,
              v_oppe_inds[0],v_oppe_inds[1],tet_cutted_edge[0], mesh)        ){

          fill_connecting_vrts(connecting_vrts,3,
                              v_oppe_inds[0],v_oppe_inds[1],tet_cutted_edge[0]);
          break;
        }

        // Face <v_oppe_inds[0], v_oppe_inds[1], tet_cutted_edge[1]>
        // Check also if it is the right diretction: against v_start.
        if(vrt_innerSegmentCrossesInnerTriangle(v_start,v_stop,
              v_oppe_inds[0],v_oppe_inds[1],tet_cutted_edge[1], mesh)    &&
           vrtsInSameHalfSpace(tet_cutted_edge[0],v_start,
              v_oppe_inds[0],v_oppe_inds[1],tet_cutted_edge[1], mesh)        ){

          fill_connecting_vrts(connecting_vrts,3,
                              v_oppe_inds[0],v_oppe_inds[1],tet_cutted_edge[1]);
          break;
        }

        // Check if case (2):
        // (v_start, v_stop) intersect the interior of the edge opposite
        // to tet_cutted_edge.
        // Check also if it is the right diretction: against v_start.
        if(vrt_innerSegmentsCross(v_start,v_stop,
                                  v_oppe_inds[0],v_oppe_inds[1], mesh)   &&
          vrtsInSameHalfSpace(tet_cutted_edge[0],v_start,
              v_oppe_inds[0],v_oppe_inds[1],tet_cutted_edge[1], mesh)        ){

          fill_connecting_vrts(connecting_vrts,2,
                               v_oppe_inds[0],v_oppe_inds[1],UINT32_MAX);
          break;
        }


        // Check if case (3):
        // (v_start,v_stop) intersects the common face between two tetrahedra
        // and pass through one of the endpoints of v_oppe.

        // Endpoint v_oppe_inds[0].
        // Check also if it is the right diretction: against v_start.
        if(vrt_pointInInnerSegment(v_oppe_inds[0], v_start, v_stop, mesh) &&
           vrt_same_half_plane(v_oppe_inds[0],v_stop,
                               tet_cutted_edge[0],tet_cutted_edge[1], mesh)  ){

          fill_connecting_vrts(connecting_vrts,1,
                               v_oppe_inds[0],UINT32_MAX,UINT32_MAX);
          break;
        }

        // Endpoint v_oppe_inds[1].
        // Check also if it is the right diretction: against v_start.
        if(vrt_pointInInnerSegment(v_oppe_inds[1], v_start, v_stop, mesh)  &&
           vrt_same_half_plane(v_oppe_inds[1],v_stop,
                               tet_cutted_edge[0],tet_cutted_edge[1], mesh)  ){

          fill_connecting_vrts(connecting_vrts,1,
                               v_oppe_inds[1],UINT32_MAX,UINT32_MAX);
          break;
        }


        // Check if case (4):
        // (v_start,v_stop) intersects the common face between two tetrahedra
        // and cross the interior of one of the edges connecting
        // tet_cutted_edge with the opposite vertex.

        // Edge <tet_cutted_edge[0],v_oppe_inds[0]>
        // Check also if it is the right diretction: against v_start.
        if(vrt_innerSegmentsCross(v_start,v_stop,
                                tet_cutted_edge[0],v_oppe_inds[0], mesh)    &&
            vrt_same_half_plane(v_oppe_inds[0],v_stop,
                                tet_cutted_edge[0],tet_cutted_edge[1], mesh) ){

          fill_connecting_vrts(connecting_vrts,2,
                               tet_cutted_edge[0],v_oppe_inds[0],UINT32_MAX);
          break;
        }

        // Edge <tet_cutted_edge[0],v_oppe_inds[1]>
        // Check also if it is the right diretction: against v_start.
        if(vrt_innerSegmentsCross(v_start,v_stop,
                                  tet_cutted_edge[0],v_oppe_inds[1], mesh) &&
           vrt_same_half_plane(v_oppe_inds[1],v_stop,
                               tet_cutted_edge[0],tet_cutted_edge[1], mesh)  ){

          fill_connecting_vrts(connecting_vrts,2,
                               tet_cutted_edge[0],v_oppe_inds[1],UINT32_MAX);
          break;
        }

        // Edge <tet_cutted_edge[1],v_oppe_inds[0]>
        // Check also if it is the right diretction: against v_start.
        if(vrt_innerSegmentsCross(v_start,v_stop,
                                   tet_cutted_edge[1],v_oppe_inds[0],  mesh) &&
           vrt_same_half_plane(v_oppe_inds[0],v_stop,
                               tet_cutted_edge[0],tet_cutted_edge[1], mesh)  ){

          fill_connecting_vrts(connecting_vrts,2,
                               tet_cutted_edge[1],v_oppe_inds[0],UINT32_MAX);
          break;
        }

        // Edge <tet_cutted_edge[1],v_oppe_inds[1]>
        // Check also if it is the right diretction: against v_start.
        if(vrt_innerSegmentsCross(v_start,v_stop,
                                  tet_cutted_edge[1],v_oppe_inds[1],  mesh) &&
           vrt_same_half_plane(v_oppe_inds[1],v_stop,
                               tet_cutted_edge[0],tet_cutted_edge[1], mesh)   ){

          fill_connecting_vrts(connecting_vrts,2,
                               tet_cutted_edge[1],v_oppe_inds[1],UINT32_MAX);
          break;
        }

        // Case (0): intersects only tet_cutted_edge -> visit the inc_Tet[i+1].

    } // END Cycle over the tetrahedra incident in tet_cutted_edge

    // Returns the number of new tetrahedra intersecated.
    *num_intersecated_tet = num_newTetIn_incTet;

    // Returns the tetrahedron from which the side travelling continue.
    // In case (3') return the tetrahedron itself.
    *nextTet_ind = tet_ind;

    free(incTet);
    incTet = NULL;
    return newTetIn_incTet;
}

//  Input: pointer to the mesh,
//         vertices of the tetrahedron face whose interior intersects
//         the constraint: tet_crossed_face,
//         vertex index of the constraint side endponit we are
//         moving away: v_start,
//         vertex index of the constraint side endponit we are
//         moving towards: v_stop,
//         index of the constraint vertex that do not belong to the side of
//         the constraint we are travelling: other_constr_vrt,
//         array of tetrahedra marker: mark_TetIntersection,
//         pointer to array of vertex index type: connecting_vrts[4].
//         pointer to tetrahedron index type: nextTet_ind,
//         (contains the tet from which start searching).
// Output: returns 1 if the tetrahedron incident in tet_crossed_face and
//         oriented throuh v_stop it has not been visited yet, 0 otherwise.
//         by uising connecting_vrts returns the informations to continue
//         the travelling along the constraint side,
//         by using nextTet_ind returns the tetrahedron from which
//         continue the side travelling.
//         by using mark_TetIntersection marks the intersection as
//         explained at the beginning of the function insert_constraints,
uint32_t intersections_TetFacePiercedConstraintSide(TetMesh* mesh,
                                               const uint32_t* tet_crossed_face,
                                                    uint32_t v_start,
                                                    uint32_t v_stop,
                                                    uint32_t other_constr_vrt,
                                                uint32_t* mark_TetIntersection,
                                                    uint32_t* connecting_vrts,
                                                    uint64_t* nextTet_ind   ){

    // There is only one tetrahedron to consider:
    // the one (nextTet) of index nextTet_ind.
    // This tetrahedron has been alraedy marked during the previous step.
    // There are 4 cases:
    //  (1) (v_start,v_stop) intersects the interior of a fece of nextTet
    //      different from tet_crossed_face.
    //  (2) (v_start,v_stop) intersects the interior of an edge
    //      of nextTet that do not belong to tet_crossed_face.
    //  (3) (v_start,v_stop) pass through the vertex of nextTet opposite
    //      to tet_crossed_face.
    // (3') like (3), but pass through v_stop.

    uint64_t tet_ind = *nextTet_ind;
    uint32_t v_opp_ind = opposite_vertex_face(mesh, tet_ind, tet_crossed_face);
    uint64_t opp_tet_ind = adjTet_oppsiteTo_vertex(mesh, tet_ind, v_opp_ind);
    v_opp_ind = opposite_vertex_face(mesh, opp_tet_ind, tet_crossed_face);

    // Check if opp_tet has been already visited.
    uint32_t isNew=0;
    if(mark_TetIntersection[ opp_tet_ind ] == 0){
        isNew=1;
        // Mark the new tetrahedron:
        // it has an improper intersection with the constraint.
        mark_TetIntersection[ opp_tet_ind ] = INTERSECTION;
    }

    // Returns the tetrahedron from which continue travelling the side,
    // it has to be added to intersecatedTet if isNew==1.
    *nextTet_ind = opp_tet_ind;

    // Check if case (3'): v_opp is v_stop.
    if(v_opp_ind == v_stop){

      fill_connecting_vrts(connecting_vrts, 1, v_stop, UINT32_MAX, UINT32_MAX);
      return isNew;
    }


    // Check if case (1):
    // (v_start, v_stop) intersect the interior of a face different
    // from tet_crossed_face.

    // Face <tet_crossed_face[0], tet_crossed_face[1], v_opp_ind>
    if(vrt_innerSegmentCrossesInnerTriangle(v_start,v_stop,
                              tet_crossed_face[0],tet_crossed_face[1],v_opp_ind,
                                            mesh)                             ){

      fill_connecting_vrts(connecting_vrts,3,
                           tet_crossed_face[0], tet_crossed_face[1], v_opp_ind);
      return isNew;
    }

    // Face <tet_crossed_face[1], tet_crossed_face[2], v_opp_ind>
    if(vrt_innerSegmentCrossesInnerTriangle(v_start,v_stop,
                              tet_crossed_face[1],tet_crossed_face[2],v_opp_ind,
                                            mesh)                             ){

      fill_connecting_vrts(connecting_vrts,3,
                           tet_crossed_face[1], tet_crossed_face[2], v_opp_ind);
      return isNew;
    }

    // Face <tet_crossed_face[2], tet_crossed_face[0], v_opp_ind>
    if( vrt_innerSegmentCrossesInnerTriangle(v_start,v_stop,
                              tet_crossed_face[2],tet_crossed_face[0],v_opp_ind,
                                             mesh)                            ){

      fill_connecting_vrts(connecting_vrts,3,
                           tet_crossed_face[2], tet_crossed_face[0], v_opp_ind);
      return isNew;
    }


    // Check if case (2):
    // (v_start, v_stop) intersect the interior of the edge
    // that have v_opp_ind as endpoint.

    // Edge <tet_crossed_face[0], v_opp_ind>
    if( vrt_innerSegmentsCross(v_start,v_stop,
                               tet_crossed_face[0],v_opp_ind, mesh) ){

      fill_connecting_vrts(connecting_vrts, 2,
                           tet_crossed_face[0], v_opp_ind, UINT32_MAX);
      return isNew;
    }

    // Edge <tet_crossed_face[1], v_opp_ind>
    if( vrt_innerSegmentsCross(v_start,v_stop,
                               tet_crossed_face[1],v_opp_ind, mesh) ){

      fill_connecting_vrts(connecting_vrts, 2,
                           tet_crossed_face[1], v_opp_ind, UINT32_MAX);
      return isNew;
    }

    // Edge <tet_crossed_face[2], v_opp_ind>
    if( vrt_innerSegmentsCross(v_start,v_stop,
                               tet_crossed_face[2], v_opp_ind, mesh) ){

      fill_connecting_vrts(connecting_vrts, 2,
                           tet_crossed_face[2], v_opp_ind, UINT32_MAX);
      return isNew;
    }

    // Otherwise it must be case (3):
    // (v_start,v_stop) pass through the vertex of index v_opp_ind.

    fill_connecting_vrts(connecting_vrts, 1, v_opp_ind, UINT32_MAX, UINT32_MAX);
    return isNew;

}


//-----------------------------------------------------------------------------
// FUNCTIONS TO FIND PROPER & IMPROPER INTERSECTIONS IN THE CONSTRAINT INTERIOR
//-----------------------------------------------------------------------------

//  Input: the index of the tetrahedron (tet): tet_ind,
//         the vertices index of the triangle: c,
//         pointer to mesh.
// Output: 0 -> no intersection;
//         1 -> triangle intersects tetrahedron interior;(improper intersection)
//         2 -> triangle intersects tetrahedron boundary;  (proper intersection)
//         3 -> a tet-face is complely contained in the constraint. (proper but save)
// Note. It is assumed that tet does not intersects the triangle boundary.
uint32_t tet_intersects_triInterior(uint64_t tet_ind, const uint32_t* c,
                                                  const TetMesh* mesh,
                                                  uint64_t* contain_face_ID){

    uint32_t t[4], t_in_c[4], or_t_WRT_c[4];
    uint32_t num_t_in_c=0;
    extract_tetVrts(t, tet_ind, mesh);

    or_t_WRT_c[0] = vrt_signe_orient3d(t[0], c[0], c[1], c[2], mesh);
    or_t_WRT_c[1] = vrt_signe_orient3d(t[1], c[0], c[1], c[2], mesh);
    or_t_WRT_c[2] = vrt_signe_orient3d(t[2], c[0], c[1], c[2], mesh);

    t_in_c[0] = (or_t_WRT_c[0]==0 && vrt_pointInInnerTriangle(t[0], c[0],c[1],c[2], mesh));
    t_in_c[1] = (or_t_WRT_c[1]==0 && vrt_pointInInnerTriangle(t[1], c[0],c[1],c[2], mesh));
    t_in_c[2] = (or_t_WRT_c[2]==0 && vrt_pointInInnerTriangle(t[2], c[0],c[1],c[2], mesh));

    num_t_in_c = t_in_c[0] + t_in_c[1] + t_in_c[2];

    // Ghost-tetrahedron case (ghost vertex may be only in last position).
    if(t[3] == UINT32_MAX){
      // Face opposite to ghost vertex is contained in the trinagle.
      if(num_t_in_c == 3)  return 2; // It is useless for the BSPsubdivision.
      // Any other intersection with constraint-interior are not possible.
      return 0;
    }

    or_t_WRT_c[3] = vrt_signe_orient3d(t[3], c[0], c[1], c[2], mesh);
    t_in_c[3] = (or_t_WRT_c[3]==0 && vrt_pointInInnerTriangle(t[3], c[0],c[1],c[2], mesh));
    num_t_in_c += t_in_c[3];

    // 3 tetrahedron vertices belong to the triangle interior ->
    // a tetrahedron face is completely contained in to the constraint:
    // it is a proper intersection but it has to be saved to correctly divide
    // BSPfaces.
    if(num_t_in_c == 3){

      if(t_in_c[0]==0) *contain_face_ID = 0;
      else if(t_in_c[1]==0) *contain_face_ID = 1;
      else if(t_in_c[2]==0) *contain_face_ID = 2;
      else if(t_in_c[3]==0) *contain_face_ID = 3;

      return 3;
    }

    // 2 tetrahedron vertices belong to the triangle interior.
    if(num_t_in_c == 2){

      uint32_t or_out_t[2];
      uint32_t j=0;
      if(t_in_c[0] == 0)  or_out_t[j++] = or_t_WRT_c[0];
      if(t_in_c[1] == 0)  or_out_t[j++] = or_t_WRT_c[1];
      if(t_in_c[2] == 0)  or_out_t[j++] = or_t_WRT_c[2];
      if(t_in_c[3] == 0)  or_out_t[j++] = or_t_WRT_c[3];

      // The edge, whose endpoint are the 2 tet-verts taht dont lie
      // on the constraint-plane, do not cross the constraint-plane.
      if(or_out_t[0] == or_out_t[1])  return 2;

      return 1;
    }

    // One tetrahedron vertex belong to the triangle interior.
    if(num_t_in_c == 1) {

      uint32_t or_out_t[3];
      uint32_t j=0;
      if(t_in_c[0] == 0)  or_out_t[j++] = or_t_WRT_c[0];
      if(t_in_c[1] == 0)  or_out_t[j++] = or_t_WRT_c[1];
      if(t_in_c[2] == 0)  or_out_t[j++] = or_t_WRT_c[2];
      if(t_in_c[3] == 0)  or_out_t[j++] = or_t_WRT_c[3];

      // The edges, whose endpoint are 2 of the 3 tet-verts taht dont lie
      // on the constraint-plane, do not cross the constraint-plane.
      if(or_out_t[0] == or_out_t[1]  &&  or_out_t[0] == or_out_t[2] )  return 2;

      return 1;
    }


    // None of the tetrahedron vertices belong to the triangle interior.

    if( (or_t_WRT_c[0] != or_t_WRT_c[1])  &&
        vrt_innerSegmentCrossesInnerTriangle(t[0],t[1], c[0],c[1],c[2], mesh) ) return 1;

    if( (or_t_WRT_c[0] != or_t_WRT_c[2])  &&
        vrt_innerSegmentCrossesInnerTriangle(t[0],t[2], c[0],c[1],c[2], mesh) ) return 1;

    if( (or_t_WRT_c[1] != or_t_WRT_c[2])  &&
        vrt_innerSegmentCrossesInnerTriangle(t[1],t[2], c[0],c[1],c[2], mesh) ) return 1;

    if( (or_t_WRT_c[0] != or_t_WRT_c[3])  &&
        vrt_innerSegmentCrossesInnerTriangle(t[0],t[3], c[0],c[1],c[2], mesh) ) return 1;

    if( (or_t_WRT_c[1] != or_t_WRT_c[3])  &&
        vrt_innerSegmentCrossesInnerTriangle(t[1],t[3], c[0],c[1],c[2], mesh) ) return 1;

    if( (or_t_WRT_c[2] != or_t_WRT_c[3])  &&
        vrt_innerSegmentCrossesInnerTriangle(t[2],t[3], c[0],c[1],c[2], mesh) ) return 1;

    return 0;
}

//-----------------------
// COORDIANTING FUNCTIONS
//-----------------------

//  Input: pointer to the mesh,
//         the vertices of the constraint: constraint_vrts,
//         pointer to a tetrahedron index type: num_intersecatedTet,
//         pointer to array tetrahedra indices: intersecatedTet,
//         array of tetrahedra marker: mark_TetIntersection.
// Output: by using num_intersecatedTet returns the number of tetrahedra
//           intersecated by the boundary of the constraint-triangle,
//         by using intersecatedTet returns the indices of tetrahedra
//           intersecated by the boundary of the constraint-triangle,
//         by using mark_TetIntersection marks the intersection as explained
//           at the beginning of the function insert_constraints.
void intersections_constraint_sides(TetMesh* mesh, const uint32_t* constraint_vrts,
                                    uint64_t* num_intersecatedTet,
                                    uint64_t** intersecatedTet,
                                    uint32_t* mark_TetIntersection){

    // Cycle over the 3 sides of the constraints.
    for(uint32_t constr_side=0; constr_side<3; constr_side++){

        // Endpoints indices
        uint32_t v_start = constraint_vrts[ constr_side ];
        uint32_t v_stop  = constraint_vrts[ (constr_side+1)%3 ];
        uint32_t other_constr_vrt = constraint_vrts[ (constr_side+2)%3 ];

        #ifdef DEBUG
        printf("--> Analysing side (%u,%u) of the constraint.\n",v_start,v_stop);
        #endif

        // Array to store the information relative to the
        // constraint side travelling (Growing-region).
        // It has 4 elements:
        // - element_0 -> takes values in the set {1,2,3};
        // - if element_0 is 1:
        //      element_1 is the index of the vertex that intersect the
        //      constarint side in the travelling direction (v_start -> v_stop);
        // - if element_0 is 2:
        //      element_1 and element_2 are the indices of two vertices that
        //      are the endpoints of a segment (shared by two tetrahedra) whose
        //      interior is intersecated by the constarint side along travelling
        //      direction (v_start -> v_stop);
        // - if element_0 is 3,
        //      element_1, element_2 and element_3 are the indices of three
        //      vertices that define a face (of one tetrahedron) whose interior
        //      is intersecated by the constarint side along travelling
        //      direction (v_start -> v_stop);
        uint32_t connecting_vrts[4];
        uint64_t nextTet_ind;

        //-------------
        // BEGIN-phase: begin travelling along the constraint side searching
        //-------------   between tetrahedra incident in v_start.
        // Note. All tetrahedra incident in v_start intersect the constraint
        //       (in v_start). Those goes in found_tet.
        uint64_t num_found_tet = 0;
        uint64_t* found_tet =
                intersections_TetVrtOnConstraintSide(mesh, v_start, v_stop,
                                                     other_constr_vrt,
                                                     mark_TetIntersection,
                                                     connecting_vrts,
                                                     &nextTet_ind,
                                                     &num_found_tet);

        enqueueTetsArray(found_tet, num_found_tet,
                         intersecatedTet, num_intersecatedTet);
        free(found_tet);
        found_tet = NULL;

        //----------------
        // CONTINUE-phase: continue travelling along constraint side searching
        //----------------   between tetrahedra in the growing-up region.
         while ( connecting_vrts[1] != v_stop ) {

            uint64_t num_found_tet = 0;
            uint64_t* found_tet = NULL;

            uint32_t tet_edge[2], tet_face[3];

            switch (connecting_vrts[0]) {

              case 1:
                // Tetrahedra incident in connecting_vrts[1]
                // intersect the constraint side.
                found_tet = intersections_TetVrtOnConstraintSide(mesh,
                                  connecting_vrts[1], v_stop, other_constr_vrt,
                                         mark_TetIntersection, connecting_vrts,
                                                &nextTet_ind, &num_found_tet);
                break;

              case 2:
                // Tetrahedra incident in <connecting_vrts[1],connecting_vrts[2]>
                // intersect the constraint side.
                tet_edge[0]=connecting_vrts[1];
                tet_edge[1]=connecting_vrts[2];
                found_tet = intersections_TetEdgeCrossConstraintSide(mesh,
                                   tet_edge, v_start, v_stop, other_constr_vrt,
                                         mark_TetIntersection, connecting_vrts,
                                                &nextTet_ind, &num_found_tet);
                break;

              case 3:
                // Tetrahedra incident in
                // <connecting_vrts[1],connecting_vrts[2],connecting_vrts[3]>
                // intersect the constraint side.
                tet_face[0]=connecting_vrts[1];
                tet_face[1]=connecting_vrts[2];
                tet_face[2]=connecting_vrts[3];
                uint32_t to_add = intersections_TetFacePiercedConstraintSide(
                                                                          mesh,
                                   tet_face, v_start, v_stop, other_constr_vrt,
                                         mark_TetIntersection, connecting_vrts,
                                                                &nextTet_ind);
                if(to_add){
                  num_found_tet = 1;
                  found_tet = (uint64_t*) malloc(sizeof(uint64_t));
                  found_tet[0] = nextTet_ind;
                }
                break;
             }

             if(num_found_tet>0){
               enqueueTetsArray(found_tet, num_found_tet,
                                intersecatedTet, num_intersecatedTet);
               free(found_tet);
               found_tet = NULL;
             }
        }

    } // END Cycle over the 3 sides of the constraints.

}


//  Input: pointer to the mesh,
//         indices of a constraint-tiangle vertices: constraints_vrts,
//         index of a tetrahedron tet: tet_ind,
//         pointer to the array of tetrahedra marker: mark_TetIntersection.
// Output: return 1 if tet is a not-already-visited tetrahedron AND
//         intersects the interior of the constraint, 0 otherwise.
static inline uint32_t constrInterior_found(const TetMesh* mesh,
                                            const uint32_t* constraints_vrts,
                                            uint64_t tet_ind,
                                            uint32_t* mark_TetIntersection){

    // Check if the tetrahedron has been already visited.
    if( mark_TetIntersection[tet_ind] != 0 )    return 0;

    uint64_t f_ID = UINT64_MAX;
    uint32_t result = tet_intersects_triInterior(tet_ind, constraints_vrts, mesh, &f_ID);

    if(result == 0){ // Tetrahedron does not intersects the constraint interior.
      return 0;
    }
    if( result == 1 ){  // Constraint interior intersects tetrahedron interior.
      mark_TetIntersection[tet_ind] = IMPROPER_INTERSECTION_COUNTED;
    }
    if( result == 3 ){  // A tet-face (2D)overlaps with the constraint.

      if(f_ID==0) mark_TetIntersection[tet_ind] = OVERLAP2D_F0_COUNTED;
      else if(f_ID==1) mark_TetIntersection[tet_ind] = OVERLAP2D_F1_COUNTED;
      else if(f_ID==2) mark_TetIntersection[tet_ind] = OVERLAP2D_F2_COUNTED;
      else if(f_ID==3) mark_TetIntersection[tet_ind] = OVERLAP2D_F3_COUNTED;
    }
    if(result==2){ // result==2 Constraint interior intersects only tetrahedron boundary.
      mark_TetIntersection[tet_ind] = PROPER_INTERSECTION_COUNTED;
    }
    return 1;
}

//  Input: pointer to the mesh,
//         indices of a constraint-tiangle vertices: constraints_vrts,
//         index of a tetrahedron tet intersecting the constraint: bnd_tet_ind,
//         pointer to the array of tetrahedra marker: mark_TetIntersection,
//         pointer to tetrahedra index type: adj_tet_ind_return.
// Output: return 1 if exist a not-already-visited tetrahedron adjacent to tet
//         that intersects the interior of the constraint, 0 otherwise;
//         by using adj_tet_ind_return returns the index of the tetrahedron
//         if 1 is rerurned.
uint32_t constrInterior_firstStep(const TetMesh* mesh,
                                              const uint32_t* constraints_vrts,
                                                uint64_t bnd_tet_ind,
                                                uint32_t* mark_TetIntersection,
                                                uint64_t* adj_tet_ind_return){
    // Cycle over the tetrahedra adjacent to bnd_tet.
    for(uint64_t i=0; i<4; i++){

        uint64_t adj_tet_ind = mesh->tet_neigh[4*bnd_tet_ind + i]>>2;

        // #ifdef DEBUG
        // printf("Looking at tetrahedron %llu adjacent to %llu.\n",
        //        adj_tet_ind, bnd_tet_ind);
        // #endif

        if(constrInterior_found(mesh, constraints_vrts, adj_tet_ind, mark_TetIntersection) ){
          *adj_tet_ind_return = adj_tet_ind;

          // #ifdef DEBUG
          // printf("(First Step) This tetrahedron intersects ONLY the "
          //        "interior of the constraint.\n");
          // #endif

          return 1;
        }

    }

    return 0;
}

//  Input: pointer to the mesh,
//         indices of a constraint-tiangle vertices: constraints_vrts,
//         index of a tetrahedron tet intersecting the constraint: tet_ind,
//         pointer to the array of tetrahedra marker: mark_TetIntersection.
// Output: return the number of not-already-visited tetrahedron that intersects
//         the interior of the constraint in a connected region limitated by
//         the already-visited tetrahedra (all that intersects constraint
//         boundary + all the alraedy-visited of those that intersectcontraint
//         interior).
uint64_t constrInterior_count(const TetMesh* mesh,
                              const uint32_t* constraints_vrts,
                              uint64_t tet_ind,
                              uint32_t* mark_TetIntersection){
    uint64_t num_tets_intersects = 1;   // tet intersects (only)
                                        // the interior of the constraint.

    // Cycle over the tetrahedra adjacent to tet.
    for(uint64_t i=0; i<4; i++){

        uint64_t adj_tet_ind = mesh->tet_neigh[4*tet_ind + i]>>2;

        if(constrInterior_found(mesh, constraints_vrts, adj_tet_ind, mark_TetIntersection) )
          num_tets_intersects += constrInterior_count(mesh, constraints_vrts,
                                                      adj_tet_ind,
                                                      mark_TetIntersection);
    }

    return num_tets_intersects;
}

//  Input: pointer to the mesh,
//         index of a tetrahedron tet intersecting the constraint: tet_ind,
//         pointer to the array of tetrahedra marker: mark_TetIntersection,
//         pointer to the current last-position (length-1) of tet_list: pos,
//         pointer to an array of tetrahedra index type: tet_list.
// Output: by using pos returns the current last-position (length-1) of tet_list.
// Note. this recursive function enqueue to tet_list each tetrahedra adjacent
//       to tet that intersects the only the interior of the constraint.
void constrInterior_save(const TetMesh* mesh, uint64_t tet_ind,
                         uint32_t* mark_TetIntersection, uint64_t* pos,
                         uint64_t* tet_list){

    tet_list[*pos] = tet_ind;
    (*pos)++;

    if(mark_TetIntersection[tet_ind] == IMPROPER_INTERSECTION_COUNTED)
      mark_TetIntersection[tet_ind] = IMPROPER_INTERSECTION;
    else if(mark_TetIntersection[tet_ind] == OVERLAP2D_F0_COUNTED)
      mark_TetIntersection[tet_ind] = OVERLAP2D_F0;
    else if(mark_TetIntersection[tet_ind] == OVERLAP2D_F1_COUNTED)
      mark_TetIntersection[tet_ind] = OVERLAP2D_F1;
    else if(mark_TetIntersection[tet_ind] == OVERLAP2D_F2_COUNTED)
      mark_TetIntersection[tet_ind] = OVERLAP2D_F2;
    else if(mark_TetIntersection[tet_ind] == OVERLAP2D_F3_COUNTED)
      mark_TetIntersection[tet_ind] = OVERLAP2D_F3;
    else
      mark_TetIntersection[tet_ind] = INTERSECTION;

    for(uint64_t i=0; i<4; i++){
        uint64_t adj_tet_ind = mesh->tet_neigh[4*tet_ind + i]>>2;

      if(mark_TetIntersection[adj_tet_ind] == IMPROPER_INTERSECTION_COUNTED ||
         mark_TetIntersection[adj_tet_ind] == OVERLAP2D_F0_COUNTED          ||
         mark_TetIntersection[adj_tet_ind] == OVERLAP2D_F1_COUNTED          ||
         mark_TetIntersection[adj_tet_ind] == OVERLAP2D_F2_COUNTED          ||
         mark_TetIntersection[adj_tet_ind] == OVERLAP2D_F3_COUNTED          ||
         mark_TetIntersection[adj_tet_ind] == PROPER_INTERSECTION_COUNTED     )
          constrInterior_save(mesh, adj_tet_ind, mark_TetIntersection, pos,
                              tet_list);
    }
    return;
}

//  Input: pointer to the mesh,
//         the vertices of the constraint: constraint_vrts,
//         pointer to a tetrahedron index type: num_intersecatedTet,
//         pointer to array tetrahedra indices: intersecatedTet,
//         pinter to a tetrahedra marker: mark_TetIntersection.
// Output: by using num_intersecatedTet returns the sum between the number of
//         tetrahedra intersecated by the boundary and the number of tetrahedra
//         intersecated by the interior of the constraint-triangle,
//         by using intersecatedTet returns the indices of tetrahedra
//         intersecated by the interior of the constraint-triangle enqueued to
//         the array of the indices of tetrahedra intersecated by the boundary.
//         by using mark_TetIntersection marks the intersection as explained
//         at the beginning of the function insert_constraints.
void intersections_constraint_interior(TetMesh* mesh, const uint32_t* constraint_vrts,
                                       uint64_t* num_intersecatedTet,
                                       uint64_t** intersecatedTet,
                                       uint32_t* mark_TetIntersection){

    uint64_t num_bnd_tets = *num_intersecatedTet;
    uint64_t bnd_tet;

    for(uint64_t side_tet_ind=0; side_tet_ind< num_bnd_tets; side_tet_ind++){

        bnd_tet = (*intersecatedTet)[side_tet_ind];

        // First-step: consider a tetrahedron on the boundary (bnd_tet) and
        //             search if exists any adjacent tetrahedron which
        //             intersects the constraint.
        //             If NOT pass to next tetrahedron visited on the boundary.
        uint64_t adjIN_tet;

        if(constrInterior_firstStep(mesh, constraint_vrts, bnd_tet,
                                    mark_TetIntersection, &adjIN_tet) == 0 )
            continue;

        // RecursiveFun-call:
        // explore a connected region of the interior of the constraint -> COUNT
        uint64_t num = constrInterior_count(mesh, constraint_vrts, adjIN_tet,
                                            mark_TetIntersection);

        if(num>1){

          uint64_t* interiorConstrInterct_tet = (uint64_t*) malloc(sizeof(uint64_t) * num);
          // RecursiveFun-call:
          // explore a connected region of the constraint interior -> SAVE
          uint64_t pos = 0;
          constrInterior_save(mesh, adjIN_tet, mark_TetIntersection, &pos,
                              interiorConstrInterct_tet);

          enqueueTetsArray(interiorConstrInterct_tet, num,
                          intersecatedTet, num_intersecatedTet);
          free(interiorConstrInterct_tet);
          interiorConstrInterct_tet = NULL;
        }
        else{
          if(mark_TetIntersection[adjIN_tet] == IMPROPER_INTERSECTION_COUNTED)
            mark_TetIntersection[adjIN_tet] = IMPROPER_INTERSECTION;
          else if(mark_TetIntersection[adjIN_tet] == OVERLAP2D_F0_COUNTED)
            mark_TetIntersection[adjIN_tet] = OVERLAP2D_F0;
          else if(mark_TetIntersection[adjIN_tet] == OVERLAP2D_F1_COUNTED)
            mark_TetIntersection[adjIN_tet] = OVERLAP2D_F1;
          else if(mark_TetIntersection[adjIN_tet] == OVERLAP2D_F2_COUNTED)
            mark_TetIntersection[adjIN_tet] = OVERLAP2D_F2;
          else if(mark_TetIntersection[adjIN_tet] == OVERLAP2D_F3_COUNTED)
            mark_TetIntersection[adjIN_tet] = OVERLAP2D_F3;
          else
            mark_TetIntersection[adjIN_tet] = INTERSECTION;

          enqueueTetsArray(&adjIN_tet, 1, intersecatedTet, num_intersecatedTet);
        }
    }
}

//
//
static inline void compile_map_innerInt(uint32_t tri_ind, uint64_t tet_ind,
                                        uint32_t* num_map, uint32_t** map  ){
  num_map[tet_ind]++;  // Number of intersections increases for the tetrahedron.

  if(num_map[tet_ind] == 1) // 1st intersection for the tetrahedron.
    map[tet_ind] = (uint32_t*) malloc( sizeof(uint32_t) );
  else            // New intersection has to be added to the existing ones.
    map[tet_ind] = (uint32_t*) realloc( map[tet_ind], sizeof(uint32_t)*num_map[tet_ind] );

  map [ tet_ind ] [ num_map[tet_ind]-1 ] = tri_ind;
}

//  Input:
// Output:
static inline void compile_map_fi(uint32_t tri_ind, uint64_t tet_ind,
                                     uint32_t* num_map_fi, uint32_t** map_fi){

    num_map_fi[tet_ind]++;   // The number of overlapping constraints increases
                             // for the tetrahedron tet_ind face_i.

    if(num_map_fi[tet_ind] == 1)   // It is the first overlapping constraint.
      map_fi[tet_ind] = (uint32_t*) malloc( sizeof(uint32_t) );
    else           // The new overlapping constraints has to be added
                   // to the existing ones.
      map_fi[tet_ind] = (uint32_t*) realloc( map_fi[tet_ind],
                                          sizeof(uint32_t)*num_map_fi[tet_ind] );


    map_fi [ tet_ind ] [ num_map_fi[tet_ind]-1 ] = tri_ind;
}

//  Input:
// Output:
static inline void compile_tetfFaces_map(uint32_t tri_ind, uint64_t tet_face_ind,
                                         uint32_t* num_map_f0, uint32_t** map_f0,
                                         uint32_t* num_map_f1, uint32_t** map_f1,
                                         uint32_t* num_map_f2, uint32_t** map_f2,
                                         uint32_t* num_map_f3, uint32_t** map_f3  ){
  uint64_t i = tet_face_ind%4;        //tet_face_ind%4
  uint64_t tet_ind = tet_face_ind/4; //tet_face_ind/4
  if(i==0)      compile_map_fi(tri_ind, tet_ind, num_map_f0, map_f0);
  else if(i==1) compile_map_fi(tri_ind, tet_ind, num_map_f1, map_f1);
  else if(i==2) compile_map_fi(tri_ind, tet_ind, num_map_f2, map_f2);
  else if(i==3) compile_map_fi(tri_ind, tet_ind, num_map_f3, map_f3);
}

//  Input: index of the constraint tri: tri_ind,
//         number of tetrahedra that intersect the constraint tri: n,
//         array of indices of the tetrahedra that intersect the
//         constraint tri: tets,
//         array of marker to distiguish between general intersection and
//         intersection that have to be mapped,
//         array, i-th elemet counts the number of constraints intersecated
//         by i-th tetrahedron: num_map,
//         array, i-th elemet points to the array listing the constraints
//         intersecated by i-th tetrahedron: map.
// Output: by using map returns the map updated with information of tets,
//         by using num_map returns the num_map updated with
//         information of tets.
// Note1. Update means that, at the lists of map relative to tetrahedra in tets
//        is added tri_ind, consequently the relative num_map are incrementated.
// Note2. The entries of the array of marker have to be set to 0 after their
//        iformation have been used.
void compile_maps(uint32_t tri_ind, uint64_t n,
                               const uint64_t* tets,
                               uint32_t* mark_TetIntersection,
                               uint32_t* num_map, uint32_t** map,
                               uint32_t* num_map_f0, uint32_t** map_f0,
                               uint32_t* num_map_f1, uint32_t** map_f1,
                               uint32_t* num_map_f2, uint32_t** map_f2,
                               uint32_t* num_map_f3, uint32_t** map_f3 ){

  for(uint64_t i=0; i<n; i++){
      uint64_t tet_ind = tets[i];

      if(mark_TetIntersection[tet_ind] == IMPROPER_INTERSECTION)
        compile_map_innerInt(tri_ind, tet_ind, num_map, map);
      else if(mark_TetIntersection[tet_ind] == OVERLAP2D_F0)
        compile_map_fi(tri_ind, tet_ind, num_map_f0, map_f0);
      else if(mark_TetIntersection[tet_ind] == OVERLAP2D_F1)
        compile_map_fi(tri_ind, tet_ind, num_map_f1, map_f1);
      else if(mark_TetIntersection[tet_ind] == OVERLAP2D_F2)
        compile_map_fi(tri_ind, tet_ind, num_map_f2, map_f2);
      else if(mark_TetIntersection[tet_ind] == OVERLAP2D_F3)
        compile_map_fi(tri_ind, tet_ind, num_map_f3, map_f3);

      mark_TetIntersection[tet_ind] = 0;  // Reset the tetrahedron marker.
  }
}

/***********************************/
/** Constraints insertion GENERAL **/
/***********************************/
// This function explores the mesh and traces the intersections between the
// constraints-triangles and the mesh-tetrahdra.
// This procedure is essential in order to perfor later the BSP subdivision.
// Some "maps" are created here to memory:
//  1- what are the constraints that IMPROPERLY intersects a certain tetrahedon,
//     unless those that are coplanar with a tetrahedon-face.
//     We use 1 map: it has as many elements as the number of tet_
//     If one constraint IMPROPERLY intersects the i-th tetrahedon (whithout
//     having a coplanar tetrahedon-face), than the constraint_ind is
//     saved in map[i].
//  2- what are the constraints that intersects a certain tetrahedon-face, and
//     the intersection have a non-zero area.
//     To this end 4 maps are created: map_f0, map_f1, map_f2, map_f3,
//     each map have as many elements as the number of tetrahedra, and refers
//     to the tetrahedron-face opposite to vertex 0,1,2 or 3 respectivelly.
//     If the face j of the i-th tetrahedron is coplanar with a constraint and
//     their intersection have a non-zero area, than the constraint_ind is
//     saved in map_fj[i].
// Note. By IMPROPER we mean that the intersection between a constraint and a
//       tetrahedron is not PROPER.
//       A PROPER intersection occours when the intersection is a
//       sub-simplex of the tetrahedron, i.e.
//        - ONLY a face of the tetrahedron,
//        - ONLY an edge of the tetrahedron,
//        - ONLY a vertex of the tetrahedron.

void insert_constraints(TetMesh* mesh, constraints_t* constraints,
                        uint32_t* num_map, uint32_t** map,
                        uint32_t* num_map_f0, uint32_t** map_f0,
                        uint32_t* num_map_f1, uint32_t** map_f1,
                        uint32_t* num_map_f2, uint32_t** map_f2,
                        uint32_t* num_map_f3, uint32_t** map_f3  ){

    // We will cycle over constraints using an array of marker to mark the
    // tetrahedra that intersect a constraint.
    // We will distinguish between:
    // - no intersection (after visit) or not already visited (before visit)
    //      -> marked as 0,
    // - generic intersection (proper or improper,
    //    unless those of [The very trivial case] - see below)
    //      -> marked as INTERSECTION (1),
    // - improper intersection
    //      -> marked as IMPROPER_INTERSECTION (2).
    // - faces (2D)overlapping with constraints, wherever its the face opposite
    //   to vertex j=0,1,2 or 3...
    //      -> marked as OVERLAP2D_Fj
    uint32_t* mark_TetIntersection = (uint32_t*) calloc(mesh->tet_num, sizeof(uint32_t));

    // Search interections on each constraint.
    for(uint32_t tri_ind=0; tri_ind<constraints->num_triangles; tri_ind++){

      uint32_t v[3]; // vertices of the constraint-triangle.
      uint32_t tri_ID = 3*tri_ind;
      v[0] = constraints->tri_vertices[tri_ID    ];
      v[1] = constraints->tri_vertices[tri_ID + 1];
      v[2] = constraints->tri_vertices[tri_ID + 2];

      // ---STEP 0--- [The very trivial case]
      // Check if the constraint-triangle is a face of a tetrahedron
      // incident in v0. In this case there is a proper intersection,
      // but there is a (2D)overlapping.
      uint64_t tet_face_ind = UINT64_MAX;
      if(triangle_in_VT(v[0], v[1], v[2], mesh, &tet_face_ind)==1 ){
         const uint64_t ng = mesh->tet_neigh[tet_face_ind];
         compile_tetfFaces_map(tri_ind, ng,
                               num_map_f0, map_f0, num_map_f1, map_f1,
                               num_map_f2, map_f2, num_map_f3, map_f3);
         compile_tetfFaces_map(tri_ind, tet_face_ind,
                               num_map_f0, map_f0, num_map_f1, map_f1,
                               num_map_f2, map_f2, num_map_f3, map_f3  );

         continue;
      }

      // We use an array to save the indices all tetrahedra that intersect
      // the constraint, properly or improperly.
      uint64_t num_intersecatedTet = 0;
      uint64_t* intersecatedTet = NULL;

      // ---STEP 1--- [Intersections with the BOUNDARY of the constraint]
      intersections_constraint_sides(mesh, v, &num_intersecatedTet,
                                              &intersecatedTet,
                                              mark_TetIntersection);

      // ---STEP 2--- [Search for improper intersections]
      find_improperIntersection(v, intersecatedTet, num_intersecatedTet,
                                   mark_TetIntersection, mesh);

      // ---STEP 3--- [Intersections with the constraint INTERIOR]
      intersections_constraint_interior(mesh, v, &num_intersecatedTet,
                                                 &intersecatedTet,
                                                 mark_TetIntersection);

      // ---STEP 4--- [Fill intersection map & reset mark_TetIntersection]
      compile_maps(tri_ind, num_intersecatedTet, intersecatedTet, mark_TetIntersection,
                   num_map, map, num_map_f0, map_f0, num_map_f1, map_f1,
                   num_map_f2, map_f2, num_map_f3, map_f3);
      free(intersecatedTet);
      intersecatedTet = NULL;
    }

    free(mark_TetIntersection);
    mark_TetIntersection = NULL;
}
