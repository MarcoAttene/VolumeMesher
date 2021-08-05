#include <iostream>
#include <fstream>
#include <set>
#include <time.h>
#include <algorithm>
#include "BSP.h"
#include "delaunay.h"

#define OPPOSITE_SIGNE(a,b) (a<0 && b>0) || (a>0 && b<0)
#define MIN_VECT_ELEM(v,n,it,i_min)  it=1; i_min=0; do{ if(v[it]<v[i_min]) i_min=it; it++; }while(it<n)
#define SAME_EDGE_ENDPTS(e1,e2,E1,E2)  ((e1==E1 && e2==E2) || (e1==E2 && e2==E1))
#define CONSECUTIVE_EDGES(v1,v2,u1,u2) (v1==u1 || v1==u2 || v2==u1 || v2==u2) // assuming <u1,u2> != <v1,v2>
#define FIND_VECT_POS(e,v,pos) pos=0;  do{ if(v[pos]==e) break; pos++; }while(pos<v.size()) // !! assumes e in v. !!
#define REMOVE_ELEM_VECT(e,v) v.erase(std::find(v.begin(), v.end(), e))
#define IS_GHOST_CELL(c) (c==UINT64_MAX)
#define IS_GHOST_TET(t) (mesh->tet_node[4*t+3]==UINT32_MAX)

//----------------
// General purpose
//----------------

// Input: vector of uint64_t type elements: vect,
//        number of elements to shift: num_shift.
// Output: by using vect returns the original vector shifted-down (i.e.
//         according to increasing order) of num_shift position.
// EX. vect = {3,7,2,5,2} num_shift = 2 -> vect = {5,2,3,7,2}.
inline void UINT64_vect_down_shift(vector<uint64_t>& vect, uint64_t num_shift){
  vector<uint64_t> tmp(num_shift, UINT64_MAX);
  uint64_t shift_start_pos = vect.size() - num_shift;
  for(uint64_t pos=0; pos<num_shift; pos++)
     tmp[pos] = vect[shift_start_pos + pos];
  for(uint64_t pos=shift_start_pos-1; pos>=0; pos--){
     vect[pos + num_shift] = vect[pos];
     if(pos==0) break; // Otherwise pos(uint64_t) becomes negative, i.e. error!!
   }
  for(uint64_t pos=0; pos<num_shift; pos++)  vect[pos] = tmp[pos];
}

//  Input: vector of uint64_t type elements: vect,
//        number of elements to shift: num_shift.
// Output: by using vect returns the original vector shifted-up (i.e.
//         according to decreasing order) of num_shift position.
// EX. vect = {3,7,2,5,2} num_shift = 2 -> vect = {3,7,2,5,2}.
inline void UINT64_vect_up_shift(vector<uint64_t>& vect, uint64_t num_shift){
  vector<uint64_t> tmp(num_shift, UINT64_MAX);
  uint64_t pos;
  for(pos=0; pos<vect.size(); pos++){
    if(pos<num_shift) tmp[pos] = vect[pos];
    else vect[pos - num_shift] = vect[pos];
  }
  uint64_t refill_pos = vect.size() - num_shift;
  for(pos=0; pos<num_shift; pos++)
    vect[refill_pos + pos] = tmp[pos];
}

//  Input: the endpoints of two CONSECUTIVE edges u = <u0,u1>, v=<v0,v1>
// Output: the common endpoint.
inline uint32_t consecEdges_common_endpt(uint32_t u0, uint32_t u1,
                                         uint32_t v0, uint32_t v1){
    if(u0 == v0 || u0 == v1) return u0;
    if(u1 == v0 || u1 == v1) return u1;
    // If goes here, something wrong.
    printf("\n[BSP.cpp]consecEdges_common_endpt: ERROR no match.\n");
    return UINT32_MAX;
}

//  Input: the endpoints of an edge u = <u0,u1>,
//         the common endpoint between an edge v CONSECUTIVE to u,
//         and u itself: w.
// Output: the endpoint of u different from v_comm.
inline uint32_t other_edge_endpt(uint32_t u0, uint32_t u1, uint32_t w){
    if(w == u0) return u1;
    if(w == u1) return u0;
    // If goes here, something wrong.
    printf("\n[BSP.cpp]other_edge_endpt_ind: ERROR no match.\n");
    return UINT32_MAX;

}

//----------------
// BSPedge Methods
//----------------

BSPedge BSPedge::split(uint32_t new_point){
    BSPedge e;
    if(meshVertices[2]==UINT32_MAX){
        e.meshVertices[0] = meshVertices[0];
        e.meshVertices[1] = meshVertices[1];
        e.meshVertices[2] = UINT32_MAX;
    }
    else{
        e.meshVertices[0] = meshVertices[0];
        e.meshVertices[1] = meshVertices[1];
        e.meshVertices[2] = meshVertices[2];
        e.meshVertices[3] = meshVertices[3];
        e.meshVertices[4] = meshVertices[4];
        e.meshVertices[5] = meshVertices[5];
    }
    e.vertices[0] = vertices[0];
    e.vertices[1] = new_point;
    vertices[0] = new_point;
    e.conn_face_0 = conn_face_0;
    return e;
}


//----------------
// BSPface Methods
//----------------
inline void BSPface::exchange_conn_cell(uint64_t cell,uint64_t newCell){

    #ifdef DEBUG_BSP
      if(conn_cells[0]!=cell && conn_cells[1]!=cell)
         printf("\n[BSP.h]BSPface::exchange_conn_cell: ERROR no macth for cell "
                "#%llu with this face.\n",cell);
    #endif

    if(conn_cells[0]==cell)   conn_cells[0]=newCell;
    else                      conn_cells[1]=newCell;
}

inline void BSPface::removeEdge(uint64_t edge_face_ind){
  // Special case: edge is the last element of edges.
  if(edge_face_ind == edges.size() -1) edges.pop_back();
  // General case: erase it, but later.. at the moment mark it
  else edges[edge_face_ind] = UINT64_MAX;
}

//----------------
// BSPcell Methods
//----------------

inline void BSPcell::removeFace(uint64_t face_pos){
  if(face_pos != faces.size()-1)  faces[face_pos] = faces[faces.size()-1];
  faces.pop_back();
}

//--------------------
// BSPcomplex Methods
//--------------------

//  Input: the index of a BSPedge: edge,
//         the index of a BSPface: face.
// Output: nothing.
// Note. BSPedge is attached to BSPface and viceversa.
// Note. this function is used ONLY during the conversion from Delaunay mesh to
//       BSPcomplex.
inline void BSPcomplex::assigne_edge_to_face(uint64_t edge, uint64_t face){
  faces[face].edges.push_back(edge);
  //edges[edge].conn_faces.push_back(face);
  edges[edge].conn_face_0 = face;
}

//  Input: the index of 2 cells w.r.t. the vector cells: c1, c2,
// Output: the index of the fece, w.r.t. the vector faces, that lies between
//         the two input-cells.
// Note. It is assumed that both the input-cells are bounded by the input-face.
uint64_t BSPcomplex::faceSharedWithCell(uint64_t c1, uint64_t c2){

  BSPcell& cell1 = cells[c1];
  for(uint64_t i=0; i<cell1.faces.size(); i++)
    if(faces[ cell1.faces[i] ].conn_cells[0] == c2 ||
       faces[ cell1.faces[i] ].conn_cells[1] == c2   )
       return cell1.faces[i];

  #ifdef DEBUG_BSP
  printf("\n[BSP.cpp]BSPcomplex::faceSharedWithCell: ERROR no common face "
         "between cell #%llu and cell #%llu.\n", c1, c2);
  #endif
  return UINT64_MAX;  // Never reached.
}

//  Input: a BSPcell: cell.
// Output: returns the number of edges of the BSPcell.
uint64_t BSPcomplex::count_cellEdges(const BSPcell& cell){
    uint64_t num_halfedges = 0;
    for(const uint64_t f : cell.faces)  num_halfedges += faces[f].edges.size();
    return num_halfedges / 2;
}

//  Input: a BSPcell: cell.
// Output: returns the number of edges of the BSPcell.
uint32_t BSPcomplex::count_cellVertices(const BSPcell& cell,
                                               uint64_t* num_cellEdges){
  if((*num_cellEdges) == UINT64_MAX)  (*num_cellEdges)=count_cellEdges(cell);
  // Euler formula: num_cellVrts = num_cellEdges + 2 - num_cellFaces
  return (uint32_t)((*num_cellEdges) + 2 - cell.faces.size());
}

//  Input: a BSPcell: cell,
//         vector of type edge index: cell_edges.
// Output: by using cell_edges returns the indices of
//         the edges (w.r.t. vector edges) of the BSPcell.
void BSPcomplex::list_cellEdges(BSPcell& cell, vector<uint64_t>& cell_edges){

    uint64_t edge_ind, ce_ind=0;
    for(uint64_t f=0; f<cell.faces.size(); f++){
        BSPface& face = faces[ cell.faces[f] ];
        for(uint64_t e=0; e<face.edges.size(); e++){
            edge_ind = face.edges[e];
            if(edge_visit[edge_ind]==0){
              cell_edges[ce_ind++] = edge_ind;
              edge_visit[edge_ind]=1;
            }
        }
    }

    // Reset edge_visit
    for(uint64_t e=0; e<cell_edges.size(); e++) edge_visit[ cell_edges[e] ]=0;
}

//  Input: a BSPcell: cell,
//         the number of edges of the BSPcell: num_cellEdges,
//         vector of type vertex index: cell_vrts.
// Output: by using cell_vrts returns the indices of
//         the vertices (w.r.t. vector vertices) of the BSPcell.
void BSPcomplex::list_cellVertices(BSPcell& cell, uint64_t num_cellEdges,
                                   vector<uint32_t>& cell_vrts){

    vector<uint64_t> cell_edges(num_cellEdges, UINT64_MAX);
    list_cellEdges(cell, cell_edges);

    uint32_t v, cv_ind=0;
    for(uint64_t e=0; e<cell_edges.size(); e++){
        BSPedge& edge = edges[ cell_edges[e] ];
        v = edge.vertices[0];
        if(vrts_visit[v]==0){
          cell_vrts[cv_ind++] = v;
          vrts_visit[v]=1;
        }
        v = edge.vertices[1];
        if(vrts_visit[v]==0){
          cell_vrts[cv_ind++] = v;
          vrts_visit[v]=1;
        }
    }

    // Reset vrts_visit
    for(uint32_t u=0; u<cell_vrts.size(); u++) vrts_visit[ cell_vrts[u] ]=0;
}

//  Input: a BSPface: face,
//         vector of type vertex index: face_vrts.
// Output: by using face_vrts returns the indices of
//         the vertices (w.r.t. vector vertices) of the BSPface.
void BSPcomplex::list_faceVertices(BSPface& face, vector<uint32_t>& face_vrts){

    uint32_t fv_ind=0;

    // Add both endpoints of first edge.
    BSPedge& edge0 = edges[ face.edges[0] ];
    uint32_t e0 = edge0.vertices[0];
    uint32_t e1 = edge0.vertices[1];
    face_vrts[fv_ind++] = e0;
    face_vrts[fv_ind++] = e1;

    // Find the common endpoint between first and second edge.
    BSPedge& edge1 = edges[ face.edges[1] ];
    uint32_t link_vrt = consecEdges_common_endpt(e0, e1, edge1.vertices[0],
                                                         edge1.vertices[1]);
    if(link_vrt==e0){
      face_vrts[0] = e1;
      face_vrts[1] = e0;
    }

    // Walk the boundary and add the endpoint != link_vrt.
    for(uint64_t e=1; e<face.edges.size()-1; e++){
      BSPedge& edge = edges[ face.edges[e] ];

      if(link_vrt == edge.vertices[0]) link_vrt = edge.vertices[1];
      else                             link_vrt = edge.vertices[0];

      face_vrts[fv_ind++] = link_vrt;
    }
}

//  Input: a BSPcell: cell,
//         vector of edges indices type: cell_edges,
//         vector of vertices indices type: cell_vrts.
// Output: by using cell_edges returns the indices (w.r.t. vector edges)
//         of cell's edges,
//         by using cell_vrts returns the indices (w.r.t. vector vertices)
//         of cell's vertices,
void BSPcomplex::fill_cell_locDS(BSPcell& cell, vector<uint64_t>& cell_edges,
                                 vector<uint32_t>& cell_vrts){

  uint64_t edge_ind, ce_ind=0, cv_ind=0;
  uint32_t e0, e1;
  for(uint64_t f=0; f<cell.faces.size(); f++){
    BSPface& face = faces[ cell.faces[f] ];
    for(uint64_t e=0; e<face.edges.size(); e++){
        edge_ind = face.edges[e];
        BSPedge& edge = edges[edge_ind];

        if(edge_visit[edge_ind]==0){
          // Fill cell_edges.
          cell_edges[ce_ind++] = edge_ind;
          edge_visit[edge_ind]=1;

          e0 = edge.vertices[0];
          e1 = edge.vertices[1];

          // Fill cell_vrts.
          if(vrts_visit[e0]==0){
            cell_vrts[cv_ind++] = e0;
            vrts_visit[e0]=1;
          }
          if(vrts_visit[e1]==0){
            cell_vrts[cv_ind++] = e1;
            vrts_visit[e1]=1;
          }
        }

    }
  }

  // Reset edge_visit and vrts_visit.
  for(uint64_t e=0; e<cell_edges.size(); e++) edge_visit[cell_edges[e]]=0;
  for(uint64_t v=0; v<cell_vrts.size(); v++) vrts_visit[cell_vrts[v]]=0;
}

//  Input: a BSPface: face,
//         indices of two edge endpoints (w.r.t. vector vertices): u, v,
// Output: if exists returns the index of a BSPedge with endpoints u and v
//         that belong to faces[face], otherwise return UINT64_MAX.
inline uint64_t BSPcomplex::find_face_edge(const BSPface& face,
                                           uint32_t v, uint32_t u){

  #ifdef DEBUG_BSP
  if(face.edges.size()==0)
    printf("\n[BSP.cpp]find_face_edge: ERROR face have no edges.\n");
  #endif

  uint64_t edge_ind;
  for(uint64_t e=0; e<face.edges.size(); e++){
    edge_ind = face.edges[e];
    BSPedge& edge = edges[edge_ind];
    if(SAME_EDGE_ENDPTS(u, v, edge.vertices[0], edge.vertices[1]) )
      return edge_ind;
  }

  #ifdef DEBUG_BSP
  printf("[BSP.cpp]BSPcomplex::find_face_edge: "
         "WARNING no match for edge <%u,%u> with this face\n", u, v);
  #endif

  return UINT64_MAX;  // Never reached if <u,v> belong to the BSPcell.
}

//
//
uint64_t BSPcomplex::count_cellFaces_inc_cellVrt(const BSPcell& cell, uint32_t v){
  uint64_t k=0;
  for(const uint64_t fi : cell.faces)
      for(const uint64_t ei : faces[fi].edges)
          if(edges[ei].vertices[0] == v || edges[ei].vertices[1] == v) k++;
  return k / 2;
}

//
//
void BSPcomplex::cell_VFrelation(const BSPcell& cell, uint32_t v,
                                 vector<uint64_t>& v_incFaces_ind){
  uint64_t k=0;
  for(const uint64_t fi : cell.faces)
      for(const uint64_t ei : faces[fi].edges)
          if(edges[ei].vertices[0] == v || edges[ei].vertices[1] == v)
          {
            v_incFaces_ind[k++] = fi;
            break;
          }
}

//
//
void BSPcomplex::COMPL_cell_VFrelation(const BSPcell& cell, uint32_t v,
                                      vector<uint64_t>& v_NOT_incFaces_ind){
    uint64_t k = 0;
    for(const uint64_t fi : cell.faces){
        bool has_edge = false;
        for(const uint64_t ei : faces[fi].edges)
            if(edges[ei].vertices[0] == v || edges[ei].vertices[1] == v){
              has_edge = true;
              break;
            }
        if(!has_edge) v_NOT_incFaces_ind[k++] = fi;
    }
}

//
//
bool BSPcomplex::is_virtual(uint32_t constr_ind){
  return (constr_ind >= first_virtual_constraint);
}

// -- Geometric Predicates ---------------------------------------------------

// Returns TRUE if v either is a vertex of the plane p0, p1, p2, or it was built
// by intersecting such a plane with other simplexes.
bool isVertexBuiltFromPlane(const genericPoint* v,
                            const explicitPoint3D *p0,
                            const explicitPoint3D *p1,
                            const explicitPoint3D *p2){
    if(v->isExplicit3D()){
      return ((v == p0) || (v == p1) || (v == p2));
    }
    else if(v->isLPI()){
          const implicitPoint3D_LPI& lpi = v->toLPI();
          return (
            (((p0 == &lpi.P()) && (p1 == &lpi.Q())) || ((p1 == &lpi.P()) && (p0 == &lpi.Q()))) ||
            (((p1 == &lpi.P()) && (p2 == &lpi.Q())) || ((p2 == &lpi.P()) && (p1 == &lpi.Q()))) ||
            (((p2 == &lpi.P()) && (p0 == &lpi.Q())) || ((p0 == &lpi.P()) && (p2 == &lpi.Q()))) ||
            ((p0 == &lpi.R()) && (p1 == &lpi.S()) && (p2 == &lpi.T()))
            );
          }
          else{
            // TPI
            const implicitPoint3D_TPI& tpi = v->toTPI();
            return (((p0 == &tpi.V1()) && (p1 == &tpi.V2()) && (p2 == &tpi.V3())) ||
                   ((p0 == &tpi.W1()) && (p1 == &tpi.W2()) && (p2 == &tpi.W3())) ||
                   ((p0 == &tpi.U1()) && (p1 == &tpi.U2()) && (p2 == &tpi.U3())));
          }
}


//  Input: vector of vertices indices (w.r.t. vector vertices): vrts_inds,
//         3 indices of vertices that define a plane: plane_pt0, plane_pt1,
//         plane_pt2.
// Output: by using the global vector vrts_orBin returns the orientations of
//         each point of vrts_inds w.r.t. the plane for
//         {plane_pt0, plane_pt1, plane_pt2}.
void BSPcomplex::vrts_orient_wrtPlane(const vector<uint32_t>& vrts_inds,
                uint32_t plane_pt0, uint32_t plane_pt1, uint32_t plane_pt2,
                uint32_t count){

   const explicitPoint3D& p0 = vertices[ plane_pt0 ]->toExplicit3D();
   const explicitPoint3D& p1 = vertices[ plane_pt1 ]->toExplicit3D();
   const explicitPoint3D& p2 = vertices[ plane_pt2 ]->toExplicit3D();

   for(uint32_t v=0; v<vrts_inds.size(); v++){
      genericPoint* vrt = vertices[vrts_inds[v]];
      if(isVertexBuiltFromPlane(vrt, &p0, &p1, &p2)) vrts_orBin[vrts_inds[v]]=0;
      else{
        vrts_orBin[vrts_inds[v]] = genericPoint::orient3D(*vrt, p0, p2, p1);

      }
    }
}


//  Input:
// Output:
inline void BSPcomplex::count_vrt_orBin(const vector<uint32_t>& inds,
                                uint32_t* pos, uint32_t* neg, uint32_t* zero){

  (*pos)=0;
  (*neg)=0;
  (*zero)=0;

  for(uint32_t i=0; i<inds.size(); i++)
    if(vrts_orBin[inds[i]]==0)       (*zero)++;
    else if(vrts_orBin[inds[i]]==1)  (*pos)++;
    else if(vrts_orBin[inds[i]]==-1) (*neg)++;
    #ifdef DEBUG_BSP
    else printf("[BSP.cpp]BSPcomplex::count_vrt_orBin: ERROR "
                "wrong access to vrt_orBin[%u]\n", inds[i]);
    #endif
}

//  Input: a BSPedge: edge,
//         vector of the inices of the vertices of the BSPcell to which the
//         BSPedge belong to: cell_vrts.
// Output: true if the constraint intersects the edge interior,
//         false otherwise.
inline bool BSPcomplex::constraint_innerIntersects_edge(const BSPedge& e,
                                            const vector<uint32_t>& cell_vrts){
  return OPPOSITE_SIGNE(vrts_orBin[ e.vertices[0] ], vrts_orBin[ e.vertices[1] ]);
}

//  Input: vector of the inices of the vertices of the BSPface: face_vrts.
// Output: true if the constraint intersects the face interior,
//         false otherwise.
inline bool BSPcomplex::constraint_innerIntersects_face(
                                             const vector<uint32_t>& face_vrts){

  // Face vertices disposition w.r.t. constraint-plane.
  uint32_t vrtsOVER, vrtsUNDER, vrtsON;
  count_vrt_orBin(face_vrts, &vrtsOVER, &vrtsUNDER, &vrtsON);

  // Check if faces[face_ind] has to be splitted
  // (at least 2 vertices with non-zero opposite cell_vrts_orient)
  return (vrtsUNDER > 0 && vrtsOVER > 0);

}

//
//
bool BSPcomplex::coplanar_constraint_innerIntersects_face(const std::vector<uint64_t>& fedges,
                                                          const uint32_t tri[3], int xyz){
    // ASSUMPTION: No face vertices are in the interior of the constraint

    // Process: intersection must be bound by at least three unaligned points on the boundary of tri
    // These points must be on either two different tri_edges,
    // or one of a vertex and two on the opposite edge, or on three vertices

    int mask = 0;

    const BSPedge& edge0 = edges[fedges.back()];
    const BSPedge& edge1 = edges[fedges[0]];
    uint32_t vid_0 = consecEdges_common_endpt(edge0.vertices[0], edge0.vertices[1],
                                              edge1.vertices[0], edge1.vertices[1]);

    // 1) Cerca i tre vertici sul bordo della faccia
    // Per ogni tri_vertex
    //  cerca vertici coincidenti in faccia: se trovi attiva campo in mask e passa al vertice successivo
    //  se campo in mask non è gia attivo, cerca edge della faccia che lo contengono:
    //     se trovi attiva campo in mask e passa al vertice successivo
    // Se mask == 7 torna true
    for (int i = 0; i < 3; i++)
    {
        uint32_t vid = vid_0;
        for (uint64_t e = 0; e < fedges.size(); e++) {
            const BSPedge& edge = edges[fedges[e]];
            if (vid == edge.vertices[0]) vid = edge.vertices[1];
            else vid = edge.vertices[0];
            if (tri[i] == vid) { mask |= (1 << i); break; }
        }
    }
    if (mask == 7) return true;

    for (int i = 0; i < 3; i++) if (!(mask & (1 << i)))
    {
        for (uint64_t e = 0; e < fedges.size(); e++) {
            const BSPedge& edge = edges[fedges[e]];
            const genericPoint* ev1 = vertices[edge.vertices[0]];
            const genericPoint* ev2 = vertices[edge.vertices[1]];
            if (genericPoint::pointInInnerSegment(*vertices[tri[i]], *ev1, *ev2, xyz)) { mask |= (1 << i); break; }
        }
    }
    if (mask == 7) return true;

    // 2) Cerca vertici della faccia su edge del constraint
    //  Per ogni tri_edge
    //   cerca vertici in faccia che stanno innerEdge: se trovi attiva campi in mask e edge successivo
    //   se campi in mask non sono gia attivi, cerca edge della faccia che lo innerCrossano: se trovi attiva campi in mask e passa a edge successivo
    // Se mask == 7 torna true altrimenti false
    for (int i = 0; i < 3; i++)
    {
        uint32_t ti0 = (i + 1) % 3;
        uint32_t ti1 = (i + 2) % 3;
        if ((mask & (1 << ti0)) && (mask & (1 << ti1))) continue;
        uint32_t vid = vid_0;
        for (uint64_t e = 0; e < fedges.size(); e++) {
            const BSPedge& edge = edges[fedges[e]];
            if (vid == edge.vertices[0]) vid = edge.vertices[1];
            else vid = edge.vertices[0];
            const genericPoint* ev1 = vertices[tri[ti0]];
            const genericPoint* ev2 = vertices[tri[ti1]];
            if (genericPoint::pointInInnerSegment(*vertices[vid], *ev1, *ev2, xyz)) { mask |= (1 << ti0); mask |= (1 << ti1); break; }
        }
    }
    if (mask == 7) return true;

    for (int i = 0; i < 3; i++)
    {
        uint32_t ti0 = (i + 1) % 3;
        uint32_t ti1 = (i + 2) % 3;
        if ((mask & (1 << ti0)) && (mask & (1 << ti1))) continue;
        for (uint64_t e = 0; e < fedges.size(); e++) {
            const BSPedge& edge = edges[fedges[e]];
            const genericPoint* ev1 = vertices[edge.vertices[0]];
            const genericPoint* ev2 = vertices[edge.vertices[1]];
            const genericPoint* fv1 = vertices[tri[ti0]];
            const genericPoint* fv2 = vertices[tri[ti1]];
            if (genericPoint::innerSegmentsCross(*ev1, *ev2, *fv1, *fv2, xyz)) return true;
        }
    }

    return false;
}

// Returns 1 if point is on the boundary of the triangle, 2 if it is in the interior, 0 otherwise
int localizedPointInTriangle(const genericPoint& P, const genericPoint& A,
                             const genericPoint& B, const genericPoint& C, int xyz)
{
    int o1, o2, o3;
    if (xyz == 2)
    {
        o1 = genericPoint::orient2Dxy(P, A, B);
        o2 = genericPoint::orient2Dxy(P, B, C);
        o3 = genericPoint::orient2Dxy(P, C, A);
    }
    else if (xyz == 0)
    {
        o1 = genericPoint::orient2Dyz(P, A, B);
        o2 = genericPoint::orient2Dyz(P, B, C);
        o3 = genericPoint::orient2Dyz(P, C, A);
    }
    else
    {
        o1 = genericPoint::orient2Dzx(P, A, B);
        o2 = genericPoint::orient2Dzx(P, B, C);
        o3 = genericPoint::orient2Dzx(P, C, A);
    }
    return ((o1 >= 0 && o2 >= 0 && o3 >= 0) || (o1 <= 0 && o2 <= 0 && o3 <= 0))  +
           ((o1 > 0 && o2 > 0 && o3 > 0) ||  (o1 < 0 && o2 < 0 && o3 < 0)   );
}

//-Upload Delaunay triangolation----------------------

//  Input: pointer to mesh,
//         vector of tetrahedra index type: new_order.
// Output: by using new_order returns new cells indexing (the same as Delaunay
//         tetrahedra indexing, but without ghost-tets).
uint64_t BSPcomplex::removing_ghost_tets(const TetMesh* mesh,
                                         vector<uint64_t>& new_order){
  // new_order have as many element as mesh->tet_num,
  // new_order[i-th tet] =
  //    i - (num of ghost_tet between 0 and i) IF i-th tet is non-ghost
  //    UINT64_MAX                             IF i-th tet is ghost
  uint64_t ghost_tet_count = 0;
  for(uint64_t tet_ind=0; tet_ind<mesh->tet_num; tet_ind++){
      if(IS_GHOST_TET(tet_ind)){
        new_order[tet_ind] = UINT64_MAX;
        ghost_tet_count++;
      }
      else  new_order[tet_ind] = tet_ind - ghost_tet_count;
  }
  return mesh->tet_num - ghost_tet_count;
}

//  Input: pointer to mesh,
//         the 2 endpoints of a tetrahedron edge: e0, e1,
//         the index of the tetrahedron (tet) to which edge belongs: tet_ind,
//         new cells indexing (the same as tetrahedra indexing, but without
//         ghost-tets): new_order.
// Output: the index of the edge (w.r.t. the vector edges) of the
//         edge <endpt0, endpt1>.
// Note. tetrahedra are added in crescent index order.
uint64_t BSPcomplex::add_tetEdge(const TetMesh* mesh, uint32_t e0, uint32_t e1,
                                 uint64_t tet_ind,
                                 const vector<uint64_t>& new_order){
    // Tetrahedra incident in <endpt0,endpt1>
    uint32_t edge_ends[2] = { e0, e1 };
    uint64_t num_incTet;
    uint64_t* incTet = mesh->ETrelation(edge_ends, tet_ind, &num_incTet);
    // Note. ETrelation can return ghost-tet.

    uint64_t min = UINT64_MAX;
    for (uint32_t i = 0; i < num_incTet; i++)
        if (!IS_GHOST_TET(incTet[i]) && incTet[i] < min) min = incTet[i];

    free(incTet);

    if(min == tet_ind){
      // None of the tetrahedra incident in <e0,e1> has been visited,
      // furthermore the edge <e0,e1> do not belongs to a face of tet_ind
      // (as a consequence of the constructor BSPcomplex).
      // A new BSPedge has to be created.
      edges.push_back( BSPedge(e0,e1, e0,e1) );
      return edges.size()-1;
    }

    uint64_t edge_ind;
    for(uint64_t f=0; f<4; f++){        // Note. also with f<3 sholud work.
      BSPface& face = faces[ cells[ new_order[min] ].faces[f] ];
      for(uint64_t e=0; e<3; e++){
        edge_ind = face.edges[e];
        BSPedge& edge = edges[edge_ind];
        if(SAME_EDGE_ENDPTS(e0, e1, edge.vertices[0], edge.vertices[1]) ) break;
      }
      BSPedge& edge = edges[edge_ind];
      if(SAME_EDGE_ENDPTS(e0, e1, edge.vertices[0], edge.vertices[1]) ) break;
    }
    return edge_ind;
}

//  Input: three face (triangle) vertices: v0, v1, v2,
//         index of a BSPcell (w.r.t. vector cells) owning the face: cell_ind,
//         index of a BSPcell (w.r.t. vector cells) faced at the previous one
//         throught the face, if it does not exists (i.e. the previous cell is
//         a complex boundary cell) then it is UINT64_MAX: adjCell_ind.
// Output: returns the index (w.r.t. the vector faces) of the new BSPface.
inline uint64_t BSPcomplex::add_tetFace(uint32_t v0, uint32_t v1, uint32_t v2,
                                      uint64_t cell_ind, uint64_t adjCell_ind){
  // Note. adjCell==UINT64_MAX -> convex-hull face.
  faces.push_back( BSPface(v0,v1,v2, cell_ind,adjCell_ind) );
  uint64_t face_ind = faces.size() -1;
  cells[cell_ind].faces.push_back(face_ind);
  if(!IS_GHOST_CELL(adjCell_ind))  cells[adjCell_ind].faces.push_back(face_ind);

  return face_ind;
}

//  Input: the index of a Delaunay mesh tetrahedron: tet_ind,
//         the index of a Delaunay mesh tetrahedron adjacent to tet: adjTet_ind,
//         the index of the BSPcell corresponding to the adjacent tetrahedron,
//         (in the case it is a ghost-tet it is UINT64_MAX): adjCell_ind.
// Output: returns true if the face between the two input tetrahedra has not
//         been turned into a BSPface yet, false otherwise.
inline bool BSPcomplex::tet_face_isNew(uint64_t tet_ind, uint64_t adjTet_ind,
                                       uint64_t adjCell_ind){
  // The face is the common one between tetrahedra indexed as tet and adjTet.
  // Assuming that the tetrahedron indexed as tet has not been visited yet.
  // Check if:
  // - adjTet has not been visited yet -> new face.
  // - adjTet has been already visited -> face already exists.(not new)
  // - adjTet is ghost (i.e. adjCell is UINT64_MAX by new_order) -> new face.
  return ( adjTet_ind > tet_ind || IS_GHOST_CELL(adjCell_ind) );
}

//
//
inline void BSPcomplex::fill_face_colour(uint64_t tet_ind, uint64_t face_ind,
                                         const uint32_t** map_fi,
                                         const uint32_t* num_map_fi){

  if(num_map_fi[tet_ind] == 0) faces[face_ind].colour = WHITE;
  else{
    // Count non-virtual constraints.
    uint32_t n = 0;
    for(uint32_t cc=0; cc<num_map_fi[tet_ind]; cc++)
      if(!is_virtual(map_fi[tet_ind][cc])) n++;

    if(n == 0) faces[face_ind].colour = WHITE;
    else{
      faces[face_ind].colour = GREY;
      faces[face_ind].coplanar_constraints.resize(n);
      uint32_t pos = 0;
      for(uint32_t cc=0; cc<num_map_fi[tet_ind]; cc++)
        if(!is_virtual(map_fi[tet_ind][cc]))
          faces[face_ind].coplanar_constraints[pos++] = map_fi[tet_ind][cc];
    }

  }
}


//bool faceHasCorrectOrientation(BSPcomplex* cpx, uint64_t f_id)
//{
//    const BSPface& f = cpx->faces[f_id];
//    const uint64_t c_id = f.conn_cells[0];
//    BSPcell& c = cpx->cells[c_id];
//    const uint64_t e0_id = f.edges[0];
//    const uint64_t e1_id = f.edges[1];
//    const uint64_t e2_id = f.edges[2];
//    const BSPedge& e0 = cpx->edges[e0_id];
//    const BSPedge& e1 = cpx->edges[e1_id];
//    const BSPedge& e2 = cpx->edges[e2_id];
//    const uint32_t v0_id = cpx->getFaceVertex(f, 0);
//    const uint32_t v1_id = cpx->getFaceVertex(f, 1);
//    genericPoint* v0 = cpx->vertices[v0_id];
//    genericPoint* v1 = cpx->vertices[v1_id];
//    genericPoint* v2;
//
//    size_t i;
//    for (i = 2; i < f.edges.size(); i++)
//    {
//        const uint32_t v2_id = cpx->getFaceVertex(f, i);
//        v2 = cpx->vertices[v2_id];
//        if (genericPoint::misaligned(*v0, *v1, *v2)) break;
//    }
//    if (i == f.edges.size()) ip_error("Degenerate face\n");
//
//    uint64_t num_cellEdges = UINT64_MAX;
//    uint32_t num_cellVrts = cpx->count_cellVertices(c, &num_cellEdges);
//    vector<uint32_t> cell_vrts(num_cellVrts, UINT32_MAX);
//    cpx->list_cellVertices(c, num_cellEdges, cell_vrts);
//    for (uint32_t vi : cell_vrts) if (!cpx->faceHasVertex(f, vi))
//    {
//        genericPoint* ov = cpx->vertices[vi];
//
//        int ori = genericPoint::orient3D(*ov, *v0, *v1, *v2);
//        if (ori == 0) continue;
//        return (ori < 0);
//    }
//    ip_error("Degenerate cell\n");
//}
//
//// Returns the index of the v_ind'th vertex in f.
//// Vertex 'i' is the common vertex between edge 'i' and edge 'i+1'%num_edges
//uint32_t BSPcomplex::getFaceVertex(const BSPface& f, uint32_t v_ind)
//{
//    const BSPedge& e0 = edges[f.edges[v_ind]];
//    const BSPedge& e1 = edges[f.edges[(v_ind + 1) % (f.edges.size())]];
//    return consecEdges_common_endpt(e0.vertices[0], e0.vertices[1], e1.vertices[0], e1.vertices[1]);
//}
//
//bool BSPcomplex::faceHasVertex(const BSPface& f, uint32_t v_ind)
//{
//    for (uint64_t eid : f.edges)
//    {
//        const BSPedge& e = edges[eid];
//        if (e.vertices[0] == v_ind || e.vertices[1] == v_ind) return true;
//    }
//    return false;
//}

//
// Fills the data scruture with the information of the Delauany mesh.
BSPcomplex::BSPcomplex(const TetMesh* mesh, const constraints_t* _constraints,
                       const uint32_t** map, const uint32_t* num_map,
                       const uint32_t** map_f0, const uint32_t* num_map_f0,
                       const uint32_t** map_f1, const uint32_t* num_map_f1,
                       const uint32_t** map_f2, const uint32_t* num_map_f2,
                       const uint32_t** map_f3, const uint32_t* num_map_f3 ){

  // Uploading the vertices of the mesh
  vertices.resize(mesh->num_vertices);
  for(uint32_t vrt=0; vrt<mesh->num_vertices; vrt++)
    vertices[vrt] = new explicitPoint3D(mesh->vertices[vrt].coord[0],
                                            mesh->vertices[vrt].coord[1],
                                            mesh->vertices[vrt].coord[2]);

  // Initialize vrts_orBin:
  // since orient3D can be -1, 0 or 1 all elements are set to 2.
  vrts_orBin.resize(mesh->num_vertices, 2);

  // Uploading the constraints (the last num_virtual_triangles constraints are virtual.)
  first_virtual_constraint = _constraints->num_triangles - _constraints->num_virtual_triangles;
  constraints_vrts.resize(3*_constraints->num_triangles);
  constraint_group.resize(_constraints->num_triangles);
  for(uint32_t i=0; i < _constraints->num_triangles; i++){
    constraints_vrts[3*i   ] = _constraints->tri_vertices[3*i   ];
    constraints_vrts[3*i +1] = _constraints->tri_vertices[3*i +1];
    constraints_vrts[3*i +2] = _constraints->tri_vertices[3*i +2];
    constraint_group[i] = _constraints->constr_group[i];
  }

  // Establish new tetrahedtra-(cell) indexing: only non-ghost cell are indexed.
  vector<uint64_t> new_order(mesh->tet_num, UINT64_MAX);
  uint64_t cell_num = removing_ghost_tets(mesh, new_order);

  // Creating as many empty cells as the number of non-ghost tet_
  cells.resize(cell_num);
  edges.reserve(cell_num + mesh->num_vertices);
  faces.reserve(cell_num * 2);


  // Loading the cells, creating faces and eadges of the BSP:
  // cells -> the non-ghost tetrahedra in the mesh,
  // faces -> the faces of the non-ghost tetrahedra in the mesh,
  // edges -> the edges of the non-ghost tetrahedra in the mesh.
  for(uint64_t tet_ind=0; tet_ind<mesh->tet_num; tet_ind++){
    // Here each BSPcell is a non-ghost tetrahedron of the mesh:
    // consider a tetrahedron (tet) whose index is tet_ind.
    uint64_t cell_ind = new_order[tet_ind];
    if( IS_GHOST_CELL(cell_ind) )   continue; // Skip ghost-tet.

    // Create BSPcells from tetrahedra by following increasing indexing:
    // all non-ghost tetrahedra which have index lower than tet_ind
    // have been already turned into BSP cells.

    // Constraints improperly intersecated by tet.
    if(num_map[tet_ind]>0){
      cells[cell_ind].constraints.resize( num_map[tet_ind] );
      for(uint32_t i=0; i<num_map[tet_ind]; i++)
        cells[cell_ind].constraints[i] = map[tet_ind][i];
    }

    // Adding BSPface and BSPedges to create a BSPcell conformed to tetrahedron.
    uint32_t v[4]; // Indices of tet vertices.
    v[0] = mesh->tet_node[4 * tet_ind];
    v[1] = mesh->tet_node[4 * tet_ind + 1];
    v[2] = mesh->tet_node[4 * tet_ind + 2];
    v[3] = mesh->tet_node[4 * tet_ind + 3];

    uint64_t face_ind, adjCell_ind, adjTet_ind;
    uint64_t tet_edge[6];
    // --- face <v0,v1,v2> -----------------------------
    adjTet_ind = mesh->tet_neigh[4 * tet_ind + 3] >> 2;
    adjCell_ind = new_order[adjTet_ind];
    if (tet_face_isNew(tet_ind, adjTet_ind, adjCell_ind)) {
        face_ind = add_tetFace(v[0], v[1], v[2], cell_ind, adjCell_ind);
        // At most three edges may have to be created <v0,v1>, <v1,v2>, <v2,v0>.
        tet_edge[0] = add_tetEdge(mesh, v[0], v[1], tet_ind, new_order);
        assigne_edge_to_face(tet_edge[0], face_ind);
        tet_edge[2] = add_tetEdge(mesh, v[2], v[0], tet_ind, new_order);
        assigne_edge_to_face(tet_edge[2], face_ind);
        tet_edge[1] = add_tetEdge(mesh, v[1], v[2], tet_ind, new_order);
        assigne_edge_to_face(tet_edge[1], face_ind);
        // Color and coplanar-constraints
        fill_face_colour(tet_ind, face_ind, map_f3, num_map_f3);
    }
    else {
        face_ind = faceSharedWithCell(cell_ind, adjCell_ind);
        BSPface& face = faces[face_ind];
        tet_edge[0] = find_face_edge(face, v[0], v[1]);
        tet_edge[1] = find_face_edge(face, v[1], v[2]);
        tet_edge[2] = find_face_edge(face, v[2], v[0]);
    }
    // --- face <v3,v0,v1> -----------------------------
    adjTet_ind = mesh->tet_neigh[4 * tet_ind + 2] >> 2;
    adjCell_ind = new_order[adjTet_ind];
    if (tet_face_isNew(tet_ind, adjTet_ind, adjCell_ind)) {
        face_ind = add_tetFace(v[3], v[0], v[1], cell_ind, adjCell_ind);
        // At most two edges may have to be created <v3,v0>, <v1,v3>.
        tet_edge[4] = add_tetEdge(mesh, v[1], v[3], tet_ind, new_order);
        assigne_edge_to_face(tet_edge[4], face_ind);
        tet_edge[3] = add_tetEdge(mesh, v[3], v[0], tet_ind, new_order);
        assigne_edge_to_face(tet_edge[3], face_ind);
        // <v0,v1> is tet_edge[0].
        assigne_edge_to_face(tet_edge[0], face_ind);
        // Color and coplanar-constraints
        fill_face_colour(tet_ind, face_ind, map_f2, num_map_f2);
    }
    else {
        face_ind = faceSharedWithCell(cell_ind, adjCell_ind);
        BSPface& face = faces[face_ind];
        tet_edge[3] = find_face_edge(face, v[3], v[0]);
        tet_edge[4] = find_face_edge(face, v[1], v[3]);
    }
    // --- face <v2,v3,v0> -----------------------------
    adjTet_ind = mesh->tet_neigh[4 * tet_ind + 1] >> 2;
    adjCell_ind = new_order[adjTet_ind];
    if (tet_face_isNew(tet_ind, adjTet_ind, adjCell_ind)) {
        face_ind = add_tetFace(v[2], v[3], v[0], cell_ind, adjCell_ind);
        // At most one edges may have to be created <v2,v3>.
        tet_edge[5] = add_tetEdge(mesh, v[2], v[3], tet_ind, new_order);
        assigne_edge_to_face(tet_edge[5], face_ind);
        // <v3,v0> is tet_edge[3], <v0,v2> is tet_edge[2].
        assigne_edge_to_face(tet_edge[2], face_ind);
        assigne_edge_to_face(tet_edge[3], face_ind);
        // Color and coplanar-constraints
        fill_face_colour(tet_ind, face_ind, map_f1, num_map_f1);
    }
    else {
        face_ind = faceSharedWithCell(cell_ind, adjCell_ind);
        tet_edge[5] = find_face_edge(faces[face_ind], v[2], v[3]);
    }
    // --- face <v1,v2,v3> -----------------------------
    adjTet_ind = mesh->tet_neigh[4 * tet_ind] >> 2;
    adjCell_ind = new_order[adjTet_ind];
    if (tet_face_isNew(tet_ind, adjTet_ind, adjCell_ind)) {
        face_ind = add_tetFace(v[1], v[2], v[3], cell_ind, adjCell_ind);
        // No edges have to be created.
        // <v1,v2> is tet_edge[1], <v2,v3> is tet_edge[5], <v3,v1> is tet_edge[4].
        assigne_edge_to_face(tet_edge[1], face_ind);
        assigne_edge_to_face(tet_edge[5], face_ind);
        assigne_edge_to_face(tet_edge[4], face_ind);
        // Color and coplanar-constraints
        fill_face_colour(tet_ind, face_ind, map_f0, num_map_f0);
    }
  }

  // Initialize visit-flag vectors: all the values are set to upper-limit.
  vrts_visit.resize(vertices.size(), 0);
  edge_visit.resize(edges.size(), 0);


  // Verify that conn_cell[0] is below the face for every face
  //for (size_t fid = 0; fid < faces.size(); fid++)
  //    if (!faceHasCorrectOrientation(this, fid))
  //        ip_error("Wrong orientation\n");
}

//-BSP subdivision----------------------

//  Input: index of a BSPedge (w.r.t. vector face.edges): edge_face_ind,
//         index of two BSPfaces (w.r.t. vector faces): face_ind, newFace_ind.
// Output: nothing.
// Note. edge is removed from faces[face_ind] and assigned to faces[newFace_ind].
inline void BSPcomplex::move_edge(uint64_t edge_face_ind, uint64_t face_ind,
                                  uint64_t newFace_ind){
   uint64_t edge_ind = faces[face_ind].edges[edge_face_ind];
   edges[edge_ind].conn_face_0 = newFace_ind;
   faces[newFace_ind].edges.push_back(edge_ind);
   faces[face_ind].removeEdge(edge_face_ind);
}

//  Input: index of a BSPface (w.r.t. vector cells[cell].faces): face_cell_ind,
//         index of two BSPcells (w.r.t. vector cells): cell_ind, newCell_ind.
// Output: nothing.
// Note. face is removed from cells[cell_ind] and assigned to cells[newCell_ind].
inline void BSPcomplex::move_face(uint64_t face_cell_ind, uint64_t cell_ind,
                                  uint64_t newCell_ind){
   uint64_t face_ind = cells[cell_ind].faces[face_cell_ind];
   faces[face_ind].exchange_conn_cell(cell_ind, newCell_ind);
   cells[newCell_ind].faces.push_back(face_ind);
   cells[cell_ind].removeFace(face_cell_ind);
}

//  Input: index of a constraint (w.r.t. vector cells[cell].constraints):
//         constr_cell_ind,
//         index of a BSPcells (w.r.t. vector cells): cell_ind.
// Output: nothing.
inline void BSPcomplex::remove_constraint(uint32_t constr_cell_ind,
                                          uint64_t cell_ind){
  BSPcell& cell = cells[cell_ind];
  size_t last_ind = cell.constraints.size()-1;
  if(constr_cell_ind != last_ind)
    cell.constraints[constr_cell_ind] = cell.constraints[last_ind];

  cell.constraints.pop_back();
}

// This function is used to identify the edges that have been marked as to
// remove (with UINT64_MAX) from a certain faces[face].edge vector
inline bool remove_edge_from_face(uint64_t edge_face_ind){
  return edge_face_ind == UINT64_MAX;
}

//  Input: index of the splitted BSPface (w.r.t. vector faces) that is going
//         to be turned into the downSubface: face_ind,
//         index of the upSubface (w.r.t. vector faces): newFace_ind.
// Output: nothing.
void BSPcomplex::edgesPartition(uint64_t face_ind, uint64_t newFace_ind){
  // Partition of the edges between upSubface and downSubface.
  BSPface& face = faces[face_ind];

  // Search the first edge having first vertex ==0 and second vertex <0
  // (in face order). This will be the first in upFace.
  uint64_t e;
  for (e = 0; e < face.edges.size(); e++) {
      const uint64_t edge_ind = face.edges[e];
      const uint64_t next_edge_ind = face.edges[((e + 1) == face.edges.size()) ? (0) : (e + 1)];
      const BSPedge& edge = edges[edge_ind];
      const BSPedge& nedge = edges[next_edge_ind];
      const uint32_t comm_vert = (edge.vertices[0] == nedge.vertices[0] || edge.vertices[0] == nedge.vertices[1]) ? 0 : 1;
      if (vrts_orBin[edge.vertices[comm_vert]] < 0 && vrts_orBin[edge.vertices[!comm_vert]] == 0) break;
  }
  if (e == face.edges.size()) ip_error("pippo1\n");

  std::rotate(face.edges.begin(), face.edges.begin() + e, face.edges.end());

  // Search the last edge belonging to upFace
  for (e = 1; e < face.edges.size(); e++) {
      const uint64_t edge_ind = face.edges[e];
      const BSPedge& edge = edges[edge_ind];
      if (vrts_orBin[edge.vertices[0]] == 0 || vrts_orBin[edge.vertices[1]] == 0) break;
  }
  if (e == face.edges.size()) ip_error("pippo2\n");
  e++;
  // Move tail edges to new face
  faces[newFace_ind].edges.assign(face.edges.begin() + e, face.edges.end());
  face.edges.resize(e);
  for (uint64_t ei : faces[newFace_ind].edges) edges[ei].conn_face_0 = newFace_ind;

  //std::move(face.edges.begin() + e, face.edges.end(), faces[newFace_ind].edges.begin());

  //for(uint64_t e=0; e<face.edges.size(); e++){
  //  const uint64_t edge_ind = face.edges[e];
  //  BSPedge& edge = edges[edge_ind];

  //  if(vrts_orBin[ edge.vertices[0] ]>0 || vrts_orBin[ edge.vertices[1] ] >0 )
  //      move_edge(e, face_ind, newFace_ind);
  //}

  //// Remove all edges of face.edges marked as UINT64_MAX
  //face.edges.erase(
  //  std::remove_if(face.edges.begin(), face.edges.end(), remove_edge_from_face),
  //  face.edges.end() );
}

//  Input: index of the splitted BSPcell (w.r.t. vector cells) that is going
//         to be turned into the downSubcell: cell_ind,
//         index of the upSubcell (w.r.t. vector cells): newCell_ind,
//         vector of the indices of the vertices (w.r.t. vector vertices)
//         of the BSPcell to which the BSPface belong to: cell_vrts.
// Output: nothing.
void BSPcomplex::facesPartition(uint64_t cell_ind, uint64_t newCell_ind,
                                const vector<uint32_t>& cell_vrts){

  // Faces whose indices (w.r.t. vector faces) are listed in
  // cells[cell_ind].faces have to be partitioned between
  // upSubcell (i.e. cells[newCell]) and downSubcell (i.e. cells[cell]).
  BSPcell& cell = cells[cell_ind];
  uint64_t num_faces = cell.faces.size();
  uint64_t face_ind;
  for(uint64_t f = 0; f<num_faces; f++){
     face_ind = cell.faces[f];
     BSPface& face = faces[face_ind];

     vector<uint32_t> face_vrts(face.edges.size(), UINT32_MAX);
     list_faceVertices(face, face_vrts);

     // Face vertices disposition w.r.t. constraint-plane.
     uint32_t vrtsOVER, vrtsUNDER, vrtsON;
     count_vrt_orBin(face_vrts, &vrtsOVER, &vrtsUNDER, &vrtsON);

     #ifdef DEBUG_BSP_DEEP
     print_BSPface_edges(edges, face_ind, face.edges);
     print_vrt_orBin(vrts_orBin, face_vrts);
     #endif

     // It is impossible that a face has:
     // - all vertices ON the constraint-plane,
     // - two (or more) vertices on opposite sides w.r.t. the constraint-plane.

     #ifdef DEBUG_BSP
     if(vrtsOVER==0 && vrtsUNDER==0)
        printf("\n[BSP.cpp]BSPcomplex::facesPartition: ERROR face have "
               "all vertices on the constraint plane.\n");
     if(vrtsOVER>0 && vrtsUNDER>0)
        printf("\n[BSP.cpp]BSPcomplex::facesPartition: ERROR face intersects "
               "the constraint plane.\n");
     #endif

     // IF one of the face vertices is OVER the constraint-plane (vrtsOVER>0),
     // the face is assigned to the upSubcell (i.e. cells[newCell]).
     if(vrtsOVER > 0){
       move_face(f, cell_ind, newCell_ind);
       f--;
       num_faces--;
     }
     // OTHERWISE
     // faces[face_ind] goes to down-subcell (i.e. remain to cells[cell_ind])
     // indeed, all its vertices have non-positive cell_vrts_orient.

     #ifdef DEBUG_BSP_DEEP
     if(vrtsOVER > 0)
      printf("\n\tface #%llu goes to up-subcell (cell #%llu).\n",
              face_ind, newCell_ind);
     else
      printf("\n\tface #%llu goes to down-subcell (cell #%llu).\n",
              face_ind, cell_ind);
     #endif
 }
}

//  Input: index of the constraint that has divided the original BSPcell
//         cells[down_cell] in to the current sub-cells cells[up_cell] and
//         the sub-cells cells[down_cell]: ref_constr,
//         index of the down sub-cell (w.r.t. vector cells): down_cell_ind,
//         index of the up sub-cell (w.r.t. vector cells): up_cell_ind,
//         vector of the indices of the vertices (w.r.t. vector vertices)
//         of the down sub-cell (before splitting): cell_vrts,
// Output: nothing.
void BSPcomplex::constraintsPartition(uint32_t ref_constr,
                                      uint64_t down_cell_ind,
                                      uint64_t up_cell_ind,
                                      const vector<uint32_t>& cell_vrts){

  // At this point down sub-cell has all the constraints,
  // while up sub-cell has none.
  BSPcell& down_cell = cells[down_cell_ind];
  BSPcell& up_cell = cells[up_cell_ind];

  #ifdef DEBUG_BSP_DEEP
  printf("\tCurrent constraint is: ");
  print_constraint(constraints_vrts, ref_constr);
  vector<uint32_t> vrts_to_print;
  for(uint32_t v=0; v<cell_vrts.size(); v++)
    if(vrts_orBin[cell_vrts[v]] >= 0) vrts_to_print.push_back(cell_vrts[v]);
  print_BSPcell_vrts(vrts_to_print, up_cell_ind);
  vrts_to_print.clear();
  for(uint32_t v=0; v<cell_vrts.size(); v++)
    if(vrts_orBin[cell_vrts[v]] <= 0) vrts_to_print.push_back(cell_vrts[v]);
  print_BSPcell_vrts(vrts_to_print, down_cell_ind);
  printf("\tConstraints intersecting (original)cell #%llu are:\n", down_cell_ind);
  for(uint32_t c=0; c<down_cell.constraints.size(); c++){
      printf("\t");
      print_constraint(constraints_vrts, down_cell.constraints[c]);
  }
  printf("\n");
  #endif

  // 3 vertices of the ref_constraint seen as triangle (k0,k1,k2).
  uint32_t ref_constr_ID = 3*ref_constr;
  uint32_t k0 = constraints_vrts[ref_constr_ID   ];
  uint32_t k1 = constraints_vrts[ref_constr_ID +1];
  uint32_t k2 = constraints_vrts[ref_constr_ID +2];

  uint64_t num_constr = down_cell.constraints.size();
  uint32_t constr, constr_ID;
  vector<uint32_t> constr_vrts(3, UINT32_MAX);
  for(uint32_t c=0; c<num_constr; c++){
      constr = down_cell.constraints[c];
      constr_ID = 3*constr;
      constr_vrts[0] = constraints_vrts[constr_ID   ];
      constr_vrts[1] = constraints_vrts[constr_ID +1];
      constr_vrts[2] = constraints_vrts[constr_ID +2];

      // commFace_vrts disposition w.r.t. constr vertices.
      vrts_orient_wrtPlane(constr_vrts, k0, k1, k2, 2);
      uint32_t vrtsOVER, vrtsUNDER, vrtsON;
      count_vrt_orBin(constr_vrts, &vrtsOVER, &vrtsUNDER, &vrtsON);

      #ifdef DEBUG_BSP_DEEP
      printf("\n\t orient3d of constraint %u vertices w.r.t. plane for ",
             down_cell.constraints[c]);
      print_constraint(constraints_vrts, ref_constr);
      printf("\n");
      print_vrt_orBin(vrts_orBin, constr_vrts);
      printf("\n");
      #endif

      // If constr and commFace_vrts define the same plane, remove constr since
      // the cut will not produce further cell-split.
      if(vrtsOVER == 0 && vrtsUNDER == 0){
        remove_constraint(c, down_cell_ind);
        c--;
        num_constr--;

        #ifdef DEBUG_BSP_DEEP
        printf("\tConsraints %u and %u define the same plane, constraint %u have "
               "been removed.\n", ref_constr, constr, constr);
        if(vrtsON==0)
          printf("\n[BSP.cpp]BSPcomplex::constraintsPartition: ERROR vertices "
                 "of consraints %u have an undefined position wrt constraint %u.\n",
                 constr, ref_constr);
        #endif

        continue;  // jump to next constraint.
      }

      const bool up = (vrtsOVER > 0);
      const bool down = (vrtsUNDER > 0);

      #ifdef DEBUG_BSP
      if(!up && !down)
        printf("[BSP.cpp]BSPcomplex::constraintsPartition: WARNING constraint "
               "#%u will be removed from down sub-cell #%llu, there are no "
               "intersection with both sub-cells.\n", constr, down_cell_ind );
      #endif

      if(up)  up_cell.constraints.push_back(constr);
      if(!down){
        remove_constraint(c, down_cell_ind);
        c--;
        num_constr--;
      }
      // OTHERWISE
      // the constraint indexed as constr goes to down-subcell (i.e. remain
      // to cells[down_cell]).

      #ifdef DEBUG_BSP_DEEP
      if(up && down)  printf("\tconstraint #%u goes to both subcell (cell "
                             "#%llu and #%llu).\n", constr, up_cell_ind, down_cell_ind);

      if(up && !down) printf("\tconstraint #%u goes to up-subcell (cell "
                             "#%llu).\n", constr, up_cell_ind);
      if(!up && down) printf("\tconstraint #%u goes to down-subcell (cell "
                             "#%llu).\n", constr, down_cell_ind);
      #endif
  }
}

//
//
void BSPcomplex::add_edgeToOrdFaceEdges(BSPface& face, uint64_t newEdge_ind){

  BSPedge& newEdge = edges[newEdge_ind];
  uint64_t edge_ind, num_faceEdges = face.edges.size();
  uint32_t n0, n1, e0, e1;
  n0 = newEdge.vertices[0];
  n1 = newEdge.vertices[1];

  for(uint64_t e=0; e<num_faceEdges; e++){
    edge_ind = face.edges[e];
    e0 = edges[edge_ind].vertices[0];
    e1 = edges[edge_ind].vertices[1];

    if(CONSECUTIVE_EDGES(n0, n1, e0, e1) ){

      // Special case: edge is the first vector element
      if(e==0){
        BSPedge& cons = edges[ face.edges.back() ];
        if(CONSECUTIVE_EDGES(n0, n1, cons.vertices[0], cons.vertices[1]) )
          face.edges.push_back(newEdge_ind);
        else{ // cons is faces[face].edges[1]
          face.edges.push_back(edge_ind);
          face.edges[0] = newEdge_ind;
        }

        return;
      }

      // Special case: edge is the last vector element
      if(e == num_faceEdges-1){
        BSPedge& cons = edges[ face.edges[0] ];
        if(CONSECUTIVE_EDGES(n0, n1, cons.vertices[0], cons.vertices[1]) )
          face.edges.push_back(newEdge_ind);
        else{ // cons is faces[face].edges.size() -2
          face.edges.push_back(edge_ind);
          face.edges[num_faceEdges-1] = newEdge_ind;
        }

        return;
      }

      // General case
      BSPedge& cons = edges[ face.edges[e+1] ];
      if(CONSECUTIVE_EDGES(n0, n1, cons.vertices[0], cons.vertices[1]) )
        face.edges.insert(face.edges.begin() + e+1, newEdge_ind);
      else  // cons = edges[ faces[face].edges[e-1] ]
        face.edges.insert(face.edges.begin() + e, newEdge_ind);

      return;
    }
  }
}

//  Input: index of the constraint splitting the BSPface: constr,
//         index of the splitted BSPface (w.r.t. vector faces) that is going
//         to be turned into the downSubface: face_ind,
//         index of the upSubface (w.r.t. vector cells): newFace_ind,
//         vector of the indices of the 2 vertices (w.r.t. vector vertices)
//         of the original BSPface faces[face] that lie
//         on the constraint-plane: endpts.
// Output: nothing.
void BSPcomplex::add_commonEdge(uint32_t constr, uint64_t face_ind,
                                uint64_t newFace_ind, const uint32_t* endpts){
  BSPface& face = faces[face_ind];
  BSPface& newFace = faces[newFace_ind];
  // The new faces "face" and "newFace", originated by splitting the original
  // "face" with the constraint-plane (constr), have a common edge.
  // The endpoint of the common edge are those that have vrt_orBin=0, i.e.
  // endpts[0], endpts[1].

  // New edge originated by the face-constraint intersection:
  // it is defined by the intersection of two planes (face and constraint)
  uint32_t constr_ID = 3*constr;
  edges.push_back( BSPedge(endpts[0], endpts[1],
                           face.meshVertices[0],
                           face.meshVertices[1],
                           face.meshVertices[2],
                           constraints_vrts[constr_ID   ],
                           constraints_vrts[constr_ID +1],
                           constraints_vrts[constr_ID +2] ) );
  uint64_t commEdge_ind = edges.size() -1;
  BSPedge& commEdge = edges[commEdge_ind];
  //commEdge.conn_faces.push_back(newFace_ind);
  //commEdge.conn_faces.push_back(face_ind);
  commEdge.conn_face_0 = face_ind;

  // Add an element to global vector edge_visit.
  edge_visit.push_back(0);

  // Add the common edge to face and newFace.
  //add_edgeToOrdFaceEdges(face, commEdge_ind);
  //add_edgeToOrdFaceEdges(newFace, commEdge_ind);
  face.edges.push_back(commEdge_ind);
  newFace.edges.push_back(commEdge_ind);
}

//  Input: an empty BSPface, created by intersecating a BSPcell with a
//         constraint: face,
//         the indices of all BSPedges (w.r.t. vector edges) that bound the
//         BSPface: edges_ind.
// Output: nothing.
void BSPcomplex::add_edges_toCommFaceEdges(BSPface& face,
                                       const vector<uint64_t>& edges_ind){
  // Find face vertices: since the face boundary is closed,
  //                     there are as many vertices as are the edges.
  vector<uint32_t> face_vrts(edges_ind.size(), UINT32_MAX);
  uint64_t edge_ind;
  uint32_t e0, e1, fv_ind=0;
  for(uint64_t e=0; e<edges_ind.size(); e++){
    BSPedge& edge = edges[ edges_ind[e] ];
    e0 = edge.vertices[0];
    if(vrts_visit[e0]==0){
      face_vrts[fv_ind++] = e0;
      vrts_visit[e0]=1;
    }
    e1 = edge.vertices[1];
    if(vrts_visit[e1]==0){
      face_vrts[fv_ind++] = e1;
      vrts_visit[e1]=1;
    }
  }

  //Relate each face vertex with its two incident edges.
  vector<uint64_t> rel_VE(2*face_vrts.size(), UINT64_MAX);

  // Set vrts_visit of face_vrts indices in order to memory positions.
  for(uint32_t u=0; u<face_vrts.size(); u++)
    vrts_visit[ face_vrts[u] ]=UINT32_MAX;

  // Fill rel_VE
  uint32_t pos=0;
  for(uint64_t e=0; e<edges_ind.size(); e++){
    edge_ind = edges_ind[e];
    BSPedge& edge = edges[edge_ind];
    e0 = edge.vertices[0];
    e1 = edge.vertices[1];
    if( vrts_visit[e0] == UINT32_MAX ){
      rel_VE[2*pos]=edge_ind;
      vrts_visit[e0]=pos++;
    }
    else  rel_VE[ 2*vrts_visit[e0] +1 ] = edge_ind;

    if( vrts_visit[e1] == UINT32_MAX ){
      rel_VE[2*pos]=edge_ind;
      vrts_visit[e1]=pos++;
    }
    else  rel_VE[ 2*vrts_visit[e1] +1 ] = edge_ind;

  }

  // Fill faces[face_ind].edges
  uint32_t num_ins_vrts=0, next_vrt=face_vrts[0];
  edge_ind = rel_VE[ 2*vrts_visit[next_vrt] ];

  face.edges.resize(edges_ind.size());
  while( num_ins_vrts < face_vrts.size() ){
    if(edge_visit[edge_ind]==0){
      edge_visit[edge_ind]=1;
      face.edges[num_ins_vrts++] = edge_ind;

      e0 = edges[edge_ind].vertices[0];
      e1 = edges[edge_ind].vertices[1];
      if( next_vrt == e0 ) next_vrt = e1;
      else next_vrt = e0;
      edge_ind = rel_VE[ 2*vrts_visit[next_vrt] ];
    }
    else{
      if(edge_ind == rel_VE[ 2*vrts_visit[next_vrt] ])
        edge_ind = rel_VE[ 2*vrts_visit[next_vrt] +1 ];
      else edge_ind = rel_VE[ 2*vrts_visit[next_vrt] ];
    }
  }

  // Reset vrts_visit and edge_visit
  for(uint32_t u=0; u<face_vrts.size(); u++) vrts_visit[ face_vrts[u] ]=0;
  for(uint64_t u=0; u<edges_ind.size(); u++) edge_visit[ edges_ind[u] ]=0;
}

//  Input: index of the constraint splitting the BSPcell: constr,
//         index of the splitted BSPcell (w.r.t. vector cells) that is going
//         to be turned into the downSubcell: cell_ind,
//         index of the upSubcell (w.r.t. vector cells): newCell_ind,
//         vector of the indices of the vertices (w.r.t. vector vertices)
//         of the BSPcell to which the BSPface belong to: cell_vrts.
// Output: nothing.
void BSPcomplex::add_commonFace(uint32_t constr,
                                uint64_t cell_ind, uint64_t newCell_ind,
                                const vector<uint32_t>& cell_vrts,
                                const vector<uint64_t>& cell_edges){
  // Common face between up-subcell and down-subcell: the edge of that face
  // are those of cells[cell_ind] that have vrts_orBin = 0.
  uint32_t constr_ID = 3*constr;
  COLOUR_T colour = GREY;
  if(is_virtual(constr)) colour = WHITE;
  faces.push_back( BSPface(constraints_vrts[constr_ID],
                           constraints_vrts[constr_ID+1],
                           constraints_vrts[constr_ID+2],
                           cell_ind, newCell_ind, colour )   );
  uint64_t face_ind = faces.size() - 1;

  uint64_t edge_ind, num_commFace_edges=0;
  // Count edges whose endpoints are both on the constraint-plane.
  for(uint64_t e=0; e<cell_edges.size(); e++){
      BSPedge& edge = edges[ cell_edges[e] ];
      if(vrts_orBin[ edge.vertices[0] ]==0 &&
         vrts_orBin[ edge.vertices[1] ]==0    )
          num_commFace_edges++;
  }
  // Fill a vector with those edges.
  vector<uint64_t> commFace_edges(num_commFace_edges, UINT64_MAX);
  num_commFace_edges = 0;
  for(uint64_t e=0; e<cell_edges.size(); e++){
      edge_ind = cell_edges[e];
      BSPedge& edge = edges[edge_ind];
      if(vrts_orBin[ edge.vertices[0] ]==0 &&
         vrts_orBin[ edge.vertices[1] ]==0    )
          commFace_edges[num_commFace_edges++] = edge_ind;
  }

  add_edges_toCommFaceEdges(faces[face_ind], commFace_edges);

  //for (uint64_t e = 0; e < commFace_edges.size(); e++)
  //    edges[commFace_edges[e]].conn_faces.push_back(face_ind);
  for (uint64_t e = 0; e < commFace_edges.size(); e++)
      edges[commFace_edges[e]].conn_face_0 = face_ind;

  cells[cell_ind].faces.push_back(face_ind);
  cells[newCell_ind].faces.push_back(face_ind);

  fixCommonFaceOrientation(face_ind);
}

inline bool edges_ShareCommonPlanes(const BSPedge& a, const BSPedge& b)
{
    const uint32_t* mva = a.meshVertices;
    const uint32_t* mvb = b.meshVertices;
    return (mva[0] == mvb[0] && mva[1] == mvb[1] && mva[2] == mvb[2] && mva[3] == mvb[3] && mva[4] == mvb[4] && mva[5] == mvb[5]);
}

void BSPcomplex::fixCommonFaceOrientation(uint64_t cf_id)
{
    BSPface& f = faces[cf_id];
    const uint32_t* mv = f.meshVertices;
    double mvc[9]; // Coords of the original input triangle
    vertices[mv[0]]->getApproxXYZCoordinates(mvc[0], mvc[1], mvc[2]);
    vertices[mv[1]]->getApproxXYZCoordinates(mvc[3], mvc[4], mvc[5]);
    vertices[mv[2]]->getApproxXYZCoordinates(mvc[6], mvc[7], mvc[8]);
    int xyz = genericPoint::maxComponentInTriangleNormal(mvc[0], mvc[1], mvc[2],
        mvc[3], mvc[4], mvc[5], mvc[6], mvc[7], mvc[8]);

    int ori0;
    if (xyz == 2) ori0 = orient2d(mvc[0], mvc[1], mvc[3], mvc[4], mvc[6], mvc[7]);
    else if (xyz == 0) ori0 = orient2d(mvc[1], mvc[2], mvc[4], mvc[5], mvc[7], mvc[8]);
    else ori0 = orient2d(mvc[2], mvc[0], mvc[5], mvc[3], mvc[8], mvc[6]);

    const BSPedge& edge0 = edges[f.edges.back()];
    const BSPedge& edge1 = edges[f.edges[0]];
    uint32_t vid = consecEdges_common_endpt(edge0.vertices[0], edge0.vertices[1], edge1.vertices[0], edge1.vertices[1]);

    genericPoint* v0 = vertices[vid];
    if (vid == edge1.vertices[0]) vid = edge1.vertices[1];
    else vid = edge1.vertices[0];
    genericPoint* v1 = vertices[vid];
    for (uint64_t e = 1; e < f.edges.size(); e++) {
        const BSPedge& edge = edges[f.edges[e]];
        if (vid == edge.vertices[0]) vid = edge.vertices[1];
        else vid = edge.vertices[0];
        if (edges_ShareCommonPlanes(edge1, edge)) continue;
        genericPoint* v2 = vertices[vid];
        int ori = ori0*genericPoint::orient2D(*v0, *v1, *v2, xyz);
        if (ori < 0) return;
        if (ori > 0) {
            if (f.conn_cells[1] == UINT64_MAX)
            {
                printf("Mmh... this should not happen\n");
                return;
            }
            std::swap(f.conn_cells[0], f.conn_cells[1]);
            return;
        }
    }
    ip_error("Degenerate face\n");
}

//  Input: a BSPedge: edge,
//         index of a constraint (w.r.t. vector constraints_vrt): constr.
// Output: returns the index of the new vertex (w.r.t. vector vertices).
// Note: a LPI (Line Plane Intersection) vertex can be generated only as
//       intesection beween a constraint-plane and an edge that is part of the
//       Delaunay triangulation (or a pice of Delaunay edge).
uint32_t BSPcomplex::add_LPIvrt(const BSPedge& edge, uint32_t constr){

  uint32_t constr_ID = 3*constr;
  genericPoint* c0 = vertices[ constraints_vrts[constr_ID  ] ];
  genericPoint* c1 = vertices[ constraints_vrts[constr_ID+1] ];
  genericPoint* c2 = vertices[ constraints_vrts[constr_ID+2] ];
  genericPoint* e0 = vertices[ edge.meshVertices[0] ];
  genericPoint* e1 = vertices[ edge.meshVertices[1] ];

  #ifdef DEBUG_BSP
  if(!(e0->isExplicit3D()) ||
     !(e1->isExplicit3D()) ||
     !(c0->isExplicit3D()) ||
     !(c1->isExplicit3D()) ||
     !(c2->isExplicit3D())    )
    printf("\n[BSP.cpp]BSPcomplex::add_LPIvrt: ERROR explicitPoint3D are "
           "expected (LPI).\n");
  #endif

  vertices.push_back(
        new implicitPoint3D_LPI(
            e0->toExplicit3D(), e1->toExplicit3D(),
            c0->toExplicit3D(), c1->toExplicit3D(), c2->toExplicit3D() ) );

  // Add new element to global vectors vrts_orBin and vrts_visit.
  vrts_orBin.push_back(2);
  vrts_visit.push_back(0);

  return (uint32_t)(vertices.size() -1);

}

// Return TRUE if the intersection of the sets {p,q,r} and {t,u,v} has at least two elements.
// Fill 'e' with the intersecting elements.
inline bool twoEqualVertices(uint32_t p, uint32_t q, uint32_t r,
                             uint32_t s, uint32_t t, uint32_t u, uint32_t *e){
    int i = 0;
    if (p == s || p == t || p == u) e[i++] = p;
    if (q == s || q == t || q == u) e[i++] = q;
    if (r == s || r == t || r == u) e[i++] = r;
    return (i >= 2);
}

//  Input: a BSPedge: edge,
//         index of a constraint (w.r.t. vector constraints_vrt): constr.
// Output: returns the index of the new vertex (w.r.t. vector vertices).
// Note: a TPI (Triple Plane Intersection) vertex can be generated only as
//       intesection beween a constraint-plane and an edge that is defined as
//       the intersection between 3 constraint-planes (i.e. not included in
//       the original Delaunay triangulation).
uint32_t BSPcomplex::add_TPIvrt(const BSPedge& edge, uint32_t constr){

  const uint32_t constr_ID = 3*constr;

  const uint32_t ic0 = constraints_vrts[constr_ID];
  const uint32_t ic1 = constraints_vrts[constr_ID + 1];
  const uint32_t ic2 = constraints_vrts[constr_ID + 2];
  const uint32_t ie0 = edge.meshVertices[0];
  const uint32_t ie1 = edge.meshVertices[1];
  const uint32_t ie2 = edge.meshVertices[2];
  const uint32_t ie3 = edge.meshVertices[3];
  const uint32_t ie4 = edge.meshVertices[4];
  const uint32_t ie5 = edge.meshVertices[5];

  // If two of the three triangles share two vertices -> create an LPI
  uint32_t comm[3];
  if(twoEqualVertices(ic0, ic1, ic2, ie0, ie1, ie2, comm))
      vertices.push_back(new implicitPoint3D_LPI(
          vertices[comm[0]]->toExplicit3D(),
          vertices[comm[1]]->toExplicit3D(),
          vertices[ie3]->toExplicit3D(),
          vertices[ie4]->toExplicit3D(),
          vertices[ie5]->toExplicit3D()       ) );
  else if(twoEqualVertices(ic0, ic1, ic2, ie3, ie4, ie5, comm))
          vertices.push_back(new implicitPoint3D_LPI(
              vertices[comm[0]]->toExplicit3D(),
              vertices[comm[1]]->toExplicit3D(),
              vertices[ie0]->toExplicit3D(),
              vertices[ie1]->toExplicit3D(),
              vertices[ie2]->toExplicit3D()     ) );
  else if(twoEqualVertices(ie3, ie4, ie5, ie0, ie1, ie2, comm))
          vertices.push_back(new implicitPoint3D_LPI(
              vertices[comm[0]]->toExplicit3D(),
              vertices[comm[1]]->toExplicit3D(),
              vertices[ic0]->toExplicit3D(),
              vertices[ic1]->toExplicit3D(),
              vertices[ic2]->toExplicit3D()         ) );
  else{
    vertices.push_back( new implicitPoint3D_TPI( vertices[ie0]->toExplicit3D(),
                                                 vertices[ie1]->toExplicit3D(),
                                                 vertices[ie2]->toExplicit3D(),
                                                 vertices[ie3]->toExplicit3D(),
                                                 vertices[ie4]->toExplicit3D(),
                                                 vertices[ie5]->toExplicit3D(),
                                                 vertices[ic0]->toExplicit3D(),
                                                 vertices[ic1]->toExplicit3D(),
                                                 vertices[ic2]->toExplicit3D() ));
  }

  // Add new element to global vectors vrts_orBin and vrts_visit.
  vrts_orBin.push_back(2);
  vrts_visit.push_back(0);

  return (uint32_t)(vertices.size() -1);
}



// 2 funzioni
// 1 - 

// Given a cell 'c', one of its faces 'f0', and one of f0 edges 'e0'
// return the face in 'c' that shares 'e0' with 'f0'
uint64_t BSPcomplex::getOppositeEdgeFace(const uint64_t e0, const uint64_t f0, const uint64_t c)
{
    const std::vector<uint64_t>& cfaces = cells[c].faces;
    for (uint64_t fid : cfaces) if (fid != f0)
    {
        const std::vector<uint64_t>& fedges = faces[fid].edges;
        for (uint64_t e : fedges) if (e == e0) return fid;
    }
    return UINT64_MAX;
}

static inline uint64_t oppositeCellId(const uint64_t c_id, const BSPface& f)
{
    if (f.conn_cells[0] == c_id) return f.conn_cells[1];
    else return f.conn_cells[0];
}

void BSPcomplex::makeEFrelation(const uint64_t e_id, std::vector<uint64_t>& ef)
{
    const BSPedge& e = edges[e_id];
    uint64_t f = e.conn_face_0;
    uint64_t c = faces[f].conn_cells[0];
    ef.push_back(f);

    for (;;)
    {
        f = getOppositeEdgeFace(e_id, f, c);
        if (f == e.conn_face_0) return;
        else ef.push_back(f);
        c = oppositeCellId(c, faces[f]);
        if (c == UINT64_MAX) break;
    }

    f = e.conn_face_0;
    if ((c = faces[f].conn_cells[1]) == UINT64_MAX) return;

    for (;;)
    {
        f = getOppositeEdgeFace(e_id, f, c);
        ef.push_back(f);
        c = oppositeCellId(c, faces[f]);
        if (c == UINT64_MAX) return;
    }
}
//  Input: a BSPedge: edge,
//         index of a constraint that intersects the edge interior: constr.
// Output: the index of the sub-edge (w.r.t. vector edges) originated by the
//         intersection between the edge and the constraint, whose endpoints
//         have non-negative orient3D w.r.t. the constraint plane.
void BSPcomplex::splitEdge(uint64_t edge_id, uint32_t constr){

  BSPedge& edge = edges[edge_id];

  std::vector<uint64_t> ef;
  makeEFrelation(edge_id, ef);

  // Add the point in which the edge intersects the constraint-plane
  // to the vector vertices.
  // edge is a Delauany mesh edge -> Line Plane Intersection.
  // edge is a 2 constraints intersection -> Triple Plane Intersection.
  uint32_t new_point;
  if(edge.meshVertices[2]==UINT32_MAX) new_point = add_LPIvrt(edge, constr);
  else                                 new_point = add_TPIvrt(edge, constr);

  // Split the edge <e0,e1> -> <e0,new_point> + <new_point,e1>
  // edge <- <new_point,e1>
  // new_edge = edges[edges.size()-1] <- <e0,new_point>.
  edges.push_back( edge.split(new_point) );
  uint64_t new_edge_ind = edges.size()-1;
  // Note. new edge is created with the same conn_faces of the old edge.

  // Add new element to edge_visit
  edge_visit.push_back(0);

  // Add new edge to all conn_faces of old edge.
  //for(uint64_t f=0; f< edges[edge_id].conn_faces.size(); f++)
  //  add_edgeToOrdFaceEdges(faces[edges[edge_id].conn_faces[f] ], new_edge_ind);

  for (uint64_t f : ef) add_edgeToOrdFaceEdges(faces[f], new_edge_ind);
}

//  Input: index of a BSPface: face_ind,
//         index of a constraint that intersects the face interior: constr,
//         index of the BSPcell to which the BSPface belongs to: cell_ind,
//         a vector with the indices of the face vertices: face_vrts.
// Output: nothing.
void BSPcomplex::splitFace(uint64_t face_ind, uint32_t constr,
                          uint64_t cell_ind, const vector<uint32_t>& face_vrts){

  #ifdef DEBUG_BSP_DEEP
  printf("\n\tDividing face #%llu with constraint %u.\n", face_ind, constr);
  #endif

  BSPface& face = faces[face_ind];

  // The face faces[face] is divided in two subfaces:
  // the up-subface and the down-subface.
  // The edges of faces[face] have to be partitioned between them.
  // up-subface inherits edges whose vertices have non-negative vrt_orBin.
  // down-subface inherits edges whose vertices have non-positive vrt_orBin.

  // Face faces[face], after splitting, becomes the down-subface.
  // New BSPface is the up-subface of the original faces[face].
  faces.push_back( BSPface(face.meshVertices[0],
                           face.meshVertices[1],
                           face.meshVertices[2],
                           face.conn_cells[0],
                           face.conn_cells[1],
                           face.colour,
                           face.coplanar_constraints) );
  uint64_t newFace_ind = faces.size() -1;
  // Note. the edge of the new face (i.e. up-subface) will be assigned later.

  #ifdef DEBUG_BSP_DEEP
  printf("\tSplit face %llu -> up-subface (face #%llu) + down-subface "
         "(face #%llu)\n", face_ind, newFace_ind, face_ind);
  #endif

  // Add the new face to its adjacent cell (the same of faces[face]).
  // If it is a convex-hull face (i.e. conn_cells[1] = UINT64_MAX),
  // there is no adjacent cell.

  cells[faces[face_ind].conn_cells[0] ].faces.push_back(newFace_ind);
  if( !IS_GHOST_CELL(faces[face_ind].conn_cells[1]) )
    cells[faces[face_ind].conn_cells[1] ].faces.push_back(newFace_ind);

  // Find face vertices that have vrts_orBin=0.
  uint32_t zero_vrts[2];
  uint32_t pos=0;
  for(uint32_t v=0; v<face_vrts.size(); v++)
    if(vrts_orBin[ face_vrts[v] ]==0)  zero_vrts[pos++]=face_vrts[v];

  #ifdef DEBUG_BSP
  if(pos!=2)
    printf("\n[BSP.cpp]BSPcomplex::splitFace: ERROR there must be exactly "
           "2 (not %u) face-vertices on the constraint plane.\n", pos);
  #endif

  // Partition of the edges between up-subface and down-subface.
  edgesPartition(face_ind, newFace_ind);

  // The new faces faces[face_ind] and faces[newFace_ind] have a common edge.
  add_commonEdge(constr, face_ind, newFace_ind, zero_vrts);


  //if (!faceHasCorrectOrientation(this, face_ind)) ip_error("Wrong f1\n");
  //if (!faceHasCorrectOrientation(this, newFace_ind)) ip_error("Wrong f2\n");
}

//
//
void BSPcomplex::find_coplanar_constraints(uint64_t cell_ind, uint32_t constr,
                                               vector<uint32_t>& coplanar_c){
  BSPcell& cell = cells[cell_ind];
  uint32_t constr_ID = 3*constr;
  uint32_t c0 = constraints_vrts[constr_ID  ];
  uint32_t c1 = constraints_vrts[constr_ID+1];
  uint32_t c2 = constraints_vrts[constr_ID+2];

  // Count coplanar constraints.
  uint32_t num_coplanar = 0;
  vector<uint32_t> k_vrts(3, UINT32_MAX);
  for(uint32_t k=0; k < cell.constraints.size(); k++){
      if (is_virtual(cell.constraints[k])) continue;
      uint32_t kID = 3*cell.constraints[k];
      k_vrts[0] = constraints_vrts[kID  ];
      k_vrts[1] = constraints_vrts[kID+1];
      k_vrts[2] = constraints_vrts[kID+2];
      vrts_orient_wrtPlane(k_vrts, c0, c1, c2, 0);
      if(vrts_orBin[ k_vrts[0] ]==0 &&
         vrts_orBin[ k_vrts[1] ]==0 &&
         vrts_orBin[ k_vrts[2] ]==0    ){

        #ifdef DEBUG_BSP_DEEP
        printf("cell #%llu: constraints %u and %u are aligned.\n",
                cell_ind, constr_ID/3, kID/3);
        #endif

        num_coplanar++;
      }
  }
  if (!is_virtual(constr)) num_coplanar++;

  // Save coplanar constraints in a new vector and remove them from cell's list.
  if(num_coplanar > 0){
    coplanar_c.resize(num_coplanar, UINT32_MAX);
    uint32_t pos=0;
    for(uint32_t k=0; k < cell.constraints.size(); k++){
        if (is_virtual(cell.constraints[k])) continue;
        uint32_t kID = 3*cell.constraints[k];
        k_vrts[0] = constraints_vrts[kID  ];
        k_vrts[1] = constraints_vrts[kID+1];
        k_vrts[2] = constraints_vrts[kID+2];

        if(vrts_orBin[ k_vrts[0] ]==0 &&
           vrts_orBin[ k_vrts[1] ]==0 &&
           vrts_orBin[ k_vrts[2] ]==0    ){

          coplanar_c[pos++] = cell.constraints[k];
          cell.constraints[k] = cell.constraints.back();
          cell.constraints.pop_back();
          k--;
        }
    }
    if (!is_virtual(constr)) coplanar_c[pos++] = constr;
  }

}

//  Input: index of a BSPcell: cell_ind.
// Output: nothing.
// Note. it is assumed that the BScell is (only) convex,
//       strictly convexity is not guaranteed.
void BSPcomplex::splitCell(uint64_t cell_ind){

  // Extract the last contraint that intersect the cell and remove it from
  // the list. The cell will be splitted by that constraint.
  BSPcell& cell = cells[cell_ind];
  uint32_t constr = cell.constraints.back();
  uint32_t constr_ID = 3*constr;
  uint32_t c0 = constraints_vrts[constr_ID  ];
  uint32_t c1 = constraints_vrts[constr_ID+1];
  uint32_t c2 = constraints_vrts[constr_ID+2];

  cell.constraints.pop_back();

  // Search for coplanar constraints.
  vector<uint32_t> coplanar_constr;
  find_coplanar_constraints(cell_ind, constr, coplanar_constr);

  // Distinguish between two mutually exclusive cases:
  // CASE. NO SPLIT: only cell boundary elements (face or edge) lie on the
  //                 constraint-plane, while the other elements belong to
  //                 the same half-space w.r.t. the constraint-plane.
  // CASE. SPLIT-INTERIOR: constraint-plane pass through the cell interior.
  // Note. it is not possible that only a vertex lies on the constraint-plane.

  // Create a local richer data structure to avoid multilple extracions of cell
  // edges and vertices.
  uint64_t num_cellEdges = count_cellEdges(cell);
  uint64_t num_cellVrts = count_cellVertices(cell, &num_cellEdges);
  vector<uint64_t> cell_edges(num_cellEdges, UINT64_MAX);
  vector<uint32_t> cell_vrts(num_cellVrts, UINT32_MAX);
  fill_cell_locDS(cell, cell_edges, cell_vrts);

  // Compute the orientation of cell vertices w.r.t. the constraint plane.
  vrts_orient_wrtPlane(cell_vrts, c0, c1, c2, 1);

  // Analysis of cell vertices disposition w.r.t. constraint-plane.
  uint32_t vrtsON, vrtsOVER, vrtsUNDER;
  count_vrt_orBin(cell_vrts, &vrtsOVER, &vrtsUNDER, &vrtsON);

  // CASE. NO SPLIT:
  // at least two cell_vrts_or are 0, the other (!=0) have the same signe.
  if(vrtsUNDER==0 || vrtsOVER==0) return;

  // (else) CASE. SPLIT-INTERIOR:
  // at least two cell_vrts_or have opposite signe.

  // Split cell edges whose endpoints have opposite cell_vrts_or signe.

  for(uint64_t e=0; e<cell_edges.size(); e++){
      BSPedge& edge = edges[cell_edges[e]];
      if(constraint_innerIntersects_edge(edge, cell_vrts)){
        splitEdge(cell_edges[e], constr);  // Here the edge is splitted.
        // Add the new vertex to cell_vrts, compute its orient3D w.r.t. the
        // constraint-plane and add it to cell_vrts_or.
        uint32_t new_vrt = (uint32_t)(vertices.size() -1);
        cell_vrts.push_back(new_vrt);
        cell_edges.push_back(edges.size() -1);
        vrts_orBin[new_vrt] = 0;

        #ifdef DEBUG_BSP_DEEP
        vector<uint32_t> vrt_to_print;
        vrt_to_print.push_back(new_vrt);
        printf("\n");
        print_vrt_orBin(vrts_orBin, vrt_to_print);
        #endif
      }
  }
  // Now all the BSPcell edges are over or under the constraint.

  // Split cell-faces that have at lesat two vertices with opposite
  // cell_vrts_or signe.

  uint64_t num_faces = cell.faces.size();
  for(uint64_t f=0; f<num_faces; f++){
    uint64_t face_ind = cell.faces[f];
    BSPface& face = faces[face_ind];
    vector<uint32_t> face_vrts(face.edges.size(), UINT32_MAX);
    list_faceVertices(face, face_vrts);

    if(constraint_innerIntersects_face(face_vrts) ){

      // Here the face is splitted.
      // The intersection between face and constraint-plane is a new edge.
      splitFace(face_ind, constr, cell_ind, face_vrts);

      // Add new edge (the face-slpitting one) to cell_edges
      cell_edges.push_back(edges.size()-1);
    }
  }

  // The cell is divided in two subcells: the up-subcell and the down-subcell.
  // The up-subcell have vertices with non-negative cell_vrts_or.
  // The down-subcell have vertices with non-positive cell_vrts_or.
  //
  // The cell faces are partitioned between the two subcells.
  // A cell face intesecting the constraint-plane is splitted in sub-faces.
  //
  // A new face is created: the common face between up-subcell and down-subcell.
  //
  // Conventionally the down-subcell replace cell in the vector cells,
  // while up-subcell is appended to the vector cells.

  // up-subcell
  cells.push_back( BSPcell() );
  uint64_t newCell_ind = cells.size() - 1;

  facesPartition(cell_ind, newCell_ind, cell_vrts);
  // Add common face between up-subcell and down-subcell.
  add_commonFace(constr, cell_ind, newCell_ind, cell_vrts, cell_edges);

  // If there are coplanar constraints (to the one used for splitting) those
  // constraints have to be aded to the common face.
  uint32_t num_coplanar_constr = (uint32_t)coplanar_constr.size();
  faces.back().coplanar_constraints.resize(1 + num_coplanar_constr);
  faces.back().coplanar_constraints[0] = constr;
  if(num_coplanar_constr > 0)
    for(uint32_t cc=0; cc<num_coplanar_constr; cc++)
      faces.back().coplanar_constraints[1+cc] = coplanar_constr[cc];

  // Constraints that have to be partitioned between up-subcell
  // and down-subcell are: cells[cell].constarints.

  constraintsPartition(constr, cell_ind, newCell_ind, cell_vrts);
}

//--Complex tetrahedralization-------------------------

//  Input: the index of a BSPface w.r.t. vector faces: face_ind.
// Output: nothing.
// Note. The last 2 edges of the vector face[face_ind].edges are used to
//       detach a triangualr face from the face faces[face_ind].
//       This 2 last edges are replaced in the vector face[face_ind].edges by
//       one new edge: the one that closes the detached triangualr face.
void BSPcomplex::triangle_detach(uint64_t face_ind){
  BSPface& face = faces[face_ind];
  uint64_t num_face_edges = face.edges.size();

  // A triangle will be created by using:
  // - edges[face.edges.size()-1] and edges[face.edges.size()-2],
  // - introducing a new_edge to close the triangle.
  uint64_t s_01_ind = face.edges[num_face_edges-1];
//  BSPedge& s_01 = edges[s_01_ind];
  uint64_t s_12_ind = face.edges[num_face_edges-2];
//  BSPedge& s_12 = edges[s_12_ind];
  uint32_t t1 = consecEdges_common_endpt(edges[s_01_ind].vertices[0], edges[s_01_ind].vertices[1],
      edges[s_12_ind].vertices[0], edges[s_12_ind].vertices[1]);
  uint32_t t0 = other_edge_endpt(edges[s_01_ind].vertices[0], edges[s_01_ind].vertices[1], t1);
  uint32_t t2 = other_edge_endpt(edges[s_12_ind].vertices[0], edges[s_12_ind].vertices[1], t1);

  // Connect t2 and t0 with a new edge.
  edges.push_back(BSPedge());
  uint64_t s_20_ind = edges.size() -1;
  BSPedge& s_20 = edges.back();
  s_20.vertices[0] = t2;
  s_20.vertices[1] = t0;
  // meshVertices -> all UINT32_MAX.
  s_20.meshVertices[0] = UINT32_MAX;
  s_20.meshVertices[1] = UINT32_MAX;
  s_20.meshVertices[2] = UINT32_MAX;
  s_20.meshVertices[3] = UINT32_MAX;
  s_20.meshVertices[4] = UINT32_MAX;
  s_20.meshVertices[5] = UINT32_MAX;
  //s_20.conn_faces.push_back(face_ind);

  // The other conn_faces element is the new triangle that is going to be
  // created below.

  // Add an element to edge_visit
  edge_visit.push_back(0);

  // Remove triangle <t0,t1,t2> from current BSPface and save it as a new one.
  face.edges.pop_back();
  face.edges[ face.edges.size()-1 ] = s_20_ind;
  uint64_t i=0;
  //REMOVE_ELEM_VECT(face_ind, edges[s_01_ind].conn_faces);
  //REMOVE_ELEM_VECT(face_ind, edges[s_12_ind].conn_faces);

  // New face-triangle is created.
  faces.push_back( BSPface(face.meshVertices[0],
                           face.meshVertices[1],
                           face.meshVertices[2],
                           face.conn_cells[0],
                           face.conn_cells[1],
                           face.colour           ) );
  uint64_t new_face_ind = faces.size() -1;
  BSPface& new_face = faces.back();
  new_face.edges.push_back(s_20_ind);
  new_face.edges.push_back(s_12_ind);
  new_face.edges.push_back(s_01_ind);

  cells[faces[face_ind].conn_cells[0] ].faces.push_back(new_face_ind);
  if( !IS_GHOST_CELL(faces[face_ind].conn_cells[1]) )
     cells[faces[face_ind].conn_cells[1] ].faces.push_back(new_face_ind);

  //edges[s_01_ind].conn_faces.push_back(new_face_ind);
  //edges[s_12_ind].conn_faces.push_back(new_face_ind);
  //s_20.conn_faces.push_back(new_face_ind);
  edges[s_01_ind].conn_face_0 = new_face_ind;
  edges[s_12_ind].conn_face_0 = new_face_ind;
  s_20.conn_face_0 = new_face_ind;

  #ifdef DEBUG_BSP_DEEP
  printf("detached triangle %llu -> <%u,%u,%u>\n", new_face_ind, t0, t1, t2);
  triangular_BSPface_isDegenerate(faces, edges, vertices, new_face_ind);
  #endif

}

//  Input:
// Output:
// Note. assumes that aligned face-edges are sub-edges of the same original edge.
bool BSPcomplex::aligned_face_edges(uint64_t fe0, uint64_t fe1, const BSPface& face){
  BSPedge& edge0 = edges[ face.edges[fe0] ];
  BSPedge& edge1 = edges[ face.edges[fe1] ];
  if(edge0.meshVertices[0] == UINT32_MAX) return false;
  if(edge0.meshVertices[0] != edge1.meshVertices[0]) return false;
  if(edge0.meshVertices[1] != edge1.meshVertices[1]) return false;
  if(edge0.meshVertices[2] != edge1.meshVertices[2]) return false;
  if(edge0.meshVertices[2] == UINT32_MAX) return true;
  if(edge0.meshVertices[3] != edge1.meshVertices[3]) return false;
  if(edge0.meshVertices[4] != edge1.meshVertices[4]) return false;
  if(edge0.meshVertices[5] != edge1.meshVertices[5]) return false;
  return true;
}

//
//
//  Input: the index of a BSPface: face_ind.
// Output: nothing.
void BSPcomplex::triangulateFace(uint64_t face_ind){

  // edges indices, of BSPface faces[face_ind], are listed in the vector
  // faces[face_ind].edges ordered as walking face-boundary clockwise
  // (or counterclockwise).
  BSPface& face = faces[face_ind];
  uint64_t num_face_edges = face.edges.size();

  while(num_face_edges > 3){

    // Check if last two edges are not-aligned.
    if(!aligned_face_edges(num_face_edges-1, num_face_edges-2, faces[face_ind]) ){

      // To remove the triangle with the vertices of the last two edges
      // one must be sure that the remaining face does not become degenerate.

      if(!aligned_face_edges(0, num_face_edges-3, faces[face_ind]) ){

        triangle_detach(face_ind);
        num_face_edges--;
      }
      else UINT64_vect_down_shift(faces[face_ind].edges, 1);

    }
    else UINT64_vect_down_shift(faces[face_ind].edges, 1);

  }
}

//  Input: vertices indices of a BSPelement: vrts.
// Output: nothing.
// Note. a point representing the baricenter (or its approximation) is created
//       and added to the vector vertices.
void BSPcomplex::computeBaricenter(const vector<uint32_t>& vrts){
    double cx, cy, cz;
    double sum_x=0.0, sum_y=0.0, sum_z=0.0;
    uint32_t np = 0;
    for (const uint32_t v : vrts)
        if (vertices[v]->getApproxXYZCoordinates(cx, cy, cz))
        {
            sum_x += cx;
            sum_y += cy;
            sum_z += cz;
            np++;
            break; // This line should be commented to have an actual barycenter !!!!!!
        }

    vertices.push_back(new explicitPoint3D(sum_x / np, sum_y / np, sum_z / np));
    vrts_visit.push_back(0);
}

//
//
inline uint64_t BSPcomplex::triFace_oppEdge(const BSPface& face, uint32_t v){
  uint64_t edge_ind = face.edges[0];
  uint32_t e0 = edges[ edge_ind ].vertices[0];
  uint32_t e1 = edges[ edge_ind ].vertices[1];
  if(e0 == v || e1 == v){
    edge_ind = face.edges[1];
    e0 = edges[ edge_ind ].vertices[0];
    e1 = edges[ edge_ind ].vertices[1];
    if(e0 == v || e1 == v){
      edge_ind = face.edges[2];
      e0 = edges[ edge_ind ].vertices[0];
      e1 = edges[ edge_ind ].vertices[1];
    }
  }
  return edge_ind;
}

//
//
uint64_t BSPcomplex::triFace_shareEdge(const BSPcell& cell, uint64_t face_ind,
                                       uint64_t vOppEdge_ind){
  for(uint64_t f=0; f<cell.faces.size(); f++){
    uint64_t adj_face_ind = cell.faces[f];
    if( (faces[ adj_face_ind ].edges[0] == vOppEdge_ind ||
         faces[ adj_face_ind ].edges[1] == vOppEdge_ind ||
         faces[ adj_face_ind ].edges[2] == vOppEdge_ind   ) &&
        face_ind != adj_face_ind                               )
      return adj_face_ind;
  }

  // never reached
  printf("[BSP.cpp]BSPcomplex::triFace_shareEdge: ERROR return UINT64_MAX.\n");
  return UINT64_MAX;
}

//
//
bool BSPcomplex::cell_is_tetrahedrizable_from_v(const BSPcell& cell, uint32_t v){
  uint64_t num_incFaces = count_cellFaces_inc_cellVrt(cell, v);
  vector<uint64_t> v_incFaces(num_incFaces, UINT64_MAX);
  cell_VFrelation(cell, v, v_incFaces);

  //bool return_zero = false;

  for(uint64_t f=0; f<num_incFaces; f++){
    //return_zero = false;
    BSPface& face = faces[ v_incFaces[f] ];
    uint64_t vOppEdge_ind = triFace_oppEdge(face, v);
    uint64_t faceShareEdge_ind = triFace_shareEdge(cell, v_incFaces[f], vOppEdge_ind);
    BSPface& oppFace = faces[faceShareEdge_ind];
    if(oppFace.meshVertices[0] == face.meshVertices[0] &&
       oppFace.meshVertices[1] == face.meshVertices[1] &&
       oppFace.meshVertices[2] == face.meshVertices[2]   ) //return_zero = true;
       return false;
  }

  return true;
}


//
//
void BSPcomplex::makeTetrahedra()
{
    uint64_t tet_num = 0; // total number of tetrahedra in which the cell will
                          // be decomposed.
    std::vector<uint32_t> decomposition_type(cells.size(), 0);
    // Possible decoposition techniques are:
    // 0 - cell is a tetrahedron -> no decomposition.
    // 1 - cell can be decomposed in tetrahedra by connecting a vertex with
    //     the other vertices that do not belong to its link.
    // 2 - cell is decomposed by connecting its baricenter with all the cell's
    //     vertices.
    std::vector<uint32_t> decomposition_vrt(cells.size(), UINT32_MAX);
    // decomposition_vrt is:
    // - UINT32_MAX for a tetrhedron,
    // - a cell vertex for decomposition type 1,
    // - the cell baricenter for decomposition type 2.

    for(uint64_t cell_i = 0; cell_i < cells.size(); cell_i++) {
        BSPcell& cell = cells[cell_i];
        if (cell.place != INTERNAL_A) continue;

        // If cell has more than 4 faces -> chose between types 1 and 2
        if(cell.faces.size() > 4) {

          uint64_t num_cellEdges = UINT64_MAX;
          uint32_t num_cellVrts = count_cellVertices(cell, &num_cellEdges);
          vector<uint32_t> cell_vrts(num_cellVrts, UINT32_MAX);
          list_cellVertices(cell, num_cellEdges, cell_vrts);

          // Check if cell is tetrahedralizable from a vertex.
          bool needs_barycenter = true;
          for(uint32_t cv=0; cv<cell_vrts.size(); cv++)
            if(cell_is_tetrahedrizable_from_v(cell, cell_vrts[cv]) ){
              decomposition_type[cell_i] = 1;
              decomposition_vrt[cell_i] = cell_vrts[cv];
              tet_num += cell.faces.size() - count_cellFaces_inc_cellVrt(cell, cell_vrts[cv]);
              needs_barycenter = false;
              break;
            }

          if(needs_barycenter){ // Cell need baricenter
            decomposition_type[cell_i] = 2;
            computeBaricenter(cell_vrts);
            decomposition_vrt[cell_i] = ((uint32_t)vertices.size() - 1);
            tet_num += cell.faces.size();
          }

        }
        else  tet_num++; // The cell is a tet.
    }

    final_tets.reserve(tet_num * 4);

    // Make tets
    for(uint64_t cell_i = 0; cell_i < cells.size(); cell_i++) {
        BSPcell& cell = cells[cell_i];
        if (cell.place != INTERNAL_A) continue;
        if(decomposition_type[cell_i] == 0){ // Simple tet
        vector<uint32_t> cell_vrts(4, UINT32_MAX);
        list_cellVertices(cells[cell_i], 6, cell_vrts);
        final_tets.insert(final_tets.end(), cell_vrts.begin(), cell_vrts.end());
      }
      else if(decomposition_type[cell_i] == 1){ // Tetrahedralizable from vertex
        uint32_t v = decomposition_vrt[cell_i];
        uint64_t num_incFaces = count_cellFaces_inc_cellVrt(cells[cell_i], v);
        uint64_t num_NOT_incFaces = cells[cell_i].faces.size() - num_incFaces;
        vector<uint64_t> v_NOT_incFaces(num_NOT_incFaces, UINT64_MAX);
        COMPL_cell_VFrelation(cells[cell_i], v, v_NOT_incFaces);
        for(uint64_t face_i : v_NOT_incFaces){
          // Simple triangle
          vector<uint32_t> face_vrts(3, UINT32_MAX);
          list_faceVertices(faces[face_i], face_vrts);
          final_tets.insert(final_tets.end(), face_vrts.begin(), face_vrts.end());
          final_tets.push_back(v);
        }
      }
      else{ // Uses cell barycenter
        for(uint64_t face_i : cells[cell_i].faces) {
          // Simple triangle
          vector<uint32_t> face_vrts(3, UINT32_MAX);
          list_faceVertices(faces[face_i], face_vrts);
          final_tets.insert(final_tets.end(), face_vrts.begin(), face_vrts.end());
          final_tets.push_back(decomposition_vrt[cell_i]);
        }
      }
    }

}




//-Decide colour of GREY faces--------------------------------------------------

//
//
int BSPcomplex::face_dominant_normal_component(const BSPface& face) {
    const uint32_t* mv = face.meshVertices;
    double mvc[9]; // Coords of the original input triangle
    vertices[mv[0]]->getApproxXYZCoordinates(mvc[0], mvc[1], mvc[2]);
    vertices[mv[1]]->getApproxXYZCoordinates(mvc[3], mvc[4], mvc[5]);
    vertices[mv[2]]->getApproxXYZCoordinates(mvc[6], mvc[7], mvc[8]);
    return genericPoint::maxComponentInTriangleNormal(mvc[0], mvc[1], mvc[2],
        mvc[3], mvc[4], mvc[5],
        mvc[6], mvc[7], mvc[8]);
}

//
//
void BSPcomplex::get_approx_faceBaricenterCoord(const BSPface& face, double* bar){
  bar[0] = 0;
  bar[1] = 0;
  bar[2] = 0;
  double tp[3];
  const BSPedge& edge0 = edges[face.edges.back()];
  const BSPedge& edge1 = edges[face.edges[0]];
  uint32_t vid = consecEdges_common_endpt(edge0.vertices[0], edge0.vertices[1],
                                          edge1.vertices[0], edge1.vertices[1]);
  for(uint64_t e = 0; e < face.edges.size(); e++) {
      const BSPedge& edge = edges[face.edges[e]];
      if (vid == edge.vertices[0]) vid = edge.vertices[1];
      else vid = edge.vertices[0];
      vertices[vid]->getApproxXYZCoordinates(tp[0], tp[1], tp[2]);
      bar[0] += tp[0];
      bar[1] += tp[1];
      bar[2] += tp[2];
  }
  bar[0] /= face.edges.size();
  bar[1] /= face.edges.size();
  bar[2] /= face.edges.size();
}

//
//
bool BSPcomplex::is_baricenter_inFace(const BSPface& face,
                                      const explicitPoint3D& face_center,
                                      int max_normComp){
  const BSPedge& edge0 = edges[face.edges.back()];
  const BSPedge& edge1 = edges[face.edges[0]];
  uint32_t vid = consecEdges_common_endpt(edge0.vertices[0], edge0.vertices[1],
                                          edge1.vertices[0], edge1.vertices[1]);

  int oro = 0;
  uint64_t e;
  for(e = 0; e < face.edges.size(); e++) {
      const BSPedge& edge = edges[face.edges[e]];
      const uint32_t pvid = vid;
      if (vid == edge.vertices[0]) vid = edge.vertices[1];
      else vid = edge.vertices[0];

      const int ao = genericPoint::orient2D(face_center, *vertices[pvid], *vertices[vid], max_normComp);
      if (ao == 0) break;
      if (ao != oro)
      {
          if (oro) break;
          else oro = ao;
      }
  }

  if(e == face.edges.size() ) return true;
  return false;
}

//
//
inline bool faceColour_matches_constrGroup(CONSTR_GROUP_T group, bool f_blackA, bool f_blackB){
  return (group == CONSTR_A && f_blackA) || (group == CONSTR_B && f_blackB);
}

//
//
COLOUR_T BSPcomplex::blackAB_or_white(uint64_t face_ind, bool two_input){
  const BSPface& face = faces[face_ind];

  // Get dominant normal component
  int xyz = face_dominant_normal_component(face);

  // Calculate approximated face barycenter
  double p[3];
  get_approx_faceBaricenterCoord(face, p);
  const explicitPoint3D face_center(p[0], p[1], p[2]);

  // Needed only for two input case.
  bool face_is_blackA = false;
  bool face_is_blackB = false;

  // Check whether the barycenter is indeed inside the face (might be not due to approximation)
  if ( is_baricenter_inFace(face, face_center, xyz) ) // Barycenter is inside the face: just check that it is in one of the constraints too
  {
      for (uint32_t c = 0; c < face.coplanar_constraints.size(); c++) {
          const uint32_t constr = face.coplanar_constraints[c];
          const CONSTR_GROUP_T c_group = constraint_group[constr];

          if( two_input && faceColour_matches_constrGroup(c_group, face_is_blackA, face_is_blackB) ) continue;

          const uint32_t constr_ID = 3 * constr;
          const genericPoint* c0 = vertices[ constraints_vrts[constr_ID   ] ];
          const genericPoint* c1 = vertices[ constraints_vrts[constr_ID +1] ];
          const genericPoint* c2 = vertices[ constraints_vrts[constr_ID +2] ];

          if (genericPoint::pointInTriangle(face_center, *c0, *c1, *c2, xyz)){

            if(!two_input) return BLACK_A;

            if(c_group == CONSTR_A)      face_is_blackA = true;
            else if(c_group == CONSTR_B) face_is_blackB = true;
            if(face_is_blackA && face_is_blackB) return BLACK_AB;
          }
      }

      if(two_input){
        if(face_is_blackA) return BLACK_A;
        else if(face_is_blackB) return BLACK_B;
      }

      return WHITE;
  }
  else // Barycenter is not inside the face: revert to slow version
  {
      const BSPedge& edge0 = edges[face.edges.back()];
      const BSPedge& edge1 = edges[face.edges[0]];
      uint32_t vid = consecEdges_common_endpt(edge0.vertices[0], edge0.vertices[1],
                                              edge1.vertices[0], edge1.vertices[1]);

      for (uint64_t e = 0; e < face.edges.size(); e++) {
          const BSPedge& edge = edges[face.edges[e]];
          if (vid == edge.vertices[0]) vid = edge.vertices[1];
          else vid = edge.vertices[0];

          uint32_t out_from_all = 0;
          const genericPoint* face_pt = vertices[vid];
          for (uint32_t c = 0; c < face.coplanar_constraints.size(); c++) {
              const uint32_t constr = face.coplanar_constraints[c];
              const CONSTR_GROUP_T c_group = constraint_group[constr];

              if( two_input && faceColour_matches_constrGroup(c_group, face_is_blackA, face_is_blackB) ) continue;

              const uint32_t constr_ID = 3 * constr;
              const uint32_t vid1 = constraints_vrts[constr_ID];
              const uint32_t vid2 = constraints_vrts[constr_ID + 1];
              const uint32_t vid3 = constraints_vrts[constr_ID + 2];
              const genericPoint* c0 = vertices[vid1];
              const genericPoint* c1 = vertices[vid2];
              const genericPoint* c2 = vertices[vid3];

              if (vid == vid1 || vid == vid2 || vid == vid3) break; // point on boundary.
              int lpt = localizedPointInTriangle(*face_pt, *c0, *c1, *c2, xyz);
              // lpt = 1 -> point on boundary, lpt = 2 -> point in interior, otherwise lpt = 0.

              if (lpt == 2){

                if(!two_input) return BLACK_A;

                if(c_group == CONSTR_A) face_is_blackA = true;
                else if(c_group == CONSTR_B) face_is_blackB = true;
                if(face_is_blackA && face_is_blackB) return BLACK_AB;
              }

              if (lpt)  break;

              out_from_all++;
          }
          if (out_from_all == face.coplanar_constraints.size()) return WHITE;
      }

      // All face vertices are on the boundary of some coplanar constraints.
      for (uint32_t c = 0; c < face.coplanar_constraints.size(); c++) {

          const uint32_t constr = face.coplanar_constraints[c];
          const CONSTR_GROUP_T c_group = constraint_group[ constr ];
          if( two_input && faceColour_matches_constrGroup(c_group, face_is_blackA, face_is_blackB) ) continue;

          const uint32_t *constraint = constraints_vrts.data() + 3*constr;
          if (coplanar_constraint_innerIntersects_face(face.edges, constraint, xyz)){

            if(!two_input) return BLACK_A;

            if(c_group == CONSTR_A) face_is_blackA = true;
            else if(c_group == CONSTR_B) face_is_blackB = true;
            if(face_is_blackA && face_is_blackB) return BLACK_AB;

          }
      }

      if(face_is_blackA) return BLACK_A;
      else if(face_is_blackB) return BLACK_B;

      return WHITE;
  }

}


//----------------
// BSP subdivision
//----------------

//
//

void BSPcomplex::extractSkinTriMesh(const char* filename, const char bool_opcode,
    double** coords, uint32_t* npts, uint32_t** tri_idx, uint32_t* ntri)
{
    const uint64_t num_faces = faces.size();
    for (uint64_t face_ind = 0; face_ind < num_faces; face_ind++) triangulateFace(face_ind);

    if (bool_opcode == 'U') { // Union
        for (BSPcell& cell : cells)
            cell.place = (cell.place == INTERNAL_A || cell.place == INTERNAL_B || cell.place == INTERNAL_AB) ? (INTERNAL_A) : (EXTERNAL);
    }

    else if (bool_opcode == 'I') { // Intersection
        for (BSPcell& cell : cells)
            cell.place = (cell.place == INTERNAL_AB) ? (INTERNAL_A) : (EXTERNAL);
    }

    else if (bool_opcode == 'D') { // Difference A\B
        for (BSPcell& cell : cells)
            cell.place = (cell.place == INTERNAL_A && !(cell.place == INTERNAL_AB)) ? (INTERNAL_A) : (EXTERNAL);
    }

    // Set "internal" depending on bool_opcode and find border faces to save
    vector<uint64_t> mark(faces.size(), 0);
    for (BSPcell& cell : cells)
        if (cell.place == INTERNAL_A)
            for (uint64_t fi = 0; fi < cell.faces.size(); fi++)
                mark[cell.faces[fi]]++;

    for (size_t i = 0; i < vrts_visit.size(); i++) vrts_visit[i] = 0;
    for (size_t i = 0; i < edges.size(); i++) edge_visit[i] = 0;

    uint64_t num_border_faces = 0;
    for (uint64_t f_i = 0; f_i < faces.size(); f_i++) if (mark[f_i] == 1)
    {
        num_border_faces++;
        for (uint64_t eid : faces[f_i].edges) edge_visit[eid] = 1;
    }
    for (size_t i = 0; i < edges.size(); i++) if (edge_visit[i])
        vrts_visit[edges[i].vertices[0]] = vrts_visit[edges[i].vertices[1]] = 1;

    std::vector<uint32_t> vmap(vertices.size(), 0);
    size_t num_v = 0;
    for (size_t i = 0; i < vrts_visit.size(); i++)
    {
        vmap[i] = (uint32_t)num_v;
        if (vrts_visit[i]) num_v++;
    }

    *coords = (double*)malloc(sizeof(double) * 3 * num_v);
    *tri_idx = (uint32_t*)malloc(sizeof(uint32_t) * 3 * num_border_faces);
    *npts = num_v;
    *ntri = num_border_faces;

    // Store vertex coordinates
    for (uint32_t v = 0, i = 0; v < vertices.size(); v++) if (vrts_visit[v])
    {
        vertices[v]->getApproxXYZCoordinates((*coords)[i], (*coords)[i+1], (*coords)[i+2], true);
        i += 3;
    }

    // Print border faces
    for (uint32_t f_i = 0, i = 0; f_i < faces.size(); f_i++)
        if (mark[f_i] == 1) {
            BSPface& face = faces[f_i];
            vector<uint32_t> face_vrts(face.edges.size(), UINT32_MAX);
            list_faceVertices(face, face_vrts);

            if (cells[face.conn_cells[0]].place == INTERNAL_A)
                for (uint32_t v = (uint32_t)face_vrts.size(); v > 0; v--) (*tri_idx)[i++] = vmap[face_vrts[v - 1]];
            else
                for (uint32_t v = 0; v < face_vrts.size(); v++) (*tri_idx)[i++] = vmap[face_vrts[v]];
        }
}


void BSPcomplex::saveSkin(const char *filename, const char bool_opcode, bool triangulate) {

    if (triangulate)
    {
        const uint64_t num_faces = faces.size();
        for (uint64_t face_ind = 0; face_ind < num_faces; face_ind++) triangulateFace(face_ind);
    }

    if (bool_opcode == 'U') { // Union
        for (BSPcell& cell : cells)
            cell.place = (cell.place == INTERNAL_A || cell.place == INTERNAL_B || cell.place == INTERNAL_AB) ? (INTERNAL_A) : (EXTERNAL);
    }

    else if (bool_opcode == 'I') { // Intersection
        for (BSPcell& cell : cells)
            cell.place = (cell.place == INTERNAL_AB) ? (INTERNAL_A) : (EXTERNAL);
    }

    else if (bool_opcode == 'D') { // Difference A\B
        for (BSPcell& cell : cells)
            cell.place = (cell.place == INTERNAL_A && !(cell.place == INTERNAL_AB)) ? (INTERNAL_A) : (EXTERNAL);
    }

    // Set "internal" depending on bool_opcode and find border faces to save
    vector<uint64_t> mark(faces.size(), 0);
    for (BSPcell& cell : cells)
            if (cell.place == INTERNAL_A)
                for (uint64_t fi = 0; fi < cell.faces.size(); fi++)
                    mark[cell.faces[fi]]++;

    for (size_t i = 0; i < vrts_visit.size(); i++) vrts_visit[i] = 0;
    for (size_t i = 0; i < edges.size(); i++) edge_visit[i] = 0;

    uint64_t num_border_faces = 0;
    for (uint64_t f_i = 0; f_i < faces.size(); f_i++) if (mark[f_i] == 1)
    {
         num_border_faces++;
         for (uint64_t eid : faces[f_i].edges) edge_visit[eid] = 1;
    }
    for (size_t i = 0; i < edges.size(); i++) if (edge_visit[i])
        vrts_visit[edges[i].vertices[0]] = vrts_visit[edges[i].vertices[1]] = 1;

    std::vector<uint32_t> vmap(vertices.size(), 0);
    size_t num_v = 0;
    for (size_t i = 0; i < vrts_visit.size(); i++)
    {
        vmap[i] = (uint32_t)num_v;
        if (vrts_visit[i]) num_v++;
    }

    ofstream f(filename);

    if (!f)
        ip_error("BSPcomplex::saveSkin: cannot open the file.\n");

    f << "OFF\n";
    f << num_v << " ";
    f << num_border_faces << " ";
    f << "0\n";

    // Print vertices coordinates
    for (uint32_t v = 0; v < vertices.size(); v++) if (vrts_visit[v])
        f << (*(vertices)[v]) << "\n";

    // Print border faces
    for (uint32_t f_i = 0; f_i < faces.size(); f_i++)
        if (mark[f_i] == 1) {
            BSPface& face = faces[f_i];
            vector<uint32_t> face_vrts(face.edges.size(), UINT32_MAX);
            list_faceVertices(face, face_vrts);
            f << face_vrts.size();

            if (cells[face.conn_cells[0]].place == INTERNAL_A)
                for (uint32_t v = (uint32_t)face_vrts.size(); v > 0; v--) f << " " << vmap[face_vrts[v - 1]];
            else
                for (uint32_t v = 0; v < face_vrts.size(); v++) f << " " << vmap[face_vrts[v]];

            f << "\n";
        }

    f.close();

}


//
//
void BSPcomplex::saveMesh(const char* filename, const char bool_opcode, bool tetrahedrize)
{
    ofstream f(filename);

    if (!f) ip_error("\nBSPcomplex::[BSP.cpp]saveTetMesh: FATAL ERROR cannot open the file.\n");

    const uint64_t num_faces = faces.size();

    if (tetrahedrize)
        for (uint64_t face_ind = 0; face_ind < num_faces; face_ind++) triangulateFace(face_ind);

    if (bool_opcode == 'U') { // Union
        for (BSPcell& cell : cells)
            cell.place = (cell.place == INTERNAL_A || cell.place == INTERNAL_B || cell.place == INTERNAL_AB) ? (INTERNAL_A) : (EXTERNAL);
    }

    else if (bool_opcode == 'I') { // Intersection
        for (BSPcell& cell : cells)
            cell.place = (cell.place == INTERNAL_AB) ? (INTERNAL_A) : (EXTERNAL);
    }

    else if (bool_opcode == 'D') { // Difference A\B
        for (BSPcell& cell : cells)
            cell.place = (cell.place == INTERNAL_A && !(cell.place == INTERNAL_AB)) ? (INTERNAL_A) : (EXTERNAL);
    }

    if (tetrahedrize)
    {
        makeTetrahedra();

        for (size_t i = 0; i < vrts_visit.size(); i++) vrts_visit[i] = 0;
        for (uint32_t t = 0; t < final_tets.size(); t++) vrts_visit[final_tets[t]] = 1;
        uint32_t final_numver = 0;
        for (size_t i = 0; i < vrts_visit.size(); i++) if (vrts_visit[i]) final_numver++;

        std::vector<uint32_t> vmap(vertices.size(), 0);
        size_t num_v = 0;
        for (size_t i = 0; i < vrts_visit.size(); i++)
        {
            vmap[i] = (uint32_t)num_v;
            if (vrts_visit[i]) num_v++;
        }

        f << final_numver << " vertices\n";
        f << final_tets.size() / 4 << " tets\n";

        // Print vertices coordinates
        for (uint32_t v = 0; v < vertices.size(); v++) if (vrts_visit[v]) f << (*vertices[v]) << "\n";

        // Print tets
        for (uint32_t t = 0; t < final_tets.size(); t += 4)
            f << "4 " << vmap[final_tets[t]] << " "
            << vmap[final_tets[t + 1]] << " "
            << vmap[final_tets[t + 2]] << " "
            << vmap[final_tets[t + 3]] << "\n";
    }
    else
    {
        size_t internal_cell_num = 0;
        std::vector<uint32_t> face_visit(faces.size(), 0);
        for (size_t i = 0; i < vrts_visit.size(); i++) vrts_visit[i] = 0;
        for (size_t i = 0; i < edge_visit.size(); i++) edge_visit[i] = 0;
        for (BSPcell& cell : cells) if (cell.place == INTERNAL_A)
        {
            internal_cell_num++;
            for (uint64_t f : cell.faces) face_visit[f] = 1;
        }
        for (size_t f = 0; f < faces.size(); f++) if (face_visit[f]) for (uint64_t e : faces[f].edges) edge_visit[e] = 1;
        for (size_t e = 0; e < edges.size(); e++) if (edge_visit[e]) vrts_visit[edges[e].vertices[0]] = vrts_visit[edges[e].vertices[1]] = 1;
        uint32_t final_numver = 0;
        for (size_t i = 0; i < vrts_visit.size(); i++) if (vrts_visit[i]) final_numver++;

        std::vector<uint32_t> vmap(vertices.size(), 0);
        size_t num_v = 0;
        for (size_t i = 0; i < vrts_visit.size(); i++)
        {
            vmap[i] = (uint32_t)num_v;
            if (vrts_visit[i]) num_v++;
        }

        f << final_numver << " vertices\n";
        f << internal_cell_num << " cells\n";

        // Print vertices coordinates
        for (uint32_t v = 0; v < vertices.size(); v++) if (vrts_visit[v]) f << (*vertices[v]) << "\n";

        // Print cells
        for (size_t v = 0; v < vertices.size(); v++) vrts_visit[v] = 0;
        for (size_t e = 0; e < edges.size(); e++) edge_visit[e] = 0;
        for (BSPcell& cell : cells) if (cell.place == INTERNAL_A)
        {
            uint64_t numce = count_cellEdges(cell);
            uint32_t numcv = count_cellVertices(cell, &numce);
            std::vector<uint32_t> cell_vrts(numcv);
            list_cellVertices(cell, numce, cell_vrts);
            f << cell_vrts.size() << " ";
            for (uint32_t i : cell_vrts) f << vmap[i] << " ";
            f << "\n";
        }
    }

    f.close();

    final_tets.clear();
}


void BSPcomplex::saveBlackFaces(const char* filename, bool triangulate) {
    ofstream f(filename);

    if (!f)
        ip_error("BSPcomplex::saveBlackFaces: cannot open the file.\n");

    if (triangulate)
    {
        const uint64_t num_faces = faces.size();
        for (uint64_t face_ind = 0; face_ind < num_faces; face_ind++) triangulateFace(face_ind);
    }

    for (size_t i = 0; i < vrts_visit.size(); i++) vrts_visit[i] = 0;
    for (size_t i = 0; i < edges.size(); i++) edge_visit[i] = 0;

    uint64_t num_border_faces = 0;
    for (uint64_t f_i = 0; f_i < faces.size(); f_i++) if (faces[f_i].colour != WHITE)
    {
        num_border_faces++;
        for (uint64_t eid : faces[f_i].edges) edge_visit[eid] = 1;
    }
    for (size_t i = 0; i < edges.size(); i++) if (edge_visit[i])
        vrts_visit[edges[i].vertices[0]] = vrts_visit[edges[i].vertices[1]] = 1;

    std::vector<uint32_t> vmap(vertices.size(), 0);
    size_t num_v = 0;
    for (size_t i = 0; i < vrts_visit.size(); i++)
    {
        vmap[i] = (uint32_t)num_v;
        if (vrts_visit[i]) num_v++;
    }

    f << "OFF\n";
    f << num_v << " ";
    f << num_border_faces << " ";
    f << "0\n";

    // Print vertices coordinates
    for (uint32_t v = 0; v < vertices.size(); v++) if (vrts_visit[v])
        f << (*(vertices)[v]) << "\n";

    // Print border faces
    for (uint32_t f_i = 0; f_i < faces.size(); f_i++)
        if (faces[f_i].colour != WHITE) {
            BSPface& face = faces[f_i];
            vector<uint32_t> face_vrts(face.edges.size(), UINT32_MAX);
            list_faceVertices(face, face_vrts);
            f << face_vrts.size();

            for (uint32_t v = 0; v < face_vrts.size(); v++) f << " " << vmap[face_vrts[v]];

            f << "\n";
        }

    f.close();
}


size_t BSPcomplex::getStructureSize() const
{
    size_t tot = 0;

    // Size for points
    for (genericPoint* p : vertices)
    {
        if (p->isExplicit3D()) tot += sizeof(explicitPoint3D);
        else if (p->isLPI()) tot += sizeof(implicitPoint3D_LPI);
        else tot += sizeof(implicitPoint3D_TPI);
    }

    // Size of the array of pointers
    tot += vertices.size() * sizeof(genericPoint*);

    // Size of edge objects
    tot += edges.size() * sizeof(BSPedge);

    // Size of face objects (including pointed arrays)
    for (const BSPface& f : faces) tot += f.getSize();

    // Size of cell objects (including pointed arrays)
    for (const BSPcell& c : cells) tot += c.getSize();

    // And all the other vectors use by the structure...
    tot += sizeof(uint32_t) * constraints_vrts.size();
    tot += sizeof(CONSTR_GROUP_T) * constraint_group.size();
    tot += sizeof(uint32_t) * final_tets.size();
    tot += sizeof(char) * vrts_orBin.size();
    tot += sizeof(uint32_t) * vrts_visit.size();
    tot += sizeof(uint64_t) * edge_visit.size();
    tot += sizeof(BSPcomplex);

    return tot;
}

