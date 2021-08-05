#ifndef BSP_h
#define BSP_h

#include <stdio.h>
#include <list>
#include <vector>
#include <ctype.h>
#include <utility>
#include "delaunay.h"
#include "conforming_mesh.h"
#include "implicit_point.h"

// BSPface colours:
// WHITE -> no 2D-intersection with constraints.
// GREY -> partially or uncertain 2D-intersection with constraints.
// BLACK_A -> fully contained in one constraint of input A.
// BLACK_B -> fully contained in one constraint of input B. (two input case)
// BLACK_B -> fully contained in one constraint of both input. (two input case)
#define WHITE 0
#define BLACK_A 1
#define BLACK_B 2
#define BLACK_AB 3
#define GREY 4
#define COLOUR_T uint32_t

// Cell places
#define UNDEFINED 4
#define INTERNAL_A 1
#define INTERNAL_B 2
#define INTERNAL_AB 3
#define EXTERNAL 0

// To manage with constraints edges
#define ENDPTS_T pair<uint32_t,uint32_t>
#define NEW_ENDPTS make_pair(UINT32_MAX, UINT32_MAX)

class BSPedge{  // The edge of a BSPcell.
public:
    uint32_t meshVertices[6];         // 2 vertices of the mesh-tet-edge which
                                      // contains BSPedge (and 4 UINT32_MAX
                                      // in the 3rd to 5th position)
                                      // OR
                                      // 6 vertices of the 2 mesh-triangles
                                      // whose intersection defines the edge.

    uint32_t vertices[2];             // BSPvertices belonging to BSPedge.

    uint64_t conn_face_0;             // One of the incident faces

    BSPedge(){}

    BSPedge(uint32_t v1, uint32_t v2, uint32_t m_v1, uint32_t m_v2){
        vertices[0] = v1;
        vertices[1] = v2;
        meshVertices[0] = m_v1;
        meshVertices[1] = m_v2;
        meshVertices[2] = UINT32_MAX;
        meshVertices[3] = UINT32_MAX;
        meshVertices[4] = UINT32_MAX;
        meshVertices[5] = UINT32_MAX;
    }

    BSPedge(uint32_t v1, uint32_t v2,
            uint32_t m_t1v1, uint32_t m_t1v2, uint32_t m_t1v3,
            uint32_t m_t2v1, uint32_t m_t2v2, uint32_t m_t2v3 ){
        vertices[0] = v1;
        vertices[1] = v2;
        meshVertices[0] = m_t1v1;
        meshVertices[1] = m_t1v2;
        meshVertices[2] = m_t1v3;
        meshVertices[3] = m_t2v1;
        meshVertices[4] = m_t2v2;
        meshVertices[5] = m_t2v3;
    }

    BSPedge split(uint32_t new_point);
};

class BSPface{  // The face of a BSPcell.
public:
    std::vector<uint64_t> edges; // BSPedges bounding the face:
                                 //   the position (i) of an edge in the vector
                                 //   is such that the previous (i-1) and next
                                 //   (i+i) edge are consecutive by walking the
                                 //   face boundary.
    uint64_t conn_cells[2];      // The two cells that share this face.
    std::vector<uint32_t> coplanar_constraints;
    uint32_t meshVertices[3];    // 3 vertices of the mesh-tet-face
                                 // which contains BSPedge.
    COLOUR_T colour;

    BSPface(){}

    // Face (general) constructor
    BSPface(uint32_t m_v1, uint32_t m_v2, uint32_t m_v3,
            uint64_t c1, uint64_t c2){

        meshVertices[0] = m_v1;
        meshVertices[1] = m_v2;
        meshVertices[2] = m_v3;

        conn_cells[0] = c1;
        conn_cells[1] = c2;
    }

    size_t getSize() const
    {
        return sizeof(BSPface) + sizeof(uint64_t) * edges.size() + sizeof(uint32_t) *
            coplanar_constraints.size();
    }

    BSPface(uint32_t m_v1, uint32_t m_v2, uint32_t m_v3,
            uint64_t c1, uint64_t c2, COLOUR_T _colour,
            const vector<uint32_t>& constraints){

        meshVertices[0] = m_v1;
        meshVertices[1] = m_v2;
        meshVertices[2] = m_v3;

        conn_cells[0] = c1;
        conn_cells[1] = c2;

        colour = _colour;
        coplanar_constraints.assign(constraints.begin(), constraints.end());
    }

    // Common face (between splitted cell sub-cells) constructor
    BSPface(uint32_t m_v1, uint32_t m_v2, uint32_t m_v3,
            uint64_t c1, uint64_t c2, COLOUR_T _colour){

        meshVertices[0] = m_v1;
        meshVertices[1] = m_v2;
        meshVertices[2] = m_v3;

        conn_cells[0] = c1;
        conn_cells[1] = c2;

        colour = _colour;
    }

    inline void removeEdge(uint64_t edge);
    inline void exchange_conn_cell(uint64_t cell, uint64_t newCell);
};

class BSPcell{  // A convex polyhedron defined by the intersection of a
                // mesh-tet and a certain number of constraints.
public:
    std::vector<uint64_t> faces;       // BSPfaces bounding the BSPcell.
    std::vector<uint32_t> constraints; // Constraints that intersect
                                       // the BSPcell.
    uint32_t place = UNDEFINED;  // Internal, external or undefined (see macros)
                                 // w.r.t. constraints surface.

    BSPcell(){}

    BSPcell(const vector<uint64_t>& cell_faces){
        faces.assign(cell_faces.begin(), cell_faces.end());
    }

    BSPcell(const vector<uint64_t>& cell_faces,
            const vector<uint32_t>& intersect_constrs){
        faces.assign(cell_faces.begin(), cell_faces.end());
        constraints.assign(intersect_constrs.begin(), intersect_constrs.end());
    }

    inline void removeFace(uint64_t face);

    size_t getSize() const
    {
        return sizeof(BSPcell) + sizeof(uint64_t) * faces.size() + sizeof(uint32_t) *
            constraints.size();
    }

};

class BSPcomplex{
    public:

    std::vector<genericPoint*> vertices; // mesh vertices + new vertices.
                                         // new vertices -> intersections
                                         // between tetrahedra and constraints.
    std::vector<BSPedge> edges;
    std::vector<BSPface> faces;
    std::vector<BSPcell> cells;
    std::vector<uint32_t> constraints_vrts; // The constraint-triangles are
                                            // constraints_vrts.size()/3 .
    std::vector<CONSTR_GROUP_T> constraint_group;
    uint32_t first_virtual_constraint;

    std::vector<uint32_t> final_tets;       // Simple vector storing the tetrahedra
											// (only used when saving a tet-mesh)


    // Supporting vectors
    std::vector<char> vrts_orBin; // Used to "cache" vertex orientations
                                 //  w.r.t. some plane.
                                 // (same length of vertices)
    std::vector<uint32_t> vrts_visit; // To flag visited vertices when needed
                                      // (same length of vertices)
    std::vector<uint64_t> edge_visit; // To flag visited edges when needed
                                      // (same length of edges)

    BSPcomplex(const TetMesh* mesh, const constraints_t* constraints,
               const uint32_t** map, const uint32_t* num_map,
               const uint32_t** map_f0, const uint32_t* num_map_f0,
               const uint32_t** map_f1, const uint32_t* num_map_f1,
               const uint32_t** map_f2, const uint32_t* num_map_f2,
               const uint32_t** map_f3, const uint32_t* num_map_f3 );

    ~BSPcomplex() { for (genericPoint* v : vertices) delete v; }

    // Save the faces representing the input constraints
    void saveBlackFaces(const char* filename, bool triangulate = false);

    // Save the faces that separate in and out
    void saveSkin(const char* filename, const char bool_opcode, bool triangulate = false);

    // Save the mesh
    void saveMesh(const char* filename, const char bool_opcode, bool tetrahedrize = false);

    // Makes a triangle mesh out of the skin faces
    void extractSkinTriMesh(const char* filename, const char bool_opcode, 
        double** coords, uint32_t* npts, uint32_t** tri_idx, uint32_t* ntri);

    // Return the overall number of bytes occupied by the BSP structure
    size_t getStructureSize() const;

    // Complex elements relations
    inline void assigne_edge_to_face(uint64_t edge, uint64_t face);

    // Explore the complex
    //uint32_t getFaceVertex(const BSPface& f, uint32_t v_ind);
    //bool faceHasVertex(const BSPface& f, uint32_t v_ind);
    uint64_t faceSharedWithCell(uint64_t c1, uint64_t c2);
    uint64_t count_cellEdges(const BSPcell& cell);
    uint32_t count_cellVertices(const BSPcell& cell,
                                       uint64_t* num_cellEdges);
    void list_cellEdges(BSPcell& cell, vector<uint64_t>& cell_edges);
    void list_cellVertices(BSPcell& cell, uint64_t num_cellEdges,
                           vector<uint32_t>& cell_vrts);
    void list_faceVertices(BSPface& face, vector<uint32_t>& face_vrts);
    void fill_cell_locDS(BSPcell& cell, vector<uint64_t>& cell_edges,
                         vector<uint32_t>& cell_vrts);
    inline uint64_t find_face_edge(const BSPface& face, uint32_t v, uint32_t u);
    uint64_t count_cellFaces_inc_cellVrt(const BSPcell& cell, uint32_t v);
    void cell_VFrelation(const BSPcell& cell, uint32_t v,
                         vector<uint64_t>& v_incFaces_ind);
    void COMPL_cell_VFrelation(const BSPcell& cell, uint32_t v,
                                vector<uint64_t>& v_NOT_incFaces_ind);
    bool is_virtual(uint32_t constr_ind);

    uint64_t getOppositeEdgeFace(const uint64_t e0, const uint64_t f0, const uint64_t c);
    void makeEFrelation(const uint64_t e_id, std::vector<uint64_t>& ef);

    // Geometric predicates
    void vrts_orient_wrtPlane(const vector<uint32_t>& vrts_inds,
                    uint32_t plane_pt0, uint32_t plane_pt1, uint32_t plane_pt2, uint32_t count);
    inline void count_vrt_orBin(const vector<uint32_t>& inds,
                                uint32_t* pos, uint32_t* neg, uint32_t* zero);
    inline bool constraint_innerIntersects_edge(const BSPedge& edge,
                                         const vector<uint32_t>& cell_vrts);
    inline bool constraint_innerIntersects_face(const vector<uint32_t>& face_vrts);
    bool coplanar_constraint_innerIntersects_face(const vector<uint64_t>& face_edges,
                                                  const uint32_t constraint[3],
                                                  const int dominant_normal_comp      );

    // Upload Delaunay triangolation
    uint64_t removing_ghost_tets(const TetMesh* mesh, vector<uint64_t>& new_order);
    uint64_t add_tetEdge(const TetMesh* mesh, uint32_t e0, uint32_t e1,
                         uint64_t tet_ind, const vector<uint64_t>& new_order);
    inline uint64_t add_tetFace(uint32_t v0, uint32_t v1, uint32_t v2,
                         uint64_t cell_ind, uint64_t adjCell_ind);
    inline bool tet_face_isNew(uint64_t tet_ind, uint64_t adjTet_ind,
                               uint64_t adjCell_ind);
    inline void fill_face_colour(uint64_t tet_ind, uint64_t face_ind,
                          const uint32_t** map_fi, const uint32_t* num_map_fi);

    // BSPsubdivision
    inline void move_edge(uint64_t edge_face_ind,
                          uint64_t face_ind, uint64_t newFace_ind);
    inline void move_face(uint64_t face_cell_ind, uint64_t cell_ind,
                          uint64_t newCell_ind);
    inline void remove_constraint(uint32_t constr_cell_ind, uint64_t cell_ind);
    void edgesPartition(uint64_t face_ind, uint64_t newFace_ind);
    void facesPartition(uint64_t cell_ind, uint64_t newCell_ind,
                        const vector<uint32_t>& cell_vrts);
    void constraintsPartition(uint32_t ref_constr,
                              uint64_t down_cell_ind, uint64_t up_cell_ind,
                              const vector<uint32_t>& cell_vrts);
    void add_edgeToOrdFaceEdges(BSPface& face, uint64_t newEdge_ind);
    void add_commonEdge(uint32_t constr, uint64_t face_ind, uint64_t newFace_ind,
                        const uint32_t* endpts);
    void add_edges_toCommFaceEdges(BSPface& face,
                                   const vector<uint64_t>& edges_ind);
    void add_commonFace(uint32_t constr,
                        uint64_t cell_ind, uint64_t newCell_ind,
                        const vector<uint32_t>& cell_vrts,
                        const vector<uint64_t>& cell_edges);
    void fixCommonFaceOrientation(uint64_t cf_id);
    uint32_t add_LPIvrt(const BSPedge& edge, uint32_t constr);
    uint32_t add_TPIvrt(const BSPedge& edge, uint32_t constr);
    void splitEdge(uint64_t edge_ind, uint32_t constr);
    void splitFace(uint64_t face_ind, uint32_t constr, uint64_t cell_ind,
                   const vector<uint32_t>& face_vrts);
    void splitCell(uint64_t cell_ind);
    void find_coplanar_constraints(uint64_t cell_ind, uint32_t constr,
                                              vector<uint32_t>& coplanar_c);

    // Decide colour of GREY faces
    int face_dominant_normal_component(const BSPface& face);
    void get_approx_faceBaricenterCoord(const BSPface& face, double* bar);
    bool is_baricenter_inFace(const BSPface& face,
                             const explicitPoint3D& face_center, int max_normComp);
    COLOUR_T blackAB_or_white(uint64_t face_ind, bool two_input);

    // Interior-exterior constraint surface
    void constraintsSurface_complexPartition(bool two_files=false);
    void markInternalCells(uint32_t skin_colour, uint32_t internal_label, const std::vector<double>& face_costs);

    // Tetrahedralization
    void triangle_detach(uint64_t face_ind);
    bool aligned_face_edges(uint64_t fe0, uint64_t fe1, const BSPface& face);
    void triangulateFace(uint64_t face_ind);
    void computeBaricenter(const vector<uint32_t>& vrts);
    inline uint64_t triFace_oppEdge(const BSPface& face, uint32_t v);
    uint64_t triFace_shareEdge(const BSPcell& cell, uint64_t face_ind,
                                           uint64_t vOppEdge_ind);
    bool cell_is_tetrahedrizable_from_v(const BSPcell& cell, uint32_t v);
    void makeTetrahedra();
};

/// <summary>
/// Main function - Create a polyhedral mesh out of the input
/// Input may be made of either one or two models to be combined into a boolean composition
/// </summary>
/// <param name="fileA_name">Name of first model</param>
/// <param name="coords_A">Serialized coordinates of first model vertices</param>
/// <param name="npts_A">Number of first model vertices</param>
/// <param name="tri_idx_A">Serialized indexes of first model triangles</param>
/// <param name="ntri_A">Number of first model triangles</param>
/// <param name="fileB_name">Name of second model</param>
/// <param name="coords_B">Serialized coordinates of second model vertices</param>
/// <param name="npts_B">Number of second model vertices</param>
/// <param name="tri_idx_B">Serialized indexes of second model triangles</param>
/// <param name="ntri_B">Number of second model triangles</param>
/// <param name="bool_opcode">Boolean operation (0 = no op, U = union, D = difference, I = intersection</param>
/// <param name="free_mem">Input models memory is released</param>
/// <param name="verbose">Print useful info during the process</param>
/// <param name="logging">Append info to mesh_generator.log</param>
/// <returns>The resulting BSPcomplex structure</returns>
BSPcomplex* makePolyhedralMesh(
    const char* fileA_name, double* coords_A, uint32_t npts_A, uint32_t* tri_idx_A, uint32_t ntri_A,
    const char* fileB_name = NULL, double* coords_B = NULL, uint32_t npts_B = 0, uint32_t* tri_idx_B = NULL, uint32_t ntri_B = 0,
    char bool_opcode = '0',
    bool free_mem = false,
    bool verbose = false, bool logging = false);

#endif /* BSP_h */
