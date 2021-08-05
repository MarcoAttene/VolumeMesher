#include <iostream>
#include <fstream>
#include <algorithm>
#include "BSP.h"

#include "graph_cut/GCoptimization.h"
#include "graph_cut/LinkedBlockList.cpp"
#include "graph_cut/GCoptimization.cpp"

// Label = 1 means INTERNAL
// Label = 0 means EXTERNAL

// Take first coplanar constraint associated to this face
// and return TRUE if the cell vertices are 'below' such constraint.
bool isFirstConnCellBelowFace(BSPface& f, BSPcomplex* cpx)
{
    const uint32_t* cid = cpx->constraints_vrts.data() + f.coplanar_constraints[0] * 3;
    const genericPoint* pv1 = cpx->vertices[cid[0]];
    const genericPoint* pv2 = cpx->vertices[cid[1]];
    const genericPoint* pv3 = cpx->vertices[cid[2]];

    BSPcell& cell = cpx->cells[f.conn_cells[0]];
    uint64_t num_cellEdges = UINT64_MAX;
    uint32_t num_cellVrts = cpx->count_cellVertices(cell, &num_cellEdges);
    vector<uint32_t> cell_vrts(num_cellVrts, UINT32_MAX);
    cpx->list_cellVertices(cell, num_cellEdges, cell_vrts);

    for (uint64_t ei : f.edges) cpx->vrts_visit[cpx->edges[ei].vertices[0]] = cpx->vrts_visit[cpx->edges[ei].vertices[1]] = 1;

    for (uint32_t vi : cell_vrts) if (!cpx->vrts_visit[vi])
    {
        const genericPoint* cv = cpx->vertices[vi];
        const int o = genericPoint::orient3D(*cv, *pv1, *pv2, *pv3);
        if (o)
        {
            for (uint64_t ei : f.edges) cpx->vrts_visit[cpx->edges[ei].vertices[0]] = cpx->vrts_visit[cpx->edges[ei].vertices[1]] = 0;
            return (o > 0);
        }
    }

    ip_error("Degenerate cell\n");
    return false;
}


#define DETERMINANT3X3(a11, a12, a13, a21, a22, a23, a31, a32, a33) ((a11)*((a22)*(a33) - (a23)*(a32)) - (a12)*((a21)*(a33) - (a23)*(a31)) + (a13)*((a21)*(a32) - (a22)*(a31)))

double approxFaceArea(BSPface& face, BSPcomplex* cpx, const std::vector<double>& approxCoords)
{
    std::vector<uint32_t> vs(face.edges.size(), 0);
    cpx->list_faceVertices(face, vs);

    const double* acp = approxCoords.data();

    const double* tv0, * tv1, * tv2;
    double a = 0.0;
    tv0 = acp + vs[0] * 3;
    for (size_t i = 2; i < vs.size(); i++)
    {
        tv1 = acp + vs[i - 1] * 3;
        tv2 = acp + vs[i    ] * 3;
        a += DETERMINANT3X3(tv0[0], tv0[1], tv0[2], tv1[0], tv1[1], tv1[2], tv2[0], tv2[1], tv2[2]);
    }

    return fabs(a);
}


// Returns TRUE if face is part of the skin according to skin_colour
inline bool isSkinFace(const BSPface& face, uint32_t skin_colour)
{
    return face.colour & skin_colour;
}

inline void setInternalCell(BSPcell& c, uint32_t internal_label)
{
    c.place |= internal_label;
}

inline void setExternalCell(BSPcell& c)
{
    c.place = EXTERNAL;
}


void BSPcomplex::markInternalCells(uint32_t skin_colour, uint32_t internal_label, const std::vector<double>& face_areas)
{
    // Allocate dual graph: num cells + 1 to account for the external "ghost" cell
    GCoptimizationGeneralGraph gc((GCoptimization::SiteID)cells.size() + 1, 2);

    // gc is the dual graph of the cell complex
    // - a node in gc corresponds to a cell in the complex
    // - an arc in gc exists if two cells share a WHITE face

    // In gc, cells that share a white face are connected by an arc weighted on face area
    // The 'data cost' associated to each cell is:
    //  - total area of BLACK faces 'consistently oriented' with the cell, if label is EXTERNAL
    //  - 0, if label is INTERNAL
    //
    // The 'smooth cost' associated to each WHITE face is:
    //  - the face area, if label_A and label_B are different
    //  - 0 otherwise
    //
    // Note: a BLACK face is considered to be consistently oriented with one of its incident cells
    // if the cell is 'below' the first of its coplanar constraints.

    // evs == 1 if edge is on boundary of skin
    std::vector<uint8_t> evs(edges.size(), 0);
    for (BSPface& f : faces) if (isSkinFace(f, skin_colour))
      for (uint64_t eid : f.edges) if (evs[eid] < 2) evs[eid]++;

    // vvs == 1 if vertex is on boundary of skin
    std::vector<uint8_t> vvs(vertices.size(), 0);
    for (size_t i = 0; i < edges.size(); i++) if (evs[i] == 1)
    {
        const BSPedge& e = edges[i];
        vvs[e.vertices[0]] = vvs[e.vertices[1]] = 1;
    }

    std::vector<double> cell_costs_external(cells.size() + 1, 0.0);
    std::vector<double> cell_costs_internal(cells.size() + 1, 0.0);

    for (size_t i = 0; i < faces.size(); i++)
    {
        BSPface& f = faces[i];
        const uint64_t cell1 = f.conn_cells[0];
        const uint64_t cell2 = (f.conn_cells[1] == UINT64_MAX) ? (cells.size()) : (f.conn_cells[1]);

        if (isSkinFace(f, skin_colour))
        {
            if (isFirstConnCellBelowFace(f, this))
            {
                cell_costs_external[cell1] += face_areas[i];
                cell_costs_internal[cell2] += face_areas[i];
            }
            else
            {
                cell_costs_external[cell2] += face_areas[i];
                cell_costs_internal[cell1] += face_areas[i];
            }
        }
    }

    for (size_t i = 0; i < faces.size(); i++)
    {
        BSPface& f = faces[i];
        const uint64_t cell1 = f.conn_cells[0];
        if (!isSkinFace(f, skin_colour))
        {
            // 'w' is an additional weight for arcs that promotes the cut of arcs corresp. to faces having
            // all their vertices on the boundary of the input surface (i.e. promote hole filling)
            double w = 0.1;
            for (uint64_t eid : f.edges) if (vvs[edges[eid].vertices[0]] == 0 || vvs[edges[eid].vertices[1]] == 0) { w = 1.0; break; }
            const uint64_t cell2 = (f.conn_cells[1] == UINT64_MAX) ? (cells.size()) : (f.conn_cells[1]);
            gc.setNeighbors((GCoptimization::SiteID)cell1, (GCoptimization::SiteID)cell2, face_areas[i]*w);
        }
    }

    const double int_weight = 0.1; // Internal cell penalization less than external to avoid artifacts at intersections
    for (size_t i = 0; i < cells.size(); i++) gc.setDataCost((GCoptimization::SiteID)i, 0, cell_costs_external[i]);
    for (size_t i = 0; i < cells.size(); i++) gc.setDataCost((GCoptimization::SiteID)i, 1, cell_costs_internal[i] * int_weight);
    gc.setDataCost((GCoptimization::SiteID)cells.size(), 1, 1.0); // Ghost cell must be external

    // Run graph cut algorithm
    // I.e., label all the cells so that the total data cost + smooth cost is minimized
    gc.swap();

    for (size_t i = 0; i < cells.size(); i++)
        if (gc.whatLabel((GCoptimization::SiteID)i)) setInternalCell(cells[i], internal_label);
}

void BSPcomplex::constraintsSurface_complexPartition(bool two_files)
{
    // Make all cells external
    for (size_t i = 0; i < cells.size(); i++) setExternalCell(cells[i]);

    // Clear vrts_visit for use in isFirstConnCellBelowFace()
    for (size_t i = 0; i < vertices.size(); i++) vrts_visit[i] = 0;

    // Precalculate approximate vertex coordinates for use in approxFaceArea()
    std::vector<double> approxCoords(vertices.size() * 3);
    for (size_t i = 0; i < vertices.size(); i++)
        vertices[i]->getApproxXYZCoordinates(approxCoords[i * 3], approxCoords[i * 3 + 1], approxCoords[i * 3 + 2]);

    // Precalculate approximate face areas for use in markInternalCells()
    std::vector<double> face_areas(faces.size(), 0.0);
    double tot_face_area = 0.0;
    for (size_t i = 0; i < faces.size(); i++)
    {
        face_areas[i] = approxFaceArea(faces[i], this, approxCoords);
        tot_face_area += face_areas[i];
    }

    // Normalize all areas to avoid overflows in graphcut
    for (size_t i = 0; i < faces.size(); i++) face_areas[i] /= tot_face_area;

    if (two_files)
    {
        markInternalCells(BLACK_A, INTERNAL_A, face_areas);
        markInternalCells(BLACK_B, INTERNAL_B, face_areas);
    }
    else markInternalCells(BLACK_A, INTERNAL_A, face_areas);
}

