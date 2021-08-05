#include "BSP.h"

void read_OFF_file(const char* filename,
    double** vertices_p, uint32_t* npts,
    uint32_t** tri_vertices_p, uint32_t* ntri, bool verbose) {

    FILE* file = fopen(filename, "r");
    if (file == NULL)
        ip_error("read_nodes_and_constraints: FATAL ERROR "
            "cannot open input file.\n");

    // Check OFF mark (1st line).
    char file_ext_read[3];
    char file_ext_target[] = { 'O','F','F' };
    if (fscanf(file, "%3c", file_ext_read) == 0)
        ip_error("read_nodes_and_constraints: FATAL ERROR "
            "cannot read 1st line of input file\n");

    for (uint32_t i = 0; i < 3; i++)
        if (file_ext_read[i] != file_ext_target[i])
            ip_error("read_nodes_and_constraints: FATAL ERROR "
                "1st line of input file is different from OFF\n");

    // Reading number of points and triangles.
    if (fscanf(file, " %d %d %*d ", npts, ntri) == 0)
        ip_error("read_nodes_and_constraints: FATAL ERROR 2st line of "
            "input file do not contanins point and triangles numbers.\n");

    if (verbose) printf("file %s contains %d vertices and %d constraints (triangles)\n",
        filename, *npts, *ntri);

    // Reading points coordinates.
    *vertices_p = (double*)malloc(sizeof(double) * 3 * (*npts));
    *tri_vertices_p = (uint32_t*)malloc(sizeof(uint32_t) * 3 * (*ntri));

    for (uint32_t i = 0; i < (*npts); i++) {
        if (fscanf(file, " %lf %lf %lf ",
            (*vertices_p) + (i * 3), (*vertices_p) + (i * 3 + 1), (*vertices_p) + (i * 3 + 2)) == 0)
            ip_error("error reading input file\n");
    }

    uint32_t nv;
    for (uint32_t i = 0; i < (*ntri); i++) {
        if (fscanf(file, " %u %u %u %u ", &nv,
            (*tri_vertices_p) + (i * 3), (*tri_vertices_p) + (i * 3 + 1), (*tri_vertices_p) + (i * 3 + 2)) == 0)
            ip_error("error reading input file\n");
        if (nv != 3) ip_error("Non-triangular faces not supported\n");
    }
    fclose(file);
}

/// <summary>
/// Main function
/// </summary>
/// <param name="argc"></param>
/// <param name="argv"></param>
/// <returns></returns>
int main(int argc, char** argv)
{
#ifndef DEBUG
    if (argc < 2) {
        printf("\nUsage: mesh_generator [-v | -l | -s | -b | -t] inputfile_A.off [bool_opcode inputfile_B.off]\n\n"
            "Defines the volume enclosed by the input OFF file(s) and saves a volume mesh to 'volume.msh'\n\n"
            "Command line arguments:\n"
            "-v = verbose mode\n"
            "-l = logging mode (appends a line to mesh_generator.log)\n"
            "-s = save the mesh bounding surface to 'skin.off'\n"
            "-b = save the subdivided constraints to 'black_faces.off'\n"
            "-t = triangulate/tetrahedrize output\n"
            "bool_opcode: {U, I, D}\n"
            "  U -> union (AuB),\n"
            "  I -> intersection (A^B),\n"
            "  D -> difference (A\\B)\n\n"
            "Example:\n"
            "mesh_generator ant.off U pig.off\n");
        return 0;
    }

    bool triangulate = false;
    bool verbose = false;
    bool logging = false;
    bool surfmesh = false;
    bool blackfaces = false;
    char* fileA_name = NULL;
    char* fileB_name = NULL;
    char bool_opcode = '0';

    for (int i = 1; i < argc; i++)
    {
        if (argv[i][0] == '-')
        {
            if (argv[i][1] == 'v') verbose = true;
            else if (argv[i][1] == 't') triangulate = true;
            else if (argv[i][1] == 'l') logging = true;
            else if (argv[i][1] == 'b') blackfaces = true;
            else if (argv[i][1] == 's') surfmesh = true;
            else ip_error("Unknown option\n");
        }
        else if (fileA_name == NULL) fileA_name = argv[i];
        else if (bool_opcode == '0') bool_opcode = argv[i][0];
        else if (fileB_name == NULL) fileB_name = argv[i];
        else ip_error("Too many args passed\n");
    }

    bool two_input = (bool_opcode != '0');

    if (verbose)
    {
        if (fileB_name == NULL) {
            printf("\nResolve auto-intersections and/or repair.\n\n");
            printf("Loading %s\n", fileA_name);
        }
        else {
            printf("\nBoolean operator: ");
            if (bool_opcode == 'U') printf(" union.\n\n");
            else if (bool_opcode == 'I') printf(" intersection.\n\n");
            else if (bool_opcode == 'D') printf(" difference.\n\n");
            else { printf("INVALID\n\n"); return 0; }
            printf("Loading %s and %s.\n\n", fileA_name, fileB_name);
        }
    }
#else
    bool triangulate = false;
    bool verbose = true;
    bool logging = false;
    bool surfmesh = false;
    bool blackfaces = false;
    char* fileA_name = "whatsoever";
    char* fileB_name = NULL;
    char bool_opcode = '0';
    bool two_input = (bool_opcode != '0');
#endif

    double* coords_A, * coords_B = NULL;
    uint32_t ncoords_A, ncoords_B;
    uint32_t* tri_idx_A, * tri_idx_B;
    uint32_t ntriidx_A, ntriidx_B;

    read_OFF_file(fileA_name, &coords_A, &ncoords_A, &tri_idx_A, &ntriidx_A, verbose);
    if (two_input) read_OFF_file(fileB_name, &coords_B, &ncoords_B, &tri_idx_B, &ntriidx_B, verbose);

    BSPcomplex* complex = makePolyhedralMesh(
        fileA_name, coords_A, ncoords_A, tri_idx_A, ntriidx_A,
        fileB_name, coords_B, ncoords_B, tri_idx_B, ntriidx_B,
        bool_opcode, true, verbose, logging
        );

    printf("Writing output file...\n");
    if (blackfaces) complex->saveBlackFaces("black_faces.off", triangulate);
    else if (surfmesh) complex->saveSkin("skin.off", bool_opcode, triangulate);
    else complex->saveMesh((triangulate)?("volume.tet"):("volume.msh"), bool_opcode, triangulate);
    printf("Done.\n");

    return 0;
}
