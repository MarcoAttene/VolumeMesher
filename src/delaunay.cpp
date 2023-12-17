#include "delaunay.h"
#include <float.h>
#include <iostream>
#include <fstream>

TetMesh::TetMesh() : 
    vertices(NULL), num_vertices(0), 
    tet_node(NULL), tet_neigh(NULL), tet_subdet(NULL), tet_num(0), tet_size(0), tet_num_vertices(0), mark_tetrahedra(NULL), 
    Del_size_tmp(1024), Del_num_tmp(0), Del_tmp(NULL), Del_size_deleted(1024), Del_num_deleted(0), Del_deleted(NULL), Del_buffer(NULL)
{
}

TetMesh::~TetMesh() {
    free(vertices);
    free(tet_node);
    free(tet_neigh);
    free(tet_subdet);
    free(mark_tetrahedra);
}

void TetMesh::reserve(uint64_t ntet){
  ntet+=tet_num;
  tet_neigh = (uint64_t *)realloc(tet_neigh, ntet*4*sizeof(uint64_t));
  tet_subdet = (double *)realloc(tet_subdet, ntet*4*sizeof(double));
  tet_node = (uint32_t *)realloc(tet_node, ntet*4*sizeof(uint32_t));
  tet_size = ntet;
  mark_tetrahedra = (uint32_t *)realloc(mark_tetrahedra, ntet*sizeof(uint32_t));    // Mod.1
}

void TetMesh::init(){
  uint32_t n = num_vertices;

  // Calculate static filters
  double c;
  double minx = DBL_MAX, miny = DBL_MAX, minz = DBL_MAX;
  double maxx = -DBL_MAX, maxy = -DBL_MAX, maxz = -DBL_MAX;

  for (uint32_t i = 0; i < n; i++)
  {
      const double* cp = vertices[i].coord;
      c = *(cp++);
      if (c < minx) minx = c;
      if (c > maxx) maxx = c;
      c = *(cp++);
      if (c < miny) miny = c;
      if (c > maxy) maxy = c;
      c = *cp;
      if (c < minz) minz = c;
      if (c > maxz) maxz = c;
  }
  maxx = fabs(maxx); maxy = fabs(maxy); maxz = fabs(maxz);
  minx = fabs(minx); miny = fabs(miny); minz = fabs(minz);
  if (minx > maxx) maxx = minx;
  if (miny > maxy) maxy = miny;
  if (minz > maxz) maxz = minz;
  if (maxx > maxz) std::swap(maxx, maxz);
  if (maxy > maxz) std::swap(maxy, maxz); 
  else if (maxy < maxx) std::swap(maxx, maxy); 

  o3d_static_filter = 5.1107127829973299e-15 * maxx * maxy * maxz;
  isp_static_filter = 1.2466136531027298e-13 * maxx * maxy * maxz * (maxz * maxz);


  // Find non-coplanar vertices (we assume that no coincident vertices exist)
  double ori=0.0;
  uint32_t i=0, j=1, k=2, l=3;

  for (; ori==0.0 && k<n-1; k++)
      for (l = k + 1; ori == 0.0 && l < n; l++)
          ori = orient3d(vertices[i].coord,
              vertices[j].coord,
              vertices[k].coord,
              vertices[l].coord);
  l--; k--;

  if(ori==0.0)
    ip_error("Input vertices do not define a volume.\n");

  std::swap(vertices[k], vertices[2]); k=2;
  std::swap(vertices[l], vertices[3]); l=3;

  if(ori<0.0) std::swap(i, j); // Tets must have positive volume

  const uint32_t base_tet[] = { l, k, j, i, l, j, k, INFINITE_VERTEX, l, k, i, INFINITE_VERTEX, l, i, j, INFINITE_VERTEX, k, j, i, INFINITE_VERTEX };
  const uint64_t base_neigh[] = { 19, 15, 11, 7, 18, 10, 13, 3, 17, 14, 5, 2, 16, 6, 9, 1, 12, 8, 4, 0 };

  reserve(num_vertices * 10);
  std::memcpy(tet_node, base_tet, 20 * sizeof(uint32_t));
  std::memcpy(tet_neigh, base_neigh, 20 * sizeof(uint64_t));

  compute_subDet(0);
  compute_subDet(4);
  compute_subDet(8);
  compute_subDet(12);
  compute_subDet(16);

  tet_num = 5;
  tet_num_vertices = 4; // not counting the ghost vertex

  // mark the tetrahedra as non-visited
  for(i=0; i<tet_num; i++) mark_tetrahedra[i]=0;

  // set the vertex-(one_of_the)incident-tetrahedron relation
  for(i=0; i<tet_num_vertices; i++) vertices[i].inc_tet = 0;
}



void TetMesh::tetrahedrize() {
  init();
  allocTmpStruct(num_vertices);

  uint64_t ct = 0;
  for (uint32_t i=4; i< num_vertices; i++)
  {
      ct = searchTetrahedron(ct, i);

      deleteInSphereTets(ct, i);
      tetrahedrizeHole(&ct);

      uint64_t Tet2update = ct;
      if(tet_node[Tet2update+3]==INFINITE_VERTEX)
        Tet2update=tet_neigh[Tet2update+3];
      
      vertices[i].inc_tet = Tet2update>>2;
  }

  tet_num_vertices = num_vertices;

  removeDelTets();
  releaseTmpStruct();

  mark_tetrahedra = (uint32_t *)realloc(mark_tetrahedra,(tet_num)*sizeof(uint32_t));

  for(uint64_t i=0; i<tet_num; i++) mark_tetrahedra[i]=0;
}


void TetMesh::saveTET(const char* filename)
{
    ofstream f(filename);

    if (!f) ip_error("\nTetMesh::saveTET: FATAL ERROR cannot open the file.\n");

    f << num_vertices << " vertices\n";

    uint64_t ngnt = 0;
    for (uint64_t i = 0; i < tet_num; i++) if (tet_node[i * 4 + 3] != INFINITE_VERTEX) ngnt++;

    f << ngnt << " tets\n";
    for (uint32_t i = 0; i < num_vertices; i++)
        f << vertices[i].coord[0] << " " << vertices[i].coord[1] << " " << vertices[i].coord[2] << "\n";
    for (uint64_t i = 0; i < tet_num; i++) if (tet_node[i*4 + 3] != INFINITE_VERTEX)
        f << "4 " << tet_node[i * 4] << " " << tet_node[i * 4 + 1] << " " << tet_node[i * 4 + 2] << " " << tet_node[i * 4 + 3] << "\n";
    
    f.close();
}

void TetMesh::removeDelTets(){
  uint64_t j;
  for (uint64_t i=0; i<Del_num_deleted; i++)
  {
    uint64_t to_delete = Del_deleted[i];

    tet_num--;
    uint64_t lastTet = tet_num*4;

    if(tet_subdet[lastTet + 3]==-1.0)
    {
      for (j=i; j<Del_num_deleted; j++)
        if(Del_deleted[j]==lastTet) break;

      Del_deleted[j] = Del_deleted[i];
    }
    else {
      for (j=0; j<4; j++)
      {
        tet_node[to_delete+j] = tet_node[lastTet+j];
        tet_subdet[to_delete+j] = tet_subdet[lastTet+j];

        uint64_t neigh = tet_neigh[lastTet+j];
        tet_neigh[to_delete+j] = neigh;
        tet_neigh[neigh] = to_delete+j;

        if(tet_node[lastTet + j] != INFINITE_VERTEX &&
           vertices[ tet_node[lastTet+j] ].inc_tet == lastTet>>2)
              vertices[ tet_node[lastTet+j] ].inc_tet = to_delete>>2;
      }
    }
  }

  Del_num_deleted = 0;
}


uint64_t TetMesh::searchTetrahedron(uint64_t tet, const uint32_t v_id)
{
    if (tet_node[tet + 3] == INFINITE_VERTEX)
        tet = getNeighbor(tet, 3);

    const double* vc = vertices[v_id].coord;

    for (uint64_t f0 = 4;;) {
        const uint32_t* Node = tet_node + tet;
        if (Node[3] == INFINITE_VERTEX) return tet;

        const uint64_t* Neigh = tet_neigh + tet;
        uint64_t i = 0;
        for (; i < 4; i++)
        {
            const double* a = vertices[Node[(i + 1) & 3]].coord;
            const double* b = vertices[Node[(i & 2) ^ 3]].coord;
            const double* c = vertices[Node[(i + 3) & 2]].coord;

            if (i != f0 && orient3d(a, b, c, vc) < 0.0) {
                tet = getIthNeighbor(Neigh, i);
                f0 = Neigh[i] & 3;
                break;
            }
        }

        if (i == 4) return tet;
    }
}


static inline double symbolicPerturbation(uint32_t indices[5], 
	const double* i, const double* j, const double* k, const double* l, const double* m)
{
  const double *pt[5] = {i,j,k,l,m};

  int swaps = 0;
  int n = 5;
  int count;
  do {
    count = 0;
    n = n - 1;
    for (int i = 0; i < n; i++) {
      if (indices[i] > indices[i+1]) {
        std::swap(pt[i], pt[i+1]);
        std::swap(indices[i], indices[i+1]);
        count++;
      }
    }
    swaps += count;
  } while (count);

  double oriA = orient3d(pt[1], pt[2], pt[3], pt[4]);
  if (oriA != 0.0) {
    if ((swaps % 2) != 0) oriA = -oriA;
    return oriA;
  }

  double oriB = -orient3d(pt[0], pt[2], pt[3], pt[4]);
  if ((swaps % 2) != 0) oriB = -oriB;
  return oriB;
}


double TetMesh::vertexInTetSphere(uint64_t tet, uint32_t v_id){

  const uint32_t* Node = tet_node + tet;
  const double* SubDet = tet_subdet + tet;

  const double* a = vertices[Node[0]].coord;
  const double* e = vertices[v_id].coord;

  const double aex = e[0] - a[0];
  const double aey = e[1] - a[1];
  const double aez = e[2] - a[2];

  if(Node[3]==INFINITE_VERTEX){
    double det = aex*SubDet[0]+aey*SubDet[1]+aez*SubDet[2];
    if (fabs(det) > o3d_static_filter) return det;

    const double* b = vertices[Node[1]].coord;
    const double* c = vertices[Node[2]].coord;

    det = orient3d(a,b,c,e);
    if(det!=0.0) return det;

    const uint32_t oppositeNode = tet_node[tet_neigh[tet+3]];
    const double* oppositeVertex = vertices[oppositeNode].coord;
    det = -insphere(a,b,c,oppositeVertex,e);

    if (det == 0.0) {
      uint32_t nn[5] = {Node[0],Node[1],Node[2],oppositeNode,v_id};
      det = -symbolicPerturbation (nn, a,b,c,oppositeVertex,e);
      if (det == 0.0) ip_error("Symbolic perturbation failed\n");
    }
    return det;
  }

  const double aer = aex*aex + aey*aey + aez*aez;
  double det = aex*SubDet[0]-aey*SubDet[1]+aez*SubDet[2]-aer*SubDet[3];
  if (fabs(det) > isp_static_filter) return det;

  const double* b = vertices[Node[1]].coord;
  const double* c = vertices[Node[2]].coord;
  const double* d = vertices[Node[3]].coord;

  det = insphere(a,b,c,d,e);

  if (det == 0.0) {
    uint32_t nn[5] = {Node[0],Node[1],Node[2],Node[3],v_id};
    det = symbolicPerturbation (nn, a,b,c,d,e);
    if (det == 0.0) ip_error("Symbolic perturbation failed! This should never happen :-(\n");
  }
  return det;
}

void TetMesh::bnd_push(uint32_t v_id,
              uint32_t node1, uint32_t node2,
              uint32_t node3, uint64_t bnd)
{
  uint64_t n = Del_num_tmp;
  Del_tmp[n].node[0] = v_id;
  Del_tmp[n].node[1] = node1;
  Del_tmp[n].node[2] = node2;
  Del_tmp[n].node[3] = node3;
  Del_tmp[n].bnd = bnd;
  Del_num_tmp++;
}


void TetMesh::deleteInSphereTets(uint64_t tet, const uint32_t v_id)
{
  uint64_t start;
  Del_deleted[Del_num_deleted++] = tet;
  tet_subdet[tet + 3] = -1.0;

  uint64_t first_deleted_pos = Del_num_deleted-1;

  for(start=Del_num_deleted-1; start < Del_num_deleted; start++){
    uint64_t tet = Del_deleted[start];
    uint64_t* Neigh = tet_neigh + tet;
    uint32_t* Node = tet_node + tet;

    if(Del_num_tmp + 4 > Del_size_tmp){
        Del_tmp = (DelTmp *)realloc(Del_tmp, 2*Del_num_tmp*sizeof(Del_tmp[0]));
      Del_size_tmp = 2*Del_num_tmp;
    }

    if(Del_num_deleted + 4 > Del_size_deleted){
      Del_deleted = (uint64_t *)realloc(Del_deleted, 2*Del_num_deleted*sizeof(uint64_t));
      Del_size_deleted = 2*Del_num_deleted;
    }

    uint64_t neigh = getIthNeighbor(Neigh, 0);
    if(tet_subdet[neigh + 3]!=-1.0){
      if(vertexInTetSphere(neigh, v_id)<0.0){
        bnd_push(v_id, Node[1], Node[2], Node[3], Neigh[0]);
      }
      else{
        Del_deleted[Del_num_deleted++] = neigh;
        tet_subdet[neigh + 3] = -1.0;
      }
    }

    neigh = getIthNeighbor(Neigh, 1);
    if(tet_subdet[neigh + 3]!=-1.0){
      if(vertexInTetSphere(neigh, v_id)<0.0){
        bnd_push(v_id, Node[2], Node[0], Node[3], Neigh[1]);
      }
      else{
        Del_deleted[Del_num_deleted++] = neigh;
        tet_subdet[neigh + 3] = -1.0;
      }
    }

    neigh = getIthNeighbor(Neigh, 2);
    if(tet_subdet[neigh + 3]!=-1.0){
      if(vertexInTetSphere(neigh, v_id)<0.0){
        bnd_push(v_id, Node[0], Node[1], Node[3], Neigh[2]);
      }
      else{
        Del_deleted[Del_num_deleted++] = neigh;
        tet_subdet[neigh + 3] = -1.0;
      }
    }

    neigh = getIthNeighbor(Neigh, 3);
    if(tet_subdet[neigh + 3]!=-1.0){
      if(vertexInTetSphere(neigh, v_id)<0.0){
        if(Node[1]<Node[2])
          bnd_push(v_id, Node[0], Node[2], Node[1], Neigh[3]);
        else
          bnd_push(v_id, Node[1], Node[0], Node[2], Neigh[3]);
      }
      else{
        Del_deleted[Del_num_deleted++] = neigh;
        tet_subdet[neigh + 3] = -1.0;
      }
    }
  }
}


void TetMesh::tetrahedrizeHole(uint64_t* tet){
  uint64_t clength = Del_num_deleted;
  uint64_t blength = Del_num_tmp;

  if(blength > clength){
    if(blength>Del_size_deleted){
        Del_deleted = (uint64_t *)realloc(Del_deleted, 2*blength*sizeof(uint64_t));
      Del_size_deleted = 2*blength;
    }

    uint64_t i;
    for (i=clength; i<blength; i++){
      Del_deleted[i] = tet_num<<2;
      tet_num++;
    }

    clength = blength;

    if(tet_num > tet_size) reserve(tet_num);
  }

  uint64_t start = clength - blength;

  for (uint64_t i=0; i<blength; i++)
  {
    const uint64_t tet = Del_deleted[i + start];
    uint32_t* Node = tet_node + tet;

    Node[0] = Del_tmp[i].node[0];
    Node[1] = Del_tmp[i].node[1];
    Node[2] = Del_tmp[i].node[2];
    Node[3] = Del_tmp[i].node[3];

    uint64_t bnd = Del_tmp[i].bnd;
    tet_neigh[tet] = bnd;
    tet_neigh[bnd] = tet;
    Del_tmp[i].bnd = tet;

    compute_subDet(tet);

    // set all the incidence vertex-tetrahedron relation of the vertices
    // of the current tetrahedron to the tetrahedron itself
    if(tet_node[tet+3]!=INFINITE_VERTEX)            // Mod.1
      for(uint32_t j=0; j<4; j++)
          vertices[tet_node[tet + j]].inc_tet = tet>>2;
  }

  uint64_t tlength = 0;
  const uint64_t middle = blength * 3 / 2;

  uint64_t* Tmp = (uint64_t*)Del_tmp;
  const unsigned index[4] = { 2,3,1,2 };

  uint64_t i;
  for (i = 0; i < blength; i++)
  {
      uint64_t tet = Del_deleted[start + i];
      const uint32_t* const Node = tet_node + tet;

      uint64_t j;
      for (j = 0; j < 3; j++)
      {
          uint64_t key = ((uint64_t)Node[index[j]] << 32) + Node[index[j + 1]];
          tet++;

          uint64_t k;
          for (k = 0; k < tlength; k++)
          {
              if (Tmp[k] == key)
                  break;
          }

          if (k == tlength) {
              Tmp[tlength] = (key >> 32) + (key << 32);
              Tmp[middle + tlength] = tet;
              tlength++;
          }
          else {
              uint64_t pairValue = Tmp[middle + k];
              tet_neigh[tet] = pairValue;
              tet_neigh[pairValue] = tet;
              tlength--;
              if (k < tlength) {
                  Tmp[k] = Tmp[tlength];
                  Tmp[middle + k] = Tmp[middle + tlength];
              }
          }
      }
  }

  Del_num_tmp = 0;
  Del_num_deleted = start;


  *tet = Del_deleted[start];
}


void TetMesh::compute_subDet(const uint64_t tet)
{
  double ab[4],ac[4];

  uint32_t* Node = tet_node + tet;
  double* SubDet = tet_subdet + tet;

  double* a = vertices[Node[0]].coord;
  double* b = vertices[Node[1]].coord;
  double* c = vertices[Node[2]].coord;

  if(Node[3]!=INFINITE_VERTEX){
    double* d = vertices[Node[3]].coord;
    double ad[4];

    unsigned i;
    for (i=0; i<3; i++)
    {
      ab[i]=b[i]-a[i]; 
      ac[i]=c[i]-a[i]; 
      ad[i]=d[i]-a[i]; 
    }

    ab[3] = ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2];
    ac[3] = ac[0]*ac[0] + ac[1]*ac[1] + ac[2]*ac[2];
    ad[3] = ad[0]*ad[0] + ad[1]*ad[1] + ad[2]*ad[2];

    double cd12 = ac[2]*ad[1] - ac[1]*ad[2];
    double db12 = ad[2]*ab[1] - ad[1]*ab[2];
    double bc12 = ab[2]*ac[1] - ab[1]*ac[2];

    double cd30 = ac[0]*ad[3] - ac[3]*ad[0];
    double db30 = ad[0]*ab[3] - ad[3]*ab[0];
    double bc30 = ab[0]*ac[3] - ab[3]*ac[0];

    SubDet[0] = ab[3]*cd12 + ac[3]*db12 + ad[3]*bc12;
    SubDet[1] = ab[2]*cd30 + ac[2]*db30 + ad[2]*bc30;
    SubDet[2] = ab[1]*cd30 + ac[1]*db30 + ad[1]*bc30;
    SubDet[3] = ab[0]*cd12 + ac[0]*db12 + ad[0]*bc12;
  }
  else {
    unsigned i;
    for (i=0; i<3; i++)
    {
      ab[i]=b[i]-a[i];
      ac[i]=c[i]-a[i];
    }

    SubDet[0] = ac[1]*ab[2] - ac[2]*ab[1];
    SubDet[1] = ac[2]*ab[0] - ac[0]*ab[2];
    SubDet[2] = ac[0]*ab[1] - ac[1]*ab[0];
    SubDet[3] = SubDet[0]*SubDet[0] + SubDet[1]*SubDet[1] + SubDet[2]*SubDet[2];
  }
}




//---------------------------------
// incident tetrahedra in a vertex
//---------------------------------

// - Recursive function -
//  Input: tetrahedra (tet) index: tet_ind,
//         the index of the vertex in which incidences are searched: central_vertex_ind,
//         pointer to mesh,
//         pointer to a tetraheron index type: tot_incTet.
// Output: by using tot_incTet increments of one unity the number of incident tetrahedra each
//         time one of tetrahedra sharing a face with tet
//         - is INCIDENT in central_vertex_ind,
//         - is NOT a GHOST TETRHEDRON,
//         - is NOT been already VISITED (mark_tetrahedra=0).
//         If tot_incTet is incremented the relative tetrahedra is passed to this function...
void count_incTet(uint64_t tet_ind, const uint32_t central_vertex_ind, TetMesh* mesh, uint64_t* tot_incTet){   // Mod.1
    for(uint32_t i=0; i<4; i++)
    {
        if(mesh->tet_node[4*tet_ind+i] != central_vertex_ind)
        {
            uint64_t neigh_tet_ind = mesh->tet_neigh[4*tet_ind+i]>>2;
            // ghost vertex is always in the last slot of the tetrahedron vertices.

            if(mesh->mark_tetrahedra[neigh_tet_ind]==1                        ||
               mesh->tet_node[4* neigh_tet_ind +3] == INFINITE_VERTEX   )
                continue;

            mesh->mark_tetrahedra[ neigh_tet_ind ]=1;
            (*tot_incTet)++;
            count_incTet(neigh_tet_ind, central_vertex_ind, mesh, tot_incTet);
        }
    }
}


// - Recursive function -
//  Input: tetrahedra (tet) index: tet_ind,
//         pointer to mesh,
//         pointer to a tetraheron index type: incTet,
//         pointer to a tetraheron index type: pos.
// Output: by using incTet returns the indices of the tetrahedra icident in a common vertex,
//         those tetrahedra have been counted by the function count_incTet and marked as 1.
//         Face neighbour of tet is added to incTet if its marker is equal to 1.
//         pos stores the number of elements of incTet.
void save_incTet(uint64_t tet_ind, TetMesh* mesh, uint64_t* incTet, uint64_t* pos){      // Mod.1
    for(uint32_t i=0; i<4; i++)
    {
        if( mesh->mark_tetrahedra[ mesh->tet_neigh[4*tet_ind+i]>>2 ] == 1 )
        {
            (*pos)++;
            uint64_t neigh_tet_ind = mesh->tet_neigh[4*tet_ind+i]>>2;
            incTet[*pos]=neigh_tet_ind;
            mesh->mark_tetrahedra[neigh_tet_ind]=0;

            save_incTet(neigh_tet_ind, mesh, incTet, pos);
        }
    }

}

//  Input: the index of the vertex in which incidences are searched: central_vertex_ind,
//         pointer to mesh,
//         pointer to a tetraheron index type: num_incTet.
// Output: by using num_incTet returns the number of non-ghost tetrahedra incident in central_vertex,
//         returns an array containing the indices of non-ghost tetrahedra incident in central_vertex.
uint64_t* TetMesh::incident_tetrahedra(const uint32_t central_vertex_ind, uint64_t* num_incTet)     // Mod.1
{
    uint64_t tet_ind = vertices[central_vertex_ind].inc_tet;

   // count the number of incident tetrahedra
    uint64_t l_incTet=1;
    mark_tetrahedra[tet_ind]=1;
    count_incTet(tet_ind, central_vertex_ind, this, &l_incTet);

    // collect the indices of incident tetrahedra
    uint64_t* incTet = (uint64_t*) malloc(sizeof(uint64_t)*l_incTet);
    incTet[0]=tet_ind;
    mark_tetrahedra[tet_ind]=0;
    uint64_t pos=0;
    save_incTet(tet_ind, this, incTet, &pos);

    * num_incTet = l_incTet;
    return incTet;
}


//---------------------------------
// incident tetrahedra at an edge
//---------------------------------


//  Input: pointer to the mesh,
//         2 vertices of an edge: (edge_ends[0],edge_ends[1]),
//         index of a tetrhedron: tet0_ind,
//         pointer of tetrahedra index type: length_related_tet.
// Output: by num_related_tet returns the number of tetrahedra that share with first_tet the edge,
//         returns the array containing the indices of those tet_
//  Note. The tetrahedron must have (edge_ends[0],edge_ends[1]) as its side.
uint64_t* TetMesh::ETrelation(const uint32_t* edge_ends, const uint64_t tet0_ind, uint64_t* length_related_tet) const
{
    uint64_t prec_tet_ind, curr_tet_ind;
    uint64_t move_dir, next_move_dir;
    uint32_t v_ind, other_tet_vrts_ind[2];
    uint64_t num_related_tet = 1; // counts first_tet

    uint64_t i=0;
    for(uint64_t j=0; j<4; j++){
        v_ind = tet_node[4*tet0_ind + j];
        if(v_ind != edge_ends[0] && v_ind != edge_ends[1] ){
            other_tet_vrts_ind[i++] = v_ind;
            move_dir = j;
        }
    }

    // move_dir is the face-ID of tet0->tet1, and vrt-ID opposite to that face [wrt tet0]
    // Move from tet0 to tet1
    prec_tet_ind = tet0_ind;
    curr_tet_ind = tet_neigh[ 4*tet0_ind + move_dir] >> 2; // tet1

    // other_tet_vrts_ind[1] is the index of the vertex opposite to the face between tet0 and tet1,
    // other_tet_vrts_ind[0] is the index of the common vertex (not on the edge) between tet0 and tet1,
    // other_tet_vrts_ind[0] is the index of the vertex opposite to the face between tet1 and tet2
    for(uint64_t j=0; j<4; j++)
        if(tet_node[4*curr_tet_ind + j] == other_tet_vrts_ind[0]){
            next_move_dir = j;     //is the face-ID of tet1 -> tet2, and vrt-ID opposite to that face [wrt tet1]
            break;
        }

    // at begininng we have: prec_tet=tet0, curr_tet=tet1, next_tet=tet2.
    while(curr_tet_ind != tet0_ind){

        // Conut curr_tet.
        num_related_tet++;

        // Find the vrt-ID of the common vertex (not on the edge) between curr_tet and next_tet [wrt curr_tet].
        uint64_t tmp = tet_neigh[ 4*prec_tet_ind + move_dir ] & 3;
        // Find the index of the vertex of curr_tet opposite to the face between curr_tet and prec_tet.
        v_ind = tet_node[ 4*curr_tet_ind + tmp];

        move_dir = next_move_dir; // to move from curr_tet to next_tet.
        prec_tet_ind = curr_tet_ind; // new prec_tet is curr_tet
        // Move from curr_tet to next_tet.
        curr_tet_ind = tet_neigh[ 4* curr_tet_ind + move_dir ] >> 2; // new curr_tet is nex_tet

        // Find the face-ID of the common face between new curr_tet and new prec_tet [wrt new curr_tet].
        for(uint64_t j=0; j<4; j++)
            if(tet_node[ 4* curr_tet_ind + j ] == v_ind)
                next_move_dir = j;
    }

    uint64_t* related_tet = (uint64_t*) malloc(sizeof(uint64_t) * num_related_tet);

    for(uint64_t i=0; i<num_related_tet; i++){

        // add curr_tet
        related_tet[i] = curr_tet_ind;

        // Move from curr_tet to next_tet
        uint64_t tmp = tet_neigh[ 4*prec_tet_ind + move_dir ] & 3;
        v_ind = tet_node[ 4*curr_tet_ind + tmp];
        move_dir = next_move_dir;

        prec_tet_ind = curr_tet_ind;
        curr_tet_ind = tet_neigh[ 4* curr_tet_ind + move_dir ] >> 2;

        for(uint64_t j=0; j<4; j++)
            if(tet_node[ 4* curr_tet_ind + j ] == v_ind)
                next_move_dir = j;

    }

    *length_related_tet = num_related_tet;
    return related_tet;
}







void TetMesh::allocTmpStruct(uint32_t num_vertices)
{
    Del_tmp = (DelTmp*) malloc(Del_size_tmp * sizeof(Del_tmp[0]));
    Del_deleted = (uint64_t*)malloc(Del_size_deleted * sizeof(uint64_t));
    Del_buffer = (uint32_t*)malloc((num_vertices + 1) * sizeof(uint64_t));
}

void TetMesh::releaseTmpStruct()
{
    free(Del_buffer);
    free(Del_deleted);
    free(Del_tmp);
}
