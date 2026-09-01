#pragma once

// Uncommenting the following macro definition makes the code use modified parts of hxt_SeqDel (Copyright (C) 2018 Célestin Marot).
// hxt_SeqDel is a sequential Delaunay triangulator hosted at https://git.immc.ucl.ac.be/hextreme/hxt_seqdel as of 2020.
// hxt_SeqDel is GPL licensed, meaning that if you uncomment the following line you accept the terms of the GPL license for
// the whole code which uses this library.
// If you need to use this code under the less restrictive LGPL license, please comment the following line.
// This will make the code slightly slower.
#define USE_MAROTS_METHOD

#include "delaunay.h"

typedef basicVec3d pointType;
typedef TetMesh_t<pointType> TetMesh;
