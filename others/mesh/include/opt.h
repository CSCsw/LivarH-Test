#ifndef FEASNEWT_H
#define FEASNEWT_H

#include "mesh.h"
extern "C" {
int optMesh(Mesh *m, double conv_tol, int precond);
}
#endif

