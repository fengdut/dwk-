#ifndef READINPUT_H
#define READINPUT_H
#include"tokamak.h"
#include"mode.h"

int read_tokamak(char* filename,Tokamak *ptok,Grid *pgrid,Slowing *pslowing,Mode *mode);
#endif
