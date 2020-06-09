#include <iostream>

#include "AMRIO.H"
#include "FArrayBox.H"
#include "LevelData.H"

#ifndef READ_HDF5_H
#define READ_HDF5_H

void Read_HDF5(Vector<DisjointBoxLayout> &a_grids, const std::string filename)
{
    Vector<LevelData<FArrayBox>*> in_ld;
    Vector<string> in_names;
    Box in_domain;
    Real in_dx, in_dt, in_time;
    Vector <int> in_refRatio;
    int in_numLevels;
    

    // Now we read
    int in_status = ReadAMRHierarchyHDF5(filename, a_grids, in_ld, in_names, in_domain, in_dx,
        in_dt, in_time, in_refRatio, in_numLevels);
}

#endif