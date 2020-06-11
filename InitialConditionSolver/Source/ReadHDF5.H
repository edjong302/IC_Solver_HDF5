#include <iostream>

#include "AMRIO.H"
#include "CH_HDF5.H"
#include "DisjointBoxLayout.H"
#include "FArrayBox.H"
#include "GRChomboUserVariables.hpp"
#include "LevelData.H"
#include "LoadBalance.H"
#include "MultigridUserVariables.hpp"
#include "NodeAMRIO.H"

#ifndef READ_HDF5_H
#define READ_HDF5_H

void Read_HDF5(Vector<DisjointBoxLayout> &aa_grids, std::string filename)
{
    aa_grids.resize(2);
    HDF5Handle handle(filename, HDF5Handle::OPEN_RDONLY);

    for (int level = 0; level < aa_grids.size(); level++)
    {
        HDF5HeaderData level_header;
        level_header.readFromFile(handle);
        char level_str[20];
        sprintf(level_str, "%d", level);
        const std::string label = std::string("level_") + level_str;
        handle.setGroup(label);

        // Read boxes and processor assignments to create DisjointBoxLayout
        // for this level
        ProblemDomain a_domain = level_header.m_box["prob_domain"];
        Vector<Box> boxes;
        int error = read(handle, boxes);
        Vector<int> procIDs;
        LoadBalance(procIDs, boxes);
        DisjointBoxLayout layout(boxes, procIDs, a_domain);
        aa_grids[level] = layout;
    }
    

    handle.close();
}

#endif