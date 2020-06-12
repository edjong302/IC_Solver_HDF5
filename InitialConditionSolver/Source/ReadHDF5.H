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
        // Set the HDF5 group and load level header
        char level_str[20];
        sprintf(level_str, "%d", level);
        const std::string label = std::string("level_") + level_str;
        handle.setGroup(label);
        level_header.readFromFile(handle);

        // Check periodicity of file
        bool periodic_bool[3];
        for (int dir = 0; dir < SpaceDim; ++dir)
        {
            char dir_str[20];
            sprintf(dir_str, "%d", dir);
            const std::string periodic_label =
                std::string("is_periodic_") + dir_str;
            if (level_header.m_int[periodic_label] == true)
            {
                periodic_bool[dir] = true;
                pout() << periodic_bool[dir] << endl;
            }
        }

        // Read boxes and assign processors to create DisjointBoxLayout
        // for this level
        ProblemDomain a_domain(level_header.m_box["prob_domain"], periodic_bool);
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