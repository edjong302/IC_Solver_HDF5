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

void Read_grid_from_HDF5(Vector<DisjointBoxLayout> &aa_grids,
               const PoissonParameters &a_params
               )
{
    HDF5Handle handle(a_params.read_from_file, HDF5Handle::OPEN_RDONLY);
    HDF5HeaderData header;
    header.readFromFile(handle);
    int num_levels = a_params.numLevels;
    pout() << "numlevels is " << num_levels << endl;
    aa_grids.resize(num_levels);

    for (int level = 0; level < num_levels; ++level)
    {
        HDF5HeaderData level_header;
        // Set the HDF5 group and load level header
        char level_str[20];
        sprintf(level_str, "%d", level);
        const std::string label = std::string("level_") + level_str;
        handle.setGroup(label);
        level_header.readFromFile(handle);

        // Check and write periodicity of file
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
            }
        }
        // For now we need the loop below because the final output bdry
        // conditions are always periodic.
        if (a_params.periodic[0] == 0)
        {
            for (int dir = 0; dir < SpaceDim; ++dir){periodic_bool[dir] = false;}
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
        pout() << "Closing handle\n";
    }
    handle.close();
    pout() << "Closed handle\n";
}

void transfer_data(LevelData<FArrayBox> &a_grchombo_vars,
                   LevelData<FArrayBox> &a_multigrid_vars)
{
    DataIterator dit = a_grchombo_vars.dataIterator();
    // Here comes loop over boxes on a given processor
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &grchombo_vars_box = a_grchombo_vars[dit()];
        FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];

        Box this_box = grchombo_vars_box.box();
        BoxIterator bit(this_box);
        // Here comes loop over all grid point in this box
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            // I want to copy the conformal factor and the 
            // extrinsic curvature info
            multigrid_vars_box(iv, c_A11_0) = 
                grchombo_vars_box(iv, c_A11);
            multigrid_vars_box(iv, c_A12_0) = 
                grchombo_vars_box(iv, c_A12);
            multigrid_vars_box(iv, c_A13_0) = 
                grchombo_vars_box(iv, c_A13);
            multigrid_vars_box(iv, c_A22_0) = 
                grchombo_vars_box(iv, c_A22);
            multigrid_vars_box(iv, c_A23_0) = 
                grchombo_vars_box(iv, c_A23);
            multigrid_vars_box(iv, c_A33_0) = 
                grchombo_vars_box(iv, c_A33);
            multigrid_vars_box(iv, c_psi_reg) = 
                pow(grchombo_vars_box(iv, c_chi), -.25);
            multigrid_vars_box(iv, c_phi_0) = 
                grchombo_vars_box(iv, c_phi);
        }
    }

}

void Read_vars_from_HDF5(LevelData<FArrayBox> &a_multigrid_vars,
                         const Vector<DisjointBoxLayout> aa_grids,
                         const PoissonParameters a_params, 
                         const IntVect a_ghosts, 
                         int level
                        )
{
    LevelData<FArrayBox> *gr_chombo_vars;
    gr_chombo_vars = new LevelData<FArrayBox>(aa_grids[level], NUM_GRCHOMBO_VARS, a_ghosts);

    HDF5Handle handle(a_params.read_from_file, HDF5Handle::OPEN_RDONLY);
    HDF5HeaderData level_header;
    // Set the HDF5 group and load level header
    char level_str[20];
    sprintf(level_str, "%d", level);
    const std::string label = std::string("level_") + level_str;
    handle.setGroup(label);
    level_header.readFromFile(handle);

    read(handle, *gr_chombo_vars, "data", aa_grids[level], Interval(), false);

    

    transfer_data(*gr_chombo_vars, a_multigrid_vars);

    // DataIterator dit = a_multigrid_vars.dataIterator();
    // for (dit.begin(); dit.ok(); ++dit)
    // {
    //     FArrayBox &fabox = gr_chombo_vars_variable[dit()];
    //     Box this_box = fabox.box();
    //     BoxIterator bit(this_box);
    //     for (bit.begin(); bit.ok(); ++bit)
    //     {
    //         IntVect iv = bit();
    //         if (iv[0] == 12 && iv[1] == 12)
    //         {
    //             pout() << fabox(iv, c_phi_0) <<
    //             " at " << iv[0] << " " << iv[1] << " " << iv[2] << endl;
    //         }
    //         if ((iv[0] < 0 || iv[1] < 0) || iv[2] < 0)
    //         {
    //             pout() << 
    //                 " at " << 
    //                 iv[0] << " " << iv[1] << " " << iv[2] << endl;
    //         } 
    //     }
    // }

    handle.close();

    // Clean up storage
    delete gr_chombo_vars;
    gr_chombo_vars = NULL;

}

#endif