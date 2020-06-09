#include <iostream>

#include "AMRIO.H"
#include "FArrayBox.H"
#include "GRChomboUserVariables.hpp"
#include "LevelData.H"
#include "MultigridUserVariables.hpp"
#include "NodeAMRIO.H"

#ifndef READ_HDF5_H
#define READ_HDF5_H

void Read_grid_from_HDF5(Vector<DisjointBoxLayout> &aa_grids, 
               //Vector<LevelData<FArrayBox> *> &a_multigrid_vars, 
               const std::string filename, 
               const IntVect i_ghosts)
{
    Vector<LevelData<FArrayBox>*> in_ld;
    Vector<string> in_names;
    Box in_domain;
    Real in_dx, in_dt, in_time;
    Vector <int> in_refRatio;
    int in_numLevels;
    
    // Now we read
    int in_status = ReadAMRHierarchyHDF5(filename, aa_grids, in_ld, in_names, in_domain, in_dx,
        in_dt, in_time, in_refRatio, in_numLevels);
}

void Read_vars_from_HDF5(Vector<DisjointBoxLayout> aa_grids, 
                         Vector<LevelData<FArrayBox> *> &a_multigrid_vars, 
                         const std::string filename, 
                         const IntVect i_ghosts
                        )
{
    Vector<DisjointBoxLayout> in_grids;
    Vector<LevelData<FArrayBox>*> in_ld;
    Vector<string> in_names;
    Box in_domain;
    Real in_dx, in_dt, in_time;
    Vector <int> in_refRatio;
    int in_numLevels;

    // Now we read
    int in_status = ReadAMRHierarchyHDF5(filename, in_grids, in_ld, in_names, 
        in_domain, in_dx, in_dt, in_time, in_refRatio, in_numLevels);

    for (int ilev = 0; ilev < 2; ilev++)
    {
        a_multigrid_vars[ilev] = new LevelData<FArrayBox>(aa_grids[ilev], NUM_MULTIGRID_VARS, i_ghosts);
    }

    // Now we copy the the parameters from the input file to multigrid_vars. This is a 
    // bit messy since they need to end up in the order given in MultigridUserVariables.hpp
    Interval temp_out;
    Interval temp_in;
    for (int level = 0; level < aa_grids.size(); level++)
    {
        // Multigrid Variables
        temp_out.define(c_chi, c_chi); // variables that are copied from file, in this line just chi
        temp_in.define(c_psi_reg, c_psi_reg); // variables that are coming in for the solver to use
        in_ld[level]->copyTo(temp_out, *a_multigrid_vars[level], temp_in);
        // temp_out.define(c_h11, c_h33); // the metric variables
        // temp_in.define(c_h11_0, c_h33_0);
        // in_ld[level]->copyTo(temp_out, *a_multigrid_vars[level], temp_in);
        temp_out.define(c_A11, c_A33); // traceless extrinsic curvature
        temp_in.define(c_A11_0, c_A33_0);
        in_ld[level]->copyTo(temp_out, *a_multigrid_vars[level], temp_in);
        temp_out.define(c_phi, c_phi); // matter field
        temp_in.define(c_phi_0, c_phi_0);
        in_ld[level]->copyTo(temp_out, *a_multigrid_vars[level], temp_in);
    }

}

#endif