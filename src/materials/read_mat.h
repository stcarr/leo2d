/*
 * File:   read_mat.h
 * Author:  Stephen Carr
 *
 * Created on May 21, 2018, 5:46 PM
 */


#ifndef read_mat_h
#define read_mat_h

#include <array>
#include <string>

#include "materials/loaded_mat.h"
#include "geom/sdata.h"


/**
 * This module is intended as an interface between our data library about 2D materials
 * and the computational part of this program.
 *
 * In particular, it provides a list of materials, and for each one its geometrical
 * characteristics (lattice basis, search cutoff radius for the hamiltonian elements)
 * as well as a tight-binding model (number of orbitals and hopping terms.)
 *
 * The intralayer term is computed as a function of orbital indices as well
 * as the vector between lattice sites in integer (lattice) coordinates.
 * The interlayer term is computed as a function of orbital indices as well as real-space
 * vector between site positions, as well as some angular information as necessary.
 */

namespace ReadMat {
   /*********************************************************************/
   /* Methods for returning each material's geometry and overall useful */
   /* spatial quantities related to tight-binding:                      */
   /*********************************************************************/
   /**
   * First, the lattice basis
   *
   * !!! CONVENTION !!!
   * The returned basis is given as an array of columns, so to construct the corresponding
   * matrix, the i-th row and j-th column entry is basis[j][i].
   */
    std::array<std::array<double, 2>, 2>  getLattice(LoadedMat& mat, Sdata& sheet_data);
    /**
     * The number of orbitals in the tight-binding model implemented in this library
     * for this material.
     */
    int n_orbitals(LoadedMat& mat,  Sdata& sheet_data);
    /**
     * The cut-off radius for nonzero terms in the intra-layer Hamiltonian.
     */
    double intra_search_radius(LoadedMat& mat);
    /**
     * The cut-off radius for nonzero terms in the inter-layer Hamiltonian.
     * This depends only on the broader category of materials (graphene, TMDCs...)
     * at this point.
     */
    double inter_search_radius(LoadedMat& mat);
    /**
     * The vertical position of a given orbital in a given layered material, provided the
     * unit cell is centered around the z=0 plane.
     */
    double orbital_pos(LoadedMat& mat, int idx, int dim, Sdata& sheet_data);

    /***********************************************************************************/
    /* Methods for returning tight-binding hopping terms and associated nonzero checks */
    /***********************************************************************************/
    /**
     * Returns an intra-layer tight-binding hopping term for a one dimensional material,
     * for a given pair of orbitals and a lattice vector in grid coordinates computed
     * for the material mat.
     *
     * !!! CONVENTION !!!
     * The lattice vector goes FROM the 'row' orbital site TO the 'column' orbital site.
     */
    double intralayer_term(int orbital_row, int orbital_col,
                    std::array<int, 2>& vector,
                    LoadedMat& mat, Sdata& sheet_data);
    /**
     * Returns whether a given intra-layer tight-binding hopping term for a two dimensional material
     * is nonzero, for a given pair of orbitals and a lattice hopping vector in grid coordinates computed
     * for the material mat.
     * Calling this method is cheaper than the corresponding calculation of the hopping term and
     * should be used when building the minimal static sparsity pattern.
     *
     * !!! CONVENTION !!!
     * The lattice vector goes FROM the 'row' orbital site TO the 'column' orbital site.
     */

    /**
     * Returns an inter-layer tight-binding hopping term between a pair of two dimensional
     * materials, for a given pair of orbitals and a real-space hopping vector
     * in cartesian coordinates.
     *
     * The two angle variables are the twist angles for the corresponding layer, with respect
     * to the un-rotated lattice as given in this library in the corresponding header file.
     *
     * !!! CONVENTION !!!
     * Rotation angles use a COUNTERCLOCKWISE rotation convention.
     *
     * !!! CONVENTION !!!
     * The real-space vector goes FROM the 'row' orbital site TO the 'column' orbital site.
     */

     // has xy symm for TMDC loaded mats baked in.
     double interlayer_term_xy_sym(int orbital_row, int orbital_col,
                     std::array<double, 3>& vector,
                     double angle_row, double angle_col,
                     LoadedMat& mat,  Sdata& sheet_data1, Sdata& sheet_data2);

     double interlayer_term(int orbital_row, int orbital_col,
                    std::array<double, 3>& vector,
                    double angle_row, double angle_col,
                    LoadedMat& mat, Sdata& sheet_data1, Sdata& sheet_data2);

   /**
    * Also a method to read all the info for a loaded material from file
    */

    LoadedMat loadMat(std::string filename, std::vector<Sdata> sdata);

}   /* End namespace ReadMat */
#endif
