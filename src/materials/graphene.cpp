/*
 * File:   graphene.cpp
 * Author:  Paul Cazeaux
 *
 * Created on June 12, 2017, 9:00 AM
 */

#include "materials/graphene.h"

using Graphene::Atom;
using Graphene::Orbital;

double Coupling::Intralayer::graphene(
                        const Orbital orbit_row, const Orbital orbit_col,
                        const std::array<int, 2>& vector)
{

    /**
     * Compute hexagonal homogeneous coordinates from grid coordinates.
     * We use the following system:
     *
     *                A                   A                 |             (-1,2,-1)            (1,1,-2)
     *                             (1/2, sqrt(3)/2)         |
     *                          .                           |                        (0,1,-1)
     *                                                      |
     *                B                   B                 |              (-1,1,0)            (1,0,-1)
     *                                                      |
     *                                                      |
     *      A                   A                    A      |   (-2,1,1)             (0,0,0)              (2,-1,-1)
     *                        (0,0)                (1,0)    |
     *                                                      |
     *                .                   .                 |              (-1,0,1)            (1,-1,0)
     *                                                      |
     *                          B                           |                        (0,-1,1)
     *                                                      |
     *                A                   A                 |             (-1,-1,2)            (1,-2,1)
     *
     * In this system, the distance between a point and the center is still the euclidean distance, on all three coordinates.
     */
    std::array<int, 2> hom_vec {{   2 * vector[0] +     vector[1],
                                   -1 * vector[0] +     vector[1] }};
       // redundant 3rd coordinate: -1 * vector[0] - 2 * vector[1]

    /* Shift the arrow vector by the orbital coordinates */

    if (atom(orbit_col) == Atom::A && atom(orbit_row) == Atom::B)
    {
        hom_vec[0] -=  1;
        // hom_vec[2] -= -1;
    }
    else if (atom(orbit_col) == Atom::B && atom(orbit_row) == Atom::A)
    {
        hom_vec[0] +=  1;
        // hom_vec[2] += -1;
    }

    /**
     * Compute the distance to the origin:
     * we use the identity r = x^2 + y^2 + z^2 = x^2 + y^2 + (x+y)^2 = 2 * (x * (x+y) + y^2)
     */
    int r = hom_vec[0] * (hom_vec[0] + hom_vec[1]) + hom_vec[1]*hom_vec[1];

    const double onsite = 0.3504;
	  //const double t_arr[9] = {0, -2.8922, 0, 0, 0, 0, 0, 0, 0}; // only nearest neighbors (debug/testing)!
    const double t_arr[9] = {onsite, -2.8922, 0.2425, -0.2656, 0.0235, 0.0524, -0.0209, -0.0148, -0.0211};
    switch (r)
    {
        case 0:
            return t_arr[0];
        case 1:
            return t_arr[1];
        case 3:
            return t_arr[2];
        case 4:
            return t_arr[3];
        case 7:
            return t_arr[4];
        case 9:
            return t_arr[5];
        case 12:
            return t_arr[6];
        case 13:
            return t_arr[7];
        case 16:
            return t_arr[8];
        default:
            return 0.;
    }
}

double Coupling::Interlayer::C_to_C(const Orbital orbit_row, const Orbital orbit_col,
                    std::array<double, 3> vector,
                    const double theta_row, const double theta_col)
{
    /* Shift the arrow vector by the orbital coordinates */

        /**
         * !!!!!!!!    Counterclockwise rotation angles theta_row and theta_col !!!!!!!!
         *                      ToDo: discuss with Stephen to harmonize
         */

    /*
    if (atom(orbit_row) == Atom::B)
    {
        vector[0] -=  std::cos(theta_row) * Graphene::atom_pos.at (Atom::B)[0] + std::sin(theta_row) * Graphene::atom_pos.at (Atom::B)[1];
        vector[1] -= -std::sin(theta_row) * Graphene::atom_pos.at (Atom::B)[0] + std::cos(theta_row) * Graphene::atom_pos.at (Atom::B)[1];
    }
    if (atom(orbit_col) == Atom::B)
    {
        vector[0] +=  std::cos(theta_col) * Graphene::atom_pos.at (Atom::B)[0] + std::sin(theta_col) * Graphene::atom_pos.at (Atom::B)[1];
        vector[1] += -std::sin(theta_col) * Graphene::atom_pos.at (Atom::B)[0] + std::cos(theta_col) * Graphene::atom_pos.at (Atom::B)[1];
    }
    */

    // can uncomment to turn of interlayer coupling in graphene
    //return 0.0;
    if (abs(vector[2])>5){
        //printf("The layer is decoupled");
        return 0;
    }
    
    double r = std::sqrt(vector[0]*vector[0] + vector[1]*vector[1]);
    if (r < Graphene::inter_cutoff_radius)
    {
        // An angularly symmetric form, for testing purposes
        //return 0.310*std::exp(-1.75*(r/Graphene::a)*(r/Graphene::a));

        double ac = std::atan2(vector[1], vector[0]);

        // theta21 (angle to bond on sheet 1)

        double theta21 = ac - theta_row;
        if (orbit_row == Orbital::B_pz)
            theta21 += numbers::PI_6;
        else // (orbit_row == Orbital::A_pz)
            theta21 -= numbers::PI_6;

        // theta12 (angle to bond on sheet 2)

        double theta12 = ac - theta_col + numbers::PI;
        if (orbit_col == Orbital::B_pz)
            theta12 += numbers::PI_6;
        else // (orbit2_col == Orbital::A_pz)
            theta12 -= numbers::PI_6;

        double rs = r/Graphene::a;
        double z_eq = 3.35; // Equilbrium AB height from our DFT Calculations
    		double z_eps = (fabs(vector[2]) - z_eq)/z_eq;
        

        // New coefficients (with compression dependence fitting)
        // /*
    		double lambda_0 =  0.310 - 1.882*z_eps + 7.741*z_eps*z_eps;
    		double xi_0     =  1.750 - 1.618*z_eps + 1.848*z_eps*z_eps;
    		double kappa_0  =  1.990 + 1.007*z_eps + 2.427*z_eps*z_eps;

    		double lambda_3 = -0.068 + 0.399*z_eps - 1.739*z_eps*z_eps;
    		double xi_3     =  3.286 - 0.914*z_eps + 12.011*z_eps*z_eps;
    		double x_3      =  0.500 + 0.322*z_eps + 0.908*z_eps*z_eps;

    		double lambda_6 = -0.008 + 0.046*z_eps - 0.183*z_eps*z_eps;
    		double xi_6     =  2.272 - 0.721*z_eps - 4.414*z_eps*z_eps;
    		double x_6      =  1.217 + 0.027*z_eps - 0.658*z_eps*z_eps;
    		double kappa_6  =  1.562 - 0.371*z_eps - 0.134*z_eps*z_eps;
        // */
        // Debugging the new coefficients (with compression dependence fitting)
        /*
    		double lambda_0 =  0.310 + 1.882*z_eps + 7.741*z_eps*z_eps;
    		double xi_0     =  1.750 + 1.618*z_eps + 1.848*z_eps*z_eps;
    		double kappa_0  =  1.989 + 1.006*z_eps +  2.426*z_eps*z_eps;

    		double lambda_3 = -0.068 + 0.399*z_eps + 1.739*z_eps*z_eps;
    		double xi_3     =  3.286 - 0.916*z_eps + 11.996*z_eps*z_eps;
    		double x_3      =  0.500 + 0.322*z_eps + 0.906*z_eps*z_eps;

    		double lambda_6 = -0.008 + 0.046*z_eps - 0.183*z_eps*z_eps;
    		double xi_6     =  2.272 + 0.719*z_eps - 4.415*z_eps*z_eps;
    		double x_6      =  1.5206 + 0.027*z_eps - 0.658*z_eps*z_eps;
    		double kappa_6  =  1.561 + 0.371*z_eps - 0.134*z_eps*z_eps;
        */
        // Old coefficients (From: S. Fang, Phys. Rev. B 93, 235153 – Published 27 June 2016)
        /*
        double lambda_0 =  0.3155;
    		double xi_0     =  1.7543;
    		double kappa_0  =  2.0010;

    		double lambda_3 = -0.0688;
    		double xi_3     =  3.4692;
    		double x_3      =  0.5212;

    		double lambda_6 = -0.0083;
    		double xi_6     =  2.8764;
    		double x_6      =  1.5206;
    		double kappa_6  =  1.5731;
        */

        double V0 = lambda_0 *         std::exp(-xi_0*(rs    )*(rs    )) * std::cos(kappa_0*rs);
        double V3 = lambda_3 * rs*rs * std::exp(-xi_3*(rs-x_3)*(rs-x_3))                       ;
        double V6 = lambda_6 *         std::exp(-xi_6*(rs-x_6)*(rs-x_6)) * std::sin(kappa_6*rs);

    		// Values from model w/o compression fitting
    		/*
            double V0 = .3155 * std::exp(-1.7543*rs*rs) * std::cos(2.001*rs);
            double V3 = -.0688 * rs*rs * std::exp(-3.4692*(rs-.5212)*(rs-.5212));
            double V6 = -.0083 * std::exp(-2.8764*(rs-1.5206)*(rs-1.5206)) * std::sin(1.5731*rs);
    		*/

        double t = V0 + V3*(std::cos(3*theta12)+std::cos(3*theta21)) + V6*(std::cos(6*theta12)+std::cos(6*theta21));

        /**
         * Smooth cutoff
         */
        double d = Graphene::inter_cutoff_radius - r;
        if (d < 1. && d > 0.)
        {
            t *= std::exp(1. - 1./(d*d) );
        }
        return t;
    }
    else
        return 0.;
}

bool IsNonZero::Intralayer::graphene(
                        const Orbital orbit_row, const Orbital orbit_col,
                        const std::array<int, 2>& vector)
{
    /**
     * Compute hexagonal homogeneous coordinates from grid coordinates.
     * We use the following system:
     *
     *                A                   A                 |             (-1,2,-1)            (1,1,-2)
     *                             (1/2, sqrt(3)/2)         |
     *                          .                           |                        (0,1,-1)
     *                                                      |
     *                B                   B                 |              (-1,1,0)            (1,0,-1)
     *                                                      |
     *                                                      |
     *      A                   A                    A      |   (-2,1,1)             (0,0,0)              (2,-1,-1)
     *                        (0,0)                (1,0)    |
     *                                                      |
     *                .                   .                 |              (-1,0,1)            (1,-1,0)
     *                                                      |
     *                          B                           |                        (0,-1,1)
     *                                                      |
     *                A                   A                 |             (-1,-1,2)            (1,-2,1)
     *
     * In this system, the distance between a point and the center is still the euclidean distance, on all three coordinates.
     */
    std::array<int, 2> hom_vec  {{   2 * vector[0] +     vector[1],
                                    -1 * vector[0] +     vector[1] }};
       // redundant 3rd coordinate: -1 * vector[0] - 2 * vector[1]

    /* Shift the arrow vector by the orbital coordinates */
    if (atom(orbit_col) == Atom::A && atom(orbit_row) == Atom::B)
    {
        hom_vec[0] -=  1;
        // hom_vec[2] -= -1;
    }
    else if (atom(orbit_col) == Atom::B && atom(orbit_row) == Atom::A)
    {
        hom_vec[0] +=  1;
        // hom_vec[2] += -1;
    }
    /**
     * Compute the distance to the origin:
     * we use the identity r = x^2 + y^2 + z^2 = x^2 + y^2 + (x+y)^2 = 2 * (x * (x+y) + y^2)
     */
    int r = hom_vec[0] * (hom_vec[0] + hom_vec[1]) + hom_vec[1]*hom_vec[1];

    switch (r)
    {
        case 0:
        case 1:
        case 3:
        case 4:
        case 7:
        case 9:
        case 12:
        case 13:
        case 16:
            return true;
        default:
            return false;
    }
}

bool IsNonZero::Interlayer::C_to_C(const Orbital orbit_row, const Orbital orbit_col,
                    std::array<double, 3> vector,
                    const double theta_row, const double theta_col)
{

    /* Shift the arrow vector by the orbital coordinates */
    /*
    if (atom(orbit_row) == Atom::B)
    {
        vector[0] -=  std::cos(theta_row) * Graphene::atom_pos.at (Atom::B)[0] + std::sin(theta_row) * Graphene::atom_pos.at (Atom::B)[1];
        vector[1] -= -std::sin(theta_row) * Graphene::atom_pos.at (Atom::B)[0] + std::cos(theta_row) * Graphene::atom_pos.at (Atom::B)[1];
    }
    if (atom(orbit_col) == Atom::B)
    {
        vector[0] +=  std::cos(theta_col) * Graphene::atom_pos.at (Atom::B)[0] + std::sin(theta_col) * Graphene::atom_pos.at (Atom::B)[1];
        vector[1] += -std::sin(theta_col) * Graphene::atom_pos.at (Atom::B)[0] + std::cos(theta_col) * Graphene::atom_pos.at (Atom::B)[1];
    }
    */

    double r2 = vector[0]*vector[0] + vector[1]*vector[1];
    return (r2 < Graphene::inter_cutoff_radius * Graphene::inter_cutoff_radius);
}
