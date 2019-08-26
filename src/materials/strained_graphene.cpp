
/*
 * File:   strained_graphene.h
 * Author:  Paul Cazeaux
 *
 * Created on June 12, 2017, 9:00 AM
 */

#include "materials/strained_graphene.h"
#include "tools/numbers.h"

using Graphene::Atom;
using Graphene::Orbital;
using namespace numbers;

double Coupling::Intralayer::strained_graphene(
        const Orbital orbit_row, const Orbital orbit_col,
        const std::array<int, 2>& vector)
{

    /* Compute hexagonal homogeneous coordinates from grid coordinates.
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
    /* Compute the distance to the origin:
     * we use the identity r = x^2 + y^2 + z^2 = x^2 + y^2 + (x+y)^2 = 2 * (x * (x+y) + y^2)
     */
    int r = hom_vec[0] * (hom_vec[0] + hom_vec[1]) + hom_vec[1]*hom_vec[1];


	// The t_0, or unperturbed, hopping parameters
    //const double t_arr[4] = {-3.6134, -2.8219,  0.2543, -0.1803}; // has vacuum-defined onsite energy
    const double t_arr[4] = {0.7833, -2.8219,  0.2543, -0.1803};
  	// the alpha, or scalar-like strain (u_xx + u_yy), scaling parameters
  	//const double a_arr[4] = {-4.8782,  4.0066, -0.4633,  0.6236};
  	// the beta, or vector-like strain ( <u_xx - u_yy, -2u_xy> ), scaling parameters
  	//const double b_arr[4] = { 0.0000, -3.0868,  0.8017,  0.4793};
    switch (r)
    {
        case 0: {
          return t_arr[0];
			}
        case 1: {
          return t_arr[1];

			}
        case 3: {
          return t_arr[2];

			}
        case 4: {
          return t_arr[3];

			}
        default:
            return 0.;
    }
}


double Coupling::Intralayer::strained_graphene(
        const Orbital orbit_row, const Orbital orbit_col,
        const std::array<int, 2>& vector, const std::vector< std::vector<double> >& strain)
{

    /* Compute hexagonal homogeneous coordinates from grid coordinates.
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
    /* Compute the distance to the origin:
     * we use the identity r = x^2 + y^2 + z^2 = x^2 + y^2 + (x+y)^2 = 2 * (x * (x+y) + y^2)
     */
    int r = hom_vec[0] * (hom_vec[0] + hom_vec[1]) + hom_vec[1]*hom_vec[1];


	// The t_0, or unperturbed, hopping parameters
    //const double t_arr[4] = {-3.6134, -2.8219,  0.2543, -0.1803}; // has vacuum-defined onsite energy
    const double t_arr[4] = {0.7833, -2.8219,  0.2543, -0.1803};
	// the alpha, or scalar-like strain (u_xx + u_yy), scaling parameters
	const double a_arr[4] = {-4.8782,  4.0066, -0.4633,  0.6236};
	// the beta, or vector-like strain ( <u_xx - u_yy, -2u_xy> ), scaling parameters
	const double b_arr[4] = { 0.0000, -3.0868,  0.8017,  0.4793};
    switch (r)
    {
        case 0: {
			// onsite, only need scalar strain
			double scalar_strain = strain[0][0] + strain[1][1];
            return (t_arr[0] + a_arr[0]*scalar_strain);
			}
        case 1: {
			// nearest neighbour, need to rotate by -pi/2

			std::vector< std::vector<double> > strain_rot;
			strain_rot.resize(2);
			strain_rot[0].resize(2);
			strain_rot[1].resize(2);

			double strain_theta = -PI_2;

			strain_rot[0][0] =  strain[0][0]*cos(strain_theta)*cos(strain_theta) +
								strain[1][1]*sin(strain_theta)*sin(strain_theta) +
								strain[0][1]*sin(strain_theta)*cos(strain_theta);
			strain_rot[1][1] =  strain[0][0]*sin(strain_theta)*sin(strain_theta) +
								strain[1][1]*cos(strain_theta)*cos(strain_theta) +
							 -2*strain[0][1]*sin(strain_theta)*cos(strain_theta);
			strain_rot[0][1] =  (strain[1][1] - strain[0][0])*sin(strain_theta)*cos(strain_theta) +
								 strain[0][1]*(cos(strain_theta)*cos(strain_theta) - sin(strain_theta)*sin(strain_theta));
			strain_rot[1][0] = strain_rot[0][1];

			double scalar_strain = strain_rot[0][0] + strain_rot[1][1];
			double vector_strain_y = strain_rot[0][0] - strain_rot[1][1];
            return (t_arr[1] + a_arr[1]*scalar_strain + b_arr[1]*vector_strain_y);

			}
        case 3: {
			// 2nd nearest neighbor, no rotation necessary!
			double scalar_strain = strain[0][0] + strain[1][1];
			double vector_strain_y = strain[0][0] - strain[1][1];
            return (t_arr[2] + a_arr[2]*scalar_strain + b_arr[2]*vector_strain_y);

			}
        case 4: {
			// 3rd nearest neighbor, need to rotate by +pi/2

			std::vector< std::vector<double> > strain_rot;
			strain_rot.resize(2);
			strain_rot[0].resize(2);
			strain_rot[1].resize(2);

			double strain_theta = PI_2;

			strain_rot[0][0] =  strain[0][0]*cos(strain_theta)*cos(strain_theta) +
								strain[1][1]*sin(strain_theta)*sin(strain_theta) +
								strain[0][1]*sin(strain_theta)*cos(strain_theta);
			strain_rot[1][1] =  strain[0][0]*sin(strain_theta)*sin(strain_theta) +
								strain[1][1]*cos(strain_theta)*cos(strain_theta) +
							 -2*strain[0][1]*sin(strain_theta)*cos(strain_theta);
			strain_rot[0][1] =  (strain[1][1] - strain[0][0])*sin(strain_theta)*cos(strain_theta) +
								 strain[0][1]*(cos(strain_theta)*cos(strain_theta) - sin(strain_theta)*sin(strain_theta));
			strain_rot[1][0] = strain_rot[0][1];

			double scalar_strain = strain_rot[0][0] + strain_rot[1][1];
			double vector_strain_y = strain_rot[0][0] - strain_rot[1][1];
            return (t_arr[3] + a_arr[3]*scalar_strain - b_arr[3]*vector_strain_y);

			}
        default:
            return 0.;
    }
}

double Coupling::Intralayer::strained_graphene(
        const Orbital orbit_row, const Orbital orbit_col,
        const std::array<int, 2>& vector, const double r_hop)
{

    /* Compute hexagonal homogeneous coordinates from grid coordinates.
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
    /* Compute the distance to the origin:
     * we use the identity r = x^2 + y^2 + z^2 = x^2 + y^2 + (x+y)^2 = 2 * (x * (x+y) + y^2)
     */
    int r = hom_vec[0] * (hom_vec[0] + hom_vec[1]) + hom_vec[1]*hom_vec[1];


	// The t_0, or unperturbed, hopping parameters
    //const double t_arr[4] = {-3.6134, -2.8219,  0.2543, -0.1803}; // has vacuum-defined onsite energy
    const double t_arr[4] = {0.7833, -2.8219,  0.2543, -0.1803};
  	// the alpha, or scalar-like strain (u_xx + u_yy), scaling parameters
  	const double a_arr[4] = {-4.8782,  4.0066, -0.4633,  0.6236};
  	// the beta, or vector-like strain ( <u_xx - u_yy, -2u_xy> ), scaling parameters
  	const double b_arr[4] = { 0.0000, -3.0868,  0.8017,  0.4793};
    double lattice_a = 2.4768;

    switch (r)
    {
        case 0: {
            return t_arr[0];
			}
        case 1: {
			// nearest neighbour
        double scalar_strain = (r_hop-lattice_a/sqrt(3.0))/(lattice_a/sqrt(3.0));
        return (t_arr[1] + a_arr[1]*scalar_strain);

			}
        case 3: {
			// 2nd nearest neighbor
			      double scalar_strain = (r_hop-lattice_a)/(lattice_a);
            return (t_arr[2] + a_arr[2]*scalar_strain);

			}
        case 4: {
			// 3rd nearest neighbor
            double scalar_strain = (r_hop-2.0*lattice_a/sqrt(3.0))/(2.0*lattice_a/sqrt(3.0));
            return (t_arr[3] + a_arr[3]*scalar_strain);

			}
        default:
            return 0.;
    }
}

double Coupling::Intralayer::strained_graphene(
        const Orbital orbit_row, const Orbital orbit_col,
        const std::array<int, 2>& vector, const std::vector<double> hop_dir,
        const std::vector< std::vector<double> >& strain, const double theta_here)
{

    /* Compute hexagonal homogeneous coordinates from grid coordinates.
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
    /* Compute the distance to the origin:
     * we use the identity r = x^2 + y^2 + z^2 = x^2 + y^2 + (x+y)^2 = 2 * (x * (x+y) + y^2)
     */
    int r = hom_vec[0] * (hom_vec[0] + hom_vec[1]) + hom_vec[1]*hom_vec[1];

    double onedim = strain[0][0] + strain[1][1]; // u_xx + u_yy
    std::vector<double> twodim;
    twodim.resize(2);
    twodim[0] = strain[0][0] - strain[1][1]; // u_xx - u_yy;
    twodim[0] = strain[0][1] + strain[1][0]; // u_xy + u_yx;

    // Rotate hop_dir clockwise by theta_here to get dir in local basis

    // TODO: Should compute hop_dir from the disp_vector, not from relaxed positions!
    // e.g. Need a formula for bonding_angle as a function of home_vec;
    std::vector<double> w;
    w.resize(2);
    w[0] =  cos(theta_here)*hop_dir[0] + sin(theta_here)*hop_dir[1];
    w[1] = -sin(theta_here)*hop_dir[0] + cos(theta_here)*hop_dir[1];


	// The t_0, or unperturbed, hopping parameters
    //const double t_arr[4] = {-3.6134, -2.8219,  0.2543, -0.1803}; // has vacuum-defined onsite energy
    const double t_arr[4] = {0.7833, -2.8219,  0.2543, -0.1803};
  	// the alpha, or scalar-like strain (u_xx + u_yy), scaling parameters
  	const double a_arr[4] = {-4.8782,  4.0066, -0.4633,  0.6236};
  	// the beta, or vector-like strain ( <u_xx - u_yy, -2u_xy> ), scaling parameters
  	const double b_arr[4] = { 0.0000, -3.0868,  0.8017,  0.4793};
    double lattice_a = 2.4768;

    switch (r)
    {
        case 0: {
            return t_arr[0] + a_arr[0]*onedim;
			}
        case 1: {
			// nearest neighbour
        return t_arr[1] + a_arr[1]*onedim + b_arr[1]*(w[0]*twodim[0] + w[1]*twodim[1]);

			}
        case 3: {
			// 2nd nearest neighbor
        // NOTE: w is rotated 90 degrees CCW relative to bonding dir for this term:
        return t_arr[2] + a_arr[2]*onedim + b_arr[2]*(-w[1]*twodim[0] + w[0]*twodim[1]);

			}
        case 4: {
			// 3rd nearest neighbor
        return t_arr[3] + a_arr[3]*onedim + b_arr[3]*(w[0]*twodim[0] + w[1]*twodim[1]);


			}
        default:
            return 0.;
    }
}

bool IsNonZero::Intralayer::strained_graphene(
        const Orbital orbit_row, const Orbital orbit_col,
        const std::array<int, 2>& vector)
{
    /* Compute hexagonal homogeneous coordinates from grid coordinates.
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
    /* Compute the distance to the origin:
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
