#ifndef _MAGNETICFIELDCALCULATOR_HPP_
#define _MAGNETICFIELDCALCULATOR_HPP_

#include "AbstractFunctionalCalculator.hpp"
#include "ChastePoint.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include "Hdf5DataReader.hpp"

/**
 * Calculates the magnetic field in cT (10^-2 T) from electric potential
 * quasistatically, using Ohm's law and Biot-Savart-law.
 * 
 * @tparam ELEMENT_DIM 
 * @tparam SPACE_DIM   Only 3D is currently supported.
 * @tparam PROBLEM_DIM Number of electrical potentials in the solution.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class MagneticFieldCalculator : public AbstractFunctionalCalculator<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>
{
private:

    friend class TestMagneticFieldCalculator;

    Hdf5DataReader* mpDataReader; /**< An HDF5 reader from which to get the solution*/
    unsigned mNumberOfNodes; /**< Number of nodes in the mesh (got from the data reader)*/
    unsigned mNumTimeSteps;/**< Number of time steps in the simulation (got from the data reader)*/
    AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& mrMesh;/**< A mesh used by the calculator*/
    ChastePoint<SPACE_DIM> mObservationPoint;
    std::string mVariableName;/**< the variable for which we want to calculate the pseudo ecg, defaults to "V"*/
    unsigned mTimestepStride; /**< The number of timesteps in a stride (so that we don't have to compute all the ECGs).  This defaults to 1.*/
    double mConductivity; /**< The conductivity of the domain (mS/cm)*/

    /**
     * @return the integrand
     * The magnetic field is computed using the Biot-Savart law:
     * 
     *  B(r) = \int \int \int_V grad(solution) cross (r-r') / |r-r'|^3
     * 
     * @param rX The point in space
     * @param rU The electric potential as a vector, u(i) = u_i
     * @param rGradU The gradient of the electric potential as a matrix, rGradU(i, j) = d(u_i)/d(X_j)
     */
    double GetIntegrand(ChastePoint<SPACE_DIM>& rX,
                        c_vector<double,PROBLEM_DIM>& rU,
                        c_matrix<double,PROBLEM_DIM,SPACE_DIM>& rGradU);

public:

    /**
     * Calculates the magnetic field and returns its value at the given time step.
     *
     * @param timeStep the time step where we want to calculate the magnetic field
     * @return the magnetic field at the given time step.
     *
     */
    double ComputeMagneticFieldAtOneTimeStep (unsigned timeStep);


    /**
     * The parameters just copied from PseudoEcgCalculator. 
     * 
     */
    MagneticFieldCalculator(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh,
                            const ChastePoint<SPACE_DIM>& rObservationPoint,
                            const FileFinder& rDirectory,
                            const std::string& rHdf5FileName,
                            const std::string& rVariableName = "V",
                            unsigned timestepStride = 1);

    /**
     * Destructor
     */
    ~MagneticFieldCalculator();
};


/**
 * Calculates the cross product (u x v) of two 3D-vectors u and v.
 * 
 * @param u the first vector
 * @param v the second vector
 * @return c_vector<double, 3> 
 */
c_vector<double, 3> CrossProduct(const c_vector<double, 3>& u, const c_vector<double, 3>& v);




#endif /* _MAGNETICFIELDCALCULATOR_ */