#ifndef _MAGNETICFIELDCALCULATOR_HPP_
#define _MAGNETICFIELDCALCULATOR_HPP_

#include "AbstractFunctionalCalculator.hpp"
#include "ChastePoint.hpp"
#include "AbstractTetrahedralMesh.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class MagneticFieldCalculator : public AbstractFunctionalCalculator<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>
{
private:

    friend class TestMagneticFieldCalculator;

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
     * The parameters just copied from PseudoEcgCalculator. 
     * 
     */
    MagneticFieldCalculator(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh,
                            const ChastePoint<SPACE_DIM>& rObservationPoint,
                            const FileFinder& rDirectory,
                            const std::string& rHdf5FileName,
                            const std::string& rVariableName = "B",
                            unsigned timestepStride = 1);

    /**
     * Destructor
     */
    ~MagneticFieldCalculator();
};

#endif /* _MAGNETICFIELDCALCULATOR_ */