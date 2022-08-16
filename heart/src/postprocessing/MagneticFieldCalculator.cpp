#include "MagneticFieldCalculator.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
double MagneticFieldCalculator<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM> ::GetIntegrand(ChastePoint<SPACE_DIM>& rX,
                                    c_vector<double,PROBLEM_DIM>& rU,
                                    c_matrix<double,PROBLEM_DIM,SPACE_DIM>& rGradU)
{
    return 0.0;
}
