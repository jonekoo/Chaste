#ifndef TESTSIMPLEVOLUMECALCULATOR_HPP_
#define TESTSIMPLEVOLUMECALCULATOR_HPP_

#include <cxxtest/TestSuite.h>

#include "AbstractFunctionalCalculator.hpp"
#include "TrianglesMeshReader.hpp"
#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"

template<unsigned DIM>
class SimpleVolumeCalculator : public AbstractFunctionalCalculator<DIM, DIM, 1>
{
    double GetIntegrand(ChastePoint<DIM>& rX,
                        c_vector<double,1>& rU,
                        c_matrix<double,1,DIM>& rGradU)
    {
        return 1.0;
    }
};


class TestSimpleVolumeCalculator : public CxxTest::TestSuite
{
public: 

    void TestWithVolumeCalculator()
    {
        TrianglesMeshReader<2,2> reader("mesh/test/data/square_128_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);
        DistributedTetrahedralMesh<2,2> distributed_mesh;
        distributed_mesh.ConstructFromMeshReader(reader);

        SimpleVolumeCalculator<2> volume_calculator;

        Vec vec = PetscTools::CreateAndSetVec(mesh.GetNumNodes(), 0.0);

        double result = volume_calculator.Calculate(mesh, vec);
        TS_ASSERT_DELTA(result, mesh.GetVolume(), 1e-12);
        double distributed_result = volume_calculator.Calculate(distributed_mesh, vec);
        TS_ASSERT_DELTA(result, distributed_result, 1e-12);

        PetscTools::Destroy(vec);
    }
};
#endif /*TESTSIMPLEVOLUMECALCULATOR_HPP_*/