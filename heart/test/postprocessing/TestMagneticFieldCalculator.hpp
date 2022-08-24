#ifndef _TESTMAGNETICFIELDCALCULATOR_HPP_
#define _TESTMAGNETICFIELDCALCULATOR_HPP_

#include "MagneticFieldCalculator.hpp"
#include "TrianglesMeshReader.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "Hdf5DataWriter.hpp"
#include "PetscVecTools.hpp"

class TestMagneticFieldCalculator : public CxxTest::TestSuite
{

public:


    /**
     * Create a linear gradient of electric potential along the z-direction
     * in a cube. Compute magnetic field outside the cube. 
     * 
     */
    void TestCubic()
    {
        // Read cubic mesh with 2mm sides
        TrianglesMeshReader<3,3> reader("mesh/test/data/cube_2mm_152_elements");
        DistributedTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);
        // Correct result is expressed in cT (10^-2 T), which is the unit of
        // magnetic field flux density in Chaste units

        const double correct_magnetic_field_x=2.47533e-8;  
        unsigned number_nodes = mesh.GetNumNodes();

        DistributedVectorFactory factory(number_nodes);
        Hdf5DataWriter writer(factory, "hdf5", "linear_V", false);
        writer.DefineFixedDimension(number_nodes);
        
        int node_id = writer.DefineVariable("Node", "dimensionless");
        int V_id = writer.DefineVariable("V", "mV");
        
        writer.DefineUnlimitedDimension("Time", "msec");
        writer.EndDefineMode();
        
        
        Vec petsc_data_1 = factory.CreateVec();
        DistributedVector distributed_vector_1 = factory.CreateDistributedVector(petsc_data_1);

        Vec petsc_data_2 = factory.CreateVec();
        DistributedVector distributed_vector_2 = factory.CreateDistributedVector(petsc_data_2);


        // Create a result with an electric potential that has a linear
        // gradient along the z-axis from V = 10 mV to V = -10 mV. Chaste
        // uses mV for voltage and cm for length, so 
        //
        // k = (-10mV - 10mV)/(0.2cm - 0cm) = -100mV/cm
        // V = k*z -10mV
        //
        const double k = -100;
        unsigned time_step = 0;

        for (DistributedTetrahedralMesh<3,3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
             node_iter != mesh.GetNodeIteratorEnd();
             ++node_iter)
        {
            distributed_vector_1[node_iter -> GetIndex()] = node_iter -> GetIndex(); 
            c_vector<double, 3> pos = node_iter -> rGetLocation();
            distributed_vector_2[node_iter -> GetIndex()] = pos[2] * k - 10;
        }
        distributed_vector_1.Restore();
        distributed_vector_2.Restore();
    
        writer.PutVector(node_id, petsc_data_1);
        writer.PutVector(V_id, petsc_data_2);
        writer.PutUnlimitedVariable(time_step);
        writer.AdvanceAlongUnlimitedDimension();
        writer.Close();

        PetscTools::Destroy(petsc_data_1);
        PetscTools::Destroy(petsc_data_2);

        // Compute the magnetic field at x = 1 mm, y = 0, z = -2 mm.
        ChastePoint<3> reference_point(0.1, 0, -0.2);
        MagneticFieldCalculator<3,3,1> calculator(mesh, reference_point,
                                       FileFinder("hdf5", RelativeTo::ChasteTestOutput),
                                       "linear_V", "V");
        double magnetic_field_x = calculator.ComputeMagneticFieldAtOneTimeStep(0);

        printf("%g", magnetic_field_x);
        TS_ASSERT_DELTA(magnetic_field_x, correct_magnetic_field_x, 1e-12);
    }

    void TestCrossProduct()
    {
        c_vector<double, 3> x_vec, y_vec, z_vec;
        x_vec[0] = 1.; x_vec[1] = 0.; x_vec[2] = 0.;
        y_vec[0] = 0.; y_vec[1] = 2.; y_vec[2] = 0.;
        z_vec[0] = 0.; z_vec[1] = 0.; z_vec[2] = 3.;
        c_vector<double, 3> result = CrossProduct(x_vec, y_vec);
        TS_ASSERT_DELTA(result[0], 0., 1e-6);
        TS_ASSERT_DELTA(result[1], 0., 1e-6);
        TS_ASSERT_DELTA(result[2], 2., 1e-6);
        result = CrossProduct(z_vec, y_vec);
        TS_ASSERT_DELTA(result[0], -6, 1e-6);
        TS_ASSERT_DELTA(result[1], 0, 1e-6);
        TS_ASSERT_DELTA(result[2], 0, 1e-6);
        result = CrossProduct(x_vec, z_vec);
        TS_ASSERT_DELTA(result[0], 0, 1e-6);
        TS_ASSERT_DELTA(result[1], -3, 1e-6);
        TS_ASSERT_DELTA(result[2], 0, 1e-6); 
    }
};


#endif /* _TESTMAGNETICFIELDCALCULATOR_HPP_ */