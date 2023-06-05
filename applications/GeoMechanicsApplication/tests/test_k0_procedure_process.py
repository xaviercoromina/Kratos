import os

import KratosMultiphysics                as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

class KratosGeoMechanicsK0ProcedureProcessTests(KratosUnittest.TestCase):
    """
    This class contains tests which check the result of a K0 Procedure Process
    """

    def test_k0_procedure_k0_nc(self):
        """
        Test to check if CAUCHY_STRESS_XX is correctly derived from CAUCHY_STRESS_YY using K0_NC
        """

        test_name = os.path.join("test_k0_procedure_process", "test_k0_procedure_k0_nc")
        file_path = test_helper.get_file_path(test_name)

        # run simulation
        simulation = test_helper.run_kratos(file_path)

        # retrieve Cauchy stress tensor
        cauchy_stresses = test_helper.get_on_integration_points(simulation,Kratos.CAUCHY_STRESS_TENSOR)

        # compare cauchy_stress_xx = k0_nc * cauchy_stress_yy, cauchy_stress_xy = 0.
        k0_nc = 0.6
        sig_integrationpoint1_element1 = cauchy_stresses[0][0]
        sig_yy = sig_integrationpoint1_element1[1,1]
        sig_xx = sig_integrationpoint1_element1[0,0]
        self.assertAlmostEqual( sig_xx, k0_nc*sig_yy )
        sig_xy = sig_integrationpoint1_element1[0,1]
        self.assertEqual( sig_xy, 0.0 )

    def test_k0_procedure_k0_nc_layers(self):
        """
        Test to check if CAUCHY_STRESS_XX is correctly derived from CAUCHY_STRESS_YY using K0_NC,
        even when materials are stacked in layers
        """

        test_name = os.path.join("test_k0_procedure_process", "test_k0_procedure_k0_nc_layers")
        file_path = test_helper.get_file_path(test_name)

        # run simulation
        simulation = test_helper.run_kratos(file_path)

        # retrieve Cauchy stress tensor
        cauchy_stresses = test_helper.get_on_integration_points(simulation,Kratos.CAUCHY_STRESS_TENSOR)

        # compare top layer cauchy_stress_xx = k0_nc * cauchy_stress_yy, cauchy_stress_xy = 0.
        k0_nc = 0.6
        sig_integrationpoint1_element6 = cauchy_stresses[5][0]
        sig_yy = sig_integrationpoint1_element6[1,1]
        sig_xx = sig_integrationpoint1_element6[0,0]
        self.assertAlmostEqual( sig_xx, k0_nc*sig_yy )
        sig_xy = sig_integrationpoint1_element6[0,1]
        self.assertEqual( sig_xy, 0.0 )

        # compare bottom layer cauchy_stress_xx = k0_nc * cauchy_stress_yy, cauchy_stress_xy = 0.
        k0_nc = 0.8
        sig_integrationpoint1_element1 = cauchy_stresses[0][0]
        sig_yy = sig_integrationpoint1_element1[1,1]
        sig_xx = sig_integrationpoint1_element1[0,0]
        self.assertAlmostEqual( sig_xx, k0_nc*sig_yy )
        sig_xy = sig_integrationpoint1_element1[0,1]
        self.assertEqual( sig_xy, 0.0 )

    def test_k0_procedure_k0_nc_skew_layers(self):
        """
        Test to check if CAUCHY_STRESS_XX is correctly derived from CAUCHY_STRESS_YY using K0_NC,
        even when materials are stacked in skewed layers
        """

        test_name = os.path.join("test_k0_procedure_process", "test_k0_procedure_k0_nc_skew_layers")
        file_path = test_helper.get_file_path(test_name)

        # run simulation
        simulation = test_helper.run_kratos(file_path)

        # retrieve Cauchy stress tensor
        cauchy_stresses = test_helper.get_on_integration_points(simulation,Kratos.CAUCHY_STRESS_TENSOR)

        # compare top layer left cauchy_stress_xx = k0_nc * cauchy_stress_yy, cauchy_stress_xy = 0.
        k0_nc = 0.6
        sig_integrationpoint1_element3 = cauchy_stresses[2][0]
        sig_yy = sig_integrationpoint1_element3[1,1]
        sig_xx = sig_integrationpoint1_element3[0,0]
        self.assertAlmostEqual( sig_xx, k0_nc*sig_yy )
        sig_xy = sig_integrationpoint1_element3[0,1]
        self.assertEqual( sig_xy, 0.0 )

        # compare middle layer right cauchy_stress_xx = k0_nc * cauchy_stress_yy, cauchy_stress_xy = 0.
        k0_nc = 0.7
        sig_integrationpoint1_element80 = cauchy_stresses[79][0]
        sig_yy = sig_integrationpoint1_element80[1,1]
        sig_xx = sig_integrationpoint1_element80[0,0]
        self.assertAlmostEqual( sig_xx, k0_nc*sig_yy )
        sig_xy = sig_integrationpoint1_element80[0,1]
        self.assertEqual( sig_xy, 0.0 )

        # compare bottom layer left cauchy_stress_xx = k0_nc * cauchy_stress_yy, cauchy_stress_xy = 0.
        k0_nc = 0.8
        sig_integrationpoint1_element117 = cauchy_stresses[116][0]
        sig_yy = sig_integrationpoint1_element117[1,1]
        sig_xx = sig_integrationpoint1_element117[0,0]
        self.assertAlmostEqual( sig_xx, k0_nc*sig_yy )
        sig_xy = sig_integrationpoint1_element117[0,1]
        self.assertEqual( sig_xy, 0.0 )

    def test_k0_procedure_k0_nc_ocr(self):
        """
        Test to check if CAUCHY_STRESS_XX is correctly derived from CAUCHY_STRESS_YY using K0_NC and OCR
        """

        test_name = os.path.join("test_k0_procedure_process", "test_k0_procedure_k0_nc_ocr")
        file_path = test_helper.get_file_path(test_name)

        # run simulation
        simulation = test_helper.run_kratos(file_path)

        # retrieve Cauchy stress tensor
        cauchy_stresses = test_helper.get_on_integration_points(simulation,Kratos.CAUCHY_STRESS_TENSOR)

        # compare cauchy_stress_xx = k0 * cauchy_stress_yy, cauchy_stress_xy = 0.
        k0_nc      = 0.6
        poisson_ur = 0.
        ocr        = 1.5
        k0 = k0_nc * ocr + ( poisson_ur / ( 1.0 - poisson_ur ) ) * ( ocr - 1.0 )
        sig_integrationpoint1_element1 = cauchy_stresses[0][0]
        sig_yy = sig_integrationpoint1_element1[1,1]
        sig_xx = sig_integrationpoint1_element1[0,0]
        self.assertAlmostEqual( sig_xx, k0*sig_yy )
        sig_xy = sig_integrationpoint1_element1[0,1]
        self.assertEqual( sig_xy, 0.0 )

    def test_k0_procedure_k0_umat(self):
        """
        Test to check if CAUCHY_STRESS_XX is correctly derived from CAUCHY_STRESS_YY using K0_NC = 1 - sin( PHI ),
        with PHI from UMAT material parameters
        """

        test_name = os.path.join("test_k0_procedure_process", "test_k0_procedure_k0_umat")
        file_path = test_helper.get_file_path(test_name)

        # run simulation
        simulation = test_helper.run_kratos(file_path)

        # retrieve Cauchy stress tensor
        cauchy_stresses = test_helper.get_on_integration_points(simulation,Kratos.CAUCHY_STRESS_TENSOR)

        # compare cauchy_stress_xx = k0_nc * cauchy_stress_yy, cauchy_stress_xy = 0. k0_nc = 1 - sin( 30 degrees )
        k0_nc = 0.5
        sig_integrationpoint1_element1 = cauchy_stresses[0][0]
        sig_yy = sig_integrationpoint1_element1[1,1]
        sig_xx = sig_integrationpoint1_element1[0,0]
        self.assertAlmostEqual( sig_xx, k0_nc*sig_yy )
        sig_xy = sig_integrationpoint1_element1[0,1]
        self.assertEqual( sig_xy, 0.0 )

if __name__ == '__main__':
    KratosUnittest.main()
