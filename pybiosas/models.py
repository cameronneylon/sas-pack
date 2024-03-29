# Central repository for model information. Could also contain
# additional models built up from components

models = {
           'cylinder' : {
                         'library_name':'CylinderModel',
                         'model_name'  :'CylinderModel',
                         'test_data'   :'test_data_cylinder.xml',
                         'test_params' :[{"value" : 20,
                                          "paramname" : "radius"},
                                         {"value" : 50,
                                          "paramname" : "length"},
                                         {'paramname' : 'sldSolv',
                                          'value'     : 1e-5,
                                          'fixed'     : True},
                                         {'paramname' : 'scale',
                                          'value'     : 0.01,
                                          'fixed'     : True}],
                         'exp_vals'    :[{'paramname' : 'scale',
                                          'value'     : 0.01},
                                         {'paramname' : 'radius',
                                          'value'     : 30.0},
                                         {'paramname' : 'length',
                                          'value'     : 100.0},
                                         {'paramname' : 'sldCyl',
                                          'value'     : 2e-6},
                                         {'paramname' : 'sldSolv',
                                          'value'     : 1e-5},
                                         {'paramname' : 'background',
                                          'value'     : 0.2}]

                         },
                         
           'sphere'   : {
                         'library_name':'SphereModel',
                         'model_name'  :'SphereModel',
                         'test_data'   :'test_data_sphere.xml',
                         'test_params' : [{"value": 60.0,
                                           "paramname": "radius"},
                                          {'paramname' : 'sldSolv',
                                           'value'     : 1e-5,
                                           'fixed'     : True},
                                          {'paramname' : 'scale',
                                           'value'     : 0.01,
                                           'fixed'     : True}],
                         'exp_vals'    :[{'paramname' : 'scale',
                                          'value'     : 0.01},
                                         {'paramname' : 'radius',
                                          'value'     : 40.0},
                                         {'paramname' : 'sldSph',
                                          'value'     : 2e-6},
                                         {'paramname' : 'sldSolv',
                                          'value'     : 1e-5},
                                         {'paramname' : 'background',
                                          'value'     : 0.1}]
                         },
                         
           'ellipse'  : {
                         'library_name':'EllipsoidModel',
                         'model_name'  :'EllipsoidModel',
                         'test_data'   :'test_data_ellipse.xml',
                         'test_params' : [{"value" : 40,
                                           "paramname" : "radius_a"},
                                          {"value" : 100,
                                           "paramname" : "radius_b"},
                                          {'paramname' : 'sldSolv',
                                           'value'     : 1e-5,
                                           'fixed'     : True},
                                          {'paramname' : 'scale',
                                           'value'     : 0.01,
                                           'fixed'     : True}],
                         'exp_vals'    :[{'paramname' : 'scale',
                                          'value'     : 0.01},
                                         {'paramname' : 'radius_a',
                                          'value'     : 30.0},
                                         {'paramname' : 'radius_b',
                                          'value'     : 200.0},
                                         {'paramname' : 'sldEll',
                                          'value'     : 2e-6},
                                         {'paramname' : 'sldSolv',
                                          'value'     : 1e-5},
                                         {'paramname' : 'background',
                                          'value'     : 0.001}]
                         },
           'coreShellBicelle' : {
                         'library_name':'CoreShellBicelleModel',
                         'model_name'  :'CoreShellBicelleModel',
                         'test_data'   :'test_data_coreshellbicelle.xml',
                         'test_params' : [{'paramname' : 'solvent_sld',
                                           'value'     : 1e-5,
                                           'fixed'     : True},
                                          {'paramname' : 'scale',
                                           'value'     : 0.1,
                                           'fixed'     : True},
                                          {'paramname' : 'core_sld',
                                           'value'     : 2e-6,
                                           'fixed'     : True},
                                          {'paramname' : 'radius',
                                           'value'     : 100},
                                          {'paramname' : 'length',
                                           'value'     : 20},
                                          {'paramname' : 'rim_sld',
                                           'value'     : 4e-6}],
                         'exp_vals'    :[{'paramname' : 'scale',
                                          'value'     : 0.1},
                                         {'paramname' : 'radius',
                                          'value'     : 100},
                                         {'paramname' : 'rim_thick',
                                          'value'     : 12},
                                         {'paramname' : 'face_thick',
                                          'value'     : 12},
                                         {'paramname' : 'length',
                                          'value'     : 20},
                                         {'paramname' : 'core_sld',
                                          'value'     : 2e-6},
                                         {'paramname' : 'face_sld',
                                          'value'     : 1e-6},
                                         {'paramname' : 'rim_sld',
                                          'value'     : 4e-6},
                                         {'paramname' : 'solvent_sld',
                                          'value'     : 1e-5},
                                         {'paramname' : 'background',
                                          'value'     : 0.04}]
                                  },
            'ellipticalCylinder' : {
                         'library_name':'EllipticalCylinderModel',
                         'model_name'  :'EllipticalCylinderModel',
                         'test_data'   : '',
                         'test_params' : [{'paramname' : 'r_minor',
                                           'fixed'     : False,
                                           'value'     : 20.0},
                                          {'paramname' : 'scale',
                                           'fixed'     : True,
                                           'value'     : 1.0},
                                          {'paramname' : 'r_ratio',
                                           'fixed'     : False,
                                           'value'     : 1.5},
                                          {'paramname' : 'length',
                                           'fixed'     : False,
                                           'value'     : 400.0},
                                          {'paramname' : 'sldCyl',
                                           'fixed'     : True,
                                           'value'     : 4e-6},
                                          {'paramname' : 'sldSolv',
                                           'fixed'     : True,
                                           'value'     : 1e-6},
                                          {'paramname' : 'background',
                                           'fixed'     : False,
                                           'value'     : 0.0}],
                         'exp_vals'    : [{'paramname' : 'r_minor',
                                           'value'     : 20.0},
                                          {'paramname' : 'scale',
                                           'value'     : 1.0},
                                          {'paramname' : 'r_ratio',
                                           'value'     : 1.5},
                                          {'paramname' : 'length',
                                           'value'     : 400.0},
                                          {'paramname' : 'sldCyl',
                                           'value'     : 4e-6},
                                          {'paramname' : 'sldSolv',
                                           'value'     : 1e-6},
                                          {'paramname' : 'background',
                                           'value'     : 0.0}]
                                   }
                                          

        }
