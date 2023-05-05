from ossdbs.type_check import TypeChecker


setting = {
            'CaseGrounding': {
                'Active': False,
                'Current[A]': 0,
                'Voltage[V]': 0.0
                },
            'CurrentControled': False,
            'DielectricModel': {
                'Type': 'ColeCole4',
                'PathToCustomParameters': ''
                },
            'Electrodes': [{
                            'Name': 'Electrode',
                            'PathToCustomParameters': '',
                            'Rotation[Degrees]': 120,
                            'Direction': {'x[mm]': 0,
                                          'y[mm]': 0,
                                          'z[mm]': 1,
                                          },
                            'TipPosition': {'x[mm]': 0,
                                            'y[mm]': 0,
                                            'z[mm]': 0,
                                            },
                            'Contacts': [{"Contact_ID": 1,
                                          "Active": True,
                                          "Current[A]": 0.0,
                                          "Voltage[V]": 1.0,
                                          "Floating": False,
                                          "SurfaceImpedance[Ωm]": {"real": 0.0,
                                                                   "imag": 0.0
                                                                   }
                                          },
                                         ],
                        }
                    ],
            'Contacts': {
                'MaxMeshSizeHeight': 0.1
                },
            'EncapsulatingLayer': {
                'Thickness[mm]': 0.0,
                'Material': 'Blood',
                'MaxMeshSizeHeight': 0.5
                },
            'EQSMode': False,
            'Floating': {
                'Active': False,
                'FloatingImpedance': False
            },
            'MaterialDistribution': {
                'MRIPath': './input_files/segmask.nii',
                'DiffusionTensorActive': False,
                'DTIPath': '',
            },
            'Mesh': {
                'LoadMesh': False,
                'LoadPath': '',
                'MeshElementOrder': 2,
                'MeshingHypothesis': {
                    'Type': 'Default',
                    'MaxMeshSizeHeight': 0.0
                    },
                'SaveMesh': False
            },
            'OutputPath': 'result',
            'BrainRegion': {
                'Center': {'x[mm]': 0, 'y[mm]': 0, 'z[mm]': 0},
                'Dimension': {'x[mm]': 0, 'y[mm]': 0, 'z[mm]': 8}
                },
            'RegionOfInterest': {
                'Center': {'x[mm]': 0, 'y[mm]': 0, 'z[mm]': 0},
                'Dimension': {'x[mm]': 0, 'y[mm]': 0, 'z[mm]': 8}
                },
            'Solver': {
                'Type': 'CG',
                'Preconditioner': 'bddc',
                'PrintRates': False,
                'MaximumSteps': 10000,
                'Precision': 1e-12
            },
            'SpectrumMode': 'FullSpectrum',
            'StimulationSignal': {
                'Type': 'Rectangle',
                'Frequency[Hz]': 130.0,
                'PulseWidth[µs]': 60.0,
                'PulseTopWidth[µs]': 0.0,
                'CounterPulseWidth[µs]': 0.0,
                'InterPulseWidth[µs]': 0.0
            },
            'PointModel': {
                'Pathway': {
                    'Active': False,
                    'FileName': ''
                    },
                'Lattice': {
                    'Center': {'x[mm]': 0, 'y[mm]': 0, 'z[mm]': 0},
                    'Direction': {'x[mm]': 0, 'y[mm]': 0, 'z[mm]': 1},
                    'PointDistance[mm]': 0.1,
                    'Shape': {'x': 10, 'y': 10, 'z': 10}
                }
                }
            }


def test_type_check():
    TypeChecker().check(setting)
