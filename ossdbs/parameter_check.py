
class ParamterExminer:

    TYPES = {'DiffusionTensorImage': {'Path': str},
             'Electrodes': [{"Name": str,
                             "Rotation[Degrees]": float,
                             "Direction": {"x[mm]": float,
                                           "y[mm]": float,
                                           "z[mm]": float},
                             "TipPosition": {"x[mm]": float,
                                             "y[mm]": float,
                                             "z[mm]": float},
                            "Contact_1": {
                                "Active": bool,
                                "Current[A]": float,
                                "Voltage[V]": float,
                                "Floating": bool,
                                "SurfaceImpedance[Ωm]": {"real": float,
                                                         "imag": float},
                                }
                             }],
             'EQSMode': bool,
             'Floating': {'Active': bool, 'FloatingImpedance': bool},
             'MaterialDistribution': {
                    'MaterialCoding': dict,
                    'MRIPath': str},
             'Mesh': {
                'LoadMesh': bool,
                'LoadPath': str,
                'MeshElementOrder': int,
                'MeshRefinement': str,
                'SavePath': str,
                },
                'OutputPath': str,
                'RegionOfInterest': {
                    'Center': {'x[mm]': float, 'y[mm]': float, 'z[mm]': float},
                    'Shape': {'x[mm]': float, 'y[mm]': float, 'z[mm]': float},
                },
                'Solver': {'Type': str,
                            'Preconditioner': str,
                            'PrintRates': bool,
                            'MaximumSteps': int,
                            'Precision': float,
                            },
                'SpectrumMode': str,
                'StimulationSignal': {
                            'Type': str,
                            'Frequency[Hz]': float,
                            'PulseWidth[µs]': float,
                            'PulseTopWidth[µs]': float,
                            'CounterPulseWidth[µs]': float,
                            'InterPulseWidth[µs]': float,
                        },
                'Points': str}

    def check_input(self, input: dict) -> None:
        pass
