

class Settings:

    CUSTOM_SETTING = {'BrainRegion':
                      {'Center': {'x[mm]': 0, 'y[mm]': 0, 'z[mm]': 0},
                       'Dimension': {'x[mm]': 0, 'y[mm]': 0, 'z[mm]': 8},
                       'Shape': "Ellipsoid"
                       },
                      }

    ELECTRODE_SETTING = {'Name': 'BostonScientificVerciseDirected',
                         'Rotation[Degrees]': 0.0,
                         'Direction': {'x[mm]': 0.0,
                                       'y[mm]': 0.0,
                                       'z[mm]': 1.0},
                         'TipPosition': {'x[mm]': 0.0,
                                         'y[mm]': 0.0,
                                         'z[mm]': 0.0},
                         'Contacts': [],
                         'EncapsulationLayer': {}
                         }

    CONTACT_SETTING = {'Contact_ID': 0,
                       'Active': False,
                       'Current[A]': 0.0,
                       'Voltage[V]': 0.0,
                       'Floating': False,
                       'SurfaceImpedance[Ohmm]': {'real': 0.0, 'imag': 0.0},
                       'MaxMeshSize': 1e6,
                       'MaxMeshSizeEdge': 1e6
                       }

    ENCAPSULATION_SETTING = {"Thickness[mm]": 0.0,
                             "Material": "",
                             "DielectricModel": "",
                             "DielectricParameters": {},
                             "MaxMeshSize": 1e6
                             }

    SETTING = {'DielectricModel':
               {'Type': 'ColeCole4'
                },
               'Electrodes': [],
               'EQSMode': False,
               'FEMOrder': 2,
               'MaterialDistribution': {
                   'MRIPath': './input_files/segmask.nii',
                   'DiffusionTensorActive': False,
                   'DTIPath': '',
               },
               'Surfaces': [],
               'Mesh': {
                   'LoadMesh': False,
                   'LoadPath': '',
                   'MeshingHypothesis':
                   {'Type': 'Default',
                    },
                   'SaveMesh': False,
                   'SavePath': 'mesh'
               },
               'MeshSize': {},
               'OutputPath': 'result',
               "ComputeImpedance": False,
               "SaveImpedance": False,
               "ExportVTK": False,
               'Solver': {
                   'Type': 'CG',
                   'Preconditioner': 'bddc',
                   'PrintRates': False,
                   'MaximumSteps': 10000,
                   'Precision': 1e-12
               },
               'StimulationSignal':
               {'Type': 'Rectangle',
                'ListOfFrequencies': [130.],
                'Frequency[Hz]': 130.0,
                'PulseWidth[us]': 60.0,
                'PulseTopWidth[us]': 0.0,
                'CounterPulseWidth[us]': 0.0,
                'InterPulseWidth[us]': 0.0,
                'SpectrumMode': 'FullSpectrum',
                'CutoffFrequency': 1e8,
                'CurrentControlled': False,
                },
               'PointModel':
               {'Pathway':
                {'Active': False,
                 'FileName': ''
                 },
                'Lattice':
                {'Center': {'x[mm]': 0, 'y[mm]': 0, 'z[mm]': 0},
                 'Direction': {'x[mm]': 0, 'y[mm]': 0, 'z[mm]': 1},
                 'PointDistance[mm]': 0.1,
                 'Shape': {'x': 10, 'y': 10, 'z': 10}
                 }
                }
               }

    def __init__(self, partial_settings: dict) -> None:
        self._partial_settings = partial_settings

    def complete_settings(self) -> dict:
        # TODO needs to be reworked
        settings = self.SETTING.copy()
        self._update(settings, self._partial_settings)
        self._update_electrodes(settings)
        return settings

    def _update(self, target: dict, settings: dict) -> dict:
        for key in [key for key in target.keys() if key in settings.keys()]:
            if isinstance(target[key], dict):
                self._update(target[key], settings[key])
            else:
                target[key] = settings[key]

    def _update_electrodes(self, settings: dict) -> None:
        for index, electrode in enumerate(settings['Electrodes']):
            for key, value in self.ELECTRODE_SETTING.items():
                if key not in electrode:
                    settings['Electrodes'][index][key] = value
                if key == 'Contacts':
                    contacts = settings['Electrodes'][index][key]
                    self._update_contacts(contacts)
                    settings['Electrodes'][index][key] = contacts

    def _update_contacts(self, contacts: list) -> None:
        for index, contact in enumerate(contacts):
            for key, value in self.CONTACT_SETTING.items():
                if key not in contact:
                    contacts[index][key] = value
