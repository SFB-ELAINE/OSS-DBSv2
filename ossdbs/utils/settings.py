from .materials import MATERIALS


class Settings:

    CUSTOM_SETTING = {'BrainRegion': {'Center': {'x[mm]': 0, 'y[mm]': 0, 'z[mm]': 0},
                                      'Dimension': {'x[mm]': 40, 'y[mm]': 40, 'z[mm]': 40},
                                      'Shape': 'Ellipsoid'
                                      }
                      }

    ELECTRODE_SETTING = {'Name': 'BostonScientificVerciseDirected',
                         'Rotation[Degrees]': 0.0,
                         'Direction': {'x[mm]': 0.0, 'y[mm]': 0.0, 'z[mm]': 1.0},
                         'TipPosition': {'x[mm]': 0.0, 'y[mm]': 0.0, 'z[mm]': 0.0},
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

    ENCAPSULATION_SETTING = {'Thickness[mm]': 0.0,
                             'Material': 'Gray matter',
                             'DielectricModel': 'ColeCole4',
                             'DielectricParameters': None,
                             'MaxMeshSize': 1e6
                             }

    SETTING = {'Electrodes': [],
               'Surfaces': [],
               'MaterialDistribution':
               {'MRIPath': '',
                'MRIMapping': MATERIALS,
                'DiffusionTensorActive': False,
                'DTIPath': ''
                },
               'DielectricModel': {'Type': 'ColeCole4',
                                   'CustomParameters': None
                                   },
               'Mesh': {'LoadMesh': False,
                        'LoadPath': '',
                        'MeshingHypothesis':
                        {'Type': 'Default',
                         'MaxMeshSize': 1e6
                         },
                        'MeshSize':
                        {'Edges': {},
                         'Faces': {},
                         'Volumes': {}
                         },
                        'SaveMesh': False,
                        'SavePath': 'mesh'
                        },
               'EQSMode': False,
               'FEMOrder': 2,
               'ComputeImpedance': False,
               'StimulationSignal':
               {'Type': 'Rectangle',
                'ListOfFrequencies': [130.0],
                'Frequency[Hz]': 130.0,
                'PulseWidth[us]': 60.0,
                'PulseTopWidth[us]': 0.0,
                'CounterPulseWidth[us]': 0.0,
                'InterPulseWidth[us]': 0.0,
                'SpectrumMode': 'FullSpectrum',
                'CutoffFrequency': 1e8,
                'CurrentControlled': False
                },
               'Solver': {'Type': 'CG',
                          'Preconditioner': 'bddc',
                          'PreconditionerKwargs': {},
                          'PrintRates': False,
                          'MaximumSteps': 10000,
                          'Precision': 1e-12
                          },
               'PointModel':
               {'Pathway': {'Active': False,
                            'FileName': ''
                            },
                'Lattice': {'Center': {'x[mm]': 0, 'y[mm]': 0, 'z[mm]': 0},
                            'Shape': {'x': 10, 'y': 10, 'z': 10},
                            'Direction': {'x[mm]': 0, 'y[mm]': 0, 'z[mm]': 1},
                            'PointDistance[mm]': 0.1
                            },
                'VoxelLattice': {
                    'Active': False,
                    'Shape': {'x': 10, 'y': 10, 'z': 10},
                    }
                },
               'OutputPath': 'Results',
               'SaveImpedance': False,
               'ExportVTK': False,
               'TemplateSpace': True,
               "ModelSide": 0,
               'CalcAxonActivation': False,
               'ActivationThresholdVTA': 0
               }

    def __init__(self, partial_settings: dict) -> None:
        self._partial_settings = partial_settings

    def complete_settings(self) -> dict:
        settings = self.CUSTOM_SETTING.copy()
        settings.update(self.SETTING.copy())
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
            for key, value in self.ENCAPSULATION_SETTING.items():
                if key not in settings['Electrodes'][index]['EncapsulationLayer']:
                    settings['Electrodes'][index]['EncapsulationLayer'][key] = value

    def _update_contacts(self, contacts: list) -> None:
        for index, contact in enumerate(contacts):
            for key, value in self.CONTACT_SETTING.items():
                if key not in contact:
                    contacts[index][key] = value
