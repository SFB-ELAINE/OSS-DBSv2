

class Settings:

    SETTING = {
                'CaseGrounding': {
                    'Active': False,
                    'Current[A]': 0.0,
                    'Voltage[V]': 0.0
                    },
                'CurrentControled': False,
                'DielectricModel': {
                    'Type': 'ColeCole4',
                    'PathToCustomParameters': ''
                    },
                'Electrodes': [],
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

    ELECTRODE_SETTING = {'Name': 'BostonScientificVerciseDirected',
                         'PathToCustomParameters': '',
                         'Rotation[Degrees]': 0.0,
                         'Direction': {'x[mm]': 0.0,
                                       'y[mm]': 0.0,
                                       'z[mm]': 1.0},
                         'TipPosition': {'x[mm]': 0.0,
                                         'y[mm]': 0.0,
                                         'z[mm]': 0.0},
                         'Contacts': [],
                         }

    CONTACT_SETTING = {'Contact_ID': 0,
                       'Active': False,
                       'Current[A]': 0.0,
                       'Voltage[V]': 0.0,
                       'Floating': False,
                       'SurfaceImpedance[Ωm]': {'real': 0.0, 'imag': 0.0}
                       }

    def __init__(self, partial_settings: dict) -> None:
        self.__partial_settings = partial_settings

    def complete_settings(self) -> dict:
        settings = self.SETTING.copy()
        self.__update(settings, self.__partial_settings)
        self.__update_electrodes(settings)
        return settings

    def __update(self, target: dict, settings: dict) -> dict:
        for key in [key for key in target.keys() if key in settings.keys()]:
            if isinstance(target[key], dict):
                self.__update(target[key], settings[key])
            else:
                target[key] = settings[key]

    def __update_electrodes(self, settings: dict) -> None:
        for index, electrode in enumerate(settings['Electrodes']):
            for key, value in self.ELECTRODE_SETTING.items():
                if key not in electrode:
                    settings['Electrodes'][index][key] = value
                if key == 'Contacts':
                    contacts = settings['Electrodes'][index][key]
                    self.__update_contacts(contacts)
                    settings['Electrodes'][index][key] = contacts

    def __update_contacts(self, contacts: list) -> None:
        for index, contact in enumerate(contacts):
            for key, value in self.CONTACT_SETTING.items():
                if key not in contact:
                    contacts[index][key] = value
