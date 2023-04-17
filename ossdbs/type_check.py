

class Settings:

    TYPES = {
             'CaseGrounding': {
                 'Active': bool,
                 'Current[A]': float,
                 'Voltage[V]': float
                 },
             'CurrentControled': bool,
             'DielectricModel': {
                 'Type': str,
                 'PathToCustomParameters': str
                 },
             'Electrodes': list,
             'Contacts': {
                 'MaxMeshSizeHeight': float
                 },
             'EncapsulatingLayer': {
                 'Thickness[mm]': float,
                 'Material': str,
                 'MaxMeshSizeHeight': float
                 },
             'EQSMode': bool,
             'Floating': {
                 'Active': bool,
                 'FloatingImpedance': bool
             },
             'MaterialDistribution': {
                 'MRIPath': str,
                 'DiffusionTensorActive': bool,
                 'DTIPath': str,
             },
             'Mesh': {
                 'LoadMesh': bool,
                 'LoadPath': str,
                 'MeshElementOrder': int,
                 'MeshingHypothesis': {
                     'Type': str,
                     'MaxMeshSizeHeight': float
                     },
                 'SaveMesh': bool
             },
             'OutputPath': str,
             'RegionOfInterest': {
                 'Center': {'x[mm]': float,
                            'y[mm]': float,
                            'z[mm]': float},
                 'Dimension': {'x[mm]': float,
                               'y[mm]': float,
                               'z[mm]': float}
                 },
             'Solver': {
                 'Type': str,
                 'Preconditioner': str,
                 'PrintRates': bool,
                 'MaximumSteps': int,
                 'Precision': float
             },
             'SpectrumMode': str,
             'StimulationSignal': {
                 'Type': str,
                 'Frequency[Hz]': float,
                 'PulseWidth[µs]': float,
                 'PulseTopWidth[µs]': float,
                 'CounterPulseWidth[µs]': float,
                 'InterPulseWidth[µs]': float
             },
             'PointModel': {
                 'Pathway': {
                     'Active': bool,
                     'FileName': str
                     },
                 'Lattice': {
                     'Center': {'x[mm]': float,
                                'y[mm]': float,
                                'z[mm]': float},
                     'Direction': {'x[mm]': float,
                                   'y[mm]': float,
                                   'z[mm]': float},
                     'PointDistance[mm]': float,
                     'Shape': {'x': int,
                               'y': int,
                               'z': int}
                 }
             }
            }

    ELECTRODE_SETTING = {'Name': str,
                         'PathToCustomParameters': str,
                         'Rotation[Degrees]': float,
                         'Direction': {'x[mm]': float,
                                       'y[mm]': float,
                                       'z[mm]': float},
                         'TipPosition': {'x[mm]': float,
                                         'y[mm]': float,
                                         'z[mm]': float},
                         'Contacts': list,
                         }

    CONTACT_SETTING = {'Contact_ID': int,
                       'Active': bool,
                       'Current[A]': float,
                       'Voltage[V]': float,
                       'Floating': bool,
                       'SurfaceImpedance[Ωm]': {'real': float,
                                                'imag': float}
                       }

    def __init__(self, partial_settings: dict) -> None:
        self.__partial_settings = partial_settings

    def complete_settings(self) -> None:
        settings = self.TYPES.copy()
        self.__check(settings, self.__partial_settings)
        self.__check_electrodes(settings)

    def __check(self, target: dict, settings: dict) -> dict:
        for key in [key for key in target.keys() if key in settings.keys()]:
            if isinstance(target[key], dict):
                try:
                    self.__check(target[key], settings[key])
                except TypeError as e:
                    message = '[{type}]'.format(type=target[key]) + e.message
                    raise TypeError(message)
            else:
                if not isinstance(settings[key], target[key]):
                    message = '[{key}] is not of instance {type}'.format(
                                                key=settings[key], type=[key])
                    raise TypeError(message)

    def __check_contacts(self, contacts: list) -> None:
        for index, contact in enumerate(contacts):
            try:
                self.__check(self.CONTACT_SETTING, contact)
            except TypeError as e:
                message = 'Contact {index}'.format(index=index) + e.message
                raise TypeError(message)

    def __check_electrodes(self, settings: dict) -> None:
        for index, electrode in enumerate(settings['Electrodes']):
            for key in electrode.keys():
                try:
                    if not isinstance(electrode[key],
                                      self.ELECTRODE_SETTING[key]):
                        pass
                    if key == 'Contacts':
                        self.__check_contacts(electrode['Contacts'])
                except TypeError as e:
                    message = 'Electrode {index}'.format(index=index) \
                              + e.message
                    raise TypeError(message)
