

class TypeChecker:

    TYPES = {
             'CaseGrounding': {
                 'Active': bool,
                 'Current[A]': (int, float),
                 'Voltage[V]': (int, float)
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
                 'Thickness[mm]': (int, float),
                 'Material': str,
                 'MaxMeshSizeHeight': (int, float)
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
                     'MaxMeshSizeHeight': (int, float)
                     },
                 'SaveMesh': bool
             },
             'OutputPath': str,
             'RegionOfInterest': {
                 'Center': {'x[mm]': (int, float),
                            'y[mm]': (int, float),
                            'z[mm]': (int, float)},
                 'Dimension': {'x[mm]': (int, float),
                               'y[mm]': (int, float),
                               'z[mm]': (int, float)}
                 },
             'Solver': {
                 'Type': str,
                 'Preconditioner': str,
                 'PrintRates': bool,
                 'MaximumSteps': int,
                 'Precision': (int, float)
             },
             'SpectrumMode': str,
             'StimulationSignal': {
                 'Type': str,
                 'Frequency[Hz]': (int, float),
                 'PulseWidth[µs]': (int, float),
                 'PulseTopWidth[µs]': (int, float),
                 'CounterPulseWidth[µs]': (int, float),
                 'InterPulseWidth[µs]': (int, float)
             },
             'PointModel': {
                 'Pathway': {
                     'Active': bool,
                     'FileName': str
                     },
                 'Lattice': {
                     'Center': {'x[mm]': (int, float),
                                'y[mm]': (int, float),
                                'z[mm]': (int, float)},
                     'Direction': {'x[mm]': (int, float),
                                   'y[mm]': (int, float),
                                   'z[mm]': (int, float)},
                     'PointDistance[mm]': (int, float),
                     'Shape': {'x': int,
                               'y': int,
                               'z': int}
                 }
             }
            }

    ELECTRODE_SETTING = {'Name': str,
                         'PathToCustomParameters': str,
                         'Rotation[Degrees]': (int, float),
                         'Direction': {'x[mm]': (int, float),
                                       'y[mm]': (int, float),
                                       'z[mm]': (int, float)},
                         'TipPosition': {'x[mm]': (int, float),
                                         'y[mm]': (int, float),
                                         'z[mm]': (int, float)},
                         'Contacts': list,
                         }

    CONTACT_SETTING = {'Contact_ID': int,
                       'Active': bool,
                       'Current[A]': (int, float),
                       'Voltage[V]': (int, float),
                       'Floating': bool,
                       'SurfaceImpedance[Ωm]': {'real': (int, float),
                                                'imag': (int, float)}
                       }

    @classmethod
    def check(cls, settings: dict) -> None:
        cls.__check(cls.TYPES, settings)
        cls.__check_electrodes(settings)

    @classmethod
    def __check(cls, target: dict, settings: dict) -> dict:
        for key in [key for key in target.keys() if key in settings.keys()]:
            if isinstance(target[key], dict):
                try:
                    cls.__check(target[key], settings[key])
                except TypeError as e:
                    message = '[\'{key}\']'.format(key=key) + str(e)
                    raise TypeError(message) from None
            else:
                if not isinstance(settings[key], target[key]):
                    message = '[\'{key}\'] is not of instance {type}'.format(
                                                key=key, type=target[key])
                    raise TypeError(message)

    @classmethod
    def __check_contacts(cls, contacts: list) -> None:
        for index, contact in enumerate(contacts):
            try:
                cls.__check(cls.CONTACT_SETTING, contact)
            except TypeError as e:
                message = '[\'Contacts\'][{index}]'.format(index=index)
                raise TypeError(message + str(e)) from None

    @classmethod
    def __check_electrodes(cls, settings: dict) -> None:
        for index, electrode in enumerate(settings['Electrodes']):
            try:
                cls.__check(cls.ELECTRODE_SETTING, electrode)
                cls.__check_contacts(electrode['Contacts'])
            except TypeError as e:
                message = '[\'Electrodes\'][{index}]'.format(index=index)
                raise TypeError(message + str(e)) from None
