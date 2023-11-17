class TypeChecker:
    TYPES = {
        "DielectricModel": {
            "Type": str,
        },
        "Electrodes": list,
        "FailFlag": str,
        "ActivationThresholdVTA": (type(None), float),
        "Contacts": {"MaxMeshSize": float, "MaxMeshSizeEdge": float},
        "EncapsulationLayer": {
            "Thickness[mm]": (int, float),
            "Material": str,
            "MaxMeshSize": (int, float),
        },
        "EQSMode": bool,
        "FEMOrder": int,
        "MaterialDistribution": {
            "MRIPath": str,
            "DiffusionTensorActive": bool,
            "DTIPath": str,
        },
        "Mesh": {
            "LoadMesh": bool,
            "LoadPath": str,
            "MeshElementOrder": int,
            "MeshingHypothesis": {"Type": str},
            "SaveMesh": bool,
            "SavePath": str,
        },
        "OutputPath": str,
        "BrainRegion": {
            "Center": {
                "x[mm]": (int, float),
                "y[mm]": (int, float),
                "z[mm]": (int, float),
            },
            "Dimension": {
                "x[mm]": (int, float),
                "y[mm]": (int, float),
                "z[mm]": (int, float),
            },
            "Shape": str,
        },
        "Solver": {
            "Type": str,
            "Preconditioner": str,
            "PreconditionerKwargs": dict,
            "PrintRates": bool,
            "MaximumSteps": int,
            "Precision": (int, float),
        },
        "StimulationSignal": {
            "Type": str,
            "Frequency[Hz]": (int, float),
            "PulseWidth[us]": (int, float),
            "PulseTopWidth[us]": (int, float),
            "CounterPulseWidth[us]": (int, float),
            "InterPulseWidth[us]": (int, float),
            "SpectrumMode": str,
            "CutoffFrequency": float,
        },
        "PointModel": {
            "Pathway": {"Active": bool, "FileName": str},
            "Lattice": {
                "Center": {
                    "x[mm]": (int, float),
                    "y[mm]": (int, float),
                    "z[mm]": (int, float),
                },
                "Direction": {
                    "x[mm]": (int, float),
                    "y[mm]": (int, float),
                    "z[mm]": (int, float),
                },
                "CollapseVTA": bool,
                "PointDistance[mm]": (int, float),
                "Shape": {"x": int, "y": int, "z": int},
            },
        },
    }

    ELECTRODE_SETTING = {
        "Name": str,
        "Rotation[Degrees]": (int, float),
        "Direction": {
            "x[mm]": (int, float),
            "y[mm]": (int, float),
            "z[mm]": (int, float),
        },
        "TipPosition": {
            "x[mm]": (int, float),
            "y[mm]": (int, float),
            "z[mm]": (int, float),
        },
        "Contacts": list,
    }

    CONTACT_SETTING = {
        "Contact_ID": int,
        "Active": bool,
        "Current[A]": (int, float),
        "Voltage[V]": (int, float),
        "Floating": bool,
        "SurfaceImpedance[Ohmm]": {"real": (int, float), "imag": (int, float)},
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
                    message = f"['{key}']" + str(e)
                    raise TypeError(message) from None
            else:
                if not isinstance(settings[key], target[key]):
                    message = f"['{key}'] is not of instance {target[key]}"
                    raise TypeError(message)

    @classmethod
    def __check_contacts(cls, contacts: list) -> None:
        for index, contact in enumerate(contacts):
            try:
                cls.__check(cls.CONTACT_SETTING, contact)
            except TypeError as e:
                message = f"['Contacts'][{index}]"
                raise TypeError(message + str(e)) from None

    @classmethod
    def __check_electrodes(cls, settings: dict) -> None:
        for index, electrode in enumerate(settings["Electrodes"]):
            try:
                cls.__check(cls.ELECTRODE_SETTING, electrode)
                cls.__check_contacts(electrode["Contacts"])
            except TypeError as e:
                message = f"['Electrodes'][{index}]"
                raise TypeError(message + str(e)) from None
