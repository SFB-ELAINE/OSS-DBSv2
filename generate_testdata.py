from ossdbs.electrodes import AbbottStjudeDirected6172, \
                              AbbottStjudeDirected6173, \
                              AbbottStjudeActiveTip6142_6145, \
                              AbbottStjudeActiveTip6146_6149, \
                              AbbottStjudeActiveTip6146_6149, \
                              BostonScientificVerciseDirected, \
                              BostonScientificVercise, \
                              Medtronic3387, \
                              Medtronic3389, \
                              Medtronic3391, \
                              PINSMedicalL301, \
                              PINSMedicalL302, \
                              PINSMedicalL303, \
                              MicroProbesCustomRodent, \
                              MicroProbesSNEX_100               \
                            
from tests.geometry_converter import GeometryConverter  
import json
import numpy as np
import netgen





if __name__ == "__main__":
    FILE_PREFIX=""
    classes = [ (AbbottStjudeActiveTip6142_6145,"AbbottStjudeActiveTip6142_6145"), 
                (AbbottStjudeActiveTip6146_6149,"AbbottStjudeActiveTip6146_6149"), 
                (AbbottStjudeActiveTip6146_6149,"AbbottStjudeActiveTip6146_6149"),
                (BostonScientificVercise,"BostonScientificVercise"), 
                (Medtronic3387,"Medtronic3387"), 
                (Medtronic3389,"Medtronic3389"),
                (Medtronic3391,"Medtronic3391"),
                (PINSMedicalL301,"PINSMedicalL301"),
                (PINSMedicalL302,"PINSMedicalL302" ),
                (PINSMedicalL303,"PINSMedicalL303"),
                (MicroProbesCustomRodent,"MicroProbesCustomRodent"),
                (MicroProbesSNEX_100,"MicroProbesSNEX_100")]
    DIRECTED_classes=[(BostonScientificVerciseDirected,"BostonScientificVerciseDirected"),
                      (AbbottStjudeDirected6172,"AbbotStjudeDirected6172"),
                      (AbbottStjudeDirected6173,"AbbotStjudeDirected6173")]
    DEFAULT_TEST_DATA=[
        ((0.0, (0, 0, 0), (0, 0, 1)),  '_0.json'),
        ((0.0, (0, 0, 0), (0, 0, 0)),  '_0.json'),
        ((3.0, (0, 0, 0), (0, 0, 1)),  '_0.json'),
        ((0.0, (1, -2, 3), (0, 0, 1)),  '_1.json'),
        ((0.0, (1, -2, 3), (0, 0, 0)),  '_1.json'),
        ((3.0, (1, -2, 3), (0, 0, 0)),  '_1.json'),
        ((0.0, (1, -2, 3), (2.0, 0, 1.0)), '_2.json'),
        ((0.0, (1, -2, 3), (2.0/3.0, 0, 1.0/3.0)), '_2.json')
    ]
    DIRECTED_TEST_DATA=[
        ((0.0, (0, 0, 0), (0, 0, 1)),  '_0.json'),
        ((0.0, (0, 0, 0), (0, 0, 0)),  '_0.json'),
        ((3.0, (0, 0, 0), (0, 0, 1)),  '_1.json'),
        ((0.0, (1, -2, 3), (0, 0, 1)),  '_2.json'),
        ((0.0, (1, -2, 3), (0, 0, 0)),  '_2.json'),
        ((3.0, (1, -2, 3), (0, 0, 0)),  '_3.json'),
        ((0.0, (1, -2, 3), (2.0, 0, 1.0)), '_4.json'),
        ((0.0, (1, -2, 3), (2.0/3.0, 0, 1.0/3.0)), '_4.json')
    ]
    DIRECTORY="tests/test_data/"
    for class_,prefix in DIRECTED_classes: 
        for electrode_parameters,suffix in DIRECTED_TEST_DATA:
            with open(DIRECTORY+prefix+suffix,"w") as file:
                rotation, translation, direction = electrode_parameters
                electrode=class_(rotation,direction,translation).generate_geometry()
                dictionary = GeometryConverter(electrode).to_dictionary()
                json.dump(dictionary, file)
    for class_,prefix in classes: 
        for electrode_parameters,suffix in DEFAULT_TEST_DATA:
            with open(DIRECTORY+prefix+suffix,"w") as file:
                rotation, translation, direction = electrode_parameters
                electrode=class_(rotation,direction,translation).generate_geometry()
                dictionary = GeometryConverter(electrode).to_dictionary()
                json.dump(dictionary, file)
                


