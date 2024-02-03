#! /bin/bash

cd input_case1
echo $PWD
ossdbs input_homogeneous.json
ossdbs input_inhomogeneous.json
cd ..

cd input_case2
echo $PWD
ossdbs input_custom_electrode.json
ossdbs input_custom_material.json
cd ..


cd input_case3
echo $PWD
ossdbs input_case_grounding.json
cd ..

cd input_case4
echo $PWD
ossdbs input_current_controlled.json
ossdbs input_multi_current.json
cd ..

cd input_case5
echo $PWD
ossdbs input_stimulation_signal.json
cd ..

cd input_case6
echo $PWD
ossdbs input_floating.json
cd ..

cd input_case7
echo $PWD
ossdbs input_vta.json
ossdbs input_vta_out_of_core.json
cd ..

cd input_case8
echo $PWD
ossdbs input_pathway.json
ossdbs input_pathway_out_of_core.json
cd ..
