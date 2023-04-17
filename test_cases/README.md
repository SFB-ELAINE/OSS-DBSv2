# Test cases

Go to the respective directory and run then

`ossdbs name_of_input_file.json` 

## Input case 1

TODO: documentation

## Input case 2

Does not work
Error:
```
importing NGSolve-6.2.2204
Traceback (most recent call last):
  File "/home/julius/Documents/gitlab/oss-dbsv2/venv/bin/ossdbs", line 33, in <module>
    sys.exit(load_entry_point('ossdbs', 'console_scripts', 'ossdbs')())
  File "/home/julius/Documents/gitlab/oss-dbsv2/ossdbs/main.py", line 21, in main
    point_analysis(input)
  File "/home/julius/Documents/gitlab/oss-dbsv2/ossdbs/point_analysis/analysis.py", line 35, in point_analysis
    electrodes = ElectrodesFactory.create(settings['Electrodes'])
  File "/home/julius/Documents/gitlab/oss-dbsv2/ossdbs/factories/electrodes_factory.py", line 20, in create
    electrode = cls.__create_electrode(index, input_par)
  File "/home/julius/Documents/gitlab/oss-dbsv2/ossdbs/factories/electrodes_factory.py", line 57, in __create_electrode
    model_parameters = cls.read_json(path)
AttributeError: 'int' object has no attribute 'read_json'
```

## Input case 3

Does not work
Error:

```
importing NGSolve-6.2.2204
Traceback (most recent call last):
  File "/home/julius/Documents/gitlab/oss-dbsv2/venv/bin/ossdbs", line 33, in <module>
    sys.exit(load_entry_point('ossdbs', 'console_scripts', 'ossdbs')())
  File "/home/julius/Documents/gitlab/oss-dbsv2/ossdbs/main.py", line 21, in main
    point_analysis(input)
  File "/home/julius/Documents/gitlab/oss-dbsv2/ossdbs/point_analysis/analysis.py", line 35, in point_analysis
    electrodes = ElectrodesFactory.create(settings['Electrodes'])
  File "/home/julius/Documents/gitlab/oss-dbsv2/ossdbs/factories/electrodes_factory.py", line 20, in create
    electrode = cls.__create_electrode(index, input_par)
  File "/home/julius/Documents/gitlab/oss-dbsv2/ossdbs/factories/electrodes_factory.py", line 57, in __create_electrode
    model_parameters = cls.read_json(path)
AttributeError: 'int' object has no attribute 'read_json'
```

## Input case 4

TODO: documentation

## Input case 5

TODO: documentation

## Input case 6 

TODO: documentation
