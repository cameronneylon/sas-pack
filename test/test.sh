# Fit a model to a file of test data given a dictionary of starting parameters on the command line.
# Parameters are REQUIRED to be a dictionary that matches a subset of the parameters of the relevant
# SansView model (i.e. that returned by model.params)
echo 'Testing a fit to a model with parameters passed on command line'
python ../pybiosas/modelling.py fit -m cylinder -d testdata.xml -p '[{"fixed": true, "value": 2.0, "paramname": "scale"}, {"value": 40.0, "paramname": "radius"}, {"fixed": false, "value": 200.0, "paramname": "length"}, {"fixed": true, "value": 6e-06, "paramname": "sldSolv"}]'

# Fit a model to a file of test data given a file of starting parameters
echo 'Testing a fit to a model with a parameter file'
python ../pybiosas/modelling.py fit -m cylinder -d testdata.xml -p test.json 

# Fit a model to a file of test data given a file of starting parameters, write to new directory/file
echo 'Testing a fit to a model with a parameter file'
python ../pybiosas/modelling.py fit -m cylinder -d testdata.xml -p test.json -o "test/testout.json"

# Fit a model to a file of test data given a file of starting parameters, write to new directory without filename
echo 'Testing a fit to a model with a parameter file'
python ../pybiosas/modelling.py fit -m cylinder -d testdata.xml -p test.json -o "test/"



