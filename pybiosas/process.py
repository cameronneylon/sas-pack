import json
import os
import os.path
import optparse
import models
REGISTERED_MODELS = models.models

def compare_chi2(first_fit, second_fit):
    diff = first_fit[0]*1e10 - second_fit[0]*1e10
    if diff < 0: # Put chi2 as first parameter
        return -1
    elif diff == 0:
        return 0

    else:
        return 1

def approx_equal(x, y, tol=1e-18, rel=1e-7):
    if tol is rel is None:
        raise TypeError('cannot specify both absolute and relative errors are None')
    tests = []
    if tol is not None: tests.append(tol)
    if rel is not None: tests.append(rel*abs(x))
    assert tests
    return abs(x - y) <= max(tests)

parser = optparse.OptionParser()
parser.add_option('-p', '--path', dest = 'path', type = str,
                  help = "Path to a folder of data to process")

(options, args) = parser.parse_args()
if options.path:
    directory = options.path
else:
    directory = raw_input("Path to data?")

fit_list = []

for file in os.listdir(directory):
    if file[len(file)-4:len(file)] == 'json':

        fit = []
        path = os.path.join(directory, file)
        print directory, file
        with open(path, 'r') as f:
            output = json.load(f)

        if 'fit' in output:
            model = output['model']
            params = [p['paramname'] for p in REGISTERED_MODELS[model]['exp_vals']]
            params.insert(0, 'chi2')

            for param in params:
                param_value = float(output['fit'][param]['value'])
                fit.append(param_value)
            fit.append(file)
            fit_list.append(fit)

fit_list.sort(compare_chi2)

print params
for fit in fit_list:
    print fit

concatenated = [fit_list[0]]
count = 1
i=0
while i < len(fit_list)-1:
    if approx_equal(fit_list[i][0], fit_list[i+1][0], rel=1e-6):
        count+=1
    else:
        fit_list[i].append(count)
        concatenated.append(fit_list[i+1])
        count = 1
    i+=1

print '\n\n'
print params
for fit in concatenated:
    print fit
    



        
        

