# PyBioSas.modelling: A python library for interacting with SansView 
# implementation of the NIST SANS models at command line level
#
# Public Domain Waiver:
# To the extent possible under law, Cameron Neylon has waived all 
# copyright and related or neighboring rights to lablogpost.py
# This work is published from United Kingdom.
#
# See http://creativecommons.org/publicdomain/zero/1.0/
#
# Dependencies: The application requires a range of modules from the
# Python 2.6 standard library including json, sys, argparse, datetime,
# os.path. Python is Copyright 2001-2012 Python Software Foundation
# and used here under the PSF Licence for Python 2.6.x
#
# In addition it requires non-standard library elements including:
#     * Sansview (and all its dependencies)
#     * Numpy and Scipy
#     * pybiosas.sas_utils
#
# In principle these should all be installed for you if you've used
# pip or easy_install to pull this package from PyPi

import optparse
import json
import sys
import datetime
import os
import os.path
try:
    import pybiosas.sas_utils
    import pybiosas.models
except ImportError:
    import sas_utils
    import models
import scipy.optimize
import copy
import numpy as np


class ApplicationRun():
    """Command line app for running refinements and model calculations

    """

    def __init__(self):
        """Initialisation method for the app class"""

        self.parser = None
        self.command = None
        self.model = None
        self.parameters = None
        self.dataset = None
        self.datain = None
        self.outpath = None


        self.__init_parser()
        self._raw_args, command = self.parser.parse_args()
        self.args = vars(self._raw_args)
        self.args['command'] = command[0]
        
    def __init_parser(self):

        self.parser = optparse.OptionParser()
        
        # TODO System for monitoring and registering models automatically
        # The available models are just hard coded at the moment
        self._registered_models = pybiosas.models.models

        self.parser.add_option('-d', '--dataset', type = str, dest='dataset',
                                     help = "The dataset to fit in SasXML format")

        models = []
        for model in self._registered_models.iterkeys():
            models.append(model)
        self.parser.add_option('-m', '--model', type=str, dest='model',
                                 help = ("""The model to fitted or calculated.
                                 Available models are""" +
                                     str(models)))
        
        self.parser.add_option('-p', '--parameters', type = str, dest='parameters',
                                 help = """The paramaters, as either a json
                          file or a list of dictionaries with structure as defined
                          for the parinfo option of mpfit.""")

        self.parser.add_option('-q', '--q_vals', '--q_values', type = str,
                                 dest = 'q_vals',
                                 help = """A list of q values for calculating
                                 intensities from a model. Provide as a string
                                 or file containing either the list of q values
                                 [1st, 2nd...last] or [1st, last, step] where
                                 step is optional""",
                                 default = "[0.0, 0.5, 0.0005]")

        self.parser.add_option('-o', '--outpath', type = str, dest='outpath',
                                 help = """A path to a directory for the output
                                 files. If it is desired to name the output file
                                 then terminate the path with that filename.
                                 Otherwise terminate the path with a '/'.
                                 Default is to place output files in
                                 same directory as the input parameter file. The
                                 default output filename is 'outfile.json'. If
                                 parameters are passed as a string on the
                                 command line then outpath defaults to the
                                 current working directory.""")



    def execute(self):
        """Generate Model Wrapper and Execute."""

        self.model_instance = ModelWrapper(self.args)
        self.model_instance.execute()

    def write_out(self):
        self.model_instance.write()
        

class ModelWrapper:
    def __init__(self, args):

        self.outfile = 'sansmodel_output'
        self.args = args
        self.fitsuccess = False
        self.traceback = None
        self.__distribute_args()
        self._registered_models = pybiosas.models.models


    def __distribute_args(self):
        """Initialiser to distribute commandline args to internal variables"""

        for key in self.args.iterkeys():
            self.__dict__[key] = self.args[key]
        if not 'q_vals' in self.args:
            self.q_vals = None

    def calculate(self):
        """Calculate values of i for given model and q values

        The function first sets up the list of q values, either from an existing
        imported list of values, or from a list of three values (start, stop,
        numpts) or two values (start, stop, assumed 100 points). The function then
        sets the appropriate parameters for the model and evaluates the model for
        each value in q, setting self.i_vals_out and self.q_vals_out in preparation
        for writing out the results.
        """
        
        q_vals_list = self.q_vals
        if len(q_vals_list) == 2:
            q_vals = np.arange(q_vals_list[0], q_vals_list[1], 100)

        elif len(q_vals_list) == 3:
            q_vals = np.arange(q_vals_list[0], q_vals_list[1],
                                    q_vals_list[2]).tolist()

        else:
            q_vals = q_vals_list
            
        # Calculate i for each value of q
        i_vals_out = map(self.__model_func.run, q_vals)

        self.i_vals_out = i_vals_out
        self.q_vals_out = q_vals
        return True
    
    def execute(self):
        """Routine to execute the calculation or fit"""

        self.__model_func = self.__model_importer()
        self.__load_files_from_args() # load data from files

        if self.command == 'fit':
            print "Fitting"
            self.fit()
            print "Fitted:", self.fitsuccess
            print "Calculating"
            self.calculate()
            print "Calculated"
            
        elif self.command == ('calc' or 'calculate'):
            self.calculate()

        else:
            raise NotImplementedError, "arguments not passing properly"

    def fit(self):
        """Run the fit process for the given model

        The parameters list is set up from those parameters that are not set as
        fixed in self.parameters. As these parameters have already been set in
        the model in the __load_args function they do not need to be set again here.
        If this library is being used in scripts it might be appropriate to reset
        parameters for the model here just to be safe.
        """
        
        parameters=[]
        self.q_vals = self.datain.q
        for par in self.parameters:
            if not par.get('fixed', False):
                parameters.append(Parameter(self.__model_func, par['paramname'],
                                                               value=par['value']))
            
        def f(params):
            i=0
            for p in parameters:
                p.set(params[i])
                i+=1

            residuals=[]
            for j in range(len(self.datain.q)):
                residuals.append(self.datain.i[j] -
                             self.__model_func.run(self.datain.q[j]))

            return residuals

        def chi2(params):
            sum = 0
            res = f(params)
            for item in res:
                sum += item * item
            return sum

        p = [param() for param in parameters]
        out, self.cov_x, self.fit_info, self.mesg, success = scipy.optimize.leastsq(f, p, 
                                                                 full_output=1,
                                                                  maxfev = 1000*len(p))
        # Calculate chi squared
        if len(parameters) > 1:
            self.chisqr = chi2(out)
        elif len(parameters) == 1:
            self.chisqr = chi2([out])
        
        # Update the main parameter list at self.parameters with finalised values
        paramlist = []
        for p in self.parameters:
            paramlist.append(p['paramname'])

        for p in parameters:
            self.parameters[paramlist.index(p.get_name())]['value'] = p.get()

        if self.cov_x != None:
            self.fitsuccess = True


    def write(self):
        outdict = {'model'            : self.model,
                   'data_out'         : {'q'       : json.dumps(self.q_vals_out),
                                         'i'       : json.dumps(self.i_vals_out),
                                         'units'   : 'A^-1'},
                   'run'              : {'command' : self.command,
                                         'date'    : str(datetime.date.today()),
                                         'time'    : str(datetime.time())},
                   'parameters_in'    : self.parameters_in}

        if self.dataset:
            outdict['dataset'] = {'q_in' : json.dumps(self.datain.q),
                                  'i_in' : json.dumps(self.datain.i)}
            
        if (self.fitsuccess and (self.command == 'fit')):
            outdict['fit'] = {'chi2'           : {'value' : self.chisqr},
                              'cov_x'          : {'value' : str(self.cov_x)}}

            for param in self.parameters:
                outdict['fit'][param['paramname']] = param

        if os.path.isdir(self.outpath):
            self.outpath = os.path.join(self.outpath, "sansmodel_output.json")
            
        path, filename = os.path.split(self.outpath)
        print "Path:", path, "Filename:", filename
             
        if not os.path.exists(path):
            os.mkdir(path)

        f = open(self.outpath, 'w')
        json.dump(outdict, f)
        f.close()


    def __load_files_from_args(self):
        """Load files in based on command arguments

        """

        if self.dataset:
            try:
                print "trying", self.dataset
                if os.path.splitext(self.dataset)[1] == '.xml':
                    self.datain = pybiosas.sas_utils.loadsasxml(self.dataset)
                else:
                    self.datain = pybiosas.sas_utils.load_two_column_data(self.dataset, rows_to_skip=1)

                self.datain.err = None
    
            except OSError:
                errmsg = "Unable to load file: " + self.dataset
                raise InputError, errmsg
                return False

        # Load and parse the parameters
        if self.parameters:
            if os.path.isfile(self.parameters):
                print "self.outpath:", self.outpath
                if not self.outpath:
                    if os.path.dirname(self.parameters):
                        self.outpath = os.path.dirname(self.parameters)
                    else:
                        self.outpath = os.getcwd()
                    
                f = open(self.parameters, 'r')
                self.parameters = json.load(f)
                f.close()
            else:
                if not self.outpath:
                    self.outpath = os.getcwd()
                self.parameters = json.loads(self.parameters)

            self.parameters_in = copy.deepcopy(self.parameters)
            
        # Check we have all the parameters for model including those left as defaults
        input_param_list = []
        for param in self.parameters:
            input_param_list.append(param['paramname'])
        for paramname in self.__model_func.details.keys():
            if (paramname not in input_param_list) and (
                    paramname not in self.__model_func.orientation_params):
                self.parameters.append({'paramname' : paramname,
                                        'value' :
                                       self.__model_func.getParam(paramname)
                                        })

        # Set the parameters for the model
        for parameter in self.parameters:
            self.__model_func.setParam(parameter['paramname'], parameter['value'])
        
        if (self.q_vals and os.path.isfile(self.q_vals)):
            try:
                f = open(self.q_vals, 'r')
                q_in = json.load(f)
            except OSError:
                errmsg = "Failed to load q_vals in: " + self.q_vals
                raise InputError, errmsg
                return False
            self.q_vals = q_in

        elif self.q_vals:
            q_in = json.loads(self.q_vals)
            self.q_vals = q_in

        else:
            pass

        
            
    def __model_importer(self):
        """Function for importing correct model library and returning model function

        Because we don't know either the module name or the model function name
        until runtime (as it is selected by the user) we need to look it up in
        the dictionary of registered model (self._registered_models) and then
        use __import__ so that we can pass the import command the string with
        the module location. Having identified the library we then want to pass
        the actual model calculation function back to the data model. The
        function is available from the module __dict__ and if we know the model
        we know the function name, so can create a pointer to this function and
        return it to the main app execution thread.
        """
        
        model_location = self._registered_models[self.model]['library_name']
        library_location = 'sans.models.' + model_location
        __import__(library_location)
        model_func = getattr(sys.modules[library_location],
                             self._registered_models[self.model]['model_name'])()
        return model_func

class Parameter:
    """
    Convenience class to handle model parameters for a fit
    """

    def __init__(self, model, name, value=None, fixed=False):
            self.model = model
            self.name = name
            self.fixedparam = fixed
            if not value == None:
                self.model.setParam(self.name, value)
           
    def set(self, value):
        """
            Set the value of the parameter
        """
        self.model.setParam(self.name, value)

    def get(self):
        """ 
            Return the current value of the parameter
        """
        return self.model.getParam(self.name)

    def get_name(self):
        return self.name

    def fixed(self):
        """
            Return the value of the "fixedparam" internal variable

            self.fixedparam is used to toggle whether this parameter should be fixed
            in the subsequent fit. If self.fixed==True then the parameter is not
            included in the list for the fit.
        """

        return self.fixedparam

    def set_fixed(self, fix=True):
        """
            Set the parameter as being fixed for fit.

            Default after being set is True.
        """
        self.fixedparam = fix

    def __call__(self):
        """ 
            Return the current value of the parameter
        """
        return self.model.getParam(self.name)
        
if __name__ == '__main__':
    run = ApplicationRun()
    run.parser.parse_args()
    run.execute()
    run.write_out()
