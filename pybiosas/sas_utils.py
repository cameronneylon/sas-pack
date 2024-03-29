import numpy as np

class SasData(object):
    """Root class for data object for holding 1-d Q versus I SAS data.

    Class has a series of methods for initialising and
    doing basic operations on SAS data. The root class
    is very simple and only contains the Q and I lists
    providing simple addition, multiplication, length
    and string operations.
    """

    def __init__(self, q_vals, i_vals):
        """Initializing the SasData object.

        Just passes two lists to the initializing object which
        are then stored internally. Possibly an argument for 
        including the units of Q in the root object
        """
        
        assert type(q_vals) == list
        assert type(i_vals) == list
        assert len(q_vals) == len(i_vals), 'q and i not the same length'
        self.q = q_vals
        self.i = i_vals


    def __len__(self):
        return len(self.q)


    def __add__(self, other):
        """Test whether other is SasData or float/int and add together.

        The basic add operation covers two specific cases. The first is
        where two SasData objects are added together where the wish is
        for the intensities of both to be combined. The second cases is
        when adding (or more likely subtracting) a numeric value (float 
        or int) which is handled separately. Probably both of these could
        be written neater and faster.

        Inputs must required are a SasData object and either a SasData object or
        an int or float. Returns a new SasData object.
        """

        # addition of two SasData objects
        if isinstance(other, SasData):
            assert len(self) == len(other), 'datasets not the same length'
	    assert all(self.q) == all(other.q), 'q values not the same'

            # initialize a SasData object with zero for all intensities
	    out = SasData(self.q, ([0] * len(self)))
	    for j in range(len(self)):
                out.i[j] = self.i[j] + other.i[j]

	    return out

        # addition of a float or int to  SasData objects
        elif type(other) is int or type(other) is float:
            out = SasData(self.q, ([0]*len(self)))
	    for j in range(len(self)):
	        out.i[j] = self.i[j] + other

	    return out

        else:
            return NotImplemented

    def __sub__(self, other):
        """Subtraction method for SasData Objects.

        See documentation for __add__ which is essentially identical
        """

        # addition of two SasData objects
        if isinstance(other, SasData):
            assert len(self) == len(other), 'datasets not the same length'
	    assert all(self.q) == all(other.q), 'q values not the same'

            # initialize a SasData object with zero for all intensities
	    out = SasData(self.q, ([0] * len(self)))
	    for j in range(len(self)):
                out.i[j] = self.i[j] - other.i[j]

	    return out

        # addition of a float or int to  SasData objects
        elif type(other) is int or type(other) is float:
            out = SasData(self.q, ([0]*len(self)))
	    for j in range(len(self)):
	        out.i[j] = self.i[j] - other

	    return out

        else:
            return NotImplemented

    def __mul__(self, other):
        """Basic mutplication function for SasData objects.

        Requires a SasData object and an int or a float. Returns a new
        SasData object.
        """

        assert isinstance(self, SasData)
        try:
            assert type(other) is int or type(other) is float
        except:
            return NotImplemented

        # initialize a SasData object with zero for all intensities 
        out = SasData(self.q, (([0] * len(self)))) 
        for j in range(len(self)):
            out.i[j] = self.i[j] * other
	return out


class ExpSasData(SasData):
    """Derived class for experimental SAS data obects.

    The ExpSasData class is derived from the SasDatq class. The ExpSasData
    class adds support for masks over the data to remove parts of the
    experimental pattern. The mask and the masked data can be stored with
    the experimental data. 
    """

    def __init__(self, q, i):
        """Initialization routine adds additional SasData object for the 
        masked data at self.masked"""

        SasData.__init__(self, q, i)
        self.id = ''
        self.instrument = ''
        self.mask = []
        self.masked = SasData([],[])

    def set_id(self, id):
        try:
            assert (type(id) == float) or (type(id) == str)
        except AssertionError:
            raise TypeError

        self.id = str(id)

    def get_id(self):
        return self.id

    def set_instrument(self, instrument):
        try:
            assert type(instrument) == str
            assert instrument in ['loq', 'sans2d', 'i22', 'I22']
        except AssertionError:
            raise TypeError

        self.instrument = instrument

    def get_instrument(self):
        return self.instrument     


##################################################
#
# Loaders for various SAS data formats to SasData objects
#
##################################################

def load_file():
    """A generic loader that will call specific loaders.

    The function will call tkFileDialog.askopenfilename() to 
    get a file and then will attempt to match it against known file
    types and call the correct one.
    """

    # Use tkfd to open an OS specific file dialog and get the path to the file
    path = tkfd.askopenfilename()

    DLS_I22_recogniser = 'Created at DLS-I22'
    sasxml_recogniser = 'cansas1d/1.0'

    file = open(path, 'r').readlines()

    for j in range(0,5):
        if DLS_I22_recogniser in file[j]:
            sas_data_object = load_two_column_data(path, 3)
            break
        elif sasxml_recogniser in file[j]:
            sas_data_object = loadsasxml(path)
            break

        else:
            pass

    sas_data_object = load_two_column_data(path,0)

    return sas_data_object
    

def load_two_column_data(file, rows_to_skip=0):
    """Loader for i22 two column data files.

    i22 files currently have two columns with three lines of text
    at the top. This just does a quick and dirty load of a the file
    into a SasData object. Currently setup to be called in the form
    data = load_two_column_data('file')
    """

    data = np.loadtxt(file, skiprows = rows_to_skip)

    data_q = data[:,0]
    data_i = data[:,1]
    return ExpSasData(list(data_q), list(data_i))

import xml.etree.ElementTree as ET

def loadsasxml(file):
    """Loaded for SASxml 1.0 format data.

    The loader uses xml.etree.ElementTree to parse the xml file and
    then searches through the file to find the {cansas1d/1.0}Q and
    {cansas1d/1,0}I tags and then extract the text attribute from each 
    of these. The list is then converted from text to floats and the Q
    and I lists passed to a new SasData object. Currently nothing else
    from the sas xml folder is loaded.
    """

    # Check that file is a sasxml file
    # assert (first line of file is what it should be) is True

    # Parse the xml file and find the root element
    tree = ET.parse(file)
    elem = tree.getroot()

    # return a list of all the <Q> tags and get the Q values
    q_tags = elem.getiterator("{cansas1d/1.0}Q")

    q_list = []

    for elements in q_tags:
        q_list.append(float(elements.text)) # need to convert text to float

    # then do the same for the <I> tags and values
    i_tags = elem.getiterator("{cansas1d/1.0}I")
    i_list = []

    for elements in i_tags:
        i_list.append(float(elements.text))

    # check everything is ok with q_list and i_list
    assert len(q_list) == len(i_list), 'different number of q and i values?'
    assert len(q_list) != 0, 'appear to be no q values'
    assert len(i_list) != 0, 'appear to be no i values'
    assert q_list[0] < q_list[-1], 'q values not in order?'

    # generate and return a SasData object
    return ExpSasData(q_list, i_list)


    

