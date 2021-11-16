# import matplotlib

## SET BACKEND
import matplotlib

# matplotlib.use('PS')

import matplotlib.pyplot as plt

import numpy as np
import sys, os

def readtxt(txtfile="untitled.txt", epsilon=0.001, ui=False):
    """
	reads the txtfile <txtfile> and
	cuts off values lower than <epsilon> to epsilon (this is practical for logy plots)
	ui = 

	returns:
	header, data
	
	where
	header is a list of str (col entries, like 'time', ...)
	data is a dict with the keys = header and the values of the tsd files
	"""

    TXT_FILE = txtfile
    file = open(TXT_FILE, "r")
    file_str = file.read()
    tuples = file_str.strip().split("\n")

    # generate the header
    if ui:
        # remove the '#' from the header
        header_str = tuples[0][1:].replace("[", "").replace("]", "")
    else:
        header_str = tuples[0].replace("[", "").replace("]", "")
    header = header_str.strip().split("\t")
    # remove the "" at the beginning and end of the strings
    header = ["time" if title == "Time" else title for title in header]

    # check if time in header
    assert( "time" in  header ),  f"TSDReader: no time/Time in header: {tuples[0]}"

    # extract the data
    data = {}
    for c in range(len(header)):
        key = header[c]

        data[key] = []
        for i in range(1, len(tuples)):
            filteredstr = tuples[i].split("\t")[c].strip().strip(")")
            if not filteredstr == "":
                data[key] += [max(float(filteredstr), epsilon)]

    return header, data


def readtsd(tsdfile="run-1.tsd", epsilon=0.001):
    """
	reads the tsdfile <tsdfile> and
	cuts off values lower than <epsilon> to epsilon (this is practical for logy plots)

	returns:
	header, data
	
	where
	header is a list of str (col entries, like 'time', ...)
	data is a dict with the keys = header and the values of the tsd files
	"""

    TSD_FILE = tsdfile
    file = open(TSD_FILE, "r")
    file_str = file.read()
    tuples = file_str.split("),(")

    # generate the header
    # remove the '((' from the header
    header_str = tuples[0][2:]
    header = header_str.split(",")
    # remove the "" at the beginning and end of the strings
    header = [title[1:-1] for title in header]

    # extract the data
    data = {}
    for c in range(len(header)):
        key = header[c]

        data[key] = []
        for i in range(1, len(tuples)):
            data[key] += [
                max(float(tuples[i].split(",")[c].strip().strip(")")), epsilon)
            ]

    return header, data


def writetsd(header, data, tsdfile="new-1.tsd"):

    file = open(tsdfile, "w")
    mystr = "("
    mystr += str(tuple(header)).replace(" ", "").replace("'", '"')

    for i in range(len(data["time"])):
        mystr += ","
        d = []
        for key in data.keys():
            d += [data[key][i]]
        mystr += str(tuple(d)).replace(" ", "")

    mystr += ")"

    file.write(mystr)
    file.close()


def writetxt(header, data, txtfile="unititled.txt"):

    file = open(txtfile, "w")
    if not os.path.isfile(txtfile):
        print('TSDReader cant find file!', file=sys.stderr) 

    mystr = ""

    # first line
    for title in header:
        if title == "time":
            mystr += "Time"
        else:
            mystr += "[" + title + "]"
        mystr += "\t"
    mystr += "\n"

    # remaining lines
    for t in range(len(data["time"])):
        for key in data.keys():
            mystr += str(data[key][t])
            mystr += "\t"
        mystr += "\n"

    file.write(mystr)
    file.close()


def appenddata(data1, data2, timeoffset):

    data2_shifted = data2
    data2_shifted["time"] = [t + timeoffset for t in data2_shifted["time"]]

    data = data1
    for key in data:
        data[key] += data2_shifted[key]

    return data
