"""
Purpose: processing txt file output from real time qPCR run (SDS 2.4) into 
input for R or DART-PCR
"""

import csv

def getInput():
    """
    Get file names
    """
    infile_name = input("Enter input file name: ")
    print()
    outfile_name = input("Enter output file name: ")
    print()
    print("There are two types of input/output options right now.")
    print("1. input = results table .txt file / output = R-friendly .csv file inlcuding columns for sample name, reporter, and Ct")
    print("2. input = multicomponent .txt file / output = .csv file for DART-PCR input")
    response = input("Which input/output option do you want (1 or 2)? ")
    return infile_name, outfile_name, response


def readData(infile_name):
    """
    Reads txt file and stores as a list of lists (lines)
    
    The first (0) item in each list containing data is a single integer (e.g. 1)
    0 = well number
    1 = sample label (e.g. D1-S1-1)
    2 = reporter name (e.g. GAPDH)
    5 = Ct
    
    item 8 in the bigger list (line number 8) contains the headers
    """
    lines = []
    infile = open(infile_name, "r")
    for l in infile:
        line = l.rstrip("\n").split("\t")
        lines.append(line)
    infile.close()
    return lines


def whichOption(lines, response):
    """
    Direct the lines to the correct selectData function depending on selected 
    option (1 or 2) and get correct sublines output
    
    This essentially makes it so that I don't have to write separate main
    code for each of these options, just different selectData functions
    """
    if response == "1":
        selectData1(lines)
    elif response == "2":
        selectData2(lines)
    else:
        print("That was an invalid option. Please enter 1 or 2 next time.")


def selectData1(lines):
    """
    Select correct elements of results table .txt file (option 1) and prepare 
    output for getOutput function
    
    Used well number to "call" rows of data and ignore other lines
    make new list of lists only containing the data I want:
        -sample label
        -reporter name
        -Ct
    """
    sublines = []
    sublines.append(["sample_id", "gene", "Ct"])
    for l in lines:
        subline = []
        # this try/except loop bypasses the issue of trying to convert strings
        # into integers, since eval wasn't working for some reason
        try:
            int(l[0])
        except ValueError:
            pass
        else:
            l[0] = int(l[0])
        if isinstance(l[0], int) == True:
            subline.append(l[1])
            subline.append(l[2])
            subline.append(l[5])
            sublines.append(subline)
        else:
            pass
    return sublines

# will need to add something in here to account for undetermined and/or
# outlier data to remove
# i would also like to do averages of remaining technical replicates so there's
# only one measure per individual per primer set
# and then calculate the ratio between 


def selectData2a(lines):
    """
    Select correct elements of multicomponent .txt file (option 2) and prepare
    output for getOutput function
    
    Implementation notes: 
    
    The ThermoFisher qPCR tech folks said that if I want fluorescence by cycle
    number per well/sample, then I need to use the data in the SYBR column of 
    the multicomponent file exported from SDS 2.4. I also need to select 
    either the average of the values at the 72 degrees stage or just one 
    (e.g. the last one) consistently for each sample
    
    For each well/sample and cycle, I want to take the average of the three 
    values in the SYBR column at ~72 degrees (Temp) for each cycle
    
    The output that I want is one row for each cycle and one column for each
    well/sample
    
    I think I can do this relatively simply if I use nested for loops where I 
    have a for loop iterating through the sample numbers (1-384, easier to just
    do all and filter out the wells I didn't actually use later; there might 
    be a way of making sure this part is more automated like doing max well 
    number and add one for the range function to make sure all wells are 
    included no matter how many there are...) and then another loop inside that
    iterating through the lines to find the rows at ~72 degrees. Then I need 
    to store the average SYBR measurements with well number and cycle number.
    
    That maintains the technical replicates as separate samples, so some
    manual filtering will be needed after making this output
    
    There are also some other data in here that might interfere with my ability
    to just use well numbers as a filter, so I may have to use in combo with
    temp data (like Well == 1 and Temp == 72). I also have to use the 
    cycle column (0, 1, 2), because there are 72 temps that I don't want in the 
    dissociation (2) phase (0 is start up, 1 is cycles 1-40).
    
    Apparently it's not 3 lines for all wells. Some of the later wells
    (e.g. 384) in the dataset that I'm using have 3 lines at the 
    elongation phase, so I guess it's a good thing I'm doing it by temp
    and not number of lines, because that should automatically mean
    that all rows at ~72 degrees are included regardless of number of rows.
    """
    # I really should turn max_well_n and max_cycle_n parts into separate functions...
    
    # first need to determine what the max well number is
    # this is column zero for each line in lines
    well_n = list()
    for line in lines:
        try:
            int(line[0])
        except ValueError:
            pass
        except TypeError: 
            pass
        else:
            well_n.append(int(line[0]))
    max_well_n = max(well_n)
    
    # then need to determine what the max cycle number is
    # this is column 5 (repeats) for each line in lines
    cycle_n = list()
    for line in lines:
        try:
            int(line[5])
        except ValueError:
            pass
        except TypeError: 
            pass
        # needed to add this because some lines don't have column 5
        except IndexError:
            pass
        else:
            cycle_n.append(int(line[5]))
    max_cycle_n = max(cycle_n)
    print(max_cycle_n)
    
    # this part is to iterate through each well/sample up to the max number
    # of wells in your file (which may or may not be 384) and retrieve
    # fluorescence data from the correct column/rows (3 rows in SYBR column)
    # and cycle number for output purposes
    
    # make empty sublines list for output
    sublines = list()
    # make empty header list and then quickly iterate through well numbers
    # to add them to header list
    header = ["cycle"]
    for well in range(1, max_well_n + 1):
        header.append(well)
    # add header list to sublines
    sublines.append(header)
    # iterate through wells to add fluorescence measurements
    for well in range(1, max_well_n + 1):
        for line in lines:
            # ignore lines that don't start with a valid well number (not headers)
            try:
                int(line[0])
            except ValueError: 
                pass
            else:
                # for each row matching well number in elongation phase of cycles
                # line 0 = well, line 3 = cycle (section containing repeats), 
                # line 2 = temp (elongation phase at 72 degrees)
                if eval(line[0]) == well and line[3] == "1" and float(line[2]) >= 71 and float(line[2]) <= 73:
                    # if I print line here, then what I should get is 40*3 lines per well (3 lines per cycle per well)
                    # what I want to do at this point is, for each well, is calculate the average across the three or four
                    # rows for each cycle and well then store this number along with cycle number and well number
                    # in the correct format for getOutput
                    # I should use a list that I add to and then take the average of for the fluorescence measurements
                    # SYBR is line[6], well number is line[0], and cycle number is line[5]
                    # that means that I need to reset the list after each well number or cycle number by well number?
                    # maybe I have to iterate through each cycle number within each well for the averages
                    for cycle in range(1, max_cycle_n + 1):
                        # empty list that will be filled with fluorescence measurements to be averaged per sample and cycle
                        flo = list()
                        # for each row/line matching cycle number
                        # int(line[5]) might be an issue if it's still looking at the lines without a column 5
                        if int(line[5]) == cycle:
                            # add fluorescence measurement to flo list
                            flo.append(line[6])
                        else:
                            pass
                else:
                    pass
    

# need to make sure I get the right sublines output from this function
# which means I need a header row that is ("cycle", 1, 2, 3...) where the
# numbers are sample numbers and then all other rows (lists within list)
# start with the cycle number and then have the fluorescence/SYBR measurements
# for each sample after that
# because each row is based on cycle number, not well number, maybe it would 
# make more sense to do cycle number first then well number


# making new version of selectData2 function below to try this
def selectData2(lines):
    """
    Select correct elements of multicomponent .txt file (option 2) and prepare
    output for getOutput function
    
    Implementation notes: 
    
    The ThermoFisher qPCR tech folks said that if I want fluorescence by cycle
    number per well/sample, then I need to use the data in the SYBR column of 
    the multicomponent file exported from SDS 2.4. I also need to select 
    either the average of the values at the 72 degrees stage or just one 
    (e.g. the last one) consistently for each sample
    
    For each well/sample and cycle, I want to take the average of the three 
    values in the SYBR column at ~72 degrees (Temp) for each cycle
    
    The output that I want is one row for each cycle and one column for each
    well/sample
    
    I think I can do this relatively simply if I use nested for loops where I 
    have a for loop iterating through the sample numbers (1-384, easier to just
    do all and filter out the wells I didn't actually use later; there might 
    be a way of making sure this part is more automated like doing max well 
    number and add one for the range function to make sure all wells are 
    included no matter how many there are...) and then another loop inside that
    iterating through the lines to find the rows at ~72 degrees. Then I need 
    to store the average SYBR measurements with well number and cycle number.
    
    That maintains the technical replicates as separate samples, so some
    manual filtering will be needed after making this output
    
    There are also some other data in here that might interfere with my ability
    to just use well numbers as a filter, so I may have to use in combo with
    temp data (like Well == 1 and Temp == 72). I also have to use the 
    cycle column (0, 1, 2), because there are 72 temps that I don't want in the 
    dissociation (2) phase (0 is start up, 1 is cycles 1-40).
    
    Apparently it's not 3 lines for all wells. Some of the later wells
    (e.g. 384) in the dataset that I'm using have 3 lines at the 
    elongation phase, so I guess it's a good thing I'm doing it by temp
    and not number of lines, because that should automatically mean
    that all rows at ~72 degrees are included regardless of number of rows.
    """
    # I really should turn max_well_n and max_cycle_n parts into separate functions...
    
    # first need to determine what the max well number is
    # this is column zero for each line in lines
    well_n = list()
    for line in lines:
        try:
            int(line[0])
        except ValueError:
            pass
        except TypeError: 
            pass
        else:
            well_n.append(int(line[0]))
    max_well_n = max(well_n)
    
    # then need to determine what the max cycle number is
    # this is column 5 (repeats) for each line in lines
    cycle_n = list()
    for line in lines:
        try:
            int(line[5])
        except ValueError:
            pass
        except TypeError: 
            pass
        # needed to add this because some lines don't have column 5
        except IndexError:
            pass
        else:
            cycle_n.append(int(line[5]))
    max_cycle_n = max(cycle_n)
    print(max_cycle_n)
    
    # this part is to iterate through each well/sample up to the max number
    # of wells in your file (which may or may not be 384) and retrieve
    # fluorescence data from the correct column/rows (3 rows in SYBR column)
    # and cycle number for output purposes
    
    # make empty sublines list for output
    sublines = list()
    # make empty header list and then quickly iterate through well numbers
    # to add them to header list
    header = ["cycle"]
    for well in range(1, max_well_n + 1):
        header.append(well)
    # add header list to sublines
    sublines.append(header)
    # iterate through wells to add fluorescence measurements
    for well in range(1, max_well_n + 1):
        for line in lines:
            # ignore lines that don't start with a valid well number (not headers)
            try:
                int(line[0])
            except ValueError: 
                pass
            else:
                # for each row matching well number in elongation phase of cycles
                # line 0 = well, line 3 = cycle (section containing repeats), 
                # line 2 = temp (elongation phase at 72 degrees)
                if eval(line[0]) == well and line[3] == "1" and float(line[2]) >= 71 and float(line[2]) <= 73:
                    # if I print line here, then what I should get is 40*3 lines per well (3 lines per cycle per well)
                    # what I want to do at this point is, for each well, is calculate the average across the three or four
                    # rows for each cycle and well then store this number along with cycle number and well number
                    # in the correct format for getOutput
                    # I should use a list that I add to and then take the average of for the fluorescence measurements
                    # SYBR is line[6], well number is line[0], and cycle number is line[5]
                    # that means that I need to reset the list after each well number or cycle number by well number?
                    # maybe I have to iterate through each cycle number within each well for the averages
                    for cycle in range(1, max_cycle_n + 1):
                        # empty list that will be filled with fluorescence measurements to be averaged per sample and cycle
                        flo = list()
                        # for each row/line matching cycle number
                        # int(line[5]) might be an issue if it's still looking at the lines without a column 5
                        if int(line[5]) == cycle:
                            # add fluorescence measurement to flo list
                            flo.append(line[6])
                        else:
                            pass
                else:
                    pass
    



def getOutput(outfile_name, sublines):
    outfile = open(outfile_name, "w")
    with outfile:
        writer = csv.writer(outfile)
        writer.writerows(sublines)


def main():
    infile_name, outfile_name, response = getInput()
    lines = readData(infile_name)
    whichOption(lines, response)
    #print(well_n.__dict__)
    #print(max_well_n)
    #print(well_n)
    #sublines = selectData(lines)
    #getOutput(outfile_name, sublines)

if __name__ == "__main__":
    main()


