# import numpy as np
# import matplotlib.pyplot as plt


class TimeBlock:
    def __init__(self,name,level,time,percentage):
        self.name = name
        self.level = int(level)
        self.percentage = float(percentage)
        self.time = float(time)
        self.children = []
        self.parent = []

    def AddChild(self,block):
        self.children.append(block)
        print("chile "+block.name+" added to "+self.name)

    def AddParent(self,block):
        self.parent = block
        print("parent "+block.name+" added to "+self.name)

def ReadProfiler(folder,profName):
    # opent the profiler data
    filename = folder +"/"+ profName + "_time.csv"
    #initialize the list of first level blocks
    profile = []
    current = None # = TimeBlock("root",0,0.0,0.0)
    # Open and read the profiler file
    with open(filename,"r") as file:
        filecont = file.read()
        data = filecont.split("\n")
        # read every line -> this is a profiler entry
        for s in data:
            entry = s.split(";",100)
            # if the entry has some interesting numbers
            if(len(entry) > 1):
                level = int(entry[1])                
                # create a new block
                # name,level, max_time, glob_percent, mean_time_per_count, max_count
                block = TimeBlock(entry[0],entry[1],entry[2],entry[3])

                if(level == 1):
                    # if the block is level 1, simply create a new 
                    profile.append(block)
                else:
                    # go back in time
                    for i in range(level,current.level+1):
                        current = current.parent
                    # add the child - parent link
                    current.AddChild(block)
                    block.AddParent(current)
                # in any case, the current block is the new boss
                current = block
                print("current is now block "+block.name)

        # close the file
        file.close()
    return profile


profile = ReadProfiler("/scratch/ucl/tfl/tgillis/murphy_2261/prof","Navier-Stokes_128ranks");