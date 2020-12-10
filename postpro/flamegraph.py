import argparse
from readprof import ReadProfiler
from matplotlib import pyplot as plt
from matplotlib import patches as patches

#------------------------------------------------------------------------------
# add some argument parsing
# parser = argparse.ArgumentParser(description='compute the flame graph for the result of a MURPHY profiler')
# parser.add_argument('--folder', nargs=1, help='the folder in which the profiler is stored')
# parser.add_argument('--sum', dest='accumulate', action='store_const',
#                     const=sum, default=max,
#                     help='sum the integers (default: find the max)')
# args = parser.parse_args()
# print(args.accumulate(args.integers))

def FoldBlock(file,block,string):
    # increment the string for the children
    child_string = string + block.name + ";"
    
    # loop and count how much is left
    total_child = 0
    for child in block.children:
        FoldBlock(file,child,child_string)
        total_child = total_child + child.ToMiliSecond()
    
    # print my name
    my_time = max(0,block.ToMiliSecond()-total_child)
    file.write(string + block.name + " " + str(my_time)+"\n")


#------------------------------------------------------------------------------
# read the file
# folder = "./"
# profname = "Navier-Stokes_128ranks"
folder = "../prof"
profname = "Navier-Stokes_3ranks"
profile = ReadProfiler(folder,profname)

# open the file
file = open(profname+".fold","w+")

# start the folding
for block in profile:
    print(block.name)
    FoldBlock(file,block,"")


file.close()


