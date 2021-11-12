import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
import tikzplotlib
from readprof import ReadProfiler

cmap10 = matplotlib.cm.get_cmap('tab10')
cmap20 = matplotlib.cm.get_cmap('tab20')
cmap20b = matplotlib.cm.get_cmap('tab20b')
cmap20c = matplotlib.cm.get_cmap('tab20c')

col = dict()
col['sten'] = cmap20c(0)
col['ghost'] = cmap20b(4)
col['comp'] = cmap20c(2)
col['adapt'] = cmap20c(4)
col['aghost'] = cmap20c(7)
col['acrit'] = cmap20c(6)
col['aset'] = cmap20c(5)

col_pie = dict()
col_pie['sten'] = cmap20c(3)
col_pie['ghost'] = cmap20b(7)
col_pie['adapt'] = cmap20c(7)
col_bar = dict()
col_bar['sten'] = cmap20c(1)
col_bar['ghost'] = cmap20b(5)
col_bar['adapt'] = cmap20c(5)



def Block2Dict(dict, block, name):
    # add myself
    dict[name + block.name] = block.time
    # go to my children
    child_name = name + block.name + ";"
    for child in block.children:
        Block2Dict(dict, child, child_name)


# folder = "/Users/tgillis/Desktop/weak-nic5/"
# folder = "/Volumes/Thomas_MIT_backup/2021_10_ercap_scaling/weak_nic5_2021-10-01-0637-5189f74c"
folder = "/Users/tgillis/Desktop/weak_nic5_2021-10-01-0637-5189f74c"
folder = "/Volumes/Thomas_MIT_backup/2021_10_ercap_scaling/weak_nic5_2021-10-01-0637-5189f74c"

overleaf = "/Users/tgillis/Dropbox/research/Overleaf/2021_10_ercap/figures"

ranks = np.arange(0, 37, 4)*32
ranks[0] = 32
nodes = ranks/64
nodes[0] = 1

time_dodt = np.zeros(len(ranks))
time_stencil = np.zeros(len(ranks))
time_stencil_ghost = np.zeros(len(ranks))
time_stencil_inner = np.zeros(len(ranks))
time_stencil_outer = np.zeros(len(ranks))
time_comp = np.zeros(len(ranks))
time_ghost = np.zeros(len(ranks))
time_rma_wait = np.zeros(len(ranks))
time_p4est = np.zeros(len(ranks))
time_adapt = np.zeros(len(ranks))
time_adapt_ginit = np.zeros(len(ranks))
time_adapt_setup = np.zeros(len(ranks))
time_adapt_smooth = np.zeros(len(ranks))
time_adapt_criterion = np.zeros(len(ranks))
time_adapt_ghost = np.zeros(len(ranks))
time_adapt_setupghost = np.zeros(len(ranks))
time_adapt_destroyghost = np.zeros(len(ranks))
time_adapt_misc = np.zeros(len(ranks))

time_ghost_RMA_PS= np.zeros(len(ranks))
time_ghost_RMA_getput= np.zeros(len(ranks))
time_ghost_RMA_CW= np.zeros(len(ranks))
time_ghost_RMA_comp= np.zeros(len(ranks))
time_adapt_partition= np.zeros(len(ranks))
time_adapt_adapt= np.zeros(len(ranks))
time_adapt_misc= np.zeros(len(ranks))

for irun in np.arange(len(ranks)):
    # # get the folder
    # print(folder+"/murphy_weak_" + str(nodes[irun]) + "_*")
    # runs = glob.glob(folder+"/prof/WeakScalability_" + str(ranks[irun]) + "ranks*")
    # assert(len(runs) == 1)
    # run = runs[0]
    rank = ranks[irun]

    # read the profile
    # print("reading from run: "+run)
    profile = ReadProfiler(folder + "/prof", "WeakScalability_" + str(rank) + "ranks_w40")

    # get my dict out of it
    mydic = {}
    for block in profile:
        Block2Dict(mydic, block, "")

    # print(mydic)
    ghost_tag_1 = "run;do dt;rk3;rhs;stencil;ghost;pullghost post (2,3);"
    ghost_tag_2 = "run;do dt;rk3;rhs;stencil;ghost;pullghost wait (2,3);"
    time_ghost_RMA_PS[irun] = mydic[ghost_tag_1+"(02) RMA - post start"] + mydic[ghost_tag_2+"(08) RMA - post start"]
    time_ghost_RMA_getput[irun] = mydic[ghost_tag_1+"(03) RMA - get"] + mydic[ghost_tag_2+"(09) RMA - put"]
    time_ghost_RMA_CW[irun] = mydic[ghost_tag_2+"(05) RMA - complete"] + mydic[ghost_tag_2+"(06) RMA - wait"] + mydic[ghost_tag_2+"(10) RMA complete"] + mydic[ghost_tag_2+"(11) RMA wait"]  
    time_ghost_RMA_comp[irun] = mydic[ghost_tag_1+"(01) push to window"]+mydic[ghost_tag_1+"(04) computation"] + mydic[ghost_tag_2+"(07) computation"]+ mydic[ghost_tag_2+"(12) pull from window"] + + mydic[ghost_tag_2+"(13) computation"]
    

    # print(mydic["Navier-Stokes run;do dt"])
    time_dodt[irun] = mydic["run;do dt"]
    time_stencil[irun] = mydic["run;do dt;rk3;rhs;stencil"]
    time_stencil_ghost[irun] = mydic["run;do dt;rk3;rhs;stencil;ghost"]
    time_stencil_inner[irun] = mydic["run;do dt;rk3;rhs;stencil;inner"]
    time_stencil_outer[irun] = mydic["run;do dt;rk3;rhs;stencil;outer"]
    time_ghost[irun] = mydic["run;do dt;rk3;rhs;stencil;ghost"]
    time_comp[irun] = mydic["run;do dt;rk3;rhs;stencil;inner"]+mydic["run;do dt;rk3;rhs;stencil;outer"]
    
    # time_rma_wait[irun] = mydic["run;do dt;rk3;rhs;stencil;ghost;pullghost wait;RMA wait"]
    time_adapt[irun] = mydic["run;adapt;adaptation"]
    time_adapt_criterion[irun] = mydic["run;adapt;adaptation;criterion"]
    time_adapt_ghost[irun] = mydic["run;adapt;adaptation;ghost for criterion"]
    time_adapt_partition[irun] = mydic["run;adapt;adaptation;partition comm"] +  mydic["run;adapt;adaptation;partition init"]
    # time_adapt_adapt[irun] = mydic["run;adapt;adaptation;solve dependency"]
    time_adapt_setupghost[irun] = mydic["run;adapt;adaptation;setup mesh and ghost"]
    time_adapt_destroyghost[irun] = mydic["run;adapt;adaptation;destroy mesh and ghost"]
    time_adapt_smooth[irun] = mydic["run;adapt;adaptation;smooth jump"] + mydic["run;adapt;adaptation;solve dependency"]
    time_adapt_misc[irun] = time_adapt[irun] - (time_adapt_smooth[irun] + time_adapt_criterion[irun] +time_adapt_ghost[irun] + time_adapt_partition[irun]+time_adapt_setupghost[irun]+time_adapt_destroyghost[irun] )
    # print(time_adapt[irun])
    # print(time_adapt_criterion[irun] + time_adapt_ghost[irun] + time_adapt_partition[irun] + time_adapt_adapt[irun]+time_adapt_setupghost[irun]+time_adapt_destroyghost[irun])
    # time_adapt_misc[irun] = time_adapt[irun] -(time_adapt_criterion[irun] + time_adapt_ghost[irun] + time_adapt_partition[irun] + time_adapt_adapt[irun]+time_adapt_setupghost[irun]+time_adapt_destroyghost[irun])
    # time_adapt_ginit[irun] = mydic["run;adapt;adaptation;ghost_init"]
    # time_p4est[irun] = mydic["run;adapt;adaptation;p4est balance"] + mydic["run;adapt;adaptation;p4est refine"] + mydic["run;adapt;adaptation;p4est coarsen"]


# fig, ax = plt.subplots(num="pie stencil")
# time_
# ax.pie()

nbars = 7
width = 0.5/nbars*32.0*4.0


#==============================================================================================================================================================
expl = 0.075
pie_rank_list = [0, 1, 9]
pad = 1.01

for pie_rank in pie_rank_list:

    print("--------------------------------------")
    print(ranks[pie_rank])
    print("--------------------------------------")
    fig, ax = plt.subplots(num="pie stencil"+str(ranks[pie_rank]))
    # ratio_stencil = np.array([time_stencil_ghost[pie_rank],time_stencil_inner[pie_rank],time_stencil_outer[pie_rank]])/time_stencil[pie_rank]
    ratio_stencil = np.array([time_stencil_ghost[pie_rank],time_stencil_inner[pie_rank],time_stencil_outer[pie_rank]])
    print("--------------------------------------")
    print(np.sum(ratio_stencil))
    print(time_stencil[pie_rank])
    print("--------------------------------------")
    data = ratio_stencil
    def absolute_value(val):
        a  = data[ np.abs(data - val/100.*data.sum()).argmin() ]
        return "{:2.2f} sec".format(a)
    cols = col_pie['sten']
    coles= col['sten']
    patches = ax.pie(ratio_stencil,explode=(expl,expl,expl),labels=("ghost","inner","outer"),normalize=True,
            autopct=absolute_value,
            colors=(cols,cols,cols,cols),
            wedgeprops={"edgecolor":coles,"linewidth":2,'antialiased': True},
            textprops={'fontsize': 14},
            pctdistance=0.65)[0]
    patches[1].set_hatch('/')
    patches[2].set_hatch('/')
    # ax.bar(xaxis - 1.0 * width, time_stencil,width= width,  label="ghost",color=col['sten'])
    # ax.bar(xaxis - 0.0 * width, time_ghost,width= width,   label="stencil - ghosting",color=col['ghost'])
    # ax.bar(xaxis + 1.0 * width, time_comp, width= width,  label="stencil - computation",color=col['comp'])
    plt.tight_layout(pad=pad)
    plt.savefig(overleaf+"/pie_stencil"+str(ranks[pie_rank])+".eps")



    # display time
    fig, ax = plt.subplots(num="pie ghost"+str(ranks[pie_rank]))
    # ratio_ghost = np.array([time_ghost_RMA_PS[pie_rank],time_ghost_RMA_getput[pie_rank],time_ghost_RMA_CW[pie_rank],time_ghost_RMA_comp[pie_rank]])/time_ghost[pie_rank]
    ratio_ghost = np.array([time_ghost_RMA_PS[pie_rank],time_ghost_RMA_getput[pie_rank],time_ghost_RMA_CW[pie_rank],time_ghost_RMA_comp[pie_rank]])
    print("--------------------------------------")
    print(np.sum(ratio_ghost))
    print(time_ghost[pie_rank])
    print("--------------------------------------")
    data = ratio_ghost
    def absolute_value(val):
        a  = data[ np.abs(data - val/100.*data.sum()).argmin() ]
        return "{:2.2f} sec".format(a)
    cols = col_pie['ghost']
    coles= col['ghost']
    patches = ax.pie(ratio_ghost,explode=(expl,expl,expl,expl),labels=("MPI PS","MPI Get/ MPI Put","MPI CW","computation"),normalize=True,
            autopct=absolute_value,
            colors=(cols,cols,cols,cols),
            wedgeprops={"edgecolor":coles,"linewidth":2,'antialiased': True},
            textprops={'fontsize': 14},
            pctdistance=0.65)[0]
    patches[3].set_hatch('/')
    plt.tight_layout(pad=pad)
    plt.savefig(overleaf+"/pie_ghost"+str(ranks[pie_rank])+".eps")


    # display time
    fig, ax = plt.subplots(num="pie adapt"+str(ranks[pie_rank]))
    # ratio_adapt = np.array([time_adapt_criterion[pie_rank],time_adapt_ghost[pie_rank],time_adapt_partition[pie_rank],time_adapt_adapt[pie_rank],time_adapt_misc[pie_rank],time_adapt_setupghost[pie_rank],time_adapt_destroyghost[pie_rank]])/time_adapt[pie_rank]
    ratio_adapt = np.array([time_adapt_criterion[pie_rank],time_adapt_smooth[pie_rank],time_adapt_ghost[pie_rank],time_adapt_partition[pie_rank], time_adapt_misc[pie_rank],time_adapt_setupghost[pie_rank]+time_adapt_destroyghost[pie_rank]])
    # print(np.sum(ratio_adapt))
    data = ratio_adapt
    def absolute_value(val):
        a  = data[ np.abs(data - val/100.*data.sum()).argmin() ]
        return "{:2.2f} sec".format(a)
    cols = col_pie['adapt']
    coles= col['adapt']
    patches = ax.pie(ratio_adapt,explode=(expl,expl,expl,expl,expl,expl),labels=("criterion","smooth +\n interp.","ghost","partition","sync","reset ghost"),normalize=True,
            autopct=absolute_value,
            colors=(cols,cols,cols,cols,cols,cols),
            wedgeprops={"edgecolor":coles,"linewidth":2,'antialiased': True},
            textprops={'fontsize': 14},
            pctdistance=0.65)[0]
    patches[0].set_hatch('/')
    patches[1].set_hatch('/')
    print("--------------------------------------")
    print(np.sum(ratio_adapt))
    print(time_adapt[pie_rank])
    print("--------------------------------------")
    plt.tight_layout(pad=pad)
    plt.savefig(overleaf+"/pie_adapt"+str(ranks[pie_rank])+".eps")



#==============================================================================================================================================================
# pie_rank = 0
# print("--------------------------------------")
# print(ranks[pie_rank])
# print("--------------------------------------")
# fig, ax = plt.subplots(num="pie stencil"+str(ranks[pie_rank]))
# # ratio_stencil = np.array([time_stencil_ghost[pie_rank],time_stencil_inner[pie_rank],time_stencil_outer[pie_rank]])/time_stencil[pie_rank]
# ratio_stencil = np.array([time_stencil_ghost[pie_rank],time_stencil_inner[pie_rank],time_stencil_outer[pie_rank]])
# print("--------------------------------------")
# print(np.sum(ratio_stencil))
# print(time_stencil[pie_rank])
# print("--------------------------------------")
# data = ratio_stencil
# def absolute_value(val):
#     a  = data[ np.abs(data - val/100.*data.sum()).argmin() ]
#     return "{:2.2f} sec".format(a)
# cols = col_pie['sten']
# coles= col['sten']
# patches = ax.pie(ratio_stencil,explode=(expl,expl,expl),labels=("ghost","inner","outer"),normalize=True,
#         autopct=absolute_value,
#         colors=(cols,cols,cols,cols),
#         wedgeprops={"edgecolor":coles,"linewidth":2,'antialiased': True},
#         textprops={'fontsize': 14},
#         pctdistance=0.65)[0]
# patches[1].set_hatch('/')
# patches[2].set_hatch('/')
# # ax.bar(xaxis - 1.0 * width, time_stencil,width= width,  label="ghost",color=col['sten'])
# # ax.bar(xaxis - 0.0 * width, time_ghost,width= width,   label="stencil - ghosting",color=col['ghost'])
# # ax.bar(xaxis + 1.0 * width, time_comp, width= width,  label="stencil - computation",color=col['comp'])
# plt.savefig(overleaf+"/pie_stencil"+str(ranks[pie_rank])+".eps")



# # display time
# fig, ax = plt.subplots(num="pie ghost"+str(ranks[pie_rank]))
# # ratio_ghost = np.array([time_ghost_RMA_PS[pie_rank],time_ghost_RMA_getput[pie_rank],time_ghost_RMA_CW[pie_rank],time_ghost_RMA_comp[pie_rank]])/time_ghost[pie_rank]
# ratio_ghost = np.array([time_ghost_RMA_PS[pie_rank],time_ghost_RMA_getput[pie_rank],time_ghost_RMA_CW[pie_rank],time_ghost_RMA_comp[pie_rank]])
# print("--------------------------------------")
# print(np.sum(ratio_ghost))
# print(time_ghost[pie_rank])
# print("--------------------------------------")
# data = ratio_ghost
# def absolute_value(val):
#     a  = data[ np.abs(data - val/100.*data.sum()).argmin() ]
#     return "{:2.2f} sec".format(a)
    
# cols = col_pie['ghost']
# coles= col['ghost']
# patches = ax.pie(ratio_ghost,explode=(expl,expl,expl,expl),labels=("MPI PS","MPI Get/ MPI Put","MPI CW","computation"),normalize=True,
#         autopct=absolute_value,
#         colors=(cols,cols,cols,cols),
#         wedgeprops={"edgecolor":coles,"linewidth":2,'antialiased': True},
#         textprops={'fontsize': 14},
#         pctdistance=0.65)[0]
# patches[3].set_hatch('/')
# # patches[2].set_hatch('/')

# plt.savefig(overleaf+"/pie_ghost"+str(ranks[pie_rank])+".eps")


# # display time
# fig, ax = plt.subplots(num="pie adapt"+str(ranks[pie_rank]))
# # ratio_adapt = np.array([time_adapt_criterion[pie_rank],time_adapt_ghost[pie_rank],time_adapt_partition[pie_rank],time_adapt_adapt[pie_rank],time_adapt_misc[pie_rank],time_adapt_setupghost[pie_rank],time_adapt_destroyghost[pie_rank]])/time_adapt[pie_rank]
# ratio_adapt = np.array([time_adapt_criterion[pie_rank],time_adapt_smooth[pie_rank],time_adapt_ghost[pie_rank],time_adapt_partition[pie_rank], time_adapt_misc[pie_rank],time_adapt_setupghost[pie_rank]+time_adapt_destroyghost[pie_rank]])
# # print(np.sum(ratio_adapt))
# cols = col_pie['adapt']
# coles= col['adapt']
# data = ratio_adapt
# def absolute_value(val):
#     a  = data[ np.abs(data - val/100.*data.sum()).argmin() ]
#     return "{:2.2f} sec".format(a)
# patches = ax.pie(ratio_adapt,explode=(expl,expl,expl,expl,expl,expl),labels=("criterion","smooth +\n interpolation","ghost","partition","sync","reset ghost"),normalize=True,
#         autopct=absolute_value,
#         colors=(cols,cols,cols,cols,cols,cols),
#         wedgeprops={"edgecolor":coles,"linewidth":2,'antialiased': True},
#         textprops={'fontsize': 14},
#         pctdistance=0.65)[0]
# patches[0].set_hatch('/')
# patches[1].set_hatch('/')
# print("--------------------------------------")
# print(np.sum(ratio_adapt))
# print(time_adapt[pie_rank])
# print("--------------------------------------")
# plt.savefig(overleaf+"/pie_adapt"+str(ranks[pie_rank])+".eps")


#==============================================================================================================================================================

# ax.bar(xaxis - 2.0 * width, time_ghost,width= width,  label="stencil",color=col['sten'])
# ax.bar(xaxis - 1.0 * width, time_ghost_RMA_PS,width= width,   label="stencil - ghosting",color=col['ghost'])
# ax.bar(xaxis + 0.0 * width, time_ghost_RMA_getput, width= width,  label="stencil - computation",color=col['comp'])
# ax.bar(xaxis + 1.0 * width, time_ghost_RMA_CW, width= width, label="adapt",color=col['adapt'])
# ax.bar(xaxis + 2.0 * width, time_ghost_RMA_comp, width= width, label="adapt - setup",color=col['aset'])
# # ax.bar(xaxis + 2.0 * width, time_adapt_criterion, width= width,  label="adapt - criterion",color=col['acrit'])
# # ax.bar(xaxis + 3.0 * width, time_adapt_ghost, width= width, label="adapt - ghost",color=col['aghost'])
xaxis = ranks

fig, ax = plt.subplots(num="barplot")
# ax.plot(xaxis, time_stencil, '.-', label="stencil",color=col['sten'])
# ax.plot(xaxis, time_adapt, '.-', label='adaptation',color=col['adapt'])
# ax.plot(xaxis, time_ghost, '.-', label='ghost',color=col['ghost'])
ax.plot(xaxis, time_stencil, '.-', label="stencil",color=col_bar['sten'])
ax.plot(xaxis, time_adapt, '.-', label='adaptation',color=col_bar['adapt'])
ax.plot(xaxis, time_ghost, '.-', label='ghost',color=col_bar['ghost'])

plt.xticks(xaxis)
plt.xlim([min(xaxis), max(xaxis)])
# plt.xlim([min(xaxis), max(xaxis)])
plt.xticks(xaxis)
# plt.title("weak scaling - NIC 5")
# plt.legend(bbox_to_anchor=(0.5,1.0),loc='lower center',ncol = int(nbars/2+1) )
plt.grid()
plt.legend()
plt.ylabel("[sec]")
plt.xlabel("ranks")
# plt.yscale("log")
tikzplotlib.save(overleaf+"/time.tex",axis_width='\\textwidth',axis_height='0.4\\textwidth',strict=True)

id_ref = 1
# eff_dodt = fact / np.multiply(time_dodt[id_ref:],ranks[id_ref:])
eff_dodt = time_dodt[id_ref]/time_dodt
eff_stencil = time_stencil[id_ref]/time_stencil
eff_comp = time_comp[id_ref]/time_comp
eff_ghost = time_ghost[id_ref]/time_ghost
eff_adapt = time_adapt[id_ref]/time_adapt
# eff_adapt_setup = time_adapt_setup[id_ref]/time_adapt_setup
# eff_adapt_smooth = time_adapt_smooth[id_ref]/time_adapt_smooth
# eff_adapt_criterion = time_adapt_criterion[id_ref]/time_adapt_criterion
# eff_adapt_ghost = time_adapt_ghost[id_ref]/time_adapt_ghost
# eff_p4est = time_p4est[id_ref]/time_p4est
# eff_wait = time_rma_wait[id_ref]/time_rma_wait
# eff_ginit = time_adapt_ginit[id_ref]/time_adapt_ginit


fig, ax = plt.subplots(num="eff")
# ax.bar(xaxis - 3 * width, eff_stencil, label="stencil")
# ax.bar(xaxis - 2 * width, eff_ghost,  label="stencil - ghosting")
# ax.bar(xaxis - 1 * width, eff_comp,  label="stencil - computation")
# ax.bar(xaxis + 0 * width, eff_adapt, label="adapt")
# ax.bar(xaxis + 1 * width, eff_adapt_setup, label="adapt - setup")
# ax.bar(xaxis + 2 * width, eff_adapt_criterion,  label="adapt - criterion")
# ax.bar(xaxis + 3 * width, eff_adapt_ghost, label="adapt - ghost")

# ax.plot(xaxis[id_ref:], eff_stencil[id_ref:], '.-', label="stencil",color=col['sten'])
# ax.plot(xaxis[id_ref:], eff_ghost[id_ref:], '.-', label="ghost",color=col['ghost'])
# ax.plot(xaxis[id_ref:], eff_adapt[id_ref:], '.-', label="adaptation",color=col['adapt'])
ax.plot(xaxis[id_ref:], eff_stencil[id_ref:], '.-', label="stencil",color=col_bar['sten'])
ax.plot(xaxis[id_ref:], eff_ghost[id_ref:], '.-', label="ghost",color=col_bar['ghost'])
ax.plot(xaxis[id_ref:], eff_adapt[id_ref:], '.-', label="adaptation",color=col_bar['adapt'])

plt.xticks(xaxis[id_ref:])
# plt.ylim([min(eff_dodt),min(eff_dodt)])
plt.ylim([0.75,1.01])
plt.xlim([min(xaxis[id_ref:]), max(xaxis[id_ref:])])
plt.legend()
plt.grid()
plt.xlabel("ranks")
plt.ylabel("$\\eta_{w}$")

tikzplotlib.save(overleaf+"/weak.tex",axis_width='\\textwidth',axis_height='0.4\\textwidth',strict=True)


# ax = plt.subplot(2, 2, 4)
# ax.plot(xaxis[id_ref:], eff_adapt[id_ref:], '.--', label="adapt",color=col['adapt'])
# ax.plot(xaxis[id_ref:], eff_adapt_setup[id_ref:], '.-', label="adapt - setup",color=col['aset'])
# ax.plot(xaxis[id_ref:], eff_adapt_criterion[id_ref:], '.-', label="adapt - criterion",color=col['acrit'])
# ax.plot(xaxis[id_ref:], eff_adapt_ghost[id_ref:], '.-', label="adapt - ghost",color=col['aghost'])

# # ax.plot(nodes, eff_wait, '.-', label="RMA wait")
# # ax.plot(nodes, eff_p4est, '.-', label="adapt p4est")
# # ax.plot(nodes, eff_ginit, '.-', label="adapt ghost init")

# # ax.plot(ranks[id_ref:],eff_ghost,'.-',label="ghost")
# # ax.plot(ranks[id_ref:],eff_rma_wait,'.-',label="RMA wait")
# plt.xticks(xaxis[id_ref:])
# # plt.ylim([min(eff_dodt),min(eff_dodt)])
# plt.ylim([0.75,1.01])
# plt.xlim([min(xaxis[id_ref:]), max(xaxis[id_ref:])])
# plt.legend()
# plt.grid()
# plt.xlabel("ranks")
# plt.ylabel("$\eta$")

plt.show()
