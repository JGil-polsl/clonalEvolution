'''
    Cellular/Microbial Clonal Evolution simulations basing on Gillespie algorithm.
    Copyright (C) 2022 by Jarosław Gil

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''

import matplotlib.pyplot as plt
import numpy as np
import scipy as sc
import time
import os
import copy
import pandas as pd
from pathlib import Path  
from threading import Thread

import clonalEvolution.clonal_evolution_binned_loop as CEBL
import clonalEvolution.clonal_evolution_clone_matrix_loop as CECML
import clonalEvolution.clonal_evolution_loop as CEL
import clonalEvolution.wmean as wm

end = False
    
def waveMut_ar(iMuts, iClones, cln):
    muts = iMuts
    clones = iClones
    m_u = []
    c_u = []
    
    mim = int(min(muts))
    mxm = int(max(muts))
    m_u = [x for x in range(mim, mxm+1)]
    if cln:
        for i in clones:
            if i not in c_u:
                c_u.append(i)
    else:
        c_u = [0]
    TAB = np.zeros((mxm-mim+1,len(c_u)+1))
    c_u.sort()
    t = [(muts[i], clones[i]) for i in range(len(muts))]
    for m,c in t:
        if cln:
            ci = c_u.index(c)
        else:
            ci = 0
        xx = TAB[m-mim][ci+1] 
        TAB[m-mim][ci+1] = xx + 1
        TAB[m-mim][0] = m
    
    return muts, clones, m_u, c_u, TAB

def saveMatrixData(df, file_localization, file_name, iter_outer):
    names = df['Mutation matrix']
    for i in names:
        sc.sparse.save_npz(i[0], i[1])
    df['Mutation matrix'] = [x[0] for x in names]
    filepath = Path(file_localization + '/' + file_name + '_' + str(iter_outer).replace('.','_') + ".mtx")      
    filepath.parent.mkdir(parents=True, exist_ok=True)  
    df.to_csv(filepath)     

def saveToFile(df, file_localization, file_name, iter_outer, ending):    
    filepath = Path(file_localization + '/' + file_name + '_' + str(iter_outer).replace('.','_') + ending)      
    filepath.parent.mkdir(parents=True, exist_ok=True)  
    df.to_csv(filepath)  

def mutationWavePlot(iMuts, iClones=None, cln=0, name="", local="", f_num=0, ret=0):
    mutWave = []
    # iMuts = [row[0] for row in iPop]
    # iClones = [row[2] for row in iPop]
    
    MUT, CLONE, UNIMT, UNICL, TAB = waveMut_ar(iMuts, iClones, cln)
    sr = wm(MUT, np.ones(len(MUT)))
    mutWave = pd.DataFrame(TAB) 
    
    if ret:
        return mutWave
    else:       
        mutWave.columns = [str(x) for x in mutWave.columns]
        a = ["mut_num"]
        b = [str(x) + " " + str(mutWave[str(UNICL.index(x)+1)].sum()) for x in UNICL]
        a.extend(b)
        mutWave.columns = a
        ax = mutWave.plot.bar(x="mut_num", stacked=True, title=("Mutation wave, wmean %.3f" % (sr)))
        lgd = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize="small")   
        
        fig = ax.get_figure()
        try:
            os.makedirs(local + "/Figures/" + str(f_num), exist_ok=True) 
        except OSError as error:
            print(error)
        finally:
            if os.path.exists(local + "/Figures/" + str(f_num) + '/' + name + "_mutWave.jpg"):
                os.remove(local + "/Figures/" + str(f_num) + '/' + name + "_mutWave.jpg")
            fig.savefig(local + "/Figures/" + str(f_num) + '/' + name + "_mutWave.jpg", bbox_extra_artists=(lgd,), bbox_inches='tight')
            plt.close(fig)
    
def VAFPlot(iMutations, iClones=None, cln=0, name="", local="", f_num=0):
    VAF = []
    clones = []
    num = []
    temp = []
    
    if cln:
        for x in iClones:
            if x not in clones:
                clones.append(x)
                num.append(1)
                VAF.append([])
                temp.append([])
            else:
                num[clones.index(x)] = num[clones.index(x)] + 1
                
        for x in range(len(iClones)):
            idx = clones.index(iClones[x])
            for i in iMutations[x]:
                if i not in temp[idx]:
                    temp[idx].append(i)
                    VAF[idx].append(1)
                else:
                    VAF[idx][temp[idx].index(i)] = VAF[idx][temp[idx].index(i)] + 1
        
        for x in range(len(VAF)):
            VAF[x] = [f/num[x] for f in VAF[x]]
        VAF = pd.DataFrame(VAF).T
            
    else:
        pop = len(iMutations)
        for x in iMutations:
            for i in x:
                if i not in temp:
                    temp.append(i)
                    VAF.append(1)
                else:
                    VAF[temp.index(i)] = VAF[temp.index(i)] + 1
        VAF = [x/pop for x in VAF]
        VAF = pd.DataFrame(VAF).T
        
    ax = VAF.plot.hist(stacked=True)
    fig = ax.get_figure()
    try:
        os.makedirs(local + "/Figures/" + str(f_num), exist_ok=True) 
    except OSError as error:
        print(error)
    finally:
        if os.path.exists(local + "/Figures/" + str(f_num) + '/' + name + "_VAF.jpg"):
            os.remove(local + "/Figures/" + str(f_num) + '/' + name + "_VAF.jpg")
        fig.savefig(local + "/Figures/" + str(f_num) + '/' + name + "_VAF.jpg")
        plt.close(fig)

def fitnessWavePlot(iProp, iClones=None, cln=0, name="", local="", f_num=0):
    fitnessWave = []    
    if cln:
        fitnessWave = pd.DataFrame(fitnessWave)
    else:
        fitnessWave = pd.DataFrame(iProp)
        
    ax = fitnessWave.plot.hist(stacked=True)
    fig = ax.get_figure()
    try:
        os.makedirs(local + "/Figures/" + str(f_num), exist_ok=True) 
    except OSError as error:
        print(error)
    finally:
        if os.path.exists(local + "/Figures/" + str(f_num) + '/' + name + "_fitness.jpg"):
            os.remove(local + "/Figures/" + str(f_num) + '/' + name + "_fitness.jpg")
        fig.savefig(local + "/Figures/" + str(f_num) + '/' + name + "_fitness.jpg")
        plt.close(fig)

def mullerPlot():
    muller = []
    return muller

def command(queue_data, q, select, iPop, ID, iter_outer, iter_inner, tau, skip, cycle, file_localization, file_name):
    global end
    if(queue_data[0] == '1' and queue_data[1] == str(ID)):
        if(queue_data[2] == "exit"):
            print("exit")
            end = True
        elif(queue_data[2] == "size"):
            if select == 0:
                q.put(['0', str(ID), str(len(iPop))])
            elif select == 1:
                q.put(['0', str(ID), str(sum([row[1] for row in iPop]))])
        elif(queue_data[2] == "time"):
            q.put(['0', str(ID), str((iter_outer-1)*skip + (iter_inner%cycle)*tau)])
        elif(queue_data[2] == "save"):
            if file_localization == "" or file_name == "":
                q.put(['0', str(ID), "no path or file, save not active"])
            else:
                if select == 0:
                    tsf = Thread(target=saveToFile, args=(pd.DataFrame({
                                                                'Clone': [row[2] for row in iPop],
                                                                'Mutations': [row[0] for row in iPop],
                                                                'Fitness': [row[1] for row in iPop],
                                                                'Mutations_ID': [row[3] for row in iPop],
                                                                'Mutations_effect': [row[4] for row in iPop],
                                                            }),  
                                                          file_localization, file_name + '_single', iter_outer - 1 + (iter_inner%cycle)/cycle, ".csv"))
                    tsf.start()
                elif select == 1:
                    tsf = Thread(target=saveToFile, args=(pd.DataFrame({
                                                                'Clone number': [row[0] for row in iPop],
                                                                'Cells number': [row[1] for row in iPop],
                                                                'Mean fitness': [row[2] for row in iPop],
                                                                'Mean mutation number': [row[3] for row in iPop],
                                                                'Driver mutation list': [row[4] for row in iPop],
                                                                'Passener mutation list': [row[5] for row in iPop],
                                                                'Previous clone number': [row[6] for row in iPop]
                                                            }),  
                                                          file_localization, file_name + '_binned', iter_outer - 1 + (iter_inner%cycle)/cycle, ".txt"))
                    tsf.start()
        elif(queue_data[2] == "plot"):        
            if select == 0:
                MUT, CLONE, UNIMT, UNICL, TAB = waveMut_ar([row[0] for row in iPop], [row[2] for row in iPop], 1)
                sr = wm(MUT, np.ones(len(MUT)))
                mutWave = pd.DataFrame(TAB) 
                
                mutWave.columns = [str(x) for x in mutWave.columns]
                a = ["mut_num"]
                b = [str(x) + " " + str(mutWave[str(UNICL.index(x)+1)].sum()) for x in UNICL]
                a.extend(b)
                mutWave.columns = a
                
                q.put(["-1", str(ID), mutWave]) 
            elif select == 1:
                df = pd.DataFrame({
                        'Clone': [row[0] for row in iPop],
                        'Cells number': [row[1] for row in iPop],
                    })
                q.put(["-2", str(ID), df])
            elif select == 2:
                df = pd.DataFrame({
                        'Clone': [row[0] for row in iPop],
                        'Cells number': [row[1] for row in iPop],
                    })
                q.put(["-2", str(ID), df])
        elif(queue_data[2] == 'clone'):
            if select == 0:
                _CLONES = []
                for i in iPop:
                    if i[2] not in _CLONES:
                        _CLONES.append(i[2])
                q.put(['0', str(ID), str(len(_CLONES))])
            elif select == 1:
                q.put(['0', str(ID), str(len(iPop))])
    else:
        q.put(queue_data)    

def plotter(plots, iPop, file_name, file_localization, iter_outer, select):
    if select == 0:
        if plots & 1:
            tsf = Thread(target=mutationWavePlot, args=([row[0] for row in iPop], 
                                                        [row[2] for row in iPop], 
                                                        1, file_name, file_localization, iter_outer, 0))
            tsf.start()
        if plots & 2:
            tsf = Thread(target=fitnessWavePlot, args=([row[1] for row in iPop], 
                                                       [row[2] for row in iPop],
                                                        0, file_name, file_localization, iter_outer))
            tsf.start()            
        if plots & 4:
            tsf = Thread(target=VAFPlot, args=([row[3] for row in iPop], 
                                               [row[2] for row in iPop], 
                                               1, file_name, file_localization, iter_outer))
            tsf.start()        
        if plots & 8:
            _MP = [] ## TODO muller plot (for now out from simulation)        
        if plots & 16: 
            tsf = Thread(target=saveToFile, args=(pd.DataFrame({
                                                        'Clone': [row[2] for row in iPop],
                                                        'Mutations': [row[0] for row in iPop],
                                                        'Fitness': [row[1] for row in iPop],
                                                        'Mutations_ID': [row[3] for row in iPop],
                                                        'Mutations_effect': [row[4] for row in iPop],
                                                        'Parent clone': [row[5] for row in iPop]
                                                    }),  
                                                  file_localization, file_name + "_single", iter_outer, '.csv'))
            tsf.start()
            
    elif select == 1:
        if plots & 1:
            _MW = [] ## TODO mutation wave plot save for binned version
        if plots & 2:
            _FW = [] ## TODO fitness wave plot save for binned version
        if plots & 4:
            _VP = [] ## TODO VAF plot save for binned version
        if plots & 8:
            _MP = [] ## TODO muller plot (for now out from simulation)
        if plots & 16:               
            tsf = Thread(target=saveToFile, args=(pd.DataFrame({
                                                        'Clone number': [row[0] for row in iPop],
                                                        'Cells number': [row[1] for row in iPop],
                                                        'Mean fitness': [row[2] for row in iPop],
                                                        'Mean mutation number': [row[3] for row in iPop],
                                                        'Driver mutation list': [row[4] for row in iPop],
                                                        'Passener mutation list': [row[5] for row in iPop],
                                                        'Previous clone number': [row[6] for row in iPop]
                                                    }),  
                                                  file_localization, file_name + '_binned', iter_outer, ".txt"))
            tsf.start()  
    
    elif select == 2:
        if plots & 16:
            names = []
            for clone in iPop:
                s_path = file_localization + '/Data/' + str(iter_outer).replace('.','_') + '/'
                try:
                    os.makedirs(s_path, exist_ok=True) 
                except OSError as error:
                    print(error)
                finally:
                    s_path = s_path + file_name + '_clone_data_' + str(clone[0]) + '.npz'
                names.append([s_path, clone[4]])
            
            tsf = Thread(target=saveMatrixData, args=(pd.DataFrame({
                                                            'Clone number': [row[0] for row in iPop],
                                                            'Cells number': [row[1] for row in iPop],
                                                            'Driver mutation list': [row[2] for row in iPop],
                                                            'Uniqal passenger mutation list': [row[3] for row in iPop],
                                                            'Mutation matrix': [row for row in names],
                                                            'Clone fitness': [row[5] for row in iPop],
                                                            'Previous clone number': [row[6] for row in iPop]
                                                        }), file_localization, file_name, iter_outer))
            tsf.start()

def clonalEvolutionMainLoop(iPop, params, file_name="", file_description="", file_localization="", plots=0, t_iter=0, q=None, ID=0, select=0):
    """
    Main simulation loop
        iPop: population array where row is in form of: 
            SINGLE CELL ALGORITH
                Mutation number  
                Fitness
                Clone number
                Mutations list
                Mutation effect list
                Parent clone
            BINNED ALGORITH
                Clone number
                Cells number
                Mean fitness
                Mean mutation number
                Driver mutation list
                Passener mutation list
                Previous clone number
            MATRIX ALGORITHM
                Clone number
                Driver mutation list
                Uniqal passener mutation list
                Mutation matrix
                Previous clone number
        params: simulation parameters with structure 
        [initial population size, population capacity, simulation steps number, tau step [s], one cycle time [s], 
         mutation probability (list) adequate index for mutation effect, mutation effect (list), number of threads]
        file_name: save file name, data will be saved to file_name.csv, figures: file_name_typeofplot.jpg
        file_localization: path to file, figures will be saved in path: file_localization/Figures/Cycle_number
        plots: binary interpreted value defines which plots will be generated: 
            1 - mutation wave (TODO for binned version)
            2 - fitness wave (TODO for binned version)
            4 - mutation histogram (VAF) (TODO for binned version)
            8 - muller plot (TODO)
            16 - ack to save data
        t_iter - value of starting iteration (for resume t_iter != 0)
        q - queue to comunicate with simulation in thread
        ID - simulation ID
        select - 0 normal, 1 binned
    """
    global end
    cap = params[1]
    steps = params[2]
    tau = params[3]
    skip = params[4]
    mut_effect = params[6]
    mut_prob = params[5]   
    threads = params[7]

    iter_inner = 0 + 1*(t_iter>0)
    begin = 0 + 1*(t_iter==0)
    iter_outer = t_iter
    resume = 0 + 1*(t_iter>0)
    cycle = round(skip/tau)
    
    t = time.time()
    
    if select == 0 and begin and plots & 16:
        os.makedirs(file_localization, exist_ok=True) 
        FILE = open(file_localization + '/' + file_name + "_muller_plot_single_data.txt", 'w')
        FILE.write("clone, cells, previous clone")
        FILE.write('\n')
        FILE.close()
    elif select == 1 and begin and plots & 16:        
        os.makedirs(file_localization, exist_ok=True) 
        FILE = open(file_localization + '/' + file_name + "_muller_plot_binned_data.txt", 'w')
        FILE.write("clone, cells, previous clone")
        FILE.write('\n')
        FILE.close()
    elif select == 2 and begin and plots & 16:        
        os.makedirs(file_localization, exist_ok=True) 
        FILE = open(file_localization + '/' + file_name + "_muller_plot_matrix_data.txt", 'w')
        FILE.write("clone, cells, previous clone")
        FILE.write('\n')
        FILE.close()


    while 1:
        
        ## Data for muller plot
        if iter_inner % (cycle/skip) == 0 and plots & 16:
            if select == 0:
                FILE = open(file_localization + '/' + file_name + "_muller_plot_single_data.txt", 'a')
                _CLONES = []
                _CELLS = []
                _LAST = []
                for i in iPop:
                    if i[2] not in _CLONES:
                        _CLONES.append(i[2])
                        _CELLS.append(1)
                        _LAST.append(i[5])
                    else:
                        _CELLS[_CLONES.index(i[2])] = _CELLS[_CLONES.index(i[2])] + 1
                FILE.write("".join(["(%.0f, %.0f, %.0f)" % (_CLONES[i], _CELLS[i], _LAST[i]) for i in range(len(_CLONES))]))
                FILE.write('\n')
                FILE.close()
            elif select == 1:
                FILE = open(file_localization + '/' + file_name + "_muller_plot_binned_data.txt", 'a')
                FILE.write("".join(["(%.0f, %.0f, %.0f)" % (row[0], row[1], row[6]) for row in iPop]))
                FILE.write('\n')
                FILE.close()
            elif select == 2:
                FILE = open(file_localization + '/' + file_name + "_muller_plot_matrix_data.txt", 'a')
                FILE.write("".join(["(%.0f, %.0f, %.0f)" % (row[0], row[1], row[6]) for row in iPop]))
                FILE.write('\n')
                FILE.close()
        
        ## Commands from GUI
        if q != None:
            if not q.empty():
                queue_data = q.get()
                command(queue_data, q, select, iPop, ID, iter_outer, iter_inner, tau, skip, cycle, file_localization, file_name)

        if end:
            break
        
        ## Save data, plots
        if iter_inner % cycle == 0 or begin:
            begin = 0
            t = time.time() - t  
            print(str(ID) + ':' + str(t))
            
            plotter(plots, iPop, file_name, file_localization, iter_outer, select)                                 
            
            iter_outer = iter_outer + 1            
            t = time.time()

        print_time = not iter_inner % round(cycle/skip)        
        # if iter_outer % steps == 0:
        #     print("pointer")
        
        if select == 0:            
            iPop = CEL.clonalEvolutionLoop(iPop, cap, tau, mut_prob, mut_effect, resume, q, threads)
        elif select == 1:
            iPop = CEBL.clonalEvolutionBinnedLoop(iPop, cap, tau, mut_prob, mut_effect, resume, q, threads, print_time)
        elif select == 2:
            iPop = CECML.clonalEvolutionCloneMatrixLoop(iPop, cap, tau, mut_prob, mut_effect, resume, q, threads, print_time)
            
        resume = 0
        iter_inner = iter_inner + 1
