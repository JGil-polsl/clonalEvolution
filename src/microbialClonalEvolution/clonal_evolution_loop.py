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

import numpy as np
import math
from threading import Thread
import copy

MUTATION_ID = 1
CLONE_ID = 1
tCLONE = 0
semafor_1 = False
semafor_2 = False
semafor_3 = False

def division(i, c, m_d, m_p, iPop, mut_effect, mut_prob):
    global MUTATION_ID, CLONE_ID, semafor_1, semafor_2, semafor_3
    if m_d[i] or m_p[i]:
        idx = 2
        if m_d[i]:
            idx = 0
        elif m_p[i]:
            idx = 1
        else:
            idx = 2
        eff = mut_effect[idx]    
        
        t_m = c[3]          
        t_e = c[4]
        
        while semafor_1:
            continue
        semafor_1 = True
        t_m.append(MUTATION_ID)
        MUTATION_ID = MUTATION_ID + 1
        semafor_1 = False
        t_e.append(eff)
        
        mut = c[0] + 1
        prop = c[1]
        if eff > 0:
            prop = prop*(1+eff)
        else:
            prop = prop/(1-eff)
        clone = 0
        parent = 0
        if eff > 0:          
            while semafor_2:
                continue
            semafor_2 = True
            clone = CLONE_ID
            CLONE_ID = CLONE_ID + 1
            semafor_2 = False
            parent = c[2]
        else:
            parent = c[5]
            clone = c[2]
        try:
            iPop.append([mut, prop, clone, t_m, t_e, parent])
        except:
            print("cell dead")
    else:
        try:
            iPop.append(c)
        except:
            print("cell dead")

def death(val, iPop):
    iPop.remove(val)

def th_div(pdv, pdv_idx, c, m_d, m_p, iPop, mut_effect, mut_prob):
    t = pdv[c[0]:c[1]]
    for i in range(len(t)):
        division(pdv_idx[i], t[i], m_d, m_p, iPop, mut_effect, mut_prob)
    
def th_dth(pdt, c, iPop):
    t = pdt[c[0]:c[1]]
    for i in reversed(t):
        death(i, iPop)
        
def prob_dvd(pdv, pdv_idx, idx_dvd, popSize, iPop, tau, THREADS):
    pdvx = np.random.exponential(1,popSize)
    pdvx = np.array([pdvx[i]/iPop[i][1] for i in range(popSize)])
    pdvx = np.where(pdvx < tau)[0] 
    pdv_idx.extend(pdvx)
    pdv.extend([copy.deepcopy(iPop[pdvx[x]]) for x in range(len(pdvx))])
    
    th = math.ceil(len(pdvx)/THREADS)
    idx_dvdx = [(int(x*th),int((x+1)*th)) for x in range(THREADS)]
    idx_dvdx[len(idx_dvdx)-1] = (idx_dvdx[len(idx_dvdx) - 1][0], len(pdvx))
    idx_dvd.extend(idx_dvdx)
    
def prob_dth(pdt, idx_dth, popSize, iPop, tau, mdt, THREADS):
    pdtx = np.random.exponential(1, popSize)
    pdtx = np.array([pdtx[i]/(mdt) for i in range(popSize)])
    pdtx = np.where(pdtx < tau)[0]    
    pdt.extend([iPop[pdtx[x]] for x in range(len(pdtx))])

    th = math.ceil(len(pdtx)/THREADS)
    idx_dthx = [(int(x*th),int((x+1)*th)) for x in range(THREADS)]
    idx_dthx[len(idx_dthx)-1] = (idx_dthx[len(idx_dthx) - 1][0], len(pdtx))  
    idx_dth.extend(idx_dthx)
    
def clonalEvolutionLoop(iPop, cap, tau, mut_prob, mut_effect, resume, q, THREADS):
    """    
    Description:
        One cycle to update population - tau loop iteration method
        Prameters:
            iPop: population matrix where row is in form of:
                Mutation number  
                Fitness
                Clone number
                Mutations list
                Mutation effect list
                Parent clone
            cap: population capacity
            tau: tau step
            mut_prob: list in form of: [driver mutation probability, passenger mutation probability]
            mut_effect: list in form of: [driver mutation effect, passenger muatation effect]
            resume: acknowledge to resume simulation
            q: common queue
            THREADS: threads number used in simulation            
    """
    global MUTATION_ID, CLONE_ID, tCLONE
    if resume:
        CLONE_ID = max([row[2] for row in iPop]) + 1
        MUTATION_ID = 0
        i = [max(row[3]) for row in iPop]
        MUTATION_ID = max(MUTATION_ID, max(i))
        MUTATION_ID = MUTATION_ID + 1    

    popSize = len(iPop)
    mdt = popSize/cap
    # mdt = math.log(1 + (math.e - 1)*popSize/cap)

    m_d = np.random.binomial(1, mut_prob[0], popSize)
    m_p = np.random.binomial(1, mut_prob[1], popSize)
    
    pdt = []
    pdv = []
    pdv_idx = []
    
    idx_dth = []
    idx_dvd = []
    
    dth = Thread(target=prob_dth, args=(pdt, idx_dth, popSize, iPop, tau, mdt, THREADS))
    dth.start()
    dvd = Thread(target=prob_dvd, args=(pdv, pdv_idx, idx_dvd, popSize, iPop, tau, THREADS))
    dvd.start()
    
    dth.join()
    dvd.join()
    
    divides = []

    for i in range(len(idx_dvd)):
        divides.append(Thread(target=th_div, args=(pdv, pdv_idx, idx_dvd[i], m_d, m_p, iPop, mut_effect, mut_prob)))
        divides[i].start()
    
    deaths = []
  
    for i in range(len(idx_dth)):
        deaths.append(Thread(target=th_dth, args=(pdt, idx_dth[i], iPop)))
        deaths[i].start()
        
    for i in divides:
        i.join()    
    for i in deaths:
        i.join()
    
    if tCLONE != CLONE_ID:
        tCLONE = CLONE_ID        
    return iPop