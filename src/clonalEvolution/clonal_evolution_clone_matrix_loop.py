import numpy as np
import pandas as pd
import scipy as sc
import math
import time
import copy

MUTATION_ID = 1
CLONE_ID = 1
tx = [0,0,0,0]
ty = [0,0,0,0]

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
def calculateFitness(cells, passenger, driver, m_effect):
    p = sum(sum(passenger.toarray()))
    d = len(driver)*cells
    return (1 + m_effect[0])**d/(1 - m_effect[1])**p

def deleteZeroColumn(clone):
    mask = np.ones(clone[4]._shape[1], dtype=bool)
    for i in range(clone[4]._shape[1]):
        _sum = (clone[4][:,i]).count_nonzero()
        if _sum == 0:
            mask[i] = False
    clone[4] = clone[4][:,mask]
    clone[3] = np.array(clone[3])[mask].tolist()

def dying(clone, death):
    if death > 0:
        d = np.random.choice(range(clone[1]), size=death, replace=False)
        mask = np.ones(clone[1], dtype=bool)
        mask[d] = False
        clone[4] = clone[4][mask,:]
        clone[1] = clone[1] - death

def division(clone, divide):
    if divide > 0:
        d = np.random.choice(range(clone[1]), size=divide, replace=False)
        mask = np.zeros(clone[1], dtype=bool)
        mask[d] = True
        clone[4] = sc.sparse.vstack([clone[4], copy.deepcopy(clone[4][mask,:])]).tocsr()
        clone[1] = clone[1] + divide

def newClone(clone, driver, iPop, mut_effect):
    global CLONE_ID, MUTATION_ID
    if driver > 0:
        for x in range(CLONE_ID, CLONE_ID + driver, 1):
            d = np.random.choice(range(clone[1]))
            driv = copy.deepcopy(clone[2])
            driv.append(MUTATION_ID)
            iPop.append([x, 1, 
                        driv, 
                        copy.deepcopy(clone[3]), 
                        clone[4].getrow(d), 
                        clone[5]*(1+mut_effect[0]),
                        clone[0]])
            MUTATION_ID = MUTATION_ID + 1
        CLONE_ID = CLONE_ID + driver

def newMutation(clone, passenger, mut_effect):
    global MUTATION_ID
    if passenger > 0:
        f_n = clone[5]/(1-mut_effect[1])*passenger
        f_o = clone[5]*clone[1]
        clone[5] = (f_n + f_o)/(clone[1]+passenger)
        for x in range(MUTATION_ID, MUTATION_ID + passenger, 1):
            d = np.random.choice(range(clone[1]))
            clone[3].append(x)
            clone[4]._shape = (clone[4]._shape[0], clone[4]._shape[1] + 1)
            clone[4] = sc.sparse.vstack([clone[4], clone[4].getrow(d)]).tocsr()
            clone[4][clone[4]._shape[0]-1,clone[4]._shape[1]-1] = 1
            clone[1] = clone[1] + 1
        MUTATION_ID = MUTATION_ID + passenger

def clonalEvolutionCloneMatrixLoop(iPop, cap, tau, mut_prob, mut_effect, resume, q, threads, print_time):
    global tx, ty, CLONE_ID, MUTATION_ID
    """
    Assumption:
        Mutation matrix in compressed form. 1 means mutation (column) occurs in cell (row)
        
    Description:
        One cycle to update population - tau loop binned method
        Prameters:
            iPop: population matrix where row is in form of:
                Clone number
                Cell number
                Driver mutation list
                Uniqal passenger mutation list
                Mutation matrix
                Clone fitness
                Previous clone number
            cap: population capacity
            tau: tau step
            mut_prob: list in form of: [driver mutation probability, passenger mutation probability]
            mut_effect: list in form of: [driver mutation effect, passenger muatation effect]
            resume: acknowledge to resume simulation
            q: common queue
            THREADS: threads number used in simulation         
    """   
    popSize = sum([row[1] for row in iPop])
    mdt = popSize/cap
    for i in iPop:
        pdv = 1 - math.exp(-tau*i[5])
        pdt = 1 - math.exp(-tau*mdt)
        
        pdv = (1-pdt)*pdv
        pdm_d = (mut_prob[0])*(1 - mut_prob[1])
        pdm_p = (1 - mut_prob[0])*(mut_prob[1])
        pdr = (1 - pdt)*(1 - pdv)*(1 - pdm_d)*(1 - pdm_p)
        
        r = np.random.multinomial(i[1], [pdt, pdv, pdr])
        
        death = r[0]
        divide = r[1]
        
        m = np.random.multinomial(divide, [pdm_p, pdm_d, (1-pdm_p)*(1-pdm_d)])

        m_p = m[0]        
        m_d = m[1]        

        time_t = time.time()             
        ## dying cells
        dying(i, death)
          
        tx[0] = tx[0] + (time.time() - time_t) 
        ty[0] = ty[0] + 1
        time_t = time.time()     
        ## new clones
        newClone(i, m_d, iPop, mut_effect)
        
        tx[1] = tx[1] + (time.time() - time_t) 
        ty[1] = ty[1] + 1
        time_t = time.time()     
        ## division
        division(i, divide - m_p - m_d)
        
        tx[2] = tx[2] + (time.time() - time_t) 
        ty[2] = ty[2] + 1
        time_t = time.time()     
        ## mean fitness
        newMutation(i, m_p, mut_effect)
        
        tx[3] = tx[3] + (time.time() - time_t) 
        ty[3] = ty[3] + 1
        time_t = time.time()     
        
        try:
            if print_time:
                print("Clone ID: %i, Dying: %.2f, newClone: %.2f, Division: %.2f, newMutation: %.2f" % (CLONE_ID,tx[0]/ty[0],tx[1]/ty[1],tx[2]/ty[2],tx[3]/ty[3]))
                tx = [0,0,0,0]
                ty = [0,0,0,0]
                ##delete mutation with 0 occurences
                deleteZeroColumn(i)
        except ZeroDivisionError:
            continue
        
    for i in iPop:
        if i[1] == 0:
            iPop.remove(i)
        
    return iPop
    