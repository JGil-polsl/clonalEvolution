import numpy as np
import pandas as pd
import scipy as sc
import math
import time
import copy
import threading
from threading import Thread

MUTATION_ID = 1
CLONE_ID = 1
tx = [0,0,0,0]
ty = [0,0,0,0]
sem_mutations = True
sem_clone = True

sgn = lambda x: -1 * (x < 0) + 1 * (x > 0) 

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
    clone[4] = clone[4][death,:]
    clone[1] = clone[1] - len(death[death == False])
    clone[5] = clone[5][death]

def division(clone, divide):
    clone[4] = sc.sparse.vstack([clone[4], copy.deepcopy(clone[4][divide,:])]).tocsr()
    clone[1] = clone[1] + len(divide[divide])
    clone[5] = np.append(clone[5], clone[5][divide])

def newClone(clone, driver, iPop, mut_effect):
    global CLONE_ID, MUTATION_ID, sem_clone, sem_mutations
    new_d = clone[4][driver,:]
    fit = clone[5][driver]
    mutations = []
    clones = []
    
    while True:
        if sem_mutations:
            sem_mutations = False
            mutations = [x for x in range(MUTATION_ID, MUTATION_ID + len(fit), 1)]
            MUTATION_ID = MUTATION_ID + len(fit)
            sem_mutations = True
            break
    while True:
        if sem_clone:
            sem_clone = False
            clones = [x for x in range(CLONE_ID, CLONE_ID + len(fit), 1)]
            CLONE_ID = CLONE_ID + len(fit)
            sem_clone = True
            break
        
    if len(fit) == 0:
        return
    idx = 0
    for i in new_d:
        driv = copy.deepcopy(clone[2])
        driv.append(mutations[idx])
        iPop.append([
            clones[idx], 1,
            driv,
            copy.deepcopy(clone[3]),
            copy.deepcopy(i),
            np.array([fit[idx]*(1+mut_effect[0])]),
            clone[0]])
        idx = idx + 1

def newMutation(clone, passenger, mut_effect):
    global MUTATION_ID, sem_mutations
    fit = clone[5][passenger]
    if len(fit) == 0:
        return
    mutations = []
    
    while True:
        if sem_mutations:
            sem_mutations = False
            mutations = [x for x in range(MUTATION_ID, MUTATION_ID + len(fit), 1)]
            MUTATION_ID = MUTATION_ID + len(fit)
            sem_mutations = True
            break
        
    clone[3].extend(mutations)
    (a,b) = clone[4]._shape
    clone[4] = sc.sparse.vstack([clone[4], copy.deepcopy(clone[4][passenger,:])]).tocsr()
    clone[4]._shape = (clone[4]._shape[0], clone[4]._shape[1] + len(fit))
    (c,d) = clone[4]._shape
    clone[4][range(a,c,1), range(b,d,1)] = 1
    clone[1] = clone[1] + len(fit)
    clone[5] = np.append(clone[5], fit * (1 + abs(mut_effect[1]))**sgn(mut_effect[1]))

def oneCloneCycle(i, iPop, tau, mdt, mut_prob, mut_effect, print_time):
    death = np.random.exponential(1, i[1])/mdt
    death = np.where(death < tau, False, True)

    time_t = time.time()             
    ## dying cells
    dying(i, death)    
    tx[1] = tx[1] + (time.time() - time_t) 
    ty[1] = ty[1] + 1
    
    divide = np.random.exponential(1, i[1])/i[5]     
    divide = np.where(divide < tau, True, False)
    m_d = np.random.binomial(1, mut_prob[0], i[1])
    m_d = np.array(m_d, dtype=bool)
    m_p = np.random.binomial(1, mut_prob[1], i[1])
    m_p = np.array(m_p, dtype=bool)
    m_d = np.logical_and(m_d, divide)
    m_p = np.logical_and(m_p, divide)
    divide[m_d] = False
    divide[m_p] = False
    
    time_t = time.time()     
    ## new clones
    newClone(i, m_d, iPop, mut_effect)
    
    tx[2] = tx[2] + (time.time() - time_t) 
    ty[2] = ty[2] + 1
    time_t = time.time()     
    ## mean fitness
    newMutation(i, m_p, mut_effect)
    
    tx[3] = tx[3] + (time.time() - time_t) 
    ty[3] = ty[3] + 1
    
    divide = np.append(divide, np.logical_not(m_p[m_p]))
    
    time_t = time.time()     
    ## division
    division(i, divide)
    
    tx[0] = tx[0] + (time.time() - time_t) 
    ty[0] = ty[0] + 1    
    
    if print_time:
        # if i[1] > 0:
        #     print("Thread: %i, Clone ID: %i, Population: %i, Mutations: %i, Mean mutation number: %i, Fitness: %.3f" % (threading.get_ident(), i[0], i[1], i[4].count_nonzero(), i[4].count_nonzero()/i[1], np.mean(i[5])))
        # else:
        #     print("Clone ID: %i, Population: %i" % (i[0], i[1]))
                ##delete mutation with 0 occurences
        deleteZeroColumn(i)
        
def cloneCycles(i, iPop, tau, mdt, mut_prob, mut_effect, print_time):
    for x in range(i[0], i[1], 1):
        oneCloneCycle(iPop[x], iPop, tau, mdt, mut_prob, mut_effect, print_time)

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
    if resume:
        max_cl = 0
        max_mut = 0
        for x in iPop:
            if x[0] > max_cl:
                max_cl = x[0]
            try:
                if max(x[3]) > max_mut:
                    max_mut = max(x[3])
            except ValueError:
                continue
                
            try:
                if max(x[2]) > max_mut:
                    max_mut = max(x[2])
            except ValueError:
                continue
            
        CLONE_ID = max_cl + 1
        MUTATION_ID = max_mut + 1
    
    popSize = sum([row[1] for row in iPop])
    mdt = popSize/cap
    
    th_ids = []
    pop_len = math.ceil(len(iPop)/threads)
    th_ids = [x for x in range(0, len(iPop), pop_len)]
    th_ids.append(len(iPop))
    t = []
    for i in range(len(th_ids) - 1):
        t.append([th_ids[i], th_ids[i+1]])
    
    develope = []
    for i in t:
        develope.append(Thread(target=cloneCycles, args=(i, iPop, tau, mdt, mut_prob, mut_effect, print_time)))
        # cloneCycles(i, iPop, tau, mdt, mut_prob, mut_effect, print_time)
        develope[len(develope) - 1].start()        
        
    for i in develope:
        i.join()
        
    for i in iPop:
        if i[1] == 0:
            iPop.remove(i)
        
    return iPop
    