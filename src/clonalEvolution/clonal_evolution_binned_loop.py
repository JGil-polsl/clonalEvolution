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
import time
import numpy as np
import math
import random
from threading import Thread

MUTATION_ID = 1
CLONE_ID = 1
tCLONE = 0
semafor_1 = False
semafor_2 = False
semafor_3 = False

tx = [0,0,0,0]
ty = [0,0,0,0]

#idea: kalkulowac na raz wszystkie komorki umierajace i dzielace sie (mutacje) sum poisson, prob and random select
    
def division(i, divide, m_list): # copy mutation in passanger mutation list
    if divide > 0 and i[5]:
        passanger_divide = []
        if m_list:
            uni = np.unique(i[5]) # select only unique mutations ids (two same mutations cannot happen in one cell)
            prob = np.bincount(i[5])
            prob = prob[prob!=0]/len(i[5])
            for x in range(divide):
                mut = sum(np.random.poisson(i[3], divide))#math.floor(i[3]) + np.random.binomial(1,i[3]-math.floor(i[3]),1)[0] ## calculate number of mutations duplicate
                new_muts = []
                if mut >= len(uni): ## check if duplicate is greater than mutation list length
                    mut = len(uni)
                try:
                    new_muts = np.random.choice(uni, size=mut, replace=False, p=prob) ## select random mutations to duplicate
                    # new_muts = np.random.choice(i[5], size=mut)
                except IndexError as error:
                    continue
                except ValueError:
                    print("no mutation")
                    continue
                passanger_divide.extend(new_muts)
            i[5].extend(passanger_divide) ## add duplicates to mutation list
            try:
                i[3] = round((len(i[5]))/i[1],7) ## calculate new mean mutation number
            except:
                i[3] = 0
        else:
            mut = divide*math.floor(i[3]) + sum(np.random.binomial(1,i[3]-math.floor(i[3]),divide))
            i[3] = round(((i[1]-divide)*i[3] + mut)/i[1],7)
    return

def dying(i, death, m_list): # delete dying mutations from passanger mutation list
    if death > 0 and i[5]: ## check if list consider all dying cells to have maksimum mutate number 
        if m_list:
            uni = np.unique(i[5]) # select only unique mutations ids (two same mutations cannot happen in one cell)
            freq = np.bincount(i[5])
            freq = freq[freq != 0]
            for y in range(death):
                mut = np.random.poisson(i[3])#math.floor(i[3]) + np.random.binomial(1,i[3]-math.floor(i[3]),1)[0] ## calculate mutation number to erease
                prob = freq/len(i[5])
                if mut >= len(uni): ## check if mutation to erease is greater than mutation list length
                    mut = len(uni)
                try:
                    passenger_death = np.random.choice(uni, size=mut, replace=False, p=prob) ## select random mutation from mutation list
                    # passenger_death = np.random.choice(i[5], size=mut)
                    for x in passenger_death:
                        i[5].remove(x) ## remove mutation
                        freq[uni == x] = freq[uni == x] - 1
                except IndexError as error:
                    continue
                except ValueError:
                    print("no mutation")
                    continue
            try:
                i[3] = round((len(i[5]))/i[1],7) ## calculate new mean mutation number
            except:
                i[3] = 0
        else:
            mut = death*math.floor(i[3]) + sum(np.random.binomial(1,i[3]-math.floor(i[3]),death))
            i[3] = round(((i[1]+death)*i[3] - mut)/i[1],7)
    return

def newMutation(i, m_p, mut_effect, m_list):
    global MUTATION_ID, semafor_1
    if m_p > 0:
        if m_list:
            f_m = i[2]/(1-mut_effect[1])*m_p ## calculate mutated cells fitness
            f_c = i[2]*(i[1] - m_p) ## calculate population fitness
            i[2] = (f_m + f_c)/i[1] ## calculate mean fitness number
            for x in range(MUTATION_ID, MUTATION_ID + m_p):
                i[5].append(x) ## add new mutation id to list
            while semafor_1:
                continue
            semafor_1 = True
            MUTATION_ID = MUTATION_ID + m_p ## update mutation ID
            semafor_1 = False
            try:
                i[3] = round((len(i[5]))/i[1],7) ## calculate new mean mutation number
            except:
                i[3] = 0
        else:
            i[3] = round(((i[1]-m_p)*i[3] + m_p)/i[1],7)
    return

def newClone(i, m_d, iPop, mut_effect, m_list):
    global CLONE_ID, MUTATION_ID, semafor_1, semafor_2
    for x in range(CLONE_ID, CLONE_ID + m_d):
        if m_list:
            mut = np.random.poisson(i[3])#math.floor(i[3]) + np.random.binomial(1,i[3]-math.floor(i[3]),1)[0]  ## calulate mutation number to copy to new clone
            driver = i[4].copy() ## copy list of driver mutations (all are included in new clone)              
            driver.append(MUTATION_ID) ## add new mutation ID to driver mutation list
            while semafor_1:
                continue
            semafor_1 = True
            MUTATION_ID = MUTATION_ID + 1 ## update mutation ID
            semafor_1 = False
            passenger = [] 
            uni = np.unique(i[5]) # select only unique mutations ids (two same mutations cannot happen in one cell)
            prob = np.bincount(i[5])
            prob = prob[prob!=0]/len(i[5])
            if mut >= len(uni):
                mut = len(uni)
            try:
                passenger = np.random.choice(uni, size=mut, replace=False, p=prob) ## select random mutations to copy
            except IndexError as error:
                iPop.append([x, 1, i[2]*(1+mut_effect[0]), 0, driver, [], i[0]]) ## append new clone to population
                continue
            except ValueError:
                iPop.append([x, 1, i[2]*(1+mut_effect[0]), 0, driver, [], i[0]]) ## append new clone to population
                continue
            iPop.append([x, 1, i[2]*(1+mut_effect[0]), len(passenger), driver, passenger.tolist(), i[0]]) ## append new clone to population
        else:
            iPop.append([x, 1, i[2]*(1+mut_effect[0]), i[3], i[4].copy(), [], i[0]])
    while semafor_2:
        continue
    semafor_2 = True
    CLONE_ID = CLONE_ID + m_d ## update clone ID
    semafor_2 = False

def main_loop(iPop, index, mdt, popSize, tau, mut_prob, mut_effect, print_time):
    global tx, ty, CLONE_ID
    m_list = True
    x = [l for l in range(index[0],index[1])]
    for l in x:
        if l >= len(iPop):
            continue
        pdv = 1 - math.exp(-tau*iPop[l][2])
        pdt = 1 - math.exp(-tau*mdt)
        
        pdv = (1-pdt)*pdv
        pdm_d = (mut_prob[0])*(1 - mut_prob[1])
        pdm_p = (1 - mut_prob[0])*(mut_prob[1])
        pdr = (1 - pdt)*(1 - pdv)*(1 - pdm_d)*(1 - pdm_p)
        
        r = np.random.multinomial(iPop[l][1], [pdt, pdv, pdr])
        
        death = r[0]
        divide = r[1]
        
        m = np.random.multinomial(divide, [pdm_p, pdm_d, (1-pdm_p)*(1-pdm_d)])

        m_p = m[0]        
        m_d = m[1]
        
        ## clone size
        iPop[l][1] = iPop[l][1] - death + divide
        

        time_t = time.time()             
        ## dying cells
        dying(iPop[l], death, m_list)
          
        tx[0] = tx[0] + (time.time() - time_t) 
        ty[0] = ty[0] + 1
        time_t = time.time()     
        ## new clones
        newClone(iPop[l], m_d, iPop, mut_effect, m_list)
        
        tx[1] = tx[1] + (time.time() - time_t) 
        ty[1] = ty[1] + 1
        time_t = time.time()     
        ## division
        division(iPop[l], divide, m_list)
        
        tx[2] = tx[2] + (time.time() - time_t) 
        ty[2] = ty[2] + 1
        time_t = time.time()     
        ## mean fitness
        newMutation(iPop[l], m_p, mut_effect, m_list)
        
        tx[3] = tx[3] + (time.time() - time_t) 
        ty[3] = ty[3] + 1
        time_t = time.time()     
        
        try:
            if print_time:
                print("Clone ID: %i, Dying: %.2f, newClone: %.2f, Division: %.2f, newMutation: %.2f" % (CLONE_ID,tx[0]/ty[0],tx[1]/ty[1],tx[2]/ty[2],tx[3]/ty[3]))
                tx = [0,0,0,0]
                ty = [0,0,0,0]
        except ZeroDivisionError:
            continue

def clonalEvolutionBinnedLoop(iPop, cap, tau, mut_prob, mut_effect, resume, q, THREADS, print_time):
    """
    Assumption:
        Clone after obtaining mean mutation number than one have no one cell without mutation (same for 2,3,4 etc.)
        
    Description:
        One cycle to update population - tau loop binned method
        Prameters:
            iPop: population matrix where row is in form of:
                Clone number
                Cells number
                Mean fitness
                Mean mutation number
                Driver mutation list
                Passener mutation list
                Previous clone number
            cap: population capacity
            tau: tau step
            mut_prob: list in form of: [driver mutation probability, passenger mutation probability]
            mut_effect: list in form of: [driver mutation effect, passenger muatation effect]
            resume: acknowledge to resume simulation
            q: common queue
            THREADS: threads number used in simulation         
    """
    global CLONE_ID, MUTATION_ID
    popSize = sum([row[1] for row in iPop])
    mdt = popSize/cap
    # mdt = math.log(1 + (math.e - 1)*popSize/cap)
    
    if resume:
        m_max = 0
        c_max = 0
        for row in iPop:
            if row[4] and row[5]:
                if(max(row[4]) > m_max):
                    m_max = max(row[4])
                if(max(row[5]) > m_max):
                    m_max = max(row[5])
            if row[0] > c_max:
                c_max = row[0]
        MUTATION_ID = m_max + 1
        CLONE_ID = c_max + 1
    
    up = len(iPop)
    
    th = math.ceil(up/THREADS)
    idx = [(int(x*th),int((x+1)*th)) for x in range(THREADS)]
    idx[len(idx)-1] = (idx[len(idx) - 1][0], up)  
    
    develop = []
    
    for i in range(len(idx)):
        develop.append(Thread(target=main_loop, args=(iPop, idx[i], mdt, popSize, tau, mut_prob, mut_effect, print_time)))
        develop[i].start()
        # main_loop(iPop, idx[i], mdt, popSize, tau, mut_prob, mut_effect)
        
    for i in develop:
        i.join()
        
    for i in iPop:
        if i[1] == 0:
            iPop.remove(i)
            
    return iPop