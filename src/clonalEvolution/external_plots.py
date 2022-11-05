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
import pandas as pd
import numpy as np
import os
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import scipy as sc

def clonePlot(paths_in, paths_out, ids):
    prop = np.transpose([paths_in, ids])
    for path_in, _id in prop:
        if not path_in.endswith('.txt'):
            print("Wrong file")
            return
        
        path_out = paths_out + '/Figures/' + _id
        if not paths_out.endswith('/'):
            path_out = path_out + '/'
            
        df = []
        if path_in.endswith('.csv'):
            df = pd.read_csv(path_in)
        elif path_in.endswith('.txt'):
            df = loadFile(path_in)
        elif path_in.endswith('.mtx'):
            df = pd.read_csv(path_in)
            
            mm = []
            for row in df.iterrows():
                mm.append(sc.sparse.load_npz(row[1]['Mutation matrix']))
                
            df['Mutation matrix'] = mm

            for i in range(len(df)):
                try:
                    df.at[i,'Driver mutation list'] = [int(x.strip('[]')) for x in df.loc[i,'Driver mutation list'].split(',')]
                except ValueError:
                    df.at[i,'Driver mutation list'] = []
                    
                try:
                    df.at[i,'Uniqal passenger mutation list'] = [int(x.strip('[]')) for x in df.loc[i,'Uniqal passenger mutation list'].split(',')]
                except ValueError:
                    df.at[i,'Uniqal passenger mutation list'] = []
        else:
            break 
        
def mutWavePlot(paths_in, paths_out, ids):
    prop = np.transpose([paths_in, ids])
    for path_in, _id in prop:
        if not path_in.endswith('.txt'):
            print("Wrong file")
            return
        
        path_out = paths_out + '/Figures/' + _id
        if not paths_out.endswith('/'):
            path_out = path_out + '/'
            
        df = []
        if path_in.endswith('.csv'):
            df = pd.read_csv(path_in)
        elif path_in.endswith('.txt'):
            df = loadFile(path_in)
        elif path_in.endswith('.mtx'):
            df = pd.read_csv(path_in)
            
            mm = []
            for row in df.iterrows():
                mm.append(sc.sparse.load_npz(row[1]['Mutation matrix']))
                
            df['Mutation matrix'] = mm

            for i in range(len(df)):
                try:
                    df.at[i,'Driver mutation list'] = [int(x.strip('[]')) for x in df.loc[i,'Driver mutation list'].split(',')]
                except ValueError:
                    df.at[i,'Driver mutation list'] = []
                    
                try:
                    df.at[i,'Uniqal passenger mutation list'] = [int(x.strip('[]')) for x in df.loc[i,'Uniqal passenger mutation list'].split(',')]
                except ValueError:
                    df.at[i,'Uniqal passenger mutation list'] = []
        else:
            break        

def fitWavePlot(paths_in, paths_out, ids):
    prop = np.transpose([paths_in, ids])
    for path_in, _id in prop:
        if not path_in.endswith('.txt'):
            print("Wrong file")
            return
        
        path_out = paths_out + '/Figures/' + _id
        if not paths_out.endswith('/'):
            path_out = path_out + '/'
            
        df = []
        if path_in.endswith('.csv'):
            df = pd.read_csv(path_in)
        elif path_in.endswith('.txt'):
            df = loadFile(path_in)
        elif path_in.endswith('.mtx'):
            df = pd.read_csv(path_in)
            
            mm = []
            for row in df.iterrows():
                mm.append(sc.sparse.load_npz(row[1]['Mutation matrix']))
                
            df['Mutation matrix'] = mm

            for i in range(len(df)):
                try:
                    df.at[i,'Driver mutation list'] = [int(x.strip('[]')) for x in df.loc[i,'Driver mutation list'].split(',')]
                except ValueError:
                    df.at[i,'Driver mutation list'] = []
                    
                try:
                    df.at[i,'Uniqal passenger mutation list'] = [int(x.strip('[]')) for x in df.loc[i,'Uniqal passenger mutation list'].split(',')]
                except ValueError:
                    df.at[i,'Uniqal passenger mutation list'] = []
        else:
            break 

def mullerPlot(path_in, path_out):   
    '''
    path_in - absolute path to file with .txt extension (simulation data is saved into txt file)
    path_out - absolute path to figures save localization
    '''
    
    if not path_in.endswith('.txt'):
        print("Wrong file")
        return
    
    if not path_out.endswith('/'):
        path_out = path_out + '/'
    
    _F = open(path_in, 'r')
    _L = _F.readlines()
    
    clones = []
    ret = []
    first = True
    for line in _L:
        if first:
            first = False
            continue
        a = []
        x = line.split(')')
        for i in x:
            if i == "\n":
                break
            i = i.strip('()')
            i = i.split(',')
            a.append(list(map(lambda w : int(w), i)))
            if int(i[0]) not in clones:
                clones.append(int(i[0]))
        ret.append(a)
        
    Population = []
    
    for i in range(len(ret)):
        temp = [i]
        temp.extend([0 for x in range(len(clones))])
        
        Population.append(temp)
        for x in ret[i]:
            Population[i][clones.index(x[0])+1] = x[1]
    
    ## Generation number, clone number, cells count
    Population = pd.DataFrame(Population)
    temp = ["Generation"]
    temp.extend([str(x) for x in clones])
    Population.columns = temp
    ax = Population.plot.bar(x="Generation", stacked=True, legend=False, width=1, figsize=(40,20))
    
    fig = ax.get_figure()
    try:
        os.makedirs(path_out, exist_ok=True) 
    except OSError as error:
        print(error)
    finally:
        if os.path.exists(path_out + "muller_plot.jpg"):
            os.remove(path_out + "muller_plot.jpg")
        fig.savefig(path_out + "muller_plot.jpg")
        plt.close(fig)

def VAFdecode_ar(iMutations, pop):
    VAF = []
    temp = []
    for x in iMutations:
        if x == 0:
            continue
        if x not in temp:
            temp.append(x)
            VAF.append(1)
        else:
            w = temp.index(x)
            VAF[w] = VAF[w] + 1            
        
    return [x/pop for x in VAF]

def loadFile(path):   
    _f = open(path)
    _l = _f.readlines()
    df = []
    columns = []
    
    for line in _l:
        if not columns:
            columns = ((line.strip(',\n')).split(','))
            continue
        
        line = line.replace('[', ':')
        line = line.replace(']', ':')
        line = line.split(':')
        df.append([])
        idx = 0
        for i in range(len(line)):
            try:
                if i == 0:
                    t = line[i].split(',')
                    idx = int(t[0])
                    df[idx].append(int(t[1]))
                    df[idx].append(int(t[2]))
                    df[idx].append(float(t[3]))
                    df[idx].append(float(t[4]))
                elif i == 1:
                    t = line[i].strip("[]\",")
                    t = t.split(',')
                    x = []
                    for a in t:
                        if not a:
                            break
                        a = a.strip("[]\"")
                        x.append(int(a))
                    df[idx].append(x)
                elif i == 3:
                    t = line[i].strip("[]")
                    t = t.split(',')
                    x = []
                    for a in t:
                        if not a:
                            break
                        x.append(int(a))
                    df[idx].append(x)
                elif i == 4:
                    t = line[i].strip('\",')
                    df[idx].append(int(t))
            except:
                df.pop()
                break
    df = pd.DataFrame(df)
    df.columns = columns
    _f.close()
    return df

def binnedHist(paths_in, paths_out, ids):
    '''
    path_in - absolute path to file with .txt extension (binned simulation data is saved into txt file)
    path_out - absolute path to figures save localization
    '''
    prop = np.transpose([paths_in, ids])
    for path_in, _id in prop:
        if not path_in.endswith('.txt'):
            print("Wrong file")
            return
        
        path_out = paths_out + '/Figures/' + _id
        if not paths_out.endswith('/'):
            path_out = path_out + '/'
        
        df = loadFile(path_in)
        
        pass_m_l = []
        driv_m_l = []
        
        for row in df.iterrows():
            pass_m_l.extend(row[1]["Passener mutation list"])
            driv_m_l.extend(row[1]["Driver mutation list"])
        
        pass_m_l = np.unique(pass_m_l)
        driv_m_l = np.unique(driv_m_l)
        
        ##TODO whole population VAF
        vaf_data = []
        uniq_muts = pass_m_l.tolist()
        uniq_muts.extend(driv_m_l.tolist())
        
        freq_p = np.zeros([len(df), len(uniq_muts)]).tolist()
        clones = []
        row_idx = 0
        for row in df.iterrows():
            x = row[1]['Passener mutation list']
            for mut in x:
                freq_p[row_idx][uniq_muts.index(mut)] = freq_p[row_idx][uniq_muts.index(mut)]+1
        
            x = row[1]['Driver mutation list']
            for mut in x:
                freq_p[row_idx][uniq_muts.index(mut)] = freq_p[row_idx][uniq_muts.index(mut)]+row[1]['Cells number']
            row_idx = row_idx + 1
                    
        popSize = sum(df['Cells number'])
        
        freq = []
        for row in range(len(freq_p)):
            freq.append([])
            freq_t = np.array(freq_p[row])/popSize
            freq[row] = freq_t
        t = np.sum(freq, axis=0)
    
        for row in range(len(freq)):
            freq[row] = np.around(freq[row]/t,decimals=3).tolist()
            
        freq.append([])
        freq[len(freq)-1] = t.tolist()    
    
        freq = np.array(freq)
        freq = freq.T
        freq = freq[np.argsort(freq[:,len(freq[1,:])-1])]
        freq = freq.tolist()
        
        bar_x = [x/200 for x in range(0,201)]
        bar_plot = np.zeros([len(freq[0]) - 1,len(bar_x)])
        idx_bar = 1
        for i in freq:
            while i[len(i)-1] > bar_x[2*idx_bar]:
                idx_bar = idx_bar + 1
            for t in range(len(i)-1):
                bar_plot[t,2*idx_bar-1] = bar_plot[t,2*idx_bar-1] + i[t]
    
        bar_plot = bar_plot.tolist()
        bar_plot.append(bar_x)
        bar_plot = np.array(bar_plot)
        bar_plot = bar_plot.T
        
        bar_plot = pd.DataFrame(bar_plot)
        temp = []
        for row in df.iterrows():
            temp.append(str(row[1]['Clone number']))
        temp.append('id')
        bar_plot.columns = temp
    
        for x in range(0,2):
    
            ax = bar_plot.plot.bar(x='id', stacked=True, legend=False, width = 2, figsize=(40,20))
            if x == 0:
                ax.set_ylim(0,100)
            ax.set_xlabel("VAF", labelpad=50, fontdict={'fontsize':50})
            ax.set_ylabel("Frequency", labelpad=50, fontdict={'fontsize':50})
            ax.set_title("Population VAF, Population: %i" % popSize, pad=50, fontdict={'fontsize':70})
            ax.xaxis.set_major_locator(ticker.MultipleLocator(20))
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(2))
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(40) 
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(40) 
            # freq_plot = pd.Series(freq[len(freq)-1])
            # ax = freq_plot.plot.hist(bins=[x/100 for x in range(0,101)], ylim=(0,50), xlim=(-0.1,1.1), title=("Clone ID: %s, Population: %i" % ('None', popSize)))
            fig = ax.get_figure()
            try:
                os.makedirs(path_out, exist_ok=True) 
            except OSError as error:
                print(error)
            finally:
                if os.path.exists(path_out + "clone_mutations_VAF_%i.jpg"%x):
                    os.remove(path_out + "clone_mutations_VAF_%i.jpg"%x)
                fig.savefig(path_out + "clone_mutations_VAF_%i.jpg"%x)
                plt.close(fig)
        
    
        b = 0
        for row in df.iterrows():
            if row[1]['Cells number'] > 1:#0.05*popSize: 
                w = row[1]['Driver mutation list']
                x = row[1]['Passener mutation list']
                
                freq = []
                uniq = []
                for i in x:
                    if i not in uniq:
                        uniq.append(i)
                        freq.append(1)
                    else:
                        freq[uniq.index(i)] = freq[uniq.index(i)] + 1
                
                for i in w:
                    freq.append(row[1]['Cells number'])
                    
                for i in range(len(freq)):
                    freq[i] = freq[i]/row[1]['Cells number']
                
                freq = pd.Series(freq)
                ax = freq.plot.hist(bins=100, ylim=(0,50), xlim=(0,1), figsize=(40,20))
                ax.set_xlabel("VAF", labelpad=50, fontdict={'fontsize':50})
                ax.set_ylabel("Frequency", labelpad=50, fontdict={'fontsize':50})
                ax.set_title("Clone ID: %i, Population: %i" % (int(row[1]['Clone number']), int(row[1]['Cells number'])), pad=50, fontdict={'fontsize':70})
                for tick in ax.xaxis.get_major_ticks():
                    tick.label.set_fontsize(40) 
                for tick in ax.yaxis.get_major_ticks():
                    tick.label.set_fontsize(40) 
                
                fig = ax.get_figure()
                try:
                    os.makedirs(path_out, exist_ok=True) 
                except OSError as error:
                    print(error)
                finally:
                    if os.path.exists(path_out + "clone_mutations_VAF_ID_" + str(int(row[1]['Clone number'])) +".jpg"):
                        os.remove(path_out + "clone_mutations_VAF_ID_" + str(int(row[1]['Clone number'])) +".jpg")
                    fig.savefig(path_out + "clone_mutations_VAF_ID_" + str(int(row[1]['Clone number'])) +".jpg")
                    plt.close(fig)
            else:
                continue

def matrixHist(paths_in, paths_out, ids):
    '''
    path_in - absolute path to file with .txt extension (binned simulation data is saved into txt file)
    path_out - absolute path to figures save localization
    '''
    prop = np.transpose([paths_in, ids])
    for path_in, _id in prop:    
        if not path_in.endswith('.mtx'):
            print("Wrong file")
            return
        
        path_out = paths_out + '/Figures/' + _id
        if not path_out.endswith('/'):
            path_out = path_out + '/'
    
        df = pd.read_csv(path_in)
        
        mm = []
        for row in df.iterrows():
            mm.append(sc.sparse.load_npz(row[1]['Mutation matrix']))
            
        df['Mutation matrix'] = mm

        for i in range(len(df)):
            try:
                df.at[i,'Driver mutation list'] = [int(x.strip('[]')) for x in df.loc[i,'Driver mutation list'].split(',')]
            except ValueError:
                df.at[i,'Driver mutation list'] = []
                
            try:
                df.at[i,'Uniqal passenger mutation list'] = [int(x.strip('[]')) for x in df.loc[i,'Uniqal passenger mutation list'].split(',')]
            except ValueError:
                df.at[i,'Uniqal passenger mutation list'] = []
                
        popSize = sum(df['Cells number'])
        
        pIDs = []
        for i in df['Uniqal passenger mutation list']:
            pIDs.extend(i)
        for i in df['Driver mutation list']:
            pIDs.extend(i)
        pIDs = np.unique(pIDs).tolist()
        pFreq = np.array([np.zeros(len(pIDs)) for x in range(len(df)+1)])
        
        for i in range(len(pFreq)-1):
            muts = df['Uniqal passenger mutation list'][i] 
            driv = df['Driver mutation list'][i]
            freqs = df['Mutation matrix'][i]
            freqs = [freqs[:,x].count_nonzero() for x in range(freqs._shape[1])]
            muts = np.array(muts)[0:len(freqs)].tolist()
            for m in muts:
                a = freqs[muts.index(m)]     
                pFreq[i,pIDs.index(m)] = a    
            for m in driv:
                a = df['Cells number'][i]
                pFreq[i,pIDs.index(m)] = a
        
        for i in range(len(pFreq[0,:])):
            pFreq[len(pFreq[:,0])-1,i] = sum(pFreq[:,i])
        
        pFreq = pFreq.T
        pFreq = pFreq[pFreq[:,1] != 0]
        pFreq = pFreq.T
        if not np.any(pFreq):
            print(_id + ": no value to plot")
            continue
        for i in range(len(pFreq[:,0])-1):
            for j in range(len(pFreq[0,:])):
                pFreq[i,j] = round(pFreq[i,j]/pFreq[len(pFreq[:,0])-1,j],2)
            
        pFreq[len(pFreq[:,0])-1,:] = pFreq[len(pFreq[:,0])-1,:]/popSize    
        pFreq = pFreq.T
        pFreq = pFreq[np.argsort(pFreq[:,len(pFreq[0,:])-1])]
        pFreq = pFreq.tolist()
        
        t_freq = [x/200 for x in range(0,201)]
        freq = np.zeros((len(pFreq[0]),len(t_freq)))
        freq[len(freq[:,0])-1,:] = t_freq
        idx = 1
        for i in pFreq:
            while i[len(i)-1] > t_freq[2*idx]:
                idx = idx + 1
                if idx == 101:
                    break
            else:
                if i[len(i)-1] == 0:
                    continue
                for x in range(len(i)-1):
                    freq[x,2*idx - 1] = freq[x,2*idx - 1] + i[x]
        
        pl = pd.DataFrame(freq.T)
        names = [str(x) for x in df['Clone number']]
        names.append('sum')
        pl.columns = names
        
        for x in range(0,2):
            ax = pl.plot.bar(x='sum', stacked=True, width = 2, figsize=(40,20))
            ax.legend(prop={'size':30})
            if x == 0:
                ax.set_ylim(0,100)
            ax.set_xlabel("VAF", labelpad=50, fontdict={'fontsize':50})
            ax.set_ylabel("Frequency", labelpad=50, fontdict={'fontsize':50})
            ax.set_title("Population VAF, Population: %i" % popSize, pad=50, fontdict={'fontsize':70})
            ax.xaxis.set_major_locator(ticker.MultipleLocator(20))
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(2))
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(40) 
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(40) 
            fig = ax.get_figure()
            try:
                os.makedirs(path_out, exist_ok=True) 
            except OSError as error:
                print(error)
            finally:
                if os.path.exists(path_out + "clone_mutations_VAF_%i.jpg"%x):
                    os.remove(path_out + "clone_mutations_VAF_%i.jpg"%x)
                fig.savefig(path_out + "clone_mutations_VAF_%i.jpg"%x)
                plt.close(fig)
        
        for clone in df.iterrows():
            freqs = clone[1]['Mutation matrix']        
            _sum = [freqs[:,x].count_nonzero() for x in range(freqs._shape[1])]
            freqs = [1 for x in range(freqs._shape[1])]
            
            pl = np.zeros(len(freqs)).tolist()
            for i in clone[1]['Driver mutation list']:
                freqs.append(0)
                pl.append(1)
                _sum.append(clone[1]['Cells number'])
                
            freq = np.zeros((3,len(pl)))
            freq[0,:] = np.array(freqs)
            freq[1,:] = np.array(pl)
            freq[2,:] = np.array(_sum)/clone[1]['Cells number']  
            
            freq = freq.T
            freq = freq[np.argsort(freq[:,len(freq[0,:])-1])]
            freq = freq.tolist()
            
            t_freq = [x/200 for x in range(0,201)]
            pfreq = np.zeros((len(freq[0]),len(t_freq)))
            pfreq[len(pfreq[:,0])-1,:] = t_freq
            idx = 1
            for i in freq:
                while i[len(i)-1] > t_freq[2*idx]:
                    idx = idx + 1
                    if idx == 101:
                        break
                else:
                    if i[len(i)-1] == 0:
                        continue
                    for x in range(len(i)-1):
                        pfreq[x,2*idx - 1] = pfreq[x,2*idx - 1] + i[x]
            
            freq = pd.DataFrame(pfreq.T)
            freq.columns = ["passenger","driver","sum"]
            ax = freq.plot.bar(x='sum', stacked=True, xlim=(0,1), width=2, figsize=(40,20))
            ax.legend(prop={'size':30})
            ax.set_xlabel("VAF", labelpad=50, fontdict={'fontsize':50})
            ax.set_ylabel("Frequency", labelpad=50, fontdict={'fontsize':50})
            ax.set_title("Population VAF, Clone ID: %i, Population: %i" % (clone[1]['Clone number'], clone[1]['Cells number']), pad=50, fontdict={'fontsize':70})
            ax.xaxis.set_major_locator(ticker.MultipleLocator(20))
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(2))
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(40) 
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(40)
            
            fig = ax.get_figure()
            try:
                os.makedirs(path_out, exist_ok=True) 
            except OSError as error:
                print(error)
            finally:
                if os.path.exists(path_out + "clone_mutations_VAF_ID_" + str(int(clone[1]['Clone number'])) +".jpg"):
                    os.remove(path_out + "clone_mutations_VAF_ID_" + str(int(clone[1]['Clone number'])) +".jpg")
                fig.savefig(path_out + "clone_mutations_VAF_ID_" + str(int(clone[1]['Clone number'])) +".jpg")
                plt.close(fig)           

def singleHist(paths_in, paths_out, ids):
    '''
    path_in - absolute path to file with .txt extension (binned simulation data is saved into txt file)
    path_out - absolute path to figures save localization
    '''
    prop = np.transpose([paths_in, ids])
    for path_in, _id in prop:    
        if not path_in.endswith('.csv'):
            print("Wrong file")
            return
        
        path_out = paths_out + '/Figures/' + _id
        if not path_out.endswith('/'):
            path_out = path_out + '/'
        
        filepath = Path(path_in)  
        filepath.parent.mkdir(parents=True, exist_ok=True)
        
        df = pd.read_csv(filepath)
        
        b = 0
        _max = len(df)
        perc = -1
        for row in df.iterrows():
            w = row[1]['Mutations_ID']
            w = w.split(', ')
            retVal = []
            for i in w:
                try:
                    _t = int(i.strip('[]'))
                    retVal.append(_t)
                except ValueError:
                    retVal = []
            df['Mutations_ID'].loc[b] = retVal.copy()
            b = b+1
            if perc != round(b/_max * 100, 1):
                perc = round(b/_max * 100, 1)
                print("to int transform %.1f %%" % perc)
        
        VAF = []
        VAF_v = []
        clones = []
        clones_f = []
        
        b = 0
        for i in df[['Clone']].iterrows():
            if int(i[1]['Clone']) not in clones:
                clones.append(int(i[1]['Clone']))
                clones_f.append(1)
            else:
                x = clones.index(int(i[1]['Clone']))
                clones_f[x] = clones_f[x] + 1
            b = b+1
            if perc != round(b/_max * 100, 1):
                perc = round(b/_max * 100, 1)
                print("clone search %.1f %%" % perc)
                
        clones_f = np.array(clones_f)
        clones = np.array(clones)
        idx = np.where(clones_f > 10)[0]
        clones = clones[idx].tolist()
        mutations = [[] for x in clones]
        freq = clones_f[idx]
        
        _max = _max * len(clones)
        b = 0
        for x in clones:
            w = []
            for row in df[['Clone','Mutations_ID']].iterrows():
                if int(row[1]['Clone']) == x:
                    w.extend(row[1]['Mutations_ID'].copy())
                b = b + 1
                if perc != round(b/_max * 100, 1):
                    perc = round(b/_max * 100, 1)
                    print("clone mutations %.1f %%" % perc)
        
            mutations[clones.index(x)] = w
            
        b = 0
        whole_pop = []
        _max = len(df)
        for row in df[['Mutations_ID']].iterrows():
            whole_pop.extend(row[1]['Mutations_ID'].copy())
            b = b + 1
            if perc != round(b/_max * 100, 1):
                perc = round(b/_max * 100, 1)
                print("clone mutations %.1f %%" % perc)
        
            
        VAF = VAFdecode_ar(whole_pop, len(df))
        VAF = pd.DataFrame(VAF)
        try:
            ax = VAF.plot.hist(bins=[x/100 for x in range(0,101)], ylim=(0,50), xlim=(-0.1,1.1), title="Population VAF, cells = %i" % (len(df)))
            ax.set_xlabel("VAF")
            ax.set_ylabel("Mutations")
        except TypeError:
            print(_id + ": no numeric data to plot")
            continue
        
        fig = ax.get_figure()
        try:
            os.makedirs(path_out, exist_ok=True) 
        except OSError as error:
            print(error)
        finally:
            if os.path.exists(path_out + "population_mutations_VAF.jpg"):
                os.remove(path_out + "population_mutations_VAF.jpg")
            fig.savefig(path_out + "population_mutations_VAF.jpg")
            plt.close(fig)
        
        for i in range(len(clones)):
            VAF = VAFdecode_ar(mutations[i], freq[i])
            VAF = pd.DataFrame(VAF)
            try:
                ax = VAF.plot.hist(bins=[x/100 for x in range(0,101)], ylim=(0,50), xlim=(-0.1,1.1), title="Clone VAF, clone id = %i, clones = %i" % (clones[i], freq[i]))
                ax.set_xlabel("VAF")
                ax.set_ylabel("Mutations")
            except TypeError:
                print(_id + ": no numeric data to plot")
                continue  
            
            fig = ax.get_figure()
            try:
                os.makedirs(path_out, exist_ok=True) 
            except OSError as error:
                print(error)
            finally:
                if os.path.exists(path_out + "clone_mutations_VAF_ID_" + str(clones[i]) +".jpg"):
                    os.remove(path_out + "clone_mutations_VAF_ID_" + str(clones[i]) +".jpg")
                fig.savefig(path_out + "clone_mutations_VAF_ID_" + str(clones[i]) +".jpg")
                plt.close(fig)
                print("saving files %.1f %%" % (i/len(clones) * 100))
            
# if __name__ == "__main__":
#     matrixHist("E:\Simulations\Clonal Evolution\Mutation matrix aproach\mm_400.mtx", 
#                "E:\Simulations\Clonal Evolution\Mutation matrix aproach\Figures\\400")