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

import PyQt5.QtWidgets as qtWidget
import PyQt5.QtGui as qtGui
from PyQt5.QtCore import Qt
import sys
import pandas as pd
import numpy as np
import re
import scipy as sc

disclaimer = '''Cellular/Microbial Clonal Evolution simulations basing on Gillespie algorithm. Copyright (C) 2022 by Jarosław Gil. 
    This program comes with ABSOLUTELY NO WARRANTY;.
    This is free software, and you are welcome to redistribute it
    under certain conditions;'''


from pathlib import Path 
import time
import matplotlib.pyplot as plt
from threading import Thread
from multiprocessing import Process, Queue
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import copy

import clonalEvolution.external_plots as external_plots
from clonalEvolution.clonal_evolution_init import clonalEvolutionMainLoop 

class mainFormat(qtWidget.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.ID = 1
        self.idx_s = 0
        self.idx_r = 0
        self.th_r = []
        self.r_ID = []
        self.th_s = []
        self.s_ID = []
        self.q = Queue()
        self._th_muller = []
        self._th_hist = []
        self._th_mw = []
        self._th_fw = []
        self._th_cp = []
        self.createMainView()
        self._monitor = True
        self.th_monitor = Thread(target=self.monitor, args=(self.q, self.status))
        self.th_monitor.start()
        
    def about(self):
        global disclaimer
        self.showDialog(disclaimer, 'About')
    
    def plot(self, ID, mutWave, select):
        self._fig.clear()
        ax = self._fig.add_subplot(111)     
        lab = []
        if select == 0:
            last = []
            first = True
            for x in mutWave:
                if x == "mut_num":
                    lab = mutWave['mut_num'].tolist()
                else:  
                    if first:
                        last = np.array(mutWave[x].tolist())
                        ax.bar(lab, last)
                        first = False
                    else:
                        data = np.array(mutWave[x].tolist())
                        ax.bar(lab, data, bottom=last)
                        last = last + data
            ax.set_title('Mutation Wave, ID: %s' % ID)
            for i in lab:
                i = str(int(i))
            ax.set_xticks(mutWave['mut_num'], lab)
            leg = []
            first = True
            for x in mutWave:
                if first:
                    first = False
                    continue
                txt = "%s: %.0f" % (x, sum(mutWave[x].tolist()))
                leg.append(txt)
            print(leg)
            ax.legend(labels=leg, loc='center left', bbox_to_anchor=(1, 0.5), fontsize="small")  
            self._canvas.draw()
        elif select == 1:
            ax.bar([x for x in range(len(mutWave['Clone']))], mutWave['Cells number'])
            ax.set_title("Population: %i, Clones: %i" % (sum(mutWave['Cells number']), len(mutWave['Clone'])))
            ax.set_xticks([x for x in range(len(mutWave['Clone']))], [str(x) for x in mutWave['Clone']])
            self._canvas.draw()
            self.status.setText(ID + ": PLOTTED")
    
    def monitor(self, q, status):
        while True:
            if not q.empty():
                data = q.get()
                if data[0] == 'exit':
                    print('exit monitor')
                    break
                elif data[0] == '0':
                    status.setText(data[1] + ": " + data[2])
                elif data[0] == "-1":
                    status.setText(data[1] + ": plot")
                    Thread(target=self.plot, args=(data[1],data[2], 0)).start()
                elif data[0] == "-2":
                    status.setText(data[1] + ": plot")
                    Thread(target=self.plot, args=(data[1],data[2], 1)).start()
                elif data[0] == 'ended':
                    ID = data[1]
                    try:
                        temp = list(map(lambda x: (x.split(',')[0]), self.s_ID))
                        index = temp.index(str(ID))
                        self.s_ID.remove(self.s_ID[index])
                        self.th_s[index].join()
                        self.th_s.remove(self.th_s[index])
                        self.idx_s = self.idx_s - 1
                    except:
                        temp = list(map(lambda x: (x.split(',')[0]), self.r_ID))
                        index = temp.index(str(ID))
                        self.r_ID.remove(self.r_ID[index])
                        self.th_r[index].join()
                        self.th_r.remove(self.th_r[index])
                        self.idx_r = self.idx_r - 1
                else:
                    q.put(data)
            if self._th_muller:
                if self._th_muller.is_alive():
                    self._th_muller.join()                    
                    self.status.setText("muller plot done")
                    self._th_muller = None
            if self._th_hist:
                if self._th_hist.is_alive():
                    self._th_hist.join()                    
                    self.status.setText("histograms done")
                    self._th_hist = None
    
    def closeEvent(self, event):
        
        for i in self.s_ID:
            self.q.put(['1', str(i).split(',')[0], "exit"])
        for i in self.r_ID:
            self.q.put(['1', str(i).split(',')[0], "exit"])
        
        for i in self.th_r:
            i.join()
        for i in self.th_s:
            i.join()
            
        if self._th_muller:
            self._th_muller[0].join()
        if self._th_hist:
            self._th_hist[0].join()
            
        plt.close(self._fig)
        self.q.put(['exit'])
        self.th_monitor.join()
        event.accept()  
    
    def update(self):
        i = self.tabs.indexOf(self._tabUI)
        self.tabs.removeTab(i)
        self._tabUI = self.threadTabUI()
        self.tabs.addTab(self._tabUI, "Active threads")
    
    def stopSim(self):
        msg = self._msg.text()
        msg = msg.split(',')
        ID = msg[0]
        msg = msg[1]
        self.q.put(['1',ID,msg])
        index = 0
        if msg == "exit":
            try:
                temp = list(map(lambda x: (x.split(',')[0]), self.s_ID))
                index = temp.index(str(ID))
                self.s_ID.remove(self.s_ID[index])
                self.th_s[index].join()
                self.th_s.remove(self.th_s[index])
                self.idx_s = self.idx_s - 1
            except:
                temp = list(map(lambda x: (x.split(',')[0]), self.r_ID))
                index = temp.index(str(ID))
                self.r_ID.remove(self.r_ID[index])
                self.th_r[index].join()
                self.th_r.remove(self.th_r[index])
                self.idx_r = self.idx_r - 1
            
            self.status.setText("Stopped")
            if not self.th_r and not self.th_s:
                self.status.setStyleSheet("background-color: red")
    
    def createCanvasView(self):
        self.canvasTab = qtWidget.QWidget()
        _canvasTab = qtWidget.QVBoxLayout()
        
        self._fig = plt.figure()
        self._canvas = FigureCanvas(self._fig)
        self._toolbar = NavigationToolbar(self._canvas, self)
        
        _canvasTab.addWidget(self._toolbar)
        _canvasTab.addWidget(self._canvas)
        
        self.canvasTab.setLayout(_canvasTab)
        return self.canvasTab  
    
    def createMainView(self):
        layout = qtWidget.QVBoxLayout()
        # Create the tab widget with two tabs
        self.tabs = qtWidget.QTabWidget()
        self.tabs.addTab(self.generalTabUI(), "General")
        self.tabs.addTab(self.generatedTabUI(), "Generated")
        self.tabs.addTab(self.parametersTabUI(), "Parameters")
        self.tabs.addTab(self.createCanvasView(), "Plot")
        self._tabUI = self.threadTabUI()
        self.tabs.addTab(self._tabUI, "Active threads")
        layout.addWidget(self.tabs)
        layout.addStretch()
        
        self._single = qtWidget.QRadioButton("Checked: single cell version")
        self._single.setChecked(True)
        self._binned = qtWidget.QRadioButton("Checked: binned version")
        self._matrix = qtWidget.QRadioButton("Checked: matrix version")
        layout.addWidget(self._single)
        layout.addWidget(self._binned)
        layout.addWidget(self._matrix)
        
        start = qtWidget.QPushButton(self)
        start.setText("Start simulation")
        start.clicked.connect(self.simStart)
        layout.addWidget(start)
        
        self.status = qtWidget.QLabel()
        self.status.setStyleSheet("background-color: red")
        self.status.setText("Stopped")
        self.status.setAlignment(Qt.AlignCenter)
        self.status.setFont(qtGui.QFont('Consolas', 15))
        layout.addWidget(self.status)
        
        send_msg = qtWidget.QPushButton(self)
        send_msg.setText("Send message")
        send_msg.clicked.connect(self.stopSim)
        layout.addWidget(send_msg)
        
        layout.addWidget(qtWidget.QLabel("message format: thread_ID,command (without space)"))
        layout.addWidget(qtWidget.QLabel("commands: exit - stop simulation, size - get population size, time - get simulation time"))
        layout.addWidget(qtWidget.QLabel("plot - shows plot (not binned - mutation wave, binned - clonal plot), save - save to file, clone - number of clones"))
        
        self._msg = qtWidget.QLineEdit()
        self._msg.returnPressed.connect(self.stopSim)
        layout.addWidget(self._msg)
        
        _about = qtWidget.QPushButton(self)
        _about.setText("about")
        _about.clicked.connect(self.about)
        layout.addWidget(_about)
        
        self.setLayout(layout)
        
    def mullerPlotAction(self):
        try:
            if self._file_path.text() == "":
                raise Exception()
            if self._th_muller:
                raise NameError()
        except NameError:
            self.showDialog("Plotting already", "Alert")
            return
        except:
            self.showDialog("Enter save path", "Alert")
            return
        fname = qtWidget.QFileDialog.getOpenFileName(self, 'Open file', "Z://","CSV files (*.txt)")[0]  
        if fname == "":
            self.showDialog("No file selected!", "Alert")
            return
        self._th_muller = (Thread(target=external_plots.mullerPlot, args=(fname, self._file_path.text() + "/Figures")))
        # external_plots.mullerPlot(fname, self._file_path.text() + "/Figures")
        self._th_muller.start()
        self.status.setText("muller plot ongoing")
        
    def cloneHistAction(self):
        try:
            if self._file_path.text() == "":
                raise Exception()
            if self._th_hist:
                raise NameError()
        except NameError:
            self.showDialog("Plotting already", "Alert")
            return
        except:
            self.showDialog("Enter save path", "Alert")
            return
        fname = qtWidget.QFileDialog.getOpenFileNames(None, 'Open file', "Z://","Single data files (*.csv);; Binned data files (*.txt);; Matrix data files (*.mtx)")[0] 
        if fname == "":
            self.showDialog("No file selected!", "Alert")
            return
        name_id = []
        for i in fname:
            name_id.append(re.findall('\d+', i)[-1])
        if all(map(lambda x: x.endswith('.txt'), fname)):
            self._th_hist = (Thread(target=external_plots.binnedHist, args=(fname, self._file_path.text(), name_id)))
            # external_plots.binnedHist(fname, self._file_path.text(), name_id)
            self._th_hist.start()
            self.status.setText("mutations histograms ongoing")
        elif all(map(lambda x: x.endswith('.csv'), fname)):
            self._th_hist = (Thread(target=external_plots.singleHist, args=(fname, self._file_path.text(), name_id)))
            # external_plots.singleHist(fname, self._file_path.text(), name_id)
            self._th_hist.start()
            self.status.setText("mutations histograms ongoing")
        elif all(map(lambda x: x.endswith('.mtx'), fname)):
            self._th_hist = (Thread(target=external_plots.matrixHist, args=(fname, self._file_path.text(), name_id)))
            # external_plots.matrixHist(fname, self._file_path.text(), name_id)
            self._th_hist.start()
            self.status.setText("mutations histograms ongoing")
        else:
            self.showDialog("Wrong file/files extension. Use only one kind at time", "Alert")
            
    def mutWaveAction(self):
        try:
            if self._file_path.text() == "":
                raise Exception()
            if self._th_mw:
                raise NameError()
        except NameError:
            self.showDialog("Plotting already", "Alert")
            return
        except:
            self.showDialog("Enter save path", "Alert")
            return
        
        fname = qtWidget.QFileDialog.getOpenFileNames(None, 'Open file', "Z://","Single data files (*.csv);; Binned data files (*.txt);; Matrix data files (*.mtx)")[0] 
        if fname == "":
            self.showDialog("No file selected!", "Alert")
            return
        name_id = []
        for i in fname:
            name_id.append(re.findall('\d+', i)[-1])
            
        # self._th_mw = (Thread(target=external_plots.mutWavePlot, args=(fname, self._file_path.text(), name_id)))
        external_plots.mutWavePlot(fname, self._file_path.text(), name_id)
        # self._th_mw.start()
        # self.status.setText("mutation wave ongoing")
    
    def fitWaveAction(self):
        try:
            if self._file_path.text() == "":
                raise Exception()
            if self._th_fw:
                raise NameError()
        except NameError:
            self.showDialog("Plotting already", "Alert")
            return
        except:
            self.showDialog("Enter save path", "Alert")
            return
        
        fname = qtWidget.QFileDialog.getOpenFileNames(None, 'Open file', "Z://","Single data files (*.csv);; Binned data files (*.txt);; Matrix data files (*.mtx)")[0] 
        if fname == "":
            self.showDialog("No file selected!", "Alert")
            return
        name_id = []
        for i in fname:
            name_id.append(re.findall('\d+', i)[-1])
            
        # self._th_fw = (Thread(target=external_plots.fitWavePlot, args=(fname, self._file_path.text(), name_id)))
        external_plots.fitWavePlot(fname, self._file_path.text(), name_id)
        # self._th_fw.start()
        # self.status.setText("mutation wave ongoing")
        
    def clonePlotAction(self):
        try:
            if self._file_path.text() == "":
                raise Exception()
            if self._th_cp:
                raise NameError()
        except NameError:
            self.showDialog("Plotting already", "Alert")
            return
        except:
            self.showDialog("Enter save path", "Alert")
            return
        
        fname = qtWidget.QFileDialog.getOpenFileNames(None, 'Open file', "Z://","Single data files (*.csv);; Binned data files (*.txt);; Matrix data files (*.mtx)")[0] 
        if fname == "":
            self.showDialog("No file selected!", "Alert")
            return
        name_id = []
        for i in fname:
            name_id.append(re.findall('\d+', i)[-1])
        
        # self._th_cp = (Thread(target=external_plots.clonePlot, args=(fname, self._file_path.text(), name_id)))
        external_plots.clonePlot(fname, self._file_path.text(), name_id)
        # self._th_cp.start()
        # self.status.setText("mutation wave ongoing")
    
    def threadTabUI(self):
        self.threadTab = qtWidget.QWidget()
        _threadTab = qtWidget.QFormLayout()
        
        _refresh = qtWidget.QPushButton("Refresh")
        _refresh.clicked.connect(self.update)
        
        _threadTab.addRow("ID", _refresh)
        _threadTab.addRow("Started", qtWidget.QLabel())
        for i in self.s_ID:
            _threadTab.addRow(qtWidget.QLabel("ID: " + str(i)), qtWidget.QLabel())
        _threadTab.addRow(qtWidget.QLabel("Resumed"), qtWidget.QLabel())
        for i in self.r_ID:
            _threadTab.addRow(qtWidget.QLabel("ID: " + str(i)), qtWidget.QLabel()) 
        if self._th_muller:
            _threadTab.addRow(qtWidget.QLabel("muller plot"), qtWidget.QLabel()) 
        if self._th_hist:
            _threadTab.addRow(qtWidget.QLabel("histograms plot"), qtWidget.QLabel()) 

        self.threadTab.setLayout(_threadTab)
        return self.threadTab        
        
    def generalTabUI(self):
        """Create the General page UI."""
        generalTab = qtWidget.QWidget()
        _generalTab = qtWidget.QGridLayout()
        
        self._file_name = qtWidget.QLineEdit()
        _file_name_label = qtWidget.QLabel()
        _file_name_label.setText('File name')
        self._file_path = qtWidget.QLineEdit()
        _file_path_label = qtWidget.QLabel()
        _file_path_label.setText('File path')
        _file_path_button = qtWidget.QPushButton(self)
        _file_path_button.setText('Select path')
        _file_path_button.clicked.connect(self.selectPath)
        self._file_desc = qtWidget.QLineEdit()
        _file_desc_label = qtWidget.QLabel("TODO")
        _file_desc_label.setText('Description')
        
        _muller_plot_button = qtWidget.QPushButton(self)
        _muller_plot_button.setText('Clone diagram - muller plot')
        _muller_plot_button.clicked.connect(self.mullerPlotAction)
        _clone_hist_button = qtWidget.QPushButton(self)
        _clone_hist_button.setText('Mutations histograms')
        _clone_hist_button.clicked.connect(self.cloneHistAction)
        _mut_wave = qtWidget.QPushButton(self)
        _mut_wave.setText('Mutation wave plot')
        _mut_wave.clicked.connect(self.mutWaveAction)
        _fit_wave = qtWidget.QPushButton(self)
        _fit_wave.setText('Fitness wave plot')
        _fit_wave.clicked.connect(self.fitWaveAction)
        _clone_plot = qtWidget.QPushButton(self)
        _clone_plot.setText('Clone plot')
        _clone_plot.clicked.connect(self.clonePlotAction)
        
        _file_instr_fn = qtWidget.QLabel()
        _file_instr_fn.setText('File name - one cycle file name will appear with cycle number')
        _file_instr_fp = qtWidget.QLabel()
        _file_instr_fp.setText('File path - archive directory')
        _file_instr_fd = qtWidget.QLabel()
        _file_instr_fd.setText('Description - most important information to add to file')
        
        row = 0
        _generalTab.addWidget(_file_name_label, row, 0)
        _generalTab.addWidget(self._file_name, row, 1, 1, 2)
        row = row + 1
        _generalTab.addWidget(_file_path_label, row, 0)
        _generalTab.addWidget(self._file_path, row, 1)
        _generalTab.addWidget(_file_path_button, row, 2)
        row = row + 1
        _generalTab.addWidget(_file_desc_label, row, 0, 1, 3)
        row = row + 1
        _generalTab.addWidget(self._file_desc, row, 0, 1, 3)
        row = row + 1
        _generalTab.addWidget(_muller_plot_button, row, 0, 1, 3)
        row = row + 1
        _generalTab.addWidget(_clone_hist_button, row, 0, 1, 3)
        row = row + 1
        _generalTab.addWidget(_mut_wave, row, 0, 1, 3)
        row = row + 1
        _generalTab.addWidget(_fit_wave, row, 0, 1, 3)
        row = row + 1
        _generalTab.addWidget(_clone_plot, row, 0, 1, 3)
        row = row + 1
        _generalTab.addWidget(_file_instr_fd, row, 0, 1, 3)
        row = row + 1
        _generalTab.addWidget(_file_instr_fn, row, 0, 1, 3)
        row = row + 1
        _generalTab.addWidget(_file_instr_fp, row, 0, 1, 3)
        row = row + 1
        
        generalTab.setLayout(_generalTab)
        return generalTab
    
    def generatedTabUI(self):
        generatedTab = qtWidget.QWidget()
        _generatedTab = qtWidget.QGridLayout()
        
        _file_instr = qtWidget.QLabel()
        _file_instr.setText('Select figures to generate (and save optionally in file directory)')
        self._mw = qtWidget.QCheckBox("Mutation Wave")
        self._fw = qtWidget.QCheckBox("Fitness Wave")
        self._cp = qtWidget.QCheckBox("Clone plot")
        self._sf = qtWidget.QCheckBox("Save figures")
        self._sx = qtWidget.QCheckBox("Save files")
        
        _generatedTab.addWidget(_file_instr)
        _generatedTab.addWidget(self._mw)
        _generatedTab.addWidget(self._fw)
        _generatedTab.addWidget(self._cp)
        _generatedTab.addWidget(self._sf)
        _generatedTab.addWidget(self._sx)
        
        generatedTab.setLayout(_generatedTab)
        return generatedTab
    
    def parametersTabUI(self):
        parametersTab = qtWidget.QWidget()
        _parametersTab = qtWidget.QGridLayout()
        
        self._population = qtWidget.QLineEdit(str(10**6))
        self._population.setValidator(qtGui.QIntValidator(0, 10**7))
        
        self._capacity = qtWidget.QLineEdit(str(10**6))
        self._capacity.setValidator(qtGui.QIntValidator(0, 10**7))
        
        self._steps = qtWidget.QLineEdit(str(100))
        self._steps.setValidator(qtGui.QIntValidator(0, 10**5))
        
        self._tau = qtWidget.QLineEdit(str(0.005))     
        self._tau.setValidator(qtGui.QDoubleValidator(0.000, 1.000, 3))
        
        self._skip = qtWidget.QLineEdit(str(10))
        self._skip.setValidator(qtGui.QDoubleValidator(0.1, 100.0, 1))
        
        self._mut_effect = qtWidget.QLineEdit(".1, -.001")
        self._mut_effect.setInputMask("[#.#####, #.#####]")
        self._mut_prob = qtWidget.QLineEdit("0.0001, 0.05")        
        self._mut_prob.setInputMask("[#.#####, #.#####]")
        
        self._threads = qtWidget.QLineEdit(str(16))
        self._threads.setValidator(qtGui.QIntValidator(1, 2**10))
        
        _load_params = qtWidget.QPushButton(self)
        _load_params.setText('Load parameters')
        _load_params.clicked.connect(self.loadParams)
        
        _resume_button = qtWidget.QPushButton(self)
        _resume_button.setText('Resume simulation - Load File')
        _resume_button.clicked.connect(self.resume)
        
        _save_params = qtWidget.QPushButton(self)
        _save_params.setText('Save Parameters')
        _save_params.clicked.connect(self.saveParams)
        
        _calculate_crit = qtWidget.QPushButton(self)
        _calculate_crit.setText('Calculate critical population')
        _calculate_crit.clicked.connect(self.critPop)
        
        _file_instr_1 = qtWidget.QLabel()
        _file_instr_1.setText('Simulation time = steps * cycle time')
        _file_instr_2 = qtWidget.QLabel()
        _file_instr_2.setText('Mutation effect and probability are consistent lists')
        
        _parametersTab.addWidget(qtWidget.QLabel("Initial population"), 0, 0)
        _parametersTab.addWidget(self._population, 0, 1)
        _parametersTab.addWidget(qtWidget.QLabel("min: 0, max: 10^7"), 0, 2)
        
        _parametersTab.addWidget(qtWidget.QLabel("Enviroment capacity"), 1, 0)
        _parametersTab.addWidget(self._capacity, 1, 1)
        _parametersTab.addWidget(qtWidget.QLabel("min: 0, max: 10^7"), 1, 2)
        
        _parametersTab.addWidget(qtWidget.QLabel("Simulation steps"), 2, 0)
        _parametersTab.addWidget(self._steps, 2, 1)
        _parametersTab.addWidget(qtWidget.QLabel("min: 0, max: 10^5"), 2, 2)
        
        _parametersTab.addWidget(qtWidget.QLabel("Tau step"), 3, 0)
        _parametersTab.addWidget(self._tau, 3, 1)
        _parametersTab.addWidget(qtWidget.QLabel("min: 0, max: 1, step: 0.001"), 3, 2)
        
        _parametersTab.addWidget(qtWidget.QLabel("Cycle time"), 4, 0)
        _parametersTab.addWidget(self._skip, 4, 1)
        _parametersTab.addWidget(qtWidget.QLabel("min: 0.1, max: 100.0, step: 0.1"), 4, 2)
        
        _parametersTab.addWidget(qtWidget.QLabel("Mutation effect (list)"), 5, 0)
        _parametersTab.addWidget(self._mut_effect, 5, 1)
        _parametersTab.addWidget(qtWidget.QLabel("min: 0.00001, max: 0.99999"), 5, 2)
        
        _parametersTab.addWidget(qtWidget.QLabel("Mutation probability (list)"), 6, 0)
        _parametersTab.addWidget(self._mut_prob, 6, 1)
        _parametersTab.addWidget(qtWidget.QLabel("min: 0.00001, max: 0.99999"), 6, 2)  
        
        _parametersTab.addWidget(qtWidget.QLabel("Threads used in simulation"), 7, 0)
        _parametersTab.addWidget(self._threads, 7, 1)
        _parametersTab.addWidget(qtWidget.QLabel("min: 1, max: 1024"), 7, 2) 
        
        _parametersTab.addWidget(_load_params, 8, 0, 1, 3)
        _parametersTab.addWidget(_resume_button, 9, 0, 1, 3)
        _parametersTab.addWidget(_save_params, 10, 0, 1, 3)
        _parametersTab.addWidget(_calculate_crit, 11, 0, 1, 3)
        
        parametersTab.setLayout(_parametersTab)
        return parametersTab
    
    def critPop(self):        
        try:
            mut_prob = self._mut_prob.text().split(',')
            mut_prob = [float(x.strip('[]')) for x in mut_prob]
            mut_effect = self._mut_effect.text().split(',')
            mut_effect = [float(x.strip('[]')) for x in mut_effect]
        except:
            self.showDialog("Type correct parameters","Alert")
            return        
        
        self.status.setText("Critical Population: " + str(round((mut_prob[1]/mut_prob[0])*(mut_effect[1]/mut_effect[0]**2),0)))
    
    def selectPath(self):
        dir_path = qtWidget.QFileDialog.getExistingDirectory(self,"Choose Directory","Z:/")
        self._file_path.setText(dir_path)

    def saveParams(self):
        try:
            pop = int(self._population.text())
            cap = int(self._capacity.text())
            steps = int(self._steps.text())
            tau = float(self._tau.text())
            skip = float(self._skip.text())
            mut_prob = self._mut_prob.text().split(',')
            mut_prob = [float(x) for x in mut_prob]
            mut_effect = self._mut_effect.text().split(',')
            mut_effect = [float(x) for x in mut_effect]
            threads = int(self._threads.text())
        except:
            self.showDialog("Type correct parameters","Alert")
            return
        
        try:
            if self._file_name.text() == "":
                raise Exception()
            if self._file_path.text() == "":
                raise Exception()
        except:
            self.showDialog("Enter save localization and name", "Alert")
            return
            
        dfp = pd.DataFrame({
                'pop': pop,
                'cap': cap,
                'steps': steps,
                'tau': tau,
                'skip': skip,
                'mut_prob': [mut_prob],
                'mut_effect': [mut_effect],
                'threads': threads
            }, index=[0])
        
        filepath = Path(self._file_path.text() + '/params'  + ".csv")  
        filepath.parent.mkdir(parents=True, exist_ok=True)  
        dfp.to_csv(filepath)  

    def simStart(self):
        try:
            pop = int(self._population.text())
            cap = int(self._capacity.text())
            steps = int(self._steps.text())
            tau = float(self._tau.text())
            skip = float(self._skip.text())
            mut_prob = self._mut_prob.text().split(',')
            mut_prob = [float(x) for x in mut_prob]
            mut_effect = self._mut_effect.text().split(',')
            mut_effect = [float(x) for x in mut_effect]
            threads = int(self._threads.text())
        except:
            self.showDialog("Type correct parameters","Alert")
            return
        
        try:
            if self._file_name.text() == "" and self._sx.isChecked():
                raise Exception()
            if self._file_path.text() == "" and self._sx.isChecked():
                raise Exception()
        except:
            self.showDialog("Enter save localization and name", "Alert")
            return
        
        plots = 1*self._mw.isChecked() + \
                2*self._fw.isChecked() + \
                4*self._cp.isChecked() + \
                8*self._sf.isChecked() + \
                16*self._sx.isChecked()
        
        if self._binned.isChecked():
            iPop = [[0, pop, 1, 0, [], [], 0]]
            
            self.status.setText("Started")
            self.status.setStyleSheet("background-color: green")
            time.sleep(1)

            self.th_s.append(Process(target=clonalEvolutionMainLoop, args=(iPop, 
                  copy.deepcopy([pop, cap, steps, tau, skip, mut_prob, mut_effect, threads]), 
                  self._file_name.text(), 
                  self._file_desc.text(), 
                  self._file_path.text(), plots, 0, self.q, self.ID, 1)))
            self.s_ID.append(str(self.ID) + ", binned")
            
        elif self._single.isChecked():        
            iMuts = np.zeros(pop, dtype=np.int64).tolist()
            iProp = np.ones(pop).tolist()
            iClones = np.zeros(pop, dtype=np.int64).tolist()
            iMutations = [[] for x in range(len(iMuts))]
            iEffect =[[] for x in range(len(iMuts))]
            
            self.status.setText("Started")
            self.status.setStyleSheet("background-color: green")
            time.sleep(1)
    
            self.th_s.append(Process(target=clonalEvolutionMainLoop, args=(np.array([copy.deepcopy(iMuts),
                                                        copy.deepcopy(iProp), 
                                                        copy.deepcopy(iClones), 
                                                        copy.deepcopy(iMutations), 
                                                        copy.deepcopy(iEffect),
                                                        copy.deepcopy(iClones)]).T.tolist(), 
                                                        copy.deepcopy([pop, cap, steps, tau, skip, mut_prob, mut_effect, threads]), 
                                                        self._file_name.text(), 
                                                        self._file_desc.text(), 
                                                        self._file_path.text(), plots, 0, self.q, self.ID, 0)))               
            self.s_ID.append(str(self.ID) + ", single")
        elif self._matrix.isChecked():
            iPop = [[0, pop, [], [0], sc.sparse.csr_matrix(np.array([[0] for x in range(pop)])), 1, 0]]
            
            self.status.setText("Started")
            self.status.setStyleSheet("background-color: green")
            time.sleep(1)
            
            self.th_s.append(Process(target=clonalEvolutionMainLoop, args=(
                iPop, 
                copy.deepcopy([pop, cap, steps, tau, skip, mut_prob, mut_effect, threads]),
                self._file_name.text(), 
                self._file_desc.text(), 
                self._file_path.text(), plots, 0, self.q, self.ID, 2)))
            
            self.s_ID.append(str(self.ID) + ", matrix")
            
        self.ID = self.ID + 1
        self.th_s[self.idx_s].start()
        self.idx_s = self.idx_s + 1
       
    def resume(self):
        self.showPromptParams()
        fname = qtWidget.QFileDialog.getOpenFileName(self, 'Open file', "Z://","Single data files (*.csv);; Binned data files (*.txt);; Matrix data files (*.mtx)")[0]       
        if fname == "":
            self.showDialog("No file selected!", "Alert")
            return
        self.loadPaths(fname)
        
        try:
            pop = int(self._population.text())
            cap = int(self._capacity.text())
            steps = int(self._steps.text())
            tau = float(self._tau.text())
            skip = float(self._skip.text())
            mut_prob = self._mut_prob.text().split(',')
            mut_prob = [float(x.strip('[]')) for x in mut_prob]
            mut_effect = self._mut_effect.text().split(',')
            mut_effect = [float(x.strip('[]')) for x in mut_effect]
            threads = int(self._threads.text())
        except:
            self.showDialog("Type correct parameters","Alert")
            return            
        
        plots = 1*self._mw.isChecked() + \
                2*self._fw.isChecked() + \
                4*self._cp.isChecked() + \
                8*self._sf.isChecked() + \
                16*self._sx.isChecked()
        
        if self._binned.isChecked():
            df = external_plots.loadFile(fname)
            iClone = df['Clone number'].tolist()
            iCells = df['Cells number'].tolist()
            iFit = df['Mean fitness'].tolist()
            iMut = df['Mean mutation number'].tolist()
            iDriv = df['Driver mutation list'].tolist()
            iPass = df['Passener mutation list'].tolist()
            iParent = df['Previous clone number'].tolist()
            
            del df
            
            self.status.setText("Started")
            self.status.setStyleSheet("background-color: green")
            time.sleep(1)            
            self.th_r.append(Process(target=clonalEvolutionMainLoop, args=(np.array([copy.deepcopy(iClone),
                                                                  copy.deepcopy(iCells),
                                                                  copy.deepcopy(iFit),
                                                                  copy.deepcopy(iMut),
                                                                  copy.deepcopy(iDriv),
                                                                  copy.deepcopy(iPass),
                                                                  copy.deepcopy(iParent)]).T.tolist(), 
                                                        copy.deepcopy([pop, cap, steps, tau, skip, mut_prob, mut_effect, threads]), 
                                                        self._file_name.text(), 
                                                        self._file_desc.text(), 
                                                        self._file_path.text(), plots, self._last_cycle + 1, self.q, self.ID, 1)))
            self.r_ID.append(str(self.ID) + ", binned")
        elif self._single.isChecked():
            df = pd.read_csv(fname)
            iMuts = df['Mutations'].tolist()
            iClones = df['Clone'].tolist()
            iProp = df['Fitness'].tolist()
            iMutations = df['Mutations_ID'].tolist()
            iEffect = df['Mutations_effect'].tolist()
            iParent = df['Parent clone'].tolist()
            
            del df
            
            iMutations = [x.split(',') for x in iMutations]
            for x in range(len(iMutations)):
                for i in range(len(iMutations[x])):
                    iMutations[x][i] = float(iMutations[x][i].strip('[]'))
                
            iEffect = [x.split(',') for x in iEffect]
            for x in range(len(iEffect)):
                for i in range(len(iEffect[x])):
                    iEffect[x][i] = float(iEffect[x][i].strip('[]'))
                    
            self.status.setText("Started")
            self.status.setStyleSheet("background-color: green")
            time.sleep(1)
            self.th_r.append(Process(target=clonalEvolutionMainLoop, args=(np.array([copy.deepcopy(iMuts),
                                                        copy.deepcopy(iProp), 
                                                        copy.deepcopy(iClones), 
                                                        copy.deepcopy(iMutations), 
                                                        copy.deepcopy(iEffect),
                                                        copy.deepcopy(iParent)]).T.tolist(), 
                                                       copy.deepcopy(self.params).append(threads), 
                                                       self._file_name.text(), 
                                                       self._file_desc.text(), 
                                                       self._file_path.text(), plots, self._last_cycle + 1, self.q, self.ID, 0)))            
            self.r_ID.append(str(self.ID) + ", single")
            
        elif self._matrix.isChecked():
            df = pd.read_csv(fname)
            df = df.drop('Unnamed: 0', axis=1)
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
                    
            iPop = df.to_numpy().tolist()
            
            self.status.setText("Started")
            self.status.setStyleSheet("background-color: green")
            time.sleep(1)
            
            self.th_r.append(Process(target=clonalEvolutionMainLoop, args=(iPop, copy.deepcopy([pop, cap, steps, tau, skip, mut_prob, mut_effect, threads]),
                                    self._file_name.text(), 
                                    self._file_desc.text(), 
                                    self._file_path.text(), plots, self._last_cycle + 1, self.q, self.ID, 2)))
            self.r_ID.append(str(self.ID) + ", matrix")
            
        self.ID = self.ID + 1
        self.th_r[self.idx_r].start()
        self.idx_r = self.idx_r + 1     
        
    def showDialog(self, text, title):
        msgBox = qtWidget.QMessageBox()
        msgBox.setIcon(qtWidget.QMessageBox.Information)
        msgBox.setText(text)
        msgBox.setWindowTitle(title)
        msgBox.setStandardButtons(qtWidget.QMessageBox.Ok)
         
        returnValue = msgBox.exec()
        if returnValue == qtWidget.QMessageBox.Ok:
           print('OK clicked')
           
    def showPromptParams(self):
        msgBox = qtWidget.QMessageBox()
        msgBox.setIcon(qtWidget.QMessageBox.Information)
        msgBox.setText("Have you change simulation parameters as continiuing simulation?")
        msgBox.setWindowTitle("Parameters")
        msgBox.setStandardButtons(qtWidget.QMessageBox.Yes | qtWidget.QMessageBox.No)       
        returnValue = msgBox.exec()
        if returnValue == qtWidget.QMessageBox.No:
           self.loadParams()
           
    def loadParams(self):
        fname = qtWidget.QFileDialog.getOpenFileName(self, 'Load simulation parameters', "Z://","CSV files (*.csv)")[0]        
        if fname == "":
            self.showDialog("No file selected!", "Alert")
            return
        df = pd.read_csv(fname)
        
        pop = df['pop'][0]
        cap = df['cap'][0]
        steps = df['steps'][0]
        tau = df['tau'][0]
        skip = df['skip'][0]
        xx = df['mut_prob'][0]
        xx = xx.split(',')
        xx = [float(x.strip('[]')) for x in xx]
        mut_prob = ""
        for i in xx:
            mut_prob = mut_prob + str(i) + ','
        mut_prob = mut_prob.rstrip(',')
        xx = df['mut_effect'][0]       
        xx = xx.split(',')
        xx = [float(x.strip('[]')) for x in xx]
        mut_effect = ""
        for i in xx:
            mut_effect = mut_effect + str(i) + ','
        mut_effect = mut_effect.rstrip(',')
        threads = df['threads'][0]
        
        # self._population.setReadOnly(True)
        # self._capacity.setReadOnly(True)
        # self._steps.setReadOnly(True)
        # self._tau.setReadOnly(True)
        # self._skip.setReadOnly(True)
        self._mut_prob.setInputMask("")
        # self._mut_prob.setReadOnly(True)
        self._mut_effect.setInputMask("")
        # self._mut_effect.setReadOnly(True)
        
        self.params = [pop, cap, steps, tau, skip, mut_prob, mut_effect]
        self._population.setText(str(pop))
        self._capacity.setText(str(cap))
        self._steps.setText(str(steps))
        self._tau.setText(str(tau))
        self._skip.setText(str(skip))
        self._mut_prob.setText(str(mut_prob))
        self._mut_effect.setText(str(mut_effect))
        self._threads.setText(str(threads))
        
    def loadPaths(self, fname):        
        # self.showDialog("Select file same as resume simulation", "Info")
        xname = fname.split('/')
        
        ##TODO interpret only last number!!!
        
        path = ""
        name = ""
        self._last_cycle = 0
        for i in xname:
            if xname.index(i) == len(xname)-1:
                t = i.split('_')
                t = t[len(t)-1]
                self._last_cycle = int(''.join(x for x in t if x.isdigit()))
                if i.endswith('.csv'):
                    name = i.rstrip('_' + str(self._last_cycle) + '.csv')
                elif i.endswith('.txt'):
                    name = i.rstrip('_' + str(self._last_cycle) + '.txt')
            else:
                path = path + i + '/'
        
        if name.endswith("_binned"):
            name = name.rstrip("_binned")
        else:
            name = name.rstrip("_single")
        self._file_name.setText(name)
        self._file_path.setText(path)  
        self._file_desc.setText(str(self._last_cycle) + " - last cycle, starting with +1")                    

def run():
    app = qtWidget.QApplication(sys.argv)
    win = mainFormat()
    win.show()
    ret = app.exec_()
    sys.exit(ret)
    
if __name__ == "__main__":
    run()