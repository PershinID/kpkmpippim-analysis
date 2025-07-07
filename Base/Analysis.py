from enum import Enum
from copy import copy
from collections import OrderedDict

from math import pi, ceil, floor
from functools import reduce

import logging
from logging import getLogger, StreamHandler, FileHandler, Formatter

import sys

import ROOT
from ROOT import TFile
from ROOT import TH1I, TH1F, TH2F, TGraph, TCanvas
from ROOT import gROOT

from .Container import Container

EHists = Enum('EHists', ['TH1I', 'TH1F', 'TH2F'])

class CutDispatcher:
    def __init__(self, cuts_available, *, n_entries_full = 0):
        self.CutsAvailable = cuts_available

        ## Current progress
        self.CutsCurrent = {}
        self.Checklist = self.createChecklist(n_entries_full)

        ## Progress done
        self.CutsCounter = 0
        self.CutsDone = OrderedDict()
        self.CutsDone.update({'no_cut': n_entries_full}) ## {name: n_entries_selected}
        self.GraphEntriesSelected = None
        self.GraphEntriesPercentage = None
        self.GraphCutPercentage = None
        
        
    ## Operations with available cuts
    def isCutAvailable(self, cut_name):
        return cut_name in self.CutsAvailable

    ## Operations with current cuts
    def isCutCurrent(self, cut_name):
        return cut_name in self.CutsCurrent

    def isCutInversed(self, cut_name):
        return cut_name[0] == '!'
    
    def addCut(self, cut_name, is_inversed = False):
        if not self.isCutAvailable(cut_name): raise ValueError(f"Cut blueprint not found: {cut_name}")
        if self.isCutCurrent(cut_name): return
        if self.isCutDone(cut_name): return
        cut_name_ = cut_name if not is_inversed else '!' + cut_name
        self.CutsCurrent.update({cut_name_: self.CutsAvailable[cut_name]})

    def addCutNew(self, cut_name, cut_func, is_inversed = False):
        if self.isCutAvailable(cut_name): raise ValueError(f"Cut blueprint already exists: {cut_name}")
        if self.isCutCurrent(cut_name): return
        if self.isCutDone(cut_name): return
        cut_name_ = cut_name if not is_inversed else '!' + cut_name
        self.CutsCurrent.update({cut_name_: cut_func})         

    def getCutsCurrent(self):
        return self.CutsCurrent
        
    ## Operations with checklist
    def createChecklist(self, n_entries_full):
        checklist = bytearray(b'\xff' * (n_entries_full // 8))
        if n_entries_full % 8 != 0:
            checklist += ( (1 << 8) - (1 << (8 - n_entries_full % 8)) ).to_bytes(1, 'big')
        return checklist
    
    def deleteEntryFromChecklist(self, n_entry):
        byte_index, bit_index = (n_entry // 8), 0b10000000 >> (n_entry % 8)
        self.Checklist[byte_index] &= (0b11111111 - bit_index)

    def checkEntryInChecklist(self, n_entry):
        byte_index, bit_index = (n_entry // 8), 0b10000000 >> (n_entry % 8)
        return (self.Checklist[byte_index] & bit_index) != 0

    def getEntriesSelected(self):
        event_count = 0
        for byte in self.Checklist:
            event_count += bin(byte).count('1')
        return event_count

    ## Operations with done cuts
    def isCutDone(self, cut_name):
        return cut_name in self.CutsDone

    def update(self):
        self.CutsDone.update({self.formFullCutName(): self.getEntriesSelected()})
        self.CutsCounter += 1
        self.CutsCurrent = {}
    
    def buildGraphs(self):
        n_entries_full = self.CutsDone['no_cut']
        n_entries_prev = n_entries_full
        
        n_cuts = len(self.CutsDone)
        self.GraphEntriesSelected = TGraph(n_cuts)
        self.GraphEntriesSelected.SetName('g_entries_selected')
        self.GraphEntriesSelected.SetTitle('Number of selected events')
        self.GraphEntriesPercentage = TGraph(n_cuts)
        self.GraphEntriesPercentage.SetName('g_entries_percentage')
        self.GraphEntriesPercentage.SetTitle('Percentage of selected events')
        self.GraphCutPercentage = TGraph(n_cuts)
        self.GraphCutPercentage.SetName('g_cut_percentage')
        self.GraphCutPercentage.SetTitle('Percentage of selected events by cut')
        ges_xaxis = self.GraphEntriesSelected.GetXaxis()
        ges_xaxis.SetTitle('Cut')
        gep_xaxis = self.GraphEntriesPercentage.GetXaxis()
        gep_xaxis.SetTitle('Cut')
        gcp_xaxis = self.GraphCutPercentage.GetXaxis()
        gcp_xaxis.SetTitle('Cut')
        for i, cut_item in enumerate(self.CutsDone.items()):
            self.GraphEntriesSelected.SetPoint(i, i, cut_item[1])
            self.GraphEntriesPercentage.SetPoint(i, i, cut_item[1] / n_entries_full)
            self.GraphCutPercentage.SetPoint(i, i, cut_item[1] / n_entries_prev if n_entries_prev != 0 else 0.)
            n_entries_prev = cut_item[1]


    def formDirectoryName(self):
        if not self.CutsCurrent: return str(self.CutsCounter)
        return str(self.CutsCounter) + '_' + '_'.join(self.CutsCurrent)

    def formFullCutName(self):
        if not self.CutsCurrent: return 'no_cut'
        return '_'.join(self.CutsCurrent)
        

class HistogramDispatcher:
    def __init__(self, histograms_available):
        self.HistogramsAvailable = histograms_available
        self.HistogramsCurrent = {}
        self.HistogramsBuilt = {}

    def isHistogramAvailable(self, hist_name):
        return hist_name in self.HistogramsAvailable
        
    def isHistogramCurrent(self, hist_name):
        return hist_name in self.HistogramsCurrent
    
    def addHistogram(self, hist_name):
        if not self.isHistogramAvailable(hist_name): raise ValueError(f"Histogram blueprint not found: {hist_name}")
        if self.isHistogramCurrent(hist_name): return
        self.buildHistogram(hist_name)
        
    def buildHistogram(self, hist_name):
        if hist_name in self.HistogramsBuilt:
            self.HistogramsBuilt[hist_name][0].Reset()
            self.HistogramsCurrent[hist_name] = self.HistogramsBuilt[hist_name]
            return
        
        hist_blueprint = self.HistogramsAvailable[hist_name]
        hist_args = hist_blueprint['args']
            
        if hist_blueprint['type'] is EHists.TH1I:
            self.HistogramsBuilt[hist_name] = (
                TH1I(
                    hist_name,
                    ';'.join([hist_args['title'], hist_args['x-axis-title']]),
                    hist_args['x-axis-nbins'], *hist_args['x-axis-range']
                ),
                hist_args['x-variable'],
            )

        if hist_blueprint['type'] is EHists.TH1F:
            self.HistogramsBuilt[hist_name] = (
                TH1F(
                    hist_name,
                    ';'.join([hist_args['title'], hist_args['x-axis-title']]),
                    hist_args['x-axis-nbins'], *hist_args['x-axis-range']
                ),
                hist_args['x-variable'],
            )

        if hist_blueprint['type'] is EHists.TH2F:
            self.HistogramsBuilt[hist_name] = (
                TH2F(
                    hist_name,
                    ';'.join([hist_args['title'], hist_args['x-axis-title'], hist_args['y-axis-title']]),
                    hist_args['x-axis-nbins'], *hist_args['x-axis-range'],
                    hist_args['y-axis-nbins'], *hist_args['y-axis-range']
                ),
                hist_args['x-variable'], hist_args['y-variable'],
            )

        self.HistogramsCurrent[hist_name] = self.HistogramsBuilt[hist_name]


    def getHistogramsCurrent(self):
        return self.HistogramsCurrent
        
    def getHistogramsCurrentList(self):
        return self.HistogramsCurrent.keys()

    def clearHistogramsCurrent(self):
        self.HistogramsCurrent = {}
    

class Analysis:
    def __init__(self, path = None, *, logname = "analysis", logpath = None, logmode = 'w'):
        self.InputContainers = []
        self.OutputContainers = []
        self.AnalysisFile = TFile.Open(path, 'recreate') if path else None
        self.Variables = {}

        self.Logger = getLogger(logname)
        self.Logger.setLevel(logging.INFO)
        self.LogHandler = StreamHandler(sys.stdout) if logpath == None else FileHandler(logpath, logmode)
        self.LogHandler.setFormatter(
            Formatter(
                fmt = "%(asctime)s: %(name)s: %(message)s",
                datefmt = "%d.%m.%Y %H:%M:%S"
            )
        )
        self.Logger.addHandler(self.LogHandler)
        
        self.CutDispatcher = None
        self.HistogramDispatcher = None
        

    ## Operations with cuts
    def addCut(self, cut_name, is_inversed = False):
        self.CutDispatcher.addCut(cut_name, is_inversed)

    def addCutNew(self, cut_name, cut_func, is_inversed = False):
        self.CutDispatcher.addCutNew(cut_name, cut_func, is_inversed)
        
    ## Operations with histograms
    def addHistogram(self, hist_name):
        if not self.AnalysisFile: return
        self.HistogramDispatcher.addHistogram(hist_name)

    def getHistogramsCurrent(self):
        if not self.AnalysisFile: return {}
        return self.HistogramDispatcher.getHistogramsCurrent()
        
    def getHistogramsCurrentList(self):
        if not self.AnalysisFile: return []
        return self.HistogramDispatcher.getHistogramsCurrentList()

    def saveHistograms(self, hists, directory_name):
        if not self.AnalysisFile: return
        if not self.AnalysisFile.GetDirectory(directory_name):
            self.AnalysisFile.mkdir(directory_name)
        self.AnalysisFile.GetDirectory(directory_name).cd()
        for hist in hists.values(): hist[0].Write()
        self.AnalysisFile.Save()
        self.AnalysisFile.cd()

    def clearHistogramsCurrent(self):
        if not self.AnalysisFile: return
        return self.HistogramDispatcher.clearHistogramsCurrent()
        
    ## Operations with entries
    def getEntry(self, n_entry):
        for container in self.InputContainers:
            container.getEntry(n_entry)

    def getEntries(self):
        return self.InputContainers[0].getEntries()
            
    def fillEntry(self):
        for container in self.OutputContainers:
            container.fillEntry()

    def calculateEntry(self):
        pass

    ## Processing loop
    def loop(self, *, directory_name = None):
        if not directory_name:
            directory_name = self.CutDispatcher.formDirectoryName() 

        ## hists: {'hist_name': (THist, Variable1, Variable2, ...), ...}
        hists = self.getHistogramsCurrent()
        
        cut_set_name = self.CutDispatcher.formFullCutName()
        self.Logger.info(f"Starting '{cut_set_name}' cut")
        n_entries_prev = self.CutDispatcher.getEntriesSelected()
        for n_entry in range(self.getEntries()):
            if not self.CutDispatcher.checkEntryInChecklist(n_entry): continue
            self.getEntry(n_entry)

            ## Executing cuts
            for cut_name, cut_func in self.CutDispatcher.getCutsCurrent().items():
                is_inversed = self.CutDispatcher.isCutInversed(cut_name)
                if not is_inversed and cut_func():
                    self.CutDispatcher.deleteEntryFromChecklist(n_entry)
                    break
                elif is_inversed and not cut_func():
                    self.CutDispatcher.deleteEntryFromChecklist(n_entry)
                    break
                else: pass
                
            else:
                ## Filling histograms
                for hist in hists.values():
                    ## Example: (MomentumVariable, ThetaVariable, ...)
                    ## -> ([momentum1, momentum2, ..., momentumN], [theta1, theta2, ..., thetaN], ...)
                    ## -> ((momentum1, theta1, ...), (momentum2, theta2, ...), ..., (momentumN, thetaN, ...))
                    variables_values_by_tuples = list(zip(
                        *map(lambda a: [a.Content] if a.Sizes == (1,) else a.Content, hist[1:])
                    ))
                    for values in variables_values_by_tuples: hist[0].Fill(*values)

        ## Logging cut results, clearing CutsCurrent and HistogramsCurrent
        self.Logger.info(f"'{cut_set_name}' cut finished. {self.CutDispatcher.getEntriesSelected()} entries out of {n_entries_prev} selected")
        self.CutDispatcher.update()
        self.clearHistogramsCurrent()

        ## Saving histograms to directory
        self.saveHistograms(hists, directory_name)
        
        
    def dumpToFile(self):
        if self.AnalysisFile:
            self.AnalysisFile.cd()
            self.CutDispatcher.buildGraphs()
            self.CutDispatcher.GraphEntriesSelected.Write()
            self.CutDispatcher.GraphEntriesPercentage.Write()
            self.CutDispatcher.GraphCutPercentage.Write()
            self.AnalysisFile.Save()

        if self.OutputContainers:
            for n_entry in range(self.getEntries()):
                if not self.CutDispatcher.checkEntryInChecklist(n_entry): continue
                
                self.getEntry(n_entry)
                self.calculateEntry()
                self.fillEntry()

            for container in self.OutputContainers:
                container.dumpToFile()

                
    def close(self):
        if self.AnalysisFile:
            self.AnalysisFile.Close()
            del self.HistogramDispatcher
        for container in self.InputContainers: container.close()
        for container in self.OutputContainers: container.close()
        del self.CutDispatcher
        self.Logger.removeHandler(self.LogHandler)
        self.LogHandler.close()
        del self.Logger


if __name__ == '__main__':
    analysis = Analysis('test.root')
    analysis.close()

    analysis = Analysis()
    analysis.close()
