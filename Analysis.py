from enum import Enum
from copy import copy

from math import pi, ceil, floor
from functools import reduce

import logging
from logging import getLogger, StreamHandler, FileHandler, Formatter

import sys

import ROOT
from ROOT import TFile
from ROOT import TH1I, TH1F, TH2F, TGraph, TCanvas
from ROOT import gROOT

from Container import Container

EHists = Enum('EHists', ['TH1I', 'TH1F', 'TH2F'])

class CutDispatcher:
    def __init__(self, cuts_available, *, n_entries_full = 0):
        self.CutsAvailable = cuts_available

        ## Current progress
        self.CutsCurrent = []
        self.Checklist = self.createChecklist(n_entries_full)

        ## Progress done
        self.CutsCounter = 0
        self.CutsDone = {'no_cut': n_entries_full} ## {name: n_entries_selected}
        self.GraphEntriesSelected = None
        self.GraphEntriesPercentage = None
        self.GraphCutPercentage = None
        
        
    ## Operations with available cuts
    def isCutAvailable(self, cut_name):
        return cut_name in self.CutsAvailable

    ## Operations with current cuts
    def isCutCurrent(self, cut_name):
        return cut_name in self.CutsCurrent

    def addCut(self, cut_name, is_inversed = False):
        if not self.isCutAvailable(cut_name): raise ValueError(f"Cut blueprint not found: {cut_name}")
        if self.isCutCurrent(cut_name): return
        if self.isCutDone(cut_name): return
        self.CutsCurrent.append(cut_name)

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
            self.GraphCutPercentage.SetPoint(i, i, cut_item[1] / n_entries_prev)
            n_entries_prev = cut_item[1]


class HistogramDispatcher:
    def __init__(self, histograms_available):
        self.HistogramsAvailable = histograms_available
        self.HistogramsCurrent = []
        self.HistogramsBuilt = {}

    def isHistogramAvailable(self, hist_name):
        return hist_name in self.HistogramsAvailable
        
    def isHistogramCurrent(self, hist_name):
        return hist_name in self.HistogramsCurrent
    
    def addHistogram(self, hist_name):
        if not self.isHistogramAvailable(hist_name): raise ValueError(f"Histogram blueprint not found: {hist_name}")
        if self.isHistogramCurrent(hist_name): return
        self.HistogramsCurrent.append(hist_name)

    def buildHistogram(self, hist_name):
        if hist_name in self.HistogramsBuilt:
            self.HistogramsBuilt[hist_name][0].Reset()
            return self.HistogramsBuilt[hist_name]
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
            return self.HistogramsBuilt[hist_name]

        if hist_blueprint['type'] is EHists.TH1F:
            self.HistogramsBuilt[hist_name] = (
                TH1F(
                    hist_name,
                    ';'.join([hist_args['title'], hist_args['x-axis-title']]),
                    hist_args['x-axis-nbins'], *hist_args['x-axis-range']
                ),
                hist_args['x-variable'],
            )
            return self.HistogramsBuilt[hist_name]

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
            return self.HistogramsBuilt[hist_name]


class Analysis:
    def __init__(self, path, mode = 'recreate', *, logname = "analysis", logpath = None, filemode = 'w'):
        self.InputContainer = None
        self.OutputContainer = None
        self.isOutputDumped = False
        self.AnalysisFile = TFile.Open(path, mode)
        self.Variables = {}

        self.Logger = getLogger(logname)
        self.Logger.setLevel(logging.INFO)
        log_handler = StreamHandler(sys.stdout) if logpath == None else FileHandler(logpath, filemode)
        log_handler.setFormatter(
            Formatter(
                fmt = "%(asctime)s: %(name)s: %(message)s",
                datefmt = "%d.%m.%Y %H:%M:%S"
            )
        )
        self.Logger.addHandler(log_handler)
        
        self.CutDispatcher = None
        self.HistogramDispatcher = None
        

    ## Operations with cuts
    def addCut(self, cut_name, is_inversed = False):
        self.CutDispatcher.addCut(cut_name, is_inversed)

    def getCutsCurrent(self):
        return self.CutDispatcher.CutsCurrent

    ## Operations with histograms
    def addHistogram(self, hist_name):
        self.HistogramDispatcher.addHistogram(hist_name)

    def getHistogramsCurrent(self):
        return self.HistogramDispatcher.HistogramsCurrent

    ## Operations with entries
    def getEntry(self, n_entry):
        self.InputContainer.getEntry(n_entry)

    def fillEntry(self):
        self.OutputContainer.fillEntry()

    def calculateEntry(self):
        pass

    ## Processing loop
    def loop(self, *, directory_name = None):
        if directory_name is None:
            directory_name = str(self.CutDispatcher.CutsCounter) + '_' + '_'.join(self.CutDispatcher.CutsCurrent)
        if directory_name == '':
            raise ValueError('Directory not chosen')

        hists = {}
        for hist_name in self.HistogramDispatcher.HistogramsCurrent:
            hists.update({
                hist_name: self.HistogramDispatcher.buildHistogram(hist_name)
            })
        
        cut_set_name = '_'.join(self.CutDispatcher.CutsCurrent)
        self.Logger.info(f"Starting '{cut_set_name}' cut")
        n_entries_prev = self.CutDispatcher.getEntriesSelected()
        for n_entry in range(self.InputContainer.getEntries()):
            if not self.CutDispatcher.checkEntryInChecklist(n_entry): continue
            self.getEntry(n_entry)

            ## Executing cuts
            for cut_name in self.CutDispatcher.CutsCurrent:
                cut_func = self.CutDispatcher.CutsAvailable[cut_name]
                is_inversed = True if cut_name[0] == '!' else False
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
        self.CutDispatcher.CutsDone.update({cut_set_name: self.CutDispatcher.getEntriesSelected()})
        self.CutDispatcher.CutsCurrent.clear()
        self.HistogramDispatcher.HistogramsCurrent.clear()
        self.CutDispatcher.CutsCounter += 1

        ## Saving histograms to directory
        self.AnalysisFile.mkdir(directory_name)
        self.AnalysisFile.GetDirectory(directory_name).cd()
        for hist in hists.values(): hist[0].Write()
        self.AnalysisFile.Save()
        self.AnalysisFile.cd()
        
        
    def dumpToFile(self):
        if self.isOutputDumped: return

        self.AnalysisFile.cd()
        self.CutDispatcher.buildGraphs()
        self.CutDispatcher.GraphEntriesSelected.Write()
        self.CutDispatcher.GraphEntriesPercentage.Write()
        self.CutDispatcher.GraphCutPercentage.Write()
        
        self.AnalysisFile.Save()
        
        for n_entry in range(self.InputContainer.getEntries()):
            if not self.CutDispatcher.checkEntryInChecklist(n_entry): continue

            self.getEntry(n_entry)
            self.calculateEntry()
            self.fillEntry()
        self.OutputContainer.dumpToFile()
        self.isOutputDumped = True


if __name__ == '__main__':
    cut_dispatcher = CutDispatcher('test', {}, n_entries_full = 128)
    print(cut_dispatcher.Checklist)
    print(cut_dispatcher.getEntriesSelected())
    cut_dispatcher.deleteEntryFromChecklist(5)
    cut_dispatcher.deleteEntryFromChecklist(10)
    cut_dispatcher.deleteEntryFromChecklist(15)
    print(cut_dispatcher.Checklist)
    print(cut_dispatcher.getEntriesSelected())
