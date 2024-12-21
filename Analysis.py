from enum import Enum
from math import pi, ceil, floor
from copy import copy

import ROOT
from ROOT import TFile
from ROOT import TH1I, TH1F, TH2F

from Container import Container

EHists = Enum('EHists', ['TH1I', 'TH1F', 'TH2F'])
    
class Analysis:
    def __init__(self, path):
        self.InputContainer = None
        self.OutputContainer = None
        self.AnalysisFile = TFile.Open(path, 'recreate')
        self.EntriesChecklist = None
        self.Cuts = {}
        self.Histograms = {}
        self.Variables = {}

    def initializeEntriesChecklist(self):
        self.EntriesChecklist = bytearray(b'\xff' * ceil(self.InputContainer.getEntries() / 8))
        
    def cut(self, cut_name, is_inversed = False):
        if cut_name not in self.Cuts: raise ValueError("Cut not found")
        cut_func = self.Cuts[cut_name]

        for n_entry in range(self.InputContainer.getEntries()):
            if not self._checkEntryInChecklist(n_entry): continue

            self.getEntry(n_entry)
            if not is_inversed and cut_func():
                self._deleteEntryFromChecklist(n_entry)
            elif is_inversed and not cut_func():
                self._deleteEntryFromChecklist(n_entry)
            else: pass

        self.AnalysisFile.mkdir(cut_name)


    def _deleteEntryFromChecklist(self, n_entry):
        ##Finds byte from checklist and zero proper bit in this byte
        byte_index, bit_index = floor(n_entry / 8), 0b10000000 >> (n_entry % 8)
        self.EntriesChecklist[byte_index] &= (0b11111111 - bit_index)

    def _checkEntryInChecklist(self, n_entry):
        byte_index, bit_index = floor(n_entry / 8), 0b10000000 >> (n_entry % 8)
        return (self.EntriesChecklist[byte_index] & bit_index) != 0

    def getEntriesSelected(self):
        event_count = 0
        for byte in self.EntriesChecklist:
            event_count += bin(byte).count('1')
        return event_count
    
    def dumpToFile(self):
        for n_entry in range(self.InputContainer.getEntries()):
            if not self._checkEntryInChecklist(n_entry): continue

            self.getEntry(n_entry)
            self.fillEntry()
        self.OutputContainer.dumpToFile()


    def getEntry(self, n_entry):
        self.InputContainer.getEntry(n_entry)

    def fillEntry(self):
        self.OutputContainer.fillEntry()
        
    def createHistogram(self, hist_name, directory_name):
        if hist_name not in self.Histograms.keys(): raise ValueError("Histogram blueprint not found")

        hist_dict = self.Histograms[hist_name]
        hist_type = hist_dict['type']
        hist_args = hist_dict['args']
        hist_opts = hist_dict.get('opts', None)
        if hist_type is EHists.TH1I:
            hist = TH1I(
                hist_name,
                ';'.join([hist_args['title'], hist_args['x-axis-title']]),
                hist_args['x-axis-nbins'], *hist_args['x-axis-range'],
            )
            hist_vars = (hist_args['x-variable'],)

        elif hist_type is EHists.TH1F:
            hist = TH1F(
                hist_name,
                ';'.join([hist_args['title'], hist_args['x-axis-title']]),
                hist_args['x-axis-nbins'], *hist_args['x-axis-range'],
            )
            hist_vars = (hist_args['x-variable'],)

        elif hist_type is EHists.TH2F:
            hist = TH2F(
                hist_name,
                ';'.join([hist_args['title'], hist_args['x-axis-title'], hist_args['y-axis-title']]),
                hist_args['x-axis-nbins'], *hist_args['x-axis-range'],
                hist_args['y-axis-nbins'], *hist_args['y-axis-range'],
            )
            hist_vars = (hist_args['x-variable'], hist_args['y-variable'],)

        ##Filling the histogram
        for n_entry in range(0, self.InputContainer.getEntries()):
            if not self._checkEntryInChecklist(n_entry): continue

            self.getEntry(n_entry)
            vars_array = list(zip(
                *map(lambda a: a.getContentList(), hist_vars)
            ))
            for variables in vars_array: hist.Fill(*variables)

        self.AnalysisFile.GetDirectory(directory_name).cd()
        hist.Write()
        self.AnalysisFile.Save()
