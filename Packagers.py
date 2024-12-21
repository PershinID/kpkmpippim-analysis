import os

import ROOT
from ROOT import TFile

class Packager:
    def __init__(self, package_path, list_energies, data_folder_path):
        self.dataPathTemplate = "hists_e{energy}.root"
        self.dataFolder = data_folder_path
        self.energies = list_energies
        self.packageFile = TFile.Open(package_path, 'update')
        self.isStructured = False

    def createStructure(self):
        if os.listdir(self.dataFolder) == []: raise OSError("Histfiles not found")

        ##Creating directories structure
        file_for_struct = TFile.Open(self.dataFolder + os.listdir(self.dataFolder)[0], "read")
        self.hists_to_merge = {
            dir_name: [
                '/'.join(dir_name, hist.GetName())
                for hist in file_for_struct.GetDirectory(dir_name).GetListOfKeys()
                if hist.GetClassName() in ("TH1I", "TH1F", "TH2F")
            ]
            for dir_name in [
                directory.GetName()
                for directory in file_for_struct.GetListOfKeys()
                if directory.GetClassName() == "TDirectoryFile"
            ]
        }
        for dir_name in self.hists_to_merge.keys():
            self.packageFile.mkdir(dir_name)
        self.isStructured = True
        file_for_struct.Close()


    def mergeHists(self):
        pass


    def joinHistsIntoDir(self):
        if not self.isStructured: self.createStructure()

        for dir_name in self.hists_to_merge.keys():
            hist_dir_names = self.hists_to_merge[dir_name]
            for hist_dir_name in hist_dir_names.keys():
                self.packageFile.mkdir(hist_dir_name)

        for energy in self.energies:
            hists_file = TFile.Open(self.dataFolder + self.dataPathTemplate.format(energy = energy), "read")
            for dir_name in self.hists_to_merge.keys():
                for hist_dir_name in self.hists_to_merge[dir_name].keys():
                    self.packageFile.GetDirectory(hist_dir_name).cd()
                    histogram = hists_file.Get(hist_dir_name)
                    histogram.Write()
            hists_file.Close()
