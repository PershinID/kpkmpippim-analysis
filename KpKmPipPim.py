from argparse import ArgumentParser
from enum import Enum

import ROOT
from ROOT import TFile, TTree

from Container import Container
from CMD3ContainerV9 import CMD3ContainerV9
from PreliminaryContainer import PreliminaryContainer

from PreliminaryAnalysis import PreliminaryAnalysis


if __name__ == '__main__':
    ##CMD3 --> Preliminary
    ##Path to datafiles
    path_cmd3data = "root://cmd//scan2020/scan2020_tr_ph_fc_e950_v9.root"
    path_prelim = "/store11/idpershin/prelim/scan2020_tr_ph_fc_e950_v9.root"

    ##Path to histograms
    path_hists = "/store11/idpershin/prelim/hists2020_tr_ph_fc_e950_v9.root"

    prelim_analysis = PreliminaryAnalysis(path_hists, path_cmd3data, path_prelim)
    prelim_analysis.cut('nt')
    prelim_analysis.cut('tcharge')
    prelim_analysis.createHistogram('h_tptot', 'tcharge')
    prelim_analysis.createHistogram('h_tnhit', 'tcharge')
    prelim_analysis.createHistogram('h_tth', 'tcharge')
    prelim_analysis.createHistogram('h_tphi', 'tcharge')
    prelim_analysis.createHistogram('h_tptot_tdedx', 'tcharge')
    prelim_analysis.cut('tptot')
    prelim_analysis.createHistogram('h_tnhit', 'tptot')
    prelim_analysis.createHistogram('h_tth', 'tptot')
    prelim_analysis.createHistogram('h_tphi', 'tptot')
    prelim_analysis.createHistogram('h_tptot_tdedx', 'tptot')
    prelim_analysis.cut('tnhit')
    prelim_analysis.createHistogram('h_tth', 'tnhit')
    prelim_analysis.createHistogram('h_tphi', 'tnhit')
    prelim_analysis.createHistogram('h_tptot_tdedx', 'tnhit')
    prelim_analysis.cut('tth')
    prelim_analysis.createHistogram('h_tphi', 'tth')
    prelim_analysis.createHistogram('h_tptot_tdedx', 'tth')
    prelim_analysis.createHistogram('h_TotalP_DeltaE', 'tth')
    prelim_analysis.cut('trho')
    prelim_analysis.cut('tz')
    prelim_analysis.cut('DeltaE')
    prelim_analysis.cut('TotalP')
    prelim_analysis.createHistogram('h_TotalP_DeltaE', 'TotalP')
    prelim_analysis.dumpToFile()

    
