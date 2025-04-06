#!/usr/bin/python3.6
import os.path

from argparse import ArgumentParser
from datetime import date

import ROOT
from ROOT import TFile, TTree

from PreliminaryAnalysis import PreliminaryAnalysis

if __name__ == '__main__':
    ##Parsing input arguments
    parser = ArgumentParser()
    parser.add_argument('--version', default = 'v9', choices = ['v9'], help = 'Version of CMD-3 data tree')
    parser.add_argument('--year', choices = ['2019', '2020', '2021', '2022', '2023'])
    parser.add_argument('--energy')
    parser.add_argument('--target-dir')
    parser.add_argument('--input-file')
    args = parser.parse_args()
    if not os.path.exists(args.target_dir): raise OSError(f"Target directory not existing: {args.target_dir}")
    if not os.path.exists(args.input_file): raise OSError(f"Input file not existing: {args.input_file}")

    ## CMD3 --> Preliminary
    ## Path to datafiles
    path_cmd3data = args.input_file
    path_prelim = os.path.normpath(f"{args.target_dir}/cut{args.year}_tr_ph_fc_e{args.energy}_{args.version}.root")

    ## Path to histograms
    path_hists = os.path.normpath(f"{args.target_dir}/hists{args.year}_tr_ph_fc_e{args.energy}_{args.version}.root")
    ## Path to log
    path_log = os.path.normpath(f"{args.target_dir}/preliminary_cut_{args.year}_e{args.energy}_{args.version}_{date.today().strftime('%d-%m-%Y')}.log")

    prelim_analysis = PreliminaryAnalysis(path_hists, path_cmd3data, path_prelim, path_log)
    prelim_analysis.addCut('nt')
    prelim_analysis.addCut('tcharge')
    prelim_analysis.addHistogram('h_tptot')
    prelim_analysis.addHistogram('h_tnhit')
    prelim_analysis.addHistogram('h_tth')
    prelim_analysis.addHistogram('h_tphi')
    prelim_analysis.addHistogram('h_tptot_tdedx')
    prelim_analysis.addHistogram('h_KpKmPipPimLklhd')
    prelim_analysis.loop()

    prelim_analysis.addCut('trho')
    prelim_analysis.addCut('tz')
    prelim_analysis.addCut('tptot')
    prelim_analysis.addHistogram('h_tnhit')
    prelim_analysis.addHistogram('h_tth')
    prelim_analysis.addHistogram('h_tphi')
    prelim_analysis.addHistogram('h_tptot_tdedx')
    prelim_analysis.addHistogram('h_KpKmPipPimLklhd')
    prelim_analysis.loop()
    
    prelim_analysis.addCut('tnhit')
    prelim_analysis.addHistogram('h_tth')
    prelim_analysis.addHistogram('h_tphi')
    prelim_analysis.addHistogram('h_tptot_tdedx')
    prelim_analysis.addHistogram('h_KpKmPipPimLklhd')
    prelim_analysis.loop()
    
    prelim_analysis.addCut('tth')
    prelim_analysis.addHistogram('h_tth')
    prelim_analysis.addHistogram('h_tptot_tdedx')
    prelim_analysis.addHistogram('h_KpKmPipPimLklhd')
    prelim_analysis.addHistogram('h_TotalP_DeltaE')
    prelim_analysis.loop()

    prelim_analysis.addCut('nph')
    prelim_analysis.addHistogram('h_tptot_tdedx')
    prelim_analysis.addHistogram('h_KpKmPipPimLklhd')
    prelim_analysis.addHistogram('h_TotalP_DeltaE')
    prelim_analysis.loop()
    
    prelim_analysis.addCut('KpKmPipPimLklhd')
    prelim_analysis.addHistogram('h_tptot_tdedx')
    prelim_analysis.addHistogram('h_KpKmPipPimLklhd')
    prelim_analysis.addHistogram('h_TotalP_DeltaE')
    prelim_analysis.loop()
    
    prelim_analysis.dumpToFile()
