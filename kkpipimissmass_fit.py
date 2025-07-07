from math import sqrt, pi

import json

from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor

import ROOT
from ROOT import TH1F, TF1, TFile, TCanvas
from ROOT import gROOT

def process_single(version, year, elabel_data):
    energy_point = elabel_data["scan-energy-point"]

    input_path = f"/home/idpershin/analysis/kpkmpippim/pershin_cut_w_TotalP-DeltaE-2/hists{year}_tr_ph_fc_e{energy_point}_{version}.root"
    output_path = f"/home/idpershin/analysis/kpkmpippim/pershin_cut_w_TotalP-DeltaE-2/kkpipimissmass_fit/y{year}_e{energy_point}_{version}.root"

    f_missmass = TFile(input_path, "read")
    h_missmass = f_missmass.GetDirectory("2_TotalP-DeltaE-2").Get("h_KKPiPiMissMass2")
    h_missmass.Rebin(5)
    func_fit = TF1("func_fit", "gausn(0) + gausn(3) + gausn(6) + gausn(9)", -15_000., 35_000.)
    func_fit.SetNpx(1000)
    
    func_fit.SetParameter(1, -260.)
    func_fit.SetParLimits(1, -240., -280.)
    func_fit.SetParameter(2, 315.)
    func_fit.SetParLimits(2, 300., 330.)
    func_fit.SetParameter(4, -750.)
    func_fit.SetParLimits(4, -730., -770.)
    func_fit.SetParameter(5, 1_050.)
    func_fit.SetParLimits(5, 1_000., 1_100.)
    func_fit.SetParameter(7, -700.)
    func_fit.SetParLimits(7, -650., -750.)
    func_fit.SetParameter(8, 4_500.)
    func_fit.SetParLimits(8, 3_500., 5_500.)
    func_fit.SetParameter(10, 19_000.)
    func_fit.SetParLimits(10, 18_000., 20_000.)
    func_fit.SetParameter(11, 7_700.)
    func_fit.SetParLimits(11, 7_200., 8_200.)
    h_missmass.Fit(func_fit, "", "", -15_000., 35_000.)

    bin_width = h_missmass.GetBinWidth(4000)
    n_sig = 0
    n_sig_err2 = 0.
    for i in range(3):
        norm_factor = func_fit.GetParameter(i * 3)
        norm_factor_err = func_fit.GetParError(i * 3)
        n_sig += int(norm_factor / bin_width)
        n_sig_err2 += (norm_factor_err / bin_width) ** 2
    n_sig_err = int(sqrt(n_sig_err2))
    elabel_data["signal-events-number"] = n_sig
    elabel_data["signal-events-number-error"] = n_sig_err
    print("Signal events number:", n_sig, "\nSignal events number error:", n_sig_err)
    
    norm_factor = func_fit.GetParameter(9)
    norm_factor_err = func_fit.GetParError(9)
    n_bkg = int(norm_factor / bin_width)
    n_bkg_err = int(norm_factor_err / bin_width)
    elabel_data["background-events-number"] = n_bkg
    elabel_data["background-events-number-error"] = n_bkg_err
    print("Background events number:", n_bkg, "\nBackground events number error:", n_bkg_err)

    f_fit = TFile(output_path, "recreate")
    gROOT.SetBatch(ROOT.kTRUE)
    c_missmass = TCanvas("c_missmass", "c_missmass")
    h_missmass.Draw()
    func_fit.Draw("lsame")
    c_missmass.Write()
    f_fit.Save()
    c_missmass.Close()
    gROOT.SetBatch(ROOT.kFALSE)
    

def process_all(version, version_data):
    with ProcessPoolExecutor(max_workers = 10) as executor:
        futures = []
        for year in version_data["years"]:
            for elabel in version_data["years"][year]["elabels"]:
                elabel_data = version_data["years"][year]["elabels"][elabel]

                futures.append(
                    executor.submit(
                        process_single,
                        version,
                        year,
                        elabel_data,
                    )
                )
        wait(futures) 
    
if __name__ == '__main__':
    ##Parsing input arguments
    parser = ArgumentParser()
    parser.add_argument('--version', default = 'v9', choices = ['v9'], help = 'Version of CMD-3 data tree')
    subparsers = parser.add_subparsers(dest = 'mode')
    
    parser_single = subparsers.add_parser('single', help = "Process one energy point")
    parser_single.add_argument('--year', choices = ['2019', '2020', '2021', '2022', '2023'], required = True)
    parser_single.add_argument('--elabel', required = True)

    parser_all = subparsers.add_parser('all', help = "Process all available energy points")
    args = parser.parse_args()

    fit_data_path = "/home/idpershin/analysis/kpkmpippim/data_kkpipimissmass_fit.json"
    fit_data = json.load(open(fit_data_path, "r"))
    
    if args.mode == 'single':
        elabel_data  = fit_data[args.version]["years"][args.year]["elabels"][args.elabel]
        process_single(args.version, args.year, elabel_data)

    if args.mode == 'all':
        version_data  = fit_data[args.version]
        process_all(args.version, version_data)

    json.dump(fit_data, open(fit_data_path, "w"), indent = 4)
