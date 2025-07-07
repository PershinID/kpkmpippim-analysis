from math import pi, sqrt, inf

from ctypes import c_float
from copy import deepcopy
from collections import Counter
from itertools import permutations

from argparse import ArgumentParser
from datetime import date
import os

import ROOT
from ROOT import gInterpreter
from ROOT import TLorentzVector

from Base.Variable import Variable, TriggerVariable
from Base.PhysicalConstants import m_pi, m_K
from Base.Analysis import Analysis, CutDispatcher, HistogramDispatcher, EHists

from Containers.CMD3ContainerV9 import CMD3ContainerV9
from Containers.PreliminaryContainer import PreliminaryContainer

class PreliminaryAnalysis(Analysis):
    def __init__(self, input_path, *, analysis_path = None, output_path = None, log_path = None, is_sim = False):
        if not analysis_path and not output_path: raise ValueError("Analysis path and output path cannot be None at the same time")
        Analysis.__init__(self, analysis_path, logname = "preliminary_analysis", logpath = log_path)
        
        ## Likelihood calculation inclusion
        gInterpreter.ProcessLine('#include "/spoolA/idpershin/analysis/kpkmpippim/python/k_pi_dedx_v9_2025_par.h"')
        self.isSimulated = is_sim

        size_variables = { # needed to define variables in self.Variables
            "nt":               Variable("nt",                  "as", max_value = 10),
            "ntlxe":            Variable("ntlxe",               "as", max_value = 10),
            "nph":              Variable("nph",                 "as", max_value = 46),
        }
        self.Variables = {
            ## CMD3 branches
            "emeas":            Variable("emeas",               "f"),
            "demeas":           Variable("demeas",              "f"),
            "xbeam":            Variable("xbeam",               "f"),
            "ybeam":            Variable("ybeam",               "f"),
            "runnum":           Variable("runnum",              "i"),
            "evnum":            Variable("evnum",               "i"),
            "ecaltot":          Variable("ecaltot",             "f"),
            "ecalneu":          Variable("ecalneu",             "f"),
            "psumch":           Variable("psumch",              "f"),
            "psumnu":           Variable("psumnu",              "f"),
            "nv_total":         Variable("nv_total",            "i"),
            "nt_total":         Variable("nt_total",            "i"),
            "ntlxe_total":      Variable("ntlxe_total",         "i"),
            "z0":               Variable("z0",                  "f"),
            "tnhit":            Variable("tnhit",               "i", sizes = (size_variables["nt"],)),
            "tlength":          Variable("tlength",             "f", sizes = (size_variables["nt"],)),
            "tphi":             Variable("tphi",                "f", sizes = (size_variables["nt"],)),
            "tth":              Variable("tth",                 "f", sizes = (size_variables["nt"],)),
            "tptot":            Variable("tptot",               "f", sizes = (size_variables["nt"],)),
            "tphiv":            Variable("tphiv",               "f", sizes = (size_variables["nt"],)),
            "tthv":             Variable("tthv",                "f", sizes = (size_variables["nt"],)),
            "tptotv":           Variable("tptotv",              "f", sizes = (size_variables["nt"],)),
            "trho":             Variable("trho",                "f", sizes = (size_variables["nt"],)),
            "tz":               Variable("tz",                  "f", sizes = (size_variables["nt"],)),
            "tdedx":            Variable("tdedx",               "f", sizes = (size_variables["nt"],)),
            "tchi2r":           Variable("tchi2r",              "f", sizes = (size_variables["nt"],)),
            "tchi2z":           Variable("tchi2z",              "f", sizes = (size_variables["nt"],)),
            "tchi2ndf":         Variable("tchi2ndf",            "f", sizes = (size_variables["nt"],)),
            "tt0":              Variable("tt0",                 "f", sizes = (size_variables["nt"],)),
            "tant":             Variable("tant",                "f", sizes = (size_variables["nt"],)),
            "tcharge":          Variable("tcharge",             "i", sizes = (size_variables["nt"],)),
            "ten":              Variable("ten",                 "f", sizes = (size_variables["nt"],)),
            "tfc":              Variable("tfc",                 "f", sizes = (size_variables["nt"],)),
            "tenlxe":           Variable("tenlxe",              "f", sizes = (size_variables["nt"],)),
            "tlengthlxe":       Variable("tlengthlxe",          "f", sizes = (size_variables["nt"],)),
            "tenslxe_layers":   Variable("tenslxe_layers",      "f", sizes = (size_variables["nt"], 14,)),
            "tencsi":           Variable("tencsi",              "f", sizes = (size_variables["nt"],)),
            "tenbgo":           Variable("tenbgo",              "f", sizes = (size_variables["nt"],)),
            "tclth":            Variable("tclth",               "f", sizes = (size_variables["nt"],)),
            "tclphi":           Variable("tclphi",              "f", sizes = (size_variables["nt"],)),
            "terr":             Variable("terr",                "f", sizes = (size_variables["nt"], 3, 3,)),
            "terr0":            Variable("terr0",               "f", sizes = (size_variables["nt"], 6, 6,)),
            "tindlxe":          Variable("tindlxe",             "i", sizes = (size_variables["nt"],)),
            "txyzatcl":         Variable("txyzatcl",            "f", sizes = (size_variables["nt"], 3,)),
            "txyzatlxe":        Variable("txyzatlxe",           "f", sizes = (size_variables["nt"], 3,)),
            "tenconv":          Variable("tenconv",             "i", sizes = (size_variables["nt"],)),
            "ntlxelayers":      Variable("ntlxelayers",         "i", sizes = (size_variables["nt"],)),
            "tlxenhit":         Variable("tlxenhit",            "i", sizes = (size_variables["nt"],)),
            "tlxelength":       Variable("tlxelength",          "f", sizes = (size_variables["nt"],)),
            "tlxededx":         Variable("tlxededx",            "f", sizes = (size_variables["nt"],)),
            "tlxeir":           Variable("tlxeir",              "f", sizes = (size_variables["nt"],)),
            "tlxeitheta":       Variable("tlxeitheta",          "f", sizes = (size_variables["nt"],)),
            "tlxeiphi":         Variable("tlxeiphi",            "f", sizes = (size_variables["nt"],)),
            "tlxevtheta":       Variable("tlxevtheta",          "f", sizes = (size_variables["nt"],)),
            "tlxevphi":         Variable("tlxevphi",            "f", sizes = (size_variables["nt"],)),
            "tlxechi2":         Variable("tlxechi2",            "f", sizes = (size_variables["nt"],)),
            "tlxesen":          Variable("tlxesen",             "f", sizes = (size_variables["nt"],)),
            "tlxesen_layers":   Variable("tlxesen_layers",      "f", sizes = (size_variables["nt"], 14,)),
            "finalstate_id":    Variable("finalstate_id",       "i"),
            "nph_total":        Variable("nph_total",           "i"),
            "phen":             Variable("phen",                "f", sizes = (size_variables["nph"],)),
            "phth":             Variable("phth",                "f", sizes = (size_variables["nph"],)),
            "phphi":            Variable("phphi",               "f", sizes = (size_variables["nph"],)),
            "phrho":            Variable("phrho",               "f", sizes = (size_variables["nph"],)),
            "phen0":            Variable("phen0",               "f", sizes = (size_variables["nph"],)),
            "phth0":            Variable("phth0",               "f", sizes = (size_variables["nph"],)),
            "phphi0":           Variable("phphi0",              "f", sizes = (size_variables["nph"],)),
            "phlxe":            Variable("phlxe",               "f", sizes = (size_variables["nph"],)),
            "phslxe_layers":    Variable("phslxe_layers",       "f", sizes = (size_variables["nph"], 14,)),
            "pherr":            Variable("pherr",               "f", sizes = (size_variables["nph"], 3,)),
            "phcsi":            Variable("phcsi",               "f", sizes = (size_variables["nph"],)),
            "phbgo":            Variable("phbgo",               "f", sizes = (size_variables["nph"],)),
            "phflag":           Variable("phflag",              "i", sizes = (size_variables["nph"],)),
            "phconv":           Variable("phconv",              "i", sizes = (size_variables["nph"],)),
            "phfc":             Variable("phfc",                "i", sizes = (size_variables["nph"],)),

            "DeltaE":               TriggerVariable("DeltaE",          "f", self.calculateDeltaETotalP),
            "TotalP":               TriggerVariable("TotalP",          "f", self.calculateDeltaETotalP),
            "KpKmPipPimLklhd":              TriggerVariable("KpKmPipPimLklhd", "f", self.calculateKpKmPipPimLklhd),
            "KpTrackIndex":                 TriggerVariable("KpTrackIndex",    "B", self.calculateKpKmPipPimLklhd), ## [K+, K-, pi+, pi-]
            "KmTrackIndex":                 TriggerVariable("KmTrackIndex",    "B", self.calculateKpKmPipPimLklhd),
            "PipTrackIndex":                TriggerVariable("PipTrackIndex",   "B", self.calculateKpKmPipPimLklhd),
            "PimTrackIndex":                TriggerVariable("PimTrackIndex",   "B", self.calculateKpKmPipPimLklhd),
        }
        self.Variables.update(size_variables)

        self.InputContainers.append(CMD3ContainerV9(input_path, self.Variables))
        if output_path: self.OutputContainers.append(PreliminaryContainer(output_path, "recreate", self.Variables))

        cuts_available = {
            'nt': lambda: self.Variables["nt"].Content != 4,
            'tcharge': lambda: sum(self.Variables["tcharge"].Content) != 0,
            'tnhit': lambda: len(list(filter(
                lambda x: x <= 9,
                self.Variables["tnhit"].Content
            ))) != 0,
            'tptot': lambda: len(list(filter(
                lambda x: x < 50.,
                self.Variables["tptot"].Content
            ))) != 0,
            'tth': lambda: len(list(filter(
                lambda x: x < 0.9 or x > pi - 0.9,
                self.Variables["tth"].Content
            ))) != 0,
            'trho': lambda: len(list(filter(
                lambda x: abs(x) > 0.4,
                self.Variables["trho"].Content
            ))) != 0,
            'tz': lambda: len(list(filter(
                lambda x: abs(x) > 10.0,
                self.Variables["tz"].Content
            ))) != 0,

            'nph': lambda: self.Variables["nph"].Content > 1,
            'phth': lambda: len(list(filter(
                lambda x: x > 0.9 and x < pi - 0.9,
                self.Variables["phth"].Content
            ))) != 0,

            'finalstate_id': lambda: self.Variables["finalstate_id"].Content == 12, # K+K-pi+pi- multihad code = 12
        }
        self.CutDispatcher = CutDispatcher(cuts_available, n_entries_full = self.getEntries())

        histograms_available = {
            'h_tz': {
                'type': EHists.TH1F,
                'args': {
                    'title': "Track Z",
                    'x-axis-title': "z_{track}, cm",
                    'x-axis-nbins': 400,
                    'x-axis-range': (-20., 20.),
                    'x-variable': self.Variables["tz"],
                },
            },
            'h_trho': {
                'type': EHists.TH1F,
                'args': {
                    'title': "Track #rho",
                    'x-axis-title': "#rho_{track}, cm",
                    'x-axis-nbins': 4000,
                    'x-axis-range': (-2., 2.),
                    'x-variable': self.Variables["trho"],
                },
            },
            'h_tth': {
                'type': EHists.TH1F,
                'args': {
                    'title': "Track #theta",
                    'x-axis-title': "#theta_{track}, rad",
                    'x-axis-nbins': 600,
                    'x-axis-range': (0., pi),
                    'x-variable': self.Variables["tth"],
                },
            },
            'h_tphi': {
                'type': EHists.TH1F,
                'args': {
                    'title': "Track #phi",
                    'x-axis-title': "#phi_{track}, rad",
                    'x-axis-nbins': 600,
                    'x-axis-range': (0., 2. * pi),
                    'x-variable': self.Variables["tphi"],
                },
            },
            'h_tnhit': {
                'type': EHists.TH1I,
                'args': {
                    'title': "Hits on track",
                    'x-axis-title': "N_{hits}",
                    'x-axis-nbins': 35,
                    'x-axis-range': (5, 40),
                    'x-variable': self.Variables["tnhit"],
                },
            },
            'h_tptot': {
                'type': EHists.TH1F,
                'args': {
                    'title': "Track momentum",
                    'x-axis-title': "P_{track}, MeV/c",
                    'x-axis-nbins': 2000,
                    'x-axis-range': (0., 2000.),
                    'x-variable': self.Variables["tptot"],
                },
            },
            'h_tptot_tdedx': {
                'type': EHists.TH2F,
                'args': {
                    'title': "Track momentum vs Ionization losses",
                    'x-axis-title': "P_{track}, MeV/c",
                    'x-axis-nbins': 1000,
                    'x-axis-range': (0., 1000.),
                    'x-variable': self.Variables["tptot"],
                    'y-axis-title': "#frac{dE}{dx}_{track}, c.u.",
                    'y-axis-nbins': 500,
                    'y-axis-range': (0., 50000.),
                    'y-variable': self.Variables["tdedx"],
                },
            },
            'h_KpKmPipPimLklhd': {
                'type': EHists.TH1F,
                'args': {
                    'title': "K+K-#pi+#pi- hypothesis likelihood by ionization losses",
                    'x-axis-title': "Likelihood",
                    'x-axis-nbins': 150,
                    'x-axis-range': (-15., 0.),
                    'x-variable': self.Variables["KpKmPipPimLklhd"],
                },
            },
            'h_TotalP_DeltaE': {
                'type': EHists.TH2F,
                'args': {
                    'title': "Total momentum vs #DeltaE",
                    'x-axis-title': "P_{total}, MeV/c",
                    'x-axis-nbins': 1000,
                    'x-axis-range': (0., 1000.),
                    'x-variable': self.Variables["TotalP"],
                    'y-axis-title': "#DeltaE, MeV",
                    'y-axis-nbins': 1200,
                    'y-axis-range': (-800., 400.),
                    'y-variable': self.Variables["DeltaE"],
                },
            },
            'h_nph': {
                'type': EHists.TH1I,
                'args': {
                    'title': "Photons number",
                    'x-axis-title': "N_{photons}",
                    'x-axis-nbins': 10,
                    'x-axis-range': (0, 10),
                    'x-variable': self.Variables["nph"],
                },
            },
            'h_phen': {
                'type': EHists.TH1F,
                'args': {
                    'title': "Photon energy",
                    'x-axis-title': "E_{photon}, MeV",
                    'x-axis-nbins': 2000,
                    'x-axis-range': (0., 2000.),
                    'x-variable': self.Variables["phen"],
                },
            },
            'h_phth': {
                'type': EHists.TH1F,
                'args': {
                    'title': "Photon #theta",
                    'x-axis-title': "#theta_{photon}, rad",
                    'x-axis-nbins': 600,
                    'x-axis-range': (0., pi),
                    'x-variable': self.Variables["phth"],
                },
            },
            'h_finalstate_id': {
                'type': EHists.TH1I,
                'args': {
                    'title': "Final state ID",
                    'x-axis-title': "ID_{final state}",
                    'x-axis-nbins': 32,
                    'x-axis-range': (0, 32),
                    'x-variable': self.Variables["finalstate_id"],
                },
            },
        }
        if analysis_path: self.HistogramDispatcher = HistogramDispatcher(histograms_available)

        
    def calculateKpKmPipPimLklhd(self):
        ## not used in calculation, needed only for compatibility
        LklhdParsType = c_float * 12
        pars = LklhdParsType(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        ## Calculating likelihood functions in K and pi hypotheses for every track
        lklhds_by_tracks = ([-inf, -inf], [-inf, -inf], [-inf, -inf], [-inf, -inf]) ## tracks likelihoods (lklhd index = track index in tr_ph)
        for (track_lklhds, p, dedx) in zip(lklhds_by_tracks, self.Variables["tptot"].Content, self.Variables["tdedx"].Content):
            ## likelihood in pi and K hypotheses
            track_lklhds[0] = ROOT.test_k(p, dedx, self.Variables["runnum"].Content, pars, self.isSimulated, False) ## pi
            track_lklhds[1] = ROOT.test_k(p, dedx, self.Variables["runnum"].Content, pars, self.isSimulated, True) ## K

        ## Permutating tracks indices so that the charges is arranged as (+, -, +, -)
        tracks_indices = (0, 1, 2, 3)
        for perm in permutations(tracks_indices):
            prmttd_tcharge = tuple(self.Variables["tcharge"].Content[i] for i in perm)
            if prmttd_tcharge == (+1, -1, +1, -1):
                tracks_indices = perm
                break

        ## Calculating total likelihood by swapping a '+'-charged pair of tracks and a '-'-charged one
        max_total_lklhd = -inf
        max_total_lklhd_indices = tracks_indices
        tracks_permutations = (
            (0, 1, 2, 3),
            (2, 1, 0, 3),
            (0, 3, 2, 1),
            (2, 3, 0, 1)
        )
        tracks_permutations = tuple(map(
            lambda perm: tuple(tracks_indices[i] for i in perm),
            tracks_permutations
        ))
#        print("Current entry:", self.InputContainers[0].CurrentEntry)
        for perm in tracks_permutations:
#            print("Permutation:", perm)
#            print("Lklhds:", lklhds_by_tracks[ perm[0] ][1], lklhds_by_tracks[ perm[1] ][1], lklhds_by_tracks[ perm[2] ][0], lklhds_by_tracks[ perm[3] ][0])
            curr_total_lklhd = lklhds_by_tracks[ perm[0] ][1] + lklhds_by_tracks[ perm[1] ][1] + lklhds_by_tracks[ perm[2] ][0] + lklhds_by_tracks[ perm[3] ][0]
#            print("Sum lklhd:", curr_total_lklhd)
            if max_total_lklhd < curr_total_lklhd:
                max_total_lklhd = curr_total_lklhd
                max_total_lklhd_indices = perm
#        print("Max lklhd:", max_total_lklhd)
        ## As a result: max_total_lklhd_indices = [K+ index, K- index, pi+ index, pi- index]

        self.Variables["KpKmPipPimLklhd"].Content = max_total_lklhd
        self.Variables["KpTrackIndex"].Content  = max_total_lklhd_indices[0]
        self.Variables["KmTrackIndex"].Content  = max_total_lklhd_indices[1]
        self.Variables["PipTrackIndex"].Content = max_total_lklhd_indices[2]
        self.Variables["PimTrackIndex"].Content = max_total_lklhd_indices[3]

        
    def calculateDeltaETotalP(self):
        tp_lorentz = TLorentzVector(0.0, 0.0, 0.0, 0.0)
        tp_auxil = TLorentzVector(0.0, 0.0, 0.0, 0.0)
        
        for (p, th, phi) in zip(self.Variables["tptot"].Content,
                                self.Variables["tth"].Content,
                                self.Variables["tphi"].Content):
            tp_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
            tp_auxil.SetRho(p)
            tp_auxil.SetTheta(th)
            tp_auxil.SetPhi(phi)
            tp_auxil.SetE(sqrt(p * p + m_pi * m_pi))
            tp_lorentz += tp_auxil

        self.Variables["DeltaE"].Content = tp_lorentz.E() - 2.0 * self.Variables["emeas"].Content
        self.Variables["TotalP"].Content = tp_lorentz.Rho()


    def calculateEntry(self):
        self.calculateDeltaETotalP()
        self.calculateKpKmPipPimLklhd()

        
if __name__ == '__main__':
    ##Parsing input arguments
    parser = ArgumentParser()
    parser.add_argument('--version', default = 'v9', choices = ['v9'], help = 'Version of CMD-3 data tree')
    parser.add_argument('--year', choices = ['2019', '2020', '2021', '2022', '2023'], required = True)
    parser.add_argument('--energy', required = True)
    parser.add_argument('--input-path', required = True)
    parser.add_argument('--output-dir')
    parser.add_argument('--hists-dir')
    parser.add_argument('--is-multihad', action = 'store_true')
    parser.add_argument('--is-sim', action = 'store_true')
    args = parser.parse_args()

    if args.is_sim and args.is_multihad:
        raise ValueError("Is_sim and is_multihad flags exist at the same time")
    if not os.path.exists(args.input_path): raise OSError(f"Input path does not exist: {args.input_dir}")
    if not os.path.exists(args.output_dir): raise OSError(f"Output directory does not exist: {args.output_dir}")
    if not os.path.exists(args.hists_dir): raise OSError(f"Histograms directory does not exist: {args.hists_dir}")

    prefix = ''
    if args.is_sim:
        prefix = 'sim_'
    if args.is_multihad:
        prefix = 'multihad_'

    output_path = f"{args.output_dir}/{prefix}prelim_cut{args.year}_tr_ph_fc_e{args.energy}_{args.version}.root" if args.output_dir else None
    hists_path = f"{args.hists_dir}/{prefix}prelim_hists{args.year}_tr_ph_fc_e{args.energy}_{args.version}.root" if args.hists_dir else None
    log_path = f"{args.output_dir}/{prefix}prelim_cut_y{args.year}_e{args.energy}_{args.version}.log"

    ## CMD3 --> Preliminary
    prelim_analysis = PreliminaryAnalysis(args.input_path, analysis_path = hists_path, output_path = output_path, log_path = log_path, is_sim = args.is_sim or args.is_multihad)
    if args.is_multihad:
        prelim_analysis.addCut('finalstate_id')
        prelim_analysis.loop()
    
    prelim_analysis.addCut('nt')
    prelim_analysis.addHistogram('h_tz')
    prelim_analysis.addHistogram('h_trho')
    prelim_analysis.addHistogram('h_tnhit')
    prelim_analysis.addHistogram('h_tptot')
    prelim_analysis.addHistogram('h_tth')
    prelim_analysis.addHistogram('h_tptot_tdedx')
    prelim_analysis.addHistogram('h_phen')
    prelim_analysis.addHistogram('h_nph')
    prelim_analysis.addHistogram('h_phth')
    if args.is_multihad: prelim_analysis.addHistogram('h_finalstate_id')
    prelim_analysis.loop()

    prelim_analysis.addCut('tcharge')
    prelim_analysis.addHistogram('h_tz')
    prelim_analysis.addHistogram('h_trho')
    prelim_analysis.addHistogram('h_tnhit')
    prelim_analysis.addHistogram('h_tptot')
    prelim_analysis.addHistogram('h_tth')
    prelim_analysis.addHistogram('h_tptot_tdedx')
    prelim_analysis.addHistogram('h_KpKmPipPimLklhd')
    prelim_analysis.addHistogram('h_TotalP_DeltaE')
    prelim_analysis.addHistogram('h_phen')
    prelim_analysis.addHistogram('h_nph')
    prelim_analysis.addHistogram('h_phth')
    if args.is_multihad: prelim_analysis.addHistogram('h_finalstate_id')
    prelim_analysis.loop()
    
    prelim_analysis.addCut('tz')
    prelim_analysis.addHistogram('h_trho')
    prelim_analysis.addHistogram('h_tnhit')
    prelim_analysis.addHistogram('h_tptot')
    prelim_analysis.addHistogram('h_tth')
    prelim_analysis.addHistogram('h_tptot_tdedx')
    prelim_analysis.addHistogram('h_KpKmPipPimLklhd')
    prelim_analysis.addHistogram('h_TotalP_DeltaE')
    prelim_analysis.addHistogram('h_phen')
    prelim_analysis.addHistogram('h_nph')
    prelim_analysis.addHistogram('h_phth')
    if args.is_multihad: prelim_analysis.addHistogram('h_finalstate_id')
    prelim_analysis.loop()
    
    prelim_analysis.addCut('trho')
    prelim_analysis.addHistogram('h_tnhit')
    prelim_analysis.addHistogram('h_tptot')
    prelim_analysis.addHistogram('h_tth')
    prelim_analysis.addHistogram('h_tptot_tdedx')
    prelim_analysis.addHistogram('h_KpKmPipPimLklhd')
    prelim_analysis.addHistogram('h_TotalP_DeltaE')
    prelim_analysis.addHistogram('h_phen')
    prelim_analysis.addHistogram('h_nph')
    prelim_analysis.addHistogram('h_phth')
    if args.is_multihad: prelim_analysis.addHistogram('h_finalstate_id')
    prelim_analysis.loop()
    
    prelim_analysis.addCut('tnhit')
    prelim_analysis.addHistogram('h_tptot')
    prelim_analysis.addHistogram('h_tth')
    prelim_analysis.addHistogram('h_tptot_tdedx')
    prelim_analysis.addHistogram('h_KpKmPipPimLklhd')
    prelim_analysis.addHistogram('h_TotalP_DeltaE')
    prelim_analysis.addHistogram('h_phen')
    prelim_analysis.addHistogram('h_nph')
    prelim_analysis.addHistogram('h_phth')
    if args.is_multihad: prelim_analysis.addHistogram('h_finalstate_id')
    prelim_analysis.loop()
    
    prelim_analysis.addCut('tth')
    prelim_analysis.addHistogram('h_tptot')
    prelim_analysis.addHistogram('h_tptot_tdedx')
    prelim_analysis.addHistogram('h_KpKmPipPimLklhd')
    prelim_analysis.addHistogram('h_TotalP_DeltaE')
    prelim_analysis.addHistogram('h_phen')
    prelim_analysis.addHistogram('h_nph')
    prelim_analysis.addHistogram('h_phth')
    if args.is_multihad: prelim_analysis.addHistogram('h_finalstate_id')
    prelim_analysis.loop()

    prelim_analysis.addCut('tptot')
    prelim_analysis.addHistogram('h_tptot_tdedx')
    prelim_analysis.addHistogram('h_KpKmPipPimLklhd')
    prelim_analysis.addHistogram('h_TotalP_DeltaE')
    prelim_analysis.addHistogram('h_phen')
    prelim_analysis.addHistogram('h_nph')
    prelim_analysis.addHistogram('h_phth')
    if args.is_multihad: prelim_analysis.addHistogram('h_finalstate_id')
    prelim_analysis.loop()
    
    prelim_analysis.dumpToFile()
