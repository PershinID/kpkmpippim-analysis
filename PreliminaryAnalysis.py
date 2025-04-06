from math import pi, sqrt, inf

from ctypes import c_float
from copy import deepcopy
from collections import Counter
from itertools import permutations

import ROOT
from ROOT import gInterpreter
from ROOT import TLorentzVector

from Variable import Variable, TriggerVariable

from PhysicalConstants import m_pi

from CMD3ContainerV9 import CMD3ContainerV9
from PreliminaryContainer import PreliminaryContainer

from Analysis import Analysis, CutDispatcher, HistogramDispatcher, EHists

class PreliminaryAnalysis(Analysis):
    def __init__(self, path: str, input_path: str, output_path: str, log_path):
        Analysis.__init__(self, path, logname = "preliminary_analysis", logpath = log_path)

        ## Likelihood calculation inclusion
        gInterpreter.ProcessLine('#include "k_pi_dedx_v9_2025_par.h"')
        
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

            "KpTrackIndex":         TriggerVariable("KpTrackIndex",  "B", self.calculateKpKmPipPimLklhd), ## [K+, K-, pi+, pi-]
            "KmTrackIndex":         TriggerVariable("KmTrackIndex",  "B", self.calculateKpKmPipPimLklhd),
            "PipTrackIndex":        TriggerVariable("PipTrackIndex", "B", self.calculateKpKmPipPimLklhd),
            "PimTrackIndex":        TriggerVariable("PimTrackIndex", "B", self.calculateKpKmPipPimLklhd),
            "DeltaE":               TriggerVariable("DeltaE", "f", self.calculateDeltaETotalP),
            "TotalP":               TriggerVariable("TotalP", "f", self.calculateDeltaETotalP),
            "KpKmPipPimLklhd":      TriggerVariable("KpKmPipPimLklhd", "f", self.calculateKpKmPipPimLklhd),
        }
        self.Variables.update(size_variables)

        self.InputContainer = CMD3ContainerV9(input_path, self.Variables)
        self.OutputContainer = PreliminaryContainer(output_path, "recreate", self.Variables)

        cuts_available = {
            'nt': lambda: self.Variables["nt"].Content != 4,
            'nph': lambda: self.Variables["nph"].Content > 1,
            'tcharge': lambda: sum(self.Variables["tcharge"].Content) != 0,
            'tnhit': lambda: len(list(filter(
                lambda x: x <= 9,
                self.Variables["tnhit"].Content
            ))) != 0,
            'tptot': lambda: len(list(filter(
                lambda x: x < 50,
                self.Variables["tptot"].Content
            ))) != 0,
            'tth': lambda: len(list(filter(
                lambda x: x < 0.9 or x > pi - 0.9,
                self.Variables["tth"].Content
            ))) != 0,
            'trho': lambda: len(list(filter(
                lambda x: abs(x) > 0.2,
                self.Variables["trho"].Content
            ))) != 0,
            'tz': lambda: len(list(filter(
                lambda x: abs(x) > 10.0,
                self.Variables["tz"].Content
            ))) != 0,
            'DeltaE': lambda: (
                self.Variables["DeltaE"].Content < -700.0 or
                self.Variables["DeltaE"].Content > -300.0
            ),
            'TotalP': lambda: (
                self.Variables["TotalP"].Content > 100.0
            ),
            'TotalP-DeltaE': lambda: (
                self.Variables["DeltaE"].Content > -200.0 - 1.5 * self.Variables["DeltaE"].Content
            ),
            'KpKmPipPimLklhd': lambda: (
                self.Variables["KpKmPipPimLklhd"].Content < -3.0
            ),
        }
        self.CutDispatcher = CutDispatcher(cuts_available, n_entries_full = self.InputContainer.getEntries())

        histograms_available = {
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
                    'x-axis-range': (0.0, 2000.0),
                    'x-variable': self.Variables["tptot"],
                },
            },
            'h_tth': {
                'type': EHists.TH1F,
                'args': {
                    'title': "Track #theta coordinate",
                    'x-axis-title': "#theta_{track}, rad",
                    'x-axis-nbins': 600,
                    'x-axis-range': (0.0, pi),
                    'x-variable': self.Variables["tth"],
                },
            },
            'h_tphi': {
                'type': EHists.TH1F,
                'args': {
                    'title': "Track #phi",
                    'x-axis-title': "#phi_{track}, rad",
                    'x-axis-nbins': 600,
                    'x-axis-range': (0.0, 2 * pi),
                    'x-variable': self.Variables["tphi"],
                },
            },
            'h_tptot_tdedx': {
                'type': EHists.TH2F,
                'args': {
                    'title': "Track momentum vs ionization losses",
                    'x-axis-title': "P_{track}, MeV/c",
                    'x-axis-nbins': 1000,
                    'x-axis-range': (0.0, 1000.0),
                    'x-variable': self.Variables["tptot"],
                    'y-axis-title': "#frac{dE}{dx}_{track}, a.u.",
                    'y-axis-nbins': 500,
                    'y-axis-range': (0.0, 50000.0),
                    'y-variable': self.Variables["tdedx"],
                },
            },
            'h_TotalP_DeltaE': {
                'type': EHists.TH2F,
                'args': {
                    'title': "TotalP vs DeltaE",
                    'x-axis-title': "P_{total}, MeV/c",
                    'x-axis-nbins': 1000,
                    'x-axis-range': (0.0, 1000.0),
                    'x-variable': self.Variables["TotalP"],
                    'y-axis-title': "#DeltaE, MeV",
                    'y-axis-nbins': 1200,
                    'y-axis-range': (-800.0, 400.0),
                    'y-variable': self.Variables["DeltaE"],
                },
            },
            'h_KpKmPipPimLklhd': {
                'type': EHists.TH1F,
                'args': {
                    'title': "Likelihood in K+K-#pi+#pi- hypothesis",
                    'x-axis-title': "Likelihood",
                    'x-axis-nbins': 100,
                    'x-axis-range': (-10.0, 0.0),
                    'x-variable': self.Variables["KpKmPipPimLklhd"],
                },
            },
            'h_nph': {
                'type': EHists.TH1I,
                'args': {
                    'title': "Photons number",
                    'x-axis-title': "N_{photons}",
                    'x-axis-nbins': 10,
                    'x-axis-range': (0, 10),
                    'x-variable': self.Variables["KpKmPipPimLklhd"],
                },
            },
            'h_PipPimPipPimMissMass': {
                'type': EHists.TH1F,
                'args': {
                    'title': "#pi+#pi-#pi+#pi- missing mass",
                    'x-axis-title': "M_{#pi+#pi-#pi+#pi-,miss}",
                    'x-axis-nbins': 1000,
                    'x-axis-range': (100.0, 1100.0),
                    'x-variable': self.Variables["PipPimPipPimMissMass"],
                },
            },
        }
        self.HistogramDispatcher = HistogramDispatcher(histograms_available)

        
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


    def calculateKpKmPipPimLklhd(self):
        ## not used in calculation, needed only for compatibility
        LklhdParsType = c_float * 12
        pars = LklhdParsType(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        ## Calculating likelihood functions in K and pi hypotheses for every track
        lklhds_by_tracks = ([-inf, -inf], [-inf, -inf], [-inf, -inf], [-inf, -inf]) ## tracks likelihoods (lklhd index = track index in tr_ph)
        for (track_lklhds, p, dedx) in zip(lklhds_by_tracks, self.Variables["tptot"].Content, self.Variables["tdedx"].Content):
            ## likelihood in pi and K hypotheses
            track_lklhds[0] = ROOT.test_k(p, dedx, self.Variables["runnum"].Content, pars, False, False) ## pi
            track_lklhds[1] = ROOT.test_k(p, dedx, self.Variables["runnum"].Content, pars, False, True) ## K

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
            (1, 0, 2, 3),
            (0, 1, 3, 2),
            (1, 0, 3, 2)
        )
        tracks_permutations = tuple(map(
            lambda perm: tuple(tracks_indices[i] for i in perm),
            tracks_permutations
        ))
        for perm in tracks_permutations:
            curr_total_lklhd = lklhds_by_tracks[ perm[0] ][1] + lklhds_by_tracks[ perm[1] ][1] + lklhds_by_tracks[ perm[2] ][0] + lklhds_by_tracks[ perm[3] ][0]
            if max_total_lklhd < curr_total_lklhd:
                max_total_lklhd = curr_total_lklhd
                max_total_lklhd_indices = perm                
        ## As a result: max_total_lklhd_indices = [K+ index, K- index, pi+ index, pi- index]

        self.Variables["KpKmPipPimLklhd"].Content = max_total_lklhd
        self.Variables["KpTrackIndex"].Content  = max_total_lklhd_indices[0]
        self.Variables["KmTrackIndex"].Content  = max_total_lklhd_indices[1]
        self.Variables["PipTrackIndex"].Content = max_total_lklhd_indices[2]
        self.Variables["PimTrackIndex"].Content = max_total_lklhd_indices[3]


    def calculatePipPimPipPimMissMass(self):
        p4_0     = TLorentzVector(0.0, 0.0, 0.0, 2.0 * self.Variables["emeas"].Content)
        p4_auxil = TLorentzVector(0.0, 0.0, 0.0, 0.0)
        
        for (p, th, phi) in zip(self.Variables["tptot"].Content,
                                self.Variables["tth"].Content,
                                self.Variables["tphi"].Content):
            p4_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
            p4_auxil.SetRho(p)
            p4_auxil.SetTheta(th)
            p4_auxil.SetPhi(phi)
            p4_auxil.SetE(sqrt(p * p + m_pi * m_pi))
            p4_0 -= p4_auxil
        
        self.Variables["PipPimPipPimMissMass"].Content = p4_0.M()

    def calculateEntry(self):
        self.calculateDeltaETotalP()
        self.calculateKpKmPipPimLklhd()
        self.calculatePipPimPipPimMissMass()

        
if __name__ == '__main__':
    test_analysis = PreliminaryAnalysis('test_hists.root', '/spoolA/idpershin/datasets/dataset_e950.root', 'test_data.root', 'test_log.log')
    test_analysis.addCut('DeltaE')
    print(test_analysis.CutDispatcher.isCutCurrent('DeltaE'))
    test_analysis.getEntry(5)
    print(test_analysis.Variables["DeltaE"].Content, test_analysis.Variables["TotalP"].Content)
