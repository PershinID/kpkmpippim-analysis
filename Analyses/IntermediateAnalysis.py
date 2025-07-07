from math import pi, sqrt, inf
from ctypes import c_float
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

from Containers.PreliminaryContainer import PreliminaryContainer

class IntermediateAnalysis(Analysis):
    def __init__(self, input_path, *, analysis_path = None, output_path = None, log_path = None):
        Analysis.__init__(self, analysis_path, logname = "final_cut", logpath = log_path)

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
            "phflag":           Variable("phflag",              "f", sizes = (size_variables["nph"],)),
            "phconv":           Variable("phconv",              "f", sizes = (size_variables["nph"],)),
            "phfc":             Variable("phfc",                "f", sizes = (size_variables["nph"],)),

            "DeltaE":                       Variable("DeltaE",                       "f"),
            "TotalP":                       Variable("TotalP",                       "f"),
            
            "KpKmPipPimLklhd":              Variable("KpKmPipPimLklhd", "f"),
            "KpTrackIndex":                 Variable("KpTrackIndex",    "B"),
            "KmTrackIndex":                 Variable("KmTrackIndex",    "B"),
            "PipTrackIndex":                Variable("PipTrackIndex",   "B"),
            "PimTrackIndex":                Variable("PimTrackIndex",   "B"),
            
            "PiPiPiPiMissMass2":            TriggerVariable("PiPiPiPiMissMass2",      "f", self.calculatePiPiPiPiMissMass2),
            "KPiPiPiMissMass2":            TriggerVariable("KPiPiPiMissMass2",      "f", self.calculateKPiPiPiMissMass2),
            "PiPiPiMissMass2":              TriggerVariable("PiPiPiMissMass2",        "f", self.calculatePiPiPiMissMass2),
            "KPiPiMissMass2":               TriggerVariable("KPiPiMissMass2",         "f", self.calculateKPiPiMissMass2),
            "KKPiPiMissMass2":              TriggerVariable("KKPiPiMissMass2",        "f", self.calculateKKPiPiMissMass2),
            "KKMissMass":                   TriggerVariable("KKMissMass",             "f", self.calculateKKMissMass),
            "PiPiMissMass":                 TriggerVariable("PiPiMissMass",           "f", self.calculatePiPiMissMass),

            "DeltaEKKPiPi":                 TriggerVariable("DeltaKKPiPi",            "f", self.calculateDeltaEKKPiPi),
        }
        self.Variables.update(size_variables)

        self.InputContainers = [PreliminaryContainer(input_path, "read", self.Variables),]
        if output_path: self.OutputContainers = [PreliminaryContainer(output_path, "recreate", self.Variables),]

        cuts_available = {
            'nph': lambda: self.Variables["nph"].Content > 1,
            'phth': lambda: len(list(filter(
                lambda x: x > 0.9 and x < pi - 0.9,
                self.Variables["phth"].Content
            ))) != 0,
            'phen': lambda: len(list(filter(
                lambda x: x < 50.,
                self.Variables["phen"].Content
            ))) != 0,
            'KpKmPipPimLklhd': lambda: (
                self.Variables["KpKmPipPimLklhd"].Content < -3.
            ),
            'TotalP-DeltaE': lambda: (
                self.Variables["DeltaE"].Content > -150. - 1. * self.Variables["TotalP"].Content
            ),
            'TotalP': lambda: self.Variables["TotalP"].Content > 80.,
            'DeltaEKKPiPi': lambda: abs(self.Variables["DeltaEKKPiPi"].Content) > 80.,
            'finalstate_id': lambda: self.Variables["finalstate_id"].Content == 12,
        }
        self.CutDispatcher = CutDispatcher(cuts_available, n_entries_full = self.getEntries())

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
                    'x-axis-range': (0., 2000.),
                    'x-variable': self.Variables["tptot"],
                },
            },
            'h_tth': {
                'type': EHists.TH1F,
                'args': {
                    'title': "Track #theta coordinate",
                    'x-axis-title': "#theta_{track}, rad",
                    'x-axis-nbins': 300,
                    'x-axis-range': (0., pi),
                    'x-variable': self.Variables["tth"],
                },
            },
            'h_trho': {
                'type': EHists.TH1F,
                'args': {
                    'title': "Track #rho coordinate",
                    'x-axis-title': "#rho_{track}, cm",
                    'x-axis-nbins': 400,
                    'x-axis-range': (-20., 20.),
                    'x-variable': self.Variables["trho"],
                },
            },
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
            'h_tptot_tdedx': {
                'type': EHists.TH2F,
                'args': {
                    'title': "Ionization losses vs track momentum",
                    'x-axis-title': "P_{track}, MeV/c",
                    'x-axis-nbins': 1000,
                    'x-axis-range': (0., 1000.),
                    'x-variable': self.Variables["tptot"],
                    'y-axis-title': "#frac{dE}{dx}_{track}, a.u.",
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
                    'title': "#DeltaE vs total momentum",
                    'x-axis-title': "P_{total}, MeV/c",
                    'x-axis-nbins': 1000,
                    'x-axis-range': (0., 1000.),
                    'x-variable': self.Variables["TotalP"],
                    'y-axis-title': "#DeltaE, MeV",
                    'y-axis-nbins': 1200,
                    'y-axis-range': (-1100., 400.),
                    'y-variable': self.Variables["DeltaE"],
                },
            },
            'h_TotalP_DeltaEKKPiPi': {
                'type': EHists.TH2F,
                'args': {
                    'title': "#DeltaE_{KK#pi#pi} vs total momentum",
                    'x-axis-title': "P_{total}, MeV/c",
                    'x-axis-nbins': 1000,
                    'x-axis-range': (0., 1000.),
                    'x-variable': self.Variables["TotalP"],
                    'y-axis-title': "#DeltaE_{KK#pi#pi}, MeV",
                    'y-axis-nbins': 1500,
                    'y-axis-range': (-700., 800.),
                    'y-variable': self.Variables["DeltaEKKPiPi"],
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
                    'x-axis-nbins': 300,
                    'x-axis-range': (0., pi),
                    'x-variable': self.Variables["phth"],
                },
            },
            'h_PiPiPiMissMass2': {
                'type': EHists.TH1F,
                'args': {
                    'title': "#pi#pi#pi missing mass squared",
                    'x-axis-title': "m_{#pi#pi#pi,miss}^{2}, #frac{MeV^{2}}{c^{4}}",
                    'x-axis-nbins': 3000,
                    'x-axis-range': (-1_000_000., 2_000_000.),
                    'x-variable': self.Variables["PiPiPiMissMass2"],
                },
            },
            'h_KPiPiMissMass2': {
                'type': EHists.TH1F,
                'args': {
                    'title': "K#pi#pi missing mass squared",
                    'x-axis-title': "m_{K#pi#pi,miss}^{2}, #frac{MeV^{2}}{c^{4}}",
                    'x-axis-nbins': 2500,
                    'x-axis-range': (-1_000_000., 1_500_000.),
                    'x-variable': self.Variables["KPiPiMissMass2"],
                },
            },
            'h_PiPiPiPiMissMass2': {
                'type': EHists.TH1F,
                'args': {
                    'title': "#pi^{+}#pi^{-}#pi^{+}#pi^{-} missing mass squared",
                    'x-axis-title': "m_{#pi+#pi-#pi+#pi-,miss}^{2}, #frac{MeV^{2}}{c^{4}}",
                    'x-axis-nbins': 1000,
                    'x-axis-range': (-200_000., 800_000.),
                    'x-variable': self.Variables["PiPiPiPiMissMass2"],
                },
            },
            'h_KPiPiPiMissMass2': {
                'type': EHists.TH1F,
                'args': {
                    'title': "K^{+}#pi^{-}#pi^{+}#pi^{-} missing mass squared",
                    'x-axis-title': "m_{K+#pi-#pi+#pi-,miss}^{2}, #frac{MeV^{2}}{c^{4}}",
                    'x-axis-nbins': 1000,
                    'x-axis-range': (-200_000., 800_000.),
                    'x-variable': self.Variables["KPiPiPiMissMass2"],
                },
            },
            'h_KKPiPiMissMass2': {
                'type': EHists.TH1F,
                'args': {
                    'title': "K^{+}K^{-}#pi^{+}#pi^{-} missing mass squared",
                    'x-axis-title': "m_{K^{+}K^{-}#pi^{+}#pi^{-},miss}, #frac{MeV^{2}}{c^{4}}",
                    'x-axis-nbins': 8000,
                    'x-axis-range': (-400_000., 400_000.),
                    'x-variable': self.Variables["KKPiPiMissMass2"],
                },
            },
            'h_KKMissMass': {
                'type': EHists.TH1F,
                'args': {
                    'title': "K^{+}K^{-} missing mass squared",
                    'x-axis-title': "m_{K^{+}K^{-},miss}, #frac{MeV^{2}}{c^{4}}",
                    'x-axis-nbins': 1000,
                    'x-axis-range': (0., 1000.),
                    'x-variable': self.Variables["KKMissMass"],
                }
            },
            'h_PiPiMissMass': {
                'type': EHists.TH1F,
                'args': {
                    'title': "#pi^{+}#pi^{-} missing mass squared",
                    'x-axis-title': "m_{#pi^{+}#pi^{-},miss}, #frac{MeV^{2}}{c^{4}}",
                    'x-axis-nbins': 3000,
                    'x-axis-range': (0., 3000.),
                    'x-variable': self.Variables["PiPiMissMass"],
                }
            },
            'h_finalstate_id': {
                'type': EHists.TH1I,
                'args': {
                    'title': "Final state ID",
                    'x-axis-title': "ID_{final state}",
                    'x-axis-nbins': 34,
                    'x-axis-range': (0, 34),
                    'x-variable': self.Variables["finalstate_id"],
                },
            },
        }
        self.HistogramDispatcher = HistogramDispatcher(histograms_available)


    def calculatePiPiPiPiMissMass2(self):
        p4_0     = TLorentzVector(0.0, 0.0, 0.0, 2.0 * self.Variables["emeas"].Content)
        p4_auxil = TLorentzVector(0.0, 0.0, 0.0, 0.0)
        
        for (p, th, phi) in zip(self.Variables["tptot"].Content,
                                self.Variables["tth"].Content,
                                self.Variables["tphi"].Content):
            p4_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
            p4_auxil.SetRho(p)
            p4_auxil.SetTheta(th)
            p4_auxil.SetPhi(phi)
            p4_auxil.SetE(sqrt(p ** 2 + m_pi ** 2))
            p4_0 -= p4_auxil
        
        self.Variables["PiPiPiPiMissMass2"].Content = p4_0.M2()


    def calculateKPiPiPiMissMass2(self):
        p4_0     = TLorentzVector(0.0, 0.0, 0.0, 2.0 * self.Variables["emeas"].Content)
        p4_auxil = TLorentzVector(0.0, 0.0, 0.0, 0.0)
        
        for i, (p, th, phi) in enumerate(zip(self.Variables["tptot"].Content,
                                             self.Variables["tth"].Content,
                                             self.Variables["tphi"].Content)):
            mass = m_pi
            if i == self.Variables["KpTrackIndex"].Content: mass = m_K
            
            p4_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
            p4_auxil.SetRho(p)
            p4_auxil.SetTheta(th)
            p4_auxil.SetPhi(phi)
            p4_auxil.SetE(sqrt(p ** 2 + mass ** 2))
            p4_0 -= p4_auxil
            
        self.Variables["KPiPiPiMissMass2"].Content = p4_0.M2()


    def calculatePiPiPiMissMass2(self):
        p4_0     = TLorentzVector(0.0, 0.0, 0.0, 2.0 * self.Variables["emeas"].Content)
        p4_auxil = TLorentzVector(0.0, 0.0, 0.0, 0.0)

        for i, (p, th, phi) in enumerate(zip(self.Variables["tptot"].Content,
                                             self.Variables["tth"].Content,
                                             self.Variables["tphi"].Content)):
            if i >= 3: break
            p4_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
            p4_auxil.SetRho(p)
            p4_auxil.SetTheta(th)
            p4_auxil.SetPhi(phi)
            p4_auxil.SetE(sqrt(p ** 2 + m_pi ** 2))
            p4_0 -= p4_auxil

        self.Variables["PiPiPiMissMass2"].Content = p4_0.M2()


    def calculateKPiPiMissMass2(self):
        p4_0     = TLorentzVector(0.0, 0.0, 0.0, 2.0 * self.Variables["emeas"].Content)
        p4_auxil = TLorentzVector(0.0, 0.0, 0.0, 0.0)

        index = self.Variables["PipTrackIndex"].Content
        p4_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
        p4_auxil.SetRho(self.Variables["tptot"].Content[index])
        p4_auxil.SetTheta(self.Variables["tth"].Content[index])
        p4_auxil.SetPhi(self.Variables["tphi"].Content[index])
        p4_auxil.SetE(sqrt(self.Variables["tptot"].Content[index] ** 2 + m_pi ** 2))
        p4_0 -= p4_auxil

        index = self.Variables["PimTrackIndex"].Content
        p4_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
        p4_auxil.SetRho(self.Variables["tptot"].Content[index])
        p4_auxil.SetTheta(self.Variables["tth"].Content[index])
        p4_auxil.SetPhi(self.Variables["tphi"].Content[index])
        p4_auxil.SetE(sqrt(self.Variables["tptot"].Content[index] ** 2 + m_pi ** 2))
        p4_0 -= p4_auxil

        index = self.Variables["KpTrackIndex"].Content
        p4_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
        p4_auxil.SetRho(self.Variables["tptot"].Content[index])
        p4_auxil.SetTheta(self.Variables["tth"].Content[index])
        p4_auxil.SetPhi(self.Variables["tphi"].Content[index])
        p4_auxil.SetE(sqrt(self.Variables["tptot"].Content[index] ** 2 + m_K ** 2))
        p4_0 -= p4_auxil
        
        self.Variables["KPiPiMissMass2"].Content = p4_0.M2()


    def calculateKKPiPiMissMass2(self):
        p4_0     = TLorentzVector(0.0, 0.0, 0.0, 2.0 * self.Variables["emeas"].Content)
        p4_auxil = TLorentzVector(0.0, 0.0, 0.0, 0.0)
        
        index = self.Variables["PipTrackIndex"].Content
        p4_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
        p4_auxil.SetRho(self.Variables["tptot"].Content[index])
        p4_auxil.SetTheta(self.Variables["tth"].Content[index])
        p4_auxil.SetPhi(self.Variables["tphi"].Content[index])
        p4_auxil.SetE(sqrt(self.Variables["tptot"].Content[index] ** 2 + m_pi ** 2))
        p4_0 -= p4_auxil

        index = self.Variables["PimTrackIndex"].Content
        p4_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
        p4_auxil.SetRho(self.Variables["tptot"].Content[index])
        p4_auxil.SetTheta(self.Variables["tth"].Content[index])
        p4_auxil.SetPhi(self.Variables["tphi"].Content[index])
        p4_auxil.SetE(sqrt(self.Variables["tptot"].Content[index] ** 2 + m_pi ** 2))
        p4_0 -= p4_auxil

        index = self.Variables["KpTrackIndex"].Content
        p4_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
        p4_auxil.SetRho(self.Variables["tptot"].Content[index])
        p4_auxil.SetTheta(self.Variables["tth"].Content[index])
        p4_auxil.SetPhi(self.Variables["tphi"].Content[index])
        p4_auxil.SetE(sqrt(self.Variables["tptot"].Content[index] ** 2 + m_K ** 2))
        p4_0 -= p4_auxil
        
        index = self.Variables["KmTrackIndex"].Content
        p4_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
        p4_auxil.SetRho(self.Variables["tptot"].Content[index])
        p4_auxil.SetTheta(self.Variables["tth"].Content[index])
        p4_auxil.SetPhi(self.Variables["tphi"].Content[index])
        p4_auxil.SetE(sqrt(self.Variables["tptot"].Content[index] ** 2 + m_K ** 2))
        p4_0 -= p4_auxil
        
        self.Variables["KKPiPiMissMass2"].Content = p4_0.M2()


    def calculateKKMissMass(self):
        p4_0     = TLorentzVector(0.0, 0.0, 0.0, 2.0 * self.Variables["emeas"].Content)
        p4_auxil = TLorentzVector(0.0, 0.0, 0.0, 0.0)
        
        index = self.Variables["KpTrackIndex"].Content
        p4_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
        p4_auxil.SetRho(self.Variables["tptot"].Content[index])
        p4_auxil.SetTheta(self.Variables["tth"].Content[index])
        p4_auxil.SetPhi(self.Variables["tphi"].Content[index])
        p4_auxil.SetE(sqrt(self.Variables["tptot"].Content[index] ** 2 + m_K * m_K))
        p4_0 -= p4_auxil
        
        index = self.Variables["KmTrackIndex"].Content
        p4_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
        p4_auxil.SetRho(self.Variables["tptot"].Content[index])
        p4_auxil.SetTheta(self.Variables["tth"].Content[index])
        p4_auxil.SetPhi(self.Variables["tphi"].Content[index])
        p4_auxil.SetE(sqrt(self.Variables["tptot"].Content[index] ** 2 + m_K * m_K))
        p4_0 -= p4_auxil
        
        self.Variables["KKMissMass"].Content = p4_0.M()


    def calculatePiPiMissMass(self):
        p4_0     = TLorentzVector(0.0, 0.0, 0.0, 2.0 * self.Variables["emeas"].Content)
        p4_auxil = TLorentzVector(0.0, 0.0, 0.0, 0.0)
        
        index = self.Variables["PipTrackIndex"].Content
        p4_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
        p4_auxil.SetRho(self.Variables["tptot"].Content[index])
        p4_auxil.SetTheta(self.Variables["tth"].Content[index])
        p4_auxil.SetPhi(self.Variables["tphi"].Content[index])
        p4_auxil.SetE(sqrt(self.Variables["tptot"].Content[index] ** 2 + m_pi ** 2))
        p4_0 -= p4_auxil
        
        index = self.Variables["PimTrackIndex"].Content
        p4_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
        p4_auxil.SetRho(self.Variables["tptot"].Content[index])
        p4_auxil.SetTheta(self.Variables["tth"].Content[index])
        p4_auxil.SetPhi(self.Variables["tphi"].Content[index])
        p4_auxil.SetE(sqrt(self.Variables["tptot"].Content[index] ** 2 + m_pi ** 2))
        p4_0 -= p4_auxil
        
        self.Variables["PiPiMissMass"].Content = p4_0.M()
        
        
    def calculateDeltaEKKPiPi(self):
        E_0 = 0.

        index = self.Variables["PipTrackIndex"].Content
        E_0 += sqrt(self.Variables["tptot"].Content[index] ** 2 + m_pi ** 2)
        index = self.Variables["PimTrackIndex"].Content
        E_0 += sqrt(self.Variables["tptot"].Content[index] ** 2 + m_pi ** 2)
        index = self.Variables["KpTrackIndex"].Content
        E_0 += sqrt(self.Variables["tptot"].Content[index] ** 2 + m_K ** 2)
        index = self.Variables["KmTrackIndex"].Content
        E_0 += sqrt(self.Variables["tptot"].Content[index] ** 2 + m_K ** 2)

        self.Variables["DeltaEKKPiPi"].Content = E_0 - 2.0 * self.Variables["emeas"].Content


if __name__ == '__main__':
    input_path = '/store11/idpershin/kpkmpippim/prelim_cuts_new/prelim_cut2019_tr_ph_fc_e978_v9.root'
    output_path = 'test_cut.root'
    hists_path = 'test_hists.root'
    log_path = 'test_cut.log'
    
    inter_analysis = IntermediateAnalysis(input_path, analysis_path = hists_path, output_path = output_path, log_path = log_path)
