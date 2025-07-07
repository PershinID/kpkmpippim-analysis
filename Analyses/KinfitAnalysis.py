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
from Containers.Kinfit2K2PiContainer import Kinfit2K2PiContainer
from Containers.Kinfit4PiContainer   import Kinfit4PiContainer

class KinfitAnalysis(Analysis):
    def __init__(self, input_path, kf_2k2pi_path, kf_4pi_path, analysis_path, log_path = None):
        Analysis.__init__(self, analysis_path, logname = "final_cut", logpath = log_path)

        size_variables = { # needed to define variables in self.Variables
            "nt":               Variable("nt",                  "as", max_value = 10),
        }
        self.Variables = {
            ## CMD3 branches
            "emeas":            Variable("emeas",               "f"),
            "demeas":           Variable("demeas",              "f"),
            "xbeam":            Variable("xbeam",               "f"),
            "ybeam":            Variable("ybeam",               "f"),
            "runnum":           Variable("runnum",              "i"),
            "evnum":            Variable("evnum",               "i"),
            "tphi":             Variable("tphi",                "f", sizes = (size_variables["nt"],)),
            "tth":              Variable("tth",                 "f", sizes = (size_variables["nt"],)),
            "tptot":            Variable("tptot",               "f", sizes = (size_variables["nt"],)),
            "trho":             Variable("trho",                "f", sizes = (size_variables["nt"],)),
            "tz":               Variable("tz",                  "f", sizes = (size_variables["nt"],)),
            "tdedx":            Variable("tdedx",               "f", sizes = (size_variables["nt"],)),
            "tcharge":          Variable("tcharge",             "i", sizes = (size_variables["nt"],)),
            "ten":              Variable("ten",                 "f", sizes = (size_variables["nt"],)),
            "terr":             Variable("terr",                "f", sizes = (size_variables["nt"], 3, 3,)),
            "terr0":            Variable("terr0",               "f", sizes = (size_variables["nt"], 6, 6,)),
            "finalstate_id":    Variable("finalstate_id",       "i"),

            "KpKmPipPimLklhd":              Variable("KpKmPipPimLklhd", "f"),
            "KpTrackIndex":                 Variable("KpTrackIndex",    "B"),
            "KmTrackIndex":                 Variable("KmTrackIndex",    "B"),
            "PipTrackIndex":                Variable("PipTrackIndex",   "B"),
            "PimTrackIndex":                Variable("PimTrackIndex",   "B"),

            "KpKmPipPimKinfitIsConverged":   Variable("KpKmPipPimKinfitIsConverged",   "B"),
            "KpKmPipPimKinfitChi2":          Variable("KpKmPipPimKinfitChi2",          "f"),
            "KpKmPipPimKinfitTrackMomenta":  Variable("KpKmPipPimKinfitTrackMomenta",  "f", sizes = (size_variables["nt"],)),
            "KpKmPipPimKinfitTrackThetas":   Variable("KpKmPipPimKinfitTrackThetas",   "f", sizes = (size_variables["nt"],)),
            "KpKmPipPimKinfitTrackPhis":     Variable("KpKmPipPimKinfitTrackPhis",     "f", sizes = (size_variables["nt"],)),
            "KpKmPipPimKinfitTrackEnergies": Variable("KpKmPipPimKinfitTrackEnergies", "f", sizes = (size_variables["nt"],)),
            "KpKmPipPimKinfitTrackIndices":  Variable("KpKmPipPimKinfitTrackIndices",  "B", sizes = (size_variables["nt"],)),

            "PipPimPipPimKinfitIsConverged":   Variable("PipPimPipPimKinfitIsConverged",   "B"),
            "PipPimPipPimKinfitChi2":          Variable("PipPimPipPimKinfitChi2",          "f"),
            "PipPimPipPimKinfitTrackMomenta":  Variable("PipPimPipPimKinfitTrackMomenta",  "f", sizes = (size_variables["nt"],)),
            "PipPimPipPimKinfitTrackThetas":   Variable("PipPimPipPimKinfitTrackThetas",   "f", sizes = (size_variables["nt"],)),
            "PipPimPipPimKinfitTrackPhis":     Variable("PipPimPipPimKinfitTrackPhis",     "f", sizes = (size_variables["nt"],)),
            "PipPimPipPimKinfitTrackEnergies": Variable("PipPimPipPimKinfitTrackEnergies", "f", sizes = (size_variables["nt"],)),
            "PipPimPipPimKinfitTrackIndices":  Variable("PipPimPipPimKinfitTrackIndices",  "B", sizes = (size_variables["nt"],)),

            "PiPiPiPiMissMass2":            TriggerVariable("PiPiPiPiMissMass2",      "f", self.calculatePiPiPiPiMissMass2),
            "PiPiPiMissMass2":              TriggerVariable("PiPiPiMissMass2",        "f", self.calculatePiPiPiMissMass2),
            "KPiPiMissMass2":               TriggerVariable("KPiPiMissMass2",         "f", self.calculateKPiPiMissMass2),
            "KKPiPiMissMass2":              TriggerVariable("KKPiPiMissMass2",        "f", self.calculateKKPiPiMissMass2),
            "KKMissMass":                   TriggerVariable("KKMissMass",             "f", self.calculateKKMissMass),
            "PiPiMissMass":                 TriggerVariable("PiPiMissMass",           "f", self.calculatePiPiMissMass),
            "DeltaEKKPiPi":                 TriggerVariable("DeltaKKPiPi",            "f", self.calculateDeltaEKKPiPi),
            "DeltaE":                       TriggerVariable("DeltaE",                 "f", self.calculateDeltaETotalP),
            "TotalP":                       TriggerVariable("TotalP",                 "f", self.calculateDeltaETotalP),

            "PiPiPiPiMissMass2KF":            TriggerVariable("PiPiPiPiMissMass2KF",      "f", self.calculatePiPiPiPiMissMass2KF),
            "PiPiPiMissMass2KF":              TriggerVariable("PiPiPiMissMass2KF",        "f", self.calculatePiPiPiMissMass2KF),
            "KPiPiMissMass2KF":               TriggerVariable("KPiPiMissMass2KF",         "f", self.calculateKPiPiMissMass2KF),
            "KKPiPiMissMass2KF":              TriggerVariable("KKPiPiMissMass2KF",        "f", self.calculateKKPiPiMissMass2KF),
            "KKMissMassKF":                   TriggerVariable("KKMissMassKF",             "f", self.calculateKKMissMassKF),
            "PiPiMissMassKF":                 TriggerVariable("PiPiMissMassKF",           "f", self.calculatePiPiMissMassKF),
            "DeltaEKKPiPiKF":                 TriggerVariable("DeltaKKPiPiKF",            "f", self.calculateDeltaEKKPiPiKF),
            "DeltaEKF":                       TriggerVariable("DeltaEKF",                 "f", self.calculateDeltaETotalPKF),
            "TotalPKF":                       TriggerVariable("TotalPKF",                 "f", self.calculateDeltaETotalPKF),
        }
        self.Variables.update(size_variables)

        self.InputContainers = [PreliminaryContainer(input_path, "read", self.Variables), Kinfit2K2PiContainer(kf_2k2pi_path, self.Variables), Kinfit4PiContainer(kf_4pi_path, self.Variables),]

        cuts_available = {
            'KpKmPipPimKinfitChi2': lambda: self.Variables["KpKmPipPimKinfitChi2"].Content > 200.,
            'PipPimPipPimKinfitChi2': lambda: self.Variables["PipPimPipPimKinfitChi2"].Content < 1500.,
            'TotalP-DeltaE': lambda: (
                self.Variables["DeltaE"].Content > -150. - 1. * self.Variables["TotalP"].Content
            ),
            'TotalP': lambda: self.Variables["TotalP"].Content > 80.,
            'DeltaEKKPiPi': lambda: abs(self.Variables["DeltaEKKPiPi"].Content) > 80.,
        }
        self.CutDispatcher = CutDispatcher(cuts_available, n_entries_full = self.getEntries())

        histograms_available = {
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
                    'y-axis-nbins': 1500,
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
            'h_KpKmPipPimKinfitChi2': {
                'type': EHists.TH1F,
                'args': {
                    'title': "K^{+}K^{-}#pi^{+}#pi^{-} chi squared",
                    'x-axis-title': "#chi^{2}_{K^{+}K^{-}#pi^{+}#pi^{-}}",
                    'x-axis-nbins': 200000,
                    'x-axis-range': (0., 200000.),
                    'x-variable': self.Variables["KpKmPipPimKinfitChi2"],
                }
            },
            'h_PipPimPipPimKinfitChi2': {
                'type': EHists.TH1F,
                'args': {
                    'title': "#pi^{+}#pi^{-}#pi^{+}#pi^{-} chi squared",
                    'x-axis-title': "#chi^{2}_{#pi^{+}#pi^{-}#pi^{+}#pi^{-}}",
                    'x-axis-nbins': 200000,
                    'x-axis-range': (0., 200000.),
                    'x-variable': self.Variables["PipPimPipPimKinfitChi2"],
                }
            },
            'h_TotalP_DeltaE_KF': {
                'type': EHists.TH2F,
                'args': {
                    'title': "#DeltaE vs total momentum (kinfit)",
                    'x-axis-title': "P_{total}, MeV/c",
                    'x-axis-nbins': 1000,
                    'x-axis-range': (0., 1000.),
                    'x-variable': self.Variables["TotalPKF"],
                    'y-axis-title': "#DeltaE, MeV",
                    'y-axis-nbins': 1500,
                    'y-axis-range': (-1100., 400.),
                    'y-variable': self.Variables["DeltaEKF"],
                },
            },
            'h_TotalP_DeltaEKKPiPi_KF': {
                'type': EHists.TH2F,
                'args': {
                    'title': "#DeltaE_{KK#pi#pi} vs total momentum (kinfit)",
                    'x-axis-title': "P_{total}, MeV/c",
                    'x-axis-nbins': 1000,
                    'x-axis-range': (0., 1000.),
                    'x-variable': self.Variables["TotalPKF"],
                    'y-axis-title': "#DeltaE_{KK#pi#pi}, MeV",
                    'y-axis-nbins': 1500,
                    'y-axis-range': (-700., 800.),
                    'y-variable': self.Variables["DeltaEKKPiPiKF"],
                },
            },
            'h_PiPiPiMissMass2_KF': {
                'type': EHists.TH1F,
                'args': {
                    'title': "#pi#pi#pi missing mass squared (kinfit)",
                    'x-axis-title': "m_{#pi#pi#pi,miss}^{2}, #frac{MeV^{2}}{c^{4}}",
                    'x-axis-nbins': 3000,
                    'x-axis-range': (-1_000_000., 2_000_000.),
                    'x-variable': self.Variables["PiPiPiMissMass2KF"],
                },
            },
            'h_KPiPiMissMass2_KF': {
                'type': EHists.TH1F,
                'args': {
                    'title': "K#pi#pi missing mass squared (kinfit)",
                    'x-axis-title': "m_{K#pi#pi,miss}^{2}, #frac{MeV^{2}}{c^{4}}",
                    'x-axis-nbins': 2500,
                    'x-axis-range': (-1_000_000., 1_500_000.),
                    'x-variable': self.Variables["KPiPiMissMass2KF"],
                },
            },
            'h_PiPiPiPiMissMass2_KF': {
                'type': EHists.TH1F,
                'args': {
                    'title': "#pi^{+}#pi^{-}#pi^{+}#pi^{-} missing mass squared (kinfit)",
                    'x-axis-title': "m_{#pi+#pi-#pi+#pi-,miss}^{2}, #frac{MeV^{2}}{c^{4}}",
                    'x-axis-nbins': 1000,
                    'x-axis-range': (-200_000., 800_000.),
                    'x-variable': self.Variables["PiPiPiPiMissMass2KF"],
                },
            },
            'h_KKPiPiMissMass2_KF': {
                'type': EHists.TH1F,
                'args': {
                    'title': "K^{+}K^{-}#pi^{+}#pi^{-} missing mass squared (kinfit)",
                    'x-axis-title': "m_{K^{+}K^{-}#pi^{+}#pi^{-},miss}, #frac{MeV^{2}}{c^{4}}",
                    'x-axis-nbins': 8000,
                    'x-axis-range': (-400_000., 400_000.),
                    'x-variable': self.Variables["KKPiPiMissMass2KF"],
                },
            },
            'h_KKMissMass_KF': {
                'type': EHists.TH1F,
                'args': {
                    'title': "K^{+}K^{-} missing mass squared (kinfit)",
                    'x-axis-title': "m_{K^{+}K^{-},miss}, #frac{MeV^{2}}{c^{4}}",
                    'x-axis-nbins': 1000,
                    'x-axis-range': (0., 1000.),
                    'x-variable': self.Variables["KKMissMassKF"],
                }
            },
            'h_PiPiMissMass_KF': {
                'type': EHists.TH1F,
                'args': {
                    'title': "#pi^{+}#pi^{-} missing mass squared (kinfit)",
                    'x-axis-title': "m_{#pi^{+}#pi^{-},miss}, #frac{MeV^{2}}{c^{4}}",
                    'x-axis-nbins': 300000,
                    'x-axis-range': (0., 3000.),
                    'x-variable': self.Variables["PiPiMissMassKF"],
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


    def calculateDeltaETotalP(self):
        p4_0     = TLorentzVector(0.0, 0.0, 0.0, 0.0)
        p4_auxil = TLorentzVector(0.0, 0.0, 0.0, 0.0)
        
        for (p, th, phi) in zip(self.Variables["tptot"].Content,
                                self.Variables["tth"].Content,
                                self.Variables["tphi"].Content):
            p4_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
            p4_auxil.SetRho(p)
            p4_auxil.SetTheta(th)
            p4_auxil.SetPhi(phi)
            p4_auxil.SetE(sqrt(p * p + m_pi * m_pi))
            p4_0 += p4_auxil

        self.Variables["DeltaE"].Content = p4_0.E() - 2.0 * self.Variables["emeas"].Content
        self.Variables["TotalP"].Content = p4_0.Rho()


    def calculatePiPiPiPiMissMass2KF(self):
        p4_0     = TLorentzVector(0.0, 0.0, 0.0, 2.0 * self.Variables["emeas"].Content)
        p4_auxil = TLorentzVector(0.0, 0.0, 0.0, 0.0)
        
        for (p, th, phi) in zip(self.Variables["KpKmPipPimKinfitTrackMomenta"].Content,
                                self.Variables["KpKmPipPimKinfitTrackThetas"].Content,
                                self.Variables["KpKmPipPimKinfitTrackPhis"].Content):
            p4_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
            p4_auxil.SetRho(p)
            p4_auxil.SetTheta(th)
            p4_auxil.SetPhi(phi)
            p4_auxil.SetE(sqrt(p ** 2 + m_pi ** 2))
            p4_0 -= p4_auxil
        
        self.Variables["PiPiPiPiMissMass2KF"].Content = p4_0.M2()


    def calculatePiPiPiMissMass2KF(self):
        p4_0     = TLorentzVector(0.0, 0.0, 0.0, 2.0 * self.Variables["emeas"].Content)
        p4_auxil = TLorentzVector(0.0, 0.0, 0.0, 0.0)

        for i, (p, th, phi) in enumerate(zip(self.Variables["KpKmPipPimKinfitTrackMomenta"].Content,
                                             self.Variables["KpKmPipPimKinfitTrackThetas"].Content,
                                             self.Variables["KpKmPipPimKinfitTrackPhis"].Content)):
            if i >= 3: break
            p4_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
            p4_auxil.SetRho(p)
            p4_auxil.SetTheta(th)
            p4_auxil.SetPhi(phi)
            p4_auxil.SetE(sqrt(p ** 2 + m_pi ** 2))
            p4_0 -= p4_auxil

        self.Variables["PiPiPiMissMass2KF"].Content = p4_0.M2()


    def calculateKPiPiMissMass2KF(self):
        p4_0     = TLorentzVector(0.0, 0.0, 0.0, 2.0 * self.Variables["emeas"].Content)
        p4_auxil = TLorentzVector(0.0, 0.0, 0.0, 0.0)

        for i, (p, th, phi) in enumerate(zip(self.Variables["KpKmPipPimKinfitTrackMomenta"].Content,
                                             self.Variables["KpKmPipPimKinfitTrackThetas"].Content,
                                             self.Variables["KpKmPipPimKinfitTrackPhis"].Content)):
            mass = m_pi
            if i == 0: mass = m_K
            if i == 1: continue
            
            p4_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
            p4_auxil.SetRho(p)
            p4_auxil.SetTheta(th)
            p4_auxil.SetPhi(phi)
            p4_auxil.SetE(sqrt(p ** 2 + mass ** 2))
            p4_0 -= p4_auxil
        
        self.Variables["KPiPiMissMass2KF"].Content = p4_0.M2()


    def calculateKKPiPiMissMass2KF(self):
        p4_0     = TLorentzVector(0.0, 0.0, 0.0, 2.0 * self.Variables["emeas"].Content)
        p4_auxil = TLorentzVector(0.0, 0.0, 0.0, 0.0)

        for i, (p, th, phi) in enumerate(zip(self.Variables["KpKmPipPimKinfitTrackMomenta"].Content,
                                             self.Variables["KpKmPipPimKinfitTrackThetas"].Content,
                                             self.Variables["KpKmPipPimKinfitTrackPhis"].Content)):
            mass = m_pi
            if i in (0, 1): mass = m_K
            
            p4_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
            p4_auxil.SetRho(p)
            p4_auxil.SetTheta(th)
            p4_auxil.SetPhi(phi)
            p4_auxil.SetE(sqrt(p ** 2 + mass ** 2))
            p4_0 -= p4_auxil
        
        self.Variables["KKPiPiMissMass2KF"].Content = p4_0.M2()


    def calculateKKMissMassKF(self):
        p4_0     = TLorentzVector(0.0, 0.0, 0.0, 2.0 * self.Variables["emeas"].Content)
        p4_auxil = TLorentzVector(0.0, 0.0, 0.0, 0.0)

        for i, (p, th, phi) in enumerate(zip(self.Variables["KpKmPipPimKinfitTrackMomenta"].Content,
                                             self.Variables["KpKmPipPimKinfitTrackThetas"].Content,
                                             self.Variables["KpKmPipPimKinfitTrackPhis"].Content)):
            if i in (2, 3): continue
            
            p4_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
            p4_auxil.SetRho(p)
            p4_auxil.SetTheta(th)
            p4_auxil.SetPhi(phi)
            p4_auxil.SetE(sqrt(p ** 2 + m_K ** 2))
            p4_0 -= p4_auxil
        
        self.Variables["KKMissMassKF"].Content = p4_0.M()


    def calculatePiPiMissMassKF(self):
        p4_0     = TLorentzVector(0.0, 0.0, 0.0, 2.0 * self.Variables["emeas"].Content)
        p4_auxil = TLorentzVector(0.0, 0.0, 0.0, 0.0)

        for i, (p, th, phi) in enumerate(zip(self.Variables["KpKmPipPimKinfitTrackMomenta"].Content,
                                             self.Variables["KpKmPipPimKinfitTrackThetas"].Content,
                                             self.Variables["KpKmPipPimKinfitTrackPhis"].Content)):
            if i in (0, 1): continue
            
            p4_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
            p4_auxil.SetRho(p)
            p4_auxil.SetTheta(th)
            p4_auxil.SetPhi(phi)
            p4_auxil.SetE(sqrt(p ** 2 + m_pi ** 2))
            p4_0 -= p4_auxil
        
        self.Variables["PiPiMissMassKF"].Content = p4_0.M()
        
        
    def calculateDeltaEKKPiPiKF(self):
        E_0 = 0.

        for i, p in enumerate(self.Variables["KpKmPipPimKinfitTrackMomenta"].Content):
            mass = m_pi
            if i in (0, 1): mass = m_K
            
            E_0 += sqrt(p ** 2 + mass ** 2)

        self.Variables["DeltaEKKPiPiKF"].Content = E_0 - 2.0 * self.Variables["emeas"].Content


    def calculateDeltaETotalPKF(self):
        p4_0     = TLorentzVector(0.0, 0.0, 0.0, 0.0)
        p4_auxil = TLorentzVector(0.0, 0.0, 0.0, 0.0)

        for (p, th, phi) in zip(self.Variables["KpKmPipPimKinfitTrackMomenta"].Content,
                                self.Variables["KpKmPipPimKinfitTrackThetas"].Content,
                                self.Variables["KpKmPipPimKinfitTrackPhis"].Content):
            p4_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
            p4_auxil.SetRho(p)
            p4_auxil.SetTheta(th)
            p4_auxil.SetPhi(phi)
            p4_auxil.SetE(sqrt(p ** 2 + m_pi ** 2))
            p4_0 += p4_auxil

        self.Variables["DeltaEKF"].Content = p4_0.E() - 2.0 * self.Variables["emeas"].Content
        self.Variables["TotalPKF"].Content = p4_0.Rho()
