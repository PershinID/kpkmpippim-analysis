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
from Containers.Kinfit4PiContainer import Kinfit4PiContainer
from Containers.FinalContainer import FinalContainer

class FinalAnalysis(Analysis):
    def __init__(self, input_path, kinfit_2k2pi_path, kinfit_4pi_path, *, analysis_path = None, output_path = None, log_path = None):
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
            
            "PipPimPipPimMissMass":         TriggerVariable("PipPimPipPimMissMass",  "f", self.calculatePipPimPipPimMissMass),
            "PiPiPiMissMass":               TriggerVariable("PiPiPiMissMass",        "f", self.calculatePiPiPiMissMass),
            "KPiPiMissMass":                TriggerVariable("KPiPiMissMass",         "f", self.calculateKPiPiMissMass),

            "DeltaEKKPiPi":                 TriggerVariable("DeltaKKPiPi",           "f", self.calculateDeltaEKKPiPi),
            
            "PipPimPipPimKinfitChi2":       Variable("PipPimPipPimKinfitChi2",       "f"),
            "PipPimPipPimKinfitPip0Track":  Variable("PipPimPipPimKinfitPip0Track",  "f", sizes = (4,)),
            "PipPimPipPimKinfitPim0Track":  Variable("PipPimPipPimKinfitPim0Track",  "f", sizes = (4,)),
            "PipPimPipPimKinfitPip1Track":  Variable("PipPimPipPimKinfitPip1Track",  "f", sizes = (4,)),
            "PipPimPipPimKinfitPim1Track":  Variable("PipPimPipPimKinfitPim1Track",  "f", sizes = (4,)),
            
            "KpKmPipPimKinfitChi2":         Variable("KpKmPipPimKinfitChi2",           "f"),
            "KpKmPipPimKinfitKpTrack":      Variable("KpKmPipPimKinfitKpTrack",        "f", sizes = (4,)),
            "KpKmPipPimKinfitKmTrack":      Variable("KpKmPipPimKinfitKmTrack",        "f", sizes = (4,)),
            "KpKmPipPimKinfitPipTrack":     Variable("KpKmPipPimKinfitPipTrack",       "f", sizes = (4,)),
            "KpKmPipPimKinfitPimTrack":     Variable("KpKmPipPimKinfitPimTrack",       "f", sizes = (4,)),
            "DeltaEKF":                     TriggerVariable("DeltaEKF",                "f", self.calculateDeltaETotalPKF),
            "TotalPKF":                     TriggerVariable("TotalPKF",                "f", self.calculateDeltaETotalPKF),
            "PipPimPipPimMissMassKF":       TriggerVariable("PipPimPipPimMissMassKF",  "f", self.calculatePipPimPipPimMissMassKF),
            "PiPiPiMissMassKF":             TriggerVariable("PiPiPiMissMassKF",        "f", self.calculatePiPiPiMissMassKF),
            "KPiPiMissMassKF":              TriggerVariable("KPiPiMissMassKF",         "f", self.calculateKPiPiMissMassKF),
        }
        self.Variables.update(size_variables)

        self.InputContainers = [
            PreliminaryContainer(input_path, 'read', self.Variables),
            Kinfit2K2PiContainer(kinfit_2k2pi_path, self.Variables),
            Kinfit4PiContainer(kinfit_4pi_path, self.Variables),
        ]
        self.OutputContainers = [FinalContainer(output_path, "recreate", self.Variables)]

        cuts_available = {
            'KpKmPipPimLklhd': lambda: (
                self.Variables["KpKmPipPimLklhd"].Content < -3.
            ),
            'TotalP-DeltaE': lambda: (
                self.Variables["DeltaE"].Content > -50. - 1. * self.Variables["TotalP"].Content
            ),
            'TotalP': lambda: self.Variables["TotalP"].Content > 80.,
            'DeltaEKKPiPi': lambda: abs(self.Variables["DeltaEKKPiPi"].Content) > 80.,
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
            'h_trho': {
                'type': EHists.TH1F,
                'args': {
                    'title': "Track #rho coordinate",
                    'x-axis-title': "#rho_{track}, cm",
                    'x-axis-nbins': 400,
                    'x-axis-range': (-20.0, 20.0),
                    'x-variable': self.Variables["trho"],
                },
            },
            'h_tz': {
                'type': EHists.TH1F,
                'args': {
                    'title': "Track Z",
                    'x-axis-title': "z_{track}, cm",
                    'x-axis-nbins': 400,
                    'x-axis-range': (-20.0, 20.0),
                    'x-variable': self.Variables["tz"],
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
            'h_KpKmPipPimLklhd': {
                'type': EHists.TH1F,
                'args': {
                    'title': "K+K-#pi+#pi- hypothesis likelihood by ionization losses",
                    'x-axis-title': "Likelihood",
                    'x-axis-nbins': 100,
                    'x-axis-range': (-10.0, 0.0),
                    'x-variable': self.Variables["KpKmPipPimLklhd"],
                },
            },
            'h_TotalP_DeltaE': {
                'type': EHists.TH2F,
                'args': {
                    'title': "Total momentum vs #DeltaE",
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
            'h_TotalP_DeltaEKKPiPi': {
                'type': EHists.TH2F,
                'args': {
                    'title': "Total momentum vs #DeltaE_{KK#pi#pi}",
                    'x-axis-title': "P_{total}, MeV/c",
                    'x-axis-nbins': 1000,
                    'x-axis-range': (0.0, 1000.0),
                    'x-variable': self.Variables["TotalP"],
                    'y-axis-title': "#DeltaE_{KK#pi#pi}, MeV",
                    'y-axis-nbins': 1600,
                    'y-axis-range': (-800.0, 800.0),
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
                    'x-variable': self.Variables["KpKmPipPimLklhd"],
                },
            },
            'h_phen': {
                'type': EHists.TH1F,
                'args': {
                    'title': "Photon energy",
                    'x-axis-title': "E_{photons}, MeV",
                    'x-axis-nbins': 2000,
                    'x-axis-range': (0., 2000.),
                    'x-variable': self.Variables["phen"],
                },
            },
            'h_PiPiPiMissMass': {
                'type': EHists.TH1F,
                'args': {
                    'title': "#pi#pi#pi missing mass squared",
                    'x-axis-title': "m_{#pi#pi#pi,miss}, MeV^{2}/c^{4}",
                    'x-axis-nbins': 300,
                    'x-axis-range': (-3000000., 3000000.),
                    'x-variable': self.Variables["PiPiPiMissMass"],
                },
            },
            'h_KPiPiMissMass': {
                'type': EHists.TH1F,
                'args': {
                    'title': "K#pi#pi missing mass squared",
                    'x-axis-title': "m_{K#pi#pi,miss}, MeV^{2}/c^{4}",
                    'x-axis-nbins': 300,
                    'x-axis-range': (-3000000., 3000000.),
                    'x-variable': self.Variables["KPiPiMissMass"],
                },
            },
            'h_PipPimPipPimMissMass': {
                'type': EHists.TH1F,
                'args': {
                    'title': "#pi+#pi-#pi+#pi- missing mass",
                    'x-axis-title': "M_{#pi+#pi-#pi+#pi-,miss}",
                    'x-axis-nbins': 300,
                    'x-axis-range': (0.0, 3000000.0),
                    'x-variable': self.Variables["PipPimPipPimMissMass"],
                },
            },
            'h_PipPimPipPimKinfitChi2': {
                'type': EHists.TH1F,
                'args': {
                    'title': "#pi+#pi-#pi+#pi- kinfit hypothesis #chi^{2}",
                    'x-axis-title': "#chi^{2}",
                    'x-axis-nbins': 200,
                    'x-axis-range': (0.0, 10000.0),
                    'x-variable': self.Variables["PipPimPipPimKinfitChi2"],
                },
            },
            'h_KpKmPipPimKinfitChi2': {
                'type': EHists.TH1F,
                'args': {
                    'title': "K+K-#pi+#pi- kinfit hypothesis #chi^{2}",
                    'x-axis-title': "#chi^{2}",
                    'x-axis-nbins': 200,
                    'x-axis-range': (0.0, 10000.0),
                    'x-variable': self.Variables["KpKmPipPimKinfitChi2"],
                },
            },
            'h_PiPiPiMissMass_KF': {
                'type': EHists.TH1F,
                'args': {
                    'title': "#pi#pi#pi missing mass squared with kinfit",
                    'x-axis-title': "m_{#pi#pi#pi,miss,KF}, MeV^{2}/c^{4}",
                    'x-axis-nbins': 300,
                    'x-axis-range': (-3000000., 3000000.),
                    'x-variable': self.Variables["PiPiPiMissMassKF"],
                },
            },
            'h_KPiPiMissMass_KF': {
                'type': EHists.TH1F,
                'args': {
                    'title': "K#pi#pi missing mass squared with kinfit",
                    'x-axis-title': "m_{K#pi#pi,miss,KF}, MeV^{2}/c^{4}",
                    'x-axis-nbins': 300,
                    'x-axis-range': (-3000., 3000.),
                    'x-variable': self.Variables["KPiPiMissMassKF"],
                },
            },
            'h_PipPimPipPimMissMass_KF': {
                'type': EHists.TH1F,
                'args': {
                    'title': "#pi+#pi-#pi+#pi- missing mass squared with kinfit",
                    'x-axis-title': "M_{#pi+#pi-#pi+#pi-,miss,KF}",
                    'x-axis-nbins': 300,
                    'x-axis-range': (0.0, 3000000.0),
                    'x-variable': self.Variables["PipPimPipPimMissMassKF"],
                },
            },
            'h_TotalP_DeltaE_KF': {
                'type': EHists.TH2F,
                'args': {
                    'title': "Total momentum vs #DeltaE with kinfit",
                    'x-axis-title': "P_{total,KF}, MeV/c",
                    'x-axis-nbins': 1000,
                    'x-axis-range': (0.0, 1000.0),
                    'x-variable': self.Variables["TotalPKF"],
                    'y-axis-title': "#DeltaE_{KF}, MeV",
                    'y-axis-nbins': 2000,
                    'y-axis-range': (-2000.0, 2000.0),
                    'y-variable': self.Variables["DeltaEKF"],
                },
            },
            'h_finalstate_id': {
                'type': EHists.TH1I,
                'args': {
                    'title': "Final state ID",
                    'x-axis-title': "ID_{final state}",
                    'x-axis-nbins': 30,
                    'x-axis-range': (0, 30),
                    'x-variable': self.Variables["finalstate_id"],
                },
            },
        }
        self.HistogramDispatcher = HistogramDispatcher(histograms_available)

        
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
        
        self.Variables["PipPimPipPimMissMass"].Content = p4_0.M2()


    def calculatePiPiPiMissMass(self):
        tp_lorentz = TLorentzVector(0.0, 0.0, 0.0, 2.0 * self.Variables["emeas"].Content)
        tp_auxil = TLorentzVector(0.0, 0.0, 0.0, 0.0)

        for i, (p, th, phi) in enumerate(zip(self.Variables["tptot"].Content,
                                             self.Variables["tth"].Content,
                                             self.Variables["tphi"].Content)):
            if i > 2: break
            tp_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
            tp_auxil.SetRho(p)
            tp_auxil.SetTheta(th)
            tp_auxil.SetPhi(phi)
            tp_auxil.SetE(sqrt(p * p + m_pi * m_pi))
            tp_lorentz -= tp_auxil

        self.Variables["PiPiPiMissMass"].Content = tp_lorentz.M2()


    def calculateKPiPiMissMass(self):
        tp_lorentz = TLorentzVector(0.0, 0.0, 0.0, 2.0 * self.Variables["emeas"].Content)
        tp_auxil = TLorentzVector(0.0, 0.0, 0.0, 0.0)

        index = self.Variables["PipTrackIndex"].Content
        tp_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
        tp_auxil.SetRho(self.Variables["tptot"].Content[index])
        tp_auxil.SetTheta(self.Variables["tth"].Content[index])
        tp_auxil.SetPhi(self.Variables["tphi"].Content[index])
        tp_auxil.SetE(sqrt(self.Variables["tptot"].Content[index] ** 2 + m_pi * m_pi))
        tp_lorentz -= tp_auxil

        index = self.Variables["PimTrackIndex"].Content
        tp_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
        tp_auxil.SetRho(self.Variables["tptot"].Content[index])
        tp_auxil.SetTheta(self.Variables["tth"].Content[index])
        tp_auxil.SetPhi(self.Variables["tphi"].Content[index])
        tp_auxil.SetE(sqrt(self.Variables["tptot"].Content[index] ** 2 + m_pi * m_pi))
        tp_lorentz -= tp_auxil

        index = self.Variables["KpTrackIndex"].Content
        tp_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
        tp_auxil.SetRho(self.Variables["tptot"].Content[index])
        tp_auxil.SetTheta(self.Variables["tth"].Content[index])
        tp_auxil.SetPhi(self.Variables["tphi"].Content[index])
        tp_auxil.SetE(sqrt(self.Variables["tptot"].Content[index] ** 2 + m_K * m_K))
        tp_lorentz -= tp_auxil
        
        self.Variables["KPiPiMissMass"].Content = tp_lorentz.M2()


    def calculateDeltaEKKPiPi(self):
        tp_lorentz = TLorentzVector(0.0, 0.0, 0.0, 0.0)
        tp_auxil = TLorentzVector(0.0, 0.0, 0.0, 0.0)

        index = self.Variables["PipTrackIndex"].Content
        tp_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
        tp_auxil.SetRho(self.Variables["tptot"].Content[index])
        tp_auxil.SetTheta(self.Variables["tth"].Content[index])
        tp_auxil.SetPhi(self.Variables["tphi"].Content[index])
        tp_auxil.SetE(sqrt(self.Variables["tptot"].Content[index] ** 2 + m_pi * m_pi))
        tp_lorentz += tp_auxil

        index = self.Variables["PimTrackIndex"].Content
        tp_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
        tp_auxil.SetRho(self.Variables["tptot"].Content[index])
        tp_auxil.SetTheta(self.Variables["tth"].Content[index])
        tp_auxil.SetPhi(self.Variables["tphi"].Content[index])
        tp_auxil.SetE(sqrt(self.Variables["tptot"].Content[index] ** 2 + m_pi * m_pi))
        tp_lorentz += tp_auxil

        index = self.Variables["KpTrackIndex"].Content
        tp_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
        tp_auxil.SetRho(self.Variables["tptot"].Content[index])
        tp_auxil.SetTheta(self.Variables["tth"].Content[index])
        tp_auxil.SetPhi(self.Variables["tphi"].Content[index])
        tp_auxil.SetE(sqrt(self.Variables["tptot"].Content[index] ** 2 + m_K * m_K))
        tp_lorentz += tp_auxil
        
        index = self.Variables["KmTrackIndex"].Content
        tp_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
        tp_auxil.SetRho(self.Variables["tptot"].Content[index])
        tp_auxil.SetTheta(self.Variables["tth"].Content[index])
        tp_auxil.SetPhi(self.Variables["tphi"].Content[index])
        tp_auxil.SetE(sqrt(self.Variables["tptot"].Content[index] ** 2 + m_K * m_K))
        tp_lorentz += tp_auxil

        self.Variables["DeltaE"].Content = tp_lorentz.E() - 2.0 * self.Variables["emeas"].Content
        

    def calculateDeltaETotalPKF(self):
        tp_lorentz = TLorentzVector(0.0, 0.0, 0.0, 0.0)
        tp_auxil = TLorentzVector(0.0, 0.0, 0.0, 0.0)
        
        tp_auxil.SetPxPyPzE(*self.Variables["KpKmPipPimKinfitKpTrack"].Content)
        tp_auxil.SetE(sqrt(tp_auxil.P() ** 2 + m_pi * m_pi))
        tp_lorentz += tp_auxil
        tp_auxil.SetPxPyPzE(*self.Variables["KpKmPipPimKinfitKmTrack"].Content)
        tp_auxil.SetE(sqrt(tp_auxil.P() ** 2 + m_pi * m_pi))
        tp_lorentz += tp_auxil
        tp_auxil.SetPxPyPzE(*self.Variables["KpKmPipPimKinfitPipTrack"].Content)
        tp_auxil.SetE(sqrt(tp_auxil.P() ** 2 + m_pi * m_pi))
        tp_lorentz += tp_auxil
        tp_auxil.SetPxPyPzE(*self.Variables["KpKmPipPimKinfitPipTrack"].Content)
        tp_auxil.SetE(sqrt(tp_auxil.P() ** 2 + m_pi * m_pi))
        tp_lorentz += tp_auxil

        self.Variables["DeltaEKF"].Content = tp_lorentz.E() - 2.0 * self.Variables["emeas"].Content
        self.Variables["TotalPKF"].Content = tp_lorentz.Rho()


    def calculatePipPimPipPimMissMassKF(self):
        tp_lorentz = TLorentzVector(0.0, 0.0, 0.0, 2. * self.Variables["emeas"].Content)
        tp_auxil = TLorentzVector(0.0, 0.0, 0.0, 0.0)

        tp_auxil.SetPxPyPzE(*self.Variables["KpKmPipPimKinfitKpTrack"].Content)
        tp_auxil.SetE(sqrt(tp_auxil.P() ** 2 + m_pi * m_pi))
        tp_lorentz -= tp_auxil
        tp_auxil.SetXYZT(*self.Variables["KpKmPipPimKinfitKmTrack"].Content)
        tp_auxil.SetE(sqrt(tp_auxil.P() ** 2 + m_pi * m_pi))
        tp_lorentz -= tp_auxil
        tp_auxil.SetXYZT(*self.Variables["KpKmPipPimKinfitPipTrack"].Content)
        tp_auxil.SetE(sqrt(tp_auxil.P() ** 2 + m_pi * m_pi))
        tp_lorentz -= tp_auxil
        tp_auxil.SetXYZT(*self.Variables["KpKmPipPimKinfitPipTrack"].Content)
        tp_auxil.SetE(sqrt(tp_auxil.P() ** 2 + m_pi * m_pi))
        tp_lorentz -= tp_auxil
        
        self.Variables["PipPimPipPimMissMassKF"].Content = tp_lorentz.M2()

        
    def calculatePiPiPiMissMassKF(self):
        tp_lorentz = TLorentzVector(0.0, 0.0, 0.0, 2. * self.Variables["emeas"].Content)
        tp_auxil = TLorentzVector(0.0, 0.0, 0.0, 0.0)

        tp_auxil.SetPxPyPzE(*self.Variables["KpKmPipPimKinfitKpTrack"].Content)
        tp_auxil.SetE(sqrt(tp_auxil.P() ** 2 + m_pi * m_pi))
        tp_lorentz -= tp_auxil
        tp_auxil.SetXYZT(*self.Variables["KpKmPipPimKinfitPimTrack"].Content)
        tp_auxil.SetE(sqrt(tp_auxil.P() ** 2 + m_pi * m_pi))
        tp_lorentz -= tp_auxil
        tp_auxil.SetXYZT(*self.Variables["KpKmPipPimKinfitPipTrack"].Content)
        tp_auxil.SetE(sqrt(tp_auxil.P() ** 2 + m_pi * m_pi))
        tp_lorentz -= tp_auxil
        
        self.Variables["PiPiPiMissMassKF"].Content = tp_lorentz.M2()


    def calculateKPiPiMissMassKF(self):
        tp_lorentz = TLorentzVector(0.0, 0.0, 0.0, 2. * self.Variables["emeas"].Content)
        tp_auxil = TLorentzVector(0.0, 0.0, 0.0, 0.0)

        tp_auxil.SetPxPyPzE(*self.Variables["KpKmPipPimKinfitKpTrack"].Content)
        tp_auxil.SetE(sqrt(tp_auxil.P() ** 2 + m_K * m_K))
        tp_lorentz -= tp_auxil
        tp_auxil.SetXYZT(*self.Variables["KpKmPipPimKinfitPipTrack"].Content)
        tp_auxil.SetE(sqrt(tp_auxil.P() ** 2 + m_pi * m_pi))
        tp_lorentz -= tp_auxil
        tp_auxil.SetXYZT(*self.Variables["KpKmPipPimKinfitPimTrack"].Content)
        tp_auxil.SetE(sqrt(tp_auxil.P() ** 2 + m_pi * m_pi))
        tp_lorentz -= tp_auxil
        
        self.Variables["KPiPiMissMassKF"].Content = tp_lorentz.M2()


    def calculateEntry(self):
        pass
        
        
if __name__ == '__main__':
    ##Parsing input arguments
    parser = ArgumentParser()
    parser.add_argument('--version', default = 'v9', choices = ['v9'], help = 'Version of CMD-3 data tree')
    parser.add_argument('--year', choices = ['2019', '2020', '2021', '2022', '2023'], required = True)
    parser.add_argument('--energy', required = True)
    parser.add_argument('--input-dir', required = True)
    parser.add_argument('--output-dir')
    parser.add_argument('--hists-dir')
    parser.add_argument('--is-multihad', action = 'store_true')
    parser.add_argument('--is-sim', action = 'store_true')
    args = parser.parse_args()

    if args.is_sim and args.is_multihad:
        raise ValueError("Is_sim and is_multihad flags exist at the same time")
    if not os.path.exists(args.input_dir): raise OSError(f"Input directory does not exist: {args.input_dir}")
    if not os.path.exists(args.output_dir): raise OSError(f"Output directory does not exist: {args.output_dir}")
    if not os.path.exists(args.hists_dir): raise OSError(f"Histograms directory does not exist: {args.hists_dir}")

    prefix = ''
    if args.is_sim:
        prefix = 'sim_'
    if args.is_multihad:
        prefix = 'multihad_'
        
    input_path = f"{args.input_dir}/{prefix}prelim_cut{args.year}_tr_ph_fc_e{args.energy}_{args.version}.root"
    kinfit_2k2pi_path = f"{args.input_dir}/{prefix}kinfit_2k2pi{args.year}_tr_ph_fc_e{args.energy}_{args.version}.root"
    kinfit_4pi_path = f"{args.input_dir}/{prefix}kinfit_4pi{args.year}_tr_ph_fc_e{args.energy}_{args.version}.root"
    output_path = f"{args.output_dir}/{prefix}final_cut{args.year}_tr_ph_fc_e{args.energy}_{args.version}.root" if args.output_dir else None
    hists_path = f"{args.hists_dir}/{prefix}final_hists{args.year}_tr_ph_fc_e{args.energy}_{args.version}.root" if args.hists_dir else None
    log_path = f"{args.output_dir}/{prefix}final_cut_y{args.year}_e{args.energy}_{args.version}.log"

    if not os.path.exists(input_path): raise OSError(f"Input path does not exist: {input_path}")
    if not os.path.exists(kinfit_2k2pi_path): raise OSError(f"2K2Pi kinfit path does not exist: {kinfit_2k2pi_path}")
    if not os.path.exists(kinfit_4pi_path): raise OSError(f"4Pi kinfit path does not exist: {kinfit_4pi_path}")
    
    ## Preliminary, Kinfit --> Final
    final_analysis = FinalKpKmPipPimAnalysis(
        input_path, kinfit_2k2pi_path, kinfit_4pi_path,
        analysis_path = hists_path, output_path = output_path, log_path = log_path,
    )
    final_analysis.addHistogram('h_TotalP_DeltaE')
    final_analysis.addHistogram('h_KpKmPipPimLklhd')
    final_analysis.addHistogram('h_KpKmPipPimKinfitChi2')
    final_analysis.addHistogram('h_PipPimPipPimKinfitChi2')
    final_analysis.addHistogram('h_PipPimPipPimMissMass')
    final_analysis.addHistogram('h_PiPiPiMissMass')
    final_analysis.addHistogram('h_KPiPiMissMass')
    final_analysis.addHistogram('h_TotalP_DeltaE_KF')
    final_analysis.addHistogram('h_PipPimPipPimMissMass_KF')
    final_analysis.addHistogram('h_PiPiPiMissMass_KF')
    final_analysis.addHistogram('h_KPiPiMissMass_KF')
    if args.is_multihad: final_analysis.addHistogram('h_finalstate_id')
    final_analysis.loop()

    final_analysis.addCut('KpKmPipPimLklhd')
    final_analysis.addHistogram('h_TotalP_DeltaE')
    final_analysis.addHistogram('h_KpKmPipPimLklhd')
    final_analysis.addHistogram('h_KpKmPipPimKinfitChi2')
    final_analysis.addHistogram('h_PipPimPipPimKinfitChi2')
    final_analysis.addHistogram('h_PipPimPipPimMissMass')
    final_analysis.addHistogram('h_PiPiPiMissMass')
    final_analysis.addHistogram('h_KPiPiMissMass')
    final_analysis.addHistogram('h_TotalP_DeltaE_KF')
    final_analysis.addHistogram('h_PipPimPipPimMissMass_KF')
    final_analysis.addHistogram('h_PiPiPiMissMass_KF')
    final_analysis.addHistogram('h_KPiPiMissMass_KF')
    if args.is_multihad: final_analysis.addHistogram('h_finalstate_id')
    final_analysis.loop()
    
    final_analysis.addCut('TotalP-DeltaE')
    final_analysis.addHistogram('h_TotalP_DeltaE')
    final_analysis.addHistogram('h_KpKmPipPimLklhd')
    final_analysis.addHistogram('h_KpKmPipPimKinfitChi2')
    final_analysis.addHistogram('h_PipPimPipPimKinfitChi2')
    final_analysis.addHistogram('h_PipPimPipPimMissMass')
    final_analysis.addHistogram('h_PiPiPiMissMass')
    final_analysis.addHistogram('h_KPiPiMissMass')
    final_analysis.addHistogram('h_TotalP_DeltaE_KF')
    final_analysis.addHistogram('h_PipPimPipPimMissMass_KF')
    final_analysis.addHistogram('h_PiPiPiMissMass_KF')
    final_analysis.addHistogram('h_KPiPiMissMass_KF')
    if args.is_multihad: final_analysis.addHistogram('h_finalstate_id')
    final_analysis.loop()
    
    final_analysis.dumpToFile()
