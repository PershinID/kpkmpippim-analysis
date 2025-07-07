from math import pi, sqrt, inf, acos
from functools import partialmethod

from argparse import ArgumentParser
from datetime import date
import os

import ROOT
from ROOT import TLorentzVector, TVector3

from Base.Variable import Variable, TriggerVariable
from Base.PhysicalConstants import m_pi, m_K
from Base.Analysis import Analysis, CutDispatcher, HistogramDispatcher, EHists

from Containers.FinalContainer import FinalContainer

class DynamicsAnalysis(Analysis):
    def __init__(self, input_path, analysis_path, *, log_path = None):
        Analysis.__init__(self, analysis_path, logname = "dynamics", logpath = log_path)

        self.Variables = {
            "emeas":                        Variable("emeas",               "f"),
            "demeas":                       Variable("demeas",              "f"),
            "xbeam":                        Variable("xbeam",               "f"),
            "ybeam":                        Variable("ybeam",               "f"),
            "runnum":                       Variable("runnum",              "i"),
            "evnum":                        Variable("evnum",               "i"),
            "tphi":                         Variable("tphi",                "f", sizes = (4,)),
            "tth":                          Variable("tth",                 "f", sizes = (4,)),
            "tptot":                        Variable("tptot",               "f", sizes = (4,)),
            "KpTrackIndex":                 Variable("KpTrackIndex",        "B"),
            "KmTrackIndex":                 Variable("KmTrackIndex",        "B"),
            "PipTrackIndex":                Variable("PipTrackIndex",       "B"),
            "PimTrackIndex":                Variable("PimTrackIndex",       "B"),

            "KpPimInvarMass":               TriggerVariable("KpPimInvarMass",  "f", self.calculateKpPimInvarMass),
            "KmPipInvarMass":               TriggerVariable("KmPipInvarMass",  "f", self.calculateKmPipInvarMass),
            "KpKmInvarMass":                TriggerVariable("KpKmInvarMass",   "f", self.calculateKpKmInvarMass),
            "PipPimInvarMass":              TriggerVariable("PipPimInvarMass", "f", self.calculatePipPimInvarMass),
            "KpPimAngle":                   TriggerVariable("KpPimAngle",      "f", self.calculateKpPimAngle),
            "KmPipAngle":                   TriggerVariable("KmPipAngle",      "f", self.calculateKmPipAngle),
            "KpKmAngle":                    TriggerVariable("KpKmAngle",       "f", self.calculateKpKmAngle),
            "PipPimAngle":                  TriggerVariable("PipPimAngle",     "f", self.calculatePipPimAngle),
        }

        self.InputContainers = [
            FinalContainer(input_path, 'read', self.Variables),
        ]

        self.CutDispatcher = CutDispatcher({}, n_entries_full = self.getEntries())
        
        histograms_available = {
            'h_KpPimInvarMass_KmPipInvarMass': {
                'type': EHists.TH2F,
                'args': {
                    'title': "K^{-}#pi^{+} invariant mass vs K^{+}#pi^{-} invariant mass",
                    'x-axis-title': "M_{K^{+}#pi^{-},invar}, MeV/c^{2}",
                    'x-axis-nbins': 200,
                    'x-axis-range': (0., 2000.),
                    'x-variable': self.Variables["KpPimInvarMass"],
                    'y-axis-title': "M_{K^{-}#pi^{+},invar}, MeV/c^{2}",
                    'y-axis-nbins': 200,
                    'y-axis-range': (0., 2000.),
                    'y-variable': self.Variables["KmPipInvarMass"],
                },
            },
            'h_KpKmInvarMass': {
                'type': EHists.TH1F,
                'args': {
                    'title': "K^{+}K^{-} invariant mass",
                    'x-axis-title': "M_{K^{+}K^{-},invar}, MeV/c^{2}",
                    'x-axis-range': (0., 2000.),
                    'x-axis-nbins': 200,
                    'x-variable': self.Variables["KpKmInvarMass"],
                },
            },
            'h_PipPimInvarMass': {
                'type': EHists.TH1F,
                'args': {
                    'title': "#pi^{+}#pi^{-} invariant mass",
                    'x-axis-title': "M_{#pi^{+}#pi^{-},invar}, MeV/c^{2}",
                    'x-axis-range': (0., 2000.),
                    'x-axis-nbins': 200,
                    'x-variable': self.Variables["PipPimInvarMass"],
                },
            },
            'h_PipPimInvarMass_KpKmInvarMass': {
                'type': EHists.TH2F,
                'args': {
                    'title': "K^{+}K^{-} invariant mass vs #pi^{+}#pi^{-} invariant mass",
                    'x-axis-title': "M_{#pi^{+}#pi^{-},invar}, MeV/c^{2}",
                    'x-axis-nbins': 200,
                    'x-axis-range': (0., 2000.),
                    'x-variable': self.Variables["PipPimInvarMass"],
                    'y-axis-title': "M_{K^{+}K^{-},invar}, MeV/c^{2}",
                    'y-axis-nbins': 200,
                    'y-axis-range': (0., 2000.),
                    'y-variable': self.Variables["KpKmInvarMass"],
                },
            },
            
            'h_KpPimAngle_KmPipAngle': {
                'type': EHists.TH2F,
                'args': {
                    'title': "K^{-}#pi^{+} angle vs K^{+}#pi^{-} angle",
                    'x-axis-title': "#alpha_{K^{+}#pi^{-}}, rad",
                    'x-axis-nbins': 300,
                    'x-axis-range': (0., pi),
                    'x-variable': self.Variables["KpPimAngle"],
                    'y-axis-title': "#alpha_{K^{-}#pi^{+}}, rad",
                    'y-axis-nbins': 300,
                    'y-axis-range': (0., pi),
                    'y-variable': self.Variables["KmPipAngle"],
                },
            },
            'h_KpKmAngle': {
                'type': EHists.TH1F,
                'args': {
                    'title': "K^{+}K^{-} angle",
                    'x-axis-title': "#alpha_{K^{+}K^{-}}, rad",
                    'x-axis-range': (0., pi),
                    'x-axis-nbins': 300,
                    'x-variable': self.Variables["KpKmAngle"],
                },
            },
            'h_PipPimAngle': {
                'type': EHists.TH1F,
                'args': {
                    'title': "#pi^{+}#pi^{-} angle",
                    'x-axis-title': "#alpha_{#pi^{+}#pi^{-}}, rad",
                    'x-axis-range': (0., pi),
                    'x-axis-nbins': 300,
                    'x-variable': self.Variables["PipPimAngle"],
                },
            },
        }
        self.HistogramDispatcher = HistogramDispatcher(histograms_available)


    def calculateTwoParticlesInvarMass(self, particle_name_1, particle_name_2):
        p4       = TLorentzVector(0.0, 0.0, 0.0, 0.0)
        p4_auxil = TLorentzVector(0.0, 0.0, 0.0, 0.0)

        index = self.Variables[f"{particle_name_1}TrackIndex"].Content
        if particle_name_1 in ("Kp", "Km"):
            mass = m_K
        elif particle_name_1 in ("Pip", "Pim"):
            mass = m_pi
        p4_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
        p4_auxil.SetRho(self.Variables["tptot"].Content[index])
        p4_auxil.SetTheta(self.Variables["tth"].Content[index])
        p4_auxil.SetPhi(self.Variables["tphi"].Content[index])
        p4_auxil.SetE(sqrt(self.Variables["tptot"].Content[index] ** 2 + mass ** 2))
        p4 += p4_auxil

        index = self.Variables[f"{particle_name_2}TrackIndex"].Content
        if particle_name_2 in ("Kp", "Km"):
            mass = m_K
        elif particle_name_2 in ("Pip", "Pim"):
            mass = m_pi
        p4_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
        p4_auxil.SetRho(self.Variables["tptot"].Content[index])
        p4_auxil.SetTheta(self.Variables["tth"].Content[index])
        p4_auxil.SetPhi(self.Variables["tphi"].Content[index])
        p4_auxil.SetE(sqrt(self.Variables["tptot"].Content[index] ** 2 + mass ** 2))
        p4 += p4_auxil

        self.Variables[f"{particle_name_1}{particle_name_2}InvarMass"].Content = p4.M()

        
    calculateKpPimInvarMass  = partialmethod(calculateTwoParticlesInvarMass, "Kp", "Pim")
    calculateKmPipInvarMass  = partialmethod(calculateTwoParticlesInvarMass, "Km", "Pip")
    calculateKpKmInvarMass   = partialmethod(calculateTwoParticlesInvarMass, "Kp", "Km")
    calculatePipPimInvarMass = partialmethod(calculateTwoParticlesInvarMass, "Pip", "Pim")

    def calculateTwoParticlesAngle(self, particle_name_1, particle_name_2):
        p3_1 = TVector3(0.0, 0.0, 0.0)
        p3_2 = TVector3(0.0, 0.0, 0.0)

        index = self.Variables[f"{particle_name_1}TrackIndex"].Content
        p3_1.SetXYZ(1.0, 0.0, 0.0)
        p3_1.SetMag(self.Variables["tptot"].Content[index])
        p3_1.SetTheta(self.Variables["tth"].Content[index])
        p3_1.SetPhi(self.Variables["tphi"].Content[index])

        index = self.Variables[f"{particle_name_2}TrackIndex"].Content
        p3_2.SetXYZ(1.0, 0.0, 0.0)
        p3_2.SetMag(self.Variables["tptot"].Content[index])
        p3_2.SetTheta(self.Variables["tth"].Content[index])
        p3_2.SetPhi(self.Variables["tphi"].Content[index])

        self.Variables[f"{particle_name_1}{particle_name_2}Angle"].Content = acos(p3_1 * p3_2/(p3_1.Mag() * p3_2.Mag()))


    calculateKpPimAngle  = partialmethod(calculateTwoParticlesAngle, "Kp", "Pim")
    calculateKmPipAngle  = partialmethod(calculateTwoParticlesAngle, "Km", "Pip")
    calculateKpKmAngle   = partialmethod(calculateTwoParticlesAngle, "Kp", "Km")
    calculatePipPimAngle = partialmethod(calculateTwoParticlesAngle, "Pip", "Pim")
