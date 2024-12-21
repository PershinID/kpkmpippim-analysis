from math import pi, sqrt

from ROOT import TLorentzVector

from Variable import Variable

from PhysicalConstants import m_pi, m_K

from CMD3ContainerV9 import CMD3ContainerV9
from PreliminaryContainer import PreliminaryContainer

from Analysis import Analysis, EHists


class PreliminaryAnalysis(Analysis):
    def __init__(self, path: str, input_path: str, output_path: str):
        Analysis.__init__(self, path)
        
        size_variables = { # needed to define variables in self.Variables
            "nt_total":         Variable("nt_total",            "as", max_value = 10),
            "nv_total":         Variable("nv_total",            "as", max_value = 10),
            "nt":               Variable("nt",                  "as", max_value = 10),
            "ntlxe_total":      Variable("ntlxe_total",         "as", max_value = 10),
            "ntlxe":            Variable("ntlxe",               "as", max_value = 10),
        }
        self.Variables = {
            ##CMD3 branches
            "emeas":            Variable("emeas",               "f"),
            "demeas":           Variable("demeas",              "f"),
            "runnum":           Variable("runnum",              "i"),
            "evnum":            Variable("evnum",               "i"),
            "ecaltot":          Variable("ecaltot",             "f"),
            "ecalneu":          Variable("ecalneu",             "f"),
            "psumch":           Variable("psumch",              "f"),
            "psumnu":           Variable("psumnu",              "f"),
            "tnhit":            Variable("tnhit",               "i", sizes = (size_variables["nt"],)),
            "tlength":          Variable("tlength",             "f", sizes = (size_variables["nt"],)),
            "tphi":             Variable("tphi",                "f", sizes = (size_variables["nt"],)),
            "tth":              Variable("tth",                 "f", sizes = (size_variables["nt"],)),
            "tptot":            Variable("tptot",               "f", sizes = (size_variables["nt"],)),
            "tphiv":            Variable("tphiv",               "f", sizes = (size_variables["nt"],)),
            "tthv":             Variable("tthv",                "f", sizes = (size_variables["nt"],)),
            "tptotv":           Variable("tptotv",              "f", sizes = (size_variables["nt"],)),
            "trho":             Variable("trho",                "f", sizes = (size_variables["nt"],)),
            "tdedx":            Variable("tdedx",               "f", sizes = (size_variables["nt"],)),
            "tz":               Variable("tz",                  "f", sizes = (size_variables["nt"],)),
            "tchi2r":           Variable("tchi2r",              "f", sizes = (size_variables["nt"],)),
            "tchi2z":           Variable("tchi2z",              "f", sizes = (size_variables["nt"],)),
            "tchi2ndf":         Variable("tchi2ndf",            "f", sizes = (size_variables["nt"],)),
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
            "tindlxe":          Variable("tindlxe",             "f", sizes = (size_variables["nt"],)),
            "txyzatcl":         Variable("txyzatcl",            "f", sizes = (size_variables["nt"], 3,)),
            "txyzatlxe":        Variable("txyzatlxe",           "f", sizes = (size_variables["nt"], 3,)),
            "tenconv":          Variable("tenconv",             "f", sizes = (size_variables["nt"],)),
            "ntlxelayers":      Variable("ntlxelayers",         "f", sizes = (size_variables["nt"],)),
            "tlxenhit":         Variable("tlxenhit",            "f", sizes = (size_variables["nt"],)),
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

            ##User branches
            "DeltaE":           Variable("DeltaE",              "f"),
            "TotalP":           Variable("TotalP",              "f"),
        }
        self.Variables.update(size_variables) # enclose dict in itself

        self.InputContainer = CMD3ContainerV9(input_path, self.Variables)
        self.initializeEntriesChecklist()
        self.OutputContainer = PreliminaryContainer(output_path, "recreate", self.Variables)

        self.Cuts = {
            'nt': lambda: self.Variables["nt"].getContent() != 4,
            'tcharge': lambda: sum(self.Variables["tcharge"].getContent()) != 0,
            'tnhit': lambda: len(list(filter(
                lambda x: x <= 9,
                self.Variables["tnhit"].getContent()
            ))) != 0,
            'tptot': lambda: len(list(filter(
                lambda x: x < 50,
                self.Variables["tptot"].getContent()
            ))) != 0,
            'tth': lambda: len(list(filter(
                lambda x: x < 0.9 or x > pi - 0.9,
                self.Variables["tth"].getContent()
            ))) != 0,
            'trho': lambda: len(list(filter(
                lambda x: abs(x) > 0.4,
                self.Variables["trho"].getContent()
            ))) != 0,
            'tz': lambda: len(list(filter(
                lambda x: abs(x) > 10.0,
                self.Variables["tz"].getContent()
            ))) != 0,
            'DeltaE': lambda: (
                self.Variables["DeltaE"].getContent() < -700.0 or
                self.Variables["DeltaE"].getContent() > -300.0
            ),
            'TotalP': lambda: (
                self.Variables["TotalP"].getContent() < 100.0
            ),
        }

        self.Histograms = {
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
                'type': EHists.TH1I,
                'args': {
                    'title': "Track momentum",
                    'x-axis-title': "P_{track}, MeV/c",
                    'x-axis-nbins': 2000,
                    'x-axis-range': (0.0, 2000.0),
                    'x-variable': self.Variables["tptot"],
                },
            },
            'h_tth': {
                'type': EHists.TH1I,
                'args': {
                    'title': "Track #theta coordinate",
                    'x-axis-title': "#theta_{track}, rad",
                    'x-axis-nbins': 600,
                    'x-axis-range': (0.0, pi),
                    'x-variable': self.Variables["tth"],
                },
            },
            'h_tphi': {
                'type': EHists.TH1I,
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
                    'x-axis-title': "#frac{dE}{dx}_{track}, a.u.",
                    'x-axis-nbins': 500,
                    'x-axis-range': (0.0, 50.0),
                    'x-variable': self.Variables["tdedx"],
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
                }
            },
        }

    def getEntry(self, n_entry):
        self.InputContainer.getEntry(n_entry)

        ##User variables calculation
        ##DeltaE and TotalP calculation
        de_content = self.Variables["DeltaE"].Content
        tp_content = self.Variables["TotalP"].Content
        
        tp_lorentz = TLorentzVector(0.0, 0.0, 0.0, 0.0)
        tp_auxil = TLorentzVector(0.0, 0.0, 0.0, 0.0)
        
        for (p, th, phi) in zip(self.Variables["tptot"].getContent(),
                                self.Variables["tth"].getContent(),
                                self.Variables["tphi"].getContent()):
            tp_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
            tp_auxil.SetRho(p)
            tp_auxil.SetTheta(th)
            tp_auxil.SetPhi(phi)
            tp_auxil.SetE(sqrt(p * p + m_pi * m_pi))
            tp_lorentz += tp_auxil

        de_content[0] = tp_lorentz.E() - 2.0 * self.Variables["emeas"].getContent()
        tp_content[0] = tp_lorentz.Rho()


if __name__ == '__main__':
    test_analysis = PreliminaryAnalysis('test_hists.root', '/spoolA/idpershin/datasets/dataset_e950.root',
                                                           'test_data.root')
    test_analysis.getEntry(5)
    print(test_analysis.Variables["DeltaE"].getContent(), ' ', test_analysis.Variables["TotalP"].getContent())
