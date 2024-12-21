from pprint import pprint
from typing import Tuple

from ROOT import TFile, TTree

from Variable import Variable

def prod(iterable, *, start = 1):
    if len(iterable) is 0: return start

    _prod = start
    for number in iterable: _prod *= number
    return _prod


class Container:
    """
    Base class for tr_ph containment.
    Derived classes deliver:
    - branches set (self.dictBranches)
    - methods FillEntry() (because it has to provide some calculations)
    """

    def __init__(self, path: str, mode: str, branches: dict):
        if mode in ("read", "new", "recreate", "update"):
            self.Mode = mode
        else:
            raise ValueError(f"Wrong mode: '{mode}'")
        self.CurrEntry = -1

        self.ContainerFile = TFile.Open(path, self.Mode)
        if self.Mode in ("read", "update"):
            self.Tree = self.ContainerFile.Get("tr_ph")
        elif self.Mode in ("new", "recreate"):
            self.Tree = TTree("tr_ph", "")

        self.Branches = branches
        for branch_name, variable in self.Branches.items():
            self.addBranch(branch_name, variable["variable"])


    def getEntries(self) -> int:
        return self.Tree.GetEntries()

    def getEntry(self, entry: int):
        self.Tree.GetEntry(entry)
        self.CurrEntry = entry

    def fillEntry(self):
        self.Tree.Fill()
        self.CurrEntry += 1

    def showEntry(self):
        print(f"Current entry: {self.currEntry}")
        for branch_name, variable in self.Branches.items():
            print(f"{branch_name}: ",
                  f"{variable.getContent()}")

    """
    Adds a new variable to contain tr_ph data and a new branch when the tree is newly created
    Doesn't throw an exception when the tree is read-only
    """
    def addBranch(self, branch_name: str, variable: Variable):
        if self.Mode == "read" and not self.Tree.FindBranch(branch_name):
            raise ValueError(f"Branch not found: {branch_name}")

        if (self.Mode == "read" or (self.Mode == "update" and self.Tree.FindBranch(branch_name))):
            self.Tree.SetBranchAddress(branch_name, variable.Content)
        elif ((self.Mode in ("new", "recreate")) or (self.Mode == "update" and not self.Tree.FindBranch(branch_name))):
            self.Tree.Branch(branch_name, variable.Content, variable.getBranchSignature(branch_name))

    def dumpToFile(self):
        if self.Mode == "read":
            raise AttributeError("Cannot write into read-only file: 'self.ContainerFile.GetName()'")
        else:
            self.ContainerFile.cd()
            self.Tree.Write()
            self.ContainerFile.Save()
