from pprint import pprint
from typing import Tuple

from ROOT import TFile, TTree

from .Variable import Variable

def prod(iterable, *, start = 1):
    if len(iterable) == 0: return start

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

    def __init__(self, path: str, mode: str, branches_variables):
        if mode in ("read", "new", "recreate", "update"):
            self.Mode = mode
        else:
            raise ValueError(f"Wrong mode: '{mode}'")
        self.CurrentEntry = -1

        self.ContainerFile = TFile.Open(path, self.Mode)
        if self.Mode in ("read", "update"):
            self.Tree = self.ContainerFile.Get("tr_ph")
        elif self.Mode in ("new", "recreate"):
            self.Tree = TTree("tr_ph", "tr_ph")

        self.Variables = branches_variables
        for branch_name, variable in self.Variables.items():
            self.addBranch(branch_name, variable)


    def getEntries(self) -> int:
        return self.Tree.GetEntries()

    def getEntry(self, entry: int):
        self.Tree.GetEntry(entry)
        self.CurrentEntry = entry

    def getBranchSignature(self, branch_name):
        variable = self.Variables[branch_name]
        if variable.Sizes == (1,):
            return f"{branch_name}/{variable.getTypecode()}"
        named_sizes = tuple(map(
            lambda size: size if not isinstance(size, Variable) else size.Name,
            variable.Sizes
        ))
        return f"{branch_name}{''.join([f'[{size}]' for size in named_sizes])}/{variable.getTypecode()}"

    def fillEntry(self):
        self.Tree.Fill()
        self.CurrentEntry += 1

    def showEntry(self):
        print(f"Current entry: {self.CurrentEntry}")
        for branch_name, variable in self.Variables.items():
            print(f"{branch_name}: ",
                  f"{variable.Content}")

    """
    Adds a new variable to contain tr_ph data and a new branch when the tree is newly created
    Doesn't throw an exception when the tree is read-only
    """
    def addBranch(self, branch_name: str, variable: Variable):
        if self.Mode == "read" and not self.Tree.FindBranch(branch_name):
            raise ValueError(f"Branch not found: {branch_name}")

        if (self.Mode == "read" or (self.Mode == "update" and self.Tree.FindBranch(branch_name))):
            self.Tree.SetBranchAddress(branch_name, variable.getArray())
        elif ((self.Mode in ("new", "recreate")) or (self.Mode == "update" and not self.Tree.FindBranch(branch_name))):
            self.Tree.Branch(branch_name, variable.getArray(), self.getBranchSignature(branch_name))

    def dumpToFile(self):
        if self.Mode == "read":
            raise AttributeError("Cannot write into read-only file: 'self.ContainerFile.GetName()'")
        else:
            self.ContainerFile.cd()
            self.Tree.Write()
            self.ContainerFile.Save()
            return self.Tree

    def close(self):
        if self.Mode in ("new", "recreate"): del self.Tree
        self.ContainerFile.Close()
