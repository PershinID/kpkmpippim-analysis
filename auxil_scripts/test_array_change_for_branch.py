from array import array

from ROOT import TFile, TTree

if __name__ == "__main__":
    t_example = TTree("mytree", "mytree")
    n_example = 10
    a1_example = array('i', [0] * n_example)
    a2_example = array('i', [0] * n_example)
    t_example.Branch("myvar", a1_example, "myvar[10]/I")
    t_example.SetBranchAddress("myvar", a2_example)
    del a1_example
