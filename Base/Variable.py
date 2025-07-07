from array import array
from copy import copy

from operator import mul
from functools import reduce
from math import pi, sqrt
import random

from ROOT import TLorentzVector

from .PhysicalConstants import m_pi

def prod(iterable, *, start = 1):
    return reduce(mul, iterable, start)


class Variable:
    """
    Class for a humane access to data in TTree branch
    Supports:
    - plain data: int, float (typecode 'i', 'f', 'b' respectively);
    - sizes for arrays in branches (typecode 'as') (int type used);
    - arrays of variable sizes (variable sizes first, then fixed sizes)
    """
    ## Correspondence between Variable typecodes (Python array typecodes except array size) and TBranch typecodes 
    Typecodes = {
        'as': 'I', ## array size -> signed int (Int_t)
        'i': 'I',  ## signed int -> signed int (Int_t)
        'f': 'F',  ## float -> float (Float_t)
        'b': 'B',  ## signed char -> signed char (Char_t)
        'B': 'b',  ## unsigned char -> unsigned char (UChar_t)
        'h': 'S',  ## signed short -> signed short (Short_t)
        'H': 's',  ## unsigned short -> unsigned short (UShort_t)
    }
    
    def __init__(self, name, typecode, *, sizes = (1,), max_value = None):
        if typecode == 'as':
            if not max_value: raise ValueError(f"'max_value' must be entered: {name}")
            if sizes[0] != 1: raise ValueError(f"'sizes' must be equal to 1: {name}")

            self.MaxValue = max_value
            self._Content = array('i', [0])
        elif typecode in Variable.Typecodes:
            self.MaxSizes = tuple(map(
                lambda size: size if not isinstance(size, Variable) else size.MaxValue,
                sizes
            ))
            self._Content = array(typecode, [0] * prod(self.MaxSizes))

        self.Typecode = typecode
        self.Sizes = sizes
        self.Name = name
                    

    def __index__(self):
        if self.Typecode == 'as':
            return self._Content[0]
        else:
            raise ValueError(f"Only array-size variables support arithmetic operations: {self.Name}")

    def __int__(self):
        if self.Typecode == 'as':
            return self._Content[0]
        else:
            raise ValueError(f"Only array-size variables support arithmetic operations: {self.Name}")

    def getTypecode(self):
        return Variable.Typecodes[self.Typecode]
        
    def getContent(self):
        if self.Sizes == (1,):
            return self._Content[0]
        else:
            auxil_array = []
            sizes = list(map(int, self.Sizes))
            full_size = prod(sizes)
            main_array = self._Content.tolist()[:full_size]
            for size in reversed(sizes[1:]):
                full_size //= size
                for i in range(full_size):
                    auxil_array.append(main_array[i * size : (i + 1) * size])
                main_array = auxil_array
                auxil_array = []
            return main_array

    def setContent(self, value):
        if self.Sizes == (1,):
            self._Content[0] = value
        else:
            ## check if value has sizes smaller than or equal to variable's ones
            value_tmp = value
            sizes = list(map(int, self.Sizes))
            for size in sizes:
                if len(value_tmp) != size: raise ValueError(f"Value has sizes other than variable's ones: value size = {len(value_tmp)}, variable size = {size}")
                value_tmp = value_tmp[0]

            ## flatten array
            auxil_array = []
            full_size = 1
            main_array = value
            for size in sizes[:-1]:
                full_size *= size
                for i in range(full_size):
                    auxil_array += main_array[i]
                main_array = auxil_array
                auxil_array = []

            ## fulfill array with zeros
            ## [ i ] [ j ] [ k ] [ l ]
            ##   |     |     |  *  |
            ##   |     |     ---------> deep_size
            ##   |  *  |
            ##   ----------> shallow_size
            shallow_size = prod(sizes)
            deep_size = 1
            for size, max_size in reversed(list(zip(sizes, self.MaxSizes))):
                n_zeroes = deep_size * (max_size - size)
                deep_size *= size
                shallow_size //= size
                if size == max_size:
                    continue

                for i in range(shallow_size, 0, -1):
                    main_array[i * deep_size : i * deep_size] = [0] * n_zeroes
                deep_size //= size
                deep_size *= max_size
                
            self._Content = array(self.Typecode, main_array)
        
    Content = property(getContent, setContent)
        
    ## Methods for dealing with ROOT TTree branches
    def getArray(self):
        return self._Content
                

class TriggerVariable(Variable):
    def __init__(self, name, typecode, trig_func, *, sizes = (1,)):
        if typecode == 'as':
            raise ValueError('Lazy variable cannot be array-size (typecode "as"): {name}')
        Variable.__init__(self, name, typecode, sizes = sizes)
        self.TrigFunc = trig_func

    def getContent(self):
        self.TrigFunc()
        return Variable.getContent(self)

    Content = property(getContent, Variable.setContent)


if __name__ == '__main__':
    nt_exam = Variable('nt_exam', 'as', max_value = 10)
    txyzcorr_exam = Variable('txyzcorr_exam', 'f', sizes = (nt_exam, 3, 3,))
    nt_exam.Content = 6
    print(txyzcorr_exam.Content)
    new_txyzcorr = 6 * [[[1.,2.,3.],
                         [4.,5.,6.],
                         [7.,8.,9.]]]
    txyzcorr_exam.Content = new_txyzcorr
    print(txyzcorr_exam.Content)

    emeas_exam = Variable('emeas_exam', 'f')
    emeas_exam.Content = 1000.0
    tptot_exam = Variable('tptot_exam', 'f', sizes = (nt_exam,))
    tth_exam = Variable('tth_exam', 'f', sizes = (nt_exam,))
    tphi_exam = Variable('tphi_exam', 'f', sizes = (nt_exam,))
    tptot_exam.Content = [random.uniform(20., 30.) for i in range(nt_exam.Content)]
    tth_exam.Content = [random.uniform(0., pi) for i in range(nt_exam.Content)]
    tphi_exam.Content = [random.uniform(0., 2.0 * pi) for i in range(nt_exam.Content)]
    print(tphi_exam.Content)

    tnhit_exam = Variable('tnhit_exam', 'i', sizes = (4,))
    tnhit_exam.Content = [0,1,2,3]
    print(tnhit_exam.Content)
    
    def delta_E_func():
        tp_lorentz = TLorentzVector(0.0, 0.0, 0.0, 0.0)
        tp_auxil = TLorentzVector(0.0, 0.0, 0.0, 0.0)
        
        for (p, th, phi) in zip(tptot_exam.Content,
                                tth_exam.Content,
                                tphi_exam.Content):
            tp_auxil.SetPxPyPzE(1.0, 0.0, 0.0, 0.0)
            tp_auxil.SetRho(p)
            tp_auxil.SetTheta(th)
            tp_auxil.SetPhi(phi)
            tp_auxil.SetE(sqrt(p * p + m_pi * m_pi))
            tp_lorentz += tp_auxil

        de_exam.Content = tp_lorentz.E() - 2.0 * emeas_exam.Content

    de_exam = TriggerVariable('de_exam', 'f', delta_E_func)
    delta_E_func.isTriggered = False
    print(de_exam.Content)
    print(delta_E_func.isTriggered)
    print(de_exam.Content)
