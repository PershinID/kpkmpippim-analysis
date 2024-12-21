from array import array
from copy import copy

def prod(iterable, *, start = 1):
    if len(iterable) is 0: return start

    _prod = start
    for number in iterable: _prod *= int(number)
    return _prod


class Variable:
    """
    Class for a humane access to data in TTree branch
    Supports:
    - plain data: int, float (typecode 'i', 'f' respectively);
    - sizes for arrays in branches (typecode 'as') (int type used);
    - arrays of variable sizes (variable sizes first, then fixed sizes)
    """
    def __init__(self, name, typecode, *, sizes = (1,), max_value = None):
        if typecode == 'as':
            if not max_value: raise ValueError(f"'max_value' must be entered: {name}")
            if sizes[0] != 1: raise ValueError(f"'sizes' must be equal to 1: {name}")

            self.Content = array('i', [0])
            self.MaxValue = max_value
        elif typecode in ('i', 'f'):
            self.MaxSizes = tuple(map(
                lambda size: size if not isinstance(size, Variable) else size.MaxValue,
                sizes
            ))
            self.Content = array(typecode, [0] * prod(self.MaxSizes))

        self.Typecode = typecode
        self.Sizes = sizes
        self.Name = name
                    

    def __index__(self):
        if self.Typecode == 'as':
            return self.Content[0]
        else:
            raise ValueError(f"Only array-size variables support arithmetic operations: {self.Name}")

    def __int__(self):
        if self.Typecode == 'as':
            return self.Content[0]
        else:
            raise ValueError(f"Only array-size variables support arithmetic operations: {self.Name}")    

    def getContent(self):
        if self.Sizes == (1,):
            return self.Content[0]
        else:
            return self.getContentList()

    def getContentList(self):
        auxil_array = []
        max_size = prod(self.Sizes)
        main_array = self.Content.tolist()[:max_size]
        for size in reversed(self.Sizes[1:]):
            max_size //= size
            for i in range(max_size):
                auxil_array.append(main_array[i * size : (i + 1) * size])
            main_array = auxil_array
            auxil_array = []
        return main_array


    def getBranchSignature(self, branch_name):
        if self.Typecode == 'as':
            return f"{branch_name}/I"
        elif self.Sizes[0] == 1:
            return f"{branch_name}/{self.Typecode.capitalize()}"
        else:
            named_sizes = list(map(
                lambda size: size if not isinstance(size, Variable) else size.Name,
                self.Sizes
            ))
            return f"{branch_name}{''.join([f'[{size}]' for size in named_sizes])}/{self.Typecode.capitalize()}"


if __name__ == '__main__':
    nt_exam = Variable('nt_exam', 'as', max_value = 10)
    txyzcorr_exam = Variable('txyzcorr_exam', 'f', sizes = (nt_exam, 3, 3,))
    nt_exam._content[0] = 10
    print(txyzcorr_exam.getContent())
    print(txyzcorr_exam.getBranchSignature())
