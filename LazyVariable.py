from Variable import Variable

class LazyVariable(Variable):
    def __init__(self, name, typecode, def_func, *, sizes = (1,)):
        if typecode == 'as':
            raise ValueError('Lazy variable cannot be array-size (typecode "as"): {name}')
        Variable.__init__(self, name, typecode, sizes)
        self.DefFunc = def_func
        self.CalcFlag = False

    def getContent(self):
        self.DefFunc()
        if self.Sizes == (1,):
            return self.Content[0]
        else:
            return self.getContentList()
        self.CalcFlag = True

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

        
