from abc import abstractmethod

class Model:

    # returns a list of strings
    @abstractmethod
    def flag(self):
        pass

    # returns a list of reactions
    @abstractmethod
    def build(self):
        pass

class Mortal(Model):

    def __init__(self):
        self.speciesnames = ['X']

    def flag(self):
        return ['']

    def build(self):
        return [f'{species} -> ' for species in self.speciesnames]

class Eater(Model):

    def __init__(self):
        self.speciesnames = ['X']

    def flag(self):
        return ['']

    def build(self):
        return [f'{species} + R -> {species} + {species}' for species in self.speciesnames]

class Ager(Model):

    def __init__(self):
        self.speciesnames = ['young', 'old']

    def flag(self):
        return ['young', 'old']

    def build(self):
        return ['young -> old']



class Bacteria(Mortal, Eater, Ager):

    def __init__(self):
        Mortal.__init__(self)
        Eater.__init__(self)
        Ager.__init__(self)

        self.speciesnames = [f'Bacteria{f}' for f in self.flag()]

    def flag(self):
        return [a+b+c for a in Mortal.flag(self) for b in Eater.flag(self) for c in Ager.flag(self)]

    def build(self):
        # TODO: this is plus
        ret = []
        ret += Mortal.build(self)
        ret += Eater.build(self)
        ret += Ager.build(self)
        return ret

def simulate(model):
    pass

def main():
    b = Bacteria()
    print(b.flag())
    print(b.speciesnames)
    print(b.build())
    simulate(Bacteria)
#    Bacteria = Mortal * Eater * Ager + Phage
#    NewClass = Ager + Mortal

if __name__ == '__main__':
    main()
