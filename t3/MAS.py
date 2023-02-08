class SAM(): # successive aproximations method
    def __init__(self, phi, x, iter):
        self.phi = phi
        self.x = x
        self.iter = iter

    def find_root(self):
        x = self.x
        for i in range (0, self.iter):
            x = eval(self.phi)
        return x
