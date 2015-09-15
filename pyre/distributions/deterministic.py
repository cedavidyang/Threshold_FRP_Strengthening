import numpy as np

class DummyObject(object):
    def __init__(self, value):
        self.value = value
    def rvs(self, size = None):
        return self.value * np.ones(size)

class Deterministic(object):
    def __init__(self, name, mu, sigma=None):
        self.name = name
        self.value = mu
        self.rv = DummyObject(mu)