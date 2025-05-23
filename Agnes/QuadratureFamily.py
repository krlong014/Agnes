from abc import ABC, abstractmethod
import numpy as np
import pprint

class QuadratureFamily(ABC):
  def __init__(self, name='Quad', order=2):
    self._order = order
    self._name = name

  @abstractmethod
  def getRule(self, dim):
    pass

  def name(self):
    return self._name

  def order(self):
    return self._order

  def test(self, dim, tol=1.0e-14):
    return self.getRule(dim).test(tol)
  
  

