import numpy as np

class P1Basis:

  @staticmethod
  def phi(x, y):
    return np.array([1.0-x-y, x, y])
  
  @staticmethod
  def grad_phi(x, y):
    return np.array([[-1.0, 1.0, 0.0],
                     [-1.0, 0.0, 1.0]])
  

  