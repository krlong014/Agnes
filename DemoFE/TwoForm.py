from abc import ABC, abstractmethod


class TwoForm(ABC):
  '''
  Abstract interface for two forms
  '''

  @abstractmethod
  def localMat(self, tri):
    '''
    Interface for computing a local matrix on a triangle. 
    '''
    pass


class QuadratureOneForm(OneForm):
  '''
  Abstract interface for one form to be computed by quadrature 
  '''

  @abstractmethod
  def localVec(self, tri):
    '''
    Interface for computing a local vector on a triangle. 
    '''
    pass
