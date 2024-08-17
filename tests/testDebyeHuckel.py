import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from Agnes import *


class DHProb:
  def __init__(self, m=1, n=1, beta=1):
    self.m = m
    self.n = n
    self.beta = beta

  def rhsFunc(self, x, y, u=None):
    pi = np.pi
    m = self.m
    n = self.n
    return np.cos(m*pi*x)*np.cos(n*pi*y)

  def uExact(self, x, y):
    pi = np.pi
    m = self.m
    n = self.n
    beta = self.beta
    A0 = 1.0/((m**2 + n**2)*pi**2 + beta**2)

    return A0 * self.rhsFunc(x, y)

def test_DebyeHuckel():
  '''Solve the Debye-Huckel equation with homogeneous Neumann BCs'''

  n = 128
  beta = 1.0
  prob = DHProb(beta=beta)

  print('Building mesh...')
  mesh = meshRectangle(nx=n, ny=n, ax=0, bx=1, ay=0, by=1)
  ds = DiscreteSpace(mesh, 1)

  twoForms = (
    LaplacianTwoForm(),
    MassTwoForm(coeff=beta*beta)
  )

  quad = GaussRule(2)
  
  
  
  oneForms = (
    VarCoeffOneForm(quad, prob.rhsFunc), 
    )

  print('Building matrix and vector...')
  assembler = Assembler(ds, twoForms, oneForms)

  (A, b) = assembler.assemble()

  print('Solving equation...')
  solnVec = spla.spsolve(A, b)

  u = DiscreteFunction(ds, 'u')
  u.setVector(solnVec)

  print('Writing solution...')

  writer = VTKWriter('DH-{}.vtu'.format(n))
  writer.addMesh(mesh)
  writer.addField('u', solnVec, 0)
  writer.write()





if __name__=='__main__':

  test_DebyeHuckel()



