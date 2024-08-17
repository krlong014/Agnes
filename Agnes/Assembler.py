import numpy as np
import scipy.sparse as sp
from collections.abc import Iterable
from .Triangle import Triangle
from .LoadableMesh import LoadableMesh
from .DiscreteSpace import DiscreteSpace
from .OneForm import OneForm
from .TwoForm import TwoForm

class Assembler:
  '''
  Assembler drives the matrix/vector assembly process
  '''

  def __init__(self, ds : DiscreteSpace,
               twoForms : Iterable, 
               oneForms : Iterable):
    
    self._ds = ds
    self._oneForms = oneForms
    self._twoForms = twoForms

    
  def assemble(self, buildMat=True, buildVec=True):

    N = self._ds.numDofs()
    print('N=', N)
    A = None
    b = None
    if buildMat:
      A = sp.dok_matrix((N,N))
    if buildVec:
      b = np.zeros(N)

    # Indicate whether debugging output is desired
    debug = False

    # Loop over elements
    mesh = self._ds.mesh()
    for ie, eVerts in enumerate(mesh.elems):

      if debug: print('processing element #%d' % ie, " vertices=", eVerts)
      # Get coordinates of the three elements
      va = mesh.verts[eVerts[0]]
      vb = mesh.verts[eVerts[1]]
      vc = mesh.verts[eVerts[2]]

      # Create a triangle object with the three vertices
      T = Triangle(va, vb, vc)

      if buildMat:
        for aForm in self._twoForms:

          A_loc = aForm.localMat(T, eVerts)

          testID = aForm.testID()
          unkID = aForm.unkID()

          testDofs, unkDofs = self._ds.getDofs(eVerts, (testID, unkID))
          
          for i in range(0,3):
            r = testDofs[i]
            for j in range(0,3):
              c = unkDofs[j]
              A[r,c] += A_loc[i,j]

      if buildVec:
        for bForm in self._oneForms:
          b_loc = bForm.localVec(T, eVerts)
          testID = bForm.testID()

          testDofs, = self._ds.getDofs(eVerts, testID)

          for i in range(0,3):
            r = testDofs[i]
            b[r] += b_loc[i]

    A = A.tocsr()

    return (A,b)
            


