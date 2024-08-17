from .LoadableMesh import LoadableMesh
from PyUtils import NamedObject 
from copy import deepcopy
import numpy as np




class DiscreteSpace:
  def __init__(self, mesh, numFuncs):
    self._mesh = mesh
    self._nf = numFuncs
    self._nDofs = numFuncs * len(mesh.verts)

  def numFuncs(self):
    return self._nf
  
  def numDofs(self):
    return self._nDofs

  def getDofs(self, elemID, funcIDs):

    verts = self._mesh.elems[elemID]

    if isinstance(funcIDs, int):
      return (self._dofsForFunc(verts, funcIDs),)      
    elif isinstance(funcIDs, (list, tuple)):
      rtn = []
      for f in funcIDs:
        rtn.append(self._dofsForFunc(verts,f))
      return rtn
  
  def _dofsForFunc(self, verts, fid):

    rtn = (
      self._nf * verts[0] + fid,
      self._nf * verts[1] + fid,
      self._nf * verts[2] + fid
    )

    return rtn
  
  def getDof(self, nodeID, funcID):
    return self._nf * nodeID + funcID


class _DiscFuncElem(NamedObject):

  def __init__(self, master, funcIndex : int, name : str):

    super().__init__(name=name)

    self._master = master
    self._funcIndex = funcIndex

  def setValue(self, node, val):
    dof = self._master._ds.getDof(node, self._funcIndex)
    self._master._vec[dof] = val

  def interpolate(self, nodes, phiAtQuadPts):

    dofs = self._master._ds._dofsForFunc(nodes, self._funcIndex)
    nodalVals = self._master.getValues(dofs)

    return phiAtQuadPts @ nodalVals
  
  def copyVecSlice(self):
    nDofs = self._master._ds.numDofs()
    nFuncs = self._master._ds.numFuncs()
    nNodes = nDofs // nFuncs

    v = self._master.getVector().reshape((nNodes, nFuncs))
    return deepcopy(v[:,self._funcIndex])
                                         

class _DiscFunc(NamedObject):

  def __init__(self, ds, name : str, vec = None):
    super().__init__(name=name)

    self._ds = ds

    if vec==None:
      self._vec = np.zeros(ds.numDofs())
    else:
      self._vec = vec

    self._elems = []

    for i in range(self._ds.numFuncs()):
      dfe = _DiscFuncElem(self, i, name + '[{}]'.format(i))
      print('dfe has type ', type(dfe))
      self._elems.append(dfe)

  def getVector(self):
    return self._vec
  
  def getValues(self, dofs):
    return np.array([self._vec.take(dofs),])
  
  def setVector(self, vec):
    if len(vec) != len(self._vec):
      raise ValueError('Bad vector size in _DiscFunc.setVector: '
                       'expected size {}, got {}'.format(len(self._vec),
                                                         len(vec)))

  def __getitem__(self, i : int):
    if i<0 or i>=self._ds.numFuncs():
      raise IndexError('_DiscFunc component index {} out of range '
                       '[{},{}).'.format(i, 0, self._ds.numFuncs()))
    return self._elems[i]
  
  def __deepcopy__(self):
    '''
    deepcopy() makes a deep copy of the function's vector, a shallow
    copy of everything else.
    '''
    return _DiscFunc(self._ds, self.name(), deepcopy(self._vec))
  

  



  
def DiscreteFunction(ds : DiscreteSpace, name : str):
  return _DiscFunc(ds, name)

def silly(f):
  assert(isinstance(f, _DiscFuncElem))

def test_DS():

  from LoadableMesh import TwoElemSquare
  from VTKWriter import VTKWriter
  from P1Basis import P1Basis
  from Quadrature import GaussRule
  from Triangle import Triangle
  mesh = TwoElemSquare()

  ds = DiscreteSpace(mesh, 2)

  U = DiscreteFunction(ds, 'u')
  u = U[0]
  v = U[1]

  assert(isinstance(u, _DiscFuncElem))
  assert(isinstance(v, _DiscFuncElem))

  print('v has type{}'.format(type(v)))

  for i,vert in enumerate(mesh.verts):
    u.setValue(i, i)
    v.setValue(i, np.sqrt(i))

  print('Writing...')
  wr = VTKWriter('DSTest.vtu')
  wr.addMesh(mesh)
  wr.addField('u', U.getVector(), 0)
  wr.addField('v', U.getVector(), 1)
  wr.write()

  quad = GaussRule(3)
  phiAtQP = quad.evalFunc(P1Basis.phi, 3)

  print('\nphi at ref quad pts')
  print(phiAtQP)

  print('u on triangle (0,1,2)')
  uVals = u.interpolate((0,1,2), phiAtQP)
  print(uVals)

  (A, B, C) = (mesh.verts[0], mesh.verts[1], mesh.verts[2])
  tri = Triangle(A, B, C)
  print('Triangle: ')
  for vtx in (A,B,C):
    print('\t{}'.format(vtx))

  XY = tri.ref_to_phys(quad.X())

  for xy, f in zip(XY, uVals):
    fInt = 0*(1-xy[0]) + 1*(xy[0]-xy[1]) + 2*xy[1]
    print('{:12.5g} {:12.5g}   {:12.5g} {:12.5}'.format(xy[0], xy[1],
                                                         f, fInt))


  

if __name__=='__main__':
  test_DS()