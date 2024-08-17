from .LoadableMesh import LoadableMesh
import numpy as np
import itertools as it

def meshRectangle(nx:int=4, ny:int=4, 
                  ax:float=0.0, bx:float=1.0, 
                  ay:float=0.0, by:float=1.0):
  '''
  Mesh the rectangle $[a_x,a_y] \\otimes [b_x,b_y]$, with nx, ny intervals
  in each direction. There are therefore 2 nx ny triangular elements. 

  '''

  X = np.linspace(ax, bx, nx+1)
  Y = np.linspace(ay, by, ny+1)  

  mesh = LoadableMesh()
  indexPairToVertMap = {}

  # Put the vertices into the mesh
  v = 0
  for ix,x in enumerate(X):
    for iy,y in enumerate(Y):
      mesh.addVertex((x,y))
      indexPairToVertMap[(ix,iy)] = v
      v += 1

  
  for ix in range(nx):
    for iy in range(ny):
      isEven = (((ix+iy) % 2)==0)
      corners = (indexPairToVertMap[(ix,iy)],
                 indexPairToVertMap[(ix+1,iy)],
                 indexPairToVertMap[(ix+1,iy+1)],
                 indexPairToVertMap[(ix,iy+1)])

      if isEven:
        T1 = (0,1,2)
        T2 = (0,2,3)
      else:
        T1 = (0,1,3)
        T2 = (1,2,3)


      edges = list(it.combinations(T1, 2)) + list(it.combinations(T2,2))
      #print('({},{}), even={}, edges are {}'.format(ix, iy, isEven,edges))
      
      for e in edges:
        mesh.addSide(corners[e[0]], corners[e[1]], 0)
      

      mesh.addElem(corners[T1[0]], corners[T1[1]], corners[T1[2]])
      mesh.addElem(corners[T2[0]], corners[T2[1]], corners[T2[2]])
      
      
  return mesh


  



