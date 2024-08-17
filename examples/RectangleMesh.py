import numpy as np
from Agnes import *

if __name__ == '__main__':

  n = 4
  mesh = meshRectangle(nx=n, ny=n)

  w = VTKWriter('testRect.vtu')
  w.addMesh(mesh)
  w.write()
