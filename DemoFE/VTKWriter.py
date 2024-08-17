import numpy as np
from .LoadableMesh import *
from PyUtils.XMLHeader import *
from .DiscreteSpace import _DiscFuncElem

class VTKWriter:

    def __init__(self, file : str):
        
           
        self.file = open(file, 'w')
        self.fields = {}
        self.mesh = None

    def addMesh(self, mesh):
        self.mesh = mesh

    def addField(self, name, f, funcIndex=0):

      print('type of f is ', type(f))
      print('type of self is ', type(self))
        
      if isinstance(f, np.ndarray):
        if self.mesh == None:
          raise RuntimeError('VTKWriter must add mesh before'
                             ' adding fields')
        nNodes = len(self.mesh.verts)
        if len(f) % nNodes != 0:
          raise ValueError('Vector size isn\'t an integer '
                           'multiple of # mesh nodes')
        fVec = f.reshape((nNodes, len(f)//nNodes))[:,funcIndex]
                    
        self.fields[name] = fVec

      elif isinstance(f, _DiscFuncElem):
        self.addField(name, f.copyVecSlice())

      else:
        raise TypeError('Argument f should be a numpy array '
                        'or a _DiscFuncElem; found type {}'\
                          .format(type(f)))
        
      

    def write(self):

        head = XMLHeader('VTKFile')
        head.addAttribute('type', 'UnstructuredGrid')
        head.addAttribute('version', '0.1')
        head.writeHeader(self.file)

        ug = XMLHeader('UnstructuredGrid')
        ug.writeHeader(self.file)

        pc = XMLHeader('Piece')
        pc.addAttribute('NumberOfPoints', len(self.mesh.verts))
        pc.addAttribute('NumberOfCells', len(self.mesh.elems))

        pc.writeHeader(self.file)

        self.writePoints()
        self.writeCells()
        self.writePointData()
        self.writeCellData()

        pc.writeFooter(self.file)

        ug.writeFooter(self.file)
        head.writeFooter(self.file)

    def writePoints(self):

        pts = XMLHeader('Points')
        pts.writeHeader(self.file)

        data = XMLHeader('DataArray')
        data.addAttribute('NumberOfComponents', '3')
        data.addAttribute('type', 'Float32')
        data.addAttribute('format', 'ascii')
        data.writeHeader(self.file)

        for p in self.mesh.verts:
            self.file.write('%g %g 0.0\n' % p)

        data.writeFooter(self.file)
        pts.writeFooter(self.file)

    def writeCells(self):

        cells = XMLHeader('Cells')
        cells.writeHeader(self.file)

        conn = XMLHeader('DataArray')
        conn.addAttribute('Name', 'connectivity')
        conn.addAttribute('type', 'Int32')
        conn.addAttribute('format', 'ascii')
        conn.writeHeader(self.file)

        for e in self.mesh.elems:
            self.file.write('%d %d %d\n' % e)

        conn.writeFooter(self.file)

        offsets = XMLHeader('DataArray')
        offsets.addAttribute('type', 'Int32')
        offsets.addAttribute('Name', 'offsets')
        offsets.addAttribute('format', 'ascii')
        offsets.writeHeader(self.file)

        count=0
        numCells = len(self.mesh.elems)
        for c in range(numCells):
            count += 3
            self.file.write('%d\n' % count)


        offsets.writeFooter(self.file)

        types = XMLHeader('DataArray')
        types.addAttribute('type', 'UInt8')
        types.addAttribute('Name', 'types')
        types.addAttribute('format', 'ascii')
        types.writeHeader(self.file)

        for c in range(numCells):
            self.file.write('5\n') # code for triangle elements

        types.writeFooter(self.file)
        cells.writeFooter(self.file)

    def writePointData(self):

        pd = XMLHeader('PointData')
        if len(self.fields)>0:
            pd.addAttribute('Scalars', list(self.fields.keys())[0])
            #print(self.fields.keys())
        pd.writeHeader(self.file)

        for name,field in self.fields.items():

            xml = XMLHeader('DataArray');
            xml.addAttribute('type', 'Float32')
            xml.addAttribute('Name', name)
            xml.addAttribute('format', 'ascii')
            xml.writeHeader(self.file)

            for f in field:
                self.file.write('%g\n' % f)

            xml.writeFooter(self.file)


        pd.writeFooter(self.file)

    def writeCellData(self):

        cd = XMLHeader('CellData')
        cd.writeHeader(self.file)
        cd.writeFooter(self.file)



if __name__=='__main__':

    from .TriangleMeshReader import *
    import numpy as np

    reader = TriangleMeshReader('../Geometry/oneHole.1')
    mesh = reader.getMesh()

    vec = np.zeros(len(mesh.verts))
    for i,p in enumerate(mesh.verts):
        vec[i] = p[0]*p[1]



    
    writer = VTKWriter('test.vtu')
    writer.addMesh(mesh)
    writer.addField('test', vec)
    writer.write()
