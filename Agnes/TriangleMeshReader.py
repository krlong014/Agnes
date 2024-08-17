

# --------------------------------------------------------------------------
# Reader for constructing a mesh from data stored in the file formats
# output from Jonathan Shewchuk's Triangle mesh generator. See the
# Triangle documentation at http://www.cs.cmu.edu/~quake/triangle.html.
#
# The reader expects to find three files:
# --- vertices in a .node file
# --- edges in a .edge file
# --- elements in a .ele file
#
# Note: by default, Triangle does not produce a .edge file. Use the "-e"
# option when running Triangle to produce this information.
#
# If no boundary markers are provided, all edges will be given the label "0".
# If edges have boundary markers, those markers are used as labels.
#
#
# Katharine Long, Sep 2020
# For Math 5344
# --------------------------------------------------------------------------

from .LoadableMesh import LoadableMesh

class TriangleMeshReader:
    # Create a reader object to read a mesh from the files
    # filename.node, filename.edge, and filename.ele
    def __init__(self, filename):
        self.filename = filename
        self.offset = 0

    # Call this function to read the mesh and return it to the user
    def getMesh(self):
        mesh = LoadableMesh()
        self.readVerts(mesh)
        self.readSides(mesh)
        self.readElems(mesh)
        return mesh

    # ---- Functions past this point are for internal use

    # Tokenize a line after stripping comments
    def tokenize(self, line):
        # Strip comments
        h = line.find('#');
        if h>=0:
            line = line[0:h]

        # strip trailing whitespace
        line = line.rstrip()

        if len(line)==0:
            return []
        return line.split()

    # Read the vertices from the .node file
    def readVerts(self, mesh):

        #print('reading verts')

        with open('%s.node' % self.filename) as f:

            lineCount = 0

            for line in f:
                toks = self.tokenize(line)
                if len(toks)==0:
                    continue

                # If this is our first line of content, the number of nodes
                # is the first entry.
                if lineCount==0:
                    nNodes = int(toks[0])
                    #print('nVerts = %d' % nNodes)
                    lineCount +=1
                    continue

                # If this is our first node, use its index to determine whether
                # indices are numbered starting from 0 or from 1.
                if lineCount==1:
                    self.offset = int(toks[0])

                # Add the vertex to the mesh
                v = (float(toks[1 ]), float(toks[2]))
                mesh.addVertex(v)

                lineCount += 1


    # Read the sides from the .edge file
    def readSides(self, mesh):

        #print('reading edges')

        with open('%s.edge' % self.filename) as f:

            first = True

            for line in f:
                toks = self.tokenize(line)
                if len(toks)==0:
                    continue

                if first:
                    nSides = int(toks[0])
                    #print('nSides = %d' % nSides)
                    first = False
                    continue

                s = [int(toks[1])-self.offset, int(toks[2])-self.offset]
                if len(toks)>=4:
                    label = int(toks[3])
                else:
                    label = 0
                mesh.addSide(s[0], s[1], label)



    # Read the elements from the .ele file
    def readElems(self, mesh):

        #print('reading elems')

        with open('%s.ele' % self.filename) as f:

            first = True

            for line in f:
                toks = self.tokenize(line)
                if len(toks)==0:
                    continue

                if first:
                    nElems = int(toks[0])
                    #print('nElems = %d' % nElems)
                    first = False
                    continue

                a = int(toks[1])-self.offset
                b = int(toks[2])-self.offset
                c = int(toks[3])-self.offset
                mesh.addElem(a,b,c)


# ------------------------------------------------------------------------
#
# Produce a set of files for a simple test mesh

def writeExampleTriangleFiles():
    vertexData = """
    4 2 0 0
    0 0 0
    1 1 0
    2 1 1
    3 0 1
    """

    sideData = """
    5 1
    0 0 1 1
    1 1 2 2
    2 2 3 1
    3 0 3 2
    4 0 2 0
    """

    elemData = """
    2 3 0
    0 0 1 2
    1 0 2 3
    """


    tmpName = 'triExample'
    with open('%s.node' % tmpName, 'w') as f:
        f.write(vertexData)

    with open('%s.edge' % tmpName, 'w') as f:
        f.write(sideData)

    with open('%s.ele' % tmpName, 'w') as f:
        f.write(elemData)

    return tmpName


# ---------------------------------------------------------------------------
# Test code

if __name__=='__main__':

    tmpName = writeExampleTriangleFiles()
    print('Writing test mesh to file %s' % tmpName)

    reader = TriangleMeshReader(tmpName)
    mesh = reader.getMesh()
    mesh.dump()
