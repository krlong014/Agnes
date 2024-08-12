from abc import ABC, abstractmethod

class BaseMesh:
    # Constructor has nothing to do
    def __init__(self):
        pass

    # Iterable view of elements
    @abstractmethod
    def elements(self):
        pass

    # Iterable view of sides having a given label
    @abstractmethod
    def elements(self):
        pass


    # Look up the coordinates of a vertex
    @abstractmethod
    def getVertexCoords(self, vertIndex):
        pass

    # Look up the indices of a cell's facets of a specified dimension
    @abstractmethod
    def getFacets(self, cellDim, cellIndex, facetDim):
        pass

    # Look up the elements connected to a specified cell
    @abstractmethod
    def getConnectedElems(self, cellDim, cellIndex):
        pass



    def refToPhysMap(self, points):
