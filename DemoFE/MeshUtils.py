import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from .TriangleMeshReader import *

# --------------------------------------------------------------------------
# Some functions to work with fields defined on triangular meshes.
#
# (*) evalOnMesh(mesh, func) -- produces a vector of function values evaluated
#     at the mesh vertices.
#
# (*) hMesh(mesh) -- computes the average side length in the mesh
#
# (*) showField(mesh, vec) -- makes a density plot of a function
#     discretized on a mesh
#
# (*) hLocal(mesh) -- produces a vector containing the average length
#                     of the edges connected to each vertex
#
#
# Katharine Long, Sep 2020
# For Math 5344
# --------------------------------------------------------------------------

# --- evalOnMesh -- produce a vector of function values at mesh vertices
#
# Input:
# (*) Argument "mesh" should have an iterable attribute verts, each item
#     in which should be an (x,y) pair
# (*) Argument "func" should be callable with a point given as an (x,y)
#     pair.
def evalOnMesh(mesh, func):

    fVals = np.array([func(v) for v in mesh.verts])

    return fVals

# --- hMesh -- find the average edge length in the mesh
#
# Input:
# (*) Argument "mesh" should have an iterable, indexable attribute verts,
#     each item in which should be an (x,y) pair, and an iterable
#     attribute sides, each item is which a pair of vertex indices.
# Output:
# (*) The average edge length
#
def hMesh(mesh):
    hAvg = 0.0
    for s in mesh.sides:
        A = np.array(mesh.verts[s[0]])
        B = np.array(mesh.verts[s[1]])
        r = A-B
        hEdge = np.sqrt(np.dot(r,r))
        hAvg += hEdge
    hAvg /= len(mesh.sides)
    return hAvg

# --- showField -- do a density plot of a function discretized on a mesh
#
# Input:
# (*) Argument "mesh" should have an iterable attribute "verts" containing
#     (x,y) coordinates of the vertex positions, and an iterable attribute
#     "elems" listing the vertex index triplets for each triangle.
# (*) Argument "vec" is the data to be plotted
# (*) colorBar: whether or not to plot a colorbar
# (*) meshLines: whether or not to show the mesh edges
# (*) contours: number of contour levels
# (*) linewidth: size of the edge lines
# (*) linecolor: color for edges; use the color codes from plot()
# (*) aspect: aspect ratio for plot
def showField(mesh, vec, meshLines=True, colorBar=True, contours=128,
    linewidth=0.4, linecolor='k', aspect='equal'):

    XY = np.array(mesh.verts)

    X = XY[:,0]
    Y = XY[:,1]

    tri = mtri.Triangulation(X, Y, mesh.elems)
    interp = mtri.LinearTriInterpolator(tri, vec)

    fig, ax = plt.subplots()
    cnt = ax.tricontourf(tri, vec,levels=contours)
    if meshLines:
        ax.triplot(tri,'%s-'%linecolor, linewidth=linewidth)
    if colorBar:
        fig.colorbar(cnt)

    ax.set_aspect(aspect)

    return fig

# --- hLocal -- find the average edge length connected to each vertex
#
# Input:
# (*) Argument "mesh" should have an iterable, indexable attribute verts,
#     each item in which should be an (x,y) pair, and an iterable
#     attribute sides, each item is which a pair of vertex indices.
def hLocal(mesh):

    hAvg = np.zeros(len(mesh.verts))
    vertEdgeCount = np.zeros(len(mesh.verts))

    for s in mesh.sides:
        a = s[0]
        b = s[1]
        A = np.array(mesh.verts[a])
        B = np.array(mesh.verts[b])
        r = A-B
        hEdge = np.sqrt(np.dot(r,r))
        hAvg[a] += hEdge
        hAvg[b] += hEdge
        vertEdgeCount[a] += 1
        vertEdgeCount[b] += 1

    hAvg /= vertEdgeCount
    return hAvg

## ---- Test

if __name__=='__main__':

    level = 5
    reader = TriangleMeshReader('./Meshes/triExample.%d' % level)
    mesh = reader.getMesh()

    pi = np.pi

    F = np.array([np.sin(2*pi*p[0])*np.cos(3.0*pi*p[1]) for p in mesh.verts])
    #F = hLocal(mesh)


    fig = showField(mesh, F, meshLines=True, linewidth=0.08, linecolor='k')
    plt.savefig('fig.png')
    #plt.show()
