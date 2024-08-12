import matplotlib.pyplot as plt
import numpy as np
from math import sqrt
from LoadableMesh import *

class MPLMeshViewer:

    def __init__(self, fontSize=16, vertRad=0.06,
            colors=['black','blue','red','seagreen', 'orangered','dodgerblue',
                'forestgreen']):
        self.fontSize = fontSize
        self.vertRad = vertRad
        self.colors = colors

    def show(self, mesh):
        ax = plt.axes()

        for i,v in enumerate(mesh.verts):
            self.addVert(ax, i, v)

        for i,e in enumerate(mesh.elems):
            self.addElem(ax, mesh, i, e)

        ax.set_aspect('equal')
        ax.autoscale_view()

        plt.show()


    def addVert(self, ax, index, xy):
        circle = plt.Circle(xy, radius=self.vertRad, fc='lightgray', ec='black')
        ax.add_patch(circle)
        ax.annotate(index.__str__(), xy, ha='center',
            va='center', fontsize=self.fontSize)

    def addElem(self, ax, mesh, index, abc):
        A = np.array(mesh.verts[abc[0]])
        B = np.array(mesh.verts[abc[1]])
        C = np.array(mesh.verts[abc[2]])

        pairs = [[abc[0], abc[1]], [abc[1], abc[2]], [abc[2], abc[0]]]

        for p in pairs:
            p.sort()
            q = tuple(p)
            label = mesh.getSideLabel(q)
            col = self.colors[label % len(self.colors)]

            P1 = np.array(mesh.verts[q[0]])
            P2 = np.array(mesh.verts[q[1]])
            P12 = P2-P1
            PUnit = P12/sqrt(np.dot(P12,P12))
            p1 = P1 + self.vertRad*PUnit
            p2 = P2 - self.vertRad*PUnit
            line = plt.Line2D((p1[0],p2[0]), (p1[1], p2[1]), color=col)
            ax.add_line(line)

        # put label at centroid of element
        centroid = (A+B+C)/3.0

        ax.annotate(index.__str__(), centroid, ha='center',
            va='center', fontsize=self.fontSize)


if __name__=='__main__':

    from TriangleMeshReader import *
    reader = TriangleMeshReader('../Geometry/oneHole.1')
    mesh = reader.getMesh()

    viewer = MPLMeshViewer(vertRad=0.05, fontSize=8)
    viewer.show(mesh)
