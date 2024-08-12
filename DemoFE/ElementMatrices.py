import numpy as np

class Triangle:
    def __init__(self, A, B, C):
        self.A = A
        self.B = B
        self.C = C
        self.Jt = np.array([ [B[0]-A[0], B[1]-A[1]], [C[0]-A[0], C[1]-A[1]]])
        self.detJ = self.Jt[0,0]*self.Jt[1,1]-self.Jt[0,1]*self.Jt[1,0]
        self.JtInv = np.array([[self.Jt[1,1], -self.Jt[0,1]],
            [-self.Jt[1,0],self.Jt[0,0]]])/self.detJ


class LocalLaplacian:
    gradPhi = np.array([[-1,1,0],[-1,0,1]])

    def localMat(self, tri):
        grad = np.matmul(tri.JtInv, LocalLaplacian.gradPhi)
        print('grad=\n', grad)
        print('grad.t=\n', np.transpose(grad))
        return 0.5 * np.abs(tri.detJ) * np.matmul(np.transpose(grad), grad)

class LocalMass:
    M = np.array([[2,1,1],[1,2,1],[1,1,2]])/4.0

    def localMat(self, tri):
        return 0.5 * np.abs(tri.detJ) * LocalMass.M

class LocalLoad:

    def localVec(self, tri, func):
        funcAB = func(0.5*(tri.A + tri.B))
        funcBC = func(0.5*(tri.B + tri.C))
        funcAC = func(0.5*(tri.A + tri.C))

        return 0.5/6.0 * np.abs(detJ) * np.array(
            [
                funcAB + funcAC,
                funcAB + funcBC,
                funcAC + funcBC
            ]
        )


if __name__=='__main__':

    T = Triangle((1,1), (2,1.5),(1.5,2))
    print('Jt=\n', T.Jt)
    print('det(J)=', T.detJ)

    print('Jt*inv(Jt)=\n', np.matmul(T.Jt, T.JtInv))

    lapl = LocalLaplacian()
    print('lapl=\n', lapl.localMat(T))

    mass = LocalMass()
    print('mass=\n', mass.localMat(T))

    def loadFunc(xy):
        return 1

    load = LocalLoad()
    print('load=\n', load.localVec(T, loadFunc))
