import numpy as np
import numpy.linalg as npla

def makeStiffnessMatrix(N, h):
    '''
    Form the stiffness matrix for the negative Laplacian operator on
    a 1D uniform mesh on [0,1] with homogeneous Dirichlet BCs and
    N internal nodes.
    '''
    K = np.zeros([N,N])

    for i in range(N): # i in (0,1,2,...,N-1)
        K[i,i] = 2/h
        if i>0:
            K[i,i-1] = -1/h
        if i<N-1:
            K[i,i+1] = -1/h

    return K

def makeLoadVector(N, h, func):
    '''Integrate func(x)phi(x) by quadrature'''

    b = np.zeros(N)

    # Loop over elements (numbered 0 through N)
    for e in range(0, N+1):
        localB = localLoadVector(e, h, func)
        #print('elem={}, local vec={}'.format(e,localB))
        if e==0:
            #print('entry {} going into b_0'.format(localB[1]))
            b[0] -= localB[1]
        elif e==N:
            #print('entry {} going into b_(N-1)'.format(localB[0]))
            b[N-1] -= localB[0]
        else:
            #print('entries {} going into b_{}, b_{}'.format(localB, e-1,e))
            b[e-1:e+1] -= localB
        #print('global load vector is {}'.format(b))

    return b


def localLoadVector(e, h, func):

    squiggle1 = -1/np.sqrt(3)
    squiggle2 = 1/np.sqrt(3)

    x_e = e*h
    x1 = x_e + (squiggle1+1)*h/2
    x2 = x_e + (squiggle2+1)*h/2

    phi_e_at_x1 = (x1-x_e)/h
    phi_e_at_x2 = (x2-x_e)/h

    phi_e_plus_1_at_x1 = ((x_e+h)-x1)/h
    phi_e_plus_1_at_x2 = ((x_e+h)-x2)/h

    f1 = func(x1)
    f2 = func(x2)
    #print('f1={}, f2={}'.format(f1,f2))
    #print('phi_e={},{}'.format(phi_e_at_x1, phi_e_at_x2))
    #print('h={}'.format(h))


    b0 = h/2*(phi_e_at_x1*f1 + phi_e_at_x2*f2)
    b1 = h/2*(phi_e_plus_1_at_x1*f1 + phi_e_plus_1_at_x2*f2)
    localB = np.array([b0, b1])
    #print('in localLoadVector: localB = {}'.format(localB))
    return localB


def integrateFunction(N, h, func):
    '''Integrate the function func(x) using 2-point Gaussian quadrature on
    each element.'''

    rtn = 0.0

    for e in range(0, N+1):
        squiggle1 = -1/np.sqrt(3)
        squiggle2 = 1/np.sqrt(3)

        x_e = e*h
        x1 = x_e + (squiggle1+1)*h/2
        x2 = x_e + (squiggle2+1)*h/2

        rtn += h/2*(func(x1) + func(x2))

    return rtn

class InterpolatedFunction:

    def __init__(self, h, vec):
        self.h = h
        self.vec = vec

    def __call__(self, x):
        h = self.h
        k = int(np.floor(x/h))

        if k < 0:
            return self.vec[0]
        if k >= len(self.vec)-1:
            return self.vec[len(self.vec)-1]

        x_k = h*k
        x_k1 = h*(k+1)

        # This is vec[k+1]*phi_{k+1}(x) + vec[k]*phi_k(x)
        return self.vec[k+1]*(x-x_k)/h + self.vec[k]*(x_k1-x)/h

def fConst(x):
    return 1

if __name__=='__main__':

    N = 127
    h = 1/(N+1)

    K = makeStiffnessMatrix(N, h)

    #print('{}-by-{} matrix, h={}, K=\n{}'.format(N, N, h, K))

    b = makeLoadVector(N, h, fConst)

    #print('load vector, h={}, b=\n{}'.format(h, b))

    #print('-'*60)

    # Solve system
    uTrunc = npla.solve(K, b)

    # Put boundary values back in the solution
    uVec = np.zeros(N+2)
    uVec[1:N+1] = uTrunc
    #print('augmented solution = ', uVec)

    # Create a function that interpolates solution values between the nodes
    uSoln = InterpolatedFunction(h, uVec)

    import matplotlib.pyplot as plt

    #plt.plot(X, U, 'o-')
    #plt.show()


    def uEx(x):
        return 0.5*x*(x-1)

    def errSq(x):
        return (uEx(x) - uSoln(x))**2


    errNorm = np.sqrt(integrateFunction(N, h, errSq))

    print('h={}, ||err||={}'.format(h, errNorm))


    X = np.linspace(0,1,10)
    U = np.array([uSoln(x) for x in X]) # list comprehension to form array
    UEx = np.array([uEx(x) for x in X]) # list comprehension to form array

    #plt.plot(X, U, 'o-r', X, UEx, 'o-b')
    #plt.show()
