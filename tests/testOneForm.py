from DemoFE import *
import numpy as np
import numpy.linalg as npla

def test_ConstOneForm():

  print('testing constant one form')
  tri = Triangle((-1/4,1/8),(5/4,0),(1,4/3))
  c = np.pi

  bExact = np.array((c,c,c)) 
  bExact = 21.0/64.0 * bExact

  oneForm = ConstCoeffOneForm(c)
  b = oneForm.localVec(tri, (0,1,2))

  error = npla.norm(b - bExact)/npla.norm(bExact)

  print('error = {:12.5g}'.format(error))

  tol = 1.0e-14

  if error <= tol:
    print('test_ConstOneForm() PASSED!')
  else:
    print('test_ConstOneForm() FAILED!')

  assert(error <= tol)



def test_VarOneForm():
  print('testing variable one form')

  tri = Triangle((-1/4,1/8),(5/4,0),(1,4/3))

  def f(x, y, discFuncs):
      return np.exp(x + y)


  # Exact soln for comparison, computed in Mathematica
  exactB = np.array((0.842898747945908,
                    1.1466748135442901,
                    1.5345857052411658))
  # Upper bound on error, computed in Mathematica
  def bound(n):
    return np.exp(2.91972240007496 - 2.22659575688739 * n)


  allGood = True

  for n in range(1,7):

    print('n={}'.format(n))
    quad = GaussRule(n)

    oneForm = VarCoeffOneForm(quad, f)

    b = oneForm.localVec(tri, (0,1,2))

    print('\tlocal vector = ', b)

    err = exactB - b
    relErr = np.linalg.norm(err)/np.linalg.norm(exactB)
    stat = 'PASSED'
    if relErr >= bound(n):
      allGood = False
      stat = 'FAIL'
    print('error norm = {:12.5g}, status={}'.format(relErr, stat))

  if allGood:
    print('------- All cases PASSED! ')
  else:
    print('------- FAILURES DETECTED! ')

  assert(allGood)

    

if __name__=='__main__':

  test_ConstOneForm()

  test_VarOneForm()


  


  


