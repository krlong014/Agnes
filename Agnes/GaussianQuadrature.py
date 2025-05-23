import numpy as np
from QuadratureRule import QuadratureRule
from QuadratureFamily import QuadratureFamily

class GaussianQuadratureFamily(QuadratureFamily):

	def __init__(self, order):
		super().__init__(name='GaussianQuadrature', order=order)

	def getRule(self, dim):

		if dim==1:
			# If in 1D, the order is p=2n-1; solve for n to find n=(p+1)/2.
			# No 1D Gauss rule can have even order, so to get order of accuracy
			# p=2k we have to use the rule for p=2k+1. 
			p = self.order()
			if p%2 == 0:
				p = p+1
			nPts = (p+1)//2
			assert(nPts in _gauss1D.keys())

			rule = _gauss1D[nPts]
		
		elif dim==2:
			# We've already tabulated 2D rules by order, so simply look up 
			# the appropriate entry from the dictionary
			assert(self.order() in _gauss2D.keys())
			rule = _gauss2D[self.order()]
		else:
			raise RuntimeError('No quadrature rules available for dim={}'.format(dim))
		
		return QuadratureRule(rule['w'], rule['x'], self.order(), dim,
												name='Gauss(p={},d={})'.format(self.order(),dim))


def main():
	
	maxP1D = 9
	maxP2D = 6
	tol = 1.0e-13

	print('Testing 1D Gauss rules')
	for p in range(1, maxP1D+1):
		gauss = GaussianQuadratureFamily(p)
		passed, err = gauss.test(1, tol)
		if passed:
			stat = 'passed'
		else:
			stat = 'FAILED'

		print('\tGauss 1D(p={}) {} with error {:10.5g}'.format(p,stat,err))



	print('Testing 2D Gauss rules')

	for p in range(1, maxP2D+1):
		gauss = GaussianQuadratureFamily(p)
		passed, err = gauss.test(2, tol)
		if passed:
			stat = 'passed'
		else:
			stat = 'FAILED'

		print('\tGauss 2D(p={}) {} with error {:10.5g}'.format(p,stat,err))
			
			
	

# ####################################################################
#
#         Tables of Gauss points
#
# ####################################################################


# --------------------------------------------------------------------
# The 1D Gauss-Legendre points through order 9 can be computed exactly.
# In finite element applications that's usually all we'll need. These 
# points are for the reference element [-1,1]
# --------------------------------------------------------------------

_gauss1D = {
	1 : {  # 1-point, order p=1
		'x' : (0,),
		'w' : (2,)
	},

	2 : { # 2-point, order p=3
		'x' : (-1/np.sqrt(3), 1/np.sqrt(3)),
		'w' : (1,1)
	},

	3 : { # 3-point, order p = 5
		'x' : (-np.sqrt(3/5), 0, np.sqrt(3/5)),
		'w' : (5/9, 8/9, 5/9) 
	},

	4 : { # 4-point, order p=7
		'x' : (
				-np.sqrt(3/7 + 2/7*np.sqrt(6/5)),
				-np.sqrt(3/7 - 2/7*np.sqrt(6/5)),
				np.sqrt(3/7 - 2/7*np.sqrt(6/5)),
				np.sqrt(3/7 + 2/7*np.sqrt(6/5))
			),
		'w' : (
				(18 - np.sqrt(30))/36, 
				(18 + np.sqrt(30))/36,
				(18 + np.sqrt(30))/36,
				(18 - np.sqrt(30))/36
			)
	},

	5 : { # 5-point, order p=9
		'x' : (
			-1/3*np.sqrt(5 + 2*np.sqrt(10/7)),
			-1/3*np.sqrt(5 - 2*np.sqrt(10/7)),
			0,
			1/3*np.sqrt(5 - 2*np.sqrt(10/7)),
			1/3*np.sqrt(5 + 2*np.sqrt(10/7))
		),
		'w' : (
			(322 - 13*np.sqrt(70))/900,
			(322 + 13*np.sqrt(70))/900,
			128/225,
			(322 + 13*np.sqrt(70))/900,
			(322 - 13*np.sqrt(70))/900
		)
	}
}


# --------------------------------------------------------------------
# Tables of 2D Gauss points on triangles from Strang & Fix, obtained
# from John Burkhardt's tabulation of quadrature rules. These rules
# are for the reference triangle ((0,0), (1,0), (0,1)).
# --------------------------------------------------------------------

_gauss2D = {
	1: {
		'x' : ((1/3, 1/3),),
		'w' : (1.0,)
	},

	2: {
		'x' : ((2/3, 1/6),(1/6, 2/3), (1/6,1/6)),
		'w' : (1/3,1/3,1/3)    
	},

	3: {
		'x': [(0.659027622374092, 0.231933368553031),
					(0.659027622374092, 0.109039009072877),
					(0.231933368553031, 0.659027622374092),
					(0.231933368553031, 0.109039009072877),
					(0.109039009072877, 0.659027622374092),
					(0.109039009072877, 0.231933368553031)],
		'w': [0.16666666666666666,
					0.16666666666666666,
					0.16666666666666666,
					0.16666666666666666,
					0.16666666666666666,
					0.16666666666666666]
		},

	4: {
		'x': [(0.816847572980459, 0.091576213509771),
					(0.091576213509771, 0.816847572980459),
					(0.091576213509771, 0.091576213509771),
					(0.10810301816807, 0.445948490915965),
					(0.445948490915965, 0.10810301816807),
					(0.445948490915965, 0.445948490915965)],
		'w': [0.109951743655322,
					0.109951743655322,
					0.109951743655322,
					0.223381589678011,
					0.223381589678011,
					0.223381589678011]
		},
	5: {
			'x': [(0.3333333333333333, 0.3333333333333333),
					(0.7974269853530872, 0.10128650732345633),
					(0.10128650732345633, 0.7974269853530872),
					(0.10128650732345633, 0.10128650732345633),
					(0.05971587178976981, 0.47014206410511505),
					(0.47014206410511505, 0.05971587178976981),
					(0.47014206410511505, 0.47014206410511505)],
			'w': [0.225,
					0.12593918054482717,
					0.12593918054482717,
					0.12593918054482717,
					0.13239415278850616,
					0.13239415278850616,
					0.13239415278850616]
		},
	6: {'x': [(0.873821971016996, 0.063089014491502),
					(0.063089014491502, 0.873821971016996),
					(0.063089014491502, 0.063089014491502),
					(0.501426509658179, 0.24928674517091),
					(0.24928674517091, 0.501426509658179),
					(0.24928674517091, 0.24928674517091),
					(0.636502499121399, 0.310352451033785),
					(0.636502499121399, 0.053145049844816),
					(0.310352451033785, 0.636502499121399),
					(0.310352451033785, 0.053145049844816),
					(0.053145049844816, 0.636502499121399),
					(0.053145049844816, 0.310352451033785)],
		'w': [0.050844906370207,
					0.050844906370207,
					0.050844906370207,
					0.116786275726379,
					0.116786275726379,
					0.116786275726379,
					0.082851075618374,
					0.082851075618374,
					0.082851075618374,
					0.082851075618374,
					0.082851075618374,
					0.082851075618374]}
	}



if __name__=='__main__':

	main()

