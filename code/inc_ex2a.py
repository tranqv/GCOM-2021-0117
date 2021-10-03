#
#   Example 4.2
#
#   Revised manuscript: 
#       "On mild solutions of the p-Laplacian fractional Langevin equations with
#       anti-periodic type boundary conditions"
#
#   The exact solution u(t) = ( ln t )^(2.34)
#



import numpy as np 
from   scipy.special import gamma

# epsilon:
c_epsilon = 1.0e-10 

# a, b:
c_a = 1
c_b = 5


#  
#   For p, q > 1 such that 
#   1/p + 1/q = 1
#  

# p, q: 
c_p = 2
c_q = 1.0 / ( 1.0 - 1.0/c_p )

# alpha, beta, lambda, eta1, eta2:
c_alpha = 0.78
c_beta  = 0.89
c_lambda = 1.23 
c_eta1 = 1.0
c_eta2 = 0.0 

# sigma:
c_sigma = 0.6 

# r, k, theta:
c1_r = 2.34
c1_k = 0.67 
c1_theta = 0.56


#print ( 'r - alpha - beta  = ', c1_r - c_alpha - c_beta ) 
#print ( 'theta * r - sigma = ', c1_theta * c1_r - c_sigma )
#print ( '            sigma = ', c_sigma )

#
#

c2_gam1 = gamma( c1_r + 1 )
c2_gam2 = gamma( c1_r - c_alpha + 1 - c_beta ) 
c2_gam3 = gamma( c1_r + 1 - c_beta )

c3_gam1 = c2_gam1 / c2_gam2
c3_gam2 = ( c2_gam1 / c2_gam3 ) * c_lambda 


#
# Uexact(t) = ( ln t )^(2.34):
def f_uexact (t) :
    return abs( np.log(t) )**c1_r 
#

# psi(t) = ln(t):
def f_psi(t):
    return np.log(t)
#

# psi'(t) = 1/t:
def f_dpsi(t):
    return 1.0/t
#

# f(t,psi,u):
def f_f ( s, u ):
    t = np.log( s ) 
    f = c1_k * abs(u)**c1_theta / ( abs(t)**c_sigma ) + \
        c3_gam1 * abs(t)**( c1_r - c_alpha - c_beta ) + \
        c3_gam2 * abs(t)**( c1_r - c_beta ) - \
        c1_k    * abs(t)**( c1_theta*c1_r - c_sigma ) 
    return f
#

#
# Some further checking:
def check_bc ( N, U ):
    print ( "\t\tChecking BC: = 0 = ", c_eta1*U[0] + c_eta2*U[N-1] )
    return None 
#
#
