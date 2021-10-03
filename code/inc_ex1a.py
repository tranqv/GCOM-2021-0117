#
#   Example 4.1 (a) 
#
#   Revised manuscript: 
#       "On mild solutions of the p-Laplacian fractional Langevin equations with
#       anti-periodic type boundary conditions"
#
#   The exact solution u(t) = 0.
#

import numpy as np 


# epsilon = 10^(-10):
c_epsilon = 1.0e-10 

# a=0, b=1:
c_a = 0
c_b = 1


#  
#   For p, q > 1 such that 
#   1/p + 1/q = 1
#  

# p = 5/3, q = 1/( 1 - 1/p ):
c_p = 5.0 / 3.0
c_q = 1.0 / ( 1.0 - 1.0/c_p )

# alpha = 5/8, beta = 4/5, lambda = 1/2, eta1 = eta2 = 1:
c_alpha = 5.0 / 8
c_beta  = 4.0 / 5
c_lambda = 1.0 / 2 
c_eta1 = 1.0
c_eta2 = 1.0 

# sigma = 1/3:
c_sigma = 1.0/3

#
# U_exact = 0:
def f_uexact ( t ) :
    return np.zeros( len(t) ) 

#
# psi(t) = t:
def f_psi(t):
    return t
#

# psi'(t) = 1:
def f_dpsi(t):
    return 1
#

# f(t,psi,u):
def f_f ( t, u ):
    return  ( abs(t)**(-1.0/3) / 10 ) * abs(u) / ( 1 + abs(u) )
#


#
# Some further checking:
def check_bc ( N, U ):
    print ( "\t\tChecking BC: U(0) + U(1) = 0 = ", c_eta1*U[0] + c_eta2*U[N-1] )
    return None 
#
#
