#
#   Example 4.3 (a) 
#       without random noise to input 
#
#   Revised manuscript: 
#       "On mild solutions of the p-Laplacian fractional Langevin equations with
#       anti-periodic type boundary conditions"
#
#   Exact solution is unknown in advance.
#

import numpy as np 

# epsilon:
c_epsilon = 1.0e-10 

# a, b:
c_a = 1
c_b = np.exp(1.0) 


#  
#   For p, q > 1 such that 
#   1/p + 1/q = 1
#  

# p, q: 
c_p = 9.0 / 5.0
c_q = 1.0 / ( 1.0 - 1.0/c_p )

# alpha, beta, lambda, eta1, eta2:
c_alpha = 3.0 / 5
c_beta  = 2.0 / 3
c_lambda = 1.0 / 5 
c_eta1 = 1.0
c_eta2 = 2.0 

# sigma:
c_sigma = 0.5 

#
# psi(t):
def f_psi(t):
    return np.log(t)
#

# psi'(t):
def f_dpsi(t):
    return 1/t
#

# f(t,psi,u):
def f_f ( t, u ):
    psi = abs(f_psi(t))
    return 2*abs(psi)**(-0.5)  + abs(psi)**(1.0/3) * abs(u)**(4.0/7)
#

#
# Some further checking:
def check_bc ( N, U ):
    Ua = U[0] 
    Ub = U[N-1]
    Bc = 2*Ub + Ua
    print ( "Checking BC: U(1) + 2U(e) = 0 = ", Bc )
    return None 
#
#
