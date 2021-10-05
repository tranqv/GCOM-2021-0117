#!/usr/bin/env python
# coding: utf-8

#
#   Subroutines for Numerical Examples of the revised manuscript: 
#
#       "On mild solutions of the p-Laplacian fractional Langevin equations with
#       anti-periodic type boundary conditions"
#
#   Usage:
#   from inc_sub2 import * 
#

import numpy as np 
from   mittag_leffler import ml 
from   scipy.special import gamma
import matplotlib.pyplot as plt


MachEps = np.finfo(float).eps

#
# Default values of local parameters:

loc_epsilon = 1.0e-10 
loc_a = 0
loc_b = 10
loc_p = 0.9      
loc_q = 1.0 / ( 1.0 - 1.0/loc_p )
loc_alpha = 0.777
loc_beta  = 0.888
loc_lambda = 0.555 
loc_eta1 = 1.0
loc_eta2 = 1.0 
loc_sigma = 0.333


def parameter_show () :
    print ( 'epsilon =', loc_epsilon )
    print ( 'a       =', loc_a )
    print ( 'b       =', loc_b )
    print ( 'p       =', loc_p )
    print ( 'q       =', loc_q )
    print ( 'alpha   =', loc_alpha )
    print ( 'beta    =', loc_beta )
    print ( 'lambda  =', loc_lambda )
    print ( 'eta1    =', loc_eta1 )
    print ( 'eta2    =', loc_eta2 )
    print ( 'sigma   =', loc_sigma )
    return None 

def parameter_sync ( 
    epsilon, a, b, p, q, alpha, beta, lam, eta1, eta2, sigma ):

    import inc_sub2 as xx 

    xx.loc_epsilon = epsilon 
    xx.loc_a = a
    xx.loc_b = b
    xx.loc_p = p
    xx.loc_q = q
    xx.loc_alpha = alpha 
    xx.loc_beta  = beta 
    xx.loc_lambda = lam
    xx.loc_eta1 = eta1
    xx.loc_eta2 = eta2
    xx.loc_sigma = sigma
    print ( "Parameters were synchronized. Continue." )
    return None 


def f_phi( p, s ):
    return s*abs(s)**(p-2)
def f_phip( s ):
    return f_phi( loc_p, s )
def f_phiq( s ):
    return f_phi( loc_q, s )
#

#
# Mittag-Leffler function:
def f_E ( alpha, beta, t ):
    return ml ( complex(t,0.0), alpha, beta )
#

#
# I_{a+}^{alpha,psi} f (t) 
# = 
# int_{a}^{t} psi'(tau) (psi(t)-psi(tau))**(alpha-1) f(tau,u(tau)) dtau 
#                                                          -----
#                                                          u
#             ----------------------------------------------------
#             f_F = F( t, tau )
#
def f_F ( alpha, t, tau, u, psi, dpsi, f ):
    tmp = abs( psi(t) - psi(tau) )
    return dpsi(tau) * tmp**(alpha-1) * f(tau,u) 
#   

def f_Eps ( alpha, t, tau, lam, psi, dpsi ):
    tmp = abs( psi(t) - psi(tau) )
    return dpsi(tau) * tmp**(alpha-1) * \
           f_E ( alpha, alpha, complex(-lam * tmp**alpha, 0.0) )
#

def f_g ( a, b, eta1, eta2, alpha, lam, psi, t ):
    return -eta2*f_E( alpha, 1.0, complex(-lam*abs(psi(t)-psi(a))**alpha,0.0) )/( \
    eta1  + eta2*f_E( alpha, 1.0, complex(-lam*abs(psi(b)-psi(a))**alpha,0.0) ) )
#

def get_sigma ( sigma, afac=0.9 ):
    return afac*sigma + ( 1.0 - afac )
#

def gene_meshes ( ns, a, b ):
    from numpy   import zeros, linspace  
    tk   = linspace ( a, b, ns )
    xik  = zeros ( ns  )
    dt = ( b - a )/( ns - 1 ) 
    xik[0] = tk[0]
    for k in range( 1, ns ):
        xik[k] = ( tk[k-1] + tk[k] )/2
#
    print ( '\n*** N = %10d' % ( ns ) ) 
    print ( '        dt = %14.7E\n' % ( dt )  )
#
    return dt, tk, xik
#

def calc_G ( N, t, psi ):
    G = np.zeros( N )
    for k in range(0,N):
        G[k] = f_g ( 
                loc_a, loc_b, loc_eta1, loc_eta2, 
                loc_alpha, loc_lambda, psi, t[k] ) 
    return G
#

def calc_E ( N, t, xi, psi, dpsi ):
    a_E = np.zeros( [N,N] )
    for k in range(1,N):
        for j in range(1,k+1):
            a_E[k,j] = f_Eps ( 
                        loc_alpha, t[k], xi[j], 
                        loc_lambda, psi, dpsi )
    return a_E
#


def f_errint ( ti, t1, u1, t2, u2 ):

    from numpy import sqrt 
    from scipy.interpolate import Akima1DInterpolator as scipy_akima
#   from scipy.interpolate import interp1d            as scipy_inter

#   f1_int = scipy_inter ( t1, u1, kind='cubic' )
    f1_int = scipy_akima ( t1, u1 )
    v1_int = f1_int ( ti )

#   f2_int = scipy_inter ( t2, u2, kind='cubic' )
    f2_int = scipy_akima ( t2, u2 )
    v2_int = f2_int ( ti )

    e0 = max( abs(v1_int - v2_int) )
    e2 = sqrt( sum( abs(v1_int - v2_int)**2 ) / len(ti) )

    return e0, e2
#


def f_interp ( t1, u1, t2 ):

    from scipy.interpolate import Akima1DInterpolator as scipy_akima
#   from scipy.interpolate import interp1d            as scipy_inter

#   f1_int = scipy_inter ( t1, u1, kind='cubic' )
    f1_int = scipy_akima ( t1, u1 )
    
    return f1_int ( t2 )
#


def f_errest ( u1, u2 ):

    from numpy import sqrt 

    e0 = max( abs(u1 - u2) )
    e2 = sqrt( sum( abs(u1 - u2)**2 ) / len(u1) )

    return e0, e2
#



def f_Dfrac (  t, u, alpha, psi ):
#
#   W(t) 
#   = ^C D_{a+}^{psi,alpha} u(t) 
#   =   1/gamma(1-alpha) * 
#       int_{a}^{t} ( psi(t)-psi(s) )**(-alpha) u'(s) ds
#                   ---------------------------
#                           v(s) 
#
#   for t = midpoints in ( a, b ), a=t[0], b=t[N-1]
#
#   Input: t[k], u[k]
#   Ouput: w[i]
#
#
    from numpy import zeros 
    from scipy.special import gamma

    N = len (t) 

    if ( N < 1 ):
        print ( "ERROR: no data" )
        return None 

    a = t[0]
    b = t[N-1]
    
    dt = t[1] - t[0]

    du = zeros( N )
    xi = zeros( N )
    wk = zeros( N )

    gx = 1.0 / gamma ( 1 - alpha )

    for i in range ( 1, N ):
        du[i] = ( u[i] - u[i-1] ) 
        xi[i] = ( t[i] + t[i-1] ) / 2 

    for k in range ( 1, N ):
        stmp = 0.0 
        for i in range ( 1, k+1 ):
            stmp = stmp + du[i]*abs( psi(t[k]) - psi(xi[i]) )**(-alpha) 
        wk[k] = stmp 
        
    return xi, wk 
#


def intshift ( num, wid ):
    st = str(num)
    ls = len(st)   
    for i in range(0,wid-ls):
        st = '0' + st 
    return st 
#

def make_paths():

    import os

#   get current directory
    here = os.getcwd()
#   prints parent directory
    ppat = os.path.abspath( os.path.join( here, os.pardir )) 

    OUTDIR = 'xout'
    OUTDIR = os.path.join( ppat, OUTDIR ) 

    FIGDIR = 'figs'
    FIGDIR = os.path.join( ppat, FIGDIR ) 

    if not os.path.exists( OUTDIR ):
        os.makedirs ( OUTDIR )
        print("Directory %s was created" %( OUTDIR ) )

    if not os.path.exists( FIGDIR ):
        os.makedirs ( FIGDIR )
        print("Directory %s was created" %( FIGDIR ) )

    return OUTDIR, FIGDIR
#




