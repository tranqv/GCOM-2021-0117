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
#   Applying Richardson extrapolation method to approximate f'(x) 
#
#   + For x well inside the domain:
#
def deri_d1f_c ( \
    f, x, eps=1.0e-7, maxit=20, delta=0.2, report=False, dref=0 ):
#
    if ( maxit == 1 ):
        return ( f(x+delta) - f(x-delta) )/delta/2
    df = 0
    d2 = 0
    d1 = 0
    dx = delta 
    if ( report ):
        print ( "   x = %12.5e" % ( x ) )
        print ( "dref = %12.5e" % ( dref) )
        print ( "%4s %12s %12s %12s %12s" %( \
                "i", "dx =", "d2-dref =", "df-dref =", "errt =" ) ) 
    for i in range( maxit ):
        df_old = df 
        d1 = d2
        d2 = ( f(x+dx) - f(x-dx) )/dx/2
        df = ( 4*d2 - d1 )/ 3 
        errt = abs( df - df_old )/( max( abs(df), 1.0 ) )
        if ( report ):
            print ( "%4d %12.5e %12.5e %12.5e %12.5e" %( \
                    i, dx, d2-dref, df-dref, errt ) ) 
        if ( errt < eps ):
            break 
        dx = dx / 2 
    return df 
#
#   + For x near the left boundary:
#
def deri_d1f_l ( \
    f, x, eps=1.0e-7, maxit=20, delta=1.0e-1, report=False, dref=0 ):
#
    if ( maxit == 1 ):
        return (-3*f(x) + 4*f(x+delta) - f(x+2*delta))/delta/2
    df = 0
    d2 = 0
    d1 = 0
    dx = delta 
    if ( report ):
        print ( "   x = %12.5e" % ( x ) )
        print ( "dref = %12.5e" % ( dref) )
        print ( "%4s %12s %12s %12s %12s" %( \
                "i", "dx =", "d2-dref =", "df-dref =", "errt =" ) ) 
    for i in range( maxit ):
        df_old = df 
        d1 = d2
        d2 = ( 3*(f(x+dx) - f(x)) + (f(x+dx) - f(x+2*dx)) )/dx/2 
        df = ( 4*d2 - d1 )/ 3 
        errt = abs( df - df_old )/( max( abs(df), 1.0 ) )
        if ( report ):
            print ( "%4d %12.5e %12.5e %12.5e %12.5e" %( \
                    i, dx, d2-dref, df-dref, errt ) ) 
        if ( errt < eps ):
            break 
        dx = dx / 2 
    return df 
#
#   + For x near the right boundary:
#
def deri_d1f_r ( \
    f, x, eps=1.0e-7, maxit=20, delta=0.2, report=False, dref=0 ):
#
    if ( maxit == 1 ):
        return -(-3*f(x) + 4*f(x-delta) - f(x-2*delta) )/delta/2
    df = 0
    d2 = 0
    d1 = 0
    dx = delta 
    if ( report ):
        print ( "   x = %12.5e" % ( x ) )
        print ( "dref = %12.5e" % ( dref) )
        print ( "%4s %12s %12s %12s %12s" %( \
                "i", "dx =", "d2-dref =", "df-dref =", "errt =" ) ) 
    for i in range( maxit ):
        df_old = df 
        d1 = d2
        d2 =-( 3*(f(x-dx) - f(x)) + (f(x-dx)- f(x-2*dx)) )/dx/2 
        df = ( 4*d2 - d1 )/ 3 
        errt = abs( df - df_old )/( max( abs(df), 1.0 ) )
        if ( report ):
            print ( "%4d %12.5e %12.5e %12.5e %12.5e" %( \
                    i, dx, d2-dref, df-dref, errt ) ) 
        if ( errt < eps ):
            break 
        dx = dx / 2 
    return df 
#

#
# Mittag-Leffler function:
def f_E ( alpha, beta, t ):
    return ml ( complex(t,0.0), alpha, beta )
#
# and its derivative:
def f_dE ( alpha, beta, t, opt=0 ):    
    def ftmp (s):
        return f_E ( alpha, beta, s )
    dE = 0.0 
    if ( opt == 0 ): # t well inside the domain (default mode)
        dE = deri_d1f_c ( ftmp, t, eps=loc_epsilon/10, delta=0.1, report=False )
    elif( opt > 0 ): # t at the left boundary 
        dE = deri_d1f_l ( ftmp, t, eps=loc_epsilon/10, delta=0.1, report=False )
    else:            # t at the right boundary
        dE = deri_d1f_r ( ftmp, t, eps=loc_epsilon/10, delta=0.1, report=False )
    return dE
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

#############################################################################
#
#   Advanced checking 
#
#
#   ^C D_{a+}^{psi,alpha} u(t) 
#
#   =   1/gamma(1-alpha) * 
#       int_{a}^{t} ( psi(t)-psi(s) )**(-alpha) u'(s) ds
#                   ---------------------------
#                           v(s) 
#
#   We wish to check:
#
#   +   eta1* u(a)                       + eta2* u(b)                       = 0
#   +   eta1* ^C D_{a+}^{psi,alpha} u(a) + eta2* ^C D_{a+}^{psi,alpha} u(b) = 0
#
#   But, 
#   ^C D_{a+}^{psi,alpha} u(a) = ??? 
#
#   =   1/gamma(1-alpha) * 
#       lim_{ d -> 0 }
#       int_{a}^{a+d} ( psi(a+d)-psi(s) )**(-alpha) u'(s) ds
#
#   =   1/gamma(1-alpha) * 
#       lim_{ d -> 0 }
#       d * ( psi(a+d)-psi(a+d/2) )**(-alpha) * ( [u(a+d) - u(a)]/d +  O(d) ) 
#
#   =   1/gamma(1-alpha) * 
#       lim_{ d -> 0 }
#       ( psi(a+d)-psi(a+d/2) )**(-alpha) * ( [u(a+d) - u(a)] + d*O(d) )
#
#   Linear interpolation:
#
#   u(a+d) = u[0]*(dt-d)/dt  + u[1]*d/dt
#
#   ^C D_{a+}^{psi,alpha} u(a) = ??? 
#   =   1/gamma(1-alpha) * 
#       ( psi(a+d)-psi(a+d/2) )**(-alpha) * 
#       [ u[0]*(dt-d)/dt  + u[1]*d/dt - u[0] ]
#   =   1/gamma(1-alpha) * 
#       ( psi(a+d)-psi(a+d/2) )**(-alpha) * 
#       [ u[1] - u[0] ]*d/dt
#   ???
#
#
#

#
def f_checkBC ( t, u, alpha, eta1, eta2, psi ):

    from numpy import zeros 
    from scipy.special import gamma

    N = len(t)

    a = t[0]
    b = t[N-1]

    bc1 = eta1*u[0] + eta2*u[N-1]

    gam = 1.0 / gamma( 1 - alpha )

    dt = t[1] - t[0]

    Sb = 0 
    for i in range(1,N):
        s  = ( t[i-1] + t[i] )/2
        v  = abs( psi(b) - psi(s) )**(-alpha)
        du = u[i] - u[i-1]
        Sb = Sb + du*v
    dub = Sb * gam 


#
#   ^C D_{a+}^{psi,alpha} u(a) = ??? 
#   =   1/gamma(1-alpha) * 
#       ( psi(a+d)-psi(a+d/2) )**(-alpha) * 
#       [ u[0]*(dt-d)/dt  + u[1]*d/dt - u[0] ]
#   =   1/gamma(1-alpha) * 
#       ( psi(a+d)-psi(a+d/2) )**(-alpha) * 
#       [ u[1] - u[0] ]*d/dt
#
#   d   = dt    #1.0e-9
#   dua = abs( psi(a+d) - psi(a+d/2) )**(-alpha) *( u[1] - u[0] )*d/dt 

    d   = dt   
    du  = ( -3*u[0] + 4*u[1] - u[2] )/dt/2
    dua = abs( psi(a+d) - psi(a+d/2) )**(-alpha) * du * d 
    dua = dua * gam 


#
    bc2 = eta1*dua + eta2*dub 

    return bc1, bc2 




def f_Dfrac (  t, u, alpha, psi ):
#
#   Advanced checkings ...
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




