#!/usr/bin/env python
# coding: utf-8

#
#   Numerical Examples (Main program for Examples 1a, 1b, 2, 3a)
#
#   Revised manuscript: 
#       "On mild solutions of the p-Laplacian fractional Langevin equations with
#       anti-periodic type boundary conditions"
#

import os
import numpy as np 
from   scipy.special import gamma
import scipy.interpolate as si
import scipy.linalg as la
from   scipy.interpolate import Akima1DInterpolator as scipy_akima

#

from mittag_leffler import ml 
from inc_post   import * 
from inc_sub2   import * 

################################################################################
#
# Select the numerical test by setting fprefix as follows:
#   Example 4.1 (a)     fprefix = "ex01a"  with psi = t 
#   Example 4.1 (b)     fprefix = "ex01b"  with psi = e^(5t)
#   Example 4.2         fprefix = "ex02a"
#   Example 4.3 (a)     fprefix = "ex03a"  without noise to input 
#*  Example 4.3 (b)     fprefix = "ex03b"  with noise to input (SEPARATED)
#

#fprefix = "ex01a"
#fprefix = "ex01b"
#fprefix = "ex02a"
fprefix = "ex03a"
#fprefix = "ex03b"      # SEPARATED. Use main_ex3b.py


#
# Define the number of N_l = 10 * 2^(l) for l=0,lmax-1: 
#   lmax = 5    To check if the code works properly.
#   lmax = 11   To run fully for the manuscript. But it takes time! 

lmax = 11

#

#
# To define some important paths for output:

OUTDIR, FIGDIR = make_paths()


# Define the largest number of the fixed-point interation steps:
M = 200 


################################################################################

if   fprefix == "ex01a":
    from inc_ex1a  import * 
    HaveExact = True  
#    
elif fprefix == "ex01b":
    from inc_ex1b  import * 
    HaveExact = True  
#    
elif fprefix == "ex02a":
    from inc_ex2   import * 
    HaveExact = True  
#
elif fprefix == "ex03a":
    from inc_ex3a  import * 
    HaveExact = False 
    furef = "ex03a"
    luref = lmax-1
    Nuref = int(10 * int(2)**luref)
#
#elif fprefix == "ex03b":
#   from inc_ex3b  import * 
#   HaveExact = False 
#   furef = "ex03a"
#   luref = lmax-1
#   Nuref = int(10 * int(2)**luref)
#

print ( '\nfprefix: %s\n' %( fprefix ) )

N_int = 500 + 2 

T_int = np.linspace( c_a, c_b, N_int )

t_int = T_int[1:N_int-1].copy()

print ( '\nBefore the setting (default values):' )
parameter_show () 

parameter_sync ( \
    c_epsilon, c_a, c_b, c_p, c_q, c_alpha, c_beta, \
    c_lambda, c_eta1, c_eta2, c_sigma )

print ( '\nAfter the setting:' )
parameter_show () 

coef_phip = -f_phip(c_eta2) / ( f_phip(c_eta1) + f_phip(c_eta2) )


if ( HaveExact ) :
#   to save errors |U_N - Uexa| in the max norm:
    o_Err0 = np.zeros( lmax )
#   to save errors |U_N - Uexa| in the l2- norm:
    o_Err2 = np.zeros( lmax )
#


# to save errors in the max norm:
a_Err0 = np.zeros( lmax )

# to save errors in the l2- norm:
a_Err2 = np.zeros( lmax )

# to save array of the mesh sizes N:
a_N    = np.zeros( lmax )

for l in range(0,lmax):
    a_N[l] = int( 10 * int(2)**l )

print ( a_N )


# contraction rate estimate: 
a_cont = np.zeros(M)

#
# Convergence rate estimating:
#
for l in range ( 0, lmax ):

    N = int(a_N[l])
#
#   tau = zeros ( N + 1 )
#   s   = zeros ( N + 1 )
#   t   = zeros ( N     )
#   xi  = zeros ( N     )
#
    dt, t, xi = gene_meshes ( N, c_a, c_b )

    coef = dt / gamma( c_beta )

    print ( "*** l = %3d,  N_l = %7d "%( l, N ) )
    print ( "a, b, sigma, dt =", c_a, c_b, c_sigma, dt )

    print ( "Procedure P: (1)")

    print ( "Calculating G ...")

    a_G = calc_G ( N, t, f_psi )

#   print ( "G =\n", a_G )

    print ( "Calculating E ...")

#   a_E = np.zeros( [N,N] )
#   for k in range(1,N):
#       for j in range(1,k+1):
#           a_E[k,j] = f_Eps ( c_alpha, t[k], xi[j], c_lambda, f_psi, f_dpsi )

    a_E = calc_E ( N, t, xi, f_psi, f_dpsi )

#   print ( "E =\n" )
#   for k in range(1,N):
#       print ( "k:", k, a_E[k,1:k+1] )

    if ( HaveExact ) :  
        Uexa = np.zeros(N)
        Uexa = f_uexact ( t )
#

    print ( "Procedure P: (2)")

    print ( "+ Initializing guess ...")

    Uold = np.random.rand( N )*20 - 10

    Unew = np.zeros(N)

    Utau = np.zeros(N+1)

    a_V = np.zeros ( N )
    a_R = np.zeros ( N )
    a_W = np.zeros ( N )
    a_Y = np.zeros ( N )
    a_H = np.zeros( [N,N] )

#
#   contraction rate estimate: 
    M1 = 0 

    for m in range(1,M+1):

        print ( "\n*   Iteration m =", m )

        print ( "Procedure P: (2) (a)")

        print ( "+ U^{m-1} (xi_j): Interpolating ...")
#
#<< Eq. (52):
#
#       f_utau = si.interp1d ( t, Uold, kind='cubic' )
#       f_utau = scipy_akima ( t, Uold )
#       Utau = f_utau ( xi )
#
#       Linear interpolation is fine enough:

        Utau = np.zeros(N)
        for k in range( 1, N ):
            Utau[k] = ( Uold[k-1] + Uold[k] )/2
#
        print ( "+ I_k: Calculating ...")

        a_V[0] = 0 
        for k in range(1,N):
            Ik = 0
            for i in range(1,k+1):
                Ik = Ik + f_F ( c_beta, t[k], xi[i], Utau[i], f_psi, f_dpsi, f_f )
            a_V[k] = Ik * coef 

        print ( "Procedure P: (2) (b)")

        print ( "+ R_i: Interpolating ...")

#       f_Vt = si.interp1d ( t, a_V, kind='cubic' )
#       f_Vt = scipy_akima ( t, a_V )
#       a_R = f_Vt ( xi )
#
#       Linear interpolation is fine enough:

        a_R = np.zeros ( N )
        for k in range(1,N):
            a_R[k] = ( a_V[k-1] + a_V[k] )/2
#>> 

#<< Eq. (51):

        print ( "+ W_i: Calculating ...")

        for j in range (1,N):
            a_W[j] = a_R[j] + coef_phip*a_V[N-1]
#>> 

        print ( "Procedure P: (2) (c)")

#<< Eq. (50):

        print ( "+ H_{i,k}^{m-1}: Calculating ...")

        for j in range (1,N):       
            a_Y[j] = f_phiq ( a_W[j] )
        for k in range(1,N):
            for j in range(1,k+1):
                a_H[k,j] = a_E[k,j] * a_Y[j]
#>>

        print ( "Procedure P: (2) (d)")

#<< Eq (49): 

        print ( "+ U_{0}^{m}: Calculating ...")

        stmp = 0 
        for k in range(1,N):
            stmp = stmp + a_H[N-1,k]
        u0g0 = stmp * dt 
        Unew[0] = a_G[0] * u0g0
#
        print ( "+ U_{k}^{m}: Calculating ...")

        for k in range(1,N):
            stmp = 0 
            for j in range(1,k+1):
                stmp = stmp + a_H[k,j]
            Unew[k] = u0g0*a_G[k] + stmp*dt 
#>>

        print ( "Procedure P: (2) (e)")

        aerr = max( abs( Unew - Uold )  ) 
        unom = max( abs( Unew ) ) 

        rerr = aerr / max( unom, 1.0e-20 ) 

        print ( "+ Contraction err = | U^{m} - U^{m-1} |:"  )
        print ( "\t Ab.err = %16.9e \t Re.err = %16.9e" %( aerr, rerr )  )

#       a_cont[m] = aerr 
        a_cont[m] = max( aerr, MachEps ) 

        if ( aerr <= c_epsilon * unom ):
            M1 = m
            print ( "\nCONVERGED! Exit Procedure P." )
            break 

        Uold = Unew.copy()
#

#   Checking some aspects:

    check_bc ( N, Unew )

    ftmp = os.path.join( FIGDIR, fprefix + str('-ContR') + intshift(l,3) )
    show_ContRate ( M1, a_cont, fname=ftmp , report=False  )

    funame = fprefix + str('-Umild') + intshift(l,3) 

    ftmp = os.path.join( FIGDIR, funame )
    print ( "Plotting Unew ... \n\t" + ftmp  )
    show_MildSol ( t, Unew, fname=ftmp, report=False )

    ftmp = os.path.join( OUTDIR, funame )
    print ( "Writing Unew to files ... \n\t" + ftmp + " .cvs, .txt" )
    w_MildSol ( t, Unew, os.path.join( OUTDIR, ftmp ), ext='csv' )
    w_MildSol ( t, Unew, os.path.join( OUTDIR, ftmp ), ext='txt' )
#

    if ( HaveExact ) :  
        funame = os.path.join( FIGDIR,fprefix + str('-CompU') + intshift(l,3) )
        plot_CompareSol ( t, Unew, Uexa, fname=funame, report=False )
#

    if ( l == 0 ):

        t_sav = np.zeros(N)
        t_sav = t.copy()
        U_sav = np.zeros(N)
        U_sav = Unew.copy()
        a_Err0[l] = max( abs( Unew ) )
        a_Err2[l] = np.sqrt( sum(abs(Unew)**2) / len(Unew) )

    else:
#
#       Computing Err[l] for l = 1,2,...
#  

        e0, e2 = f_errint ( t_int, t_sav, U_sav, t, Unew )

        a_Err2[l] = e2 

        a_Err0[l] = max( e0, MachEps )

        print ( "\n"  + "-" * 70  )
        print ( "l = %2d, N_l = %6d, Err = | U_{N_{l}} - U_{N_{l-1}} |:"%( l, N ) )
        print ( "\tMaxErr = %16.9e, L2-Err = %16.9e"%( e0, e2 ) )
        print ( "-" * 70  )
#        
        t_sav = None 
        t_sav = np.zeros(N)
        t_sav = t.copy()

        U_sav = None 
        U_sav = np.zeros(N)
        U_sav = Unew.copy()
#

#** Compare to the exact solution:

    if ( HaveExact ) :  

        e0, e2 = f_errest ( Unew, Uexa )

        print ( \
            "\nl = %3d, Max|Uexact-U_N| = %16.9e, L2 |Uexact-U_N| = %16.9e"%( \
             l, e0, e2 ) )
        print ( "*" * 80 + "\n" )

        o_Err0[l] = e0
        o_Err2[l] = e2 

#
#
#
    t   = None 
    tau = None 
    s   = None 
    xi  = None
    Uold = None
    Unew = None 
    Utau = None 
    a_V = None 
    a_R = None 
    a_W = None 
    a_Y = None
    a_H = None
    a_G = None
    a_E = None 
    Uexa = None 
#
#

#*********************************************************************
#  |Unew - Uold|:  i 

ftmp = os.path.join( FIGDIR, fprefix + str('-iConvR0') )
show_ConvRate ( lmax, a_N, a_Err0, fname=ftmp, report=False, alog=True ) 

ftmp = os.path.join( FIGDIR, fprefix + str('-iConvR0') )
show_ConvRate ( lmax, a_N, a_Err0, fname=ftmp, report=False, alog=False ) 


print ( '\nThe error ||U_{N2} - U_{N1}|| in the max norm:' )
w_ConvRate ( lmax, a_N, a_Err0, fname=None )

ftmp = os.path.join( OUTDIR, fprefix + str('-iConvR0') )
w_ConvRate ( lmax, a_N, a_Err0, fname=ftmp )

#

ftmp = os.path.join( FIGDIR, fprefix + str('-iConvR2') )
show_ConvRate ( lmax, a_N, a_Err2, fname=ftmp, report=False, alog=True ) 

ftmp = os.path.join( FIGDIR, fprefix + str('-iConvR2') )
show_ConvRate ( lmax, a_N, a_Err2, fname=ftmp, report=False, alog=False ) 

print ( '\nThe error ||U_{N2} - U_{N1}|| in the l2- norm:' )
w_ConvRate ( lmax, a_N, a_Err2, fname=None )

ftmp = os.path.join( OUTDIR, fprefix + str('-iConvR2') )
w_ConvRate ( lmax, a_N, a_Err2, fname=ftmp )


#*********************************************************************
#  |Ucal - Uexa|: e
#
if ( HaveExact ) :     

    ftmp = os.path.join( FIGDIR, fprefix + str('-eConvR0') )
    plot_ConvRate ( lmax, a_N, o_Err0, fname=ftmp, report=False, alog=True  ) 
    plot_ConvRate ( lmax, a_N, o_Err0, fname=ftmp, report=False, alog=False ) 

    print ( '\nThe error ||U_{calc} - U_{exact}|| in the max norm:' )
    w_ConvRate ( lmax, a_N, o_Err0, fname=None )

    ftmp = os.path.join( OUTDIR, fprefix + str('-eConvR0') )
    w_ConvRate ( lmax, a_N, o_Err0, fname=ftmp )
#
    ftmp = os.path.join( FIGDIR, fprefix + str('-eConvR2') )
    plot_ConvRate ( lmax, a_N, o_Err2, fname=ftmp, report=False, alog=True  ) 
    plot_ConvRate ( lmax, a_N, o_Err2, fname=ftmp, report=False, alog=False ) 

    print ( '\nThe error ||U_{calc} - U_{exact}|| in the l2- norm:' )
    w_ConvRate ( lmax, a_N, o_Err2, fname=None )

    ftmp = os.path.join( OUTDIR, fprefix + str('-eConvR2') )
    w_ConvRate ( lmax, a_N, o_Err2, fname=ftmp )

#
#
