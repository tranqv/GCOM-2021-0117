#!/usr/bin/env python
# coding: utf-8

#
#   Numerical Examples (Example 3b & postprocessing)
#
#   Revised manuscript: 
#       "On mild solutions of the p-Laplacian fractional Langevin equations with
#       anti-periodic type boundary conditions"
#

import os
import numpy as np 
from   scipy.special import gamma
import scipy.interpolate as si
#import scipy.linalg as la
from   scipy.interpolate import Akima1DInterpolator as scipy_akima

#

from mittag_leffler import ml 
from inc_post   import * 
from inc_sub2   import * 

################################################################################
#
# Select the numerical test by setting fprefix as follows:
#   Example 4.3 (b)     fprefix = "ex03b"  with noise to input 
#

fprefix = "ex03b"

#
# For N_l = 10 * 2^(l) for l=0,lmax-1 
#   luref = 3    To check if the code works properly. luref <= lmax-1
#   luref = 7    To run for the manuscript. But it takes time! 

luref = 3


#
# To define some important paths for output:

OUTDIR, FIGDIR = make_paths()


# Define the largest number of the fixed-point interation steps:
M = 200 


# the max number of the noise levels:
nmax = 10 

################################################################################

from inc_ex3b  import * 

HaveExact = False 
furef = "ex03a"
Nuref = int( 10 * int(2)**luref )

#


N_int = 500 + 2 

T_int = np.linspace( c_a, c_b, N_int )

t_int = T_int[1:N_int-1].copy()

o_Err0 = np.zeros( nmax )
o_Err2 = np.zeros( nmax )

#

# Noise amplitudes:
a_N = np.zeros( nmax )
for n in range(0,nmax):
    a_N[n] = 2.0**( -n ) / 20       # The amplitudes in the manuscript.
#   a_N[n] = 2.0**( -n ) / 10       # Larger apmplitudes. 
#




# to save errors in the max norm:
a_Err0 = np.zeros( nmax )

# to save errors in the l2- norm:
a_Err2 = np.zeros( nmax )


N = int( 10 * int(2)**luref )

dt, t, xi = gene_meshes ( N, c_a, c_b )

print ( "N, a, b, sigma, dt =", N, c_a, c_b, c_sigma, dt )

# to back up the orginal parameters:
s_alpha  = c_alpha  
s_beta   = c_beta   
s_lambda = c_lambda 
s_eta1   = c_eta1   
s_eta2   = c_eta2   


if ( HaveExact ) :  
    Uref = f_uexact ( t )
else:
#   Pick up the reference:
    funame = os.path.join( OUTDIR, furef + str('-Umild') + intshift(luref,3) )
    print ( 'Loading U_N as solution without noise ...\n\t' + funame ) 
    tE, UE = load_U ( luref, funame )
    Uref = UE.copy() 
    tE = None 
    uE = None 
#




Utmp = np.zeros(N) # pointer for the swapping procedure 
Uold = np.zeros(N)
Unew = np.zeros(N)

t_sav = np.zeros(N)
U_sav = np.zeros(N)

Utau = np.zeros(N+1)

a_V = np.zeros ( N )
a_R = np.zeros ( N )
a_W = np.zeros ( N )
a_Y = np.zeros ( N )
a_H = np.zeros( [N,N] )

#   contraction rate estimate: 
a_cont = np.zeros(M)


#
# Convergence rate estimating w.r.t the noise amplitudes:
#
for n in range ( 0, nmax ):

    noise_amplitude = a_N[n]

    print ( "*** n =", n )
    print ( '\n*** Noise amplitude = %14.7E\n' %( noise_amplitude ) )

    print ( '\nBefore the setting (default values):' )

    parameter_show () 

#
#   generating random noise to the parameters:
#
    n_s = 6
    x_rand   = np.random.rand( n_s ) - 0.5*np.ones( n_s ) 
    x_sign   = np.sign( x_rand ) 
    d_noise  = noise_amplitude * x_sign[0]
    c_alpha  = s_alpha  * ( 1 + noise_amplitude * x_sign[1] ) 
    c_beta   = s_beta   * ( 1 + noise_amplitude * x_sign[2] )
    c_lambda = s_lambda * ( 1 + noise_amplitude * x_sign[3] ) 
    c_eta1   = s_eta1   * ( 1 + noise_amplitude * x_sign[4] )  
    c_eta2   = s_eta2   * ( 1 + noise_amplitude * x_sign[5] ) 

    parameter_sync ( \
        c_epsilon, c_a, c_b, c_p, c_q, \
        c_alpha, c_beta, c_lambda, c_eta1, c_eta2, \
        c_sigma )

    print ( '\nAfter the setting:' )
    parameter_show () 

    print ( "d_noise = ", d_noise )

    coef_phip = -f_phip(c_eta2) / ( f_phip(c_eta1) + f_phip(c_eta2) )

    coef = dt / gamma( c_beta )

    print ( "Procedure P: (1)")

    print ( "+ Calculating G ..." )

    a_G = calc_G ( N, t, f_psi )

#   print ( "G =\n", a_G )

    print ( "+ Calculating E ...")

    a_E = calc_E ( N, t, xi, f_psi, f_dpsi )

    print ( "Procedure P: (2)")

    print ( "+ Initializing guess ...")

    Uini = np.random.rand( N )*20 - 10

    Uold = Uini.copy()

#
#   contraction rate estimate: 
    M1 = 0 

    for m in range(0,M):

        print ( "\n*   Iteration m =", m+1 )

        print ( "Procedure P: (2) (a)")

        print ( "+ U^{m-1} (xi_j): Interpolating ...")

#       f_utau = si.interp1d ( t, Uold, kind='cubic' )
#       f_utau = scipy_akima ( t, Uold )
#       Utau = f_utau ( xi )
#
#       Linear interpolation is fine enough:

        Utau[0] = 0
        for k in range( 1, N ):
            Utau[k] = ( Uold[k-1] + Uold[k] )/2

        print ( "+ I_k: Calculating ...")

        a_V[0] = 0 
        for k in range(1,N):
            Ik = 0
            for i in range(1,k+1):
                Ik = Ik + f_F ( c_beta, t[k], xi[i], Utau[i], f_psi, f_dpsi, f_f )
            a_V[k] = Ik * coef 

        print ( "Procedure P: (2) (b)")

        print ( "+ R_i: Interpolating ...")

#
#       f_Vt = si.interp1d ( t, a_V, kind='cubic' )
#       f_Vt = scipy_akima ( t, a_V )
#       a_R = f_Vt ( xi )
#
#       Linear interpolation is fine enough:
        a_R[0] = 0
        for k in range(1,N):
            a_R[k] = ( a_V[k-1] + a_V[k] )/2

        print ( "+ W_i: Calculating ...")

        for j in range (1,N):
            a_W[j] = a_R[j] + coef_phip*a_V[N-1]

        print ( "Procedure P: (2) (c)")

        print ( "+ H_{i,k}^{m-1}: Calculating ...")

        for j in range (1,N):       
            a_Y[j] = f_phiq ( a_W[j] )
        for k in range(1,N):
            for j in range(1,k+1):
                a_H[k,j] = a_E[k,j] * a_Y[j]
#
        print ( "Procedure P: (2) (d)")

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
#
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

    ftmp = os.path.join( FIGDIR, fprefix + str('-CTR-n') + intshift(n,3) )
    show_ContRate ( M1, a_cont, fname=ftmp , report=False )

    funame = fprefix + str('-U---n') + intshift(n,3)  

    ftmp = os.path.join( FIGDIR, funame )
    print ( "Plotting U_N^delta ... \n\t" + ftmp )
    show_MildSol ( t, Unew, fname=ftmp, report=False  )

    ftmp = os.path.join( OUTDIR, funame )
    print ( "Writing U_N^delta to files ... \n\t" + ftmp + " .cvs, .txt" )
    w_MildSol ( t, Unew, ftmp, ext='csv' )
    w_MildSol ( t, Unew, ftmp, ext='txt' )
#
#
    funame = os.path.join( FIGDIR, fprefix + str('-UeU-n') + intshift(n,3) )
    plot_Sols_3 ( t, Unew, t, Uref, fname=funame, report=False )
#

    if ( n == 0 ):

        t_sav = t.copy()
        U_sav = Unew.copy()
        a_Err0[n] = max( abs( Unew ) )
        a_Err2[n] = np.sqrt( sum(abs(Unew)**2) / len(Unew) )

    else:
#
#       Computing Err[n] for l = 1,2,...
#  
        e0, e2 = f_errint ( t_int, t_sav, U_sav, t, Unew )
        a_Err2[n] = e2 
        a_Err0[n] = max( e0, MachEps )

        print ( "\n"  + "-" * 70  )
        print ( "n = %2d, Err = | U_N^{delta_{n}} - U_N^{delta_{n-1}} |:"%( n ) )
        print ( "\tMaxErr = %16.9e, L2-Err = %16.9e"%( e0, e2 ) )
        print ( "-" * 70  )
#        
        t_sav = t.copy()
        U_sav = Unew.copy()
#
#
#** Compare to the exact solution:

#   if ( HaveExact ) :  

    e0, e2 = f_errest ( Unew, Uref )

    print ( "\n"  + "*" * 80  )
    print ( \
        "n = %3d, Max|U_N^delta - U_N| = %16.9e, L2 |U_N^delta - U_N| = %16.9e"% ( \
        n, e0, e2) )
    print ( "*" * 80 + "\n" )

    o_Err0[n] = e0
    o_Err2[n] = e2 
#
    a_G = None
    a_E = None 
    Uini = None 
#
#

#*********************************************************************

ftmp = os.path.join( FIGDIR, fprefix + str('-iCVR0') )
show_ConvRate ( nmax, a_N, a_Err0, fname=ftmp, report=False, alog=True  ) 
show_ConvRate ( nmax, a_N, a_Err0, fname=ftmp, report=False, alog=False ) 

print ( '\nThe error | U_N^{delta_{n}} - U_N^{delta_{n-1}} | in the max norm:' )
w_ConvRate ( nmax, a_N, a_Err0, fname=None )
ftmp = os.path.join( OUTDIR, fprefix + str('-iCVR0') )
w_ConvRate ( nmax, a_N, a_Err0, fname=ftmp )

#

ftmp = os.path.join( FIGDIR, fprefix + str('-iCVR2') )
show_ConvRate ( nmax, a_N, a_Err2, fname=ftmp, report=False, alog=True  ) 
show_ConvRate ( nmax, a_N, a_Err2, fname=ftmp, report=False, alog=False ) 

print ( '\nThe error | U_N^{delta_{n}} - U_N^{delta_{n-1}} | in the l2- norm:' )
w_ConvRate ( nmax, a_N, a_Err2, fname=None )
ftmp = os.path.join( OUTDIR, fprefix + str('-iCVR2') )
w_ConvRate ( nmax, a_N, a_Err2, fname=fprefix + str('-iCVR2') )


#*********************************************************************

ftmp = os.path.join( FIGDIR, fprefix + str('-eCVR0') )
plot_CVR0 ( nmax, a_N, o_Err0, fname=ftmp, report=False, alog=True  ) 
plot_CVR0 ( nmax, a_N, o_Err0, fname=ftmp, report=False, alog=False ) 

print ( '\nThe error | U_N^delta - U_N | in the max norm:' )
w_ConvRate ( nmax, a_N, o_Err0, fname=None )
ftmp = os.path.join( OUTDIR, fprefix + str('-eCVR0') )
w_ConvRate ( nmax, a_N, o_Err0, fname=ftmp )

#

ftmp = os.path.join( FIGDIR, fprefix + str('-eCVR2') )
plot_CVR0 ( nmax, a_N, o_Err2, fname=ftmp, report=False, alog=True  ) 
plot_CVR0 ( nmax, a_N, o_Err2, fname=ftmp, report=False, alog=False ) 

print ( '\nThe error | U_N^delta - U_N | in the l2- norm:' )
w_ConvRate ( nmax, a_N, o_Err2, fname=None )
ftmp = os.path.join( OUTDIR, fprefix + str('-eCVR2') )
w_ConvRate ( nmax, a_N, o_Err2, fname=ftmp )

#
#
