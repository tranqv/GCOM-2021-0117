#!/usr/bin/env python
# coding: utf-8

import numpy as np 
from   scipy.special import gamma
import matplotlib.pyplot as plt
from   mittag_leffler import ml 

MachEps = np.finfo(float).eps

def cm2in (x):
    return x/2.54

#

def show_UV ( t, Unew, a_V, fname=None, report=False ):

    plt.rcParams['figure.figsize'] = [ cm2in(10) , cm2in(10) ]

    fig, axs = plt.subplots(1,2)

    axs[0].plot( t, Unew,
        color='red', linestyle='dashed', linewidth=2.0,
        label=r'\mathbf{U}' )
    axs[0].grid(True)
    axs[0].set_xlabel( r'$t$' )
    axs[0].set_ylabel( r'$U$' )
#   axs[0].set_title( "$U$ v.s. $t$" )

    axs[0].legend()

    axs[1].plot( t, a_V, 
#                    color='blue', linestyle='dashed', linewidth=2.0, 
                     color='blue', linewidth=2.0, 
                     label=r'I_{a+}^{\beta,\psi} f (t,\psi,\mathbf{U})' )
    axs[1].grid(True)
    axs[1].set_xlabel( r'$t$' )
    axs[1].set_ylabel( r'$I$' )
#   axs[1].set_title( "$I$ v.s. $t$" )
    axs[1].legend()

    if ( report ):
        plt.show()
    else:
        if ( fname is not None ):
            fig.savefig ( fname + '.jpg' , bbox_inches = 'tight' )
            fig.savefig ( fname + '.pdf' , bbox_inches = 'tight' )

    plt.close('all')
    return None 
#


def show_MildSol ( t, Unew, fname=None, report=True ):

    plt.rcParams['figure.figsize'] = [ cm2in(10) , cm2in(10) ]

    fig, axs = plt.subplots()

    axs.plot( t, Unew,
#       color='red', linestyle='dashed', linewidth=2.0,
#       color='black', linestyle='dashed', linewidth=2.0,
        color='black', linewidth=2.0,
        label=r'$ \mathbf{U} $' )
    axs.grid(True)
    axs.set_xlabel( r'$t$' )
#   axs.set_ylabel( r'$u$' )
#   axs.set_title( "Mild solution" )

    axs.legend()

    if ( report ): 
        plt.show()
    else:
        if ( fname is not None ):
            fig.savefig ( fname + '.jpg' , bbox_inches = 'tight' )
            fig.savefig ( fname + '.pdf' , bbox_inches = 'tight' )

    plt.close('all')
    return None 
#


#
#  Plotting  
#  U_calc., U_exac
#
def plot_CompareSol ( t, Unew, Uexa, fname=None, report=True ):

    plt.rcParams['figure.figsize'] = [ cm2in(10) , cm2in(10) ]

    fig, axs = plt.subplots()

    axs.plot( t, Unew,
#       color='red', linestyle='dashed', linewidth=2.0,
#       color='black', linestyle='dashed', linewidth=2.0,
        color='black', linestyle='dashed', linewidth=1.5,
        label=r'$ \mathbf{U}_{ \rm calc. } $' )

    axs.plot( t, Uexa,
#       color='red', linestyle='dashed', linewidth=2.0,
#       color='black', linestyle='dashed', linewidth=2.0,
        color='black', linewidth=1.5,
        label=r'$ \mathbf{U}_{ \rm exact } $' )


    axs.grid(True)
    axs.set_xlabel( r'$t$' )
#   axs.set_ylabel( r'$u$' )
#   axs.set_title( "Mild solution" )

    axs.legend()

    if ( report ): 
        plt.show()
    else:
        if ( fname is not None ):
            fig.savefig ( fname + '.jpg' , bbox_inches = 'tight' )
            fig.savefig ( fname + '.pdf' , bbox_inches = 'tight' )

    plt.close('all')
    return None 










def show_ContRate ( M1, a_cont, fname=None, report=True ):

    plt.rcParams['figure.figsize'] = [ cm2in(10) , cm2in(10) ]

    vm = np.zeros ( M1 )
    ve = np.zeros ( M1 )
    for i in range(M1):
        vm[i] = i+1
        ve[i] = np.log( max( a_cont[i], MachEps ) )

    fig, axs = plt.subplots()

    axs.plot( vm, ve,
#       color='red', linestyle='dashed', linewidth=0.5,
        color='black', linestyle='dashed', linewidth=1.0,
        marker='s', markerfacecolor='darkgreen', markersize=8,
        label=r'$\log \left\Vert \mathbf{U}^{m}-\mathbf{U}^{m-1} \right\Vert$' )
    axs.grid(True)
    axs.set_xlabel( r'$m$' )
#   axs.set_title( r'$\log \left\| \mathbf{U}^{m} - \mathbf{U}^{m-1} \right\|$' ) 
    axs.legend()

    if ( report ): 
        plt.show()
    else:
        if ( fname is not None ):
            fig.savefig ( fname + '.jpg' , bbox_inches = 'tight' )
            fig.savefig ( fname + '.pdf' , bbox_inches = 'tight' )

    plt.close('all')

    return None 
#



#
# | UN1 - UN2 |:
#
def show_ConvRate ( lmax, a_N, a_Err, fname=None, report=True, alog=True ):

    plt.rcParams['figure.figsize'] = [ cm2in(10) , cm2in(10) ]

    if ( alog ): 

        fig, axs = plt.subplots()

        tmp = np.log ( np.fmax( a_Err, MachEps ) ) 

        axs.plot( np.log ( a_N[1:lmax] ), tmp[1:lmax],
#           color='red', linestyle='dashed', linewidth=0.5,
            color='black', linestyle='dashed', linewidth=1.0,
            marker='o', markerfacecolor='blue', markersize=8,
            label=r'$\log \left\Vert\mathbf{U}_{N_l}-\mathbf{U}_{N_{l-1}} \right\Vert$' )
        axs.grid(True)
        axs.set_xlabel( r'$\log N_l$' )
#       axs.set_ylabel( "$\log ||U_{N_l}-U_{N_{l-1}}||_{F}$ " )
#       axs.set_title( "$\log ||U_{N_l}-U_{N_{l-1}}||$ v.s. $\log N_l$" )
        axs.legend()

        if ( report ): 
            plt.show()
        else:
            if ( fname is not None ):
                fig.savefig ( fname + '_log.jpg' , bbox_inches = 'tight' )
                fig.savefig ( fname + '_log.pdf' , bbox_inches = 'tight' )

    else:

        fig, axs = plt.subplots()

        axs.plot( a_N[1:lmax], a_Err[1:lmax] ,
#           color='blue', linestyle='dashed', linewidth=0.5,
            color='black', linestyle='dashed', linewidth=1.0,
            marker='o', markerfacecolor='red', markersize=8,
            label=r'$\left\Vert\mathbf{U}_{N_l}-\mathbf{U}_{N_{l-1}}\right\Vert$' )
        axs.grid(True)
        axs.set_xlabel( r'$N_l$' )
#       axs.set_ylabel( "$||U_{N_l}-U_{N_{l-1}}||$ " )
#       axs.set_title( "$||U_{N_l}-U_{N_{l-1}}||$ v.s. $N_l$" )
        axs.legend()

        if ( report ): 
            plt.show()
        else:
            if ( fname is not None ):
                fig.savefig ( fname + '.jpg' , bbox_inches = 'tight' )
                fig.savefig ( fname + '.pdf' , bbox_inches = 'tight' )
#
    plt.close('all')
    return None 
#



#
# | Unew - Uexa |:
#
def plot_CVR0 ( lmax, a_N, a_Err, fname=None, report=True, alog=True ):

    plt.rcParams['figure.figsize'] = [ cm2in(10) , cm2in(10) ]

    if ( alog ): 

        fig, axs = plt.subplots()

        tmp = np.log ( np.fmax( a_Err, MachEps ) ) 

        axs.plot( np.log ( a_N ), tmp ,
#           color='red', linestyle='dashed', linewidth=0.5,
            color='black', linestyle='dashed', linewidth=1.0,
            marker='o', markerfacecolor='red', markersize=8,
            label=r'$\log \left\Vert\mathbf{U}_{N}^{\delta}-\mathbf{U}_{N} \right\Vert$' )
        axs.grid(True)
        axs.set_xlabel( r'$\log \delta $' )
#       axs.set_ylabel( "$\log ||U_{N_l}-U_{N_{l-1}}||_{F}$ " )
#       axs.set_title( "$\log ||U_{N_l}-U_{N_{l-1}}||$ v.s. $\log N_l$" )
        axs.legend()

        if ( report ): 
            plt.show()

        if ( fname is not None ):
            fig.savefig ( fname + '_log.jpg' , bbox_inches = 'tight' )
            fig.savefig ( fname + '_log.pdf' , bbox_inches = 'tight' )

    else:

        fig, axs = plt.subplots()

        axs.plot( a_N, a_Err ,
#           color='blue', linestyle='dashed', linewidth=0.5,
            color='black', linestyle='dashed', linewidth=1.0,
            marker='o', markerfacecolor='red', markersize=8,
            label=r'$\left\Vert\mathbf{U}_{N}^{\delta}-\mathbf{U}_{N} \right\Vert$' )
        axs.grid(True)
        axs.set_xlabel( r'$\delta$' )
#       axs.set_ylabel( "$||U_{N_l}-U_{N_{l-1}}||$ " )
#       axs.set_title( "$||U_{N_l}-U_{N_{l-1}}||$ v.s. $N_l$" )
        axs.legend()

        if ( report ): 
            plt.show()
        else:
            if ( fname is not None ):
                fig.savefig ( fname + '.jpg' , bbox_inches = 'tight' )
                fig.savefig ( fname + '.pdf' , bbox_inches = 'tight' )
#
    plt.close('all')
    return None 
#



#
# | Unew - Uexa |:
#
def plot_ConvRate ( lmax, a_N, a_Err, fname=None, report=True, alog=True ):

    plt.rcParams['figure.figsize'] = [ cm2in(10) , cm2in(10) ]

    if ( alog ): 

        fig, axs = plt.subplots()

        tmp = np.log ( np.fmax( a_Err, MachEps ) ) 

        axs.plot( np.log ( a_N[1:lmax] ), tmp[1:lmax] ,
#           color='red', linestyle='dashed', linewidth=0.5,
            color='black', linestyle='dashed', linewidth=1.0,
            marker='o', markerfacecolor='blue', markersize=8,
            label=r'$\log \left\Vert\mathbf{U}_{\rm calc.}-\mathbf{U}_{\rm exact} \right\Vert$' )
        axs.grid(True)
        axs.set_xlabel( r'$\log N_l$' )
#       axs.set_ylabel( "$\log ||U_{N_l}-U_{N_{l-1}}||_{F}$ " )
#       axs.set_title( "$\log ||U_{N_l}-U_{N_{l-1}}||$ v.s. $\log N_l$" )
        axs.legend()

        if ( report ): 
            plt.show()

        if ( fname is not None ):
            fig.savefig ( fname + '_log.jpg' , bbox_inches = 'tight' )
            fig.savefig ( fname + '_log.pdf' , bbox_inches = 'tight' )

    else:

        fig, axs = plt.subplots()

        axs.plot( a_N[1:lmax], a_Err[1:lmax] ,
#           color='blue', linestyle='dashed', linewidth=0.5,
            color='black', linestyle='dashed', linewidth=1.0,
            marker='o', markerfacecolor='red', markersize=8,
            label=r'$\left\Vert\mathbf{U}_{\rm calc.}-\mathbf{U}_{\rm exact} \right\Vert$' )
        axs.grid(True)
        axs.set_xlabel( r'$N_l$' )
#       axs.set_ylabel( "$||U_{N_l}-U_{N_{l-1}}||$ " )
#       axs.set_title( "$||U_{N_l}-U_{N_{l-1}}||$ v.s. $N_l$" )
        axs.legend()

        if ( report ): 
            plt.show()
        else:
            if ( fname is not None ):
                fig.savefig ( fname + '.jpg' , bbox_inches = 'tight' )
                fig.savefig ( fname + '.pdf' , bbox_inches = 'tight' )
#
    plt.close('all')
    return None 
#





#
def w_ConvRate ( lmax, a_N, a_Err, fname=None ):

    from numpy   import zeros 

    r = zeros(lmax)

    for l in range (1,lmax):
        r[l] = np.log( a_Err[l]/a_Err[l-1] ) / np.log( a_N[l-1]/a_N[l] )
    r[0] = 0

    if ( fname is not None ) :

        f = open( fname + ".csv", "w" )
        for l in range (0,lmax):
            txt = '%d,%.9E,%.9E,%.9E,%.9E,%.9E\n' % ( l, a_N[l], a_Err[l], 
                    np.log(a_N[l]),  np.log(a_Err[l]), r[l]  )
            f.write ( txt )
        f.close()

    else:

        for l in range (0,lmax):
            txt = '%d, %.9E, %.9E, %.9E, %.9E, %.9E' % ( l, a_N[l], a_Err[l], 
                    np.log(a_N[l]),  np.log(a_Err[l]), r[l]  )

            print ( txt ) 

    return None 
#
#



def w_MildSol ( t, U, fname, ext='csv' ):

    if ( ext == 'csv' ) :
        f = open( fname + ".csv", "w" )
        N = len(t)
        for k in range (0,N):
            f.write ( '%d,%.15E,%.15E\n' % ( k, t[k], U[k] )  )
        f.close()

    if ( ext == 'txt' ) :
        f = open( fname + ".txt", "w" )
        N = len(t)
        for k in range (0,N):
            f.write ( '%10d  %23.15E  %23.15E\n' % ( k, t[k], U[k] )  )
        f.close()

    return None 



def r_MildSol ( fname, ext='csv', report=False ):
#
#   Relying on 
#       st.rstrip()
#       st.split()
#       st.split(',')
#
#
    if ( ext == 'csv' ) :

        f = open( fname + ".csv", "r" )

        contents = f.readlines() 

        N = len( contents )

        if ( report ):
            print ( "N = ", N  )

        t = np.zeros (N)
        U = np.zeros (N)

        k = 0

        for line in contents :

            line = line.rstrip()
            cols = line.split(',')

            t[k] = float( cols[1] ) 
            U[k] = float( cols[2] ) 

#           print ( k, t[k], U[k] )

            k = k + 1 

        f.close()

        return t, U 

    if ( ext == 'txt' ) :

        f = open( fname + ".txt", "r" )

        contents = f.readlines() 

        N = len( contents )

        if ( report ):
            print ( "N = ", N  )

        t = np.zeros (N)
        U = np.zeros (N)

        k = 0

        for line in contents :

            line = line.rstrip()
            cols = line.split()

            t[k] = float( cols[1] ) 
            U[k] = float( cols[2] ) 

#           print ( k, t[k], U[k] )

            k = k + 1 

        f.close()

        return t, U 


    return None 

###############################################################################

def load_U ( l, funame, report=False ):
    if ( report ): 
        print ( "Reading " + funame )
#   return r_MildSol ( OUTDIR + funame, ext='csv' )
    return r_MildSol ( funame, ext='txt' )
#


def plot_Sols_prime ( t, U, te, Ue, labelU="???", labelUe="???", \
    fname=None, report=True ):
#
#  Plotting  
#  U_calc., U_exac
#
    plt.rcParams['figure.figsize'] = [ cm2in(10) , cm2in(10) ]

    fig, axs = plt.subplots()

    axs.plot( te, Ue,
#       color='red', linestyle='dashed', linewidth=2.0,
#       color='black', linestyle='dashed', linewidth=2.0,
        color='black', linewidth=1.5,
        label=labelUe )

    axs.plot( t, U,
#       color='red', linestyle='dashed', linewidth=2.0,
#       color='black', linestyle='dashed', linewidth=2.0,
        color='black', linestyle='dashed', linewidth=1.5,
        label=labelU )

    axs.grid(True)
    axs.set_xlabel( r'$t$' )
#   axs.set_ylabel( r'$u$' )
#   axs.set_title( "Mild solution" )

    axs.legend()

    if ( report ): 
        plt.show()
    else:
        if ( fname is not None ):
            fig.savefig ( fname + '.jpg' , bbox_inches = 'tight' )
            fig.savefig ( fname + '.pdf' , bbox_inches = 'tight' )

    plt.close('all')
    return None 
#


def plot_Sols ( t, U, te, Ue, fname=None, report=True ):
#
#  Plotting  
#  U_calc., U_exac
#

    plot_Sols_prime ( t, U, te, Ue, 
        labelU=r'$ \mathbf{U}_{ N_l } $', 
        labelUe=r'$ \mathbf{U}_{ N_{\rm max} } $',
        fname=fname, report=report )

    return None 
#



#
#   have exact solution 
#
def plot_Sols_1 ( t, U, te, Ue, fname=None, report=True ):
#
#  Plotting  
#  U_calc., U_exac
#

    plot_Sols_prime ( t, U, te, Ue, 
        labelU=r'$ \mathbf{U}_{ N_l } $', 
        labelUe=r'$ \mathbf{U}_{ \rm exact } $',
        fname=fname, report=report )

    return None 
#

#
#   have no exact solution 
#
def plot_Sols_2 ( t, U, te, Ue, fname=None, report=True ):
#
#  Plotting  
#  U_calc., U_exac
#

    plot_Sols_prime ( t, U, te, Ue, 
        labelU=r'$ \mathbf{U}_{ N_l } $', 
        labelUe=r'$ \mathbf{U}_{ N_{\rm max} } $',
        fname=fname, report=report )

    return None 
#

#
#   plotting solution with and without noise 
#
def plot_Sols_3 ( t, U, te, Ue, fname=None, report=True ):
#
#  Plotting  
#  U_N, U_N^delta
#

    plot_Sols_prime ( t, U, te, Ue, 
        labelU=r'$ \mathbf{U}_{ N }^{ \delta } $', 
        labelUe=r'$ \mathbf{U}_{ N } $',
        fname=fname, report=report )

    return None 
#







##################
#
#  basic tool
#
def plot_Err_prime ( a_N, a_Err, \
    label="???", fname=None, report=True, alog=True ):

    plt.rcParams['figure.figsize'] = [ cm2in(10) , cm2in(10) ]

    if ( alog ): 

        fig, axs = plt.subplots()

        tmp = np.log ( np.fmax( a_Err, MachEps ) ) 

        axs.plot( np.log(a_N), tmp,
#           color='red', linestyle='dashed', linewidth=1.0,
            color='black', linestyle='dashed', linewidth=1.0,
            marker='o', markerfacecolor='blue', markersize=8,
            label=label ) 
        axs.grid(True)
        axs.set_xlabel( r'$\log N_l$' )
#       axs.set_ylabel( "$\log ||U_{N_l}-U_{N_{l-1}}||_{F}$ " )
#       axs.set_title( "$\log ||U_{N_l}-U_{N_{l-1}}||$ v.s. $\log N_l$" )
        axs.legend()

        if ( report ): 
            plt.show()
        else:
            if ( fname is not None ):
                fig.savefig ( fname + '_log.jpg' , bbox_inches = 'tight' )
                fig.savefig ( fname + '_log.pdf' , bbox_inches = 'tight' )

    else:

        fig, axs = plt.subplots()

        axs.plot( a_N, a_Err ,
#           color='blue', linestyle='dashed', linewidth=1.0,
            color='black', linestyle='dashed', linewidth=1.0,
            marker='o', markerfacecolor='red', markersize=8,
            label=label )
        axs.grid(True)
        axs.set_xlabel( r'$N_l$' )
#       axs.set_ylabel( "$||U_{N_l}-U_{N_{l-1}}||$ " )
#       axs.set_title( "$||U_{N_l}-U_{N_{l-1}}||$ v.s. $N_l$" )
        axs.legend()

        if ( report ): 
            plt.show()
        else:
            if ( fname is not None ):
                fig.savefig ( fname + '.jpg' , bbox_inches = 'tight' )
                fig.savefig ( fname + '.pdf' , bbox_inches = 'tight' )
#
    plt.close('all')
    return None 
#

#
# | U_Nl - U_exact |:
#
def plot_Err_1 ( lmax, a_N, a_Err, fname=None, report=True, alog=True ):

    if ( alog ): 

        plot_Err_prime ( a_N, a_Err, 
        label=r'$\log\left\Vert\mathbf{U}_{N_l} - \mathbf{U}_{\rm exact}\right\Vert$',
        fname=fname, report=report, alog=True )

    else:

        plot_Err_prime ( a_N, a_Err, 
        label=r'$\left\Vert\mathbf{U}_{N_l}-\mathbf{U}_{N_{l-1}}\right\Vert$', 
        fname=fname, report=report, alog=False )
#
    plt.close('all')
    return None 
#

#
# | UN - UNmax |:
#
def plot_Err_2 ( lmax, a_N, a_Err, fname=None, report=True, alog=True ):

    if ( alog ): 

        plot_Err_prime ( a_N, a_Err, 
        label=r'$\log\left\Vert\mathbf{U}_{N_l} - \mathbf{U}_{N_{\rm max}}\right\Vert$',
        fname=fname, report=report, alog=True )

    else:

        plot_Err_prime ( a_N, a_Err, 
        label=r'$\left\Vert\mathbf{U}_{N_l} - \mathbf{U}_{N_{\rm max}}\right\Vert$', 
        fname=fname, report=report, alog=False )
#
    plt.close('all')
    return None 
#




#################
#
#   legacy  
#
def plot_Err ( lmax, a_N, a_Err, fname=None, report=True, alog=True ):

    if ( alog ): 

        plot_Err_prime ( a_N, a_Err, 
        label=r'$\log\left\Vert\mathbf{U}_{N_l} - \mathbf{U}_{N_{\rm max}}\right\Vert$',
        fname=fname, report=report, alog=True )

    else:

        plot_Err_prime ( a_N, a_Err, 
        label=r'$\left\Vert\mathbf{U}_{N_l} - \mathbf{U}_{N_{\rm max}}\right\Vert$', 
        fname=fname, report=report, alog=False )
#
    plt.close('all')
    return None 
#


def plot_dU ( xi, dU, dUref=None, fname=None, report=True ):
#
#  Plotting  
#  U_calc., U_exac
#
    plt.rcParams['figure.figsize'] = [ cm2in(10) , cm2in(10) ]

    fig, axs = plt.subplots()

    axs.plot( xi, dU,
#       color='red', linestyle='dashed', linewidth=2.0,
#       color='black', linestyle='dashed', linewidth=2.0,
        color='black', linewidth=1.5,
        label=r'$ ^C D^{ \psi, \alpha }_{a+} \mathbf{U}_{N} $' )

    if ( dUref is not None ):
        axs.plot( xi, dUref,
#       color='red', linestyle='dashed', linewidth=2.0,
#       color='black', linestyle='dashed', linewidth=2.0,
        color='black', linestyle='dashed', linewidth=1.5,
        label=r'$ ^C D^{ \psi, \alpha }_{a+} \mathbf{U}_{\rm exact} $' )


    axs.grid(True)
    axs.set_xlabel( r'$t$' )
#   axs.set_ylabel( r'$u$' )
#   axs.set_title( "Mild solution" )

    axs.legend()

    if ( report ): 
        plt.show()
    else:
        if ( fname is not None ):
            fig.savefig ( fname + '.jpg' , bbox_inches = 'tight' )
            fig.savefig ( fname + '.pdf' , bbox_inches = 'tight' )

    plt.close('all')
    return None 
#
#
