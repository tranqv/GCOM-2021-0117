import os
 
# get current directory

here = os.getcwd()
print( "Current Directory:", here )

print( "os.pardir:",  os.pardir ) 


# prints parent directory

ppat = os.path.abspath( os.path.join( here, os.pardir )) 

print( "Parent  Directory:", ppat )


OUTDIR = 'xout'
OUTDIR = os.path.join( ppat, OUTDIR ) 

print( "outdir :", OUTDIR )


