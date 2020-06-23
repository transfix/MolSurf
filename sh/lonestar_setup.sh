UPDATE:

module load boost
module load cmake

export BOOST_LIBRARYDIR=$TACC_BOOST_LIB
export BOOST_INCLUDEDIR=$TACC_BOOST_INC




-------------------------------------

# this doesn't seemt to run correctly as a script but these are the commands that I need

# don't forget to disable building the CoarseGrain directory... nfft doesn't seem to be available.

module load Linux
module load boost
module load cmake

setenv CPPFLAGS "-I/usr/include -I$TACC_BOOST_INC/boost-1_39"
setenv LDFLAGS "-L/usr/lib -L$TACC_BOOST_INC"
