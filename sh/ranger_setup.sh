module swap pgi intel
module load cmake
module load boost

setenv BOOST_INCLUDEDIR $TACC_BOOST_INC
setenv BOOST_LIBRARYDIR $TACC_BOOST_LIB
