CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

ADD_SUBDIRECTORY (ARM)
ADD_SUBDIRECTORY (OptimizerLib)

# this is redundant in MolSurf but
# TexMol needs this to work...
ADD_SUBDIRECTORY (CoarseGrain)

