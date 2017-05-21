#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun May 21 17:36:04 2017

@author: andrew


Function to calculate Ewald-Kornfeld representation energy for a dipole in
a lattice of dipoles.

See: "Dipolar rotor-rotor interactions in a difluorobenzene molecular rotor 
crystal", Horanksy et al, Phys. Rev. B, 74, 054306 (2006).

"""

from numpy import dot, sqrt
from numpy.linalg import norm
from scipy.special import erfc
from scipy import exp, pi

#==============================================================================
# PSEUDOCODE FIRST
#==============================================================================
def ewaldEnergy():
    
    realSpaceEnergy = realSpaceEnergy() 
    reciprocalSpaceEnergy = reciprocalSpaceEnergy()
    selfEnergy = selfEnergy()
    
    ....