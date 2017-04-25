# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 17:58:45 2017

@author: andrew



PyPerCarlo (Python Perovskite Monte Carlo): HyPerCarlo gone pythonic!

"""
import configuration as cfg
import lattice
        
if __name__ == "__main__":
    conf = cfg.Configuration()
    #conf.loadJSON('./input/input.json')
    test = {"settings": [
    {   
        "dimension": 2,
        "latticeXdim": 10, 
        "latticeYdim": 10, 
        "periodicBoundaryConditions": "True",
        "model": "dipole-dipole",
        "degOfFreedom": 6,
        "runMode": "anneal",
        "Tmin": 2.0,
        "Tmax": 0.0,
        "initialiseLatticeType": "ferro"
    }]}

    conf.createJSON(test, './input/test.json')
    
    lat = lattice.Lattice(lengthX=2, lengthY=2)
    lat.initialiseLattice('ising-antiferro')
    
    #lat2 = lattice.Lattice(lengthX=3, lengthY=2, lengthZ = 4)
    