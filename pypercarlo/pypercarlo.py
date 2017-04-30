# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 17:58:45 2017

@author: andrew



PyPerCarlo (Python Perovskite Monte Carlo): HyPerCarlo gone pythonic!

"""
import configuration as cfg
import lattice
    
"""
Will put everything that is in main into a def main(): func
soon. 

"""
import csv

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
    
    lat = lattice.Lattice(lengthX=10, lengthY=10)
    lat.initialiseFerro()
    lat.visualise('./matt/ferro.pdf')
    lat.savecsv('./matt/ferro.csv')
    
    lat.initialiseColAntiFerro()
    lat.visualise('./matt/colAntiFerro.pdf')
    lat.savecsv('./matt/colAntiFerro.csv')
    
    lat.initialiseAntiFerro()
    lat.visualise('./matt/antiFerro.pdf')
    lat.savecsv('./matt/antiFerro.csv')
    
    """
    data = []
    with open('test.csv') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            data.append(row)
            
    print data
    """
