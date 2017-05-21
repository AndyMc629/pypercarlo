#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun May 14 11:14:30 2017

@author: andrew

main program for running sims.
"""
# pyalps modules
import pyalps
import pyalps.plot
# my modules
import simulation
import analyzer 
# standard modules
import matplotlib.pyplot as plt
import numpy as np


if __name__ == '__main__':
    
    #==============================================================================
    # SET UP AND RUN SIMULATION
    #==============================================================================
    Lmax = 20    # Linear lattice size, assume cubic
    N = 5000    # of simulation steps

    print '# Lmax:', Lmax, 'N:', N
    model = 'dipole-dipole6'

    # Scan beta range [0,1] in steps of 0.1
    #for beta in [0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.]:
    #for beta in [0., .2, .4, .6, .8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]: 
    #for beta in [2.0, 2.5, 3.0, 3.5, 4.0]:    
    for beta in np.arange(1.0, 8.0, 1.0):
        #for l in range(4,Lmax,2):
            l=Lmax
            print '-----------'
            print 'beta =', beta
            print 'l =', l
            sim = simulation.Simulation(beta,l,model)
            sim.run(N,4*N)
            sim.save('data/'+str(model)+'.L_'+str(l)+'beta_'+str(beta)+'.h5')
    
    print('Accepted moves:', sim.accepted)
    
    
    analysis = analyzer.Analyzer(model=model)
    
    #==============================================================================
    # DATA ANALYSIS     
    #==============================================================================
    #how to calculate the Binder Ratio within Python:
    resultsDir = 'data/'
    dataLocationPattern = 'data/'+str(model)
    
    infiles=pyalps.getResultFiles(pattern=dataLocationPattern)

    data = pyalps.loadMeasurements(pyalps.getResultFiles(pattern=dataLocationPattern+'*'),['E','m^2', 'm^4'])
    m2 = pyalps.collectXY(data,x='BETA',y='m^2',foreach=['L'])
    m4 = pyalps.collectXY(data,x='BETA',y='m^4',foreach=['L'])
    E = pyalps.collectXY(data,x='BETA',y='E',foreach=['L'])
    
    m2plot = []
    eplot = []
    u4=[]
    for i in range(len(m2)):
        d = pyalps.DataSet()
        d.propsylabel='U4'
        d.props = m2[i].props
        d.x= m2[i].x
        d.y = m4[i].y/m2[i].y/m2[i].y
        u4.append(d)
        
        e = pyalps.DataSet()
        e.plopsylabel='E/site'
        e.props = E[i].props
        e.x = E[i].x
        e.y = E[i].y
        eplot.append(e)
    
    plt.figure()
    pyalps.plot.plot(u4)
    plt.ylim([1,3.5])
    plt.xlabel('Inverse Temperature $\\beta$')
    plt.ylabel('Binder Cumulant U4 $g$')
    plt.title('2D '+str(model)+' model')
    plt.legend()
    plt.savefig(str(resultsDir)+str(model)+'Lmax_'+str(Lmax)+'Binder.pdf')
    plt.close()
    
    plt.figure()
    pyalps.plot.plot(eplot)
    plt.xlabel('$\\beta$')
    plt.ylabel('E/site')
    plt.title('2D '+str(model)+' model')
    plt.legend()
    plt.savefig(str(resultsDir)+str(model)+'Lmax_'+str(Lmax)+'Energy.pdf')
    plt.close()
    
    plt.figure()
    pyalps.plot.plot(E)
    plt.show()
    
    sim.spinsToFile('finalSpins.csv')
    