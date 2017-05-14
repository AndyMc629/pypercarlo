#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun May 14 11:14:30 2017

@author: andrew

main program for running sims.
"""
import pyalps
import simulation

import matplotlib.pyplot as plt
import pyalps.plot

if __name__ == '__main__':
    L = 4    # Linear lattice size
    N = 5000    # of simulation steps

    print '# L:', L, 'N:', N

    # Scan beta range [0,1] in steps of 0.1
    for beta in [0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.]:
        for l in [4,6,8]:
            print '-----------'
            print 'beta =', beta
            sim = simulation.Simulation(beta,l)
            sim.run(N/2,N)
            sim.save('ising.L_'+str(l)+'beta_'+str(beta)+'.h5')
    
    #how to calculate the Binder Ratio within Python:
    infiles=pyalps.getResultFiles(pattern='ising.L')

    data = pyalps.loadMeasurements(pyalps.getResultFiles(pattern='ising.L*'),['E','m^2', 'm^4'])
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
    plt.title('2D Ising model')
    plt.legend()
    plt.savefig('Binder.pdf')
    plt.close()
    
    plt.figure()
    pyalps.plot.plot(eplot)
    plt.xlabel('$\\beta$')
    plt.ylabel('E/site')
    plt.title('2D Ising model')
    plt.legend()
    plt.savefig('Energy.pdf')
    plt.close()
    print('Accepted moves:', sim.accepted)