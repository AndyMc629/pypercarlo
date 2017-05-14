#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun May 14 10:48:15 2017

@author: andrew
"""

import pyalps.hdf5 as hdf5
import sys, time, traceback, getopt

import ising

if __name__ == '__main__':

    sim = ising.sim({
        'L': 100,
        'THERMALIZATION': 100,
        'SWEEPS': 1000,
        'T': 2
    })

    outfile = 'test'
    limit = 1000
    
    if limit == 0:
        sim.run(lambda: False)
    else:
        start = time.time()
        sim.run(lambda: time.time() > start + float(limit))

    with hdf5.archive(outfile[0:outfile.rfind('.h5')] + '.clone0.h5', 'w') as ar:
        ar['/'] = sim

    results = sim.collectResults() # TODO: how should we do that?
    for key, value in results.iteritems():
        print "{}: {}".format(key, value)

    with hdf5.archive(outfile, 'w') as ar:
        ar['/parameters'] = sim.parameters
        ar['/simulation/results'] = results