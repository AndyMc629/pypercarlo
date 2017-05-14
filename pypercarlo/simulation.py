#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun May 14 11:09:23 2017

@author: andrew

Attempt at building a generic simulation class using the ALPS library.

"""
import math
import pyalps
import pyalps.alea as alpsalea
import pyalps.pytools as alpstools

import pyalps.plot

import spin

class Simulation:
    # Seed random number generator: self.rng() will give a random float from the interval [0,1)
    rng = alpstools.rng(42)
    
    def __init__(self,beta,L):
        self.L = L
        self.beta = beta
        
        # Init exponential map
        self.exp_table = dict()
        for E in range(-8,10,4): 
          self.exp_table[E] = math.exp(2*beta*E)
          #self.exp_table[E] = math.exp(-beta*E) #math.exp(beta*E)
          
        # Init random spin configuration
        self.spins = [ [spin.Spin(sx=0.0,sy=0.0,sz=2*self.randint(2)-1) for j in range(L)] for i in range(L) ]
        
        # Init observables
        self.energy = alpsalea.RealObservable('E')
        self.magnetization = alpsalea.RealObservable('m')
        self.abs_magnetization = alpsalea.RealObservable('|m|')
        self.magnetization_2 = alpsalea.RealObservable('m^2')
        self.magnetization_4 = alpsalea.RealObservable('m^4')
        self.accepted = 0
        
    def save(self, filename):
        pyalps.save_parameters(filename, {'L':self.L, 'BETA':self.beta, 'SWEEPS':self.n, 'THERMALIZATION':self.ntherm})
        self.abs_magnetization.save(filename)
        self.energy.save(filename)
        self.magnetization.save(filename)
        self.magnetization_2.save(filename)
        self.magnetization_4.save(filename)
        
    def run(self,ntherm,n):
        # Thermalize for ntherm steps
        self.n = n
        self.ntherm = ntherm
        while ntherm > 0:
            self.step()
            ntherm = ntherm-1
            
        # Run n steps
        while n > 0:
            self.step()
            self.measure()
            n = n-1
            
        # Print observables
        print '|m|:\t', self.abs_magnetization.mean, '+-', self.abs_magnetization.error, ',\t tau =', self.abs_magnetization.tau
        print 'E:\t', self.energy.mean, '+-', self.energy.error, ',\t tau =', self.energy.tau
        print 'm:\t', self.magnetization.mean, '+-', self.magnetization.error, ',\t tau =', self.magnetization.tau
    

    def energyLocal(self, spin, i, j):
        e = self.spins[(i-1+self.L)%self.L][j].sz +\
        self.spins[(i+1)%self.L][j].sz +\
        self.spins[i][(j-1+self.L)%self.L].sz +\
        self.spins[i][(j+1)%self.L].sz
        e *= -spin.sz
        return e
    
    def flip(self, spinToFlip):
        return spin.Spin(sx=spinToFlip.sx, sy=spinToFlip.sy, sz=-spinToFlip.sz)
        
    def step(self):
        for s in range(self.L*self.L):
            # Pick random site k=(i,j)
            i = self.randint(self.L)
            j = self.randint(self.L)
            
            newSpin = self.flip(self.spins[i][j])
            
            de = self.energyLocal(newSpin, i, j) - self.energyLocal(self.spins[i][j], i, j)
            #de = self.deltaEnergy(i, j)
            
            # Flip s_k with probability exp(2 beta e)
            #if e > 0 or self.rng() < self.exp_table[e]:
            #    self.spins[i][j] = -self.spins[i][j]
            if de < 0 or self.rng() < math.exp(-self.beta*de): #de<=0 breaks the binder plot ...
                self.spins[i][j] = newSpin
                self.accepted += 1
                
    def measure(self):
        E = 0.    # energy
        M = 0.    # magnetization
        for i in range(self.L):
            for j in range(self.L):
                E -= self.spins[i][j].sz * (self.spins[(i+1)%self.L][j].sz + self.spins[i][(j+1)%self.L].sz)
                M += self.spins[i][j].sz
                
        # Add sample to observables
        self.energy << E/(self.L*self.L)
        self.magnetization << M/(self.L*self.L)
        self.abs_magnetization << abs(M)/(self.L*self.L)
        self.magnetization_2 << (M/(self.L*self.L))*(M/(self.L*self.L))
        self.magnetization_4 << (M/(self.L*self.L))*(M/(self.L*self.L))*(M/(self.L*self.L))*(M/(self.L*self.L))
        
    # Random int from the interval [0,max)
    def randint(self,max):
        return int(max*self.rng())
