# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 19:42:33 2017

@author: andrew

Lattice class, which represents the physical system.
"""
import math

class spin(object):
    def __init__(self, sx=0.0, sy=0.0, sz=0.0):
        self.sx = sx
        self.sy = sy
        self.sz = sz
    
    def magnitude(self):
        return math.sqrt( 
        (self.sx*self.sx)+\
        (self.sy*self.sy)+\
        (self.sz*self.sz)
        )
    
    def asList(self):
        return [self.sx, self.sy, self.sz]
    
        
class Lattice(object):
    def __init__(self, lattice = None, lengthX=0, lengthY=0, lengthZ=0):
        
        dimcheck = 0        
        if lengthX>0:
            dimcheck+=1
            if lengthY>0:
                dimcheck+=1
                if lengthZ>0:
                    dimcheck+=1
                    
        # spatial dimension
        if dimcheck>3 or dimcheck<1: 
            raise ValueError("Lattice: 'dimension' must be equal to 1,2 or 3.")
        else:
            self.dimension = dimcheck
            
        # lattice dimensions
        if not (isinstance(lengthX, int) 
        and isinstance(lengthY, int) 
        and isinstance(lengthZ, int)):
                raise ValueError("Lattice: one or more dimensions of the lattice\
                                    are not integers.")
        else:
            self.lengthX = lengthX
            self.lengthY = lengthY
            self.lengthZ = lengthZ
            if self.dimension == 1:
               self.volume = self.lengthX
            elif self.dimension == 2:
                self.volume = self.lengthX*self.lengthY
            elif self.dimension == 3:
                self.volume = self.lengthX*self.lengthY*self.lengthZ
                
        # default lattice    
        if lattice is None:
            defaultSpin = spin()
            self.lattice = [defaultSpin]*self.volume            
    
    
    def initialiseLattice(self, keyword):
        if keyword == 'ising-ferro':
            for spin in self.lattice:
                spin.sz = 1
        elif keyword == 'ising-antiferro':
            up = 1.0            
            for spin in self.lattice:
                up *= -1.0
                spin.sz = up
                print spin.sz
                
    def printLattice(self):
        if self.volume<100:        
            for spin in self.lattice:
                print spin.asList()
                
    
                
                
        
            
            
            
        
                
        
