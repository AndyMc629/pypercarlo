# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 19:42:33 2017

@author: andrew

Lattice class, which represents the physical system.
"""
import math
import matplotlib.pyplot as plt
import numpy as np
import csv

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
    
class latticeDataReader(object):
    def __init__(self, path):
        self.path = path
        self.data = []
        self.lx = 0
        self.ly = 0
        self.lz = 0
        with open(self.path) as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            for row in reader:
                self.data.append(row)
                if(row[0]>self.lx):
                    self.lx = row[0]
                if(row[1]>self.ly):
                    self.ly = row[1]
                if(row[2]>self.lz):
                    self.lz = row[2]
    
    def returnLattice(self):
        dimension = len(self.data[0])-3 # always have coord's plus 3d vector.
        if dimension == 3:
            lattice = Lattice(lengthX = self.lx, 
                              lengthY = self.ly, 
                              lengthZ = self.lz)  
            
        return lattice
            
            
class Lattice(object):
    def __init__(self,  lengthX=0, lengthY=0, lengthZ=0, lattice = None):
        
        if(lengthX==0 and lengthY==0 and lengthZ==0):
            raise ValueError("must supply at least one dimension length or data path!")
    
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
            self.lattice = []
            for i in range(0, self.volume): #volume = int
                self.lattice.append(spin())            
    
    def getCoordinates(self, idx):
        x = idx%self.lengthX
        y = int(math.floor((idx-x)/self.lengthX))
        return y, x # C++ 2d array coordinates.
    """
    def initialiseFromData(self, data, dim): # data = list of lists here
        if dim == 2:
            for i in range(0, len(data)):
                for j in range(0, len(data[i])):
    """                
    
    def initialiseFerro(self):
        for spin in self.lattice:
            spin.sy = 1.0
    
    def initialiseAntiFerro(self):
        if self.dimension == 1:
            for i in range(0, len(self.lattice), 2):
                self.lattice[i].sy = -1.0
            for i in range(1, len(self.lattice), 2):
                self.lattice[i].sy = 1.0
        
        elif self.dimension == 2:
            for y in range(0, self.lengthY):
                for x in range(0, self.lengthX):
                    self.lattice[y*self.lengthX+x].sy = (-1)**(x+y) 
    
    def initialiseColAntiFerro(self):
        if self.dimension == 1:
            raise ValueError("cannot have colAntiFerro for 1D lattice.")
        elif self.dimension == 2:
            for y in range(0, self.lengthY):
                for x in range(0, self.lengthX):
                    self.lattice[y*self.lengthX+x].sy = (-1)**y
            
            
    def printLattice(self):
        if self.volume<100:        
            for i in range(0, len(self.lattice)):
                print self.lattice[i].asList()
                
    def savecsv(self, path=None):
        if path==None:
            raise TypeError("save() requires path (string) as arg.")
        elif isinstance(path, str):
            saveFile = open(path, "w")
            for idx, spin in enumerate(self.lattice):
                y, x = self.getCoordinates(idx)
                # y and x are in correct order, use C++ array notation.
                saveFile.write(
                               str(y)+","+\
                               str(x)+","+\
                               str(self.lattice[idx].sx)+","+\
                               str(self.lattice[idx].sy)+","+\
                               str(self.lattice[idx].sz)+"\n"
                               )
            saveFile.close()
            
    def visualise(self, path=None):
        
        if self.dimension == 2:
            
            X, Y = np.mgrid[0:self.lengthX, 0:self.lengthY]
            
            U = np.zeros((self.lengthX,self.lengthY))
            V = np.zeros((self.lengthX,self.lengthY))
            
            for i in range(0, self.volume):
                U[self.getCoordinates(i)] = self.lattice[i].sx
                V[self.getCoordinates(i)] = self.lattice[i].sy
            
            #plt.axes([0.025, 0.025, 0.95, 0.95])
            plt.axis('equal')
            
            plt.quiver(X, Y, U, V, 0.5, alpha=.5)
            plt.quiver(X, Y, U, V, 0.5, edgecolor='k', facecolor='None', linewidth=.5)
            
            plt.xlim(-1, self.lengthX)
            plt.xticks(())
            plt.ylim(-1, self.lengthY)
            plt.yticks(())
            
            if path == None:
                plt.show()
            else:
                plt.savefig(path)
            plt.close()
        
    
        
            
            
            
        
                
                
        
            
            
            
        
                
        
