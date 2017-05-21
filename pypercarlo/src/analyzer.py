#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun May 21 11:45:15 2017

@author: andrew

Data analysis and plotting wrapper class

"""
import pyalps.plot
import matplotlib.pyplot as plt

class Analyzer(object):
    
    def __init__(self, dataLocationPattern=None, model=None):
        self.dataLocationPattern = dataLocationPattern
        self.resultFiles = None
        self.data = None
        self.model = model
        
    def getFiles(self, dataLocationPattern):
        if self.dataLocationPattern is None:
            self.dataLocationPattern = dataLocationPattern
            
        self.resultFiles = pyalps.getResultFiles(pattern=dataLocationPattern)
    
    def printObservables(self):
        print pyalps.loadObservableList(self.resultFiles)
    
    def loadMeasurements(self, measurements):
        self.data = pyalps.loadMeasurements(self.resultFiles, measurements)
    
    def plotMeasurement(self, xaxis, measurement):
        plotData = pyalps.collectXY(self.data, xaxis, measurement)
        plt.figure()
        pyalps.plot.plot(plotData)
        plt.title(self.model)
        plt.show()