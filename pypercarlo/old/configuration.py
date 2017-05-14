# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 18:01:22 2017

@author: andrew

Configuration class.

"""
import json

class Configuration(object):
    """ 
    A class to read in, manipulate and store configuration information 
    for pypercarlo.    
    """
    def __init__(self):
        self.settings = None
    
    def loadJSON(self, pathToFile):
        with open(pathToFile, 'r') as f:        
            jsonString = json.load(f)
        if 'settings' in jsonString.keys():        
            self.settings =  jsonString['settings'][0]
        else:
            raise ValueError("Input json file must contain 'settings' key value pair.")
        f.close()
        
    def createJSON(self, settingsDict, pathToFile):
        with open(pathToFile, 'w') as f:
            json.dump(settingsDict, f, sort_keys=True, indent=2)
        f.close()              
        