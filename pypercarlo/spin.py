#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun May 14 15:13:59 2017

@author: andrew

Spin class (no pun intended).
"""
import math

class Spin(object):
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
