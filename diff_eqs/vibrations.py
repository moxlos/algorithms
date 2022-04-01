#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: lefteris

@subject: oscilations
"""

import sympy as sym
from matplotlib import pyplot as plt
import numpy as np


class Oscilations:
    def __init__(self, F="0"):
        #Initialize symbolic variables
        #self.per_person_cost = per_person_cost
        self.t = sym.Symbol("t")
        self.k = sym.Symbol("k")
        self.m = sym.Symbol("m")
        self.gam = sym.Symbol("gam")
        self.u_0 = sym.Symbol("u_0")
        self.du_0 = sym.Symbol("du_0")
        self.u = sym.Function("u")
        self.F = sym.Derivative(F,'t',0).doit()
        self.du = sym.diff(self.u(self.t), self.t)
    
    def get_equation(self, k=None, m=None,gam=None):
        if k is None:
            k = self.k
        if gam is None:
            gam = self.gam
        if m is None:
            m = self.m
        #TODO: Add forced vibrations
        return sym.Eq(m*sym.Derivative(self.u(self.t), self.t,2) + 
                      gam *sym.Derivative(self.u(self.t), self.t) +
                      k* self.u(self.t), self.F)
    
    def get_solution(self,u_0=None,
                     du_0=None, k=None, m=None,gam=None):
        if k is None:
            k = self.k
        if gam is None:
            gam = self.gam
        if m is None:
            m = self.m
        if du_0 is None:
            du_0  = self.du_0 
        if u_0 is None:
            u_0  = self.u_0 
            
        eq = self.get_equation(k=k, m=m,gam=gam)
        return sym.dsolve(eq, self.u(self.t), ics={self.u(0): u_0,
                                                   self.du.subs(self.t,0): du_0})
    
    def plot_solution(self, u_0,du_0,k,m,gam, N =10):
        eq = self.get_equation(k=k, m=m,gam=gam)
        sol = self.get_solution(u_0=u_0,du_0=du_0,
                     m=m,k=k,gam=gam)
        
        func = sym.lambdify(self.t, sol.rhs,'numpy')
        xvals = np.arange(0,N,.1)
        yvals = func(xvals)

        # make figure
        fig, ax = plt.subplots(1,1)
        ax.axhline(0, color='black')
        ax.plot(xvals, yvals)
        ax.set_xlabel('t')
        ax.set_ylabel('u(t)')
        plt.show()
    






