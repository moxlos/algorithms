#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: lefteris

@subject: disease dynamics
"""


import sympy as sym

class Disease:
    def __init__(self, per_person_cost):
        #Initialize symbolic variables
        self.per_person_cost = per_person_cost
        self.t = sym.Symbol("t")
        self.k = sym.Symbol("k")
        self.I_0 = sym.Symbol("I_0")
        self.I = sym.Function("I")
    
    def get_equation(self, k=None):
        if k is None:
            k = self.k
        return sym.Eq(sym.Derivative(self.I(self.t), self.t), - k * self.I(self.t))
    def get_solution(self,I_0=None, k=None):
        if k is None:
            k = self.k
        if I_0 is None:
            I_0 = self.I_0
        eq = self.get_equation(k=k)
        return sym.dsolve(eq, self.I(self.t), ics={self.I(0): I_0})
    def get_cost(self,I_0=None,k=None, cure_cost=0):
        
        if k is None:
            k = self.k
        if I_0 is None:
            I_0 = self.I_0

        I_sol = self.get_solution(I_0=I_0, k=k)
        area = sym.integrate(I_sol.rhs, (self.t, 0, sym.oo))
        productivity_cost = area * self.per_person_cost
        total_cost_of_cure = cure_cost * I_0
        return productivity_cost + total_cost_of_cure
