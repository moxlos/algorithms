#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: lefteris

@subject: supply-demand optimization class
"""


import pandas as pd
import numpy as np
import pulp


class SupplyDemand:
    def __init__(self,C, D,S,c):
        #All arrays
        self.C = C #transportation cost
        self.D = D #demand cost
        self.S = S #supply cost
        self.c = c #oporational cost
        self.nVar = c.shape[0]
        self.mVar = D.shape[0]
    
        
    def get_variables(self):
        xshape = (range(self.nVar), range(self.mVar))
        x = pulp.LpVariable.dicts("X", xshape,lowBound = 0,  cat='Continuous')
        y = pulp.LpVariable.dicts("Y", range(self.nVar), cat=pulp.LpBinary)
        return x, y

    def get_solution(self):
        x, y = self.get_variables()
        
        prob = pulp.LpProblem("distribution_opt", pulp.LpMinimize)
        
        objective_function = sum([ self.C[0,m_idx] * x[0][m_idx] for m_idx in range(self.mVar)]) + self.c[0] * y[0]
        for n_idx in range(1,self.nVar):
            objective_function += sum([ self.C[n_idx,m_idx] * x[n_idx][m_idx] for m_idx in range(self.mVar)]) + self.c[n_idx] * y[n_idx]
        
        prob += objective_function
        
        #Constraints
        for n_idx in range(self.nVar):
            prob += sum(x[n_idx][m_idx] for m_idx in range(self.mVar)) <= self.S[n_idx] * y[n_idx]
        for m_idx in range(self.mVar):
            prob += sum(x[n_idx][m_idx] for n_idx in range(self.nVar)) == self.D[m_idx]
    
    
        prob.solve(pulp.apis.PULP_CBC_CMD(msg=False))
        return prob


if __name__ == "__main__":
    df = pd.read_csv("data/demand.csv", sep='\s+')
    C = np.array(df.iloc[:-1,1:-2])
    D = np.array(df.iloc[-1,1:-2])
    S = np.array(df.iloc[:-1,-2])
    c = np.array(df.iloc[:-1,-1])
    df
    
    supp = SupplyDemand(C,D,S,c)
    prob = supp.get_solution()
    prob
    
    X = []
    Y = []
    for v in prob.variables():
        if v.name[0] == 'X':
            X.append(v.varValue)
        else:
            Y.append(v.varValue)
            
    X = np.reshape(np.array(X), (S.shape[0], D.shape[0]))
    Y = np.array(Y)
    X