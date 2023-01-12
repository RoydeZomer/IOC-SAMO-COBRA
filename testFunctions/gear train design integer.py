# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 14:03:37 2022

@author: r.dewinter
"""


Minimize f1(x) = abs( 6.931 − x3*x1/x4*x2 )
Minimize f2(~x) = max(x1, x2, x3, x4),
Subject to f1(~x)/6.931 ≤ 0.5,
12 ≤ x1, x2, x3, x4 ≤ 60,
all xi’s are integers