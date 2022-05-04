#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  4 17:13:05 2022

@author: matthew
"""

import matplotlib.pyplot as plt

from shapely.geometry import Polygon
from shapely.ops import unary_union

test1 = Polygon([(1,1), (1,2), (2,2), (2,1), (1,1)])
test2 = Polygon([(0,0), (0,1.5), (1.5,1.5), (1.5,0), (0,0)])
    
test3 = unary_union([test1, test2])


#plot it
fig, axes = plt.subplots(1,2)
xs, ys = test1.exterior.xy
axes[0].fill(xs, ys, alpha=0.5, fc='r', ec='none')
xs, ys = test2.exterior.xy
axes[0].fill(xs, ys, alpha=0.5, fc='y', ec='none')

xs, ys = test3.exterior.xy
axes[1].fill(xs, ys, alpha=0.5, fc='y', ec='none')