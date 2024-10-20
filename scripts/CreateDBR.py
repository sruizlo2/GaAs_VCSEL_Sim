# -*- coding: utf-8 -*-
"""
Created on Sat Oct 19 18:56:11 2024

@author: user
"""

def Create_DBR(n_pairs, index_in, index_1, index_2, thickness, index_out=1):
    # Entry medium
    layers = [(index_in, 0)]
    # Alteranting layers
    for i in range(n_pairs):
        layer_1 = (index_1, thickness)
        layer_2 = (index_2, thickness)
        layers.append(layer_1)
        layers.append(layer_2)
    # Exit medium
    layers.append((index_out, 0))
    
    return layers