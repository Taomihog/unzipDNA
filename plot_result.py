# -*- coding: utf-8 -*-
"""
Created on Sun Aug 13 20:49:16 2023

@author: taomihog
"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

pyResult = pd.read_csv("standard.txt",delimiter='\t')
cppResult = pd.read_csv("out_exampleData.csv")
print(pyResult.columns)
print(cppResult.columns)

fig,ax = plt.subplots()
ax.plot(pyResult['# zMean'],pyResult['FMean'], "-b", label = "by original python code (may not be 100% correct)")
ax.plot(cppResult['DNA extension (nm)'],cppResult['average force (pN)'], "-.r", label = "by this program")
plt.ylim(0,30)
plt.xlim(600,5500)
plt.show()




        
            
        