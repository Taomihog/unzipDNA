# -*- coding: utf-8 -*-
"""
Created on Sun Aug 13 20:49:16 2023

@author: xiang
"""
import matplotlib.pyplot as plt
import pandas as pd

pyResult = pd.read_csv("standard.txt",delimiter='\t')
cppResult = pd.read_csv("out_exampleData.txt")
print(pyResult.columns)
print(cppResult.columns)

fig,ax = plt.subplots()
ax.plot(pyResult['# zMean'],pyResult['FMean'], "-b", label = "standard curve calculated from original python code")
ax.plot(cppResult['extension (nm)'],cppResult['average force (pN)'], "-r", label = "curve by this program")
plt.show()