# -*- coding: utf-8 -*-
"""
Created on Sat May 15 01:35:57 2021

@author: Jitao Zhang
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


df = pd.read_excel(r"E:\work in CNGB\DNA_storage\image_processing\Image_processing\Bmp\lena512_16_16_matrix\size-test\psnr_diff_mapping.xlsx",
                   sheet_name="Summary").iloc[:20,1:6]
#fig = plt.figure(figsize=(8,8))
#df.plot.figsize(20,20)



df.plot.box(figsize=(10,8),fontsize=14,)
#plt.title("match rate of images in genome 4^{}".format(kmer),fontsize=25)
plt.xlabel("kmer",fontsize=20)
plt.ylabel("match rate",fontsize=20)
plt.grid(linestyle="--",alpha=0.3)
plt.ylim([22,32])
#plt.plot([3.6,4.4],[1.0,1.0],"--",color = "red",label="The threshold at which a file can be compressed")
#plt.plot([4.6,5.4],[0.5,0.5],"--",color = "red",label="The threshold at which\na file can be compressed")
#plt.plot([5.6,6.4],[0.25,0.25],"--",color = "red")
#plt.plot([6.6,7.4],[0.17,0.17],"--",color = "red")
 
#plt.legend(fontsize=14)
#average_df = pd.read_excel(r"E:\work in CNGB\DNA_storage\image_processing\Image_processing\Bmp\lena512_16_16_matrix\size-test\psnr_diff_mapping.xlsx",
                   #sheet_name="Summary").iloc[20:21,1:6]
#plt.scatter(average_df.columns.values, average_df.values)

#plt.scatter([7,8,9,10,11], [25,25,25,25,25])
plt.show()
#plt.savefig("match rate of images in genome 4^{}".format(kmer))
 