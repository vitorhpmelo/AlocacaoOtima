#%%
import numpy as np


Hfat=np.array([[1,2,3],[1,2,3],[1,2,3]])

#%%
i=0
j=2

Hfat[:,[j,i]]=Hfat[:,[i,j]]
# %%
