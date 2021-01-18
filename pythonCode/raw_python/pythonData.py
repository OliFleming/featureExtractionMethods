#!/usr/bin/env python
# coding: utf-8

# In[41]:


import numpy as np
import matplotlib.pyplot as plt
from getzcfeat import *


# In[42]:


data = np.genfromtxt('smallnoise_sin.csv', delimiter=',')
t = data[:,0]
s = data[:,1]
np.matrix.view(data)


# In[43]:


get_ipython().run_line_magic('matplotlib', 'inline')
plt.plot(t,s)


# In[44]:


size = np.shape(s)
print(size)


# In[45]:


result = getzcfeat(s,0.01,2000,2000)
print(result)


# In[ ]:




