#!/usr/bin/env python
# coding: utf-8

# In[96]:


import ginga
from ginga.util import iqcalc
import numpy as np
from ginga.web.pgw import ipg


# In[100]:


use_opencv = False

server = ipg.make_server(host='localhost', port=8701, use_opencv=use_opencv)
server.start(no_ioloop=True)
v1 = server.get_viewer('v1')
v1.open()
v1.load('bgc.fits')
