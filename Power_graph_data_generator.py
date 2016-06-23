# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 15:04:14 2016

@author: patrick
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 10:50:52 2016

Purpose: Generate likelihood ratio tests for Model 1 versus Model 2 (consistent effect) and Model 2 versus Model 3 (interaction effect)

"""

from scipy.stats import *
import csv
import math
import time
import timeit
from random import *
import numpy as np
import matplotlib.pyplot as plt

  
ime13r=19.0
ime13a=21.0
iml13r=21.0
iml13a=19.0

ime14r=20.0
ime14a=20.0
iml14r=20.0
iml14a=20.0


for i in range(0,22):
    if i>0 and i<22:
        ime13r+=1.0
        iml13a+=1.0
        ime13a-=1.0
        iml13r-=1.0
                                                             
        pime13=ime13r/(ime13r+ime13a)
        piml13=iml13r/(iml13r+iml13a)

        vime13=(0.025+(1/(ime13r+ime13a)))
        xime13=2*math.asin(pime13**0.5)
        viml13=(0.025+(1/(iml13r+iml13a)))
        ximl13=2*math.asin(piml13**0.5)
            
        pime14=0.5
        piml14=0.5
        vime14=(0.025+(1/(ime14r+ime14a)))
        xime14=2*math.asin(pime14**0.5)
        viml14=(0.025+(1/(iml14r+iml14a)))
        ximl14=2*math.asin(piml14**0.5)
                
            
        """Model 5 = Year and Pop effect.  No differences among bulks.""" 
        ll_Y_im13=0.0
        ll_Y_im14=0.0
        ll_Y_im=0.0

        xim13=[xime13,ximl13]
        xim14=[xime14,ximl14]
        vim13=[vime13,viml13]
        vim14=[vime14,viml14]
        wim13=[1/vime13,1/viml13]
        wim14=[1/vime14,1/viml14]
        uim13=sum([xx*(ww/sum(wim13)) for xx,ww in zip(xim13,wim13)])
        uim14=sum([xx*(ww/sum(wim14)) for xx,ww in zip(xim14,wim14)])
        
        df_B_im=2              
                
        for j,xx in enumerate(xim13):
            ll_Y_im13+=(-(xx-uim13)**2)/(2*vim13[j])
        ll_Y_im+=ll_Y_im13
        
        for j,xx in enumerate(xim14):
            ll_Y_im14+=(-(xx-uim14)**2)/(2*vim14[j])
        ll_Y_im+=ll_Y_im14
                
                
        """Model 6b = Year effect and consistent effect of bulk"""
        ll_YCB_im=0.0
        uim13 = -((-viml13*xime13 - vime13*ximl13 - vime14*ximl13 - viml14*ximl13 + viml13*xime14 - viml13*ximl14)/(vime13 + viml13 + vime14 + viml14))
        uim14 = -((viml14*xime13 - viml14*ximl13 - viml14*xime14 - vime13*ximl14 - viml13*ximl14 - vime14*ximl14)/(vime13 + viml13 + vime14 + viml14))                
        aim = -((-vime14*xime13 - viml14*xime13 + vime14*ximl13 + viml14*ximl13 - vime13*xime14 - viml13*xime14 + 
        vime13*ximl14 + viml13*ximl14)/(vime13 + viml13 + vime14 + viml14))

        df_YCB_im=1                  

        ll_YCB_im+=(-(xime13-(uim13+aim))**2)/(2*vime13)
        ll_YCB_im+=(-(xime14-(uim14+aim))**2)/(2*vime14)
        ll_YCB_im+=(-(ximl13-(uim13))**2)/(2*viml13)
        ll_YCB_im+=(-(ximl14-(uim14))**2)/(2*viml14)
        
        lrt_B_im_M12=-2*(ll_Y_im-ll_YCB_im)#Compare model 1 (no difference among early versus late) to model 2 (consistent difference across populations: one parameter describes difference in both populations)
        lrt_B_im_M13=-2*(ll_Y_im)
        lrt_B_im_M23=-2*(ll_YCB_im)

        
        print "1",lrt_B_im_M12,lrt_B_im_M23,pime13-piml13,pime14-piml14


 