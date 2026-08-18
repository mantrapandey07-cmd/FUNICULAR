import numpy as np
from qpsolvers import solve_qp
import math
w,l=10,10
init_x,init_y=w/2,l/2
x,y=w/2,l/2
target_x,target_y=0,0
W= np.array([-x/(math.sqrt(x**2+(l-y)**2)),(w-x)/(math.sqrt((w-x)**2+(l-y)**2)),(w-x)/(math.sqrt((w-x)**2+y**2)),-x/(math.sqrt(x**2+y**2))]
            ,[(l-y)/(math.sqrt(x**2+(l-y)**2)),(l-y)/(math.sqrt((w-x)**2+(l-y)**2)),-y/(math.sqrt((w-x)**2+y**2)),-y/(math.sqrt(x**2+y**2))],
            [0,0,0,0])
dist=math.sqrt((target_x-x)^2+(target_x-y)^2)
f=((target_x-x)/(dist,target_y-y)/dist,0) if (() and ())
lb = np.array([1.0, 1.0, 1.0, 1.0])
ub = np.array([100.0, 100.0, 100.0, 100.0])
tensions=solve_qp(np.eye(4),np.zeros(4),None,None,W,f,lb,ub,solver="osqp")

