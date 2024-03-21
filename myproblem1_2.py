import numpy as np
from numpy import *
import geatpy as ea
# from gaussian_ei import gaussian_ei as EI
import warnings
from gekko import GEKKO
from transfer2 import inv_Transfer_ea



from sklearn.metrics import mean_squared_error

def objective(m,Param, model, NIND,I0,X0,S0,z,light,Y_,t_):
    
    y_model = np.zeros(shape = [NIND,t_.shape[0]])
    error = np.zeros(t_.shape[0])
    E = np.zeros(NIND)
    t0 = t_.tolist()
    # y_=np.zeros(NIND)
    for i in range(NIND):
        m = inv_Transfer_ea(m, Param, i,light = light)
        model0 = model(m=m,I0=I0,X0=X0,S0=S0,z=z,light = light)
        yX,yS = model0.simulation(t_ = t_)
        for v,k in enumerate(t_):
            y_model[i][v] = np.array(yX[v])
        for j in range(t_.shape[0]):
            error[j] = np.sqrt((Y_[j] - y_model[i][j])**2)/(t_.shape[0])#np.abs(Y_[j] - y_model[i][j])/Y_[j]
        E[i] = np.abs(np.sum(error))
    
    return E



class MyProblem1(ea.Problem):
    def __init__(self,m,model,NIND,Y_,t_,X0,I0,S0,z,light):

        self.m = m
        self.model = model
        self.NIND = NIND
        self.X_exp = Y_
        self.t_exp = t_
        # self.delta = delta
        self.X0 = X0
        self.S0 = S0
        self.I0 = I0
        self.z  = z
        self.light = light
        
        name = 'Myproblem'
        M = 1
        maxormins = [1] # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        # Dim = len(self.space.transform(x)[0])
        Dim = 8
        
        
        #parameter = [mum,kd,kr,tau,delt,k,me,xmax,delta,Yxs,ks,ki]  
        # parameter = {'mum':0.45354,
        #               'kd':2.99e-04,  # 损伤速率 
        #               'kr':6.8e-03,   # 光和系统修复速率 s^-1
        #               'tau':0.25,     # 周转率 s
        #               'delt':0.047,    # m^2.(µmol)^−1
        #               'k':8.7e-6,     # 连接接受能量和生长速率的关系
        #               'me':1.389e-7,   # 呼吸速率
        #               'xmax':33.9261,
        #               'delta':0.1,
        #               'Yxs':0.1,
        #               'ks':0.1,
        #               'ki':0.1
        #               }
                 
        #          }
        varTypes = [0]*Dim #[0,0,0,0,0,0]
        # lb = [-10,1e-6, 1e-6,1e-2,1e-4,1e-3,1e-3, 0.1, 1e-5,1e-8,1e-8,1e-8]
        # ub = [100, 0.1, 0.1, 100, 100, 100,  1, 100, 1.0,100.0,100,100]
        lb = [1e-6, 1e-6, 20,  1e-6,1e-6, 20,  1e-6, 1e-6]
        ub = [1, 1,   100, 10,  1.0,  100, 10,   1.0]
        lbin = [0]*Dim
        ubin = [1]*Dim
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self,name,M,maxormins,Dim,varTypes,lb,ub,lbin,ubin)
        
    # 定义目标函数
    
    def aimFunc(self,Param):
        
        # mum = Param.Phen[:,[0]]
        # ks = Param.Phen[:,[1]]
        # ki = Param.Phen[:,[2]]
        # I0 = Param.Phen[:,[3]]
        # Ka = Param.Phen[:,[4]]
        # xmax = Param.Phen[:,[5]]

               
        lenx = len(Param.Phen)
        f = objective(self.m,Param, self.model, lenx,self.I0,self.X0,self.S0,self.z,self.light,self.X_exp,self.t_exp)
        Param.ObjV = f.reshape(lenx,1)#.reshape(lenx,1)
        
        

        # CV1 = np.array([np.sum(X[:,3:7])-1 ,-np.sum(X[:,3:7])+1])
        # CV2 = np.array([np.sum(X[:,7:10])-1 ,-np.sum(X[:,7:10])+1])
        # CV = np.hstack([np.abs(np.sum(X[:,3:7])-1),np.abs(np.sum(X[:,7:10])-1)])
        # CV = np.hstack([np.abs(X4%5),np.abs(X5%2)])
        #X.CV = CV

        return Param.ObjV#,X.CV