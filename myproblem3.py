import numpy as np
import geatpy as ea
from scipy.integrate import odeint
from scipy.integrate import solve_ivp

def objective(Param,model, NIND,I0,X0,S0,tf,z,X_exp,t_,lenx):
    
    # model0 = model(parameters = None)
    y_model = np.zeros(shape = [lenx,t_.shape[0]])
    error = np.zeros(t_.shape[0])
    E = np.zeros(lenx)
    # y_=np.zeros(NIND)
    for i in range(lenx): 
        model.parameter = Param.Phen[i]
        # X ,N = odeint(model.Simulation, [X0,S0], t_)
        X = solve_ivp(model.Simulation,[0,tf], [X0,S0],t_eval=t_ )
        for v,k in enumerate(t_):
            y_model[i][v] = np.array(X.y[0][v])
        for j in range(t_.shape[0]):
            error[j] = np.sqrt((X_exp[j] - y_model[i][j])**2)/(t_.shape[0])#np.abs(Y_[j] - y_model[i][j])/Y_[j]
        E[i] = np.abs(np.sum(error))
    
    return E

class MyProblem3(ea.Problem):
    def __init__(self, model, parameters,param_domain,NIND,I0,X0,S0,tf,z,X_exp,t_,ligh):
        self.model = model
        self.NIND = NIND
        self.X_exp = X_exp
        self.t_exp = t_
        self.X0 = X0
        self.S0 = S0
        self.tf = tf
        self.I0 = I0
        self.z  = z
        self.parameter = parameters
        self.domain = param_domain
        self.ligh = ligh
        name = 'Myproblem'
        M = 1
        maxormins = [1] # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）

        
        Dim = len(self.parameter)
        # parameter = {"xmax":0.5, # 环境最大浓度 0
        #              "um":0.75,  # 最大比增长率 1
        #              "Ka":0.0114, # 光衰减系数  2
        #              "KIi":20,    # 光抑制常数   3
        #              "KIs":2.0 ,   # 光饱和常数   4
        #              "ud":0.15,   # 比光合自养衰减率 5
        #              "KNi":20,    # 营养抑制常数   6
        #              "KNs":2.0,    # 营养饱和常数   7
        #              "Yx/N":200    # N元素吸收速率   8 mg/L-1
        #             }
        varTypes = [0]*Dim #[0,0,0,0,0,0]
        lb = [0, 1e-6, 1e-6,20,  1e-6,1e-6, 20,  1.0, 1e-6]
        ub = [10, 3,    1.0, 250, 100,  1.0,  800, 400,  210.0]
        lbin = [0]*Dim
        ubin = [1]*Dim
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self,name,M,maxormins,Dim,varTypes,lb,ub,lbin,ubin)
        
    # 定义目标函数
    
    def aimFunc(self,Param):
        
        lenx = len(Param.Phen)
        f = objective(Param,self.model,self.NIND,self.I0,self.X0,self.S0,self.tf,self.z,self.X_exp,self.t_exp,lenx)
        Param.ObjV = f.reshape(lenx,1) #.reshape(lenx,1)
        
        

        # CV1 = np.array([np.sum(X[:,3:7])-1 ,-np.sum(X[:,3:7])+1])
        # CV2 = np.array([np.sum(X[:,7:10])-1 ,-np.sum(X[:,7:10])+1])
        # CV = np.hstack([np.abs(np.sum(X[:,3:7])-1),np.abs(np.sum(X[:,7:10])-1)])
        # CV = np.hstack([np.abs(X4%5),np.abs(X5%2)])
        #X.CV = CV

        return Param.ObjV#,X.CV

class MyProblem3_dark(ea.Problem):
    def __init__(self, model, NIND,X0,S0,tf,z,X_exp,t_):
        self.model = model
        self.NIND = NIND
        self.X_exp = X_exp
        self.t_exp = t_
        self.X0 = X0
        self.S0 = S0
        self.tf = tf
        self.I0 = I0
        self.z  = z
        name = 'Myproblem'
        M = 1
        maxormins = [1] # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        # Dim = len(self.space.transform(x)[0])
        Dim = 9
        # parameter = {"xmax":0.5, # 环境最大浓度 0
        #              "um":0.75,  # 最大比增长率 1
        #              "Ka":0.0114, # 光衰减系数  2
        #              "KIi":20,    # 光抑制常数   3
        #              "KIs":2.0 ,   # 光饱和常数   4
        #              "ud":0.15,   # 比光合自养衰减率 5
        #              "KNi":20,    # 营养抑制常数   6
        #              "KNs":2.0,    # 营养饱和常数   7
        #              "Yx/N":200    # N元素吸收速率   8 mg/L-1
        #             }
        varTypes = [0]*Dim #[0,0,0,0,0,0]
        lb = [0, 1e-6, 1e-6,20,  1e-6,1e-6, 20,  1.0, 1e-6]
        ub = [10, 3,    1.0, 250, 100,  1.0,  800, 400,  210.0]
        lbin = [0]*Dim
        ubin = [1]*Dim
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self,name,M,maxormins,Dim,varTypes,lb,ub,lbin,ubin)
        
    # 定义目标函数
    
    def aimFunc(self,Param):
        
        lenx = len(Param.Phen)
        f = objective(Param,self.model,self.NIND,self.I0,self.X0,self.S0,self.tf,self.z,self.X_exp,self.t_exp,lenx)
        Param.ObjV = f.reshape(lenx,1) #.reshape(lenx,1)
        
        

        # CV1 = np.array([np.sum(X[:,3:7])-1 ,-np.sum(X[:,3:7])+1])
        # CV2 = np.array([np.sum(X[:,7:10])-1 ,-np.sum(X[:,7:10])+1])
        # CV = np.hstack([np.abs(np.sum(X[:,3:7])-1),np.abs(np.sum(X[:,7:10])-1)])
        # CV = np.hstack([np.abs(X4%5),np.abs(X5%2)])
        #X.CV = CV

        return Param.ObjV#,X.CV




