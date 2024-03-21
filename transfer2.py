import numpy as np
from gekko import  GEKKO

def Transfer_ea(m):
    Param = m.m._parameters.copy()
    # Const = m._constants.copy()
    if isinstance(Param[0].VALUE.value,float ):
        xmax = [i.VALUE.value for i in Param if i.NAME == 'xmax'][0]      
        mum =  [i.VALUE.value for i in Param if i.NAME == 'mum'][0]
        ki  = [i.VALUE.value for i in Param if i.NAME == 'ki'][0] # 光抑制常数 
        ks  = [i.VALUE.value for i in Param if i.NAME == 'ks'][0] # 光饱和常数
        
        ud = [i.VALUE.value for i in Param if i.NAME == 'ud'][0]
        kNi = [i.VALUE.value for i in Param if i.NAME == 'kni'][0]# 营养抑制常数
        kNs = [i.VALUE.value for i in Param if i.NAME == 'kns'][0]# 营养饱和常数
       
        Yxn = [i.VALUE.value for i in Param if i.NAME == 'yxn'][0]
        # kd = [i.VALUE.value for i in Param if i.NAME == 'kd'][0]
        
        # kr = [i.VALUE.value for i in Param if i.NAME == 'kr'][0]
        
        # tau = [i.VALUE.value for i in Param if i.NAME == 'tau'][0]
        
        # delt = [i.VALUE.value for i in Param if i.NAME == 'delt'][0]
        
        # k = [i.VALUE.value for i in Param if i.NAME == 'k'][0]
        
        # R = [i.VALUE.value for i in Param if i.NAME == 'R'][0]
        
        
        
        # delta = [i.VALUE.value for i in Param if i.NAME == 'delta'][0]
        
        # Yxs = [i.VALUE.value for i in Param if i.NAME == 'Yxs'][0]
       
        
        parameter = [xmax,mum,ki,ks,ud,kNi,kNs,Yxn]
    
    if isinstance(Param[0].VALUE.value,np.ndarray ):
        xmax = [i.VALUE.value for i in Param if i.NAME == 'xmax'][0][0]      
        mum =  [i.VALUE.value for i in Param if i.NAME == 'mum'][0][0]
        ki  = [i.VALUE.value for i in Param if i.NAME == 'ki'][0][0] # 光抑制常数 
        ks  = [i.VALUE.value for i in Param if i.NAME == 'ks'][0][0] # 光饱和常数
        
        ud = [i.VALUE.value for i in Param if i.NAME == 'ud'][0][0]
        kNi = [i.VALUE.value for i in Param if i.NAME == 'kni'][0][0]# 营养抑制常数
        kNs = [i.VALUE.value for i in Param if i.NAME == 'kns'][0][0]# 营养饱和常数
       
        Yxn = [i.VALUE.value for i in Param if i.NAME == 'yxn'][0][0]
        
        # kd = [i.VALUE.value for i in Param if i.NAME == 'kd'][0][0]
        
        # kr = [i.VALUE.value for i in Param if i.NAME == 'kr'][0][0]
        
        # tau = [i.VALUE.value for i in Param if i.NAME == 'tau'][0][0]
        
        # delt = [i.VALUE.value for i in Param if i.NAME == 'delt'][0][0]
        
        # k = [i.VALUE.value for i in Param if i.NAME == 'k'][0][0]
        
        # R = [i.VALUE.value for i in Param if i.NAME == 'R'][0][0]
        
        
        # xmax = [i.VALUE.value for i in Param if i.NAME == 'xmax'][0][0]
        
        # delta = [i.VALUE.value for i in Param if i.NAME == 'delta'][0][0]
        
        # Yxs = [i.VALUE.value for i in Param if i.NAME == 'Yxs'][0][0]
       
        
        parameter = [xmax,mum,ki,ks,ud,kNi,kNs,Yxn]
      
    
    return np.array(parameter)

# def inv_Transfer_ea(m,Param,NIND):
    
#     parameter = np.zeros(shape = [NIND,6])
#     # ks  = np.zeros_like(mum)
#     # ki  = np.zeros_like(mum)
#     # I0  = np.zeros_like(mum)
#     # Ka  = np.zeros_like(mum)
#     # xmax  = np.zeros_like(mum)
    
#     for i in range(NIND):
#         for j in range(6):
#             parameter[i][j] = m.Param(value = Param.Phen[[i][:]])

def inv_Transfer_ea(m,Param,i,light):
    
    # parameter = np.zeros(6)
    m._parameters.clear()
    m._equations.clear()
    m._intermediates.clear()
    m._inter_equations.clear()
    m._variables.clear()
    m._constants.clear()
    # m.clear()
    if light == 1:
        xmax = m.Param(value = Param.Phen[i][0], name = 'xmax')
        mum  = m.Param(value = Param.Phen[i][1], name = 'mum')
        ki = m.Param(value = Param.Phen[i][2], name = 'ki')
        ks = m.Param(value = Param.Phen[i][3], name = 'ks')
        ud = m.Param(value = Param.Phen[i][4], name = 'ud')
        kNi = m.Param(value = Param.Phen[i][5], name = 'kni')
        kNs = m.Param(value = Param.Phen[i][6], name = 'kns')
        Yxn = m.Param(value = Param.Phen[i][7], name = 'yxn')
    #     delta = m.Param(value = Param.Phen[i][8], name = 'delta')
    #     Yxs = m.Param(value = Param.Phen[i][9], name = 'Yxs')
    #     ks = m.Param(value = Param.Phen[i][10], name = 'ks')
    #     ki = m.Param(value = Param.Phen[i][11], name = 'ki')
    if light == 1:
        xmax = m.Param(value = Param.Phen[i][0], name = 'xmax')
        mum  = m.Param(value = Param.Phen[i][1], name = 'mum')
        ud = m.Param(value = Param.Phen[i][2], name = 'ud')
        kNi = m.Param(value = Param.Phen[i][3], name = 'kni')
        kNs = m.Param(value = Param.Phen[i][4], name = 'kns')
        Yxn = m.Param(value = Param.Phen[i][5], name = 'yxn')

    return m

def inv_Transfer_gekko(m,light=1):
    # 把gekko Fv转化为字典
    Param = m._parameters

    
    if light == 1:
        # parameter = {'mum':Param[0],
        #               'kd':Param[1],
        #               'kr':Param[2],
        #               'tau':Param[3],
        #               'delt':Param[4],
        #               'k':Param[5],
        #               'R':Param[6],
        #               'xmax':Param[7],
        #               'delta':Param[8],
        #               'Yxs':Param[9]
        #             }
        parameter = {'mum':Param[0],
                      'kd':Param[1],  # 损伤速率 
                      'kr':Param[2],   # 光和系统修复速率 s^-1
                      'tau':Param[3],     # 周转率 s
                      'delt':Param[4],    # m^2.(µmol)^−1
                      'k':Param[5],     # 连接接受能量和生长速率的关系
                      'me':Param[6],   # 呼吸速率
                      'xmax':Param[7],
                      'delta':Param[8],
                      'Yxs':Param[9],
                      'ks':Param[10],
                      'ki':Param[11]
                      }
    if light == -1:
        parameter = {'mum':Param[0],
                      
                      'xmax':Param[7],
                      'delta':Param[8],
                      'Yxs':Param[9]
                    }
    return parameter
    
