import numpy as np
from gekko import GEKKO


def Initial(m,parameter,light=1):
    m._parameters.clear()
    m._equations.clear()
    m._intermediates.clear()
    m._inter_equations.clear()
    m._variables.clear()
    m._constants.clear()
    # m.clear()
    if light == 1:
        xmax = m.Param(value = parameter['xmax'], name = 'xmax')
        mum = m.Param(value = parameter['mum'], name = 'mum')
        ks = m.Param(value = parameter['ks'], name = 'ks')
        ki = m.Param(value = parameter['ki'], name = 'ki')
        ud = m.Param(value = parameter['ud'], name = 'ud')
        
        # kr = m.Param(value = parameter['kr'], name ='kr')
        
        # tau = m.Param(value = parameter['tau'], name ='tau')
        
        # delt = m.Param(value = parameter['delt'], name ='delt')
        
        # k = m.Param(value = parameter['k'], name ='k')
        
        # me = m.Param(value = parameter['me'], name ='me')
        
        
        
        # delta = m.Param(value = parameter['delta'], name = 'delta')
        
        # alpha = m.Param(value = parameter['Yxs'], name = 'Yxs')
        
        
        kNi = m.Param(value = parameter['kNi'], name = 'kni')
        kNs = m.Param(value = parameter['kNs'], name = 'kns')
        Yxs = m.Param(value = parameter['Yx/N'], name = 'yxn')
       
        # tau = self.m.Param(value = self.parameter['tau'], name ='tau')
       
        # L = m.Const(value = parameter['L'], name ='L')
        
        # z = m.Const(value = parameter['z'], name ='z')
    if light == -1:
        xmax = m.Param(value = parameter['xmax'], name = 'xmax')
        mum = m.Param(value = parameter['mum'], name = 'mum')
        
        delta = m.Param(value = parameter['delta'], name = 'delta')
        
        alpha = m.Param(value = parameter['Yxs'], name = 'Yxs')
       

    return m#, Param

def Initial_gekko(m,parameter=None,light=1):
    m._parameters.clear()
    m._equations.clear()
    m._intermediates.clear()
    m._inter_equations.clear()
    m._variables.clear()
    m._constants.clear()
    
    if light == 1:
        # mum = m.FV(value = parameter['mum'], name = 'mum')
        
        # kd = m.FV(value = parameter['kd'], name = 'kd')
        
        # kr = m.FV(value = parameter['kr'], name ='kr')
        
        # tau = m.FV(value = parameter['tau'], name ='tau')
        
        # delt = m.FV(value = parameter['delt'], name ='delt')
        
        # k = m.FV(value = parameter['k'], name ='k')
        
        # R = m.FV(value = parameter['R'], name ='R')
        
        # xmax = m.FV(value = parameter['xmax'], name = 'xmax')
        
        # delta = m.FV(value = parameter['delta'], name = 'delta')
        
        # alpha = m.FV(value = parameter['alpha'], name = 'alpha')
        
       
        mum = m.FV(value = parameter['mum'],lb=-10,ub=100, name = 'mum');mum.STATUS = 1;mum.FSTATUS = 0
        
        kd = m.FV(value = parameter['kd'],lb=1e-6,ub=1e-1,name ='kd');kd.STATUS = 1;kd.FSTATUS = 0
        
        kr = m.FV(value = parameter['kr'],lb=1e-6,ub=1e-1,name ='kr');kr.STATUS = 1;kr.FSTATUS = 0
        
        tau = m.FV(value = parameter['tau'],lb=1e-2,ub=100,name ='tau');tau.STATUS = 1;tau.FSTATUS = 0
        
        delt = m.FV(value = parameter['delt'],lb=1e-4,ub=100,name ='delt');delt.STATUS = 1;delt.FSTATUS = 0
        
        k = m.FV(value = parameter['k'],lb=1e-3,ub=100,name ='k');k.STATUS = 1;k.FSTATUS = 0
        
        me = m.FV(value = parameter['me'],lb=1e-3,ub=1,name ='me');me.STATUS = 1;me.FSTATUS = 0
        
        xmax = m.FV(value = parameter['xmax'],lb=0,name = 'xmax');xmax.STATUS = 1;xmax.FSTATUS = 0
        
        delta = m.FV(value = parameter['delta'],name = 'delta');delta.STATUS = 1;delta.FSTATUS = 0
        
        alpha = m.FV(value = parameter['Yxs'],lb=1e-8,name = 'Yxs');alpha.STATUS = 1;alpha.FSTATUS = 0
        ks = m.FV(value = parameter['ks'],lb=1e-8,name = 'ks');alpha.STATUS = 1;alpha.FSTATUS = 0
        ki = m.FV(value = parameter['ki'],lb=1e-8,name = 'ki');alpha.STATUS = 1;alpha.FSTATUS = 0
    if light == -1:
        mum = m.FV(name = 'mum')
        
        xmax = m.FV( name = 'xmax')
        
        delta = m.FV(name = 'delta')
        
        alpha = m.FV(name = 'Yxs')
    # if parameter != None:
    #     if light == 1:
    #         mum = m.FV(name = 'mum')
            
    #         kd = m.FV(name = 'kd')
            
    #         kr = m.FV(name ='kr')
            
    #         tau = m.FV(name ='tau')
            
    #         delt = m.FV( name ='delt')
            
    #         k = m.FV(name ='k')
            
    #         R = m.FV(name ='R')
            
    #         xmax = m.FV(name = 'xmax')
            
    #         delta = m.FV(name = 'delta')
            
    #         alpha = m.FV(name = 'alpha')
           
    #         # tau = self.m.Param(value = self.parameter['tau'], name ='tau')
           
    #         # L = m.Const(value = parameter['L'], name ='L')
            
    #         # z = m.Const(value = parameter['z'], name ='z')
    #     if light == -1:
    #         mum = m.FV(name = 'mum')
            
    #         xmax = m.FV( name = 'xmax')
            
    #         delta = m.FV(name = 'delta')
            
    #         alpha = m.FV(name = 'alpha')
       

    return m

