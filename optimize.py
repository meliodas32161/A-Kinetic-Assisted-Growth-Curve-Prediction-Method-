import geatpy as ea
from myproblem3 import MyProblem3
def param_estimation(model,NIND,I0,X0,S0,tf,z,X_exp,t):
    problem = MyProblem3(model,NIND,I0,X0,S0,tf,z,X_exp,t)
    
    algorithm = ea.soea_EGA_templet(problem,
                                ea.Population(Encoding='RI',NIND=NIND),
                                MAXGEN=100,  # 最大进化代数
                                logTras=1
                                )  
    # algorithm.mutOper.F = 0.5  # 差分进化中的参数F
    # algorithm.recOper.XOVR = 0.7  # 重组概率
    res = ea.optimize(algorithm, verbose=False, drawing=1, outputMsg=False, drawLog=False, saveFlag=False)
    
    return res