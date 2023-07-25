import numpy as np
from scipy.spatial import ConvexHull
from fractions import Fraction
np.set_printoptions(formatter = {'all':lambda x : str(Fraction(x).limit_denominator())})

def GE(A):
    B = np.copy(A)
    rw,cl = B.shape
    curR = 0
    curC = 0
    while(curR<rw and curC<cl):
        #選一個第curC行中，該列為非零數的列
        #為了減少接下來的誤差，我們會挑第curC行中的絕對值的最大值
        maxEc = abs(B[curR,curC])
        maxRow = curR
        for k in range(curR+1,rw):
            if abs(B[k,curC])> maxEc:
                maxEc = abs(B[k,curC])
                maxRow = k
        #把選到的那一列(第maxRow行)與第curR列互換
        B[[curR, maxRow]] = B[[maxRow, curR]]
        #print('c',B)
        #如果找不到非零(B[curR,curC]!=0),就換到下一行繼續動作
        if(B[curR,curC]!=0):
            #把第curR列從左算的第一個非零數變成1(該位置稱leading)
            for j in range(curC+1,cl):#每一項除leading
                B[curR,j] = (1/B[curR,curC])*B[curR,j]
            B[curR,curC] = 1 #leading變成1(自己除自己)
            #把curR列curC行下面的數字都變成0
            for i in range(curR+1,rw):
                #第i行 = -B[i,curC]*第curR列+第i列
                for j in range(curC+1,cl):
                    #為了減少誤差, 使用Fraction().limit_denominator()讓數字變模糊
                    B[i,j] =Fraction( B[i,j]-B[i,curC]*B[curR,j]).limit_denominator()
                B[i,curC] = 0
            #換到下一列繼續動作
            curR+=1  
        #換到下一行繼續動作
        curC+=1
    return B

def GE_solve(A,b):
    r = A.shape[0] #r為A的行數
    c = A.shape[1] #c為A的列數
    CL = np.mat(np.zeros([r,c+1])) #建立一個矩陣存放[A|b]
    CL[:,:c] = A
    CL[:,c] = b

    #print("CL",CL)
    #print(GE(CL))
    CL = GE(CL) #對[A|b]做高斯消去法
    #print("gCL",CL)
    sol =np.zeros([c,1])#準備好一個cx1矩陣存放解
    #用back-substitution求解
    for i in range(c-1,-1,-1):
        sol[i,0]=(CL[i,c])
        for j in range(c-1,i,-1):       
            sol[i,0] = sol[i,0] - CL[i,j]*sol[j,0]

    return (sol)