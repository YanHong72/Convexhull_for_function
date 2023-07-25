import numpy as np
import numpy.random as rn
from scipy.spatial import ConvexHull
from fractions import Fraction


def hullfordim1(A):
    Max = Fraction(A[0,0])#最大值
    indexOfMax = 0 #最大值的索引值
    Min = Fraction(A[0,0])#最小值
    indexOfMin = 0 #最小值的索引值
    #尋找最大值和最小值的索引值
    for i in range(A.shape[0]):
        if Max < Fraction(A[i,0]):
            Max = Fraction(A[i,0])
            indexOfMax = i
        if Min > Fraction(A[i,0]):
            Min = Fraction(A[i,0])
            indexOfMin = i
    #回傳索引值
    return [indexOfMin,indexOfMax,]
def UseVerticesToCheck(A): 
    n,m = A.shape  
    if m == 1: #特殊情況:如果m=1,只要找最大最小值的位置即可
        vt = hullfordim1(A)
    else: #找哪一列是頂點
        vt = ConvexHull(A).vertices
    #print(vt)
    if n-1 in vt:#確定n是否為頂點
        print("提示: p為頂點")
        return False
    else:
        print("提示: p不為頂點")
        return True