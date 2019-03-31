#!/usr/bin/python
# -*-coding:UTF-8 -*-
# File:1.2.py
# Author:Yan Ruyu
# Data:2019/3/18

# HMM.py
# Using Vertibi algorithm

import numpy as np
import math

def main():
    #read 观测集合
    input_file = input("inputFile:")

    #读取观测序列
    with open(input_file) as f:
        a=f.readlines()
        aa = ''.join(a[1:])   #把第一行略去
        myfile=aa.upper()
        #print len(myfile)
        seq=myfile.replace('A','0').replace('C', '1').replace('G', '2').replace('T', '3').replace(' ','').replace('\n','').replace('\r','')
        #print type(seq)
        #print len(seq)
        #print seq

        ls= list(map(int, seq))

        #print ls
        #new_ls = filter(None, ls)  #去除空格
        #print len(new_ls)
    f.close()


    global N,M

#读取概率信息文件
    input_DNA = input("inputDNA:")
    with open(input_DNA) as g:
        HMM_file=g.read().split()
        #print HMM_file
        N = int( HMM_file[0])    #用map函数没用 map（int，HMM_file）
        M =int( HMM_file[1])

        B = [[0 for i in range(M)] for j in range(N)]      #2*4 N*M  [[0] * M] * N 不行
        A = [[0 for i in range(N)] for j in range(N)]

# 状态集合（自己改）
        Q = ('state1', 'state2')

        # 观测集合V
        VV= HMM_file[2].replace('A', '0').replace('C', '1').replace('G', '2').replace('T', '3')
        V= list(VV)       #string to List


        # 转移概率: A  发射概率:B
        x =3 + N
        for i in range(0,N):                     #range:2
            for j in range(N+M):                 #range:6
                if(j<N):
                    A[i][j]=float(HMM_file[x])
                    x+=1
                    #print A
                    #print A[i][j]
                else:
                    B[i][j-N]=float(HMM_file[x])
                    x+=1
        #print A,B

        # 初始概率 N=2
        PI = [[0]*1]*N
        for i in range(N-1):
            PI[0][i]=float( HMM_file[i+3])
        #print PI
    g.close()

    P, hidden_states = Viterbi(A,B,PI,V,Q,ls)


#输出state2段数
    out_index = 1
    list_hidden = list(hidden_states)
    state2_num=[]
    for i in range(1,len(list_hidden)):
        if list_hidden[i] is not list_hidden[i-1]:
            print('%d  %d  %s' %(out_index ,i-1 ,list_hidden[i-1]))
            state2_num.append(list_hidden[i-1])
            out_index = i
        else :
            if i==len(list_hidden)-1:
                print('%d  %d  %s'  %(out_index, i+1, list_hidden[i]))
                state2_num.append(list_hidden[i])
            else:
                continue


    key=0
    for i in range(len(state2_num)):
        if state2_num[i] is Q[1]:
            key+=1

    print ('num of %s is %d'  %(Q[1],key))


#viterbi函数
def Viterbi(A, B, PI, V, Q, obs):                   #A为转移概率 B发射概率  PI 初始概率 V观测集合 Q状态集合  obs观测序列

    N = len(Q)                                      # Q有几个字符组 Q=2
    T = len(obs)                                    #obs有几个字符组
    #print T
    delta =  np.array([[0] * N] * T, dtype=np.float64)   #delta T*2
    phi = np.array([[0] * N] * T, dtype=np.int64)       #phi  T*2
                                                   # 初始化
    for i in range(0,N):
        delta[0][i] = math.log(PI[i][0],2) +math.log( B[i][obs[0]],2)  #检索V中有没有obs[0],并返回它的位置
        phi[0][i] = 0
    #print delta[0][0],delta[0][1]             #array*column 行乘列

# 主体递归计算
    for i in range(1, T):
        for j in range(N):              #N=2
            tmp = [delta[i-1, k]+math.log(A[k][j],2) for k in range(N)]    #1*2  Sj,L
            delta[i][j] = max(tmp)+ math.log(B[j][obs[i]],2)   #当前最终值
            #print delta
            phi[i][j] = tmp.index(max(tmp))
            #print phi
    # 最终的概率及节点
    P = max(delta[T-1, :])
    I = int(np.argmax(delta[T-1, :]))  #np.argmax表示输出最大数的位置


# 最优路径path  backtrace回退过程
    path = [I]
    #print path
    for i in reversed(range(1, T)):
        end = path[-1]
        path.append(phi[i][end])
        #print len(path)

    hidden_states = [Q[i] for i in reversed(path)] #reversed 表示列表数据的反转

    return P, hidden_states


main()