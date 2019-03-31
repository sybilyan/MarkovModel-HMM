#!/usr/bin/python
# -*-coding:UTF-8 -*-
# File:1.1.py
# Author:Yan Ruyu
# Data:2019/3/17
import numpy as np
import math

input_file = input("inputFile:")
inside_file = input("insideFile:")
outside_file=input("outsideFile:")

with open(inside_file) as f:
    read_data = f.read()
    a = read_data.split()
    aa = map(float, a)
    inside0 = np.array(aa)
    inside = inside0.reshape(4,4)
    #print (inside)
f.closed

with open(outside_file) as g:
    read_data = g.read()
    b = read_data.split()  # 按空格分开字符串
    bb = map(float, b)  # 将字符串化为浮点
    print bb
    outside0 = np.array(bb)  # 矩阵排列浮点数
    outside = outside0.reshape(4,4)  # 得到要求矩阵
    #print (outside)
g.closed

with open(input_file) as text:
    text_word = text.read()
    new_text_word = text_word.replace('A','0').replace('C', '1').replace('G', '2').replace('T', '3')
    s = new_text_word.split()
    #print(s)
text.closed

for x in range(len(s)):
    current_word = s[x]
    P = 0
    Q = 0
    for index in range(len(current_word)-1):
         i = current_word[index]
         j = current_word[index+1]
         i = map(int, i)
         j = map(int, j)
         in_key = inside[i,j]
         out_key = outside[i,j]
         P = float (P) + math.log(in_key,2)
         Q = float (Q) + math.log(out_key,2)
         log_ratio=P-Q
    if P > Q:
        print(" %f  inside" %(P-Q))
    else:
        print(" %f  outside" %(P-Q))

