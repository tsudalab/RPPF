import numpy as np
import csv

import pandas as pd
import linecache
import random
import itertools
from scipy import integrate

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import seaborn as sns


#############################
# parameters
#############################

# parameter in augmented weighted Tchebycheff
rho = 0.1

# number of proposals
num_top = 10

# window of weights
wind = 0.01

# value of Hstar
Hstar = 0.05

#############################

def load_data():

    name = []
    y = []

    with open('data.csv', 'r') as f:
        header = next(csv.reader(f))
        reader = csv.reader(f)
        
        for line in reader:
            name.append(line[0])

            each_y = []
            for i in range(len(line)-1):
                each_y.append(float(line[i+1]))

            y.append(each_y)
    
    return name, y

name, y = load_data()


datanum=len(name)

print("number of data =", datanum)
print("number of objectives =", len(y[0]))



if len(y[0]) == 2:

    ###################
    # two objectives
    ###################

    def Multi_func(id_can):
        return y[id_can][0], y[id_can][1]

    E1=[]
    E2=[]

    for i in range(datanum):

        E1.append(Multi_func(i)[0])
        E2.append(Multi_func(i)[1])

    # pareto solutions
    pareto_actions = []

    for i in range(datanum):

        num_strong = 0

        for j in range(datanum):

            if E1[i] > E1[j] and E2[i] > E2[j]:
                num_strong = num_strong +1

        if num_strong == 0:
            pareto_actions.append(i)


    #min-max normalization for each objective function
    E1 = np.array(E1)
    E2 = np.array(E2)

    E1_std=(E1-np.min(E1))/(np.max(E1)-np.min(E1))
    E2_std=(E2-np.min(E2))/(np.max(E2)-np.min(E2))


    #definition of alpha
    a = []

    for i in range(int(1.0/wind)+1):
        a.append(np.round(i*wind,5))

    weights = []

    for i in range(len(a)):

        weights.append([a[i],1-a[i]])

    weights = np.array(weights)



    ###################################
    ##### free energy evaluations #####
    ###################################

    FT_list=[]
    pareto_list = []
    H_all = []

    for ind_weight in range(len(weights)):

        #augmented weighted Tchebycheff
        H=[]
        for i in range(len(E1)):
            H.append(max([weights[ind_weight,0] * E1_std[i], weights[ind_weight,1] * E2_std[i]])
            + rho * (E1_std[i] + E2_std[i]))

        Emin=min(H)
        Emax=max(H)
        
        pareto_list.append(np.argmin(H))
        H_all.append((H-Emin)/(Emax-Emin))

        #print("****************")

        #print('weights=', weights[ind_weight])
        #print("")
        #print("pareto solution")
        #print("id =", np.argmin(H))


        #calculation of MIPS score

        T = - Hstar / np.log(0.5)
        
        Z = 0.0
        for i in range(datanum):
            Ene = (H[i] - Emin) / (Emax - Emin)
            Z = Z + np.exp(- Ene / T)

        if Z == 0.0:
            F = 0.0
        else:
            F = - T*(np.log(Z))    

        FT_list.append(F)

    

if len(y[0]) == 3:

    ###################
    # three objectives
    ###################

    def Multi_func(id_can):
        return y[id_can][0], y[id_can][1], y[id_can][2]


    E1=[]
    E2=[]
    E3=[]

    for i in range(datanum):

        E1.append(Multi_func(i)[0])
        E2.append(Multi_func(i)[1])
        E3.append(Multi_func(i)[2])

    # pareto solutions
    pareto_actions = []

    for i in range(datanum):

        num_strong = 0

        for j in range(datanum):

            if E1[i] > E1[j] and E2[i] > E2[j] and E3[i] > E3[j]:
                num_strong = num_strong +1

        if num_strong == 0:
            pareto_actions.append(i)

    #min-max normalization for each objective function
    E1 = np.array(E1)
    E2 = np.array(E2)
    E3 = np.array(E3)

    E1_std=(E1-np.min(E1))/(np.max(E1)-np.min(E1))
    E2_std=(E2-np.min(E2))/(np.max(E2)-np.min(E2))
    E3_std=(E3-np.min(E3))/(np.max(E3)-np.min(E3))


    #definition of alpha
    weights = []

    window_width = wind*100

    add_max = int(100/window_width)

    point = [1,0,0]

    weights.append(point)


    for i in range(add_max):

        a = add_max-i-1

        for j in range(add_max-a+1):
            b = j
            c = add_max-a-b

            point = [a*window_width/100,b*window_width/100,c*window_width/100]
                    
            weights.append(point)

    weights = np.array(weights)



    ###################################
    ##### free energy evaluations #####
    ###################################

    FT_list=[]
    pareto_list = []
    H_all = []

    for ind_weight in range(len(weights)):

        #augmented weighted Tchebycheff
        H=[]
        for i in range(len(E1)):
            H.append(max([weights[ind_weight,0] * E1_std[i], weights[ind_weight,1] * E2_std[i], weights[ind_weight,2] * E3_std[i]])
            + rho * (E1_std[i] + E2_std[i] + E3_std[i]))

        Emin=min(H)
        Emax=max(H)
        
        pareto_list.append(np.argmin(H))
        H_all.append((H-Emin)/(Emax-Emin))

        #calculation of MIPS score
        T = - Hstar / np.log(0.5)
        
        Z = 0.0
        for i in range(datanum):
            Ene = (H[i] - Emin) / (Emax - Emin)
            Z = Z + np.exp(- Ene / T)

        if Z == 0.0:
            F = 0.0
        else:
            F = - T*(np.log(Z))    

        FT_list.append(F)




#####################
##### opt value #####
#####################

arg_index = np.array(FT_list).argsort()[::-1]

ranking_index = []
ranking_MIPS = []

for i in range(len(arg_index)):

    if pareto_list[arg_index[i]] not in ranking_index:

        ranking_index.append(pareto_list[arg_index[i]])
        ranking_MIPS.append(FT_list[arg_index[i]])


print("Ranking of Pareto solutions based on MIPS score")


for i in range(len(ranking_index)):

    print(i+1, ",", name[ranking_index[i]], ",", "MIPS score =", ranking_MIPS[i], ",", "properties = ", y[ranking_index[i]])



