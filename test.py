import matplotlib.pyplot as plt
import math
import numpy as np

# ro = p|r><r| + (1-p)|r_(perpendicular)><r_(perpendicular)|
# Purity = Tr[ro^2] = 2p^2 - 2p + 1

def p_entropy_func(p):
    if p == 0 or p == 1:
        return 0
    return -p * np.log2(p) - (1-p) * np.log2(1-p)

def p_search(H_ro):
    error = 0.001
    p_list = [i/100000 for i in range(0, 50000)]
    init = 0
    fin = len(p_list) - 1
    while fin >= init:
        mid = init + (fin - init) // 2
        p = p_list[mid]
        entropy_diff = p_entropy_func(p) - H_ro
        if abs(entropy_diff) < error:
            return p
        elif entropy_diff > 0:
            fin = mid - 1
        else:
            init = mid + 1
    return -1

def purity_calculator(H_ro):
    if H_ro == 0:
        return 1, 0.999
    else:
        p = p_search(H_ro)
        return 2*(p**2) - 2*p +1, p

def p_calculator(H_ro):
    if H_ro == 0:
        return 0.999
    elif H_ro == 0.3:
        p = 0.053239
        #p = 0.946761
    elif H_ro == 0.6:
        #p = 0.853898
        p = 0.146102
    elif H_ro == 0.9:
        p = 0.316019
        # p = 0.683981
    return p

def lower_bound_eq14(P, H_ro, c):
    binary_entropy = []
    for elem in c:
        value = (math.sqrt(2*P-1)*(2*math.sqrt(elem)-1)+1)/2
        binary_entropy.append( -value*math.log2(value)-(1-value)*math.log2(1-value)-H_ro)
    return binary_entropy

def lower_bound_eq4(H_ro, c):
    lower_bound = []
    for elem in c:
        lower_bound.append(-(1-H_ro)*math.log2(elem))
    return lower_bound

def lower_bound_eq12(H_ro, c):
    lower_bound = []
    for elem in c:
        value = -math.log2(elem) - H_ro
        if value > 0:
            lower_bound.append(value)
        else:
            lower_bound.append(0)
    return lower_bound

def lower_bound_eq13(H_ro, c):
    binary_entropy = []
    for elem in c:
        value = (1 + math.sqrt(2 * elem - 1))/2
        bn_entr = -value*math.log2(value)-(1-value)*math.log2(1-value) - 2*H_ro
        if bn_entr > 0:
            binary_entropy.append( bn_entr)
        else:
            binary_entropy.append(0)
    return binary_entropy

def entropy_calculator(x, p):
    value = x*p + (1-x)*(1-p)
    entropy = - value * math.log2(value) - (1-value) * math.log2(1-value)
    entropy_p = -p * math.log2(p) - (1-p) * math.log2(1-p)
    return entropy - entropy_p


def numerical_optimum(c, p):
    num_optim = []
    for elem in c:
        gamma = math.acos(2*elem - 1)
        fmin = 5
        alfa_list = [i/100 for i in range(int((gamma/2)*100), int(gamma*100)+1)]
        for alfa in alfa_list:
            x1 = (math.cos(alfa) + 1)/2
            x2 = (math.cos(gamma - alfa)+1)/2
            entropy_1 = entropy_calculator(x1, p)
            entropy_2 = entropy_calculator(x2, p)
            f = entropy_1 + entropy_2
            if f < fmin:
                fmin = f 
        num_optim.append(fmin)   
    return num_optim 


def l1_entropy(H_ro, P, c):
    l1 = []
    for elem in c:
        value = 2 * math.sqrt(2 * P - 1) * math.sqrt(elem * (1 - elem))
        l1.append(value)
    return l1

def mod(data):
    if data >= 0 :
        return data 
    return (-1*data)

fig_1, axs_1 = plt.subplots(2, 2, figsize=(10, 8))  # 2 rows, 2 columns
axs_1 = axs_1.flatten()

#c = [i/1000 for i in range(500, 1000)]
c = [i/1000 for i in range(500, 1000)]
h_ro = [0, 0.3, 0.6, 0.9]
ip = 0
for H_ro in h_ro:
    
    P, p = purity_calculator(H_ro)
    lb4 = lower_bound_eq4(H_ro, c)
    lb12 = lower_bound_eq12(H_ro, c)
    lb13 = lower_bound_eq13(H_ro, c)
    lb14 = lower_bound_eq14(P, H_ro, c)
    no = numerical_optimum(c, p)

    axs_1[ip]
    axs_1[ip].plot(c, lb4, label = "Eq. (4)")
    axs_1[ip].plot(c, lb12, label = "Eq. (12)")
    axs_1[ip].plot(c, lb13, label = "Eq. (13)")
    axs_1[ip].plot(c, lb14, label = "Eq. (14)")
    axs_1[ip].plot(c, no, label = "Numerical optimum")

    axs_1[ip].set_xlabel('c')
    axs_1[ip].set_ylabel('Lower bound')
    axs_1[ip].set_title('H(ro)=' + str(H_ro))
    axs_1[ip].legend()
    ip = ip + 1

fig_1.savefig('plots/H(rho)'  + '.png', bbox_inches='tight')


lb4 =[]
lb12 = []
lb13 = []
lb14 = []
no = []

num_rows = 500
num_cols = 500

hros = [i/500 for i in range(0, num_cols)]
i =0 

for H_ro in hros:

    P, p = purity_calculator(H_ro)

    lb4 = lb4 + (lower_bound_eq4(H_ro, c))
    lb12 = lb12 + (lower_bound_eq12(H_ro, c))
    lb13 = lb13 + (lower_bound_eq13(H_ro, c))
    lb14 = lb14 + (lower_bound_eq14(P, H_ro, c))
    no = no + (numerical_optimum(c, p))
    i = i + 1
    print(f'\r Calculating purity {i//5} % \n', end=' ', flush=True) 

print("Done")

difflb4_12 =[ mod(x - y) for x, y in zip(lb4, lb12)]
difflb4_13 = [mod(x - y) for x, y in zip(lb4, lb13)]
difflb4_14 = ([mod(x - y) for x, y in zip(lb4, lb14)])
difflb12_13 = ([mod(x - y) for x, y in zip(lb12, lb13)])
difflb12_14 = [mod(x - y) for x, y in zip(lb12, lb14)]
difflb13_14 = [mod(x - y) for x, y in zip(lb13, lb14)]

difflbno_4 = [mod(x - y) for x, y in zip(no, lb4)]
difflbno_12 = [mod(x - y) for x, y in zip(no, lb12)]
difflbno_13 = [mod(x - y) for x, y in zip(no, lb13)]
difflbno_14 = [mod(x - y) for x, y in zip(no, lb14)]

list_vals_plot = [(lb4,"eq4"), (lb12,"eq12"), (lb13,"eq13"), (lb14,"eq14"), (difflb4_12,"diff eq4 - eq12"), (difflb4_13,"diff eq4 - eq13"), (difflb4_14,"diff eq4 - eq14"), (difflb12_13,"diff eq12 -eq13"), (difflb12_14,"diff eq12 - eq14"), (difflb13_14,"diff eq13 -eq14")]
list_no_plot = [(no,"Numerical optimum"), (difflbno_4,"diff optimum - eq4"), (difflbno_12,"diff optimum - 12"), (difflbno_13,"diff optimum - eq13"), (difflbno_14,"diff optimum - eq14")]

# For the equation heatmaps
first_four_data_arrays = list_vals_plot[:4]
all_values_first_four = [value for array, _ in first_four_data_arrays for value in array]

overall_min_first_four = min(all_values_first_four)
overall_max_first_four = max(all_values_first_four)


# For the diff heatmaps
last_six_data_arrays = list_vals_plot[-6:]
all_values_last_six = [value for array, _ in last_six_data_arrays for value in array]

overall_min_last_six = min(all_values_last_six)
overall_max_last_six = max(all_values_last_six)


# For the diff heatmaps

all_values_no_plots = [value for array, _ in list_no_plot for value in array]

overall_min_no_plots = min(all_values_no_plots)
overall_max_no_plots = max(all_values_no_plots)

# For the diff heatmaps
data_no_plots414 = [(difflbno_4,"diff optimum - eq4"),(difflbno_14,"diff optimum - eq14")]
all_values_no_plots414 = [value for array, _ in data_no_plots414 for value in array]


overall_min_no_plots414 = min(all_values_no_plots414)
overall_max_no_plots414 = max(all_values_no_plots414)




fig_2, axs_2 = plt.subplots(2, 5, figsize=(8*5, 6*2))

axs_2 = axs_2.flatten()

for i, (data,name) in enumerate(list_vals_plot):
    
    #res = [ele for ele in data for i in range(num_cols)]
    matrix = np.array(data).reshape(num_rows, num_cols)
    if i >= 5 :
        im = axs_2[i].imshow(matrix, cmap='gray', interpolation='nearest',  vmin=overall_min_last_six, vmax=overall_max_last_six)
    else:
        im = axs_2[i].imshow(matrix, cmap='gray', interpolation='nearest',  vmin=overall_min_first_four, vmax=overall_max_first_four)
    fig_2.colorbar(im, ax=axs_2[i])
    axs_2[i].set_xlabel('H(rho)')
    axs_2[i].set_ylabel('c')
    axs_2[i].set_title(f"Heatmap {name}")

    x_ticks = np.linspace(0, num_cols , num=10)
    x_ticks_normalized = (x_ticks - x_ticks.min()) / (x_ticks.max() - x_ticks.min())
    
    axs_2[i].set_xticks(x_ticks)
    axs_2[i].set_xticklabels([f'{xt:.2f}' for xt in x_ticks_normalized])

    # Set the desired number of y-ticks
    y_ticks = np.linspace(0, num_rows - 1, num=5)
    
    # Normalize the y-tick labels to be between 0 and 1
    y_ticks_normalized = (y_ticks - y_ticks.min()) / (y_ticks.max() - y_ticks.min())
    
    # Set the ticks and the normalized labels
    axs_2[i].set_yticks(y_ticks)
    axs_2[i].set_yticklabels([f'{0.5 + yt/2:.2f}' for yt in y_ticks_normalized])

fig_2.suptitle('Heatmaps for eq4, eq12, eq13, eq14 and the difference between eachother, with H(rho) = 0.6 ', fontsize=16)

fig_2.savefig('plots/Heatmaps'  + '.png')





fig_3, axs_3 = plt.subplots(1, 5, figsize=(8*5, 6))

axs_3 = axs_3.flatten()

for i, (data,name) in enumerate(list_no_plot):
    
    #res = [ele for ele in data for i in range(num_cols)]
    matrix = np.array(data).reshape(num_rows, num_cols)
    
    im = axs_3[i].imshow(matrix, cmap='gray', interpolation='nearest',  vmin=overall_min_no_plots, vmax=overall_max_no_plots)
   
    fig_3.colorbar(im, ax=axs_3[i])
    axs_3[i].set_xlabel('H(rho)')
    axs_3[i].set_ylabel('c')
    axs_3[i].set_title(f"Heatmap {name}")

    x_ticks = np.linspace(0, num_cols  , num=10)
    x_ticks_normalized = (x_ticks - x_ticks.min()) / (x_ticks.max() - x_ticks.min())
    
    axs_3[i].set_xticks(x_ticks)
    axs_3[i].set_xticklabels([f'{xt:.2f}' for xt in x_ticks_normalized])

    # Set the desired number of y-ticks
    y_ticks = np.linspace(0, num_rows - 1, num=5)
    
    # Normalize the y-tick labels to be between 0 and 1
    y_ticks_normalized = (y_ticks - y_ticks.min()) / (y_ticks.max() - y_ticks.min())
    
    # Set the ticks and the normalized labels
    axs_3[i].set_yticks(y_ticks)
    axs_3[i].set_yticklabels([f'{0.5 + yt/2:.2f}' for yt in y_ticks_normalized])

fig_3.suptitle('Heatmaps for numerical optimum and the difference between numerical optimum and eq4, e12, eq13 and eq14, with H(rho) = 0.6 ', fontsize=16)

fig_3.savefig('plots/HeatmapsNumericalOptimum'  + '.png')



fig_4, axs_4 = plt.subplots(1, 2, figsize=(8*2, 6))

axs_4 = axs_4.flatten()

for i, (data,name) in enumerate(data_no_plots414):
    
    #res = [ele for ele in data for i in range(num_cols)]
    matrix = np.array(data).reshape(num_rows, num_cols)
    
    im = axs_4[i].imshow(matrix, cmap='gray', interpolation='nearest',  vmin=overall_min_no_plots414, vmax=overall_max_no_plots414)
   
    fig_4.colorbar(im, ax=axs_4[i])
    axs_4[i].set_xlabel('H(rho)')
    axs_4[i].set_ylabel('c')

    axs_4[i].set_title(f"Heatmap {name}")

    x_ticks = np.linspace(0, num_cols , num=10)
    x_ticks_normalized = (x_ticks - x_ticks.min()) / (x_ticks.max() - x_ticks.min())

    axs_4[i].set_xticks(x_ticks)
    axs_4[i].set_xticklabels([f'{xt:.2f}' for xt in x_ticks_normalized])

    

    # Set the desired number of y-ticks
    y_ticks = np.linspace(0, num_rows - 1, num=5)
    
    # Normalize the y-tick labels to be between 0 and 1
    y_ticks_normalized = (y_ticks - y_ticks.min()) / (y_ticks.max() - y_ticks.min())
    
    # Set the ticks and the normalized labels
    axs_4[i].set_yticks(y_ticks)
    axs_4[i].set_yticklabels([f'{0.5 + yt/2:.2f}' for yt in y_ticks_normalized])

fig_4.suptitle('Heatmaps for comparisons only between best 2 eq (eq4 and eq14) differences with the numerical optimum, with H(rho) = 0.6 ', fontsize=16)

fig_4.savefig('plots/HeatmapOptimal-4VSHeatmapOptimal-14'  + '.png')


print("***************************************************")
print("*         PLOTS ARE SAVED IN /plots FOLDER         *")
print("***************************************************")

