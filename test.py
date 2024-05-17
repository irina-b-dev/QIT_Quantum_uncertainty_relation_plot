import matplotlib.pyplot as plt
import math
import numpy as np

# ro = p|r><r| + (1-p)|r_(perpendicular)><r_(perpendicular)|
# Purity = Tr[ro^2] = 2p^2 - 2p + 1

def p_entropy_func(p):
    if p == 0 or p == 1:
        return 0
    return -p * math.log2(p) - (1-p) * math.log2(1-p)

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
        binary_entropy.append(-value*math.log2(value)-(1-value)*math.log2(1-value)-H_ro)
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
            binary_entropy.append(bn_entr)
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



fig_1, axs_1 = plt.subplots(2, 2, figsize=(10, 8))  # 2 rows, 2 columns

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

    axs_1[ip%2,ip//2]
    axs_1[ip%2,ip//2].plot(c, lb4, label = "Eq. (4)")
    axs_1[ip%2,ip//2].plot(c, lb12, label = "Eq. (12)")
    axs_1[ip%2,ip//2].plot(c, lb13, label = "Eq. (13)")
    axs_1[ip%2,ip//2].plot(c, lb14, label = "Eq. (14)")
    axs_1[ip%2,ip//2].plot(c, no, label = "Numerical optimum")

    axs_1[ip%2,ip//2].set_xlabel('c')
    axs_1[ip%2,ip//2].set_ylabel('Lower bound')
    axs_1[ip%2,ip//2].set_title('H(ro)=' + str(H_ro))
    axs_1[ip%2,ip//2].legend()
    ip = ip + 1

fig_1.savefig('plots/H(rho)'  + '.png', bbox_inches='tight')

H_ro = 0.6
P, p = purity_calculator(H_ro)
num_rows = 25
num_cols = 20
lb4 = lower_bound_eq4(H_ro, c)
lb12 = lower_bound_eq12(H_ro, c)
lb13 = lower_bound_eq13(H_ro, c)
lb14 = lower_bound_eq14(P, H_ro, c)
difflb4_12 = [x - y for x, y in zip(lb4, lb12)]
difflb4_13 = [x - y for x, y in zip(lb4, lb13)]
difflb4_14 = [x - y for x, y in zip(lb4, lb14)]
difflb12_13 = [x - y for x, y in zip(lb12, lb13)]
difflb12_14 = [x - y for x, y in zip(lb12, lb14)]
difflb13_14 = [x - y for x, y in zip(lb13, lb14)]

list_vals_plot = [lb4, lb12, lb13, lb14, difflb4_12, difflb4_13, difflb4_14, difflb12_13, difflb12_14, difflb13_14]

fig_2, axs_2 = plt.subplots(2, 5, figsize=(8*5, 6*2))

for i, data in enumerate([lb4, lb12, lb13, lb14, difflb4_12, difflb4_13, difflb4_14, difflb12_13, difflb12_14, difflb13_14]):
    matrix = np.array(data).reshape(num_rows, num_cols)
    im = axs_2[i%2,i%5].imshow(matrix, cmap='gray', interpolation='nearest')
    fig_2.colorbar(im, ax=axs_2[i%2,i%5])
    axs_2[i%2,i%5].set_title(f"Heatmap {i+1}")

fig_2.savefig('plots/Heatmaps'  + '.png')


print("***************************************************")
print("*         PLOTS ARE SAVED IN /plots FOLDER         *")
print("***************************************************")

