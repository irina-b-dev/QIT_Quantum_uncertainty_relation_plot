import matplotlib.pyplot as plt
import math
import numpy as np

# ro = p|r><r| + (1-p)|r_(perpendicular)><r_(perpendicular)|
# Purity = Tr[ro^2] = 2p^2 - 2p + 1


def purity_calculator(H_ro):
    if H_ro == 0:
        return 1
    elif H_ro == 0.3:
        p = 0.053239
        #p = 0.946761
    elif H_ro == 0.6:
        #p = 0.853898
        p = 0.146102
    elif H_ro == 0.9:
        p = 0.316019
        # p = 0.683981
    return 2*(p**2) - 2*p +1

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

def numerical_optimum(H_ro, c):
    p = p_calculator(H_ro)
    binary_entropy = []
    for elem in c:
        value = p * math.sqrt(elem) + (1-p) * (1-math.sqrt(elem))
        bn_entr1 = -value*math.log2(value)-(1-value)*math.log2(1-value)
        bn_entr2 = -p*math.log2(p)-(1-p)*math.log2(1-p)
        binary_entropy.append(bn_entr1-bn_entr2)
    return binary_entropy

def l1_entropy(H_ro, P, c):
    l1 = []
    for elem in c:
        value = 2 * math.sqrt(2 * P - 1) * math.sqrt(elem * (1 - elem))
        l1.append(value)
    return l1


c = [i/1000 for i in range(500, 1000)]
h_ro = [0, 0.3, 0.6, 0.9]
for H_ro in h_ro:
    P = purity_calculator(H_ro)
    lb4 = lower_bound_eq4(H_ro, c)
    lb12 = lower_bound_eq12(H_ro, c)
    lb13 = lower_bound_eq13(H_ro, c)
    lb14 = lower_bound_eq14(P, H_ro, c)
    #no = numerical_optimum(H_ro, c)

    plt.plot(c, lb4, label = "Eq. (4)")
    plt.plot(c, lb12, label = "Eq. (12)")
    plt.plot(c, lb13, label = "Eq. (13)")
    plt.plot(c, lb14, label = "Eq. (14)")
    #plt.plot(c, no, label = "Numerical optimum")

    plt.xlabel('c')
    plt.ylabel('Lower bound')
    plt.title('H(ro)=' + str(H_ro))
    plt.legend()
    plt.show()

H_ro = 0.6
P = purity_calculator(H_ro)
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

for i, data in enumerate([lb4, lb12, lb13, lb14, difflb4_12, difflb4_13, difflb4_14, difflb12_13, difflb12_14, difflb13_14]):
    matrix = np.array(data).reshape(num_rows, num_cols)
    plt.figure(figsize=(8, 6))
    plt.imshow(matrix, cmap='gray', interpolation='nearest')
    plt.colorbar()
    plt.title(f"Heatmap {i+1}")
    plt.show()

