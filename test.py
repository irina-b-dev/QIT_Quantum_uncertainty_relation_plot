import matplotlib.pyplot as plt
import math


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


c = [i/1000 for i in range(500, 1000)]
h_ro = [0, 0.3, 0.6, 0.9]
for H_ro in h_ro:
    P = purity_calculator(H_ro)
    lb4 = lower_bound_eq4(H_ro, c)
    lb12 = lower_bound_eq12(H_ro, c)
    lb13 = lower_bound_eq13(H_ro, c)
    lb14 = lower_bound_eq14(P, H_ro, c)

    plt.plot(c, lb4, label = "Eq. (4)")
    plt.plot(c, lb12, label = "Eq. (12)")
    plt.plot(c, lb13, label = "Eq. (13)")
    plt.plot(c, lb14, label = "Eq. (14)")

    plt.xlabel('c')
    plt.ylabel('Lower bound')
    plt.title('H(ro)=' + str(H_ro))
    plt.legend()
    plt.show()