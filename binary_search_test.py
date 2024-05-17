import matplotlib.pyplot as plt
import math
import numpy as np

def p_entropy_func(p):
    if p == 0 or p == 1:
        return 0
    return -p * math.log2(p) - (1-p) * math.log2(1-p)

def p_search(H_ro):
    error = 0.001
    p_list = [i/100000 for i in range(1, 50000)]
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
print(p_search(0.3))
