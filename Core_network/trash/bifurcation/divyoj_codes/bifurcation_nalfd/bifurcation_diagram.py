import matplotlib.pyplot as plt
import numpy as np

path_to_file="./matlab_output/steady_state_pparg_01.txt"
k_pparg=[]
y=[]
with open (path_to_file) as file:
    for line in file:
        a=line[:-1].split("\t")
        k_pparg.append(a[0])
        steady_state=a[1:-1]
        y.append(tuple(steady_state))


x =k_pparg

for xe, ye in zip(x, y):
    plt.scatter([xe] * len(ye), ye,s=6, c="g",alpha=0.5)    
plt.yticks([])
plt.xticks(np.arange(0, 5, 1.0))
plt.show()
    