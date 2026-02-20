import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.patches import Patch
import math
import seaborn as sns

file=open("archivo.txt", "r")
phi_l_2=[]
phi_l=[]

for line in file:
    temp=[]
    temp.extend(map(float, line.split(',')))
    phi_l_2.append(temp[0])
    phi_l.append(temp[1])
        
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Primera gráfica
ax1.plot([i*2+1 for i in range(len(phi_l_2))], phi_l_2, color='blue')
ax1.set_title(r'$\phi_{l/2}$, temperatura en la mitad de la barra')
ax1.set_xlabel('N nodos')
ax1.set_ylabel(r'$\phi$')
#ax1.set_ylim(0.0, 7)

# Segunda gráfica
ax2.plot([i*2+1 for i in range(len(phi_l))], phi_l, color='red')
ax2.set_title(r'$\phi_l$, temperatura en un extremo de la barra')
ax2.set_xlabel('N nodos')
ax2.set_ylabel(r'$\phi$')




file=open("archivo_temps.txt", "r")

temps=[]
for line in file:
        temp=[]
        
        temp.extend(map(float, line.split(',')))
        temps.append(temp)

def set_node_ticks(ax, n_nodes, L):
    xticks = np.linspace(0, L, 4)
    labels = np.linspace(0, n_nodes, 4, dtype=int)

    ax.set_xticks(xticks)
    ax.set_xticklabels(labels)

fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, figsize=(8, 4), sharex=False)
fig.suptitle(r'Valores de $\phi$ para diferente cantidad de nodos')

axes=[ax1, ax2, ax3, ax4, ax5]
indices=[0, 10, 250, 375, 500]
titles=[f"{indices[0]} nodos", 
        f"{indices[1]} nodos", 
        f"{indices[2]} nodos", 
        f"{indices[3]} nodos", 
        f"{indices[4]} nodos"]

i=0
for ax, idx in zip(axes, indices):
    im=ax.imshow([temps[idx]], cmap='hot', interpolation='nearest', extent=[0, 10.0, 0, 0.3])
    ax.set_title(titles[i], loc="left", fontsize=8)
    ax.set_yticks([])
    i+=1

    set_node_ticks(ax, idx+3, 10.0)

cbar = fig.colorbar(im, ax=axes, orientation='vertical')
cbar.set_label(r"$\phi$")

plt.show()