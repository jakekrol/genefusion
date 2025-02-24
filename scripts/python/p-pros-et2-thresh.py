#!/usr/bin/env python3

import math
import matplotlib.pyplot as plt

MIN_REPORTED_FREQ=0.2
MAX_REPORTED_FREQ=0.6
f = '/data/jake/genefusion/data/2024_11_10-fusions-pcawg-combined/25_02_24-sample_erg_tmprss2/f.dist'
d=[]
with open(f) as f:
    for line in f:
        d.append(int(line.strip()))
# print(d)

n=len(d)
t=list(range(0,101))

ys=[]
for i in t:
    y=0
    for j in d:
        if j>=i:
            y+=1
    ys.append(y)

fig,(ax1,ax2) = plt.subplots(1,2)
ax1.plot(t,ys, color = 'r')
ax1.set_xlabel('Threshold')
ax1.set_ylabel('Number of samples above threshold')
ys_freq = [y/n for y in ys]
ax2.plot(t,ys_freq, color = 'b')
ax2.set_xlabel('Threshold')
ax2.set_ylabel('Frequency of samples above threshold')
ax2.axhline(y=MIN_REPORTED_FREQ, color='g', linestyle='--', label='Min reported pop. freq')
ax2.axhline(y=MAX_REPORTED_FREQ, color='orange', linestyle='--', label='Max reported pop. freq')
ax2.legend()
plt.subplots_adjust(wspace=0.4)  # Increase the width space between subplots
plt.savefig('test2.png')

# Sample data
# import numpy as np
# x = np.arange(0, 10, 0.1)
# y1 = np.sin(x)
# y2 = np.cos(x)

# fig, ax1 = plt.subplots()

# # Plot on the primary y-axis
# ax1.plot(x, y1, 'g-')
# ax1.set_xlabel('X data')
# ax1.set_ylabel('Sin', color='g')

# # Create a secondary y-axis
# ax2 = ax1.twinx()
# ax2.plot(x, y2, 'b-')
# ax2.set_ylabel('Cos', color='b')

# plt.title('Sin and Cos on different y-axes')
# plt.savefig('test.png')
