#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

FONTSIZE=15
total_samples=10

# sample score
x=list(range(0, total_samples +1))
y=[i/total_samples for i in x]
xlabel='# of supporting samples'
ylabel='h(#samples/ #total_samples)'
plt.plot(x,y)
plt.text(
    0.4 * total_samples,
    0.2,
    '#total_samples={}'.format(total_samples),
    fontsize=FONTSIZE
)
plt.xlabel(xlabel, fontsize=FONTSIZE)
plt.ylabel(ylabel, fontsize=FONTSIZE)
plt.title('Sample score (tumor or normal)', fontsize=FONTSIZE +2)
# make ticklabels larger
plt.xticks(fontsize=FONTSIZE -3)
plt.yticks(fontsize=FONTSIZE -3)
plt.savefig('sample_score.png')
plt.close()

# read score tumor
def g_tumor(reads,m):
    if reads <= (2*m):
        return (m-abs(reads - m))/m
    else:
        return 0
cov=50
m=cov*total_samples
x=np.linspace(0, 3*m, 100)
y=[g_tumor(r,m) for r in x]
xlabel='# of tumor supporting reads'
ylabel='g_tumor(#reads)'
plt.plot(x,y)
plt.text(
    750,
    0.6,
    '#total_samples={}\n#avg_coverage={}\nm={}'.format(total_samples, cov, m),
    fontsize=FONTSIZE
)
plt.xlabel(xlabel, fontsize=FONTSIZE)
plt.ylabel(ylabel, fontsize=FONTSIZE)
plt.title('Tumor read score', fontsize=FONTSIZE +2)
plt.xticks(fontsize=FONTSIZE -3)
plt.yticks(fontsize=FONTSIZE -3)
plt.savefig('tumor_read_score.png')
plt.close()


# read score normal
x=np.linspace(0, 3*m, 100)
y=[i/m for i in x]
xlabel='# of normal supporting reads'
ylabel='g_normal(#reads)'
plt.plot(x,y)
plt.text(
    800,
    0.6,
    '#total_samples={}\n#avg_coverage={}\nm={}'.format(total_samples, cov, m),
    fontsize=FONTSIZE
)
plt.xlabel(xlabel, fontsize=FONTSIZE)
plt.ylabel(ylabel, fontsize=FONTSIZE)
plt.title('Normal read score', fontsize=FONTSIZE +2)
plt.xticks(fontsize=FONTSIZE -3)
plt.yticks(fontsize=FONTSIZE -3)
plt.savefig('normal_read_score.png')
    

plt.close()