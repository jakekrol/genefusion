#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

# Simple normal distributions with low variance
dist1 = norm(loc=-0.7, scale=0.2)
dist2 = norm(loc=0, scale=0.2)
dist3 = norm(loc=0.7, scale=0.2)

x = np.linspace(-1, 1, 500)
pdf1 = dist1.pdf(x)
pdf2 = dist2.pdf(x)
pdf3 = dist3.pdf(x)

# Plot
plt.figure(figsize=(10, 6))
plt.fill_between(x, pdf1, alpha=0.3, color='tab:blue', label='Recurrent normal')
plt.fill_between(x, pdf2, alpha=0.3, color='tab:orange', label='Random gene pairs')
plt.fill_between(x, pdf3, alpha=0.3, color='tab:green', label='Recurrent tumor')
plt.plot(x, pdf1, linewidth=2.5, color='tab:blue')
plt.plot(x, pdf2, linewidth=2.5, color='tab:orange')
plt.plot(x, pdf3, linewidth=2.5, color='tab:green')
plt.xlim(-1, 1)
plt.xlabel('Value', fontsize=12)
plt.ylabel('Density', fontsize=12)
plt.legend(fontsize=11, framealpha=0.95)
plt.title('Score function evaluation', fontsize=14)
plt.tight_layout()
plt.savefig('ideal_score_eval.png')