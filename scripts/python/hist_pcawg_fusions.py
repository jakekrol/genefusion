#!/usr/bin/env python3
import argparse
import matplotlib.pyplot as plt
import numpy as np

def main():
    parser = argparse.ArgumentParser(description='Plot overlaid histograms of two distributions')
    parser.add_argument('--pcawg', required=True, help='PCAWG fusions file (line delimited values)')
    parser.add_argument('--null', required=True, help='Null fusions file (line delimited values)')
    parser.add_argument('--output', required=True, help='Output PNG file')
    parser.add_argument('--title', default='Distribution Comparison', help='Plot title')
    parser.add_argument('--xlabel', default='Value', help='X-axis label')
    parser.add_argument('--ylabel', default='Frequency', help='Y-axis label')
    parser.add_argument('--fontsize', type=int, default=12, help='Base font size')
    parser.add_argument('--bins', type=int, default=30, help='Number of histogram bins')
    parser.add_argument('--vlines', default='', help='Vertical lines: "x1:color1:label1,x2:color2:label2,..."')
    parser.add_argument('--density', action='store_true', help='Use density instead of counts')
    
    args = parser.parse_args()
    
    # Read data files
    with open(args.pcawg) as f:
        pcawg_data = [float(line.strip()) for line in f if line.strip()]
    
    with open(args.null) as f:
        null_data = [float(line.strip()) for line in f if line.strip()]
    
    # Create plot
    fig, ax = plt.subplots()
    
    # Plot histograms with better visibility for overlapping regions
    if args.density:
        # For proportions (y-values 0-1), use weights instead of density=True
        null_weights = np.ones(len(null_data)) / len(null_data)
        pcawg_weights = np.ones(len(pcawg_data)) / len(pcawg_data)
        ax.hist(null_data, bins=args.bins, alpha=0.5, color='lightblue', 
                edgecolor='blue', linewidth=1, label='Null fusions', weights=null_weights)
        ax.hist(pcawg_data, bins=args.bins, histtype='step', color='red', 
                linewidth=3, label='PCAWG fusion calls', weights=pcawg_weights)
    else:
        # Regular count mode
        ax.hist(null_data, bins=args.bins, alpha=0.5, color='lightblue', 
                edgecolor='blue', linewidth=1, label='Null fusions')
        ax.hist(pcawg_data, bins=args.bins, histtype='step', color='red', 
                linewidth=3, label='PCAWG fusion calls')
    
    # Set labels and title with font sizes
    ax.set_title(args.title, fontsize=args.fontsize)
    ax.set_xlabel(args.xlabel, fontsize=int(args.fontsize * 0.8))
    ax.set_ylabel(args.ylabel, fontsize=int(args.fontsize * 0.8))
    ax.tick_params(labelsize=int(args.fontsize * 0.7))
    
    # Add vertical lines if specified
    if args.vlines:
        for vline in args.vlines.split(','):
            parts = vline.strip().split(':')
            if len(parts) == 2:
                x, color = float(parts[0]), parts[1]
                ax.axvline(x, color=color, linestyle='--', linewidth=2)
            elif len(parts) == 3:
                x, color, label = float(parts[0]), parts[1], parts[2]
                ax.axvline(x, color=color, linestyle='--', linewidth=2, label=label)
    
    # Add legend (after vertical lines so they appear in legend)
    ax.legend(fontsize=int(args.fontsize * 0.7))
    
    # Save plot
    plt.tight_layout()
    plt.savefig(args.output, dpi=300, bbox_inches='tight')

if __name__ == '__main__':
    main()
