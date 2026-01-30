#!/usr/bin/env python3
import argparse
import matplotlib.pyplot as plt
import numpy as np

def main():
    parser = argparse.ArgumentParser(description='Create a scatter plot with colored points')
    parser.add_argument('--point_size', type=float, default=1, help='Point size (default: 1)')
    parser.add_argument('--xlabel', default='X', help='X-axis label')
    parser.add_argument('--ylabel', default='Y', help='Y-axis label') 
    parser.add_argument('--title', default='', help='Plot title')
    parser.add_argument('--fontsize', type=int, default=12, help='Base font size')
    parser.add_argument('--input', required=True, help='Input file (3 columns: x, y, color)')
    parser.add_argument('--output', required=True, help='Output PNG filename')
    parser.add_argument('--annotate', default='', help='Manual annotations: "x1:y1:color1:size1:label1,x2:y2:color2:size2:label2,..."')
    
    args = parser.parse_args()
    
    # Read data with mixed types (numeric x,y and string colors)
    x, y, color_values = [], [], []
    with open(args.input, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                x.append(float(parts[0]))
                y.append(float(parts[1]))
                color_values.append(parts[2])  # Keep as string
    
    x, y = np.array(x), np.array(y)
    
    # Handle discrete colors and markers (assume 2 categories)
    unique_colors = np.unique(color_values)
    colors = ['red', 'blue']
    markers = ['o', 's']  # circle and square
    color_map = {val: colors[i] for i, val in enumerate(unique_colors)}
    marker_map = {val: markers[i] for i, val in enumerate(unique_colors)}
    
    # Create plot
    fig, ax = plt.subplots()
    
    # Plot each category separately with different markers
    # Plot non-pcawg categories first, then pcawg on top
    categories_ordered = [cat for cat in unique_colors if cat != 'pcawg'] + [cat for cat in unique_colors if cat == 'pcawg']
    
    for category in categories_ordered:
        mask = np.array([val == category for val in color_values])
        zorder = 5 if category == 'pcawg' else 1  # pcawg on top
        ax.scatter(x[mask], y[mask], 
                  c=color_map[category], 
                  marker=marker_map[category],
                  s=args.point_size,
                  alpha=1,
                  label=category,
                  zorder=zorder)
    
    # Set labels and title with font sizes
    if args.title:
        ax.set_title(args.title, fontsize=args.fontsize)
    ax.set_xlabel(args.xlabel, fontsize=int(args.fontsize * 0.8))
    ax.set_ylabel(args.ylabel, fontsize=int(args.fontsize * 0.8))
    ax.tick_params(labelsize=int(args.fontsize * 0.7))
    
    # Add legend (already created by scatter with label parameter)
    legend_handles, legend_labels = ax.get_legend_handles_labels()
    
    # Add manual annotations if specified
    annotation_info = {}  # {color: (label, marker_info)}
    if args.annotate:
        for annotation in args.annotate.split(','):
            parts = annotation.strip().split(':')
            if len(parts) == 3:  # Old format: x:y:color
                ax_x, ax_y, ax_color = float(parts[0]), float(parts[1]), parts[2]
                ax_size = args.point_size * 3
                ax_label = f'annotation ({ax_color})'
            elif len(parts) == 4:  # Format: x:y:color:size
                ax_x, ax_y, ax_color, ax_size = float(parts[0]), float(parts[1]), parts[2], float(parts[3])
                ax_label = f'annotation ({ax_color})'
            elif len(parts) == 5:  # New format: x:y:color:size:label
                ax_x, ax_y, ax_color, ax_size, ax_label = float(parts[0]), float(parts[1]), parts[2], float(parts[3]), parts[4]
            else:
                continue
            ax.scatter(ax_x, ax_y, c=ax_color, s=ax_size, marker='*', 
                      edgecolors='black', linewidths=1, zorder=10)
            annotation_info[ax_color] = ax_label
    
    # Add annotation colors to legend
    for color in sorted(annotation_info.keys()):
        legend_handles.append(plt.Line2D([0], [0], marker='*', color='w', 
                                        markerfacecolor=color, markeredgecolor='black',
                                        markersize=8, label=annotation_info[color]))
    
    ax.legend(handles=legend_handles, fontsize=int(args.fontsize * 0.7))
    
    # Save plot
    plt.tight_layout()
    plt.savefig(args.output, dpi=300, bbox_inches='tight')

if __name__ == '__main__':
    main()
