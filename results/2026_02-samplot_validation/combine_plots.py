#!/usr/bin/env python3

import os, sys
import glob
import argparse
import pandas as pd
from PIL import Image, ImageDraw, ImageFont

parser=argparse.ArgumentParser(description="Combine plots into a single PDF")
parser.add_argument('-i', '--input_dir', help='Input directory with plots', required=True)
parser.add_argument('-o', '--output_file', help='Output PDF file', required=True)
args = parser.parse_args()

def fname2data(x):
    tissue = x.split('_')[0]
    genepair=x.split('_')[1]
    leftgene=genepair.split('--')[0]
    rightgene=genepair.split('--')[1]
    svtype=x.split('_')[2]
    total_evidence=x.split('_')[3].split('.')[0]
    # rm ev prefix from evidence string
    if total_evidence.startswith('ev'):
        total_evidence = total_evidence[2:]
    return tissue, leftgene, rightgene, svtype, total_evidence, x

data = [ fname2data(x) for x in os.listdir(args.input_dir) if x.endswith('.png') ]
df = pd.DataFrame(data, columns=['tissue','leftgene','rightgene','svtype','total_evidence','filename'])
df = df.sort_values(by=['tissue','svtype','total_evidence'], ascending=[True, True, False])

# Add filename text to each image and combine into PDF
images_with_text = []
for i, filename in enumerate(df['filename']):
    print(f"Processing {filename} ({i+1}/{len(df)})")
    img_path = f"{args.input_dir}/{filename}"
    img = Image.open(img_path)
    
    # Convert to RGB if needed (for PDF compatibility)
    if img.mode != 'RGB':
        img = img.convert('RGB')
    
    # Create a drawing context
    draw = ImageDraw.Draw(img)
    
    # Try to use a large font; fall back to default if not available
    try:
        font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 48)
    except:
        try:
            font = ImageFont.truetype("Arial.ttf", 48)
        except:
            font = ImageFont.load_default()
    
    # Add text at the top with background for readability
    text = filename
    bbox = draw.textbbox((0, 0), text, font=font)
    text_width = bbox[2] - bbox[0]
    text_height = bbox[3] - bbox[1]
    
    # Position text at top center
    x = (img.width - text_width) // 2
    y = 10
    
    # Draw white background rectangle for text
    padding = 10
    draw.rectangle([x - padding, y - padding, x + text_width + padding, y + text_height + padding], 
                   fill='white', outline='black', width=2)
    
    # Draw the filename text in black
    draw.text((x, y), text, fill='black', font=font)
    
    images_with_text.append(img)

# Save all images as a single PDF
images_with_text[0].save(args.output_file, save_all=True, append_images=images_with_text[1:], resolution=100.0)