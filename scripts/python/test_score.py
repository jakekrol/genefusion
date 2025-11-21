#!/usr/bin/env python3



def score(
    # reads
    reads_dna_normal=0,
    reads_dna_tumor=0,
    reads_rna_normal=0,
    reads_rna_tumor=0,
    reads_onekg=0,

    # sample counts
    samples_dna_normal=0,
    samples_dna_tumor=0,
    samples_rna_normal=0,
    samples_rna_tumor=0,
    samples_onekg=0,
    # pop size
    pop_size_dna_normal=0,
    pop_size_dna_tumor=0,
    pop_size_rna_normal=0,
    pop_size_rna_tumor=0,
    pop_size_dna_onekg=2535,
    # hyperparameters
    w_tumor=0.5,
    w_dna=0.5,
    w_read=0.5,
    upper_factor=50 # max expected read count per sample
):
    # Fast exit if all inputs are zero
    if not (reads_dna_normal or reads_dna_tumor or reads_rna_normal or reads_rna_tumor or reads_onekg or
            samples_dna_normal or samples_dna_tumor or samples_rna_normal or samples_rna_tumor or samples_onekg):
        return 0.0

    # Precompute weights
    w_rna = 1.0 - w_dna
    w_normal = 1.0 - w_tumor
    w_sample = 1.0 - w_read

    score_val = 0.0
    
    # Inline all calculations to avoid function call overhead
    # Process each modality directly without loops or list creation
    
    # DNA Normal
    if pop_size_dna_normal > 0:
        m = upper_factor * pop_size_dna_normal
        # Inline read_support
        g = (m - abs(reads_dna_normal - m)) / m if reads_dna_normal <= 2 * m else 0.0
        # Inline sample_support
        h = samples_dna_normal / pop_size_dna_normal
        # Inline modality_score
        f = ((w_read * g) + (w_sample * h))
        score_val -= w_normal * w_dna * f
    
    # DNA Tumor
    if pop_size_dna_tumor > 0:
        m = upper_factor * pop_size_dna_tumor
        g = (m - abs(reads_dna_tumor - m)) / m if reads_dna_tumor <= 2 * m else 0.0
        h = samples_dna_tumor / pop_size_dna_tumor
        f = ((w_read * g) + (w_sample * h))
        score_val += w_tumor * w_dna * f
    
    # RNA Normal
    if pop_size_rna_normal > 0:
        m = upper_factor * pop_size_rna_normal
        g = (m - abs(reads_rna_normal - m)) / m if reads_rna_normal <= 2 * m else 0.0
        h = samples_rna_normal / pop_size_rna_normal
        f = ((w_read * g) + (w_sample * h))
        score_val -= w_normal * w_rna * f
    
    # RNA Tumor
    if pop_size_rna_tumor > 0:
        m = upper_factor * pop_size_rna_tumor
        g = (m - abs(reads_rna_tumor - m)) / m if reads_rna_tumor <= 2 * m else 0.0
        h = samples_rna_tumor / pop_size_rna_tumor
        f = ((w_read * g) + (w_sample * h))
        score_val += w_tumor * w_rna * f
    
    # 1KG (DNA, normal population)
    if pop_size_dna_onekg > 0:
        m = upper_factor * pop_size_dna_onekg
        g = (m - abs(reads_onekg - m)) / m if reads_onekg <= 2 * m else 0.0
        h = samples_onekg / pop_size_dna_onekg
        f = ((w_read * g) + (w_sample * h))
        score_val -= w_normal * w_dna * f

    return score_val
# example
print('scoring x')
x = score(
    reads_dna_normal=0,
    reads_dna_tumor=5000,
    reads_rna_normal=0,
    reads_rna_tumor=0,  
    reads_onekg=0,
    samples_dna_normal=0,
    samples_dna_tumor=100,
    samples_rna_normal=0,
    samples_rna_tumor=0,
    samples_onekg=0,
    pop_size_dna_normal=0,
    pop_size_dna_tumor=100,
    pop_size_rna_normal=0,
    pop_size_rna_tumor=0,
    w_dna=1.0,
    w_tumor=1.0,
    upper_factor=50
)
print('score max tumor dna evidence (w_dna=1,w_tumor=1,u=50):', x)

print('scoring y')
y = score(
    reads_dna_normal=5000,
    reads_dna_tumor=5000,
    reads_rna_normal=0,
    reads_rna_tumor=0,  
    reads_onekg=0,
    samples_dna_normal=100,
    samples_dna_tumor=100,
    samples_rna_normal=0,
    samples_rna_tumor=0,
    samples_onekg=0,
    pop_size_dna_normal=100,
    pop_size_dna_tumor=100,
    pop_size_rna_normal=0,
    pop_size_rna_tumor=0,
    w_dna=1.0,
    w_tumor=0.5,
    w_read=0.5,
    upper_factor=50
)
print('score equal tumor and normal dna evidence (w_dna=1, w_tumor=0.5, w_read=0.5, u=50)', y)

# Test case 3: Max normal DNA evidence
# Expected: -1.0
print('\nscoring z')
z = score(
    reads_dna_normal=5000,
    reads_dna_tumor=0,
    samples_dna_normal=100,
    samples_dna_tumor=0,
    pop_size_dna_normal=100,
    pop_size_dna_tumor=0,
    pop_size_dna_onekg=0,
    w_dna=1.0,
    w_tumor=0.0,
    upper_factor=50
)
print('score max normal dna evidence (w_dna=1, w_tumor=0, u=50):', z, '(expected: -1.0)')

# Test case 4: Max RNA tumor evidence
# Expected: 1.0
print('\nscoring a')
a = score(
    reads_rna_tumor=5000,
    samples_rna_tumor=100,
    pop_size_rna_tumor=100,
    pop_size_dna_onekg=0,
    w_dna=0.0,
    w_tumor=1.0,
    upper_factor=50
)
print('score max tumor rna evidence (w_dna=0, w_tumor=1, u=50):', a)

# Test case 5: Balanced DNA and RNA tumor evidence
print('\nscoring b')
b = score(
    reads_dna_tumor=2500,
    reads_rna_tumor=2500,
    samples_dna_tumor=50,
    samples_rna_tumor=50,
    pop_size_dna_tumor=50,
    pop_size_rna_tumor=50,
    pop_size_dna_onekg=0,
    w_dna=0.5,
    w_tumor=1.0,
    upper_factor=50
)
print('score balanced dna/rna tumor (w_dna=0.5, w_tumor=1, u=50):', b)

# Test case 6: 1KG population evidence (should be negative)
print('\nscoring c')
c = score(
    reads_onekg=126800,
    samples_onekg=2536,
    pop_size_dna_onekg=2536,
    w_tumor=0.5,
    upper_factor=50
)
print('score max 1kg evidence (default weights, u=50):', c)

# Test case 7: Read-weighted vs sample-weighted
print('\nscoring d')
d = score(
    reads_dna_tumor=5000,
    samples_dna_tumor=100,
    pop_size_dna_tumor=100,
    pop_size_dna_onekg=0,
    w_dna=1.0,
    w_tumor=1.0,
    w_read=1.0,  # Only reads matter
    upper_factor=50
)
print('score tumor dna with w_read=1.0 (only reads):', d)

print('\nscoring e')
e = score(
    reads_dna_tumor=100,
    samples_dna_tumor=100,
    pop_size_dna_tumor=100,
    pop_size_dna_onekg=0,
    w_dna=1.0,
    w_tumor=1.0,
    w_read=0.0,  # Only samples matter
    upper_factor=50
)
print('score tumor dna with w_read=0.0 (only samples):', e)

# Test case 8: No evidence (should be 0.0)
print('\nscoring f')
f = score()
print('score no evidence (all defaults):', f)