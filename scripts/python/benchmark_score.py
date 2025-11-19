#!/usr/bin/env python3

import numpy as np
import numba
import time

# Original Python function
def score_python(
    reads_dna_normal=0,
    reads_dna_tumor=0,
    reads_rna_normal=0,
    reads_rna_tumor=0,
    reads_onekg=0,
    samples_dna_normal=0,
    samples_dna_tumor=0,
    samples_rna_normal=0,
    samples_rna_tumor=0,
    samples_onekg=0,
    pop_size_dna_normal=0,
    pop_size_dna_tumor=0,
    pop_size_rna_normal=0,
    pop_size_rna_tumor=0,
    pop_size_dna_onekg=2536,
    w_tumor=0.5,
    w_dna=0.5,
    w_read=0.5,
    upper_factor=50
):
    if not (reads_dna_normal or reads_dna_tumor or reads_rna_normal or reads_rna_tumor or reads_onekg or
            samples_dna_normal or samples_dna_tumor or samples_rna_normal or samples_rna_tumor or samples_onekg):
        return 0.0

    w_rna = 1.0 - w_dna
    w_normal = 1.0 - w_tumor
    w_sample = 1.0 - w_read

    score_val = 0.0
    
    # DNA Normal
    if pop_size_dna_normal > 0:
        m = upper_factor * pop_size_dna_normal
        g = (m - abs(reads_dna_normal - m)) / m if reads_dna_normal <= 2 * m else 0.0
        h = samples_dna_normal / pop_size_dna_normal
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
    
    # 1KG
    if pop_size_dna_onekg > 0:
        m = upper_factor * pop_size_dna_onekg
        g = (m - abs(reads_onekg - m)) / m if reads_onekg <= 2 * m else 0.0
        h = samples_onekg / pop_size_dna_onekg
        f = ((w_read * g) + (w_sample * h))
        score_val -= w_normal * w_dna * f

    return score_val


# Numba JIT-compiled version
@numba.jit(nopython=True)
def score_numba(
    reads_dna_normal,
    reads_dna_tumor,
    reads_rna_normal,
    reads_rna_tumor,
    reads_onekg,
    samples_dna_normal,
    samples_dna_tumor,
    samples_rna_normal,
    samples_rna_tumor,
    samples_onekg,
    pop_size_dna_normal,
    pop_size_dna_tumor,
    pop_size_rna_normal,
    pop_size_rna_tumor,
    pop_size_dna_onekg,
    w_tumor,
    w_dna,
    w_read,
    upper_factor
):
    if not (reads_dna_normal or reads_dna_tumor or reads_rna_normal or reads_rna_tumor or reads_onekg or
            samples_dna_normal or samples_dna_tumor or samples_rna_normal or samples_rna_tumor or samples_onekg):
        return 0.0

    w_rna = 1.0 - w_dna
    w_normal = 1.0 - w_tumor
    w_sample = 1.0 - w_read

    score_val = 0.0
    
    # DNA Normal
    if pop_size_dna_normal > 0:
        m = upper_factor * pop_size_dna_normal
        g = (m - abs(reads_dna_normal - m)) / m if reads_dna_normal <= 2 * m else 0.0
        h = samples_dna_normal / pop_size_dna_normal
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
    
    # 1KG
    if pop_size_dna_onekg > 0:
        m = upper_factor * pop_size_dna_onekg
        g = (m - abs(reads_onekg - m)) / m if reads_onekg <= 2 * m else 0.0
        h = samples_onekg / pop_size_dna_onekg
        f = ((w_read * g) + (w_sample * h))
        score_val -= w_normal * w_dna * f

    return score_val


# Vectorized version for numba
@numba.jit(nopython=True, parallel=True)
def score_numba_vectorized(
    reads_dna_normal_arr,
    reads_dna_tumor_arr,
    reads_rna_normal_arr,
    reads_rna_tumor_arr,
    reads_onekg_arr,
    samples_dna_normal_arr,
    samples_dna_tumor_arr,
    samples_rna_normal_arr,
    samples_rna_tumor_arr,
    samples_onekg_arr,
    pop_size_dna_normal_arr,
    pop_size_dna_tumor_arr,
    pop_size_rna_normal_arr,
    pop_size_rna_tumor_arr,
    pop_size_dna_onekg_arr,
    w_tumor_arr,
    w_dna_arr,
    w_read_arr,
    upper_factor_arr
):
    n = len(reads_dna_normal_arr)
    results = np.empty(n, dtype=np.float64)
    
    for i in numba.prange(n):
        results[i] = score_numba(
            reads_dna_normal_arr[i],
            reads_dna_tumor_arr[i],
            reads_rna_normal_arr[i],
            reads_rna_tumor_arr[i],
            reads_onekg_arr[i],
            samples_dna_normal_arr[i],
            samples_dna_tumor_arr[i],
            samples_rna_normal_arr[i],
            samples_rna_tumor_arr[i],
            samples_onekg_arr[i],
            pop_size_dna_normal_arr[i],
            pop_size_dna_tumor_arr[i],
            pop_size_rna_normal_arr[i],
            pop_size_rna_tumor_arr[i],
            pop_size_dna_onekg_arr[i],
            w_tumor_arr[i],
            w_dna_arr[i],
            w_read_arr[i],
            upper_factor_arr[i]
        )
    
    return results


def benchmark():
    print("=" * 80)
    print("PERFORMANCE BENCHMARK: score function")
    print("=" * 80)
    
    n = 10000
    print(f"\nGenerating {n} random test inputs...")
    
    # Generate random test data with realistic constraints
    np.random.seed(42)
    
    # First generate population sizes
    pop_size_dna_normal = np.random.randint(1, 200, n).astype(np.float64)
    pop_size_dna_tumor = np.random.randint(1, 200, n).astype(np.float64)
    pop_size_rna_normal = np.random.randint(1, 100, n).astype(np.float64)
    pop_size_rna_tumor = np.random.randint(1, 100, n).astype(np.float64)
    pop_size_dna_onekg = np.full(n, 2536.0)
    
    # Samples must be <= population size (realistic constraint)
    samples_dna_normal = np.random.randint(0, 1, n).astype(np.float64)
    samples_dna_tumor = np.random.randint(0, 1, n).astype(np.float64)
    samples_rna_normal = np.random.randint(0, 1, n).astype(np.float64)
    samples_rna_tumor = np.random.randint(0, 1, n).astype(np.float64)
    samples_onekg = np.random.randint(0, 2537, n).astype(np.float64)
    
    for i in range(n):
        samples_dna_normal[i] = np.random.randint(0, int(pop_size_dna_normal[i]) + 1)
        samples_dna_tumor[i] = np.random.randint(0, int(pop_size_dna_tumor[i]) + 1)
        samples_rna_normal[i] = np.random.randint(0, int(pop_size_rna_normal[i]) + 1)
        samples_rna_tumor[i] = np.random.randint(0, int(pop_size_rna_tumor[i]) + 1)
        samples_onekg[i] = np.random.randint(0, int(pop_size_dna_onekg[i]) + 1)
    
    # Reads must be >= samples (realistic constraint)
    reads_dna_normal = (samples_dna_normal + np.random.randint(0, 10000, n)).astype(np.float64)
    reads_dna_tumor = (samples_dna_tumor + np.random.randint(0, 10000, n)).astype(np.float64)
    reads_rna_normal = (samples_rna_normal + np.random.randint(0, 5000, n)).astype(np.float64)
    reads_rna_tumor = (samples_rna_tumor + np.random.randint(0, 5000, n)).astype(np.float64)
    reads_onekg = (samples_onekg + np.random.randint(0, 200000, n)).astype(np.float64)
    
    # Hyperparameters
    w_tumor = np.random.uniform(0.1, 0.9, n)
    w_dna = np.random.uniform(0.1, 0.9, n)
    w_read = np.random.uniform(0.1, 0.9, n)
    upper_factor = np.random.randint(10, 100, n).astype(np.float64)  # Positive integers
    
    # Warm-up: compile numba functions
    print("\nWarming up numba JIT compilation...")
    _ = score_numba(
        reads_dna_normal[0], reads_dna_tumor[0], reads_rna_normal[0], reads_rna_tumor[0],
        reads_onekg[0], samples_dna_normal[0], samples_dna_tumor[0], samples_rna_normal[0],
        samples_rna_tumor[0], samples_onekg[0], pop_size_dna_normal[0], pop_size_dna_tumor[0],
        pop_size_rna_normal[0], pop_size_rna_tumor[0], pop_size_dna_onekg[0],
        w_tumor[0], w_dna[0], w_read[0], upper_factor[0]
    )
    _ = score_numba_vectorized(
        reads_dna_normal[:10], reads_dna_tumor[:10], reads_rna_normal[:10], reads_rna_tumor[:10],
        reads_onekg[:10], samples_dna_normal[:10], samples_dna_tumor[:10], samples_rna_normal[:10],
        samples_rna_tumor[:10], samples_onekg[:10], pop_size_dna_normal[:10], pop_size_dna_tumor[:10],
        pop_size_rna_normal[:10], pop_size_rna_tumor[:10], pop_size_dna_onekg[:10],
        w_tumor[:10], w_dna[:10], w_read[:10], upper_factor[:10]
    )
    
    print("\n" + "-" * 80)
    print("BENCHMARK RESULTS")
    print("-" * 80)
    
    # 1. Pure Python loop
    print("\n1. Pure Python (loop)...")
    start = time.time()
    results_python = []
    for i in range(n):
        results_python.append(score_python(
            reads_dna_normal[i], reads_dna_tumor[i], reads_rna_normal[i], reads_rna_tumor[i],
            reads_onekg[i], samples_dna_normal[i], samples_dna_tumor[i], samples_rna_normal[i],
            samples_rna_tumor[i], samples_onekg[i], pop_size_dna_normal[i], pop_size_dna_tumor[i],
            pop_size_rna_normal[i], pop_size_rna_tumor[i], pop_size_dna_onekg[i],
            w_tumor[i], w_dna[i], w_read[i], upper_factor[i]
        ))
    time_python = time.time() - start
    results_python = np.array(results_python)
    print(f"   Time: {time_python:.4f} seconds")
    print(f"   Rate: {n/time_python:.0f} calls/sec")
    
    # 2. np.vectorize
    print("\n2. np.vectorize...")
    score_vectorized = np.vectorize(score_python)
    start = time.time()
    results_vectorize = score_vectorized(
        reads_dna_normal, reads_dna_tumor, reads_rna_normal, reads_rna_tumor,
        reads_onekg, samples_dna_normal, samples_dna_tumor, samples_rna_normal,
        samples_rna_tumor, samples_onekg, pop_size_dna_normal, pop_size_dna_tumor,
        pop_size_rna_normal, pop_size_rna_tumor, pop_size_dna_onekg,
        w_tumor, w_dna, w_read, upper_factor
    )
    time_vectorize = time.time() - start
    print(f"   Time: {time_vectorize:.4f} seconds")
    print(f"   Rate: {n/time_vectorize:.0f} calls/sec")
    print(f"   Speedup vs Python: {time_python/time_vectorize:.2f}x")
    
    # 3. Numba JIT (loop)
    print("\n3. Numba JIT (loop)...")
    start = time.time()
    results_numba = []
    for i in range(n):
        results_numba.append(score_numba(
            reads_dna_normal[i], reads_dna_tumor[i], reads_rna_normal[i], reads_rna_tumor[i],
            reads_onekg[i], samples_dna_normal[i], samples_dna_tumor[i], samples_rna_normal[i],
            samples_rna_tumor[i], samples_onekg[i], pop_size_dna_normal[i], pop_size_dna_tumor[i],
            pop_size_rna_normal[i], pop_size_rna_tumor[i], pop_size_dna_onekg[i],
            w_tumor[i], w_dna[i], w_read[i], upper_factor[i]
        ))
    time_numba = time.time() - start
    results_numba = np.array(results_numba)
    print(f"   Time: {time_numba:.4f} seconds")
    print(f"   Rate: {n/time_numba:.0f} calls/sec")
    print(f"   Speedup vs Python: {time_python/time_numba:.2f}x")
    
    # 4. Numba vectorized (parallel)
    print("\n4. Numba vectorized (parallel)...")
    start = time.time()
    results_numba_vec = score_numba_vectorized(
        reads_dna_normal, reads_dna_tumor, reads_rna_normal, reads_rna_tumor,
        reads_onekg, samples_dna_normal, samples_dna_tumor, samples_rna_normal,
        samples_rna_tumor, samples_onekg, pop_size_dna_normal, pop_size_dna_tumor,
        pop_size_rna_normal, pop_size_rna_tumor, pop_size_dna_onekg,
        w_tumor, w_dna, w_read, upper_factor
    )
    time_numba_vec = time.time() - start
    print(f"   Time: {time_numba_vec:.4f} seconds")
    print(f"   Rate: {n/time_numba_vec:.0f} calls/sec")
    print(f"   Speedup vs Python: {time_python/time_numba_vec:.2f}x")
    print(f"   Speedup vs Numba loop: {time_numba/time_numba_vec:.2f}x")
    
    # Verify results match
    print("\n" + "-" * 80)
    print("VERIFICATION")
    print("-" * 80)
    print(f"Max difference (Python vs np.vectorize): {np.max(np.abs(results_python - results_vectorize)):.2e}")
    print(f"Max difference (Python vs Numba): {np.max(np.abs(results_python - results_numba)):.2e}")
    print(f"Max difference (Python vs Numba vectorized): {np.max(np.abs(results_python - results_numba_vec)):.2e}")
    
    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"{'Method':<30} {'Time (s)':<12} {'Speedup':<12} {'Rate (calls/s)':<15}")
    print("-" * 80)
    print(f"{'Pure Python':<30} {time_python:<12.4f} {'1.00x':<12} {n/time_python:<15.0f}")
    print(f"{'np.vectorize':<30} {time_vectorize:<12.4f} {f'{time_python/time_vectorize:.2f}x':<12} {n/time_vectorize:<15.0f}")
    print(f"{'Numba JIT':<30} {time_numba:<12.4f} {f'{time_python/time_numba:.2f}x':<12} {n/time_numba:<15.0f}")
    print(f"{'Numba Vectorized (parallel)':<30} {time_numba_vec:<12.4f} {f'{time_python/time_numba_vec:.2f}x':<12} {n/time_numba_vec:<15.0f}")
    print("=" * 80)


if __name__ == '__main__':
    benchmark()
