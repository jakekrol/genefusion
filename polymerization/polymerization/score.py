import numba
import numpy as np
import math
import pandas as pd

def normalize_evidence_columns(
	df,
	column_map
):
	df_norm = df.copy()
	# validation
	for key_column_in in column_map.keys():
		assert key_column_in in df_norm.columns, f"Expected column '{key_column_in}' not found in DataFrame."
	for key, sub_dict in column_map.items():
		assert 'evidence_type' in sub_dict, f"Missing 'evidence_type' in column_map for key '{key}'."
		assert 'specimen' in sub_dict, f"Missing 'specimen' in column_map for key '{key}'."
		assert 'total_samples' in sub_dict, f"Missing 'total_samples' in column_map for key '{key}'."

	for key_column_in, sub_dict in column_map.items():
		# get parameters
		evidence_type = sub_dict['evidence_type']
		specimen = sub_dict['specimen']
		total_samples = sub_dict['total_samples']
		if "upper_bound" in sub_dict.keys():
			upper_bound = sub_dict['upper_bound']
		# apply normalization
		if evidence_type == 'sample':
			df_norm[key_column_in] = sample_score_vectorized(df_norm[key_column_in].values, total_samples)
		if evidence_type == 'read':
			if specimen == 'normal':
				df_norm[key_column_in] = read_score_normal_vectorized(df[key_column_in].values, total_samples, upper_bound)
			elif specimen == 'tumor':
				df_norm[key_column_in] = read_score_tumor_vectorized(df[key_column_in].values, total_samples, upper_bound)
			else:
				raise ValueError(f"Invalid specimen type '{specimen}' for column '{key_column_in}'. Expected 'normal' or 'tumor'.")

	return df_norm


@numba.jit(nopython=True)
def read_score_tumor(
	reads,
	total_samples,
	upper_bound
):
	n = reads.shape[0]
	out = np.empty(n, dtype=np.float64)
	for i in range(n):
		total_upper_bound = upper_bound[i] * total_samples[i]
		if total_upper_bound <= 0:
			out[i] = 0.0
		elif reads[i] <= 2.0 * total_upper_bound:
			out[i] = (total_upper_bound - abs(reads[i] - total_upper_bound)) / total_upper_bound
		else:
			out[i] = 0.0
	return out


def read_score_tumor_vectorized(reads, total_samples_const, upper_bound_const):
	"""
	Vectorized NumPy version of `read_score_tumor`.

	Parameters
	- reads: 1-D array-like of read counts (column).
	- total_samples_const: scalar used as the total_samples normalization constant for the column.
	- upper_bound_const: scalar upper bound used for normalization for the column.

	Returns
	- NumPy array of scores with the same shape as `reads`.

	This operates column-wise (broadcasting the two constants) and avoids Python loops.
	"""
	total_upper_bound = upper_bound_const * total_samples_const
	out = np.zeros_like(reads, dtype=np.float64)
	if total_upper_bound <= 0:
		return out
	# if reads in [0, 2 * upper bound], then score is (upper_bound - abs(reads - upper_bound)) / upper_bound, else score is 0.0
	# the former expression is flipped, shifted, and scaled absolute value function peaking at the upper_bound
	mask = reads <= 2.0 * total_upper_bound
	out[mask] = (total_upper_bound - np.abs(reads[mask] - total_upper_bound)) / total_upper_bound
	return out

@numba.jit(nopython=True)
def read_score_normal(
    reads,
    total_samples,
    upper_bound
):
	n = reads.shape[0]
	out = np.empty(n, dtype=np.float64)
	for i in range(n):
		total_upper_bound = upper_bound[i] * total_samples[i]
		if total_upper_bound <= 0:
			out[i] = 0.0
		elif reads[i] <= total_upper_bound:
			out[i] = reads[i] / total_upper_bound
		else:
			out[i] = 1.0
	return out

def read_score_normal_vectorized(reads, total_samples_const, upper_bound_const):
	"""
	Vectorized NumPy version of `read_score_normal`.

	Parameters
	- reads: 1-D array-like of read counts (column).
	- total_samples_const: scalar used as the total_samples normalization constant for the column.
	- upper_bound_const: scalar upper bound used for normalization for the column.

	Returns
	- NumPy array of scores with the same shape as `reads`.

	This operates column-wise (broadcasting the two constants) and avoids Python loops.
	"""
	total_upper_bound = upper_bound_const * total_samples_const
	out = np.zeros_like(reads, dtype=np.float64)
	if total_upper_bound <= 0:
		return out
	# if reads <= upper bound, then score is reads / upper bound, else score is 1.0
	mask = reads <= total_upper_bound
	out[mask] = reads[mask] / total_upper_bound
	out[~mask] = 1.0
	return out

@numba.jit(nopython=True)
def sample_score(
	samples,
	total_samples
):
	n = samples.shape[0]
	out = np.empty(n, dtype=np.float64)
	for i in range(n):
		if total_samples[i] <= 0:
			out[i] = 0.0
		else:
			out[i] = samples[i] / total_samples[i]
	return out

def sample_score_vectorized(samples, total_samples_const):
	"""
	Vectorized NumPy version of `sample_score`.

	Parameters
	- samples: 1-D array-like of sample counts (column).
	- total_samples_const: scalar used as the total_samples normalization constant for the column.

	Returns
	- NumPy array of scores with the same shape as `samples`.

	This operates column-wise (broadcasting the total_samples_const) and avoids Python loops.
	"""
	out = np.zeros_like(samples, dtype=np.float64)
	mask = total_samples_const > 0
	out[mask] = samples[mask] / total_samples_const
	return out




    
	


# Numba JIT-compiled version
@numba.jit(nopython=True)
def score_numba(
	# columns are num_supporting_reads, num_supporting_samples, total_samples, and upper_bound
	tumor_matrix, # T x 4
	normal_matrix, # N x 4
	w_normal=0.5
):
	w_tumor = 1.0 - w_normal
	# total the sample counts of sub-(tumor/normal) populations for weighted averaging
	total_samples_tumor = np.sum(tumor_matrix[:, 2])  # scalar
	total_samples_normal = np.sum(normal_matrix[:, 2])  # scalar
	### tumor
	## read
	tumor_read_scores = read_score_tumor(
    	tumor_matrix[:, 0],
		tumor_matrix[:, 2],
		tumor_matrix[:, 3]
	) # T x 1
	## sample
	tumor_sample_scores = sample_score(
		tumor_matrix[:, 1],
		tumor_matrix[:, 2]
	) # T x 1
	## combine read and sample scores
	tumor_score = tumor_read_scores + tumor_sample_scores  # T x 1
	## normalize and weighted average
	# div 2 bc read and sample score are both in [0,1], normalized score is 1
	# multiply by fraction sub-population/population for weighting
	if total_samples_tumor > 0:
		normalization_factors_tumor = w_tumor * 0.5 * (tumor_matrix[:, 2] / total_samples_tumor)  # T x 1
	else:
		normalization_factors_tumor = np.zeros(tumor_matrix.shape[0], dtype=np.float64)
	# elementwise multiply normalization term and read+sample scores
	tumor_score_weighted = tumor_score * normalization_factors_tumor  # T x 1
	# finally collapse into weighted average (weights are sub-population/population fractions)
	tumor_score_final = np.sum(tumor_score_weighted)
	# normal
	normal_read_scores = read_score_normal(
		normal_matrix[:, 0],
		normal_matrix[:, 2],
		normal_matrix[:, 3]
	) # N x 1
	normal_sample_scores = sample_score(
		normal_matrix[:, 1],
		normal_matrix[:, 2]
	) # N x 1
	normal_score = normal_read_scores + normal_sample_scores  # N x 1
	if total_samples_normal > 0:
		normalization_factors_normal = w_normal * 0.5 * (normal_matrix[:, 2] / total_samples_normal)  # N x 1
	else:
		normalization_factors_normal = np.zeros(normal_matrix.shape[0], dtype=np.float64)
	normal_score_weighted = normal_score * normalization_factors_normal  # N x 1
	normal_score_final = np.sum(normal_score_weighted)
	final_score = tumor_score_final - normal_score_final
	return final_score


@numba.jit(nopython=True, parallel=True)
def score_numba_batched(
	tumor_batch,   # M x T x 4
	normal_batch,  # M x N x 4
	w_normal=0.5
):
	m = tumor_batch.shape[0]
	out = np.empty(m, dtype=np.float64)
	for i in numba.prange(m):
		out[i] = score_numba(tumor_batch[i], normal_batch[i], w_normal=w_normal)
	return out



