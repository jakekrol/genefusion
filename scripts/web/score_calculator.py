#!/usr/bin/env python3
import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from genefusion.genefusion import score_numba

st.set_page_config(page_title="Gene Fusion Score Calculator", layout="wide")

st.title("Gene Fusion Score Calculator")

# Initialize session state for parameters
if 'params' not in st.session_state:
    st.session_state.params = {
        "reads_dna_normal": 0,
        "reads_dna_tumor": 5000,
        "reads_rna_normal": 0,
        "reads_rna_tumor": 0,
        "reads_onekg": 0,
        "samples_dna_normal": 0,
        "samples_dna_tumor": 100,
        "samples_rna_normal": 0,
        "samples_rna_tumor": 0,
        "samples_onekg": 0,
        "pop_size_dna_normal": 0,
        "pop_size_dna_tumor": 100,
        "pop_size_rna_normal": 0,
        "pop_size_rna_tumor": 0,
        "pop_size_dna_onekg": 0,
        "w_tumor": 1.0,
        "w_dna": 1.0,
        "w_read": 0.5,
        "cov_dna_tumor": 50.0,
        "cov_rna_tumor": 0.0,
        "cov_dna_normal": 0.0,
        "cov_rna_normal": 0.0,
        "cov_onekg": 0.0
    }

# Input Section
st.header("Input Parameters")

col1, col2, col3 = st.columns(3)

with col1:
    with st.expander("📊 Data (Reads & Samples)", expanded=True):
        st.subheader("Reads")
        reads_dna_normal = st.number_input("reads_dna_normal", min_value=0, value=st.session_state.params["reads_dna_normal"], step=1)
        reads_dna_tumor = st.number_input("reads_dna_tumor", min_value=0, value=st.session_state.params["reads_dna_tumor"], step=1)
        reads_rna_normal = st.number_input("reads_rna_normal", min_value=0, value=st.session_state.params["reads_rna_normal"], step=1)
        reads_rna_tumor = st.number_input("reads_rna_tumor", min_value=0, value=st.session_state.params["reads_rna_tumor"], step=1)
        reads_onekg = st.number_input("reads_onekg", min_value=0, value=st.session_state.params["reads_onekg"], step=1)
        
        st.subheader("Samples")
        samples_dna_normal = st.number_input("samples_dna_normal", min_value=0, value=st.session_state.params["samples_dna_normal"], step=1)
        samples_dna_tumor = st.number_input("samples_dna_tumor", min_value=0, value=st.session_state.params["samples_dna_tumor"], step=1)
        samples_rna_normal = st.number_input("samples_rna_normal", min_value=0, value=st.session_state.params["samples_rna_normal"], step=1)
        samples_rna_tumor = st.number_input("samples_rna_tumor", min_value=0, value=st.session_state.params["samples_rna_tumor"], step=1)
        samples_onekg = st.number_input("samples_onekg", min_value=0, value=st.session_state.params["samples_onekg"], step=1)

with col2:
    with st.expander("⚙️ Parameters (Coverage & Pop Size)", expanded=True):
        st.subheader("Coverage")
        cov_dna_tumor = st.number_input("cov_dna_tumor", min_value=0.0, value=st.session_state.params["cov_dna_tumor"], step=1.0)
        cov_rna_tumor = st.number_input("cov_rna_tumor", min_value=0.0, value=st.session_state.params["cov_rna_tumor"], step=1.0)
        cov_dna_normal = st.number_input("cov_dna_normal", min_value=0.0, value=st.session_state.params["cov_dna_normal"], step=1.0)
        cov_rna_normal = st.number_input("cov_rna_normal", min_value=0.0, value=st.session_state.params["cov_rna_normal"], step=1.0)
        cov_onekg = st.number_input("cov_onekg", min_value=0.0, value=st.session_state.params["cov_onekg"], step=1.0)
        
        st.subheader("Population Size")
        pop_size_dna_normal = st.number_input("pop_size_dna_normal", min_value=0, value=st.session_state.params["pop_size_dna_normal"], step=1)
        pop_size_dna_tumor = st.number_input("pop_size_dna_tumor", min_value=0, value=st.session_state.params["pop_size_dna_tumor"], step=1)
        pop_size_rna_normal = st.number_input("pop_size_rna_normal", min_value=0, value=st.session_state.params["pop_size_rna_normal"], step=1)
        pop_size_rna_tumor = st.number_input("pop_size_rna_tumor", min_value=0, value=st.session_state.params["pop_size_rna_tumor"], step=1)
        pop_size_dna_onekg = st.number_input("pop_size_dna_onekg", min_value=0, value=st.session_state.params["pop_size_dna_onekg"], step=1)

with col3:
    with st.expander("🎛️ Hyperparameters (Weights)", expanded=True):
        w_tumor = st.slider("w_tumor", min_value=0.0, max_value=1.0, value=st.session_state.params["w_tumor"], step=0.01)
        w_dna = st.slider("w_dna", min_value=0.0, max_value=1.0, value=st.session_state.params["w_dna"], step=0.01)
        w_read = st.slider("w_read", min_value=0.0, max_value=1.0, value=st.session_state.params["w_read"], step=0.01)

# Update session state
st.session_state.params.update({
    "reads_dna_normal": reads_dna_normal,
    "reads_dna_tumor": reads_dna_tumor,
    "reads_rna_normal": reads_rna_normal,
    "reads_rna_tumor": reads_rna_tumor,
    "reads_onekg": reads_onekg,
    "samples_dna_normal": samples_dna_normal,
    "samples_dna_tumor": samples_dna_tumor,
    "samples_rna_normal": samples_rna_normal,
    "samples_rna_tumor": samples_rna_tumor,
    "samples_onekg": samples_onekg,
    "pop_size_dna_normal": pop_size_dna_normal,
    "pop_size_dna_tumor": pop_size_dna_tumor,
    "pop_size_rna_normal": pop_size_rna_normal,
    "pop_size_rna_tumor": pop_size_rna_tumor,
    "pop_size_dna_onekg": pop_size_dna_onekg,
    "w_tumor": w_tumor,
    "w_dna": w_dna,
    "w_read": w_read,
    "cov_dna_tumor": cov_dna_tumor,
    "cov_rna_tumor": cov_rna_tumor,
    "cov_dna_normal": cov_dna_normal,
    "cov_rna_normal": cov_rna_normal,
    "cov_onekg": cov_onekg
})

st.divider()

# Compute Button
if st.button("🔬 Compute Score", type="primary", use_container_width=True):
    with st.spinner("Computing score..."):
        score = score_numba(**st.session_state.params)
        st.session_state.current_score = score

# Output Section
st.header("Output")
if 'current_score' in st.session_state:
    st.success(f"**Score: {st.session_state.current_score:.6f}**")
else:
    st.info("Click 'Compute Score' to calculate the score")

st.divider()

# Sensitivity Analysis Section
st.header("Sensitivity Analysis")

col_graph1, col_graph2 = st.columns([1, 2])

with col_graph1:
    st.subheader("Parameter Selection")
    
    # Create parameter groups
    param_groups = {
        "Data - Reads": ["reads_dna_normal", "reads_dna_tumor", "reads_rna_normal", "reads_rna_tumor", "reads_onekg"],
        "Data - Samples": ["samples_dna_normal", "samples_dna_tumor", "samples_rna_normal", "samples_rna_tumor", "samples_onekg"],
        "Parameters - Coverage": ["cov_dna_tumor", "cov_rna_tumor", "cov_dna_normal", "cov_rna_normal", "cov_onekg"],
        "Parameters - Pop Size": ["pop_size_dna_normal", "pop_size_dna_tumor", "pop_size_rna_normal", "pop_size_rna_tumor", "pop_size_dna_onekg"],
        "Hyperparameters": ["w_tumor", "w_dna", "w_read"]
    }
    
    # Flatten for selection
    all_params = []
    for group, params in param_groups.items():
        all_params.extend(params)
    
    selected_param = st.selectbox("Select parameter to vary", all_params)
    
    # Get current value
    current_val = st.session_state.params[selected_param]
    
    # Determine if integer or float
    is_weight = selected_param in ["w_tumor", "w_dna", "w_read"]
    is_coverage = "cov_" in selected_param
    is_float = is_weight or is_coverage
    
    if is_weight:
        min_val = st.number_input("Min value", min_value=0.0, max_value=1.0, value=0.0, step=0.01)
        max_val = st.number_input("Max value", min_value=0.0, max_value=1.0, value=1.0, step=0.01)
        num_points = st.slider("Number of points", min_value=10, max_value=100, value=50)
    elif is_float:
        min_val = st.number_input("Min value", min_value=0.0, value=0.0, step=1.0)
        max_val = st.number_input("Max value", min_value=0.0, value=float(max(100, current_val * 2)), step=1.0)
        num_points = st.slider("Number of points", min_value=10, max_value=100, value=50)
    else:
        min_val = st.number_input("Min value", min_value=0, value=0, step=1)
        max_val = st.number_input("Max value", min_value=0, value=max(100, int(current_val * 2)), step=1)
        num_points = st.slider("Number of points", min_value=10, max_value=100, value=50)
    
    compute_graph = st.button("📈 Generate Graph", use_container_width=True)

with col_graph2:
    st.subheader("Score vs Parameter")
    
    if compute_graph and min_val < max_val:
        with st.spinner("Generating sensitivity analysis..."):
            # Create range of values
            if is_float:
                param_values = np.linspace(min_val, max_val, num_points)
            else:
                param_values = np.linspace(min_val, max_val, num_points, dtype=int)
            
            scores = []
            
            # Calculate scores
            for val in param_values:
                temp_params = st.session_state.params.copy()
                temp_params[selected_param] = float(val) if is_float else int(val)
                score = score_numba(**temp_params)
                scores.append(score)
            
            # Create plot
            fig = go.Figure()
            fig.add_trace(go.Scatter(
                x=param_values,
                y=scores,
                mode='lines+markers',
                name='Score',
                line=dict(color='#1f77b4', width=2),
                marker=dict(size=4)
            ))
            
            # Add vertical line at current value
            fig.add_vline(
                x=current_val,
                line_dash="dash",
                line_color="red",
                annotation_text=f"Current: {current_val}",
                annotation_position="top"
            )
            
            fig.update_layout(
                xaxis_title=selected_param,
                yaxis_title="Score",
                hovermode='x unified',
                height=500
            )
            
            st.plotly_chart(fig, use_container_width=True)
            
            # Show statistics
            st.metric("Maximum Score", f"{max(scores):.6f}", f"at {selected_param}={param_values[np.argmax(scores)]:.2f}")
            st.metric("Minimum Score", f"{min(scores):.6f}", f"at {selected_param}={param_values[np.argmin(scores)]:.2f}")
    elif compute_graph:
        st.warning("Min value must be less than Max value")
    else:
        st.info("Configure parameters and click 'Generate Graph' to see sensitivity analysis")
