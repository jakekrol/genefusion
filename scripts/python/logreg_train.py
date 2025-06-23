#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, roc_auc_score, average_precision_score
import statsmodels.api as sm
import pickle 

def parse_arguments():
    parser = argparse.ArgumentParser(description="Logistic Regression Training Script")
    parser.add_argument("-i", "--input", type=str, required=True, help="Path to the input tabular file (CSV/TSV).")
    parser.add_argument("-y", "--label-idx", type=int, required=True, help="Index of the label column (1-indexed).")
    parser.add_argument("-r", "--regularization", type=str, choices=["l1", "l2", "elasticnet"], default="l2",
                        help="Type of regularization to use (default: l2).")
    parser.add_argument("--alpha", type=float, default=1.0, help="Regularization strength (default: 1.0).")
    parser.add_argument("--l1-ratio", type=float, default=0.5,
                        help="ElasticNet mixing parameter (only used if regularization is elasticnet).")
    parser.add_argument("-C", "--regularization-strength", type=float, default=1.0,
                        help="Inverse of regularization strength (C) (default: 1.0).")
    parser.add_argument("-t", "--type", type=str, choices=["model", "stats"], required=True,
                        help="Specify the output type: 'model' to save the model or 'stats' to save the statistics.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output file path for the model or statistics.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility (default: 42).")
    return parser.parse_args()

def load_data(input_path, label_idx):
    data = pd.read_csv(input_path, sep='\t')
    
    # Check for NA values
    if data.isna().any().any():
        na_rows = data[data.isna().any(axis=1)].index.tolist()
        na_columns = data.columns[data.isna().any()].tolist()
        print(f"Warning: The dataset contains missing values.")
        print(f"Rows with NA values: {na_rows}")
        print(f"Columns with NA values: {na_columns}")
    
    y = data.iloc[:, label_idx]
    X = data.drop(data.columns[label_idx], axis=1)
    return X, y

def train_logistic_regression(X, y, regularization, alpha, l1_ratio, regularization_strength,seed):
    if regularization == "elasticnet":
        solver = "saga"
        penalty = "elasticnet"
    elif regularization == "l1":
        solver = "saga"
        penalty = "l1"
    else:  # l2
        solver = "lbfgs"
        penalty = "l2"

    model = LogisticRegression(
        penalty=penalty,
        C=regularization_strength,
        solver=solver,
        l1_ratio=l1_ratio if regularization == "elasticnet" else None,
        max_iter=1000,
        random_state=seed
    )
    model.fit(X, y)
    return model

def calculate_coefficient_statistics(X, y, model):
    X_with_const = sm.add_constant(X)
    logit_model = sm.Logit(y, X_with_const)
    result = logit_model.fit(disp=False)
    coef_stats = pd.DataFrame({
        "Coefficient": result.params,
        "P-Value": result.pvalues
    })
    return coef_stats

def evaluate_model(X, y, model):
    y_pred = model.predict(X)
    y_pred_proba = model.predict_proba(X)[:, 1]
    accuracy = accuracy_score(y, y_pred)
    auroc = roc_auc_score(y, y_pred_proba)
    auprc = average_precision_score(y, y_pred_proba)
    return accuracy, auroc, auprc

def main():
    args = parse_arguments()
    X, y = load_data(args.input, args.label_idx-1)
    print("Label column: ", y.name)
    print(y.head())
    print("Feature columns: ", X.columns.tolist())
    print(X.head())
    
    model = train_logistic_regression(X, y, args.regularization, args.alpha, args.l1_ratio, args.regularization_strength, args.seed)
    # print the model coefficients and bias
    print("\nModel Coefficients:")
    for feature, coef in zip(X.columns, model.coef_[0]):
        print(f"{feature}: {coef:.4f}")
    print(f"Bias: {model.intercept_[0]:.4f}")
    
    if args.type == "model":
        # Save the trained model to the specified file
        with open(args.output, "wb") as model_file:
            pickle.dump(model, model_file)
        print(f"\nTrained model saved to '{args.output}'.")
    elif args.type == "stats":
        coef_stats = calculate_coefficient_statistics(X, y, model)
        coef_stats.to_csv(args.output, index=False,sep='\t')
        print(f"\nCoefficient statistics saved to '{args.output}'.")
    
    # Evaluate the model and print performance metrics
    accuracy, auroc, auprc = evaluate_model(X, y, model)
    print("\nPerformance Metrics:")
    print(f"Accuracy: {accuracy:.4f}")
    print(f"AUROC: {auroc:.4f}")
    print(f"AUPRC: {auprc:.4f}")

if __name__ == "__main__":
    main()