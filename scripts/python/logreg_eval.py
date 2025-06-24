#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve, precision_recall_curve
from sklearn.preprocessing import StandardScaler

def parse_args():
    parser = argparse.ArgumentParser(description="Thorough logistic regression evaluation with ElasticNet")
    parser.add_argument("-i", "--input", required=True, help="Input TSV file")
    parser.add_argument("-y", "--label_col", type=int, required=True, help="1-indexed label column (subtracts 1 internally)")
    parser.add_argument("-k", "--kfolds", type=int, default=5, help="Number of CV folds (default: 5)")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file with predicted probabilities")
    parser.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    return parser.parse_args()

def random_hyperparams(seed):
    rng = np.random.RandomState(seed)
    # Always include bounds
    C_vals = [0.0001, 100]
    l1_vals = [0.0, 1.0]
    # 50 random samples
    C_rand = rng.uniform(np.log10(0.0001), np.log10(100), 50) # sample from log space to ensure wide range
    C_rand = np.power(10, C_rand) # convert back to linear scale
    l1_rand = rng.uniform(0, 1, 50)
    C_vals += list(C_rand)
    l1_vals += list(l1_rand)
    return list(zip(C_vals, l1_vals))

def plot_hp_heatmap(C_list, l1_list, auc_list, outpath):
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.colors import LogNorm

    plt.figure(figsize=(6,5))
    sc = plt.scatter(l1_list, C_list, c=auc_list, cmap='viridis', s=60, edgecolor='k', norm=None)
    plt.yscale('log')
    plt.xlabel("l1_ratio")
    plt.ylabel("C")
    plt.title("AUROC for each (C, l1_ratio)")
    cbar = plt.colorbar(sc)
    cbar.set_label("AUROC")
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()

def plot_roc(y_true, y_pred, outpath, title):
    fpr, tpr, _ = roc_curve(y_true, y_pred)
    auc = roc_auc_score(y_true, y_pred)
    plt.figure()
    plt.plot(fpr, tpr, label=f"AUROC = {auc:.4f}")
    plt.plot([0,1],[0,1],'k--')
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()

def plot_pr(y_true, y_pred, outpath, title):
    prec, rec, _ = precision_recall_curve(y_true, y_pred)
    auc = average_precision_score(y_true, y_pred)
    plt.figure()
    plt.plot(rec, prec, label=f"AUPRC = {auc:.4f}")
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()

def plot_coefs(coefs, feat_names, outpath, title):
    plt.figure(figsize=(max(6, len(feat_names)*0.7), 4))
    plt.bar(feat_names, coefs)
    plt.xticks(rotation=45, ha='right')
    plt.ylabel("Coefficient")
    plt.hlines(0, xmin=0, xmax=len(feat_names)-1, colors='k', linestyles='dashed')
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()

def plot_pred_prob_histograms(df, outpath):
    import matplotlib.pyplot as plt

    plt.figure(figsize=(6,4))
    bins = np.linspace(0, 1, 30)
    plt.hist(df[df['label'] == 1]['pred_prob'], bins=bins, color='red', alpha=0.5, label='label=1', density=True)
    plt.hist(df[df['label'] == 0]['pred_prob'], bins=bins, color='blue', alpha=0.5, label='label=0', density=True)
    plt.xlabel("Predicted Probability")
    plt.ylabel("Density")
    plt.title("Predicted Probability Distributions")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()

def main():
    args = parse_args()
    np.random.seed(args.seed)

    # Load data
    df = pd.read_csv(args.input, sep="\t")
    y = df.iloc[:, args.label_col - 1].values
    X = df.drop(df.columns[args.label_col - 1], axis=1).values
    feat_names = df.drop(df.columns[args.label_col - 1], axis=1).columns.tolist()

    # Standardize features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Hyperparameter search
    hp_list = random_hyperparams(args.seed)
    best_auc = -np.inf
    best_hp = None
    best_model = None
    best_cv_pred = None
    best_cv_true = None

    hp_Cs, hp_l1s, hp_aucs = [], [], []
    print("Searching hyperparameters...")
    for idx, (C, l1) in enumerate(hp_list):
        cv = StratifiedKFold(n_splits=args.kfolds, shuffle=True, random_state=args.seed)
        cv_pred = []
        cv_true = []
        for train_idx, test_idx in cv.split(X_scaled, y):
            model = LogisticRegression(
                penalty="elasticnet", solver="saga", l1_ratio=l1, C=C, max_iter=10000, random_state=args.seed
            )
            model.fit(X_scaled[train_idx], y[train_idx])
            prob = model.predict_proba(X_scaled[test_idx])[:,1]
            cv_pred.append(prob)
            cv_true.append(y[test_idx])
        cv_pred = np.concatenate(cv_pred)
        cv_true = np.concatenate(cv_true)
        auc = roc_auc_score(cv_true, cv_pred)
        hp_Cs.append(C)
        hp_l1s.append(l1)
        hp_aucs.append(auc)
        if auc > best_auc:
            best_auc = auc
            best_hp = (C, l1)
            best_cv_pred = cv_pred
            best_cv_true = cv_true
        print(f"HP {idx+1}/{len(hp_list)}: C={C:.5g}, l1_ratio={l1:.3f}, AUROC={auc:.4f}")

    print(f"\nBest hyperparameters: C={best_hp[0]:.5g}, l1_ratio={best_hp[1]:.3f}, CV AUROC={best_auc:.4f}")

    # Final model on full data
    final_model = LogisticRegression(
        penalty="elasticnet", solver="saga", l1_ratio=best_hp[1], C=best_hp[0], max_iter=10000, random_state=args.seed
    )
    final_model.fit(X_scaled, y)
    full_pred = final_model.predict_proba(X_scaled)[:,1]
    full_auc = roc_auc_score(y, full_pred)
    full_auprc = average_precision_score(y, full_pred)
    print(f"Full-data AUROC: {full_auc:.4f}")
    print(f"Full-data AUPRC: {full_auprc:.4f}")

    # Plots
    plot_roc(y, full_pred, args.output + ".full_roc.png", "Full Data ROC")
    plot_pr(y, full_pred, args.output + ".full_pr.png", "Full Data PR")
    plot_roc(best_cv_true, best_cv_pred, args.output + ".cv_roc.png", "CV ROC (best HP)")
    plot_pr(best_cv_true, best_cv_pred, args.output + ".cv_pr.png", "CV PR (best HP)")
    plot_hp_heatmap(hp_Cs, hp_l1s, hp_aucs, args.output + ".hp_heatmap.png")

    # Print and save AUROC/AUPRC
    cv_auprc = average_precision_score(best_cv_true, best_cv_pred)
    print(f"CV AUROC (best HP): {best_auc:.4f}")
    print(f"CV AUPRC (best HP): {cv_auprc:.4f}")

    # Save D* with predicted probabilities
    df_out = df.copy()
    df_out["pred_prob"] = full_pred
    df_out.to_csv(args.output, sep="\t", index=False)

    # Barplot of coefficients
    plot_coefs(final_model.coef_[0], feat_names, args.output + ".coefs.png", "Final Model Coefficients")

    # Overlayed histogram of predicted probabilities
    plot_pred_prob_histograms(df_out, args.output + ".pred_prob_hist.png")

    # Output summary table
    summary = pd.DataFrame([{
        "full_data_AUROC": full_auc,
        "full_data_AUPRC": full_auprc,
        "cv_best_AUROC": best_auc,
        "cv_best_AUPRC": cv_auprc,
        "best_C": best_hp[0],
        "best_l1_ratio": best_hp[1]
    }])
    summary.to_csv(args.output + ".summary.tsv", sep="\t", index=False)

main()
