import warnings, platform
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, sklearn
from sklearn.pipeline import Pipeline
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import RepeatedStratifiedKFold, cross_val_score

print("python", platform.python_version(), "| sklearn", sklearn.__version__,
      "| numpy", np.__version__, "| pandas", pd.__version__)
meta = pd.read_csv("data/ml_meta.tsv", sep="\t")
X = pd.read_csv("data/ml_expr.tsv", sep="\t", index_col=0)
assert list(X.index) == list(meta["sample"]), "sample order mismatch"
y = meta["label"].to_numpy(); Xv = X.to_numpy(dtype=float)
C = meta[["age", "sex_male", "pmi", "ph", "rin", "batch2", "batch3"]].to_numpy(dtype=float)
n = Xv.shape[1]
print("samples", Xv.shape[0], "| probes", n, "| AD", int(y.sum()), "| control", int((y == 0).sum()), flush=True)

K = 200
cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=1)
def rf():
    return RandomForestClassifier(n_estimators=500, max_depth=10, min_samples_split=5,
                                  class_weight="balanced", n_jobs=-1, random_state=1)
def expr_block():
    return Pipeline([("imp", SimpleImputer(strategy="median")), ("sel", SelectKBest(f_classif, k=K))])

ct = ColumnTransformer([("expr", expr_block(), slice(0, n)),
                        ("cov", "passthrough", slice(n, n + C.shape[1]))])
models = {
 "A_expression_only (selection inside folds)": (Xv, Pipeline([("sel", expr_block()), ("rf", rf())])),
 "B_covariates_only (logistic regression)": (C, Pipeline([("sc", StandardScaler()), ("lr", LogisticRegression(max_iter=1000))])),
 "C_expression+covariates (selection inside folds)": (np.hstack([Xv, C]), Pipeline([("ct", ct), ("rf", rf())])),
}
rows = []
def report(name, s):
    rows.extend((name, i, float(v)) for i, v in enumerate(s))
    print(f"{name}: mean AUC {s.mean():.3f} | SD across 50 folds {s.std():.3f} | range {s.min():.3f}-{s.max():.3f}", flush=True)

for name, (M, pipe) in models.items():
    report(name, cross_val_score(pipe, M, y, cv=cv, scoring="roc_auc"))

# Leaky reference: probes chosen on ALL samples before CV (the original design). Do not report as a result.
Xleak = SelectKBest(f_classif, k=K).fit_transform(SimpleImputer(strategy="median").fit_transform(Xv), y)
report("LEAKY_REFERENCE (selection before CV; do not report)", cross_val_score(rf(), Xleak, y, cv=cv, scoring="roc_auc"))

pd.DataFrame(rows, columns=["model", "fold", "auc"]).to_csv("results/cv_fold_aucs.tsv", sep="\t", index=False)
