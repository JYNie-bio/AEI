from sklearn.model_selection import GridSearchCV, KFold
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn import metrics
from xgboost import XGBClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
#from sklearn.inspection import permutation_importance
import math
import numpy as np
import pandas as pd
import os
import sys
from sklearn.metrics import make_scorer, accuracy_score, roc_auc_score, f1_score, recall_score, precision_score, cohen_kappa_score, balanced_accuracy_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

arguments = sys.argv
new_working_directory = "/lustre/home/niejingyi2023/project/bjCDC/ARG_MIC/"
os.chdir(new_working_directory)


traindata=pd.read_csv("RS_811/14880file_ARGanno_Protein_result_copynumber_addPhenotype_" + arguments[1] + "_train.txt",sep='\t', header=0,index_col=0)
verify_testdata=pd.read_csv("RS_811/14880file_ARGanno_Protein_result_copynumber_addPhenotype_" + arguments[1] + "_verify_test.txt",sep='\t', header=0,index_col=0)
X_train=traindata.iloc[:, :-1]
y_train=traindata.iloc[:, -1]
y_train=y_train.replace({'R': 1, 'S': 0})
X_verify_test=verify_testdata.iloc[:, :-1]
y_verify_test=verify_testdata.iloc[:, -1]
y_verify_test=y_verify_test.replace({'R': 1, 'S': 0})


skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=1)

def calculate_metrics(y_true, y_pred, y_pred_proba):
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    metrics = {
        'Precision': tp / (tp + fp) if (tp + fp) != 0 else 0,
        'Recall': tp / (tp + fn) if (tp + fn) != 0 else 0,
        'Specificity': tn / (fp + tn) if (fp + tn) != 0 else 0,
        'F1_Score': f1_score(y_true, y_pred),
        'Accuracy': (tp + tn) / (tp + fn + fp + tn) if (tp + fn + fp + tn) != 0 else 0,
        'Balanced_Accuracy': balanced_accuracy_score(y_true, y_pred),
        'AUC': roc_auc_score(y_true, y_pred_proba),
        'Major_Error': fp / (fp + tn) if (fp + tn) != 0 else 0,
        'Very_Major_Error': fn / (tp + fn) if (tp + fn) != 0 else 0,
        'NPV': tn / (tn + fn) if (tn + fn) != 0 else 0,
        'PPV': tp / (tp + fp) if (tp + fp) != 0 else 0,
        'Prevalence': np.mean(y_true)
    }
    return metrics


all_metrics = []

for fold, (train_idx, val_idx) in enumerate(skf.split(X_train, y_train)):
    X_fold_train, X_fold_val = X_train.iloc[train_idx], X_train.iloc[val_idx]
    y_fold_train, y_fold_val = y_train.iloc[train_idx], y_train.iloc[val_idx]    
    fold_train_df = pd.concat([X_fold_train, y_fold_train], axis=1)
    fold_val_df = pd.concat([X_fold_val, y_fold_val], axis=1)
    fold_train_df.to_csv(f"AllFeature/XGBoost_ARG_copynumber_fold_{fold+1}_train_samples_10fold_"+ arguments[1] +".txt", sep='\t', index=True, header=True)
    fold_val_df.to_csv(f"AllFeature/XGBoost_ARG_copynumber_fold_{fold+1}_val_samples_10fold_"+ arguments[1] +".txt", sep='\t', index=True, header=True)

    model = XGBClassifier(random_state=2)
    model.fit(X_fold_train, y_fold_train)
    
    y_val_pred = model.predict(X_fold_val)
    y_verify_test_pred = model.predict(X_verify_test)
    
    y_val_pred_proba = model.predict_proba(X_fold_val)[:, 1]
    y_verify_test_pred_proba = model.predict_proba(X_verify_test)[:, 1]

    val_metrics = calculate_metrics(y_fold_val, y_val_pred, y_val_pred_proba)
    verify_test_metrics = calculate_metrics(y_verify_test, y_verify_test_pred, y_verify_test_pred_proba)
    
    fold_results = pd.DataFrame([val_metrics, verify_test_metrics], 
                                index=[f'XGBoost_Fold_{fold+1}_val', f'XGBoost_Fold_{fold+1}_verify_test'])
    
    all_metrics.append(fold_results)

all_metrics_df = pd.concat(all_metrics)
all_metrics_df.to_csv("AllFeature/XGBoost_10fold_crossvalidation_default_hyperparameter_ARG_copynumber_"+ arguments[1] +".txt", sep='\t', index=True, header=True)


