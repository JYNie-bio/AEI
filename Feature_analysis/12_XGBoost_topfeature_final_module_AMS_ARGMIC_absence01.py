from sklearn.metrics import confusion_matrix
from sklearn import metrics
from xgboost import XGBClassifier
import numpy as np
import pandas as pd
import os
import sys
import re
from sklearn.metrics import make_scorer, accuracy_score, roc_auc_score, f1_score, recall_score, precision_score, cohen_kappa_score, balanced_accuracy_score

arguments = sys.argv
new_working_directory = "~/Feature_analysis/"
os.chdir(new_working_directory)


traindata=pd.read_csv("~/ARG_MIC/RS_811/14880file_ARGanno_Protein_result_absence01_addPhenotype_AMS_train.txt",sep='\t', header=0,index_col=0)
X_train=traindata.iloc[:, :-1]
y_train=traindata.iloc[:, -1]
y_train=y_train.replace({'R': 1, 'S': 0})

verifytestdata=pd.read_csv("~/ARG_MIC/RS_811/14880file_ARGanno_Protein_result_absence01_addPhenotype_AMS_verify_test.txt",sep='\t', header=0,index_col=0)
X_verify_test=verifytestdata.iloc[:, :-1]
y_verify_test=verifytestdata.iloc[:, -1]
y_verify_test=y_verify_test.replace({'R': 1, 'S': 0})


df_sorted=pd.read_csv("~/ARG_MIC/AllFeature/XGBoost_10fold_crossvalidation_default_hyperparameter_ARG_absence01_AMS_fold_6_feature_importances_verifytest.txt",sep='\t', header=0,index_col=0)
df_filtered = df_sorted[df_sorted['Importance'] > 0]
df_filtered_top_index0 = df_filtered.head(20).index
X_train2=X_train.iloc[:,df_filtered_top_index0]
X_verify_test2=X_verify_test.iloc[:,df_filtered_top_index0]

final_selected_features_index = X_train2.columns.tolist() + ['Group']
traindata2=traindata.loc[:,final_selected_features_index]
verifytestdata2=verifytestdata.loc[:,final_selected_features_index]

old_col_pattern = "cryptic_aminoglycoside_nucleotidyltransferase_ANT\(3\''\)/ANT\(9\)"
traindata2.rename(columns={col: re.sub(old_col_pattern, "ANT(3'')-Ia_family_aminoglycoside_nucleotidyltransferase_AadA", col)
                              for col in traindata2.columns}, inplace=True)
verifytestdata2.rename(columns={col: re.sub(old_col_pattern, "ANT(3'')-Ia_family_aminoglycoside_nucleotidyltransferase_AadA", col)
                              for col in verifytestdata2.columns}, inplace=True)

traindata2.to_csv("~/Feature_analysis/use_topfeature/XGBoost_AMS_topfeature_finalmodel_traindata.txt", sep='\t', index=True, header=True)
verifytestdata2.to_csv("~/Feature_analysis/use_topfeature/XGBoost_AMS_topfeature_finalmodel_verifytestdata.txt", sep='\t', index=True, header=True)



old_col_pattern = "cryptic_aminoglycoside_nucleotidyltransferase_ANT\(3\''\)/ANT\(9\)"
X_train2.rename(columns={col: re.sub(old_col_pattern, "ANT(3'')-Ia_family_aminoglycoside_nucleotidyltransferase_AadA", col)
                              for col in X_train2.columns}, inplace=True)
X_verify_test2.rename(columns={col: re.sub(old_col_pattern, "ANT(3'')-Ia_family_aminoglycoside_nucleotidyltransferase_AadA", col)
                              for col in X_verify_test2.columns}, inplace=True)


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


model = XGBClassifier(learning_rate=0.1,max_depth=3,n_estimators=500,random_state=2)
model.fit(X_train2, y_train)
y_verify_test_pred = model.predict(X_verify_test2)
y_verify_test_pred_proba = model.predict_proba(X_verify_test2)[:, 1]
verify_test_metrics = calculate_metrics(y_verify_test, y_verify_test_pred, y_verify_test_pred_proba)
df = pd.DataFrame([verify_test_metrics],index=['XGBoost_verify_test'])
df.to_csv("~/Feature_analysis/use_topfeature/XGBoost_AMS_topfeature_finalmodel_evaluation_indicator_results.txt", sep='\t', index=True, header=True)


importances = model.feature_importances_
df0 = pd.DataFrame({
    'Feature': X_train2.columns,
    'Importance': importances
})
df0_sorted = df0.sort_values(by='Importance', ascending=False)
df0_sorted.to_csv("~/Feature_analysis/use_topfeature/XGBoost_AMS_topfeature_finalmodel_feature_importances.txt", sep='\t', index=True, header=True)


model = XGBClassifier(random_state=2)
model.fit(X_train2, y_train)
y_verify_test_pred = model.predict(X_verify_test2)
y_verify_test_pred_proba = model.predict_proba(X_verify_test2)[:, 1]
verify_test_metrics = calculate_metrics(y_verify_test, y_verify_test_pred, y_verify_test_pred_proba)
df = pd.DataFrame([verify_test_metrics],index=['XGBoost_verify_test'])
df.to_csv("~/Feature_analysis/use_topfeature2/XGBoost_AMS_topfeature_finalmodel_evaluation_indicator_results.txt", sep='\t', index=True, header=True)

importances = model.feature_importances_
df0 = pd.DataFrame({
    'Feature': X_train2.columns,
    'Importance': importances
})
df0_sorted = df0.sort_values(by='Importance', ascending=False)
df0_sorted.to_csv("~/Feature_analysis/use_topfeature2/XGBoost_AMS_topfeature_finalmodel_feature_importances.txt", sep='\t', index=True, header=True)

