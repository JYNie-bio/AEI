from pandas.api.types import CategoricalDtype
from sklearn.model_selection import GridSearchCV, KFold
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn import metrics
from xgboost import XGBClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
import math
import numpy as np
import pandas as pd
import os
import sys
from sklearn.metrics import make_scorer, accuracy_score, roc_auc_score, f1_score, recall_score, precision_score, cohen_kappa_score, balanced_accuracy_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

arguments = sys.argv
new_working_directory = "~/ARG_with_SNP/"
os.chdir(new_working_directory)


traindata=pd.read_csv("RS_811/14880file_ARGanno_copynumber_SNPmatrix959090_absence01_addPhenotype_" + arguments[1] + "_train.txt",sep='\t', header=0,index_col=0, keep_default_na=False, dtype=CategoricalDtype())
verifytestdata=pd.read_csv("RS_811/14880file_ARGanno_copynumber_SNPmatrix959090_absence01_addPhenotype_" + arguments[1] + "_verify_test.txt",sep='\t', header=0,index_col=0, keep_default_na=False, dtype=CategoricalDtype())
X_train=traindata.iloc[:, :-1]
numeric_columns = list(range(0, 908))
category_columns = list(range(909, X_train.shape[1]))
X_train.iloc[:, numeric_columns] = X_train.iloc[:, numeric_columns].apply(pd.to_numeric, errors='coerce')
X_train.iloc[:, numeric_columns] = X_train.iloc[:, numeric_columns].astype('int')
for col in X_train.columns[category_columns]:
    X_train[col] = X_train[col].cat.set_categories(['Yes', 'No'])

y_train=traindata.iloc[:, -1]
y_train=y_train.replace({'R': 1, 'S': 0})
y_train=y_train.astype('int')
X_verify_test=verifytestdata.iloc[:, :-1]
X_verify_test.iloc[:, numeric_columns] = X_verify_test.iloc[:, numeric_columns].apply(pd.to_numeric, errors='coerce')
X_verify_test.iloc[:, numeric_columns] = X_verify_test.iloc[:, numeric_columns].astype('int')
for col in X_verify_test.columns[category_columns]:
    X_verify_test[col] = X_verify_test[col].cat.set_categories(['Yes', 'No'])

y_verify_test=verifytestdata.iloc[:, -1]
y_verify_test=y_verify_test.replace({'R': 1, 'S': 0})
y_verify_test=y_verify_test.astype('int')


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
X_fold_train, X_fold_val = None, None
y_fold_train, y_fold_val = None, None
strnum=int(arguments[2])-1

for fold, (train_idx, val_idx) in enumerate(skf.split(X_train, y_train)):
    if fold==strnum:
        X_fold_train, X_fold_val = X_train.iloc[train_idx], X_train.iloc[val_idx]
        y_fold_train, y_fold_val = y_train.iloc[train_idx], y_train.iloc[val_idx]

if X_fold_train is None or y_fold_train is None:
    print("Not Found Fold or undefined training data")

model = XGBClassifier(random_state=2, enable_categorical=True)
model.fit(X_fold_train, y_fold_train)
importances = model.feature_importances_
df0 = pd.DataFrame({
    'Feature': X_fold_train.columns,
    'Importance': importances
})
df0_sorted = df0.sort_values(by='Importance', ascending=False)
df0_filtered = df0_sorted[df0_sorted['Importance'] > 0]
df0_sorted.to_csv("AllFeature/XGBoost_10fold_crossvalidation_default_hyperparameter_ARGanno_copynumber_SNPmatrix959090_absence01_"+ arguments[1] +"_fold_" + arguments[2] + "_feature_importances_category_and_combine_verifytest.txt", sep='\t', index=True, header=True)


y_val_pred = model.predict(X_fold_val)
y_verify_test_pred = model.predict(X_verify_test)
y_val_pred_proba = model.predict_proba(X_fold_val)[:, 1]
y_verify_test_pred_proba = model.predict_proba(X_verify_test)[:, 1]
val_metrics = calculate_metrics(y_fold_val, y_val_pred, y_val_pred_proba)
verify_test_metrics = calculate_metrics(y_verify_test, y_verify_test_pred, y_verify_test_pred_proba)
df = pd.DataFrame([val_metrics, verify_test_metrics],index=['XGBoost_val_featureall', 'XGBoost_verify_test_featureall'])
df.to_csv("AllFeature/XGBoost_10fold_crossvalidation_default_hyperparameter_ARGanno_copynumber_SNPmatrix959090_absence01_"+ arguments[1] +"_fold_" + arguments[2] + "_topfeature_category_and_combine_verifytest.txt", sep='\t', index=True, header=True)


df0_filtered_top=df0_filtered.head(5)
X_fold_train2=X_fold_train.iloc[:,df0_filtered_top.index]
X_fold_val2=X_fold_val.iloc[:,df0_filtered_top.index]
X_verify_test2=X_verify_test.iloc[:,df0_filtered_top.index]
X_train2=X_train.iloc[:,df0_filtered_top.index]
model = XGBClassifier(random_state=2, enable_categorical=True)
model.fit(X_fold_train2, y_fold_train)
y_val_pred = model.predict(X_fold_val2)
y_verify_test_pred = model.predict(X_verify_test2)
y_val_pred_proba = model.predict_proba(X_fold_val2)[:, 1]
y_verify_test_pred_proba = model.predict_proba(X_verify_test2)[:, 1]
val_metrics = calculate_metrics(y_fold_val, y_val_pred, y_val_pred_proba)
verify_test_metrics = calculate_metrics(y_verify_test, y_verify_test_pred, y_verify_test_pred_proba)
df5 = pd.DataFrame([val_metrics, verify_test_metrics],index=['XGBoost_val_feature5', 'XGBoost_verify_test_feature5'])
combined_df = pd.concat([df, df5], axis=0)
combined_df.to_csv("AllFeature/XGBoost_10fold_crossvalidation_default_hyperparameter_ARGanno_copynumber_SNPmatrix959090_absence01_"+ arguments[1] +"_fold_" + arguments[2] + "_topfeature_category_and_combine_verifytest.txt", sep='\t', index=True, header=True)

cv = KFold(n_splits=10, shuffle=True, random_state=42)
def prevalence_score(y_true, y_pred):
    return sum(y_true) / len(y_true)
def specificity_score(y_true, y_pred):
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    return tn/(fp+tn)
def sensitivity_score(y_true, y_pred):
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    return tp/(tp+fn)
def recall_score(y_true, y_pred):
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    return tp/(tp+fn)
def precision_score(y_true, y_pred):
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    return tp/(tp+fp)
def pos_pred_value_score(y_true, y_pred):
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    return tp / (tp + fp)
def neg_pred_value_score(y_true, y_pred):
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    return tn / (tn + fn)
scoring = {'AUC': make_scorer(roc_auc_score),
           'Balanced_AUC': make_scorer(roc_auc_score, multi_class='ovr', average='weighted'),
           'Accuracy': make_scorer(accuracy_score),
           'Balanced_Accuracy': make_scorer(balanced_accuracy_score),
           'F1_Score': make_scorer(f1_score),
           'Sensitivity': make_scorer(sensitivity_score),
           'specificity': make_scorer(specificity_score),
           'Kappa': make_scorer(cohen_kappa_score),
           'Pos_Pred_Value': make_scorer(pos_pred_value_score),
           'Neg_Pred_Value': make_scorer(neg_pred_value_score),
           'Precision': make_scorer(precision_score),
           'Recall': make_scorer(recall_score),
           'Prevalence': make_scorer(prevalence_score)}

xgb_params = {'n_estimators': [100,200,300,500,700,1000,1500],
              'max_depth': [3,4,5,6,7,8,9,10],
              'learning_rate': [0.01, 0.1, 0.3]
        }
xgb_model = XGBClassifier(random_state=2, enable_categorical=True)
xgb_grid = GridSearchCV(xgb_model, param_grid=xgb_params, scoring=scoring,refit='Accuracy', cv=cv)
xgb_grid.fit(X_train2, y_train)
print("XGBoost Best Parameters:", xgb_grid.best_params_)
print("XGBoost Best Accuracy:", xgb_grid.best_score_)
results_XGBoost = pd.DataFrame(xgb_grid.cv_results_)
results_XGBoost.to_csv("AllFeature/XGBoost_10fold_crossvalidation_default_hyperparameter_ARGanno_copynumber_SNPmatrix959090_absence01_"+ arguments[1] +"_fold_" + arguments[2] + "_top5feature_category_and_combine_verifytest.txt", sep='\t', index=False, header=True)


df0_filtered_top=df0_filtered.head(10)
X_fold_train2=X_fold_train.iloc[:,df0_filtered_top.index]
X_fold_val2=X_fold_val.iloc[:,df0_filtered_top.index]
X_verify_test2=X_verify_test.iloc[:,df0_filtered_top.index]
X_train2=X_train.iloc[:,df0_filtered_top.index]
model = XGBClassifier(random_state=2, enable_categorical=True)
model.fit(X_fold_train2, y_fold_train)
y_val_pred = model.predict(X_fold_val2)
y_verify_test_pred = model.predict(X_verify_test2)
y_val_pred_proba = model.predict_proba(X_fold_val2)[:, 1]
y_verify_test_pred_proba = model.predict_proba(X_verify_test2)[:, 1]
val_metrics = calculate_metrics(y_fold_val, y_val_pred, y_val_pred_proba)
verify_test_metrics = calculate_metrics(y_verify_test, y_verify_test_pred, y_verify_test_pred_proba)
df10 = pd.DataFrame([val_metrics, verify_test_metrics],index=['XGBoost_val_feature10', 'XGBoost_verify_test_feature10'])
combined_df = pd.concat([combined_df, df10], axis=0)
combined_df.to_csv("AllFeature/XGBoost_10fold_crossvalidation_default_hyperparameter_ARGanno_copynumber_SNPmatrix959090_absence01_"+ arguments[1] +"_fold_" + arguments[2] + "_topfeature_category_and_combine_verifytest.txt", sep='\t', index=True, header=True)

xgb_model = XGBClassifier(random_state=2, enable_categorical=True)
xgb_grid = GridSearchCV(xgb_model, param_grid=xgb_params, scoring=scoring,refit='Accuracy', cv=cv)
xgb_grid.fit(X_train2, y_train)
print("XGBoost Best Parameters:", xgb_grid.best_params_)
print("XGBoost Best Accuracy:", xgb_grid.best_score_)
results_XGBoost = pd.DataFrame(xgb_grid.cv_results_)
results_XGBoost.to_csv("AllFeature/XGBoost_10fold_crossvalidation_default_hyperparameter_ARGanno_copynumber_SNPmatrix959090_absence01_"+ arguments[1] +"_fold_" + arguments[2] + "_top10feature_category_and_combine_verifytest.txt", sep='\t', index=False, header=True)

df0_filtered_top=df0_filtered.head(20)
X_fold_train2=X_fold_train.iloc[:,df0_filtered_top.index]
X_fold_val2=X_fold_val.iloc[:,df0_filtered_top.index]
X_verify_test2=X_verify_test.iloc[:,df0_filtered_top.index]
X_train2=X_train.iloc[:,df0_filtered_top.index]
model = XGBClassifier(random_state=2, enable_categorical=True)
model.fit(X_fold_train2, y_fold_train)
y_val_pred = model.predict(X_fold_val2)
y_verify_test_pred = model.predict(X_verify_test2)
y_val_pred_proba = model.predict_proba(X_fold_val2)[:, 1]
y_verify_test_pred_proba = model.predict_proba(X_verify_test2)[:, 1]
val_metrics = calculate_metrics(y_fold_val, y_val_pred, y_val_pred_proba)
verify_test_metrics = calculate_metrics(y_verify_test, y_verify_test_pred, y_verify_test_pred_proba)
df20 = pd.DataFrame([val_metrics, verify_test_metrics],index=['XGBoost_val_feature20', 'XGBoost_verify_test_feature20'])
combined_df = pd.concat([combined_df, df20], axis=0)
combined_df.to_csv("AllFeature/XGBoost_10fold_crossvalidation_default_hyperparameter_ARGanno_copynumber_SNPmatrix959090_absence01_"+ arguments[1] +"_fold_" + arguments[2] + "_topfeature_category_and_combine_verifytest.txt", sep='\t', index=True, header=True)

xgb_model = XGBClassifier(random_state=2, enable_categorical=True)
xgb_grid = GridSearchCV(xgb_model, param_grid=xgb_params, scoring=scoring,refit='Accuracy', cv=cv)
xgb_grid.fit(X_train2, y_train)
print("XGBoost Best Parameters:", xgb_grid.best_params_)
print("XGBoost Best Accuracy:", xgb_grid.best_score_)
results_XGBoost = pd.DataFrame(xgb_grid.cv_results_)
results_XGBoost.to_csv("AllFeature/XGBoost_10fold_crossvalidation_default_hyperparameter_ARGanno_copynumber_SNPmatrix959090_absence01_"+ arguments[1] +"_fold_" + arguments[2] + "_top20feature_category_and_combine_verifytest.txt", sep='\t', index=False, header=True)

df0_filtered_top=df0_filtered.head(30)
X_fold_train2=X_fold_train.iloc[:,df0_filtered_top.index]
X_fold_val2=X_fold_val.iloc[:,df0_filtered_top.index]
X_verify_test2=X_verify_test.iloc[:,df0_filtered_top.index]
X_train2=X_train.iloc[:,df0_filtered_top.index]
model = XGBClassifier(random_state=2, enable_categorical=True)
model.fit(X_fold_train2, y_fold_train)
y_val_pred = model.predict(X_fold_val2)
y_verify_test_pred = model.predict(X_verify_test2)
y_val_pred_proba = model.predict_proba(X_fold_val2)[:, 1]
y_verify_test_pred_proba = model.predict_proba(X_verify_test2)[:, 1]
val_metrics = calculate_metrics(y_fold_val, y_val_pred, y_val_pred_proba)
verify_test_metrics = calculate_metrics(y_verify_test, y_verify_test_pred, y_verify_test_pred_proba)
df30 = pd.DataFrame([val_metrics, verify_test_metrics],index=['XGBoost_val_feature30', 'XGBoost_verify_test_feature30'])
combined_df = pd.concat([combined_df, df30], axis=0)
combined_df.to_csv("AllFeature/XGBoost_10fold_crossvalidation_default_hyperparameter_ARGanno_copynumber_SNPmatrix959090_absence01_"+ arguments[1] +"_fold_" + arguments[2] + "_topfeature_category_and_combine_verifytest.txt", sep='\t', index=True, header=True)

xgb_model = XGBClassifier(random_state=2, enable_categorical=True)
xgb_grid = GridSearchCV(xgb_model, param_grid=xgb_params, scoring=scoring,refit='Accuracy', cv=cv)
xgb_grid.fit(X_train2, y_train)
print("XGBoost Best Parameters:", xgb_grid.best_params_)
print("XGBoost Best Accuracy:", xgb_grid.best_score_)
results_XGBoost = pd.DataFrame(xgb_grid.cv_results_)
results_XGBoost.to_csv("AllFeature/XGBoost_10fold_crossvalidation_default_hyperparameter_ARGanno_copynumber_SNPmatrix959090_absence01_"+ arguments[1] +"_fold_" + arguments[2] + "_top30feature_category_and_combine_verifytest.txt", sep='\t', index=False, header=True)

df0_filtered_top=df0_filtered.head(40)
X_fold_train2=X_fold_train.iloc[:,df0_filtered_top.index]
X_fold_val2=X_fold_val.iloc[:,df0_filtered_top.index]
X_verify_test2=X_verify_test.iloc[:,df0_filtered_top.index]
X_train2=X_train.iloc[:,df0_filtered_top.index]
model = XGBClassifier(random_state=2, enable_categorical=True)
model.fit(X_fold_train2, y_fold_train)
y_val_pred = model.predict(X_fold_val2)
y_verify_test_pred = model.predict(X_verify_test2)
y_val_pred_proba = model.predict_proba(X_fold_val2)[:, 1]
y_verify_test_pred_proba = model.predict_proba(X_verify_test2)[:, 1]
val_metrics = calculate_metrics(y_fold_val, y_val_pred, y_val_pred_proba)
verify_test_metrics = calculate_metrics(y_verify_test, y_verify_test_pred, y_verify_test_pred_proba)
df40 = pd.DataFrame([val_metrics, verify_test_metrics],index=['XGBoost_val_feature40', 'XGBoost_verify_test_feature40'])
combined_df = pd.concat([combined_df, df40], axis=0)
combined_df.to_csv("AllFeature/XGBoost_10fold_crossvalidation_default_hyperparameter_ARGanno_copynumber_SNPmatrix959090_absence01_"+ arguments[1] +"_fold_" + arguments[2] + "_topfeature_category_and_combine_verifytest.txt", sep='\t', index=True, header=True)

xgb_model = XGBClassifier(random_state=2, enable_categorical=True)
xgb_grid = GridSearchCV(xgb_model, param_grid=xgb_params, scoring=scoring,refit='Accuracy', cv=cv)
xgb_grid.fit(X_train2, y_train)
print("XGBoost Best Parameters:", xgb_grid.best_params_)
print("XGBoost Best Accuracy:", xgb_grid.best_score_)
results_XGBoost = pd.DataFrame(xgb_grid.cv_results_)
results_XGBoost.to_csv("AllFeature/XGBoost_10fold_crossvalidation_default_hyperparameter_ARGanno_copynumber_SNPmatrix959090_absence01_"+ arguments[1] +"_fold_" + arguments[2] + "_top40feature_category_and_combine_verifytest.txt", sep='\t', index=False, header=True)

df0_filtered_top=df0_filtered.head(50)
X_fold_train2=X_fold_train.iloc[:,df0_filtered_top.index]
X_fold_val2=X_fold_val.iloc[:,df0_filtered_top.index]
X_verify_test2=X_verify_test.iloc[:,df0_filtered_top.index]
X_train2=X_train.iloc[:,df0_filtered_top.index]
model = XGBClassifier(random_state=2, enable_categorical=True)
model.fit(X_fold_train2, y_fold_train)
y_val_pred = model.predict(X_fold_val2)
y_verify_test_pred = model.predict(X_verify_test2)
y_val_pred_proba = model.predict_proba(X_fold_val2)[:, 1]
y_verify_test_pred_proba = model.predict_proba(X_verify_test2)[:, 1]
val_metrics = calculate_metrics(y_fold_val, y_val_pred, y_val_pred_proba)
verify_test_metrics = calculate_metrics(y_verify_test, y_verify_test_pred, y_verify_test_pred_proba)
df50 = pd.DataFrame([val_metrics, verify_test_metrics],index=['XGBoost_val_feature50', 'XGBoost_verify_test_feature50'])
combined_df = pd.concat([combined_df, df50], axis=0)
combined_df.to_csv("AllFeature/XGBoost_10fold_crossvalidation_default_hyperparameter_ARGanno_copynumber_SNPmatrix959090_absence01_"+ arguments[1] +"_fold_" + arguments[2] + "_topfeature_category_and_combine_verifytest.txt", sep='\t', index=True, header=True)

xgb_model = XGBClassifier(random_state=2, enable_categorical=True)
xgb_grid = GridSearchCV(xgb_model, param_grid=xgb_params, scoring=scoring,refit='Accuracy', cv=cv)
xgb_grid.fit(X_train2, y_train)
print("XGBoost Best Parameters:", xgb_grid.best_params_)
print("XGBoost Best Accuracy:", xgb_grid.best_score_)
results_XGBoost = pd.DataFrame(xgb_grid.cv_results_)
results_XGBoost.to_csv("AllFeature/XGBoost_10fold_crossvalidation_default_hyperparameter_ARGanno_copynumber_SNPmatrix959090_absence01_"+ arguments[1] +"_fold_" + arguments[2] + "_top50feature_category_and_combine_verifytest.txt", sep='\t', index=False, header=True)

