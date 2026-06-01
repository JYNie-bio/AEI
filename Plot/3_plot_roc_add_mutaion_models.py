from pandas.api.types import CategoricalDtype
from sklearn.model_selection import GridSearchCV, KFold
from sklearn.model_selection import StratifiedKFold
from xgboost import XGBClassifier
import pandas as pd
import os
import re
from sklearn.metrics import roc_auc_score

new_working_directory = "/lustre/home/niejingyi2023/project/bjCDC/Feature_analysis"
os.chdir(new_working_directory)


import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc

modelAMP = XGBClassifier(learning_rate=0.01,max_depth=6,n_estimators=500,random_state=2)
modelAMS = XGBClassifier(learning_rate=0.1,max_depth=3,n_estimators=500,random_state=2)
modelCFZ = XGBClassifier(learning_rate=0.01,max_depth=4,n_estimators=1000,random_state=2, enable_categorical=True)
modelCHL = XGBClassifier(learning_rate=0.01,max_depth=3,n_estimators=200,random_state=2)
modelCIP = XGBClassifier(learning_rate=0.1,max_depth=4,n_estimators=500,random_state=2, enable_categorical=True)
modelCTX = XGBClassifier(learning_rate=0.01,max_depth=3,n_estimators=1000,random_state=2)
modelGEN = XGBClassifier(learning_rate=0.3,max_depth=3,n_estimators=500,random_state=2, enable_categorical=True)
modelSXT = XGBClassifier(learning_rate=0.3,max_depth=5,n_estimators=200,random_state=2, enable_categorical=True)
modelTET = XGBClassifier(learning_rate=0.01,max_depth=3,n_estimators=100,random_state=2)
modelAZM = XGBClassifier(learning_rate=0.1,max_depth=8,n_estimators=100,random_state=2, enable_categorical=True)
modelCAZ = XGBClassifier(learning_rate=0.01,max_depth=8,n_estimators=1000,random_state=2, enable_categorical=True)
modelCT = XGBClassifier(learning_rate=0.01,max_depth=4,n_estimators=1000,random_state=2, enable_categorical=True)
modelNAL = XGBClassifier(random_state=2, enable_categorical=True)


AMP_fold_traindata=pd.read_csv("use_topfeature/XGBoost_AMP_topfeature_finalmodel_traindata.txt",sep='\t', header=0,index_col=0)
X_fold_train_AMP=AMP_fold_traindata.iloc[:, :-1]
y_fold_train_AMP=AMP_fold_traindata.iloc[:, -1]
y_fold_train_AMP=y_fold_train_AMP.replace({'R': 1, 'S': 0})
AMP_fold_verifytestdata=pd.read_csv("use_topfeature/XGBoost_AMP_topfeature_finalmodel_verifytestdata.txt",sep='\t', header=0,index_col=0)
X_fold_verifytest_AMP=AMP_fold_verifytestdata.iloc[:, :-1]
y_fold_verifytest_AMP=AMP_fold_verifytestdata.iloc[:, -1]
y_fold_verifytest_AMP=y_fold_verifytest_AMP.replace({'R': 1, 'S': 0})
AMS_fold_traindata=pd.read_csv("use_topfeature/XGBoost_AMS_topfeature_finalmodel_traindata.txt",sep='\t', header=0,index_col=0)
X_fold_train_AMS=AMS_fold_traindata.iloc[:, :-1]
y_fold_train_AMS=AMS_fold_traindata.iloc[:, -1]
y_fold_train_AMS=y_fold_train_AMS.replace({'R': 1, 'S': 0})
AMS_fold_verifytestdata=pd.read_csv("use_topfeature/XGBoost_AMS_topfeature_finalmodel_verifytestdata.txt",sep='\t', header=0,index_col=0)
X_fold_verifytest_AMS=AMS_fold_verifytestdata.iloc[:, :-1]
y_fold_verifytest_AMS=AMS_fold_verifytestdata.iloc[:, -1]
y_fold_verifytest_AMS=y_fold_verifytest_AMS.replace({'R': 1, 'S': 0})

CFZ_fold_traindata=pd.read_csv("use_topfeature/XGBoost_CFZ_topfeature_finalmodel_traindata.txt",sep='\t', header=0,index_col=0, keep_default_na=False, dtype=CategoricalDtype())
X_fold_train_CFZ=CFZ_fold_traindata.iloc[:, :-1]
for col in X_fold_train_CFZ.columns:
    if "|" in col:
        X_fold_train_CFZ[col] = X_fold_train_CFZ[col].cat.set_categories(['Yes', 'No'])
    else:
        X_fold_train_CFZ[col] = X_fold_train_CFZ[col].apply(pd.to_numeric, errors='coerce')
        X_fold_train_CFZ[col] = X_fold_train_CFZ[col].astype('int')

y_fold_train_CFZ=CFZ_fold_traindata.iloc[:, -1]
y_fold_train_CFZ=y_fold_train_CFZ.replace({'R': 1, 'S': 0})
y_fold_train_CFZ=y_fold_train_CFZ.astype('int')

CFZ_fold_verifytestdata=pd.read_csv("use_topfeature/XGBoost_CFZ_topfeature_finalmodel_verifytestdata.txt",sep='\t', header=0,index_col=0, keep_default_na=False, dtype=CategoricalDtype())
X_fold_verifytest_CFZ=CFZ_fold_verifytestdata.iloc[:, :-1]
for col in X_fold_verifytest_CFZ.columns:
    if "|" in col:
        X_fold_verifytest_CFZ[col] = X_fold_verifytest_CFZ[col].cat.set_categories(['Yes', 'No'])
    else:
        X_fold_verifytest_CFZ[col] = X_fold_verifytest_CFZ[col].apply(pd.to_numeric, errors='coerce')
        X_fold_verifytest_CFZ[col] = X_fold_verifytest_CFZ[col].astype('int')

y_fold_verifytest_CFZ=CFZ_fold_verifytestdata.iloc[:, -1]
y_fold_verifytest_CFZ=y_fold_verifytest_CFZ.replace({'R': 1, 'S': 0})
y_fold_verifytest_CFZ=y_fold_verifytest_CFZ.astype('int')

CHL_fold_traindata=pd.read_csv("use_topfeature/XGBoost_CHL_topfeature_finalmodel_traindata.txt",sep='\t', header=0,index_col=0)
X_fold_train_CHL=CHL_fold_traindata.iloc[:, :-1]
y_fold_train_CHL=CHL_fold_traindata.iloc[:, -1]
y_fold_train_CHL=y_fold_train_CHL.replace({'R': 1, 'S': 0})
CHL_fold_verifytestdata=pd.read_csv("use_topfeature/XGBoost_CHL_topfeature_finalmodel_verifytestdata.txt",sep='\t', header=0,index_col=0)
X_fold_verifytest_CHL=CHL_fold_verifytestdata.iloc[:, :-1]
y_fold_verifytest_CHL=CHL_fold_verifytestdata.iloc[:, -1]
y_fold_verifytest_CHL=y_fold_verifytest_CHL.replace({'R': 1, 'S': 0})
CIP_fold_traindata=pd.read_csv("use_topfeature/XGBoost_CIP_topfeature_finalmodel_traindata.txt",sep='\t', header=0,index_col=0, keep_default_na=False, dtype=CategoricalDtype())
X_fold_train_CIP=CIP_fold_traindata.iloc[:, :-1]
for col in X_fold_train_CIP.columns:
    if "|" in col:
        X_fold_train_CIP[col] = X_fold_train_CIP[col].cat.set_categories(['Yes', 'No', 'NA'])
    else:
        X_fold_train_CIP[col] = X_fold_train_CIP[col].cat.set_categories(['0', '1'])

y_fold_train_CIP=CIP_fold_traindata.iloc[:, -1]
y_fold_train_CIP=y_fold_train_CIP.replace({'R': 1, 'S': 0})
y_fold_train_CIP=y_fold_train_CIP.astype('int')

CIP_fold_verifytestdata=pd.read_csv("use_topfeature/XGBoost_CIP_topfeature_finalmodel_verifytestdata.txt",sep='\t', header=0,index_col=0, keep_default_na=False, dtype=CategoricalDtype())
X_fold_verifytest_CIP=CIP_fold_verifytestdata.iloc[:, :-1]
for col in X_fold_verifytest_CIP.columns:
    if "|" in col:
        X_fold_verifytest_CIP[col] = X_fold_verifytest_CIP[col].cat.set_categories(['Yes', 'No', 'NA'])
    else:
        X_fold_verifytest_CIP[col] = X_fold_verifytest_CIP[col].cat.set_categories(['0', '1'])

y_fold_verifytest_CIP=CIP_fold_verifytestdata.iloc[:, -1]
y_fold_verifytest_CIP=y_fold_verifytest_CIP.replace({'R': 1, 'S': 0})
y_fold_verifytest_CIP=y_fold_verifytest_CIP.astype('int')

CTX_fold_traindata=pd.read_csv("use_topfeature/XGBoost_CTX_topfeature_finalmodel_traindata.txt",sep='\t', header=0,index_col=0)
X_fold_train_CTX=CTX_fold_traindata.iloc[:, :-1]
y_fold_train_CTX=CTX_fold_traindata.iloc[:, -1]
y_fold_train_CTX=y_fold_train_CTX.replace({'R': 1, 'S': 0})
CTX_fold_verifytestdata=pd.read_csv("use_topfeature/XGBoost_CTX_topfeature_finalmodel_verifytestdata.txt",sep='\t', header=0,index_col=0)
X_fold_verifytest_CTX=CTX_fold_verifytestdata.iloc[:, :-1]
y_fold_verifytest_CTX=CTX_fold_verifytestdata.iloc[:, -1]
y_fold_verifytest_CTX=y_fold_verifytest_CTX.replace({'R': 1, 'S': 0})

GEN_fold_traindata=pd.read_csv("use_topfeature/XGBoost_GEN_topfeature_finalmodel_traindata.txt",sep='\t', header=0,index_col=0, keep_default_na=False, dtype=CategoricalDtype())
X_fold_train_GEN=GEN_fold_traindata.iloc[:, :-1]
for col in X_fold_train_GEN.columns:
    if "|" in col:
        X_fold_train_GEN[col] = X_fold_train_GEN[col].cat.set_categories(['Yes', 'No', 'NA'])
    else:
        X_fold_train_GEN[col] = X_fold_train_GEN[col].cat.set_categories(['0', '1'])

y_fold_train_GEN=GEN_fold_traindata.iloc[:, -1]
y_fold_train_GEN=y_fold_train_GEN.replace({'R': 1, 'S': 0})
y_fold_train_GEN=y_fold_train_GEN.astype('int')

GEN_fold_verifytestdata=pd.read_csv("use_topfeature/XGBoost_GEN_topfeature_finalmodel_verifytestdata.txt",sep='\t', header=0,index_col=0, keep_default_na=False, dtype=CategoricalDtype())
X_fold_verifytest_GEN=GEN_fold_verifytestdata.iloc[:, :-1]
for col in X_fold_verifytest_GEN.columns:
    if "|" in col:
        X_fold_verifytest_GEN[col] = X_fold_verifytest_GEN[col].cat.set_categories(['Yes', 'No', 'NA'])
    else:
        X_fold_verifytest_GEN[col] = X_fold_verifytest_GEN[col].cat.set_categories(['0', '1'])

y_fold_verifytest_GEN=GEN_fold_verifytestdata.iloc[:, -1]
y_fold_verifytest_GEN=y_fold_verifytest_GEN.replace({'R': 1, 'S': 0})
y_fold_verifytest_GEN=y_fold_verifytest_GEN.astype('int')

SXT_fold_traindata=pd.read_csv("use_topfeature/XGBoost_SXT_topfeature_finalmodel_traindata.txt",sep='\t', header=0,index_col=0, keep_default_na=False, dtype=CategoricalDtype())
X_fold_train_SXT=SXT_fold_traindata.iloc[:, :-1]
for col in X_fold_train_SXT.columns:
    if pd.api.types.is_categorical_dtype(X_fold_train_SXT[col]):
        X_fold_train_SXT[col] = X_fold_train_SXT[col].cat.set_categories(['0', '1'])

y_fold_train_SXT=SXT_fold_traindata.iloc[:, -1]
y_fold_train_SXT=y_fold_train_SXT.replace({'R': 1, 'S': 0})
y_fold_train_SXT=y_fold_train_SXT.astype('int')

SXT_fold_verifytestdata=pd.read_csv("use_topfeature/XGBoost_SXT_topfeature_finalmodel_verifytestdata.txt",sep='\t', header=0,index_col=0, keep_default_na=False, dtype=CategoricalDtype())
X_fold_verifytest_SXT=SXT_fold_verifytestdata.iloc[:, :-1]
for col in X_fold_verifytest_SXT.columns:
    if pd.api.types.is_categorical_dtype(X_fold_verifytest_SXT[col]):
        X_fold_verifytest_SXT[col] = X_fold_verifytest_SXT[col].cat.set_categories(['0', '1'])

y_fold_verifytest_SXT=SXT_fold_verifytestdata.iloc[:, -1]
y_fold_verifytest_SXT=y_fold_verifytest_SXT.replace({'R': 1, 'S': 0})
y_fold_verifytest_SXT=y_fold_verifytest_SXT.astype('int')

TET_fold_traindata=pd.read_csv("use_topfeature/XGBoost_TET_topfeature_finalmodel_traindata.txt",sep='\t', header=0,index_col=0)
X_fold_train_TET=TET_fold_traindata.iloc[:, :-1]
y_fold_train_TET=TET_fold_traindata.iloc[:, -1]
y_fold_train_TET=y_fold_train_TET.replace({'R': 1, 'S': 0})
TET_fold_verifytestdata=pd.read_csv("use_topfeature/XGBoost_TET_topfeature_finalmodel_verifytestdata.txt",sep='\t', header=0,index_col=0)
X_fold_verifytest_TET=TET_fold_verifytestdata.iloc[:, :-1]
y_fold_verifytest_TET=TET_fold_verifytestdata.iloc[:, -1]
y_fold_verifytest_TET=y_fold_verifytest_TET.replace({'R': 1, 'S': 0})

AZM_fold_traindata=pd.read_csv("use_topfeature/XGBoost_AZM_topfeature_finalmodel_traindata.txt",sep='\t', header=0,index_col=0, keep_default_na=False, dtype=CategoricalDtype())
X_fold_train_AZM=AZM_fold_traindata.iloc[:, :-1]
for col in X_fold_train_AZM.columns:
    if "|" in col:
        X_fold_train_AZM[col] = X_fold_train_AZM[col].cat.set_categories(['Yes', 'No', 'NA'])
    else:
        X_fold_train_AZM[col] = X_fold_train_AZM[col].cat.set_categories(['0', '1'])

y_fold_train_AZM=AZM_fold_traindata.iloc[:, -1]
y_fold_train_AZM=y_fold_train_AZM.replace({'R': 1, 'S': 0})
y_fold_train_AZM=y_fold_train_AZM.astype('int')

AZM_fold_verifytestdata=pd.read_csv("use_topfeature/XGBoost_AZM_topfeature_finalmodel_verifytestdata.txt",sep='\t', header=0,index_col=0, keep_default_na=False, dtype=CategoricalDtype())
X_fold_verifytest_AZM=AZM_fold_verifytestdata.iloc[:, :-1]
for col in X_fold_verifytest_AZM.columns:
    if "|" in col:
        X_fold_verifytest_AZM[col] = X_fold_verifytest_AZM[col].cat.set_categories(['Yes', 'No', 'NA'])
    else:
        X_fold_verifytest_AZM[col] = X_fold_verifytest_AZM[col].cat.set_categories(['0', '1'])

y_fold_verifytest_AZM=AZM_fold_verifytestdata.iloc[:, -1]
y_fold_verifytest_AZM=y_fold_verifytest_AZM.replace({'R': 1, 'S': 0})
y_fold_verifytest_AZM=y_fold_verifytest_AZM.astype('int')

CAZ_fold_traindata=pd.read_csv("use_topfeature/XGBoost_CAZ_topfeature_finalmodel_traindata.txt",sep='\t', header=0,index_col=0, keep_default_na=False, dtype=CategoricalDtype())
X_fold_train_CAZ=CAZ_fold_traindata.iloc[:, :-1]
for col in X_fold_train_CAZ.columns:
    if pd.api.types.is_categorical_dtype(X_fold_train_CAZ[col]):
        X_fold_train_CAZ[col] = X_fold_train_CAZ[col].cat.set_categories(['0', '1'])

y_fold_train_CAZ=CAZ_fold_traindata.iloc[:, -1]
y_fold_train_CAZ=y_fold_train_CAZ.replace({'R': 1, 'S': 0})
y_fold_train_CAZ=y_fold_train_CAZ.astype('int')

CAZ_fold_verifytestdata=pd.read_csv("use_topfeature/XGBoost_CAZ_topfeature_finalmodel_verifytestdata.txt",sep='\t', header=0,index_col=0, keep_default_na=False, dtype=CategoricalDtype())
X_fold_verifytest_CAZ=CAZ_fold_verifytestdata.iloc[:, :-1]
for col in X_fold_verifytest_CAZ.columns:
    if pd.api.types.is_categorical_dtype(X_fold_verifytest_CAZ[col]):
        X_fold_verifytest_CAZ[col] = X_fold_verifytest_CAZ[col].cat.set_categories(['0', '1'])

y_fold_verifytest_CAZ=CAZ_fold_verifytestdata.iloc[:, -1]
y_fold_verifytest_CAZ=y_fold_verifytest_CAZ.replace({'R': 1, 'S': 0})
y_fold_verifytest_CAZ=y_fold_verifytest_CAZ.astype('int')

CT_fold_traindata=pd.read_csv("use_topfeature/XGBoost_CT_topfeature_finalmodel_traindata.txt",sep='\t', header=0,index_col=0, keep_default_na=False, dtype=CategoricalDtype())
X_fold_train_CT=CT_fold_traindata.iloc[:, :-1]
for col in X_fold_train_CT.columns:
    if pd.api.types.is_categorical_dtype(X_fold_train_CT[col]):
        X_fold_train_CT[col] = X_fold_train_CT[col].cat.set_categories(['0', '1'])

y_fold_train_CT=CT_fold_traindata.iloc[:, -1]
y_fold_train_CT=y_fold_train_CT.replace({'R': 1, 'S': 0})
y_fold_train_CT=y_fold_train_CT.astype('int')

CT_fold_verifytestdata=pd.read_csv("use_topfeature/XGBoost_CT_topfeature_finalmodel_verifytestdata.txt",sep='\t', header=0,index_col=0, keep_default_na=False, dtype=CategoricalDtype())
X_fold_verifytest_CT=CT_fold_verifytestdata.iloc[:, :-1]
for col in X_fold_verifytest_CT.columns:
    if pd.api.types.is_categorical_dtype(X_fold_verifytest_CT[col]):
        X_fold_verifytest_CT[col] = X_fold_verifytest_CT[col].cat.set_categories(['0', '1'])

y_fold_verifytest_CT=CT_fold_verifytestdata.iloc[:, -1]
y_fold_verifytest_CT=y_fold_verifytest_CT.replace({'R': 1, 'S': 0})
y_fold_verifytest_CT=y_fold_verifytest_CT.astype('int')

NAL_fold_traindata=pd.read_csv("use_topfeature/XGBoost_NAL_topfeature_finalmodel_traindata.txt",sep='\t', header=0,index_col=0, keep_default_na=False, dtype=CategoricalDtype())
X_fold_train_NAL=NAL_fold_traindata.iloc[:, :-1]
for col in X_fold_train_NAL.columns:
    if "|" in col:
        X_fold_train_NAL[col] = X_fold_train_NAL[col].cat.set_categories(['Yes', 'No', 'NA'])
    else:
        X_fold_train_NAL[col] = X_fold_train_NAL[col].cat.set_categories(['0', '1'])

y_fold_train_NAL=NAL_fold_traindata.iloc[:, -1]
y_fold_train_NAL=y_fold_train_NAL.replace({'R': 1, 'S': 0})
y_fold_train_NAL=y_fold_train_NAL.astype('int')

NAL_fold_verifytestdata=pd.read_csv("use_topfeature/XGBoost_NAL_topfeature_finalmodel_verifytestdata.txt",sep='\t', header=0,index_col=0, keep_default_na=False, dtype=CategoricalDtype())
X_fold_verifytest_NAL=NAL_fold_verifytestdata.iloc[:, :-1]
for col in X_fold_verifytest_NAL.columns:
    if "|" in col:
        X_fold_verifytest_NAL[col] = X_fold_verifytest_NAL[col].cat.set_categories(['Yes', 'No', 'NA'])
    else:
        X_fold_verifytest_NAL[col] = X_fold_verifytest_NAL[col].cat.set_categories(['0', '1'])

y_fold_verifytest_NAL=NAL_fold_verifytestdata.iloc[:, -1]
y_fold_verifytest_NAL=y_fold_verifytest_NAL.replace({'R': 1, 'S': 0})
y_fold_verifytest_NAL=y_fold_verifytest_NAL.astype('int')


def plot_roc_curves(models, output_path):
    plt.figure()

    for model_name, (model, X_test_train, y_test_train, X_test_verifytest, y_test_verifytest) in models.items():
        model.fit(X_test_train, y_test_train)
        y_pred_proba_verifytest = model.predict_proba(X_test_verifytest)[:, 1]
        fpr, tpr, _ = roc_curve(y_test_verifytest, y_pred_proba_verifytest)
        roc_auc = roc_auc_score(y_test_verifytest, y_pred_proba_verifytest)
        plt.plot(fpr, tpr, label=f'{model_name} (AUC = {roc_auc:.2f})')

    plt.plot([0, 1], [0, 1], 'k--', lw=2)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve')
    plt.legend(loc="lower right")
    plt.savefig(output_path, format='pdf')
    plt.show()

# 示例：假设有多个模型，且每个模型有不同的输入数据
models = {
    'AMP': (modelAMP, X_fold_train_AMP, y_fold_train_AMP, X_fold_verifytest_AMP, y_fold_verifytest_AMP),
    'AMS': (modelAMS, X_fold_train_AMS, y_fold_train_AMS, X_fold_verifytest_AMS, y_fold_verifytest_AMS),
    'AZM': (modelAZM, X_fold_train_AZM, y_fold_train_AZM, X_fold_verifytest_AZM, y_fold_verifytest_AZM),
    'CAZ': (modelCAZ, X_fold_train_CAZ, y_fold_train_CAZ, X_fold_verifytest_CAZ, y_fold_verifytest_CAZ),
    'CFZ': (modelCFZ, X_fold_train_CFZ, y_fold_train_CFZ, X_fold_verifytest_CFZ, y_fold_verifytest_CFZ),
    'CHL': (modelCHL, X_fold_train_CHL, y_fold_train_CHL, X_fold_verifytest_CHL, y_fold_verifytest_CHL),
    'CIP': (modelCIP, X_fold_train_CIP, y_fold_train_CIP, X_fold_verifytest_CIP, y_fold_verifytest_CIP),
    'CT': (modelCT, X_fold_train_CT, y_fold_train_CT, X_fold_verifytest_CT, y_fold_verifytest_CT),
    'CTX': (modelCTX, X_fold_train_CTX, y_fold_train_CTX, X_fold_verifytest_CTX, y_fold_verifytest_CTX),
    'GEN': (modelGEN, X_fold_train_GEN, y_fold_train_GEN, X_fold_verifytest_GEN, y_fold_verifytest_GEN),
    'NAL': (modelNAL, X_fold_train_NAL, y_fold_train_NAL, X_fold_verifytest_NAL, y_fold_verifytest_NAL),
    'SXT': (modelSXT, X_fold_train_SXT, y_fold_train_SXT, X_fold_verifytest_SXT, y_fold_verifytest_SXT),
    'TET': (modelTET, X_fold_train_TET, y_fold_train_TET, X_fold_verifytest_TET, y_fold_verifytest_TET)
}

plot_roc_curves(models, 'roc_curve.pdf')

