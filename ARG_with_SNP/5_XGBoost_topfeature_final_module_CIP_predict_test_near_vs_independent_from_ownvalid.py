# 导入必要的库
#from sklearn.model_selection import GridSearchCV, KFold
#from sklearn.model_selection import StratifiedKFold
#from sklearn.metrics import accuracy_score
from pandas.api.types import CategoricalDtype
from sklearn.metrics import confusion_matrix
from sklearn import metrics
from xgboost import XGBClassifier
#from sklearn.ensemble import RandomForestClassifier
#from sklearn.svm import SVC
#from sklearn.inspection import permutation_importance
#import math
import numpy as np
import pandas as pd ##读文件
import os ##改工作目录
import sys ##传参
from sklearn.metrics import make_scorer, accuracy_score, roc_auc_score, f1_score, recall_score, precision_score, cohen_kappa_score, balanced_accuracy_score
#from sklearn.preprocessing import StandardScaler
#from sklearn.pipeline import Pipeline

# 传参
arguments = sys.argv
# 设置新的工作目录路径
new_working_directory = "~/ARG_with_SNP/"
# 更改工作目录
os.chdir(new_working_directory)

# 定义评估指标计算函数
def calculate_metrics(y_true, y_pred, y_pred_proba):
    # 计算混淆矩阵组件
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    # 计算指标
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

## 训练集
traindata=pd.read_csv("~/Feature_analysis/use_topfeature/XGBoost_CIP_topfeature_finalmodel_traindata.txt",sep='\t', header=0,index_col=0, keep_default_na=False, dtype=CategoricalDtype())
# 提取特征矩阵 X 和目标变量 y
X_train=traindata.iloc[:, :-1]# 选择除最后一列之外的所有列作为特征
for col in X_train.columns:
    if "|" in col:  # 判断列名中是否包含 "|"
        X_train[col] = X_train[col].cat.set_categories(['Yes', 'No', 'NA'])
    else:
        X_train[col] = X_train[col].cat.set_categories(['0', '1'])

y_train=traindata.iloc[:, -1]# 选择最后一列作为目标变量
# 将字符串标签转换为整数标签
y_train=y_train.replace({'R': 1, 'S': 0})
y_train=y_train.astype('int')


# 模型
model = XGBClassifier(learning_rate=0.1,max_depth=4,n_estimators=500,random_state=2, enable_categorical=True)
model.fit(X_train, y_train)

####### 独立的测试集cgmlst_near_inOwnValid #########
verifytestdata=pd.read_csv("RS_811/14880file_ARGanno_absence01_SNPmatrix959090_YesNoNA_addPhenotype_CIP_verify_test_cgmlst_near_inOwnValid.txt",sep='\t', header=0,index_col=0, keep_default_na=False, dtype=CategoricalDtype())
y_verify_test=verifytestdata.iloc[:, -1]
y_verify_test=y_verify_test.replace({'R': 1, 'S': 0})
y_verify_test=y_verify_test.astype('int')

verifytestdata2 = verifytestdata.loc[:, X_train.columns.tolist()]
for col in verifytestdata2.columns:
    if "|" in col:  # 判断列名中是否包含 "|"
        verifytestdata2[col] = verifytestdata2[col].cat.set_categories(['Yes', 'No', 'NA'])
    else:
        verifytestdata2[col] = verifytestdata2[col].cat.set_categories(['0', '1'])


# 预测独立的测试集和验证集数据
y_verify_test_pred = model.predict(verifytestdata2)
y_verify_test_pred_proba = model.predict_proba(verifytestdata2)[:, 1]
# 计算验证集和独立测试集的指标
verify_test_metrics = calculate_metrics(y_verify_test, y_verify_test_pred, y_verify_test_pred_proba)
# 将结果放入DataFrame
df = pd.DataFrame([verify_test_metrics],index=['XGBoost_test_cgmlst_near_inOwnValid'])
# 保存文件
df.to_csv("RS_811/14880file_ARGanno_absence01_SNPmatrix959090_YesNoNA_addPhenotype_CIP_verify_test_cgmlst_near_inOwnValid_finalmodel_XGBoostresult.txt", sep='\t', index=True, header=True)

### 输出测试集用matrix文件（真实值 + 预测值 + 预测概率）
output_df = verifytestdata2.copy()
output_df['y_true'] = y_verify_test          # 真实标签 R=1, S=0
output_df['y_pred'] = y_verify_test_pred    # 模型预测结果
output_df['y_pred_prob'] = y_verify_test_pred_proba  # 预测为1（R）的概率
output_df.to_csv("RS_811/14880file_ARGanno_absence01_SNPmatrix959090_YesNoNA_addPhenotype_CIP_verify_test_cgmlst_near_inOwnValid_finalmodel_usematrix.txt",sep="\t",header=True,index=True)

####### 独立的测试集cgmlst_independent_inOwnValid #########
verifytestdata=pd.read_csv("RS_811/14880file_ARGanno_absence01_SNPmatrix959090_YesNoNA_addPhenotype_CIP_verify_test_cgmlst_independent_inOwnValid.txt",sep='\t', header=0,index_col=0, keep_default_na=False, dtype=CategoricalDtype())
y_verify_test=verifytestdata.iloc[:, -1]
y_verify_test=y_verify_test.replace({'R': 1, 'S': 0})
y_verify_test=y_verify_test.astype('int')

verifytestdata2 = verifytestdata.loc[:, X_train.columns.tolist()]
for col in verifytestdata2.columns:
    if "|" in col:  # 判断列名中是否包含 "|"
        verifytestdata2[col] = verifytestdata2[col].cat.set_categories(['Yes', 'No', 'NA'])
    else:
        verifytestdata2[col] = verifytestdata2[col].cat.set_categories(['0', '1'])


# 预测独立的测试集和验证集数据
y_verify_test_pred = model.predict(verifytestdata2)
y_verify_test_pred_proba = model.predict_proba(verifytestdata2)[:, 1]
# 计算验证集和独立测试集的指标
verify_test_metrics = calculate_metrics(y_verify_test, y_verify_test_pred, y_verify_test_pred_proba)
# 将结果放入DataFrame
df = pd.DataFrame([verify_test_metrics],index=['XGBoost_test_cgmlst_independent_inOwnValid'])
# 保存文件
df.to_csv("RS_811/14880file_ARGanno_absence01_SNPmatrix959090_YesNoNA_addPhenotype_CIP_verify_test_cgmlst_independent_inOwnValid_finalmodel_XGBoostresult.txt", sep='\t', index=True, header=True)

### 输出测试集用matrix文件（真实值 + 预测值 + 预测概率）
output_df = verifytestdata2.copy()
output_df['y_true'] = y_verify_test          # 真实标签 R=1, S=0
output_df['y_pred'] = y_verify_test_pred    # 模型预测结果
output_df['y_pred_prob'] = y_verify_test_pred_proba  # 预测为1（R）的概率
output_df.to_csv("RS_811/14880file_ARGanno_absence01_SNPmatrix959090_YesNoNA_addPhenotype_CIP_verify_test_cgmlst_independent_inOwnValid_finalmodel_usematrix.txt",sep="\t",header=True,index=True)

