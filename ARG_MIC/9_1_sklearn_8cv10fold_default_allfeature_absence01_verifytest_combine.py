# 导入必要的库
from sklearn.model_selection import GridSearchCV, KFold
#from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn import metrics
from xgboost import XGBClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
#from sklearn.inspection import permutation_importance
import math
import numpy as np
import pandas as pd ##读文件
import os ##改工作目录
import sys ##传参
from sklearn.metrics import make_scorer, accuracy_score, roc_auc_score, f1_score, recall_score, precision_score, cohen_kappa_score, balanced_accuracy_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

# 传参
arguments = sys.argv
# 设置新的工作目录路径
new_working_directory = "~/ARG_MIC/"
# 更改工作目录
os.chdir(new_working_directory)


## 外部的总测试集
verifytestdata=pd.read_csv("RS_811/14880file_ARGanno_Protein_result_absence01_addPhenotype_" + arguments[1] + "_verify_test_stcgmlst.txt",sep='\t', header=0,index_col=0)
# 将字符串标签转换为整数标签
X_verify_test=verifytestdata.iloc[:, :-1]
y_verify_test=verifytestdata.iloc[:, -1]
y_verify_test=y_verify_test.replace({'R': 1, 'S': 0})


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


# 初始化存储结果的列表
all_metrics = []

for fold in range(1,11):
    ## 10倍交叉验证中的训练集验证集
    fold_traindata=pd.read_csv(f"AllFeature/XGBoost_ARG_absence01_fold_{fold}_train_samples_10fold_" + arguments[1] + "_stcgmlst.txt",sep='\t', header=0,index_col=0)
    fold_valdata=pd.read_csv(f"AllFeature/XGBoost_ARG_absence01_fold_{fold}_val_samples_10fold_" + arguments[1] + "_stcgmlst.txt",sep='\t', header=0,index_col=0)
    # 提取特征矩阵 X 和目标变量 y
    X_fold_train=fold_traindata.iloc[:, :-1]# 选择除最后一列之外的所有列作为特征
    y_fold_train=fold_traindata.iloc[:, -1]# 选择最后一列作为目标变量
    X_fold_val=fold_valdata.iloc[:, :-1]# 选择除最后一列之外的所有列作为特征
    y_fold_val=fold_valdata.iloc[:, -1]# 选择最后一列作为目标变量

    # 初始化并训练模型--XGBoost
    model = XGBClassifier(random_state=2)
    model.fit(X_fold_train, y_fold_train)

    # 预测验证集数据
    y_val_pred = model.predict(X_fold_val)
    # 预测独立的测试集和验证集数据
    y_verify_test_pred = model.predict(X_verify_test)
    y_val_pred_proba = model.predict_proba(X_fold_val)[:, 1]
    y_verify_test_pred_proba = model.predict_proba(X_verify_test)[:, 1]
    # 计算验证集和独立测试集的指标
    val_metrics = calculate_metrics(y_fold_val, y_val_pred, y_val_pred_proba)
    verify_test_metrics = calculate_metrics(y_verify_test, y_verify_test_pred, y_verify_test_pred_proba)
    # 将结果放入DataFrame
    fold_results = pd.DataFrame([val_metrics, verify_test_metrics],index=[f'XGBoost_Fold_{fold}_val', f'XGBoost_Fold_{fold}_verify_test'])

    all_metrics.append(fold_results)

# 将所有结果合并为一个DataFrame
all_metrics_df = pd.concat(all_metrics)
# 保存文件
all_metrics_df.to_csv("AllFeature/XGBoost_10fold_crossvalidation_default_hyperparameter_ARG_absence01_"+ arguments[1] +"_verifytest_stcgmlst.txt", sep='\t', index=True, header=True)

