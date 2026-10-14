# 导入必要的库
from pandas.api.types import CategoricalDtype
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
import pandas as pd ##读文件
import os ##改工作目录
import sys ##传参
from sklearn.metrics import make_scorer, accuracy_score, roc_auc_score, f1_score, recall_score, precision_score, cohen_kappa_score, balanced_accuracy_score

# 传参
arguments = sys.argv
# 设置新的工作目录路径
new_working_directory = "~/ARG_with_SNP/"
# 更改工作目录
os.chdir(new_working_directory)


## 训练集验证集
traindata=pd.read_csv("RS_811/14880file_ARGanno_absence01_SNPmatrix959090_YesNoNA_addPhenotype_" + arguments[1] + "_train_stcgmlst.txt",sep='\t', header=0,index_col=0, keep_default_na=False, dtype=CategoricalDtype())
verify_testdata=pd.read_csv("RS_811/14880file_ARGanno_absence01_SNPmatrix959090_YesNoNA_addPhenotype_" + arguments[1] + "_verify_test_stcgmlst.txt",sep='\t', header=0,index_col=0, keep_default_na=False, dtype=CategoricalDtype())
# 提取特征矩阵 X 和目标变量 y
X_train=traindata.iloc[:, :-1]# 选择除最后一列之外的所有列作为特征
# 获取类别型列的索引
category_columns1 = list(range(0, 897)) # 第1到第897列
category_columns2 = list(range(898, X_train.shape[1])) # 第898列到最后列
# 将ARGanno_absence01指定01索引
for col in X_train.columns[category_columns1]:
    X_train[col] = X_train[col].cat.set_categories(['0', '1'])

# 将SNPmatrix959090_YesNoNA指定YesNoNA索引
for col in X_train.columns[category_columns2]:
    X_train[col] = X_train[col].cat.set_categories(['Yes', 'No', 'NA'])

y_train=traindata.iloc[:, -1]# 选择最后一列作为目标变量
# 将字符串标签转换为整数标签
y_train=y_train.replace({'R': 1, 'S': 0})
y_train=y_train.astype('int')
#for col in X_train.columns:
#       X_train[col] = X_train[col].astype('category')
X_verify_test=verify_testdata.iloc[:, :-1]
for col in X_verify_test.columns[category_columns1]:
    X_verify_test[col] = X_verify_test[col].cat.set_categories(['0', '1'])

for col in X_verify_test.columns[category_columns2]:
    X_verify_test[col] = X_verify_test[col].cat.set_categories(['Yes', 'No', 'NA'])

y_verify_test=verify_testdata.iloc[:, -1]
y_verify_test=y_verify_test.replace({'R': 1, 'S': 0})
y_verify_test=y_verify_test.astype('int')
#for col in X_verify.columns:
#        X_verify[col] = X_verify[col].astype('category')


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
    fold_traindata=pd.read_csv(f"AllFeature_YesNoNA/XGBoost_ARGanno_absence01_SNPmatrix959090_YesNoNA_fold_{fold}_train_samples_10fold_" + arguments[1] + "_stcgmlst.txt",sep='\t', header=0,index_col=0, keep_default_na=False, dtype=CategoricalDtype())
    fold_valdata=pd.read_csv(f"AllFeature_YesNoNA/XGBoost_ARGanno_absence01_SNPmatrix959090_YesNoNA_fold_{fold}_val_samples_10fold_" + arguments[1] + "_stcgmlst.txt",sep='\t', header=0,index_col=0, keep_default_na=False, dtype=CategoricalDtype())

    # 提取特征矩阵 X 和目标变量 y
    X_fold_train=fold_traindata.iloc[:, :-1]# 选择除最后一列之外的所有列作为特征
    # 获取类别型列的索引
    category_columns1 = list(range(0, 897)) # 第1到第897列
    category_columns2 = list(range(898, X_fold_train.shape[1])) # 第898列到最后列
    # 将ARGanno_absence01指定01索引
    for col in X_fold_train.columns[category_columns1]:
        X_fold_train[col] = X_fold_train[col].cat.set_categories(['0', '1'])
    
    # 将SNPmatrix959090_YesNoNA指定YesNoNA索引
    for col in X_fold_train.columns[category_columns2]:
        X_fold_train[col] = X_fold_train[col].cat.set_categories(['Yes', 'No', 'NA'])
    
    y_fold_train=fold_traindata.iloc[:, -1]# 选择最后一列作为目标变量
    # 将字符串标签转换为整数标签
    y_fold_train=y_fold_train.replace({'R': 1, 'S': 0})
    y_fold_train=y_fold_train.astype('int')
    
    X_fold_val=fold_valdata.iloc[:, :-1]
    for col in X_fold_val.columns[category_columns1]:
        X_fold_val[col] = X_fold_val[col].cat.set_categories(['0', '1'])
    
    for col in X_fold_val.columns[category_columns2]:
        X_fold_val[col] = X_fold_val[col].cat.set_categories(['Yes', 'No', 'NA'])
    
    y_fold_val=fold_valdata.iloc[:, -1]
    y_fold_val=y_fold_val.replace({'R': 1, 'S': 0})
    y_fold_val=y_fold_val.astype('int')



    # 初始化并训练模型--XGBoost
    model = XGBClassifier(random_state=2, enable_categorical=True)
    model.fit(X_fold_train, y_fold_train)
    
    # 预测验证集数据
    y_val_pred = model.predict(X_fold_val)
    # 预测独立的测试集和验证集数据
    y_verify_test_pred = model.predict(X_verify_test)
    y_val_pred_proba = model.predict_proba(X_fold_val)[:, 1]
    y_verify_test_pred_proba = model.predict_proba(X_verify_test)[:, 1]
    
    # 计算每折的验证集、独立验证集、测试集指标
    val_metrics = calculate_metrics(y_fold_val, y_val_pred, y_val_pred_proba)
    verify_test_metrics = calculate_metrics(y_verify_test, y_verify_test_pred, y_verify_test_pred_proba)
    
    # 将结果存入DataFrame并标记折数
    fold_results = pd.DataFrame([val_metrics, verify_test_metrics], 
                                index=[f'XGBoost_Fold_{fold+1}_val', f'XGBoost_Fold_{fold+1}_verify_test'])
    
    all_metrics.append(fold_results)

# 将所有结果合并为一个DataFrame
all_metrics_df = pd.concat(all_metrics)

# 保存文件
all_metrics_df.to_csv("AllFeature_YesNoNA/XGBoost_10fold_crossvalidation_default_hyperparameter_ARGanno_absence01_SNPmatrix959090_YesNoNA_"+ arguments[1] +"_stcgmlst.txt", sep='\t', index=True, header=True)





############### RandomForest模型 #############
#model1 = RandomForestClassifier(random_state=2)
#model1.fit(X_train,y_train)
###### 特征重要性
##importances1 = model1.feature_importances_
##print("Feature importances:", importances1)
### 设置一个阈值
##threshold1 = 0.15
### 选择重要性大于阈值的特征
##selected_features1 = np.where(importances1 > threshold1)[0]
##selected_X_train = X_train[:, selected_features1]
##print("Selected feature indices:", selected_features1)
#
## 预测验证集数据
#y_verify_test_pred1 = model1.predict(X_verify_test)
## 预测测试集数据
#y_test_pred1 = model1.predict(X_test)
#
## 计算验证集和独立测试集的指标
#verify_test_metrics1 = calculate_metrics(y_verify_test, y_verify_test_pred1)
#test_metrics1 = calculate_metrics(y_test, y_test_pred1)
## 将两组指标放入DataFrame
#df1 = pd.DataFrame([verify_test_metrics1, test_metrics1], index=['RandomForest_verify_test', 'RandomForest_test'])
## 合并XGBoost和SVM结果
#combined_df = pd.concat([df, df1], axis=0)
## 保存文件
#combined_df.to_csv("AllFeature_YesNoNA/default_hyperparameter_ARGanno_absence01_SNPmatrix959090_YesNoNA_"+ arguments[1] +".csv", sep='\t', index=True, header=True)
#
#
############### SVM模型 #############
#model2 = SVC(random_state=3)
#model2.fit(X_train,y_train)
#####随机打乱（去掉）一个特征后看跟基线的差距，差距大的说明特征重要
##result = permutation_importance(model, X_train, y_train, n_repeats=5, random_state=42, scoring='balanced_accuracy')
##importances = result.importances_mean
###从大到小排列的索引值np.argsort(-abs(importances))
##
##### 输出各特征的重要性值
##np.sav_stcgmlst.txt("Feature_select/SVM_rbfC50Gamma0.01_permutationImportance_"+ arguments +".csv", importances,delimiter=',')
##importances.to_csv("Feature_select/SVM_rbfC50Gamma0.01_permutationImportance_"+ arguments[1] +".csv", index=False)
##
#### 打印特征的重要性
###for i in results.importances_mean.argsort()[::-1]:
###    if results.importances_mean[i] - 2 * results.importances_std[i] > 0:
###        print(f"{iris.feature_names[i]}: {results.importances_mean[i]:.3f} +/- {results.importances_std[i]:.3f}")
##
###importances.to_csv("Feature_select/SVM_rbfC50Gamma0.01_permutationImportance_"+ arguments[1] +".csv", index=False)
###results_svm = pd.DataFrame(svm_grid.cv_results_)
###results_svm.to_csv("SVM_grid_search_results2_"+ arguments[1] +".csv", index=False)
#
## 预测验证集数据
#y_verify_test_pred2 = model2.predict(X_verify_test)
## 预测测试集数据
#y_test_pred2 = model2.predict(X_test)
#
## 计算验证集和独立测试集的指标
#verify_test_metrics2 = calculate_metrics(y_verify_test, y_verify_test_pred2)
#test_metrics2 = calculate_metrics(y_test, y_test_pred2)
## 将两组指标放入DataFrame
#df2 = pd.DataFrame([verify_test_metrics2, test_metrics2], index=['SVM_verify_test', 'SVM_test'])
## 合并XGBoost和SVM结果
#combined_df2 = pd.concat([combined_df, df2], axis=0)
## 保存文件
#combined_df2.to_csv("AllFeature_YesNoNA/default_hyperparameter_ARGanno_absence01_SNPmatrix959090_YesNoNA_"+ arguments[1] +".csv", sep='\t', index=True, header=True)
#
