# 导入必要的库
#from sklearn.model_selection import GridSearchCV, KFold
#from sklearn.model_selection import StratifiedKFold
#from sklearn.metrics import accuracy_score
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
new_working_directory = "~/Feature_analysis/"
# 更改工作目录
os.chdir(new_working_directory)


## 训练集
traindata=pd.read_csv("~/ARG_MIC/RS_811/14880file_ARGanno_Protein_result_absence01_addPhenotype_AMP_train_stcgmlst.txt",sep='\t', header=0,index_col=0)
# 提取特征矩阵 X 和目标变量 y
X_train=traindata.iloc[:, :-1]# 选择除最后一列之外的所有列作为特征
y_train=traindata.iloc[:, -1]# 选择最后一列作为目标变量
# 将字符串标签转换为整数标签
y_train=y_train.replace({'R': 1, 'S': 0})

## 测试集
verifytestdata=pd.read_csv("~/ARG_MIC/RS_811/14880file_ARGanno_Protein_result_absence01_addPhenotype_AMP_verify_test_stcgmlst.txt",sep='\t', header=0,index_col=0)
# 将字符串标签转换为整数标签
X_verify_test=verifytestdata.iloc[:, :-1]
y_verify_test=verifytestdata.iloc[:, -1]
y_verify_test=y_verify_test.replace({'R': 1, 'S': 0})


### selectFeature
df_sorted=pd.read_csv("~/ARG_MIC/AllFeature/XGBoost_10fold_crossvalidation_default_hyperparameter_ARG_absence01_AMP_fold_3_feature_importances_verifytest_stcgmlst.txt",sep='\t', header=0,index_col=0)
df_filtered = df_sorted[df_sorted['Importance'] > 0]
df_filtered_top_index0 = df_filtered.head(30).index
X_train2=X_train.iloc[:,df_filtered_top_index0]
X_verify_test2=X_verify_test.iloc[:,df_filtered_top_index0]

#final_selected_features_index = df_filtered.head(30)['Feature'].tolist() + ['Group']
#traindata2=traindata.iloc[:,traindata.columns.isin(final_selected_features_index)]
#verifytestdata2=verifytestdata.iloc[:,verifytestdata.columns.isin(final_selected_features_index)]
final_selected_features_index = X_train2.columns.tolist() + ['Group']
traindata2=traindata.loc[:,final_selected_features_index]
verifytestdata2=verifytestdata.loc[:,final_selected_features_index]

# 输出最终用模型的训练集测试集
traindata2.to_csv("~/Feature_analysis/use_topfeature/XGBoost_AMP_topfeature_finalmodel_traindata_stcgmlst.txt", sep='\t', index=True, header=True)
verifytestdata2.to_csv("~/Feature_analysis/use_topfeature/XGBoost_AMP_topfeature_finalmodel_verifytestdata_stcgmlst.txt", sep='\t', index=True, header=True)



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


# 初始化并训练模型--XGBoost
model = XGBClassifier(learning_rate=0.01,max_depth=9,n_estimators=700,random_state=2)
model.fit(X_train2, y_train)
# 预测测试集数据
y_verify_test_pred = model.predict(X_verify_test2)
y_verify_test_pred_proba = model.predict_proba(X_verify_test2)[:, 1]
# 计算验证集和独立测试集的指标
verify_test_metrics = calculate_metrics(y_verify_test, y_verify_test_pred, y_verify_test_pred_proba)
# 将结果放入DataFrame
df = pd.DataFrame([verify_test_metrics],index=['XGBoost_verify_test'])
# 保存文件
df.to_csv("~/Feature_analysis/use_topfeature/XGBoost_AMP_topfeature_finalmodel_evaluation_indicator_results_stcgmlst.txt", sep='\t', index=True, header=True)


#### 输出各特征重要性
importances = model.feature_importances_
# 输出各特征的重要性值
df0 = pd.DataFrame({
    'Feature': X_train2.columns,
    'Importance': importances
})
df0_sorted = df0.sort_values(by='Importance', ascending=False)
df0_sorted.to_csv("~/Feature_analysis/use_topfeature/XGBoost_AMP_topfeature_finalmodel_feature_importances_stcgmlst.txt", sep='\t', index=True, header=True)


###################################################### 默认参数 ###############################
# 初始化并训练模型--XGBoost
model = XGBClassifier(random_state=2)
model.fit(X_train2, y_train)
# 预测测试集数据
y_verify_test_pred = model.predict(X_verify_test2)
y_verify_test_pred_proba = model.predict_proba(X_verify_test2)[:, 1]
# 计算验证集和独立测试集的指标
verify_test_metrics = calculate_metrics(y_verify_test, y_verify_test_pred, y_verify_test_pred_proba)
# 将结果放入DataFrame
df = pd.DataFrame([verify_test_metrics],index=['XGBoost_verify_test'])
# 保存文件
df.to_csv("~/Feature_analysis/use_topfeature2/XGBoost_AMP_topfeature_finalmodel_evaluation_indicator_results_stcgmlst.txt", sep='\t', index=True, header=True)

#### 输出各特征重要性
importances = model.feature_importances_
# 输出各特征的重要性值
df0 = pd.DataFrame({
    'Feature': X_train2.columns,
    'Importance': importances
})
df0_sorted = df0.sort_values(by='Importance', ascending=False)
df0_sorted.to_csv("~/Feature_analysis/use_topfeature2/XGBoost_AMP_topfeature_finalmodel_feature_importances_stcgmlst.txt", sep='\t', index=True, header=True)

