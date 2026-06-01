# 导入必要的库
from pandas.api.types import CategoricalDtype
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
import pandas as pd ##读文件
import os ##改工作目录
import sys ##传参
from sklearn.metrics import make_scorer, accuracy_score, roc_auc_score, f1_score, recall_score, precision_score, cohen_kappa_score, balanced_accuracy_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

# 传参
arguments = sys.argv
# 设置新的工作目录路径
new_working_directory = "/lustre/home/niejingyi2023/project/bjCDC/ARG_with_SNP/"
# 更改工作目录
os.chdir(new_working_directory)

traindata=pd.read_csv("RS_811/14880file_ARGanno_absence01_SNPmatrix959090_absence01_addPhenotype_" + arguments[1] + "_train.txt",sep='\t', header=0,index_col=0, keep_default_na=False, dtype=CategoricalDtype())
verifytestdata=pd.read_csv("RS_811/14880file_ARGanno_absence01_SNPmatrix959090_absence01_addPhenotype_" + arguments[1] + "_verify_test.txt",sep='\t', header=0,index_col=0, keep_default_na=False, dtype=CategoricalDtype())

# 提取特征矩阵 X 和目标变量 y
X_train=traindata.iloc[:, :-1]# 选择除最后一列之外的所有列作为特征
for col in X_train.columns:
    if pd.api.types.is_categorical_dtype(X_train[col]):
        X_train[col] = X_train[col].cat.set_categories(['0', '1'])

y_train=traindata.iloc[:, -1]# 选择最后一列作为目标变量
# 将字符串标签转换为整数标签
y_train=y_train.replace({'R': 1, 'S': 0})
y_train=y_train.astype('int')
#for col in X_train.columns:
#       X_train[col] = X_train[col].astype('category')
X_verify_test=verifytestdata.iloc[:, :-1]
for col in X_verify_test.columns:
    if pd.api.types.is_categorical_dtype(X_verify_test[col]):
        X_verify_test[col] = X_verify_test[col].cat.set_categories(['0', '1'])

y_verify_test=verifytestdata.iloc[:, -1]
y_verify_test=y_verify_test.replace({'R': 1, 'S': 0})
y_verify_test=y_verify_test.astype('int')
#for col in X_verify.columns:
#        X_verify[col] = X_verify[col].astype('category')


# 10折交叉验证
skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=1)

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
X_fold_train, X_fold_val = None, None
y_fold_train, y_fold_val = None, None
strnum=int(arguments[2])-1

# 交叉验证
for fold, (train_idx, val_idx) in enumerate(skf.split(X_train, y_train)):
    # 获取当前折的训练和验证数据
    if fold==strnum:
        X_fold_train, X_fold_val = X_train.iloc[train_idx], X_train.iloc[val_idx]
        y_fold_train, y_fold_val = y_train.iloc[train_idx], y_train.iloc[val_idx]

if X_fold_train is None or y_fold_train is None:
    print("Not Found Fold or undefined training data")

# 初始化并训练模型--XGBoost
model = XGBClassifier(random_state=2, enable_categorical=True)
model.fit(X_fold_train, y_fold_train)
#### 特征重要性
importances = model.feature_importances_
# 输出各特征的重要性值
df0 = pd.DataFrame({
    'Feature': X_fold_train.columns,
    'Importance': importances
})
df0_sorted = df0.sort_values(by='Importance', ascending=False)
df0_filtered = df0_sorted[df0_sorted['Importance'] > 0]
df0_sorted.to_csv("AllFeature/XGBoost_10fold_crossvalidation_default_hyperparameter_ARGanno_absence01_SNPmatrix959090_absence01_"+ arguments[1] +"_fold_" + arguments[2] + "_feature_importances_category_and_combine_verifytest.txt", sep='\t', index=True, header=True)


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
df = pd.DataFrame([val_metrics, verify_test_metrics],index=['XGBoost_val_featureall', 'XGBoost_verify_test_featureall'])
# 保存文件
df.to_csv("AllFeature/XGBoost_10fold_crossvalidation_default_hyperparameter_ARGanno_absence01_SNPmatrix959090_absence01_"+ arguments[1] +"_fold_" + arguments[2] + "_topfeature_category_and_combine_verifytest.txt", sep='\t', index=True, header=True)


########################### 取top5特征训练模型 ###################
df0_filtered_top=df0_filtered.head(5)
## refilt data
X_fold_train2=X_fold_train.iloc[:,df0_filtered_top.index]
X_fold_val2=X_fold_val.iloc[:,df0_filtered_top.index]
X_verify_test2=X_verify_test.iloc[:,df0_filtered_top.index]
X_train2=X_train.iloc[:,df0_filtered_top.index]
## 使用特征过滤后数据重新建模，不需要再重新交叉验证了
# 初始化并训练模型--XGBoost
model = XGBClassifier(random_state=2, enable_categorical=True)
model.fit(X_fold_train2, y_fold_train)
# 预测验证集数据
y_val_pred = model.predict(X_fold_val2)
# 预测独立的测试集和验证集数据
y_verify_test_pred = model.predict(X_verify_test2)
y_val_pred_proba = model.predict_proba(X_fold_val2)[:, 1]
y_verify_test_pred_proba = model.predict_proba(X_verify_test2)[:, 1]
# 计算验证集和独立测试集的指标
val_metrics = calculate_metrics(y_fold_val, y_val_pred, y_val_pred_proba)
verify_test_metrics = calculate_metrics(y_verify_test, y_verify_test_pred, y_verify_test_pred_proba)
# 将结果放入DataFrame
df5 = pd.DataFrame([val_metrics, verify_test_metrics],index=['XGBoost_val_feature5', 'XGBoost_verify_test_feature5'])
# 合并各类特征数目模型结果
combined_df = pd.concat([df, df5], axis=0)
# 保存文件
combined_df.to_csv("AllFeature/XGBoost_10fold_crossvalidation_default_hyperparameter_ARGanno_absence01_SNPmatrix959090_absence01_"+ arguments[1] +"_fold_" + arguments[2] + "_topfeature_category_and_combine_verifytest.txt", sep='\t', index=True, header=True)

######## 寻找最优超参数 #######
# 定义十倍交叉验证
cv = KFold(n_splits=10, shuffle=True, random_state=42)
# 自定义Prevalence评估函数
def prevalence_score(y_true, y_pred):
    return sum(y_true) / len(y_true)
# 定义评分函数
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
# 定义评估指标
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
           'Prevalence': make_scorer(prevalence_score)}  # 使用自定义的Prevalence评估函数
# 定义XGBoost模型及参数
xgb_params = {'n_estimators': [100,200,300,500,700,1000,1500],
              'max_depth': [3,4,5,6,7,8,9,10],
              'learning_rate': [0.01, 0.1, 0.3]#,
        }
xgb_model = XGBClassifier(random_state=2, enable_categorical=True)
# 使用GridSearchCV进行交叉验证和参数搜索
xgb_grid = GridSearchCV(xgb_model, param_grid=xgb_params, scoring=scoring,refit='Accuracy', cv=cv)
xgb_grid.fit(X_train2, y_train)
# 打印最佳参数和交叉验证分数
print("XGBoost Best Parameters:", xgb_grid.best_params_)
print("XGBoost Best Accuracy:", xgb_grid.best_score_)
results_XGBoost = pd.DataFrame(xgb_grid.cv_results_)
results_XGBoost.to_csv("AllFeature/XGBoost_10fold_crossvalidation_default_hyperparameter_ARGanno_absence01_SNPmatrix959090_absence01_"+ arguments[1] +"_fold_" + arguments[2] + "_top5feature_category_and_combine_verifytest.txt", sep='\t', index=False, header=True)


########################### 取top10特征训练模型 ###################
df0_filtered_top=df0_filtered.head(10)
## refilt data
X_fold_train2=X_fold_train.iloc[:,df0_filtered_top.index]
X_fold_val2=X_fold_val.iloc[:,df0_filtered_top.index]
X_verify_test2=X_verify_test.iloc[:,df0_filtered_top.index]
X_train2=X_train.iloc[:,df0_filtered_top.index]
## 使用特征过滤后数据重新建模，不需要再重新交叉验证了
# 初始化并训练模型--XGBoost
model = XGBClassifier(random_state=2, enable_categorical=True)
model.fit(X_fold_train2, y_fold_train)
# 预测验证集数据
y_val_pred = model.predict(X_fold_val2)
# 预测独立的测试集和验证集数据
y_verify_test_pred = model.predict(X_verify_test2)
y_val_pred_proba = model.predict_proba(X_fold_val2)[:, 1]
y_verify_test_pred_proba = model.predict_proba(X_verify_test2)[:, 1]
# 计算验证集和独立测试集的指标
val_metrics = calculate_metrics(y_fold_val, y_val_pred, y_val_pred_proba)
verify_test_metrics = calculate_metrics(y_verify_test, y_verify_test_pred, y_verify_test_pred_proba)
# 将结果放入DataFrame
df10 = pd.DataFrame([val_metrics, verify_test_metrics],index=['XGBoost_val_feature10', 'XGBoost_verify_test_feature10'])
# 合并各类特征数目模型结果
combined_df = pd.concat([combined_df, df10], axis=0)
# 保存文件
combined_df.to_csv("AllFeature/XGBoost_10fold_crossvalidation_default_hyperparameter_ARGanno_absence01_SNPmatrix959090_absence01_"+ arguments[1] +"_fold_" + arguments[2] + "_topfeature_category_and_combine_verifytest.txt", sep='\t', index=True, header=True)

######## 寻找最优超参数 #######
xgb_model = XGBClassifier(random_state=2, enable_categorical=True)
# 使用GridSearchCV进行交叉验证和参数搜索
xgb_grid = GridSearchCV(xgb_model, param_grid=xgb_params, scoring=scoring,refit='Accuracy', cv=cv)
xgb_grid.fit(X_train2, y_train)
# 打印最佳参数和交叉验证分数
print("XGBoost Best Parameters:", xgb_grid.best_params_)
print("XGBoost Best Accuracy:", xgb_grid.best_score_)
results_XGBoost = pd.DataFrame(xgb_grid.cv_results_)
results_XGBoost.to_csv("AllFeature/XGBoost_10fold_crossvalidation_default_hyperparameter_ARGanno_absence01_SNPmatrix959090_absence01_"+ arguments[1] +"_fold_" + arguments[2] + "_top10feature_category_and_combine_verifytest.txt", sep='\t', index=False, header=True)

########################### 取top20特征训练模型 ###################
df0_filtered_top=df0_filtered.head(20)
## refilt data
X_fold_train2=X_fold_train.iloc[:,df0_filtered_top.index]
X_fold_val2=X_fold_val.iloc[:,df0_filtered_top.index]
X_verify_test2=X_verify_test.iloc[:,df0_filtered_top.index]
X_train2=X_train.iloc[:,df0_filtered_top.index]
## 使用特征过滤后数据重新建模，不需要再重新交叉验证了
# 初始化并训练模型--XGBoost
model = XGBClassifier(random_state=2, enable_categorical=True)
model.fit(X_fold_train2, y_fold_train)
# 预测验证集数据
y_val_pred = model.predict(X_fold_val2)
# 预测独立的测试集和验证集数据
y_verify_test_pred = model.predict(X_verify_test2)
y_val_pred_proba = model.predict_proba(X_fold_val2)[:, 1]
y_verify_test_pred_proba = model.predict_proba(X_verify_test2)[:, 1]
# 计算验证集和独立测试集的指标
val_metrics = calculate_metrics(y_fold_val, y_val_pred, y_val_pred_proba)
verify_test_metrics = calculate_metrics(y_verify_test, y_verify_test_pred, y_verify_test_pred_proba)
# 将结果放入DataFrame
df20 = pd.DataFrame([val_metrics, verify_test_metrics],index=['XGBoost_val_feature20', 'XGBoost_verify_test_feature20'])
# 合并各类特征数目模型结果
combined_df = pd.concat([combined_df, df20], axis=0)
# 保存文件
combined_df.to_csv("AllFeature/XGBoost_10fold_crossvalidation_default_hyperparameter_ARGanno_absence01_SNPmatrix959090_absence01_"+ arguments[1] +"_fold_" + arguments[2] + "_topfeature_category_and_combine_verifytest.txt", sep='\t', index=True, header=True)

######## 寻找最优超参数 #######
xgb_model = XGBClassifier(random_state=2, enable_categorical=True)
# 使用GridSearchCV进行交叉验证和参数搜索
xgb_grid = GridSearchCV(xgb_model, param_grid=xgb_params, scoring=scoring,refit='Accuracy', cv=cv)
xgb_grid.fit(X_train2, y_train)
# 打印最佳参数和交叉验证分数
print("XGBoost Best Parameters:", xgb_grid.best_params_)
print("XGBoost Best Accuracy:", xgb_grid.best_score_)
results_XGBoost = pd.DataFrame(xgb_grid.cv_results_)
results_XGBoost.to_csv("AllFeature/XGBoost_10fold_crossvalidation_default_hyperparameter_ARGanno_absence01_SNPmatrix959090_absence01_"+ arguments[1] +"_fold_" + arguments[2] + "_top20feature_category_and_combine_verifytest.txt", sep='\t', index=False, header=True)

########################### 取top30特征训练模型 ###################
df0_filtered_top=df0_filtered.head(30)
## refilt data
X_fold_train2=X_fold_train.iloc[:,df0_filtered_top.index]
X_fold_val2=X_fold_val.iloc[:,df0_filtered_top.index]
X_verify_test2=X_verify_test.iloc[:,df0_filtered_top.index]
X_train2=X_train.iloc[:,df0_filtered_top.index]
## 使用特征过滤后数据重新建模，不需要再重新交叉验证了
# 初始化并训练模型--XGBoost
model = XGBClassifier(random_state=2, enable_categorical=True)
model.fit(X_fold_train2, y_fold_train)
# 预测验证集数据
y_val_pred = model.predict(X_fold_val2)
# 预测独立的测试集和验证集数据
y_verify_test_pred = model.predict(X_verify_test2)
y_val_pred_proba = model.predict_proba(X_fold_val2)[:, 1]
y_verify_test_pred_proba = model.predict_proba(X_verify_test2)[:, 1]
# 计算验证集和独立测试集的指标
val_metrics = calculate_metrics(y_fold_val, y_val_pred, y_val_pred_proba)
verify_test_metrics = calculate_metrics(y_verify_test, y_verify_test_pred, y_verify_test_pred_proba)
# 将结果放入DataFrame
df30 = pd.DataFrame([val_metrics, verify_test_metrics],index=['XGBoost_val_feature30', 'XGBoost_verify_test_feature30'])
# 合并各类特征数目模型结果
combined_df = pd.concat([combined_df, df30], axis=0)
# 保存文件
combined_df.to_csv("AllFeature/XGBoost_10fold_crossvalidation_default_hyperparameter_ARGanno_absence01_SNPmatrix959090_absence01_"+ arguments[1] +"_fold_" + arguments[2] + "_topfeature_category_and_combine_verifytest.txt", sep='\t', index=True, header=True)

######## 寻找最优超参数 #######
xgb_model = XGBClassifier(random_state=2, enable_categorical=True)
# 使用GridSearchCV进行交叉验证和参数搜索
xgb_grid = GridSearchCV(xgb_model, param_grid=xgb_params, scoring=scoring,refit='Accuracy', cv=cv)
xgb_grid.fit(X_train2, y_train)
# 打印最佳参数和交叉验证分数
print("XGBoost Best Parameters:", xgb_grid.best_params_)
print("XGBoost Best Accuracy:", xgb_grid.best_score_)
results_XGBoost = pd.DataFrame(xgb_grid.cv_results_)
results_XGBoost.to_csv("AllFeature/XGBoost_10fold_crossvalidation_default_hyperparameter_ARGanno_absence01_SNPmatrix959090_absence01_"+ arguments[1] +"_fold_" + arguments[2] + "_top30feature_category_and_combine_verifytest.txt", sep='\t', index=False, header=True)

########################### 取top40特征训练模型 ###################
df0_filtered_top=df0_filtered.head(40)
## refilt data
X_fold_train2=X_fold_train.iloc[:,df0_filtered_top.index]
X_fold_val2=X_fold_val.iloc[:,df0_filtered_top.index]
X_verify_test2=X_verify_test.iloc[:,df0_filtered_top.index]
X_train2=X_train.iloc[:,df0_filtered_top.index]
## 使用特征过滤后数据重新建模，不需要再重新交叉验证了
# 初始化并训练模型--XGBoost
model = XGBClassifier(random_state=2, enable_categorical=True)
model.fit(X_fold_train2, y_fold_train)
# 预测验证集数据
y_val_pred = model.predict(X_fold_val2)
# 预测独立的测试集和验证集数据
y_verify_test_pred = model.predict(X_verify_test2)
y_val_pred_proba = model.predict_proba(X_fold_val2)[:, 1]
y_verify_test_pred_proba = model.predict_proba(X_verify_test2)[:, 1]
# 计算验证集和独立测试集的指标
val_metrics = calculate_metrics(y_fold_val, y_val_pred, y_val_pred_proba)
verify_test_metrics = calculate_metrics(y_verify_test, y_verify_test_pred, y_verify_test_pred_proba)
# 将结果放入DataFrame
df40 = pd.DataFrame([val_metrics, verify_test_metrics],index=['XGBoost_val_feature40', 'XGBoost_verify_test_feature40'])
# 合并各类特征数目模型结果
combined_df = pd.concat([combined_df, df40], axis=0)
# 保存文件
combined_df.to_csv("AllFeature/XGBoost_10fold_crossvalidation_default_hyperparameter_ARGanno_absence01_SNPmatrix959090_absence01_"+ arguments[1] +"_fold_" + arguments[2] + "_topfeature_category_and_combine_verifytest.txt", sep='\t', index=True, header=True)

######## 寻找最优超参数 #######
xgb_model = XGBClassifier(random_state=2, enable_categorical=True)
# 使用GridSearchCV进行交叉验证和参数搜索
xgb_grid = GridSearchCV(xgb_model, param_grid=xgb_params, scoring=scoring,refit='Accuracy', cv=cv)
xgb_grid.fit(X_train2, y_train)
# 打印最佳参数和交叉验证分数
print("XGBoost Best Parameters:", xgb_grid.best_params_)
print("XGBoost Best Accuracy:", xgb_grid.best_score_)
results_XGBoost = pd.DataFrame(xgb_grid.cv_results_)
results_XGBoost.to_csv("AllFeature/XGBoost_10fold_crossvalidation_default_hyperparameter_ARGanno_absence01_SNPmatrix959090_absence01_"+ arguments[1] +"_fold_" + arguments[2] + "_top40feature_category_and_combine_verifytest.txt", sep='\t', index=False, header=True)

########################### 取top50特征训练模型 ###################
df0_filtered_top=df0_filtered.head(50)
## refilt data
X_fold_train2=X_fold_train.iloc[:,df0_filtered_top.index]
X_fold_val2=X_fold_val.iloc[:,df0_filtered_top.index]
X_verify_test2=X_verify_test.iloc[:,df0_filtered_top.index]
X_train2=X_train.iloc[:,df0_filtered_top.index]
## 使用特征过滤后数据重新建模，不需要再重新交叉验证了
# 初始化并训练模型--XGBoost
model = XGBClassifier(random_state=2, enable_categorical=True)
model.fit(X_fold_train2, y_fold_train)
# 预测验证集数据
y_val_pred = model.predict(X_fold_val2)
# 预测独立的测试集和验证集数据
y_verify_test_pred = model.predict(X_verify_test2)
y_val_pred_proba = model.predict_proba(X_fold_val2)[:, 1]
y_verify_test_pred_proba = model.predict_proba(X_verify_test2)[:, 1]
# 计算验证集和独立测试集的指标
val_metrics = calculate_metrics(y_fold_val, y_val_pred, y_val_pred_proba)
verify_test_metrics = calculate_metrics(y_verify_test, y_verify_test_pred, y_verify_test_pred_proba)
# 将结果放入DataFrame
df50 = pd.DataFrame([val_metrics, verify_test_metrics],index=['XGBoost_val_feature50', 'XGBoost_verify_test_feature50'])
# 合并各类特征数目模型结果
combined_df = pd.concat([combined_df, df50], axis=0)
# 保存文件
combined_df.to_csv("AllFeature/XGBoost_10fold_crossvalidation_default_hyperparameter_ARGanno_absence01_SNPmatrix959090_absence01_"+ arguments[1] +"_fold_" + arguments[2] + "_topfeature_category_and_combine_verifytest.txt", sep='\t', index=True, header=True)

######## 寻找最优超参数 #######
xgb_model = XGBClassifier(random_state=2, enable_categorical=True)
# 使用GridSearchCV进行交叉验证和参数搜索
xgb_grid = GridSearchCV(xgb_model, param_grid=xgb_params, scoring=scoring,refit='Accuracy', cv=cv)
xgb_grid.fit(X_train2, y_train)
# 打印最佳参数和交叉验证分数
print("XGBoost Best Parameters:", xgb_grid.best_params_)
print("XGBoost Best Accuracy:", xgb_grid.best_score_)
results_XGBoost = pd.DataFrame(xgb_grid.cv_results_)
results_XGBoost.to_csv("AllFeature/XGBoost_10fold_crossvalidation_default_hyperparameter_ARGanno_absence01_SNPmatrix959090_absence01_"+ arguments[1] +"_fold_" + arguments[2] + "_top50feature_category_and_combine_verifytest.txt", sep='\t', index=False, header=True)

