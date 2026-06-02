from sklearn.metrics import confusion_matrix
from sklearn import metrics
from xgboost import XGBClassifier
import numpy as np
import pandas as pd 
import os
import sys
from sklearn.metrics import make_scorer, accuracy_score, roc_auc_score, f1_score, recall_score, precision_score, cohen_kappa_score, balanced_accuracy_score

arguments = sys.argv
new_working_directory = "~/Feature_analysis/"
os.chdir(new_working_directory)


traindata=pd.read_csv("use_topfeature/XGBoost_TET_topfeature_finalmodel_traindata.txt",sep='\t', header=0,index_col=0)
X_train=traindata.iloc[:, :-1]
y_train=traindata.iloc[:, -1]
y_train=y_train.replace({'R': 1, 'S': 0})

model = XGBClassifier(learning_rate=0.01,max_depth=3,n_estimators=100,random_state=2)
model.fit(X_train, y_train)

predictdata=pd.read_csv("predict_use_data/Salmonella_enterica_all_checkmPass_20231128_prokkapath_info_haveyear_havecountry_2014_arg_last_use_result_last_TET_topfeature_matrix.txt",sep='\t', header=0,index_col=0)
predictdata2 = predictdata.loc[:, X_train.columns.tolist()]
data_pred = model.predict(predictdata2)
data_pred2 = np.where(data_pred == 1, 'R', data_pred)
data_pred2 = np.where(data_pred2 == '0', 'S', data_pred2)
df = pd.DataFrame([predictdata2.index, data_pred2], index=['Genome', 'TET_predict_phenotype'])
df_transposed=df.T
df_transposed.to_csv("predict_use_data/Salmonella_enterica_all_checkmPass_20231128_prokkapath_info_haveyear_havecountry_2014_arg_last_use_result_last_TET_topfeature_predict_result.txt", sep='\t', index=False, header=True)
