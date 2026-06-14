from collections import Counter
import numpy as np
import pandas as pd
from lightgbm import LGBMClassifier
import lightgbm
print(lightgbm.__version__)
import shap
import os
pd.options.mode.chained_assignment = None
import warnings
warnings.filterwarnings("ignore")

target = 'death'
food = 'Food_resid_z_new'

# Set your data path and output path here

food_df = pd.read_csv(dpath + '/data/Food/' + food + '.csv')
food_f_lst = food_df.columns.tolist()
food_f_lst.remove('eid')
target_df = pd.read_csv(dpath + '/data/death/death.csv', usecols = ['eid', 'death', 'death_duration'])
reg_df = pd.read_csv(dpath + '/data/Region_code.csv', usecols = ['eid', 'Region_code'])

mydf = pd.merge(target_df, reg_df, how = 'inner', on = ['eid'])
mydf = pd.merge(mydf, food_df, how = 'inner', on = ['eid'])

my_params = {
    'n_estimators': 1000,
    'learning_rate': 0.01,
    'max_depth': -1,
    'num_leaves': 31,
    'subsample': 0.8,
    'colsample_bytree': 1,
    'reg_alpha': 0.1,
    'reg_lambda': 1.0,
    'min_child_samples': 20,
    'objective': 'binary',
    'metric': 'auc'
}

def normal_imp(mydict):
    mysum = sum(mydict.values())
    mykeys = mydict.keys()
    for key in mykeys:
        mydict[key] = mydict[key]/mysum
    return mydict

tg_imp_cv = Counter()
tc_imp_cv = Counter()
shap_imp_cv = np.zeros(len(food_f_lst))
fold_id_lst = [i for i in range(10)]

for fold_id in fold_id_lst:
    print('fold ' + str(fold_id))
    train_idx = mydf['Region_code'].index[mydf['Region_code'] != fold_id]
    test_idx = mydf['Region_code'].index[mydf['Region_code'] == fold_id]
    X_train, X_test = mydf.iloc[train_idx][food_f_lst], mydf.iloc[test_idx][food_f_lst]
    y_train, y_test = mydf.iloc[train_idx][target], mydf.iloc[test_idx][target]
    my_lgb = LGBMClassifier(objective = 'binary', metric = 'auc', is_unbalance = True, verbosity = -1, seed = 2025)
    my_lgb.set_params(**my_params)
    my_lgb.fit(X_train, y_train)
    totalgain_imp = my_lgb.booster_.feature_importance(importance_type='gain')
    totalgain_imp = dict(zip(my_lgb.booster_.feature_name(), totalgain_imp.tolist()))
    totalcover_imp = my_lgb.booster_.feature_importance(importance_type='split')
    totalcover_imp = dict(zip(my_lgb.booster_.feature_name(), totalcover_imp.tolist()))
    tg_imp_cv += Counter(normal_imp(totalgain_imp))
    tc_imp_cv += Counter(normal_imp(totalcover_imp))
    explainer = shap.TreeExplainer(my_lgb)
    shap_values = explainer.shap_values(X_test)
    shap_values = np.abs(np.mean(shap_values, axis=0))
    shap_values /= shap_values.sum()
    shap_imp_cv += shap_values

shap_imp_df = pd.DataFrame({'Food': food_f_lst,
                            'ShapValues_cv': shap_imp_cv/10})
shap_imp_df.sort_values(by = 'ShapValues_cv', ascending = False, inplace = True)
shap_imp_df['Food'] = shap_imp_df['Food'].str.replace(' ', '_')

tg_imp_cv = normal_imp(tg_imp_cv)
tg_imp_df = pd.DataFrame({'Food': list(tg_imp_cv.keys()),
                          'TotalGain_cv': list(tg_imp_cv.values())})

tc_imp_cv = normal_imp(tc_imp_cv)
tc_imp_df = pd.DataFrame({'Food': list(tc_imp_cv.keys()),
                          'TotalCover_cv': list(tc_imp_cv.values())})

my_imp_df = pd.merge(left = shap_imp_df, right = tg_imp_df, how = 'left', on = ['Food'])
my_imp_df = pd.merge(left = my_imp_df, right = tc_imp_df, how = 'left', on = ['Food'])
my_imp_df['Ensemble'] = (my_imp_df['ShapValues_cv'] + my_imp_df['TotalGain_cv'] + my_imp_df['TotalCover_cv'])/3
my_imp_df.sort_values(by = 'TotalGain_cv', ascending = False, inplace = True)

outfile = out_path + 'Food_Importance_new.csv'
my_imp_df.to_csv(outfile, index = False)

print('finished')
