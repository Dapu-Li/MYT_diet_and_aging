import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, log_loss
from lightgbm import LGBMClassifier
import scipy.stats as stats
pd.options.mode.chained_assignment = None
import os

target = 'death'
food = 'Food_resid_z_new'

# Set your data path and output path here
ImpMethod = 'TotalGain'


food_df = pd.read_csv(dpath + '/data/Food/' + food + '.csv')
target_df = pd.read_csv(dpath + '/data/death/death.csv', usecols = ['eid', 'death', 'death_duration'])
reg_df = pd.read_csv(dpath + '/data/Region_code.csv', usecols = ['eid', 'Region_code'])

mydf = pd.merge(target_df, reg_df, how = 'inner', on = ['eid'])
mydf = pd.merge(mydf, food_df, how = 'inner', on = ['eid'])

my_f_df = pd.read_csv(out_path + '/Food_Importance_new.csv')
my_f_df.sort_values(by = ImpMethod + '_cv', ascending=False, inplace = True)
food_f_lst = my_f_df.Food.tolist()
food_f_lst = [f.replace('_', ' ') for f in my_f_df.Food.tolist()]


fold_id_lst = [i for i in range(10)]


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

y_test_full = np.zeros(shape = [1,1])
for fold_id in fold_id_lst:
    test_idx = mydf['Region_code'].index[mydf['Region_code'] == fold_id]
    y_test_full = np.concatenate([y_test_full, np.expand_dims(mydf.iloc[test_idx][target], -1)])

y_pred_full_prev = y_test_full
tmp_f, loss_cv_lst= [], []

for f in food_f_lst:
    tmp_f.append(f)
    my_X = mydf[tmp_f]
    loss_cv = []
    y_pred_full = np.zeros(shape = [1,1])
    for fold_id in fold_id_lst:
        train_idx = mydf['Region_code'].index[mydf['Region_code'] != fold_id]
        test_idx = mydf['Region_code'].index[mydf['Region_code'] == fold_id]
        X_train, X_test = mydf.iloc[train_idx][tmp_f], mydf.iloc[test_idx][tmp_f]
        y_train, y_test = mydf.iloc[train_idx][target], mydf.iloc[test_idx][target]
        my_lgb = LGBMClassifier(objective='binary', metric='auc', is_unbalance=True,  n_jobs=4, verbosity=-1, seed=2023)
        my_lgb.set_params(**my_params)
        my_lgb.fit(X_train, y_train)
        y_pred_prob = my_lgb.predict_proba(X_test)[:, 1]
        loss_cv.append(np.round(log_loss(y_test, y_pred_prob), 3))
        y_pred_full = np.concatenate([y_pred_full, np.expand_dims(y_pred_prob, -1)])
    y_pred_full_prev = y_pred_full
    loss_all = log_loss(y_test_full, y_pred_full)
    tmp_out = np.array([np.round(np.mean(loss_cv), 5), np.round(np.std(loss_cv), 5), np.round(loss_all, 5)] + loss_cv)
    loss_cv_lst.append(tmp_out)
    print((f, np.round(np.mean(loss_cv), 5), np.round(loss_all, 5)))

loss_df = pd.DataFrame(loss_cv_lst, columns = ['loss_cv_mean', 'loss_cv_sd', 'loss_all'] + ['loss_' + str(i) for i in range(10)])

loss_df = pd.concat((pd.DataFrame({'Food':tmp_f}), loss_df), axis = 1)
outfile = out_path + '/AccLOSS_' + ImpMethod + '_new.csv'
loss_df.to_csv(outfile, index = False)

print('finished')
