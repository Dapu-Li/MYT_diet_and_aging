import numpy as np
import pandas as pd
from lightgbm import LGBMClassifier
import scipy.stats as stats
import os
import matplotlib.pyplot as plt
import seaborn as sns
pd.options.mode.chained_assignment = None

target = 'death'
food = 'Food_resid_z'

# Set your data path and output path here
ImpMethod = 'TotalGain'

food_imp_df = pd.read_csv(out_path + '/Food_Importance.csv')
food_imp_df['Food'] = food_imp_df['Food'].str.replace('_', ' ')
food_imp_df.rename(columns = {ImpMethod + '_cv': 'Food_imp'}, inplace = True)
food_auc_df = pd.read_csv(out_path + 'AccLOSS_' + ImpMethod +  '.csv')
mydf = pd.merge(food_auc_df, food_imp_df, how = 'left', on = ['Food'])

mydf['loss_lower'] = mydf['loss_cv_mean'] - mydf['loss_cv_sd']
mydf['loss_upper'] = mydf['loss_cv_mean'] + mydf['loss_cv_sd']
mydf['food_idx'] = [i for i in range(1, len(mydf)+1)]
nb_f = 9

fig, ax = plt.subplots(figsize = (12, 6))
palette = sns.color_palette("Blues",n_colors=len(mydf))
palette.reverse()
sns.barplot(ax=ax, x = "Food", y = "Food_imp", palette=palette, data=mydf.sort_values(by="Food_imp", ascending=False))
y_imp_up_lim = round(mydf['Food_imp'].max() + 0.01, 2)
ax.set_ylim([0, y_imp_up_lim])
ax.tick_params(axis='y', labelsize=18)
ax.set_xticklabels(mydf['Food'], rotation=45, fontsize=12, horizontalalignment='right')
my_col = ['r']*nb_f + ['k']*(len(mydf)-nb_f)
my_size = [14]*nb_f + [10]*(len(mydf)-nb_f)

for ticklabel, tickcolor in zip(plt.gca().get_xticklabels(), my_col):
    ticklabel.set_color(tickcolor)

ax.set_ylabel('Food importance', fontsize=24)
ax.set_xlabel('')
ax.grid(which='minor', alpha=0.2, linestyle=':')
ax.grid(which='major', alpha=0.5,  linestyle='--')
ax.set_axisbelow(True)

ax2 = ax.twinx()
ax2.plot(np.arange(nb_f-1, len(mydf)), mydf['loss_cv_mean'][nb_f-1:], 'black', alpha = 1, marker='o')
ax2.plot(np.arange(nb_f), mydf['loss_cv_mean'][:nb_f], 'red', markersize=12, alpha = 1, marker='o')
plt.fill_between(mydf['food_idx']-1, mydf['loss_lower'], mydf['loss_upper'], color = 'tomato', alpha = 0.2)
ax2.set_ylabel('Log loss', fontsize=24)
ax2.tick_params(axis='y', labelsize=18)
y_auc_up_lim = round(mydf['loss_upper'].max() + 0.001, 2)
y_auc_low_lim = round(mydf['loss_lower'].min() - 0.001, 2)
ax2.set_ylim([y_auc_low_lim, y_auc_up_lim])

fig.tight_layout()
fig.tight_layout()
plt.xlim([-.6, len(mydf)-.2])
output_img = out_path + 'SFS_' + ImpMethod + '.png'
plt.savefig(output_img, dpi = 250)
