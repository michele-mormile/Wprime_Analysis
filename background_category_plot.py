#!/home/cmsuser/cmssw/slc6_amd64_gcc630/cms/cmssw/CMSSW_9_4_9/external/slc6_amd64_gcc630/bin/python
# coding: utf-8

import ROOT
import pandas as pd
pd.set_option('display.max_columns', 500)
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'figure.max_open_warning': 0})
plt.rcParams['figure.figsize'] = (14, 7)
import root_pandas
from sklearn.model_selection import train_test_split
import xgboost as xgb
from ML_config import fetch_configuration, standard_Wprime_cut, pre_cut
import os

def background_category(config, config_dict, tree_name):

    data = root_pandas.read_root(config_dict['file_name'], tree_name)
    data = data.query(config_dict['bin'])
    data = data.query(pre_cut())
    data = data.query(config_dict['preselection'])

    truth = 'Top_High_Truth'

    X_train_all_variables, X_test_all_variables = train_test_split(data, test_size=0.3)

    #Prendo X_train e X_test con tutte le colonne (da usare per le top_category)
    
    #Prendo X_train e X_test da usare per il ML
    X_train, X_test = X_train_all_variables[config_dict['variables']], X_test_all_variables[config_dict['variables']]

    #Prendo train e test per il ML
    y_train, y_test = X_train_all_variables[truth].values, X_test_all_variables[truth].values

    # Hyperparameters:
    n_estimators = config_dict['n_estimators' + tree_name[2:]]
    learning_rate = 0.1
    max_depth = config_dict['max_depth' + tree_name[2:]]
    min_child_weight = config_dict['min_child_weight' + tree_name[2:]]
    reg_alpha = config_dict['reg_alpha' + tree_name[2:]]

    # Early stopping
    early_stopping_rounds = 15

    # Define model
    model_bdt = xgb.XGBClassifier(n_estimators=n_estimators, learning_rate=learning_rate, max_depth=max_depth,
                                  min_child_weight=min_child_weight, reg_alpha=reg_alpha)

    # Last in list is used for early stopping
    eval_set = [(X_train, y_train), (X_test, y_test)]

    # Fit with early stopping
    model_bdt.fit(X_train, y_train, eval_metric=["logloss"],
                  eval_set=eval_set, early_stopping_rounds=early_stopping_rounds, verbose=False)

    dir_path = '/home/mmormile/CMSSW_9_4_9/src/Wprime/WprimeAnalysis/bin/Png/Background_Category_Plots/' + config + tree_name[2:] + '/'

    try:
        os.mkdir(dir_path)
    except OSError:
        pass

    variables = ['Muon_MiniIso', 'top_M', 'Muon_Dxy', 'Muon_DB', 'Muon_Dz', 'top_nu_M']
    low_high = [(0, .05), (0, 500), (-.02, .02), (-.02, .02), (-.02, .02), (0, 500)]
    colors = ['b', 'g', 'y', 'r']
    labels = ['B0', 'B1', 'B2', 'S']
    bins = 100

    for count, variable in enumerate(variables):
        fig_total, ax_total = plt.subplots()
        for i in range(4):
            ax_total.hist(data.query('Top_Category == ' + str(i))[variable].values,
                                  color=colors[i], label=labels[i], histtype='stepfilled',
                                  alpha=.4, normed=True, range=low_high[count], bins=bins)
        ax_total.set_xlabel(variable)
        ax_total.set_ylabel("Entries")
        ax_total.legend(loc='best')
        fig_total.savefig(dir_path + variable + '_total.png')

    y_pred_train = model_bdt.predict_proba(X_train)[:, 1]
    y_pred_test = model_bdt.predict_proba(X_test)[:, 1]

    train_under_thr = X_train_all_variables[y_pred_train < config_dict['threshold' + tree_name[2:]]]
    train_over_thr = X_train_all_variables[y_pred_train > config_dict['threshold' + tree_name[2:]]]
    test_under_thr = X_test_all_variables[y_pred_test < config_dict['threshold' + tree_name[2:]]]
    test_over_thr = X_test_all_variables[y_pred_test > config_dict['threshold' + tree_name[2:]]]

    labels_over_thr = ['B0 (over threshold)', 'B1 (over thresold)', 'B2 (over threshold)', 'S (over threshold)']
    labels_under_thr = ['B0 (under threshold)', 'B1 (under thresold)', 'B2 (under threshold)', 'S (under threshold)']

    for count, variable in enumerate(variables):
        fig_train, ax_train = plt.subplots()
        fig_test, ax_test = plt.subplots()
        for i in range(4):
            ax_train.hist(train_under_thr.query('Top_Category == ' + str(i))[variable][variable < low_high[count][0] 
                                 & variable > low_high[count] [1]].values,
                                 color=colors[i], label=labels_under_thr[i], histtype='stepfilled',
                                 alpha=.5, normed=False, range=low_high[count], bins=bins)
            ax_test.hist(test_under_thr.query('Top_Category == ' + str(i))[variable][variable < low_high[count][0] 
                                 & variable > low_high[count] [1]].values,
                                 color=colors[i], label=labels_under_thr[i], histtype='stepfilled',
                                 alpha=.5, normed=False, range=low_high[count], bins=bins)
        
            hist_train, bins_train = np.histogram(train_over_thr.query('Top_Category == ' + str(i))[variable].values,
                                                  bins=bins, range=low_high[count], normed=False)
            
            # scale = len(train_over_thr.query('Top_Category == ' + str(i))[variable].values) / sum(hist_train)
            # err = np.sqrt(hist_train * scale) / scale
            err = np.sqrt(hist_train)
            center = (bins_train[:-1] + bins_train[1:]) / 2
            
            ax_train.errorbar(center, hist_train, yerr=err, fmt='o', c=colors[i], label=labels_over_thr[i])
            
            hist_test, bins_test = np.histogram(test_over_thr.query('Top_Category == ' + str(i))[variable].values,
                                     bins=bins, range=low_high[count], normed=False)
            
            # scale = len(test_over_thr.query('Top_Category == ' + str(i))[variable].values) / sum(hist_test)
            # err = np.sqrt(hist_test * scale) / scale
            err = np.sqrt(hist_test)
            center = (bins_test[:-1] + bins_test[1:]) / 2
            
            ax_test.errorbar(center, hist_test, yerr=err, fmt='o', c=colors[i], label=labels_over_thr[i])

        ax_train.set_xlabel(variable)
        ax_train.set_ylabel("Entries")
        ax_train.legend(loc='best')
        fig_train.savefig(dir_path + variable + '_train.png')

        ax_test.set_xlabel(variable)
        ax_test.set_ylabel("Entries")
        ax_test.legend(loc='best')
        fig_test.savefig(dir_path + variable + '_test.png')



if __name__ == "__main__":

    config_dic = fetch_configuration()

    for config in config_dic:
        if config_dic[config]['enable_merged']:
            background_category(config, config_dic[config], 'is_merged')
        if config_dic[config]['enable_resolved']:
            background_category(config, config_dic[config], 'is_resolved')

