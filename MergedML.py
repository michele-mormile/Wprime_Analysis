#!/home/cmsuser/cmssw/slc6_amd64_gcc630/cms/cmssw/CMSSW_9_4_9/external/slc6_amd64_gcc630/bin/python
# coding: utf-8


import ROOT
import numpy as np
import pandas as pd
pd.set_option('display.max_columns', 500)
import matplotlib.pyplot as plt
plt.rcParams.update({'figure.max_open_warning': 0})
plt.rcParams['figure.figsize'] = (14, 7)
import root_pandas
from sklearn.model_selection import train_test_split
import xgboost as xgb
from sklearn.metrics import classification_report, roc_auc_score, accuracy_score, f1_score, confusion_matrix, auc, roc_curve
from ML_config import fetch_configuration, standard_Wprime_cut, pre_cut, file_dir, fetch_file_list
import os


def boosted_decision_top_tagging(config, config_dict, file_list, tree_name):

    data_dir = file_dir()
    data = pd.DataFrame()
    for file in file_list:
        file_name = data_dir + file + '_treeMerged.root'
        data = data.append(root_pandas.read_root(file_name, tree_name), ignore_index=True)
    
    data = data.query(config_dict['bin'])
    data = data.query(pre_cut())

    truth = 'Top_High_Truth'

    data_std = data.query(standard_Wprime_cut())
    
    total_events = data.shape[0]
    total_good_events = data.query('Top_High_Truth == 1').shape[0]
    selected_std_events = data_std.shape[0]
    selected_std_good_events = data_std.query('Top_High_Truth == 1').shape[0]

    if(total_good_events == 0):
        return {'configuration' : (config + tree_name[2:])}

    std_efficiency = 100*selected_std_good_events / (total_good_events * 1.0)
    std_accuracy = 100*selected_std_good_events / (selected_std_events * 1.0)
    std_bkg = 100*((selected_std_events - selected_std_good_events)*1.0) / (total_events - total_good_events)

    #Divido il data_set in train e test. Gli eventi sono quelli che userò per calcolare le efficienze allo step 0 (efficienze di ML + preselezione)
    X_train, X_test = train_test_split(data, test_size=0.3)
    
    total_events_test_step_0 = X_test.shape[0]
    true_events_test_step_0 = X_test.query('Top_High_Truth == 1').shape[0]
    
    #Prendo X_train e X_test con tutte le colonne (da usare per le top_category)
    X_train_all_variables, X_test_all_variables = X_train.query(config_dict['preselection']), X_test.query(config_dict['preselection'])
    
    #Prendo X_train e X_test da usare per il ML
    X_train, X_test = X_train_all_variables[config_dict['variables']], X_test_all_variables[config_dict['variables']]

    #Prendo train e test per il ML
    y_train, y_test = X_train_all_variables[truth].values, X_test_all_variables[truth].values

    if(np.count_nonzero(y_train) == 0 or np.count_nonzero(y_test) == 0):
        print("No data in configuration: ", config, tree_name[2:])
        return {'configuration' : (config + tree_name[2:])}
    # Hyperparameters:
    n_estimators = config_dict['n_estimators' + tree_name[2:]]
    learning_rate = 0.1
    max_depth = config_dict['max_depth' + tree_name[2:]]
    min_child_weight = config_dict['min_child_weight' + tree_name[2:]]
    reg_alpha = config_dict['reg_alpha' + tree_name[2:]]

    # Early stopping
    early_stopping_rounds = 15

    # Define model
    model_bdt = xgb.XGBClassifier(n_estimators=n_estimators, learning_rate=learning_rate,
                                  max_depth=max_depth, min_child_weight=min_child_weight,
                                  reg_alpha=reg_alpha)

    # Last in list is used for early stopping
    eval_set = [(X_train, y_train), (X_test, y_test)]

    # Fit with early stopping
    model_bdt.fit(X_train, y_train, eval_metric=["logloss"],
                  eval_set=eval_set, early_stopping_rounds=early_stopping_rounds, verbose=False)


    # make predictions for test data
    y_pred = model_bdt.predict_proba(X_test)[:, 1]
    # predictions = [round(value) for value in y_pred]

    X_test_all_variables['ML_Score'] = y_pred

    tree_path = '/home/mmormile/CMSSW_9_4_9/src/Wprime/WprimeAnalysis/bin/Root/ML_scores/'

    if len(file_list) == 1:
        tree_path += file_list[0] + '/'
        try:
            os.mkdir(tree_path)
        except OSError:
            pass

    X_test_all_variables.to_root(tree_path + config + tree_name[2:] + '.root', 'scores')
    # evaluate predictions

    fpr, tpr, thresholds = roc_curve(y_test, y_pred)

    #_std è per confrontare l'efficienza con quella dello standard cut; _knee è l'efficienza al ginocchio della ROC
    #_step0 sono le efficienze di tutto l'algoritmo di selezione, con la soglia messa uguale a quella per _knee

    th_index_std = 0

    while tpr[th_index_std] < std_efficiency/100:
        th_index_std += 1
    
    k_std = 0.
    s_std = 0.

    k_knee = 0.
    s_knee = 0.

    k_down_knee = 0.
    s_down_knee = 0.

    k_up_knee = 0.
    s_up_knee = 0.

    sel_bkg_0 = 0.
    sel_bkg_1 = 0.
    sel_bkg_2 = 0.
    sel_bkg_4 = 0.


    for count, y in enumerate(y_test):
       # print(y_pred[i])
        if y_pred[count] > thresholds[th_index_std]:
            k_std = k_std+1
            if y==1:
                s_std = s_std+1
            if X_test_all_variables['Top_Category'].values[count] == 0:
                sel_bkg_0 += 1
            if X_test_all_variables['Top_Category'].values[count] == 1:
                sel_bkg_1 += 1
            if X_test_all_variables['Top_Category'].values[count] == 2:
                sel_bkg_2 += 1
            if X_test_all_variables['Top_Category'].values[count] == 2:
                sel_bkg_4 += 1

        if y_pred[count] > config_dict['threshold' + tree_name[2:]]:
            k_knee = k_knee+1
            if y == 1:
                s_knee = s_knee+1

        if y_pred[count] > config_dict['threshold' + tree_name[2:]] -.1:
            k_down_knee = k_down_knee+1
            if y == 1:
                s_down_knee = s_down_knee+1

        if y_pred[count] > config_dict['threshold' + tree_name[2:]] +.1:
            k_up_knee = k_up_knee+1
            if y == 1:
                s_up_knee = s_up_knee+1

    b = len(y_test[y_test > .5])

    tot_bkg_0 = X_test_all_variables.query('Top_Category == 0').shape[0]
    tot_bkg_1 = X_test_all_variables.query('Top_Category == 1').shape[0]
    tot_bkg_2 = X_test_all_variables.query('Top_Category == 2').shape[0]
    tot_bkg_4 = X_test_all_variables.query('Top_Category == 4').shape[0]

    if k_std == 0 or b == 0:
        return {'configuration' : (config + tree_name[2:])}


    accuracy = 100*s_std/k_std
    efficiency = 100*s_std/b
    bkg = 100*(k_std-s_std)/(len(y_test)-b)

    efficiency_knee = 100 * s_knee/b
    bkg_knee = 100 * (k_knee - s_knee)/(len(y_test)-b)

    efficiency_down_knee = 100 * s_down_knee/b
    bkg_down_knee = 100 * (k_down_knee - s_down_knee)/(len(y_test)-b)

    efficiency_up_knee = 100 * s_up_knee/b
    bkg_up_knee = 100 * (k_up_knee - s_up_knee)/(len(y_test)-b)

    try:
        rel_bkg_0 = 100*sel_bkg_0/tot_bkg_0
    except ZeroDivisionError:
        rel_bkg_0 = -1
    try:
        abs_bkg_0 = 100*sel_bkg_0/(len(y_test) - b)
    except ZeroDivisionError:
        abs_bkg_0 = -1
    try:
        rel_bkg_1 = 100*sel_bkg_1/tot_bkg_1
    except ZeroDivisionError:
        rel_bkg_1 = -1
    try: 
        abs_bkg_1 = 100*sel_bkg_1/(len(y_test) - b)
    except ZeroDivisionError:
        abs_bkg_1 = -1
    try:
        rel_bkg_2 = 100*sel_bkg_2/tot_bkg_2
    except ZeroDivisionError:
        rel_bkg_2 = -1
    try: 
        abs_bkg_2 = 100*sel_bkg_2/(len(y_test) - b)
    except ZeroDivisionError:
        abs_bkg_2 = -1
    try:
        rel_bkg_4 = 100*sel_bkg_4/tot_bkg_4
    except ZeroDivisionError:
        rel_bkg_4 = -1
    try: 
        abs_bkg_4 = 100*sel_bkg_4/(len(y_test) - b)
    except ZeroDivisionError:
        abs_bkg_4 = -1

    eff_step_0 = 100*s_knee/true_events_test_step_0
    bkg_step_0 = 100*(k_knee-s_knee)/(total_events_test_step_0-true_events_test_step_0)

    # retrieve performance metrics
    png_dir = '/home/mmormile/CMSSW_9_4_9/src/Wprime/WprimeAnalysis/bin/Png/ML/'

    if len(file_list) == 1:
        png_dir += file_list[0] + '/'
        try:
            os.mkdir(png_dir)
        except OSError:
            pass

    png_name = png_dir + config  + tree_name[2:]

    compare_train_test(model_bdt, X_train, y_train, X_test, y_test, png_name)

    compare_train_test_background_category(clf=model_bdt,
                       X_train_all_variables=X_train_all_variables,
                       X_test_all_variables=X_test_all_variables,
                       file_name=png_name + '_background_category', variables=config_dict['variables'])

    compare_train_test_background_category(clf=model_bdt, X_train_all_variables=X_train_all_variables,
                       X_test_all_variables=X_test_all_variables, file_name=png_name+'_normalized_background',
                       variables=config_dict['variables'], normalize_background=True)

    importance_plot = xgb.plot_importance(model_bdt)
    importance_plot.figure.savefig(png_name + '_importance_plot.png')
    roc_auc = auc(fpr, tpr)
    fig_roc, ax_roc = plt.subplots(figsize=(7, 5))
    ax_roc.plot(fpr, tpr, lw=2, color='cyan', label='auc = %.3f' % (roc_auc))
    ax_roc.plot([0, 1], [0, 1], linestyle='--', lw=2, color='k', label='random chance')
    ax_roc.set_xlim([0, 1.0])
    ax_roc.set_ylim([0, 1.0])
    ax_roc.set_xlabel('false positive rate')
    ax_roc.set_ylabel('true positive rate')
    ax_roc.set_title('receiver operating curve')
    ax_roc.legend(loc="lower right")
    fig_roc.savefig(png_name + '_roc.png', format='png')

    dic_returned = {'configuration' : (config + tree_name[2:]),
                     'efficiency_standard_cut' : round(std_efficiency, 3),
                     'accuracy_standard_cut' : round(std_accuracy, 3),
                     'background_standard_cut' : round(std_bkg, 3),
                     'efficiency_machine_learning' : round(efficiency, 3),
                     'accuracy_machine_learning' : round(accuracy, 3),
                     'background_machine_learning' : round(bkg, 3),
                     'efficiency_total_selection' : round(eff_step_0, 3),
                     'background_total_selection' : round(bkg_step_0, 3),
                     'absolute_background_category_0' : round(abs_bkg_0, 3),
                     'absolute_background_category_1' : round(abs_bkg_1, 3),
                     'absolute_background_category_2' : round(abs_bkg_2, 3),
                     'absolute_background_category_4' : round(abs_bkg_4, 3),
                     'relative_background_category_0' : round(rel_bkg_0, 3),
                     'relative_background_category_1' : round(rel_bkg_1, 3),
                     'relative_background_category_2' : round(rel_bkg_2, 3),
                     'relative_background_category_4' : round(rel_bkg_4, 3),
                     'efficiency_knee' : round(efficiency_knee, 3),
                     'background_knee' : round(bkg_knee, 3),
                     'efficiency_down_knee' : round(efficiency_down_knee, 3),
                     'background_down_knee' : round(bkg_down_knee, 3),
                     'efficiency_up_knee' : round(efficiency_up_knee, 3),
                     'background_up_knee' : round(bkg_up_knee, 3),
                     'total_events' : total_events,
                     'total_true_events' : total_good_events}

    return dic_returned

def compare_train_test_background_category(clf, X_train_all_variables, X_test_all_variables, file_name, variables, normalize_background=False, bins=25):

    decisions = []
    norm = []
    for X in (X_train_all_variables, X_test_all_variables):
        d0 = clf.predict_proba(X.query('Top_Category == 0')[variables])[:, 1].ravel()
        d1 = clf.predict_proba(X.query('Top_Category == 1 & Jet_hasPromptLepton == 0')[variables])[:, 1].ravel()
        d2 = clf.predict_proba(X.query('Top_Category == 2')[variables])[:, 1].ravel()
        d3 = clf.predict_proba(X.query('Top_Category == 3')[variables])[:, 1].ravel()
        d4 = clf.predict_proba(X.query('Top_Category == 4')[variables])[:, 1].ravel()
        d5 = clf.predict_proba(X.query('Top_Category == 1 & Jet_hasPromptLepton == 1')[variables])[:, 1].ravel()
        decisions += [d0, d1, d2, d3, d4, d5]
        norm += [len(d0) + len(d1) + len(d2) + len(d3) + len(d4) + len(d5)]

    try:
        low = min(np.min(d) for d in decisions)
        high = max(np.max(d) for d in decisions)
        low_high = (low, high)
    except ValueError:
        low_high = (0,1)

    fig, ax = plt.subplots()

    weights_0 = []

    for decision in decisions[0:6]:
        weight = []
        for d in decision:
            weight.append(1./norm[0])
        weights_0.append(weight)
    
    weights_1 = []

    for decision in decisions[6:12]:
        weight = []
        for d in decision:
            weight.append(1./norm[1])
        weights_1.append(weight)
          
    colors = ['b', 'g', 'y', 'r', 'black', 'violet']

    if normalize_background:
        labels = ['B0 (train)', 'B1 (train)', 'B2 (train)', 'S (train)', 'B4 (train)', 'B1 bis (train)']

        ax.hist(decisions[0:6],
                color=colors, alpha=.6, range=low_high, bins=bins,
                label=labels, histtype='bar', normed=False, stacked=False, weights=weights_0)
    
        labels = ['B0 (test)', 'B1 (test)', 'B2 (test)', 'S (test)', 'B4 (test)', 'B1 bis (test)']
        ax.hist(decisions[6:12],
                color=colors, alpha=.4, range=low_high, bins=bins,
                label=labels, ls='dotted', histtype='bar', normed=False, stacked=False, weights=weights_1)

    else:
        labels = ['B0 (train)', 'B1 (train)', 'B2 (train)', 'S (train)', 'B4 (train)', 'B1 bis (train)']
        ax.hist(decisions[0:6],
                color=colors, alpha=.6, range=low_high, bins=bins,
                label=labels, histtype='bar', normed=not normalize_background, stacked=False)

        labels = ['B0 (test)', 'B1 (test)', 'B2 (test)', 'S (test)', 'B4 (test)', 'B1 bis (test)']
        ax.hist(decisions[6:12],
                color=colors, alpha=.4, range=low_high, bins=bins,
                label=labels, ls='dotted', histtype='bar', normed=not normalize_background, stacked=False)
    
    ax.set_xlabel("Output")
    ax.set_ylabel("Arbitrary units")
    ax.legend(loc='best')
    fig.savefig(file_name + "_test_train.png")

def compare_train_test(clf, X_train, y_train, X_test, y_test, file_name, bins=30):
    
    decisions = []
    for X,y in ((X_train, y_train), (X_test, y_test)):
        d1 = clf.predict_proba(X[y>0.5])[:,1].ravel()
        d2 = clf.predict_proba(X[y<0.5])[:,1].ravel()
        decisions += [d1, d2]
        
    low = min(np.min(d) for d in decisions)
    high = max(np.max(d) for d in decisions)
    low_high = (low, high)
    
    fig, ax = plt.subplots()

    ax.hist(decisions[0],
             color='r', alpha=0.5, range=low_high, bins=bins,
             histtype='stepfilled', normed=True,
             label='S (train)')
    ax.hist(decisions[1],
             color='b', alpha=0.5, range=low_high, bins=bins,
             histtype='stepfilled', normed=True,
             label='B (train)')

    hist, bins = np.histogram(decisions[2],
                              bins=bins, range=low_high, normed=True)
    scale = len(decisions[2]) / sum(hist)
    err = np.sqrt(hist * scale) / scale
    
    center = (bins[:-1] + bins[1:]) / 2
    ax.errorbar(center, hist, yerr=err, fmt='o', c='r', label='S (test)')
    
    hist, bins = np.histogram(decisions[3],
                              bins=bins, range=low_high, normed=True)
    scale = len(decisions[2]) / sum(hist)
    err = np.sqrt(hist * scale) / scale

    ax.errorbar(center, hist, yerr=err, fmt='o', c='b', label='B (test)')

    ax.set_xlabel("Output")
    ax.set_ylabel("Arbitrary units")
    ax.legend(loc='best')
    fig.savefig(file_name + "_test_train.png")

if __name__ == "__main__":
    config_dic = fetch_configuration()

    csv_dir = '/home/mmormile/CMSSW_9_4_9/src/Wprime/WprimeAnalysis/bin/Csv/'

    efficiencies = ['configuration',
                    'efficiency_standard_cut',
                    'accuracy_standard_cut',
                    'background_standard_cut',
                    'efficiency_machine_learning',
                    'accuracy_machine_learning',
                    'background_machine_learning',
                    'efficiency_total_selection',
                    'background_total_selection',
                    'absolute_background_category_0',
                    'absolute_background_category_1',
                    'absolute_background_category_2',
                    'absolute_background_category_4',
                    'relative_background_category_0',
                    'relative_background_category_1',
                    'relative_background_category_2',
                    'relative_background_category_4',
                    'efficiency_knee',
                    'background_knee',
                    'efficiency_down_knee',
                    'background_down_knee',
                    'efficiency_up_knee',
                    'background_up_knee',
                    'total_events',
                    'total_true_events']

    
    do_all_files = True


    if do_all_files:
        try:
            efficiencies_df = pd.read_csv(csv_dir + 'ML_efficiencies.csv')
        except IOError:
            efficiencies_df = pd.DataFrame(columns=efficiencies)

        file_list = fetch_file_list()

        for config in config_dic:
            if config_dic[config]['enable_merged']:
                efficiencies_dic = boosted_decision_top_tagging(config, config_dic[config], file_list.keys(), 'is_merged')
                efficiencies_df = efficiencies_df[efficiencies_df.configuration != (config + '_merged')]
                efficiencies_df = efficiencies_df.append(efficiencies_dic, ignore_index=True)
            if config_dic[config]['enable_resolved']:
                efficiencies_dic = boosted_decision_top_tagging(config, config_dic[config], file_list.keys(), 'is_resolved')
                efficiencies_df = efficiencies_df[efficiencies_df.configuration != (config + '_resolved')]
                efficiencies_df = efficiencies_df.append(efficiencies_dic, ignore_index=True)

        efficiencies_df.to_csv(csv_dir + 'ML_efficiencies.csv', index=False)

    do_single_files = True

    if do_single_files:

        file_list = fetch_file_list()

        for single_file in file_list:
            if file_list[single_file]:
                try:
                    efficiencies_df = pd.read_csv(csv_dir + single_file + '/ML_efficiencies.csv')
                except IOError:
                    efficiencies_df = pd.DataFrame(columns=efficiencies)


                for config in config_dic:
                    if config_dic[config]['enable_merged']:
                        efficiencies_dic = boosted_decision_top_tagging(config, config_dic[config], [single_file], 'is_merged')
                        efficiencies_df = efficiencies_df[efficiencies_df.configuration != (config + '_merged')]
                        efficiencies_df = efficiencies_df.append(efficiencies_dic, ignore_index=True)
                    if config_dic[config]['enable_resolved']:
                        efficiencies_dic = boosted_decision_top_tagging(config, config_dic[config], [single_file], 'is_resolved')
                        efficiencies_df = efficiencies_df[efficiencies_df.configuration != (config + '_resolved')]
                        efficiencies_df = efficiencies_df.append(efficiencies_dic, ignore_index=True)
                
                try:
                    os.mkdir(csv_dir + single_file + '/')
                except OSError:
                    pass

                efficiencies_df.to_csv(csv_dir + single_file + '/ML_efficiencies.csv', index=False)
