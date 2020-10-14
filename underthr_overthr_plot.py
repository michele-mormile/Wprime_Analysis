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
from tagging_jetpt_configuration import fetch_configuration, fetch_file_list
import os

def under_over_plot(config, config_dict, tree_name, file):

    dir_path = '/home/mmormile/CMSSW_9_4_9/src/Wprime/WprimeAnalysis/bin/Root/ML_scores/'
        
    png_path = '/home/mmormile/CMSSW_9_4_9/src/Wprime/WprimeAnalysis/bin/Png/Under_over_thr_plot/'

    root_path = '/home/mmormile/CMSSW_9_4_9/src/Wprime/WprimeAnalysis/bin/Root/Under_over_thr_plot/'

    outFile = ROOT.TFile.Open(root_path + config + tree_name[2:] + '.root', 'RECREATE')

    if file != 'all':
        png_path += file + '/'
        dir_path += file + '/'

    png_path += config + tree_name[2:] + '/'

    try:
        os.mkdir(png_path)
    except OSError:
        pass

    if not os.path.exists(dir_path + config + tree_name[2:] + '.root'):
        return

    df = root_pandas.read_root(dir_path + config + tree_name[2:] + '.root')

    outFile.cd()
    for variable in df.columns.values.tolist():
        if variable != 'ML_Score':
            histo_tot = do_histo(variable=variable, dataframe=df, condition='Top_High_Truth > -1', index=0)


            c = ROOT.TCanvas(variable + '_total')
            histo_tot.Draw()
            c.SaveAs(png_path + variable + '_total.png')
            c.Write()
            c.Close()

            stack = ROOT.THStack(variable + '_threshold', variable + '_threshold')

            histo_signal_over = do_histo(variable=variable, dataframe=df, condition='(Top_High_Truth > 0.5) & (ML_Score > ' + str(config_dict['threshold' + tree_name[2:]]) + ')', index=1)
            histo_signal_under = do_histo(variable=variable, dataframe=df, condition='(Top_High_Truth > 0.5) & (ML_Score < ' + str(config_dict['threshold' + tree_name[2:]]) + ')', index=2)
            histo_bkg_over = do_histo(variable=variable, dataframe=df, condition='(Top_High_Truth < 0.5) & (ML_Score > ' + str(config_dict['threshold' + tree_name[2:]]) + ')', index=3)
            histo_bkg_under = do_histo(variable=variable, dataframe=df, condition='(Top_High_Truth < 0.5) & (ML_Score < ' + str(config_dict['threshold' + tree_name[2:]]) + ')', index=4)

            try:
                histo_signal_over.Scale(1/histo_signal_over.GetEntries())
            except ZeroDivisionError:
                pass
            try:
                histo_signal_under.Scale(1/histo_signal_under.GetEntries())
            except ZeroDivisionError:
                pass
            try:
                histo_bkg_over.Scale(1/histo_bkg_over.GetEntries())
            except ZeroDivisionError:
                pass
            try:
                histo_bkg_under.Scale(1/histo_bkg_under.GetEntries())
            except ZeroDivisionError:
                pass

            histo_signal_over.SetLineColor(ROOT.kRed)
            histo_signal_under.SetLineColor(ROOT.kBlack)
            histo_bkg_over.SetLineColor(ROOT.kBlue)
            histo_bkg_under.SetLineColor(ROOT.kGreen)

            legend = ROOT.TLegend()
            legend.AddEntry(histo_signal_over, 'Signal over Threshold')
            legend.AddEntry(histo_signal_under, 'Signal under Threshold')
            legend.AddEntry(histo_bkg_over, 'Bkg over Threshold')
            legend.AddEntry(histo_bkg_under, 'Bkg under Threshold')

            d = ROOT.TCanvas(variable + '_thr')

            stack.Add(histo_signal_over)
            stack.Add(histo_signal_under)
            stack.Add(histo_bkg_over)
            stack.Add(histo_bkg_under)

            stack.Draw('nostackhist')

            # legend.Draw()

            d.SaveAs(png_path + variable + '_thresholds.png')
            d.Write()
            d.Close()
        if variable == 'ML_Score':
            for i in range(3):

                stack = ROOT.THStack()

                df_signal = df.query('Top_High_Truth == 1 & Neutrino_number == ' + str(i))
                df_background = df.query('Top_High_Truth == 0 & Neutrino_number == ' + str(i))

                histo_signal = ROOT.TH1F("ML_Score_Signal_Neutrino_" + str(i), "ML_Score_Signal_Neutrino_" + str(i), 30, 0, 1)
                histo_background = ROOT.TH1F("ML_Score_Background_Neutrino_" + str(i), "ML_Score_Background_Neutrino_" + str(i), 30, 0, 1)

                for value in df_signal['ML_Score'].values:
                    histo_signal.Fill(value)
                for value in df_background['ML_Score'].values:
                    histo_background.Fill(value)

                histo_signal.SetLineColor(ROOT.kRed)
                histo_background.SetLineColor(ROOT.kBlue)

                stack.Add(histo_signal)
                stack.Add(histo_background)

                legend = ROOT.TLegend(0.4,0.8,0.6,1)
                legend.AddEntry(histo_signal, "Signal")
                legend.AddEntry(histo_background, "Background")

                c = ROOT.TCanvas()
                stack.Draw('nostackhist')
                legend.Draw()
                c.SaveAs(png_path + "ML_Score_Neutrino_" + str(i) + '.png')
                c.Write()
                c.Close()

    
    outFile.Close()
            

def do_histo(variable, dataframe, condition, index): #index is necessary to avoid memory leaks

    median = dataframe[variable].median()
    std_dev = dataframe[variable].std()

    xmin = median - 3*std_dev
    xmax = median + 3*std_dev

    histo = ROOT.TH1F(variable + '_' + str(index), variable + '_' + str(index), 100, xmin, xmax)

    for value in dataframe.query(condition)[variable].values:
        histo.Fill(value)

    return histo

if __name__ == "__main__":

    ROOT.gROOT.SetBatch(ROOT.kTRUE)

    configuration_dictionary = fetch_configuration()
    file_list = fetch_file_list()
    do_all = True
    do_single = True

    for configuration in configuration_dictionary:
        if do_all:
            if configuration_dictionary[configuration]['enable_merged']:
                    under_over_plot(configuration, configuration_dictionary[configuration], 'is_merged', 'all')
            if configuration_dictionary[configuration]['enable_resolved']:
                    under_over_plot(configuration, configuration_dictionary[configuration], 'is_resolved', 'all')
        if do_single:
            for file in file_list:
                if configuration_dictionary[configuration]['enable_merged']:
                        under_over_plot(configuration, configuration_dictionary[configuration], 'is_merged', file)
                if configuration_dictionary[configuration]['enable_resolved']:
                        under_over_plot(configuration, configuration_dictionary[configuration], 'is_resolved', file)
