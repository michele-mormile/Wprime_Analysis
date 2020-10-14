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
from tagging_jetpt_configuration import fetch_configuration, pre_cut, fetch_file_list
from math import sqrt, cos

def histo_scatter(config, tree_name, file):

    dir_path = '/home/mmormile/CMSSW_9_4_9/src/Wprime/WprimeAnalysis/bin/Root/ML_scores/'
        
    png_path = '/home/mmormile/CMSSW_9_4_9/src/Wprime/WprimeAnalysis/bin/Png/Scatter_plot/'

    root_path = '/home/mmormile/CMSSW_9_4_9/src/Wprime/WprimeAnalysis/bin/Root/Scatter_plot/'
    


    if file != 'all':
        png_path += file + '/'
        dir_path += file + '/'
        root_path += file + '/'
    outFile = ROOT.TFile.Open(root_path + 'scatter_plot' + tree_name[2:] + '.root', 'RECREATE')

    df = pd.DataFrame()
    for configuration in config:
        if config[configuration]['enable' + tree_name[2:]]:
            try:
                df = df.append(root_pandas.read_root(dir_path + configuration + tree_name[2:] + '.root'), ignore_index=True)
            except IOError:
                pass

    titles_conditions = [('total', 'Top_High_Truth > -1'), ('signal', 'Top_High_Truth > .5'), ('background', 'Top_High_Truth < .5')] 

    outFile.cd()

    muon_mass = 0.10566

    df = df.query('Muon_E*Muon_E - 0.0111640356 > Muon_Pt*Muon_Pt')

    df['Muon_P'] = df.apply(lambda row: sqrt(row['Muon_E']**2 - muon_mass**2), axis=1)
    df['Muon_Pz'] = df.apply(lambda row: sqrt(row['Muon_P']**2 - row['Muon_Pt']**2), axis=1)
    df['top_P'] = df.apply(lambda row: sqrt(row['top_E']**2 - row['top_M']**2), axis=1)

    df = df.query('top_P > top_Pt')

    df['top_Pz'] = df.apply(lambda row: sqrt(row['top_P']**2 - row['top_Pt']**2), axis=1)

    df['top_M_no_MET'] = df.apply(lambda row: sqrt(row['top_M']**2 + 2*row['top_E']*row['Muon_P'] - 
                                  2*row['top_Pt']*row['Muon_Pt']*cos(row['top_Phi']-row['MET_Phi']) - 2*row['top_Pz']*row['Muon_Pz']) , axis=1)
      
    for title_condition in titles_conditions:
        histo_create(dataframe=df, varx='Muon_Pt', vary='ML_Score', png_dir=png_path, tree_name=tree_name, xbins=50, xmin=0, xmax=3000, ybins=30, ymin=0, ymax=1,
                     title=title_condition[0], condition=title_condition[1])
        histo_create(dataframe=df, varx='top_Pt', vary='top_M', png_dir=png_path, tree_name=tree_name, xbins=50, xmin=0, xmax=3000, ybins=100, ymin=0, ymax=500,
                     title=title_condition[0], condition=title_condition[1])
        histo_create(dataframe=df, varx='Muon_Pt', vary='top_M', png_dir=png_path, tree_name=tree_name, xbins=50, xmin=0, xmax=3000, ybins=100, ymin=0, ymax=500,
                     title=title_condition[0], condition=title_condition[1])
        histo_create(dataframe=df, varx='top_Pt', vary='ML_Score', png_dir=png_path, tree_name=tree_name, xbins=50, xmin=0, xmax=3000, ybins=30, ymin=0, ymax=1,
                     title=title_condition[0], condition=title_condition[1])
        histo_create(dataframe=df, varx='top_Pt', vary='top_nu_M', png_dir=png_path, tree_name=tree_name, xbins=50, xmin=0, xmax=3000, ybins=100, ymin=0, ymax=500,
                     title=title_condition[0], condition=title_condition[1])
        histo_create(dataframe=df, varx='top_Pt', vary='top_nu_Pt', png_dir=png_path, tree_name=tree_name, xbins=50, xmin=0, xmax=3000, ybins=50, ymin=0, ymax=3000,
                     title=title_condition[0], condition=title_condition[1])
        histo_create(dataframe=df, varx='top_M', vary='top_nu_M', png_dir=png_path, tree_name=tree_name, xbins=100, xmin=0, xmax=500, ybins=100, ymin=0, ymax=500,
                     title=title_condition[0], condition=title_condition[1])
        histo_create(dataframe=df, varx='top_Pt', vary='rate_top_M_pT', png_dir=png_path, tree_name=tree_name, xbins=50, xmin=0, xmax=3000, ybins=100, ymin=0, ymax=1.5,
                     title=title_condition[0], condition=title_condition[1])
        histo_create(dataframe=df, varx='top_Pt', vary='top_M_corr', png_dir=png_path, tree_name=tree_name, xbins=500, xmin=0, xmax=3000, ybins=100, ymin=0, ymax=3000,
                     title=title_condition[0], condition=title_condition[1])
        histo_create(dataframe=df, varx='top_Pt', vary='top_M_no_MET', png_dir=png_path, tree_name=tree_name, xbins=500, xmin=0, xmax=3000, ybins=100, ymin=0, ymax=3000,
                     title=title_condition[0], condition=title_condition[1])
        histo_create(dataframe=df, varx='top_Pt', vary='top_Mt', png_dir=png_path, tree_name=tree_name, xbins=500, xmin=0, xmax=3000, ybins=100, ymin=0, ymax=3000,
                     title=title_condition[0], condition=title_condition[1])


    outFile.Close()

def histo_create(dataframe, varx, vary, png_dir, tree_name, xbins, xmin, xmax, ybins, ymin, ymax, condition, title):
    

    histo = ROOT.TH2F("histo_" + varx + '_' + vary + tree_name[2:] + '_' + title, "histo_" + varx + '_' + vary + tree_name[2:] + '_' + title,
                      xbins, xmin, xmax, ybins, ymin, ymax)

    histo.GetXaxis().SetTitle(varx)
    histo.GetYaxis().SetTitle(vary)
    x_axis = dataframe.query(condition)[varx].values
    y_axis = dataframe.query(condition)[vary].values

    for i, value in enumerate(x_axis):
        histo.Fill(value, y_axis[i])

    histo.Write()
    c = ROOT.TCanvas()
    histo.Draw('colz')
    c.SaveAs(png_dir +  varx + '_' + vary + tree_name[2:] + '_' + title +'.png')



if __name__ == "__main__":
    
    configuration_dictionary = fetch_configuration()
    file_list = fetch_file_list()
    ROOT.gStyle.SetOptStat(000000000)
    ROOT.gROOT.SetBatch(ROOT.kTRUE)

    do_all = True
    do_single = True

    if do_all:
            histo_scatter(configuration_dictionary, 'is_merged', 'all')
            histo_scatter(configuration_dictionary, 'is_resolved', 'all')
    if do_single:
        for file in file_list:
            if file_list[file]:  
                histo_scatter(configuration_dictionary, 'is_merged', file)
                histo_scatter(configuration_dictionary, 'is_resolved', file)
