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

def histo_muon_pt_scatter(config, tree_name, histo, file_name):

    dir_path = '/home/mmormile/CMSSW_9_4_9/src/Wprime/WprimeAnalysis/bin/Root/ML_scores/' + file_name + '/'

    df = root_pandas.read_root(dir_path + config + tree_name[2:] + '.root')

    muon_pt = df['Muon_Pt'].values
    score = df['ML_Score'].values

    for i, value in enumerate(muon_pt):
        histo.Fill(value, score[i])

    # y_pred = df['ML_Score'].values
    # y = df['Top_High_Truth'].values



if __name__ == "__main__":
    
    configuration_dictionary = fetch_configuration()
    file_list = fetch_file_list()
    ROOT.gStyle.SetOptStat(000000000)

    outFile = ROOT.TFile.Open('/home/mmormile/CMSSW_9_4_9/src/Wprime/WprimeAnalysis/bin/Root/B1_distribution.root', 'RECREATE')


    png_dir = '/home/mmormile/CMSSW_9_4_9/src/Wprime/WprimeAnalysis/bin/Png/B1_distribution/'

    for single_file in file_list:
        if file_list[single_file]:
            histo_muon_pt = ROOT.TH2F("histo_muon_pt_" + single_file, "histo_muon_pt_" + single_file, 100, 0, 500, 30, 0, 1)
            histo_muon_pt.GetXaxis().SetTitle('Muon Pt')
            histo_muon_pt.GetYaxis().SetTitle('ML Score')
            
            for configuration in configuration_dictionary:
                if configuration_dictionary[configuration]['enable_merged']:
                    try:
                        histo_muon_pt_scatter(configuration, 'is_merged', histo_muon_pt, single_file)
                    except IOError:
                        pass
                if configuration_dictionary[configuration]['enable_resolved']:
                    try:
                        histo_muon_pt_scatter(configuration, 'is_resolved', histo_muon_pt, single_file)
                    except IOError:
                        pass
            
            canvas = ROOT.TCanvas()
            histo_muon_pt.Draw("colz")
            canvas.SaveAs(png_dir + single_file + '.png')
            outFile.cd()
            histo_muon_pt.Write()

    outFile.Close()
