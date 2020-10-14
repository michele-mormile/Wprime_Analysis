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
from tagging_jetpt_configuration import fetch_configuration, pre_cut
from array import array

def histo_jet_pt(config, config_dict, tree_name, histo_selected_signal, histo_total_signal, histo_selected_background, histo_total_background):

    dir_path = '/home/mmormile/CMSSW_9_4_9/src/Wprime/WprimeAnalysis/bin/Root/ML_scores/'

    df = root_pandas.read_root(dir_path + config + tree_name[2:] + '.root')
    df = df.query(pre_cut())
    df = df.query('not ((Top_Category == 1) and (Jet_hasPromptLepton == 1))') # toggle to avoid pathological B1 bis


    # y_pred = df['ML_Score'].values
    # y = df['Top_High_Truth'].values

    tightness = ['L', 'M', 'T']

    total_signal_events = df.query('Top_High_Truth > 0.5')
    total_background_events = df.query('Top_High_Truth < 0.5')
    threshold_increment = .2

    df_muons = df.query('(Muon_isTight == 1) & (Muon_MiniIso<0.1)')

    selected_signal_events = []
    selected_background_events = []

    selected_signal_events.append(df.query('(ML_Score > ' + str(config_dict['threshold' + tree_name[2:]]) + ') & (Top_High_Truth > 0.5)'))
    selected_background_events.append(df.query('(ML_Score > ' + str(config_dict['threshold' + tree_name[2:]]) + ') & (Top_High_Truth < 0.5)'))
    selected_signal_events.append(df.query('(ML_Score > ' + str(config_dict['threshold' + tree_name[2:]] + threshold_increment) + ') & (Top_High_Truth > 0.5)'))
    selected_background_events.append(df.query('(ML_Score > ' + str(config_dict['threshold' + tree_name[2:]] + threshold_increment) + ') & (Top_High_Truth < 0.5)'))
    for tight in tightness:
        selected_signal_events.append(df_muons.query('(Jet_isDeepCSV' + tight + ' > .5) & (Top_High_Truth > 0.5)'))
        selected_background_events.append(df_muons.query('(Jet_isDeepCSV' + tight + ' > .5) & (Top_High_Truth < 0.5)'))
    
    for index, dataframe in enumerate(selected_signal_events):
        for i in range(dataframe.shape[0]):
            histo_selected_signal[index].Fill(dataframe['Jet_Pt'].values[i])
                
    for index, dataframe in enumerate(selected_background_events):
        for i in range(dataframe.shape[0]):
            histo_selected_background[index].Fill(dataframe['Jet_Pt'].values[i])
    
    for i in range(total_signal_events.shape[0]):
        histo_total_signal.Fill(total_signal_events['Jet_Pt'].values[i])
                
    for i in range(total_background_events.shape[0]):
        histo_total_background.Fill(total_background_events['Jet_Pt'].values[i])


if __name__ == "__main__":
    
    configuration_dictionary = fetch_configuration()
    # first_bin = 0
    # last_bin = 2000
    # bins = 40

    extremes = [i for i in range(0, 1000, 100)]

    extremes.append(1200)
    extremes.append(1400)
    extremes.append(2000)
    extremes.append(2500)
    extremes.append(3000)
    
    nbins = len(extremes) - 1

    extremes = array('d', extremes)

    ROOT.gStyle.SetOptStat(000000000)

    selection_criteria = ['ML', 'ML_Tight', 'Loose', 'Medium', 'Tight']
    colors = [ROOT.kRed, ROOT.kMagenta, ROOT.kBlue, ROOT.kGreen, ROOT.kBlack]

    #Histos for signal efficiency
    histo_selected_signal_merged = [ROOT.TH1F('selected_signal_merged_' + criterion, 'selected_signal_merged_' + criterion, nbins, extremes)
                             for criterion in selection_criteria]
    histo_total_signal_merged = ROOT.TH1F('total_signal_merged', 'total_signal_merged', nbins, extremes)

    histo_selected_signal_resolved = [ROOT.TH1F('selected_signal_resolved_' + criterion, 'selected_signal_resolved_' + criterion, nbins, extremes)
                             for criterion in selection_criteria]
    histo_total_signal_resolved = ROOT.TH1F('total_signal_resolved', 'total_signal_resolved', nbins, extremes)

    #Histos for background efficiency
    histo_selected_background_merged = [ROOT.TH1F('selected_background_merged_' + criterion, 'selected_background_merged_' + criterion, nbins, extremes)
                             for criterion in selection_criteria]
    histo_total_background_merged = ROOT.TH1F('total_background_merged', 'total_background_merged', nbins, extremes)

    histo_selected_background_resolved = [ROOT.TH1F('selected_background_resolved_' + criterion, 'selected_background_resolved_' + criterion, nbins, extremes)
                             for criterion in selection_criteria]
    histo_total_background_resolved = ROOT.TH1F('total_background_resolved', 'total_background_resolved', nbins, extremes)

    for histo_list in [histo_selected_signal_merged, histo_selected_signal_resolved]:
        for histo in histo_list:
            histo.Sumw2()

    histo_total_signal_merged.Sumw2()
    histo_total_signal_resolved.Sumw2()

    for config in configuration_dictionary:
        if configuration_dictionary[config]['enable_merged']:
            histo_jet_pt(config, configuration_dictionary[config], 'is_merged', histo_selected_signal_merged, histo_total_signal_merged,
                         histo_selected_background_merged, histo_total_background_merged)
        if configuration_dictionary[config]['enable_resolved']:
            histo_jet_pt(config, configuration_dictionary[config], 'is_resolved', histo_selected_signal_resolved, histo_total_signal_resolved,
                         histo_selected_background_resolved, histo_total_background_resolved)


    histo_efficiency_signal_merged = []
    histo_efficiency_signal_resolved = []

    histo_efficiency_background_merged = []
    histo_efficiency_background_resolved = []


    for count, histo in enumerate(histo_selected_signal_merged):
        histo_efficiency_signal_merged.append(ROOT.TEfficiency(histo, histo_total_signal_merged))
        histo_efficiency_signal_merged[count].SetTitle('Efficiency_signal_merged_' + selection_criteria[count] + '; Jet pT; #epsilon')
        histo_efficiency_signal_merged[count].SetLineColor(colors[count])

    for count, histo in enumerate(histo_selected_signal_resolved):
        histo_efficiency_signal_resolved.append(ROOT.TEfficiency(histo, histo_total_signal_resolved))
        histo_efficiency_signal_resolved[count].SetTitle('Efficiency_signal_resolved_' + selection_criteria[count] + '; Jet pT; #epsilon')
        histo_efficiency_signal_resolved[count].SetLineColor(colors[count])

    for count, histo in enumerate(histo_selected_background_merged):
        histo_efficiency_background_merged.append(ROOT.TEfficiency(histo, histo_total_background_merged))
        histo_efficiency_background_merged[count].SetTitle('Efficiency_background_merged_' + selection_criteria[count] + '; Jet pT; #epsilon')
        histo_efficiency_background_merged[count].SetLineColor(colors[count])

    for count, histo in enumerate(histo_selected_background_resolved):
        histo_efficiency_background_resolved.append(ROOT.TEfficiency(histo, histo_total_background_resolved))
        histo_efficiency_background_resolved[count].SetTitle('Efficiency_background_resolved_' + selection_criteria[count] + '; Jet pT; #epsilon')
        histo_efficiency_background_resolved[count].SetLineColor(colors[count])


    c_signal_merged = ROOT.TCanvas()
    histo_efficiency_signal_merged[0].Draw('AP')
    for histo in histo_efficiency_signal_merged[1:5]:
        histo.Draw("Psame")

    legend_signal_merged = ROOT.TLegend(.85, .85, 1, 1)
    for count, histo in enumerate(histo_efficiency_signal_merged):
        legend_signal_merged.AddEntry(histo, selection_criteria[count])
    legend_signal_merged.Draw()

    c_signal_merged.SaveAs('/home/mmormile/CMSSW_9_4_9/src/Wprime/WprimeAnalysis/bin/Png/merged_jet_signal_efficiency.png')


    c_signal_resolved = ROOT.TCanvas()
    histo_efficiency_signal_resolved[0].Draw('AP')
    for histo in histo_efficiency_signal_resolved[1:5]:
        histo.Draw("Psame")

    legend_signal_resolved = ROOT.TLegend(.85, .85, 1, 1)
    for count, histo in enumerate(histo_efficiency_signal_resolved):
        legend_signal_resolved.AddEntry(histo, selection_criteria[count])
    legend_signal_resolved.Draw()

    c_signal_resolved.SaveAs('/home/mmormile/CMSSW_9_4_9/src/Wprime/WprimeAnalysis/bin/Png/resolved_jet_signal_efficiency.png')
    

    c_background_merged = ROOT.TCanvas()
    histo_efficiency_background_merged[0].Draw('AP')
    c_background_merged.Update()
    graph = histo_efficiency_background_merged[0].GetPaintedGraph()
    graph.SetMaximum(.4)
    graph.SetMinimum(0)
    c_background_merged.Update()

    for histo in histo_efficiency_background_merged[1:5]:
        histo.Draw("Psame")

    legend_background_merged = ROOT.TLegend(.85, .85, 1, 1)
    for count, histo in enumerate(histo_efficiency_background_merged):
        legend_background_merged.AddEntry(histo, selection_criteria[count])
    legend_background_merged.Draw()

    c_background_merged.SaveAs('/home/mmormile/CMSSW_9_4_9/src/Wprime/WprimeAnalysis/bin/Png/merged_jet_background_efficiency.png')


    c_background_resolved = ROOT.TCanvas()
    histo_efficiency_background_resolved[0].Draw('AP')
    c_background_resolved.Update()
    graph = histo_efficiency_background_resolved[0].GetPaintedGraph()
    graph.SetMaximum(.4)
    graph.SetMinimum(0)
    c_background_resolved.Update()

    for histo in histo_efficiency_background_resolved[1:5]:
        histo.Draw("Psame")

    legend_background_resolved = ROOT.TLegend(.85, .85, 1, 1)
    for count, histo in enumerate(histo_efficiency_background_resolved):
        legend_background_resolved.AddEntry(histo, selection_criteria[count])
    legend_background_resolved.Draw()

    c_background_resolved.SaveAs('/home/mmormile/CMSSW_9_4_9/src/Wprime/WprimeAnalysis/bin/Png/resolved_jet_background_efficiency.png')

    outFile = ROOT.TFile.Open('/home/mmormile/CMSSW_9_4_9/src/Wprime/WprimeAnalysis/bin/Root/jetpt_efficiencies.root', 'RECREATE')

    for histo_list in [histo_efficiency_signal_merged, histo_efficiency_signal_resolved,
        histo_efficiency_background_merged, histo_efficiency_background_resolved,
        histo_selected_signal_merged, histo_selected_signal_resolved,
        histo_selected_background_merged, histo_selected_background_resolved]:
        for histo in histo_list:
            histo.Write()
    
    histo_total_signal_merged.Write()
    histo_total_signal_resolved.Write()

    histo_total_background_merged.Write()
    histo_total_background_resolved.Write()

    outFile.Close()
