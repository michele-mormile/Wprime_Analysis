def fetch_configuration():
    

    variables_jet_score = [
                               'Muon_Boosted_Pt',
                               'Muon_DB',
                               'Jet_Boosted_Pt',
                               'Muon_Boosted_Phi',
                               'Muon_Boosted_Eta',
                               'Jet_Boosted_Phi',
                               'Jet_Boosted_Eta',
                               'Muon_DBerr',
                               'Muon_Dz',
                               'Muon_DB_fract',
                               'Muon_MiniIso',
                               'Muon_Pt_Rel',
                               'Muon_Dxy',
                               'top_M',
                               'top_nu_Pt',
                               'top_nu_Phi',
                               'top_nu_Eta',
                               'top_nu_M',
                               'Costheta',
                               'Muon_IsGlobalMuon',
                               'Muon_NumberMatchedStations',
                               'Muon_NumberOfValidTrackerHits',
                               'Jet_Pt',
                               'Jet_E',
                               'Jet_Phi',
                               'Jet_DCSV'
                               ]
    variables_jet_score_less_variable = [
                                        'Muon_DB',
                                        'Jet_Boosted_Phi',
                                        'Muon_DBerr',
                                        'Muon_Dz',
                                        'Muon_DB_fract',
                                        'Muon_MiniIso',
                                        'Muon_Pt_Rel',
                                        'Muon_Dxy',
                                        'top_M',
                                        'top_nu_Pt',
                                        'top_nu_M',
                                        'Costheta',
                                        'Muon_NumberMatchedStations',
                                        'Jet_Pt',
                                        'Jet_E',
                                        'Jet_DCSV'
                               ]
    variables_jet_score_rate = [
                               'Muon_Boosted_Pt',
                               'Muon_DB',
                               'Jet_Boosted_Pt',
                               'Muon_Boosted_Phi',
                               'Muon_Boosted_Eta',
                               'Jet_Boosted_Phi',
                               'Jet_Boosted_Eta',
                               'Muon_DBerr',
                               'Muon_Dz',
                               'Muon_DB_fract',
                               'Muon_MiniIso',
                               'Muon_Pt_Rel',
                               'Muon_Dxy',
                               'rate_top_M_pT',
                               'top_nu_Pt',
                               'top_nu_Phi',
                               'top_nu_Eta',
                               'top_nu_M',
                               'Costheta',
                               'Muon_IsGlobalMuon',
                               'Muon_NumberMatchedStations',
                               'Muon_NumberOfValidTrackerHits',
                               'Jet_Pt',
                               'Jet_E',
                               'Jet_Phi',
                               'Jet_DCSV'
                               ]
    variables_jet_score_less_variable_rate = [
                                        'Muon_DB',
                                        'Jet_Boosted_Phi',
                                        'Muon_DBerr',
                                        'Muon_Dz',
                                        'Muon_DB_fract',
                                        'Muon_MiniIso',
                                        'Muon_Pt_Rel',
                                        'Muon_Dxy',
                                        'rate_top_M_pT',
                                        'top_nu_Pt',
                                        'top_nu_M',
                                        'Costheta',
                                        'Muon_NumberMatchedStations',
                                        'Jet_Pt',
                                        'Jet_E',
                                        'Jet_DCSV'
                               ]

    variables_jet_score_M_corr = [
                               'Muon_Boosted_Pt',
                               'Muon_DB',
                               'Jet_Boosted_Pt',
                               'Muon_Boosted_Phi',
                               'Muon_Boosted_Eta',
                               'Jet_Boosted_Phi',
                               'Jet_Boosted_Eta',
                               'Muon_DBerr',
                               'Muon_Dz',
                               'Muon_DB_fract',
                               'Muon_MiniIso',
                               'Muon_Pt_Rel',
                               'Muon_Dxy',
                               'top_M_corr',
                               'top_nu_Pt',
                               'top_nu_Phi',
                               'top_nu_Eta',
                               'top_nu_M',
                               'Costheta',
                               'Muon_IsGlobalMuon',
                               'Muon_NumberMatchedStations',
                               'Muon_NumberOfValidTrackerHits',
                               'Jet_Pt',
                               'Jet_E',
                               'Jet_Phi',
                               'Jet_DCSV'
                               ]
    variables_jet_score_less_variable_M_corr = [
                                        'Muon_DB',
                                        'Jet_Boosted_Phi',
                                        'Muon_DBerr',
                                        'Muon_Dz',
                                        'Muon_DB_fract',
                                        'Muon_MiniIso',
                                        'Muon_Pt_Rel',
                                        'Muon_Dxy',
                                        'top_M_corr',
                                        'top_nu_Pt',
                                        'top_nu_M',
                                        'Costheta',
                                        'Muon_NumberMatchedStations',
                                        'Jet_Pt',
                                        'Jet_E',
                                        'Jet_DCSV'
                               ]


    dic = {
        'medium_pt_jet_score': {
            'enable_merged' : True,
            'enable_resolved' : True,
            'preselection': '(Tau_High_Truth==0) & (Muon_MiniIso<.5)',
            'bin' : '(top_Pt>1000) & (top_Pt < 2000)',
            'variables': variables_jet_score,
            'threshold_merged' : .6,
            'threshold_resolved' : .4,
            'n_estimators_merged' : 130,
            'max_depth_merged' : 4,
            'min_child_weight_merged' : 5,
            'reg_alpha_merged' : .01,
            'n_estimators_resolved' : 100,
            'max_depth_resolved' : 4,
            'min_child_weight_resolved' : 4,
            'reg_alpha_resolved' : 1,
        },
        'high_pt_jet_score': {
            'enable_merged' : True,
            'enable_resolved' : True,
            'preselection': '(Tau_High_Truth==0) & (Muon_MiniIso<.5)',
            'bin' : 'top_Pt>2000',
            'variables': variables_jet_score,
            'threshold_merged' : .6,
            'threshold_resolved' : .1,
            'n_estimators_merged' : 130,
            'max_depth_merged' : 4,
            'min_child_weight_merged' : 5,
            'reg_alpha_merged' : .01,
            'n_estimators_resolved' : 100,
            'max_depth_resolved' : 4,
            'min_child_weight_resolved' : 4,
            'reg_alpha_resolved' : 1,
        },
        'low_pt_jet_score': {
            'enable_merged' : True,
            'enable_resolved' : True,
            'preselection': '(Tau_High_Truth==0) & (Muon_MiniIso<.5)',
            'bin' : 'top_Pt<1000',
            'variables': variables_jet_score,
            'threshold_merged' : .5,
            'threshold_resolved' : .2,
            'n_estimators_resolved' : 130,
            'max_depth_resolved' : 4,
            'min_child_weight_resolved' : 4,
            'reg_alpha_resolved' : .01,
            'n_estimators_merged' : 130,
            'max_depth_merged' : 4,
            'min_child_weight_merged' : 4,
            'reg_alpha_merged' : .01,
        },
    }

    return dic
def standard_Wprime_cut():
    return '(Muon_isTight == 1) & (Muon_MiniIso<0.1) & (Jet_isDeepCSVL==1)'

def pre_cut():
    return '(Muon_MiniIso<5) & (Muon_Dxy<0.1) & (Muon_Dxy>-0.1)'

def file_dir():
    return '/home/mmormile/CMSSW_9_4_9/src/Wprime/WprimeAnalysis/bin/MergedTrees/'

def fetch_file_list():

    dic = {
           'Wprime'     : True,
           'Zprime4000' : True,
           'Zprime5000' : True,
           'Zprime6000' : True,
          }

    return dic

