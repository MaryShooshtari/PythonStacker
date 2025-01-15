import json
import sys

sorting_key = ['Lumi', 'Lumi2018', 'Lumi2017', 'Lumi2016', 'Lumi1718'
               'tttt_norm', 'tttw_norm', 'tttq_norm',
               'ttw_norm', 'ttz_norm', 'tth_norm', 'WZNorm', 'VVVNorm', 'XgammaNorm', 'Other_tNorm',
               'nonPromptElNorm', 'nonPromptElNorm2018', 'nonPromptElNorm2017', 'nonPromptElNorm2016',
               'nonPromptMuNorm', 'nonPromptMuNorm2018', 'nonPromptMuNorm2017', 'nonPromptMuNorm2016',
               'ChargeMisIDNorm',
               'WZNjets_2', 'WZNjets_3', 'WZNjets_4', 'WZNjets_5', 'WZNjets_6', 'ZZNjets_0', 'ZZNjets_1', 'ZZNjets_2', 'ZZNjets_3',
               'pileup', 'prefire',
               'ElectronReco', 'ElectronIDStat', "ElectronIDStat_16PreVFP", "ElectronIDStat_16PostVFP", 'ElectronIDSyst',
               'MuonIDSyst', 'MuonIDStat', "MuonIDStat_2016PreVFP", "MuonIDStat_2016PostVFP",
               'btag_hf', 'btag_hf_stats1', 'btag_hf_stats2',
               "btag_hf_stats1_16PreVFP", "btag_hf_stats2_16PreVFP", "btag_hf_stats1_16PostVFP", "btag_hf_stats2_16PostVFP",
               'btag_lf', 'btag_lf_stats1', 'btag_lf_stats2',
               "btag_lf_stats1_16PreVFP", "btag_lf_stats2_16PreVFP", "btag_lf_stats1_16PostVFP", "btag_lf_stats2_16PostVFP",
               'btag_cferr1', 'btag_cferr2',
               'ISR', 'ISRSignal', 'FSR', 'pdfUncertainty',
               'ScaleVarFactorization', 'ScaleVarFactorizationSignal', 'ScaleVarRenormalization', 'ScaleVarRenormalizationSignal',
               'JER_1p93', 'JER_2p5',"JER_1p93_16PreVFP", "JER_2p5_16PreVFP", "JER_1p93_16PostVFP", "JER_2p5_16PostVFP",
               'METUnclustered', "METUnclustered_16PreVFP", "METUnclustered_16PostVFP",
               'HEMIssue', 'Absolute', 'Absolute_2018', 'BBEC1', 'BBEC1_2018', 'RelativeBal', 'RelativeSample_2018',
               'FlavorQCD_light', 'FlavorQCD_charm', 'FlavorQCD_bottom',
               'ttV_addBJets', 'ttW_addJets']


if __name__ == "__main__":
    with open(sys.argv[1]) as f:
        data = json.load(f)

    sorted_data = {}
    for key in sorting_key:
        if key in data:
            sorted_data[key] = data[key]

    newfilename = sys.argv[1].replace('.json', '_sorted.json')
    with open(newfilename, 'w') as f:
        json.dump(sorted_data, f, indent=4)
