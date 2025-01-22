import ROOT 
from array import array


# Input and output file names
input_dir = "/pnfs/iihe/cms/store/user/nivanden/AnalysisOutput/ReducedTuples/2024-04-23_16-46/"
output_dir = "/pnfs/iihe/cms/store/user/mshoosht/AnalysisOutput/ReducedTuples/2024-12-13_12-13/"
input_file_names = [
"Tree_ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8_MCPrompt_2018_base.root",
"Tree_ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8_MCPrompt_2016PostVFP_base.root",
"Tree_ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8_MCPrompt_2016PreVFP_base.root",
"Tree_ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8_MCPrompt_2017_base.root",
"Tree_ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8_MCPrompt_2018_base.root",
]
for in_file in input_file_names:

    print (in_file)
    output_file_name = "Tree_ttHJetToNonbb_EFT"+in_file.split("ttHJetToNonbb")[1]
    
    # Open the input and output ROOT files
    input_file = ROOT.TFile.Open(input_dir + in_file, "READ")
    if not input_file or input_file.IsZombie():
        raise RuntimeError("Error: Could not open input file.")
    output_file = ROOT.TFile.Open(output_dir + output_file_name, "RECREATE")
    
    # Copy non-TTree objects
    for key in input_file.GetListOfKeys():
        obj = key.ReadObj()
        if not isinstance(obj, ROOT.TTree):
            #print(obj.GetName())
            obj.Write()
    
    # Get the input TTree
    input_tree = input_file.Get("TTH")
    if not input_tree:
        raise RuntimeError("Error: TTree 'TTH' not found in the input file.")
    # Create a clone of the tree without the eftVariations branch
    output_tree = input_tree.CloneTree(0)  # Clone structure, no data
    output_tree.SetName("TTH_EFT")  # Rename the tree
    
    eftVariations = ROOT.std.vector("double")()
    output_tree.SetBranchAddress("eftVariations", eftVariations)
    # Fill the new tree
    for i in range(input_tree.GetEntries()):
        input_tree.GetEntry(i)
        eftVariations.clear()
        for val in [1.0, 1.0, 0.99990205, 0.99992729, 0.9999462, 1.0, 0.99997367,0.88140415, 1.0014771, 0.99987528, 0.99982935,
                    0.99984825, 0.99990205, 0.99988799, 0.88134156, 1.0013791, 1.0001748, 0.99992867, 0.99992729,
                    0.99990096, 0.88136608, 1.0014043, 0.99991938, 0.9999462, 0.99991988, 0.88135975, 1.0014233,
                    1.0, 0.99997367, 0.88140415, 1.0014771, 0.99995334, 0.88138661, 1.0014508, 0.77031229,
                    0.88287952, 1.0060596]:
            eftVariations.push_back(val) 
        output_tree.Fill()
    
    # Write and close files
    output_tree.Write()
    output_file.Close()
    input_file.Close()
    print(f"Modified ROOT file has been saved to {output_dir + output_file_name}")

#import subprocess
#
## List of files to copy
#copy_files = [
#"Tree_ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8_MCPrompt_2016PostVFP_JEC.root",
#"Tree_ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8_MCPrompt_2016PostVFP_JER.root",
#"Tree_ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8_MCPrompt_2016PostVFP_MET.root",
#"Tree_ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8_MCPrompt_2016PreVFP_JEC.root",
#"Tree_ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8_MCPrompt_2016PreVFP_JER.root",
#"Tree_ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8_MCPrompt_2016PreVFP_MET.root",
#"Tree_ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8_MCPrompt_2017_JEC.root",
#"Tree_ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8_MCPrompt_2017_JER.root",
#"Tree_ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8_MCPrompt_2017_MET.root",
#"Tree_ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8_MCPrompt_2018_HEMIssue.root",
#"Tree_ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8_MCPrompt_2018_JEC.root",
#"Tree_ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8_MCPrompt_2018_JER.root",
#"Tree_ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8_MCPrompt_2018_MET.root",
#]
## Loop over files
#for c_file in copy_files:
#    source = input_dir + c_file
#    destination = output_dir + "Tree_ttHJetToNonbb_EFT"+c_file.split("ttHJetToNonbb")[1]
#    try:
#        subprocess.run(["xrdcp", source, destination], check=True)
#        print(f"Successfully copied: {source} -> {destination}")
#    except subprocess.CalledProcessError as e:
#        print(f"Error occurred while copying {source}: {e}")
