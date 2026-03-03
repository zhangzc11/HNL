import ROOT as rt

import sys
import os


def GetKeyNames( self, dir = "" ):
        self.cd(dir)
        return [key.GetName() for key in rt.gDirectory.GetListOfKeys()]
def GetClassNames( self, dir = "" ):
        self.cd(dir)
        return [key.GetClassName() for key in rt.gDirectory.GetListOfKeys()]

rt.TFile.GetKeyNames = GetKeyNames
rt.TFile.GetClassNames = GetClassNames

if __name__ == '__main__':
    inputFile = "/pool1/lhcb/zhicaiz/mc/2024/HNL_Wmumuqq/ntuple/ntuple_20260228/42912010_5GeVCtau0ps_options_20260228_111712.root"
    if len(sys.argv) > 1:
        inputFile = sys.argv[1]
    outputFile = inputFile.replace(".root","_skim.root")
    if len(sys.argv) > 2:
        outputFile = sys.argv[2]

    isMC = False

    if "/mc/" in inputFile:
        isMC = True

    if len(sys.argv) > 3:
        if sys.argv[3] == "yes" or sys.argv[3] == "True" or sys.argv[3] == "true":
            isMC = True

    do_mc_truth_matching = False
    if isMC: 
        do_mc_truth_matching = True


    os.system("rm "+outputFile)

    cut_all = ""

    #cut_reco_common = "Mu_PT > 15000.0 && SpruceQEE_SingleHighPtMuonDecision && Hlt2QEE_SingleHighPtMuonFullDecision && Hlt1SingleHighPtMuonDecision && Mu_INECAL && Mu_INHCAL && Mu_NUTHITS > 0 && Mu_TRCHI2DOF > 0.01 && Mu_ETA > 2.0 && Mu_ETA < 4.5 && MuNuR_PT > 3000.0 && MuNuR_ETA > 2.0 && MuNuR_ETA < 4.5 && MuNuR_INECAL && MuNuR_INHCAL && Jet1_CPF > 0.1 && Jet1_MPT > 1200.0 && NuR_M < 80000.0 && W_M > 60000.0 && W_M < 100000.0 && W_MmuWmuN > 20000.0 && W_MmuWmuN < 70000.0"
    cut_reco_common = "Mu_PT > 15000.0 && SpruceQEE_SingleHighPtMuonDecision && Hlt2QEE_SingleHighPtMuonFullDecision && Hlt1SingleHighPtMuonDecision && Mu_INECAL && Mu_INHCAL && Mu_NUTHITS > 0 && Mu_TRCHI2DOF > 0.01 && Mu_ETA > 2.0 && Mu_ETA < 4.5 && MuNuR_PT > 2000.0 && MuNuR_ETA > 2.0 && MuNuR_ETA < 4.5 && MuNuR_INECAL && MuNuR_INHCAL && Jet1_CPF > 0.1 && Jet1_MPT > 1200.0 && NuR_M < 80000.0 && W_M > 20000.0 && W_M < 100000.0 && W_MmuWmuN > 10000.0 && W_MmuWmuN < 70000.0"

    cut_dict = {
        "myTupleSS1J": "Jet1_PT > 10000.0", 
        "myTupleOS1J": "Jet1_PT > 10000.0", 
        "myTupleSS2J": "Jet1_PT > 10000.0 && Jet2_PT > 10000.0 && Jet2_CPF > 0.1 && Jet2_MPT > 1200.0", 
        "myTupleOS2J": "Jet1_PT > 10000.0 && Jet2_PT > 10000.0 && Jet2_CPF > 0.1 && Jet2_MPT > 1200.0", 
        }
    for key in cut_dict:
        cut_dict[key] = cut_dict[key]+" && "+cut_reco_common

    cut_truth_matching_dict = {
        "myTupleSS1J": "abs(Mu_TRUEID)==13 && abs(MuNuR_TRUEID)==13 && abs(Mu_MC_MOTHER_ID) == 9900024 && abs(MuNuR_MC_MOTHER_ID) == 9900014 && abs(MuNuR_MC_GD_MOTHER_ID) == 9900024 && (Jet1_ConMC_GDall_MOTHER_ID_N9900014>0 || Jet1_ConTRUEID_N9900014 > 0)",
        "myTupleOS1J": "abs(Mu_TRUEID)==13 && abs(MuNuR_TRUEID)==13 && abs(Mu_MC_MOTHER_ID) == 9900024 && abs(MuNuR_MC_MOTHER_ID) == 9900014 && abs(MuNuR_MC_GD_MOTHER_ID) == 9900024 && (Jet1_ConMC_GDall_MOTHER_ID_N9900014>0 || Jet1_ConTRUEID_N9900014 > 0)",
        "myTupleSS2J": "abs(Mu_TRUEID)==13 && abs(MuNuR_TRUEID)==13 && abs(Mu_MC_MOTHER_ID) == 9900024 && abs(MuNuR_MC_MOTHER_ID) == 9900014 && abs(MuNuR_MC_GD_MOTHER_ID) == 9900024 && (Jet1_ConMC_GDall_MOTHER_ID_N9900014>0 || Jet1_ConTRUEID_N9900014 > 0)  &&  (Jet2_ConMC_GDall_MOTHER_ID_N9900014>0 || Jet2_ConTRUEID_N9900014 > 0)",
        "myTupleOS2J": "abs(Mu_TRUEID)==13 && abs(MuNuR_TRUEID)==13 && abs(Mu_MC_MOTHER_ID) == 9900024 && abs(MuNuR_MC_MOTHER_ID) == 9900014 && abs(MuNuR_MC_GD_MOTHER_ID) == 9900024 && (Jet1_ConMC_GDall_MOTHER_ID_N9900014>0 || Jet1_ConTRUEID_N9900014 > 0)  &&  (Jet2_ConMC_GDall_MOTHER_ID_N9900014>0 || Jet2_ConTRUEID_N9900014 > 0)",
    }

    if do_mc_truth_matching:
        for key in cut_dict:
            if cut_truth_matching_dict.get(key) is not None:
                cut_dict[key] = cut_dict[key]+" && "+cut_truth_matching_dict[key]

    print("input: "+inputFile)
    print("output: "+outputFile)
    print("cut_all: "+cut_all)
    print("cut_dict: ")
    print(cut_dict)

    inFile = rt.TFile(inputFile, "READ")

    keyList = inFile.GetKeyNames()
    classList = inFile.GetClassNames()

    print(keyList)
    print(classList)

    outFile = rt.TFile(outputFile, "RECREATE")
    outFile.cd()

    for j in range(0, len(keyList)):
        print(classList[j] + "   ===   " + keyList[j])
        if classList[j] == "TTree":
            print("===skimming tree: "+keyList[j])
            inFile.cd()
            inputTree = inFile.Get(keyList[j])
            print("events before cut: "+str(inputTree.GetEntries()))
            outFile.cd()
            outputTree = inputTree.CopyTree(cut_all)
            outputTree.Write()
            print("events after cut: "+str(outputTree.GetEntries()))
        if classList[j] == "TH1F" or classList[j] == "TH1D":
            inFile.cd()
            histThis = inFile.Get(keyList[j])
            outFile.cd()
            histThis_out = histThis.Clone()
            histThis_out.Write()
        if classList[j] == "TDirectoryFile":
            # get all trees in the directory
            print("========skimming directory: "+keyList[j])
            subkeyList = inFile.GetKeyNames(keyList[j])
            subclassList = inFile.GetClassNames(keyList[j])
            outFile.mkdir(keyList[j])
            cut_this = cut_all
            cut_dict_this = cut_dict.get(keyList[j])
            if cut_dict_this is not None:
                if cut_this == "":
                    cut_this = cut_dict_this
                else:
                    cut_this = cut_this + " && " + cut_dict_this
            print("cut for trees in this directory: "+cut_this)
            for subj in range(0, len(subkeyList)):
                if subclassList[subj] == "TTree":
                    print("===skimming tree: "+keyList[j]+"/"+subkeyList[subj])
                    inFile.cd(keyList[j])
                    inputTree = inFile.Get(keyList[j]+"/"+subkeyList[subj])
                    print("events before cut: "+str(inputTree.GetEntries()))
                    outFile.cd(keyList[j])
                    outputTree = inputTree.CopyTree(cut_this)
                    outputTree.Write()
                    print("events after cut: "+str(outputTree.GetEntries()))
                if subclassList[subj] == "TH1F" or subclassList[subj] == "TH1D":
                    inFile.cd(keyList[j])
                    inputHist = inFile.Get(keyList[j]+"/"+subkeyList[subj])
                    outFile.cd(keyList[j])
                    outputHist = inputHist.Clone()
                    outputHist.Write()
    outFile.Close()
