#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

void MakeTrain()
{

    TMVA::Tools::Instance();
    TFile *input_s(0);
    TFile *input_b(0);

    TString fname_s = "/home/leoperes/GammaAxion_0.1_0.2GeV_51452090_29_20231231T010656Z_ntuple.root";
    TString fname_b = "/home/leoperes/ntuples_grid_newVertex/atm_hA_BR_4ana.root";

    if (!gSystem->AccessPathName(fname_s))
    {
        input_s = TFile::Open(fname_s); // check if file in local directory exists
    }
    if (!input_s)
    {
        std::cout << "ERROR: could not open signal data file" << std::endl;
        exit(1);
    }

    if (!gSystem->AccessPathName(fname_b))
    {
        input_b = TFile::Open(fname_b); // check if file in local directory exists
    }
    if (!input_b)
    {
        std::cout << "ERROR: could not open background data file" << std::endl;
        exit(1);
    }
    std::cout << "--- TMVAClassification       : Using input signal file: " << input_s->GetName() << std::endl;
    std::cout << "--- TMVAClassification       : Using input background file: " << input_b->GetName() << std::endl;

    // Register the training and test trees

    TTree *signalTree = (TTree *)input_s->Get("ana/Atm");
    TTree *background = (TTree *)input_b->Get("ana/Atm");

    // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
    TString outfileName("TMVA_ALP_Analysis.root");
    TFile *outputFile = TFile::Open(outfileName, "RECREATE");

    TMVA::Factory *factory = new TMVA::Factory("TMVAClassification", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");
    TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");

    dataloader->AddVariable("LongestTrack", 'F');
    // dataloader->AddVariable("HighestTrackSummedADC", 'F');
    dataloader->AddVariable("PIDALongestTrack", 'F');
    dataloader->AddVariable("nTracks", 'I');
    dataloader->AddVariable("nShowers", 'I');
    // dataloader->AddVariable("TotalMomentumP", 'F');
    dataloader->AddVariable("nSpacePoints", 'I');
    dataloader->AddVariable("NPrimaryDaughters", 'I');
    // dataloader->AddVariable("HighestShowerSummedADC", 'F');
    // dataloader->AddVariable("PrimaryPDGReco", 'I');
    dataloader->AddVariable("LargeShowerOpenAngle", 'F');
    dataloader->AddVariable("LongestShower", 'F');
    dataloader->AddVariable("NHits", 'F');
    // dataloader->AddVariable("CosThetaDetTotalMom", 'F');
    // dataloader->AddVariable("CosPhiDetTotalMom", 'F');
    // dataloader->AddVariable("CosThetaSunRecoRange", 'F');
    dataloader->AddVariable("EventRecoEnergy_Charge", 'F');
    // dataloader->AddVariable("AvarageTrackLength", 'F');
    dataloader->AddVariable("CVN_NumuProbability", 'F');
    dataloader->AddVariable("CVN_0protonsProbability", 'F');

    // global event weights per tree (see below for setting event-wise weights)
    Double_t signalWeight = 1.0;
    Double_t backgroundWeight = 1.0;
    // You can add an arbitrary number of signal or background trees
    dataloader->AddSignalTree(signalTree, signalWeight);
    dataloader->AddBackgroundTree(background, backgroundWeight);

    // Apply additional cuts on the signal and background samples (can be different)
    TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
    TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

    dataloader->PrepareTrainingAndTestTree(mycuts, mycutb,
                                           "nTrain_Signal=28000:nTrain_Background=200000:SplitMode=Random:NormMode=NumEvents:!V");
    // Adaptive Boost
    factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT",
                        "!H:!V:NTrees=400:MinNodeSize=2.5%:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3:SeparationType=GiniIndex");

    // Gradient Boost
    // factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTG",
    //                "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2");
    // Bagging
    // factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTB",
    //                  "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20");

    // Decorrelation + Adaptive Boost
    // factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTD",
    //              "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate");

    // Allow Using Fisher discriminant in node splitting for (strong) linearly correlated variables
    // factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTF",
    //            "!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20");

    // Train MVAs using the set of training events
    factory->TrainAllMethods();

    // Evaluate all MVAs using the set of test events
    factory->TestAllMethods();

    // Evaluate and compare performance of all configured MVAs
    factory->EvaluateAllMethods();

    // Save the output
    outputFile->Close();

    std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
    std::cout << "==> TMVAClassification is done!" << std::endl;

    delete factory;
    delete dataloader;
    // Launch the GUI for the root macros
    if (!gROOT->IsBatch())
        TMVA::TMVAGui(outfileName);
}