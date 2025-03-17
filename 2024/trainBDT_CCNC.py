#!/usr/bin/env python
# coding: utf-8

# In[1]:


import ROOT


# In[2]:


# in order to start TMVA
ROOT.TMVA.Tools.Instance()

# note that it seems to be mandatory to have an
# output file, just passing None to TMVA::Factory(..)
# does not work. Make sure you don't overwrite an
# existing file.


# In[3]:


# open file
fFile1 = ROOT.TFile("/dune/data/users/lperes/refactored_ntuples/atm/NNBarAtm_hA_BR_dune10kt_1x2x6_54053565_2782_20220330T174136Z_gen_g4_detsim_reco_ntuple.root")
#fFile2 = ROOT.TFile("/dune/data/users/lperes/refactored_ntuples/atm/NNBarAtm_hA_BR_dune10kt_1x2x6_54053565_3304_20220330T174243Z_gen_g4_detsim_reco_ntuple.root")



# In[4]:


# get trees 
fTree1_s = fFile1.Get("ana/Atm")
fTree1_b = fFile1.Get("ana/Atm")
#fTree2_s = fFile1.Get("ana/Atm")
#fTree2_b = fFile1.Get("ana/Atm")

fout = ROOT.TFile("BDT_Atm.root","RECREATE")


# In[5]:


# define factory with options
factory = ROOT.TMVA.Factory("TMVAClassification", fout, #TMVA::Factory( "<JobName>", outputFile, "<options>" );
                            ":".join([    "!V", # Verbose flag, default=False
                                          "!Silent", # Batch mode: boolean silent flag inhibiting any output from TMVA after the creation of the factory class object (default: False)
                                          "Color", # Flag for coloured screen output (default: True, if in batch mode: False)
                                          "DrawProgressBar", # Draw progress bar to display training, testing and evaluation schedule (default: True)
                                          "Transformations=I;D;P;G,D", # List of transformations to test; formatting example: Transformations=I;D;P;U;G,D, for identity, decorrelation, PCA, Uniform and Gaussianisation followed by decorrelation transformations
                                          "AnalysisType=Classification"] # Set the analysis type (Classification, Regression, Multiclass, Auto) (default: Auto)
                                     ))


# In[6]:


dataloader = ROOT.TMVA.DataLoader("dataset")



# In[7]:


# add discriminating variables for training
dataloader.AddVariable("LongestTrack","F")
dataloader.AddVariable("HighestTrackSummedADC","F")
dataloader.AddVariable("PIDALongestTrack", "F")
dataloader.AddVariable("nTracks", "I")
dataloader.AddVariable("nShowers", "I")
dataloader.AddVariable("TotalMomentumP", "F")
dataloader.AddVariable("nSpacePoints", "I")
dataloader.AddVariable("NPrimaryDaughters", "I")
dataloader.AddVariable("HighestShowerSummedADC", "F")
dataloader.AddVariable("PrimaryPDGReco", "I")
dataloader.AddVariable("LargeShowerOpenAngle", "F")
dataloader.AddVariable("LongestShower", "F")
dataloader.AddVariable("CosThetaDetTotalMom", "F")
dataloader.AddVariable("CosPhiDetTotalMom", "F")
#dataloader.AddVariable("CosThetaSunRecoRange", "F")
dataloader.AddVariable("FracTotalChargeLongTrack", "F")
dataloader.AddVariable("AvarageTrackLength", "F")
dataloader.AddVariable("CVN_NCProbability", "F")

# In[8]:


# define signal and background trees
dataloader.AddSignalTree(fTree1_s)
#dataloader.AddSignalTree(fTree2_s)
dataloader.AddBackgroundTree(fTree1_b)
#dataloader.AddBackgroundTree(fTree2_b)



# In[9]:


# define additional cuts 
sigCut = ROOT.TCut("CCNC == 1 && PIDALongestTrack <60 ")
bgCut = ROOT.TCut("CCNC == 0 && PIDALongestTrack <60 ")

# In[10]:


# set options for trainings
dataloader.PrepareTrainingAndTestTree(sigCut, 
                                   bgCut, 
                                   ":".join(["nTrain_Signal=0",
                                             "nTrain_Background=0",
                                             "SplitMode=Random",
                                             "NormMode=NumEvents",
                                             "!V"
                                             ]))


# In[11]:


method = factory.BookMethod(dataloader, ROOT.TMVA.Types.kBDT, "BDT",
                            ":".join([ "!H",
                                       "!V",
                                       "NTrees=850",
                                       "nEventsMin=150",
                                       "MaxDepth=3",
                                       "BoostType=AdaBoost",
                                       "AdaBoostBeta=0.5",
                                       "SeparationType=GiniIndex",
                                       "nCuts=20",
                                       "PruneMethod=NoPruning",
                                       ]))




# In[12]:


# self-explaining
factory.TrainAllMethods()


# In[13]:


factory.TestAllMethods()


# In[14]:


factory.EvaluateAllMethods()


# In[ ]:




