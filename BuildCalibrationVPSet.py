#!/usr/bin/env python
import FWCore.ParameterSet.Config as cms
import os
import pickle as pkl
from sys import argv

def BuildCalibrationVPSet(rootfile_path, etaBinning, outputfile_path):

  from ROOT import TFile
  parameters = cms.VPSet()
  
  tfile = TFile(rootfile_path)

  ptBinning = []

  for binIdx in xrange(0, len(etaBinning) -1 ):

    parameter = cms.PSet()
    
    graph = tfile.Get("l1corr_eta_" + str(etaBinning[binIdx]) + "_" + str(etaBinning[binIdx + 1]))
    
    if graph == None:
      raise ValueError("Object " + "l1corr_" + str(etaBinning[binIdx]) + "_" + str(etaBinning[binIdx + 1]) + " not found in file " + rootfile_path)
    
    l1tPtCentres = [pt for pt in graph.GetX()]

    parameter.l1tPtBins = cms.vdouble([float("-inf")] + [(l1tPtCentres[i] + l1tPtCentres[i+1])/2. for i in xrange(0, len(l1tPtCentres)-1)] + [float("+inf")])
    parameter.l1tCalibrationFactors = cms.vdouble([y for y in graph.GetY()])
    parameter.etaMin = cms.double(etaBinning[binIdx])
    parameter.etaMax = cms.double(etaBinning[binIdx + 1])
    
    parameters.append(parameter)

  tfile.Close()

  with open(outputfile_path, 'w') as f:
    f.write("import FWCore.ParameterSet.Config as cms\n")
    f.write("calibration = " + parameters.__repr__())

  return


if __name__ == "__main__":
  BuildCalibrationVPSet("/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/CalibrationFactors/Calibration_TTBar_PU200_Jets_Histograms_CalibrationTree_MatchAK4GenJetWithPhase1L1TJetFromPfCandidates_3766735.0.root",
  [0, 0.435, 0.783, 1.131, 1.305, 1.479, 1.653, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.964, 3.489, 4.191, 5.191],
  "Calibration_TTBar_PU200_Jets_Histograms_CalibrationTree_MatchAK4GenJetWithPhase1L1TJetFromPfCandidates_3766735.0.py"
  )
