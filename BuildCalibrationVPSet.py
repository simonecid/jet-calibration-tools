#!/usr/bin/env python
import FWCore.ParameterSet.Config as cms
import os
import pickle as pkl
from sys import argv

def BuildCalibrationVPSet(rootfile_path, etaBinning, hfBinning, outputfile_path):

  from ROOT import TFile, TCanvas, TF1
  parameters = cms.VPSet()
  
  tfile = TFile(rootfile_path)

  ptBinning = []

  for binIdx in xrange(0, len(etaBinning) -1 ):

    parameter = cms.PSet()
    
    graph = tfile.Get("l1corr_eta_" + str(etaBinning[binIdx]) + "_" + str(etaBinning[binIdx + 1]))

    linear1 = TF1("linear1", "pol1", 80, 200)
    linear2 = TF1("linear2", "pol1", 200, 600)
    quad = TF1("quad", "pol2", 30, 80)
    fit_linear1 = graph.Fit("linear1", "Q", "", 80, 200)
    fit_linear2 = graph.Fit("linear2", "Q", "", 200, 600)
    fit_quad = graph.Fit("quad", "Q", "", 30, 80)
    
    if graph == None:
      raise ValueError("Object " + "l1corr_" + str(etaBinning[binIdx]) + "_" + str(etaBinning[binIdx + 1]) + " not found in file " + rootfile_path)
    
    l1tPtCentres = [pt for pt in graph.GetX()]
    l1tCalibrationFactorPairs = []
    for index, pt in  enumerate(l1tPtCentres):
      y = quad.Eval(30) if pt < 30 else quad.Eval(pt) if pt < 80 else linear1.Eval(pt) if pt < 200 else linear2.Eval(pt) 
      if y <= 0:
        y = 1
      if y > 2:
        print("### WARNING: CALIBRATION FACTOR > 2, ARE YOU SURE EVERYTHING IS OK?")
      graph.SetPoint(index, pt, y)
      l1tCalibrationFactorPairs.append((pt, y))

    l1tCalibrationFactorPairs = sorted(l1tCalibrationFactorPairs)
    parameter.l1tPtBins = cms.vdouble([float("-inf")] + 
      [(l1tCalibrationFactorPairs[i][0] + l1tCalibrationFactorPairs[i+1][0])/2. for i in xrange(0, len(l1tCalibrationFactorPairs)-1)] +
      [float("+inf")])
    parameter.l1tCalibrationFactors = cms.vdouble([l1tCalibrationFactorPair[1] for l1tCalibrationFactorPair in l1tCalibrationFactorPairs])
    parameter.etaMin = cms.double(etaBinning[binIdx])
    parameter.etaMax = cms.double(etaBinning[binIdx + 1])

    parameters.append(parameter)

    canvas = TCanvas()

    graph.GetYaxis().SetRangeUser(0, 3)
    # graph.GetXaxis().SetRangeUser(0, 200)
    graph.Draw("APE")
    linear2.Draw("SAME")
    linear1.Draw("SAME")
    quad.Draw("SAME")

    canvas.SaveAs("l1corr_" + str(etaBinning[binIdx]) + "_" + str(etaBinning[binIdx + 1]) + ".png")

  for binIdx in xrange(0, len(hfBinning) -1 ):
    graph = tfile.Get("l1corr_eta_" + str(hfBinning[binIdx]) + "_" + str(hfBinning[binIdx + 1]))
    l1tPtCentres = [pt for pt in graph.GetX()]

    l1tCalibrationFactorPairs = []
    for index, pt in  enumerate(l1tPtCentres):
      y = 1
      graph.SetPoint(index, pt, y)
      l1tCalibrationFactorPairs.append((pt, y))

    l1tCalibrationFactorPairs = sorted(l1tCalibrationFactorPairs)
    parameter = cms.PSet()
    parameter.etaMin = cms.double(hfBinning[binIdx])
    parameter.etaMax = cms.double(hfBinning[binIdx + 1])
    parameter.l1tPtBins = cms.vdouble([float("-inf")] + [(l1tPtCentres[i] + l1tPtCentres[i+1])/2. for i in xrange(0, len(l1tPtCentres)-1)] + [float("+inf")])
    parameter.l1tCalibrationFactors = cms.vdouble([y for y in graph.GetY()])
    parameters.append(parameter)



  tfile.Close()

  with open(outputfile_path, 'w') as f:
    f.write("import FWCore.ParameterSet.Config as cms\n")
    f.write("calibration = " + parameters.__repr__())

  return


if __name__ == "__main__":
  BuildCalibrationVPSet("ComputeUncalibratedPhase1AndAK4L1TJetsFromPfCandidates_10_0_4_MTD_5x5Jets_CalibrationFactors_MatchAK4GenJetWithAK4JetFromPfInputs.root",
  [0, 0.435, 0.783, 1.131, 1.305, 1.479, 1.653, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.964],
  [2.964, 3.489, 4.191, 5.191],
  "ComputeUncalibratedPhase1AndAK4L1TJetsFromPfCandidates_10_0_4_MTD_5x5Jets_CalibrationFactors_MatchAK4GenJetWithAK4JetFromPfInputs.py"
  )
  BuildCalibrationVPSet("ComputeUncalibratedPhase1AndAK4L1TJetsFromPfCandidates_10_0_4_MTD_5x5Jets_CalibrationFactors_MatchAK4GenJetWithPhase1L1TJetFromPfInputs.root",
  [0, 0.435, 0.783, 1.131, 1.305, 1.479, 1.653, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.964],
  [2.964, 3.489, 4.191, 5.191],
  "ComputeUncalibratedPhase1AndAK4L1TJetsFromPfCandidates_10_0_4_MTD_5x5Jets_CalibrationFactors_MatchAK4GenJetWithPhase1L1TJetFromPfInputs.py"
  )
  # BuildCalibrationVPSet("ComputeUncalibratedPhase1AndAK4L1TJetsFromPfCandidates_10_0_4_MTD_CalibrationFactors_MatchAK4GenJetWithAK4JetFromPfInputs.root",
  # [],
  # [0, 0.435, 0.783, 1.131, 1.305, 1.479, 1.653, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.964, 3.489, 4.191, 5.191],
  # "Calibration_1.py"
  # )
  # BuildCalibrationVPSet("ComputeUncalibratedPhase1AndAK4L1TJetsFromPfCandidates_10_0_4_MTD_9x9Jets_CalibrationFactors_MatchAK4GenJetWithPhase1L1TJetFromPfInputs.root",
  # [0, 0.435, 0.783, 1.131, 1.305, 1.479, 1.653, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.964],
  # [2.964, 3.489, 4.191, 5.191],
  # "ComputeUncalibratedPhase1AndAK4L1TJetsFromPfCandidates_10_0_4_MTD_9x9Jets_CalibrationFactors_MatchAK4GenJetWithPhase1L1TJetFromPfInputs.py"
  # )