#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TTreeReader.h"
#include <cmath>

float deltaPhi(float p1, float p2)
{
  //Computes delta phi, handling periodic limit conditions
  float res = p1 - p2;
  while (res > M_PI) res -= 2 * M_PI;
  while (res < -M_PI) res += 2 * M_PI;
  return res;
}

void ConvertTree(TChain & inputChain, TTree & outputTree)
{
  float genJet_pt_branch;
  float genJet_eta_branch;
  float genJet_phi_branch;
  float l1tJet_pt_branch;
  float l1tJet_eta_branch;
  float l1tJet_phi_branch;
  float deltaR2_branch;

  inputChain.SetBranchAddress("genJet_pt", &genJet_pt_branch);
  inputChain.SetBranchAddress("genJet_eta", &genJet_eta_branch);
  inputChain.SetBranchAddress("genJet_phi", &genJet_phi_branch);
  inputChain.SetBranchAddress("caloJet_pt", &l1tJet_pt_branch);
  inputChain.SetBranchAddress("caloJet_eta", &l1tJet_eta_branch);
  inputChain.SetBranchAddress("caloJet_phi", &l1tJet_phi_branch);
  inputChain.SetBranchAddress("deltaR2", &deltaR2_branch);

  // Quantities for L1 jets:
  float l1tJet_pt(-1.), l1tJet_eta(99.), l1tJet_phi(99.);
  outputTree.Branch("pt",     &l1tJet_pt,     "pt/F");
  outputTree.Branch("eta",    &l1tJet_eta,    "eta/F");
  outputTree.Branch("phi",    &l1tJet_phi,    "phi/F");
  // Quantities for reference jets (GenJet, etc):
  float genJet_pt(-1.), genJet_eta(99.), genJet_phi(99.);
  outputTree.Branch("ptRef",  &genJet_pt, "ptRef/F");
  outputTree.Branch("etaRef", &genJet_eta, "etaRef/F");
  outputTree.Branch("phiRef", &genJet_phi, "phiRef/F");

  // Quantities to describe relationship between the two:
  float out_ratio(-1.), out_ratio_inv(-1.);
  float out_dr(99.), out_deta(99.), out_dphi(99.);
  float out_ptDiff(99999.), out_resL1(99.), out_resRef(99.);

  outputTree.Branch("ptDiff", &out_ptDiff, "ptDiff/F"); // L1 - Ref
  outputTree.Branch("rsp",    &out_ratio,    "rsp/F"); // response = l1 pT/ ref jet pT
  outputTree.Branch("rsp_inv",   &out_ratio_inv,   "rsp_inv/F"); // response = ref pT/ l1 jet pT
  outputTree.Branch("dr",     &out_dr,     "dr/F");
  outputTree.Branch("deta",   &out_deta,   "deta/F"); // delta eta
  outputTree.Branch("dphi",   &out_dphi,   "dphi/F"); // delta phi
  outputTree.Branch("resL1", &out_resL1, "resL1/F"); // resolution = L1 - Ref / L1
  outputTree.Branch("resRef", &out_resRef, "resRef/F"); // resolution = L1 - Ref / Ref
  // PU quantities
  //float out_trueNumInteractions(-1.), out_numPUVertices(-1.);
  //outputTree.Branch("trueNumInteractions", &out_trueNumInteractions, "trueNumInteractions/F");
  //outputTree.Branch("numPUVertices", &out_numPUVertices, "numPUVertices/F");
  Int_t nevents = inputChain.GetEntries();
  for (int x = 0; x < nevents; x++)
  {
    inputChain.GetEvent(x);
    if (x % 1000 == 0)
      std::cout << "\r" << x << "/" << nevents;
    
    genJet_pt = genJet_pt_branch;
    genJet_eta = genJet_eta_branch;
    genJet_phi = genJet_phi_branch;
    l1tJet_pt = l1tJet_pt_branch;
    l1tJet_eta = l1tJet_eta_branch;
    l1tJet_phi = l1tJet_phi_branch;
    out_dr = sqrt(deltaR2_branch);

    out_ptDiff = l1tJet_pt - genJet_pt;
    out_ratio = l1tJet_pt / genJet_pt;
    out_ratio_inv = genJet_pt / l1tJet_pt;
    out_deta = genJet_eta - l1tJet_eta;
    out_dphi = deltaPhi(genJet_phi, l1tJet_phi);
    out_resL1 = out_ptDiff/l1tJet_pt;
    out_resRef = out_ptDiff/genJet_pt;
    
    outputTree.Fill();
  }
  std::cout << "\r" << nevents << "/" << nevents;
  
  std::cout << std::endl;
  return;

}

int main(int argc, char** argv)
{
  if (argc < 3) 
  {
    std::cout << "Not enough parameters have been passed." << std::endl;
    std::cout << "Usage: " << argv[0] << " [input_file.root] [output_file.root]" << std::endl;
  }


  const char * inputFile = argv[1];
  TFile* outputFile = new TFile(argv[2], "CREATE");

  if (!outputFile -> IsOpen()) 
  {
    std::cerr << "Error while opening the file." << std::endl;
    delete inputFile;
    delete outputFile;
    return -1;
  }
  
  TChain* MatchAK4GenJetWithPhase1L1TJetFromPfClusters = new TChain("MatchAK4GenJetWithPhase1L1TJetFromPfClusters/matchedCaloJetGenJetTree");
  TChain* MatchAK4GenJetWithPhase1L1TJetFromPfCandidates = new TChain("MatchAK4GenJetWithPhase1L1TJetFromPfCandidates/matchedCaloJetGenJetTree");
  TChain* MatchAK4GenJetWithAK4JetFromPfClusters = new TChain("MatchAK4GenJetWithAK4JetFromPfClusters/matchedCaloJetGenJetTree");
  TChain* MatchAK4GenJetWithAK4JetFromPfCandidates = new TChain("MatchAK4GenJetWithAK4JetFromPfCandidates/matchedCaloJetGenJetTree");

  MatchAK4GenJetWithPhase1L1TJetFromPfClusters -> Add(inputFile);
  MatchAK4GenJetWithPhase1L1TJetFromPfCandidates -> Add(inputFile);
  MatchAK4GenJetWithAK4JetFromPfClusters -> Add(inputFile);
  MatchAK4GenJetWithAK4JetFromPfCandidates -> Add(inputFile);

  TTree* MatchAK4GenJetWithPhase1L1TJetFromPfClusters_calibration = new TTree("MatchAK4GenJetWithPhase1L1TJetFromPfClusters", "Calibration ttree with Phase1L1TJet - PfClusters");
  TTree* MatchAK4GenJetWithPhase1L1TJetFromPfCandidates_calibration = new TTree("MatchAK4GenJetWithPhase1L1TJetFromPfCandidates", "Calibration ttree with Phase1L1TJet - PfCandidates");
  TTree* MatchAK4GenJetWithAK4JetFromPfClusters_calibration = new TTree("MatchAK4GenJetWithAK4JetFromPfClusters", "Calibration ttree with AK4Jet - PfClusters");
  TTree* MatchAK4GenJetWithAK4JetFromPfCandidates_calibration = new TTree("MatchAK4GenJetWithAK4JetFromPfCandidates", "Calibration ttree with AK4Jet - PfCandidates");

  ConvertTree(*MatchAK4GenJetWithPhase1L1TJetFromPfClusters, *MatchAK4GenJetWithPhase1L1TJetFromPfClusters_calibration);
  ConvertTree(*MatchAK4GenJetWithPhase1L1TJetFromPfCandidates, *MatchAK4GenJetWithPhase1L1TJetFromPfCandidates_calibration);
  ConvertTree(*MatchAK4GenJetWithAK4JetFromPfClusters, *MatchAK4GenJetWithAK4JetFromPfClusters_calibration);
  ConvertTree(*MatchAK4GenJetWithAK4JetFromPfCandidates, *MatchAK4GenJetWithAK4JetFromPfCandidates_calibration);

  outputFile -> cd();
  MatchAK4GenJetWithPhase1L1TJetFromPfClusters_calibration -> Write();
  MatchAK4GenJetWithPhase1L1TJetFromPfCandidates_calibration -> Write();
  MatchAK4GenJetWithAK4JetFromPfClusters_calibration -> Write();
  MatchAK4GenJetWithAK4JetFromPfCandidates_calibration -> Write();

  outputFile -> Close();
  delete outputFile;

  return 0;

}
