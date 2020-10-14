#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include <vector>
#include <assert.h>
#include <TMVA/Reader.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>
#include <sstream>
#include <string>
#include "TFileCollection.h"
#include "THashList.h"
#include "TBenchmark.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "THStack.h"

using namespace std;

void savePlots(TFile *inFile, bool mergeBool, string branchList, string globalCondition, TFile *outFile, string png_dir);

int main(int argc, char **argv)
{

  string merged_resolved(argv[1]); //"merged" for merged tree, "resolved" for resolved tree, "merged_resolved" for both

  string inputFile(argv[2]); // tree file

  string branchList(argv[3]); // list of plots to do

  string fileName(argv[4]); // name of output file

  string png_dir(argv[5]);

  TFile *inFile = TFile::Open(inputFile.c_str());

  TFile *outFile = TFile::Open(fileName.c_str(), "RECREATE");

  string ptBin[3] = {"top_Pt<1000 && (Muon_MiniIso<5) && (Muon_Dxy<0.1) && (Muon_Dxy>-0.1)", "top_Pt>1000 && top_Pt<2000 && (Muon_MiniIso<5) && (Muon_Dxy<0.1) && (Muon_Dxy>-0.1)", "top_Pt>2000 && (Muon_MiniIso<5) && (Muon_Dxy<0.1) && (Muon_Dxy>-0.1)"};

  if (merged_resolved.find("merged") == string::npos && merged_resolved.find("resolved") == string::npos)
  {
    cout << "Select resolved tree or merged tree" << endl;
    return 0;
  }

  for (int i = 0; i < 3; i++)
  {
    if (merged_resolved.find("merged") != string::npos)
      savePlots(inFile, true, branchList, ptBin[i], outFile, png_dir);

    if (merged_resolved.find("resolved") != string::npos)
      savePlots(inFile, false, branchList, ptBin[i], outFile, png_dir);
  }

  inFile->Close();
  outFile->Close();

  return 1;
}

void savePlots(TFile *inFile, bool mergeBool, string branchList, string globalCondition, TFile *outFile, string png_dir)
{

  TTree *tree;
  outFile->cd();
  TDirectory *dir;
  string merged_or_resolved;

  if (mergeBool)
  {
    merged_or_resolved = "merged";
  }

  else
  {
    merged_or_resolved = "resolved";
  }

  tree = (TTree *)inFile->Get(("is_" + merged_or_resolved).c_str());
  dir = outFile->mkdir((merged_or_resolved + "_" + globalCondition).c_str());

  ifstream in;
  in.open(branchList.c_str());

  string branch;
  string bin, start, end;

  dir->cd();

  int i = 0;

  gStyle->SetOptStat(000000000);

  while (1)
  {
    in >> branch >> bin >> start >> end;

    if (!in.good())
      break;

    tree->Draw((branch + ">>htemp" + to_string(i) + "(" + bin + "," + start + "," + end + ")").c_str(), ("Top_High_Truth < .5 && " + globalCondition).c_str());
    TH1F isto1 = *((TH1F *)gPad->GetPrimitive(("htemp" + to_string(i)).c_str()));
    isto1.SetLineColor(kBlue);

    i++;

    tree->Draw((branch + ">>htemp" + to_string(i) + "(" + bin + "," + start + "," + end + ")").c_str(), ("Top_High_Truth > .5 && " + globalCondition).c_str());
    TH1F isto2 = *((TH1F *)gPad->GetPrimitive(("htemp" + to_string(i)).c_str()));
    isto2.SetLineColor(kRed);

    i++;

    TCanvas canvas(branch.c_str());

    THStack stack;

    isto1.Scale(1/isto1.GetEntries());
    isto2.Scale(1/isto2.GetEntries());
    // isto1.DrawNormalized();

    // isto2.DrawNormalized("same");

    stack.Add(&isto1);
    stack.Add(&isto2);

    stack.Draw("histnostack");

    TLegend legend(0.8, 0.8, 1, 1);
    legend.AddEntry(&isto1, "Background");
    legend.AddEntry(&isto2, "Signal");
    legend.Draw("L");

    canvas.Write();

    string globalConditionName = globalCondition;

    if (globalCondition.find(" && ") != string::npos)
      globalConditionName.replace(globalConditionName.find(" && "), 4, "");
    if (globalCondition.find("<") != string::npos)
      globalConditionName.replace(globalConditionName.find("<"), 1, "");
    if (globalCondition.find(">") != string::npos)
      globalConditionName.replace(globalConditionName.find(">"), 1, "");

    canvas.Print(("/home/mmormile/CMSSW_9_4_9/src/Wprime/WprimeAnalysis/bin/Png/topTrueFalsePlot/" + png_dir + "/" + merged_or_resolved + "_" + branch + "_" + globalConditionName + ".png").c_str());
    isto1.Delete();
    isto2.Delete();
    stack.Delete();
    legend.Delete();

    cout << "I've done: " << branch << endl;
  }

  in.close();
}
