#include "inputParams.h"

void plot(bool doPrompt = true, bool doPbPb = true){
  if (!setSystTag(doPbPb)) return;

  gStyle->SetOptStat(0);

  string filename = "";
  string outputfile = "";
  filename = Form("%s/dataUnf/unfOutput/step%i/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin%s.root",unfPath.c_str(), 1,doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", doPbPb?nIter:nIter_pp, nBinZ_reco, nBinJet_reco, nBinZ_reco, nBinJet_reco, systTag.c_str());

  outputfile = Form("%s/dataUnf/unfOutput/finalResults/trMatrix_%s_%s%s.png",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt",noSmearing?"_noSmearing":"");


  cout <<"filename = "<<filename <<endl;    
  TFile *file = new TFile(filename.c_str());
  // get inverse transfer matrix
  TMatrixT<double> *invMatrix  = (TMatrixT<double> *)file->Get(Form("invtrmat%d",doPbPb?nIter:nIter_pp));
  
  TH2D *h2_UnfData  = (TH2D*)file->Get(Form("hReco_Iter%d",doPbPb?nIter:nIter_pp));
  TH2D *h2_MeasData  = (TH2D*)file->Get("fh2MeasData");
  
  TH2D *h_invMatrix = new TH2D(*invMatrix);
  //invMatrix->Print();
  //TH2D *h_covMatrix = new TH2D(*covMatrix);
  //h_invMatrix->Rebin2D(1,nBinJet_gen*nBinZ_gen/(nBinJet_reco*nBinZ_reco));
  //h_invMatrix->GetZaxis()->SetRangeUser(0,1);

  TLine *line0 = new TLine(4,0,4,12);
  line0->SetLineColor(kBlack);
  line0->SetLineStyle(1);
  line0->SetLineWidth(2);

  TLine *line1 = new TLine(8,0,8,12);
  line1->SetLineColor(kBlack);
  line1->SetLineStyle(1);
  line1->SetLineWidth(2);

  TLine *line2 = new TLine(12,0,12,12);
  line2->SetLineColor(kBlack);
  line2->SetLineStyle(1);
  line2->SetLineWidth(2);

  TLine *line5 = new TLine(0,4,12,4);
  line5->SetLineColor(kBlack);
  line5->SetLineStyle(1);
  line5->SetLineWidth(2);

  TLine *line6 = new TLine(0,8,12,8);
  line6->SetLineColor(kBlack);
  line6->SetLineStyle(1);
  line6->SetLineWidth(2);

  TLine *line7 = new TLine(0,12,12,12);
  line7->SetLineColor(kBlack);
  line7->SetLineStyle(1);
  line7->SetLineWidth(2);

  TCanvas *can = new TCanvas("can","can",1000,900);

  gPad->SetLeftMargin(0.1);
  gPad->SetRightMargin(0.25);

  gStyle->SetPalette(56);

  
  h_invMatrix->GetXaxis()->SetTickLength(0.0);
  h_invMatrix->GetYaxis()->SetTickLength(0.0);
  

  h_invMatrix->GetXaxis()->SetNdivisions(3, 4);
  h_invMatrix->GetYaxis()->SetNdivisions(3, 4);
  
  h_invMatrix->GetXaxis()->SetTitleSize(0);
  h_invMatrix->GetXaxis()->SetLabelSize(0);
  h_invMatrix->GetYaxis()->SetTitleSize(0);
  h_invMatrix->GetYaxis()->SetLabelSize(0);

  h_invMatrix->Draw("colz");

  line0->Draw("same");
  line1->Draw("same");
  line2->Draw("same");
  line5->Draw("same");
  line6->Draw("same");
  line7->Draw("same");

  TLatex *texf = new TLatex(-1,2.,"z^{2} unfolde");
  texf->SetTextFont(42);
  texf->SetTextColor(kBlack);
  texf->SetTextAlign(22);
  texf->SetTextSize(0.035);
  texf->SetTextAngle(90);
  texf->Draw();

  TLatex *texf0 = new TLatex(2.,-1.0,"z^{2} mesure");
  texf0->SetTextFont(42);
  texf0->SetTextColor(kBlack);
  texf0->SetTextAlign(22);
  texf0->SetTextSize(0.035);
  texf0->Draw();

  TLatex *texf1 = new TLatex(2.,12.4,"10 #minus 20");
  texf1->SetTextFont(42);
  texf1->SetTextColor(kBlack);
  texf1->SetTextAlign(22);
  texf1->SetTextSize(0.03);
  texf1->Draw();

  TLatex *texf2 = new TLatex(6.,12.4,"20 #minus 30");
  texf2->SetTextFont(42);
  texf2->SetTextColor(kBlack);
  texf2->SetTextAlign(22);
  texf2->SetTextSize(0.03);
  texf2->Draw();

  TLatex *texf3 = new TLatex(10.,12.4,"30 #minus 40");
  texf3->SetTextFont(42);
  texf3->SetTextColor(kBlack);
  texf3->SetTextAlign(22);
  texf3->SetTextSize(0.03);
  texf3->Draw();

  TLatex *texf5 = new TLatex(6.,13.,"Jet p_{T} mesure [GeV]");
  texf5->SetTextFont(42);
  texf5->SetTextColor(kBlack);
  texf5->SetTextAlign(22);
  texf5->SetTextSize(0.04);
  texf5->Draw();
    
  TLatex *texf6 = new TLatex(14.5,2.,"10 #minus 20");
  texf6->SetTextFont(42);
  texf6->SetTextColor(kBlack);
  texf6->SetTextAlign(22);
  texf6->SetTextSize(0.03);
  texf6->SetTextAngle(270);
  texf6->Draw();

  TLatex *texf7 = new TLatex(14.5,6.,"20 #minus 30");
  texf7->SetTextFont(42);
  texf7->SetTextColor(kBlack);
  texf7->SetTextAlign(22);
  texf7->SetTextSize(0.03);
  texf7->SetTextAngle(270);
  texf7->Draw();

  TLatex *texf8 = new TLatex(14.5,10.,"30 #minus 40");
  texf8->SetTextFont(42);
  texf8->SetTextColor(kBlack);
  texf8->SetTextAlign(22);
  texf8->SetTextSize(0.03);
  texf8->SetTextAngle(270);
  texf8->Draw();

  TLatex *texf12 = new TLatex(15.5,6.,"Jet p_{T} unfolde [GeV]");
  texf12->SetTextFont(42);
  texf12->SetTextColor(kBlack);
  texf12->SetTextAlign(22);
  texf12->SetTextSize(0.04);
  texf12->SetTextAngle(270);
  texf12->Draw();

  TLatex *texf13 = new TLatex(-0.3,-0.3,"0");
  texf13->SetTextFont(42);
  texf13->SetTextColor(kBlack);
  texf13->SetTextAlign(22);
  texf13->SetTextSize(0.03);
  texf13->Draw();

  TLatex *texf14 = new TLatex(-0.3,4.,"1");
  texf14->SetTextFont(42);
  texf14->SetTextColor(kBlack);
  texf14->SetTextAlign(22);
  texf14->SetTextSize(0.03);
  texf14->Draw();

  TLatex *texf15 = new TLatex(4.,-0.3,"1");
  texf15->SetTextFont(42);
  texf15->SetTextColor(kBlack);
  texf15->SetTextAlign(22);
  texf15->SetTextSize(0.03);
  texf15->Draw();
  /*
  //18.,39
  TString cmsText = "CMS";
  //TLatex *latex = new TLatex(44,39.5,cmsText);
  TLatex *latex = new TLatex(4,34.5,cmsText);
  latex->SetTextFont(61);
  latex->SetTextSize(0.055);
  latex->SetTextAlign(22);

  TString simText = "Simulation";
  //TLatex *sim = new TLatex(44,37.5,simText);
  TLatex *sim = new TLatex(12,34.2,simText);
  sim->SetTextFont(52);
  sim->SetTextSize(0.04);
  sim->SetTextAlign(22);
    
  TString extraText = "Preliminary";
  TLatex *extra = new TLatex(3,31,extraText);
  extra->SetTextFont(52);
  extra->SetTextSize(0.024);
  extra->SetTextAlign(22);

  latex->Draw();
  sim->Draw();
  //extra->Draw();


  string coll = "";
  if(doPbPb) coll = "PbPb";
  else coll = "pp";
  //TLatex *texf16 = new TLatex(1.5,35,coll.c_str());
  //if (doPbPb) texf16 = new TLatex(2.6,35,coll.c_str());
  TLatex *texf16 = new TLatex(16,1.6,coll.c_str());
  if (doPbPb) texf16 = new TLatex(15,1.9,coll.c_str());
  texf16->SetTextColor(kBlack);
  texf16->SetTextFont(42);
  texf16->SetTextAlign(22);
  texf16->SetTextSize(0.04);
  texf16->Draw();

  string jpsi = "";
  if(doPrompt) jpsi = "prompt J/#psi";
  else jpsi = "nonprompt J/#psi";

  //TLatex *texf17 = new TLatex(5.,32.75,jpsi.c_str());
  TLatex *texf17 = new TLatex(23,1.6,jpsi.c_str());
  texf17->SetTextColor(kBlack);
  texf17->SetTextFont(42);
  texf17->SetTextAlign(22);
  texf17->SetTextSize(0.04);
  texf17->Draw();

  string rapRange = "|#eta_{jet}| < 2";
  TLatex *texf18 = new TLatex(32,1.4,rapRange.c_str());
  texf18->SetTextFont(42);
  texf18->SetTextColor(kBlack);
  texf18->SetTextAlign(22);
  texf18->SetTextSize(0.04);
  texf18->Draw();
  */
  can->SaveAs(outputfile.c_str());
    
}


void plotMatrix_PAS(){
  
  //plot(true,true);
  plot(true,false);
  //plot(false,true);
  //plot(false,false);
  
}
