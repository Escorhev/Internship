#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__) || !already
#include "inputParams.h"
#endif

bool removePrelim = false;
bool plotMC = true;
bool plotTrue = true;

TGraphAsymmErrors* systUncertaintyHistAll(bool doPbPb, bool doPrompt, TH1D* nominalHist);
void systUncertaintyHistReg(bool doPbPb, bool doPrompt, TH1D* nominalHist);
void systUncertaintyHistTrStat(bool doPbPb, bool doPrompt);
//void systUncertaintyQuarkoniaSyst(bool doPbPb, bool doPrompt);
void systUncertaintyPriorSyst(bool doPbPb, bool doPrompt, TH1D* nominalHist);
void plotFinalResults_oneStyle(bool doPrompt=true);
void plotFinalResultsZ() {
  plotMC = false;
  //removePrelim = false;
  //plotFinalResults_oneStyle(true);
  removePrelim = true;
  plotFinalResults_oneStyle(true);
  cout<<"goMC"<<endl;
  plotMC = true;
  //removePrelim = false;
  //plotFinalResults_oneStyle(true);
  removePrelim = true;
  plotFinalResults_oneStyle(true);

}

void plotFinalResults_oneStyle(bool doPrompt) {
  gStyle->SetOptStat(0);
    
  string filePPName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_PP_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%s_statError.root",unfPath.c_str(),doPrompt?"prompt":"nonprompt", nIter_pp, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, nSIter_pp, "_nominal");
  
  string fileMCName = Form("%s/dataUnf/unfOutput/finalResults/MC_PP_%dz%dptBins%dz%dptMeasBin_statError.pdf",unfPath.c_str(), nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco);

  string fileOutputName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_PP%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_PbPbnIter%inSIter%i_PPnIter%inSIter%i_statError%s.pdf",unfPath.c_str(),!plotMC?"":"vsMC",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, nIter, nSIter, nIter_pp, nSIter_pp,removePrelim?"_noPreliminaryLabel":"");

  TFile* filePP = TFile::Open(filePPName.c_str());
  TFile* fileMC = TFile::Open(fileMCName.c_str());

  TH1D* histPP = histPP = (TH1D*) filePP->Get(Form("zUnfSI%d",1));
  TH1D* histMC = (TH1D*) fileMC->Get("hZTrueMC");
  /*
  TH1D *cTruth1DX = dynamic_cast<TH1D*>(hist2MC->ProjectionX("TX", 0, -1, "E"));
  TH1D *cTruth1DY = dynamic_cast<TH1D*>(hist2MC->ProjectionY("TY", 0, -1, "E"));

  TCanvas * can6 = new TCanvas("can6","can6",1200,600);
  can6->Divide(2,1);
  can6->cd(1);
  cTruth1DX->raw("EP");
  can6->cd(2);
  cTruth1DY->Draw("EP");
  can6->SaveAs("/home/llr/cms/henderson/UpsUnfolding/Unfolding/dataUnf/code/debughist2MC.png");
  */
  TCanvas * can7 = new TCanvas("can7","can7",1200,600);
  histMC->Draw("EP");
  can7->SaveAs("/home/llr/cms/henderson/UpsUnfolding/Unfolding/dataUnf/code/debughistMC.png");
  
  //double normpp = normPP;
  if (plotMC) normPP=1./histPP->Integral();
  histPP->Scale(normPP*1./z_reco_binWidth);
  histMC->Scale(1./histMC->Integral()*1./z_reco_binWidth);

  histPP->SetTitle("");
  
  TGraphAsymmErrors* graphPP = new TGraphAsymmErrors(histPP);//(nBin,x_pp,y_pp,exl_pp,exh_pp,eyl_pp,eyh_pp);
  TGraphAsymmErrors* graphPPSyst = systUncertaintyHistAll(false, doPrompt, histPP);//new TGraphAsymmErrors(histPPSyst);

  TFile* ftest = new TFile(Form("%s/dataUnf/unfOutput/finalResults/Test_PP%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_PbPbnIter%inSIter%i_PPnIter%inSIter%i_statError%s.pdf",unfPath.c_str(),plotMC?(!plotMC?"":"vsMC"):"",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, nIter, nSIter, nIter_pp, nSIter_pp,removePrelim?"_noPreliminaryLabel":""),"RECREATE");

  graphPP->SetLineColor(kBlack);
  graphPP->SetMarkerColor(kBlack);
  graphPP->Write("graphPP");
  graphPPSyst->SetLineColor(kBlack);
  graphPPSyst->SetMarkerColor(kBlack);
  graphPPSyst->Write("graphPPSyst");
  //ftest->ls();
  ftest->Close();
    
  graphPP->SetLineColor(kpp);
  graphPP->SetMarkerColor(kpp);
  graphPP->SetMarkerStyle(kFullSquare);
  graphPP->SetMarkerSize(1.5);
  graphPP->SetTitle("");
  graphPP->GetXaxis()->SetTitle("z^{2}");
  graphPP->GetYaxis()->SetTitle("#frac{d#sigma^{pp}}{dz^{2}}, #frac{1}{N_{evt}} [nb]");

  //graphPP->GetYaxis()->SetRangeUser(0, graphPP->GetMaximum()*1.5);
  
  graphPPSyst->SetLineColor(kpp);
  graphPPSyst->SetMarkerColor(kpp);
  graphPPSyst->SetMarkerStyle(kFullSquare);
  graphPPSyst->SetFillColorAlpha(kppLight, 0.75);

  histMC->SetLineColor(kred);
  histMC->SetLineWidth(2);
  histMC->SetMarkerColor(kred);
  histMC->SetMarkerStyle(kFullSquare);
  histMC->SetMarkerSize(0);

  TLegend* leg = new TLegend(0.7,0.7,0.93,0.65);
  //if (plotMC) leg = new TLegend(0.45,0.5,0.75,0.65);
  if (plotMC) leg = new TLegend(0.3,0.52,0.76,0.65);
  //if (plotMC) leg = new TLegend(0.42,0.77,0.78,0.89);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  if (!plotMC) {
    leg->AddEntry(graphPP, "pp","lp");
  }
  else {
    if (plotMC) {
    leg->AddEntry(graphPP, "Data","lp");
    leg->AddEntry(histMC, "PYTHIA 8","lp");
    }
    else {
      leg->AddEntry(graphPP, "Data","lp");
    }
  }

  TLatex *  text = new TLatex(0.2 ,0.82,"CMS");
  text->SetNDC();
  text->SetTextFont(61);
  text->SetTextSize(0.075);
  text->SetLineWidth(8);
  
  TLatex *  text0 = new TLatex(0.38 ,0.82,"Preliminary");
  text0->SetNDC();
  text0->SetTextFont(52);
  text0->SetTextSize(0.055);
  text0->SetLineWidth(2);

  TLatex *  text1 = new TLatex(0.2 , 0.76, "#Upsilon (1S)");
  text1->SetNDC();
  text1->SetTextFont(42);
  text1->SetTextSize(0.044);
  text1->SetLineWidth(2);

  TLatex *  text3 = new TLatex(0.2 , 0.71, "20 < p_{T,Jet} < 30 GeV");
  text3->SetNDC();
  text3->SetTextFont(42);
  text3->SetTextSize(0.035);
  text3->SetLineWidth(2);
  
  TLatex *  text4 = new TLatex(0.2 , 0.66, "|#eta_{Jet}| < 2");
  text4->SetNDC();
  text4->SetTextFont(42);
  text4->SetTextSize(0.035);
  text4->SetLineWidth(2);

  TLatex *  text6 = new TLatex(0.344, 0.91,"pp 302 pb^{-1} (5.02 TeV)");
  if (plotMC) text6 = new TLatex(0.57, 0.91,"pp 302 pb^{-1} (5.02 TeV)");
  text6->SetNDC();
  text6->SetTextFont(42);
  text6->SetTextSize(0.037);
  text6->SetLineWidth(2);

  TCanvas* c = new TCanvas("c","",800,800);
  c->cd();
  gPad->SetLeftMargin(0.15);

  TH1D* axisHist = new TH1D("axisHist","",6,0.,1.);
  axisHist->SetTitle("");
  axisHist->GetXaxis()->SetTitle("z^{2}");
  axisHist->GetYaxis()->SetTitle("#frac{d#sigma^{pp}}{dz^{2}} [nb]");
  //cout <<"x axis title size"<<axisHist->GetXaxis()->GetTitleSize()<<endl;
  axisHist->GetXaxis()->SetTitleSize(0.045);
  axisHist->GetYaxis()->SetTitleSize(0.04);
  if (plotMC) axisHist->GetYaxis()->SetTitleSize(0.045);
  axisHist->GetYaxis()->SetTitleOffset(1.6);
  axisHist->GetYaxis()->SetRangeUser(0, histPP->GetMaximum()*1.2);
  if (plotMC) {
    axisHist->GetYaxis()->SetTitle("1/N dN/dz^{2}");
    axisHist->GetYaxis()->SetRangeUser(0, histMC->GetMaximum()*1.2);
  }
  axisHist->Draw();
  if (plotMC) {
    histMC->Draw("same e1");
    histMC->Draw("same hist");
  }
  
  //graphPP->GetYaxis()->SetRangeUser(1,1.2);
  graphPP->GetXaxis()->SetRangeUser(0.,2.0);
  //graphPP->Draw("P");
  graphPPSyst->Draw("ZP");
  graphPPSyst->Draw("5");
  leg->Draw("same");
  text->Draw("same");
  if (!removePrelim)
    text0->Draw("same");
  text1->Draw("same");
  text3->Draw("same");
  text4->Draw("same");
  if(!plotMC)
  text6->Draw("same");

  c->SaveAs(fileOutputName.c_str());
  fileOutputName.replace(fileOutputName.find(".pdf"),4,".png");
  c->SaveAs(fileOutputName.c_str());

  fileOutputName.replace(fileOutputName.find(".png"),4,".pdf");
  fileOutputName.replace(fileOutputName.find("PP"),8,"PPZ");
  //if (plotMC) fileOutputName.replace(fileOutputName.find("vsMC"),4,"");

  Double_t x[4] = {0,0,0,0};
  Double_t y[4] = {0,0,0,0};
  Double_t yMC[4] = {0,0,0,0};
  Double_t exl[4] = {0,0,0,0};
  Double_t eyl[4] = {0,0,0,0};
  Double_t eylMC[4] = {0,0,0,0};
  Double_t exh[4] = {0,0,0,0};
  Double_t eyh[4] = {0,0,0,0};
  Double_t eyhMC[4] = {0,0,0,0};
  Double_t prec = 0;
  Double_t dist = 0;

  for (int iBin = 0; iBin < nBinZ_reco; iBin++) {
    int bin = histPP->FindBin(iBin*z_reco_binWidth+min_z+unfStart);
    if (iBin > 0)
      prec = x[iBin-1];
    if (iBin > 0)
      dist = exh[iBin-1];
    x[iBin] = (sqrt(histPP->GetBinLowEdge(bin) + histPP->GetBinWidth(bin)) + prec + dist)/2.;
    exl[iBin] = x[iBin] - dist - prec;
    //cout<<x[iBin]<<" and "<<dist<<" and "<<prec<<endl;
    exh[iBin] = exl[iBin];
    eyl[iBin] = graphPPSyst->GetErrorYlow(bin-1)*z_reco_binWidth/(exh[iBin]+exl[iBin]);
    eylMC[iBin] = histMC->GetBinError(bin)*z_reco_binWidth/(exh[iBin]+exl[iBin]);
    eyh[iBin] = graphPPSyst->GetErrorYhigh(bin-1)*z_reco_binWidth/(exh[iBin]+exl[iBin]);
    eyhMC[iBin] = histMC->GetBinError(bin)*z_reco_binWidth/(exh[iBin]+exl[iBin]);
    y[iBin] = histPP->GetBinContent(bin)*z_reco_binWidth/(exh[iBin]+exl[iBin]);
    yMC[iBin] = histMC->GetBinContent(bin)*z_reco_binWidth/(exh[iBin]+exl[iBin]);
  }
  TGraphAsymmErrors* graphZ = new TGraphAsymmErrors(nBinZ_reco,x,y,exl,exh,eyl,eyh);
  TGraphAsymmErrors* graphZMC = new TGraphAsymmErrors(nBinZ_reco,x,yMC,exl,exh,eylMC,eyhMC);

  graphZ->SetLineColor(kpp);
  graphZ->SetMarkerColor(kpp);
  graphZ->SetMarkerStyle(kFullSquare);
  graphZ->SetFillColorAlpha(kppLight, 0.75);
  graphZ->SetMarkerSize(1.5);
  graphZ->SetTitle("");
  graphZ->GetXaxis()->SetTitle("z");
  graphZ->GetYaxis()->SetTitle("#frac{d#sigma^{pp}}{dz}, #frac{1}{N_{evt}} [nb]");

  graphZMC->SetLineColor(kred);
  graphZMC->SetLineWidth(2);
  graphZMC->SetMarkerColor(kred);
  graphZMC->SetMarkerStyle(kFullSquare);
  graphZMC->SetMarkerSize(0);

  TH1D* axisHist2 = new TH1D("axisHist2","",6,0.,1.);
  axisHist2->SetTitle("");
  axisHist2->GetXaxis()->SetTitle("z");
  axisHist2->GetYaxis()->SetTitle("#frac{d#sigma^{pp}}{dz} [nb]");
  axisHist2->GetXaxis()->SetTitleSize(0.045);
  axisHist2->GetYaxis()->SetTitleSize(0.04);
  if (plotMC) axisHist2->GetYaxis()->SetTitleSize(0.045);
  axisHist2->GetYaxis()->SetTitleOffset(1.6);
  axisHist2->GetYaxis()->SetRangeUser(0, histPP->GetMaximum()*1.2*0.25/0.14);
  if (plotMC) {
    axisHist2->GetYaxis()->SetTitle("1/N dN/dz");
    axisHist2->GetYaxis()->SetRangeUser(0, histMC->GetMaximum()*1.2*0.25/0.14);
  }
  axisHist2->Draw();

  graphZ->GetXaxis()->SetRangeUser(0.,2.0);
  
  graphZ->Draw("5");
  graphZ->Draw("same ZP");
  if (plotMC) {
    graphZMC->Draw("same e1");
  }
  leg->Draw("same");
  text->Draw("same");
  if (!removePrelim)
    text0->Draw("same");
  text1->Draw("same");
  text3->Draw("same");
  text4->Draw("same");
  text6->Draw("same");

  c->SaveAs(fileOutputName.c_str());
  fileOutputName.replace(fileOutputName.find(".pdf"),4,".png");
  c->SaveAs(fileOutputName.c_str());
  fileOutputName.replace(fileOutputName.find(".png"),4,".root");
  
  TFile* fsave = new TFile(fileOutputName.c_str(),"RECREATE");
  graphPP->SetLineColor(kBlack);
  graphPP->SetMarkerColor(kBlack);
  graphPP->Write("graphPP");
  graphPPSyst->SetLineColor(kBlack);
  graphPPSyst->SetMarkerColor(kBlack);
  graphPPSyst->Write("graphPPSyst");
  //fsave->ls();
  fsave->Close();
  
}
TGraphAsymmErrors* systUncertaintyHistAll(bool doPbPb, bool doPrompt, TH1D* nominalHist) {
  gStyle->SetOptStat(0);
  systUncertaintyHistReg(false, true, nominalHist);
  systUncertaintyHistTrStat(false, true);
  //systUncertaintyQuarkoniaSyst(false, true);
  systUncertaintyPriorSyst(false, true, nominalHist);

  TCanvas* cSyst = new TCanvas("cSyst","",800,800);
  nominalHist->SetLineColor(col[0]);
  nominalHist->SetMarkerColor(col[0]);
  nominalHist->SetMarkerStyle(markerStyle[0]);
  nominalHist->SetStats(0);
  nominalHist->Draw();
  nominalHist->SetTitle("");
  nominalHist->GetXaxis()->SetTitle("z^{2}");
  nominalHist->GetYaxis()->SetTitle("#frac{d#sigma^{pp}}{dz^{2}} [nb]");

  nominalHist->GetYaxis()->SetTitleOffset(1.7);
  cSyst->cd();
  gPad->SetLeftMargin(0.15);

  TCanvas* cSystInd = new TCanvas("cSystInd","",800,800);
  cSystInd->cd();
  gPad->SetLeftMargin(0.15);

  TCanvas* cSystRat = new TCanvas("cSystRat","",800,800);
  cSystRat->cd();
  gPad->SetLeftMargin(0.15);

  double maxSyst = nominalHist->GetMaximum();
  TLegend* legSyst = new TLegend(0.7,0.6,0.9,0.9);
  //if (!doPbPb) legSyst = new TLegend(0.4,0.25,0.6,0.45);
  legSyst->SetBorderSize(0);
  legSyst->SetFillStyle(0);
  legSyst->AddEntry(nominalHist,"nominal","lp");

  TLegend* legRat = new TLegend(0.6,0.5,0.9,0.8);
  if (removePrelim)
    legRat = new TLegend(0.5,0.55,0.9,0.88);

  legRat->SetBorderSize(0);
  legRat->SetFillStyle(0);

  //string systNames [] = {"_JESSyst","_SFSyst","_nominal_centShiftSyst","_Regularization", "_TrStat","_QuarkoniaSyst","_nprPriorSyst"};
  //string systNames [] = {"_QuarkoniaSyst","_JESSyst","_SFSyst","_Regularization","_truePriorSyst","_TrStat","_nominal_centShiftSyst"};
  string systNames [] = {"_Regularization","_truePriorSyst","_TrStat"};

  //string systLabels [] = {"Jet energy scale","Jet energy resolution","Underlying event","Regularization", "Transfer matrix stat.","J/#psi signal extraction","Prior"};
  string systLabels [] = {"Regularisation","Base","Erreurs statistiques"};
  int nSyst = sizeof(systNames)/sizeof(systNames[0]);
  Double_t x[4] = {0,0,0,0};
  Double_t y[4] = {0,0,0,0};
  Double_t exl[4] = {0,0,0,0};
  Double_t eyl[4] = {0,0,0,0};
  Double_t exh[4] = {0,0,0,0};
  Double_t eyh[4] = {0,0,0,0};

  for (int i=0; i<nSyst; i++) {
    string fileUpName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%sUp_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, doPbPb?nSIter:nSIter_pp, systNames[i].c_str());
    string fileDownName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%sDown_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, doPbPb?nSIter:nSIter_pp, systNames[i].c_str());
    TFile* fileUp = TFile::Open(fileUpName.c_str());
    TFile* fileDown = TFile::Open(fileDownName.c_str());
    if (!fileUp || !fileDown) {cout<<"[WARNING] systematic file not found:"<<systNames[i]<<endl; continue;}
    TH1D* histUp = (TH1D*) fileUp->Get(Form("zUnfSI%d",doPbPb?nSIter:1));
    TH1D* histDown = (TH1D*) fileDown->Get(Form("zUnfSI%d",doPbPb?nSIter:1));
    TH1D* histUpRat = NULL;
    TH1D* histDownRat = NULL;

    if (!histUp || !histDown) {cout<<"[WARNING] systematic histogram not found:"<<systNames[i]<<endl; continue;}

    histUp->SetStats(0);
    histDown->SetStats(0);

    int nB = histUp->GetNbinsX();
    for (int iB=0;iB<nB; iB++) {
      if (histUp->GetBinCenter(iB)<unfStart) {
	histUp->SetBinContent(iB,0);
	histUp->SetBinError(iB,0);
	histDown->SetBinContent(iB,0);
	histDown->SetBinError(iB,0);
	continue;
      }
    }
    
    histUp->Scale(normPP*1./z_reco_binWidth);
    histDown->Scale(normPP*1./z_reco_binWidth);

    if (histUp->GetMaximum() > maxSyst) maxSyst = histUp->GetMaximum();
    if (histDown->GetMaximum() > maxSyst) maxSyst = histDown->GetMaximum();

    nominalHist->GetYaxis()->SetRangeUser(0,1.5*maxSyst);

    histUp->SetLineColor(col[i+1]);
    histUp->SetMarkerColor(col[i+1]);
    histUp->SetMarkerStyle(markerStyle[i+1]);
    histUp->SetLineStyle(7);
    histDown->SetLineColor(col[i+1]);
    histDown->SetMarkerColor(col[i+1]);
    histDown->SetMarkerStyle(markerStyle[i+1]);
    histDown->SetLineStyle(8);

    histUpRat = (TH1D*) histUp->Clone("histUpRat"); 
    if (systDiff) {
      histUpRat->Scale(-1);
      histUpRat->Add(nominalHist);
      histUpRat->Scale(-1);
      histUpRat->GetYaxis()->SetTitle("Incertitudes");
      histUpRat->GetYaxis()->SetRangeUser(-0.05e-3,0.12e-3);
    }
    else {
      histUpRat->Divide(nominalHist); 
      histUpRat->GetYaxis()->SetTitle("syst"); 
      histUpRat->GetYaxis()->SetRangeUser(0.6,1.4); 
    }
    histUpRat->SetTitle("");
    histUpRat->GetXaxis()->SetTitle("z^{2}");
    histUpRat->GetXaxis()->SetTitleSize(0.045);
    histUpRat->GetYaxis()->SetTitleSize(0.045); 
    histDownRat = (TH1D*) histDown->Clone("histDownRat"); 
    if (systDiff) {
      histDownRat->Scale(-1);
      histDownRat->Add(nominalHist);
      histDownRat->Scale(-1);
    }
    else 
      histDownRat->Divide(nominalHist);
    histUpRat->SetLineStyle(1);
    histUpRat->SetLineWidth(2);
    histDownRat->SetLineStyle(1);
    histDownRat->SetLineWidth(2);
    histUpRat->GetYaxis()->SetRangeUser(-0.000002, 0.0000045);
    cSyst->cd();
    histUp->Draw("same");
    histDown->Draw("same");
    cSystRat->cd();
    histUpRat->Draw("same hist");
    histDownRat->Draw("same hist");
    cSystInd->cd();
    nominalHist->Draw();
    histUp->Draw("same");
    histDown->Draw("same");
    if (!plotMC)
      cSystInd->SaveAs(Form("%s/dataUnf/unfOutput/finalResults/SystematicDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%s.pdf",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco,doPbPb?nSIter:nSIter_pp,systNames[i].c_str()));
    
    legSyst->AddEntry(histUp,Form("%s",systLabels[i].c_str()),"lp");
    legRat->AddEntry(histUpRat,Form("%s",systLabels[i].c_str()),"l");
    //legSyst->AddEntry(histDown,Form("%s Down",systNames[i].c_str()),"lp");
    for (int iBin = 0; iBin <= nBinZ_reco-1; iBin++) { //-1 for the underflow
      double errUp=0;
      double errDown=0;
      int nBin = nominalHist->FindBin(iBin*z_reco_binWidth+0.0001);
      if (nominalHist->GetBinContent(nBin) == 0) continue;

      if (histUp->GetBinContent(nBin) > nominalHist->GetBinContent(nBin)) {
	if (histDown->GetBinContent(nBin) < nominalHist->GetBinContent(nBin)) {
	  errUp = histUp->GetBinContent(nBin) - nominalHist->GetBinContent(nBin);
	  errDown = nominalHist->GetBinContent(nBin) - histDown->GetBinContent(nBin);
	}
	else {
	  errUp = histUp->GetBinContent(nBin) - nominalHist->GetBinContent(nBin);
	  errDown = histDown->GetBinContent(nBin) - nominalHist->GetBinContent(nBin);
	  errUp = max(errUp,errDown);
	  errDown=0;
	}
      }
      else {
	if (histDown->GetBinContent(nBin) < nominalHist->GetBinContent(nBin)) {
	  errUp = nominalHist->GetBinContent(nBin) - histUp->GetBinContent(nBin);
	  errDown = nominalHist->GetBinContent(nBin) - histDown->GetBinContent(nBin);
	  errDown = max(errUp,errDown);
	  errUp=0;
	}
	else {
	  errUp = nominalHist->GetBinContent(nBin) - histDown->GetBinContent(nBin);
	  errDown = nominalHist->GetBinContent(nBin) - histUp->GetBinContent(nBin);
	}
      }
      errUp = errUp/nominalHist->GetBinContent(nBin);
      errDown = errDown/nominalHist->GetBinContent(nBin);
            
      if (i==0) {
	x[iBin] = nominalHist->GetBinCenter(nBin);
	y[iBin] = nominalHist->GetBinContent(nBin);
	exl[iBin] = z_reco_binWidth/2.;
	exh[iBin] = z_reco_binWidth/2.;
	}
      eyl[iBin] = y[iBin]*sqrt(pow(errDown,2) + pow(eyl[iBin]/y[iBin],2));
      eyh[iBin] = y[iBin]*sqrt(pow(errUp,2) + pow(eyh[iBin]/y[iBin],2));
      //cout << "x = "<<x[iBin]<<", y = "<<y[iBin]<<", eyl = "<<eyl[iBin]<<", eyh = "<<eyh[iBin]<<endl;
    }
  }

  TH1D *histUpNew = new TH1D("histUpNew", "histUpNew", nBinZ_reco, 0., 1.);
  TH1D *histDownNew = new TH1D("histDownNew", "histDownNew", nBinZ_reco, 0., 1.);

  for (int i=0; i<nBinZ_reco;i++) {
    histUpNew->SetBinContent(i+1, nominalHist->GetBinContent(i+1)+eyh[i]);
    histDownNew->SetBinContent(i+1, nominalHist->GetBinContent(i+1)-eyl[i]);
  }

  cSyst->cd();
  histUpNew->Draw("same");
  histDownNew->Draw("same");
  legSyst->Draw("same");
  if(!plotMC) {
    cSyst->SaveAs(Form("%s/dataUnf/unfOutput/finalResults/SystematicDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i.pdf",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco,doPbPb?nSIter:nSIter_pp));
    cSyst->SaveAs(Form("%s/dataUnf/unfOutput/finalResults/SystematicDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i.png",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco,doPbPb?nSIter:nSIter_pp));
  }
  cSystRat->cd();
  TLine* y1line = new TLine(min_z,1,1.0,1);
  if (systDiff) y1line = new TLine(min_z,0,1.0,0);
  y1line->SetLineColor(kBlack);
  y1line->SetLineStyle(1);
  y1line->SetLineWidth(2);
  y1line->Draw("same");

  TLatex *  text = new TLatex(0.2 ,0.82,"CMS");
  text->SetNDC();
  text->SetTextFont(61);
  text->SetTextSize(0.075);
  text->SetLineWidth(8);
  text->Draw("same");  

  TLatex *  text0 = new TLatex(0.38 ,0.82,"Preliminary");
  text0->SetNDC();
  text0->SetTextFont(52);
  text0->SetTextSize(0.055);
  text0->SetLineWidth(2);
  if (!removePrelim)
    text0->Draw("same");  

  TLatex *  text1 = new TLatex(0.2 , 0.76, "#Upsilon (1S)");
  text1->SetNDC();
  text1->SetTextFont(42);
  text1->SetTextSize(0.044);
  text1->SetLineWidth(2);
  text1->Draw("same");  
  /*
  TLatex *  text2 = new TLatex(0.2 , 0.72, "p_{T,J/#psi} > 6.5 GeV");
  text2->SetNDC();
  text2->SetTextFont(42);
  text2->SetTextSize(0.035);
  text2->SetLineWidth(2);
  text2->Draw("same");  

  TLatex *  text3 = new TLatex(0.2 , 0.67, "30 < p_{T,Jet} < 40 GeV");
  text3->SetNDC();
  text3->SetTextFont(42);
  text3->SetTextSize(0.035);
  text3->SetLineWidth(2);
  text3->Draw("same");  
  */
  TLatex *  text4 = new TLatex(0.2 , 0.71, "|#eta_{Jet}| < 2");
  text4->SetNDC();
  text4->SetTextFont(42);
  text4->SetTextSize(0.035);
  text4->SetLineWidth(2);
  text4->Draw("same");  
  /*
  TLatex *  text5 = new TLatex(0.57, 0.91,"pp 302 pb^{-1} (5.02 TeV)");
  text5->SetNDC();
  text5->SetTextFont(42);
  text5->SetTextSize(0.037);
  text5->SetLineWidth(2);
  text5->Draw("same");
  */
  TGraphAsymmErrors* systHist = new TGraphAsymmErrors(nBinZ_reco,x,y,exl,exh,eyl,eyh);

  Double_t yR[4] = {0,0,0,0};
  Double_t exlR[4] = {z_reco_binWidth/2.,z_reco_binWidth/2.,z_reco_binWidth/2.,z_reco_binWidth/2.};
  TGraphAsymmErrors* systHistNorm = new TGraphAsymmErrors(nBinZ_reco,x,yR,exlR,exlR,eyl,eyh);
  systHistNorm->SetLineColor(kBlack);
  systHistNorm->SetLineWidth(2);
  systHistNorm->SetMarkerColor(kBlack);
  systHistNorm->SetMarkerStyle(kOpenCircle);
  systHistNorm->SetMarkerSize(0);
  systHistNorm->SetFillColorAlpha(kWhite, 0);
  
  systHistNorm->Draw("5");
  legRat->AddEntry(systHistNorm,"Total","l");
  legRat->Draw("same");
  if(!plotMC)
    cSystRat->SaveAs(Form("%s/dataUnf/unfOutput/finalResults/Systematic%sUncertainty_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%s.png",unfPath.c_str(),systDiff?"Diff":"Rel",doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco,doPbPb?nSIter:nSIter_pp,removePrelim?"_noPreliminaryLabel":""));
  
  return systHist;
}

void systUncertaintyHistReg(bool doPbPb, bool doPrompt, TH1D* nominalHist) {
  gStyle->SetOptStat(0);
  string fileSystName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%s_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", doPbPb?nIter_syst:nIter_syst_pp, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, doPbPb?nSIter_syst:nSIter_syst_pp, "_nominal");
  TFile* fileSyst = TFile::Open(fileSystName.c_str());
  if (!fileSyst) {cout<<"[WARNING] file "<<fileSystName<<" not found for regularization syst"<<endl; return;}
  TH1D* systHist = (TH1D*) fileSyst->Get(Form("zUnfSI%d",doPbPb?nSIter_syst:1));
  TH1D* histUp = (TH1D*) systHist->Clone("histUp");
  TH1D* histDown = (TH1D*) systHist->Clone("histDown");
  TH1D* nomHist = (TH1D*) nominalHist->Clone("nomHist");
  //systHist->Scale(doPbPb?(normPbPb*1./z_reco_binWidth):(normPP*1./z_reco_binWidth));
  nomHist->Scale(doPbPb?(z_reco_binWidth*1./normPbPb):(z_reco_binWidth*1./normPP));
  int nBin = systHist->GetNbinsX();
  for (int i=0; i<=nBin;i++) {
    double binCenter = systHist->GetBinCenter(i);
    double err = fabs(nomHist->GetBinContent(nomHist->FindBin(binCenter))-systHist->GetBinContent(i));
    histUp->SetBinContent(i, nomHist->GetBinContent(i)+err);
    histDown->SetBinContent(i, nomHist->GetBinContent(i)-err);
  }
  
  string fileUpName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%sUp_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, doPbPb?nSIter:nSIter_pp, "_Regularization");
  string fileDownName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%sDown_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, doPbPb?nSIter:nSIter_pp, "_Regularization");
  TFile* fileUp = new TFile(fileUpName.c_str(),"RECREATE");
  histUp->Write(Form("zUnfSI%d",doPbPb?nSIter:1));
  nomHist->Write("nominal");
  systHist->Write("syst");
  fileUp->Close();
  TFile* fileDown = new TFile(fileDownName.c_str(),"RECREATE");
  histDown->Write(Form("zUnfSI%d",doPbPb?nSIter:1));
  fileDown->Close();
}

void systUncertaintyHistTrStat(bool doPbPb, bool doPrompt) {
  gStyle->SetOptStat(0);
  string fileSystName = Form("%s/dataUnf/unfOutput/matrixOper/systUnc_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin%s.root",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, "_nominal");
  cout <<"reading systematics from file "<<fileSystName<<endl;
  TFile* fileSyst = TFile::Open(fileSystName.c_str());
  if (!fileSyst) {cout<<"[WARNING] file "<<fileSystName<<" not found for trM stat syst"<<endl; return;}
  TH1D* systHist = (TH1D*) fileSyst->Get("zUnf_trMatrixSyst");
  TH1D* histUp = (TH1D*) systHist->Clone("histUp");
  TH1D* histDown = (TH1D*) systHist->Clone("histDown");
  int nBin = systHist->GetNbinsX();
  for (int i=0; i<=nBin;i++) {
    histUp->SetBinContent(i, histUp->GetBinContent(i)+histUp->GetBinError(i));
    histDown->SetBinContent(i, histDown->GetBinContent(i)-histDown->GetBinError(i));
  }
  string fileUpName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%sUp_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, doPbPb?nSIter:nSIter_pp, "_TrStat");
  string fileDownName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%sDown_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, doPbPb?nSIter:nSIter_pp, "_TrStat");
  TFile* fileUp = new TFile(fileUpName.c_str(),"RECREATE");
  histUp->Write(Form("zUnfSI%d",doPbPb?nSIter:1));
  fileUp->Close();
  TFile* fileDown = new TFile(fileDownName.c_str(),"RECREATE");
  histDown->Write(Form("zUnfSI%d",doPbPb?nSIter:1));
  fileDown->Close();
}
void systUncertaintyQuarkoniaSyst(bool doPbPb, bool doPrompt) {
  gStyle->SetOptStat(0);
  string fileSystName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%s_systError.root",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, doPbPb?nSIter:nSIter_pp, "_nominal");
  TFile* fileSyst = TFile::Open(fileSystName.c_str());
  if (!fileSyst) {cout<<"[WARNING] file "<<fileSystName<<" not found for quarkonia syst"<<endl; return;}
  TH1D* systHist = (TH1D*) fileSyst->Get(Form("zUnfSI%d",doPbPb?nSIter:1));
  TH1D* histUp = (TH1D*) systHist->Clone("histUp");
  TH1D* histDown = (TH1D*) systHist->Clone("histDown");
  int nBin = systHist->GetNbinsX();
  for (int i=0; i<=nBin;i++) {
    histUp->SetBinContent(i, histUp->GetBinContent(i)+histUp->GetBinError(i));
    histDown->SetBinContent(i, histDown->GetBinContent(i)-histDown->GetBinError(i));
  }

  string fileUpName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%sUp_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, doPbPb?nSIter:nSIter_pp, "_QuarkoniaSyst");
  string fileDownName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%sDown_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, doPbPb?nSIter:nSIter_pp, "_QuarkoniaSyst");
  TFile* fileUp = new TFile(fileUpName.c_str(),"RECREATE");
  histUp->Write(Form("zUnfSI%d",doPbPb?nSIter:1));
  fileUp->Close();
  TFile* fileDown = new TFile(fileDownName.c_str(),"RECREATE");
  histDown->Write(Form("zUnfSI%d",doPbPb?nSIter:1));
  fileDown->Close();
}

void systUncertaintyPriorSyst(bool doPbPb, bool doPrompt, TH1D* nominalHist) {
  gStyle->SetOptStat(0);
  string fileSystName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%s_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, doPbPb?nSIter:nSIter_pp, "_truePriorSyst");
  TFile* fileSyst = TFile::Open(fileSystName.c_str());
  if (!fileSyst) {cout<<"[WARNING] file "<<fileSystName<<" not found for regularization syst"<<endl; return;}
  TH1D* systHist = (TH1D*) fileSyst->Get(Form("zUnfSI%d",doPbPb?nSIter:1));
  TH1D* histUp = (TH1D*) systHist->Clone("histUp");
  TH1D* histDown = (TH1D*) systHist->Clone("histDown");
  TH1D* nomHist = (TH1D*) nominalHist->Clone("nomHist");
  //systHist->Scale(doPbPb?(normPbPb*1./z_reco_binWidth):(normPP*1./z_reco_binWidth));
  nomHist->Scale(doPbPb?(z_reco_binWidth*1./normPbPb):(z_reco_binWidth*1./normPP));
  int nBin = systHist->GetNbinsX();
  for (int i=0; i<=nBin;i++) {
    double binCenter = systHist->GetBinCenter(i);
    double err = fabs(nomHist->GetBinContent(nomHist->FindBin(binCenter))-systHist->GetBinContent(i));
    histUp->SetBinContent(i, nomHist->GetBinContent(i)+err);
    histDown->SetBinContent(i, nomHist->GetBinContent(i)-err);
  }
  
  string fileUpName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%sUp_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, doPbPb?nSIter:nSIter_pp, "_truePriorSyst");
  string fileDownName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%sDown_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, doPbPb?nSIter:nSIter_pp, "_truePriorSyst");
  TFile* fileUp = new TFile(fileUpName.c_str(),"RECREATE");
  histUp->Write(Form("zUnfSI%d",doPbPb?nSIter:1));
  nomHist->Write("nominal");
  systHist->Write("syst");
  fileUp->Close();
  TFile* fileDown = new TFile(fileDownName.c_str(),"RECREATE");
  histDown->Write(Form("zUnfSI%d",doPbPb?nSIter:1));
  fileDown->Close();
}
