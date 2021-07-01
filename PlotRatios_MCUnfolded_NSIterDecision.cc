#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include "inputParams.h"
#endif

void plot(bool doPrompt = false, bool doPbPb = true){
  if (!setCaseTag()) return;
  gSystem->mkdir(Form("%s/mcUnf/plots/",unfPath.c_str()));
  gSystem->mkdir(Form("%s/mcUnf/plots/chi2Test",unfPath.c_str()));
  string filename = "";
  cout<<"nIter = "<<nIter<<endl;
  int iterMin = nSIter-3;
  int iterMax = nSIter+3;
  if (!doPbPb) {
    iterMin = nSIter_pp-5;
    iterMax = nSIter_pp;
  }
  if (iterMin<1) iterMin=1;
  
  TH1D *hZTrue = NULL;
  
  TCanvas * mycan1 = new TCanvas("mycan1","mycan1",900,500);
  mycan1->Divide(2,1);
  
  TLegend *legend = new TLegend(0.8,0.15,.95,0.9,"","brNDC");
  legend->SetHeader("");  
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  float xCoord2 = unfStart;
  TLine *line1(0x0);
  TLine *line2(0x0);

  TH1D* chi2Hist = new TH1D("chi2Hist","",150,0.5,150.5);//99,0.5,99.5);
  
  //for (int iIter = iterMin; iIter<=iterMax; iIter++)
  for (int iIter = 1; iIter<=25; iIter++)
    {
      if (!doPbPb && iIter>nSIter_pp) continue;
      filename = Form("%s/mcUnf/unfOutput/step%d/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBins%s.root", unfPath.c_str(), iIter, doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, caseTag.c_str());
      
      TFile *file = new TFile(filename.c_str());
      if (!file) {cout<<"[WARNING] file not found"<<filename<<endl; continue;}
      
      if (iIter==1) hZTrue = (TH1D*)file->Get(Form("hMTru_%d;1",midLowerId)); TH1D *hZTrueTemp = NULL;
      TH1D *hZUnf = (TH1D*) file->Get(Form("hMUnf_%d_Iter%d;1",midLowerId,nIter)); TH1D *hZUnf_temp = NULL;
      
      for (int i = 1; i <= (midUpperId-midLowerId); i++) {
	if (iIter==1) {
	  hZTrueTemp = (TH1D*)file->Get(Form("hMTru_%d;1",midLowerId+i));
	  hZTrue->Add(hZTrueTemp);
	}
	
	hZUnf_temp = (TH1D*)file->Get(Form("hMUnf_%d_Iter%d;1",midLowerId+i,nIter));
	hZUnf->Add(hZUnf_temp);
      }

      if (iIter==1)
	hZTrue->Rebin(nBinZ_gen/nBinZ_reco);
      hZUnf->Rebin(nBinZ_gen/nBinZ_reco);
  
      mycan1->cd(1);
      //gPad->SetLogy();

      gPad->SetLeftMargin(0.12);
      gPad->SetRightMargin(0.2);

      if (iIter==1) {
	hZTrue->SetStats(0);
	hZTrue->SetTitle(Form("%s %s",doPrompt?"prompt":"nonprompt",doPbPb?"PbPb":"PP"));  
	hZTrue->GetXaxis()->SetTitle("z");
	hZTrue->GetYaxis()->SetTitleOffset(1.6);
	hZTrue->GetYaxis()->SetTitle(Form("dN/dz (per %f)",z_reco_binWidth));
	hZTrue->GetYaxis()->SetRangeUser(0,1.5*hZTrue->GetMaximum());
	hZTrue->SetLineColor(col[0]);
	hZTrue->SetMarkerColor(col[0]);
	hZTrue->SetMarkerStyle(0);
	hZTrue->SetMarkerSize(markerSize[0]);
	hZTrue->SetLineWidth(lineWidth[0]);
	legend->AddEntry(hZTrue, "truth","ep");
	hZTrue->Draw("EP");
      }
      if (iIter>=iterMin && iIter<=iterMax) {
      hZUnf->SetMarkerColor(col[iIter-iterMin+1]);
      hZUnf->SetMarkerStyle(markerStyle[iIter-iterMin+1]);
      hZUnf->SetMarkerSize(markerSize[iIter-iterMin+1]);
      hZUnf->SetLineColor(col[iIter-iterMin+1]);
      hZUnf->SetLineWidth(lineWidth[iIter-iterMin+1]);
      }
      TH1D* hZUnf_chi2 = (TH1D*) hZUnf->Clone("hZUnf_chi2");
      TH1D* hZTrue_chi2 = (TH1D*) hZTrue->Clone("hZTrue_chi2");
      hZUnf_chi2->SetBinContent(1,0);
      hZUnf_chi2->SetBinError(1,0);
      hZTrue_chi2->SetBinContent(1,0);
      hZTrue_chi2->SetBinError(1,0);
      for (int iBin=1;iBin<=nBinZ_reco; iBin++) {
	hZUnf_chi2->SetBinError(iBin,0);
      }

      Double_t chi2 = hZUnf_chi2->Chi2Test(hZTrue_chi2,"WW CHI2/NDF");
      cout <<"[INFO] For SI "<<iIter<<" CHI2/NDF = "<<chi2<<endl;
      chi2Hist->SetBinContent(chi2Hist->FindBin(iIter),chi2);
      if (iIter>=iterMin && iIter<=iterMax) {
      legend->AddEntry(hZUnf, Form("unf SI%d",iIter),"ep");
      hZUnf->Draw("E1same");
      }
      
      if (iIter==iterMax) {
	legend->Draw("same");
	line1 = new TLine(xCoord2,0,xCoord2,gPad->GetUymax());
	line1->SetLineColor(kRed);
	line1->SetLineStyle(2);
	line1->SetLineWidth(1);
	line1->Draw("same");
      }
      if (iIter>=iterMin && iIter<=iterMax) {
      mycan1->Update();
  
      mycan1->cd(2);
      TH1D *hZUnf_diffTruth;
      hZUnf_diffTruth = (TH1D*) hZUnf->Clone();
      hZUnf_diffTruth->Scale(-1);
      hZUnf_diffTruth->Add(hZTrue);
      hZUnf_diffTruth->Scale(-1);
      TH1D *hZUnf_ratioTruth;
      hZUnf_ratioTruth = (TH1D*) hZUnf_diffTruth->Clone();      
      hZUnf_ratioTruth->Divide(hZTrue);
      hZUnf_ratioTruth->SetStats(0);
      hZUnf_ratioTruth->GetYaxis()->SetTitle("(true-unf/true)");
      hZUnf_ratioTruth->GetYaxis()->SetTitleOffset(1.6);
      hZUnf_ratioTruth->GetYaxis()->SetRangeUser(-0.5,0.5);
      hZUnf_ratioTruth->GetXaxis()->SetTitle("z");
      
      if (iIter==iterMin)
	hZUnf_ratioTruth->Draw("E1");
      else
	hZUnf_ratioTruth->Draw("E1same");
      
      mycan1->Update();
      }
      if(iIter==iterMax) {
	line1->Draw("same");
	line2 = new TLine(min_z,0,gPad->GetUxmax(),0);
	line2->SetLineColor(kblue);
 	line2->SetLineStyle(1);
 	line2->SetLineWidth(1);
 	line2->Draw("same");
       }
       mycan1->Update();
     }  
   mycan1->SaveAs(Form("%s/mcUnf/plots/chi2Test/unf_mc_IterComparison_%s_%s_jetR%d_ratioTruth_%dz%dptBins%dz%dptMeasBins%s_%dIter_SIterFrom%dTo%d.pdf",unfPath.c_str(), doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", (int) (jetR*10), nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, caseTag.c_str(),nIter, iterMin,iterMax));
   mycan1->SaveAs(Form("%s/mcUnf/plots/chi2Test/unf_mc_IterComparison_%s_%s_jetR%d_ratioTruth_%dz%dptBins%dz%dptMeasBins%s_%dIter_SIterFrom%dTo%d.png",unfPath.c_str(), doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", (int) (jetR*10), nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, caseTag.c_str(),nIter, iterMin,iterMax));
   TFile* fileChi2Save = new TFile(Form("%s/mcUnf/plots/chi2Test/unf_mc_chi2Hist_%s_%s_jetR%d_ratioTruth_%dz%dptBins%dz%dptMeasBins%s_%dIter_SIterFrom%dTo%d.root",unfPath.c_str(), doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", (int) (jetR*10), nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, caseTag.c_str(),nIter, iterMin,iterMax),"RECREATE");
   chi2Hist->Write();
   fileChi2Save->Close();
 }

 void plotChi2(bool doPrompt, bool doPbPb) {
   if (!setCaseTag()) return;
   gStyle->SetOptStat(0);

   int iter[] = {2,3,4,5,6,7};
   int SIMin[] = {22,22,22,22,22};
   int SIMax[] = {28,28,28,28,28};
   int SIMin_pp[] = {1,1,1,1,1,1};
   int SIMax_pp[] = {4,4,4,4,4,4};
   int nCase = sizeof(iter)/sizeof(int);

  TH1D* histAxis = new TH1D("histAxis","",150,0.5,150.5);
  if (doPbPb) 
    histAxis->GetXaxis()->SetRangeUser(0.5,40.5);
  else 
    histAxis->GetXaxis()->SetRangeUser(0.5,4.5);
  histAxis->GetXaxis()->SetTitle("SI");
  histAxis->GetYaxis()->SetTitle("#chi^{2}");
  histAxis->GetYaxis()->SetTitleOffset(1.3);
  TCanvas* c = new TCanvas("c","",900,900);
  histAxis->Draw();
  TLegend* leg = new TLegend(0.4,0.5,0.85,0.8);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  int max = 0;
  TH1D* chi2HistTemp =NULL;
  TH1D* chi2Hist[nCase];

  for(int i=0; i<nCase; i++) {
    TFile* fileChi2 = TFile::Open(Form("%s/mcUnf/plots/chi2Test/unf_mc_chi2Hist_%s_%s_jetR%d_ratioTruth_%dz%dptBins%dz%dptMeasBins%s_%dIter_SIterFrom%dTo%d.root",unfPath.c_str(), doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", (int) (jetR*10), nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, caseTag.c_str(),iter[i], doPbPb?SIMin[i]:SIMin_pp[i],doPbPb?SIMax[i]:SIMax_pp[i]));
    cout <<"found the file:"<<fileChi2->GetName()<<endl;
    chi2HistTemp = (TH1D*) fileChi2->Get("chi2Hist;1");//new TH1D("chi2Hist","",150,0.5,150.5);
    cout <<"found the histogram:"<<chi2HistTemp->GetName()<<endl;
    chi2Hist[i] = (TH1D*) chi2HistTemp->Clone(Form("chi2Hist_%d",iter[i]));
    if (chi2Hist[i]->GetMaximum()>max) max = chi2Hist[i]->GetMaximum();
    chi2Hist[i]->SetMarkerColor(col[i+1]);
    chi2Hist[i]->SetLineColor(col[i+1]);
    chi2Hist[i]->SetMarkerStyle(markerStyle[i]);
    chi2Hist[i]->Draw("same hist");
    chi2Hist[i]->Draw("same p");
    leg->AddEntry(chi2Hist[i],Form("%d iterations",iter[i]),"p");
  }
  histAxis->GetYaxis()->SetRangeUser(0,1000);
  if (doPbPb && doPrompt) {
  histAxis->GetYaxis()->SetRangeUser(10,max*10);
  c->SetLogy();
  }
  leg->Draw("same");
  c->SaveAs(Form("%s/mcUnf/plots/chi2Test/unf_mc_chi2Plot_%s_%s_jetR%d_ratioTruth_%dz%dptBins%dz%dptMeasBins%s.pdf",unfPath.c_str(), doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", (int) (jetR*10), nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, caseTag.c_str()));
  c->SaveAs(Form("%s/mcUnf/plots/chi2Test/unf_mc_chi2Plot_%s_%s_jetR%d_ratioTruth_%dz%dptBins%dz%dptMeasBins%s.png",unfPath.c_str(), doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", (int) (jetR*10), nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, caseTag.c_str()));
}

void PlotRatios_MCUnfolded_NSIterDecision(){
  //plot(bool doPrompt, bool doPbPb)

  //plot(true,true);
  //if (centShift==0 && nSIter_pp>0 && !doCent && !doPeri)  
  plot(true,false);
  /*
  plot(false,true);
  if (centShift==0 && nSIter_pp>0 && !doCent && !doPeri)
    plot(false,false);
  */
  //plotChi2(1,1);
  plotChi2(1,0);
  //plotChi2(0,1);
  //plotChi2(0,0);
}
