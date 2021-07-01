#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include "inputParams.h"
#endif

void unfold(bool doPrompt = true, bool doPbPb = true, Int_t iterMax = 8, Int_t stepNumber = 1) {
  if (!setCaseTag()) return;
  Int_t iterMin =1;
  Int_t iterDef = iterMax - 2;
  
  string testInputName = "";
  string trainInputName = "";
  string outputName = "";

  testInputName = Form("%s/mcUnf/unfInput1D/step%i/response_4D_%s_%s_Test_%diter_%dz%dptBins%dz%dptMeasBins%s.root",unfPath.c_str(), stepNumber,doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", iterMax, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, caseTag.c_str());
  trainInputName = Form("%s/mcUnf/unfInput1D/step%i/response_4D_%s_%s_Train_%diter_%dz%dptBins%dz%dptMeasBins%s.root",unfPath.c_str(),stepNumber,doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", iterMax, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, caseTag.c_str());
  outputName = Form("%s/mcUnf/unfOutput1D/step%i/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBins%s.root",unfPath.c_str(),stepNumber,doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", iterMax, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, caseTag.c_str());
  
  TFile *f_measured = new TFile(testInputName.c_str());
  f_measured->ls();
  TH1D *temporaryMeasured = (TH1D*)f_measured->Get("fh1RespDimM;1");
  TH1D *temporaryTrueInd = (TH1D*)f_measured->Get("fh1RespDimT;1");
  TH2D *hMeasured = new TH2D ("hMeasured", "hMeasured", nBinZ_reco, min_z, max_z, nBinJet_reco, min_jetpt, max_jetpt);
  TH2D *hTrueInd = new TH2D ("hTrueInd", "hTrueInd", nBinZ_gen, min_z, max_z, nBinJet_gen, min_jetpt, max_jetpt);
  
  TCanvas * can1 = new TCanvas("can1","can1",1200,600);
  temporaryMeasured->Draw("EP");
  can1->SaveAs("/home/llr/cms/henderson/UpsUnfolding/Unfolding/mcUnf/code/debugMeasured.png");
  TCanvas * can3 = new TCanvas("can3","can3",1200,600);
  temporaryTrueInd->Draw("EP");
  can3->SaveAs("/home/llr/cms/henderson/UpsUnfolding/Unfolding/mcUnf/code/debugTrueInd.png");

  for(Int_t i = 1; i<=nBinZ_reco; i++) {
    for(Int_t j = 1; j<=nBinJet_reco; j++) {
      Double_t w = temporaryMeasured->GetBinContent(i + nBinZ_reco*(j-1));
      hMeasured->Fill((i-1)*0.25+0.001, 10.001+(j-1)*10, w);
    }
  }

  for(Int_t i = 1; i<=nBinZ_gen; i++) {
    for(Int_t j = 1; j<=nBinJet_gen; j++) {
      Double_t w = temporaryTrueInd->GetBinContent(i + nBinZ_gen*(j-1));
      hTrueInd->Fill((i-1)/32.+0.001, 10.001+(j-1)*2, w);
    }
  }

  TH1D *gTruth1DX = dynamic_cast<TH1D*>(hTrueInd->ProjectionX("TX", 0, -1, "E"));
  TH1D *gTruth1DY = dynamic_cast<TH1D*>(hTrueInd->ProjectionY("TY", 0, -1, "E"));

  TCanvas * can4 = new TCanvas("can4","can4",1200,600);
  can4->Divide(2,1);
  can4->cd(1);
  gTruth1DX->Draw("EP");
  can4->cd(2);
  gTruth1DY->Draw("EP");
  can4->SaveAs("/home/llr/cms/henderson/UpsUnfolding/Unfolding/mcUnf/code/debughTrueInd.png");

  TH1D *gMeas1DX = dynamic_cast<TH1D*>(hMeasured->ProjectionX("TX", 0, -1, "E"));
  TH1D *gMeas1DY = dynamic_cast<TH1D*>(hMeasured->ProjectionY("TY", 0, -1, "E"));

  TCanvas * can5 = new TCanvas("can5","can5",1200,600);
  can5->Divide(2,1);
  can5->cd(1);
  gMeas1DX->Draw("EP");
  can5->cd(2);
  gMeas1DY->Draw("EP");
  can5->SaveAs("/home/llr/cms/henderson/UpsUnfolding/Unfolding/mcUnf/code/debughMeasured.png");


  //////// if we want to smear MC errors to match data
    if (smearMeas) {
      if (stepNumber==1) {
	cout <<"[INFO] Smearing the MC measured distribution to match data stats"<<endl;
	TFile* dataFile = TFile::Open(Form("%s/dataUnf/data_results/meas_%s_data_%s%s_statErrs.root",unfPath.c_str(), doPbPb?"PbPb":"PP",doPrompt?"prompt":"prompt",doCent?"_centBin":doPeri?"_periBin":"")); //use prompt data for prompt and nonprompt
	TRandom* rnd = new TRandom3();
	TH2D *hMeasuredData = (TH2D*) dataFile->Get("h_Meas;1");
	TFile* saveFile = new TFile(Form("%s/mcUnf/unfOutput/step%i/statSmearing_%s_%s_%diter_%dz%dptBins%dz%dptMeasBins%s.root",unfPath.c_str(),stepNumber,doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", iterMax, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, caseTag.c_str()),"RECREATE");
	saveFile->cd();
	hMeasuredData->Write("data");
	hMeasured->Write("mc");
	int nBinsX = hMeasuredData->GetNbinsX();
	int nBinsY = hMeasuredData->GetNbinsY();
	TF1* fRand = new TF1("fRand","TMath::Gaus(x,[0],[1],1)",-100,100);
	fRand->SetNpx(10000);
	for (int ix = 0; ix<=nBinsX; ix++) {
	  for (int iy = 0; iy<=nBinsY; iy++) {
	    if (hMeasuredData->GetXaxis()->GetBinCenter(ix)<min_z || hMeasuredData->GetXaxis()->GetBinCenter(ix)>max_z) continue;
	    if (hMeasuredData->GetYaxis()->GetBinCenter(iy)<min_jetpt || hMeasuredData->GetYaxis()->GetBinCenter(iy)>max_jetpt ) continue;
	    int bin = hMeasuredData->GetBin(ix,iy);
	    if (hMeasuredData->GetBinContent(bin)==0) continue;
	    double dataErr = hMeasuredData->GetBinError(bin)*1.0/hMeasuredData->GetBinContent(bin);
	    double mcVal = hMeasured->GetBinContent(bin);
	    if (dataDist) mcVal = hMeasuredData->GetBinContent(bin);
	    double mcErr = dataErr*mcVal;
	    fRand->SetRange(mcVal-6.*mcErr,mcVal+6.*mcErr);
	    fRand->SetParameter(0,mcVal);
	    fRand->SetParameter(1,mcErr);
	    double newMcVal = fRand->GetRandom();	  
	    while (newMcVal<=0.) {newMcVal = fRand->GetRandom();}
	    cout <<"dataVal = "<<hMeasuredData->GetBinContent(bin)<<", dataErr = "<<hMeasuredData->GetBinError(bin)<<", rel dataErr = "<<dataErr<<", mcVal = "<<mcVal<<", newMcVal = "<<newMcVal<<", old rel mcErr = "<<hMeasured->GetBinError(bin)*1.0/hMeasured->GetBinContent(bin)<<", set as "<<dataErr*newMcVal<<endl;
	    hMeasured->SetBinContent(bin, newMcVal);
	    hMeasured->SetBinError(bin, dataErr*newMcVal);
	    cout <<"new mc error = "<<hMeasured->GetBinError(bin)/hMeasured->GetBinContent(bin)<<endl;
	  }
	}
	hMeasured->Write("smearedMC");
	saveFile->Close();
      }
      else {
	TFile* inFile = TFile::Open(Form("%s/mcUnf/unfOutput/step1/statSmearing_%s_%s_%diter_%dz%dptBins%dz%dptMeasBins%s.root",unfPath.c_str(), doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", iterMax, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, caseTag.c_str()));
	hMeasured = (TH2D*) inFile->Get("smearedMC");
      }
  }
    
  TFile *f = new TFile(trainInputName.c_str());
  RooUnfoldResponse *resp = (RooUnfoldResponse*)f->Get("resp;1");
  TH2F *hSmear = (TH2F*)f->Get("fh2RespDimM;1");
  TH2F *hTrue = (TH2F*)f->Get("fh2RespDimT;1");

  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;//kNoError;//
  //RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kNoError;
  
  Int_t nIterTmp = iterMax - iterMin+1;
  const Int_t nIter = nIterTmp;

  TH2D * cheatReco =new  TH2D ("cheatReco","cheatReco",nBinZ_reco, min_z, max_z, nBinJet_reco, min_jetpt, max_jetpt);
  TH2D * cheatGen =new  TH2D ("cheatGen","cheatGen",nBinZ_gen, min_z, max_z, nBinJet_gen, min_jetpt, max_jetpt);

  TH2D *hReco[nIter];
  TH1D *temporaryReco[nIter];
  TH2D *hFolded[nIter];
  TH1D *temporaryFolded[nIter];
  TMatrixD covmat[nIter];
  RooUnfoldBayes unfold[nIter];

  TH1D *hJetPtTrue;
  TH1D *hJetPtMeas;
  TH1D *hJetPtUnf[nIter];
  TH1D *hJetPtFol[nIter];
  
  for(Int_t iter = iterMin; iter<=iterMax; iter++) {

    cout << "iter = " << iter << endl;
    
    // RooUnfoldBayes    unfold (resp, hMeas, iter);
    unfold[iter-iterMin] = RooUnfoldBayes(resp, temporaryMeasured, iter);

    // unfold using response matrix
    
    temporaryReco[iter-iterMin] = (TH1D*)unfold[iter-iterMin].Hreco(errorTreatment);
    unfold[iter-iterMin].Print();

    //unfold[iter-iterMin].PrintTable(cout,temporaryTrueInd,RooUnfold::kCovariance);

    //fold using just unfolded result
    temporaryFolded[iter-iterMin] = (TH1D*)resp->ApplyToTruth(temporaryReco[iter-iterMin]);

    //hFolded[iter-iterMin]->Sumw2();

    // covariance matrix problem here : resize it before assigning new matrix.

    cout << "before covmat[iter-iterMin]" << endl;
    //covmat[iter-iterMin].ResizeTo(unfold[iter-iterMin].Ereco(RooUnfold::kCovariance));
    //covmat[iter-iterMin] = unfold[iter-iterMin].Ereco(RooUnfold::kCovariance);
    //cout << "cols : "<< covmat[iter-iterMin].GetNcols() << endl;
    //cout << "raws : "<< covmat[iter-iterMin].GetNrows() << endl;
    cout << "after covmat[iter-iterMin]" << endl;
    
    //cout << "covmat[iter-iterMin](0,0) " << covmat[iter-iterMin](0,0) << endl;
    //cout << "covmat[iter-iterMin](1,20) " << covmat[iter-iterMin](1,20) << endl;

    /*
    TMatrixD errmat(ntbins,ntbins);
    for (Int_t i=0; i<ntbins; i++)
      for (Int_t j=0; j<ntbins; j++)
	errmat(i,j)= covmat(i,j)>=0 ? sqrt(covmat(i,j)) : -sqrt(-covmat(i,j));
    PrintMatrix(errmat,"","covariance matrix",10);
    */

    hReco[iter-iterMin] = cheatGen;
    hReco[iter-iterMin]->SetNameTitle(Form("hReco_Iter%d",iter), Form("hReco_Iter%d",iter));
    for(Int_t i = 1; i<=nBinZ_gen; i++) {
      for(Int_t j = 1; j<=nBinJet_gen; j++) {
	Double_t w = temporaryReco[iter-iterMin]->GetBinContent(i + nBinZ_gen*(j-1));
	hReco[iter-iterMin]->Fill((i-1)/32.+0.001, 10.001+(j-1)*2., w);
      }
    }
    
    hFolded[iter-iterMin] = cheatReco;
    hFolded[iter-iterMin]->SetNameTitle(Form("hFolded_Iter%d",iter), Form("hFolded_Iter%d",iter));

    for(Int_t i = 1; i<=nBinZ_reco; i++) {
      for(Int_t j = 1; j<=nBinJet_reco; j++) {
        Double_t w = temporaryFolded[iter-iterMin]->GetBinContent(i + nBinZ_reco*(j-1));
        hFolded[iter-iterMin]->Fill((i-1)*0.25+0.001, 10.001+(j-1)*10, w);
      }
    }

    // jet pt distributions save here
    
    hJetPtUnf[iter-iterMin] = dynamic_cast<TH1D*>(hReco[iter-iterMin]->ProjectionY(Form("hJetPtUnf_Iter%d",iter)));
    hJetPtFol[iter-iterMin] = dynamic_cast<TH1D*>(hFolded[iter-iterMin]->ProjectionY(Form("hJetPtFol_Iter%d",iter)));

    //hJetPtUnf[iter-iterMin]->SetName(Form("hJetPtUnf_Iter%d",iter));
    //hJetPtFol[iter-iterMin]->SetName(Form("hJetPtFol_Iter%d",iter));

    //hJetPtUnf[iter-iterMin]->Sumw2();
    //hJetPtFol[iter-iterMin]->Sumw2();
    
  }

  // response matrix folded with true distribution (test sample)
  
  TH1D *cTruth1DX = dynamic_cast<TH1D*>(hReco[2]->ProjectionX("TX", 0, -1, "E"));
  TH1D *cTruth1DY = dynamic_cast<TH1D*>(hReco[2]->ProjectionY("TY", 0, -1, "E"));

  TCanvas * can6 = new TCanvas("can6","can6",1200,600);
  can6->Divide(2,1);
  can6->cd(1);
  cTruth1DX->Draw("EP");
  can6->cd(2);
  cTruth1DY->Draw("EP");
  can6->SaveAs("/home/llr/cms/henderson/UpsUnfolding/Unfolding/mcUnf/code/debughReco2.png");

  TCanvas * can8 = new TCanvas("can8","can8",1200,600);
  temporaryReco[6]->Draw("EP");
  can8->SaveAs("/home/llr/cms/henderson/UpsUnfolding/Unfolding/mcUnf/code/debugReco2.png");
  
  if (stepNumber == 1) {
    
    TH1D *cTruth5DX = dynamic_cast<TH1D*>(hReco[0]->ProjectionX("TX", 0, -1, "E"));
    TH1D *cTruth5DY = dynamic_cast<TH1D*>(hReco[0]->ProjectionY("TY", 0, -1, "E"));

    TCanvas * can11 = new TCanvas("can11","can11",1200,600);
    can11->Divide(2,1);
    can11->cd(1);
    cTruth5DX->Draw("EP");
    can11->cd(2);
    cTruth5DY->Draw("EP");
    can11->SaveAs("/home/llr/cms/henderson/UpsUnfolding/Unfolding/mcUnf/code/debughReco0.png");

    TCanvas * can10 = new TCanvas("can10","can10",1200,600);
    temporaryReco[0]->Draw("EP");
    can10->SaveAs("/home/llr/cms/henderson/UpsUnfolding/Unfolding/mcUnf/code/debugReco0.png");
  }

  TH1D *cTruth2DX = dynamic_cast<TH1D*>(hFolded[0]->ProjectionX("TX", 0, -1, "E"));
  TH1D *cTruth2DY = dynamic_cast<TH1D*>(hFolded[0]->ProjectionY("TY", 0, -1, "E"));

  TCanvas * can7 = new TCanvas("can7","can7",1200,600);
  can7->Divide(2,1);
  can7->cd(1);
  cTruth2DX->Draw("EP");
  can7->cd(2);
  cTruth2DY->Draw("EP");
  can7->SaveAs("/home/llr/cms/henderson/UpsUnfolding/Unfolding/mcUnf/code/debughFolded0.png");
  
  TH2D *hPriorFolded = new TH2D ("hPriorFolded", "hPriorFolded", nBinZ_reco, min_z, max_z, nBinJet_reco, min_jetpt, max_jetpt);
  TH1D *temporaryPriorFolded = (TH1D*)resp->ApplyToTruth(temporaryTrueInd);

  for(Int_t i = 1; i<=nBinZ_reco; i++) {
    for(Int_t j = 1; j<=nBinJet_reco; j++) {
      Double_t w = temporaryPriorFolded->GetBinContent(i + nBinZ_reco*(j-1));
      hPriorFolded->Fill((i-1)*0.25+0.001, 10.001+(j-1)*10, w);
    }
  }

  TH1D *gPrior1DX = dynamic_cast<TH1D*>(hPriorFolded->ProjectionX("TX", 0, -1, "E"));
  TH1D *gPrior1DY = dynamic_cast<TH1D*>(hPriorFolded->ProjectionY("TY", 0, -1, "E"));

  TCanvas * can2 = new TCanvas("can2","can2",1200,600);
  can2->Divide(2,1);
  can2->cd(1);
  gPrior1DX->Draw("EP");
  can2->cd(2);
  gPrior1DY->Draw("EP");
  can2->SaveAs("/home/llr/cms/henderson/UpsUnfolding/Unfolding/mcUnf/code/debughPriorFolded.png");

  
  Int_t min_PriorFolded = hPriorFolded->GetYaxis()->FindBin(midLowerPt+0.00001);
  Int_t max_PriorFolded = hPriorFolded->GetYaxis()->FindBin(midUpperPt-0.00001);
  TH1D *hPriorFolded_x = dynamic_cast<TH1D*>(hPriorFolded->ProjectionX("hPriorFolded_x",min_PriorFolded,max_PriorFolded));

  // 1d projections from unfolded 2d distribution
  TH1D *hRecoP[2];
  hRecoP[0] = dynamic_cast<TH1D*>(hReco[iterDef-iterMin+2]->ProjectionX("hRecoP_x"));
  hRecoP[1] = dynamic_cast<TH1D*>(hReco[iterDef-iterMin+2]->ProjectionY("hRecoP_y"));

  // 1d projections from folded 2d distribution (obtained using new unfolded distribution)
  TH1D *hFoldedP[2];
  hFoldedP[0] = dynamic_cast<TH1D*>(hFolded[iterDef-iterMin+2]->ProjectionX("hFoldedP_x"));
  hFoldedP[1] = dynamic_cast<TH1D*>(hFolded[iterDef-iterMin+2]->ProjectionY("hFoldedP_y"));

  // 1d projections of measured z and jet pt (test sample)
  TH1D *hSP[2];
  hSP[0] = dynamic_cast<TH1D*>(hMeasured->ProjectionX("hSP_x"));
  hSP[1] = dynamic_cast<TH1D*>(hMeasured->ProjectionY("hSP_y"));

  Int_t min_Measured = hMeasured->GetYaxis()->FindBin(midLowerPt+0.00001);
  Int_t max_Measured = hMeasured->GetYaxis()->FindBin(midUpperPt-0.00001);
  TH1D *hMeasured_x = dynamic_cast<TH1D*>(hMeasured->ProjectionX("hMeasured_x",min_Measured,max_Measured));
    
  // 1d projections of true z and jet pt (test sampel)
  TH1D *hTP[2];
  hTP[0] = dynamic_cast<TH1D*>(hTrueInd->ProjectionX("hTP_x"));
  hTP[1] = dynamic_cast<TH1D*>(hTrueInd->ProjectionY("hTP_y"));

  //hSP[0]->Sumw2();
  //hFoldedP[0]->Sumw2();

  //hSP[1]->Sumw2();
  //hFoldedP[1]->Sumw2();

  //hTP[0]->Sumw2();
  //hRecoP[0]->Sumw2();

  //hTP[1]->Sumw2();
  //hRecoP[1]->Sumw2();

  hJetPtTrue = (TH1D*)hTP[1]->Clone();
  hJetPtTrue->SetName("hJetPtTrue");
  hJetPtMeas = (TH1D*)hSP[1]->Clone();
  hJetPtMeas->SetName("hJetPtMeas");

  Double_t *ptmin = new Double_t[nBinJet_gen];
  Double_t *ptmax = new Double_t[nBinJet_gen];

  for (int i=0;i<nBinJet_gen;i++) {
      ptmin[i] = min_jetpt_real+i*jetPt_gen_binWidth;
      ptmax[i] = min_jetpt_real+(i+1)*jetPt_gen_binWidth;
    }

  //Double_t ptmin[nBinJet_gen] = {20.,22.,24.,26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.};
  //Double_t ptmax[nBinJet_gen] = {22.,24.,26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.,50.};

  TH1D *hMUnf[nBinJet_gen][nIter];
  TH1D *hMFol[nBinJet_gen][nIter];
  TH1D *hMTru[nBinJet_gen];
  TH1D *hMSme[nBinJet_gen];
  TH1D *hMPri[nBinJet_gen];

  for(Int_t i = 0; i<nBinJet_gen; i++) {
    Int_t min = hTrueInd->GetYaxis()->FindBin(ptmin[i]+0.00001);
    Int_t max = hTrueInd->GetYaxis()->FindBin(ptmax[i]-0.00001);
    hMTru[i] = dynamic_cast<TH1D*>(hTrueInd->ProjectionX(Form("hMTru_%d",i),min,max));

    min = hMeasured->GetYaxis()->FindBin(ptmin[i]+0.00001);
    max = hMeasured->GetYaxis()->FindBin(ptmax[i]-0.00001);
    hMSme[i] = dynamic_cast<TH1D*>(hMeasured->ProjectionX(Form("hMSme_%d",i),min,max));
    
    for(Int_t iter = iterMin; iter<=iterMax; iter++) {
      Int_t min2 = hReco[iter-iterMin]->GetYaxis()->FindBin(ptmin[i]+0.00001);
      Int_t max2 = hReco[iter-iterMin]->GetYaxis()->FindBin(ptmax[i]-0.00001);
      hMUnf[i][iter-iterMin] = dynamic_cast<TH1D*>(hReco[iter-iterMin]->ProjectionX(Form("hMUnf_%d_Iter%d",i,iter),min2,max2));

      min2 = hFolded[iter-iterMin]->GetYaxis()->FindBin(ptmin[i]+0.00001);
      max2 = hFolded[iter-iterMin]->GetYaxis()->FindBin(ptmax[i]-0.00001);
      hMFol[i][iter-iterMin] = dynamic_cast<TH1D*>(hFolded[iter-iterMin]->ProjectionX(Form("hMFol_%d_Iter%d",i,iter),min2,max2));
    }
  }


  TFile *fout = new TFile(outputName.c_str(),"RECREATE");
  hTrueInd->Write("fh2Prior");
  hMeasured->Write("fh2Smear");
  resp->Write();
  hPriorFolded->Write();

  hJetPtTrue->Write();
  hJetPtMeas->Write();
  
  for(Int_t iter = iterMin; iter<=iterMax; iter++) {
    hReco[iter-iterMin]->Write();
    hFolded[iter-iterMin]->Write();
    covmat[iter-iterMin].Write(Form("covmat%d",iter));
    unfold[iter-iterMin].Write(Form("unfold%d",iter));
    hJetPtUnf[iter-iterMin]->Write();
    hJetPtFol[iter-iterMin]->Write();
  }


  for(Int_t i = 0; i<nBinJet_gen; i++) {
    //hMPri[i]->Write();
    hMSme[i]->Write();
    hMTru[i]->Write();
    for(Int_t iter = iterMin; iter<=iterMax; iter++) {
      if(hMUnf[i][iter-iterMin]) hMUnf[i][iter-iterMin]->Write();
      if(hMFol[i][iter-iterMin]) hMFol[i][iter-iterMin]->Write();
    }
  }

  
  fout->Write();
  fout->Close();
}

void unfoldStep1D(Int_t step = 1){
  //prompt
  if (step<=nSIter_pp && centShift == 0 && !doCent && !doPeri)
    unfold(true,false,nIter,step);
}
