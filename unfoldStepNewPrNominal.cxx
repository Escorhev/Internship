#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
//#include "inputParams.h"
#endif

void unfold(bool doPrompt = true, bool doPbPb = true, Int_t iterMax = 8, Int_t stepNumber = 1) {
  if (!setSystTag(doPbPb)) return;
  
  Int_t iterMin =1;
  Int_t iterDef = iterMax - 2;
  
  string testInputName = "";
  string trainInputName = "";
  string outputName = "";
  
  testInputName = Form("%s/dataUnf/data_results/meas_%s_data_%s%s%s.root",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt",doCent?"_centBin":doPeri?"_periBin":"",systErr?"_systErrs":"_statErrs");
  trainInputName = Form("%s/dataUnf/unfInput/step%i/response_4D_%s_%s_%diter_%dz%dptBins%dz%dptMeasBins%s.root", unfPath.c_str(), stepNumber, doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", iterMax, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, systTag.c_str());
  outputName = Form("%s/dataUnf/unfOutput/step%i/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin%s.root",unfPath.c_str(), stepNumber,doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", iterMax, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, systTag.c_str());
  
  TFile *f_measured = new TFile(testInputName.c_str());
  f_measured->ls();
  TH2D *hMeasured = (TH2D*)f_measured->Get("h_Meas;1");
       
  TFile *f = new TFile(trainInputName.c_str());
  RooUnfoldResponse *resp = (RooUnfoldResponse*)f->Get("resp;1");
  TH2F *hSmear = (TH2F*)f->Get("fh2RespDimM;1");
  TH2F *hTrue = (TH2F*)f->Get("fh2RespDimT;1");

  TH1D *sTruth1DX = dynamic_cast<TH1D*>(hSmear->ProjectionX("TX", 0, -1, "E"));
  TH1D *sTruth1DY = dynamic_cast<TH1D*>(hSmear->ProjectionY("TY", 0, -1, "E"));

  TCanvas * can4 = new TCanvas("can4","can4",1200,600);
  can4->Divide(2,1);
  can4->cd(1);
  sTruth1DX->Draw("EP");
  can4->cd(2);
  sTruth1DY->Draw("EP");
  can4->SaveAs("/home/llr/cms/henderson/UpsUnfolding/Unfolding/dataUnf/code/debughSmear.png");

  //hSmear->Sumw2();
  //hTrue->Sumw2();

  TH2D *h_trMatrix = (TH2D*)f->Get("h_trMatrix;1");
  
  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovToy;//kNoError;//
  //RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kNoError;
  
  Int_t nIterTmp = iterMax - iterMin+1;
  const Int_t nIter = nIterTmp;
  TH2D *hReco[nIter];
  TH2D *hFolded[nIter];
  TMatrixD covmat[nIter];
  RooUnfoldBayes unfold[nIter];
  //RooUnfoldInvert unfoldInv[nIter];
  TMatrixD transfmat[nIter];
  
  TH1D *hJetPtTrue;
  TH1D *hJetPtMeas;
  TH1D *hJetPtUnf[nIter];
  TH1D *hJetPtFol[nIter];

  TH1D *cTruth1DX = dynamic_cast<TH1D*>(hMeasured->ProjectionX("TX", 0, -1, "E"));
  TH1D *cTruth1DY = dynamic_cast<TH1D*>(hMeasured->ProjectionY("TY", 0, -1, "E"));

  TCanvas * can6 = new TCanvas("can6","can6",1200,600);
  can6->Divide(2,1);
  can6->cd(1);
  cTruth1DX->Draw("EP");
  can6->cd(2);
  cTruth1DY->Draw("EP");
  can6->SaveAs("/home/llr/cms/henderson/UpsUnfolding/Unfolding/dataUnf/code/debughMeasured.png");
  
  for(Int_t iter = iterMin; iter<=iterMax; iter++) {

    cout << "iter = " << iter << endl;
    
    // RooUnfoldBayes    unfold (resp, hMeas, iter);
    unfold[iter-iterMin] = RooUnfoldBayes(resp, hMeasured, iter);
    //if (matrixInv) unfoldInv[iter-iterMin] = RooUnfoldInvert(resp, hMeasured);

    // unfold using response matrix
    
    hReco[iter-iterMin] = (TH2D*)unfold[iter-iterMin].Hreco(errorTreatment);
    hReco[iter-iterMin]->SetName(Form("hReco_Iter%d",iter));
    unfold[iter-iterMin].Print();

    //if (matrixInv) hReco[iter-iterMin] = (TH2D*)unfoldInv[iter-iterMin].Hreco(errorTreatment);

    // unfold[iter-iterMin].PrintTable(cout,hTrue,RooUnfold::kCovariance);
    
    // fold using just unfolded result
    
    hFolded[iter-iterMin] = (TH2D*)resp->ApplyToTruth(hReco[iter-iterMin]);
    //hFolded[iter-iterMin]->Sumw2();
    hFolded[iter-iterMin]->SetName(Form("hFolded_Iter%d",iter));

    // covariance matrix problem here : resize it before assigning new matrix.

    cout << "before covmat[iter-iterMin]" << endl;
    covmat[iter-iterMin].ResizeTo(unfold[iter-iterMin].Ereco(errorTreatment));
    covmat[iter-iterMin] = unfold[iter-iterMin].Ereco(errorTreatment);
    //if (matrixInv) {
    //covmat[iter-iterMin].ResizeTo(unfoldInv[iter-iterMin].Ereco(errorTreatment));
    //covmat[iter-iterMin] = unfoldInv[iter-iterMin].Ereco(errorTreatment);
    //}
    cout << "cols : "<< covmat[iter-iterMin].GetNcols() << endl;
    cout << "raws : "<< covmat[iter-iterMin].GetNrows() << endl;
    cout << "after covmat[iter-iterMin]" << endl;

    transfmat[iter-iterMin].ResizeTo(unfold[iter-iterMin].UnfoldingMatrix());
    transfmat[iter-iterMin] = unfold[iter-iterMin].UnfoldingMatrix();
    
    //if (matrixInv) {
    //transfmat[iter-iterMin].ResizeTo(unfoldInv[iter-iterMin].UnfoldingMatrix());
    //transfmat[iter-iterMin] = unfoldInv[iter-iterMin].UnfoldingMatrix();
    //}
    
    //cout << "covmat[iter-iterMin](0,0) " << covmat[iter-iterMin](0,0) << endl;
    //cout << "covmat[iter-iterMin](1,20) " << covmat[iter-iterMin](1,20) << endl;

    /*
    TMatrixD errmat(ntbins,ntbins);
    for (Int_t i=0; i<ntbins; i++)
      for (Int_t j=0; j<ntbins; j++)
	errmat(i,j)= covmat(i,j)>=0 ? sqrt(covmat(i,j)) : -sqrt(-covmat(i,j));
    PrintMatrix(errmat,"","covariance matrix",10);
    */
    
    // jet pt distributions save here
    
    hJetPtUnf[iter-iterMin] = dynamic_cast<TH1D*>(hReco[iter-iterMin]->ProjectionY(Form("hJetPtUnf_Iter%d",iter)));
    hJetPtFol[iter-iterMin] = dynamic_cast<TH1D*>(hFolded[iter-iterMin]->ProjectionY(Form("hJetPtFol_Iter%d",iter)));

    //hJetPtUnf[iter-iterMin]->SetName(Form("hJetPtUnf_Iter%d",iter));
    //hJetPtFol[iter-iterMin]->SetName(Form("hJetPtFol_Iter%d",iter));

    /*
    hJetPtUnf[iter-iterMin]->Sumw2();
    hJetPtFol[iter-iterMin]->Sumw2();
    */
  }

  // response matrix folded with true distribution (test sample)
  
  TH2D *hPriorFolded = (TH2D*)resp->ApplyToTruth(hTrue);
  hPriorFolded->SetName("hPriorFolded");
  
  
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
    
  // 1d projections of true z and jet pt (mc  sample)
  TH1D *hTP[2];
  hTP[0] = dynamic_cast<TH1D*>(hTrue->ProjectionX("hTP_x"));
  hTP[1] = dynamic_cast<TH1D*>(hTrue->ProjectionY("hTP_y"));

  /*
  hSP[0]->Sumw2();
  hFoldedP[0]->Sumw2();

  hSP[1]->Sumw2();
  hFoldedP[1]->Sumw2();

  hTP[0]->Sumw2();
  hRecoP[0]->Sumw2();

  hTP[1]->Sumw2();
  hRecoP[1]->Sumw2();
  */
  
  hJetPtTrue = (TH1D*)hTP[1]->Clone();
  hJetPtTrue->SetName("hJetPtTrue");
  hJetPtMeas = (TH1D*)hSP[1]->Clone();
  hJetPtMeas->SetName("hJetPtMeas");

  /*
  const Int_t nPtBins = 3;  
  Double_t ptmin[nPtBins] = {15.,25.,35.};
  Double_t ptmax[nPtBins] = {25.,35.,45.};
  */
  
  Double_t *ptmin = new Double_t[nBinJet_gen];
  Double_t *ptmax = new Double_t[nBinJet_gen];

  for (int i=0;i<nBinJet_gen;i++) {
    ptmin[i] = min_jetpt+i*jetPt_gen_binWidth;
    ptmax[i] = min_jetpt+(i+1)*jetPt_gen_binWidth;
  }
    
  TH1D *hMUnf[nBinJet_gen][nIter];
  TH1D *hMFol[nBinJet_gen][nIter];
  TH1D *hMTru[nBinJet_gen];
  TH1D *hMSme[nBinJet_gen];
  TH1D *hMPri[nBinJet_gen];

  for(Int_t i = 0; i<nBinJet_gen; i++) {
    Int_t min = hTrue->GetYaxis()->FindBin(ptmin[i]+0.00001);
    Int_t max = hTrue->GetYaxis()->FindBin(ptmax[i]-0.00001);
    hMTru[i] = dynamic_cast<TH1D*>(hTrue->ProjectionX(Form("hMTru_%d",i),min,max));

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
  hTrue->Write("fh2TrueMC");
  hMeasured->Write("fh2MeasData");
  h_trMatrix->Write("trmatrix");
  
  resp->Write();
  hPriorFolded->Write();
  hJetPtTrue->Write();
  hJetPtMeas->Write();
  
  for(Int_t iter = iterMin; iter<=iterMax; iter++) {
    hReco[iter-iterMin]->Write();
    hFolded[iter-iterMin]->Write();
    covmat[iter-iterMin].Write(Form("covmat%d",iter));
    unfold[iter-iterMin].Write(Form("unfold%d",iter));
    //if (matrixInv) 
    //unfoldInv[iter-iterMin].Write(Form("unfoldInv%d",iter));
    transfmat[iter-iterMin].Write(Form("invtrmat%d",iter));
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
  
  TH1D *tTruth1DX = dynamic_cast<TH1D*>(hTrue->ProjectionX("TX", 0, -1, "E"));
  TH1D *tTruth1DY = dynamic_cast<TH1D*>(hTrue->ProjectionY("TY", 0, -1, "E"));

  TCanvas * can5 = new TCanvas("can5","can5",1200,600);
  can5->Divide(2,1);
  can5->cd(1);
  tTruth1DX->Draw("EP");
  can5->cd(2);
  tTruth1DY->Draw("EP");
  can5->SaveAs("/home/llr/cms/henderson/UpsUnfolding/Unfolding/dataUnf/code/debughTrue.png");

  fout->Write();
  fout->Close();
}

void unfoldStepNewPrNominal(Int_t step = 1) { //, double SF = 1.1){
  //if (!matrixInv)
  //  unfold(true,true,nIter,step);
  //if (step<=(nSIter_pp+3) && centShift==0 && !doCent && !doPeri)
  unfold(true,false,nIter,step);
}
