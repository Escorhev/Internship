#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include "inputParams.h"
#endif

void create(bool doPrompt = true, bool doPbPb = true, bool doTrain = false, Int_t stepNumber = 1){
  if (!setCaseTag()) return;
  string inputName = "";
  string outputName = "";
  string partOfOutput = "response";

  inputName = Form("%s/mcUnf/unfInput1D/step%i/unfolding_4D_%s_%s_%s_%diter_%dz%dptBins%dz%dptMeasBins%s.root", unfPath.c_str(), stepNumber, doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", doTrain?"Train":"Test", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, caseTag.c_str());

  outputName = inputName;
  int idxReplace = inputName.find("unfolding_4D");
  cout <<"idxReplace "<<idxReplace<<endl;
  outputName.replace(idxReplace,9,partOfOutput);

  cout <<"input = "<<inputName<<", output = "<<outputName<<endl;
  TFile *f = new TFile(inputName.c_str());
  f->ls();

  //Get response

  string thnSparseName = "";
  //take not normalized matrix, since we need real measured 2D distribution from it
  if (!doTrain) thnSparseName = "hs;1";
  //take normalized tr matrix 
  else if (flatPrior) thnSparseName = "hs_newJetPtNorm;1";
  else if (stepNumber==1) thnSparseName = "hs;1";
  else thnSparseName = "hs_newJetPtNorm;1";

  //if(doTrain && flatPrior) thnSparseName = "hs_newJetPtNorm;1";
  //else if(doTrain && !flatPrior) {
  //if (stepNumber==1)
  //thnSparseName = "hs;1";
  //}
  //else thnSparseName = "hs;1";
  
  THnSparseD *hn = static_cast<THnSparseD*>(f->Get(thnSparseName.c_str()));
  //hn->Sumw2();
  
  TH1F * h_z_gen = (TH1F*)f->Get("h_z_gen;1");
  TH1F * hGenZJetPtCentBin = (TH1F*)f->Get("hGenZJetPtCentBin;1");
  TH1F * hGenZJetPtLowBin = (TH1F*)f->Get("hGenZJetPtLowBin;1");
  TH1F * hGenZJetPtHighBin = (TH1F*)f->Get("hGenZJetPtHighBin;1");
  
  const Int_t ndim = 2;
  Int_t dim[ndim];
  for(Int_t i = 0; i<ndim; i++)
    dim[i] = i;

  Int_t nDim = hn->GetNdimensions();
  cout <<"nDim = " << nDim << endl;
  
  Int_t iTrue = 0;
  Int_t iDet = 1;

  TH1D *fh1Smear = dynamic_cast<TH1D*>(hn->Projection(iDet,"E"));
  TH1D *fh1Prior = dynamic_cast<TH1D*>(hn->Projection(iTrue,"E"));
  TCanvas * can3 = new TCanvas("can3","can3",1200,600);
  fh1Prior->Draw("EP");
  can3->SaveAs("/home/llr/cms/henderson/UpsUnfolding/Unfolding/mcUnf/code/debugPrior.png");


  Int_t nBin[2] = {0+nBinJet_reco*nBinZ_reco, nBinJet_gen*nBinZ_gen};
  Double_t mmin[2] = {0., 0.};
  Double_t mmax[2] = {0.+nBinJet_reco*nBinZ_reco,(Double_t) nBinJet_gen*nBinZ_gen};
  
  TH1D *fh1RespDimM = new TH1D("fh1RespDimM","fh1RespDimM",nBin[0],mmin[0],mmax[0]);
  TH1D *fh1RespDimT = new TH1D("fh1RespDimT","fh1RespDimT",nBin[1],mmin[1],mmax[1]);
  TH1D *fh1Miss = new TH1D("fh1Miss","fh1Miss",nBin[1],mmin[1],mmax[1]);
  
  /*
  //dimensions of measured axis
  TH1D *fh1RespDimM = new TH1D("fh1RespDimM","fh1RespDimM",nBinZ[0],mmin[0],mmax[0],nBinPt[0],ptmin[0],ptmax[0]);
  //dimensions of true axis
  TH1D *fh1RespDimT = new TH1D("fh1RespDimT","fh1RespDimT",nBinZ[1],mmin[1],mmax[1],nBinPt[1],ptmin[1],ptmax[1]);  
  //feed-out of response
  TH1D *fh1Miss = new TH1D("fh1Miss","fh1Miss",nBinZ[1],mmin[1],mmax[1],nBinPt[1],ptmin[1],ptmax[1]);
  cout << "fh2Smear->GetEntries() " << fh2Smear->GetEntries() << endl;
  */

  //fill detector-level distribution
  for(Int_t i = 1; i<=fh1RespDimM->GetNbinsX(); i++) {
    Double_t xlow = fh1RespDimM->GetBinLowEdge(i);
    Double_t xup = xlow + 2*(fh1RespDimM->GetBinCenter(i) - xlow);
    Int_t jlow = fh1Smear->FindBin(xlow+0.001);
    Int_t jup = fh1Smear->FindBin(xup-0.001);
    
    Double_t err = 0.;
    Double_t con = fh1Smear->IntegralAndError(jlow,jup,err);
    //cout<<"ContentDet = "<<con<<endl;
    fh1RespDimM->SetBinContent(i,con/2.5);
    fh1RespDimM->SetBinError(i,err);
  }
  Printf("Created fh1RespDimM");

  //fill particle-level distribution
  for(Int_t i = 1; i<=fh1RespDimT->GetNbinsX(); i++) {
    Double_t xlow = fh1RespDimT->GetBinLowEdge(i);
    Double_t xup = xlow + 2*(fh1RespDimT->GetBinCenter(i) - xlow);;
    Int_t jlow = fh1Prior->FindBin(xlow+0.001);
    Int_t jup = fh1Prior->FindBin(xup-0.001);

    Double_t err = 0.;
    Double_t con = fh1Prior->IntegralAndError(jlow,jup,err);
    //cout<<"ContentPart = "<<con<<endl;
    fh1RespDimT->SetBinContent(i,con/16.);
    fh1RespDimT->SetBinError(i,err);
  }
  Printf("Created fh1RespDimT");

  //response object for RooUnfold
  RooUnfoldResponse *fResponse = new RooUnfoldResponse("resp","resp");
  fResponse->Setup(fh1RespDimM,fh1RespDimT);

  //fill RooUnfoldResponse object
  
  Int_t* coord = new Int_t[nDim];
  Int_t nbin = hn->GetNbins();
  
  for(Int_t bin=0; bin<nbin; bin++) {    
    Double_t w = hn->GetBinContent(bin,coord);
    Double_t cTrue = hn->GetAxis(0)->GetBinCenter(coord[0]);
    Double_t cDet = hn->GetAxis(1)->GetBinCenter(coord[1]);
    
    //cout<<"w = "<<w<<" and cTrue = "<<cTrue<<" and cDet = "<<cDet<<endl;
    fResponse->Fill(cDet, cTrue, w);
  }

  delete [] coord;

  //TMatrixD debugMat = new TMatrixD(480, 12);
  //fResponse->PrintMatrix(debugMat);

  //Write response + 1D histos to file
  TFile *fout = new TFile(outputName.c_str(),"RECREATE");
  hn->Write("fhn");
  fResponse->Write("resp");
  fh1Smear->Write("fh1Smear");
  fh1Prior->Write("fh1Prior");
  fh1RespDimM->Write();
  fh1RespDimT->Write();
  fh1Miss->Write();

  h_z_gen->Write();
  hGenZJetPtCentBin->Write();
  hGenZJetPtLowBin->Write();
  hGenZJetPtHighBin->Write();
    
  fout->Write();
  fout->Close();
}

void createRooUnfoldResponse1D(Int_t step = 1){
  //prompt pp
  if (step<=nSIter_pp && centShift ==0 && !doCent && !doPeri){
    create(true,false,true,step);
    create(true,false,false,step);
  }
}
