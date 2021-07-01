#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
//#include "inputParams.h"
#endif

void create(bool doPrompt = true, bool doPbPb = true, Int_t stepNumber = 1) {
  if (!setSystTag(doPbPb)) return;
  
  string inputName = "";
  string outputName = "";
  string partOfOutput = "response";

  //string SF_name = "";
  //if(SF == 1.1) SF_name = "_nominal";
  //if(SF == 1.2) SF_name = "_up";
  //if(SF == 1.0) SF_name= "_down";

  cout << "step # =" << stepNumber << endl;
    
  inputName = Form("%s/dataUnf/unfInput/step%i/unfolding_4D_%s_%s_%diter_%dz%dptBins%dz%dptMeasBins%s.root", unfPath.c_str(), stepNumber, doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, systTag.c_str());

  outputName = inputName;
  int idxReplace = inputName.find("unfolding_4D");
  cout <<"idxReplace "<<idxReplace<<endl;
  outputName.replace(idxReplace,9,partOfOutput);
  cout << "outputName = " << outputName << endl;
  
  TFile *f = new TFile(inputName.c_str());
  f->ls();

  //Get response

  string thnSparseName = "";
  //take normalized tr matrix 
  thnSparseName = "hs_newJetPtNorm;1";
  if (!flatPrior && stepNumber==1) thnSparseName = "hs;1";
  
  THnSparseF *hn = static_cast<THnSparseF*>(f->Get(thnSparseName.c_str()));
  //hn->Sumw2();

  TH1D *cTruth1DX = dynamic_cast<TH1D*>(hn->Projection(0, "E"));
  TH1D *cTruth1DY = dynamic_cast<TH1D*>(hn->Projection(1, "E"));

  TCanvas * can6 = new TCanvas("can6","can6",1200,600);
  can6->Divide(2,1);
  can6->cd(1);
  cTruth1DX->Draw("EP");
  can6->cd(2);
  cTruth1DY->Draw("EP");
  can6->SaveAs("/home/llr/cms/henderson/UpsUnfolding/Unfolding/dataUnf/code/debughn.png");

  
  TH1F * h_z_gen = (TH1F*)f->Get("h_z_gen;1");
  TH1F * hGenZJetPtCentBin = (TH1F*)f->Get("hGenZJetPtCentBin;1");
  TH1F * hGenZJetPtLowBin = (TH1F*)f->Get("hGenZJetPtLowBin;1");
  TH1F * hGenZJetPtHighBin = (TH1F*)f->Get("hGenZJetPtHighBin;1");

  TH2D * h_trMatrix = (TH2D*)f->Get("h_trMatrix;1");
  
  const Int_t ndim = 4;
  Int_t dim[ndim];
  for(Int_t i = 0; i<ndim; i++)
    dim[i] = i;

  Int_t nDim = hn->GetNdimensions();
  cout <<"nDim = " << nDim << endl;
  
  Int_t iPtTrue = 0;
  Int_t iZTrue = 1;
  Int_t iPtDet = 2;
  Int_t iZDet = 3;

  TH2D *fh2Smear = dynamic_cast<TH2D*>(hn->Projection(iPtDet,iZDet,"E"));
  TH2D *fh2Prior = dynamic_cast<TH2D*>(hn->Projection(iPtTrue,iZTrue,"E"));

  Int_t nBinPt[2] = {nBinJet_reco,nBinJet_gen};
  Double_t ptmin[2] = {min_jetpt,min_jetpt};
  Double_t ptmax[2] = {max_jetpt,max_jetpt};

  Int_t nBinZ[2] = {nBinZ_reco,nBinZ_gen};
  Double_t mmin[2] = {min_z,min_z};
  Double_t mmax[2] = {max_z,max_z};
    
  //dimensions of measured axis
  TH2D *fh2RespDimM = new TH2D("fh2RespDimM","fh2RespDimM",nBinZ[0],mmin[0],mmax[0],nBinPt[0],ptmin[0],ptmax[0]);
  //dimensions of true axis
  TH2D *fh2RespDimT = new TH2D("fh2RespDimT","fh2RespDimT",nBinZ[1],mmin[1],mmax[1],nBinPt[1],ptmin[1],ptmax[1]);  
  //feed-out of response
  TH2D *fh2Miss     = new TH2D("fh2Miss","fh2Miss",nBinZ[1],mmin[1],mmax[1],nBinPt[1],ptmin[1],ptmax[1]);
  cout << "fh2Smear->GetEntries() " << fh2Smear->GetEntries() << endl;

  //fill detector-level distribution
  for(Int_t ix = 1; ix<=fh2RespDimM->GetNbinsX(); ix++) {
    Double_t xlow = fh2RespDimM->GetXaxis()->GetBinLowEdge(ix);
    Double_t xup = fh2RespDimM->GetXaxis()->GetBinUpEdge(ix);
    Int_t jxlow = fh2Smear->GetXaxis()->FindBin(xlow+0.000001);
    Int_t jxup = fh2Smear->GetXaxis()->FindBin(xup-0.000001);
    for(Int_t iy = 1; iy<=fh2RespDimM->GetNbinsY(); iy++) {
      Double_t ylow = fh2RespDimM->GetYaxis()->GetBinLowEdge(iy);
      Double_t yup = fh2RespDimM->GetYaxis()->GetBinUpEdge(iy);
      Int_t jylow = fh2Smear->GetYaxis()->FindBin(ylow+0.000001);
      Int_t jyup = fh2Smear->GetYaxis()->FindBin(yup-0.000001);

      Double_t err = 0.;
      Double_t con = fh2Smear->IntegralAndError(jxlow,jxup,jylow,jyup,err);
      //cout << "ix = " << ix << " , iy = " << iy << " jxlow = " << jxlow << " jxup = " << jxup << " jylow = " << jylow << " jyup = " << jyup << endl;
      //cout << "con = " << con << endl;
      fh2RespDimM->SetBinContent(ix,iy,con);
      fh2RespDimM->SetBinError(ix,iy,err);
    }
  }
  Printf("Created fh2RespDimM");

  //fill particle-level distribution
  for(Int_t ix = 1; ix<=fh2RespDimT->GetNbinsX(); ix++) {
    Double_t xlow = fh2RespDimT->GetXaxis()->GetBinLowEdge(ix);
    Double_t xup = fh2RespDimT->GetXaxis()->GetBinUpEdge(ix);
    Int_t jxlow = fh2Prior->GetXaxis()->FindBin(xlow+0.000001);
    Int_t jxup = fh2Prior->GetXaxis()->FindBin(xup-0.000001);
    for(Int_t iy = 1; iy<=fh2RespDimT->GetNbinsY(); iy++) {
      Double_t ylow = fh2RespDimT->GetYaxis()->GetBinLowEdge(iy);
      Double_t yup = fh2RespDimT->GetYaxis()->GetBinUpEdge(iy);
      Int_t jylow = fh2Prior->GetYaxis()->FindBin(ylow+0.000001);
      Int_t jyup = fh2Prior->GetYaxis()->FindBin(yup-0.000001);

      Double_t err = 0.;
      Double_t con = fh2Prior->IntegralAndError(jxlow,jxup,jylow,jyup,err);
      //cout << "ix = " << ix << " , iy = " << iy << endl;
      //cout << "con = " << con << endl;
      fh2RespDimT->SetBinContent(ix,iy,con);
      fh2RespDimT->SetBinError(ix,iy,err);
    }
  }
  Printf("Created fh2RespDimT");

  //response object for RooUnfold
  RooUnfoldResponse *fResponse = new RooUnfoldResponse("resp","resp");
  fResponse->Setup(fh2RespDimM,fh2RespDimT);

  //Fill RooUnfoldResponse object
  Int_t* coord = new Int_t[nDim];
  Int_t nbin = hn->GetNbins();
  
  for(Int_t bin=0; bin<nbin; bin++) {    
    Double_t w = hn->GetBinContent(bin,coord);
    Double_t pttrue = hn->GetAxis(0)->GetBinCenter(coord[0]);
    Double_t ztrue = hn->GetAxis(1)->GetBinCenter(coord[1]);
    Double_t ptdet = hn->GetAxis(2)->GetBinCenter(coord[2]);
    Double_t zdet = hn->GetAxis(3)->GetBinCenter(coord[3]);    
    if(zdet>=mmin[0] && zdet<=mmax[0]
       && ztrue>=mmin[1] && ztrue<=mmax[1]
       && ptdet>=ptmin[0] && ptdet<=ptmax[0]
       && pttrue>=ptmin[1] && pttrue<=ptmax[1]
       ){
         fResponse->Fill(zdet,ptdet,ztrue,pttrue,w);
    } 
    else {      
      fResponse->Miss(ztrue,pttrue,w);
      fh2Miss->Fill(ztrue,pttrue,w);
    }
  }

  delete [] coord;
  //Write response + 2D histos to file
  TFile *fout = new TFile(outputName.c_str(),"RECREATE");
  hn->Write("fhn");
  fResponse->Write("resp");
  fh2Smear->Write("fh2Smear");
  fh2Prior->Write("fh2Prior");
  fh2RespDimM->Write();
  fh2RespDimT->Write();
  fh2Miss->Write();

  h_trMatrix->Write();
  h_z_gen->Write();
  hGenZJetPtCentBin->Write();
  hGenZJetPtLowBin->Write();
  hGenZJetPtHighBin->Write();
    
  fout->Write();
  fout->Close();
}

void createRooUnfoldResponseNewPrNominal(Int_t step = 1) { //, double SF = 1.1){
  //if (!matrixInv)
  //  create(true,true,step);
  //if (step<=(nSIter_pp+3) && centShift==0 && !doCent && !doPeri)
  create(true,false,step);
}

  
