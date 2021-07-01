#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
//#include "inputParams.h"
#endif

void unfoldDiag(bool doPrompt = true, bool doPbPb = true, Int_t iterMax = 3) {
  if (!setSystTag(doPbPb)) return;

  Int_t iterMin =1;
  Int_t iterDef = iterMax - 2;
  
  string testInputName = "";
  string trainInputName = "";
  string outputName = "";
  
  testInputName = Form("%s/dataUnf/unfOutput/step%i/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin%s.root",unfPath.c_str(), iterMax, doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, systTag.c_str());
  trainInputName = Form("%s/dataUnf/unfInput/diag4DMatrixResponseInv.root",unfPath.c_str());
  //trainInputName = "diag4DMatrixResponseInv.root";
  outputName = Form("%s/dataUnf/unfOutput/step%i/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_Diag%s.root",unfPath.c_str(), iterMax, doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, systTag.c_str());
  
  
  TFile *f_measured = new TFile(testInputName.c_str());
  f_measured->ls();
  //unfolded as measured
  TH2D *hMeasured = (TH2D*)f_measured->Get(Form("hReco_Iter%d;1",nIter));
  TMatrixD *covmat = (TMatrixD*)f_measured->Get(Form("covmat%d;1",nIter));
  
  TFile *f = new TFile(trainInputName.c_str());
  RooUnfoldResponse *resp = (RooUnfoldResponse*)f->Get("resp;1");
  TH2F *hSmear = (TH2F*)f->Get("fh2RespDimM;1");
  TH2F *hTrue = (TH2F*)f->Get("fh2RespDimT;1");
  
  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;//kNoError;//
  
  TH2D *hReco;
  RooUnfoldBayes unfold;
  //RooUnfoldInvert unfoldInv; //for the matrix inversion in pp 
  
  TH2D *hRecoInvert;
  RooUnfoldInvert unfoldInvert; //for the errors
  
  
  unfold = RooUnfoldBayes(resp, hMeasured, 1);
  unfold.SetMeasuredCov(*covmat);

  //unfoldInv = RooUnfoldInvert(resp, hMeasured);
  //unfoldInv.SetMeasuredCov(*covmat);  
  
  hReco = (TH2D*)unfold.Hreco(errorTreatment);
  //if (matrixInv) 
  //hReco = (TH2D*)unfoldInv.Hreco(errorTreatment);
  hReco->SetName(Form("hReco_Iter%d",1));

  unfoldInvert = RooUnfoldInvert(resp, hMeasured);
  unfoldInvert.SetMeasuredCov(*covmat);

  hRecoInvert = (TH2D*)unfoldInvert.Hreco(errorTreatment);
  hRecoInvert->SetName("hReco_Invert");
  
  TFile *fout = new TFile(outputName.c_str(),"RECREATE");
  resp->Write();
  
  hReco->Write();
  unfold.Write(Form("unfold%d",1));

  hRecoInvert->Write();
  unfoldInvert.Write("unfoldInvert");
  
  fout->Write();
  fout->Close();
  cout<<"done with unfoldStepNewPrNominalDiag"<<endl;
}

void unfoldStepNewPrNominalDiag(int step){
  //if (!matrixInv)
  //  unfoldDiag(true,true,step);
  //if (step<=(nSIter_pp+3) && (centShift==0) && !doCent && !doPeri)
  unfoldDiag(true,false,step);
}
