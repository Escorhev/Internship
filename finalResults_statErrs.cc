#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
//#include "inputParams.h"
#endif

void plotRes(bool doPrompt = true, bool doPbPb = false){
  if (!setSystTag(doPbPb)) return;
  gSystem->mkdir(Form("%s/dataUnf/unfOutput/finalResults/",unfPath.c_str()));

  cout<<"test"<<endl;

  string filename1 = "";
  string filename2 = "";
  string filename3 = "";
  string filename4 = "";
  string filenameMeas = "";
  string outputfile = "";

  int iterFinal = nSIter;
  if (!doPbPb) iterFinal = nSIter_pp;
  if (iterFinal<4) iterFinal = 4;

  filename1 = Form("%s/dataUnf/unfOutput/step%i/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_Diag%s.root",unfPath.c_str(),1,doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, systTag.c_str());
  filename2 = Form("%s/dataUnf/unfOutput/step%i/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_Diag%s.root",unfPath.c_str(),2,doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, systTag.c_str());
  filename3 = Form("%s/dataUnf/unfOutput/step%i/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_Diag%s.root",unfPath.c_str(),iterFinal-1,doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, systTag.c_str());
  filename4 = Form("%s/dataUnf/unfOutput/step%i/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_Diag%s.root",unfPath.c_str(),iterFinal,doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, systTag.c_str());
  filenameMeas = Form("%s/dataUnf/unfOutput/step%i/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin%s.root",unfPath.c_str(),iterFinal,doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, systTag.c_str());
  outputfile = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%s%s.root",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, doPbPb?nSIter:nSIter_pp, systTag.c_str(), systErr?"":"_statError");

  TFile *file1 = new TFile(filename1.c_str());
  TFile *file2 = new TFile(filename2.c_str());
  TFile *file3 = new TFile(filename3.c_str());
  TFile *file4 = new TFile(filename4.c_str());

  TFile *fileMeas = new TFile(filenameMeas.c_str());
  
  TH2D *h2ZMeas;
  h2ZMeas=(TH2D*)fileMeas->Get("fh2MeasData;1");
  TH1D *hZMeas = h2ZMeas->ProjectionX("hZMeas",3,3);
  hZMeas->Draw();

  cout <<"getting the hists"<<endl;  
  TH2D *h2UnfResp1 = (TH2D*)file1->Get(Form("hReco_Iter%d;1",1));
  TH2D *h2UnfResp2 = (TH2D*)file2->Get(Form("hReco_Iter%d;1",1));
  TH2D *h2UnfResp3 = (TH2D*)file3->Get(Form("hReco_Iter%d;1",1));
  TH2D *h2UnfResp4 = (TH2D*)file4->Get(Form("hReco_Iter%d;1",1));

  cout <<"hists found"<<endl;  

  Int_t min_PriorFolded = h2UnfResp1->GetYaxis()->FindBin(midLowerPt+0.00001);
  Int_t max_PriorFolded = h2UnfResp1->GetYaxis()->FindBin(midUpperPt-0.00001);

  cout <<"min_PriorFolded = "<<min_PriorFolded<<"; max_PriorFolded = "<<max_PriorFolded<<endl;
  
  TH1D *hZUnf_SI1 = (TH1D*)h2UnfResp1->ProjectionX("hZUnf_SI1",min_PriorFolded,max_PriorFolded);
  TH1D *hZUnf_SI2 = (TH1D*)h2UnfResp2->ProjectionX("hZUnf_SI2",min_PriorFolded,max_PriorFolded);
  TH1D *hZUnf_SI3 = (TH1D*)h2UnfResp3->ProjectionX("hZUnf_SI3",min_PriorFolded,max_PriorFolded);
  TH1D *hZUnf_SI4 = (TH1D*)h2UnfResp4->ProjectionX("hZUnf_SI4",min_PriorFolded,max_PriorFolded);

  cout<<"saving"<<endl;
  TFile *file_out = new TFile(outputfile.c_str(),"RECREATE");
  hZUnf_SI1->Write(Form("zUnfSI%d",1));
  hZUnf_SI2->Write(Form("zUnfSI%d",2));
  hZUnf_SI3->Write(Form("zUnfSI%d",iterFinal-1));
  hZUnf_SI4->Write(Form("zUnfSI%d",iterFinal));
  hZMeas->Write("zMeas");
  file_out->Close();
}

void finalResults_statErrs(){
  //if (!matrixInv)
  //  plotRes(true,true);
  //if (centShift ==0 && !doCent && !doPeri)
  plotRes(true,false);
}
