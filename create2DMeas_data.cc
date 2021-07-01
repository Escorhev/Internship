#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
//#include "inputParams.h"
#endif

void create2DMeas_data(bool doPrompt = true, bool doPbPb = false){
  bool statErr = true;
  //if ((doCent&&doPeri) ||!doPbPb) {doCent = false; doPeri = false;}

  gSystem->mkdir(Form("%s/dataUnf/plots",unfPath.c_str()));
  gSystem->mkdir(Form("%s/dataUnf/data_results",unfPath.c_str()));
    
  string filename0 = "";
  string filename1 = "";
  string filename2 = "";
  
  filename0 = "/data_CMS/cms/henderson/Unfolding/dataUnf/data_results/fitResults_10_20.root";
  filename1 = "/data_CMS/cms/henderson/Unfolding/dataUnf/data_results/fitResults_20_30.root";
  filename2 = "/data_CMS/cms/henderson/Unfolding/dataUnf/data_results/fitResults_30_40.root";
  
  TFile *file0 = new TFile(filename0.c_str(),"READ");
  TFile *file1 = new TFile(filename1.c_str(),"READ");
  TFile *file2 = new TFile(filename2.c_str(),"READ");

  TH1F *h_10_20_Bin;
  TH1F *h_20_30_Bin;
  TH1F *h_30_40_Bin;

  h_10_20_Bin = (TH1F*)file0->Get("signal_nSigY1S;1");
  h_20_30_Bin = (TH1F*)file1->Get("signal_nSigY1S;1");
  h_30_40_Bin = (TH1F*)file2->Get("signal_nSigY1S;1");
  
  TH2D *h_Meas = new TH2D("h_Meas","h_Meas",nBinZ_reco,min_z,max_z,nBinJet_reco,min_jetpt,max_jetpt);
  h_Meas->Sumw2();
  
  int nbins = h_20_30_Bin->GetNbinsX();

  float _10_20_BinVals[nbins];
  float _20_30_BinVals[nbins];
  float _30_40_BinVals[nbins];

  float _10_20_BinErrs[nbins];
  float _20_30_BinErrs[nbins];
  float _30_40_BinErrs[nbins];
  
  double totIntegral = 0;
  double totErr = 0;
  for(int i = 1; i < nbins+1; i++){
    cout << Form("h_10_20_Bin : %i ", i) << h_10_20_Bin->GetBinContent(i) << endl;
    _10_20_BinVals[i-1] = h_10_20_Bin->GetBinContent(i);
    _10_20_BinErrs[i-1] = h_10_20_Bin->GetBinError(i);
    totIntegral=totIntegral+_10_20_BinVals[i-1];
    totErr=totErr+_10_20_BinErrs[i-1]*_10_20_BinErrs[i-1];
  }
  
  for(int i = 1; i < nbins+1; i++){
    cout << Form("h_20_30_Bin : %i ", i) << h_20_30_Bin->GetBinContent(i) << endl;
    _20_30_BinVals[i-1] = h_20_30_Bin->GetBinContent(i);
    _20_30_BinErrs[i-1] = h_20_30_Bin->GetBinError(i);
    totIntegral=totIntegral+_20_30_BinVals[i-1];
    totErr=totErr+_20_30_BinErrs[i-1]*_20_30_BinErrs[i-1];
  }

  for(int i = 1; i < nbins+1; i++){
    cout << Form("h_30_40_Bin : %i ", i) << h_30_40_Bin->GetBinContent(i) << endl;
    _30_40_BinVals[i-1] = h_30_40_Bin->GetBinContent(i);
    _30_40_BinErrs[i-1] = h_30_40_Bin->GetBinError(i);
    totIntegral=totIntegral+_30_40_BinVals[i-1];
    totErr=totErr+_30_40_BinErrs[i-1]*_30_40_BinErrs[i-1];
  }

  cout << "***" << endl;

  totErr = sqrt(totErr);

  cout << "***"<< endl;

  float content = 0.0;
  float errVal = 0.0;
  
  for(int icount=0; icount < nBinJet_reco; icount++){
    for(int jcount=0; jcount < nBinZ_reco; jcount++){

      if(icount == 0) {
	content = _10_20_BinVals[jcount];
	errVal = _10_20_BinErrs[jcount];
      }
      if(icount == 1) {
	content = _20_30_BinVals[jcount];
	errVal = _20_30_BinErrs[jcount];
      }
      if(icount == 2) {
	content = _30_40_BinVals[jcount];
	errVal = _30_40_BinErrs[jcount];
      }
      h_Meas->SetBinContent(jcount+1,icount+1,content);
      h_Meas->SetBinError(jcount+1,icount+1,errVal);
    }
  }

  TCanvas * can0 = new TCanvas ("can0","can0",1500,500);
  can0->Divide(3,1);

  can0->cd(1);
  h_10_20_Bin->SetTitle(Form("%d < jet p_{T} < %d",(int) (min_jetpt+0*jetPt_reco_binWidth), (int)(min_jetpt+1*jetPt_reco_binWidth)));
  h_10_20_Bin->Draw();

  can0->cd(2);
  h_20_30_Bin->SetTitle(Form("%d < jet p_{T} < %d",(int) (min_jetpt+1*jetPt_reco_binWidth), (int)(min_jetpt+2*jetPt_reco_binWidth)));
  h_20_30_Bin->Draw();

  can0->cd(3);
  h_30_40_Bin->SetTitle(Form("%d < jet p_{T} < %d",(int) (min_jetpt+2*jetPt_reco_binWidth), (int)(min_jetpt+3*jetPt_reco_binWidth)));
  h_30_40_Bin->Draw();

  can0->SaveAs(Form("%s/dataUnf/plots/data_%s_z_distr_%s_%sErr.pdf",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt",statErr?"stat":"syst"));

  TCanvas * can1 = new TCanvas ("can1","can1",600,600);
  can1->SetRightMargin(0.2);
  
  h_Meas->SetTitle(Form("%s %s",doPrompt?"prompt":"nonprompt",doPbPb?"PbPb":"PP"));
  
  h_Meas->SetStats(0);
  h_Meas->GetXaxis()->SetTitle("z");
  h_Meas->GetYaxis()->SetTitle("jet p_{T}");
  h_Meas->GetYaxis()->SetTitleOffset(1.2);
  h_Meas->Draw("TEXTcolz");
  can1->SaveAs(Form("%s/dataUnf/plots/data_%s_meas_%s_%sErr.pdf",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt",statErr?"stat":"syst"));

  
  string outputfile = "";
  outputfile = Form("%s/dataUnf/data_results/meas_%s_data_%s%s_%sErrs.root",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt",doCent?"_centBin":doPeri?"_periBin":"",statErr?"stat":"syst");
  
  TFile *file_data_meas = new TFile(outputfile.c_str(),"RECREATE");

  h_Meas->Write();
  h_10_20_Bin->Write("zMeas_10_20_Bin");
  h_20_30_Bin->Write("zMeas_20_30_Bin");
  h_30_40_Bin->Write("zMeas_30_40_Bin");

  cout <<"for "<<(doPbPb?"PbPb":"PP")<<" ,total number of events = "<<totIntegral<<", err = "<<totErr<<endl;
  if (!doPbPb) {
    totIntegral = totIntegral;
    totErr = totErr;
  }
  else if(!doCent && !doPeri) {
    totIntegral = totIntegral;
    totErr = totErr;
  }
  cout <<"for "<<(doPbPb?"PbPb":"PP")<<" ,Normalized total number of events = "<<totIntegral<<", err = "<<totErr<<endl;
}
