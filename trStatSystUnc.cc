#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
//#include "inputParams.h"
#endif

void compute(bool doPrompt = true, bool doPbPb = true){

  if (!setSystTag(doPbPb)) return;
  string filename = "";
  string outputfile = "";
  
  int stepNumber = nSIter;
  if (!doPbPb) stepNumber = nSIter_pp;
  gSystem->mkdir(Form("%s/dataUnf/unfOutput/matrixOper",unfPath.c_str()));
  
  filename = Form("%s/dataUnf/unfOutput/step1/matrixOper/matrixOperation_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin%s.root",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, systTag.c_str());
  outputfile = Form("%s/dataUnf/unfOutput/matrixOper/systUnc_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin%s.root",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, systTag.c_str());

  cout <<"reading file = "<<filename<<endl;  
  TFile *file = new TFile(filename.c_str());

  TH1D *h1_zUnf = (TH1D*)file->Get("nominalZUnf");

  h1_zUnf->Rebin(nBinZ_gen/nBinZ_reco);
  
  TH1D* h1_zUnfToys[100];

  float sums[nBinZ_reco];
  float errs[nBinZ_reco];
  float relErrs[nBinZ_reco];

  for (int i = 0; i<nBinZ_reco; i++) {
    sums[i] = 0.;
    errs[i] = 0.;
    relErrs[i] = 0.;
  }
  
  int nToys = 100;
  
  for(int i = 0; i < nToys; i++){

    h1_zUnfToys[i] = (TH1D*)file->Get(Form("zUnfSmear_toy%i",i+1));

    h1_zUnfToys[i]->Rebin(nBinZ_gen/nBinZ_reco);
    
    float nominalBinVal = 0.;
    float toyBinVal = 0.;
    
    for(int ibin = 0; ibin < h1_zUnf->GetNbinsX(); ibin++){

      nominalBinVal = h1_zUnf->GetBinContent(ibin+1);
      toyBinVal = h1_zUnfToys[i]->GetBinContent(ibin+1);

      sums[ibin] += (nominalBinVal-toyBinVal)*(nominalBinVal-toyBinVal);

      //cout << "i = " << i << " ibin = " << ibin << endl;
      //cout << "(nominalBinVal-toyBinVal)*(nominalBinVal-toyBinVal) = " << (nominalBinVal-toyBinVal)*(nominalBinVal-toyBinVal) << endl;
    }  
  }

  TH1D *h1_zUnf_newSyst = (TH1D*)h1_zUnf->Clone("h1_zUnf_newSyst");
  h1_zUnf_newSyst->Reset();
  
  for(int ibin = 0; ibin < nBinZ_reco; ibin++){

    errs[ibin] = TMath::Sqrt(sums[ibin]*1.0/nToys);
    if(h1_zUnf->GetBinContent(ibin+1)>0) relErrs[ibin] = errs[ibin]/(h1_zUnf->GetBinContent(ibin+1));
    else relErrs[ibin] = 0.;

    cout << "z bin = " << ibin+1 << " , err = " << errs[ibin] << " , rel error (in %) = " << relErrs[ibin]*100 << endl;
    h1_zUnf_newSyst->SetBinContent(ibin+1,h1_zUnf->GetBinContent(ibin+1));
    h1_zUnf_newSyst->SetBinError(ibin+1,errs[ibin]);
    
  }

  
  TCanvas * can0 = new TCanvas("can0","can0",600,600);
  h1_zUnf_newSyst->Draw("EP");
  can0->SaveAs("/home/llr/cms/henderson/UpsUnfolding/Unfolding/dataUnf/code/trMatrixSystCheck.png");
  

  TFile *outfile = new TFile(outputfile.c_str(),"RECREATE");
  h1_zUnf_newSyst->Write("zUnf_trMatrixSyst");
  outfile->Close();
  
  file->Close();
}

// syst unc from MC transfer matrix stat limitation

void trStatSystUnc(){
  //compute(bool doPrompt = true, bool doPbPb = true)
  // compute(true,true);
  //if (!doCent && !doPeri)
  compute(true,false);

  //compute(false,true);
  //compute(false,false);

}
