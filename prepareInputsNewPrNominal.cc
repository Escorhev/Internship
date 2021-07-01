#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
//#include "inputParams.h"
#endif

TF1* promptJESFunc() {
  TF1 *zFunc = new TF1("pol4","pol4",0.15,1.);
  zFunc->SetParameters(1.09912e+00,-7.43642e-01,1.49694e+00,-1.35510e+00,5.20681e-01);
  return zFunc;
}

TF1* nonpromptJESFunc(){
  TF1 *zFunc = new TF1("pol5","pol5",0.15,1.);
  zFunc->SetParameters(8.49432e-01,8.41236e-01,-1.15481e+00,-6.29518e-01,2.06411e+00,-9.58378e-01);
  return zFunc;
}


float getZWeight(float z_gen = 0){
  TF1 *zFunc = new TF1("pol6","pol6",0.18,1.);
  //if(doMid) 
  zFunc->SetParameters(11.2639,-252.785,1698.17,-4651.78,6166.36,-3973.96,1002.8);
  //else zFunc->SetParameters(53.2622,-763.157,4007.63,-9781.24,12194.8,-7581.59,1870.36);
  float wZ = zFunc->Eval(z_gen);
  return wZ;
}

void preparePrior(bool doPbPb = true, string outputName="");
void prepare(bool doPrompt = false, bool doPbPb = true, Int_t stepNumber = 1) { 
  printInput();
  if (!setSystTag(doPbPb)) return;
  setSFVal();
  
  string filename = "";
  string outputfile = "";
  string filenamePrevStep = "";
  
  gSystem->mkdir(Form("%s/dataUnf/unfInput",unfPath.c_str()));
  gSystem->mkdir(Form("%s/dataUnf/unfInput/step%i",unfPath.c_str(),stepNumber));
  gSystem->mkdir(Form("%s/dataUnf/unfOutput",unfPath.c_str()));
  gSystem->mkdir(Form("%s/dataUnf/unfOutput/step%i",unfPath.c_str(),stepNumber));

  //filename = Form("~/JpsiInJetsPbPb/Fitter/TreesForUnfolding/tree_%s_%s_NoBkg%s_AccEff_JEC.root",useSystTrM?(doPrompt?"MCJPSINOPR":"MCJPSIPR"):(doPrompt?"MCJPSIPR":"MCJPSINOPR"),doPbPb?"PbPb":"PP",Form("_jetR%d",(int)(jetR*10)));
  filename = "/data_CMS/cms/henderson/merged_HiForestAOD_MC_skimjet.root";
  outputfile = Form("%s/dataUnf/unfInput/step%i/unfolding_4D_%s_%s_%diter_%dz%dptBins%dz%dptMeasBins%s.root",unfPath.c_str(),stepNumber,doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, systTag.c_str());
  if(stepNumber > 1) filenamePrevStep = Form("%s/dataUnf/unfOutput/step%i/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin%s.root",unfPath.c_str(),stepNumber-1,doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, systTag.c_str());
  if (nprPrior && stepNumber == 1) {
    filenamePrevStep = Form("%s/dataUnf/unfInput/step%i/priorDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin%s.root",unfPath.c_str(),stepNumber,doPbPb?"PbPb":"PP","nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, systTag.c_str());
    preparePrior(doPbPb,filenamePrevStep);
  }

  cout <<"filename = "<<filename<<endl;
  cout <<"outputfile = "<<outputfile<<endl;
  cout <<"filenamePrevStep = "<<filenamePrevStep<<endl;
  TFile *file = new TFile(filename.c_str());
  TTree *t_unf = (TTree*)file->Get("onia_skimmed");

  int n_entries = t_unf->GetEntries();
  //n_entries=10;  
  Double_t z_frac_JetPtBin[nBinJet_gen][nBinZ_gen];
  
  //only for prompt fwd in steps > 1, until the stats are better in mc
  Double_t *z_frac_allJetPtBins = new Double_t[nBinZ_gen];

  // if this is step 2, 3, 4 ... take the z unfolded and use it for the prior, else put z prior to flat
  
  if(stepNumber > 1) {
    TFile *filePrevStep = new TFile(filenamePrevStep.c_str());
    //    cout << "file OK" << endl;
    TH1D *h_z_unf_jetPtBin[nBinJet_gen];
    for(int ibin = 0; ibin < nBinJet_gen; ibin++){
      h_z_unf_jetPtBin[ibin] = (TH1D*)filePrevStep->Get(Form("hMUnf_%i_Iter%d;1",ibin,nIter));
    }
    
    TH1D *h_z_allBins = (TH1D*)h_z_unf_jetPtBin[0]->Clone();
    for(int ibin = 1; ibin < nBinJet_gen; ibin++){
      h_z_allBins->Add(h_z_unf_jetPtBin[ibin]);
    }
    
    h_z_allBins->Scale(1/h_z_allBins->Integral());

    for(int ibin = 1; ibin < nBinJet_gen; ibin++){
      h_z_unf_jetPtBin[ibin]->Scale(1/h_z_unf_jetPtBin[ibin]->Integral());
    }
            
    for(int iz = 1; iz <= nBinZ_gen; iz++){
      z_frac_allJetPtBins[iz-1] = h_z_allBins->GetBinContent(iz);
      //cout << "iz = " << iz << " allJtPtBins content = " << z_frac_allJetPtBins[iz-1] << endl;  
      for(int ijtpt = 0; ijtpt < nBinJet_gen; ijtpt++){
	z_frac_JetPtBin[ijtpt][iz-1] = h_z_unf_jetPtBin[ijtpt]->GetBinContent(iz);
      }
    }    
  }

  else if(nprPrior) {
    TFile *filePrevStep = new TFile(filenamePrevStep.c_str());
    //    cout << "file OK" << endl;
    TH1D *h_z_unf_jetPtBin[nBinJet_gen];
    for(int ibin = 0; ibin < nBinJet_gen; ibin++){
      h_z_unf_jetPtBin[ibin] = (TH1D*)filePrevStep->Get(Form("zDist_%i",ibin));
    }
    
    TH1D *h_z_allBins = (TH1D*)h_z_unf_jetPtBin[0]->Clone();
    for(int ibin = 1; ibin < nBinJet_gen; ibin++){
      h_z_allBins->Add(h_z_unf_jetPtBin[ibin]);
    }
    
    h_z_allBins->Scale(1/h_z_allBins->Integral());

    for(int ibin = 1; ibin < nBinJet_gen; ibin++){
      h_z_unf_jetPtBin[ibin]->Scale(1/h_z_unf_jetPtBin[ibin]->Integral());
    }
            
    for(int iz = 1; iz <= nBinZ_gen; iz++){
      z_frac_allJetPtBins[iz-1] = h_z_allBins->GetBinContent(iz);
      //cout << "iz = " << iz << " allJtPtBins content = " << z_frac_allJetPtBins[iz-1] << endl;  
      for(int ijtpt = 0; ijtpt < nBinJet_gen; ijtpt++){
	z_frac_JetPtBin[ijtpt][iz-1] = h_z_unf_jetPtBin[ijtpt]->GetBinContent(iz);
      }
    }
  }

  else {
    for(int i = 1; i <= nBinZ_gen; i++){
      for(int ijtpt = 0; ijtpt < nBinJet_gen; ijtpt++){
	z_frac_JetPtBin[ijtpt][i-1] = 1.0;
      }
    }
  }
  
  //J/Psi variables:
  float reco_pT = 0.; t_unf->SetBranchAddress("reco_pT", &reco_pT);
  float reco_eta = 0.; t_unf->SetBranchAddress("reco_eta", &reco_eta);
  float reco_mass = 0.; t_unf->SetBranchAddress("reco_mass", &reco_mass);

  //jet reco variables:
  float jt_pt = 0.; t_unf->SetBranchAddress("jt_pt", &jt_pt);
  //float jt_pt_JEU_Up = 0.; t_unf->SetBranchAddress("jt_pt_JEU_Up", &jt_pt_JEU_Up);
  //float jt_pt_JEU_Down = 0.; t_unf->SetBranchAddress("jt_pt_JEU_Down", &jt_pt_JEU_Down);
  float jt_eta = 0.; t_unf->SetBranchAddress("jt_eta", &jt_eta);

  //centrality:
  //int centr = 0; t_unf->SetBranchAddress("centr", &centr);

  //correction variables:
  //float corr_AccEff = 0.; t_unf->SetBranchAddress("corr_AccEff", &corr_AccEff);
  float accxeff = 0.; t_unf->SetBranchAddress("accxeff", &accxeff);
  float pthat = 0.; t_unf->SetBranchAddress("pthat",&pthat);
  float finalCorr = 0.;

  //z corrected
  float z = 0.; t_unf->SetBranchAddress("z", &z);

  //jet gen variables
  float jt_ref_pt = 0.; t_unf->SetBranchAddress("jt_ref_pt", &jt_ref_pt);

  //gen z:
  float gen_z = 0.; t_unf->SetBranchAddress("gen_z", &gen_z);

  t_unf->SetBranchStatus("*",0);
  t_unf->SetBranchStatus("reco_pT",1);
  t_unf->SetBranchStatus("reco_eta",1);
  t_unf->SetBranchStatus("reco_mass",1);
  t_unf->SetBranchStatus("jt_pt",1);
  //t_unf->SetBranchStatus("jt_pt_JEU_Up",1);
  //t_unf->SetBranchStatus("jt_pt_JEU_Down",1);  
  t_unf->SetBranchStatus("jt_pt",1);
  t_unf->SetBranchStatus("jt_eta",1);
  //t_unf->SetBranchStatus("centr",1);
  t_unf->SetBranchStatus("accxeff",1);
  t_unf->SetBranchStatus("pthat",1);
  t_unf->SetBranchStatus("z",1);
  t_unf->SetBranchStatus("jt_ref_pt",1);
  t_unf->SetBranchStatus("gen_z",1);
  
  // histograms for renormalizing response matrix:
  TH2F * hRespZJetPtCentBin = new TH2F ("hRespZJetPtCentBin","z_gen vs z_raw for 30 < refpt < 40",nBinZ_gen,min_z,max_z,nBinZ_gen,min_z,max_z);
  TH2F * hRespZJetPtLowBin = new TH2F ("hRespZJetPtLowBin","z_gen vs z_raw for 20 < refpt < 30",nBinZ_gen,min_z,max_z,nBinZ_gen,min_z,max_z);
  TH2F * hRespZJetPtHighBin = new TH2F ("hRespZJetPtHighBin","z_gen vs z_raw for 40 < refpt < 50",nBinZ_gen,min_z,max_z,nBinZ_gen,min_z,max_z);

  TH1D * h_z_gen = new TH1D ("h_z_gen","z_gen for 20 < refpt < 50",nBinZ_gen,min_z,max_z);
  TH1D * hGenZJetPtCentBin = new TH1D ("hGenZJetPtCentBin","z_gen for 30 < refpt < 40",nBinZ_gen,min_z,max_z);
  TH1D * hGenZJetPtLowBin = new TH1D ("hGenZJetPtLowBin","z_gen for 20 < refpt < 30",nBinZ_gen,min_z,max_z);
  TH1D * hGenZJetPtHighBin = new TH1D ("hGenZJetPtHighBin","z_gen for 40 < refpt < 50",nBinZ_gen,min_z,max_z);
    
  TH2F * h_jetpPtGen_jetPtReco = new TH2F ("h_jetpPtGen_jetPtReco","gen vs reco corr jet pt",nBinJet_reco,min_jetpt,max_jetpt,nBinJet_reco,min_jetpt,max_jetpt);

  TF1* nonpromptJES = nonpromptJESFunc();
  TF1* promptJES = promptJESFunc();
  //cout <<"nonpromptJES function parameters: [0] =  "<<nonpromptJES->GetParameter(0)<<", [1] = "<<nonpromptJES->GetParameter(1)<<", [2] = "<<nonpromptJES->GetParameter(2)<<", [3] = "<<nonpromptJES->GetParameter(3)<<", [4] = "<<nonpromptJES->GetParameter(4)<<", [5] = "<<nonpromptJES->GetParameter(5)<<endl;
  //cout <<"promptJES function parameters: [0] =  "<<promptJES->GetParameter(0)<<", [1] = "<<promptJES->GetParameter(1)<<", [2] = "<<promptJES->GetParameter(2)<<", [3] = "<<promptJES->GetParameter(3)<<", [4] = "<<promptJES->GetParameter(4)<<endl;
  hRespZJetPtCentBin->Sumw2();
  hRespZJetPtLowBin->Sumw2();
  hRespZJetPtHighBin->Sumw2();
  h_jetpPtGen_jetPtReco->Sumw2();

  // 4D histogram creation
  
  int fDim = 4.;
  Double_t fValue[fDim];
  
  Int_t* bins_sparce = new Int_t[fDim];
  Double_t *xmin_sparce = new Double_t[fDim];
  Double_t *xmax_sparce = new Double_t[fDim];

  bins_sparce[0] = nBinJet_gen;
  xmin_sparce[0] = min_jetpt;
  xmax_sparce[0] = max_jetpt;

  bins_sparce[1] = nBinZ_gen;
  xmin_sparce[1] = min_z;
  xmax_sparce[1] = max_z;

  bins_sparce[2] = nBinJet_reco;
  xmin_sparce[2] = min_jetpt;
  xmax_sparce[2] = max_jetpt;

  bins_sparce[3] = nBinZ_reco;
  xmin_sparce[3] = min_z;
  xmax_sparce[3] = max_z;  
  
  //initial not normalized 4D  
  THnSparseF * fSparse = new THnSparseF("hs", "hs", fDim, bins_sparce, xmin_sparce, xmax_sparce);
  fSparse->Sumw2();
  fSparse->CalculateErrors();

  int nBins_X = nBinJet_reco*nBinZ_reco;
  double min_X = 0;
  double max_X = nBinJet_reco*nBinZ_reco;

  int nBins_Y = nBinJet_gen*nBinZ_gen;
  double min_Y = 0;
  double max_Y =  nBinJet_gen*nBinZ_gen;
  
  TH2D *h_trMatrix = new TH2D("h_trMatrix","h_trMatrix",nBins_X,min_X,max_X,nBins_Y,min_Y,max_Y);
  h_trMatrix->Sumw2();
    
  //4D with z gen flat
  THnSparseF * fSparse_newZNorm = new THnSparseF("hs_newZNorm", "hs_newZNorm", fDim, bins_sparce, xmin_sparce, xmax_sparce);
  fSparse_newZNorm->Sumw2();
  fSparse_newZNorm->CalculateErrors();

  //4D with z gen flat and gen pt flat
  THnSparseF * fSparse_newJetPtNorm = new THnSparseF("hs_newJetPtNorm", "hs_newJetPtNorm", fDim, bins_sparce, xmin_sparce, xmax_sparce);
  fSparse_newJetPtNorm->Sumw2();
  fSparse_newJetPtNorm->CalculateErrors();

  float jp_pt_cut = 0.;

  TRandom *rand = new TRandom();

  double c2 = (5.92454e-02)*(5.92454e-02);//(2.54002e-02)*(2.54002e-02);//0.06*0.06;
  double s2 = (9.84442e-01)*(9.84442e-01);//(9.42256e-01)*(9.42256e-01);//0.8*0.8;

  int sfBin = -1;
  if(SF<1.05) sfBin=0;
  else if(SF<1.15) sfBin=1;
  else sfBin=2;
  
  for(int nEv = 0; nEv < n_entries; nEv++) {
    if (nEv%100000==0) cout<<"processing evt "<<nEv<<"/"<<n_entries<<endl;
    t_unf->GetEntry(nEv);

    // check event content
    //if (doPbPb && centShift==-1 && (centr<min_cent+6 || centr>=max_cent+6)) continue;
    //if (doPbPb && centShift==0 && (centr<(min_cent+9) || centr>=(max_cent+9))) continue;
    //if (doPbPb && centShift==1 && (centr<(min_cent+12) || centr>=(max_cent+12))) continue;

    //if (JESsyst==+1) jt_pt = jt_pt_JEU_Up;
    //else if (JESsyst==-1) jt_pt = jt_pt_JEU_Down;

    if(reco_pT < min_jp_pt || reco_pT > max_jp_pt) continue;
    if(TMath::Abs(reco_eta) > max_jp_eta ) continue;
    //if(reco_mass < 2.6 || reco_mass > 3.5) continue;

    float absJtEta = TMath::Abs(jt_eta);
    if(absJtEta > max_jt_eta) continue;

    // get eta bin for resolution
    int etaBin = -1;
    if(absJtEta < 0.522) etaBin=0;
    else if(absJtEta < 0.783) etaBin=1;
    else if(absJtEta < 1.131) etaBin=2;
    else if(absJtEta < 1.305) etaBin=3;
    else if(absJtEta < 1.740) etaBin=4;
    else if(absJtEta < 1.930) etaBin=5;
    else etaBin=6;

    float sfVal = 0.;
    if(doPbPb) sfVal = SFvalPbPb[sfBin][etaBin];
    else sfVal = SFvalPP[sfBin][etaBin];
    
    //smear reco jet pt here, and recompute z reco before filling the matrix
    double pTres = TMath::Sqrt(c2+s2/(jt_ref_pt*(1.-gen_z)));    
    double sigmaSmear = pTres*TMath::Sqrt(sfVal*sfVal - 1.);
    double smearPt  = (1.-gen_z)*jt_ref_pt*rand->Gaus(0.,sigmaSmear);
    //cout <<"[INFO] sfVal = "<<sfVal<<", orig jtpt= "<<jt_pt<<", orig z = "<<", smearPt = "<<smearPt;
    if (!noSmearing)
      jt_pt = jt_pt+smearPt;
    if (jt_pt<reco_pT) jt_pt= reco_pT;
    z = reco_pT/jt_pt;
    //cout <<", new jtpt = "<<jt_pt<<"new z = "<<z<<endl;
    if (z<0) cout<<"[WARNING] jtpt = "<<jt_pt<<", z = "<<z<<endl;
    if(z<0) z=0;
    if (z>1) z=1;

    if(gen_z >= 1. && gen_z<1.001) gen_z = 0.99;
    if(z >= 1. && z < 1.001) z = 0.99;
    if(gen_z >=1. ) continue;
    if(z >= 1. ) continue;
    if(jt_pt < min_jetpt || jt_pt > max_jetpt) continue;
    if(jt_ref_pt < min_jetpt || jt_ref_pt > max_jetpt) continue;

    //response matrix fill
    fValue[0] = jt_ref_pt;
    fValue[1] = gen_z*gen_z;
    fValue[2] = jt_pt;
    fValue[3] = z*z;
    
    //if (nEv%100000==0) cout<<fValue[0]<<" and "<<fValue[1]<<" and "<<fValue[2]<<" and "<<fValue[3]<<endl;
  
    finalCorr = accxeff*pthat;
    fSparse->Fill(fValue,finalCorr);

    double ijt = (double)(TMath::Floor((jt_ref_pt-min_jetpt)/jetPt_gen_binWidth));
    double jjt = (double)(TMath::Floor((jt_pt-min_jetpt)/jetPt_reco_binWidth));

    double iz = (double)(TMath::Floor(gen_z*nBinZ_gen));
    double jz = (double)(TMath::Floor(z*nBinZ_reco));

    if(iz>(nBinZ_gen-1)) iz=nBinZ_gen-1;
    if(jz>(nBinZ_reco-1)) jz=nBinZ_reco-1;

    h_trMatrix->Fill(nBinZ_reco*jjt+jz+0.5,nBinZ_gen*ijt+iz+0.5,finalCorr);
  } // end of entry loop

  TH1D *cTruth1DX = dynamic_cast<TH1D*>(fSparse->Projection(0, "E"));
  TH1D *cTruth1DY = dynamic_cast<TH1D*>(fSparse->Projection(1, "E"));

  TCanvas * can6 = new TCanvas("can6","can6",1200,600);
  can6->Divide(2,1);
  can6->cd(1);
  cTruth1DX->Draw("EP");
  can6->cd(2);
  cTruth1DY->Draw("EP");
  can6->SaveAs("/home/llr/cms/henderson/UpsUnfolding/Unfolding/dataUnf/code/debugfSparse.png");


  // check 4D content

  cout << "number of bins :" << fSparse->GetNbins() <<endl;
  
  for(int iz = 0; iz < nBinZ_gen; iz++){
    for(int ijet = 0; ijet < nBinJet_gen; ijet++){

      float izMid = min_z+iz*z_gen_binWidth+z_gen_binWidth/2;
      float ijetMid = min_jetpt+ijet*jetPt_gen_binWidth+jetPt_gen_binWidth/2;

      double valIntegral = 0;
      
      double zUnfFrac = 1.;
      zUnfFrac = z_frac_JetPtBin[ijet][iz];
      
      for(int jjet = 0; jjet < nBinJet_reco; jjet++){
	for(int jz = 0; jz < nBinZ_reco; jz++){
	  float jjetMid = min_jetpt+jjet*jetPt_reco_binWidth+jetPt_reco_binWidth/2;
	  float jzMid = min_z+jz*z_reco_binWidth+z_reco_binWidth/2;
	  
	  const double x[4] = {ijetMid,izMid,jjetMid,jzMid};

	  int bin = fSparse->GetBin(x);
	  valIntegral+=fSparse->GetBinContent(bin);
	}
      }
      
      double scaleFactor = 0;
      if(valIntegral>0) scaleFactor = zUnfFrac/valIntegral;
      
      int countBin = 0;
      
      for(int jjet = 0; jjet < nBinJet_reco; jjet++){
	  for(int jz = 0; jz < nBinZ_reco; jz++){
	    float jjetMid = min_jetpt+jjet*jetPt_reco_binWidth+jetPt_reco_binWidth/2;
	    float jzMid = min_z+jz*z_reco_binWidth+z_reco_binWidth/2;
	    
	    const double x[4] = {ijetMid,izMid,jjetMid,jzMid};
	    int bin = fSparse->GetBin(x);

	    //this does not scale errors properly
	    //fSparse_newPrior->Fill(x, fSparse->GetBinContent(bin)*scaleFactor);

	    //do it bin by bin
	    
	    const double x2[4] = {ijetMid,izMid,jjetMid,jzMid};
	    int bin2 = fSparse_newZNorm->GetBin(x2);

	    double binCont = fSparse_newZNorm->GetBinContent(bin2);
	    double binErr  = fSparse_newZNorm->GetBinError(bin2);

	    double newBinCont = binCont+fSparse->GetBinContent(bin)*scaleFactor;
	    double newBinErr = binErr+fSparse->GetBinError(bin)*scaleFactor;

	    fSparse_newZNorm->SetBinContent(bin2, newBinCont);
	    fSparse_newZNorm->SetBinError(bin2, newBinErr);
	    
	  }
      }
    }
  }


  
  /*
  TH1D *gTruth1DX = dynamic_cast<TH1D*>(fSparse_newZNorm->Projection(0,"E"));
  TH1D *gTruth1DY = dynamic_cast<TH1D*>(fSparse_newZNorm->Projection(1,"E"));
  
  TCanvas * can1 = new TCanvas("can1","can1",1200,600);
  can1->Divide(2,1);
  can1->cd(1);
  gTruth1DX->Draw("EP");
  can1->cd(2);
  gTruth1DY->Draw("EP");
  can1->SaveAs("matrixTest.png");
  */
  

  /// jet pt normalization -> put it to match jet pt truth

  for(int ijet = 0; ijet < nBinJet_gen; ijet++){
      
    //float ijetMid = 10.*ijet+20;
    //float ijetMid = ijet+15.5;
    float ijetMid =  min_jetpt+ijet*jetPt_gen_binWidth+jetPt_gen_binWidth/2;
    double valIntegralOrig = 0;
    double valIntegralZNorm = 0;

    for(int iz = 0; iz < nBinZ_gen; iz++){
      //float izMid = iz*z_gen_binWidth+z_gen_binWidth/2;
      float izMid =  min_z+iz*z_gen_binWidth+z_gen_binWidth/2;
           
      for(int jjet = 0; jjet < nBinJet_reco; jjet++){
	for(int jz = 0; jz < nBinZ_reco; jz++){

	  float jjetMid = min_jetpt+jjet*jetPt_reco_binWidth+jetPt_reco_binWidth/2;;
	  float jzMid = min_z+jz*z_reco_binWidth+z_reco_binWidth/2;

	  const double x[4] = {ijetMid,izMid,jjetMid,jzMid};
	  int bin = fSparse->GetBin(x);
	  valIntegralOrig+=fSparse->GetBinContent(bin);

	  const double x2[4] = {ijetMid,izMid,jjetMid,jzMid};
	  int bin2 = fSparse_newZNorm->GetBin(x2);
	  valIntegralZNorm+=fSparse_newZNorm->GetBinContent(bin2);
	  	  
	}

      }

    }

    double scaleFactor = 0;
    //if(valIntegralZNorm>0) scaleFactor = valIntegralOrig/valIntegralZNorm;
    if(valIntegralZNorm>0) scaleFactor = valIntegralOrig/valIntegralZNorm;
    
      for(int iz = 0; iz < nBinZ_gen; iz++){
	float izMid = min_z+iz*z_gen_binWidth+z_gen_binWidth/2;
		
	for(int jjet = 0; jjet < nBinJet_reco; jjet++){
	  for(int jz = 0; jz < nBinZ_reco; jz++){

	    float jjetMid = min_jetpt+jjet*jetPt_reco_binWidth+jetPt_reco_binWidth/2;
	    float jzMid = min_z+jz*z_reco_binWidth+z_reco_binWidth/2;
	    
	    const double x3[4] = {ijetMid,izMid,jjetMid,jzMid};
	    int bin3 = fSparse_newZNorm->GetBin(x3);

	    //does not scale errors properly
	    //fSparse_newJetPtNorm->Fill(x3, fSparse_newZNorm->GetBinContent(bin3)*scaleFactor);

	    const double x4[4] = {ijetMid,izMid,jjetMid,jzMid};
	    int bin4 = fSparse_newJetPtNorm->GetBin(x4);

	    double binCont = fSparse_newJetPtNorm->GetBinContent(bin4);
	    double binErr  = fSparse_newJetPtNorm->GetBinError(bin4);

	    double newBinCont = binCont+fSparse_newZNorm->GetBinContent(bin3)*scaleFactor;
	    double newBinErr = binErr+fSparse_newZNorm->GetBinError(bin3)*scaleFactor;

	    fSparse_newJetPtNorm->SetBinContent(bin4, newBinCont);
	    fSparse_newJetPtNorm->SetBinError(bin4, newBinErr);
	    	    
	  }
	}
	
      }
        
  }


  /*
  TH2D *gTruth2D = dynamic_cast<TH2D*>(fSparse_newJetPtNorm->Projection(0,1));
  gTruth2D->Draw("colz");

  
  TH1D *gTruth1DX_jtPt = dynamic_cast<TH1D*>(fSparse_newJetPtNorm->Projection(0,"E"));
  TH1D *gTruth1DY_jtPt = dynamic_cast<TH1D*>(fSparse_newJetPtNorm->Projection(1,"E"));

  TCanvas * can3 = new TCanvas("can3","can3",1200,600);
  can3->Divide(2,1);
  can3->cd(1);
  gTruth1DX_jtPt->Draw("EP");
  can3->cd(2);
  gTruth1DY_jtPt->Draw("EP");
  can3->SaveAs("matrixTest_jtPt.png");
  */
  
  ///
  cout<< "saving the output file "<<outputfile<<endl;    
  TFile *file_forProf_varBin = new TFile(outputfile.c_str(),"RECREATE");

  h_jetpPtGen_jetPtReco->Write();
  
  fSparse->Write();
  fSparse_newZNorm->Write();
  fSparse_newJetPtNorm->Write();

  h_trMatrix->Write();
  h_jetpPtGen_jetPtReco->Write();
  h_z_gen->Write();

  hGenZJetPtCentBin->Write();
  hGenZJetPtLowBin->Write();
  hGenZJetPtHighBin->Write();
  
  hRespZJetPtCentBin->Write();
  hRespZJetPtLowBin->Write();
  hRespZJetPtHighBin->Write();

  file_forProf_varBin->Close();

}

void prepareInputsNewPrNominal(Int_t step = 1) { 
  //if (!matrixInv)
  //  prepare(true,true,step);
  //if (step<=(nSIter_pp+3) && centShift==0 && !doCent && !doPeri)
  prepare(true,false,step);
}


void preparePrior(bool doPbPb, string outputName){
  string treeFilename = Form("~/JpsiInJetsPbPb/Fitter/TreesForUnfolding/tree_%s_%s_NoBkg%s_AccEff_JEC%s.root","MCJPSINOPR",doPbPb?"PbPb":"PP",Form("_jetR%d",(int)(jetR*10)),doPbPb?"":"_updatedCorr");
  TFile *file = new TFile(treeFilename.c_str());
  TTree *t_unf = (TTree*)file->Get("treeForUnfolding");
  TH1D *h_z_jetPtBin[nBinJet_gen];

  TFile *outputFile = new TFile(outputName.c_str(),"RECREATE");
    for(int ibin = 0; ibin < nBinJet_gen; ibin++){
      float min_jetpt_bin = min_jetpt+ibin*(jetPt_gen_binWidth);
      float max_jetpt_bin = min_jetpt+(ibin+1)*(jetPt_gen_binWidth);
      h_z_jetPtBin[ibin] = new TH1D(Form("zDist_%i",ibin),Form("jt_ref_pt>%f && jt_ref_pt<%f",min_jetpt_bin,max_jetpt_bin),nBinZ_gen,min_z,max_z);
      t_unf->Draw(Form("gen_z>>zDist_%i",ibin),Form("corr_AccEff%s*corr_ptw*(jp_gen_pt>%f && jp_gen_pt<%f && fabs(jp_gen_rap)<%f && jt_ref_pt>%f && jt_ref_pt<%f && gen_z>%f && gen_z<%f && fabs(jt_ref_eta)<%f %s)",doPbPb?"":"_comp",min_jp_pt,max_jp_pt,max_jp_eta,min_jetpt_bin,max_jetpt_bin,min_z,max_z,max_jt_eta,doPbPb?Form("&& centr >%d && centr <%d",min_cent+9,max_cent+9):""));
      h_z_jetPtBin[ibin]->Write(Form("zDist_%i",ibin));
    }
    outputFile->Close();
    file->Close();
}
