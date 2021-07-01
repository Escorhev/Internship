#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include "inputParams.h"
#endif

void prepare(bool doPrompt = false, bool doPbPb = true, bool doTrain=true, Int_t stepNumber = 1){
  printInput();
  if (!setCaseTag()) return;
  
  string filename = "";
  string outputfile = "";
  string filenamePrevStep = "";

  gSystem->mkdir(Form("%s/mcUnf/unfInput1D",unfPath.c_str()));
  gSystem->mkdir(Form("%s/mcUnf/unfInput1D/step%i",unfPath.c_str(),stepNumber));
  gSystem->mkdir(Form("%s/mcUnf/unfOutput1D",unfPath.c_str()));
  gSystem->mkdir(Form("%s/mcUnf/unfOutput1D/step%i",unfPath.c_str(),stepNumber));

  //  filename = Form("~/JpsiInJetsPbPb/Fitter/TreesForUnfolding/tree_%s_%s_NoBkg%s_AccEff_JEC.root",doPrompt?"MCJPSIPR":"MCJPSINOPR",doPbPb?"PbPb":"PP",mc2015?"":Form("_jetR%d",(int)(jetR*10)));
  filename = "/data_CMS/cms/henderson/tree_MC1S_PP_jetR4_AccEff_JEC.root";
  outputfile = Form("%s/mcUnf/unfInput1D/step%i/unfolding_4D_%s_%s_%s_%diter_%dz%dptBins%dz%dptMeasBins%s.root",unfPath.c_str(), stepNumber, doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", doTrain?"Train":"Test", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, caseTag.c_str());
  
  if(stepNumber > 1) filenamePrevStep = Form("%s/mcUnf/unfOutput1D/step%i/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBins%s.root",unfPath.c_str(), stepNumber-1,doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, caseTag.c_str());

  cout <<"filename = "<<filename<<endl;
  cout <<"outputfile = "<<outputfile<<endl;
  cout <<"filenamePrevStep = "<<filenamePrevStep<<endl;

  TFile *file = new TFile(filename.c_str());
  TTree *t_unf = (TTree*)file->Get("treeForUnfolding");

  int n_entries = t_unf->GetEntries();
    
  Double_t z_frac_JetPtBin_dataMC[nBinJet_reco][nBinZ_reco];
  TH2D* mcDataWeight = NULL;
  if (dataDist) {
    if (stepNumber==1) {
      TFile *fileData = TFile::Open(Form("%s/dataUnf/data_results/meas_%s_data_%s_statErrs.root",unfPath.c_str(), doPbPb?"PbPb":"PP",doPrompt?"prompt":"prompt")); //use prompt data always 
      TH2D *hMeasuredData = (TH2D*) fileData->Get("h_Meas;1");
      TH2D *hMeasuredMC = (TH2D*) hMeasuredData->Clone("hMeasuredMC");
      mcDataWeight = (TH2D*) hMeasuredData->Clone("mcDataWeight");
      t_unf->Draw("jt_pt:z>>hMeasuredMC",Form("corr_AccEff*corr_ptw*((centr>=%d && centr<%d) && (jp_pt>%f && jp_pt<%f) && (TMath::Abs(jp_eta)<%f) && (jp_mass>2.6 && jp_mass<3.5) && (jt_pt>%f && jt_pt<%f) && (jt_ref_pt>%f && jt_ref_pt<%f) && (TMath::Abs(jt_eta)<2.4) && gen_z <=1. && z<=1)",centShift==0?10:centShift==-1?5:15, centShift==0?190:centShift==-1?185:195, min_jp_pt, max_jp_pt, max_jp_eta, min_jetpt, max_jetpt, min_jetpt, max_jetpt));
      mcDataWeight->Divide(hMeasuredMC);
    }//end of (stepNumber==1)
    else {
      TFile *fileW = TFile::Open(Form("%s/mcUnf/unfInput1D/step1/unfolding_4D_%s_%s_%s_%diter_%dz%dptBins%dz%dptMeasBins%s.root",unfPath.c_str(), doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", doTrain?"Train":"Test", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, caseTag.c_str()));
      mcDataWeight = (TH2D*) fileW->Get("mcDataWeight");
    }
  }//end of if (dataDist)
  
  Double_t z_frac_JetPtBin[nBinJet_gen][nBinZ_gen];
  
  //only for prompt fwd in steps > 1, until the stats are better in mc
  Double_t *z_frac_allJetPtBins = new Double_t[nBinZ_gen];
  
  // if this is step 2, 3, 4 ... take the z unfolded and use it for the prior, else put z prior to flat    
  if(stepNumber > 1) {
    TFile *filePrevStep = new TFile(filenamePrevStep.c_str());
    // cout << "file OK" << endl;
    TH1D *h_z_unf_jetPtBin[nBinJet_gen];

    for(int ibin = 0; ibin < nBinJet_gen; ibin++){
      h_z_unf_jetPtBin[ibin] = (TH1D*)filePrevStep->Get(Form("hMUnf_%i_Iter%d;1",ibin,nIter));
    }
    
    TH1D *h_z_allBins = (TH1D*)h_z_unf_jetPtBin[0]->Clone();
    for(int ibin = 1; ibin < nBinJet_gen; ibin++){
      h_z_allBins->Add(h_z_unf_jetPtBin[ibin]);
    }

    h_z_allBins->Scale(1/h_z_allBins->Integral());

    TCanvas * can4 = new TCanvas("can4","can4",1200,600);
    h_z_allBins->Draw("EP");
    can4->SaveAs("/home/llr/cms/henderson/UpsUnfolding/Unfolding/mcUnf/code/debugzBins.png");

    for(int ibin = 1; ibin < nBinJet_gen; ibin++){
      cout<<"intégrale = "<<h_z_unf_jetPtBin[ibin]->Integral()<<endl;
      h_z_unf_jetPtBin[ibin]->Scale(1/h_z_unf_jetPtBin[ibin]->Integral());
    }
            
    for(int iz = 1; iz <= nBinZ_gen; iz++){
      z_frac_allJetPtBins[iz-1] = h_z_allBins->GetBinContent(iz);
      //cout << "iz = " << iz << " allJtPtBins content = " << z_frac_allJetPtBins[iz-1] << endl;

      for(int ijtpt = 0; ijtpt < nBinJet_gen; ijtpt++){
	z_frac_JetPtBin[ijtpt][iz-1] = h_z_unf_jetPtBin[ijtpt]->GetBinContent(iz);
	//cout <<"ijtpt = " << ijtpt << " z_frac_JetPtBin[ijtpt][iz-1] = " << z_frac_JetPtBin[ijtpt][iz-1] << endl;
      }
    }
  }
  else{
    for(int i = 1; i <= nBinZ_gen; i++){
      for(int ijtpt = 0; ijtpt < nBinJet_gen; ijtpt++){
	z_frac_JetPtBin[ijtpt][i-1] = 1.0;
      }
    }    
  }
  
  //J/Psi variables:
  float jp_pt = 0.; t_unf->SetBranchAddress("jp_pt", &jp_pt);
  float jp_eta = 0.; t_unf->SetBranchAddress("jp_eta", &jp_eta);
  float jp_mass = 0.; t_unf->SetBranchAddress("jp_mass", &jp_mass);

  //jet reco variables:
  float jt_pt = 0.; t_unf->SetBranchAddress("jt_pt", &jt_pt);
  float jt_eta = 0.; t_unf->SetBranchAddress("jt_eta", &jt_eta);

  //centrality
  int centr = 0; t_unf->SetBranchAddress("centr", &centr);
  
  //correction variables:
  float corr_AccEff = 0;//t_unf->SetBranchAddress("corr_AccEff", &corr_AccEff);
  float corr_ptw = 0;//t_unf->SetBranchAddress("corr_ptw",&corr_ptw);
  float finalCorr = 0;
      
  //z corrected
  float z = 0.; t_unf->SetBranchAddress("z", &z);

  //jet gen variables
  float jt_ref_pt = 0.; t_unf->SetBranchAddress("jt_ref_pt", &jt_ref_pt);

  //gen z:
  float gen_z = 0.; t_unf->SetBranchAddress("gen_z", &gen_z);

  // histograms for renormalizing response matrix:
  TH2F * hRespZJetPtCentBin = new TH2F ("hRespZJetPtCentBin","z_gen vs z_raw for 30 < refpt < 40",nBinZ_gen,min_z,max_z,nBinZ_gen,min_z,max_z);
  TH2F * hRespZJetPtLowBin = new TH2F ("hRespZJetPtLowBin","z_gen vs z_raw for 20 < refpt < 30",nBinZ_gen,min_z,max_z,nBinZ_gen,min_z,max_z);
  TH2F * hRespZJetPtHighBin = new TH2F ("hRespZJetPtHighBin","z_gen vs z_raw for 40 < refpt < 50",nBinZ_gen,min_z,max_z,nBinZ_gen,min_z,max_z);

  TH1D * h_z_gen = new TH1D ("h_z_gen","z_gen for 20 < refpt < 50",nBinZ_gen,min_z,max_z);
  TH1D * hGenZJetPtCentBin = new TH1D ("hGenZJetPtCentBin","z_gen for 30 < refpt < 40",nBinZ_gen,min_z,max_z);
  TH1D * hGenZJetPtLowBin = new TH1D ("hGenZJetPtLowBin","z_gen for 20 < refpt < 30",nBinZ_gen,min_z,max_z);
  TH1D * hGenZJetPtHighBin = new TH1D ("hGenZJetPtHighBin","z_gen for 40 < refpt < 50",nBinZ_gen,min_z,max_z);

  TH2F * h_jetpPtGen_jetPtReco = new TH2F ("h_jetpPtGen_jetPtReco","gen vs reco corr jet pt",nBinJet_reco,min_jetpt,max_jetpt,nBinJet_reco,min_jetpt,max_jetpt);

  hRespZJetPtCentBin->Sumw2();
  hRespZJetPtLowBin->Sumw2();
  hRespZJetPtHighBin->Sumw2();
  h_jetpPtGen_jetPtReco->Sumw2();

  // 4D histogram creation
  
  int fDim = 2.;
  Double_t fValue[fDim];
  
  Int_t* bins_sparce = new Int_t[fDim];
  Double_t *xmin_sparce = new Double_t[fDim];
  Double_t *xmax_sparce = new Double_t[fDim];

  bins_sparce[0] = nBinJet_gen*nBinZ_gen;
  xmin_sparce[0] = 0;
  xmax_sparce[0] = nBinJet_gen*nBinZ_gen;

  bins_sparce[1] = nBinJet_reco*nBinZ_reco;
  xmin_sparce[1] = 0;
  xmax_sparce[1] = nBinJet_reco*nBinZ_reco;
  
  //initial not normalized 4D  
  THnSparseF * fSparse = new THnSparseF("hs", "hs", fDim, bins_sparce, xmin_sparce, xmax_sparce);
  fSparse->Sumw2();
  fSparse->CalculateErrors();

  THnSparseF * fSparseOrig = new THnSparseF("hs_orig", "hs_orig", fDim, bins_sparce, xmin_sparce, xmax_sparce);
  fSparseOrig->Sumw2();
  fSparseOrig->CalculateErrors();

  //4D with z gen flat
  THnSparseF * fSparse_newZNorm = new THnSparseF("hs_newZNorm", "hs_newZNorm", fDim, bins_sparce, xmin_sparce, xmax_sparce);
  fSparse_newZNorm->Sumw2();
  fSparse_newZNorm->CalculateErrors();

  //4D with z gen flat and gen pt = truth
  THnSparseF * fSparse_newJetPtNorm = new THnSparseF("hs_newJetPtNorm", "hs_newJetPtNorm", fDim, bins_sparce, xmin_sparce, xmax_sparce);
  fSparse_newJetPtNorm->Sumw2();
  fSparse_newJetPtNorm->CalculateErrors();
  
  double_t tabZ_gen[nBinZ_gen+1];
  for(int k = 0; k <= nBinZ_gen; k++){
    float fk = k;
    tabZ_gen[k] = fk/nBinZ_gen;
  }
  double_t tabJet_gen[nBinJet_gen+1];
  for(int k = 0; k <= nBinJet_gen; k++){
    float fk = k;
    tabJet_gen[k] = 10. + 30.*fk/nBinJet_gen;
  }
  double_t tabZ_reco[nBinZ_reco+1];
  for(int k = 0; k <= nBinZ_reco; k++){
    float fk = k;
    tabZ_reco[k] = fk/nBinZ_reco;
  }
  double_t tabJet_reco[nBinJet_reco+1];
  for(int k = 0; k <= nBinJet_reco; k++){
    float fk = k;
    tabJet_reco[k] = 10. + 30.*fk/nBinJet_reco;
  }

  for(int nEv = 0; nEv < n_entries; nEv++){
        
    if (!sameSample) {
      if (doPbPb) {
	if(doTrain && nEv%10 >= 8) continue;
	if(!doTrain && nEv%10 < 8) continue;
      }
      else {
	if(doTrain && nEv%2 == 1) continue;
	if(!doTrain && nEv%2 == 0) continue;
      }
    }
    
    t_unf->GetEntry(nEv);

    //if (doPbPb && (centr<10 || centr>=190)) continue;
    if (!doTrain) {
      if (doPbPb && centShift==-1 && (centr<min_cent || centr>=max_cent)) continue;
      if (doPbPb && centShift==0 && (centr<min_cent+10 || centr>=max_cent+10)) continue;
      if (doPbPb && centShift==1 && (centr<min_cent+20 || centr>=max_cent+20)) continue;
    }
    else {
      if (doPbPb && (centr<min_cent+10 || centr>=max_cent+10)) continue;
    }
    // check event content
    if(jp_pt < min_jp_pt || jp_pt > max_jp_pt) continue;    
    if(TMath::Abs(jp_eta) > max_jp_eta ) continue;
    // if(jp_mass < 2.6 || jp_mass > 3.5) continue;
    if(jt_pt < min_jetpt_real || jt_pt > max_jetpt_real) continue;
    if(jt_ref_pt < min_jetpt_real || jt_ref_pt > max_jetpt_real) continue;
    if(TMath::Abs(jt_eta) > max_jt_eta) continue;    

    double gen_z2 = gen_z*gen_z;
    if(gen_z == 1.) gen_z2 = 0.99;

    double z2 = z*z;
    if(z == 1.) z2 = 0.99;

    if(gen_z2 >=1. ) continue;
    if(z2 >= 1. ) continue;
    
    double cheatBin_gen = 0.5;
    for(int l = 0; l < nBinZ_gen; l++){
      if (tabZ_gen[l]<=gen_z2 && gen_z2<tabZ_gen[l+1])
	cheatBin_gen += l;
    }
    for(int l= 0; l < nBinJet_gen; l++){
      if (tabJet_gen[l]<=jt_ref_pt && jt_ref_pt<tabJet_gen[l+1])
	cheatBin_gen += nBinZ_gen*l;
    }  

    double cheatBin_reco = 0.5;
    for(int l= 0; l < nBinZ_reco; l++){
      if (tabZ_reco[l]<=z2 && z2<tabZ_reco[l+1])
	cheatBin_reco += l;
    }
    for(int l= 0; l < nBinJet_reco; l++){
      if (tabJet_reco[l]<=jt_pt && jt_pt<tabJet_reco[l+1])
	cheatBin_reco += nBinZ_reco*l;
    }
      
    //response matrix fill
    fValue[0] = cheatBin_gen;
    fValue[1] = cheatBin_reco;
    
      /*if (jt_ref_pt > max_jetpt-jetPt_reco_binWidth)
    fValue[0] = (max_jetpt-jetPt_reco_binWidth)+(jt_ref_pt-(max_jetpt-jetPt_reco_binWidth))*squeezeHigh;//jt_ref_pt;
    else if (jt_ref_pt < min_jetpt+jetPt_reco_binWidth)
    fValue[0] = (min_jetpt+jetPt_reco_binWidth)-((min_jetpt+jetPt_reco_binWidth)-jt_ref_pt)*squeezeLow;//jt_ref_pt;
    if (jt_pt > max_jetpt-jetPt_reco_binWidth)
	fValue[2] = (max_jetpt-jetPt_reco_binWidth)+(jt_pt-(max_jetpt-jetPt_reco_binWidth))*squeezeHigh;//jt_pt;
      else if (jt_pt < min_jetpt+jetPt_reco_binWidth)
	fValue[2] = (min_jetpt+jetPt_reco_binWidth)-((min_jetpt+jetPt_reco_binWidth)-jt_pt)*squeezeLow;//jt_pt;
      */

      //if (jt_ref_pt<midLowerPt)
      //cout <<"jt_ref_pt = "<<jt_ref_pt<<", fValue[0] = "<< fValue[0]<<", jt_pt = "<<jt_pt<<", fValue[2] = "<<fValue[2]<<endl; 
    
    //cout<<fValue[0]<<" "<<fValue[1]<<" fValue"<<endl;
    finalCorr = corr_AccEff*corr_ptw;
    fSparseOrig->Fill(fValue);
    //if (dataDist) finalCorr = finalCorr*mcDataWeight->GetBinContent(mcDataWeight->FindBin(z2,jt_pt));
    fSparse->Fill(fValue);
    /*
    if (jt_ref_pt > midUpperPt)
      fValue[0] = midUpperPt+0.5*(jt_ref_pt-midUpperPt);//jt_ref_pt;
    else if (jt_ref_pt < midLowerPt)
      fValue[0] = midLowerPt-0.5*(midLowerPt-jt_ref_pt);//jt_ref_pt;
      if (jt_pt > midUpperPt)
	fValue[2] = midUpperPt+0.5*(jt_pt-midUpperPt);//jt_pt;
      else if (jt_pt < midLowerPt)
	fValue[2] = midLowerPt-0.5*(midLowerPt-jt_pt);//jt_pt;
    */
  }

  
  // check 4D content
  cout << "number of bins :" << fSparse->GetNbins() <<endl;
  fSparse->Print();
  
  for(int iZ = 0; iZ < nBinZ_gen; iZ++){
    for(int iJ = 0; iJ < nBinJet_gen; iJ++){

      float iZp = iZ + 0.5;
      float iJp = iJ;
      double valIntegral = 0;
      //cout << "iz = " << iz << " izMid = " << izMid << endl;
      double zUnfFrac = 1.;
      zUnfFrac = z_frac_JetPtBin[iJ][iZ];

      for(int j = 0; j < nBinJet_reco*nBinZ_reco; j++){

	float jp  = j + 0.5;
	const double x[2] = {iZp + nBinZ_gen*iJp, jp};

	int bin = fSparse->GetBin(x);
	valIntegral+=fSparse->GetBinContent(bin);
      }
      
      double scaleFactor = 0;
      //if(valIntegral>0) scaleFactor = 1./valIntegral;
      if(valIntegral>0) scaleFactor = zUnfFrac/valIntegral;
      int countBin = 0;

      for(int j = 0; j < nBinJet_reco*nBinZ_reco; j++){
      
	float jp = j + 0.5;
	const double x[2] = {iZp + nBinZ_gen*iJp, jp};
	int bin = fSparse->GetBin(x);
    
	const double x2[2] = {iZp + nBinZ_gen*iJp, jp};
	int bin2 = fSparse_newZNorm->GetBin(x2);

	double binCont = fSparse_newZNorm->GetBinContent(bin2);
	double binErr  = fSparse_newZNorm->GetBinError(bin2);

	double newBinCont = binCont+fSparse->GetBinContent(bin)*scaleFactor;
	double newBinErr = binErr+fSparse->GetBinError(bin)*scaleFactor;
	//cout<<"binCont = "<<binCont<<endl;
	//cout<<"bin2 = "<<bin2<<endl;
	if(binCont>0.000001){
	  cout<<"ip, jp = "<<iZp + nBinZ_gen*iJp<<", "<<jp<<endl;
	  cout<<"old bin content ="<<binCont<<"; new bin content = "<<newBinCont<<"; sf = "<<scaleFactor<<endl;
	}
      fSparse_newZNorm->SetBinContent(bin2, newBinCont);
      fSparse_newZNorm->SetBinError(bin2, newBinErr);
	    
      }
    }
  }
  
  TH1D *gTruth1D = dynamic_cast<TH1D*>(fSparse_newZNorm->Projection(0,"E"));
  TCanvas * can1 = new TCanvas("can1","can1",1200,600);
  gTruth1D->Draw("EP");
  can1->SaveAs("matrixTest1D.png");
  
  
  /// jet pt normalization -> put it to match jet pt truth

  for(int iJ = 0; iJ < nBinJet_gen; iJ++){
    double valIntegralOrig = 0;
    double valIntegralZNorm = 0;

    for(int iZ = 0; iZ < nBinZ_gen; iZ++){
    
      float iZp = iZ + 0.5;
      float iJp = iJ;

      for(int j = 0; j < nBinJet_reco*nBinZ_reco; j++){

	float jp = j + 0.5;
        const double x[2] = {nBinZ_gen*iJp + iZp, jp};
        int bin = fSparse->GetBin(x);
	//cout<<fSparse->GetBinContent(bin)<<" and the bin "<<bin<<endl;
        valIntegralOrig+=fSparse->GetBinContent(bin);

        const double x2[2] = {nBinZ_gen*iJp + iZp, jp};
        int bin2 = fSparse_newZNorm->GetBin(x2);
        valIntegralZNorm+=fSparse_newZNorm->GetBinContent(bin2);
      }
    }

    double scaleFactor = 0;
    if(valIntegralZNorm>0) scaleFactor = valIntegralOrig/valIntegralZNorm;
    //cout<<"vizN "<<valIntegralZNorm<<endl;
    for(int iZ = 0; iZ < nBinZ_gen; iZ++){
      
      float iZp = iZ + 0.5;
      float iJp = iJ;
      
      for(int j = 0; j < nBinJet_reco; j++){

	float jp = j + 0.5;
	const double x3[2] = {nBinZ_gen*iJp + iZp, jp};
	int bin3 = fSparse_newZNorm->GetBin(x3);

	const double x4[2] = {nBinZ_gen*iJp + iZp, jp};
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
  
  fSparse_newJetPtNorm->Print();
  TH1D *gTruth1D_jtPt = dynamic_cast<TH1D*>(fSparse_newJetPtNorm->Projection(0,"E"));
  TCanvas * can3 = new TCanvas("can3","can3",1200,600);
  gTruth1D_jtPt->Draw("EP");
  can3->SaveAs("matrixTest_jtPt1D.png");

  ///                                                                                                                         
  TFile *file_forProf_varBin = new TFile(outputfile.c_str(),"RECREATE");

  h_jetpPtGen_jetPtReco->Write();

  fSparse->Write();
  fSparseOrig->Write();
  fSparse_newZNorm->Write();
  fSparse_newJetPtNorm->Write();
  h_z_gen->Write();

  hGenZJetPtCentBin->Write();
  hGenZJetPtLowBin->Write();
  hGenZJetPtHighBin->Write();
  
  hRespZJetPtCentBin->Write();
  hRespZJetPtLowBin->Write();
  hRespZJetPtHighBin->Write();

  if (dataDist)
    mcDataWeight->Write();
  
  file_forProf_varBin->Close();
}

void prepareInputs1D(Int_t step = 1){
  //prepare(bool doPrompt, bool doPbPb, bool doTrain, Int_t stepNumber)
  if (step<=nSIter_pp && centShift == 0 && !doCent && !doPeri){
    prepare(true,false,true,step);
    prepare(true,false,false,step);
  }  
}
