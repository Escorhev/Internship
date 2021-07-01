#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
//#include "inputParams.h"
#endif

void makeRebinMatrix(){

  // find the transfer matrix to take a measurement of 50 bins in z and 15 bins in pT to one with 5 bins in z and 3 bins in pT.  

  TH2D *h2dpt = new TH2D("h2dpt","h2dpt",nBinJet_gen,min_jetpt,max_jetpt,nBinJet_reco,min_jetpt,max_jetpt);
  TH2D *h2dz = new TH2D("h2dz","h2dz",nBinZ_gen,min_z,max_z,nBinZ_reco,min_z,max_z);

  int fDim = 4.;  
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

  Int_t* bins_sparce2 = new Int_t[fDim];
  Double_t *xmin_sparce2 = new Double_t[fDim];
  Double_t *xmax_sparce2 = new Double_t[fDim];

  bins_sparce2[2] = nBinJet_gen;
  xmin_sparce2[2] = min_jetpt;
  xmax_sparce2[2] = max_jetpt;

  bins_sparce2[3] = nBinZ_gen;
  xmin_sparce2[3] = min_z;
  xmax_sparce2[3] = max_z;

  bins_sparce2[0] = nBinJet_reco;
  xmin_sparce2[0] = min_jetpt;
  xmax_sparce2[0] = max_jetpt;

  bins_sparce2[1] = nBinZ_reco;
  xmin_sparce2[1] = min_z;
  xmax_sparce2[1] = max_z;
  
  THnSparseF * fSparseInv = new THnSparseF("hsInv", "hsInv", fDim, bins_sparce2, xmin_sparce2, xmax_sparce2);
  fSparseInv->Sumw2();
  fSparseInv->CalculateErrors();
  

  for(int ipt = 0; ipt < nBinJet_reco; ipt++){    
    for(int iz = 0; iz < nBinZ_reco; iz++){

      float izMid = min_z+iz*z_reco_binWidth+z_reco_binWidth/2;
      float ijetMid = min_jetpt+ipt*jetPt_reco_binWidth+jetPt_reco_binWidth/2;

      int jzlo = iz*nBinZ_gen/nBinZ_reco;
      int jzhi = (iz+1)*nBinZ_gen/nBinZ_reco;
      int jptlo = ipt*nBinJet_gen/nBinJet_reco;
      int jpthi = (ipt+1)*nBinJet_gen/nBinJet_reco;
     
      for(int jpt = jptlo; jpt < jpthi; jpt++){	
	for(int jz = jzlo; jz < jzhi; jz++){
	  
	  h2dpt->SetBinContent(jpt+1,ipt+1,1.);
	  h2dz->SetBinContent(jz+1,iz+1,1.);

	  float jzMid = min_z+jz*z_gen_binWidth+z_gen_binWidth/2;
	  float jjetMid = min_jetpt+jpt*jetPt_gen_binWidth+jetPt_gen_binWidth/2;
	  	  
	  const double x[4] = {jjetMid,jzMid,ijetMid,izMid};
	  int bin = fSparse->GetBin(x);
	  fSparse->SetBinContent(bin, 1.);

	  const double x2[4] = {ijetMid,izMid,jjetMid,jzMid};
	  int bin2 = fSparseInv->GetBin(x2);
	  fSparseInv->SetBinContent(bin2, 1.);
	  
	  
	  //cout << "bin = " << bin << " izMid = " << izMid << " jzMid = " << jzMid << endl;
	}
      }
      
    }
  }
  
  TCanvas *c1=new TCanvas("c1","c1",600,600);
  h2dpt->Draw("zcol");
  TCanvas *c2=new TCanvas("c2","c2",600,600);
  h2dz->Draw("zcol");

  TCanvas *c1_4Dpr=new TCanvas("c1_4Dpr","c1_4Dpr",600,600);
  fSparseInv->Projection(0,2)->Draw("zcol");
  TCanvas *c2_4Dpr=new TCanvas("c2_4Dpr","c2_4Dpr",600,600);
  fSparseInv->Projection(1,3)->Draw("zcol");


  TFile *fout = new TFile(Form("%s/dataUnf/unfInput/diag4DMatrixInv.root",unfPath.c_str()),"RECREATE");
  fSparse->Write("diag4DMatrix");  
  fSparseInv->Write("diag4DMatrixInv");
}
