#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
//#include "inputParams.h"
#endif

void operate(bool doPrompt = true, bool doPbPb = true, int toyNumber = 1){

  if (!setSystTag(doPbPb)) return;
  string filename = "";
  string outputfile = "";
  int stepNumber = nSIter;
  if (!doPbPb) stepNumber = nSIter_pp;

  gSystem->mkdir(Form("%s/dataUnf/unfOutput/step1/matrixOper", unfPath.c_str()));
  
  filename = Form("%s/dataUnf/unfOutput/step1/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin%s.root",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, systTag.c_str());

  outputfile = Form("%s/dataUnf/unfOutput/step1/matrixOper/matrixOperation_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin%s.root",unfPath.c_str(), doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, systTag.c_str());
  
  cout <<"output = "<<outputfile<<endl;     
  TFile *file = new TFile(filename.c_str());

  // get forward transfer matrix
  TH2D *h_matrix  = (TH2D *)file->Get("trmatrix");
  
  // 2D histo for its errors
  TH2D *h_matrixErr  = (TH2D *)h_matrix->Clone("h_matrixErr");
  
  // get inverse transfer matrix
  TMatrixT<double> *invMatrix  = (TMatrixT<double> *)file->Get(Form("invtrmat%d",nIter));

  TH2D *h2_UnfData  = (TH2D*)file->Get(Form("hReco_Iter%d;1",nIter));
  TH2D *h2_MeasData  = (TH2D*)file->Get("fh2MeasData;1");

  int midStart = h2_UnfData->GetYaxis()->FindBin(midLowerPt);
  int midEnd = h2_UnfData->GetYaxis()->FindBin(midUpperPt)-1;

  //cout <<"midStart = "<<midStart<<", midEnd"<<midEnd<<endl;
  TH1D *h1_zUnf = (TH1D*)h2_UnfData->ProjectionX("h1_zUnf",midStart,midEnd);//6,10);
  
  //  TH1D *h1_zUnf = (TH1D*)file->Get("hMUnf_1_Iter3;1");
  
  TH2D *h_invMatrix = new TH2D(*invMatrix);
  h_invMatrix->Sumw2();
  
  TH2D *h_invMatrixErr = (TH2D *) h_invMatrix->Clone("h_invMatrixErr");

  // 4D inv matrix, norm so its proj gives the truth
  TH2D *h_invMatrixOut = (TH2D*)h_invMatrix->Clone("h_invMatrixOut");
  h_invMatrixOut->Reset();
  
  //cout << "here" << endl;
  
  //normalize 2d histogram, such that each measured bin integrates to unity
  for(int jjet = 0; jjet < nBinJet_reco; jjet++){
    for(int jz = 0; jz < nBinZ_reco; jz++){

      int binY = nBinZ_reco*jjet+jz+1;
      
      float valIntegral = 0;
      float errIntegral = 0;

      for(int ijet = 0; ijet < nBinJet_gen; ijet++){
	for(int iz = 0; iz < nBinZ_gen; iz++){
	  int binX = nBinZ_gen*ijet+iz+1;
	  
	  valIntegral += h_invMatrix->GetBinContent(binY,binX);
	  errIntegral += h_invMatrixErr->GetBinError(binY,binX);
	  
	}
      }

      for(int ijet = 0; ijet < nBinJet_gen; ijet++){
	for(int iz = 0; iz < nBinZ_gen; iz++){
	  int binX = nBinZ_gen*ijet+iz+1;
	  
	  if(valIntegral>0)h_invMatrix->SetBinContent(binY,binX, h_invMatrix->GetBinContent(binY,binX)/valIntegral);
	  if(errIntegral>0)h_invMatrixErr->SetBinContent(binY,binX, h_invMatrixErr->GetBinError(binY,binX)/errIntegral);

	}
      }

    }
  }

  //multiply back measured to see if we can get unfolded
  
  TH2D *h2_UnfDataFromMatrix = (TH2D*)h2_UnfData->Clone("h2_UnfDataFromMatrix");
  h2_UnfDataFromMatrix->Reset();
    
  for(int jjet = 0; jjet < nBinJet_reco; jjet++){
    for(int jz = 0; jz < nBinZ_reco; jz++){
      int binY = nBinZ_reco*jjet+jz+1;

      float measVal = h2_MeasData->GetBinContent(jz+1,jjet+1);
      float measErr = h2_MeasData->GetBinError(jz+1,jjet+1);

      double sumErr = 0.;

      for(int ijet = 0; ijet < nBinJet_gen; ijet++){
	for(int iz = 0; iz < nBinZ_gen; iz++){
	  int binX = nBinZ_gen*ijet+iz+1;
	  
	  if(h_invMatrixErr->GetBinContent(binY,binX)>0) sumErr+=h_invMatrixErr->GetBinContent(binY,binX)*h_invMatrixErr->GetBinContent(binY,binX);
	}
      }
      if(sumErr>0) sumErr = sqrt(sumErr);

      for(int ijet = 0; ijet < nBinJet_gen; ijet++){
	for(int iz = 0; iz < nBinZ_gen; iz++){
	  int binX = nBinZ_gen*ijet+iz+1;
	  
	  float valWeight = h_invMatrix->GetBinContent(binY,binX);	  
	  float errWeight = h_invMatrixErr->GetBinContent(binY,binX);
	  
	  float inVal = h2_UnfDataFromMatrix->GetBinContent(iz+1,ijet+1);
	  float inErr = h2_UnfDataFromMatrix->GetBinError(iz+1,ijet+1);

	  h2_UnfDataFromMatrix->SetBinContent(iz+1,ijet+1,inVal+valWeight*measVal);
	  h_invMatrixOut->SetBinContent(binY,binX,valWeight*measVal);

	  float outErr = 0.;
	  if(sumErr>0)outErr=errWeight*measErr/sumErr;
	  h2_UnfDataFromMatrix->SetBinError(iz+1,ijet+1,sqrt(outErr*outErr+inErr*inErr));
	  h_invMatrixOut->SetBinError(binY,binX,outErr);
	}
      }

    }
  }


  //fix the errors in h_invMatrixOut
  for(int jjet = 0; jjet < nBinJet_reco; jjet++){
    for(int jz = 0; jz < nBinZ_reco; jz++){
      int binY = nBinZ_reco*jjet+jz+1;

      for(int ijet = 0; ijet < nBinJet_gen; ijet++){
	for(int iz = 0; iz < nBinZ_gen; iz++){
	  int binX = nBinZ_gen*ijet+iz+1;
	  
	  float matrixVal = h_matrix->GetBinContent(binY,binX);
	  float matrixErr = h_matrix->GetBinError(binY,binX);
	  //cout <<"matrixVal = "<<matrixVal<<", matrixErr = "<<endl;
	  
	  //if(matrixVal>0)cout << "matrix relative error = " << matrixErr/matrixVal << endl;
	  if(matrixErr/matrixVal>0.99) cout <<"rel error is 100% in the bin with jjet = " << jjet << " , jz = " << jz << " , ijet = " << ijet << " , iz = " << iz << endl; 
	  
	  float invMatrixVal = h_invMatrixOut->GetBinContent(binY,binX);
	  //cout << "invMatrixVal  = " << invMatrixVal << endl;
	  float invMatrixErr = 0;
	  if(matrixVal>0) invMatrixErr = matrixErr*invMatrixVal/matrixVal;

	  //if(matrixVal>0)cout << "inv matrix relative error = " << invMatrixErr/invMatrixVal << endl;
	  if(invMatrixErr/invMatrixVal>0.99) cout<<"rel error is 100% in the bin with jjet = " << jjet << " , jz = " << jz << " , ijet = " << ijet << " , iz = " <<iz << endl;

	  h_invMatrixOut->SetBinError(binY,binX,invMatrixErr);  

	}
      }

    }
  }

  //project the nominal, take the unfolded nominal jet pt bin, loop over the unfolded z bins, sum over all measured bins.

  TH1D *h_zUnfNomBin = (TH1D *)h1_zUnf->Clone("h_zUnfNomBin");
  h_zUnfNomBin->Reset();
  float zBinsVals[nBinZ_gen];  
  
  for(int iz = 0; iz < nBinZ_gen; iz++){
    //int ijet = 1;
    zBinsVals[iz] = 0;
    for(int ijet = 0; ijet < nBinJet_gen; ijet++){
      int binX = nBinZ_gen*ijet+iz+1;
          
      for(int jjet = 0; jjet < nBinJet_reco; jjet++){
	for(int jz = 0; jz < nBinZ_reco; jz++){
	  int binY = nBinZ_reco*jjet+jz+1;
	  
	  if(ijet>=midLowerId && ijet<midLowerId+nBinJet_gen/nBinJet_reco) {
	    zBinsVals[iz] += h_invMatrixOut->GetBinContent(binY,binX);
	    //cout << "bin content to add = " << h_invMatrixOut->GetBinContent(binY,binX) << endl;
	  }
	  
	}
      }
    }
    
    //cout << "iz+1 = " << iz+1 << " bin content = " << zBinsVals[iz] << endl;
    h_zUnfNomBin->SetBinContent(iz+1,zBinsVals[iz]);
    
  }

  //compare the unfolded result from the one obtained with the matrix projection
  
  for(int iz = 0; iz < nBinZ_gen; iz++){
    cout << "bin = " << iz << endl;
    cout << "projected bin content = "<< h_zUnfNomBin->GetBinContent(iz+1) << endl;
    cout << "unfolded bin content = "<< h1_zUnf->GetBinContent(iz+1) << endl;
  }

  /*
  TCanvas *can_nominal = new TCanvas("can_nominal","can_nominal",1200,600);
  can_nominal->Divide(2,1);
  can_nominal->cd(1);
  h1_zUnf->GetXaxis()->SetTitle("unfolded result");
  h1_zUnf->Draw();
  can_nominal->cd(2);
  h_zUnfNomBin->GetXaxis()->SetTitle("projection from the matrix");
  h_zUnfNomBin->Draw("EP");
  can_nominal->SaveAs("../plots/nominalBinCheck.png");
  */
  
  //smear h_invMatrixOut -> ok
  TH2D *h_invMatrixOutSmear = (TH2D*)h_invMatrixOut->Clone("h_invMatrixOutSmear");
  h_invMatrixOutSmear->Reset();
    
  TRandom3 * myRand = new TRandom3();
  myRand->SetSeed(0);

  for(int jjet = 0; jjet < nBinJet_reco; jjet++){
    for(int jz = 0; jz < nBinZ_reco; jz++){
      int binY = nBinZ_reco*jjet+jz+1;
      
      float oldVal = 0;
      float errVal = 0;
      float smearVal = 0;

      for(int ijet = 0; ijet < nBinJet_gen; ijet++){
	for(int iz = 0; iz < nBinZ_gen; iz++){
	  int binX = nBinZ_gen*ijet+iz+1;

	  oldVal = h_invMatrixOut->GetBinContent(binY,binX);
	  errVal = h_invMatrixOut->GetBinError(binY,binX);
	  smearVal = myRand->Gaus(oldVal,errVal);
	  //cout <<"oldVal = "<<oldVal<<", errVal = "<<errVal<<", smearVal = "<<smearVal<<endl;
	  if(smearVal<0) {
	    smearVal = 0;
	  }
	  
	  h_invMatrixOutSmear->SetBinContent(binY,binX,smearVal);
	  
	}
      }

    }

  }

  //project the correct axis
  
  TH1D *h_zUnfNomBinSmear = (TH1D*)h1_zUnf->Clone("h_zUnfNomBinSmear");
  h_zUnfNomBinSmear->Reset();
  float zBinsValsSmear[nBinZ_gen];

  for(int iz = 0; iz < nBinZ_gen; iz++){
    zBinsValsSmear[iz] = 0;
    // int ijet = 1; // not 1 jet pt bin, but 5,6,7,8,9
    for(int ijet = 0; ijet < nBinJet_gen; ijet++){
      
      int binX = nBinZ_gen*ijet+iz+1;

      //      zBinsValsSmear[iz] = 0;
      
      for(int jjet = 0; jjet < nBinJet_reco; jjet++){
	for(int jz = 0; jz < nBinZ_reco; jz++){
	  int binY = nBinZ_reco*jjet+jz+1;
	  // do it only for the unf jt pt 25-35 GeV
	  int nStart = midLowerId;
	  int nEnd = midLowerId+(nBinJet_gen/nBinJet_reco);
	  if (ijet>=nStart && ijet<nEnd) zBinsValsSmear[iz] += h_invMatrixOutSmear->GetBinContent(binY,binX);
	}
      }
      
    }
    
    h_zUnfNomBinSmear->SetBinContent(iz+1,zBinsValsSmear[iz]);
    
  }
  

  h_zUnfNomBinSmear->Draw("EP");
    
  
  TFile *outfile = new TFile(outputfile.c_str(),toyNumber==1?"RECREATE":"UPDATE");
  if(toyNumber == 1) {
    h1_zUnf->Write("nominalZUnf");
    h_invMatrixOutSmear->Write("invMatrixOutSmear");
  }
  h_zUnfNomBinSmear->Write(Form("zUnfSmear_toy%i",toyNumber));
  outfile->ls();
  outfile->Close();
  
}

void smearTransferMatrix(){
  for(int i = 1; i < 101; i++){
    //operate(true,true,i);
    //if (!doCent && !doPeri)
    operate(true,false,i);
    //nonprompt
    //operate(false,true,i);
    //operate(false,false,i);
  }
}
