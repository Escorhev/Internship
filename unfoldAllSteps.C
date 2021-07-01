#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__) || !already
#include "inputParams.h"
#endif

#include "create2DMeas_data.cc"
#include "PlotRatios_DataUnfolded_afterDiag.cc"
#include "finalResults_statErrs.cc"
#include "prepareInputsNewPrNominal.cc"
#include "createRooUnfoldResponseNewPrNominal.cxx"
#include "unfoldStepNewPrNominal.cxx"
#include "unfoldStepNewPrNominalDiag.cxx"
#include "plotFinalResults.cc"
#include "smearTransferMatrix.cc" 
#include "trStatSystUnc.cc"
#include "makeRebinMatrix.C"
#include "createRooUnfoldResponseDiag.cxx"

/////////////////////change in the create2DMeas_data_systErrs.cc
void unfoldAllSteps(int SetCentShift = 0, int setJESsyst = 0, double setSF = 1.1, bool setMatrixInv = false, bool setSystTrMint = false, int nIt = 5, int nSI = 4, int nSI_pp = 4, int istart=1, int istop=4) {
  JESsyst = setJESsyst;
  SF = setSF;
  centShift = SetCentShift;
  matrixInv = setMatrixInv;
  useSystTrM = setSystTrMint;
  if (matrixInv) {nIt=1; nSI_pp=1; istart=1; istop=1;}
  nIter = nIt;
  nSIter = nSI;
  nSIter_pp = nSI_pp;
  
  for (int i=istart;i<=istop;i++) {
    std::cout<<"[INFO] SI "<<i<<" begins"<<endl;
    prepareInputsNewPrNominal(i);
    createRooUnfoldResponseNewPrNominal(i);
    unfoldStepNewPrNominal(i);
    unfoldStepNewPrNominalDiag(i);
    std::cout<<"[INFO] SI "<<i<<" done"<<endl;
  }  
  std::cout<<"[INFO] plotting done"<<endl;
  
  //PlotRatios_DataUnfolded_afterDiag();
  //if (!matrixInv && (istop>=nSI))
  finalResults_statErrs();
}

void unfoldAllSystematics(int nIt = 3, int nSI = 25, int nSI_pp = 1, int istart=1, int istop=27) {
  makeRebinMatrix();
  createRooUnfoldResponseDiag();
  create2DMeas_data(1,1);
  if (!doCent && !doPeri)
  create2DMeas_data(1,0);
  unfoldAllSteps(0, 0, 1.1, false, false, nIt, nSI, nSI_pp, istart, istop); //nominal
  if (istop>nSI) {
  smearTransferMatrix();
  trStatSystUnc();
  }
  systErr = true;
  create2DMeas_data(1,1);
  if (!doCent && !doPeri)
  create2DMeas_data(1,0);
  unfoldAllSteps(0, 0, 1.1, false, false, nIt, nSI, nSI_pp, istart, istop); //systError
  systErr = false;
  unfoldAllSteps(+1, 0, 1.1, false, false, nIt, nSI, nSI_pp, istart, istop); //cent up
  unfoldAllSteps(-1, 0, 1.1, false, false, nIt, nSI, nSI_pp, istart, istop); //cent down
  unfoldAllSteps(0, +1, 1.1, false, false, nIt, nSI, nSI_pp, istart, istop); //JES up
  unfoldAllSteps(0, -1, 1.1, false, false, nIt, nSI, nSI_pp, istart, istop); //JES down
  unfoldAllSteps(0, 0, 1.2, false, false, nIt, nSI, nSI_pp, istart, istop); //SF up
  unfoldAllSteps(0, 0, 1.0, false, false, nIt, nSI, nSI_pp, istart, istop); //SF down

  nprPrior = true;
  unfoldAllSteps(0, 0, 1.1, false, false, nIt, nSI, nSI_pp, istart, istop); //prior syst
  nprPrior = false;

  unfoldAllSteps(0, 0, 1.1, false, false, doCent?nIter_syst_cent:doPeri?nIter_syst_peri:nIter_syst,doCent?nSIter_syst_cent:doPeri?nSIter_syst_peri:nSIter_syst, nSIter_syst_pp, istart, istop+5); //regularisation 
  nIter = nIt;
  nSIter = nSI;
  nSIter_pp = nSI_pp;
  plotFinalResults();
}

void unfoldAllCent() {
  doCent=true;doPeri=false;
  unfoldAllSystematics(nIter_cent, nSIter_cent, 1, 1, 30);
}

void unfoldAllPeri() {
  doCent=false;doPeri=true;
  unfoldAllSystematics(nIter_peri, nSIter_peri, 1, 1, 22);
}

