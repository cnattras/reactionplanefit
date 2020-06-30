// A CorrData emulates 2 detectors A and B producing each
// a TClonesArray of Hit objects.
// A TClonesArray  of Track objects is built with Hits objects
// of detectors A and B. Eack Track object has a TRefArray of hits.
// A TClonesArray of Jets is made with a subset of the Track objects
// also stored in a TRefArray.
// see $ROOTSYS/tutorials/jets.C for an example creating a Tree
// with CorrDatas.
      
#include "TMath.h"   
#include "TRandom.h"   
#include "CorrData.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"
#include "TColor.h"
#include "Rtypes.h"
#include <fstream>
#include <iostream>
#include "TLatex.h"
Float_t Rn[5] = {1.0,1.0,1.0,1.0,1.0};

ClassImp(CorrData)

Int_t  CorrData::nBinsInRxnPlane = 6;
//the phi range used for the fit needs to be static because it is called by some of the functions and one cannot define a TF1 with a Double_t or Float_t which is not static.  Therefore one has to edit these by hand.  It is not elegant but it works.
Double_t CorrData::phiRange =TMath::Pi()/2.0;// 1.309;//TMath::Pi()/2.0;//1.0472;//TMath::Pi()/2.0;//1.0472;
Float_t  CorrData::mypi = TMath::Pi();
Float_t CorrData::sqrt3 = TMath::Sqrt(3.0);
//Likewise, the Rn have to be edited and input by hand because they are used by a TF1 and therefore need to be static.
// Double_t CorrData::R2 = 1.0;
// Double_t CorrData::R4 = 1.0;
// Double_t CorrData::R6 = 1.0;
// Double_t CorrData::R8 = 1.0;
// Double_t CorrData::R10 = 1.0;
//1.0-1.5 GeV
// Double_t CorrData::R2 = 0.637;
// Double_t CorrData::R4 = 0.247;
// Double_t CorrData::R6 = 0.104;
// Double_t CorrData::R8 = 0.104;
// Double_t CorrData::R10 = 0.0;
//1.5-2.0 GeV
// Double_t CorrData::R2 = 0.676;
// Double_t CorrData::R4 = 0.328;
// Double_t CorrData::R6 = 0.129;
// Double_t CorrData::R8 = 0;//0.129;
// Double_t CorrData::R10 = 0.0;
//2.0-3.0 GeV
// Double_t CorrData::R2 = 0.704;
// Double_t CorrData::R4 = 0.360;
// Double_t CorrData::R6 = 0.150;
// Double_t CorrData::R8 = 0.150;
// Double_t CorrData::R10 = 0.0;
//3.0-4.0 GeV
Double_t CorrData::R2 = 0.704;
Double_t CorrData::R4 = 0.360;
Double_t CorrData::R6 = 0.150;
Double_t CorrData::R8 = 0.150;
Double_t CorrData::R10 = 0.0;
Double_t CorrData::R2FracErr = 0.01;
Double_t CorrData::R4FracErr = 0.03;
Double_t CorrData::R6FracErr = 0.03;
// Double_t CorrData::R2FracErr = 0.01;
// Double_t CorrData::R4FracErr = 0.03;
// Double_t CorrData::R6FracErr = 0.03;
//Histogram ranges - these generally should not be changed.  The input data do not run over the same ranges and therefore are shifted.  Care should be used if and when shifting these.  The phi range used in the fit needs to fit on the negative side here and the range needs to be 2Pi
Double_t CorrData::MINPHI = -TMath::Pi()/2.0;
Double_t CorrData::MAXPHI = TMath::Pi()*3.0/2.0;
//______________________________________________________________________________
CorrData::CorrData()
{
  for(Int_t i= 0;i<nBinsInRxnPlane;i++){
    hSignalPlusBkgd[i] = NULL;
    hSignal[i] = NULL;
    hBkgd[i] = NULL;
    hCalcBkgd[i] = NULL;
    hCalcBkgdInSignalPlusBackgroundRegion[i] = NULL;
    fBkgd[i] = NULL;
    fBkgdInSignalPlusBkgdRegion[i] = NULL;
    hSignalLowScale[i] = NULL;
    hCalcBkgdInSignalPlusBackgroundRegionLowScale[i] = NULL;
    fBkgdInSignalPlusBkgdRegionLowScale[i] = NULL;
    hSignalHighScale[i] = NULL;
    hCalcBkgdInSignalPlusBackgroundRegionHighScale[i] = NULL;
    fBkgdInSignalPlusBkgdRegionHighScale[i] = NULL;
    hSignalLowR[i] = NULL;
    hCalcBkgdLowR[i] = NULL;
    hCalcBkgdInSignalPlusBackgroundRegionLowR[i] = NULL;
    fBkgdLowR[i] = NULL;
    fBkgdInSignalPlusBkgdRegionLowR[i] = NULL;
    hSignalHighR[i] = NULL;
    hCalcBkgdHighR[i] = NULL;
    hCalcBkgdInSignalPlusBackgroundRegionHighR[i] = NULL;
    fBkgdHighR[i] = NULL;
    fBkgdInSignalPlusBkgdRegionHighR[i] = NULL;
    nsYield[i] = 0;
    asYield[i] = 0;
    nsRMS[i] = 0;
    asRMS[i] = 0;
    nsYieldError[i] = 0;
    asYieldError[i] = 0;
    nsRMSError[i] = 0;
    asRMSError[i] = 0;
    nsYieldLowScale[i] = 0;
    asYieldLowScale[i] = 0;
    nsRMSLowScale[i] = 0;
    asRMSLowScale[i] = 0;
    nsYieldErrorLowScale[i] = 0;
    asYieldErrorLowScale[i] = 0;
    nsRMSErrorLowScale[i] = 0;
    asRMSErrorLowScale[i] = 0;
    nsYieldHighScale[i] = 0;
    asYieldHighScale[i] = 0;
    nsRMSHighScale[i] = 0;
    asRMSHighScale[i] = 0;
    nsYieldErrorHighScale[i] = 0;
    asYieldErrorHighScale[i] = 0;
    nsRMSErrorHighScale[i] = 0;
    asRMSErrorHighScale[i] = 0;
    nsYieldLowR[i] = 0;
    asYieldLowR[i] = 0;
    nsRMSLowR[i] = 0;
    asRMSLowR[i] = 0;
    nsYieldErrorLowR[i] = 0;
    asYieldErrorLowR[i] = 0;
    nsRMSErrorLowR[i] = 0;
    asRMSErrorLowR[i] = 0;
    nsYieldHighR[i] = 0;
    asYieldHighR[i] = 0;
    nsRMSHighR[i] = 0;
    asRMSHighR[i] = 0;
    nsYieldErrorHighR[i] = 0;
    asYieldErrorHighR[i] = 0;
    nsRMSErrorHighR[i] = 0;
    asRMSErrorHighR[i] = 0;
  }
  fitFunc = NULL;
  fitFuncLowR = NULL;
  fitFuncHighR = NULL;
  hCommonBkgd = NULL;
  etaRangeSigToBkgdRatio = 1.0/0.4225;
  etaRangeSigToBkgdRatioLow = 1.0/0.4225;
  etaRangeSigToBkgdRatioHigh = 1.0/0.4225;
  cSignalPlusBackground = NULL;
  cBackground = NULL;
  nsPhiCut = 0.7;
  asPhiCut = 0.7;
  nParameters = 6;
}

//______________________________________________________________________________
CorrData::~CorrData()
{
  for(Int_t i= 0;i<nBinsInRxnPlane;i++){
    delete hSignalPlusBkgd[i];
    delete hSignal[i];
    delete hBkgd[i];
    delete hCalcBkgd[i];
    delete hCalcBkgdInSignalPlusBackgroundRegion[i];
    delete fBkgd[i];
    delete fBkgdInSignalPlusBkgdRegion[i];
    delete hSignalLowR[i];
    delete hCalcBkgdLowR[i];
    delete hCalcBkgdInSignalPlusBackgroundRegionLowR[i];
    delete fBkgdLowR[i];
    delete fBkgdInSignalPlusBkgdRegionLowR[i];
    delete hSignalHighR[i];
    delete hCalcBkgdHighR[i];
    delete hCalcBkgdInSignalPlusBackgroundRegionHighR[i];
    delete fBkgdHighR[i];
    delete fBkgdInSignalPlusBkgdRegionHighR[i];
  }
  delete fitFunc;
  delete fitFuncLowR;
  delete fitFuncHighR;
  delete hCommonBkgd;
  delete cSignalPlusBackground;
  delete cBackground;
}

//______________________________________________________________________________
TH1F *CorrData::GetSignalPlusBackgroundHistogram(Int_t bin) {
  if(bin>=0 && bin<nBinsInRxnPlane){
    return hSignalPlusBkgd[bin];
  }
  else{
    cerr<<"Warning:  Requested histogram "<<bin<<" and there are only "<<nBinsInRxnPlane<<" histograms"<<endl;
    return NULL;
  }
}


TH1F *CorrData::GetBackgroundHistogram(Int_t bin) {
  if(bin>=0 && bin<nBinsInRxnPlane){
    return hBkgd[bin];
  }
  else{
    cerr<<"Warning:  Requested histogram "<<bin<<" and there are only "<<nBinsInRxnPlane<<" histograms"<<endl;
    return NULL;
  }
}


TH1F *CorrData::GetSignalHistogram(Int_t bin) {
  if(bin>=0 && bin<nBinsInRxnPlane){
    return hSignal[bin];
  }
  else{
    cerr<<"Warning:  Requested histogram "<<bin<<" and there are only "<<nBinsInRxnPlane<<" histograms"<<endl;
    return NULL;
  }
}


Double_t CorrData::R(Int_t n) {
  switch(n){
  case 2:
    return R2;
  case 4:
    return R4;
  case 6:
    return R6;
  case 8:
    return R8;
  case 10:
    return R10;
  }
  return -1.0;
}

void CorrData::SetR(Int_t n, Float_t R) {
  switch(n){
  case 2:
    Rn[0] = R;
    break;
  case 4:
    Rn[1] = R;
    break;
  case 6:
    Rn[2] = R;
    break;
  case 8:
    Rn[3] = R;
    break;
  case 10:
    Rn[4] = R;
    break;
  }
  return;
}

void CorrData::SetR(Float_t thisR2, Float_t thisR4, Float_t thisR6, Float_t thisR8, Float_t thisR10){
    Rn[0] = thisR2;
    Rn[1] = thisR4;
    Rn[2] = thisR6;
    Rn[3] = thisR8;
    Rn[4] = thisR10;
}
void CorrData::CreateCommonBackgroundHisto(char *name){
  Int_t phiLowBin = hBkgd[0]->FindBin(-phiRange+1e-3);//Set the phi range to the actual range usable in the histogram
  Int_t phiHighBin = hBkgd[0]->FindBin(phiRange-1e-3);//Set the phi range to the actual range usable in the histogram
  Float_t realphiRange = hBkgd[0]->GetBinLowEdge(phiHighBin+1);
  cout<<"Using phi range "<<-realphiRange<<" to "<<realphiRange<<" fixed phi range is "<<phiRange<<" warning if these are not equal you will have bugs!!"<<endl;
  nbinsBkgd = (phiHighBin+1-phiLowBin)*nBinsInRxnPlane;
  phiRangeBkgd = realphiRange*nBinsInRxnPlane*2.0;
  //cout<<"nbins "<<nbinsBkgd<<" total range "<<phiRangeBkgd<<endl;
  hCommonBkgd = new TH1F(name,name,nbinsBkgd,0.0,phiRangeBkgd);
  for(int i=0;i<nBinsInRxnPlane;i++){
    for(int j=phiLowBin;j<=phiHighBin;j++){
      Float_t phi = ConvertPhiToCommonPhi(i,hBkgd[i]->GetBinCenter(j));
      Float_t content = hBkgd[i]->GetBinContent(j);
      Float_t error = hBkgd[i]->GetBinError(j);
      Int_t newbin = hCommonBkgd->FindBin(phi);
      //cout<<"old bin i "<<i<<", j "<<j<<" "<<content<<" +/- "<<error;
      hCommonBkgd->SetBinContent(newbin,content);
      hCommonBkgd->SetBinError(newbin,error);
      //cout<<" setting bin "<<newbin<<" to "<<content<<" +/- "<<error<<endl;
    }
  }
}
Double_t CorrData::ConvertPhiToCommonPhi(Int_t bin, Double_t phi){
  if(TMath::Abs(phi)>phiRange){
    cerr<<"phi "<<phi<<" out of range!"<<endl;
    return -1;
  }
  //cout<<"changing "<<phi<<" in bin "<<bin<<" to "<<phiRange*(bin*2.0+1)+phi<<endl;
  return phiRange*(bin*2.0+1)+phi;
}
Double_t CorrData::ConvertCommonPhiToPhi(Int_t bin, Double_t phi){
  Double_t newphi = phi - phiRange*(bin*2.0+1);
  return newphi;
  if(TMath::Abs(phi)>phiRange){
    cerr<<"phi "<<phi<<" out of range!"<<endl;
    return -1;
  }
  //cout<<"changing "<<phi<<" in bin "<<bin<<" to "<<phiRange*(bin*2.0+1)+phi<<endl;
  return phiRange*(bin*2.0+1)+phi;
}
Int_t CorrData::CalculateReactionPlaneBin(Float_t commonphi){
  return (Int_t) (commonphi/2.0/phiRange);
}
Double_t CorrData::vntrig(Int_t n, Int_t rxnplanebin, Double_t v2, Double_t v4, Int_t rangeoption){
  Double_t R2eff = R2;
  Double_t R4eff = R4;
  Double_t R6eff = R6;
  Double_t R8eff = R8;
  Double_t R10eff = R10;
  if(rangeoption==-1){
    R2eff = R2*(1-R2FracErr);
    R4eff = R4*(1-R4FracErr);
    R6eff = R6*(1-R6FracErr);
    R8eff = 0;
    R10eff = 0;
  }
  if(rangeoption==1){
    R2eff = R2*(1+R2FracErr);
    R4eff = R4*(1+R4FracErr);
    R6eff = R6*(1+R6FracErr);
    R8eff = R6eff;
    R10eff = R6eff;
  }
  if(n==2){
    switch(rxnplanebin){
    case 0:
      return (2.0*mypi*v2+3.0*sqrt3*R4eff*v2+4.0*R6eff*v4+6.0*R2eff*(1.0+v4))/2.0/(mypi+6.0*R2eff*v2+3*sqrt3*R4eff*v4);
      break;
    case 1:
      return (mypi*v2-2.0*R6eff*v4+3.0*(sqrt3-1)*R2eff*(1.0+v4))/(mypi+6.0*(-1.0+sqrt3)*R2eff*v2);
      break;
    case 2:
      return (2.0*mypi*v2-3.0*sqrt3*R4eff*v2-4.0*R6eff*v4-6.0*(-2.0+sqrt3)*R2eff*(1.0+v4))/2.0/(mypi-3.0*(2.0*(-2.0+sqrt3)*R2eff*v2+sqrt3*R4eff*v4));
    case 3:
      return (2.0*mypi*v2-3.0*sqrt3*R4eff*v2+4.0*R6eff*v4+6.0*(-2.0+sqrt3)*R2eff*(1.0+v4))/2.0/(mypi+6.0*(-2.0+sqrt3)*R2eff*v2-3.0*sqrt3*R4eff*v4);
      //case 2:
    case 4:
      return (mypi*v2+2.0*R6eff*v4-3.0*(-1.0+sqrt3)*R2eff*(1.0+v4))/(mypi-6.0*(-1.0+sqrt3)*R2eff*v2);
    case 5:
      return (2.0*mypi*v2+3.0*sqrt3*R4eff*v2-4.0*R6eff*v4-6.0*R2eff*(1+v4))/2.0/(mypi-6.0*R2eff*v2+3.0*sqrt3*R4eff*v4);
    }
  }
  if(n==4){
    switch(rxnplanebin){
    case 0:
      return (6.0*sqrt3*R4eff+12.0*R2eff*v2+8.0*R6eff*v2+4.0*mypi*v4+3.0*sqrt3*R8eff*v4)/4.0/(mypi+6.0*R2eff*v2+3.0*sqrt3*R4eff*v4);
      break;
    case 1:
      return (6.0*(-1.0+sqrt3)*R2eff*v2-4.0*R6eff*v2+2.0*mypi*v4-3.0*sqrt3*R8eff*v4)/2.0/(mypi+6.0*(-1.0+sqrt3)*R2eff*v2);
      break;
    case 2:
      return (-6.0*sqrt3*R4eff-12.0*(-2+sqrt3)*R2eff*v2-8.0*R6eff*v2+4*mypi*v4+3.0*sqrt3*R8eff*v4)/4.0/(mypi-3.0*(2.0*(-2.0+sqrt3)*R2eff*v2+sqrt3*R4eff*v4));
    case 3:
      return (-6.0*sqrt3*R4eff+12.0*(-2.0+sqrt3)*R2eff*v2+8.0*R6eff*v2+4.0*mypi*v4+3*sqrt3*R8eff*v4)/4.0/(mypi+6.0*(-2.0+sqrt3)*R2eff*v2-3.0*sqrt3*R4eff*v4);
    case 4:
      return (-6.0*(-1.0+sqrt3)*R2eff*v2+4.0*R6eff*v2+2.0*mypi*v4-3.0*sqrt3*R8eff*v4)/2.0/(mypi-6.0*(-1.0+sqrt3)*R2eff*v2);
    case 5:
      return (6.0*sqrt3*R4eff-12.0*R2eff*v2-8.0*R6eff*v2+4.0*mypi*v4+3.0*sqrt3*R8eff*v4)/4.0/(mypi-6.0*R2eff*v2+3.0*sqrt3*R4eff*v4);
    }
  }
  return 0.0;
}

Double_t CorrData::Background(const double * x, const double * p, Int_t rxnplanebin){
  //reordering parameters allows us to fix the jet vn by just choosing a different number of parameters :)
  //this option isn't implemented yet but we might want to.  The jet vn would have to be hard coded.
  Double_t B = p[0];
  Double_t v2assoc = p[1];
  Double_t v3 = p[2];
  Double_t v4assoc = p[3];
  Double_t v2jet = p[4];
  Double_t v4jet = p[5]; 
  Double_t v5 = p[6];
  //Double_t phi = x[0];
  Double_t phi = ConvertCommonPhiToPhi(rxnplanebin,x[0]);
  Double_t Part2  = 2.0*vntrig(2,rxnplanebin,v2jet,v4jet,0) *v2assoc*TMath::Cos(2.0*phi);
  Double_t Part3 = 2.0*v3*TMath::Cos(3*phi);
  Double_t Part4  = 2.0*vntrig(4,rxnplanebin,v2jet,v4jet,0) *v4assoc*TMath::Cos(4.0*phi);
  Double_t Part5 = 2.0*v5*TMath::Cos(5*phi);
  return B*(1+Part2+Part3+Part4+Part5);
}
Double_t CorrData::BackgroundGeneral(const double * x, const double * p, Int_t rxnplanebin){
  //reordering parameters allows us to fix the jet vn by just choosing a different number of parameters :)
  //this option isn't implemented yet but we might want to.  The jet vn would have to be hard coded.
  Double_t B = p[0];
  Double_t v2assoc = p[1];
  Double_t v3 = p[2];
  Double_t v4assoc = p[3];
  Double_t v2jet = p[4];
  Double_t v4jet = p[5]; 
  Double_t v5 = p[6];
  Double_t phi = x[0];
  //Double_t phi = ConvertCommonPhiToPhi(rxnplanebin,x[0]);
  Double_t Part2  = 2.0*vntrig(2,rxnplanebin,v2jet,v4jet,0) *v2assoc*TMath::Cos(2.0*phi);
  Double_t Part3 = 2.0*v3*TMath::Cos(3*phi);
  Double_t Part4  = 2.0*vntrig(4,rxnplanebin,v2jet,v4jet,0) *v4assoc*TMath::Cos(4.0*phi);
  Double_t Part5 = 2.0*v5*TMath::Cos(5*phi);
  return B*(1+Part2+Part3+Part4+Part5);
}
Double_t CorrData::BackgroundCombined(const double * x, const double * p){
  Int_t bin = (Int_t) (x[0]/2.0/phiRange);
  //x[0] = ConvertCommonPhiToPhi(bin,x[0]);
  return Background(x,p,bin);
}
Double_t CorrData::BackgroundLowR(const double * x, const double * p, Int_t rxnplanebin){
  //reordering parameters allows us to fix the jet vn by just choosing a different number of parameters :)
  //this option isn't implemented yet but we might want to.  The jet vn would have to be hard coded.
  Double_t B = p[0];
  Double_t v2assoc = p[1];
  Double_t v3 = p[2];
  Double_t v4assoc = p[3];
  Double_t v2jet = p[4];
  Double_t v4jet = p[5]; 
  Double_t v5 = p[6];
  //Double_t phi = x[0];
  Double_t phi = ConvertCommonPhiToPhi(rxnplanebin,x[0]);
  Double_t Part2  = 2.0*vntrig(2,rxnplanebin,v2jet,v4jet,-1) *v2assoc*TMath::Cos(2.0*phi);
  Double_t Part3 = 2.0*v3*TMath::Cos(3*phi);
  Double_t Part4  = 2.0*vntrig(4,rxnplanebin,v2jet,v4jet,-1) *v4assoc*TMath::Cos(4.0*phi);
  Double_t Part5 = 2.0*v5*TMath::Cos(5*phi);
  return B*(1+Part2+Part3+Part4+Part5);
}
Double_t CorrData::BackgroundGeneralLowR(const double * x, const double * p, Int_t rxnplanebin){
  //reordering parameters allows us to fix the jet vn by just choosing a different number of parameters :)
  //this option isn't implemented yet but we might want to.  The jet vn would have to be hard coded.
  Double_t B = p[0];
  Double_t v2assoc = p[1];
  Double_t v3 = p[2];
  Double_t v4assoc = p[3];
  Double_t v2jet = p[4];
  Double_t v4jet = p[5]; 
  Double_t v5 = p[6];
  Double_t phi = x[0];
  //Double_t phi = ConvertCommonPhiToPhi(rxnplanebin,x[0]);
  Double_t Part2  = 2.0*vntrig(2,rxnplanebin,v2jet,v4jet,-1) *v2assoc*TMath::Cos(2.0*phi);
  Double_t Part3 = 2.0*v3*TMath::Cos(3*phi);
  Double_t Part4  = 2.0*vntrig(4,rxnplanebin,v2jet,v4jet,-1) *v4assoc*TMath::Cos(4.0*phi);
  Double_t Part5 = 2.0*v5*TMath::Cos(5*phi);
  return B*(1+Part2+Part3+Part4+Part5);
}
Double_t CorrData::BackgroundCombinedLowR(const double * x, const double * p){
  Int_t bin = (Int_t) (x[0]/2.0/phiRange);
  //x[0] = ConvertCommonPhiToPhi(bin,x[0]);
  return BackgroundLowR(x,p,bin);
}

Double_t CorrData::BackgroundHighR(const double * x, const double * p, Int_t rxnplanebin){
  //reordering parameters allows us to fix the jet vn by just choosing a different number of parameters :)
  //this option isn't implemented yet but we might want to.  The jet vn would have to be hard coded.
  Double_t B = p[0];
  Double_t v2assoc = p[1];
  Double_t v3 = p[2];
  Double_t v4assoc = p[3];
  Double_t v2jet = p[4];
  Double_t v4jet = p[5]; 
  Double_t v5 = p[6];
  //Double_t phi = x[0];
  Double_t phi = ConvertCommonPhiToPhi(rxnplanebin,x[0]);
  Double_t Part2  = 2.0*vntrig(2,rxnplanebin,v2jet,v4jet,1) *v2assoc*TMath::Cos(2.0*phi);
  Double_t Part3 = 2.0*v3*TMath::Cos(3*phi);
  Double_t Part4  = 2.0*vntrig(4,rxnplanebin,v2jet,v4jet,1) *v4assoc*TMath::Cos(4.0*phi);
  Double_t Part5 = 2.0*v5*TMath::Cos(5*phi);
  return B*(1+Part2+Part3+Part4+Part5);
}
Double_t CorrData::BackgroundGeneralHighR(const double * x, const double * p, Int_t rxnplanebin){
  //reordering parameters allows us to fix the jet vn by just choosing a different number of parameters :)
  //this option isn't implemented yet but we might want to.  The jet vn would have to be hard coded.
  Double_t B = p[0];
  Double_t v2assoc = p[1];
  Double_t v3 = p[2];
  Double_t v4assoc = p[3];
  Double_t v2jet = p[4];
  Double_t v4jet = p[5]; 
  Double_t v5 = p[6];
  Double_t phi = x[0];
  //Double_t phi = ConvertCommonPhiToPhi(rxnplanebin,x[0]);
  Double_t Part2  = 2.0*vntrig(2,rxnplanebin,v2jet,v4jet,1) *v2assoc*TMath::Cos(2.0*phi);
  Double_t Part3 = 2.0*v3*TMath::Cos(3*phi);
  Double_t Part4  = 2.0*vntrig(4,rxnplanebin,v2jet,v4jet,1) *v4assoc*TMath::Cos(4.0*phi);
  Double_t Part5 = 2.0*v5*TMath::Cos(5*phi);
  return B*(1+Part2+Part3+Part4+Part5);
}
Double_t CorrData::BackgroundCombinedHighR(const double * x, const double * p){
  Int_t bin = (Int_t) (x[0]/2.0/phiRange);
  //x[0] = ConvertCommonPhiToPhi(bin,x[0]);
  return BackgroundHighR(x,p,bin);
}

void CorrData::InitializeCommonFitFunction(){
  //fitFunc = func;
   Double_t min = 0.0;
   Double_t max = phiRange*2.0*nBinsInRxnPlane;
  fitFunc = new TF1("func",BackgroundCombined, min,max , nParameters);
  fitFunc->SetParNames("B","v2assoc","v3","v4assoc","v2jet","v4jet","v5");
  Float_t B = (hCommonBkgd->GetMaximum() + hCommonBkgd->GetMinimum())/2.0;
  //cout<<"Setting B at "<<B<<endl;
  fitFunc->SetParameters(B,0.1,0.1*0.1,0.01,0.0,0.0);
  fitFuncLowR = new TF1("funcLowR",BackgroundCombinedLowR, min,max , nParameters);
  fitFuncLowR->SetParNames("B","v2assoc","v3","v4assoc","v2jetLowR","v4jetLowR","v5");
  fitFuncLowR->SetParameters(B,0.1,0.10*0.10,0.010,0.0,0.0);
  fitFuncHighR = new TF1("funcHighR",BackgroundCombinedHighR, min,max , nParameters);
  fitFuncHighR->SetParNames("B","v2assoc","v3","v4assoc","v2jetHighR","v4jetHighR","v5");
  fitFuncHighR->SetParameters(B,0.1,0.1*0.1,0.01,0.0,0.0);
//   fitFunc->SetParLimits(3,0,0.1);
//   fitFunc->SetParLimits(4,0,1);
//   fitFunc->SetParLimits(5,0,1e-1);
//   fitFuncHighR->SetParLimits(4,0,1);
//   fitFuncHighR->SetParLimits(5,0,1e-1);
//   fitFuncLowR->SetParLimits(4,0,1);
//   fitFuncLowR->SetParLimits(5,0,1e-1);

  fBkgd[0] = new TF1("fBkgd0",Background0, MINPHI,MAXPHI , nParameters);
  fBkgd[1] = new TF1("fBkgd1",Background1, MINPHI,MAXPHI , nParameters);
  fBkgd[2] = new TF1("fBkgd2",Background2, MINPHI,MAXPHI , nParameters);
  fBkgd[3] = new TF1("fBkgd3",Background3, MINPHI,MAXPHI , nParameters);
  fBkgd[4] = new TF1("fBkgd4",Background4, MINPHI,MAXPHI , nParameters);
  fBkgd[5] = new TF1("fBkgd5",Background5, MINPHI,MAXPHI , nParameters);
  fBkgdInSignalPlusBkgdRegion[0] = new TF1("fBkgdInSignalPlusBkgdRegion0",Background0, MINPHI,MAXPHI , nParameters);
  fBkgdInSignalPlusBkgdRegion[1] = new TF1("fBkgdInSignalPlusBkgdRegion1",Background1, MINPHI,MAXPHI , nParameters);
  fBkgdInSignalPlusBkgdRegion[2] = new TF1("fBkgdInSignalPlusBkgdRegion2",Background2, MINPHI,MAXPHI , nParameters);
  fBkgdInSignalPlusBkgdRegion[3] = new TF1("fBkgdInSignalPlusBkgdRegion3",Background3, MINPHI,MAXPHI , nParameters);
  fBkgdInSignalPlusBkgdRegion[4] = new TF1("fBkgdInSignalPlusBkgdRegion4",Background4, MINPHI,MAXPHI , nParameters);
  fBkgdInSignalPlusBkgdRegion[5] = new TF1("fBkgdInSignalPlusBkgdRegion5",Background5, MINPHI,MAXPHI , nParameters);
  fBkgdInSignalPlusBkgdRegionLowScale[0] = new TF1("fBkgdInSignalPlusBkgdRegion0LowScale",Background0, MINPHI,MAXPHI , nParameters);
  fBkgdInSignalPlusBkgdRegionLowScale[1] = new TF1("fBkgdInSignalPlusBkgdRegion1LowScale",Background1, MINPHI,MAXPHI , nParameters);
  fBkgdInSignalPlusBkgdRegionLowScale[2] = new TF1("fBkgdInSignalPlusBkgdRegion2LowScale",Background2, MINPHI,MAXPHI , nParameters);
  fBkgdInSignalPlusBkgdRegionLowScale[3] = new TF1("fBkgdInSignalPlusBkgdRegion3LowScale",Background3, MINPHI,MAXPHI , nParameters);
  fBkgdInSignalPlusBkgdRegionLowScale[4] = new TF1("fBkgdInSignalPlusBkgdRegion4LowScale",Background4, MINPHI,MAXPHI , nParameters);
  fBkgdInSignalPlusBkgdRegionLowScale[5] = new TF1("fBkgdInSignalPlusBkgdRegion5LowScale",Background5, MINPHI,MAXPHI , nParameters);
  fBkgdInSignalPlusBkgdRegionHighScale[0] = new TF1("fBkgdInSignalPlusBkgdRegion0HighScale",Background0, MINPHI,MAXPHI , nParameters);
  fBkgdInSignalPlusBkgdRegionHighScale[1] = new TF1("fBkgdInSignalPlusBkgdRegion1HighScale",Background1, MINPHI,MAXPHI , nParameters);
  fBkgdInSignalPlusBkgdRegionHighScale[2] = new TF1("fBkgdInSignalPlusBkgdRegion2HighScale",Background2, MINPHI,MAXPHI , nParameters);
  fBkgdInSignalPlusBkgdRegionHighScale[3] = new TF1("fBkgdInSignalPlusBkgdRegion3HighScale",Background3, MINPHI,MAXPHI , nParameters);
  fBkgdInSignalPlusBkgdRegionHighScale[4] = new TF1("fBkgdInSignalPlusBkgdRegion4HighScale",Background4, MINPHI,MAXPHI , nParameters);
  fBkgdInSignalPlusBkgdRegionHighScale[5] = new TF1("fBkgdInSignalPlusBkgdRegion5HighScale",Background5, MINPHI,MAXPHI , nParameters);
  fBkgdLowR[0] = new TF1("fBkgd0LowR",Background0LowR, MINPHI,MAXPHI , nParameters);
  fBkgdLowR[1] = new TF1("fBkgd1LowR",Background1LowR, MINPHI,MAXPHI , nParameters);
  fBkgdLowR[2] = new TF1("fBkgd2LowR",Background2LowR, MINPHI,MAXPHI , nParameters);
  fBkgdLowR[3] = new TF1("fBkgd3LowR",Background3LowR, MINPHI,MAXPHI , nParameters);
  fBkgdLowR[4] = new TF1("fBkgd4LowR",Background4LowR, MINPHI,MAXPHI , nParameters);
  fBkgdLowR[5] = new TF1("fBkgd5LowR",Background5LowR, MINPHI,MAXPHI , nParameters);
  fBkgdInSignalPlusBkgdRegionLowR[0] = new TF1("fBkgdInSignalPlusBkgdRegion0LowR",Background0LowR, MINPHI,MAXPHI , nParameters);
  fBkgdInSignalPlusBkgdRegionLowR[1] = new TF1("fBkgdInSignalPlusBkgdRegion1LowR",Background1LowR, MINPHI,MAXPHI , nParameters);
  fBkgdInSignalPlusBkgdRegionLowR[2] = new TF1("fBkgdInSignalPlusBkgdRegion2LowR",Background2LowR, MINPHI,MAXPHI , nParameters);
  fBkgdInSignalPlusBkgdRegionLowR[3] = new TF1("fBkgdInSignalPlusBkgdRegion3LowR",Background3LowR, MINPHI,MAXPHI , nParameters);
  fBkgdInSignalPlusBkgdRegionLowR[4] = new TF1("fBkgdInSignalPlusBkgdRegion4LowR",Background4LowR, MINPHI,MAXPHI , nParameters);
  fBkgdInSignalPlusBkgdRegionLowR[5] = new TF1("fBkgdInSignalPlusBkgdRegion5LowR",Background5LowR, MINPHI,MAXPHI , nParameters);
  fBkgdHighR[0] = new TF1("fBkgd0HighR",Background0HighR, MINPHI,MAXPHI , nParameters);
  fBkgdHighR[1] = new TF1("fBkgd1HighR",Background1HighR, MINPHI,MAXPHI , nParameters);
  fBkgdHighR[2] = new TF1("fBkgd2HighR",Background2HighR, MINPHI,MAXPHI , nParameters);
  fBkgdHighR[3] = new TF1("fBkgd3HighR",Background3HighR, MINPHI,MAXPHI , nParameters);
  fBkgdHighR[4] = new TF1("fBkgd4HighR",Background4HighR, MINPHI,MAXPHI , nParameters);
  fBkgdHighR[5] = new TF1("fBkgd5HighR",Background5HighR, MINPHI,MAXPHI , nParameters);
  fBkgdInSignalPlusBkgdRegionHighR[0] = new TF1("fBkgdInSignalPlusBkgdRegion0HighR",Background0HighR, MINPHI,MAXPHI , nParameters);
  fBkgdInSignalPlusBkgdRegionHighR[1] = new TF1("fBkgdInSignalPlusBkgdRegion1HighR",Background1HighR, MINPHI,MAXPHI , nParameters);
  fBkgdInSignalPlusBkgdRegionHighR[2] = new TF1("fBkgdInSignalPlusBkgdRegion2HighR",Background2HighR, MINPHI,MAXPHI , nParameters);
  fBkgdInSignalPlusBkgdRegionHighR[3] = new TF1("fBkgdInSignalPlusBkgdRegion3HighR",Background3HighR, MINPHI,MAXPHI , nParameters);
  fBkgdInSignalPlusBkgdRegionHighR[4] = new TF1("fBkgdInSignalPlusBkgdRegion4HighR",Background4HighR, MINPHI,MAXPHI , nParameters);
  fBkgdInSignalPlusBkgdRegionHighR[5] = new TF1("fBkgdInSignalPlusBkgdRegion5HighR",Background5HighR, MINPHI,MAXPHI , nParameters);
  for(Int_t i=0;i<nBinsInRxnPlane;i++){
    fBkgd[i]->SetParNames("B","v2assoc","v3","v4assoc","v2jet","v4jet");
    fBkgdInSignalPlusBkgdRegion[i]->SetParNames("B","v2assoc","v3","v4assoc","v2jet","v4jet");
    fBkgdLowR[i]->SetParNames("B","v2assoc","v3","v4assoc","v2jetLowR","v4jetLowR");
    fBkgdInSignalPlusBkgdRegionLowR[i]->SetParNames("B","v2assoc","v3","v4assoc","v2jetLowR","v4jetLowR");
    fBkgdHighR[i]->SetParNames("B","v2assoc","v3","v4assoc","v2jetHighR","v4jetHighR");
    fBkgdInSignalPlusBkgdRegionHighR[i]->SetParNames("B","v2assoc","v3","v4assoc","v2jetHighR","v4jetHighR");
  }
  return;
}
void CorrData::SetSignalHistogram(Int_t bin, TH1F *histo){
  delete hSignal[bin];
  hSignal[bin] = histo;
  delete hSignalLowR[bin];
  TString lowR = "LowR";
  TString nameLowR = histo->GetName() + lowR;
  hSignalLowR[bin] = (TH1F*) histo->Clone(nameLowR.Data());
  TString highR = "HighR";
  TString nameHighR = histo->GetName() + highR;
  hSignalHighR[bin] = (TH1F*) histo->Clone(nameHighR.Data());//Low
  TString lowScale = "LowScale";
  TString nameLowScale = histo->GetName() + lowScale;
  hSignalLowScale[bin] = (TH1F*) histo->Clone(nameLowScale.Data());
  TString highScale = "HighScale";
  TString nameHighScale = histo->GetName() + highScale;
  hSignalHighScale[bin] = (TH1F*) histo->Clone(nameHighScale.Data());//Low
return;
}
void CorrData::SetBackgroundHistogram(Int_t bin, TH1F *histo){
  delete hBkgd[bin];
  hBkgd[bin] = histo;
  return;
}

void CorrData::FitCommonBackgroundHistogram(){
  //we need to do the fits in succession because we get the error bars by getting the current TVirtualFitter
  hCommonBkgd->Fit(fitFunc);
  cerr<<"done with fit 1"<<endl;
  cout<<"Chi2 "<<fitFunc->GetChisquare()<<" NDF "<<fitFunc->GetNDF()<<" chi^2/ndf "<< fitFunc->GetChisquare() / fitFunc->GetNDF() <<endl;//hCommonBkgd
  cout<<" & ";
  cout<<Form("%2.2f & ",fitFunc->GetChisquare() / fitFunc->GetNDF());
  cout<<Form("%2.1f $\\pm$ %2.1f & ",fitFunc->GetParameter(0)*1.0, fitFunc->GetParError(0)*1.0);
  cout<<Form("& %2.1f $\\pm$ %2.1f & ",fitFunc->GetParameter(1)*100.0, fitFunc->GetParError(1)*100.0);
  cout<<Form("%2.1f $\\pm$ %2.1f & ",fitFunc->GetParameter(4)*100.0, fitFunc->GetParError(4)*100.0);
  cout<<Form("%2.0f $\\pm$ %2.0f & ",fitFunc->GetParameter(2)*10000.0, fitFunc->GetParError(2)*10000.0);
  cout<<Form("%2.1f $\\pm$ %2.1f & ",fitFunc->GetParameter(3)*100.0, fitFunc->GetParError(3)*100.0);
  cout<<Form("%2.2f $\\pm$ %2.2f \\\\ ",fitFunc->GetParameter(5)*100.0, fitFunc->GetParError(5)*100.0)<<endl;
  //the default is to give 95% confidence intervals but we want the standard 1 sigma (68%) confidence intervals
  cerr<<"Getting confidence intervals fit 1"<<endl;
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hCommonBkgd, 0.68);
  cerr<<"done with confidence intervals fit 1"<<endl;
  for(Int_t i=0;i<nBinsInRxnPlane;i++){
    for(Int_t j=0;j<6;j++){
      fBkgd[i]->SetParameter(j,fitFunc->GetParameter(j));
      fBkgdInSignalPlusBkgdRegion[i]->SetParameter(j,fitFunc->GetParameter(j));
      fBkgdInSignalPlusBkgdRegionLowScale[i]->SetParameter(j,fitFunc->GetParameter(j));
      fBkgdInSignalPlusBkgdRegionHighScale[i]->SetParameter(j,fitFunc->GetParameter(j));
      fBkgd[i]->SetParError(j,fitFunc->GetParError(j));
      fBkgdInSignalPlusBkgdRegion[i]->SetParError(j,fitFunc->GetParError(j));
      fBkgdInSignalPlusBkgdRegionLowScale[i]->SetParError(j,fitFunc->GetParError(j));
      fBkgdInSignalPlusBkgdRegionHighScale[i]->SetParError(j,fitFunc->GetParError(j));
    }
    //cerr<<"Line 619"<<endl;
    hCalcBkgd[i] = GetBackgroundErrorPW(hSignalPlusBkgd[i],fBkgd[i]);
    hCalcBkgd[i]->SetLineColor(kRed);
    hCalcBkgd[i]->SetMarkerColor(kRed);
    hCalcBkgdInSignalPlusBackgroundRegion[i] = (TH1F*) hCalcBkgd[i]->Clone(Form("%sInSignalPlusBackgroundRegion",hCalcBkgd[i]->GetName()));
    hCalcBkgdInSignalPlusBackgroundRegion[i]->Scale(etaRangeSigToBkgdRatio);
    hSignal[i]->Add(hCalcBkgdInSignalPlusBackgroundRegion[i],-1.0);
    cerr<<"Line 626->"<<endl;
    CalculateNetYield(hSignalPlusBkgd[i],fBkgdInSignalPlusBkgdRegion[i],i,0);

    cerr<<"<-Line 628"<<endl;
    hCalcBkgdInSignalPlusBackgroundRegionLowScale[i] = (TH1F*) hCalcBkgd[i]->Clone(Form("%sInSignalPlusBackgroundRegion",hCalcBkgd[i]->GetName()));
    hCalcBkgdInSignalPlusBackgroundRegionLowScale[i]->Scale(etaRangeSigToBkgdRatioLow);
    hSignalLowScale[i]->Add(hCalcBkgdInSignalPlusBackgroundRegionLowScale[i],-1.0);
    cerr<<"Line 633->"<<endl;
    CalculateNetYield(hSignalPlusBkgd[i],fBkgdInSignalPlusBkgdRegionLowScale[i],i,-2);

cerr<<"<-Line 634"<<endl; 
    hCalcBkgdInSignalPlusBackgroundRegionHighScale[i] = (TH1F*) hCalcBkgd[i]->Clone(Form("%sInSignalPlusBackgroundRegion",hCalcBkgd[i]->GetName()));
    hCalcBkgdInSignalPlusBackgroundRegionHighScale[i]->Scale(etaRangeSigToBkgdRatioHigh);
    hSignalHighScale[i]->Add(hCalcBkgdInSignalPlusBackgroundRegionHighScale[i],-1.0);
    cerr<<"Line 640->"<<endl;
    CalculateNetYield(hSignalPlusBkgd[i],fBkgdInSignalPlusBkgdRegionHighScale[i],i,2);
cerr<<"<-Line 642"<<endl;
    
    fBkgdInSignalPlusBkgdRegion[i]->SetParameter(0,etaRangeSigToBkgdRatio*fitFunc->GetParameter(0));//a function called in the previous function resets the background parameter but we need to do that to get the error bars correctly. 
    fBkgdInSignalPlusBkgdRegion[i]->SetParError(0,etaRangeSigToBkgdRatio*fitFunc->GetParError(0));
//cerr<<"Line 642"<<endl; 

  }
  //low R fits
  cerr<<"starting fit 2"<<endl;
  hCommonBkgd->Fit(fitFuncLowR);
  cerr<<"done with fit 2"<<endl;
  //the default is to give 95% confidence intervals but we want the standard 1 sigma (68%) confidence intervals
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hCommonBkgd, 0.68);
  for(Int_t i=0;i<nBinsInRxnPlane;i++){
    for(Int_t j=0;j<6;j++){
      fBkgdLowR[i]->SetParameter(j,fitFuncLowR->GetParameter(j));
      fBkgdInSignalPlusBkgdRegionLowR[i]->SetParameter(j,fitFuncLowR->GetParameter(j));
      fBkgdLowR[i]->SetParError(j,fitFuncLowR->GetParError(j));
      fBkgdInSignalPlusBkgdRegionLowR[i]->SetParError(j,fitFuncLowR->GetParError(j));
    }
    hCalcBkgdLowR[i] = GetBackgroundErrorPW(hSignalPlusBkgd[i],fBkgdLowR[i]);
    hCalcBkgdLowR[i]->SetLineColor(kRed);
    hCalcBkgdLowR[i]->SetMarkerColor(kRed);
    hCalcBkgdInSignalPlusBackgroundRegionLowR[i] = (TH1F*) hCalcBkgdLowR[i]->Clone(Form("%sInSignalPlusBackgroundRegionLowR",hCalcBkgdLowR[i]->GetName()));
    hCalcBkgdInSignalPlusBackgroundRegionLowR[i]->Scale(etaRangeSigToBkgdRatio);
    hSignalLowR[i]->Add(hCalcBkgdInSignalPlusBackgroundRegionLowR[i],-1.0);
    CalculateNetYield(hSignalPlusBkgd[i],fBkgdInSignalPlusBkgdRegionLowR[i],i,-1);
    //cout<<"test "<<nsYield[i]<<" +/- "<<nsYieldError[i]<<" +/- "<<nsYieldLowR[i]<<" - "<<nsYieldHighR[i]<<" +/- "<<nsYieldLowScale[i]<<" - "<<nsYieldHighScale[i]<<endl;//"\t";

    fBkgdInSignalPlusBkgdRegionLowR[i]->SetParameter(0,etaRangeSigToBkgdRatio*fitFuncLowR->GetParameter(0));//a function called in the previous function resets the background parameter but we need to do that to get the error bars correctly. 
    fBkgdInSignalPlusBkgdRegionLowR[i]->SetParError(0,etaRangeSigToBkgdRatio*fitFuncLowR->GetParError(0));
  }
  //high R fits
  cerr<<"starting fit 3"<<endl;
  hCommonBkgd->Fit(fitFuncHighR);
cerr<<"done with fit 3"<<endl;
    //the default is to give 95% confidence intervals but we want the standard 1 sigma (68%) confidence intervals
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hCommonBkgd, 0.68);
  for(Int_t i=0;i<nBinsInRxnPlane;i++){
    for(Int_t j=0;j<6;j++){
      fBkgdHighR[i]->SetParameter(j,fitFuncHighR->GetParameter(j));
      fBkgdInSignalPlusBkgdRegionHighR[i]->SetParameter(j,fitFuncHighR->GetParameter(j));
      fBkgdHighR[i]->SetParError(j,fitFuncHighR->GetParError(j));
      fBkgdInSignalPlusBkgdRegionHighR[i]->SetParError(j,fitFuncHighR->GetParError(j));
    }
    hCalcBkgdHighR[i] = GetBackgroundErrorPW(hSignalPlusBkgd[i],fBkgdHighR[i]);
    hCalcBkgdHighR[i]->SetLineColor(kRed);
    hCalcBkgdHighR[i]->SetMarkerColor(kRed);
    hCalcBkgdInSignalPlusBackgroundRegionHighR[i] = (TH1F*) hCalcBkgdHighR[i]->Clone(Form("%sInSignalPlusBackgroundRegionHighR",hCalcBkgdHighR[i]->GetName()));
    hCalcBkgdInSignalPlusBackgroundRegionHighR[i]->Scale(etaRangeSigToBkgdRatio);
    hSignalHighR[i]->Add(hCalcBkgdInSignalPlusBackgroundRegionHighR[i],-1.0);
    CalculateNetYield(hSignalPlusBkgd[i],fBkgdInSignalPlusBkgdRegionHighR[i],i,1);
    //cout<<"test "<<nsYield[i]<<" +/- "<<nsYieldError[i]<<" +/- "<<nsYieldLowR[i]<<" - "<<nsYieldHighR[i]<<" +/- "<<nsYieldLowScale[i]<<" - "<<nsYieldHighScale[i]<<endl;//"\t";
    fBkgdInSignalPlusBkgdRegionHighR[i]->SetParameter(0,etaRangeSigToBkgdRatio*fitFuncHighR->GetParameter(0));//a function called in the previous function resets the background parameter but we need to do that to get the error bars correctly. 
    fBkgdInSignalPlusBkgdRegionHighR[i]->SetParError(0,etaRangeSigToBkgdRatio*fitFuncHighR->GetParError(0));
  }

  //cout<<"Calculated near-side yields from dphi = "<<nsRealPhiCutLow<<" to "<<nsRealPhiCutHigh<<" and away-side yields from "<<asRealPhiCutLow<<" to "<<asRealPhiCutHigh<<endl;
  //Double_t *covMatrixPW = (TVirtualFitter::GetFitter())->GetCovarianceMatrix();
  hCommonBkgd->Draw();
  fitFuncHighR->Draw("same");
  fitFuncLowR->Draw("same");
  fitFunc->Draw("same");

  fitFuncHighR->SetLineColor(1);
return;
}


void CorrData::DrawSignalPlusBackground(){
  TCanvas *cSignalPlusBackground = new TCanvas("cSignalPlusBackground","Signal Plus Background",1000,800);
  cSignalPlusBackground->SetFillColor(10);
  cSignalPlusBackground->SetBorderMode(0);
  cSignalPlusBackground->SetBorderSize(0);
  cSignalPlusBackground->SetFrameFillColor(0);
  cSignalPlusBackground->SetFrameLineWidth(0);
  cSignalPlusBackground->Divide(3,2);
  for(Int_t i=0;i<nBinsInRxnPlane;i++){
    cSignalPlusBackground->cd(i+1);
    hSignalPlusBkgd[i]->Draw();
    fBkgdInSignalPlusBkgdRegion[i]->Draw("same");
    if(hCalcBkgdInSignalPlusBackgroundRegion[i]){
      hCalcBkgdInSignalPlusBackgroundRegion[i]->Draw("same");
    }
  }
  return;
}
void CorrData::DrawBackground(){
  TCanvas *cBackground = new TCanvas("cBackground","Background",1000,800);
  cBackground->SetFillColor(10);
  cBackground->SetBorderMode(0);
  cBackground->SetBorderSize(0);
  cBackground->SetFrameFillColor(0);
  cBackground->SetFrameLineWidth(0);
  cBackground->Divide(3,2);
  for(Int_t i=0;i<nBinsInRxnPlane;i++){
    cBackground->cd(i+1);
    hBkgd[i]->Draw();
    fBkgd[i]->Draw("same");
    fBkgdLowR[i]->Draw("same");
    fBkgdHighR[i]->Draw("same");
    fBkgdHighR[i]->SetLineColor(4);
    fBkgdLowR[i]->SetLineColor(1);
    fBkgdHighR[i]->SetLineStyle(2);
    fBkgdLowR[i]->SetLineStyle(2);
    if(hCalcBkgd[i]){
      hCalcBkgd[i]->Draw("same");
    }
  }
  return;
}
void CorrData::DrawSignal(){
  TCanvas *cSignal = new TCanvas("cSignal","Signal",1000,300);
  cSignal->SetFillColor(10);
  cSignal->SetBorderMode(0);
  cSignal->SetBorderSize(0);
  cSignal->SetFrameFillColor(0);
  cSignal->SetFrameLineWidth(0);
  //cSignal->Divide(3,2);
  Float_t fracA =0.05;// 0.2;
  Float_t fracB =0;//0.09;// 0.2;
   Float_t fractionYAxis =fracA+fracB;// 0.2;
   Float_t fractionXAxis =0.15;// 0.12;
   Float_t xwidth = (1.0-fractionYAxis)/6.0;//these are crossed because the x-axis takes up room in the y direction and vice versa
   Float_t ywidth = (1.0-fractionXAxis)/1.0;
    Float_t topmargin = 0.02;
    Float_t leftmargin = fracA/(fracA+xwidth);
    Float_t rightmargin = fracB/(fracB+xwidth);
    Float_t bottommargin = fractionXAxis/(fractionXAxis+ywidth);
   //cout<<"xwidth "<<xwidth <<" ywidth "<<ywidth<<endl;
   //Double_t xlow, Double_t ylow, Double_t xup, Double_t yup
   // pad1 pad2 pad3
   // pad4 pad5 pad6
   TPad *part1 = new TPad("part1","part1", 0.0                 ,0.0,1.0-5.0*xwidth-fracB,1.0       ,1);
   TPad *part2 = new TPad("part2","part2", 1.0-5.0*xwidth-fracB,0.0,1.0-4.0*xwidth-fracB,1.0       ,2);
   TPad *part3 = new TPad("part3","part3", 1.0-4.0*xwidth-fracB,0.0,1.0-3.0*xwidth-fracB,1.0       ,3);
   TPad *part4 = new TPad("part4","part4", 1.0-3.0*xwidth-fracB,0.0,1.0-2.0*xwidth-fracB,1.0       ,4);
   TPad *part5 = new TPad("part5","part5", 1.0-2.0*xwidth-fracB,0.0,1.0-1.0*xwidth-fracB,1.0       ,5);
   TPad *part6 = new TPad("part6","part6",1.0-1.0*xwidth-fracB,0.0,1.0                 ,1.0       ,6);
   part1->SetFillColor(0);
   part2->SetFillColor(0);
   part3->SetFillColor(0);
   part4->SetFillColor(0);
   part5->SetFillColor(0);
   part6->SetFillColor(0);
   part1->SetBorderMode(0);
   part1->SetBorderSize(0);
   part1->SetFrameFillColor(0);
   part1->SetFrameBorderMode(0);
   part2->SetBorderMode(0);
   part2->SetBorderSize(0);
   part2->SetFrameFillColor(0);
   part2->SetFrameBorderMode(0);
   part3->SetBorderMode(0);
   part3->SetBorderSize(0);
   part3->SetFrameFillColor(0);
   part3->SetFrameBorderMode(0);
   part4->SetBorderMode(0);
   part4->SetBorderSize(0);
   part4->SetFrameFillColor(0);
   part4->SetFrameBorderMode(0);
   part5->SetBorderMode(0);
   part5->SetBorderSize(0);
   part5->SetFrameFillColor(0);
   part5->SetFrameBorderMode(0);
   part6->SetBorderMode(0);
   part6->SetBorderSize(0);
   part6->SetFrameFillColor(0);
   part6->SetFrameBorderMode(0);
   part1->SetTopMargin(topmargin);
   part2->SetTopMargin(topmargin);
   part3->SetTopMargin(topmargin);
   part4->SetTopMargin(topmargin);
   part5->SetTopMargin(topmargin);
   part6->SetTopMargin(topmargin);

   part1->SetBottomMargin(bottommargin);
   part2->SetBottomMargin(bottommargin);
   part3->SetBottomMargin(bottommargin);
   part4->SetBottomMargin(bottommargin);
   part5->SetBottomMargin(bottommargin);
   part6->SetBottomMargin(bottommargin);

   part1->SetLeftMargin(leftmargin);
   part2->SetLeftMargin(0.0);
   part3->SetLeftMargin(0.00);
   part4->SetLeftMargin(0.00);
   part5->SetLeftMargin(0.00);
   part6->SetLeftMargin(0.00);

   part6->SetRightMargin(rightmargin);
   part1->SetRightMargin(0.0);
   part2->SetRightMargin(0.0);
   part3->SetRightMargin(0.0);
   part4->SetRightMargin(0.0);
   part5->SetRightMargin(0.0);
   part1->Draw();
   part2->Draw();
   part3->Draw();
   part4->Draw();
   part5->Draw();
   part6->Draw();
   //part1
   TPad *pads[] = {part1,part2,part3,part4,part5,part6};
   TLatex *tex[] = {NULL,NULL,NULL,NULL,NULL,NULL};
   TString psi1 = "0<#phi^{t}-#psi<#pi/12";
   TString psi2 = "#pi/12<#phi^{t}-#psi<2#pi/12";
   TString psi3 = "2#pi/12<#phi^{t}-#psi<3#pi/12";
   TString psi4 = "3#pi/12<#phi^{t}-#psi<4#pi/12";
   TString psi5 = "4#pi/12<#phi^{t}-#psi<5#pi/12";
   TString psi6 = "5#pi/12<#phi^{t}-#psi<6#pi/12";
   TString psis[] = {psi1,psi2,psi3,psi4,psi5,psi6};
   Float_t size[] = {0.080292,0.10219,0.10219,0.10219,0.10219,0.10219};
   //return;
   Float_t max = -1;
   Float_t min = 1.0;
   for(Int_t i=0;i<nBinsInRxnPlane;i++){
     if(hSignal[i]->GetMaximum()>max) max = hSignal[i]->GetMaximum();
     if(hSignal[i]->GetMinimum()<min) min = hSignal[i]->GetMinimum();
   }
   max = max + (max-min)*0.1;
   min = min - (max-min)*0.1;
  cout<<"Max "<< max <<" min "<<min<<endl;
  for(Int_t i=0;i<nBinsInRxnPlane;i++){
    //cSignal->cd(i+1);
    pads[i]->cd();
    hSignal[i]->SetMaximum(max);
    hSignal[i]->SetMinimum(min);
    hSignal[i]->GetXaxis()->SetTitle("#Delta#phi");
    hSignal[i]->GetXaxis()->CenterTitle();
    hSignal[i]->GetYaxis()->SetTitleOffset(1.5);
    hSignal[i]->GetXaxis()->SetTitleSize(0.13);
    hSignal[i]->GetYaxis()->SetTitleSize(0.07);
    hSignal[i]->GetXaxis()->SetLabelSize(0.07);
    hSignal[i]->GetYaxis()->SetLabelSize(0.07);
    hSignal[i]->GetXaxis()->SetLabelOffset(-0.01);
    hSignal[i]->GetXaxis()->SetLabelSize(0.08);
    hSignal[i]->GetXaxis()->SetTitleSize(0.11);
    hSignal[i]->GetXaxis()->SetTitleOffset(0.6);
    if(i>0){
      hSignal[i]->GetXaxis()->SetTitleSize(0.13);
      hSignal[i]->GetXaxis()->SetLabelSize(0.1);
      hSignal[i]->GetXaxis()->SetLabelOffset(-0.025);
      hSignal[i]->GetXaxis()->SetTitleOffset(0.5);
    }
    hSignal[i]->GetYaxis()->SetTitle("1/N_{trig}dN/d#Delta#phi");
    hSignal[i]->Draw();
//     hSignalLowR[i]->Draw("same");
//     hSignalHighR[i]->Draw("same");
//     hSignalLowR[i]->SetLineColor(2);
//     hSignalHighR[i]->SetLineColor(2);
    //hSignalLowScale[i]->Draw("same");
    //hSignalHighScale[i]->Draw("same");
    //hSignalLowScale[i]->SetLineColor(2);
    //hSignalHighScale[i]->SetLineColor(2);
    TGraph *graph = GetErrorGraph(hSignalLowScale[i],hSignalHighScale[i]);
    graph->SetFillColor(kYellow);
    graph->Draw("f same");
    TGraph *graph2 = GetErrorGraph(hSignalLowR[i],hSignalHighR[i]);
    graph2->SetFillColor(1);
    //graph2->Draw("f same");
    hSignal[i]->Draw("same");
    tex[i] = new TLatex(-1.12456,min+(max-min)*0.93,psis[i].Data());
    tex[i]->SetTextSize(size[i]);
    tex[i]->Draw();
  }
 return; 
}

//From TF1.cxx
Double_t CorrData::GradientParPW(Int_t ipar, const Double_t x, Double_t eps, TF1 *inputFuncForErrors, Double_t *fParams) {
   // Compute the gradient (derivative) wrt a parameter ipar
   // Parameters:
   // ipar - index of parameter for which the derivative is computed
   // x - point, where the derivative is computed
   // eps - if the errors of parameters have been computed, the step used in
   // numerical differentiation is eps*parameter_error.
   // if the errors have not been computed, step=eps is used
   // default value of eps = 0.01
   // Method is the same as in Derivative() function
   //
   // If a paramter is fixed, the gradient on this parameter = 0

   if (inputFuncForErrors->GetNumberFreeParameters() == 0){
     cerr<<"No free parameters!"<<endl;
     return 0; 
   }

   Double_t h;
   Double_t al, bl;
   Double_t f1, f2, g1, g2, h2, d0, d2;

   inputFuncForErrors->GetParLimits(ipar,al,bl);
   if (al*bl != 0 && al >= bl) {
      //this parameter is fixed
      cerr<<"Par "<<ipar<<"is fixed!"<<endl;
      return 0;
   }

   // check if error has been computer (is not zero)
   if (inputFuncForErrors->GetParError(ipar)!=0) h = eps*inputFuncForErrors->GetParError(ipar);
   else h=eps;

   //h = 0.01;

   // get parameter
   Double_t par0 = inputFuncForErrors->GetParameter(ipar);

   // set parameter to par + small value
   inputFuncForErrors->SetParameter(ipar,par0 + h);
   // evaluate function f1 at x  
   f1 = inputFuncForErrors->Eval(x);
   // set parameter to par - small value
   inputFuncForErrors->SetParameter(ipar,par0 - h);     
   // evaluate function f2 at x
   f2 = inputFuncForErrors->Eval(x);

   inputFuncForErrors->SetParameter(ipar,par0 + h/2);   
   g1 = inputFuncForErrors->Eval(x);
   inputFuncForErrors->SetParameter(ipar,par0 - h/2);   
   g2 = inputFuncForErrors->Eval(x);
   inputFuncForErrors->SetParameter(ipar,par0);

   //compute the central differences
   h2    = 1/(2.*h);
   d0    = f1 - f2;
   d2    = 2*(g1 - g2);

   // compute gradient
   Double_t  grad = h2*(4*d2 - d0)/3.;

   return grad;
}

//______________________________________________________________________________
void CorrData::GradientParPW(const Double_t x, Double_t *grad, TF1 *inputFuncForErrors,   Double_t *fParams, Double_t eps)
{
   // Compute the gradient wrt parameters
   // Parameters:
   // x - point, were the gradient is computed
   // grad - used to return the computed gradient, assumed to be of at least fNpar size
   // eps - if the errors of parameters have been computed, the step used in
   // numerical differentiation is eps*parameter_error.
   // if the errors have not been computed, step=eps is used
   // default value of eps = 0.01
   // Method is the same as in Derivative() function
   //
   // If a paramter is fixed, the gradient on this parameter = 0
   // 0.01
   if(eps< 1e-10 || eps > 1) {
      Warning("Derivative","parameter esp=%g out of allowed range[1e-10,1], reset to 0.01",eps);
      eps = 0.01;
   }

   Int_t fNpar = inputFuncForErrors->GetNumberFreeParameters();
   for (Int_t ipar=0; ipar<fNpar; ipar++){
     grad[ipar] = GradientParPW(ipar,x,eps,inputFuncForErrors,fParams);
   }
   return;
}
//From TF1.cxx
Double_t CorrData::DERPW(Int_t ipar, const Double_t x, Double_t eps, TF1 *inputFuncForErrors, Double_t *fParams) {
   // Compute the gradient (derivative) wrt a parameter ipar
   // Parameters:
   // ipar - index of parameter for which the derivative is computed
   // x - point, where the derivative is computed
   // eps - if the errors of parameters have been computed, the step used in
   // numerical differentiation is eps*parameter_error.
   // if the errors have not been computed, step=eps is used 
   // default value of eps = 0.01
   // Method is the same as in Derivative() function
   //
   // If a paramter is fixed, the gradient on this parameter = 0

   if (inputFuncForErrors->GetNumberFreeParameters() == 0){
     cerr<<"No free parameters!"<<endl;
     return 0; 
   }

   Double_t h;
   Double_t al, bl;
   Double_t f1, f2;

   inputFuncForErrors->GetParLimits(ipar,al,bl);
   if (al*bl != 0 && al >= bl) {
      //this parameter is fixed
      cerr<<"Par "<<ipar<<"is fixed!"<<endl;
      return 0;
   }

   // check if error has been computed (is not zero)
   if (inputFuncForErrors->GetParError(ipar)!=0) h = eps*inputFuncForErrors->GetParError(ipar);
   else h=eps;

   // get parameter
   Double_t par0 = inputFuncForErrors->GetParameter(ipar);

   // set parameter to par + small value
   inputFuncForErrors->SetParameter(ipar,par0 + h);
   // evaluate function f1 at x  
   f1 = inputFuncForErrors->Eval(x);
   // set parameter to par - small value
   inputFuncForErrors->SetParameter(ipar,par0 - h);     
   // evaluate function f2 at x
   f2 = inputFuncForErrors->Eval(x);
   inputFuncForErrors->SetParameter(ipar,par0); // need to reset back to default?

   Double_t gradX = (f1 - f2)/(2*h);

   return gradX;
}

//______________________________________________________________________________
void CorrData::DerivPW(const Double_t x, Double_t *gradX, TF1 *inputFuncForErrors, Double_t *fParams, Double_t eps){
   // Compute the gradient wrt parameters
   // Parameters:
   // x - point, were the gradient is computed
   // grad - used to return the computed gradient, assumed to be of at least fNpar size
   // eps - if the errors of parameters have been computed, the step used in
   // numerical differentiation is eps*parameter_error.
   // if the errors have not been computed, step=eps is used
   // default value of eps = 0.01
   // Method is the same as in Derivative() function
   //
   // If a paramter is fixed, the gradient on this parameter = 0
   // 0.01
   if(eps< 1e-10 || eps > 1) {
      Warning("Derivative","parameter esp=%g out of allowed range[1e-10,1], reset to 0.01",eps);
      eps = 0.01;
   }

   Int_t fNpar = inputFuncForErrors->GetNumberFreeParameters();
   for (Int_t ipar=0; ipar<fNpar; ipar++){
     gradX[ipar] = DERPW(ipar,x,eps,inputFuncForErrors,fParams);
   }
   return;
}

TH1F *CorrData::GetBackgroundErrorPW(TH1F *inputHisto, TF1 *inputFunc){
   // Name and Create histo
   TString outputName = inputHisto->GetName();
   outputName += "BackgroundError";
   TH1F *hOutput = (TH1F*) inputHisto->Clone(outputName.Data());

   // get parameter info
//   Int_t npar = inputFunc->GetNumberFreeParameters();
//   Int_t npar_real = inputFunc->GetNpar();
   const Int_t npar = fitFunc->GetNumberFreeParameters();
   const Int_t npar_real = fitFunc->GetNpar();
   Double_t   *fParams   = new Double_t[npar];
   for(int i=0;i<npar;i++){ 
     fParams[i] = fitFunc->GetParameter(i); 
   }   // loop over free parameters?

   Double_t *grad = new Double_t[npar_real];
   Double_t *gradX = new Double_t[npar_real]; // added
   Double_t *sum_vector = new Double_t[npar];
   Bool_t *parIsFixed=0;
   Double_t parLimLow, parLimHigh;

   // determine which parameters are fixed
   if (npar_real != npar){
      parIsFixed = new Bool_t[npar_real];

      // loop over real parameters
      for (Int_t ipar=0; ipar<npar_real; ipar++){
         parIsFixed[ipar]=0;

         // get parameter limits 
         fitFunc->GetParLimits(ipar,parLimLow,parLimHigh);
         if (parLimLow*parLimHigh != 0 && parLimLow >= parLimHigh) {
            //this parameter is fixed
            parIsFixed[ipar]=1;
         }
      }
   }


// ========================================================================================================================
   // Covariance matrix elements already stored here..
   //Double_t *covMat = new Double_t[npar][npar];

   //cout<<"Get student quantile.."<<endl;   
   //cout<<"NDFpw: "<<fitFunc->GetNDF()<<"   NDF: "<<inputFunc->GetNDF()<<endl;
//   Double_t t = TMath::StudentQuantile(0.5 + confidenceLimit/2, fitFunc->GetNDF());   
//   Double_t chidf = TMath::Sqrt(fitFunc->GetChisquare()/fitFunc->GetNDF());
   Int_t igrad, ifree=0;
   Double_t c=0;
   Int_t npoints = hOutput->GetXaxis()->GetNbins();
   //cout<<"tried getting student quantile and chidf - about to loop over points now"<<endl;

   //cout<<"About to loop over points..."<<endl;
  // loop over bins
  for (Int_t ipoint=1; ipoint<=npoints; ipoint++){
    c = 0;
    //the input histogram runs from -pi/2 to 3/pi/2
    const Double_t x = hOutput->GetXaxis()->GetBinCenter(ipoint);
    const Double_t xshifted = hOutput->GetXaxis()->GetBinCenter(ipoint);

    // fill gradX or grad with gradient values
    //cout<<"Evaluate df/dx at x = "<<x<<endl;
    DerivPW(x, gradX, inputFunc, fParams);
    //DerivPW(x, gradX, inputFunc, fParams);
    // loop over (ROWS) parameters to get gradient - for all parameters
    for (Int_t irow=0; irow<npar; irow++){
      sum_vector[irow]=0;

      // loop over (COLS) parameters to get gradient - for all parameters
      for (Int_t icol=0; icol<npar; icol++){
        igrad = 0;
        ifree=0;
	    
        if (parIsFixed) {//if any parameters were fixed
          //find the free parameter #icol
          while (ifree<icol+1){
            if (parIsFixed[igrad]==0) ifree++;
            igrad++;
          }
          igrad--;
          //now the [igrad] element of gradient corresponds to [icol] element of cov.matrix
        } else { igrad = icol; }

	Double_t covariance = (TVirtualFitter::GetFitter())->GetCovarianceMatrixElement(irow, icol);
        c+=gradX[icol]*gradX[irow]*covariance;
      } // loop over 'cols'
    } // loop over 'rows'

    // take square root of error
    c=TMath::Sqrt(c);

    hOutput->SetBinContent(ipoint, inputFunc->Eval(x));
    hOutput->SetBinError(ipoint, c);


  } // loop over bins
  //hOutput->SetMarkerStyle(2);
  delete[] grad;
  delete[] gradX;
  delete[] sum_vector;
  delete parIsFixed;
  delete[] fParams;
  return hOutput;
}
void CorrData::CalculateNetYield(TH1 *histo, TF1 *inputFunc, Int_t rxnPlaneBin, Int_t option){
  //cout<<"Calculating yields from "<<histo->GetName()<<" using function "<<inputFunc->GetName()<<endl;
  Double_t *covMatrixPW = (TVirtualFitter::GetFitter())->GetCovarianceMatrix();
  const Int_t npar = fitFunc->GetNumberFreeParameters();
  Double_t   *fParams   = new Double_t[npar];
  for(int i=0;i<npar;i++){ //this uses the parameters with the correct error bars but then we need to scale up the background 
    fParams[i] = fitFunc->GetParameter(i); 
  }   // loop over free parameters

  //Near Side
  // get low and high bin of range for yield
  Int_t lowBin = histo->GetXaxis()->FindBin(-nsPhiCut+1e-3);
  Int_t highBin = histo->GetXaxis()->FindBin(nsPhiCut-1e-3);
  // get bin width
  Double_t binw = histo->GetBinWidth(lowBin);
  // get real ranges of integral on histo from bin edges
  nsRealPhiCutLow = histo->GetXaxis()->GetBinLowEdge(lowBin);
  nsRealPhiCutHigh = histo->GetXaxis()->GetBinUpEdge(highBin);
  //cout<<"Ns yield covers bins "<<lowBin<<" - "<<highBin<<", phi "<<nsRealPhiCutLow<<" - "<<nsRealPhiCutHigh<<endl;
  // get yield and yield error of histo
  Double_t yielderr = 0.0;
  Double_t yield = histo->IntegralAndError(lowBin,highBin,yielderr,"width");
  // get integral and integral error of fit function of background for yield calculation
  Double_t background = etaRangeSigToBkgdRatio*inputFunc->Integral(nsRealPhiCutLow, nsRealPhiCutHigh, 1E-12);   
  Double_t backgrounderr = etaRangeSigToBkgdRatio*inputFunc->IntegralError(nsRealPhiCutLow, nsRealPhiCutHigh, fParams, covMatrixPW, 0.01);
  if(option==-2){
    background = etaRangeSigToBkgdRatioLow*inputFunc->Integral(nsRealPhiCutLow, nsRealPhiCutHigh, 1E-12);   
    backgrounderr = etaRangeSigToBkgdRatioLow*inputFunc->IntegralError(nsRealPhiCutLow, nsRealPhiCutHigh, fParams, covMatrixPW, 0.01);
    nsYieldLowScale[rxnPlaneBin] = (yield - background);
    nsYieldErrorLowScale[rxnPlaneBin] = 1.0*TMath::Sqrt(yielderr*yielderr + backgrounderr*backgrounderr);
  }
  if(option==2){
    background = etaRangeSigToBkgdRatioHigh*inputFunc->Integral(nsRealPhiCutLow, nsRealPhiCutHigh, 1E-12);   
    backgrounderr = etaRangeSigToBkgdRatioHigh*inputFunc->IntegralError(nsRealPhiCutLow, nsRealPhiCutHigh, fParams, covMatrixPW, 0.01);
    nsYieldHighScale[rxnPlaneBin] = (yield - background);
    nsYieldErrorHighScale[rxnPlaneBin] = 1.0*TMath::Sqrt(yielderr*yielderr + backgrounderr*backgrounderr);
  }
  // get difference and scale to make readable
  if(option ==0){
    nsYield[rxnPlaneBin] = (yield - background);
    nsYieldError[rxnPlaneBin] = 1.0*TMath::Sqrt(yielderr*yielderr + backgrounderr*backgrounderr);
  }
  if(option ==-1){
    nsYieldLowR[rxnPlaneBin] = (yield - background);
    nsYieldErrorLowR[rxnPlaneBin] = 1.0*TMath::Sqrt(yielderr*yielderr + backgrounderr*backgrounderr);
  }
  if(option ==1){
    nsYieldHighR[rxnPlaneBin] = (yield - background);
    nsYieldErrorHighR[rxnPlaneBin] = 1.0*TMath::Sqrt(yielderr*yielderr + backgrounderr*backgrounderr);
  }

  //Away Side
  // get low and high bin of range for yield
  lowBin = histo->GetXaxis()->FindBin(-asPhiCut+mypi+1e-3);
  highBin = histo->GetXaxis()->FindBin(asPhiCut+mypi-1e-3);
  // get bin width
  binw = histo->GetBinWidth(lowBin);
  // get real ranges of integral on histo from bin edges
  asRealPhiCutLow = histo->GetXaxis()->GetBinLowEdge(lowBin);
  asRealPhiCutHigh = histo->GetXaxis()->GetBinUpEdge(highBin);
  //cout<<"As yield covers bins "<<lowBin<<" - "<<highBin<<", phi "<<asRealPhiCutLow<<" - "<<asRealPhiCutHigh<<endl;
  // get yield and yield error of histo
  yielderr = 0.0;
  yield = histo->IntegralAndError(lowBin,highBin,yielderr,"width");
  // get integral and integral error of fit function of background for yield calculation
  background = etaRangeSigToBkgdRatio*inputFunc->Integral(nsRealPhiCutLow, nsRealPhiCutHigh, 1E-12);   
  backgrounderr = etaRangeSigToBkgdRatio*inputFunc->IntegralError(nsRealPhiCutLow, nsRealPhiCutHigh, fParams, covMatrixPW, 0.01);
  // get difference and scale to make readable
  if(option==0){
    asYield[rxnPlaneBin] = (yield - background);
    asYieldError[rxnPlaneBin] = 1.0*TMath::Sqrt(yielderr*yielderr + backgrounderr*backgrounderr);
    CalculateRMS(histo,inputFunc,nsRealPhiCutLow, nsRealPhiCutHigh, 0.0, rxnPlaneBin,nsRMS[rxnPlaneBin],nsRMSError[rxnPlaneBin] ,0);
    CalculateRMS(histo,inputFunc,asRealPhiCutLow, asRealPhiCutHigh, mypi, rxnPlaneBin,asRMS[rxnPlaneBin],asRMSError[rxnPlaneBin] ,0 );
  }
  if(option==-2){
    background = etaRangeSigToBkgdRatioLow*inputFunc->Integral(nsRealPhiCutLow, nsRealPhiCutHigh, 1E-12);   
    backgrounderr = etaRangeSigToBkgdRatioLow*inputFunc->IntegralError(nsRealPhiCutLow, nsRealPhiCutHigh, fParams, covMatrixPW, 0.01);
    asYieldLowScale[rxnPlaneBin] = (yield - background);
    asYieldErrorLowScale[rxnPlaneBin] = 1.0*TMath::Sqrt(yielderr*yielderr + backgrounderr*backgrounderr);
    CalculateRMS(histo,inputFunc,nsRealPhiCutLow, nsRealPhiCutHigh, 0.0, rxnPlaneBin,nsRMSLowScale[rxnPlaneBin],nsRMSErrorLowScale[rxnPlaneBin] ,-2);
    CalculateRMS(histo,inputFunc,asRealPhiCutLow, asRealPhiCutHigh, mypi, rxnPlaneBin,asRMSLowScale[rxnPlaneBin],asRMSErrorLowScale[rxnPlaneBin] ,-2 );
  }
  if(option==2){
    background = etaRangeSigToBkgdRatioHigh*inputFunc->Integral(nsRealPhiCutLow, nsRealPhiCutHigh, 1E-12);   
    backgrounderr = etaRangeSigToBkgdRatioHigh*inputFunc->IntegralError(nsRealPhiCutLow, nsRealPhiCutHigh, fParams, covMatrixPW, 0.01);
    asYieldHighScale[rxnPlaneBin] = (yield - background);
    asYieldErrorHighScale[rxnPlaneBin] = 1.0*TMath::Sqrt(yielderr*yielderr + backgrounderr*backgrounderr);
    CalculateRMS(histo,inputFunc,nsRealPhiCutLow, nsRealPhiCutHigh, 0.0, rxnPlaneBin,nsRMSHighScale[rxnPlaneBin],nsRMSErrorHighScale[rxnPlaneBin] ,2);
    CalculateRMS(histo,inputFunc,asRealPhiCutLow, asRealPhiCutHigh, mypi, rxnPlaneBin,asRMSHighScale[rxnPlaneBin],asRMSErrorHighScale[rxnPlaneBin] ,2 );
  }
  if(option==-1){
    asYieldLowR[rxnPlaneBin] = (yield - background);
    asYieldErrorLowR[rxnPlaneBin] = 1.0*TMath::Sqrt(yielderr*yielderr + backgrounderr*backgrounderr);
    CalculateRMS(histo,inputFunc,nsRealPhiCutLow, nsRealPhiCutHigh, 0.0, rxnPlaneBin,nsRMSLowR[rxnPlaneBin],nsRMSErrorLowR[rxnPlaneBin] ,-1);
    CalculateRMS(histo,inputFunc,asRealPhiCutLow, asRealPhiCutHigh, mypi, rxnPlaneBin,asRMSLowR[rxnPlaneBin],asRMSErrorLowR[rxnPlaneBin] ,-1 );
  }
  if(option==1){
    asYieldHighR[rxnPlaneBin] = (yield - background);
    asYieldErrorHighR[rxnPlaneBin] = 1.0*TMath::Sqrt(yielderr*yielderr + backgrounderr*backgrounderr);
    CalculateRMS(histo,inputFunc,nsRealPhiCutLow, nsRealPhiCutHigh, 0.0, rxnPlaneBin,nsRMSHighR[rxnPlaneBin],nsRMSErrorHighR[rxnPlaneBin] ,1);
    CalculateRMS(histo,inputFunc,asRealPhiCutLow, asRealPhiCutHigh, mypi, rxnPlaneBin,asRMSHighR[rxnPlaneBin],asRMSErrorHighR[rxnPlaneBin] ,1 );
  }


  //cout<<"RMS: "<<RMSSigPlusBkgdNS<<" +/- "<<RMSErrSigPlusBkgdNS<<"\t"<<RMSSigPlusBkgdAS<<" +/- "<<RMSErrSigPlusBkgdAS<<endl;

  delete[] fParams;
  return;// netyieldel;

}

void CorrData::MakeYieldTable(char *filename){
  TString nsyields = "NSYields";
  TString nsyieldname = nsyields+filename;
  TString asyields = "ASYields";
  TString asyieldname = asyields+filename;
  TString nsrms = "NSRMS";
  TString nsrmsname = nsrms+filename;
  TString asrms = "ASRMS";
  TString asrmsname = asrms+filename;
  TString corrs = "Correlations";
  TString corrsname = corrs+filename;
  ofstream myfilensy;
  ofstream myfileasy;
  ofstream myfilensr;
  ofstream myfileasr;
  ofstream myfilecorr;
  myfilensy.open (nsyieldname.Data());
  myfileasy.open (asyieldname.Data());
  myfilensr.open (nsrmsname.Data());
  myfileasr.open (asrmsname.Data());
  myfilecorr.open (corrsname.Data());
  for(Int_t i=0;i<nBinsInRxnPlane;i++){
    //cout<<nsYield[i]<<" +/- "<<nsYieldError[i]<<" +/- "<<nsYieldLowR[i]<<" - "<<nsYieldHighR[i]<<" +/- "<<nsYieldLowScale[i]<<" - "<<nsYieldHighScale[i]<<endl;//"\t";
    //cout<<asYield[i]<<" +/- "<<asYieldError[i]<<" +/- "<<asYieldLowR[i]<<" - "<<asYieldHighR[i]<<" +/- "<<asYieldLowScale[i]<<" - "<<asYieldHighScale[i]<<endl;
    myfilensy<<nsYield[i]<<"\t"<<nsYieldError[i]<<"\t"<<TMath::Abs(nsYieldLowR[i]-nsYieldHighR[i])/2<<"\t"<<TMath::Abs(nsYieldLowScale[i]-nsYieldHighScale[i])/2<<endl;
    myfileasy<<asYield[i]<<"\t"<<asYieldError[i]<<"\t"<<TMath::Abs(asYieldLowR[i]-asYieldHighR[i])/2<<"\t"<<TMath::Abs(asYieldLowScale[i]-asYieldHighScale[i])/2<<endl;
    //cout<<nsRMS[i]<<" +/- "<<nsRMSError[i]<<" +/- "<<nsRMSLowR[i]<<" - "<<nsRMSHighR[i]<<" +/- "<<nsRMSLowScale[i]<<" - "<<nsRMSHighScale[i]<<endl;
    //cout<<asRMS[i]<<" +/- "<<asRMSError[i]<<" +/- "<<asRMSLowR[i]<<" - "<<asRMSHighR[i]<<" +/- "<<asRMSLowScale[i]<<" - "<<asRMSHighScale[i]<<endl;;
    myfilensr<<nsRMS[i]<<"\t"<<nsRMSError[i]<<"\t"<<TMath::Abs(nsRMSLowR[i]-nsRMSHighR[i])/2<<"\t"<<TMath::Abs(nsRMSLowScale[i]-nsRMSHighScale[i])/2<<endl;
    myfileasr<<asRMS[i]<<"\t"<<asRMSError[i]<<"\t"<<TMath::Abs(asRMSLowR[i]-asRMSHighR[i])/2<<"\t"<<TMath::Abs(asRMSLowScale[i]-asRMSHighScale[i])/2<<endl;;
    //myfile<<endl;
  }
  myfilensy.close ();
  myfileasy.close ();
  myfilensr.close ();
  myfileasr.close ();
  for(Int_t i=1;i<=hSignal[i]->GetNbinsX();i++){
    myfilecorr<<hSignal[0]->GetBinCenter(i);
    for(Int_t j=0;j<nBinsInRxnPlane;j++){
      myfilecorr<<"\t"<<hSignal[j]->GetBinContent(i)<<"\t"<<hSignal[j]->GetBinError(i)<<"\t"<<hSignalLowR[j]->GetBinContent(i)<<"\t"<<hSignalHighR[j]->GetBinContent(i)<<"\t"<<hSignalLowScale[j]->GetBinContent(i)<<"\t"<<hSignalHighScale[j]->GetBinContent(i);
    }
    myfilecorr<<endl;
  }
  myfilecorr.close();
  //myfile.close();
}
void CorrData::CalculateRMS(TH1 *histo, TF1 *fMyBackground,Double_t lowrange, Double_t highrange, Double_t center,Int_t rxnplanebin, Double_t & RMS, Double_t & RMSError, Int_t option){
  Int_t lowBin = histo->GetXaxis()->FindBin(lowrange+1e-3);
  Int_t highBin = histo->GetXaxis()->FindBin(highrange-1e-3);  
  Double_t sum = 0.0;
  Double_t sumx = 0.0;
  Double_t sumxx = 0.0;
  Double_t part1err = 0.0;//this is square of error until we take the square root at the end!
  Double_t scale = etaRangeSigToBkgdRatio;
  if(option==-2) scale = etaRangeSigToBkgdRatioLow;
  if(option==2)scale = etaRangeSigToBkgdRatioHigh;
  for(Int_t i=lowBin;i<=highBin;i++){
    Double_t bincontent = histo->GetBinContent(i);
    Double_t bincontenterr = histo->GetBinError(i);
    Double_t bincenter = histo->GetBinCenter(i);
    Double_t binwidth = histo->GetBinWidth(i);
    Double_t mybkgd = scale*fMyBackground->Eval(bincenter);
    sum += (bincontent-mybkgd)*binwidth;
    sumx += TMath::Abs(bincontent-mybkgd)*(bincenter-center)*binwidth;
    sumxx += TMath::Abs(bincontent-mybkgd)*(bincenter-center)*(bincenter-center)*binwidth;
    //treat the bin center as a constant.  We could do a better estimate by estimating the error on the slope from the neighboring bins 
    //we also might want to do a better estimate of the x value than the bin center but we can come back to that once the basics are working
    part1err += TMath::Power(bincontenterr*(bincenter-center),2.0)*binwidth;
//       Double_t altbincenter = 0;
//       Double_t altbincentererr = 0;
//       GetBinCenter(histo,i,altbincenter,altbincentererr);
  }
  Double_t meanx = sumx/sum;
  Double_t meanxx = sumxx/sum;
  RMS = TMath::Sqrt(TMath::Abs(meanxx));

  TF1 *fBkgdTemp =NULL;
  if(option==0 ||option==2 ||option==-2 ){
    if(center<1){//near side
      switch(rxnplanebin){
      case 0:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxNS0, MINPHI,MAXPHI , nParameters);
	break;
      case 1:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxNS1, MINPHI,MAXPHI , nParameters);
	break;
      case 2:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxNS2, MINPHI,MAXPHI , nParameters);
	break;
      case 3:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxNS3, MINPHI,MAXPHI , nParameters);
	break;
      case 4:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxNS4, MINPHI,MAXPHI , nParameters);
	break;
      default:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxNS5, MINPHI,MAXPHI , nParameters);
	break;
      }
    }
    else{
      switch(rxnplanebin){
      case 0:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxAS0, MINPHI,MAXPHI , nParameters);
	break;
      case 1:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxAS1, MINPHI,MAXPHI , nParameters);
	break;
      case 2:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxAS2, MINPHI,MAXPHI , nParameters);
	break;
      case 3:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxAS3, MINPHI,MAXPHI , nParameters);
	break;
      case 4:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxAS4, MINPHI,MAXPHI , nParameters);
	break;
      default:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxAS5, MINPHI,MAXPHI , nParameters);
	break;
      }
    }
  }
  if(option==-1){
    if(center<1){//near side
      switch(rxnplanebin){
      case 0:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxNS0LowR, MINPHI,MAXPHI , nParameters);
	break;
      case 1:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxNS1LowR, MINPHI,MAXPHI , nParameters);
	break;
      case 2:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxNS2LowR, MINPHI,MAXPHI , nParameters);
	break;
      case 3:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxNS3LowR, MINPHI,MAXPHI , nParameters);
	break;
      case 4:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxNS4LowR, MINPHI,MAXPHI , nParameters);
	break;
      default:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxNS5LowR, MINPHI,MAXPHI , nParameters);
	break;
      }
    }
    else{
      switch(rxnplanebin){
      case 0:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxAS0LowR, MINPHI,MAXPHI , nParameters);
	break;
      case 1:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxAS1LowR, MINPHI,MAXPHI , nParameters);
	break;
      case 2:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxAS2LowR, MINPHI,MAXPHI , nParameters);
	break;
      case 3:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxAS3LowR, MINPHI,MAXPHI , nParameters);
	break;
      case 4:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxAS4LowR, MINPHI,MAXPHI , nParameters);
	break;
      default:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxAS5LowR, MINPHI,MAXPHI , nParameters);
	break;
      }
    }
  } 
  if(option==1){
    if(center<1){//near side
      switch(rxnplanebin){
      case 0:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxNS0HighR, MINPHI,MAXPHI , nParameters);
	break;
      case 1:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxNS1HighR, MINPHI,MAXPHI , nParameters);
	break;
      case 2:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxNS2HighR, MINPHI,MAXPHI , nParameters);
	break;
      case 3:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxNS3HighR, MINPHI,MAXPHI , nParameters);
	break;
      case 4:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxNS4HighR, MINPHI,MAXPHI , nParameters);
	break;
      default:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxNS5HighR, MINPHI,MAXPHI , nParameters);
	break;
      }
    }
    else{
      switch(rxnplanebin){
      case 0:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxAS0HighR, MINPHI,MAXPHI , nParameters);
	break;
      case 1:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxAS1HighR, MINPHI,MAXPHI , nParameters);
	break;
      case 2:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxAS2HighR, MINPHI,MAXPHI , nParameters);
	break;
      case 3:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxAS3HighR, MINPHI,MAXPHI , nParameters);
	break;
      case 4:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxAS4HighR, MINPHI,MAXPHI , nParameters);
	break;
      default:
	fBkgdTemp = new TF1("fBkgdTemp",BackgroundxxAS5HighR, MINPHI,MAXPHI , nParameters);
	break;
      }
    }
  }

  if(!fBkgdTemp){
    cerr<<"Warning: function not set!  Cannot calculate RMS"<<endl;
    return;
  }
  Double_t *covMatrixPW = (TVirtualFitter::GetFitter())->GetCovarianceMatrix();
  const Int_t npar = fitFunc->GetNumberFreeParameters();
  Double_t   *fParams   = new Double_t[npar];
  for(int i=0;i<npar;i++){ //this uses the parameters with the correct error bars but then we need to scale up the background 
    fParams[i] = fitFunc->GetParameter(i);
    fBkgdTemp->SetParameter(i, fMyBackground->GetParameter(i));
  }   // loop over free parameters

  Double_t background = scale*fBkgdTemp->Integral(lowrange, highrange, 1E-12);   
  Double_t backgrounderr = scale*fBkgdTemp->IntegralError(lowrange, highrange, fParams, covMatrixPW, 0.01);
  //The RMS is sqrt(1/yield * (sum(content*x^2) - integral(background*x^2)))
  //The yield is really a normalization so we will treat it as errorless
  //We'll write this as R = sqrt((a-b)/c)
  //dR/da = 1/2 1/R
  //the error on this term is then
  //1/2R sqrt(errb^2 + err^2)
  RMSError = 0.5/RMS*TMath::Sqrt(backgrounderr*backgrounderr+part1err);
  //cout<<"range "<<lowrange<<" - "<<highrange<<" ";
  //cout<<"background "<<background/sum<<" +/- "<<backgrounderr/sum<<" meanxx "<<meanxx<<" +/- "<<TMath::Sqrt(part1err)/sum<<" sum "<<sum<<" sumx "<<sumx<<" sumxx "<<sumxx<<" meanx "<<meanx<<" meanxx "<<meanxx<<endl;
//   Double_t bkgdtmp = fMyBackground->Integral(lowrange, highrange, fParams, 1E-12);   
//   cout<<"part1 "<<part1<<" sqrt part1/part2 "<<TMath::Sqrt(part1/part2)<<" bkgd*x*x "<<background<<" sqrt bkgd*x*x/bkgd "<<TMath::Sqrt(background/bkgdtmp)<<endl;
//   RMS = TMath::Sqrt((part1-background)/yield);
//   RMSError = TMath::Sqrt(part1err/yield);
//  delete fBkgdTemp;
    
}

void CorrData::GetBinCenter(TH1 *histo, Int_t bin, Double_t & mean, Double_t & meanerr){
  //This function is to calculate the error on the bin center because this is actually a significant effect for the RMS
  //We want to use the average dphi in the bin but the bin has a finite width so the bin center is not the best estimate for the position
  //We estimate the distribution as a straight line y=mx+b
  //The average y (calculated by int(x y(x))/int(y(x))) is
  //[m/3(x2^3-x1^3)+b/2(x2^2-x1^2)] / [b*(x2-x1)+m/2*(x2^2-x1^2)]
  //as a check, this does reduce to the center of the bin when the slope is 0
  //We will calculate the average x from the slope and the intercept from the two neighboring bins and use the difference as the error on the average x
  Double_t bincenter1 = histo->GetBinCenter(bin);
  Double_t bincenter0 = histo->GetBinCenter(bin-1);
  Double_t bincenter2 = histo->GetBinCenter(bin+1);
  Double_t y1 = histo->GetBinContent(bin);
  Double_t y0 = histo->GetBinContent(bin-1);
  Double_t y2 = histo->GetBinContent(bin+1);
  //calculate the possible slopes and intercepts
  //slope = rise/run
  Double_t mA = (y1-y0)/(bincenter1-bincenter0);
  Double_t mB = (y1-y2)/(bincenter1-bincenter2);
  //intercept b = y - m x
  Double_t bA = y1 - mA * bincenter1;
  Double_t bB = y1 - mB * bincenter1;
  //calculate high and low edges of the bin
  Double_t x1 = histo->GetBinLowEdge(bin);
  Double_t x2 = histo->GetBinLowEdge(bin+1);
  Double_t meanXA = (mA/3.0*(x2*x2*x2-x1*x1*x1)+bA/2.0*(x2*x2-x1*x1)) / (bA*(x2-x1)+mA/2.0*(x2*x2-x1*x1));
  Double_t meanXB = (mB/3.0*(x2*x2*x2-x1*x1*x1)+bB/2.0*(x2*x2-x1*x1)) / (bB*(x2-x1)+mB/2.0*(x2*x2-x1*x1));
  mean = (meanXA+meanXB)/2.0;
  meanerr = TMath::Abs(meanXA-meanXB)/2.0;
  //cout<<"bin center is "<<bincenter1<<" calculating mean as "<<mean<<" +/- "<<meanerr<<Form("(%2.3f percent of bin width)",meanerr*100.0/histo->GetBinWidth(bin))<<" bin width is "<<histo->GetBinWidth(bin)<<endl;
}
TGraph *CorrData::GetErrorGraph(TH1 *low, TH1* high){
  TGraph *graph = new TGraph();
  //set first point as low
  Int_t nPoint = 0;
  Float_t x = low->GetBinLowEdge(1);
  Float_t y = low->GetBinContent(1);
  graph->SetPoint(0,x,y);
  nPoint++;
  for(Int_t i=1;i<=low->GetNbinsX();i++){
    x = low->GetBinCenter(i);
    y = low ->GetBinContent(i);
    graph->SetPoint(nPoint,x,y);
    nPoint++;
  }
  x = low->GetBinLowEdge(low->GetNbinsX()+1);
  y = low->GetBinContent(low->GetNbinsX()+1);
  graph->SetPoint(nPoint,x,y);
  nPoint++;
  x = high->GetBinLowEdge(high->GetNbinsX()+1);
  y = high->GetBinContent(high->GetNbinsX()+1);
  graph->SetPoint(nPoint,x,y);
  nPoint++;
  for(Int_t i=high->GetNbinsX();i>=1;i--){
    x = high->GetBinCenter(i);
    y = high->GetBinContent(i);
    graph->SetPoint(nPoint,x,y);
    nPoint++;
  }
  x = high->GetBinLowEdge(1);
  y = high->GetBinContent(1);
  graph->SetPoint(nPoint,x,y);
  nPoint++;
  x = low->GetBinLowEdge(1);
  y = low->GetBinContent(1);
  graph->SetPoint(nPoint,x,y);
  nPoint++;


  return graph;
}
