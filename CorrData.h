#ifndef ROOT_CorrData
#define ROOT_CorrData

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// CorrData                                                             //
//                                                                      //
// Description of the event and track parameters                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include "TH1.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"

class CorrData : public TObject {

private:
  static Int_t nBinsInRxnPlane;
  TH1F *hSignalPlusBkgd[6];//signal+background, data
  TH1F *hSignal[6];//signal, calculated
  TH1F *hBkgd[6];//background dominated region
  TH1F *hCalcBkgd[6];//calculated background with errors
  TH1F *hCalcBkgdInSignalPlusBackgroundRegion[6];//calculated background with errors
  TF1 *fBkgd[6];//fit background function in background dominated region
  TF1 *fBkgdInSignalPlusBkgdRegion[6];//fit background in signal+background region

  TH1F *hSignalLowScale[6];//signal, calculated
  TH1F *hCalcBkgdInSignalPlusBackgroundRegionLowScale[6];//calculated background with errors
  TF1 *fBkgdInSignalPlusBkgdRegionLowScale[6];//fit background in signal+background region

  TH1F *hSignalHighScale[6];//signal, calculated
  TH1F *hCalcBkgdInSignalPlusBackgroundRegionHighScale[6];//calculated background with errors
  TF1 *fBkgdInSignalPlusBkgdRegionHighScale[6];//fit background in signal+background region
  TH1F *hSignalLowR[6];//signal, calculated
  TH1F *hCalcBkgdLowR[6];//calculated background with errors
  TH1F *hCalcBkgdInSignalPlusBackgroundRegionLowR[6];//calculated background with errors
  TF1 *fBkgdLowR[6];//fit background function in background dominated region
  TF1 *fBkgdInSignalPlusBkgdRegionLowR[6];//fit background in signal+background region
  TH1F *hSignalHighR[6];//signal, calculated
  TH1F *hCalcBkgdHighR[6];//calculated background with errors
  TH1F *hCalcBkgdInSignalPlusBackgroundRegionHighR[6];//calculated background with errors
  TF1 *fBkgdHighR[6];//fit background function in background dominated region
  TF1 *fBkgdInSignalPlusBkgdRegionHighR[6];//fit background in signal+background region
  TH1F *hCommonBkgd;
  TF1 *fitFunc;
  TF1 *fitFuncLowR;
  TF1 *fitFuncHighR;
  Float_t Rn[5];
  Float_t nsYield[6];
  Float_t asYield[6];
  Double_t nsRMS[6];
  Double_t asRMS[6];
  Float_t nsYieldError[6];
  Float_t asYieldError[6];
  Double_t nsRMSError[6];
  Double_t asRMSError[6];
  Float_t nsYieldLowScale[6];
  Float_t asYieldLowScale[6];
  Double_t nsRMSLowScale[6];
  Double_t asRMSLowScale[6];
  Float_t nsYieldErrorLowScale[6];
  Float_t asYieldErrorLowScale[6];
  Double_t nsRMSErrorLowScale[6];
  Double_t asRMSErrorLowScale[6];
  Float_t nsYieldHighScale[6];
  Float_t asYieldHighScale[6];
  Double_t nsRMSHighScale[6];
  Double_t asRMSHighScale[6];
  Float_t nsYieldErrorHighScale[6];
  Float_t asYieldErrorHighScale[6];
  Double_t nsRMSErrorHighScale[6];
  Double_t asRMSErrorHighScale[6];
  Float_t nsYieldLowR[6];
  Float_t asYieldLowR[6];
  Double_t nsRMSLowR[6];
  Double_t asRMSLowR[6];
  Float_t nsYieldErrorLowR[6];
  Float_t asYieldErrorLowR[6];
  Double_t nsRMSErrorLowR[6];
  Double_t asRMSErrorLowR[6];
  Float_t nsYieldHighR[6];
  Float_t asYieldHighR[6];
  Double_t nsRMSHighR[6];
  Double_t asRMSHighR[6];
  Float_t nsYieldErrorHighR[6];
  Float_t asYieldErrorHighR[6];
  Double_t nsRMSErrorHighR[6];
  Double_t asRMSErrorHighR[6];
  //nominal phi cuts for yield calculations
  Float_t nsPhiCut;
  Float_t asPhiCut;//from pi
  //real cuts using bin widths
  Float_t nsRealPhiCutLow;
  Float_t nsRealPhiCutHigh;
  Float_t asRealPhiCutLow;
  Float_t asRealPhiCutHigh;
  static Double_t R2;
  static Double_t R4;
  static Double_t R6;
  static Double_t R8;
  static Double_t R10;
  static Double_t R2FracErr;
  static Double_t R4FracErr;
  static Double_t R6FracErr;
  static  Double_t phiRange;
  static  Double_t MINPHI;
  static  Double_t MAXPHI;
  static Double_t IntegralTolerance;
  Int_t nbinsBkgd;
  Float_t phiRangeBkgd;
  //Constants so that code runs more efficiently
  static Float_t mypi;
  static Float_t sqrt3;
  static Double_t ConvertPhiToCommonPhi(Int_t bin,Double_t phi);
  static Double_t ConvertCommonPhiToPhi(Int_t bin,Double_t phi);
  Int_t CalculateReactionPlaneBin(Float_t commonphi);//input phi from common background histogram, get out the bin it's in relative to the reaction plane
  Double_t etaRangeSigToBkgdRatio;//Ratio of eta range in Signal+Background to eta range in Background dominated range
  Double_t etaRangeSigToBkgdRatioLow;//Ratio of eta range in Signal+Background to eta range in Background dominated range
  Double_t etaRangeSigToBkgdRatioHigh;//Ratio of eta range in Signal+Background to eta range in Background dominated range
  TCanvas *cBackground;
  TCanvas *cSignalPlusBackground;
  //The following functions are used for propagation of errors.  They are largely taken from TF1.cxx but modified for the special case we have here.
  //Here we want to be able to propagate errors from one fit function to another function.  These functions therefore needed to be modified to take a function used for a fit which is used for the error propagation.
  void GradientParPW(const Double_t x, Double_t *grad, TF1 *inputFuncForErrors,   Double_t *fParams, Double_t eps = 0.01);
  Double_t GradientParPW(Int_t ipar, const Double_t x, Double_t eps, TF1 *inputFuncForErrors,   Double_t *fParams);
  Double_t DERPW(Int_t ipar, const Double_t x, Double_t eps, TF1 *inputFuncForErrors,   Double_t *fParams);
  void DerivPW(const Double_t x, Double_t *gradX, TF1 *inputFuncForErrors,   Double_t *fParams, Double_t eps = 0.01);
  //input histo is used to get binning correct
  TH1F *GetBackgroundErrorPW(TH1F *inputHisto, TF1 *inputFunc);
  void CalculateNetYield(TH1 *histo, TF1 *inputFunc, Int_t rxnPlaneBin, Int_t option);
  void CalculateRMS(TH1 *histo, TF1 *fMyBackground,Double_t lowrange, Double_t highrange, Double_t center,Int_t rxnplanebin, Double_t & RMS, Double_t & RMSError, Int_t option);
  void GetBinCenter(TH1 *histo, Int_t bin, Double_t & mean, Double_t & meanerr);//This function is to calculate the error on the bin center because this is actually a significant effect for the RMS
  TGraph *GetErrorGraph(TH1 *low, TH1* high);
  Int_t nParameters;

public:
  CorrData();
  virtual ~CorrData();
  TH1F *GetSignalPlusBackgroundHistogram(Int_t bin);
  TH1F *GetSignalHistogram(Int_t bin);
  TH1F *GetBackgroundHistogram(Int_t bin);
  int *SetSignalPlusBackgroundHistogram(Int_t bin, TH1F *histo){delete hSignalPlusBkgd[bin];hSignalPlusBkgd[bin] = histo;return 0;};
  void SetSignalHistogram(Int_t bin, TH1F *histo);//{delete hSignal[bin];hSignal[bin] = histo;};
  void SetBackgroundHistogram(Int_t bin, TH1F *histo);//{delete hBkgd[bin];hBkgd[bin] = histo;};
  TH1F *GetCommonBackgroundHistogram(){return hCommonBkgd;}
  static Double_t R(Int_t n);
  void SetR(Int_t n, Float_t R);
  void SetR(Float_t R2, Float_t R4, Float_t R6, Float_t R8, Float_t R10);
  void SetPhiFitRange(Float_t range){phiRange = 1.0;}
  void CreateCommonBackgroundHisto(char *name);
  void SetEtaRangeSignalToBackgroundRatio(Double_t valueavg, Double_t valuelow, Double_t valuehigh){etaRangeSigToBkgdRatio = valueavg;etaRangeSigToBkgdRatioLow = valuelow;etaRangeSigToBkgdRatioHigh = valuehigh;}
  static Double_t vntrig(Int_t n, Int_t rxnplanebin, Double_t v2, Double_t v4, Int_t rangeoption);
  void InitializeCommonFitFunction();
  static Double_t Background(const double * x, const double * p, Int_t rxnplanebin);
  static Double_t BackgroundGeneral(const double * x, const double * p, Int_t rxnplanebin);
  static Double_t BackgroundCombined(const double * x, const double * p);
  static Double_t BackgroundLowR(const double * x, const double * p, Int_t rxnplanebin);
  static Double_t BackgroundGeneralLowR(const double * x, const double * p, Int_t rxnplanebin);
  static Double_t BackgroundCombinedLowR(const double * x, const double * p);
  static Double_t BackgroundHighR(const double * x, const double * p, Int_t rxnplanebin);
  static Double_t BackgroundGeneralHighR(const double * x, const double * p, Int_t rxnplanebin);
  static Double_t BackgroundCombinedHighR(const double * x, const double * p);
  static Double_t Background0(const double * x, const double * p){return BackgroundGeneral(x,p,0);};
  static Double_t Background1(const double * x, const double * p){return BackgroundGeneral(x,p,1);};
  static Double_t Background2(const double * x, const double * p){return BackgroundGeneral(x,p,2);};
  static Double_t Background3(const double * x, const double * p){return BackgroundGeneral(x,p,3);};
  static Double_t Background4(const double * x, const double * p){return BackgroundGeneral(x,p,4);};
  static Double_t Background5(const double * x, const double * p){return BackgroundGeneral(x,p,5);};
  static Double_t BackgroundxxNS0(const double * x, const double * p){return x[0]*x[0]*BackgroundGeneral(x,p,0);};
  static Double_t BackgroundxxNS1(const double * x, const double * p){return x[0]*x[0]*BackgroundGeneral(x,p,1);};
  static Double_t BackgroundxxNS2(const double * x, const double * p){return x[0]*x[0]*BackgroundGeneral(x,p,2);};
  static Double_t BackgroundxxNS3(const double * x, const double * p){return x[0]*x[0]*BackgroundGeneral(x,p,3);};
  static Double_t BackgroundxxNS4(const double * x, const double * p){return x[0]*x[0]*BackgroundGeneral(x,p,4);};
  static Double_t BackgroundxxNS5(const double * x, const double * p){return x[0]*x[0]*BackgroundGeneral(x,p,5);};
  static Double_t BackgroundxxAS0(const double * x, const double * p){return (x[0]-mypi)*(x[0]-mypi)*BackgroundGeneral(x,p,0);};
  static Double_t BackgroundxxAS1(const double * x, const double * p){return (x[0]-mypi)*(x[0]-mypi)*BackgroundGeneral(x,p,1);};
  static Double_t BackgroundxxAS2(const double * x, const double * p){return (x[0]-mypi)*(x[0]-mypi)*BackgroundGeneral(x,p,2);};
  static Double_t BackgroundxxAS3(const double * x, const double * p){return (x[0]-mypi)*(x[0]-mypi)*BackgroundGeneral(x,p,3);};
  static Double_t BackgroundxxAS4(const double * x, const double * p){return (x[0]-mypi)*(x[0]-mypi)*BackgroundGeneral(x,p,4);};
  static Double_t BackgroundxxAS5(const double * x, const double * p){return (x[0]-mypi)*(x[0]-mypi)*BackgroundGeneral(x,p,5);};
  static Double_t Background0LowR(const double * x, const double * p){return BackgroundGeneralLowR(x,p,0);};
  static Double_t Background1LowR(const double * x, const double * p){return BackgroundGeneralLowR(x,p,1);};
  static Double_t Background2LowR(const double * x, const double * p){return BackgroundGeneralLowR(x,p,2);};
  static Double_t Background3LowR(const double * x, const double * p){return BackgroundGeneralLowR(x,p,3);};
  static Double_t Background4LowR(const double * x, const double * p){return BackgroundGeneralLowR(x,p,4);};
  static Double_t Background5LowR(const double * x, const double * p){return BackgroundGeneralLowR(x,p,5);};
  static Double_t BackgroundxxNS0LowR(const double * x, const double * p){return x[0]*x[0]*BackgroundGeneralLowR(x,p,0);};
  static Double_t BackgroundxxNS1LowR(const double * x, const double * p){return x[0]*x[0]*BackgroundGeneralLowR(x,p,1);};
  static Double_t BackgroundxxNS2LowR(const double * x, const double * p){return x[0]*x[0]*BackgroundGeneralLowR(x,p,2);};
  static Double_t BackgroundxxNS3LowR(const double * x, const double * p){return x[0]*x[0]*BackgroundGeneralLowR(x,p,3);};
  static Double_t BackgroundxxNS4LowR(const double * x, const double * p){return x[0]*x[0]*BackgroundGeneralLowR(x,p,4);};
  static Double_t BackgroundxxNS5LowR(const double * x, const double * p){return x[0]*x[0]*BackgroundGeneralLowR(x,p,5);};
  static Double_t BackgroundxxAS0LowR(const double * x, const double * p){return (x[0]-mypi)*(x[0]-mypi)*BackgroundGeneralLowR(x,p,0);};
  static Double_t BackgroundxxAS1LowR(const double * x, const double * p){return (x[0]-mypi)*(x[0]-mypi)*BackgroundGeneralLowR(x,p,1);};
  static Double_t BackgroundxxAS2LowR(const double * x, const double * p){return (x[0]-mypi)*(x[0]-mypi)*BackgroundGeneralLowR(x,p,2);};
  static Double_t BackgroundxxAS3LowR(const double * x, const double * p){return (x[0]-mypi)*(x[0]-mypi)*BackgroundGeneralLowR(x,p,3);};
  static Double_t BackgroundxxAS4LowR(const double * x, const double * p){return (x[0]-mypi)*(x[0]-mypi)*BackgroundGeneralLowR(x,p,4);};
  static Double_t BackgroundxxAS5LowR(const double * x, const double * p){return (x[0]-mypi)*(x[0]-mypi)*BackgroundGeneralLowR(x,p,5);};
  static Double_t Background0HighR(const double * x, const double * p){return BackgroundGeneralHighR(x,p,0);};
  static Double_t Background1HighR(const double * x, const double * p){return BackgroundGeneralHighR(x,p,1);};
  static Double_t Background2HighR(const double * x, const double * p){return BackgroundGeneralHighR(x,p,2);};
  static Double_t Background3HighR(const double * x, const double * p){return BackgroundGeneralHighR(x,p,3);};
  static Double_t Background4HighR(const double * x, const double * p){return BackgroundGeneralHighR(x,p,4);};
  static Double_t Background5HighR(const double * x, const double * p){return BackgroundGeneralHighR(x,p,5);};
  static Double_t BackgroundxxNS0HighR(const double * x, const double * p){return x[0]*x[0]*BackgroundGeneralHighR(x,p,0);};
  static Double_t BackgroundxxNS1HighR(const double * x, const double * p){return x[0]*x[0]*BackgroundGeneralHighR(x,p,1);};
  static Double_t BackgroundxxNS2HighR(const double * x, const double * p){return x[0]*x[0]*BackgroundGeneralHighR(x,p,2);};
  static Double_t BackgroundxxNS3HighR(const double * x, const double * p){return x[0]*x[0]*BackgroundGeneralHighR(x,p,3);};
  static Double_t BackgroundxxNS4HighR(const double * x, const double * p){return x[0]*x[0]*BackgroundGeneralHighR(x,p,4);};
  static Double_t BackgroundxxNS5HighR(const double * x, const double * p){return x[0]*x[0]*BackgroundGeneralHighR(x,p,5);};
  static Double_t BackgroundxxAS0HighR(const double * x, const double * p){return (x[0]-mypi)*(x[0]-mypi)*BackgroundGeneralHighR(x,p,0);};
  static Double_t BackgroundxxAS1HighR(const double * x, const double * p){return (x[0]-mypi)*(x[0]-mypi)*BackgroundGeneralHighR(x,p,1);};
  static Double_t BackgroundxxAS2HighR(const double * x, const double * p){return (x[0]-mypi)*(x[0]-mypi)*BackgroundGeneralHighR(x,p,2);};
  static Double_t BackgroundxxAS3HighR(const double * x, const double * p){return (x[0]-mypi)*(x[0]-mypi)*BackgroundGeneralHighR(x,p,3);};
  static Double_t BackgroundxxAS4HighR(const double * x, const double * p){return (x[0]-mypi)*(x[0]-mypi)*BackgroundGeneralHighR(x,p,4);};
  static Double_t BackgroundxxAS5HighR(const double * x, const double * p){return (x[0]-mypi)*(x[0]-mypi)*BackgroundGeneralHighR(x,p,5);};
  TF1 *GetCommonBackgroundFunction(){return fitFunc;}
  void FitCommonBackgroundHistogram();
  void DrawSignalPlusBackground();
  void DrawBackground();
  void DrawSignal();
  void MakeYieldTable(char *filename);
  void SetNParameters(Int_t nPar){nParameters = nPar;}

   
   ClassDef(CorrData,1)  //Event structure
};

#endif




















