#include "CorrData.h"
//#include "CorrData.cxx"
//class CorrData;
//set ranges for histograms
Float_t MINPHI = -TMath::Pi()/2.0;
Float_t MAXPHI = TMath::Pi()*3.0/2.0;
Float_t shiftPhi(Float_t phi){
  if(phi<MINPHI){
    phi += 2.0*TMath::Pi();
  }
  if(phi>MAXPHI){
    phi -= 2.0*TMath::Pi();
  }
  return phi;
}
void ReadSTARData(char *filenamesignalplusbkgd,char *filenamebkgd,char *histoname,CorrData *mydata){
  Float_t binlowedge;
  Float_t val[6] = {0,0,0,0,0,0};
  Float_t err[6] = {0,0,0,0,0,0};
  Int_t nbins = 48;
  TH1F *hSignalPlusBkgd[6];
  TH1F *hBkgd[6];
  TH1F *hSignal[6];
  char myhistoname[200];
  for(int i=0;i<6;i++){
    sprintf(myhistoname,"SigPlusBkgd%sBin%i",histoname,i);
    hSignalPlusBkgd[i] = new TH1F(myhistoname,myhistoname,nbins,MINPHI,MAXPHI);
    sprintf(myhistoname,"Bkgd%sBin%i",histoname,i);
    hBkgd[i] = new TH1F(myhistoname,myhistoname,nbins,MINPHI,MAXPHI);
    sprintf(myhistoname,"Sig%sBin%i",histoname,i);
    hSignal[i] = new TH1F(myhistoname,myhistoname,nbins,MINPHI,MAXPHI);
  }
  ifstream filesignalplusbkgd(filenamesignalplusbkgd);
  Int_t bin = 1;
  while(!filesignalplusbkgd.eof())
    {
      filesignalplusbkgd>>binlowedge;
      for(int i=0;i<6;i++){
	filesignalplusbkgd>>val[i]>>err[i];
      }
  Int_t thisbin;
      if(bin<=nbins){//number of bins I should read in, avoids weird bug where it actually goes one past the number of lines.  Could write loop better but why?
 	binlowedge = shiftPhi(binlowedge);//shift the bin low edge to the ones in our coordinates
	thisbin = hSignalPlusBkgd[0]->FindBin(binlowedge+1e-3);
      }
      for(int i=0;i<6;i++){
	hSignalPlusBkgd[i]->SetBinContent(thisbin,val[i]);
	hSignalPlusBkgd[i]->SetBinError(thisbin,err[i]);
	hSignal[i]->SetBinContent(thisbin,val[i]);
	hSignal[i]->SetBinError(thisbin,err[i]);
      }
      bin++;
    }
  filesignalplusbkgd.close();
  ifstream filebkgd(filenamebkgd);
  bin = 1;
  while(!filebkgd.eof())
    {
      filebkgd>>binlowedge;
      for(int i=0;i<6;i++){
	filebkgd>>val[i]>>err[i];
      }
  Int_t thisbin;
      if(bin<=nbins){//number of bins I should read in, avoids weird bug where it actually goes one past the number of lines.  Could write loop better but why?
 	binlowedge = shiftPhi(binlowedge);//shift the bin low edge to the ones in our coordinates
	thisbin = hBkgd[0]->FindBin(binlowedge+1e-3);
      }
      for(int i=0;i<6;i++){
	hBkgd[i]->SetBinContent(thisbin,val[i]);
	hBkgd[i]->SetBinError(thisbin,err[i]);
      }
      bin++;
    }
  filebkgd.close();
  for(int i=0;i<6;i++){
    mydata->SetSignalPlusBackgroundHistogram(i,hSignalPlusBkgd[i]);
    mydata->SetSignalHistogram(i,hSignal[i]);
    mydata->SetBackgroundHistogram(i,hBkgd[i]);
  }

}

Float_t GetPtAssocLow(Int_t assocBin){//0.15-0.5, 0.5-1.0, 1.0-1.5, 1.5-2.0, 2.0-3.0, first bin rounded to 0.0 to make logistics easier
  switch(assocBin){
  case 0:
    return 0.0;
  case 1:
    return 0.5;
  case 2:
    return 1.0;
  case 3:
    return 1.5;
  case 4:
    return 2.0;
  case 5://define an extra bin so we can be lazy for the next function...
    return 3.0;
  case 6:
    return 4.0;
  }
  return 0.0;
}
Float_t GetPtAssocHigh(Int_t assocBin){
  return GetPtAssocLow(assocBin+1);
}
Float_t GetPtTrigLow(Int_t trigBin){
  if(trigBin==0) return 3.0;
  return 4.0;
}
Float_t GetPtTrigHigh(Int_t trigBin){
  if(trigBin==0) return 4.0;
  return 6.0;
}

void FitSTARData(Int_t trigBin = 1, Int_t assocBin = 4){
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

   gROOT->ProcessLine(".L CorrData.cxx++");
   gSystem->Load("libPhysics");
   //gROOT->ProcessLine(".L BackgroundFunction.cxx+");
   //gROOT->ProcessLine(".L BackgroundFunction.C+");
   //gROOT->LoadMacro("BackgroundFunction.C");
   gROOT->ProcessLine(".L CorrData.cxx++");
    // gROOT->LoadLibrary(CorrData.)
   CorrData *dataAssoc2Trig4 = new CorrData();
   //this line reads in the data
   Float_t assocLow = GetPtAssocLow(assocBin);
   Float_t assocHigh = GetPtAssocHigh(assocBin);
   Float_t trigLow = GetPtTrigLow(trigBin);
   Float_t trigHigh = GetPtTrigHigh(trigBin);

   char canvasname[200];
   char canvastitle[200];
   char backgroundhistoname[200];
   char datafile1[200];
   char datafile2[200];
   char histotag[200];
   char outputfilename[200];
   sprintf(canvasname,"canvasAssoc%2.1fTo%2.1fTrig%2.1fTo%2.1f",assocLow,assocHigh,trigLow,trigHigh);
   sprintf(canvastitle,"Associated %2.1f To %2.1f Trig %2.1f To %2.1f",assocLow,assocHigh,trigLow,trigHigh);
   sprintf(datafile1,"datafiles/SigPlusBkgdAssoc%2.1fTo%2.1fTrig%2.1fTo%2.1f.dat",assocLow,assocHigh,trigLow,trigHigh);
   sprintf(datafile2,"datafiles/BkgdAssoc%2.1fTo%2.1fTrig%2.1fTo%2.1f.dat",assocLow,assocHigh,trigLow,trigHigh);
   sprintf(histotag,"Assoc%2.1fTo%2.1fTrig%2.1fTo%2.1f",assocLow,assocHigh,trigLow,trigHigh);
   sprintf(backgroundhistoname,"hBkgdAssoc%2.1fTo%2.1fTrig%2.1fTo%2.1f",assocLow,assocHigh,trigLow,trigHigh);
   sprintf(outputfilename,"Assoc%2.1fTo%2.1fTrig%2.1fTo%2.1f.dat",assocLow,assocHigh,trigLow,trigHigh);

   cout<<datafile1<<endl;
   cout<<datafile2<<endl;

   TCanvas *canvasAssoc2Trig4 = new TCanvas(canvasname,canvastitle,600,400);
   //ReadSTARData("datafiles/SigPlusBkgdAssoc2.0To3.0Trig4.0To6.0.dat","datafiles/BkgdAssoc2.0To3.0Trig4.0To6.0.dat","Assoc2.0To3.0Trig4.0To6.0",dataAssoc2Trig4);
   ReadSTARData(datafile1,datafile2,histotag,dataAssoc2Trig4);
   dataAssoc2Trig4->SetEtaRangeSignalToBackgroundRatio(1.0/0.4225*1.03,1.0/0.4225*1.02,1.0/0.4225*1.04);//reasonable ranges appear to be 1.03 +/- 0.01, 100% correlated
   dataAssoc2Trig4->SetNParameters(6);
   dataAssoc2Trig4->CreateCommonBackgroundHisto(backgroundhistoname);//argument is used as histogram name
   dataAssoc2Trig4->GetCommonBackgroundHistogram()->Draw();
   dataAssoc2Trig4->InitializeCommonFitFunction();
   TF1 *fitFunc =    dataAssoc2Trig4->GetCommonBackgroundFunction();
   fitFunc->SetParLimits(5,0,1e-3);
   dataAssoc2Trig4->FitCommonBackgroundHistogram();
   //return;
   dataAssoc2Trig4->GetCommonBackgroundFunction()->Draw("same");
   dataAssoc2Trig4->DrawSignalPlusBackground();
   dataAssoc2Trig4->DrawBackground();
   dataAssoc2Trig4->DrawSignal();
   dataAssoc2Trig4->MakeYieldTable(outputfilename);
return;

}
