#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include <algorithm>
 
const Int_t np = 6;
Int_t npeaks;
Double_t fpeaks(Double_t *x, Double_t *par) {
  Double_t result = 0;
  for (Int_t p=0;p<npeaks;p++) {
    Double_t norm  = par[3*p];
    Double_t mean  = par[3*p+1];
    Double_t sigma = par[3*p+2];
    result += norm*TMath::Gaus(x[0],mean,sigma);
  }

  return result;
}


void peaksYY1(Int_t nadc=192) {

  // Int_t np=3;
  TFile *f = new TFile("/home/iris/anaIris/root_files/output01530.root");
  npeaks = np;
  Char_t var[10];
  sprintf(var,"adc%i",nadc);
  TH1F *h = gFile->FindObjectAny(var); //adc number
  if (nadc>191 && nadc<320)  
    h->GetXaxis()->SetRange(200,4096); // ignoring pedestal for YY1
  else if (nadc<192)
  h->GetXaxis()->SetRange(60,4096); //ignoring pedestal for s3
  //  TH1F *h = new TH1F("h","test",500,0,1000);
  //generate n peaks at random
    Double_t par[30];
   Int_t p;

  TCanvas *c1 = new TCanvas("c1","c1",10,10,1000,900);
  c1->Divide(1,3);
  c1->cd(1);
  
  TH1F *h2 = (TH1F*)h->Clone("h2");
  
  h->Draw();
  //Use TSpectrum to find the peak candidates
 Double_t inSig = 1.5.; // sigma of the peaks to be founds 
  TSpectrum *s = new TSpectrum(npeaks);
  Int_t nfound = s->Search(h,inSig,"new");
  printf("Found %d candidate peaks to fit\n",nfound);
  c1->Update();
  c1->cd(2);
 
  npeaks = 0;
  Float_t *xpeaks = s->GetPositionX();
  for (p=0;p<nfound;p++) {
    Float_t xp = xpeaks[p];
    Int_t bin = h->GetXaxis()->FindBin(xp);
    Float_t yp = h->GetBinContent(bin);
 
    par[3*npeaks] = yp;
    par[3*npeaks+1] = xp;
    par[3*npeaks+2] = xp/100.;
    npeaks++;
  }
  
printf("Found %d useful peaks to fitn",npeaks);
  printf("Now fitting: Be patient\n");
  TF1 *fit = new TF1("fit",fpeaks,0,1000,3*npeaks); //Function
  TVirtualFitter::Fitter(h2,3*npeaks);
  fit->SetParameters(par);
  fit->SetNpx(1000);
  h2->Fit("fit"); 

  Int_t xlow = 2*par[7]/3;
  Int_t xhigh = 4*par[7]/3;

  h->GetXaxis()->SetRange(xlow,xhigh);
  h2->GetXaxis()->SetRange(xlow,xhigh);
  
// Ordering the peaks in order of ascension.
  Int_t roi[np]={1,4,7,10,13,16}; //Re Ordering Index. used for reordering fit parameters in ascending order

  Double_t rom[np];//Re Ordered Mean.
  
  for (Int_t i = 0; i <npeaks; i++){
    rom[i] = fit->GetParameter(roi[i]);
} 
  // Use swap function from algorithm library.
 
  for(Int_t i = 0; i< np;i++) {
    for (Int_t j = 0; j< i;j++){ 
   if (rom[i] < rom[j]){
     swap(rom[i],rom[j]); // swap array values
     swap(roi[i],roi[j]); // swap array index (array for the order of the rom array)
   }//if
  }//for
  }//for

//Linear fit using triple-alpha peaks.(README)

  TF1 *fitlin = new TF1("fitlin","pol1",60,120);
  Double_t xfit[np];
   Double_t exfit[np];

   Double_t yfit[np]  = {5.1443,5.15659,5.4428,5.48556, 5.76264, 5.80477}; //Triple Alpha source energies, in MeV
   Double_t eyfit[np]  = {0.005,0.005,0.005, 0.005,0.005,0.005}; // errors

   for(Int_t i = 0; i< np;i++) {
     xfit[i] = fit->GetParameter(roi[i]);
     exfit[i] = abs(fit->GetParameter(roi[i]+1));
 
}

  
 c1->cd(3);
   TGraphErrors *gr = new TGraphErrors(6,xfit, yfit, exfit,  eyfit);
   gr->Fit("fitlin");

gr->SetTitle("Si calibration");
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);
   gr->Draw("AP");

   c1->Update();

   //write the calibration parameters to a file

   ofstream ofile;

   ofile.open("alphaResolYY11480.txt",ios::app);
   ofile << nadc << "\t"<< 2.35*5794.83*fit->GetParameter(8)/fit->GetParameter(7) << endl;
   // ofile << nadc << "\t"<< fitlin->GetParameter(0) << "\t"<< fitlin->GetParameter(1) << endl;
   ofile.close();
           
}
