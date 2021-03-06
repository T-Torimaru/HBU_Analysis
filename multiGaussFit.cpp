#include "langaus.C"
using namespace TMath;

void multiGaussFit() {
  TFile * inputFile = TFile::Open("../rootfile/VCalib_5000cycles_ET_IC_200fF_5000mV__06p12p2017__22o20o24.root");
  TTree * tr = (TTree *) inputFile->Get("tree");

  //  Float_t pedestal[4][36];
  
  Int_t adc;
  Int_t chipID;
  Int_t channel;
  Int_t HBit;
  
  Int_t nEntries = tr->GetEntries();
  Int_t k;
  
  TH1F *adcHisto[4][36];
  
  TF1 *f1;
  Double_t peakP, peakError;
  
  TCanvas *pedCanvas = new TCanvas("c1","canvas",696,500);
  char histoName[256];
  char imgName[256];
  
  TH1F *gainH;
  gainH = new TH1F(histoName,"ADC distribution;ADC",500,300,800);
  gainH->SetLineColor(1);
  gainH->SetXTitle("ADC");
  gainH->SetYTitle("Event Number");
  
  tr->SetBranchAddress("ChipID2", &chipID);
  tr->SetBranchAddress("chn", &channel);
  tr->SetBranchAddress("ADC", &adc);
  tr->SetBranchAddress("Hit_Bit", &HBit);
  
  for (int i=0;i<nEntries;i++) {
    tr->GetEntry(i);
    
    if (chipID==132 || chipID==135 || chipID==138 || chipID==141) { 
      //      if (HBit!=1) continue;
    k = chipID-132;
    k = k/3;
    gainH->Fill(adc);
    }
  }
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  
  TF1 *Gaus = new TF1("Gaus", "[3]*exp(-((x-[0])/[2])^2)+[5]*exp(-((x-[0]-[1])/[4])^2)+[7]*exp(-((x-[0]-[1]*2)/[6])^2)",300,800);
  Gaus->SetParameter(0,500);
  Gaus->SetParameter(1,50);
  Gaus->SetParameter(2,1000);

  for (int i=0;i<4;i++) {
    for (int j=0;j<36;j++) {
      adcHisto[i][j]->Fit(Gaus, "NQ", "", 300, 800);
      peakP = Gaus->GetParameter(0);
      peakError = Gaus->GetParameter(1);      
      fprintf(output,"%d %d %f %f\n", i, j, peakP, peakError);
    }
  }
 
  // for (int i=0;i<4;i++) {
  //   for (int j=0;j<36;j++) {
      adcHisto[1][0]->Draw();
      adcHisto[1][0]->Fit(Gaus,"NQ","",300,800);
      Gaus->Draw("same");
      //      cout>>Gaus->GetParameter(0)>>endl;
      //    }}
}
