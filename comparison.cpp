void comparison() {

  Int_t n,m,l;
  Float_t pea,error;
  char canvasName[256];
  Double_t effa[6][14][14],eff[6][14][14];

  ifstream ifs("../txt/shift/shift_6layer_0220_align.txt");
  while(ifs>>n>>m>>l>>pea){
      effa[n][m][l]=pea;
  }
  ifs.close();

  ifstream ifa("../txt/shift/shift_6layer_0220.txt");
  while(ifa>>n>>m>>l>>pea){
      eff[n][m][l]=pea;
  }
  ifa.close();
  
  TCanvas *c1 = new TCanvas("c1","c1",696,500);
  TH1F *efficiency = new TH1F("efficiency","Alignment Calibration",25,0.8,1.);
  TH1F *efficiency_align = new TH1F("efficiency2","efficiency2",35,0.8,1.);
  //  efficiency_align->SetFillStyle(1001);
  efficiency_align->SetLineColor(2);

  for (int layer=0;layer<6;layer++) {
    for (int j=0;j<14;j++) {
      for (int k=0;k<14;k++) {
  	efficiency->Fill(eff[layer][j][k]);
  	efficiency_align->Fill(effa[layer][j][k]);	
      }
    }
  }
  //  gStyle->SetFitStat(1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  
  c1->cd();
  efficiency->Draw();
  efficiency_align->Draw("same");
    }
