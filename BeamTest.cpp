#include "langaus.C"
using namespace TMath;

void BeamTest() {

  Int_t n,m;
  Float_t pea,error;
  char canvasName[256];

  //  TFile * inputFile = TFile::Open("../rootfile/MIP/MIPcombination.root");
  TFile * inputFile = TFile::Open("../rootfile/central_200fF_beam/central_200fF.root");
  TTree * tr = (TTree *) inputFile->Get("tree");

  // Float_t ped[12][12] = {
  //   {542.002, 544.64 , 530.871, 534.821, 538.6  , 539.584, 561.556, 549.517, 541.918, 544.732, 569.881, 564.457},
  //   {542.213, 560.895, 535.873, 537.6  , 538.676, 539.863, 535.333, 548.768, 541.527, 562.925, 564.827, 546.642},
  //   {539.369, 549.664, 547.178, 537.825, 536.627, 537.584, 559.983, 579.721, 555.046, 549.801, 570.647, 554.804},
  //   {531.222, 515.267, 544.139, 543.698, 528.689, 536.902, 558.955, 560.388, 564.95 , 562.565, 570.078, 550.418},
  //   {549.339, 536.103, 529.394, 539.566, 533.816, 539.976, 559.292, 551.643, 569.486, 557.887, 563.654, 557.262},
  //   {534.371, 550.846, 540.311, 522.177, 550.872, 540.139, 550.238, 552.114, 564.563, 565.567, 538.668, 536.28 },
  //   {454.785, 458.474, 468.254, 462.743, 472.091, 452.494, 555.2  , 549.603, 549.036, 539.652, 538.437, 553.163},
  //   {456.756, 456.505, 455.448, 468.68 , 471.075, 462.852, 561.873, 553.023, 535.437, 556.29 , 551.017, 549.05},
  //   {468.168, 459.247, 464.965, 456.937, 461.345, 474.811, 556.913, 557.434, 551.73 , 549.013, 560.464, 542.963},
  //   {457.001, 460.302, 455.084, 474.297, 461.828, 464.308, 547.468, 532.67 , 536.547, 545.375, 550.799, 553.222},
  //   {466.826, 449.56 , 466.83 , 477.023, 464.621, 460.071, 538.309, 552.867, 558.133, 539.996, 555.578, 550.251},
  //   {461.257, 465.075, 460.383, 475.312, 450.769, 468.472, 542.936, 552.445, 552.138, 538.162, 568.533, 548.085}
  // };

  // Float_t ga[12][12] = {
  //   {36.2143, 34.8   , 34.8929, 35.25  , 34.9301, 34.7738, 36.8167, 35.9   , 36.2262, 36.1071, 36.4   , 36.8   },
  //   {35.8091, 35.4286, 35.4048, 35.0362, 35.2   , 35.4405, 36.7   , 35.3929, 36.4833, 36.2145, 35.958 , 36.2145},
  //   {35.9   , 35.8357, 34.6286, 35.2121, 34.6333, 34.65  , 36.6   , 36.7857, 35.5   , 36.2145, 36.7   , 36.2145},
  //   {35.6538, 34.9758, 34.1364, 34.6429, 34.5833, 34.8571, 36     , 35.6209, 35.8833, 36.1091, 36.3571, 36.7273},
  //   {34.7262, 36.05  , 33.5   , 35.1429, 34.3214, 34.8667, 36.3   , 36.7143, 37.1   , 35.4762, 36.7   , 36.3   },
  //   {34.5   , 35.5   , 34.5868, 34.0833, 34.8   , 34.631 , 35.5714, 36.2286, 35.9667, 35.3182, 36.5394, 36.1273},
  //   {30.4   , 31.2333, 30.6429, 30.5357, 30.4   , 29.8833, 40.2857, 39.1079, 39.4   , 38.7619, 41.0352, 38.9231},
  //   {30.2857, 29.7937, 29.6786, 30.1833, 30.1905, 29.8717, 39.1079, 39.7818, 39.1079, 39.3357, 39.1079, 40.1319},
  //   {30.4167, 30.0238, 29.8788, 29.75  , 29.8455, 29.8717, 39.1079, 39.1079, 39.1079, 38.7028, 39.4286, 39.1079},
  //   {29.958 , 29.9405, 29.8717, 29.7152, 30.4182, 28.8929, 39.1079, 39.1079, 39.1079, 39.1079, 39.1079, 39.1079},
  //   {29.4593, 29.8717, 29.5   , 29.6667, 29.8667, 30.3   , 39.1079, 39.1079, 39.1079, 39.1079, 40.0857, 38.4121},
  //   {29.3333, 29.4405, 29.2857, 29.7   , 29.8717, 29.1429, 39.1079, 39.1079, 39.1079, 39.1079, 39.3714, 42.5   },
  // };

  Float_t pedestal[4][36];
  Float_t gain[4][36];

  // for (int i=0;i<4;i++) {
  //   for (int j=0;j<36;j++){
  //     if (i==0 || i==1){
  // 	//	pedestal[i][j] = ped[5-int(j/6)][11*(1-i)+5*i-j%6];
  // 	gain[i][j] = ga[5-int(j/6)][11*(1-i)+5*i-j%6];
  // 	// pedestal[i][j] = ped[11*(1-i)+6*i-j%6][6+int(j/6)];
  // 	// gain[i][j] = ga[11*(1-i)+6*i-j%6][6+int(j/6)];
  //     }
  //     else {
  // 	//	pedestal[i][j] = ped[11-int(j/6)][6*(3-i)+j%6];
  // 	gain[i][j] = ga[11-int(j/6)][6*(3-i)+j%6];
  // 	// pedestal[i][j] = ped[6*(3-i)+j%6][int(j/6)];
  // 	// gain[i][j] = ga[6*(3-i)+j%6][int(j/6)];
  //     }
  //   }
  // }

  //  cout<<gain[4][35]<<endl;

  ifstream ifs("../txt/testBeam_pedestal_calib.txt");
  //  ifstream ifs("../txt/pedestal_MIP_all2.txt");
  while(ifs>>n>>m>>pea>>error){
    // if (n==132 || n==135 || n==138 || n==141){
    //   n=n-132;
    //   n=n/3;
      pedestal[n][m]=pea;
  //   }
  }
  ifs.close();
  
  //  ifstream ifa("../txt/testBeam_Allchannel_gain600fF_chip138.txt");
  ifstream ifa("../txt/testBeam_gain_138.txt");
  while(ifa>>n>>m>>pea){
    // if (n==132 || n==135 || n==138 || n==141){
    //   n=n-132;
    //   n=n/3;
      gain[n][m]=pea;
      //    }
  }
  ifa.close();
  // for (int i=0;i<36;i++) {
  // gain[0][i]=15.5155;
  // gain[1][i]=15.087;
  // gain[2][i]=23.0949;
  // gain[3][i]=12.801;
  // }

  
  Int_t adc;
  Int_t chipID;
  Int_t channel;
  Int_t HBit;
  
  Int_t nEntries = tr->GetEntries();
  Int_t k;
  
  TCanvas *c1 = new TCanvas("c1","c1",696,500);
  TCanvas *c2 = new TCanvas("c2","c2",896,900);
  TH1F *chargeH[4][36];
  TH2F *hitMap = new TH2F("hitMap","HBU light yield map in DESY;y [mm];x [mm]",12,-180,180,12,-180,180);
  hitMap->SetMaximum(20.);
  hitMap->SetMinimum(6.3);
  
  TF1 *f1;
  Double_t peakP, peakError;
  
  char histoName[256];
  char imgName[256];
  
  for (int i=0;i<4;i++) {
    for (int j=0;j<36;j++) {
      sprintf(histoName,"chargeH%d_%d",132+3*i,j+1);
      chargeH[i][j] = new TH1F(histoName,"light yield histogram;Npe",50,0,50);
      chargeH[i][j]->SetLineColor(1);
    }
  }
  
  tr->SetBranchAddress("ChipID2", &chipID);
  tr->SetBranchAddress("chn", &channel);
  tr->SetBranchAddress("ADC", &adc);
  tr->SetBranchAddress("Hit_Bit", &HBit);
  
  for (int i=0;i<nEntries;i++) {
    tr->GetEntry(i);
    
    if (chipID==132 || chipID==135 || chipID==138 || chipID==141) { 
      if (HBit!=1) continue;
      k = chipID-132;
      k = k/3;
      chargeH[k][channel]->Fill((adc-pedestal[k][channel])/gain[k][channel]);
    }
  }
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  
  // chargeH[2][18]->Draw();
  // langaus(chargeH[2][18], &f1, &peakP, &peakError);
  // cout<<peakP<<endl;
  // langaus(chargeH[0][5], &f2, &peakP, &peakError);
  // cout<<peakP<<endl;
  // c1->cd();
  // f1->Draw("same");
  

  FILE *output = fopen("../txt/testBeamresults_200fF_central_138.txt","w");
  for (int i=0;i<4;i++) {
    for (int j=0;j<36;j++) {
      langaus(chargeH[i][j], &f1, &peakP, &peakError);
      if (i==0 || i==1) {
  	hitMap->Fill(165*(1-i)-15*i-30*(j%6),-15-30*int(j/6),peakP);
  	fprintf(output,"%d %d\t%f\n", i, j, peakP);
      }
      else {
  	hitMap->Fill(-165*(i-2)+15*(3-i)+30*(j%6),165-30*int(j/6),peakP);
  	fprintf(output,"%d %d\t%f\n", i, j, peakP);
      }
      c1->cd();
      chargeH[i][j]->Draw();
      f1->Draw("same");
      sprintf(canvasName,"../png/12_12_600fF/chip138_fit%d_%d.png",i,j);
      c1->SaveAs(canvasName);
    }
  }
  c2->cd();
  hitMap->Draw("colztext");
  fclose(output);
}
