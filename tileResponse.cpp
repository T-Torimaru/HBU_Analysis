#include"langaus.C"
using namespace TMath;

void tileResponse() {
  TFile * inputFile = TFile::Open("../rootfile/testWE.root");
  TTree * tree = (TTree *) inputFile->Get("bigtree");

  Float_t pedestal[12][12] = {
    {542.002, 544.64 , 530.871, 534.821, 538.6  , 539.584, 561.556, 549.517, 541.918, 544.732, 569.881, 564.457},
    {542.213, 560.895, 535.873, 537.6  , 538.676, 539.863, 535.333, 548.768, 541.527, 562.925, 564.827, 546.642},
    {539.369, 549.664, 547.178, 537.825, 536.627, 537.584, 559.983, 579.721, 555.046, 549.801, 570.647, 554.804},
    {531.222, 515.267, 544.139, 543.698, 528.689, 536.902, 558.955, 560.388, 564.95 , 562.565, 570.078, 550.418},
    {549.339, 536.103, 529.394, 539.566, 533.816, 539.976, 559.292, 551.643, 569.486, 557.887, 563.654, 557.262},
    {534.371, 550.846, 540.311, 522.177, 550.872, 540.139, 550.238, 552.114, 564.563, 565.567, 538.668, 536.28 },
    {454.785, 458.474, 468.254, 462.743, 472.091, 452.494, 555.2  , 549.603, 549.036, 539.652, 538.437, 553.163},
    {456.756, 456.505, 455.448, 468.68 , 471.075, 462.852, 561.873, 553.023, 535.437, 556.29 , 551.017, 549.05},
    {468.168, 459.247, 464.965, 456.937, 461.345, 474.811, 556.913, 557.434, 551.73 , 549.013, 560.464, 542.963},
    {457.001, 460.302, 455.084, 474.297, 461.828, 464.308, 547.468, 532.67 , 536.547, 545.375, 550.799, 553.222},
    {466.826, 449.56 , 466.83 , 477.023, 464.621, 460.071, 538.309, 552.867, 558.133, 539.996, 555.578, 550.251},
    {461.257, 465.075, 460.383, 475.312, 450.769, 468.472, 542.936, 552.445, 552.138, 538.162, 568.533, 548.085}
  };
  Float_t gain[12][12] = {
    {36.2143, 34.8   , 34.8929, 35.25  , 34.9301, 34.7738, 36.8167, 35.9   , 36.2262, 36.1071, 36.4   , 36.8   },
    {35.8091, 35.4286, 35.4048, 35.0362, 35.2   , 35.4405, 36.7   , 35.3929, 36.4833, 36.2145, 35.958 , 36.2145},
    {35.9   , 35.8357, 34.6286, 35.2121, 34.6333, 34.65  , 36.6   , 36.7857, 35.5   , 36.2145, 36.7   , 36.2145},
    {35.6538, 34.9758, 34.1364, 34.6429, 34.5833, 34.8571, 36     , 35.6209, 35.8833, 36.1091, 36.3571, 36.7273},
    {34.7262, 36.05  , 33.5   , 35.1429, 34.3214, 34.8667, 36.3   , 36.7143, 37.1   , 35.4762, 36.7   , 36.3   },
    {34.5   , 35.5   , 34.5868, 34.0833, 34.8   , 34.631 , 35.5714, 36.2286, 35.9667, 35.3182, 36.5394, 36.1273},
    {30.4   , 31.2333, 30.6429, 30.5357, 30.4   , 29.8833, 40.2857, 39.1079, 39.4   , 38.7619, 41.0352, 38.9231},
    {30.2857, 29.7937, 29.6786, 30.1833, 30.1905, 29.8717, 39.1079, 39.7818, 39.1079, 39.3357, 39.1079, 40.1319},
    {30.4167, 30.0238, 29.8788, 29.75  , 29.8455, 29.8717, 39.1079, 39.1079, 39.1079, 38.7028, 39.4286, 39.1079},
    {29.958 , 29.9405, 29.8717, 29.7152, 30.4182, 28.8929, 39.1079, 39.1079, 39.1079, 39.1079, 39.1079, 39.1079},
    {29.4593, 29.8717, 29.5   , 29.6667, 29.8667, 30.3   , 39.1079, 39.1079, 39.1079, 39.1079, 40.0857, 38.4121},
    {29.3333, 29.4405, 29.2857, 29.7   , 29.8717, 29.1429, 39.1079, 39.1079, 39.1079, 39.1079, 39.3714, 42.5   },
  };

  Double_t reco[2][3];
  Float_t ahcADC[1024];
  Int_t ahcHitLayer[1024];
  Int_t ahcHitI[1024];
  Int_t ahcHitJ[1024];
  Float_t ahcHitPos[1024][3];
  Int_t maxPass;

  Int_t nEntries = tree->GetEntries();
  Int_t nHits;

  TCanvas *c1 = new TCanvas("c1","c1",896,900);
  TH1F *chargeH[12][12];
  TH1F *noangleH[12][12];
  TH2F *hitMap = new TH2F("hitMap","HBU light yield map in DESY;y [mm];x [mm]", 12,-180,180, 12,-180,180);

  TF1 *f1, *f2;
  Double_t peakP, peakPError;
    
  char histoName[256];
  char imgName[256];

  for (int i=0;i<12;i++) {
    for (int j=0;j<12;j++) {
        sprintf(histoName,"chargeH%d_%d",i+7,j+7);
        chargeH[i][j] = new TH1F(histoName,"light yield histogram;Npe",50,0,50);
        chargeH[i][j]->SetLineColor(1);
        sprintf(histoName,"noangleH%d_%d",i+7,j+7);
        noangleH[i][j] = new TH1F(histoName,"light yield histogram;Npe",50,0,50);
        noangleH[i][j]->SetLineColor(2);
    }
  }

  tree->SetBranchAddress("hod1_trueRecoX",&reco[0][0]);
  tree->SetBranchAddress("hod1_trueRecoY",&reco[0][1]);
  tree->SetBranchAddress("hod1_trueRecoZ",&reco[0][2]);
  tree->SetBranchAddress("hod2_trueRecoX",&reco[1][0]);
  tree->SetBranchAddress("hod2_trueRecoY",&reco[1][1]);
  tree->SetBranchAddress("hod2_trueRecoZ",&reco[1][2]);
  tree->SetBranchAddress("ahc_nHits",&nHits);
  tree->SetBranchAddress("ahc_hitI",ahcHitI);
  tree->SetBranchAddress("ahc_hitJ",ahcHitJ);
  tree->SetBranchAddress("ahc_hitK",ahcHitLayer);
  tree->SetBranchAddress("ahc_hitPos",ahcHitPos);
  tree->SetBranchAddress("ahc_hitEnergy",ahcADC);
  tree->SetBranchAddress("hod_maxPass",&maxPass);

  for (int i=0;i<nEntries;i++) {
    tree->GetEntry(i);

    if (maxPass<=1) continue;
    //if (nHits>100) continue;

    Double_t angleTrack;
    angleTrack = (reco[1][2]-reco[0][2])/Sqrt((reco[0][0]-reco[1][0])*(reco[0][0]-reco[1][0])+(reco[0][1]-reco[1][1])*(reco[0][1]-reco[1][1])+(reco[0][2]-reco[1][2])*(reco[0][2]-reco[1][2]));


    for (int n=0;n<nHits;n++) {
      if (ahcHitLayer[n]!=3) continue;
      if (ahcHitI[n]<=6 || ahcHitJ[n]<=6 || ahcHitI[n]>18 || ahcHitJ[n]>18) continue;
        Double_t newPos[2];
        newPos[0] = ((ahcHitPos[n][2]-reco[0][2])*reco[1][0]+(reco[1][2]-ahcHitPos[n][2])*reco[0][0])/(reco[1][2]-reco[0][2]);
        newPos[1] = ((ahcHitPos[n][2]-reco[0][2])*reco[1][1]+(reco[1][2]-ahcHitPos[n][2])*reco[0][1])/(reco[1][2]-reco[0][2]);

      if (ahcHitPos[n][0]>newPos[0]+25 || ahcHitPos[n][0]<newPos[0]-25
          || ahcHitPos[n][1]>newPos[1]+25 || ahcHitPos[n][1]<newPos[1]-25)
         continue;
      //if (ahcHitI[n]!=12 || ahcHitJ[n]!=12) continue;

      chargeH[ahcHitI[n]-7][ahcHitJ[n]-7]->Fill((ahcADC[n]-pedestal[ahcHitI[n]-7][ahcHitJ[n]-7])/gain[ahcHitI[n]-7][ahcHitJ[n]-7]*angleTrack);
      noangleH[ahcHitI[n]-7][ahcHitJ[n]-7]->Fill((ahcADC[n]-pedestal[ahcHitI[n]-7][ahcHitJ[n]-7])/gain[ahcHitI[n]-7][ahcHitJ[n]-7]);
      //chargeH[ahcHitI[n]-7][ahcHitJ[n]-7]->Fill((ahcADC[n]-pedestal[ahcHitI[n]-7][ahcHitJ[n]-7])*angleTrack);
    }
  }
  gStyle->SetOptStat(0);


    langaus(chargeH[5][6], &f1, &peakP, &peakPError);
    cout<<peakP<<endl;
    langaus(noangleH[5][6], &f2, &peakP, &peakPError);
    cout<<peakP<<endl;
/*
    noangleH[5][6]->Draw();
    chargeH[5][6]->Draw("same");
    f1->SetLineColor(1);
    f1->Draw("same");
    f2->SetLineColor(2);
    f2->Draw("same");
*/


  FILE *output = fopen("../txt/DESY_true.txt","w");
  for (int i=0;i<12;i++) {
    for (int j=0;j<12;j++) {
      //sprintf(imgName,"../img/npe%d_%d.png",i+7,j+7);
      //chargeH[i][j]->Draw();
      //c1->SaveAs(imgName);
      langaus(chargeH[i][j], &f1, &peakP, &peakPError);
      //hitMap->Fill((i-5.5)*30,(j-5.5)*30,chargeH[i][j]->GetMean());
      hitMap->Fill((j-5.5)*30,(i-5.5)*30,peakP);
      fprintf(output,"%d %d\t%f\n",i+7,j+7,peakP);
    }
  }
  hitMap->Draw("colztext");
  fclose(output);
}
