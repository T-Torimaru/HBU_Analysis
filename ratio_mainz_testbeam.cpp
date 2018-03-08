using namespace TMath;

void ratio_mainz_testbeam() {

  Int_t n,m;

  Float_t pea;
  Float_t beam[4][36];
  Float_t mainz[4][36];

  ifstream ifs("../txt/testBeamresults_200fF_central_138.txt");
  while(ifs>>n>>m>>pea){
    beam[n][m]=pea;
  }
  ifs.close();

  ifstream aaa("../txt/testBeamAllMIPcentral600fF_chip138.txt");
  while(aaa>>n>>m>>pea){
    //    if (n==3 || n==6 || n==9 || n==12){
      // n=n/3;
      // n=n-1;
      mainz[n][m]=pea;
      //    }
  }
  aaa.close();

  TCanvas *c2 = new TCanvas("c2","c2",896,900);
  TH2F *hitMap = new TH2F("hitMap","Light yield ratio 200fF/600fF;y [mm];x [mm]",12,-180,180,12,-180,180);
  hitMap->SetMinimum(0.88);
  hitMap->SetMaximum(1.12);

  char histoName[256];
  char imgName[256];
  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  FILE *output = fopen("../txt/ratio_200_600_138_138.txt","w");
  for (int i=0;i<4;i++) {
    for (int j=0;j<36;j++) {
      if (i==0 || i==1) {
  	hitMap->Fill(165*(1-i)-15*i-30*(j%6),-15-30*int(j/6),beam[i][j]/mainz[i][j]);
  	fprintf(output,"%d %d\t%f\n", i, j, beam[i][j]/mainz[i][j]);
      }
      else {
  	hitMap->Fill(-165*(i-2)+15*(3-i)+30*(j%6),165-30*int(j/6),beam[i][j]/mainz[i][j]);
  	fprintf(output,"%d %d\t%f\n", i, j, beam[i][j]/mainz[i][j]);
      }
    }
  }
  hitMap->Draw("colztext");
  fclose(output);
}
