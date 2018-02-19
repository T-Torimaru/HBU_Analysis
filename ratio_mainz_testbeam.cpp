using namespace TMath;

void ratio_mainz_testbeam() {

  Int_t n,m;

  Float_t pea;
  Float_t beam[4][36];
  Float_t mainz[4][36];

  ifstream ifs("../txt/testBeamresults.txt");
  while(ifs>>n>>m>>pea){
    beam[n][m]=pea;
  }
  ifs.close();

  ifstream aaa("../txt/mainz_hitmap.txt");
  while(aaa>>n>>m>>pea){
    mainz[n][m]=pea;
  }
  aaa.close();

  TCanvas *c2 = new TCanvas("c2","c2",896,900);
  TH2F *hitMap = new TH2F("hitMap","Light yield ratio TestBeam/Mainz;y [mm];x [mm]",12,-180,180,12,-180,180);
  hitMap->SetMinimum(0.88);
  hitMap->SetMaximum(1.12);

  char histoName[256];
  char imgName[256];
  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  FILE *output = fopen("../txt/ratio_mainz_testbeam.txt","w");
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
