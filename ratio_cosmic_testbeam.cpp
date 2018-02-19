using namespace TMath;

void ratio_cosmic_testbeam() {

  Int_t n,m,k,l;

  Float_t pea;
  Float_t beam[4][36];
  Float_t cosmic[12][12];
  Float_t cosm[4][36];

  ifstream ifs("../txt/testBeamresults.txt");
  while(ifs>>n>>m>>pea){
    beam[n][m]=pea;
  }
  ifs.close();

  ifstream aaa("../txt/DESY_true.txt");
  while(aaa>>n>>m>>pea){
    cosmic[n-7][m-7]=pea;
  }
  aaa.close();
  
  for (int i=0;i<4;i++) {
    for (int j=0;j<36;j++){
      if (i==0 || i==1){
        cosm[i][j] = cosmic[5-int(j/6)][11*(1-i)+5*i-j%6];
      }
      else {
        cosm[i][j] = cosmic[11-int(j/6)][6*(3-i)+j%6];
      }
    }
  }
  

  TCanvas *c2 = new TCanvas("c2","c2",896,900);
  TH2F *hitMap = new TH2F("hitMap","Light yield ratio TestBeam/Cosmic;y [mm];x [mm]",12,-180,180,12,-180,180);
  hitMap->SetMaximum(1.12);
  //  hitMap->SetMinimum(0.88);
  hitMap->SetMinimum(0.78);
  
  char histoName[256];
  char imgName[256];
  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  FILE *output = fopen("../txt/ratio_cosmic_testbeam.txt","w");
  for (int i=0;i<4;i++) {
    for (int j=0;j<36;j++) {
      if (i==0 || i==1) {
  	hitMap->Fill(165*(1-i)-15*i-30*(j%6),-15-30*int(j/6),beam[i][j]/cosm[i][j]);
  	fprintf(output,"%d %d\t%8.1f\n", i, j, beam[i][j]/cosm[i][j]);
      }
      else {
  	hitMap->Fill(-165*(i-2)+15*(3-i)+30*(j%6),165-30*int(j/6),beam[i][j]/cosm[i][j]);
  	fprintf(output,"%d %d\t%8.1f\n", i, j, beam[i][j]/cosm[i][j]);
      }
    }
  }
  hitMap->Draw("colztext");
  fclose(output);
}
