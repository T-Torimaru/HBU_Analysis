using namespace TMath;
void event_display2(int evenum = 3370){
  TString matchfile;
  matchfile.Form("../rootfile/cosmic/18Feb20_6layer/cosmi_0220_2262.root");
  //  matchfile.Form("../rootfile/everything/testWE.root");
  TFile *fmatch = new TFile(matchfile.Data(), "read");
  TTree *tmatch = (TTree*)fmatch->Get("bigtree");
  Int_t HhitI[256],HhitJ[256],HhitK[256],Hnhit;
  Float_t HhitPos[256][3];
  Float_t  HADC[256];
  Double_t nph[2][64];
  Int_t    nX[2];
  Int_t    nY[2];
  Double_t recoX[2][8];
  Double_t recoY[2][8];
  Double_t trueReco[6];
  Int_t nMax;
  // Double_t HBU_x,HBU_y;
  // Double_t HBU_z=217.8;
  //  Double_t SquDis;

  tmatch->SetBranchAddress("ahc_nHits", &Hnhit);
  tmatch->SetBranchAddress("hod_maxPass", &nMax);
  tmatch->SetBranchAddress("ahc_hitI", HhitI);
  tmatch->SetBranchAddress("ahc_hitJ", HhitJ);
  tmatch->SetBranchAddress("ahc_hitK", HhitK);
  tmatch->SetBranchAddress("ahc_hitPos", HhitPos);
  tmatch->SetBranchAddress("ahc_hitEnergy", HADC);
  //  tmatch->SetBranchAddress("squareDis", SquDis);
  TString ephname,nxname,nyname,recxname,recyname,truexname,trueyname,truezname;
  for (int m = 0; m < 2; m++) {
    ephname.Form("hod%d_nph",m+1);
    nxname.Form("hod%d_nRecoX",m+1);
    nyname.Form("hod%d_nRecoY",m+1);
    recxname.Form("hod%d_recoX",m+1);
    recyname.Form("hod%d_recoY",m+1);
    truexname.Form("hod%d_trueRecoX",m+1);
    trueyname.Form("hod%d_trueRecoY",m+1);
    truezname.Form("hod%d_trueRecoZ",m+1);
    tmatch->SetBranchAddress(ephname,nph[m]);
    tmatch->SetBranchAddress(nxname,&nX[m]);
    tmatch->SetBranchAddress(nyname,&nY[m]);
    tmatch->SetBranchAddress(recxname,recoX[m]);
    tmatch->SetBranchAddress(recyname,recoY[m]);
    tmatch->SetBranchAddress(truexname,&trueReco[3*m+0]);
    tmatch->SetBranchAddress(trueyname,&trueReco[3*m+1]);
    tmatch->SetBranchAddress(truezname,&trueReco[3*m+2]);
  }

  for (int i=evenum;i<tmatch->GetEntries();i++) {
    tmatch->GetEntry(i);
    if (Hnhit<=2) continue;
    int flag[2]={0,0};
    for (int n=0;n<Hnhit;n++) {
      if (HhitK[n]==2) flag[0]=1;
      else flag[1]=1;
    }
    if (flag[0]==1 && flag[1]==1 ) {
      evenum=i;
      break;
    }
    else continue;
  }
  
  TCanvas *canvas1 = new TCanvas("canvas1","scintillator plate",500,500);
  TCanvas *canvas2 = new TCanvas("canvas2","scintillator plate",500,500);
  TCanvas *canvas3 = new TCanvas("canvas3","combined",900,600);
  //  TCanvas *canvas4 = new TCanvas("canvas4","Square Distance",696,500);

  /*Get event*/
  tmatch->GetEntry(evenum);
  Double_t ly[2][2][16];
  for(int ch=0;ch<16;ch++) {
    ly[0][0][ch]=nph[0][32+(ch+12)%16]+nph[0][63-(ch+12)%16];
    ly[0][1][ch]=nph[0][15-(ch+12)%16]+nph[0][16+(ch+12)%16];
    ly[1][0][ch]=nph[1][ch]+nph[1][16+(ch+12)%16];
    ly[1][1][ch]=nph[1][47-ch]+nph[1][48+ch];
  }

  cout<<"evenum: "<<evenum<<endl;
  cout<<"  hod1: "<<nX[0]<<"x"<<nY[0]<<endl;
  cout<<"  hod2: "<<nX[1]<<"x"<<nY[1]<<endl;
  cout<<"  nhit: "<<Hnhit<<endl;
  cout<<"nMax: "<<nMax<<endl;
  //  cout<<"  hitK: "<<HhitK[3]<<endl;

  TH2F *hcr[2];
  //  TH1F *squdist;
  TGraphErrors *recoG[2];
  TGraph2D *Edisplay = new TGraph2D();
  TGraph2D *Hdisplay = new TGraph2D();
  TPolyLine3D *tile[256];
  int gcount=0;
  double height[2]={0,43.3*13};

  hcr[0] = new TH2F("hcr1","CR1 ly display;x[mm];y[mm]",84,-210,210,84,-210,210);
  hcr[1] = new TH2F("hcr2","CR2 ly display;x[mm];y[mm]",84,-210,210,84,-210,210);
  //  squdist = new TH1F("squdist","Square Distance;x^2[mm];Entry",200,0,300);
  for (int m=0;m<2;m++) {
    for (int chx=0;chx<84;chx++) {
      for (int chy=0;chy<84;chy++) {
        hcr[m]->Fill((chx+0.5-42)*5,(chy+0.5-42)*5,ly[m][0][chx%16]*ly[m][1][chy%16]);
      }
    }
  }
  // Double_t fakeHBU_x,fakeHBU_y;
  // Double_t distance;
  // Double_t sigma=9. + (30/Sqrt(12.))**2;
  //   cout<<"  sigma: "<<sigma<<endl;
  // if (HhitK[3]==2){
  //   HBU_x=(trueReco[3]-trueReco[0])*(HBU_z-trueReco[2])/(trueReco[5]-trueReco[2])+trueReco[0];
  //   HBU_y=(trueReco[4]-trueReco[1])*(HBU_z-trueReco[2])/(trueReco[5]-trueReco[2])+trueReco[1];
  //   for (int i=0;i<nX[0];i++){
  //     for (int j=0;j<nX[1];j++){      
  // 	for (int k=0;k<nY[0];k++){
  // 	  for (int l=0;l<nY[1];l++){      
  // 	    fakeHBU_x=(recoX[1][i]*HBU_z+(trueReco[5]-HBU_z)*recoX[0][j])/(trueReco[5]-trueReco[2]);
  // 	    fakeHBU_y=(recoY[1][l]*HBU_z+(trueReco[5]-HBU_z)*recoY[0][k])/(trueReco[5]-trueReco[2]);
  // 	    distance=Sqrt((HhitPos[Hnhit][0]-fakeHBU_x)*(HhitPos[Hnhit][0]-fakeHBU_x)+(HhitPos[Hnhit][1]-fakeHBU_y)*(HhitPos[Hnhit][1]-fakeHBU_y));
  // 	    squdist->Fill(distance**2/sigma);
  // 	    cout<<"  distance: "<<distance<<endl;
  // 	  }
  // 	}
  //     }
  //   }
  // }
  Edisplay->SetMarkerStyle(20);
  Hdisplay->SetMarkerStyle(21);
  Hdisplay->SetMarkerColor(2);
  Edisplay->SetTitle("3d display;x[mm];y[mm];z[mm]");
  for (int m=0;m<2;m++) {
    recoG[m] = new TGraphErrors();
    recoG[m]->SetMarkerStyle(20);
    recoG[m]->SetMarkerColor(0);
    for (int x=0;x<nX[m];x++) {
      for (int y=0;y<nY[m];y++) {
        recoG[m]->SetPoint(x*nY[m]+y,recoX[m][x],recoY[m][y]);
        Edisplay->SetPoint(gcount,recoX[m][x],recoY[m][y],height[m]);
        gcount++;
      }
    }
  }
  //trueReco[2] = height[0];
  //trueReco[5] = height[1];
  TPolyLine3D *line = new TPolyLine3D(2,trueReco);
  line->SetLineColor(2);

  cout<<trueReco[0]<<" "<<trueReco[1]<<" "<<trueReco[2]<<" "<<trueReco[3]
  <<" "<<trueReco[4]<<" "<<trueReco[5]<<endl;

  for (int n=0;n<Hnhit;n++) {
/*
    if(HhitK[n]==8) {
      HhitPos[n][0]=(HhitI[n]-12.5)*30.;
      HhitPos[n][1]=(HhitJ[n]-18.5)*30.;
      HhitPos[n][2]=43.3*3;
    } else if(HhitK[n]==9) {
      HhitPos[n][0]=(HhitI[n]-12.5)*30.;
      HhitPos[n][1]=(HhitJ[n]-6.5)*30.;
      HhitPos[n][2]=43.3*3;
    } else {
      HhitPos[n][0]=(HhitI[n]-12.5)*30.;
      HhitPos[n][1]=(HhitJ[n]-6.5)*30.+75.;
      HhitPos[n][2]=43.3*12;
    }

    if(HhitK[n]==1) {
      HhitPos[n][0]+=0;
      HhitPos[n][1]+=255;
      HhitPos[n][2]+=43.3;
    } else if(HhitK[n]==2) {
      HhitPos[n][0]+=-45;
      HhitPos[n][1]+=0;
      HhitPos[n][2]+=43.3*9;
    } else if(HhitK[n]==3) {
      HhitPos[n][0]+=-30;
      HhitPos[n][1]+=138;
      HhitPos[n][2]+=43.3*3;
    }
*/
    Hdisplay->SetPoint(n,HhitPos[n][0],HhitPos[n][1],HhitPos[n][2]);
    Double_t tilecorner[15] = {
HhitPos[n][0]-15,HhitPos[n][1]-15,HhitPos[n][2],
HhitPos[n][0]-15,HhitPos[n][1]+15,HhitPos[n][2],
HhitPos[n][0]+15,HhitPos[n][1]+15,HhitPos[n][2],
HhitPos[n][0]+15,HhitPos[n][1]-15,HhitPos[n][2],
HhitPos[n][0]-15,HhitPos[n][1]-15,HhitPos[n][2]};
    tile[n] = new TPolyLine3D(5,tilecorner);
    cout<<n<<" "<<HhitPos[n][0]<<" "<<HhitPos[n][1]<<" "<<HhitPos[n][2]<<endl;
  }

  Edisplay->GetXaxis()->SetRangeUser(-210,210);
  Edisplay->GetYaxis()->SetRangeUser(-210,210);
  Edisplay->GetZaxis()->SetRangeUser(0,43.3*13);
  canvas1->cd();
  hcr[0]->Draw("colz");
  recoG[0]->Draw("Psame");
  canvas2->cd();
  hcr[1]->Draw("colz");
  recoG[1]->Draw("Psame");
  canvas3->cd();
  Edisplay->Draw("P");
  Hdisplay->Draw("P same");
  for (int n=0;n<Hnhit;n++) tile[n]->Draw("same");
  line->Draw("same");
  // canvas4->cd();
  // squdist->Draw();
}
