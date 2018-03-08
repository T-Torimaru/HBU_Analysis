void shift(){
  TString matchfile;
  matchfile.Form("../rootfile/cosmic/18Feb20_6layer/cosmi_0220_2262.root");
  TFile *fmatch = new TFile(matchfile.Data(), "read");
  TTree *tmatch = (TTree*)fmatch->Get("bigtree");
  char canvasName[256];
  char histoName[256];
  Int_t HhitI[256],HhitJ[256],HhitK[256],Hnhit;
  Float_t HhitPos[256][3];
  Float_t  HADC[256];
  Double_t nph[2][64];
  Int_t    nX[2];
  Int_t    nY[2];
  Int_t tileI,tileJ,xaxisI,yaxisJ;
  Double_t entry_center[6]={0};
  Double_t recoX[2][8];
  Double_t recoY[2][8];
  Double_t tReco[6];
  Double_t HBUPos[4][3];
  Double_t effic1[4][14][14];
  Double_t effic2[4][14][14];
  Int_t nPass;
  Int_t nEntries=tmatch->GetEntries();
  Int_t rangeI = 6;
  Int_t rangeJ = 6;
  Int_t count[6][3][3]={0};
  Double_t HBU_x[64]={0},HBU_y[64]={0},HBU_z[64]={0};
  Double_t constHBU_z[6]={131.15,174.45,217.75,261.05,304.35,347.65};
    
  tmatch->SetBranchAddress("ahc_nHits", &Hnhit);
  tmatch->SetBranchAddress("ahc_hitI", HhitI);
  tmatch->SetBranchAddress("ahc_hitJ", HhitJ);
  tmatch->SetBranchAddress("ahc_hitK", HhitK);
  tmatch->SetBranchAddress("ahc_hitPos", HhitPos);
  tmatch->SetBranchAddress("ahc_hitEnergy", HADC);
  tmatch->SetBranchAddress("hod_maxPass",&nPass);
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
    tmatch->SetBranchAddress(truexname,&tReco[3*m+0]);
    tmatch->SetBranchAddress(trueyname,&tReco[3*m+1]);
    tmatch->SetBranchAddress(truezname,&tReco[3*m+2]);
  }

  TH2F *shift[6];
  TCanvas *canv1 = new TCanvas("c1","c1",500,500);

  for (int layer=0;layer<6;layer++) {
    sprintf(histoName,"Shift_Layer%d",layer+1);
    shift[layer] = new TH2F(histoName,"Shift;y [mm];x[mm]",3,0,90,3,0,90);
    shift[layer]->SetMaximum(1.);
    shift[layer]->SetMinimum(0.);
  }

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  for (int i=0;i<nEntries;i++){
    if (i==7294) continue;
    tmatch->GetEntry(i);
    if (Hnhit>64) continue;
    Double_t tmp;
    Int_t nTmp,ntp;
    for (int j=0;j<Hnhit;j++){
      for (int k=j+1;k<Hnhit;k++){
  	if (HhitPos[j][2] > HhitPos[k][2]){
  	  tmp = HhitPos[j][2];
  	  HhitPos[j][2] = HhitPos[k][2];
  	  HhitPos[k][2] = tmp;
  	  tmp = HhitPos[j][0];
  	  HhitPos[j][0] = HhitPos[k][0];
  	  HhitPos[k][0] = tmp;
  	  tmp = HhitPos[j][1];
  	  HhitPos[j][1] = HhitPos[k][1];
  	  HhitPos[k][1] = tmp;
  	}
  	if (HhitK[j] > HhitK[k]){
  	  nTmp = HhitI[j];
  	  HhitI[j] = HhitI[k];
  	  HhitI[k] = nTmp;
  	  nTmp = HhitJ[j];
  	  HhitJ[j] = HhitJ[k];
  	  HhitJ[k] = nTmp;
  	  nTmp = HhitK[j];
  	  HhitK[j] = HhitK[k];
  	  HhitK[k] = nTmp;
  	}
      }
    }
    for (int j=0;j<Hnhit;j++){
      ntp=HhitK[j]-1;
      HBU_z[j]=constHBU_z[ntp];
      HBU_x[j]=(tReco[3]-tReco[0])*(HBU_z[j]-tReco[2])/(tReco[5]-tReco[2])+tReco[0];
      HBU_y[j]=(tReco[4]-tReco[1])*(HBU_z[j]-tReco[2])/(tReco[5]-tReco[2])+tReco[1];
      xaxisI=HhitI[j]-rangeI-6;
      yaxisJ=HhitJ[j]-rangeJ-6;

      if (   HBU_x[j]>-210.+30.*(rangeI+1)
      	     && HBU_x[j]<-210.+30.*(rangeI+2)
      	     && HBU_y[j]>-210.+30.*(rangeJ+1)
      	     && HBU_y[j]<-210.+30.*(rangeJ+2)
      	     && xaxisI < 3 && yaxisJ < 3
      	     && xaxisI >= 0 && yaxisJ >= 0
      	     ) {
	if (HhitK[j+1]!=HhitK[j]) {
	  entry_center[ntp]++;
	}
      	count[ntp][xaxisI][yaxisJ]++;
      	}
    }
      
    // if (i==20) {
    //   cout<<"Hnhit: "<<Hnhit<<endl;      
    //   for (int j=0;j<Hnhit;j++){
    // 	cout<<"HhitPosX: "<<HhitPos[j][0]<<endl;
    // 	cout<<"HhitPosY: "<<HhitPos[j][1]<<endl;
    // 	cout<<"HhitPosZ: "<<HhitPos[j][2]<<endl;
    // 	cout<<"HBUPosX: "<<HBU_x[j]<<endl;
    // 	cout<<"HBUPosY: "<<HBU_y[j]<<endl;
    // 	cout<<"HBUPosZ: "<<HBU_z[j]<<endl;
    // 	cout<<"HhitI: "<<HhitI[j]-6<<endl;
    // 	cout<<"HhitJ: "<<HhitJ[j]-6<<endl;
    // 	cout<<"HhitK: "<<HhitK[j]<<endl;
    //   }
    // }
    
  }

  FILE *output = fopen("../txt/shift/shift_6layer_0220_6_6.txt","w");
  for (int layer=0;layer<6;layer++) {
    for (int xaxis=0;xaxis<3;xaxis++) {
      for (int yaxis=0;yaxis<3;yaxis++) {
  	shift[layer]->Fill(xaxis*30+15,yaxis*30+15,count[layer][xaxis][yaxis]/entry_center[layer]);
  	fprintf(output,"%d %d %d\t%f\n",layer,xaxis,yaxis,count[layer][xaxis][yaxis]/entry_center[layer]);
      }
    }
  }
  canv1->cd();
  shift[2]->Draw("colztext");
  sprintf(canvasName,"../png/efficiency/shift/cosmic0220_%d_%d_layer3.png",rangeI+1,rangeJ+1);
  canv1->SaveAs(canvasName);
  fclose(output);

}
