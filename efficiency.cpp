void efficiency(){
  TString matchfile;
  matchfile.Form("../rootfile/cosmic/18Feb20_6layer/cosmi_0220.root");
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
  Double_t entry_center[6][14][14]={0};
  Double_t recoX[2][8];
  Double_t recoY[2][8];
  Double_t tReco[6];
  Double_t HBUPos[4][3];
  Double_t deciX,deciY;
  Int_t nPass;
  Int_t nEntries=tmatch->GetEntries();
  Int_t count[6][14][14]={0};
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

  TCanvas *canv1 = new TCanvas("c1","c1",696,500);
  TH1F *efficiency =new TH1F("Efficiency","Efficiency",30,0.8,1.);

  gStyle->SetOptStat(1111);

  for (int i=0;i<nEntries;i++){
    if (i==7294 || i==17572 || i==33381 || i==42705 || i==44465 || i==48873 || i==48979 || i==50060 || i==68587 || i==108446 || i==145457 || i==198356 || i==217244 || i==231350 || i==291337) continue;
    tmatch->GetEntry(i);
    if (Hnhit>18 || nPass<5) continue;
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

  	  nTmp = HhitI[j];
  	  HhitI[j] = HhitI[k];
  	  HhitI[k] = nTmp;
  	  nTmp = HhitJ[j];
  	  HhitJ[j] = HhitJ[k];
  	  HhitJ[k] = nTmp;
  	  nTmp = HhitK[j];
  	  HhitK[j] = HhitK[k];
  	  HhitK[k] = nTmp;

  	  tmp = HADC[j];
  	  HADC[j] = HADC[k];
  	  HADC[k] = tmp;
  	}
      }
    }
      
    if (nPass==5&&HhitK[0]!=1) continue;

    for (int j=0;j<Hnhit;j++){
      if (HADC[j]<650) continue;
      ntp=HhitK[j]-1;
      HBU_z[j]=constHBU_z[ntp];
      HBU_x[j]=(tReco[3]-tReco[0])*(HBU_z[j]-tReco[2])/(tReco[5]-tReco[2])+tReco[0];
      HBU_y[j]=(tReco[4]-tReco[1])*(HBU_z[j]-tReco[2])/(tReco[5]-tReco[2])+tReco[1];
      xaxisI=HhitI[j]-6;
      yaxisJ=HhitJ[j]-6;
      deciX=HBU_x[j]+210;
      deciY=HBU_y[j]+210;
      tileI=deciX/30;
      tileJ=deciY/30;

      // if (ntp==1) {
      // 	HBU_x[j]=HBU_x[j]+4.;
      // 	HBU_y[j]=HBU_y[j]+2.;
      // }
      //      if (ntp==2) {
	//      	HBU_x[j]=HBU_x[j]+0.35;
      // 	HBU_y[j]=HBU_y[j]+5.;
      // }
      if (HhitK[j+1]==HhitK[j]) continue;

      if ( 
	  HBU_x[j]>-205.+30.*tileI
	  && HBU_x[j]<-215.+30.*(tileI+1)
	  && HBU_y[j]>-205.+30.*tileJ
	  && HBU_y[j]<-215.+30.*(tileJ+1)
	  && tileI < 14 && tileJ < 14
	  && tileI >= 0 && tileJ >= 0
	   )
	{
	  entry_center[ntp][tileI][tileJ]++;
	  if ( 
	      tileI==xaxisI && tileJ==yaxisJ
	      && xaxisI < 14 && yaxisJ < 14
	      && xaxisI >= 0 && yaxisJ >= 0
	       )
	    {
	      count[ntp][xaxisI][yaxisJ]++;
	    }
	}
    }
    // cout<<"entry_center : "<<entry_center[5]<<endl;	 
    // cout<<"count[0][1][2] : "<<count[0][1][2]<<endl;
    // cout<<"count[0][2][1] : "<<count[0][2][1]<<endl;
    // cout<<"count[0][1][0] : "<<count[0][1][0]<<endl;
    // cout<<"count[0][0][1] : "<<count[0][0][1]<<endl;
    // cout<<"count[0][1][1] : "<<count[0][1][1]<<endl;
    // cout<<"entry_center[0]: "<<entry_center[0]<<endl;
    // cout<<""<<endl;
  }
  
  //  FILE *output = fopen("../txt/shift/shift_6layer_0220_2_2.txt","w");
  for (int layer=0;layer<6;layer++) {
    for (int xaxis=0;xaxis<14;xaxis++) {
      for (int yaxis=0;yaxis<14;yaxis++) {
  	efficiency->Fill(count[layer][xaxis][yaxis]/entry_center[layer][xaxis][yaxis]);
	//  	fprintf(output,"%d %d %d\t%f\n",layer,xaxis,yaxis,count[layer][xaxis][yaxis]/entry_center[layer]);
      }
    }
  }
  canv1->cd();
  efficiency->Draw();
  // sprintf(canvasName,"../png/efficiency/shift/cosmic0220_%d_%d_layer5.png",rangeI+1,rangeJ+1);
  // canv1->SaveAs(canvasName);
  //  fclose(output);
  cout<<"count[0][5][5]: "<<count[0][5][5]<<endl;
  cout<<"entry_center[0][5][5]: "<<entry_center[0][5][5]<<endl;
}
