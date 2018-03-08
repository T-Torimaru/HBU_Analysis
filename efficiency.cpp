void efficiency(){
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

  // TH1F *efficiencyH = new TH1F("histo","Tile Response Efficiency",20,0,1);
  // TCanvas *c1 = new TCanvas("c1","c1",696,500);

  TH2F *shift[6];
  TCanvas *canv1 = new TCanvas("c1","c1",500,500);
  // TCanvas *canv2 = new TCanvas("c1","c2",500,500);
  // TCanvas *canv3 = new TCanvas("c1","c3",500,500);
  // TCanvas *canv4 = new TCanvas("c1","c4",500,500);
  // TCanvas *canv5 = new TCanvas("c1","c5",500,500);
  // TCanvas *canv6 = new TCanvas("c1","c6",500,500);

  for (int layer=0;layer<6;layer++) {
    sprintf(histoName,"Shift_Layer%d",layer+1);
    shift[layer] = new TH2F(histoName,"Shift;y [mm];x[mm]",3,0,90,3,0,90);
    shift[layer]->SetMaximum(1.);
    shift[layer]->SetMinimum(0.);
  }

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  // Int_t modX,modY;
  // Double_t deciX,deciY;
  // Double_t valid1=2;
  // Double_t valid2=28;
  // Double_t val1,val2;
  // Double_t tileEdge[2][2];
  // bool eff[4];
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

  // cout<<"count[0][0][0]: "<<count[0][1][1]<<endl;
  // cout<<"Entry: "<<nEntries<<endl;
  // cout<<"entry_center[0]: "<<entry_center[0]<<endl;

  FILE *output = fopen("../txt//shift/shift_6layer_0220_6_6.txt","w");
  for (int layer=0;layer<6;layer++) {
    for (int xaxis=0;xaxis<3;xaxis++) {
      for (int yaxis=0;yaxis<3;yaxis++) {
  	shift[layer]->Fill(xaxis*30+15,yaxis*30+15,count[layer][xaxis][yaxis]/entry_center[layer]);
  	fprintf(output,"%d %d %d\t%f\n",layer,xaxis,yaxis,count[layer][xaxis][yaxis]/entry_center[layer]);
      }
    }
  }
  canv1->cd();
  shift[0]->Draw("colztext");
  sprintf(canvasName,"../png/efficiency/shift/cosmic0220_%d_%d.png",rangeI+1,rangeJ+1);
  canv1->SaveAs(canvasName);
  fclose(output);
  cout<<"devide: "<<count[0][1][1]/entry_center[0]<<endl;
  cout<<"count[0][2][2]: "<<count[0][2][2]<<endl;
  cout<<"entry_center[0]: "<<entry_center[0]<<endl;
  // canv2->cd();
  // shift[1]->Draw();
  // sprintf(canvasName,"../png/efficiency/shift/cosmic0220_%d_%d.png",rangeI+1,rangeJ+1);
  // canv2->SaveAs(canvasName);
  // canv3->cd();
  // shift[2]->Draw();
  // sprintf(canvasName,"../png/efficiency/shift/cosmic0220_%d_%d.png",rangeI+1,rangeJ+1);
  // canv3->SaveAs(canvasName);
  // canv4->cd();
  // shift[3]->Draw();
  // sprintf(canvasName,"../png/efficiency/shift/cosmic0220_%d_%d.png",rangeI+1,rangeJ+1);
  // canv4->SaveAs(canvasName);
  // canv5->cd();
  // shift[4]->Draw();
  // sprintf(canvasName,"../png/efficiency/shift/cosmic0220_%d_%d.png",rangeI+1,rangeJ+1);
  // canv5->SaveAs(canvasName);
  // canv6->cd();
  // shift[5]->Draw();
  // sprintf(canvasName,"../png/efficiency/shift/cosmic0220_%d_%d.png",rangeI+1,rangeJ+1);
  // canv6->SaveAs(canvasName);

  // cout<<"count : "<<count[0][1][1]<<endl;
  // cout<<"count : "<<count[1][0][1]<<endl;
    
    // Double_t count1=0;
  // Double_t count2=0;

  // for (int j=0;j<4;j++){
  //   HBU_x[j]=(tReco[3]-tReco[0])*(HBU_z[j]-tReco[2])/(tReco[5]-tReco[2])+tReco[0];
  //   HBU_y[j]=(tReco[4]-tReco[1])*(HBU_z[j]-tReco[2])/(tReco[5]-tReco[2])+tReco[1];
  //   if (HBU_x[j]<0) {
  //     modX=-HBU_x[j];
  //     tileI=6-modX/30;
  //     deciX=-HBU_x[j]-modX;
  //   }
  //   else {
  //     modX=HBU_x[j];
  //     deciX=HBU_x[j]-modX;
  //     tileI=modX/30+7;
  //   }
  //   if (HBU_y[j]<0) {
  //     modY=-HBU_y[j];
  //     tileJ=6-modY/30;
  //     deciY=-HBU_y[j]-modY;
  //   }
  //   else {
  //     modY=HBU_y[j];
  //     deciY=HBU_y[j]-modY;
  //     tileJ=modY/30+7;
  //   }
  //   modX=modX%30;
  //   modY=modY%30;
  //   val1=modX+deciX;
  //   val2=modY+deciY;
  //   if (val1<valid1 || val1>valid2 || val2<valid1 || val2>valid2 || tileI<0 || tileI>13 || tileJ<0 || tileJ>13) eff[j]=false;
  //   else {
  //     eff[j]=true;
  //     effic1[j][tileI][tileJ]=effic1[j][tileI][tileJ]+1;
  //   }
    // if (i==16) {
    //   //      cout<<effic1[j][tileI][tileJ]<<endl; 
    //   //      cout<<HBU_x[j]<<" "<<HBU_y[j]<<endl;    
    //   cout<<j<<" "<<tileI<<" "<<tileJ<<endl; 
    //    }
  // if (i==16){
  //   cout<<effic1[0][3][7]<<endl;    
  //   cout<<effic1[1][3][7]<<endl;    
  //   cout<<effic1[3][8][13]<<endl;    
  // }
  // for (int k=0;k<Hnhit;k++){
  //   if (eff[HhitK[k]-1]){
  //   tileEdge[0][0]=HhitPos[k][0]+15;
  //   tileEdge[0][1]=HhitPos[k][0]-15;
  //   tileEdge[1][0]=HhitPos[k][1]+15;
  //   tileEdge[1][1]=HhitPos[k][1]-15;
       // if (i==16){
       // 	 cout<<"HhitK : "<<HhitK[k]-1<<" "<<"HhitI : "<<HhitI[k]-6<<" "<<"HhitJ : "<<HhitJ[k]-6<<endl;
       // //       effic2[HhitK-1][HhitI-7][HhitJ-7];
       // }
    // cout<<HhitK[k]<<endl;
    // cout<<HBU_x[HhitK[k]-1]<<" "<<HBU_y[HhitK[k]-1]<<endl;
    // cout<<tileEdge[0][1]<<" "<<tileEdge[0][0]<<" "<<tileEdge[1][1]<<" "<<tileEdge[1][0]<<endl;    
  //    if (tileEdge[0][1]<HBU_x[HhitK[k]-1] && tileEdge[0][0]>HBU_x[HhitK[k]-1] && tileEdge[1][0]>HBU_y[HhitK[k]-1] && tileEdge[1][1]<HBU_y[HhitK[k]-1]) {
  //      effic2[HhitK[k]-1][HhitI[k]-6][HhitJ[k]-6]++;
  //    }
  //   }
  // }

  // for (int i=0;i<4;i++){
  //   for (int j=0;j<14;j++) {
  //     for (int k=0;k<14;k++) {
  // 	if (effic1[i][j][k]!=0){
  // 	  efficiencyH->Fill(effic2[i][j][k]/effic1[i][j][k]);
  // 	}
  //     }
  //   }
  // }
  // c1->cd();
  // efficiencyH->Draw();
  // sprintf(canvasName,"../png/efficiency_2017Dec_test.png");
  // c1->SaveAs(canvasName);
  // cout<<effic1[3][5][8]<<endl;
  // cout<<effic1[2][5][8]<<endl;
  // cout<<effic2[3][5][8]<<endl;    
  // cout<<effic2[2][5][8]<<endl;    
}
