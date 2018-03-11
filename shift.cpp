void shift(){
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
  Double_t entry_center[6]={0};
  Double_t recoX[2][8];
  Double_t recoY[2][8];
  Double_t tReco[6];
  Double_t HBUPos[4][3];
  Double_t effic1[4][14][14];
  Double_t effic2[4][14][14];
  Int_t nPass;
  Int_t nEntries=tmatch->GetEntries();
  Int_t rangeI = 2;
  Int_t rangeJ = 9;
  Int_t count[6][3][3]={0};
  Double_t HBU_x[64]={0},HBU_y[64]={0},HBU_z[64]={0},fHBU_x,fHBU_y,distance;
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
  TCanvas *canv2 = new TCanvas("c2","c2",696,500);
  TH1F *xi =new TH1F("xi^2","xi^2",100,0,100);

  for (int layer=0;layer<6;layer++) {
    sprintf(histoName,"Shift_Layer%d",layer+1);
    shift[layer] = new TH2F(histoName,"Shift;y [mm];x[mm]",3,0,90,3,0,90);
    shift[layer]->SetMaximum(1.);
    shift[layer]->SetMinimum(0.);
  }

  // gStyle->SetOptStat(0);
  // gStyle->SetPalette(1);
  //  gStyle->SetOptFit();
  gStyle->SetOptStat(1111);

  for (int i=0;i<50;i++){
    if (i==7294 || i==17572 || i==33381 || i==42705 || i==44465 || i==48873 || i==48979 || i==50060 || i==68587 || i==108446 || i==145457 || i==198356 || i==217244 || i==231350 || i==291337) continue;
    tmatch->GetEntry(i);
    if (Hnhit>18 || nPass<5) continue;
    Double_t tmp;
    Int_t nTmp,ntp;
    // for (int j=0;j<Hnhit;j++){
    //   cout<<"HhitK : "<<HhitK[j]<<endl;     
    // }
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
      //      cout<<"HhitK : "<<HhitK[j]<<endl;     
    }
    //      cout<<" "<<endl; 
      
    if (nPass==5&&HhitK[0]!=1) continue;

    for (int j=0;j<Hnhit;j++){
      if (HADC[j]<650) continue;
      ntp=HhitK[j]-1;
      HBU_z[j]=constHBU_z[ntp];
      HBU_x[j]=(tReco[3]-tReco[0])*(HBU_z[j]-tReco[2])/(tReco[5]-tReco[2])+tReco[0];
      HBU_y[j]=(tReco[4]-tReco[1])*(HBU_z[j]-tReco[2])/(tReco[5]-tReco[2])+tReco[1];
      xaxisI=HhitI[j]-rangeI-6;
      yaxisJ=HhitJ[j]-rangeJ-6;

      // if (ntp==1) {
      // 	HBU_x[j]=HBU_x[j]+4.;
      // 	HBU_y[j]=HBU_y[j]+2.;
      // }
      //      if (ntp==2) {
	//      	HBU_x[j]=HBU_x[j]+0.35;
      // 	HBU_y[j]=HBU_y[j]+5.;
      // }
      if (   HBU_x[j]>-205.+30.*(rangeI+1)
      	     && HBU_x[j]<-215.+30.*(rangeI+2)
      	     && HBU_y[j]>-205.+30.*(rangeJ+1)
      	     && HBU_y[j]<-215.+30.*(rangeJ+2)
      	     && xaxisI < 3 && yaxisJ < 3
      	     && xaxisI >= 0 && yaxisJ >= 0
      	     ) {
	if (HhitK[j+1]!=HhitK[j]) {
	  entry_center[ntp]++;
      }
      	count[ntp][xaxisI][yaxisJ]++;
      }
      // cout<<"HBU_x : "<<HBU_x[j]<<endl;
      // cout<<"HBU_y : "<<HBU_y[j]<<endl;
      // cout<<"HBU_z : "<<HBU_z[j]<<endl;
      // cout<<"HhitI : "<<HhitI[j]<<endl;
      // cout<<"HhitJ : "<<HhitJ[j]<<endl;
      // cout<<"HhitK : "<<HhitK[j]<<endl;
      // cout<<" "<<endl;
     }
    // cout<<"entry_center : "<<entry_center[5]<<endl;	 
    // cout<<"count[0][1][2] : "<<count[0][1][2]<<endl;
    // cout<<"count[0][2][1] : "<<count[0][2][1]<<endl;
    // cout<<"count[0][1][0] : "<<count[0][1][0]<<endl;
    // cout<<"count[0][0][1] : "<<count[0][0][1]<<endl;
    // cout<<"count[0][1][1] : "<<count[0][1][1]<<endl;
    // cout<<"entry_center[0]: "<<entry_center[0]<<endl;
    // cout<<""<<endl;
    if (i==25) {
      for (int n1x=0;n1x<nX[0];n1x++) {
	for (int n1y=0;n1y<nY[0];n1y++) {
	  for (int n2x=0;n2x<nX[1];n2x++) {
	    for (int n2y=0;n2y<nY[1];n2y++) {
	      fHBU_x=(recoX[1][n2x]-recoX[0][n1x])*(HBU_z[0]-tReco[2])/(tReco[5]-tReco[2])+recoX[0][n1x];
	      fHBU_y=(recoY[1][n2y]-recoY[0][n2y])*(HBU_z[0]-tReco[2])/(tReco[5]-tReco[2])+recoY[0][n2y];
	      xi->Fill(((fHBU_x-HhitPos[0][0])*(fHBU_x-HhitPos[0][0])+(fHBU_y-HhitPos[0][1])*(fHBU_y-HhitPos[0][1]))/84);

	      // if (n1x==1&&n1y==1&&n2x==1&&n2y==1){
	      // xi->Fill(((HBU_x[0]-HhitPos[0][0])*(HBU_x[0]-HhitPos[0][0])+(HBU_y[0]-HhitPos[0][1])*(HBU_y[0]-HhitPos[0][1]))/84);
	      // 	      }
	    }
	  }
	}
      }

      //      cout<<"nhits : "<<(HBU_x[0]-HhitPos[0][0])*(HBU_x[0]-HhitPos[0][0])+(HBU_y[0]-HhitPos[0][1])*(HBU_y[0]-HhitPos[0][1])<<endl;
      // cout<<"HBU_x : "<<HBU_x[0]<<endl;
      // cout<<"HBU_y : "<<HBU_y[0]<<endl;
      // cout<<"HhitPosX : "<<HhitPos[0][0]<<endl;
      // cout<<"HhitPosY : "<<HhitPos[0][1]<<endl;
      // cout<<"HhitI : "<<HhitI[0]<<endl;
      // cout<<"HhitJ : "<<HhitJ[0]<<endl;
    }
  }

  //  FILE *output = fopen("../txt/shift/shift_6layer_0220_2_2.txt","w");
  for (int layer=0;layer<6;layer++) {
    for (int xaxis=0;xaxis<3;xaxis++) {
      for (int yaxis=0;yaxis<3;yaxis++) {
  	shift[layer]->Fill(xaxis*30+15,yaxis*30+15,count[layer][xaxis][yaxis]/entry_center[layer]);
	//  	fprintf(output,"%d %d %d\t%f\n",layer,xaxis,yaxis,count[layer][xaxis][yaxis]/entry_center[layer]);
      }
    }
  }
  canv1->cd();
  shift[0]->Draw("colztext");
  // sprintf(canvasName,"../png/efficiency/shift/cosmic0220_%d_%d_layer5.png",rangeI+1,rangeJ+1);
  // canv1->SaveAs(canvasName);
  //  fclose(output);

  canv2->cd();
  xi->Draw();
  
  cout<<"count[0][1][2]: "<<count[0][1][2]<<endl;
  cout<<"count[0][2][1]: "<<count[0][2][1]<<endl;
  cout<<"count[0][1][0]: "<<count[0][1][0]<<endl;
  cout<<"count[0][0][1]: "<<count[0][0][1]<<endl;
  cout<<"count[0][1][1]: "<<count[0][1][1]<<endl;
  cout<<"entry_center[0]: "<<entry_center[0]<<endl;

	
}
