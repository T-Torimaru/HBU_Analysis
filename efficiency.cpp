void efficiency(){
  TString matchfile;
  matchfile.Form("../rootfile/cosmic/cosmi_1224_2163.root");
  TFile *fmatch = new TFile(matchfile.Data(), "read");
  TTree *tmatch = (TTree*)fmatch->Get("bigtree");
  char canvasName[256];
  Int_t HhitI[256],HhitJ[256],HhitK[256],Hnhit;
  Float_t HhitPos[256][3];
  Float_t  HADC[256];
  Double_t nph[2][64];
  Int_t    nX[2];
  Int_t    nY[2];
  Double_t recoX[2][8];
  Double_t recoY[2][8];
  Double_t tReco[6];
  Double_t HBUPos[4][3];
  Int_t maxPass;
  Int_t nEntries=tmatch->GetEntries();
  Double_t HBU_z[4]={217.75,261.05,304.35,347.65};
  Double_t HBU_x[4];
  Double_t HBU_y[4];

  tmatch->SetBranchAddress("ahc_nHits", &Hnhit);
  tmatch->SetBranchAddress("ahc_hitI", HhitI);
  tmatch->SetBranchAddress("ahc_hitJ", HhitJ);
  tmatch->SetBranchAddress("ahc_hitK", HhitK);
  tmatch->SetBranchAddress("ahc_hitPos", HhitPos);
  tmatch->SetBranchAddress("ahc_hitEnergy", HADC);
  tmatch->SetBranchAddress("hod_maxPass",&maxPass);
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

  TH1F *efficiencyH = new TH1F("histo","Tile Response Efficiency",12,0,1.2);
  TCanvas *c1 = new TCanvas("c1","c1",696,500);
  Int_t modX,modY;
  Double_t deciX,deciY;
  Double_t valid1=0.2;
  Double_t valid2=29.8;
  Double_t val1,val2;
  Double_t tileEdge[2][2];
  bool eff[4];
  for (int i=0;i<nEntries;i++){
  tmatch->GetEntry(i);
  Double_t count1=0;
  Double_t count2=0;

  for (int j=0;j<4;j++){
    HBU_x[j]=(tReco[3]-tReco[0])*(HBU_z[j]-tReco[2])/(tReco[5]-tReco[2])+tReco[0];
    HBU_y[j]=(tReco[4]-tReco[1])*(HBU_z[j]-tReco[2])/(tReco[5]-tReco[2])+tReco[1];
    if (HBU_x[j]<0) {
      modX=-HBU_x[j];
      deciX=-HBU_x[j]-modX;
    }
    else {
      modX=HBU_x[j];
      deciX=HBU_x[j]-modX;
    }
    if (HBU_y[j]<0) {
      modY=-HBU_y[j];
      deciY=-HBU_y[j]-modY;
    }
    else {
      modY=HBU_y[j];
      deciY=HBU_y[j]-modY;
    }
    modX=modX%30;
    modY=modY%30;
    val1=modX+deciX;
    val2=modY+deciY;
    if (val1<valid1 || val1>valid2 || val2<valid1 || val2>valid2) eff[j]=false;
    else {
      eff[j]=true;
      count1++;
    }
  }
  //  if (i==11){
  for (int k=0;k<Hnhit;k++){
    if (eff[HhitK[k]-1]){
    tileEdge[0][0]=HhitPos[k][0]+15;
    tileEdge[0][1]=HhitPos[k][0]-15;
    tileEdge[1][0]=HhitPos[k][1]+15;
    tileEdge[1][1]=HhitPos[k][1]-15;
    // cout<<HhitK[k]<<endl;
    // cout<<HBU_x[HhitK[k]-1]<<" "<<HBU_y[HhitK[k]-1]<<endl;
    // cout<<tileEdge[0][1]<<" "<<tileEdge[0][0]<<" "<<tileEdge[1][1]<<" "<<tileEdge[1][0]<<endl;    
     if (tileEdge[0][1]<HBU_x[HhitK[k]-1] && tileEdge[0][0]>HBU_x[HhitK[k]-1] && tileEdge[1][0]>HBU_y[HhitK[k]-1] && tileEdge[1][1]<HBU_y[HhitK[k]-1]) {
       count2++;
     }
    }
  }
  Double_t div=count2/count1;
  if (count1!=0){
  efficiencyH->Fill(div);
  }
  }
  c1->cd();
  efficiencyH->Draw();
  sprintf(canvasName,"../png/efficiency_2017Dec.png");
  c1->SaveAs(canvasName);
}

