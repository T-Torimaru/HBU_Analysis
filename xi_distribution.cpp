using namespace TMath;

void xi_distribution(){
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
  Int_t rangeI = 9;
  Int_t rangeJ = 2;
  Int_t count[6][3][3]={0};
  Double_t HBU_x[64]={0},HBU_y[64]={0},HBU_z[64]={0},fHBU_x,fHBU_y,distance,XI;
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

  TCanvas *canv2 = new TCanvas("c2","c2",696,500);
  //  TH1F *xi =new TH1F("xi^2","Xi^2(EventNumber=0-500);Xi^2;Entry",10000,0,15000);
    TH1F *xi =new TH1F("xi^2","Xi^2(EventNumber=13);Xi^2;Entry",50,0,150);

  // gStyle->SetOptStat(0);
  // gStyle->SetPalette(1);
  //  gStyle->SetOptFit();
  gStyle->SetOptStat(1111);

  for (int i=0;i<500;i++){
    if (i==7294 || i==17572 || i==33381 || i==42705 || i==44465 || i==48873 || i==48979 || i==50060 || i==68587 || i==108446 || i==145457 || i==198356 || i==217244 || i==231350 || i==291337) continue;
    tmatch->GetEntry(i);
    //    if (Hnhit>18 || nPass<5) continue;
    if (Hnhit>18) continue;
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

    for (int j=0;j<Hnhit;j++){
    //   //      if (HADC[j]<650) continue;
      ntp=HhitK[j]-1;
      HBU_z[j]=constHBU_z[ntp];
      // HBU_x[j]=(tReco[3]-tReco[0])*(HBU_z[j]-tReco[2])/(tReco[5]-tReco[2])+tReco[0];
      // HBU_y[j]=(tReco[4]-tReco[1])*(HBU_z[j]-tReco[2])/(tReco[5]-tReco[2])+tReco[1];
    //   xaxisI=HhitI[j]-rangeI-6;
    //   yaxisJ=HhitJ[j]-rangeJ-6;
    }

    if (i==13) {
      //      cout<<"recoX[0][0] : "<<endl;
      // cout<<"recoX[0][1] : "<<recoX[0][1]<<endl;
      for (int n1x=0;n1x<nX[0];n1x++) {
    	for (int n1y=0;n1y<nY[0];n1y++) {
    	  for (int n2x=0;n2x<nX[1];n2x++) {
    	    for (int n2y=0;n2y<nY[1];n2y++) {
	      XI=0;
	      for (int j=0;j<Hnhit;j++) {
    	      fHBU_x=(recoX[1][n2x]-recoX[0][n1x])*(HBU_z[j]-tReco[2])/(tReco[5]-tReco[2])+recoX[0][n1x];
    	      fHBU_y=(recoY[1][n2y]-recoY[0][n1y])*(HBU_z[j]-tReco[2])/(tReco[5]-tReco[2])+recoY[0][n1y];
	      distance=Sqrt((fHBU_x-HhitPos[j][0])*(fHBU_x-HhitPos[j][0])+(fHBU_y-HhitPos[j][1])*(fHBU_y-HhitPos[j][1]));
	      XI+=(distance**2)/84.;
	      }
    	      xi->Fill(XI);
      	    }
      	  }
	}
      }
    }
  }
  canv2->cd();
  //  canv2->SetLogx();
  xi->Draw();
  
}
