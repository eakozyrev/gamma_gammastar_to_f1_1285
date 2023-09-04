#include <TH2.h>
#include <TF1.h>
#include <TH1.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TProfile.h"
#include "TFile.h"
#include "TText.h"
#include "TMatrixTSym.h"
#include "TGraphErrors.h"
#include "TGaxis.h"
#include "TMath.h"
#include "TVector.h"
#include "TStyle.h"
#include "TBenchmark.h"
#include <iostream>
#include <fstream>
#include "RooChi2Var.h"
#include <string>
#include "TCut.h"
#include "TFile.h"
#include "TFitter.h"
#include "TTree.h"
#include "THStack.h"
#include "RooPolynomial.h"
#include <complex>
#include "TEfficiency.h"
#include <stdlib.h>
#include <dirent.h>
#include <sys/types.h>
#include "TGraphAsymmErrors.h"

using namespace std;
TFile *newfile0, *newfile1, *newfile2, *newfile3, *newfile4;
TTree *tree0, *tree1, *tree2, *tree3, *tree4, *tree;
TTree *tree_m0_gen, *tree_m1_gen;

string SELECT = "1==1";
double result[10];

double cross_total_m0 = 33.8*2.;//  pb
double cross_total_m1 = 17.4*2.;//  pb
double dcr2[6];

void open(){

  newfile0 = TFile::Open("../eef1/histograms/data.root");//data.root");
  tree0 = (TTree*)newfile0->Get("Tree");

  newfile1 = TFile::Open("../eef1/histograms/mc_m0.root");
  tree1 = (TTree*)newfile1->Get("Tree");

  newfile2 = TFile::Open("../eef1/histograms/mc_m1.root");
  tree2 = (TTree*)newfile2->Get("Tree");

  newfile3 = TFile::Open("../eef1/histograms/mc_f0_eta_m0.root");
  tree3 = (TTree*)newfile3->Get("Tree");

  newfile4 = TFile::Open("../eef1/histograms/mc_f0_eta_m1.root");
  tree4 = (TTree*)newfile4->Get("Tree");

  tree_m0_gen = (TTree*)(TFile::Open("../eef1/histograms/mc_m0.root"))->Get("Treegen");
  tree_m1_gen = (TTree*)(TFile::Open("../eef1/histograms/mc_m1.root"))->Get("Treegen");

}


double dcrs(double *e, double *par){

  double M = 1.285;
  double s = 10.58*10.58;
  double Q2cut_m2 = 0.1/0.0005/0.0005;
  double log1tau = log(s/(M*M+e[0]));
  double x = e[0]/M/M;
  double FF = 1./(1+e[0]/0.7/0.7);
  complex<double> ph = 1.+complex<double>(cos(par[0]),sin(par[0]))*(1.+x)*par[1];
  FF = (2. + x*(ph.real()*ph.real() + ph.imag()*ph.imag()))/(2.+x)/(1.+x)/(1.+x)/(1+x/0.36)/(1+x/0.36);
  double Ggg = 3./1000./1000.;
  double res = 24./137.15/137.15/pow(M,5.)*Ggg;
  res = res*FF*FF*(log(Q2cut_m2)*(log1tau - 7./4.) + log1tau*log1tau - 3.*log1tau - 3.14*3.14/6. + 23./8. + 0.5*e[0]/M/M*(log(Q2cut_m2)*(log1tau-1.5) + log1tau*log1tau - 2.5*log1tau - 3.14*3.14/6. + 19./8.));

  return res*pow(197.327,2.)*pow(10.,6.);

}


// TF1 *fcr = new TF1("fcr",dcrs,0.1,30,0);


double prob_ma0(double *s, double *par){

  double mpl = s[0];
  double mmi = s[1];
  double wide = par[2];
  TCut cut = Form((SELECT + "&& ifsig==1 && ma0p > %g && ma0p < %g && ma0m > %g && ma0m < %g").c_str(),mpl-wide,mpl+wide,mmi-wide,mmi+wide);
  double norm = 1;
  double N = 1;
  //cout << tree1->GetEntries() << " " << par[1] << " ";
  if(par[1] == 0){norm = 1.; N = tree0->GetEntries(cut);} 
  if(par[1] == 1){norm = tree1->GetEntries((TCut)SELECT.c_str()); N = tree1->GetEntries(cut);}
  if(par[1] == 2){norm = tree2->GetEntries((TCut)SELECT.c_str()); N = tree2->GetEntries(cut);}
  if(par[1] == 3){norm = tree3->GetEntries((TCut)SELECT.c_str()); N = tree3->GetEntries(cut);}
  if(par[1] == 4){norm = tree4->GetEntries((TCut)SELECT.c_str()); N = tree4->GetEntries(cut);}
  if(par[1] == 5){norm = 1.;N = 0.5*tree0->GetEntries(Form((SELECT + " && ifsig==0 && ma0p > %g && ma0p < %g && ma0m > %g && ma0m < %g").c_str(),mpl-wide,mpl+wide,mmi-wide,mmi+wide));} 
  //cout << N << " " << norm << endl;
  return par[0]*N/norm;
}


/*
open()
TF2 *u = new TF2("u",prob_ma0,0.6,1.2,0.6,1.2,3)
u->SetParameters(1,0,0.05);
u->Draw("box");
*/



double diff_cr_sec(double mode, double par0, double par1){
  if(mode == 0)return (double)tree_m0_gen->GetEntries(Form("Q2gen > %g && Q2gen < %g",par0,par1))/tree_m0_gen->GetEntries()*cross_total_m0;
  else if(mode == 1)return (double)tree_m1_gen->GetEntries(Form("Q2gen > %g && Q2gen < %g",par0,par1))/tree_m1_gen->GetEntries()*cross_total_m1;
  else return 0.;
}



double prob_cosa0(double *s, double *par){

  double mpl = s[0];
  double mmi = s[1];
  double wide = par[2];
  TCut cut = Form((SELECT + " && ifsig==1 && cosa0_p > %g && cosa0_p < %g && cosa0_m > %g && cosa0_m < %g").c_str(),mpl-wide,mpl+wide,mmi-wide,mmi+wide);
  double norm = 1;
  double N = 1;
  

  double alphamc = (double)tree1->GetEntries((SELECT + " && ifsig==0").c_str())/tree1->GetEntries((SELECT + " && ifsig==1").c_str());
    
  if(par[1] == 0){norm = 1.; N = tree0->GetEntries(cut);} 
  if(par[1] == 1){norm = tree1->GetEntries((TCut)(SELECT+"&& ifsig==1").c_str()); N = tree1->GetEntries(cut)*(1. - alphamc/2.);}
  if(par[1] == 2){norm = tree2->GetEntries((TCut)(SELECT+"&& ifsig==1").c_str()); N = tree2->GetEntries(cut)*(1. - alphamc/2.);}
  if(par[1] == 3){norm = tree3->GetEntries((TCut)(SELECT+"&& ifsig==1").c_str()); N = tree3->GetEntries(cut)*(1. - alphamc/2.);}
  if(par[1] == 4){norm = tree4->GetEntries((TCut)(SELECT+"&& ifsig==1").c_str()); N = tree4->GetEntries(cut)*(1. - alphamc/2.);}
  if(par[1] == 5){
    norm = 1.;
    N = tree0->GetEntries(Form((SELECT + " && ifsig==0 && cosa0_p > %g && cosa0_p < %g && cosa0_m > %g && cosa0_m < %g").c_str(),mpl-wide,mpl+wide,mmi-wide,mmi+wide));
    //N =N*(1. - alphamc);
    //cout << "alphamc ===================================  " << alphamc << endl;
  }
  
 
  return par[0]*N/norm;
}

void draw_bkg(){
  SELECT = Form("Q2 > %g && Q2 < %g",6.,20.);
  
  int ncos = 5;
  TH2D *hbkg = new TH2D("hbkg","hbkg",ncos,-1,1,ncos,-1,1);
  hbkg->SetXTitle("cos(#theta #pi+)");
  hbkg->SetYTitle("cos(#theta #pi-)");
  double cos_start = -1;
  double wide = (1.-cos_start)/ncos;
  
  for(int i = 0; i< ncos; i++){
    for(int j = 0; j < ncos; j++){
      double mpl = cos_start + i*wide + wide/2.;
      double mmi = cos_start + j*wide + wide/2.;
      double BG =  0.5*tree1->GetEntries(Form((SELECT + " && ifsig==0 && cosa0_p > %g && cosa0_p < %g && cosa0_m > %g && cosa0_m < %g").c_str(),mpl-wide/2.,mpl+wide/2.,mmi-wide/2.,mmi+wide/2.));
      double N = tree1->GetEntries(Form((SELECT + " && ifsig==1 && cosa0_p > %g && cosa0_p < %g && cosa0_m > %g && cosa0_m < %g").c_str(),mpl-wide/2.,mpl+wide/2.,mmi-wide/2.,mmi+wide/2.));
      hbkg->SetBinContent(i+1,j+1, N > 0 ? BG/N : 0);
      
    }
  }

  hbkg->Draw("boxtext");
  
}


double chi2(double *par, double *par1){

  double res = 0;
  int nma0 = 6;
  double ma0_start = 0.65;
  double wide_ma0 = (1.15-ma0_start)/nma0;
  int ncos = 5;
  double cos_start = -1;
  double wide_cos = (1.-cos_start)/ncos;
  for(int i = 0; i< nma0; i++){
    for(int j = 0; j < nma0; j++){
      double s[] = {ma0_start + i*wide_ma0 + wide_ma0/2., ma0_start + j*wide_ma0 + wide_ma0/2.};
      double parm[] = {1,0,wide_ma0/2.};
      double btv = prob_ma0(s,parm);
      if(btv <= 0.)continue;
      parm[1] = 1;
      double nmc = par[0]*par[1]*(1.-par[2])*prob_ma0(s,parm);
      parm[1] = 2;
      nmc = nmc + par[0]*(1.-par[1])*(1.-par[2])*prob_ma0(s,parm);
      parm[1] = 3;
      nmc = nmc + par[0]*par[1]*par[2]*prob_ma0(s,parm);
      parm[1] = 4;
      nmc = nmc + par[0]*(1.-par[1])*par[2]*prob_ma0(s,parm);
      parm[1] = 5;
      nmc = nmc + prob_ma0(s,parm);
      //result = result - log(TMath::Poisson(nmc,btv));
    }
  }
  double Bkg = 0;
  for(int i = 0; i< ncos; i++){
    for(int j = 0; j < ncos; j++){
      double s[] = {cos_start + i*wide_cos + wide_cos/2., cos_start + j*wide_cos + wide_cos/2.};
      double parm[] = {1,0,wide_cos/2.};
      double btv = prob_cosa0(s,parm);
      if(btv <= 0.)continue;
      parm[1] = 1;
      double nmc = par[0]*par[1]*(1.-par[2])*prob_cosa0(s,parm);
      parm[1] = 2;
      nmc = nmc + par[0]*(1.-par[1])*(1.-par[2])*prob_cosa0(s,parm);
      parm[1] = 3;
      nmc = nmc + par[0]*par[1]*par[2]*prob_cosa0(s,parm);
      parm[1] = 4;
      nmc = nmc + par[0]*(1.-par[1])*par[2]*prob_cosa0(s,parm);
      parm[1] = 5;
      double bkg = prob_cosa0(s,parm);
      nmc = nmc + bkg/2.;
      Bkg+=bkg;
      //cout << "bkg = " << bkg << "            Bkg = " << Bkg <<  "   par[0] = " << par[0] << "    nmc = " << nmc << endl;
      res = res - 2.*log(TMath::Poisson(nmc,btv));
    }
  }

  return res;
    
}

/*
open()
TF2 *u = new TF2("u",prob_cosa0,-1,1.,-1,1.,2)
u->SetParameters(1,0);
u->Draw("box");
*/


void minuitCMD3(int &nDim,double*gout,double&result,double par[],int flg){
  result= chi2(par,par);
}

void minim(){
  open();
  for(int i=0; i < 10; i++){result[i]=0.;}
  double parm[] = {1,0,10.};
  double s[] = {0.7,0.7};
  double N_data = prob_cosa0(s,parm);
  parm[1] = 5;
  cout << "N_data = " << N_data << endl;
  N_data = N_data - prob_cosa0(s,parm)/2.;
  cout << "N_data = " << N_data << endl;
  result[9] = N_data;
  TFitter*minimizer=new TFitter(3);
    double minimum = 0;
    double resultssss[1000];
    minimizer->SetFitMethod("loglikelihood");
    minimizer->SetParameter(0,"lambda1 ", N_data, 0.0001,10,10000.);
    minimizer->SetParameter(1,"lambda2 ", 0.3, 0.0001,0.0,1.);
    minimizer->SetParameter(2,"lambda3 ", 0., 0.00001,-1,1.);
    minimizer->FixParameter(0);
    minimizer->FixParameter(2);
    minimizer->SetFCN(minuitCMD3);
    double arglist[100];
    arglist[0] = 4000000; // number of function calls
    arglist[1] = 0.00001; // tolerance
    TFitResultPtr r = minimizer->ExecuteCommand("MINOS",arglist,2);

    result[0] = minimizer->GetParameter(1);
    result[1] = minimizer->GetParError(1);
    // cout << r->Print() << endl;
}


int script_minim(){  
  int nrun = 6;
  double en0[] = {2,4,5,6,7,10,20.,30};
  //double en0[] = {2,6,15.,30};
  ofstream stream("results/m0_fraction.dat");
  
  for(int i = 0; i < nrun; i++){
    cout << "======================================================= " << endl;
    cout << " START MINIMIZATION FOR " << en0[i] << " --  " << en0[i+1] << endl;
    SELECT = Form("Q2 > %g && Q2 < %g && hel < 0.9",en0[i],en0[i+1]);

    minim();
    stream << en0[i] << " " << en0[i+1] << " " << result[9] << " " << result[0] << " " << result[1] << endl;
  }
      return 1;
}



void draw_fraction(){

 double en[100000],den[100000],cr[100000],dcr[100000];
 double dcrl[1000],dcrh[1000];
   int nrun = 0;
   ifstream stream("results/m0_froaction/m0_fraction_0.dat");
   while(stream.eof()==0){
	stream >> en[nrun] >> den[nrun] >> cr[nrun] >> cr[nrun] >> dcr[nrun];
	den[nrun] = ( den[nrun]  - en[nrun])/2.;
	en[nrun] = en[nrun] + den[nrun];
	if(stream.eof()==1)break;
        nrun++;
   }

   vector<double> model_un[6];

   struct dirent *entry;
   DIR *dir = opendir("results/m0_froaction/");
   if (dir == NULL) {
     return;
   }
   int num = 0;
   while ((entry = readdir(dir)) != NULL) {
     string newfile = entry->d_name;
     if(newfile.find(".dat") > 1000)continue;
     ifstream stream("results/m0_froaction/"+newfile);
     cout << newfile << endl;
     num = 0;
     while(stream.eof()==0){
       double btv;
       stream >> btv >> btv >> btv >> btv; //  >> dcr[nrun];
       if(stream.eof()==1)break;
       model_un[num].push_back(cr[num] - btv);
       cout << num << " " << cr[num] << " " << cr[num] - btv << endl;
       stream >> btv;
       num++;
       }
   }
     
    
     for(int i = 0; i < 6; i++){
       cout << i << " " << dcr2[i] << " " <<  model_un[i].at(0)<< endl;
       vector<double>::iterator it; 
       double btv = 0;
       for(it = model_un[i].begin(); it != model_un[i].end(); it++)    {
	 btv += (*it)*(*it);
	 
       }
       dcr2[i] = sqrt(btv); 
     }
   
   TCanvas *s = new TCanvas();
   TH1F *frd  = s->DrawFrame(en[0]-10,-100.,en[nrun-1]+10,10000.);
   frd->SetXTitle("Q^{2}, GeV^{2}");
   frd->SetYTitle("#sigma_{m=0}/#sigma_{total}");

   
   for(int i = 0; i < nrun; i++){
     dcrh[i] = cr[i] + dcr[i] > 1 ?  1-cr[i] : dcr[i];
     dcrl[i] = cr[i] - dcr[i] < 0 ? cr[i] : dcr[i];
     dcr2[i] = cr[i] + dcr2[i] > 1 ?  1-cr[i] : dcr2[i];
     dcr2[i] = cr[i] - dcr2[i] < 0 ? cr[i] : dcr2[i];     
   }
   TGraphAsymmErrors *Cross  = new TGraphAsymmErrors(nrun,en,cr,den,den,dcrl,dcrh);
   Cross->SetMarkerColor(2);
   Cross->SetMarkerStyle(20);
   Cross->SetLineColor(2);
   Cross->SetLineWidth(2.);
   Cross->SetTitle("");
   // Cross->Draw("P");
   
   TGraphAsymmErrors *Cross2  = new TGraphAsymmErrors(nrun,en,cr,den, den, dcr2, dcr2);
   Cross2->SetMarkerColor(2);
   Cross2->SetMarkerStyle(20);
   Cross2->SetLineColor(2);
   Cross2->SetLineWidth(2.);
   Cross2->SetTitle("");
   Cross2->SetFillColor(3);
   Cross2->SetFillStyle(3001);
   Cross2->Draw("a2");
   //Cross->SetXTitle("Q^{2}, GeV^{2}");
   //Cross->SetYTitle("#sigma_{0}/#sigma_{total}");
   //Cross2->Draw("p");
   
   //   Cross->Draw("P");

   //TGraphAsymmErrors* gme = new TGraphAsymmErrors("gme", "TGraphMultiErrors Example", nrun, en, cr, den, den, dcr, dcr);
   /* gme->AddYError(5, dcr2, dcr2);
   gme->SetMarkerStyle(20);
   gme->SetLineColor(kRed);
   gme->GetAttLine(0)->SetLineColor(kRed);
   gme->GetAttLine(1)->SetLineColor(kBlue);
   gme->GetAttFill(1)->SetFillStyle(0);
 
   gme->Draw("a p s ; ; 5 s=0.5");
   */

    for(int i = 0; i < nrun; i++){
     dcr2[i] = dcr2[i]/cr[i];
   }

    Cross->Draw("P");
}



void compare_cos(string filee, string filemc0, string filemc1, double factorm1, string variable,string cutt = "1==1"){

  TCanvas *s = new TCanvas();
  TCut cut = (cutt + " && ifsig==1").c_str();
  TFile *newfile = TFile::Open((filee).c_str());
  TTree* tree = (TTree*)newfile->Get("Tree");
  TH1D *htest = new TH1D("htest","htest",1000000,-10000,10000);
  tree->Draw(("cosa0_p >> htest"),cut);
  tree->Draw(("cosa0_m >> htest"),cut);  
  double start = -10000;
  double end = 100000;
  for(int i = 1; i <= 1000000; i++){
	  if(htest->GetBinContent(i) > 0)start = htest->GetBinCenter(i);
  }
  for(int i = 1000000; i > 0; i--){
	  if(htest->GetBinContent(i) > 0)end = htest->GetBinCenter(i);
  }
  start = -1.;
  end = 1.;
  TH1D *hm2ph = new TH1D("hm2ph","hm2ph",40,start,end);
  tree->Draw(Form("cosa0_p >> hm2ph(40,%g,%g)",start,end),cut);
  tree->Draw(Form("cosa0_m >> hm2ph(40,%g,%g)",start,end),cut);
  TH1D *hm2ph_b = new TH1D("hm2ph_b","hm2ph_b",40,start,end);
  TCut cut_b = (cutt + "&& ifsig==0").c_str();
  tree->Draw(Form("cosa0_p >> hm2ph_b(40,%g,%g)",start,end),cut_b);
  tree->Draw(Form("cosa0_m >> hm2ph_b(40,%g,%g)",start,end),cut_b);
  //hm2ph->Add(hm2ph_b,-0.5);

  
  double notm = tree->GetEntries(cut) - 0.5*tree->GetEntries(cut_b);
  cout << notm << endl;
  TFile *newfilemc0 = TFile::Open(filemc0.c_str());
  TTree* treemc0 = (TTree*)newfilemc0->Get("Tree");
  TH1D *hm2ph_mc = new TH1D("hm2ph_mc","hm2ph_mc",40,start,end);
  TH1D *hm2ph_mc0 = new TH1D("hm2ph_mc0","hm2ph_mc0",40,start,end);
  treemc0->Draw(Form("cosa0_p >> hm2ph_mc0(40,%g,%g)",start,end),cut);
  treemc0->Draw(Form("cosa0_m >> hm2ph_mc0(40,%g,%g)",start,end),cut);
  TFile *newfilemc1 = TFile::Open(filemc1.c_str());
  TTree* treemc1 = (TTree*)newfilemc1->Get("Tree");
  TH1D *hm2ph_mc1 = new TH1D("hm2ph_mc1","hm2ph_mc1",40,start,end);
  treemc1->Draw(Form("cosa0_p >> hm2ph_mc1(40,%g,%g)",start,end),cut);
  treemc1->Draw(Form("cosa0_m >> hm2ph_mc1(40,%g,%g)",start,end),cut);

  hm2ph_mc0->SetNormFactor((double)notm*(1.-factorm1));
  hm2ph_mc1->SetNormFactor((double)notm*factorm1);

  hm2ph_mc->Add(hm2ph_b,0.5);
  hm2ph_mc->Add(hm2ph_mc0,1);
  hm2ph_mc->Add(hm2ph_mc1,1);
  
  hm2ph->Draw("e");
  hm2ph->SetLineColor(1);
  hm2ph->SetLineWidth(2.);
  hm2ph_mc->SetFillColor(4);
  hm2ph_mc->SetFillStyle(3002);
  hm2ph_mc->Draw("same");
  hm2ph_b->SetNormFactor(0.5*hm2ph_b->GetEntries());hm2ph_b->Draw("same");
  hm2ph_b->SetFillStyle(3001);hm2ph_b->SetLineColor(13);hm2ph_b->SetFillColor(13);

  hm2ph_mc1->SetLineColor(1);hm2ph_mc1->SetLineWidth(3);hm2ph_mc1->Draw("same");
  hm2ph_mc0->SetLineColor(2);hm2ph_mc0->SetLineWidth(3);hm2ph_mc0->Draw("same");
  hm2ph->Draw("esame");

  //hm2ph->SetAxisRange(0,50,"Y");
  hm2ph->SetXTitle(variable.c_str());

  s->SaveAs(("plots/"+(variable + cutt)+".png").c_str());
  s->Close();
  newfile->Close();
}



void compare(string filee, string filemc0, string filemc1, double factorm0, string variable,string cutt = "1==1"){

  TCanvas *s = new TCanvas();
  TCut cut = (cutt + " && ifsig==1").c_str();
  TFile *newfile = TFile::Open((filee).c_str());
  TTree* tree = (TTree*)newfile->Get("Tree");
  TH1D *htest = new TH1D("htest","htest",1000000,-10000,10000);
  tree->Draw((variable + " >> htest").c_str(),cut);
  double start = -10000;
  double end = 100000;
  for(int i = 1; i <= 1000000; i++){
	  if(htest->GetBinContent(i) > 0)start = htest->GetBinCenter(i);
  }
  for(int i = 1000000; i > 0; i--){
	  if(htest->GetBinContent(i) > 0)end = htest->GetBinCenter(i);
  }
  start = 0.6;
  end = 1.3;
  TH1D *h1;
  TH1D *hm2ph = new TH1D("hm2ph","hm2ph",40,start,end);
  tree->Draw((variable + Form(" >> h1(40,%g,%g)",start,end)).c_str(),cut);
  h1 = (TH1D*)gPad->GetPrimitive("h1");
  hm2ph->Add(h1);
  cout << "here0" << endl;
  TH1D *hm2ph_b = new TH1D("hm2ph_b","hm2ph_b",40,start,end);
  TCut cut_b = (cutt + "&& ifsig==0").c_str();
  tree->Draw((variable + Form(" >> h1(40,%g,%g)",start,end)).c_str(),cut_b);
  hm2ph_b->Add(h1);

  
  double notm = tree->GetEntries(cut) - 0.5*tree->GetEntries(cut_b);
  cout << notm << endl;
  TFile *newfilemc0 = TFile::Open(filemc0.c_str());
  TTree* treemc0 = (TTree*)newfilemc0->Get("Tree");
  TH1D *hm2ph_mc = new TH1D("hm2ph_mc","hm2ph_mc",40,start,end);
  TH1D *hm2ph_mc0 = new TH1D("hm2ph_mc0","hm2ph_mc0",40,start,end);
  treemc0->Draw((variable + Form(" >> h1(40,%g,%g)",start,end)).c_str(),cut);
  h1 = (TH1D*)gPad->GetPrimitive("h1");
  hm2ph_mc0->Add(h1);
  TFile *newfilemc1 = TFile::Open(filemc1.c_str());
  TTree* treemc1 = (TTree*)newfilemc1->Get("Tree");
  TH1D *hm2ph_mc1 = new TH1D("hm2ph_mc1","hm2ph_mc1",40,start,end);
  treemc1->Draw((variable + Form(" >> h1(40,%g,%g)",start,end)).c_str(),cut);
  h1 = (TH1D*)gPad->GetPrimitive("h1");
  hm2ph_mc1->Add(h1);

  hm2ph_mc0->SetNormFactor((double)notm*factorm0);
  hm2ph_mc1->SetNormFactor((double)notm*(1.-factorm0));

  hm2ph_mc->Add(hm2ph_b,0.5);
  hm2ph_mc->Add(hm2ph_mc0,1);
  hm2ph_mc->Add(hm2ph_mc1,1);
  
  hm2ph->Draw("e");
  hm2ph->SetLineColor(1);
  hm2ph->SetLineWidth(2.);
  hm2ph_mc->SetFillColor(4);
  hm2ph_mc->SetFillStyle(3002);
  hm2ph_mc->Draw("same");
  hm2ph_b->SetNormFactor(0.5*hm2ph_b->GetEntries());hm2ph_b->Draw("same");
  hm2ph_b->SetFillStyle(3001);hm2ph_b->SetLineColor(13);hm2ph_b->SetFillColor(13);

  hm2ph_mc1->SetLineColor(1);hm2ph_mc1->SetLineWidth(3);hm2ph_mc1->Draw("same");
  hm2ph_mc0->SetLineColor(2);hm2ph_mc0->SetLineWidth(3);hm2ph_mc0->Draw("same");
  hm2ph->Draw("esame");

  //hm2ph->SetAxisRange(0,50,"Y");
  hm2ph->SetXTitle(variable.c_str());

  s->SaveAs(("plots/"+variable + cutt+".png").c_str());
  //s->Close();
  //newfile->Close();
}



TH1D *res_hist(TFile *filee, string variable, TCut cutt = "1==1"){

  double start = -1.;
  double end = 1.;
  TH1D *hm2ph = new TH1D("hm2ph","hm2ph",40,start,end);
  TTree* treemc0 = (TTree*)filee->Get("Tree");
  treemc0->Draw(Form("cosa0_p >> hm2ph(40,%g,%g)",start,end),cutt);

  return hm2ph;
 
}

void compare_cos_all(string filee, string filemc0, string filemc1, double factorm0, string variable,string cutt = "1==1",string name_fs=""){

  TCanvas *s = new TCanvas();
  TCut cut = (cutt + " && ifsig==1 && cosa0_p > -0.99999 && cosa0_p < 0.99999 && cosa0_m > -0.99999 && cosa0_m < 0.99999").c_str();
  TFile *newfile = TFile::Open((filee).c_str());
  TTree* tree = (TTree*)newfile->Get("Tree");
  TH1D *htest = new TH1D("htest","htest",1000000,-10000,10000);
  tree->Draw(("cosa0_p >> htest"),cut);
  tree->Draw(("cosa0_m >> htest"),cut);  
  double start = -0.99999;
  double end = 0.99999;
  int nbin = 5;

  TH1D *hm2ph = new TH1D("","",nbin,start,end);
  tree->Draw(Form("cosa0_p >> h1(%i,%g,%g",nbin,start,end),cut);
  TH1D *h1 = (TH1D*)gPad->GetPrimitive("h1");
  hm2ph->Add(h1);
  tree->Draw(Form("cosa0_m >> h1(%i,%g,%g)",nbin,start,end),cut);
  h1 = (TH1D*)gPad->GetPrimitive("h1");
  hm2ph->Add(h1);

  TH1D *hm2ph_b = new TH1D("hm2ph_b","hm2ph_b",nbin,start,end);
  TCut cut_b = (cutt + "&& ifsig==0").c_str();
  tree->Draw(Form("cosa0_p >> h1(%i,%g,%g)",nbin,start,end),cut_b);
  h1 = (TH1D*)gPad->GetPrimitive("h1");
  hm2ph_b->Add(h1);
  tree->Draw(Form("cosa0_m >> h1(%i,%g,%g)",nbin,start,end),cut_b);
  h1 = (TH1D*)gPad->GetPrimitive("h1");
  hm2ph_b->Add(h1);
  
  //hm2ph->Add(hm2ph_b,-0.5);

  TH1D *hm2ph_mc = new TH1D("hm2ph_mc","hm2ph_mc",nbin,start,end);

  
  TH1D *hm2ph_mc0 = new TH1D("hm2ph_mc0","hm2ph_mc0",nbin,start,end);
  TFile *newfilemc0 = TFile::Open(filemc0.c_str());
  TTree* treemc0 = (TTree*)newfilemc0->Get("Tree");
  treemc0->Draw(Form("cosa0_p >> h1(%i,%g,%g)",nbin,start,end),cut);
  h1 = (TH1D*)gPad->GetPrimitive("h1");
  hm2ph_mc0->Add(h1);
  treemc0->Draw(Form("cosa0_m >> h1(%i,%g,%g)",nbin,start,end),cut);
  h1 = (TH1D*)gPad->GetPrimitive("h1");
  hm2ph_mc0->Add(h1);  
  TFile *newfilemc1 = TFile::Open(filemc1.c_str());
  TTree* treemc1 = (TTree*)newfilemc1->Get("Tree");
  TH1D *hm2ph_mc1 = new TH1D("hm2ph_mc1","hm2ph_mc1",nbin,start,end);
  treemc1->Draw(Form("cosa0_p >> h1(%i,%g,%g)",nbin,start,end),cut);
  h1 = (TH1D*)gPad->GetPrimitive("h1");
  hm2ph_mc1->Add(h1);
  treemc1->Draw(Form("cosa0_m >> h1(%i,%g,%g)",nbin,start,end),cut);
  h1 = (TH1D*)gPad->GetPrimitive("h1");
  hm2ph_mc1->Add(h1);


  double factor = treemc0->GetEntries(cut_b)/treemc0->GetEntries(cut);
  factor = 0.5*(1.-factor);
  double notm = hm2ph->GetEntries() - factor*hm2ph_b->GetEntries();
  cout << notm << endl;
  
  hm2ph_mc0->SetNormFactor((double)notm*factorm0);
  hm2ph_mc1->SetNormFactor((double)notm*(1.-factorm0));

  hm2ph_mc->Add(hm2ph_b,factor);
  hm2ph_mc->Add(hm2ph_mc0,1);
  hm2ph_mc->Add(hm2ph_mc1,1);

  hm2ph->Draw("e");
  hm2ph->SetLineColor(1);
  hm2ph->SetLineWidth(2.);
  hm2ph_mc->SetFillColor(4);
  hm2ph_mc->SetFillStyle(3002);
  hm2ph_mc->Draw("same");
  hm2ph_b->SetNormFactor(factor*hm2ph_b->GetEntries());hm2ph_b->Draw("same");
  hm2ph_b->SetFillStyle(3001);hm2ph_b->SetLineColor(13);hm2ph_b->SetFillColor(13);

  hm2ph_mc1->SetLineColor(1);hm2ph_mc1->SetLineWidth(3);hm2ph_mc1->Draw("same");
  hm2ph_mc0->SetLineColor(2);hm2ph_mc0->SetLineWidth(3);hm2ph_mc0->Draw("same");
  hm2ph->GetXaxis()->SetRange(0,hm2ph->GetMaximum());
  hm2ph->Draw("esame");

  //hm2ph->SetAxisRange(0,50,"Y");
  hm2ph->SetXTitle(variable.c_str());

  s->SaveAs(("plots/"+name_fs+".png").c_str());
  s->Close();
  newfile->Close();
  
}




void script_draw_cos(){

  ifstream stream_ratio("results/m0_froaction/m0_fraction_0.dat");  
  double factorm0 = 0.54;double dfactorm1;
  double btv;
  double ena,enb,Nev,dNev;
  string name_fs = "";
  int i = 0;
  while(stream_ratio.eof()==0){
    stream_ratio >> ena >> enb >> btv >>  factorm0 >> dfactorm1;
    cout << factorm0 << " " << ena << " " << enb << endl;
    //    factorm0 = 0.5;
    name_fs = Form("cosa0_fit_%i",i);
    i++;
    compare_cos_all("histograms/data.root","histograms/mc_m0.root","histograms/mc_m1.root",factorm0,"cos(#theta_{#pi^{#pm}})",Form("Q2 > %g && Q2 < %g && hel < 0.9",ena,enb),name_fs);
    //   compare("histograms/data.root","histograms/mc_m0.root","histograms/mc_m1.root",factorm1,"cosa0_m",Form("Q2 > %g && Q2 < %g",ena,enb));    
  }
}


void cross_section(){
  ifstream stream("../eef1/results/nevents.dat");
  ifstream stream_ratio("../eef1/results/m0_froaction/m0_fraction_0.dat");
  ifstream stream_rad("../eef1/results/rad.dat");
  ofstream streamof("../eef1/results/cross.dat");
  double crossmc_m0, crossmc_m1;

  ofstream stream_table_cr0("../eef1/results/table_cr0.dat");
  ofstream stream_table_cr1("../eef1/results/table_cr1.dat");
  
  double  effic_Q2(double, double, int, bool);
  double factorm0 = 0.54;double dfactorm0;
  double btv;
  while(stream.eof()==0){
    double ena,enb,Nev,dNev;
    stream >> ena >> enb >> Nev >> dNev;
    stream_ratio >> btv >> btv >> btv >>  factorm0 >> dfactorm0;
    if(stream.eof()==1)break;
    cout << "factorm0 = " << factorm0 << endl;
    double rad;
    stream_rad >> rad >> rad >> rad;
    cout << rad << endl;
    double Nm0 = factorm0*Nev;
    double dNm0 = sqrt(pow(factorm0*dNev,2.) + pow(dfactorm0*Nev,2.));
    double Nm1 = (1.-factorm0)*Nev;
    double dNm1 = sqrt(pow((1.-factorm0)*dNev,2.) + pow(dfactorm0*Nev,2.));
    TFile *newfilem0 = TFile::Open("../eef1/histograms/mc_m0.root");
    TTree* treem0 = (TTree*)newfilem0->Get("Tree");
    TTree* treem0gen = (TTree*)newfilem0->Get("Treegen");
    double partm0 = treem0gen->GetEntries(Form("Q2gen > %g && Q2gen < %g",ena,enb))/(double)treem0gen->GetEntries();
    double dpartm0 = sqrt(partm0*(1. - partm0)/(double)treem0gen->GetEntries());
    TFile *newfilem1 = TFile::Open("../eef1/histograms/mc_m1.root");
    TTree* treem1 = (TTree*)newfilem1->Get("Tree");
    TTree* treem1gen = (TTree*)newfilem1->Get("Treegen");
    double partm1 = treem1gen->GetEntries(Form("Q2gen > %g & Q2gen < %g",ena,enb))/(double)treem1gen->GetEntries();
    double dpartm1 = sqrt(partm1*(1. - partm1)/(double)treem1gen->GetEntries());
    double cross_mc_m0 = partm0*cross_total_m0*1000.;
    double dcross_mc_m0 = dpartm0/partm0*cross_mc_m0;
    double cross_mc_m1 = partm1*cross_total_m1*1000.;
    double dcross_mc_m1 = dpartm1/partm1*cross_mc_m1;
    cout << "Here1" << endl;
    double eff_m0 = effic_Q2(ena,enb,0,false);
    cout << "Here2" << endl;
    double deff_m0 = sqrt(eff_m0*(1. - eff_m0)/(double)treem0gen->GetEntries(Form("Q2gen > %g && Q2gen < %g",ena,enb)));
    cout << "Here3" << endl;
    double eff_m1 = effic_Q2(ena,enb,1,false);
    cout << "Here4" << endl;
    double deff_m1 = sqrt(eff_m1*(1.-eff_m1)/(double)treem1gen->GetEntries(Form("Q2gen > %g && Q2gen < %g",ena,enb)));
cout << "Here5" << endl;
    double eff_m0_true = effic_Q2(ena,enb,0,true);
    double eff_m1_true = effic_Q2(ena,enb,1,true);
    cout << "Here6" << endl;
    double deff_m0_true = sqrt(eff_m0_true*(1. - eff_m0_true)/(double)treem0gen->GetEntries(Form("Q2gen > %g && Q2gen < %g",ena,enb)));
    double deff_m1_true = sqrt(eff_m1_true*(1. - eff_m1_true)/(double)treem1gen->GetEntries(Form("Q2gen > %g && Q2gen < %g",ena,enb)));
    cout << "Here7" << endl;
    double cross_m0 = Nm0/eff_m0_true/469./0.348/rad;
    double dcross_m0 = dNm0/Nm0*cross_m0;
    cout << Nm0 << " " << dNm0  << " " << factorm0 << " " << Nev << " eff_m0 = " << eff_m0 <<   endl;
    double cross_m1 = Nm1/eff_m1_true/469./0.348/rad;
    double dcross_m1 = dNm1/Nm1*cross_m1;
    
    double ff_m0 = sqrt(cross_m0/cross_mc_m0);//form_factor_1((ena+enb)/2.)
    double dff_m0 = sqrt((cross_m0 +dcross_m0)/cross_mc_m0) - ff_m0;
    double ff_m1 = sqrt(cross_m1/cross_mc_m1);
    double dff_m1 = sqrt((cross_m1 + dcross_m1)/cross_mc_m1) - ff_m1;
    
    cout << ena << " " << enb << " " << cross_m0 << " " << dcross_m0 << " "  << cross_m1 << " " << dcross_m1 << " " << ff_m0 << " " << dff_m0 << " "  << ff_m1 << "  " << dff_m1 << endl;
    streamof <<  ena << " " << enb << " " << cross_m0 << " " << dcross_m0 << " "  << cross_m1 << " " << dcross_m1 << " " << ff_m0 << " " << dff_m0 << " "  << ff_m1 << "  " << dff_m1 << endl;

    stream_table_cr0 <<  ena << "$\\div$" << enb << " & " << ((int)(Nm0*10.))/10. << " $\\pm$ " << ((int)(dNm0*10.))/10. << " & " << ((int)(eff_m0*100000.))/1000. << " $\\pm$ " << ((int)(deff_m0*100000.))/1000. << " & " << ((int)(eff_m0_true*100000.))/1000. << " $\\pm$ " << ((int)(deff_m0_true*100000.))/1000. << " & " << (int)cross_mc_m0 << " $\\pm$ " << (int)dcross_mc_m0 <<  " & " << ((int)(cross_m0*100.))/100 << " $\\pm$ " << ((int)(dcross_m0*100.))/100 << " & " << ((int)(rad*100.))/100. <<  " \\\\" << endl;

    stream_table_cr1 << ena << "$\\div$" << enb << " & " <<  Nm1 << " $\\pm$ " << dNm1 << " & " << ((int)(eff_m1*100000.))/1000. << " $\\pm$ " << ((int)(deff_m1*100000.))/1000.  << " & " <<  ((int)(eff_m1_true*100000.))/1000. << " $\\pm$ " << ((int)(deff_m1_true*100000.))/1000.  << " & " << (int)cross_mc_m1 << " $\\pm$ " << (int)dcross_mc_m1 <<  " & " << ((int)(cross_m1*100.))/100 << " $\\pm$ " << ((int)(dcross_m1*100.))/100 << " & " << ((int)(rad*100.))/100. << " \\\\" << endl;
    
  }
  
}

void draw_cross_sect_diff(int mode = 0){


  double ss[]={0.02,0.1,0.4,0.9,6.};
  double cross_FF1_norad_m1[]={60.1, 54.5, 42.1, 29.7, 4.29};
  //TFF_Asym
  
  double cross_section_dQ2_TL(double *Q2, double *par);
  double cross_section_dQ2_TT(double *Q2, double *par);
  double factm0 = (0.86/1000. - 0.517/10000.)/(0.365/1000.-0.517/10000.);
  double factm1 = (0.161/100.-0.169/10000.)/(0.344/1000.-0.169/10000.);
  double en[100000],den[100000],cr0[100000],dcr0[100000],cr1[100000],dcr1[100000];
   int nrun = 0;
   double total[1000],dtotal[1000];
   ifstream stream("../eef1/results/cross.dat");
   double llll;
   double enO[1],denO[1],totalO[1],dtotalO[1];
   while(stream.eof()==0){
     stream >> en[nrun] >> den[nrun] >> cr0[nrun] >> dcr0[nrun] >> cr1[nrun] >> dcr1[nrun] >> llll >> llll >> llll >> llll;
     cr0[nrun] = cr0[nrun]/1000./(den[nrun]-en[nrun]);
     dcr0[nrun] = dcr0[nrun]/1000./(den[nrun]-en[nrun]);
     cr1[nrun] = cr1[nrun]/1000./(den[nrun]-en[nrun]);
     dcr1[nrun] = dcr1[nrun]/1000./(den[nrun]-en[nrun]);
     den[nrun] = (den[nrun] - en[nrun])/2.;
     en[nrun] = en[nrun] + den[nrun];
     total[nrun] = cr0[nrun] + cr1[nrun];
     dtotal[nrun] = sqrt(pow(dcr0[nrun],2.) + pow(dcr1[nrun],2.));
     if(stream.eof()==1)break;
     nrun++;
   }
   enO[0] = (6.+0.9)/2.;
   denO[0] = (6.-0.9)/2.;
   totalO[0] = (cr0[0] + cr0[1] +cr0[2])*factm0 + (cr1[0] + cr1[1] +cr1[2])*factm1;
   dtotalO[0] = sqrt(pow(dcr0[0]*factm0,2.)+pow(dcr0[1]*factm0,2.)+pow(dcr0[2]*factm0,2.)+pow(dcr1[0]*factm1,2.)+pow(dcr1[1]*factm1,2.)+pow(dcr1[2]*factm1,2.));
   cout << " TOTAL0 = " << totalO[0] << " + /- " << dtotalO[0] << endl;
   TCanvas *s = new TCanvas();
   TH1F *frd  = s->DrawFrame(-2.,0.001,22.,80.);
   if(mode == 1)frd  = s->DrawFrame(1.,0.001,22.,1.);
   frd->SetXTitle("Q^{2}, GeV^{2}");
   frd->SetYTitle("d#sigma/dQ^{2} (pb/GeV^{2})");

   TGraphErrors *CrossO  = new TGraphErrors(1,enO,totalO,denO,dtotalO);
   CrossO->SetMarkerColor(1);
   CrossO->SetMarkerStyle(27);
   CrossO->SetLineColor(1);
   CrossO->SetLineWidth(3.);
   CrossO->SetTitle("#sigma(0.9--6)");
   if(mode == 0)CrossO->Draw("P");
   
   TGraphErrors *Crossr  = new TGraphErrors(nrun,en,total,den,dtotal);
   Crossr->SetMarkerColor(1);
   Crossr->SetMarkerStyle(21);
   Crossr->SetLineColor(1);
   Crossr->SetLineWidth(2.);
   Crossr->SetTitle("#sigma(TT) + #sigma(ST)");
   Crossr->Draw("P");

   
   TGraphErrors *Cross  = new TGraphErrors(nrun,en,cr0,den,dcr0);
   Cross->SetMarkerColor(2);
   Cross->SetMarkerStyle(26);
   Cross->SetLineColor(2);
   Cross->SetLineWidth(2.);
   Cross->SetTitle("#sigma(TT)");
   Cross->Draw("P");

   TGraphErrors *Cross1  = new TGraphErrors(nrun,en,cr1,den,dcr1);
   Cross1->SetMarkerColor(4);
   Cross1->SetMarkerStyle(23);
   Cross1->SetLineColor(4);
   Cross1->SetLineWidth(2.);
   Cross1->SetTitle("#sigma(ST)");
   Cross1->Draw("P");

   double ener0[] = {0.02, 0.1, 0.4, 0.9};
   double ener1[] = {0.1, 0.4, 0.9, 6.0};
   double cross[] = {28.9, 57.7, 29.8, 28.2};
   double dcross[] = {8., 7.6, 4.7, 3.8};
   double ddcross[] = {2.7, 5.3, 3.8, 4.3};
   for(int i = 0; i < 4; i++){
     double btv = (ener1[i] - ener0[i])/2.;
     ener0[i] = (ener1[i] + ener0[i])/2.;
     ener1[i] = btv;
     cross[i] = cross[i]/5.4/2./btv;
     dcross[i] = dcross[i]/5.4/2./btv;
     ddcross[i] = ddcross[i]/5.4/2./btv;
     dcross[i] = sqrt(pow(dcross[i],2.) + pow(ddcross[i],2.));
   }

   
   TGraphErrors *Cross2  = new TGraphErrors(4,ener0,cross,ener1,dcross);
   Cross2->SetMarkerColor(1);
   Cross2->SetMarkerStyle(20);
   Cross2->SetLineColor(1);
   Cross2->SetLineWidth(2.);
   Cross2->SetTitle("L3");
   Cross2->SetFillColor(2);
   Cross2->SetFillStyle(3001);
   if(mode == 0)Cross2->Draw("2");

   if(mode == 0)Cross->Draw("P");
   if(mode == 0)Crossr->Draw("P");

   TF1 *u1 = new TF1("ST(L3 ff)",cross_section_dQ2_TL,0.001,20);
   u1->SetLineColor(4);
   u1->Draw("same");
   TF1 *u2 = new TF1("TT(L3 ff)",cross_section_dQ2_TT,0.001,20);
   u2->SetLineColor(2);
   u2->Draw("same");
   
   s->BuildLegend();

   
}

void draw_cross_sect(int mode = 0){
  double factm0 = (0.86/1000. - 0.517/10000.)/(0.365/1000.-0.517/10000.);
  double factm1 = (0.161/100.-0.169/10000.)/(0.344/1000.-0.169/10000.);
  double en[100000],den[100000],cr0[100000],dcr0[100000],cr1[100000],dcr1[100000];
   int nrun = 0;
   double total[1000],dtotal[1000];
   ifstream stream("../eef1/results/cross.dat");
   double llll;
   double enO[1],denO[1],totalO[1],dtotalO[1];
   while(stream.eof()==0){
     stream >> en[nrun] >> den[nrun] >> cr0[nrun] >> dcr0[nrun] >> cr1[nrun] >> dcr1[nrun] >> llll >> llll >> llll >> llll;
     cr0[nrun] = cr0[nrun]/1000.;
     dcr0[nrun] = dcr0[nrun]/1000.;
     cr1[nrun] = cr1[nrun]/1000.;
     dcr1[nrun] = dcr1[nrun]/1000.;
     double btv = (den[nrun] - en[nrun])/2.;
     en[nrun] = en[nrun]/2. + den[nrun]/2.;
     den[nrun] = btv;
     total[nrun] = cr0[nrun] + cr1[nrun];
     dtotal[nrun] = sqrt(pow(dcr0[nrun],2.) + pow(dcr1[nrun],2.));
     if(stream.eof()==1)break;
     nrun++;
   }
   enO[0] = (6.+0.9)/2.;
   denO[0] = (6.-0.9)/2.;
   totalO[0] = (cr0[0] + cr0[1] +cr0[2])*factm0 + (cr1[0] + cr1[1] +cr1[2])*factm1;
   dtotalO[0] = sqrt(pow(dcr0[0]*factm0,2.)+pow(dcr0[1]*factm0,2.)+pow(dcr0[2]*factm0,2.)+pow(dcr1[0]*factm1,2.)+pow(dcr1[1]*factm1,2.)+pow(dcr1[2]*factm1,2.));
   
   TCanvas *s = new TCanvas();
   TH1F *frd  = s->DrawFrame(-2.,0.001,22.,80.);
   if(mode == 1)frd  = s->DrawFrame(1.,0.001,22.,1.);
   frd->SetXTitle("Q^{2}, GeV^{2}");
   frd->SetYTitle("#Delta#sigma (pb)");


   double unc_cor_cr1[nrun];
   for(int i = 0; i < nrun; i++)unc_cor_cr1[i] = dcr2[i]*total[i];
   TGraphAsymmErrors *Crossdf  = new TGraphAsymmErrors(nrun,en,cr1,den,den,unc_cor_cr1,unc_cor_cr1);
   Crossdf->SetMarkerColor(4);
   Crossdf->SetMarkerStyle(23);
   Crossdf->SetLineColor(4);
   Crossdf->SetLineWidth(1.);
   Crossdf->SetTitle("corr uncer");
   Crossdf->SetMarkerColor(2);
   Crossdf->SetMarkerStyle(20);
   Crossdf->SetLineColor(2);
   Crossdf->SetLineWidth(2.);
   Crossdf->SetTitle("");
   Crossdf->SetFillColor(3);
   Crossdf->SetFillStyle(3001);
   Crossdf->Draw("a2");

   double unc_cor_cr0[nrun];
   for(int i = 0; i < nrun; i++)unc_cor_cr0[i] = dcr2[i]*cr0[i];
   TGraphAsymmErrors *Crossdf0  = new TGraphAsymmErrors(nrun,en,cr0,den,den,unc_cor_cr0,unc_cor_cr0);
   Crossdf0->SetMarkerColor(4);
   Crossdf0->SetMarkerStyle(23);
   Crossdf0->SetLineColor(4);
   Crossdf0->SetLineWidth(1.);
   Crossdf0->SetTitle("corr uncer");
   Crossdf0->SetMarkerColor(2);
   Crossdf0->SetMarkerStyle(20);
   Crossdf0->SetLineColor(2);
   Crossdf0->SetLineWidth(2.);
   Crossdf0->SetTitle("");
   Crossdf0->SetFillColor(3);
   Crossdf0->SetFillStyle(3001);
   //Crossdf0->Draw("a2 same");
   

   TGraphErrors *CrossO  = new TGraphErrors(1,enO,totalO,denO,dtotalO);
   CrossO->SetMarkerColor(1);
   CrossO->SetMarkerStyle(27);
   CrossO->SetMarkerSize(3);
   CrossO->SetLineColor(1);
   CrossO->SetLineWidth(3.);
   CrossO->SetTitle("#sigma(0.9--6)");
   if(mode == 0)CrossO->Draw("P");
   
   TGraphErrors *Crossr  = new TGraphErrors(nrun,en,total,den,dtotal);
   Crossr->SetMarkerColor(1);
   Crossr->SetMarkerStyle(21);
   Crossr->SetMarkerSize(1.3);
   Crossr->SetLineColor(1);
   Crossr->SetLineWidth(3.);
   Crossr->SetTitle("#sigma(TT) + #sigma(ST)");
   Crossr->Draw("P");

   
   TGraphErrors *Cross  = new TGraphErrors(nrun,en,cr0,den,dcr0);
   Cross->SetMarkerColor(2);
   Cross->SetMarkerStyle(26);
   Cross->SetLineColor(2);
   Cross->SetLineWidth(1.);
   Cross->SetTitle("#sigma(TT)");
   Cross->Draw("P");

   TGraphErrors *Cross1  = new TGraphErrors(nrun,en,cr1,den,dcr1);
   Cross1->SetMarkerColor(4);
   Cross1->SetMarkerStyle(23);
   Cross1->SetLineColor(4);
   Cross1->SetLineWidth(1.);
   Cross1->SetTitle("#sigma(TL)");
   Cross1->Draw("P");

   double ener0[] = {0.02, 0.1, 0.4, 0.9};
   double ener1[] = {0.1, 0.4, 0.9, 6.0};
   double cross[] = {28.9, 57.7, 29.8, 28.2};
   double dcross[] = {8., 7.6, 4.7, 3.8};
   double ddcross[] = {2.7, 5.3, 3.8, 4.3};
   for(int i = 0; i < 4; i++){
     double btv = (ener1[i] - ener0[i])/2.;
     ener0[i] = (ener1[i] + ener0[i])/2.;
     ener1[i] = btv;
     cross[i] = cross[i]/5.4;
     dcross[i] = dcross[i]/5.4;
     ddcross[i] = ddcross[i]/5.4;
     dcross[i] = sqrt(pow(dcross[i],2.) + pow(ddcross[i],2.));
   }

   
   TGraphErrors *Cross2  = new TGraphErrors(4,ener0,cross,ener1,dcross);
   Cross2->SetMarkerColor(1);
   Cross2->SetMarkerStyle(20);
   Cross2->SetLineColor(1);
   Cross2->SetLineWidth(2.);
   Cross2->SetTitle("L3");
   Cross2->SetFillColor(2);
   Cross2->SetFillStyle(3001);
   if(mode == 0)Cross2->Draw("2");

   if(mode == 0)Cross->Draw("P");
   if(mode == 0)Crossr->Draw("P");

   double crm1[10],crm0[10],dcrm1[10],enth[10],denrh[10],dcrm0[10];
   TF1 *fcr = new TF1("fcr",dcrs,0.01,30,2);
   fcr->SetParameter(2.,0.2);
   for(int i = 0; i < 4; i++){
     crm1[i] = fcr->Integral(ener0[i] - ener1[i], ener0[i] + ener1[i]);
     enth[i] = ener0[i];
     denrh[i] = ener1[i];
     dcrm0[i] = 0.;
   }

   for(int i = 4; i < 10; i++){

     crm1[i] = fcr->Integral(en[i-4] - den[i-4], en[i-4] + den[i-4]);
     enth[i] = en[i-4];
     denrh[i] = den[i-4];     dcrm0[i] = 0.;
   }   
   
   TGraphErrors *Cross4  = new TGraphErrors(10,enth,crm1,denrh,dcrm0);
   Cross4->SetMarkerColor(1);
   Cross4->SetMarkerStyle(20);
   Cross4->SetLineColor(1);
   Cross4->SetLineWidth(2.);
   Cross4->SetTitle("MC SIM m0");
   Cross4->SetFillColor(1);
   Cross4->SetFillStyle(1001);
   //   if(mode == 0)Cross4->Draw("P");

   

   s->BuildLegend();
   
}


double TFF_Fpr(double *e, double *par){

  double omega2 = 1.756*1.756;
  double M2 = 1.271*1.271;
  double Mrho2 = 0.77*0.77;
  return par[0]/(1.+e[0]/Mrho2);
  return par[0]*4.*M2*(omega2+e[0])/2./Mrho2/(e[0]+Mrho2);
    
}

double TFF_Mil(double *e, double *par){

  double M2 = 1.271*1.271*1.271*1.271;
  double Mrho2 = 0.77*0.77;
  double g2 = 2.5/10000.;
  double q = 1./1.271*0.5*(1.271*1.271+e[0]);
  double res = 2.*g2*M2/(e[0]+Mrho2)/Mrho2;
  return res;
    
}

double TFF_Asym(double *e, double *par){

  double F0 = 0.245/sqrt(2);
  double F3 = 0.238/sqrt(2);
  double F8 = 0.239/sqrt(2);
  double C0 = 2./3./sqrt(6);
  double C3 = 1./6.;
  double C8 = 1./6./sqrt(3);
  //  double res = 3.*(C0*F0 + C3*F3 + C8*F8)*pow(1.285,3.);
  double res = 3.*0.146*pow(1.285,3.); 
  res = res/pow(1.285,2.)*2.*(e[0] + 1.285*1.285)/2.;
  res = res/e[0]/e[0];
  double X = -pow(1.285,2.)/e[0];
  res = res*2./X/X*(X/(1.-X) + log(1.-X));
  return res;
    
}

double TFF_Asym_err(double *e, double *par){

  double F0 = 0.245/sqrt(2);
  double F3 = 0.238/sqrt(2);
  double F8 = 0.239/sqrt(2);
  double C0 = 2./3./sqrt(6);
  double C3 = 1./6.;
  double C8 = 1./6./sqrt(3);
  //  double res = 3.*(C0*F0 + C3*F3 + C8*F8)*pow(1.285,3.);
  double res = 3.*0.014*pow(1.285,3.); 
  res = res/pow(1.285,2.)*2.*(e[0] + 1.285*1.285)/2.;
  res = res/e[0]/e[0];
  double X = -pow(1.285,2.)/e[0];
  res = res*2./X/X*(X/(1.-X) + log(1.-X));
  return res;
}


double TFF_Asym_interp(double *e, double *par){

  double res0 = sqrt(1.4/0.52/1000000./1.285/3.14*137.*137.*12.);
  double mf1_2 = pow(1.285,2.);
  double Qp = e[0] + mf1_2;
  double exp1 = exp(-0.001*e[0]/1.285/1.285);
  double res = res0/TFF_Asym(&mf1_2,par)*TFF_Asym(&Qp,par)*(TFF_Asym(&mf1_2,par)/res0+exp1)/(1+exp1);
  return res;
    
}


double F_L3_squared(double *Q2, double *par){

  double M = 1.281;
  double L = 1.04;
  double res = Q2[0]/M/M*(1.+0.5*Q2[0]/M/M)*2/pow(1+Q2[0]/L/L,4.);

  
  return res;
  

}


double cross_section_dQ2_TT_(double *Q2, double *par){

  // PRL 59, 18 (1987)
  double alpha = 1/137.;
  double hc = 3.88*pow(10.,8.); // pb/GeV^2  197 MeV*fm
  double G = 3.5/1000/1000;
  double M = 1.281;
  double Qmax_me = 100/511./511.*pow(10.,12);
  double s = 10.8*10.8; //par[0]*par[0];
  double dzeta = (M*M + Q2[0])/s;
  double res = alpha*alpha*12.*G/pow(M,3.)*(log(Qmax_me)*(log(1/dzeta)-3./2.) + pow(log(1/dzeta),2.) - 2.5*log(1/dzeta)-3.14*3.14/6.+19./8.);
  res = hc*res*Q2[0]/pow(M,4.)*F_L3_squared(Q2,par);
  
  return res;
}


double cross_section_dQ2_TL(double *Q2, double *par){

  // PRL 59, 18 (1987)
  cout << " Hello " << endl;
  double alpha = 1/137.;

  double hc = 3.88*pow(10.,8.); // pb/GeV^2  197 MeV*fm

  double G = 3.5/1000/1000;
  double M = 1.285;
  cout << " Hello1 " << endl;
  double Qmax_me = 100/511./511.*pow(10.,12);
  cout << " Hello1 " << endl;
  double s = 10.8*10.8; //par[0]*par[0];
  cout << " Hello2 " << endl;
  double dzeta = (M*M + Q2[0])/s;
  cout << "dzeta = " << dzeta << " ";
  double res = alpha*alpha*24.*G/pow(M,3.)*(log(Qmax_me)*(log(1/dzeta)-7./4.) + pow(log(1/dzeta),2.) - 3.*log(1/dzeta)-3.14*3.14/6.+23./8.);
  cout << "res = " << res << " ";
  res = hc*res/M/M*F_L3_squared(Q2,par);
  
  return res;
}

double cross_section_dQ2(double *Q2, double *par){

  // PRD 35, 11 (1987)
  double alpha_ = 1/137.;
  double pi_ = TMath::Pi();
  double hc = 3.88*pow(10.,8.); // pb/GeV^2  197 MeV*fm
  double s = 10.8*10.8;
  double M = 1.285;
  double G = 3.5/1000/1000;
  double me_ = 0.511/1000.;
  double res = hc*pow(alpha_/2./pi_,2.)*log(s/4./me_/me_)*pi_/4.;
  res*=96.*pi_/pow(M,5.)*G*F_L3_squared(Q2,par);

  return res;
}

double cross_section_dQ2_TT(double *Q2, double *par){

  double M = 1.285;
  double s = 10.8*10.8;
  double tau_ = (M*M+Q2[0])/s;
  double res = cross_section_dQ2(Q2, par);
  res*=4.*(1+tau_)*log(1/tau_)-(1-tau_)*(7+tau_);

  return res;
}
  


double  effic_Q2(double left, double right, int m, bool FF){

  TTree *tree_effic;
  TTree *tree_effic_gen;
  if(m==0){
    tree_effic = tree1;
    tree_effic_gen = tree_m0_gen;
  }
  else if(m==1){
    tree_effic = tree2;
    tree_effic_gen = tree_m1_gen;
  }
  else {std::cerr << "File not found" << std::endl; return 0.;}
  
  double res = 0;
  double res0 = 0;
  int Nsteps=10;
  double step = (right-left)/(double)Nsteps;
  double factor;
  for(int i = 0; i < Nsteps; i++){
    double q2 = left + step*i + step/2.;
    double Neff = tree_effic->GetEntries(Form("mf1 > 1.17 && mf1 < 1.4 && Q2 > %g && Q2 < %g",q2-step/2.,q2+step/2.));
    factor = 1;
    if(FF)factor = pow(TFF_Asym_interp(&q2, &q2),2.);
    res += factor*Neff;
    res0 += factor*(double)tree_effic_gen->GetEntries(Form("Q2gen > %g && Q2gen < %g",q2-step/2.,q2+step/2.));
  }
  return res/res0;
  
}



void draw_ff(){

  double en[100000],den[100000],cr0[100000],dcr0[100000],cr1[100000],dcr1[100000];
   int nrun = 0;
   double total[1000],dtotal[1000];
   ifstream stream("../eef1/results/cross.dat");
   double llll;
   while(stream.eof()==0){
     stream >> en[nrun] >> den[nrun] >> llll >> llll >> llll >> llll >> cr0[nrun] >> dcr0[nrun] >> cr1[nrun] >> dcr1[nrun];
     if(stream.eof()==1)break;
     double btv = (den[nrun] - en[nrun])/2.;
     en[nrun] = en[nrun]/2. + den[nrun]/2.;
     den[nrun] = btv;
     //dcr1[nrun] = sqrt(/*pow(dcr0[nrun]/cr0[nrun],2.) + */pow(dcr1[nrun]/cr1[nrun],2.));
     //cr1[nrun] = cr1[nrun];///cr0[nrun];
     //dcr1[nrun] =  dcr1[nrun]*cr1[nrun];
     nrun++;
   }


   TCanvas *s = new TCanvas();
   TH1F *frd  = s->DrawFrame(-2.,1.,22.,1.);
   frd->SetXTitle("Q^{2}, GeV^{2}");
   frd->SetYTitle("form factor");

   double ax[] = {0};
   double ay[] = {0.39};
   double aex[] = {0.2};
   double aey[] = {0.04};
   TGraphErrors* gae = new TGraphErrors(1, ax, ay, aex, aey);
   gae->SetFillColor(1);
   gae->SetFillStyle(1001);
   gae->SetTitle("G(0,0)");
   gae->Draw("2");
   
   TGraphErrors *Crossr  = new TGraphErrors(nrun,en,cr0,den,dcr0);
   Crossr->SetMarkerColor(2);
   Crossr->SetMarkerStyle(20);
   Crossr->SetLineColor(2);
   Crossr->SetLineWidth(2.);
   Crossr->SetTitle("|G'-F'|");
   
   
   TGraphErrors *Crossr1  = new TGraphErrors(nrun,en,cr1,den,dcr1);
   Crossr1->SetMarkerColor(4);
   Crossr1->SetMarkerStyle(20);
   Crossr1->SetLineColor(4);
   Crossr1->SetLineWidth(2.);
   Crossr1->SetTitle("|G'|");

   TF1 *u_fpr = new TF1("Asymp",TFF_Asym,1.,40,0);
   u_fpr->SetLineColor(1);
   //   u_fpr->Draw("same");


   double er[100],err[100],der[100],derr[100];
   for(int i = 0; i < 100; i++){

     er[i] = 0.5 + i/100.*30.;
     err[i] = TFF_Asym(&er[i],&er[i]);
     der[i] = 0.1;
     derr[i] = TFF_Asym_err(&er[i],&er[i]);

   }  
   TGraphErrors* ge = new TGraphErrors(100, er, err, der, derr);
   ge->SetFillColor(3);
   ge->SetFillStyle(3001);
   ge->Draw("3");
   ge->SetLineColor(3);
   ge->SetTitle("Asymptotic");
   Crossr->Draw("P");
   Crossr1->Draw("P");
}




