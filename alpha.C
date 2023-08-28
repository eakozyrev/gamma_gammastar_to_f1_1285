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
#include "TTree.h"
#include "THStack.h"
#include "RooPolynomial.h"
#include "TEfficiency.h"

#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>

#define chi25Ccut 25.
#include <stdlib.h>

using namespace std;


#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooConstVar.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#include "RooKeysPdf.h"
#include "RooFFTConvPdf.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
using namespace RooFit;


double paramm[100];

inline bool exists_file(const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
}


double cross_total_m0 = 33.8*2.;//  pb
double cross_total_m1 = 17.4*2.;//  pb


TFile *newfile0, *newfile1, *newfile2, *newfile3, *newfile4;
TTree *tree0, *tree1, *tree2, *tree3, *tree4, *tree;


void open(){

  newfile0 = TFile::Open("../eef1/histograms/data.root");//data.root");
  tree0 = (TTree*)newfile0->Get("Tree");

  newfile1 = TFile::Open("../eef1/histograms/mc_m0.root");
  tree1 = (TTree*)newfile1->Get("Treegen");

  newfile2 = TFile::Open("../eef1/histograms/mc_m1.root");
  tree2 = (TTree*)newfile2->Get("Treegen");

  newfile3 = TFile::Open("../eef1/histograms/mc_f0_eta_m0.root");
  tree3 = (TTree*)newfile3->Get("Tree");

  newfile4 = TFile::Open("../eef1/histograms/mc_f0_eta_m1.root");
  tree4 = (TTree*)newfile4->Get("Tree");


}




double fit2g1(string filee, string filemc, double qmin, double qmax){

  // Observable
  double startt = 1.17;//startt, endd
  double endd = 1.4;
  RooRealVar x("mf1","m_{#eta2#pi}, MeV/c^{2}",startt, endd);
  RooRealVar frac("frac","frac",0.9,0.040,1.0);

  TFile *newfile = TFile::Open((""+filee).c_str());
  TTree* tree = (TTree*)newfile->Get("Tree");
  
   x.setBins(50); 
   //   RooRealVar cosa0("cosa0","cosa0",-4,4); 
   RooRealVar Q2("Q2","Q2",0,400); 
   string selection = Form("mf1 > %g && mf1 < %g && Q2 > %g && Q2 < %g",startt,endd, qmin, qmax);
   RooDataSet *data = new RooDataSet("data","data",RooArgSet(x,Q2),Import(*tree),Cut(selection.c_str()));


  //*************************************************************************
  //===========================MC============================================

  TFile *newfilemc = TFile::Open((""+filemc).c_str());
  TTree* treemc = (TTree*)newfilemc->Get("Tree");

  string selectionmc = Form("mf1 > %g && mf1 < %g && Q2 > %g-5 && Q2 < %g+2",startt,endd, qmin, qmax); 
  RooDataSet datamc("datamc","datamc",RooArgSet(x,Q2),Import(*treemc),Cut(selectionmc.c_str()));
  RooKeysPdf kest1("kest1","kest1",x,datamc,RooKeysPdf::MirrorBoth) ;
  RooRealVar mg("mg","mg",0.);//,-3,3); 
  RooRealVar sg("sg","sg",0.0003);//,0.0000001,0.1); 
  RooGaussian gauss("gauss","gauss",x,mg,sg);
  cout << "===========================================================================9" << endl;
  RooFFTConvPdf resol("lxg","landau (X) gauss",x,kest1,gauss);
  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++9" << endl;
  //**************************************************************************
  //========================== BCKGR ===========================================
  /*
  TFile *newfileBCKGR = TFile::Open(("../4pisel/histograms/"+filebkgr).c_str());
  TTree* treeBCKGR = (TTree*)newfileBCKGR->Get("Tree");
  RooDataSet dataBCKGR("dataBCKGR","dataBCKGR",RooArgSet(RooArgSet(x,helicity1,helicity2,chi25C,tthg1, tthg2, tthg3, tthg4),RooArgSet(tthpic1,tthpic2,momg1,momg2,momg3,momg4)),Import(*treeBCKGR),Cut(selection.c_str()));
  RooKeysPdf kestBCKGR("kestBCKGR","kestBCKGR",x,dataBCKGR,RooKeysPdf::MirrorBoth);
  */
  //**************************************************************************
  //==========================MODEL===========================================
  RooRealVar c0("c0","c0",1.,0.,30.); 
  RooRealVar c1("c1","c1",0.,-10.,10.); 
  RooRealVar c2("c2","c2",0.,-1,1.); 
  RooGenericPdf P("P","c0 + c1*(mf1-1.28) + c2*(mf1-1.28)*(mf1-1.28)",RooArgSet(x,c0,c1,c2)); 
  RooPolynomial pol0("pol0","pol0",x,RooArgList());
  
  RooAddPdf model("model","model",RooArgList(resol,P),frac) ;
  //model.plotOn(xframe1);
  // Construct unbinned likelihood of model w.r.t. data
  RooAbsReal* nll = model.createNLL(*data) ;

  // I n t e r a c t i v e   m i n i m i z a t i o n ,   e r r o r   a n a l y s i s
  // -------------------------------------------------------------------------------

  // Create MINUIT interface object
  //RooMinuit m(*nll) ;

  // Activate verbose logging of MINUIT parameter space stepping
  //m.setVerbose(kTRUE) ;
  //m.setVerbose(kFALSE) ;

  // Call MIGRAD to minimize the likelihood
  RooFitResult* r = model.fitTo(*data,Save());
  //m.migrad();
  cout << "-log(L) at minimum = " << r->minNll() << endl ;
  cout << " (double)frac.getVal() = " <<  (double)frac.getVal() << endl;
  
  paramm[0] = (double)frac.getVal()*tree->GetEntries(TCut(selection.c_str()));
  paramm[1] = frac.getError()*tree->GetEntries(TCut(selection.c_str()));

  paramm[2] = (double)mg.getVal();
  paramm[3] = mg.getError();

  paramm[4] = (double)c2.getVal();
  paramm[5] = c2.getError();
  
  // Print values of all parameters, that reflect values (and error estimates)
  // that are back propagated from MINUIT
  //model.getParameters(x)->Print("s") ;
  RooPlot* xframe1 = x.frame(Title(Form("%g < Q^{2} < %g GeV^{2}",qmin,qmax)));
  data->plotOn(xframe1);
  model.plotOn(xframe1); 
  model.plotOn(xframe1,Components(P),LineStyle(kDashed));
  
  TCanvas *s = new TCanvas();
  TPad*    upperPad = new TPad("upperPad", "upperPad", 0.0,0.05,1.0,1.0);
  TPad*    lowerPad = new TPad("lowerPad", "lowerPad", 0.0,0.1,1.0,0.3);
  upperPad->Draw();
  //lowerPad->Draw();
  upperPad->cd();
  //gPad->SetLeftMargin(0.15) ; 
  xframe1->GetYaxis()->SetTitleOffset(0.7);
  xframe1->SetNdivisions(8,"Y");
  xframe1->GetYaxis()->SetLabelOffset(0.007); 
  xframe1->Draw();
  
  /*
  TH1F *hmc = new TH1F("hmc","dfssddf",myh->GetNbinsX(),myh->GetBinLowEdge(1),myh->GetBinCenter(myh->GetNbinsX())+myh->GetBinWidth(1)/2.);
  model.fillHistogram(hmc,RooArgList(x),tree->GetEntries(cutt));
  myh->Divide(hmc);
  myh->SetYTitle("data/fit");
  myh->SetTitle("");
  myh->SetAxisRange(0.5,1.5,"Y");
  myh->SetAxisRange(startt,endd,"X");
  myh->SetXTitle("m_{#pi^{0}}, MeV/c^{2}");
  
  
  lowerPad->cd();
  gStyle->SetOptStat(0);
  myh->SetLineWidth(3.);
  myh->SetLineColor(1);
  gStyle->SetTextSize(1.8);
  myh->SetLabelSize(0.19,"xy");
  myh->SetTitleSize(0.19,"xy");
  myh->SetNdivisions(2,"Y");
  myh->SetTitleOffset(1.2,"x");
  myh->SetTitleOffset(0.2,"y");
  lowerPad->SetLeftMargin(0.15) ;	
  myh->Draw();
  s->SaveAs(("plots/"+filee+".png").c_str());
  s->SaveAs(("plots/"+filee+".root").c_str());
  s->Close(); 
  */
  //  newfile->Close();
  //  newfilemc->Close();
  cout << "===============================================" << endl;
  cout << "                     end                       " << endl;
  cout << "===============================================" << endl;
  s->SaveAs(Form("../eef1/figs/%g_%g.png",qmin,qmax));
  return 1.;

}


void fitall(){

  int nrun = 6;
  double en0[] = {2,4,5,6,7,10,20.};
  double en[20],den[20];
  double cr[20],dcr[20];
  ofstream stream("../eef1/results/nevents.dat");
  
  for(int i = 0; i < nrun; i++){
    en[i] = (en0[i]+en0[i+1])/2.;
    den[i] = (en0[i+1] - en0[i])/2.;
    stream << en0[i] << " " << en0[i+1] << " "; 
    fit2g1("../eef1/histograms/data.root","../eef1/histograms/mc_m0.root",en0[i],en0[i+1]);
    stream << paramm[0] << " " << paramm[1] << endl;
    cr[i] = paramm[4];
    dcr[i] = paramm[5];
  }

  TCanvas *s0 = new TCanvas();
  TH1F *frd  = s0->DrawFrame(-1,0.,20,cr[1]);
  frd->SetXTitle("Q^{2} (GeV^{2})");
  frd->SetYTitle("#sigma_{E}/E");
  TGraphErrors *Crossu  = new TGraphErrors(nrun,en,cr,den,dcr);
  Crossu->SetMarkerColor(2);
  Crossu->SetMarkerStyle(20);
  Crossu->SetLineColor(2);
  Crossu->SetLineWidth(2.);
  Crossu->SetTitle("");
  Crossu->Draw("P");


}



void effic_draw(string file, string same){
  TFile *newfilem0 = TFile::Open(file.c_str());
  TTree* treem0 = (TTree*)newfilem0->Get("Tree");
  TH1D *h_Q2 = new TH1D("h_Q2","h_Q2",40,1,20);
  h_Q2->SetXTitle("Q^{2} (GeV)^{2}");
  TCanvas s;
  treem0->Draw("Q2 >> h_Q2");
  TTree* treem0gen = (TTree*)newfilem0->Get("Treegen");
  TH1D *h_Q2gen = new TH1D("h_Q2gen","h_Q2gen",40,1,20);
  treem0gen->Draw("Q2gen >> h_Q2gen");
  s.Close();
  TEfficiency *hEfficiency = new TEfficiency(*h_Q2,*h_Q2gen);
  hEfficiency->SetTitle(" ; Q^{2} (GeV^{2}) ; #varepsilon");
  hEfficiency->Draw(same.c_str());
  //  h_Q2->Divide(h_Q2gen);
  //  h_Q2->Draw();
}

double form_factor_1(double Q2){

  return 1./(1.+Q2/0.775/0.775);

}

/*
void cross_section(){

  open();
  
  ifstream stream("results/nevents.dat");
  ifstream stream_ratio("results/ratio.dat");  
  ofstream streamof("results/cross.dat");
  ofstream streamof0("results/cross0.dat");
  ofstream streamof1("results/cross1.dat");
  ofstream stream_cr_mc("results/cross_mc.dat");

  ofstream stream_table_cr0("results/table_cr0.dat");
  ofstream stream_table_cr1("results/table_cr1.dat");
  
  double crossmc_m0, crossmc_m1;


  double factorm0 = 0.54;double dfactorm0;
  double btv;
  while(stream.eof()==0){
    double ena,enb,Nev,dNev;
    stream >> ena >> enb >> Nev >> dNev;
    stream_ratio >> btv >> btv >>  factorm0 >> dfactorm0;
    if(stream.eof()==1)break;
    double Nm0 = factorm0*Nev;
    double dNm0 = sqrt(pow(factorm0*dNev,2.) + pow(dfactorm0*Nev,2.));
    double Nm1 = (1.-factorm0)*Nev;
    double dNm1 = sqrt(pow((1.-factorm0)*dNev,2.) + pow(dfactorm0*Nev,2.));
    TTree* treem0 = (TTree*)newfile1->Get("Tree");
    TTree* treem0gen = (TTree*)newfile1->Get("Treegen");
    double partm0 = treem0gen->GetEntries(Form("Q2gen > %g && Q2gen < %g",ena,enb))/(double)treem0gen->GetEntries();
    double dpartm0 = sqrt(partm0*(1. - partm0)/(double)treem0gen->GetEntries());
    TTree* treem1 = (TTree*)newfile2->Get("Tree");
    TTree* treem1gen = (TTree*)newfile2->Get("Treegen");
    double partm1 = treem1gen->GetEntries(Form("Q2gen > %g & Q2gen < %g",ena,enb))/(double)treem1gen->GetEntries();
    double dpartm1 = sqrt(partm1*(1. - partm1)/(double)treem1gen->GetEntries());
    double cross_mc_m0 = partm0*cross_total_m0/(enb-ena)*1000.;
    double dcross_mc_m0 = dpartm0/partm0*cross_mc_m0;
    double cross_mc_m1 = partm1*cross_total_m1/(enb-ena)*1000.;
    double dcross_mc_m1 = dpartm1/partm1*cross_mc_m1;
    stream_cr_mc << ena << " " << enb << " " << cross_mc_m0 << " " << dcross_mc_m0 << " " << cross_mc_m1 << " " << dcross_mc_m1 << endl;
      
    double eff_m0 = treem0->GetEntries(Form("Q2 > %g && Q2 < %g",ena,enb))/(double)treem0gen->GetEntries(Form("Q2gen > %g && Q2gen < %g",ena,enb));
    double deff_m0 = sqrt(eff_m0*(1. - eff_m0)/(double)treem0gen->GetEntries(Form("Q2gen > %g && Q2gen < %g",ena,enb)));
    double eff_m1 = treem1->GetEntries(Form("Q2 > %g && Q2 < %g",ena,enb))/(double)treem1gen->GetEntries(Form("Q2gen > %g && Q2gen < %g",ena,enb));    
    double deff_m1 = sqrt(eff_m1*(1.-eff_m1)/(double)treem1gen->GetEntries(Form("Q2gen > %g && Q2gen < %g",ena,enb)));
    
    double cross_m0 = Nm0/eff_m0/469./0.348/(enb-ena);
    double dcross_m0 = dNm0/eff_m0/469./0.348/(enb-ena);
    double cross_m1 = Nm1/eff_m1/469./0.348/(enb-ena);
    double dcross_m1 = dNm1/eff_m1/469./0.348/(enb-ena);
    
    double ff_m0 = sqrt(cross_m0/cross_mc_m0);//form_factor_1((ena+enb)/2.)
    double dff_m0 = sqrt((cross_m0 +dcross_m0)/cross_mc_m0) - ff_m0;
    double ff_m1 = sqrt(cross_m1/cross_mc_m1);
    double dff_m1 = sqrt((cross_m1 + dcross_m1)/cross_mc_m1) - ff_m1;
    
    cout << ena << "-" << enb << " & " << ((int)(Nm0*10.))/10. << " $\\pm$ " << ((int)(dNm0*10.))/10. << " & " << ((int)(eff_m0*100000.))/1000. << " $\\pm$ " << ((int)(deff_m0*100000.))/1000. << " & " << (int)cross_mc_m0 << " $\\pm$ " << (int)dcross_mc_m0 <<  " & " << ((int)(cross_m0*100.))/100 << " $\\pm$ " << ((int)(dcross_m0*100.))/100 << " & " <<  Nm1 << " $\\pm$ " << dNm1 << " & " << ((int)(eff_m1*100000.))/1000. << " $\\pm$ " << ((int)(deff_m1*100000.))/1000.  << " & " << (int)cross_mc_m1 << " $\\pm$ " << (int)dcross_mc_m1 <<  " & " << ((int)(cross_m1*100.))/100 << " $\\pm$ " << ((int)(dcross_m1*100.))/100 << " \\\\" << endl;

    streamof << ena << " " << enb << " " << cross_m0 << " " << dcross_m0 << " "  << cross_m1 << " " << dcross_m1 << " " << ff_m0 << " " << dff_m0 << " "  << ff_m1 << "  " << dff_m1 << endl;

    streamof0 << ena << " " << enb << " " << Nm0 << " " << dNm0 << " " << eff_m0 << " " << deff_m0 << " " <<  cross_m0 << " " << dcross_m0 << " "  << ff_m0 << " " << dff_m0 << endl;

    streamof1 << ena << " " << enb << " " << Nm1 << " " << dNm1 << " " << eff_m1 << " " << deff_m1 << " " <<  cross_m1 << " " << dcross_m1 << " "  << ff_m1 << " " << dff_m1 << endl;

    stream_table_cr0 <<  ena << "$\\div$" << enb << " & " << ((int)(Nm0*10.))/10. << " $\\pm$ " << ((int)(dNm0*10.))/10. << " & " << ((int)(eff_m0*100000.))/1000. << " $\\pm$ " << ((int)(deff_m0*100000.))/1000. << " & " << (int)cross_mc_m0 << " $\\pm$ " << (int)dcross_mc_m0 <<  " & " << ((int)(cross_m0*100.))/100 << " $\\pm$ " << ((int)(dcross_m0*100.))/100 << " \\\\" << endl;
    

    stream_table_cr1 << ena << "$\\div$" << enb << " & " <<  Nm1 << " $\\pm$ " << dNm1 << " & " << ((int)(eff_m1*100000.))/1000. << " $\\pm$ " << ((int)(deff_m1*100000.))/1000.  << " & " << (int)cross_mc_m1 << " $\\pm$ " << (int)dcross_mc_m1 <<  " & " << ((int)(cross_m1*100.))/100 << " $\\pm$ " << ((int)(dcross_m1*100.))/100 << " \\\\" << endl;
    
  }
  
}


double fit_all(){

  fit2g1("histograms/data.root","histograms/mc_m0.root",0,100);
  cout << paramm[0] << " " << paramm[1] << endl;

  TFile *newfilem0 = TFile::Open("mc_m0.root");
  TTree* treem0 = (TTree*)newfilem0->Get("Tree");
  TTree* treem0gen = (TTree*)newfilem0->Get("Treegen");
  double effm0 = treem0->GetEntries()/(double)treem0gen->GetEntries();
  double deffm0 = sqrt(effm0*(1.-effm0)/(double)treem0gen->GetEntries());
  TFile *newfilem1 = TFile::Open("mc_m1.root");
  TTree* treem1 = (TTree*)newfilem1->Get("Tree");
  TTree* treem1gen = (TTree*)newfilem1->Get("Treegen");
  double effm1 = treem1->GetEntries()/(double)treem1gen->GetEntries();
  double deffm1 = sqrt(effm1*(1.-effm1)/(double)treem1gen->GetEntries());
  
  cout << effm0 << " +/- " << deffm0 << "     "  << effm1 << " +/- " << deffm1 << endl;
  

  return 1;
}
*/



void draw_cross_sect(){

  double en[100000],den[100000],cr0[100000],dcr0[100000],cr1[100000],dcr1[100000];
   int nrun = 0;
   double total[1000],dtotal[1000];
   ifstream stream("results/cross.dat");
   double llll;
   while(stream.eof()==0){
     stream >> en[nrun] >> den[nrun] >> cr0[nrun] >> dcr0[nrun] >> cr1[nrun] >> dcr1[nrun] >> llll >> llll >> llll >> llll;;
     double btv = (den[nrun] - en[nrun])/2.;
     en[nrun] = en[nrun]/2. + den[nrun]/2.;
     den[nrun] = btv;
     total[nrun] = cr0[nrun] + cr1[nrun];
     dtotal[nrun] = sqrt(pow(dcr0[nrun],2.) + pow(dcr1[nrun],2.));
     if(stream.eof()==1)break;
     nrun++;
   }
   TCanvas *s = new TCanvas();
   TH1F *frd  = s->DrawFrame(-2.,1.,22.,900000.);
   frd->SetXTitle("Q^{2}, GeV^{2}");
   frd->SetYTitle("d#sigma/dQ^{2} (fb/GeV^{2})");

   TGraphErrors *Crossr  = new TGraphErrors(nrun,en,total,den,dtotal);
   Crossr->SetMarkerColor(3);
   Crossr->SetMarkerStyle(20);
   Crossr->SetLineColor(3);
   Crossr->SetLineWidth(2.);
   Crossr->SetTitle("#sigma(TT) + #sigma(ST)");
   Crossr->Draw("P");

   
   TGraphErrors *Cross  = new TGraphErrors(nrun,en,cr0,den,dcr0);
   Cross->SetMarkerColor(2);
   Cross->SetMarkerStyle(20);
   Cross->SetLineColor(2);
   Cross->SetLineWidth(2.);
   Cross->SetTitle("#sigma(TT)");
   Cross->Draw("P");

   TGraphErrors *Cross1  = new TGraphErrors(nrun,en,cr1,den,dcr1);
   Cross1->SetMarkerColor(4);
   Cross1->SetMarkerStyle(20);
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
     cross[i] = cross[i]*1000./2./btv;
     dcross[i] = dcross[i]*1000./2./btv;
     ddcross[i] = ddcross[i]*1000./2./btv;
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
   Cross2->Draw("2");
   //   Cross2->Draw("p");   

   Cross->Draw("P");
   Crossr->Draw("P");
}


double TFF_Fpr(double *e, double *par){

  double omega2 = 1.756*1.756;
  double M2 = 1.271*1.271;
  double Mrho2 = 0.77*0.77;
  return par[0]/(1.+e[0]/Mrho2);
  return par[0]*4.*M2*(omega2+e[0])/2./Mrho2/(e[0]+Mrho2);
    
}

double TFF_Milsh(double *e, double *par){
  
  double M = 1.285;
  double M5 = pow(M,5.);
  double g2 = 2.5/10000.;
  double Mrho2 = 0.77*0.77;
  double q = (1.271*1.271 + e[0])/2./M;
  double res = g2*M5/q/(e[0]+Mrho2)/Mrho2;
  double nu = (e[0] + 1.281*1.281)/2.;
  return res/M/M*2*nu*137.;
    
}


double TFF_Milsh_F(double *e, double *par){

  double g1 = 2.9/10000.;
  double M = 1.285;
  double nu = (e[0] + 1.281*1.281)/2.;
  double q = nu/M;
  double res = g1*pow(M,3.)*e[0]/q/(e[0]+0.778*0.778)/0.778/0.778;
  res = res/e[0]/M/M*4*nu*nu*137.;
  
  return res;
					   }
//TF1 *u = new TF1("u",TFF_Milsh_F,1,30,0)
		 
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


double TFF_Asym_interp(double *e, double *par){

  double res0 = sqrt(1.4/0.52/1000000./1.285/3.14*137.*137.*12.);
  double mf1_2 = pow(1.285,2.);
  double Qp = e[0] + mf1_2;
  double res = res0/TFF_Asym(&mf1_2,par)*TFF_Asym(&Qp,par);
  return res;
    
}




// TF1 *u = new TF1("u",TFF_Asym,1,30,0)

double TFF_QM(double *e, double *par){

  double M = 1.285;
  double res = pow(M*M/(M*M+e[0]),2.)/pow(1.285,2.)*2.*1.4;
  return res;
    
}



void draw_ff(){

  double en[100000],den[100000],cr0[100000],dcr0[100000],cr1[100000],dcr1[100000];
   int nrun = 0;
   double total[1000],dtotal[1000];
   ifstream stream("results/cross.dat");
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
     
     TGraphErrors *Crossr  = new TGraphErrors(nrun,en,cr0,den,dcr0);
     Crossr->SetMarkerColor(2);
     Crossr->SetMarkerStyle(20);
     Crossr->SetLineColor(2);
     Crossr->SetLineWidth(2.);
     Crossr->SetTitle("|G'-F'|");
     Crossr->Draw("P");

     TGraphErrors *Crossr1  = new TGraphErrors(nrun,en,cr1,den,dcr1);
     Crossr1->SetMarkerColor(4);
     Crossr1->SetMarkerStyle(20);
     Crossr1->SetLineColor(4);
     Crossr1->SetLineWidth(2.);
     Crossr1->SetTitle("|G'|");
     Crossr1->Draw("P");

     TF1 *u_fpr = new TF1("Asymp",TFF_Asym,0.,40,0);
     u_fpr->Draw("same");

     TF1 *u_fpr1 = new TF1("intermediate",TFF_Asym_interp,0.,40,0);
     u_fpr1->SetLineColor(4);
     u_fpr1->Draw("same");


   
}


