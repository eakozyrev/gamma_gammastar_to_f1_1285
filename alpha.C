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
#include <iomanip>
//#include "RooChi2Var.h"
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
//#include "RooMinuit.h"
#include "RooFitResult.h"
#include "RooKeysPdf.h"
#include "RooFFTConvPdf.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"

//#include "/home/eakozyrev/Babar/fftw-3.3.10/api/fftw3.h"
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
  RooRealVar m2pi("m2pi","m2pi, MeV/c^{2}",0,1);
  RooRealVar frac("frac","frac",0.9,0.040,1.0);

  TFile *newfile = TFile::Open((""+filee).c_str());
  TTree* tree = (TTree*)newfile->Get("Tree");
  
   x.setBins(50);
   RooRealVar Q2("Q2","Q2",0,400); 
   //m2pi < 0.55 && 
   string selection = Form("mf1 > %g && mf1 < %g && Q2 > %g && Q2 < %g",startt,endd, qmin, qmax);
   RooDataSet *Data = new RooDataSet("data","data",RooArgSet(x,Q2,m2pi),Import(*tree),Cut(selection.c_str()));


  //*************************************************************************
  //===========================MC============================================

  TFile *newfilemc = TFile::Open((""+filemc).c_str());
  TTree* treemc = (TTree*)newfilemc->Get("Tree");
  //m2pi < 0.55 && 
  string selectionmc = Form("mf1 > %g && mf1 < %g && Q2 > %g-5 && Q2 < %g+2",startt,endd, qmin, qmax); 
  RooDataSet datamc("datamc","datamc",RooArgSet(x,Q2,m2pi),Import(*treemc),Cut(selectionmc.c_str()));
  RooKeysPdf kest1("kest1","kest1",x,datamc,RooKeysPdf::MirrorBoth) ;
  RooRealVar mg("mg","mg",0.,-3,3); 
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
  c2.setConstant(true);
  RooGenericPdf P("P","c0 + c1*(mf1-1.28) + c2*(mf1-1.28)*(mf1-1.28)",RooArgSet(x,c0,c1,c2)); 
  
  RooAddPdf model("model","model",RooArgList(resol,P),frac) ;
  //model.plotOn(xframe1);
  // Construct unbinned likelihood of model w.r.t. data
  RooAbsReal* nll = model.createNLL(*Data) ;

  // I n t e r a c t i v e   m i n i m i z a t i o n ,   e r r o r   a n a l y s i s
  // -------------------------------------------------------------------------------

  // Create MINUIT interface object
  //RooMinuit m(*nll) ;

  // Activate verbose logging of MINUIT parameter space stepping
  //m.setVerbose(kTRUE) ;
  //m.setVerbose(kFALSE) ;

  // Call MIGRAD to minimize the likelihood
  RooFitResult* r = model.fitTo(*Data,Save());
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
  Data->plotOn(xframe1);
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

void effic_draw_m2pi(string file, string same){
  const double massMin = 0.35;
  const double massMax = 0.75;

  TFile *newfilem0 = TFile::Open(file.c_str());
  TTree* treem0 = (TTree*)newfilem0->Get("Tree");
  TH1D *h_Q2 = new TH1D("h_Q2","h_Q2",40,0.3,0.8);
  h_Q2->SetXTitle("Q^{2} (GeV)^{2}");
  TCanvas s;
  treem0->Draw("m2pi >> h_Q2", "mf1 > 1.17 && mf1 < 1.4 && m2pi > 0.35 && m2pi < 0.75");

  TTree* treem0gen = (TTree*)newfilem0->Get("Treegen");
  TH1D *h_Q2gen = new TH1D("h_Q2gen","h_Q2gen",40,0.3,0.8);
  treem0gen->Draw("m2pigen >> h_Q2gen");
  s.Close();
  TEfficiency *hEfficiency = new TEfficiency(*h_Q2,*h_Q2gen);
  hEfficiency->SetTitle(" ; m_{2#pi} (GeV^{2}/c^{2}) ; #varepsilon");
  hEfficiency->Draw(same.c_str());


  TFile *newfileDATA = TFile::Open("../eef1/histograms/data.root");
  TTree* tree_data = (TTree*)newfileDATA->Get("Tree");
  TH1D *h_Q2_data = new TH1D("h_Q2_data","h_Q2_data",40,0.3,0.8);
  h_Q2_data->SetXTitle("Q^{2} (GeV)^{2}");
  TCanvas s2;
  tree_data->Draw("m2pi >> h_Q2_data", "mf1 > 1.17 && mf1 < 1.4 && m2pi > 0.35 && m2pi < 0.75");
  s2.Close();
  // A TEfficiency cannot be passed to TH1::Divide, so make an ordinary
  // histogram containing the same bin-by-bin MC efficiency.
  TH1D *h_efficiency_hist = (TH1D*)h_Q2->Clone("h_efficiency_hist");
  h_efficiency_hist->Divide(h_Q2gen);

  // Keep the original histograms unchanged and make corrected copies.
  TH1D *h_Q2_corrected = (TH1D*)h_Q2->Clone("h_Q2_corrected");
  h_Q2_corrected->Divide(h_efficiency_hist);

  TH1D *h_Q2_data_corrected =
      (TH1D*)h_Q2_data->Clone("h_Q2_data_corrected");
  h_Q2_data_corrected->Divide(h_efficiency_hist);

  const int firstBin = h_Q2->FindBin(massMin);
  const int lastBin = h_Q2->FindBin(massMax);
  const double mcRecoYield = h_Q2->Integral(firstBin, lastBin);
  const double mcCorrectedYield =
      h_Q2_corrected->Integral(firstBin, lastBin);
  const double dataRecoYield = h_Q2_data->Integral(firstBin, lastBin);
  const double dataCorrectedYield =
      h_Q2_data_corrected->Integral(firstBin, lastBin);

  if (mcCorrectedYield == 0. || dataRecoYield == 0.) {
    std::cerr << "Cannot estimate the efficiency-model difference: "
              << "one of the normalization yields is zero." << std::endl;
    return;
  }

  // Expected reconstructed data yield if the MC efficiency model were exact.
  const double expectedDataReco =
      dataCorrectedYield / mcCorrectedYield * mcRecoYield;
  const double difference = expectedDataReco - dataRecoYield;
  const double relativeDifference = difference / dataRecoYield * 100.;

  std::cout << std::fixed << std::setprecision(2)
            << "Efficiency-model check in " << massMin << " < m2pi < "
            << massMax << ":\n"
            << "  MC reconstructed yield       = " << mcRecoYield << '\n'
            << "  MC efficiency-corrected yield = " << mcCorrectedYield << '\n'
            << "  Data reconstructed yield      = " << dataRecoYield << '\n'
            << "  Data efficiency-corrected yield = " << dataCorrectedYield << '\n'
            << "  Expected data reconstructed yield = "
            << expectedDataReco << '\n'
            << "  Difference (expected - observed)  = " << difference << '\n'
            << "  Relative difference               = "
            << relativeDifference << "%" << std::endl;
  
  TCanvas *s3 = new TCanvas();
  h_Q2_data_corrected->Draw();

}