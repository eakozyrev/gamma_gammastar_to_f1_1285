#define yield_cxx
#include "yield.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TProfile.h"
#include "TFile.h"
#include "TMatrixTSym.h"
#include "TMath.h"
#include "TVector.h"
#include "TText.h"
#include "TLorentzVector.h"
#include "TStyle.h"
#include "TBenchmark.h"
#include <iostream>
#include "TF1.h"
#include <string>
#include "TGraphErrors.h"

void yield::Loop(std::string histFileName)
{

  if (fChain == 0)
    return;

  Long64_t nentries = fChain->GetEntriesFast();
  TFile *newfile = new TFile(histFileName.c_str(), "recreate");

  TH2D *h_Ptr_P = new TH2D("h_Ptr_P", "h_Ptr_P", 200, -1, 3, 200, -1, 9);
  TH1D *h_cos_el = new TH1D("h_cos_el", "h_cos_el", 200, -1, 1);
  TH1D *h_r = new TH1D("h_r", "h_r", 1000, -7, 7);
  TH1D *h_helicity = new TH1D("h_helicity", "h_helicity", 100, -0.1, 1.1);
  TH1D *h_m2gamma = new TH1D("h_m2gamma", "h_m2gamma", 200, 0.4, 0.7);
  TH1D *h_energy_gamma = new TH1D("h_energy_gamma", "h_energy_gamma", 600, 0, 5);
  TH1D *h_psicosthcm = new TH1D("h_psicosthcm", "h_psicosthcm", 2000, -1.1, 1.1);
  TH1D *h_mf1 = new TH1D("h_mf1", "h_mf1", 200, 0.9, 1.4);
  TH1D *h_Q2 = new TH1D("h_Q2", "h_Q2", 500, -5, 100);
  TH2D *h_Q2_Q2 = new TH2D("h_Q2_Q2", "h_Q2_Q2", 500, -2, 20, 500, -2, 20);
  h_Q2_Q2->SetXTitle("Q^{2}_electron, GeV^{2}");
  h_Q2_Q2->SetYTitle("Q^{2}_positron, GeV^{2}");
  TH2D *h_Q2_mf1 = new TH2D("h_Q2_mf1", "h_Q2_mf1", 500, -5, 100, 300, 0.8, 1.5);
  TH2D *h_Q2_Ptr = new TH2D("h_Q2_Ptr", "h_Q2_Ptr", 500, -5, 100, 700, 0, 7);
  TH1D *h_th_pi = new TH1D("h_th_pi", "h_th_pi", 1000, 0, 3.14);
  h_th_pi->SetXTitle("#theta_{#pi^{#pm}} (rad)");
  TH1D *h_th_el = new TH1D("h_th_el", "h_th_el", 1000, 0, 3.14);
  h_th_el->SetXTitle("#theta_{e^{#pm}} (rad)");
  TH1D *h_energy_el = new TH1D("h_energy_el", "h_energy_el", 2000, 0, 5);
  TH2D *h_m_eta_pi = new TH2D("h_m_eta_pi", "h_m_eta_pi", 100, 0.6, 1.2, 100, 0.6, 1.2);
  h_m_eta_pi->SetXTitle("m_{#pi^{-}#eta} (GeV/c^{2})");
  h_m_eta_pi->SetYTitle("m_{#pi^{+}#eta} (GeV/c^{2})");
  TH2D *h_m_eta_pi_b = new TH2D("h_m_eta_pi_b", "h_m_eta_pi_b", 100, 0.6, 1.2, 100, 0.6, 1.2);
  h_m_eta_pi_b->SetXTitle("m_{#pi^{-}#eta} (GeV/c^{2})");
  h_m_eta_pi_b->SetYTitle("m_{#pi^{+}#eta} (GeV/c^{2})");
  TH1D *h_m_pipeta = new TH1D("h_m_pipeta", "h_m_pipeta", 100, 0.6, 1.2);
  h_m_pipeta->SetXTitle("m_{#pi^{+}#eta} (GeV/c^{2})");
  TH1D *h_m_pimeta = new TH1D("h_m_pimeta", "h_m_pimeta", 100, 0.6, 1.2);
  h_m_pimeta->SetXTitle("m_{#pi^{-}#eta} (GeV/c^{2})");
  TH1D *h_m_pippim = new TH1D("h_m_pippim", "h_m_pippim", 200, 0.2, 1.2);
  h_m_pippim->SetXTitle("m_{#pi^{+}#pi^{-}} (GeV/c^{2})");
  TH1D *h_cosa0 = new TH1D("h_cosa0", "h_cosa0", 30, -1, 1);

  TH2D *h_cosa0_p_m = new TH2D("h_cosa0_p_m", "h_cosa0_p_m", 10, -1, 1, 10, -1, 1);
  h_cosa0_p_m->SetXTitle("cos(#pi-)");
  h_cosa0_p_m->SetYTitle("cos(#pi+)");

  TH2D *h_cosa0_p_m_b = new TH2D("h_cosa0_p_m_b", "h_cosa0_p_m_b", 10, -1, 1, 10, -1, 1);
  h_cosa0_p_m_b->SetXTitle("cos(#pi-)");
  h_cosa0_p_m_b->SetYTitle("cos(#pi+)");

  TH1D *h_piselectorsmap = new TH1D("h_piselectorsmap", "h_piselectorsmap", 10000, 0, 1000000);
  TH1D *h_tr_xy = new TH1D("h_tr_xy", "h_tr_xy", 1000, 0, 5);
  TH1D *h_tr_z = new TH1D("h_tr_z", "h_tr_z", 2000, -10, 10);
  TH1D *h_Q2_mc = new TH1D("h_Q2_mc", "h_Q2_mc", 100, 0, 100);
  TH1D *hcheck = new TH1D("hcheck", "hcheck", 200, -5, 5);

  TTree Tree("Tree", "Tree");
  double mf1, cosa0_p, cosa0_m, Q2, omega, hel, ma0m, ma0p, ifsig, m2pi;
  double Q2gen, m2pigen;
  Tree.Branch("mf1", &mf1, "mf1/D");
  Tree.Branch("cosa0_p", &cosa0_p, "cosa0_p/D");
  Tree.Branch("cosa0_m", &cosa0_m, "cosa0_m/D");
  Tree.Branch("ma0m", &ma0m, "ma0m/D");
  Tree.Branch("ma0p", &ma0p, "ma0p/D");
  Tree.Branch("m2pi", &m2pi, "m2pi/D");
  Tree.Branch("ifsig", &ifsig, "ifsig/D");
  Tree.Branch("Q2", &Q2, "Q2/D");
  Tree.Branch("omega", &omega, "omega/D");
  Tree.Branch("hel", &hel, "hel/D");
  Tree.Branch("m2pigen", &m2pigen, "m2pigen/D");

  TTree Treegen("Treegen", "Treegen");
  Treegen.Branch("Q2gen", &Q2gen, "Q2gen/D");
  Treegen.Branch("m2pigen", &m2pigen, "m2pigen/D");

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry = 0; jentry < nentries; jentry++)
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
      break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    // if (Cut(ientry) < 0) continue;

    double sqrts = sqrt(eee * eee - eepx * eepx - eepy * eepy - eepz * eepz);
    TLorentzVector Pbeam(0., 0., 0., sqrts);
    TLorentzVector P_init_el(0., 0., -sqrts / 2., sqrts / 2.);
    TLorentzVector P_init_pos(0., 0., sqrts / 2., sqrts / 2.);
    TLorentzVector PPh[100];
    TLorentzVector Psi[100], Pel[100], Ppi[100], Peta[100], P_f1[100], Psim[100];
    TLorentzVector Pel_lab[100], Ppi_lab[100];
    TLorentzVector P_eloosen, P_ploosen, P_2g_rest, P_a0p, P_a0m;
    int pip_sim{-1}, pim_sim{-1};
    for (int i = 0; i < mclen; i++)
    {
      Psim[i] = TLorentzVector(mcp3cm[i] * sqrt(1. - mccosthcm[i] * mccosthcm[i]) * cos(mcphicm[i]),
                               mcp3cm[i] * sqrt(1. - mccosthcm[i] * mccosthcm[i]) * sin(mcphicm[i]),
                               mcp3cm[i] * mccosthcm[i],
                               mcenergycm[i]);
      if (mclund[i] == 211 && mclund[mothidx[i]] == 20223) {
        pip_sim = i;
      }
      if (mclund[i] == -211 && mclund[mothidx[i]] == 20223) {
        pim_sim = i;
      }
    }
    // Psim[0] - initial electron
    // Psim[1] - initial positron
    int nel = 3;
    if (mclund[4] == 11)
      nel = 4;
    if (mclund[5] == 11)
      nel = 5;
    int npos = 3;
    if (mclund[4] == -11)
      npos = 4;
    if (mclund[5] == -11)
      npos = 5;
    // Psim[npos] - final positron
    // Psim[nel] - final electron
    // Psim[2] - f1-meson
    m2pigen = (pip_sim >= 0 && pim_sim >= 0) ? (Psim[pip_sim] + Psim[pim_sim]).M() : 0;
    Q2gen = -(Psim[0] - Psim[nel]).Dot(Psim[0] - Psim[nel]);
    double Q2gen_1 = -(Psim[1] - Psim[npos]).Dot(Psim[1] - Psim[npos]);
    h_Q2_mc->Fill(Q2gen);
    h_Q2_Q2->Fill(Q2gen, Q2gen_1);
    Q2gen = -(Psim[0] - Psim[nel]).Dot(Psim[0] - Psim[nel]);
    if (-(Psim[1] - Psim[npos]).Dot(Psim[1] - Psim[npos]) > Q2gen)
      Q2gen = -(Psim[1] - Psim[npos]).Dot(Psim[1] - Psim[npos]);


    Treegen.Fill();

    if (L3outdch == 0 && L3outemc == 0)
      continue;
    if (Bgfmultihadron == 0)
      continue;
    int npip = -1;
    int npim = -1;
    // for(int i = 4; i < 9; i++){if(mclund[i]==221)neta=i;}
    for (int i = 4; i < 9; i++)
    {
      if (mclund[i] == 211)
        npip = i;
    }
    for (int i = 4; i < 9; i++)
    {
      if (mclund[i] == -211)
        npim = i;
    }

    for (int i = 0; i < ngamma; i++)
    {
      PPh[i] = TLorentzVector(gammaenergycm[i] * sqrt(1. - gammacosthcm[i] * gammacosthcm[i]) * cos(gammaphicm[i]),
                              gammaenergycm[i] * sqrt(1. - gammacosthcm[i] * gammacosthcm[i]) * sin(gammaphicm[i]),
                              gammaenergycm[i] * gammacosthcm[i],
                              gammaenergycm[i]);
    }
    for (int i = 0; i < npsi; i++)
    {
      Psi[i] = TLorentzVector(psip3cm[i] * sqrt(1. - psicosthcm[i] * psicosthcm[i]) * cos(psiphicm[i]),
                              psip3cm[i] * sqrt(1. - psicosthcm[i] * psicosthcm[i]) * sin(psiphicm[i]),
                              psip3cm[i] * psicosthcm[i],
                              psienergycm[i]);
    }

    for (int i = 0; i < netap; i++)
    {
      P_f1[i] = TLorentzVector(etapp3cm[i] * sqrt(1. - etapcosthcm[i] * etapcosthcm[i]) * cos(etapphicm[i]),
                               etapp3cm[i] * sqrt(1. - etapcosthcm[i] * etapcosthcm[i]) * sin(etapphicm[i]),
                               etapp3cm[i] * etapcosthcm[i],
                               etapenergycm[i]);
    }

    for (int i = 0; i < nel; i++)
    {
      Pel[i] = TLorentzVector(elp3cm[i] * sqrt(1. - elcosthcm[i] * elcosthcm[i]) * cos(elphicm[i]),
                              elp3cm[i] * sqrt(1. - elcosthcm[i] * elcosthcm[i]) * sin(elphicm[i]),
                              elp3cm[i] * elcosthcm[i],
                              elenergycm[i]);
    }
    for (int i = 0; i < npi; i++)
    {
      Ppi[i] = TLorentzVector(pip3cm[i] * sqrt(1. - picosthcm[i] * picosthcm[i]) * cos(piphicm[i]),
                              pip3cm[i] * sqrt(1. - picosthcm[i] * picosthcm[i]) * sin(piphicm[i]),
                              pip3cm[i] * picosthcm[i],
                              pienergycm[i]);
    }
    for (int i = 0; i < nel; i++)
    {
      Pel_lab[i] = TLorentzVector(elp3[i] * sqrt(1. - elcosth[i] * elcosth[i]) * cos(elphi[i]),
                                  elp3[i] * sqrt(1. - elcosth[i] * elcosth[i]) * sin(elphi[i]),
                                  elp3[i] * elcosth[i],
                                  elenergy[i]);
    }
    for (int i = 0; i < npi; i++)
    {
      Ppi_lab[i] = TLorentzVector(pip3[i] * sqrt(1. - picosth[i] * picosth[i]) * cos(piphi[i]),
                                  pip3[i] * sqrt(1. - picosth[i] * picosth[i]) * sin(piphi[i]),
                                  pip3[i] * picosth[i],
                                  pienergy[i]);
    }
    for (int i = 0; i < neta; i++)
    {
      Peta[i] = TLorentzVector(etap3cm[i] * sqrt(1. - etacosthcm[i] * etacosthcm[i]) * cos(etaphicm[i]),
                               etap3cm[i] * sqrt(1. - etacosthcm[i] * etacosthcm[i]) * sin(etaphicm[i]),
                               etap3cm[i] * etacosthcm[i],
                               etaenergycm[i]);
    }

    for (int i = 0; i < npsi; i++)
    {

      int nfermion = psid1idx[i];
      double diffP = (Pbeam - Psi[i]).P();
      double diffPtr = sqrt(pow((Pbeam - Psi[i]).Px(), 2.) + pow((Pbeam - Psi[i]).Py(), 2.));
      h_Ptr_P->Fill(diffPtr, diffP);
      if (diffPtr > 0.08)
        continue;
      if (diffP < 3.4)
        continue;
      int n_etap = psid2idx[i]; // P_f1[n_etap]
      int n_eta = etapd3idx[n_etap];
      int n_pip = etapd1idx[n_etap];
      int n_pim = etapd2idx[n_etap];
      int n_gamma1 = etad1idx[n_eta];
      int n_gamma2 = etad2idx[n_eta];
      int ntrk1 = pitrkidx[n_pip];
      int ntrk2 = pitrkidx[n_pim];
      P_a0p = Peta[n_eta] + Ppi[n_pip];
      P_a0m = Peta[n_eta] + Ppi[n_pim];

      h_cos_el->Fill(cos(Pel[nfermion].Theta()));
      hel = fabs(PPh[n_gamma1].E() - PPh[n_gamma2].E()) / (PPh[n_gamma1] + PPh[n_gamma2]).P();
      //if(hel > 0.9)continue;
      P_eloosen = P_init_el - Pel[nfermion];
      P_ploosen = P_init_pos - (Pbeam - Pel[nfermion] - P_f1[n_etap]);
      if (psid1lund[i] > 0.)
        P_ploosen = P_init_pos - Pel[nfermion];
      if (psid1lund[i] > 0.)
        P_eloosen = P_init_el - (Pbeam - Pel[nfermion] - P_f1[n_etap]);
      if ((Ppi[n_pim] + Ppi[n_pip]).M() < 0.36)
        continue;
      if (piselectorsmap[pitrkidx[n_pip]] < 40000 || piselectorsmap[pitrkidx[n_pim]] < 40000)
        continue;

      h_mf1->Fill(etapmass[n_etap]);
      mf1 = etapmass[n_etap];
      h_m2gamma->Fill(etarmass[n_eta]);
      if (etarmass[n_eta] < 0.5)
        continue;
      h_psicosthcm->Fill(psicosthcm[i]);
      h_r->Fill((sqrts - Psi[i].E() - Psi[i].P()) / sqrts);
      h_helicity->Fill(hel);
      Q2 = pow((P_init_el - Pel[nfermion]).M(), 2.);
      omega = P_init_el.E() - Pel[nfermion].E();

      if (psid1lund[i] > 0.)
      {
        Q2 = pow((P_init_pos - Pel[nfermion]).M(), 2.);
        omega = P_init_pos.E() - Pel[nfermion].E();
      }
      omega = (P_eloosen + P_ploosen).M();
      h_Q2->Fill(Q2);
      h_Q2_mf1->Fill(Q2, etapmass[n_etap]);
      bool ifsignal = false;
      bool ifbckg = false;
      if (etapmass[n_etap] > 1.26 && etapmass[n_etap] < 1.3)
        ifsignal = true;
      if (etapmass[n_etap] > 1.2 && etapmass[n_etap] < 1.24)
        ifbckg = true;
      if (etapmass[n_etap] > 1.32 && etapmass[n_etap] < 1.36)
        ifbckg = true;

      h_energy_el->Fill(Pel_lab[nfermion].E());
      h_th_pi->Fill(Ppi_lab[n_pim].Theta());
      h_th_pi->Fill(Ppi_lab[n_pip].Theta());
      h_th_el->Fill(Pel_lab[nfermion].Theta());
      h_Q2_Ptr->Fill(Q2, diffPtr);

      h_tr_xy->Fill(Trkdocaxy_xy[ntrk1]);
      h_tr_xy->Fill(Trkdocaxy_xy[ntrk2]);
      h_tr_z->Fill(Trkdocaxy_z[ntrk1]);
      h_tr_z->Fill(Trkdocaxy_z[ntrk2]);
      h_piselectorsmap->Fill(eselectorsmap[eltrkidx[nfermion]]);
      h_piselectorsmap->Fill(piselectorsmap[pitrkidx[n_pim]]);
      ma0m = P_a0m.M() / etapmass[n_etap] * 1.285;
      ma0p = P_a0p.M() / etapmass[n_etap] * 1.285;
      m2pi = (Ppi[n_pip] + Ppi[n_pim]).M();

      h_energy_el->Fill(Pel_lab[nfermion].E());
      h_th_pi->Fill(Ppi_lab[n_pim].Theta());
      h_th_pi->Fill(Ppi_lab[n_pip].Theta());
      h_th_el->Fill(Pel_lab[nfermion].Theta());
      h_Q2_Ptr->Fill(Q2, diffPtr);

      h_tr_xy->Fill(Trkdocaxy_xy[ntrk1]);
      h_tr_xy->Fill(Trkdocaxy_xy[ntrk2]);
      h_tr_z->Fill(Trkdocaxy_z[ntrk1]);
      h_tr_z->Fill(Trkdocaxy_z[ntrk2]);
      h_piselectorsmap->Fill(eselectorsmap[eltrkidx[nfermion]]);
      h_piselectorsmap->Fill(piselectorsmap[pitrkidx[n_pim]]);
      ma0m = P_a0m.M() / etapmass[n_etap] * 1.285;
      ma0p = P_a0p.M() / etapmass[n_etap] * 1.285;
      ifsig = -1;
      if (ifsignal == true)
        ifsig = 1.;
      if (ifbckg == true)
        ifsig = 0.;
      if (ifsignal == true)
        h_m_eta_pi->Fill(P_a0m.M(), P_a0p.M());
      if (ifbckg == true)
        h_m_eta_pi_b->Fill(P_a0m.M(), P_a0p.M());
      if (ifsignal == true)
        h_m_pipeta->Fill(P_a0p.M());
      if (ifsignal == true)
        h_m_pimeta->Fill(P_a0m.M());
      h_m_pippim->Fill((Ppi[n_pip] + Ppi[n_pim]).M());
      h_energy_gamma->Fill(PPh[n_gamma1].E());
      h_energy_gamma->Fill(PPh[n_gamma2].E());

      double betax = P_f1[n_etap].Px() / P_f1[n_etap].E();
      double betay = P_f1[n_etap].Py() / P_f1[n_etap].E();
      double betaz = P_f1[n_etap].Pz() / P_f1[n_etap].E();
      P_a0p.Boost(-betax, -betay, -betaz);
      P_a0m.Boost(-betax, -betay, -betaz);
      Ppi[n_pim].Boost(-betax, -betay, -betaz);
      Ppi[n_pip].Boost(-betax, -betay, -betaz);

      P_eloosen.Boost(-betax, -betay, -betaz);
      P_2g_rest = P_eloosen;
      cosa0_p = (P_2g_rest.E() * Ppi[n_pim].E() - P_2g_rest.Dot(Ppi[n_pim])) / P_2g_rest.P() / Ppi[n_pim].P();
      cosa0_m = (P_2g_rest.E() * Ppi[n_pip].E() - P_2g_rest.Dot(Ppi[n_pip])) / P_2g_rest.P() / Ppi[n_pip].P();
      if (ifsignal == true)
        h_cosa0_p_m->Fill(cosa0_p, cosa0_m);
      if (ifbckg == true)
        h_cosa0_p_m_b->Fill(cosa0_p, cosa0_m);

      Tree.Fill();
    }
  }
  newfile->Write();
  h_cosa0_p_m->Draw("colz");
  h_cosa0_p_m->SaveAs((histFileName + "_h_cosa0_p_m.eps").c_str());
  h_cosa0_p_m_b->SaveAs((histFileName + "_h_cosa0_p_m_b.eps").c_str());
  h_m_eta_pi->SaveAs((histFileName + "_h_m_eta_pi.eps").c_str());
  h_m_eta_pi_b->SaveAs((histFileName + "_h_m_eta_pi_b.eps").c_str());
}

TH1D *h_signal1, *h_signal2, *h_bkg;

double f_signal(double *x, double *par)
{

  int nbins = h_signal1->GetNbinsX();
  int nb = (x[0] + 1) / 2. * nbins;
  return par[0] * (par[1] * h_signal1->GetArray()[nb + 1] / (double)h_signal1->GetEntries() + (1. - par[1]) * h_signal2->GetArray()[nb + 1] / (double)h_signal2->GetEntries()) + par[2];
}

double f_signal1(double *x, double *par)
{

  int nbins = h_signal1->GetNbinsX();
  int nb = (x[0] + 1) / 2. * nbins;
  return par[0] * (par[1] * h_signal1->GetArray()[nb + 1] / (double)h_signal1->GetEntries() + (1. - par[1]) * h_signal2->GetArray()[nb + 1] / (double)h_signal2->GetEntries() * 0.);
}

double result[5];

double f_signal2(double *x, double *par)
{

  int nbins = h_signal1->GetNbinsX();
  int nb = (x[0] + 1) / 2. * nbins;
  return par[0] * (par[1] * h_signal1->GetArray()[nb + 1] / (double)h_signal1->GetEntries() * 0. + (1. - par[1]) * h_signal2->GetArray()[nb + 1] / (double)h_signal2->GetEntries());
}

void fit(double q1 = 0, double q2 = 10)
{

  int nbins = 33;
  TFile *newfile = TFile::Open("run2.root");
  TTree *tree = (TTree *)newfile->Get("Tree");
  TH1D *h_data = new TH1D("h_data", Form("%g < Q^{2} < %g (GeV^{2})", q1, q2), nbins, -1, 1);
  h_data->SetXTitle("cos(#alpha_{#pi}*)");
  h_data->SetYTitle("yields");
  TH1D *h_bkg = new TH1D("h_bkg", Form("%g < Q^{2} < %g (GeV^{2})", q1, q2), nbins, -1, 1);
  TH1D *h_bkg1 = new TH1D("h_bkg1", Form("%g < Q^{2} < %g (GeV^{2})", q1, q2), nbins, -1, 1);
  string selec = "mf1 > 1.26 && mf1 < 1.3 && ";
  string selec_bkg1 = "mf1 > 1.18 && mf1 < 1.2 &&";
  string selec_bkg2 = "mf1 > 1.38 && mf1 < 1.4 &&";
  tree->Draw("cosa0 >> h_data", Form((selec + "Q2 > %g && Q2 < %g").c_str(), q1, q2));
  cout << h_data->GetEntries() << endl;
  tree->Draw("cosa0 >> h_bkg1", Form((selec_bkg1 + "Q2 > %g && Q2 < %g").c_str(), q1, q2));
  tree->Draw("cosa0 >> h_bkg", Form((selec_bkg2 + "Q2 > %g && Q2 < %g").c_str(), q1, q2));
  *h_bkg = *h_bkg + *h_bkg1;
  TFile *newfilem0 = TFile::Open("mc_m0.root");
  TTree *treem0 = (TTree *)newfilem0->Get("Tree");
  TH1D *h_m0 = new TH1D("h_m0", "h_m0", nbins, -1, 1);
  treem0->Draw("cosa0 >> h_m0", Form((selec + "Q2 > %g && Q2 < %g").c_str(), q1 - 5, q2 + 3));
  h_signal1 = h_m0;

  TFile *newfilem1 = TFile::Open("mc_m1.root");
  TTree *treem1 = (TTree *)newfilem1->Get("Tree");
  TH1D *h_m1 = new TH1D("h_m1", "h_m1", nbins, -1, 1);
  treem1->Draw("cosa0 >> h_m1", Form((selec + "Q2 > %g && Q2 < %g").c_str(), q1 - 5, q2 + 3));
  h_signal2 = h_m1;

  TF1 *ffit = new TF1("ffit", f_signal, -1, 1, 3);
  ffit->FixParameter(2, 0.);
  ffit->SetNpx(10000);
  ffit->SetLineWidth(4.);
  h_data->Add(h_bkg, -1);
  h_data->Fit(ffit);
  h_data->SetLineWidth(3);
  h_data->SetLineColor(1);
  h_data->SetMarkerStyle(21);
  h_data->Draw("e");
  double par1 = ffit->GetParameter(0);
  double par2 = ffit->GetParameter(1);
  result[0] = ffit->GetParameter(1);
  result[1] = ffit->GetParError(1);
  TF1 *ffit1 = new TF1("ffit1", f_signal1, -1, 1, 2);
  ffit1->SetParameters(par1, par2);
  ffit1->SetLineColor(3);
  ffit1->SetLineStyle(2);
  ffit1->SetLineWidth(4);
  ffit1->SetNpx(10000);
  ffit1->SetTitle("m = 0");
  ffit1->Draw("same");
  h_bkg->Draw("same");

  TF1 *ffit2 = new TF1("ffit2", f_signal2, -1, 1, 2);
  ffit2->SetParameters(par1, par2);
  ffit2->SetLineColor(4);
  ffit2->SetLineStyle(9);
  ffit2->SetLineWidth(4);
  ffit2->SetNpx(10000);
  ffit2->SetTitle("m = 1");
  ffit2->Draw("same");
}

void draw_Q2()
{

  double en[] = {2, 2.6, 3.2, 4., 5., 6., 7, 9, 12, 20, 30};
  double den[1000];
  double cr[100], dcr[100];

  for (int i = 0; i < 9; i++)
  {
    fit(en[i], en[i + 1]);
    cr[i] = result[0];
    dcr[i] = result[1];
    den[i] = (en[i + 1] - en[i]) / 2.;
    en[i] = (en[i] + en[i + 1]) / 2.;
  }

  TCanvas *s0 = new TCanvas();
  TH1F *frd = s0->DrawFrame(-1, 0., 20, 1.);
  frd->SetXTitle("Q^{2} (GeV^{2})");
  frd->SetYTitle("#sigma_{E}/E");
  TGraphErrors *Crossu = new TGraphErrors(9, en, cr, den, dcr);
  Crossu->SetMarkerColor(2);
  Crossu->SetMarkerStyle(20);
  Crossu->SetLineColor(2);
  Crossu->SetLineWidth(2.);
  Crossu->SetTitle("");
  Crossu->Draw("P");
}
