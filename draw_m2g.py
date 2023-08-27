from __future__ import print_function
from ROOT import TCanvas, TGraphErrors
from ROOT import gROOT, TF1, TFile, gStyle
from math import sin
from array import array


file = TFile('histograms/data.root')
h1f = gROOT.FindObject( 'h_m2gamma' )
file1 = TFile('histograms/mc_m0.root')
h1fm = gROOT.FindObject( 'h_m2gamma' )

gStyle.SetOptStat(000000)
gStyle.SetPadLeftMargin(0.13);
gStyle.SetPadBottomMargin(0.14);

h1f.SetNdivisions(6,"X")
h1fm.SetXTitle("m_{#gamma#gamma}, GeV/c^{2}")
h1fm.SetYTitle("yields")
h1fm.SetTitle("")
h1fm.SetTitleOffset(0.87,"Y")
h1f.Rebin(2)
h1fm.Rebin(2)
h1fm.SetAxisRange(0.45,0.65,"X")
h1f.SetLineColor(1)
h1f.SetLineWidth(3)
h1f.SetTitleSize(0.07,"X")
h1f.SetTitleSize(0.07,"Y")
h1f.SetLabelSize(0.07,"X")
h1f.SetLabelSize(0.07,"Y")
s1  = TCanvas();
h1fm.SetNormFactor(h1f.GetEntries())
h1fm.SetFillColor(3)
h1fm.SetLineColor(3)
h1fm.SetFillStyle(3001)
h1fm.Draw()
h1f.Draw("esame")
s1.SaveAs("plots_for_paper/P_m2g.eps")

raw_input()
