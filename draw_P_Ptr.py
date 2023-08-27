from __future__ import print_function
from ROOT import TCanvas, TGraphErrors
from ROOT import gROOT, TF1, TFile, gStyle
from math import sin
from array import array

file = TFile('histograms/data.root')
#file = TFile('histograms/mc_m0.root')

gStyle.SetOptStat(000000)
gStyle.SetPadLeftMargin(0.13);
gStyle.SetPadBottomMargin(0.14);
h1f = gROOT.FindObject( 'h_Ptr_P' )
h1f.SetNdivisions(6,"X")
h1f.SetYTitle("P, GeV/c")
h1f.SetXTitle("P_{tr}, GeV/c")
h1f.SetTitle("")
h1f.RebinX(8)
h1f.RebinY(8)
h1f.SetAxisRange(0,2,"X")
h1f.SetAxisRange(2,7,"Y")
h1f.SetLineColor(1)
h1f.SetLineWidth(3)
h1f.SetTitleSize(0.07,"X")
h1f.SetTitleSize(0.07,"Y")
h1f.SetLabelSize(0.07,"X")
h1f.SetLabelSize(0.07,"Y")
s1  = TCanvas();
h1f.Draw("box")
s1.SaveAs("plots_for_paper/P_Ptr_data.eps")
#s1.SaveAs("plots_for_paper/P_Ptr_mc.eps")

raw_input()
