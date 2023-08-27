from __future__ import print_function
from ROOT import TCanvas, TGraphErrors
from ROOT import gROOT, TF1
from math import sin
from array import array
 
f = open("results/cross.dat","r");
line = f.readlines()
x0 = array( 'd' )
dx0 = array( 'd' )
y0 = array( 'd' )
dy0 = array( 'd' )
y1 = array( 'd' )
dy1 = array( 'd' )
z = 0
for i in line:
    print(i)
    arr = i.split(' ')
    print(arr[10])
    x0.append(float(arr[0])/2.+float(arr[1])/2.)
    dx0.append(-float(arr[0])/2.+float(arr[1])/2.)
    dy1.append(float(arr[10]))
    y1.append(float(arr[8]))
    dy0.append(float(arr[7]))
    y0.append(float(arr[6]))    
    z = z + 1

 
 
c1 = TCanvas( 'c1', 'A Simple Graph Example', 200, 10, 700, 500 )
 
gr = TGraphErrors(z,x0,y0,dx0,dy0)
gr.SetLineColor(2)
gr.SetLineWidth( 4 )
gr.SetMarkerColor( 2 )
gr.SetMarkerStyle( 21 )
gr.SetTitle("|G'-xF'/#nu^{2}|")

gr1 = TGraphErrors(z,x0,y1,dx0,dy1)
gr1.SetLineColor(4)
gr1.SetLineWidth(4)
gr1.SetMarkerColor(4)
gr1.SetMarkerStyle(81)
gr1.SetTitle("|G'|")
gr1.GetXaxis().SetTitle('Q^{2} (GeV^{2})')
gr1.GetYaxis().SetTitle('formfactor')
gr1.GetXaxis().SetLabelSize(0.06)
gr1.GetXaxis().SetTitleSize(0.06)
gr1.GetYaxis().SetLabelSize(0.06)
gr1.GetYaxis().SetTitleSize(0.06)

gr1.Draw("AP")
gr.Draw("P")
# TCanvas.Update() draws the frame, after which one can change it
c1.Update()
#c1.GetFrame().SetFillColor( 21 )
#c1.GetFrame().SetBorderSize( 12 )
c1.Modified()
c1.Update()

# Defining external function:
def TFF_Mil(e,par) :
    M2 = 1.271*1.271*1.271*1.271;
    Mrho2 = 0.77*0.77;
    g2 = 2.5/10000.;
    res = 2.*g2*M2/(e[0]+Mrho2)/Mrho2;
    return par[0]*res

u_fpr = TF1("u_fpr",TFF_Mil,0,40,1);
#u_fpr.Draw("same");
gr1.Fit("u_fpr")

raw_input('Press Enter to exit')    

