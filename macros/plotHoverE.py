import ROOT
import sys
if len(sys.argv) != 3:
   print " USAGE : %s <input file > <output file >"%( sys . argv [0])
inFileName = sys.argv [1]
outFileName = sys.argv [2]
print " Reading from ", inFileName , "Dataset is ", outFileName
inFile = ROOT.TFile.Open(inFileName ,"READ") 
histos = {}
lines= {'EB_hOverE_2D' : .236, 'EB_hOverE_2D_pass' : .236, 'EE_hOverE_2D' : .0801, 'EE_hOverE_2D_pass' : .0801}
c = ROOT.TCanvas("c","c",600,600)      
c.SetBorderSize(0);   
plot2d = False
c.SetFrameBorderMode(0)                                    
ROOT.gStyle.SetTitleBorderSize(0)
ROOT.gStyle.SetOptStat(0)         
for key in inFile.GetListOfKeys():
   h = key.ReadObj()
   if h.ClassName() == 'TH1F':
      h.SetTitle(outFileName + ': ' + h.GetName())
      histos[h.GetName()]= inFile.Get(str(h.GetName()))
      if h.GetName() == 'EB_hOverE' or h.GetName()=='EE_hOverE':
         h.SetFillColor(4)
         h.SetTitle(h.GetName())
         #h.Draw('BAR')
         #fileName = "./newgraphs/"+outFileName+ "_"+h.GetTitle()  +'.png'
         #c.SaveAs(fileName)
         #ROOT.gPad.SetLogy()
         #fileName = "./newgraphs/"+outFileName+ "_"+h.GetTitle()  +'_LOGY.png'
         #c.SaveAs(fileName)  
   if  h.ClassName() == 'TH2F' and plot2d==True:
      h.SetTitle(outFileName + ': ' + h.GetName()) 
      h.GetXaxis().SetTitle("H/E") 
      h.GetYaxis().SetTitle("Loose(Rho,Esc)")
      h.Draw("COL Z CJUST")
      flatCut = lines[h.GetName()] 
      oldCutLine = ROOT.TLine(flatCut,0,flatCut,0.3) 
      oldCutLine.SetLineColor(2) 
      oldCutLine.SetLineWidth(2)
      oldCutLine.Draw()
      newCutLine = ROOT.TLine(0,0,.3,.3)
      newCutLine.SetLineColor(6)
      newCutLine.SetLineWidth(2)                     
      newCutLine.Draw()  
      fileName = "./newgraphs/"+outFileName+ "_"+h.GetTitle()  +'.png'
      print "Printing 2DHist for " + h.GetTitle() +"to"+ fileName 
      c.SaveAs(fileName)  
      #print h.GetTitle()
nmatched = inFile.Get('nmatched')
nmatched.SetTitle(outFileName+' Matched/UnMatched')
nmatched.Draw()
c.SaveAs('./newgraphs/'+outFileName+'_nmatched.png')

Ec_HoverE_2D = inFile.Get("Ec_hOverE_2D")
Ec_HoverE_2D.SetTitle(outFileName+' Esc/HoverE')
Ec_HoverE_2D.GetXaxis().SetTitle("H/E")
Ec_HoverE_2D.GetYaxis().SetTitle("Esc")
Ec_HoverE_2D.Draw("COL Z CJUST")
c.SaveAs('./newgraphs/'+outFileName+'_EcHoverE.png')

EE_Ec_HoverE_2D = inFile.Get("EE_Ec_hOverE_2D")
EE_Ec_HoverE_2D.SetTitle(outFileName+' Esc/HoverE EE')
EE_Ec_HoverE_2D.GetXaxis().SetTitle("H/E")
EE_Ec_HoverE_2D.GetYaxis().SetTitle("Esc")
EE_Ec_HoverE_2D.Draw("COL Z CJUST")
c.SaveAs('./newgraphs/'+outFileName+'_EE_EcHoverE.png')

EB_Ec_HoverE_2D = inFile.Get("EB_Ec_hOverE_2D")
EB_Ec_HoverE_2D.SetTitle(outFileName+' Esc/HoverE EB')
EB_Ec_HoverE_2D.GetXaxis().SetTitle("H/E")
EB_Ec_HoverE_2D.GetYaxis().SetTitle("Esc")
EB_Ec_HoverE_2D.Draw("COL Z CJUST")
c.SaveAs('./newgraphs/'+outFileName+'_EB_EcHoverE.png')

Ec_HoverE_2D = inFile.Get("Ec_hOverE_2D_pass")
Ec_HoverE_2D.SetTitle(outFileName+' Esc/HoverE pass')
Ec_HoverE_2D.GetXaxis().SetTitle("H/E")
Ec_HoverE_2D.GetYaxis().SetTitle("Esc")
Ec_HoverE_2D.Draw("COL Z CJUST")
c.SaveAs('./newgraphs/'+outFileName+'_EcHoverE_pass.png')

EE_Ec_HoverE_2D = inFile.Get("EE_Ec_hOverE_2D_pass")
EE_Ec_HoverE_2D.SetTitle(outFileName+' Esc/HoverE EE pass')
EE_Ec_HoverE_2D.GetXaxis().SetTitle("H/E")
EE_Ec_HoverE_2D.GetYaxis().SetTitle("Esc")

EE_Ec_HoverE_2D.Draw("COL Z CJUST")
c.SaveAs('./newgraphs/'+outFileName+'_EE_EcHoverE_pass.png')

EB_Ec_HoverE_2D = inFile.Get("EB_Ec_hOverE_2D_pass")
EB_Ec_HoverE_2D.SetTitle(outFileName+' Esc/HoverE EB pass')
EB_Ec_HoverE_2D.GetXaxis().SetTitle("H/E")
EB_Ec_HoverE_2D.GetYaxis().SetTitle("Esc")
EB_Ec_HoverE_2D.Draw("COL Z CJUST")
c.SaveAs('./newgraphs/'+outFileName+'_EB_EcHoverE_pass.png')

hOverE = inFile.Get('hOverE')
hOverE.SetTitle('H/E')
hOverE.Draw('BAR')
c.SaveAs('./newgraphs/'+outFileName+'_hOverE.png')

hOverE_pass = inFile.Get('hOverE_pass')
hOverE_pass.SetTitle(outFileName+' H/E Pass')
hOverE_pass.Draw('BAR')
c.SaveAs('./newgraphs/'+outFileName+'_hOverE_pass.png')

EE_hOverE_2D_pass= inFile.Get('EE_hOverE_2D_pass')
EE_hOverE_2D_pass.SetTitle(outFileName+' EE H/E vs Cut pass')
EE_hOverE_2D_pass.GetXaxis().SetTitle("H/E")
EE_hOverE_2D_pass.GetYaxis().SetTitle("Loose(Rho,Esc)")
EE_hOverE_2D_pass.Draw("COL Z CJUST")
flatCut = lines[EE_hOverE_2D_pass.GetName()]
oldCutLine = ROOT.TLine(flatCut, 0, flatCut, 0.3)
oldCutLine.SetLineColor(2)
oldCutLine.SetLineWidth(2)
oldCutLine.Draw()
newCutLine = ROOT.TLine(0, 0, .3, .3)
newCutLine.SetLineColor(6)
newCutLine.SetLineWidth(2)
newCutLine.Draw()
fileName = "./newgraphs/"+outFileName+"_EE_HoverE_Pass.png"
c.SaveAs(fileName)

EB_hOverE_2D_pass= inFile.Get('EB_hOverE_2D_pass')
EB_hOverE_2D_pass.SetTitle(outFileName+' EE H/E vs Cut pass')
EB_hOverE_2D_pass.GetXaxis().SetTitle("H/E")
EB_hOverE_2D_pass.GetYaxis().SetTitle("Loose(Rho,Esc)")
EB_hOverE_2D_pass.Draw("COL Z CJUST")
flatCut = lines[EB_hOverE_2D_pass.GetName()]
oldCutLine = ROOT.TLine(flatCut, 0, flatCut, 0.3)
oldCutLine.SetLineColor(2)
oldCutLine.SetLineWidth(2)
oldCutLine.Draw()
newCutLine = ROOT.TLine(0, 0, .3, .3)
newCutLine.SetLineColor(6)
newCutLine.SetLineWidth(2)
newCutLine.Draw()
fileName = "./newgraphs/"+outFileName+"_EB_HoverE_Pass.png"
c.SaveAs(fileName)