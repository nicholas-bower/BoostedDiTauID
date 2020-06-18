import ROOT
import sys
if len(sys.argv) != 3:
   print " USAGE : %s <input file > <output file >"%( sys . argv [0])
inFileName = sys.argv [1]
sample = sys.argv [2]    
print " Reading from ", inFileName , "Dataset is ", sample
inFile = ROOT.TFile.Open(inFileName ,"READ") 
histos = {}
lines= {'EB_hOverE_2D' : .298, 'EB_hOverE_2D_pass' : .298, 'EE_hOverE_2D' : .101, 'EE_hOverE_2D_pass' : .101}
c = ROOT.TCanvas("c","c",600,600)      
c.SetBorderSize(0);                                
c.SetFrameBorderMode(0)                                    
plot2d=False
ROOT.gStyle.SetTitleBorderSize(0)
ROOT.gStyle.SetOptStat(0)

sample = sys.argv [2]
outFolder = './newgraphs/'+ sample+'/'+sample
nlooseMatched = inFile.Get('matched_ePt_loose')
nMatched=inFile.Get('matched_ePt')
NElectron_matched=inFile.Get('NElectron_matched')
NElectron_matched.Draw()
c.SaveAs(outFolder+'_NElectron_matched.png')
#ePt Plots
NElectron_looseHoverE_ePt=inFile.Get('ePt_loose')
#NElectrons=inFile.Get('NElectrons')
#nE = NElectrons.GetBinContent(2)
#print(nE)
#ID Plots
NElectronID = inFile.Get('NElectronID')
NElectronID.SetTitle(sample+': Passing Electron ID')
#NElectronID.Scale(1/nE)
NElectronID.Draw()
#c.SaveAs (outFolder+"NelecID.png")

NIso = inFile.Get('NElectronIso')
NIso.SetTitle(sample+': Passing Non H/E')
#NIso.Scale(1/nE)
NIso.Draw()
#c.SaveAs(outFolder+"_NIso.png")

NOld = inFile.Get('NOld')
NOld.SetTitle(sample+': NElectrons Passing 80X H/E')
#NOld.Scale(1/nE)
NOld.Draw()
#c.SaveAs(outFolder+"_NOld.png")

effHoverE = ROOT.TEfficiency(NElectronID, NIso)
effHoverE.SetTitle(sample+': Efficiency of 94X H/E [Only]')
effHoverE.Draw()
c.SaveAs(outFolder+"_HoverE_eff.png")

#effOvN = ROOT.TEfficiency(NElectronID,NOld)
#effOvN.SetTitle(': Eff of 80X vs 94X (HoverE)')
#effOvN.Draw()
#c.SaveAs(outFolder+'_eff8v9.png')


effOldHoverE = ROOT.TEfficiency(NOld, NIso)
effOldHoverE.SetTitle(sample+': Efficiency of 80XH/E on selected')
effOldHoverE.Draw()
c.SaveAs(outFolder+'_effOldHoverE.png')

matched_NElectronID= inFile.Get('NElectronID_matched')
matched_NOld=inFile.Get('matched_nOld')
matched_NOld.Scale(.5)
c.SaveAs(outFolder+'_effblahblah.png')  
matched_Iso = inFile.Get('matched_nElectronIso')

eff_NIso_matched = ROOT.TEfficiency(matched_Iso, NElectron_matched)
eff_NIso_matched.SetTitle(sample+': Efficiency of 94X non H/E [Matched]')
eff_NIso_matched.Draw()
c.SaveAs(outFolder+'_M_NIsoEff.png')

eff_ElectronID_M= ROOT.TEfficiency(matched_NElectronID, NElectron_matched)
eff_ElectronID_M.SetTitle(sample+': Efficiency of 94X ID [Matched]')
eff_ElectronID_M.Draw()
c.SaveAs(outFolder+"_M_NID_eff.png")

eff_NOld_M=ROOT.TEfficiency(matched_NOld, NElectron_matched)
eff_NOld_M.SetTitle(sample+': Efficiency of 80X ID [Matched]')
eff_NOld_M.Draw()
c.SaveAs(outFolder+"_M_80X_eff.png")

effHoverE_M= ROOT.TEfficiency(matched_NElectronID,matched_Iso)
effHoverE_M.SetTitle(sample+': Efficiency of 94X H/E on selected[Matched]')
effHoverE_M.Draw()
c.SaveAs(outFolder+"_M_HoverE_eff.png")


effOldHoverE_M = ROOT.TEfficiency(matched_NOld, matched_Iso)
effOldHoverE_M.SetTitle(sample+': Efficiency of 80XH/E on selected[Matched]')
effOldHoverE_M.Draw()
c.SaveAs(outFolder+'_M_effOldHoverE.png')

effHoverE = ROOT.TEfficiency(NElectronID, NIso)
effHoverE.SetTitle(sample+': Efficiency of 94X H/E [Only]')
effHoverE.Draw()
c.SaveAs(outFolder+"_HoverE_eff.png")

