import ROOT,sys,os
from DataFormats.FWLite import Events, Handle
from looseElectron import *
import numpy as np

ROOT.gROOT.SetStyle('Plain')
ROOT.gROOT.SetBatch()

isCheckVar = False
isTCP = False

process = sys.argv[2]#Var or eff
sample = sys.argv[3]##Name of sample 

inputFileListName=sys.argv[1]

inputFileList=inputFileListName

#if len(sys.argv)>2:
#    outputFileDir=sys.argv[2]
#else:
#outputFileDir = "./plots/"
#outputFileDir = "root://cmseos.fnal.gov//store/user/nbower/QCD_var_plots/"
outputFileDir = "root://cmseos.fnal.gov//store/user/nbower/"+sample+"_"+process+"_plots/"
#outputFileDir = "root://cmseos.fnal.gov//store/user/nbower/TCP_var_plots/" 
outputFileName = outputFileDir+"h_"+inputFileListName.split("/")[-1].replace(".txt", "_"+process+".root")
print outputFileName
print "process is "+ process
out=ROOT.TFile.Open(outputFileName,'recreate')

handleMuon = Handle ('vector<pat::Muon>')
labelMuon = ('slimmedMuons')

handleElectron = Handle ('vector<pat::Electron>')
labelElectron = ('slimmedElectrons')

handleJet = Handle ('vector<pat::Jet>')
labelJet = ('slimmedJets')

handleTaus = Handle ('vector<pat::Tau>')
labelTaus = ('slimmedTaus')

handleElectronCleanedTaus = Handle ('vector<pat::Tau>')
labelElectronCleanedTaus = ('slimmedTausElectronCleaned', '', 'PAT')

handleMuonCleanedTaus = Handle ('vector<pat::Tau>')
labelMuonCleanedTaus = ('slimmedTausMuonCleaned', '', 'PAT')

handleVertex = Handle ('vector<reco::Vertex>')
labelVertex = ('offlineSlimmedPrimaryVertices')

handleGenParticle = Handle ('vector<reco::GenParticle>')
labelGenParticle = ('genParticles')

handleHLT = Handle ('edm::TriggerResults')
labelHLT = ('TriggerResults')

handleRho = Handle ('double')
labelRho = ('fixedGridRhoFastjetAll')

handleGenParticle = Handle ('vector<reco::GenParticle>')
labelGenParticle = ('prunedGenParticles')

handleGenInfo = Handle('GenEventInfoProduct')
labelGenInfo = ( 'generator' )

# Histograms

h={}

h['NEvent'] = ROOT.TH1F ("NEvent", "", 2, 0, 2)
h['NElectron'] = ROOT.TH1F ("NElectron", "", 6, 0, 6)
h['NElectron_matched'] = ROOT.TH1F ("NElectron_matched", "", 4, 0, 4)
h['NElectronID'] = ROOT.TH1F ("NElectronID", "", 4, 0, 4)
h['NElectronID_matched'] = ROOT.TH1F ("NElectronID_matched", "", 4, 0, 4)

h['ePt'] = ROOT.TH1F("ePt", "", 500, 0, 500)
h['ePt_veto'] = ROOT.TH1F("ePt_veto", "", 500, 0, 500)
h['ePt_loose'] = ROOT.TH1F("ePt_loose", "", 500, 0, 500)
h['ePt_medium'] = ROOT.TH1F("ePt_medium", "", 500, 0, 500)
h['ePt_tight'] = ROOT.TH1F("ePt_tight", "", 500, 0, 500)
h['matched_ePt'] = ROOT.TH1F("matched_ePt", "", 500, 0, 500)
h['matched_ePt_veto'] = ROOT.TH1F("matched_ePt_veto", "", 500, 0, 500)
h['matched_ePt_loose'] = ROOT.TH1F("matched_ePt_loose", "", 500, 0, 500)
h['matched_ePt_medium'] = ROOT.TH1F("matched_ePt_medium", "", 500, 0, 500)
h['matched_ePt_tight'] = ROOT.TH1F("matched_ePt_tight", "", 500, 0, 500)

h['gen_ePt'] = ROOT.TH1F("gen_ePt", "", 500, 0, 500)
h['matched_gen_ePt'] = ROOT.TH1F("matched_gen_ePt", "", 500, 0, 500)
h['matched_gen_ePt_veto'] = ROOT.TH1F("matched_gen_ePt_veto", "", 500, 0, 500)
h['matched_gen_ePt_loose'] = ROOT.TH1F("matched_gen_ePt_loose", "", 500, 0, 500)
h['matched_gen_ePt_medium'] = ROOT.TH1F("matched_gen_ePt_medium", "", 500, 0, 500)
h['matched_gen_ePt_tight'] = ROOT.TH1F("matched_gen_ePt_tight", "", 500, 0, 500)

h['EE_sigmaIetaIeta'] = ROOT.TH1F ("EE_sigmaIetaIeta", "", 480, 0, 60E-3)
h['EE_dEtaSeed'] = ROOT.TH1F ("EE_dEtaSeed", "", 240, -15E-3, 15E-3)
h['EE_dPhiIn'] = ROOT.TH1F ("EE_dPhiIn", "", 160, -0.2, 0.2)
h['EE_hOverE'] = ROOT.TH1F ("EE_hOverE", "", 160, 0, 0.3)
h['EE_hOverE_2D'] = ROOT.TH2F ("EE_hOverE_2D", "", 160, 0, 0.3, 160, 0, 0.3)
h['EE_hOverE_2D_pass'] = ROOT.TH2F ("EE_hOverE_2D_pass", "", 160, 0, 0.3, 160, 0, 0.3)
h['EE_relIso'] = ROOT.TH1F ("EE_relIso", "", 200, 0, 0.5)
h['EE_eInvpInv'] = ROOT.TH1F ("EE_eInvpInv", "", 160, 0, 0.2)
h['EE_missHits'] = ROOT.TH1F ("EE_missHits", "", 5, -0.5, 4.5)               
h['EE_dz'] = ROOT.TH1F ("EE_dz", "", 160,-0.2, 0.2)
h['EE_dxdy'] = ROOT.TH1F ("EE_dxdy", "", 160,-0.2, 0.2)
h['EE_Ec_hOverE_2D'] = ROOT.TH2F ("EE_Ec_hOverE_2D", "", 160, 0, 0.3, 500, 0, 500)
h['EE_Ec']= ROOT.TH1F ("EE_Ec", "", 500, 0, 500)
h['EE_Ec_pass']= ROOT.TH1F("EE_Ec_pass", "", 500, 0, 500)
h['EE_Ec_hOverE_2D_pass'] = ROOT.TH2F ("EE_Ec_hOverE_2D_pass", "", 160, 0, 0.3, 500, 0, 500)

h['hOverE_2D_pass'] = ROOT.TH2F ("hOverE_2D_pass", "", 160, 0, 0.3, 160, 0, 0.3)
h['hOverE_pass'] = ROOT.TH1F ("hOverE_pass", "", 160, 0, 0.3)
h['Ec_hOverE_2D_pass'] = ROOT.TH2F ("Ec_hOverE_2D_pass", "", 160, 0, 0.3, 500, 0, 500)

h['Ec_hOverE_2D'] = ROOT.TH2F ("Ec_hOverE_2D", "", 160, 0, 0.3, 500, 0, 500)
h['Ec']= ROOT.TH1F ("Ec", "", 500, 0, 500)
h['hOverE'] = ROOT.TH1F ("hOverE", "", 500, 0, 500)
h['Ec_pass']= ROOT.TH1F("Ec_pass", "", 500, 0, 500)

h['EB_Ec_hOverE_2D'] = ROOT.TH2F ("EB_Ec_hOverE_2D", "", 160, 0, 0.3, 500, 0, 500)
h['EB_Ec_hOverE_2D_pass'] = ROOT.TH2F ("EB_Ec_hOverE_2D_pass", "", 160, 0, 0.3, 500, 0, 500)
h['EB_Ec']= ROOT.TH1F ("EB_Ec", "", 160, 0, 0.3)
h['EB_Ec_pass']= ROOT.TH1F("EB_Ec_pass", "", 160, 0, 0.3)
h['EB_sigmaIetaIeta'] = ROOT.TH1F ("EB_sigmaIetaIeta", "", 160, 0, 20E-3)
h['EB_dEtaSeed'] = ROOT.TH1F ("EB_dEtaSeed", "", 240, -15E-3, 15E-3)
h['EB_dPhiIn'] = ROOT.TH1F ("EB_dPhiIn", "", 80, -0.1, 0.1)
h['EB_hOverE'] = ROOT.TH1F ("EB_hOverE", "", 160, 0, 0.3)
h['EB_hOverE_2D'] = ROOT.TH2F ("EB_hOverE_2D", "", 160, 0, 0.3, 160, 0, 0.3)
h['EB_hOverE_2D_pass'] = ROOT.TH2F ("EB_hOverE_2D_pass", "", 160, 0, 0.3, 160, 0, 0.3)
h['EB_relIso'] = ROOT.TH1F ("EB_relIso", "", 200, 0, 0.5)
h['EB_eInvpInv'] = ROOT.TH1F ("EB_eInvpInv", "", 160, 0, 0.3)
h['EB_missHits'] = ROOT.TH1F ("EB_missHits", "", 5, -0.5, 4.5)               
h['EB_dz'] = ROOT.TH1F ("EB_dz", "", 160,-0.2, 0.2)
h['EB_dxdy'] = ROOT.TH1F ("EB_dxdy", "", 160,-0.2, 0.2)
######################## MY Efficiency Histos###################################
h['EE_looseIso_ePt'] = ROOT.TH1F("EE_looseIso_ePt", "", 500, 0, 500) 
h['EE_nElectronIso'] = ROOT.TH1F ("EE_nElectronIso", "", 4, 0, 4)


h['EB_looseIso_ePt'] = ROOT.TH1F("EB_looseIso_ePt", "", 500, 0, 500) 
h['EB_nElectronIso'] = ROOT.TH1F ("EB_nElectronIso", "", 4, 0, 4)


h['NElectronIso'] = ROOT.TH1F ("NElectronIso", "", 4, 0, 4)
h['NElectron_looseIso_ePt'] = ROOT.TH1F("NElectron_looseIso_ePt", "", 500, 0, 500)

h['matched_nElectronIso'] = ROOT.TH1F ("matched_nElectronIso", "", 4, 0, 4)
h['matched_nElectron_looseIso_ePt'] = ROOT.TH1F("matched_nElectron_looseIso_ePt", "", 500, 0, 500)
h['matched_gen_nElectron_looseIso_ePt'] = ROOT.TH1F("matched_gen_nElectron_looseIso_ePt", "", 500, 0, 500)

h['EE_nOld'] = ROOT.TH1F ("EE_nOld", "", 4, 0, 4)
h['EB_nOld'] = ROOT.TH1F ("EB_nOld", "", 4, 0, 4)

h['NOld'] = ROOT.TH1F ("NOld", "", 4, 0, 4)
h['NOld_loose_ePt']= ROOT.TH1F('NOld_loose_ePt', "", 500, 0, 500)

h['fake_matched_NOld']= ROOT.TH1F ("fake_matched_NOld", "", 4, 0, 4)


h['matched_NOld'] = ROOT.TH1F ("matched_nOld", "", 4, 0, 4)
h['matched_NOld_loose_ePt']= ROOT.TH1F('matched_nOld_loose_ePt', "", 500, 0, 500)

h['matched_gen_nOld_loose_ePt']= ROOT.TH1F('matched_gen_nOld_loose_ePt', "", 500, 0, 500)
h['nmatched'] =ROOT.TH1F ("nmatched", "", 2,0,2)
# workingPoint={
#     'barrel':{
#         'veto':{'sigmaIEIE':0.0115,'dEtaSeed':0.00749,'dPhi':0.228,'invEinvP':0.299,'miss':2},
#         'loose':{'sigmaIEIE':0.011,'dEtaSeed':0.00477,'dPhi':0.222,'invEinvP':0.299,'miss':1},
#         'medium':{'sigmaIEIE':0.00998,'dEtaSeed':0.00311,'dPhi':0.103,'invEinvP':0.134,'miss':1},
#         'tight':{'sigmaIEIE':0.00998,'dEtaSeed':0.00308,'dPhi':0.0816,'invEinvP':0.0129,'miss':1}
#     },
#     'endcap':{
#         'veto':{'sigmaIEIE':0.037,'dEtaSeed':0.00895,'dPhi':0.213,'invEinvP':0.15,'miss':3},
#         'loose':{'sigmaIEIE':0.0314,'dEtaSeed':0.00868,'dPhi':0.213,'invEinvP':0.14,'miss':1},
#         'medium':{'sigmaIEIE':0.0298,'dEtaSeed':0.00609,'dPhi':0.045,'invEinvP':0.13,'miss':1},
#         'tight':{'sigmaIEIE':0.0292,'dEtaSeed':0.00605,'dPhi':0.0394,'invEinvP':0.0129,'miss':1}
#     }
# } #80X
workingPoint_HoverE={
    'barrel':{'veto':.247,'loose':.236, 'medium' : .0859, 'tight': .0833},
    'endcap':{'veto':.0982,'loose':.0801, 'medium' : .0604, 'tight': .0543}
}   
        
workingPoint={
    'barrel':{
        'veto':{'sigmaIEIE':0.0126,'dEtaSeed':0.00463,'dPhi':0.148,'invEinvP':0.209,'miss':2},
        'loose':{'sigmaIEIE':0.0112,'dEtaSeed':0.00377,'dPhi':0.0884,'invEinvP':0.193,'miss':1},
        'medium':{'sigmaIEIE':0.0106,'dEtaSeed':0.0032,'dPhi':0.0547,'invEinvP':0.184,'miss':1},
        'tight':{'sigmaIEIE':0.0104,'dEtaSeed':0.00255,'dPhi':0.022,'invEinvP':0.159,'miss':1}
    },
    'endcap':{
        'veto':{'sigmaIEIE':0.0457,'dEtaSeed':0.00814,'dPhi':0.19,'invEinvP':0.132,'miss':3},
        'loose':{'sigmaIEIE':0.0425,'dEtaSeed':0.00674,'dPhi':0.169,'invEinvP':0.111,'miss':1},
        'medium':{'sigmaIEIE':0.0387,'dEtaSeed':0.00632,'dPhi':0.0394,'invEinvP':0.0721,'miss':1},
        'tight':{'sigmaIEIE':0.0353,'dEtaSeed':0.00501,'dPhi':0.0236,'invEinvP':0.0197,'miss':1}
    }
}


# -----------

inputFileNames=open(inputFileList, 'r')
for inputFileName in inputFileNames:
    if sample == "TCP":
        print "root://cmseos.fnal.gov//store/user/zhangj/events/ALP/RunIISummer17DR94Premix/"+inputFileName.replace("\n","")
        inputFileName="root://cmseos.fnal.gov//store/user/zhangj/events/ALP/RunIISummer17DR94Premix/"+inputFileName.replace("\n","")
#        print inputFileName.replace("\n","")
#        inputFileName=inputFileName.replace("\n","")
    if sample == "ZToEE":
        print "root://cmxrootd.fnal.gov/"+inputFileName.replace("\n","")
        inputFileName="root://cmsxrootd.fnal.gov/"+inputFileName.replace("\n","")
    #if sample == "DYJetsToLL":
    #print inputFileName.replace("\n","")
    inputFileName=inputFileName.replace("\n","")  
    f=ROOT.TFile.Open(inputFileName)
    print inputFileName
    if not f.IsZombie():
        events=Events(inputFileName)
    else:
        print "Can't Open File: "+inputFileName
        continue
    events=Events(inputFileName)    
 #   print events
    for event in events:
        #print event
        event.getByLabel(labelVertex, handleVertex)
        vertex=handleVertex.product()
        pv = vertex[0].position()
    
        event.getByLabel(labelElectron, handleElectron)
        electrons=handleElectron.product()
        event.getByLabel(labelRho, handleRho)
        if len(handleRho.product())>0:
            rho = handleRho.product()[0]
        else:
            rho = 0
    
        event.getByLabel(labelGenParticle, handleGenParticle)
        particles=handleGenParticle.product()
    
        event.getByLabel(labelGenInfo, handleGenInfo)
        geninfo=handleGenInfo.product()
        genweight=geninfo.weight()

        h['NEvent'].Fill(0.5, 1)
        h['NEvent'].Fill(1.5, genweight)  

        genElectrons = []
        #print genweight
        for particle in particles:
            if sample == "TCP": 
                if abs(particle.pdgId())==11 and particle.isDirectHardProcessTauDecayProductFinalState(): #TCP
                    genElectrons+=[particle]
            if sample == "ZToEE": 
                if abs(particle.pdgId())==11 and particle.isHardProcess(): #ZToEE
                    genElectrons+=[particle]
            if sample == "DYJetsToLL":
                if abs(particle.pdgId())==11 and particle.isHardProcess():
                    genElectrons+=[particle]
            if sample=="QCD":
                if abs(particle.pdgId())==11:  ##### For QCD we want to see our efficiency on Fake Electrons>>electrons not associated with a gen electron
                   genElectrons+=[particle]  
                   #print "test"
        selected_genElectrons = []
        #print genElectrons
        #print "TCP: "+ str( particle.isDirectHardProcessTauDecayProductFinalState())
        #print "ZToEE: " + str(particle.isHardProcess())
        for electron in genElectrons:
            if electron.pt()>7 and abs(electron.eta())<2.5: 
                selected_genElectrons+=[electron]
                h['gen_ePt'].Fill(electron.pt(), genweight)

        selected_genElectrons.sort(key=lambda x: x.pt(), reverse=True)
#        if len(selected_genElectrons) > 0: h['gen_ePt'].Fill(selected_genElectrons[0].pt(), genweight)

        if process == "var" and process != "eff":
            selected_recoElectrons = []   
            selected_electrons = [[],[],[],[]]
            selected_IsoElectrons = [[],[],[],[]]

            for electron in electrons: 
                h['NElectron'].Fill(0.5, genweight)
                if electron.pt()>7 and abs(electron.eta())<2.5:
                    #print "selected electron"
                    selected_recoElectrons+=[electron]
                    E_c = electron.superCluster().energy()
                    h['NElectron'].Fill(1.5, genweight)
                    h['ePt'].Fill(electron.pt(), genweight)
                    #print "eta = " +str(electron.eta())
                    matched = False
                    for genElec in selected_genElectrons:
                        lE = electron
                        lGenE = selected_genElectrons[0]
                        if lE.charge() == lGenE.charge():
                            genE = ROOT.TLorentzVector()
                            genE.SetPtEtaPhiM(lGenE.pt(), lGenE.eta(), lGenE.phi(), lGenE.mass())
                            e = ROOT.TLorentzVector()
                            e.SetPtEtaPhiM(lE.pt(), lE.eta(), lE.phi(), lE.mass())
                            if genE.DeltaR(e) < 0.1:
                                matched = True
                    if matched:
                        h['nmatched'].Fill(0.5,genweight)
                    if not matched:
                        h['nmatched'].Fill(1.5, genweight)
                    if electron.isEB():
                        if sample=='QCD' and not matched:
                            h['EB_sigmaIetaIeta'].Fill(electron.full5x5_sigmaIetaIeta(), genweight)
                            h['EB_dEtaSeed'].Fill(dEtaInSeed(electron), genweight)
                            h['EB_dPhiIn'].Fill(electron.deltaPhiSuperClusterTrackAtVtx(), genweight)
                            h['EB_hOverE'].Fill(electron.hadronicOverEm(), genweight)
                            h['hOverE'].Fill(electron.hadronicOverEm(), genweight)
                            h['Ec'].Fill(E_c, genweight)
                            h['EB_Ec'].Fill(E_c, genweight)
                            h['EB_hOverE_2D'].Fill(electron.hadronicOverEm(), 0.05+1.16/E_c+0.0324*rho/E_c, genweight)
                            if electron.full5x5_sigmaIetaIeta() < workingPoint['barrel']['loose']['sigmaIEIE'] \
                            and abs(dEtaInSeed(electron)) < workingPoint['barrel']['loose']['dEtaSeed'] \
                            and abs(electron.deltaPhiSuperClusterTrackAtVtx()) < workingPoint['barrel']['loose']['dPhi'] \
                            and GsfEleEInverseMinusPInverse(electron) < workingPoint['barrel']['loose']['invEinvP'] \
                            and GsfEleMissingHitsCut(electron) <= workingPoint['barrel']['loose']['miss']:
                                h['EB_hOverE_2D_pass'].Fill(electron.hadronicOverEm(), 0.05+1.16/E_c+0.0324*rho/E_c, genweight)
                                h['hOverE_pass'].Fill(electron.hadronicOverEm(),genweight)
                                h['hOverE_2D_pass'].Fill(electron.hadronicOverEm(), 0.05+1.16/E_c+0.0324*rho/E_c, genweight)
                                h['Ec_hOverE_2D_pass'].Fill(electron.hadronicOverEm(), E_c, genweight)
                                h['Ec_pass'].Fill(E_c,genweight)
                                h['EB_Ec_pass'].Fill(E_c, genweight)
                                h['EB_Ec_hOverE_2D_pass'].Fill(electron.hadronicOverEm(), E_c, genweight)
                            h['EB_relIso'].Fill(GsfEleEffAreaPFIsoCut(electron, rho), genweight)
                            h['EB_eInvpInv'].Fill(GsfEleEInverseMinusPInverse(electron), genweight)
                            h['EB_missHits'].Fill(GsfEleMissingHitsCut(electron), genweight)
                            h['EB_dz'].Fill(electron.gsfTrack().dz(pv), genweight)
                            h['EB_dxdy'].Fill(electron.gsfTrack().dxy(pv), genweight)
                        if sample=='TCP' and matched:
                            h['EB_sigmaIetaIeta'].Fill(electron.full5x5_sigmaIetaIeta(), genweight)
                            h['EB_dEtaSeed'].Fill(dEtaInSeed(electron), genweight)
                            h['EB_dPhiIn'].Fill(electron.deltaPhiSuperClusterTrackAtVtx(), genweight)
                            h['EB_hOverE'].Fill(electron.hadronicOverEm(), genweight)
                            h['hOverE'].Fill(electron.hadronicOverEm(), genweight)
                            h['Ec'].Fill(E_c, genweight)
                            h['EB_Ec'].Fill(E_c, genweight)
                            h['EB_hOverE_2D'].Fill(electron.hadronicOverEm(), 0.05+1.16/E_c+0.0324*rho/E_c, genweight)
                            if electron.full5x5_sigmaIetaIeta() < workingPoint['barrel']['loose']['sigmaIEIE'] \
                            and abs(dEtaInSeed(electron)) < workingPoint['barrel']['loose']['dEtaSeed'] \
                            and abs(electron.deltaPhiSuperClusterTrackAtVtx()) < workingPoint['barrel']['loose']['dPhi'] \
                            and GsfEleEInverseMinusPInverse(electron) < workingPoint['barrel']['loose']['invEinvP'] \
                            and GsfEleMissingHitsCut(electron) <= workingPoint['barrel']['loose']['miss']:
                                h['EB_hOverE_2D_pass'].Fill(electron.hadronicOverEm(), 0.05+1.16/E_c+0.0324*rho/E_c, genweight)
                                h['hOverE_pass'].Fill(electron.hadronicOverEm(),genweight)
                                h['hOverE_2D_pass'].Fill(electron.hadronicOverEm(), 0.05+1.16/E_c+0.0324*rho/E_c, genweight)
                                h['Ec_hOverE_2D_pass'].Fill(electron.hadronicOverEm(), E_c, genweight)
                                h['Ec_pass'].Fill(E_c,genweight)
                                h['EB_Ec_pass'].Fill(E_c, genweight)
                                h['EB_Ec_hOverE_2D_pass'].Fill(electron.hadronicOverEm(), E_c, genweight)
                            h['EB_relIso'].Fill(GsfEleEffAreaPFIsoCut(electron, rho), genweight)
                            h['EB_eInvpInv'].Fill(GsfEleEInverseMinusPInverse(electron), genweight)
                            h['EB_missHits'].Fill(GsfEleMissingHitsCut(electron), genweight)
                            h['EB_dz'].Fill(electron.gsfTrack().dz(pv), genweight)
                            h['EB_dxdy'].Fill(electron.gsfTrack().dxy(pv), genweight)

                    if electron.isEE():
                        if sample == 'QCD' and not matched:
                            h['EE_sigmaIetaIeta'].Fill(electron.full5x5_sigmaIetaIeta(), genweight)
                            h['EE_dEtaSeed'].Fill(dEtaInSeed(electron), genweight)
                            h['EE_dPhiIn'].Fill(electron.deltaPhiSuperClusterTrackAtVtx(), genweight)
                            h['EE_hOverE'].Fill(electron.hadronicOverEm(), genweight)
                            h['hOverE'].Fill(electron.hadronicOverEm(), genweight)
                            h['Ec'].Fill(E_c, genweight)
                            h['EE_Ec'].Fill(E_c, genweight)
                            h['EE_hOverE_2D'].Fill(electron.hadronicOverEm(), 0.0441+2.54/E_c+0.183*rho/E_c, genweight)
                            if electron.full5x5_sigmaIetaIeta() < workingPoint['endcap']['loose']['sigmaIEIE'] \
                            and abs(dEtaInSeed(electron)) < workingPoint['endcap']['loose']['dEtaSeed'] \
                            and abs(electron.deltaPhiSuperClusterTrackAtVtx()) < workingPoint['endcap']['loose']['dPhi'] \
                            and GsfEleEInverseMinusPInverse(electron) < workingPoint['endcap']['loose']['invEinvP'] \
                            and GsfEleMissingHitsCut(electron) <= workingPoint['endcap']['loose']['miss']:
                                h['EE_hOverE_2D_pass'].Fill(electron.hadronicOverEm(), 0.0441+2.54/E_c+0.183*rho/E_c, genweight)
                                h['hOverE_pass'].Fill(electron.hadronicOverEm(),genweight)
                                h['hOverE_2D_pass'].Fill(electron.hadronicOverEm(), 0.0441+2.54/E_c+0.183*rho/E_c/E_c, genweight)
                                h['Ec_hOverE_2D_pass'].Fill(electron.hadronicOverEm(), E_c, genweight)
                                h['Ec_pass'].Fill(E_c,genweight)
                                h['EE_Ec_pass'].Fill(E_c, genweight)
                                h['EE_Ec_hOverE_2D_pass'].Fill(electron.hadronicOverEm(), E_c, genweight)
                            h['EE_relIso'].Fill(GsfEleEffAreaPFIsoCut(electron, rho), genweight)
                            h['EE_eInvpInv'].Fill(GsfEleEInverseMinusPInverse(electron), genweight)
                            h['EE_missHits'].Fill(GsfEleMissingHitsCut(electron), genweight)
                            h['EE_dz'].Fill(electron.gsfTrack().dz(pv), genweight)
                            h['EE_dxdy'].Fill(electron.gsfTrack().dxy(pv), genweight)
                        if sample == 'TCP' and matched:
                            h['EE_sigmaIetaIeta'].Fill(electron.full5x5_sigmaIetaIeta(), genweight)
                            h['EE_dEtaSeed'].Fill(dEtaInSeed(electron), genweight)
                            h['EE_dPhiIn'].Fill(electron.deltaPhiSuperClusterTrackAtVtx(), genweight)
                            h['EE_hOverE'].Fill(electron.hadronicOverEm(), genweight)
                            h['hOverE'].Fill(electron.hadronicOverEm(), genweight)
                            h['Ec'].Fill(E_c, genweight)
                            h['EE_Ec'].Fill(E_c, genweight)
                            h['EE_hOverE_2D'].Fill(electron.hadronicOverEm(),
                                                       0.0441 + 2.54 / E_c + 0.183 * rho / E_c, genweight)
                            if electron.full5x5_sigmaIetaIeta() < workingPoint['endcap']['loose']['sigmaIEIE'] \
                                and abs(dEtaInSeed(electron)) < workingPoint['endcap']['loose']['dEtaSeed'] \
                                and abs(electron.deltaPhiSuperClusterTrackAtVtx()) < \
                                workingPoint['endcap']['loose']['dPhi'] \
                                and GsfEleEInverseMinusPInverse(electron) < workingPoint['endcap']['loose'][
                                'invEinvP'] \
                                and GsfEleMissingHitsCut(electron) <= workingPoint['endcap']['loose']['miss']:
                                h['EE_hOverE_2D_pass'].Fill(electron.hadronicOverEm(), 0.0441 + 2.54 / E_c + 0.183 * rho / E_c, genweight)
                                h['hOverE_pass'].Fill(electron.hadronicOverEm(), genweight)
                                h['hOverE_2D_pass'].Fill(electron.hadronicOverEm(), 0.0441 + 2.54 / E_c + 0.183 * rho / E_c / E_c, genweight)
                                h['Ec_hOverE_2D_pass'].Fill(electron.hadronicOverEm(), E_c, genweight)
                                h['Ec_pass'].Fill(E_c, genweight)
                                h['EE_Ec_pass'].Fill(E_c, genweight)
                                h['EE_Ec_hOverE_2D_pass'].Fill(electron.hadronicOverEm(), E_c, genweight)
                            h['EE_relIso'].Fill(GsfEleEffAreaPFIsoCut(electron, rho), genweight)
                            h['EE_eInvpInv'].Fill(GsfEleEInverseMinusPInverse(electron), genweight)
                            h['EE_missHits'].Fill(GsfEleMissingHitsCut(electron), genweight)
                            h['EE_dz'].Fill(electron.gsfTrack().dz(pv), genweight)
                            h['EE_dxdy'].Fill(electron.gsfTrack().dxy(pv), genweight)


        if process == "eff" and process!= "var":

            selected_recoElectrons = []   
            selected_electrons = [[],[],[],[]]
            selected_IsoElectrons = [[],[],[],[]]
            selected_electrons_iso_EB=[[],[],[],[]]
            selected_electrons_iso_EE=[[],[],[],[]] 
            selected_electrons_old_EE=[[],[],[],[]]     
            selected_electrons_old_EB=[[],[],[],[]]
            selected_electrons_old=[[],[],[],[]]
            selected_electrons_iso=[[],[],[],[]]
            selected_electrons_EE = [[],[],[],[]] 
            selected_electrons_EB = [[],[],[],[]]     
            for electron in electrons: 
                h['NElectron'].Fill(0.5, genweight)
                if electron.pt()>7 and abs(electron.eta())<2.5:
                    selected_recoElectrons+=[electron]
                    E_c = electron.superCluster().energy()
                    h['NElectron'].Fill(1.5, genweight)
                    h['ePt'].Fill(electron.pt(), genweight)
                    if electron.isEB():
                        loc = "barrel"
                        for wp in list(workingPoint[loc]):
                            if wp == "veto":
                                workingPoint[loc][wp]['hOvE'] = 0.05+1.16/E_c+0.0324*rho/E_c
                                workingPoint[loc][wp]['relIso'] = 0.198+0.506/electron.pt()
                                n = 0.5
                                w = 0
                            if wp == "loose":
                                workingPoint[loc][wp]['hOvE'] = 0.05+1.16/E_c+0.0324*rho/E_c
                                workingPoint[loc][wp]['relIso'] = 0.112+0.506/electron.pt()
                                n = 1.5
                                w = 1
                            if wp == "medium":
                                workingPoint[loc][wp]['hOvE'] = 0.046+1.16/E_c+0.0324*rho/E_c
                                workingPoint[loc][wp]['relIso'] = 0.0478+0.506/electron.pt()
                                n = 2.5
                                w = 2
                            if wp == "tight":
                                workingPoint[loc][wp]['hOvE'] = 0.026+1.15/E_c+0.0324*rho/E_c
                                workingPoint[loc][wp]['relIso'] = 0.0287+0.506/electron.pt()
                                n = 3.5
                                w = 3

                            #if electron.full5x5_sigmaIetaIeta() < workingPoint[loc][wp]['sigmaIEIE'] \
                            #and abs(dEtaInSeed(electron)) < workingPoint[loc][wp]['dEtaSeed'] \
                            #and abs(electron.deltaPhiSuperClusterTrackAtVtx()) < workingPoint[loc][wp]['dPhi'] \
                            #and GsfEleEInverseMinusPInverse(electron) < workingPoint[loc][wp]['invEinvP'] \
                            #and GsfEleMissingHitsCut(electron) <= workingPoint[loc][wp]['miss'] \
                            #and electron.hadronicOverEm() < workingPoint[loc][wp]['hOvE'] \
                            #and electron.passConversionVeto() \
                            #and abs(electron.gsfTrack().dz(pv))<0.1 \
                            #and abs(electron.gsfTrack().dxy(pv))<0.05:
                                #selected_electrons[w]+=[electron]
                                #h['NElectronID'].Fill(n, genweight)
                            if electron.full5x5_sigmaIetaIeta() < workingPoint[loc][wp]['sigmaIEIE'] \
                            and abs(dEtaInSeed(electron)) < workingPoint[loc][wp]['dEtaSeed'] \
                            and abs(electron.deltaPhiSuperClusterTrackAtVtx()) < workingPoint[loc][wp]['dPhi'] \
                            and GsfEleEInverseMinusPInverse(electron) < workingPoint[loc][wp]['invEinvP'] \
                            and GsfEleMissingHitsCut(electron) <= workingPoint[loc][wp]['miss'] \
                            and electron.hadronicOverEm() < workingPoint[loc][wp]['hOvE'] \
                            and electron.passConversionVeto() \
                            and abs(electron.gsfTrack().dz(pv))<0.2 \
                            and abs(electron.gsfTrack().dxy(pv))<0.1:
                                selected_electrons[w]+=[electron]      
                                selected_electrons_EB[w]+=[electron]
                                h['NElectronID'].Fill(n, genweight)
                                ############################## Find Efficiency of the HoverE cut ONLY#############################################

                            if electron.full5x5_sigmaIetaIeta() < workingPoint[loc][wp]['sigmaIEIE'] \
                            and abs(dEtaInSeed(electron)) < workingPoint[loc][wp]['dEtaSeed'] \
                            and GsfEleEInverseMinusPInverse(electron) < workingPoint[loc][wp]['invEinvP'] \
                            and GsfEleMissingHitsCut(electron) <= workingPoint[loc][wp]['miss'] \
                            and abs(electron.deltaPhiSuperClusterTrackAtVtx()) < workingPoint[loc][wp]['dPhi']\
                            and abs(electron.gsfTrack().dz(pv))<0.2 \
                            and electron.passConversionVeto() \
                            and abs(electron.gsfTrack().dxy(pv))<0.1:
                                selected_electrons_iso_EB[w]+=[electron]
                                selected_electrons_iso[w]+=[electron]                 
                                h['NElectronIso'].Fill(n, genweight)
                                h['EB_nElectronIso'].Fill(n,genweight)
                                if electron.hadronicOverEm() < workingPoint_HoverE[loc][wp]:
                                    selected_electrons_old_EB[w]+=[electron]
                                    selected_electrons_old[w]+=[electron]
                                    h['EB_nOld'].Fill(n, genweight)
                                    h['NOld'].Fill(n, genweight)

                    if electron.isEE():
                        loc = "endcap"
                        for wp in list(workingPoint[loc]):
                            if wp == "veto":
                                workingPoint[loc][wp]['hOvE'] = 0.05+2.54/E_c+0.183*rho/E_c
                                workingPoint[loc][wp]['relIso'] = 0.203+0.963/electron.pt()
                                n = 0.5
                                w = 0
                            if wp == "loose": 
                                workingPoint[loc][wp]['hOvE'] = 0.0441+2.54/E_c+0.183*rho/E_c
                                workingPoint[loc][wp]['relIso'] = 0.108+0.963/electron.pt()
                                n = 1.5
                                w = 1
                            if wp == "medium": 
                                workingPoint[loc][wp]['hOvE'] = 0.0275+2.52/E_c+0.183*rho/E_c
                                workingPoint[loc][wp]['relIso'] = 0.0658+0.963/electron.pt()
                                n = 2.5
                                w = 2
                            if wp == "tight": 
                                workingPoint[loc][wp]['hOvE'] = 0.0188+2.06/E_c+0.183*rho/E_c
                                workingPoint[loc][wp]['relIso'] = 0.0445+0.963/electron.pt()
                                n = 3.5
                                w = 3

                            if electron.full5x5_sigmaIetaIeta() < workingPoint[loc][wp]['sigmaIEIE'] \
                            and abs(dEtaInSeed(electron)) < workingPoint[loc][wp]['dEtaSeed'] \
                            and abs(electron.deltaPhiSuperClusterTrackAtVtx()) < workingPoint[loc][wp]['dPhi'] \
                            and GsfEleEInverseMinusPInverse(electron) < workingPoint[loc][wp]['invEinvP'] \
                            and GsfEleMissingHitsCut(electron) <= workingPoint[loc][wp]['miss'] \
                            and electron.hadronicOverEm() < workingPoint[loc][wp]['hOvE'] \
                            and electron.passConversionVeto() \
                            and abs(electron.gsfTrack().dz(pv))<0.2 \
                            and abs(electron.gsfTrack().dxy(pv))<0.1:
                                selected_electrons[w]+=[electron]
                                selected_electrons_EE[w]+=[electron]
                                h['NElectronID'].Fill(n, genweight)
                            ############################## Find Efficiency of the HoverE cut ONLY#############################################    
                            if electron.full5x5_sigmaIetaIeta() < workingPoint[loc][wp]['sigmaIEIE'] \
                            and abs(electron.deltaPhiSuperClusterTrackAtVtx()) < workingPoint[loc][wp]['dPhi'] \
                            and abs(dEtaInSeed(electron)) < workingPoint[loc][wp]['dEtaSeed'] \
                            and GsfEleEInverseMinusPInverse(electron) < workingPoint[loc][wp]['invEinvP'] \
                            and GsfEleMissingHitsCut(electron) <= workingPoint[loc][wp]['miss'] \
                            and abs(electron.gsfTrack().dz(pv))<0.2 \
                            and electron.passConversionVeto() \
                            and abs(electron.gsfTrack().dxy(pv))<0.1: 
                                selected_electrons_iso_EE[w]+=[electron]
                                selected_electrons_iso[w]+=[electron]
                                h['EE_nElectronIso'].Fill(n,genweight)
                                h['NElectronIso'].Fill(n, genweight)
                                if electron.hadronicOverEm() < workingPoint_HoverE[loc][wp]:
                                    h['EE_nOld'].Fill(n, genweight)
                                    h['NOld'].Fill(n, genweight)
                                    selected_electrons_old[w]+=[electron]
                                    if wp==1:
                                        h['NOld_loose_ePt'].Fill(electron.pt(), genweight)


            for i in range(4):
                selected_electrons[i].sort(key=lambda x: x.pt(), reverse=True)
                selected_IsoElectrons[i].sort(key=lambda x: x.pt(), reverse=True)
                selected_electrons_iso[i].sort(key=lambda x: x.pt(), reverse=True)
                selected_electrons_old[i].sort(key=lambda x: x.pt(), reverse=True)
            if len(selected_electrons[0])>0: h['ePt_veto'].Fill(selected_electrons[0][0].pt(), genweight)
            if len(selected_electrons[1])>0: h['ePt_loose'].Fill(selected_electrons[1][0].pt(), genweight)
            if len(selected_electrons[2])>0: h['ePt_medium'].Fill(selected_electrons[2][0].pt(), genweight)
            if len(selected_electrons[3])>0: h['ePt_tight'].Fill(selected_electrons[3][0].pt(), genweight)
            if len(selected_electrons_iso_EE[1])>0: h['EE_looseIso_ePt'].Fill(selected_electrons_iso_EE[1][0].pt(), genweight)
            if len(selected_electrons_iso_EB[1]) >0:  h['EB_looseIso_ePt'].Fill(selected_electrons_iso_EB[1][0].pt(), genweight)
            if len(selected_electrons_iso[1]) > 0:  h['NElectron_looseIso_ePt'].Fill(selected_electrons_iso[1][0].pt(), genweight)
            selected_recoElectrons.sort(key=lambda x: x.pt(), reverse=True)
            #if len(selected_recoElectrons)>0: h['ePt'].Fill(electron.pt(), genweight)
############################## Find Efficiency of the HoverE cut ONLY#############################################
#######Seperate into EB and EE, start with lose only and ad more working points as needed#################

            for electron in selected_recoElectrons:
                matched = False
                for gen in selected_genElectrons:
                    lE = electron
                    lGenE = selected_genElectrons[0]
                    if lE.charge() == lGenE.charge():
                        genE = ROOT.TLorentzVector()
                        genE.SetPtEtaPhiM(lGenE.pt(), lGenE.eta(), lGenE.phi(), lGenE.mass())
                        e = ROOT.TLorentzVector()
                        e.SetPtEtaPhiM(lE.pt(), lE.eta(), lE.phi(), lE.mass())
                        if genE.DeltaR(e) < 0.1:
                            matched=True
                if sample=='QCD' and not matched:
                    for i in range (0,4):
                        h['NElectron_matched'].Fill(i+0.5, genweight)
                if sample == 'TCP' and matched:
                    for i in range (0,4):
                        h['NElectron_matched'].Fill(i+0.5, genweight)
            for electron in selected_electrons[0]:
                lE = electron
                matched=False
                for gen in selected_genElectrons:
                    lGenE = gen
                    if lE.charge() == lGenE.charge():
                        genE = ROOT.TLorentzVector()
                        genE.SetPtEtaPhiM(lGenE.pt(), lGenE.eta(), lGenE.phi(), lGenE.mass())
                        e = ROOT.TLorentzVector()
                        e.SetPtEtaPhiM(lE.pt(), lE.eta(), lE.phi(), lE.mass())
                        if genE.DeltaR(e) < 0.1:
                            matched = True
                if sample=='QCD' and not matched:
                    h['matched_ePt_veto'].Fill(lE.pt(), genweight)
                    h['NElectronID_matched'].Fill(0.5, genweight)
                if sample == 'TCP' and matched:
                    h['matched_ePt_veto'].Fill(lE.pt(), genweight)
                    h['NElectronID_matched'].Fill(0.5, genweight)
            for electron in selected_electrons[1]:
                lE = electron
                matched = False
                for gen in selected_genElectrons:
                    lGenE = gen
                    if lE.charge() == lGenE.charge():
                            genE = ROOT.TLorentzVector()
                            genE.SetPtEtaPhiM(lGenE.pt(), lGenE.eta(), lGenE.phi(), lGenE.mass())
                            e = ROOT.TLorentzVector()
                            e.SetPtEtaPhiM(lE.pt(), lE.eta(), lE.phi(), lE.mass())
                            if genE.DeltaR(e) < 0.1:
                                matched = True
                if sample=='QCD' and not matched:
                    h['matched_ePt_loose'].Fill(lE.pt(), genweight)
                    h['NElectronID_matched'].Fill(1.5, genweight)
                if sample == 'TCP' and matched:
                    h['matched_ePt_loose'].Fill(lE.pt(), genweight)
                    h['NElectronID_matched'].Fill(1.5, genweight)
            for electron in selected_electrons[2]:
                lE = electron
                matched = False
                for gen in selected_genElectrons:
                    lGenE = gen
                    if lE.charge() == lGenE.charge():
                        genE = ROOT.TLorentzVector()
                        genE.SetPtEtaPhiM(lGenE.pt(), lGenE.eta(), lGenE.phi(), lGenE.mass())
                        e = ROOT.TLorentzVector()
                        e.SetPtEtaPhiM(lE.pt(), lE.eta(), lE.phi(), lE.mass())
                        if genE.DeltaR(e) < 0.1:
                            matched = True
                if sample=='QCD' and not matched:
                    h['matched_ePt_medium'].Fill(lE.pt(), genweight)
                    h['NElectronID_matched'].Fill(2.5, genweight)
                if sample == 'TCP' and matched:
                    h['matched_ePt_medium'].Fill(lE.pt(), genweight)
                    h['NElectronID_matched'].Fill(2.5, genweight)
            for electron in selected_electrons[3]:
                lE = electron
                matched = False
                for gen in selected_genElectrons:
                    lGenE = gen
                    if lE.charge() == lGenE.charge():
                        genE = ROOT.TLorentzVector()
                        genE.SetPtEtaPhiM(lGenE.pt(), lGenE.eta(), lGenE.phi(), lGenE.mass())
                        e = ROOT.TLorentzVector()
                        e.SetPtEtaPhiM(lE.pt(), lE.eta(), lE.phi(), lE.mass())
                        if genE.DeltaR(e) < 0.1:
                            matched = True
                if sample=='QCD' and not matched:
                    h['matched_ePt_tight'].Fill(lE.pt(), genweight)
                    h['NElectronID_matched'].Fill(3.5, genweight)
                if sample == 'TCP' and matched:
                    h['matched_ePt_tight'].Fill(lE.pt(), genweight)
                    h['NElectronID_matched'].Fill(3.5, genweight)

#################Do matching using only the OLD Hover Ecut#####################
            for electron in selected_electrons_old[0]:
                lE = electron
                matched=False
                for gen in selected_genElectrons:
                    lGenE = gen
                    if lE.charge() == lGenE.charge():
                        genE = ROOT.TLorentzVector()
                        genE.SetPtEtaPhiM(lGenE.pt(), lGenE.eta(), lGenE.phi(), lGenE.mass())
                        e = ROOT.TLorentzVector()
                        e.SetPtEtaPhiM(lE.pt(), lE.eta(), lE.phi(), lE.mass())
                        if genE.DeltaR(e) < 0.1:
                            matched = True
                if sample=='QCD' and not matched:
                    h['matched_NOld'].Fill(0.5, genweight)
                if sample == 'TCP' and matched:
                    h['matched_NOld'].Fill(0.5, genweight)
            for electron in selected_electrons_old[1]:
                lE = electron
                matched = False
                for gen in selected_genElectrons:
                    lGenE = gen
                    if lE.charge() == lGenE.charge():
                            genE = ROOT.TLorentzVector()
                            genE.SetPtEtaPhiM(lGenE.pt(), lGenE.eta(), lGenE.phi(), lGenE.mass())
                            e = ROOT.TLorentzVector()
                            e.SetPtEtaPhiM(lE.pt(), lE.eta(), lE.phi(), lE.mass())
                            if genE.DeltaR(e) < 0.1:
                                matched = True
                if sample=='QCD' and not matched:
                    h['matched_NOld_loose_ePt'].Fill(lE.pt(), genweight)
                    h['matched_NOld'].Fill(1.5, genweight)
                if sample == 'TCP' and matched:
                    h['matched_NOld_loose_ePt'].Fill(lE.pt(), genweight)
                    h['matched_NOld'].Fill(1.5, genweight)
            for electron in selected_electrons_old[2]:
                lE = electron
                matched = False
                for gen in selected_genElectrons:
                    lGenE = gen
                    if lE.charge() == lGenE.charge():
                        genE = ROOT.TLorentzVector()
                        genE.SetPtEtaPhiM(lGenE.pt(), lGenE.eta(), lGenE.phi(), lGenE.mass())
                        e = ROOT.TLorentzVector()
                        e.SetPtEtaPhiM(lE.pt(), lE.eta(), lE.phi(), lE.mass())
                        if genE.DeltaR(e) < 0.1:
                            matched = True
                if sample=='QCD' and not matched:
                    h['matched_NOld'].Fill(2.5, genweight)
                if sample == 'TCP' and matched:
                    h['matched_NOld'].Fill(2.5, genweight)
            for electron in selected_electrons_old[3]:
                lE = electron
                matched = False
                for gen in selected_genElectrons:
                    lGenE = gen
                    if lE.charge() == lGenE.charge():
                        genE = ROOT.TLorentzVector()
                        genE.SetPtEtaPhiM(lGenE.pt(), lGenE.eta(), lGenE.phi(), lGenE.mass())
                        e = ROOT.TLorentzVector()
                        e.SetPtEtaPhiM(lE.pt(), lE.eta(), lE.phi(), lE.mass())
                        if genE.DeltaR(e) < 0.1:
                            matched = True
                if sample=='QCD' and not matched:
                    h['matched_NOld'].Fill(3.5, genweight)
                if sample == 'TCP' and matched:
                    h['matched_NOld'].Fill(3.5, genweight)


            #######################Matching the electrons the pass cut excluding HoverE
            for electron in selected_electrons_iso[0]:
                lE = electron
                matched=False
                for gen in selected_genElectrons:
                    lGenE = gen
                    if lE.charge() == lGenE.charge():
                        genE = ROOT.TLorentzVector()
                        genE.SetPtEtaPhiM(lGenE.pt(), lGenE.eta(), lGenE.phi(), lGenE.mass())
                        e = ROOT.TLorentzVector()
                        e.SetPtEtaPhiM(lE.pt(), lE.eta(), lE.phi(), lE.mass())
                        if genE.DeltaR(e) < 0.1:
                            matched = True
                if sample=='QCD' and not matched:
                    h['matched_nElectronIso'].Fill(0.5, genweight)
                if sample == 'TCP' and matched:
                    h['matched_nElectronIso'].Fill(0.5, genweight)
            for electron in selected_electrons_iso[1]:
                lE = electron
                matched = False
                for gen in selected_genElectrons:
                    lGenE = gen
                    if lE.charge() == lGenE.charge():
                            genE = ROOT.TLorentzVector()
                            genE.SetPtEtaPhiM(lGenE.pt(), lGenE.eta(), lGenE.phi(), lGenE.mass())
                            e = ROOT.TLorentzVector()
                            e.SetPtEtaPhiM(lE.pt(), lE.eta(), lE.phi(), lE.mass())
                            if genE.DeltaR(e) < 0.1:
                                matched = True
                if sample=='QCD' and not matched:
                    h['matched_nElectron_looseIso_ePt'].Fill(lE.pt(), genweight)
                    h['matched_nElectronIso'].Fill(1.5, genweight)
                if sample =='TCP' and matched:
                    h['matched_nElectron_looseIso_ePt'].Fill(lE.pt(), genweight)
                    h['matched_nElectronIso'].Fill(1.5, genweight)
            for electron in selected_electrons_iso[2]:
                lE = electron
                matched = False
                for gen in selected_genElectrons:
                    lGenE = gen
                    if lE.charge() == lGenE.charge():
                        genE = ROOT.TLorentzVector()
                        genE.SetPtEtaPhiM(lGenE.pt(), lGenE.eta(), lGenE.phi(), lGenE.mass())
                        e = ROOT.TLorentzVector()
                        e.SetPtEtaPhiM(lE.pt(), lE.eta(), lE.phi(), lE.mass())
                        if genE.DeltaR(e) < 0.1:
                            matched = True
                if sample=='QCD' and not matched:
                    h['matched_nElectronIso'].Fill(2.5, genweight)
                if sample == 'TCP' and matched:
                    h['matched_nElectronIso'].Fill(2.5, genweight)
            for electron in selected_electrons_iso[3]:
                lE = electron
                matched = False
                for gen in selected_genElectrons:
                    lGenE = gen
                    if lE.charge() == lGenE.charge():
                        genE = ROOT.TLorentzVector()
                        genE.SetPtEtaPhiM(lGenE.pt(), lGenE.eta(), lGenE.phi(), lGenE.mass())
                        e = ROOT.TLorentzVector()
                        e.SetPtEtaPhiM(lE.pt(), lE.eta(), lE.phi(), lE.mass())
                        if genE.DeltaR(e) < 0.1:
                            matched = True
                if sample=='QCD' and not matched:
                    h['matched_nElectronIso'].Fill(3.5, genweight)
                if sample =='TCP' and matched:
                    h['matched_nElectronIso'].Fill(3.5, genweight)

                            #do your fucking job nick. You're better than this

out.cd()

for key in h.keys():
    h[key].Write()

out.Close()
