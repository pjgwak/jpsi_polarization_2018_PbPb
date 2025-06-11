import os, ROOT
from Plot import Plot1D
from ROOT import TFile, TCanvas, TLegend, kGreen, kBlack, TCanvas, TH1

# set batch mode
ROOT.gROOT.SetBatch(True)

# set CMS plot style
rootlogon_path = "../../../rootlogon.C"
if os.path.exists(rootlogon_path):
    print('Use rootlogon.C')
    ROOT.gROOT.ProcessLine(f".X {rootlogon_path}")
else:
    print(f"Not exist: {rootlogon_path}")

# ====================== #
# ===== Data - Fwd ===== #
# ====================== #
# get data TTree
f = TFile.Open('OniaFlowSkim.Data2023.AlgDefault.root')
myTree = f.Get('myTree')

# rapidity common
yRegion = 'Fwd'
yName = '_pT3_6p5_y1p6_2p4'
legEntry1 = '1.6 < |y| < 2.4'
legEntry2 = '3 < #it{p}_{T} < 50'
cutBase = 'pt > 3 && pt < 50 && fabs(y) > 1.6 && fabs(y) < 2.4'


# ===== Fwd mass ===== 
baseName = f'Data2023_{yRegion}_Mass'
myCanvas = TCanvas(baseName+yName, '', 800, 600)
leg = TLegend(0.8, 0.8, 0.9, 0.9)
legEntry = [legEntry1, legEntry2, '', ''] # Up to 4 lines
myHist = ROOT.TH1D('h_'+baseName, '', 30, 2.6, 3.5)
histTitles = ';mass (GeV/c^{2});# of dimuon'
myCut = f'{cutBase}' # Functionize - this part

fwdMass = Plot1D(tree=myTree, var='mass', cut=myCut, hist=myHist, histTitles=histTitles, leg=leg, legEntry=legEntry, canvas=myCanvas)
fwdMass.draw()
fwdMass.save(fPath='figs/', isPng=True, isPdf=False)

# ===== Fwd pT ===== 
baseName = f'Data2023_{yRegion}_Pt'
myCanvas = TCanvas(baseName+yName, '', 800, 600)
leg = TLegend(0.8, 0.8, 0.9, 0.9)
legEntry = [legEntry1, legEntry2, '', ''] # Up to 4 lines
myHist = ROOT.TH1D('h_'+baseName, '', 100, 0, 50)
histTitles = ';#it{p}_{T} (GeV/c);# of dimuon'
myCut = f'{cutBase}'

fwdPt = Plot1D(tree=myTree, var='pt', cut=myCut, hist=myHist, histTitles=histTitles, leg=leg, legEntry=legEntry, canvas=myCanvas)
fwdPt.draw()
fwdPt.save(fPath='figs/', isPng=True, isPdf=False)

# ===== Fwd phi ===== 
baseName = f'Data2023_{yRegion}_Phi'
myCanvas = TCanvas(baseName+yName, '', 800, 600)
leg = TLegend(0.8, 0.8, 0.9, 0.9)
legEntry = [legEntry1, legEntry2, '', ''] # Up to 4 lines
myHist = ROOT.TH1D('h_'+baseName, '', 30, -3.141592, 3.141592)
histTitles = ';#varphi_{1} (radian);# of dimuon'
myCut = f'{cutBase}'

fwdPhi = Plot1D(tree=myTree, var='phi1', cut=myCut, hist=myHist, histTitles=histTitles, leg=leg, legEntry=legEntry, canvas=myCanvas)
fwdPhi.draw()
fwdPhi.save(fPath='figs/', isPng=True, isPdf=False)

# ===== Fwd phi HX ===== 
baseName = f'Data2023_{yRegion}_PhiHX'
myCanvas = TCanvas(baseName+yName, '', 800, 600)
leg = TLegend(0.8, 0.8, 0.9, 0.9)
legEntry = [legEntry1, legEntry2, '', ''] # Up to 4 lines
myHist = ROOT.TH1D('h_'+baseName, '', 30, -3.141592, 3.141592)
histTitles = ';#varphi_{HX} (radian);# of dimuon'
myCut = f'{cutBase}'

fwdPhiHX = Plot1D(tree=myTree, var='phiHX', cut=myCut, hist=myHist, histTitles=histTitles, leg=leg, legEntry=legEntry, canvas=myCanvas)
fwdPhiHX.draw()
fwdPhiHX.save(fPath='figs/', isPng=True, isPdf=False)

# ===== Fwd phi CS ===== 
baseName = f'Data2023_{yRegion}_PhiCS'
myCanvas = TCanvas(baseName+yName, '', 800, 600)
leg = TLegend(0.8, 0.8, 0.9, 0.9)
legEntry = [legEntry1, legEntry2, '', ''] # Up to 4 lines
myHist = ROOT.TH1D('h_'+baseName, '', 30, -3.141592, 3.141592)
histTitles = ';#varphi_{CS} (radian);# of dimuon'
myCut = f'{cutBase}'

fwdPhiCS = Plot1D(tree=myTree, var='phiCS', cut=myCut, hist=myHist, histTitles=histTitles, leg=leg, legEntry=legEntry, canvas=myCanvas)
fwdPhiCS.draw()
fwdPhiCS.save(fPath='figs/', isPng=True, isPdf=False)

# ===== Fwd cos HX ===== 
baseName = f'Data2023_{yRegion}_CosHX'
myCanvas = TCanvas(baseName+yName, '', 800, 600)
leg = TLegend(0.8, 0.8, 0.9, 0.9)
legEntry = [legEntry1, legEntry2, '', ''] # Up to 4 lines
myHist = ROOT.TH1D('h_'+baseName, '', 30, -1, 1)
histTitles = ';cos#theta_{HX};# of dimuon'
myCut = f'{cutBase}'

fwdCosHX = Plot1D(tree=myTree, var='cosHX', cut=myCut, hist=myHist, histTitles=histTitles, leg=leg, legEntry=legEntry, canvas=myCanvas)
fwdCosHX.draw()
fwdCosHX.save(fPath='figs/', isPng=True, isPdf=False)

# ===== Fwd cos CS ===== 
baseName = f'Data2023_{yRegion}_CosCS'
myCanvas = TCanvas(baseName+yName, '', 800, 600)
leg = TLegend(0.8, 0.8, 0.9, 0.9)
legEntry = [legEntry1, legEntry2, '', ''] # Up to 4 lines
myHist = ROOT.TH1D('h_'+baseName, '', 30, -1, 1)
histTitles = ';cos#theta_{CS};# of dimuon'
myCut = f'{cutBase}'

fwdCosCS = Plot1D(tree=myTree, var='cosCS', cut=myCut, hist=myHist, histTitles=histTitles, leg=leg, legEntry=legEntry, canvas=myCanvas)
fwdCosCS.draw()
fwdCosCS.save(fPath='figs/', isPng=True, isPdf=False)


# ====================== #
# ===== Data - Mid ===== #
# ====================== #
# rapidity common
yRegion = 'Mid'
yName = '_pT6p5_50_y1p6' # Functionize - this part
legEntry1 = '|y| < 1.6'
legEntry2 = '6.5 < #it{p}_{T} < 50'
cutBase = 'pt > 6.5 && pt < 50 && fabs(y) > 0 && fabs(y) < 1.6'


# ===== Mid mass ===== 
baseName = f'Data2023_{yRegion}_Mass'
myCanvas = TCanvas(baseName+yName, '', 800, 600)
leg = TLegend(0.8, 0.8, 0.9, 0.9)
legEntry = [legEntry1, legEntry2, '', ''] # Up to 4 lines
myHist = ROOT.TH1D('h_'+baseName, '', 30, 2.6, 3.5)
histTitles = ';mass (GeV/c^{2});# of dimuon'
myCut = f'{cutBase}' # Functionize - this part

midMass = Plot1D(tree=myTree, var='mass', cut=myCut, hist=myHist, histTitles=histTitles, leg=leg, legEntry=legEntry, canvas=myCanvas)
midMass.draw()
midMass.save(fPath='figs/', isPng=True, isPdf=False)

# ===== Mid pT ===== 
baseName = f'Data2023_{yRegion}_Pt'
myCanvas = TCanvas(baseName+yName, '', 800, 600)
leg = TLegend(0.8, 0.8, 0.9, 0.9)
legEntry = [legEntry1, legEntry2, '', ''] # Up to 4 lines
myHist = ROOT.TH1D('h_'+baseName, '', 100, 0, 50)
histTitles = ';#it{p}_{T} (GeV/c);# of dimuon'
myCut = f'{cutBase}'

midPt = Plot1D(tree=myTree, var='pt', cut=myCut, hist=myHist, histTitles=histTitles, leg=leg, legEntry=legEntry, canvas=myCanvas)
midPt.draw()
midPt.save(fPath='figs/', isPng=True, isPdf=False)

# ===== Mid phi ===== 
baseName = f'Data2023_{yRegion}_Phi'
myCanvas = TCanvas(baseName+yName, '', 800, 600)
leg = TLegend(0.8, 0.8, 0.9, 0.9)
legEntry = [legEntry1, legEntry2, '', ''] # Up to 4 lines
myHist = ROOT.TH1D('h_'+baseName, '', 30, -3.141592, 3.141592)
histTitles = ';#varphi_{1} (radian);# of dimuon'
myCut = f'{cutBase}'

midPhi = Plot1D(tree=myTree, var='phi1', cut=myCut, hist=myHist, histTitles=histTitles, leg=leg, legEntry=legEntry, canvas=myCanvas)
midPhi.draw()
midPhi.save(fPath='figs/', isPng=True, isPdf=False)

# ===== Mid phi HX ===== 
baseName = f'Data2023_{yRegion}_PhiHX'
myCanvas = TCanvas(baseName+yName, '', 800, 600)
leg = TLegend(0.8, 0.8, 0.9, 0.9)
legEntry = [legEntry1, legEntry2, '', ''] # Up to 4 lines
myHist = ROOT.TH1D('h_'+baseName, '', 30, -3.141592, 3.141592)
histTitles = ';#varphi_{HX} (radian);# of dimuon'
myCut = f'{cutBase}'

midPhiHX = Plot1D(tree=myTree, var='phiHX', cut=myCut, hist=myHist, histTitles=histTitles, leg=leg, legEntry=legEntry, canvas=myCanvas)
midPhiHX.draw()
midPhiHX.save(fPath='figs/', isPng=True, isPdf=False)

# ===== Mid phi CS ===== 
baseName = f'Data2023_{yRegion}_PhiCS'
myCanvas = TCanvas(baseName+yName, '', 800, 600)
leg = TLegend(0.8, 0.8, 0.9, 0.9)
legEntry = [legEntry1, legEntry2, '', ''] # Up to 4 lines
myHist = ROOT.TH1D('h_'+baseName, '', 30, -3.141592, 3.141592)
histTitles = ';#varphi_{CS} (radian);# of dimuon'
myCut = f'{cutBase}'

midPhiCS = Plot1D(tree=myTree, var='phiCS', cut=myCut, hist=myHist, histTitles=histTitles, leg=leg, legEntry=legEntry, canvas=myCanvas)
midPhiCS.draw()
midPhiCS.save(fPath='figs/', isPng=True, isPdf=False)

# ===== Mid cos HX ===== 
baseName = f'Data2023_{yRegion}_CosHX'
myCanvas = TCanvas(baseName+yName, '', 800, 600)
leg = TLegend(0.8, 0.8, 0.9, 0.9)
legEntry = [legEntry1, legEntry2, '', ''] # Up to 4 lines
myHist = ROOT.TH1D('h_'+baseName, '', 30, -1, 1)
histTitles = ';cos#theta_{HX};# of dimuon'
myCut = f'{cutBase}'

midCosHX = Plot1D(tree=myTree, var='cosHX', cut=myCut, hist=myHist, histTitles=histTitles, leg=leg, legEntry=legEntry, canvas=myCanvas)
midCosHX.draw()
midCosHX.save(fPath='figs/', isPng=True, isPdf=False)

# ===== Mid cos CS ===== 
baseName = f'Data2023_{yRegion}_CosCS'
myCanvas = TCanvas(baseName+yName, '', 800, 600)
leg = TLegend(0.8, 0.8, 0.9, 0.9)
legEntry = [legEntry1, legEntry2, '', ''] # Up to 4 lines
myHist = ROOT.TH1D('h_'+baseName, '', 30, -1, 1)
histTitles = ';cos#theta_{CS};# of dimuon'
myCut = f'{cutBase}'

midCosCS = Plot1D(tree=myTree, var='cosCS', cut=myCut, hist=myHist, histTitles=histTitles, leg=leg, legEntry=legEntry, canvas=myCanvas)
midCosCS.draw()
midCosCS.save(fPath='figs/', isPng=True, isPdf=False)


# ====================== #
# ===== Data - AllY ===== #
# ====================== #
# rapidity common
yRegion = 'AllY'
yName = '_pT6p5_50_y2p4' # Functionize - this part
legEntry1 = '|y| < 2.4'
legEntry2 = '6.5 < #it{p}_{T} < 50'
cutBase = 'pt > 6.5 && pt < 50 && fabs(y) > 0 && fabs(y) < 2.4'


# ===== AllY mass ===== 
baseName = f'Data2023_{yRegion}_Mass'
myCanvas = TCanvas(baseName+yName, '', 800, 600)
leg = TLegend(0.8, 0.8, 0.9, 0.9)
legEntry = [legEntry1, legEntry2, '', ''] # Up to 4 lines
myHist = ROOT.TH1D('h_'+baseName, '', 30, 2.6, 3.5)
histTitles = ';mass (GeV/c^{2});# of dimuon'
myCut = f'{cutBase}' # Functionize - this part

allYMass = Plot1D(tree=myTree, var='mass', cut=myCut, hist=myHist, histTitles=histTitles, leg=leg, legEntry=legEntry, canvas=myCanvas)
allYMass.draw()
allYMass.save(fPath='figs/', isPng=True, isPdf=False)

# ===== AllY pT ===== 
baseName = f'Data2023_{yRegion}_Pt'
myCanvas = TCanvas(baseName+yName, '', 800, 600)
leg = TLegend(0.8, 0.8, 0.9, 0.9)
legEntry = [legEntry1, legEntry2, '', ''] # Up to 4 lines
myHist = ROOT.TH1D('h_'+baseName, '', 100, 0, 50)
histTitles = ';#it{p}_{T} (GeV/c);# of dimuon'
myCut = f'{cutBase}'

allYPt = Plot1D(tree=myTree, var='pt', cut=myCut, hist=myHist, histTitles=histTitles, leg=leg, legEntry=legEntry, canvas=myCanvas)
allYPt.draw()
allYPt.save(fPath='figs/', isPng=True, isPdf=False)

# ===== AllY phi ===== 
baseName = f'Data2023_{yRegion}_Phi'
myCanvas = TCanvas(baseName+yName, '', 800, 600)
leg = TLegend(0.8, 0.8, 0.9, 0.9)
legEntry = [legEntry1, legEntry2, '', ''] # Up to 4 lines
myHist = ROOT.TH1D('h_'+baseName, '', 30, -3.141592, 3.141592)
histTitles = ';#varphi_{1} (radian);# of dimuon'
myCut = f'{cutBase}'

allYPhi = Plot1D(tree=myTree, var='phi1', cut=myCut, hist=myHist, histTitles=histTitles, leg=leg, legEntry=legEntry, canvas=myCanvas)
allYPhi.draw()
allYPhi.save(fPath='figs/', isPng=True, isPdf=False)

# ===== AllY phi HX ===== 
baseName = f'Data2023_{yRegion}_PhiHX'
myCanvas = TCanvas(baseName+yName, '', 800, 600)
leg = TLegend(0.8, 0.8, 0.9, 0.9)
legEntry = [legEntry1, legEntry2, '', ''] # Up to 4 lines
myHist = ROOT.TH1D('h_'+baseName, '', 30, -3.141592, 3.141592)
histTitles = ';#varphi_{HX} (radian);# of dimuon'
myCut = f'{cutBase}'

allYPhiHX = Plot1D(tree=myTree, var='phiHX', cut=myCut, hist=myHist, histTitles=histTitles, leg=leg, legEntry=legEntry, canvas=myCanvas)
allYPhiHX.draw()
allYPhiHX.save(fPath='figs/', isPng=True, isPdf=False)

# ===== AllY phi CS ===== 
baseName = f'Data2023_{yRegion}_PhiCS'
myCanvas = TCanvas(baseName+yName, '', 800, 600)
leg = TLegend(0.8, 0.8, 0.9, 0.9)
legEntry = [legEntry1, legEntry2, '', ''] # Up to 4 lines
myHist = ROOT.TH1D('h_'+baseName, '', 30, -3.141592, 3.141592)
histTitles = ';#varphi_{CS} (radian);# of dimuon'
myCut = f'{cutBase}'

allYPhiCS = Plot1D(tree=myTree, var='phiCS', cut=myCut, hist=myHist, histTitles=histTitles, leg=leg, legEntry=legEntry, canvas=myCanvas)
allYPhiCS.draw()
allYPhiCS.save(fPath='figs/', isPng=True, isPdf=False)

# ===== AllY cos HX ===== 
baseName = f'Data2023_{yRegion}_CosHX'
myCanvas = TCanvas(baseName+yName, '', 800, 600)
leg = TLegend(0.8, 0.8, 0.9, 0.9)
legEntry = [legEntry1, legEntry2, '', ''] # Up to 4 lines
myHist = ROOT.TH1D('h_'+baseName, '', 30, -1, 1)
histTitles = ';cos#theta_{HX};# of dimuon'
myCut = f'{cutBase}'

allYCosHX = Plot1D(tree=myTree, var='cosHX', cut=myCut, hist=myHist, histTitles=histTitles, leg=leg, legEntry=legEntry, canvas=myCanvas)
allYCosHX.draw()
allYCosHX.save(fPath='figs/', isPng=True, isPdf=False)

# ===== AllY cos CS ===== 
baseName = f'Data2023_{yRegion}_CosCS'
myCanvas = TCanvas(baseName+yName, '', 800, 600)
leg = TLegend(0.8, 0.8, 0.9, 0.9)
legEntry = [legEntry1, legEntry2, '', ''] # Up to 4 lines
myHist = ROOT.TH1D('h_'+baseName, '', 30, -1, 1)
histTitles = ';cos#theta_{CS};# of dimuon'
myCut = f'{cutBase}'

allYCosCS = Plot1D(tree=myTree, var='cosCS', cut=myCut, hist=myHist, histTitles=histTitles, leg=leg, legEntry=legEntry, canvas=myCanvas)
allYCosCS.draw()
allYCosCS.save(fPath='figs/', isPng=True, isPdf=False)


# cleaning
f.Close()


print('\n\n===== Finish drawMe.py =====\n')
exit(1)