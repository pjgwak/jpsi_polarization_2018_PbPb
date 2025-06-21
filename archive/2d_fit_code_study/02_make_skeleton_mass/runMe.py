import ROOT
from ROOT import gSystem
gSystem.Load('JpsiFitter.so')

from ROOT import RooWorkspace
from ROOT import ModelBuilder, JpsiFitter

# 이 부분은 좋아.
ws = RooWorkspace('ws', '')
model = ModelBuilder(ws)

myFitter = JpsiFitter(ws)
myFitter.processTree('/work/pjgwak/pol24/files_roodata/RooDataSet_miniAOD_isMC1_NP_Jpsi_cent0_180_Effw1_Accw1_PtW1_TnP1_HFNom_250221.root', 'myTree', 'mass', 'mass>2.6 && mass <3.5') # wrapper로 감싸기?
# /work/pjgwak/pol24/files_roodata/RooDataSet_miniAOD_isMC1_NP_Jpsi_cent0_180_Effw1_Accw1_PtW1_TnP1_HFNom_250221.root

model.buildMassSignal()
model.buildMassBkg()
model.buildMassModel()

# mass = ws.var("mass")
dataset = ws.data("ds_mc_temp")
pdf = ws.pdf("massModel")


mass = ws.var('mass')
mass.setRange(2.6, 3.5)
myFitter.print()
# exit(1)

mass.setRange("fit_region", 2.6, 3.5)
fit_result = pdf.fitTo(dataset, ROOT.RooFit.Save(), ROOT.RooFit.PrintLevel(-1), ROOT.RooFit.Verbose(False), ROOT.RooFit.Extended(1), ROOT.RooFit.Range("fit_region"))

c = ROOT.TCanvas("c", "Fit Plot", 800, 600)
c.SetLogy()

frame = mass.frame(ROOT.RooFit.Title("Fit of massModel to ds_mc_temp"), ROOT.RooFit.Range("fit_region"))
dataset.plotOn(frame, ROOT.RooFit.Range("fit_region"))


pdf.plotOn(frame, ROOT.RooFit.Range("fit_region"))
# ROOT.RooFit.NormRange("signal")


frame.Draw()
c.Update()

c.SaveAs("fit_plot.png")

print("\n--- Fit Results ---")
fit_result.Print()


# 일단 
# 1. input 파일 외부에서 넣어주기
# 2. cut 적용하기
# 3. pdf 종류 넣어주기 (sig, bkg)
# 4. 매개변수 외부에서 넣어주기


# = builder = #
# 변수 이름 세팅 - 유니크하게
# 나중에 변수, 함수 이름 알 수 있어야 한다.
# 내부적으로 sigModel, bkgModel 이름을 클래스가 가지고 있는 게 맞다.
# 코드 짜가려면 규칙성 있게 만들어서 나중에 이용자가 코드 짤 때 이해하기 쉽도록 할 것.
# mc mass는 build 할 때 signal만.

# = fitter = #
# 샘플 경로 전달 받아서 열기
# 샘플에 컷 적용
# 매개변수 범위 바꾸기 
# Fit option
# Drawing
# Saving


print('===== Finish runMe.py =====')