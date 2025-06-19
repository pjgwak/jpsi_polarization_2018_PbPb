## ==================================
from ROOT import TStyle, gROOT, TColor

def setStyle():
  atlasStyle  =  TStyle('ATLAS','Atlas style')

  # use plain black on white colors
  icol = 0 # WHITE
  atlasStyle.SetFrameBorderMode(icol)
  atlasStyle.SetFrameFillColor(icol)
  atlasStyle.SetCanvasBorderMode(icol)
  atlasStyle.SetCanvasColor(icol)
  atlasStyle.SetPadBorderMode(icol)
  atlasStyle.SetPadColor(icol)
  atlasStyle.SetStatColor(icol)

  # set the paper & margin sizes
  atlasStyle.SetPaperSize(20,26)

  # set margin sizes
  atlasStyle.SetPadTopMargin(0.05)

  atlasStyle.SetPadRightMargin(0.05)
  atlasStyle.SetPadBottomMargin(0.16)
  atlasStyle.SetPadLeftMargin(0.16)

  # set title offsets (for axis label)
  atlasStyle.SetTitleXOffset(1.4)
  atlasStyle.SetTitleYOffset(1.5)

  # use large fonts
  font = 43 # Helvetica  - fixed size
  tsize = 27
  atlasStyle.SetTextFont(font)

  atlasStyle.SetTextSize(tsize)
  atlasStyle.SetLabelFont(font,'x')
  atlasStyle.SetTitleFont(font,'x')
  atlasStyle.SetLabelFont(font,'y')
  atlasStyle.SetTitleFont(font,'y')
  atlasStyle.SetLabelFont(font,'z')
  atlasStyle.SetTitleFont(font,'z')

  atlasStyle.SetLabelSize(tsize,'x')
  atlasStyle.SetTitleSize(tsize,'x')
  atlasStyle.SetLabelSize(tsize,'y')
  atlasStyle.SetTitleSize(tsize,'y')
  atlasStyle.SetLabelSize(tsize,'z')
  atlasStyle.SetTitleSize(tsize,'z')

  # use bold lines and markers
  atlasStyle.SetMarkerStyle(20)
  atlasStyle.SetMarkerSize(1.2)
  atlasStyle.SetHistLineWidth(2)
  atlasStyle.SetLineStyleString(2,'[12 12]') # postscript dashes

  # get rid of X error bars
  atlasStyle.SetErrorX(0.5)
  # get rid of error bar caps
  atlasStyle.SetEndErrorSize(0.)

  # default stats
  atlasStyle.SetOptStat(0)
  atlasStyle.SetOptFit(0)
  atlasStyle.SetOptFit(0)

  # title
  atlasStyle.SetOptTitle(0)
  atlasStyle.SetTitleBorderSize(0)
  atlasStyle.SetTitleFillColor(icol)
  atlasStyle.SetTitleFont(font,'')
  atlasStyle.SetTitleStyle(0)

  # put tick marks on top and RHS of plots
  atlasStyle.SetPadTickX(1)
  atlasStyle.SetPadTickY(1)
  
  # plot title offset
  atlasStyle.SetTitleOffset(1.7,'y')
  atlasStyle.SetTitleOffset(1.7,'z')

  # legend
  atlasStyle.SetLegendFont(font)
  atlasStyle.SetLegendBorderSize(0)

  # rainbow palette
  atlasStyle.SetPalette(1)
  NCont = 255

  ## rainbow palette
  gROOT.ProcessLine("const Int_t NRGBs = 5;")
  gROOT.ProcessLine("Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };")
  gROOT.ProcessLine("Double_t red  [NRGBs] = { 0.00, 0.00, 0.87, 1.00, 0.51 };")
  gROOT.ProcessLine("Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };")
  gROOT.ProcessLine("Double_t blue [NRGBs] = { 0.51, 1.00, 0.12, 0.00, 0.00 };")
  gROOT.ProcessLine("TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, %i);" % NCont)
  
  gROOT.ProcessLine("TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, %i);" % NCont)

  atlasStyle.SetNumberContours(NCont) 

  ## set ATLAS style
  gROOT.SetStyle('ATLAS')
  gROOT.ForceStyle()

  from ROOT import gStyle
  gStyle.Print()
  print(gStyle.GetErrorX())