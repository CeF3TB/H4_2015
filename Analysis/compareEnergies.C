{

  TStyle* style = new TStyle("DrawBaseStyle", "");
  style->SetCanvasColor(0);
  style->SetPadColor(0);
  style->SetFrameFillColor(0);
  style->SetStatColor(0);
  style->SetOptStat(0);
  style->SetTitleFillColor(0);
  style->SetCanvasBorderMode(0);
  style->SetPadBorderMode(0);
  style->SetFrameBorderMode(0);
  style->SetPadBottomMargin(0.12);
  style->SetPadLeftMargin(0.12);
  style->SetPalette(1,0);

///////// pretty palette ///////////

  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  style->SetNumberContours(NCont);



  style->cd();

  //for the histos
  style->SetHistLineWidth(2);

  // For the canvas:
  style->SetCanvasBorderMode(0);
  style->SetCanvasColor(kWhite);
  style->SetCanvasDefH(600); //Height of canvas
  style->SetCanvasDefW(600); //Width of canvas
  style->SetCanvasDefX(0); //POsition on screen
  style->SetCanvasDefY(0);

  // For the Pad:
  style->SetPadBorderMode(0);
  style->SetPadColor(kWhite);
  style->SetPadGridX(false);
  style->SetPadGridY(false);
  style->SetGridColor(0);
  style->SetGridStyle(3);
  style->SetGridWidth(1);

  // For the frame:
  style->SetFrameBorderMode(0);
  style->SetFrameBorderSize(1);
  style->SetFrameFillColor(0);
  style->SetFrameFillStyle(0);
  style->SetFrameLineColor(1);
  style->SetFrameLineStyle(1);
  style->SetFrameLineWidth(1);


  // Margins:
  style->SetPadTopMargin(0.05);
  style->SetPadBottomMargin(0.15);//0.13);
  style->SetPadLeftMargin(0.15);//0.16);
  style->SetPadRightMargin(0.10);//0.02);

  // For the Global title:

  style->SetOptTitle(0);
  style->SetTitleFont(42);
  style->SetTitleColor(1);
  style->SetTitleTextColor(1);
  style->SetTitleFillColor(10);
  style->SetTitleFontSize(0.05);

  // For the axis titles:

  style->SetTitleColor(1, "XYZ");
  style->SetTitleFont(42, "XYZ");
  style->SetTitleSize(0.05, "XYZ");
  style->SetTitleXOffset(1.15);//0.9);
  style->SetTitleYOffset(1.5); // => 1.15 if exponents

  // For the axis labels:

  style->SetLabelColor(1, "XYZ");
  style->SetLabelFont(42, "XYZ");
  style->SetLabelOffset(0.007, "XYZ");
  style->SetLabelSize(0.045, "XYZ");

  // For the axis:

  style->SetAxisColor(1, "XYZ");
  style->SetStripDecimals(kTRUE);
  style->SetTickLength(0.03, "XYZ");
  style->SetNdivisions(510, "XYZ");
  style->SetPadTickX(1); // To get tick marks on the opposite side of the frame
  style->SetPadTickY(1);

  style->cd();

  TLegend* lego = new TLegend(0.6, 0.7, 0.8, 0.9);
  lego->SetTextSize(0.036);
  lego->SetFillColor(0);
  lego->SetLineColor(0);




  TFile *_file0 = TFile::Open("plots_timingPerformance_SiPM_200GeV/plotOverallPerformances_V10.root");
  int  nbins=resVsAmplitude_channel->GetN();
  float maxX=(resVsAmplitude_channel->GetX())[nbins-1];

  float  maxY=(resVsAmplitude_channel->GetY())[1]*1.5;

  TH2F* axis = new TH2F("axis","axis",100,5,maxX,100,0,maxY);


  axis->GetYaxis()->SetTitle("#sigma_{t} [ps]");
  axis->GetXaxis()->SetTitle("Amplitude [ADC]");

  axis->Draw("P");
  resVsAmplitude_channel->SetLineColor(kBlack);
  resVsAmplitude_channel->SetMarkerColor(kBlack);
  resVsAmplitude_channel->Draw("EPsame");
  lego->AddEntry(  resVsAmplitude_channel ,"200 GeV", "P");

  TFile *_file0 = TFile::Open("plots_timingPerformance_SiPM_150GeV/plotOverallPerformances_V10.root");
  resVsAmplitude_channel->SetLineColor(kViolet);
  resVsAmplitude_channel->SetMarkerColor(kViolet);
  resVsAmplitude_channel->Draw("EPsame");
  lego->AddEntry(  resVsAmplitude_channel ,"150 GeV", "P");

  TFile *_file0 = TFile::Open("plots_timingPerformance_SiPM_100GeV/plotOverallPerformances_V10.root");
  resVsAmplitude_channel->SetLineColor(kRed);
  resVsAmplitude_channel->SetMarkerColor(kRed);
  resVsAmplitude_channel->Draw("EPsame");
  lego->AddEntry(  resVsAmplitude_channel ,"100 GeV", "P");

  TFile *_file0 = TFile::Open("plots_timingPerformance_SiPM_50GeV/plotOverallPerformances_V10.root");
  resVsAmplitude_channel->SetLineColor(kBlue);
  resVsAmplitude_channel->SetMarkerColor(kBlue);
  resVsAmplitude_channel->Draw("EPsame");
  lego->AddEntry(  resVsAmplitude_channel ,"50 GeV", "P");

  TFile *_file0 = TFile::Open("plots_timingPerformance_SiPM_20GeV/plotOverallPerformances_V10.root");
  resVsAmplitude_channel->SetLineColor(kGreen+2);
  resVsAmplitude_channel->SetMarkerColor(kGreen+2);
  resVsAmplitude_channel->Draw("EPsame");
  lego->AddEntry(  resVsAmplitude_channel ,"20 GeV", "P");


  lego->Draw("same");

   std::string constDirName = "plots_timingPerformance_SiPM";
   TString dir(constDirName);

   c1->SaveAs(dir+"/resVsAmplitude_EnergyComparison.png");
   c1->SaveAs(dir+"/resVsAmplitude_EnergyComparison.pdf");
}
