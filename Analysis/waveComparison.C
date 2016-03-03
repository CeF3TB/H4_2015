{
  TFile *_file0 = TFile::Open("analysisTrees_V10/Reco_3971_waveComp.root");
  TFile *_file1 = TFile::Open("analysisTrees_V10/Reco_3990_waveComp.root");
  TFile *_file2 = TFile::Open("analysisTrees_V10/Reco_4001_waveComp.root");
  TFile *_file3 = TFile::Open("analysisTrees_V10/Reco_4026_waveComp.root");
  TGraph* wave_50= (TGraph*)_file0->Get("waveform_2");
  TGraph* wave_100= (TGraph*)_file1->Get("waveform_2");
  TGraph* wave_150= (TGraph*)_file2->Get("waveform_2");
  TGraph* wave_200= (TGraph*)_file3->Get("waveform_2");
 
  //  TH2F* histo_50=new TH2F("histo_50","histo_50",400,0,400,1100,0,1.1);
  TH1F* histo_50=new TH1F("histo_50","histo_50",400,0,400);
  TH1F* histo_100=new TH1F("histo_100","histo_100",400,0,400);
  TH1F* histo_150=new TH1F("histo_150","histo_150",400,0,400);
  TH1F* histo_200=new TH1F("histo_200","histo_200",400,0,400);

  float max_50 = TMath::MaxElement(wave_50->GetN(),wave_50->GetY());
  float max_100 = TMath::MaxElement(wave_100->GetN(),wave_100->GetY());
  float max_150 = TMath::MaxElement(wave_150->GetN(),wave_150->GetY());
  float max_200 = TMath::MaxElement(wave_200->GetN(),wave_200->GetY());

  for(int i=0; i<wave_50->GetN();++i){
    histo_50->SetBinContent(wave_50->GetX()[i]*1.e9,wave_50->GetY()[i]/max_50);
    histo_100->SetBinContent(wave_100->GetX()[i]*1.e9,wave_100->GetY()[i]/max_100);
    histo_150->SetBinContent(wave_150->GetX()[i]*1.e9,wave_150->GetY()[i]/max_150);
    histo_200->SetBinContent(wave_200->GetX()[i]*1.e9,wave_200->GetY()[i]/max_200);

    histo_50->GetXaxis()->SetTitle("time [ns]");
    histo_100->GetXaxis()->SetTitle("time [ns]");
    histo_150->GetXaxis()->SetTitle("time [ns]");
    histo_200->GetXaxis()->SetTitle("time [ns]");
  }

  histo_100->SetLineColor(kBlue);
  histo_150->SetLineColor(kRed);
  histo_200->SetLineColor(kViolet);

  TCanvas* c1=new TCanvas();
  histo_50->GetXaxis()->SetRangeUser(0,200);
  histo_50->Draw("l");
  histo_100->Draw("lsame");
  histo_150->Draw("lsame");
  histo_200->Draw("lsame");


  c1->SaveAs("wave.pdf");

  histo_50->GetXaxis()->SetRangeUser(20,60);
  histo_50->Draw("l");
  histo_100->Draw("lsame");
  histo_150->Draw("lsame");
  histo_200->Draw("lsame");
  c1->SaveAs("wave_zoom.pdf");

  c1->Clear();
  c1->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0.); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();  

  histo_50->GetXaxis()->SetRangeUser(20,200);
  histo_50->GetYaxis()->SetRangeUser(0.01,1.1);
  histo_50->Draw("l");
  histo_100->Draw("lsame");



  //histo_50->GetYaxis()->SetLabelSize(0.);
  //  TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
  //  TGaxis *axis = new TGaxis( 0, 400, 0, 1.1, 0,1.0,510,"");
  //axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  //axis->SetLabelSize(15);
  //axis->Draw();

  // lower plot will be in pad
  c1->cd();          // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0);
  //  pad2->SetBottomMargin(0.);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();

  TH1F *h3 = (TH1F*)histo_50->Clone("h3");
  h3->SetLineColor(kBlack);
  h3->SetMinimum(0.8);  // Define Y ..
  h3->SetMaximum(1.13); // .. range
  h3->Divide(histo_100);

  h3->GetXaxis()->SetTitle("time [ns]");
  h3->SetMarkerStyle(21);
  h3->Draw("p");    

  gPad->RedrawAxis();

  c1->SaveAs("wave_ratio50_100.pdf");

  //50-200 ratio
  c1->Clear();
  c1->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0.); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();  

  histo_50->GetXaxis()->SetRangeUser(20,200);
  histo_50->GetYaxis()->SetRangeUser(0.01,1.1);
  histo_50->Draw("l");
  histo_200->Draw("lsame");



  //histo_50->GetYaxis()->SetLabelSize(0.);
  //  TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
  //  TGaxis *axis = new TGaxis( 0, 400, 0, 1.1, 0,1.0,510,"");
  //axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  //axis->SetLabelSize(15);
  //axis->Draw();

  // lower plot will be in pad
  c1->cd();          // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0);
  //  pad2->SetBottomMargin(0.);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();

  TH1F *h3 = (TH1F*)histo_50->Clone("h3");
  h3->SetLineColor(kBlack);
  h3->SetMinimum(0.8);  // Define Y ..
  h3->SetMaximum(1.13); // .. range
  h3->Divide(histo_200);

  h3->GetXaxis()->SetTitle("time [ns]");
  h3->SetMarkerStyle(21);
  h3->Draw("p");    

  gPad->RedrawAxis();

  c1->SaveAs("wave_ratio50_200.pdf");


  //100-200 ratio
  c1->Clear();
  c1->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0.); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();  

  histo_100->GetXaxis()->SetRangeUser(20,200);
  histo_100->GetYaxis()->SetRangeUser(0.01,1.1);
  histo_100->Draw("l");
  histo_200->Draw("lsame");



  //histo_100->GetYaxis()->SetLabelSize(0.);
  //  TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
  //  TGaxis *axis = new TGaxis( 0, 400, 0, 1.1, 0,1.0,510,"");
  //axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  //axis->SetLabelSize(15);
  //axis->Draw();

  // lower plot will be in pad
  c1->cd();          // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0);
  //  pad2->SetBottomMargin(0.);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();

  TH1F *h3 = (TH1F*)histo_100->Clone("h3");
  h3->SetLineColor(kBlack);
  h3->SetMinimum(0.8);  // Define Y ..
  h3->SetMaximum(1.13); // .. range
  h3->Divide(histo_200);

  h3->GetXaxis()->SetTitle("time [ns]");
  h3->SetMarkerStyle(21);
  h3->Draw("p");    

  gPad->RedrawAxis();

  c1->SaveAs("wave_ratio100_200.pdf");

  //100-150 ratio
  c1->Clear();
  c1->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0.); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();  

  histo_100->GetXaxis()->SetRangeUser(20,200);
  histo_100->GetYaxis()->SetRangeUser(0.01,1.1);
  histo_100->Draw("l");
  histo_150->Draw("lsame");



  //histo_100->GetYaxis()->SetLabelSize(0.);
  //  TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
  //  TGaxis *axis = new TGaxis( 0, 400, 0, 1.1, 0,1.0,510,"");
  //axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  //axis->SetLabelSize(15);
  //axis->Draw();

  // lower plot will be in pad
  c1->cd();          // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0);
  //  pad2->SetBottomMargin(0.);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();

  TH1F *h3 = (TH1F*)histo_100->Clone("h3");
  h3->SetLineColor(kBlack);
  h3->SetMinimum(0.8);  // Define Y ..
  h3->SetMaximum(1.13); // .. range
  h3->Divide(histo_150);

  h3->GetXaxis()->SetTitle("time [ns]");
  h3->SetMarkerStyle(21);
  h3->Draw("p");    

  gPad->RedrawAxis();

  c1->SaveAs("wave_ratio100_150.pdf");


}
