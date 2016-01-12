{

  recoTree->Draw("wc_y_corr>>h2(80,-20.5,19.5)","cef3_time_at_thresh[1]>0 && cef3_time_at_thresh[1]<50 && mcp_max_amplitude>200 && cef3_maxAmpl[2]>100 && TMath::Abs(cef3_time_at_frac50[1]-mcp_time_frac50)<11" ,"");
  c1->SaveAs("plots_timing_posOpt/y_total.pdf");

  recoTree->Draw("wc_x_corr>>h2(80,-20.5,19.5)","cef3_time_at_thresh[1]>0 && cef3_time_at_thresh[1]<50 && mcp_max_amplitude>200 && cef3_maxAmpl[2]>100 && TMath::Abs(cef3_time_at_frac50[1]-mcp_time_frac50)<11" ,"");
  c1->SaveAs("plots_timing_posOpt/x_total.pdf");


  recoTree->Draw("wc_y_corr>>h2(80,-20.5,19.5)","cef3_time_at_thresh[1]>0 && cef3_time_at_thresh[1]<50 && mcp_max_amplitude>200 && cef3_maxAmpl[2]>100 && TMath::Abs(cef3_time_at_frac50[1]-mcp_time_frac50)<11 && wc_x>-20 && (cef3_time_at_frac50[1]-mcp_time_frac50)>7 && (cef3_time_at_frac50[1]-mcp_time_frac50)<8.6","");
  c1->SaveAs("plots_timing_posOpt/y_lower.pdf");
  recoTree->Draw("wc_y_corr>>h2(80,-20.5,19.5)","cef3_time_at_thresh[1]>0 && cef3_time_at_thresh[1]<50 && mcp_max_amplitude>200 && cef3_maxAmpl[2]>100 && TMath::Abs(cef3_time_at_frac50[1]-mcp_time_frac50)<11 && wc_x>-20 && (cef3_time_at_frac50[1]-mcp_time_frac50)>7 && (cef3_time_at_frac50[1]-mcp_time_frac50)>8.6","");
 c1->SaveAs("plots_timing_posOpt/y_upper.pdf");

  recoTree->Draw("wc_x_corr>>h2(80,-20.5,19.5)","cef3_time_at_thresh[1]>0 && cef3_time_at_thresh[1]<50 && mcp_max_amplitude>200 && cef3_maxAmpl[2]>100 && TMath::Abs(cef3_time_at_frac50[1]-mcp_time_frac50)<11 && wc_x>-20 && (cef3_time_at_frac50[1]-mcp_time_frac50)>7 && (cef3_time_at_frac50[1]-mcp_time_frac50)<8.6","");
  c1->SaveAs("plots_timing_posOpt/x_lower.pdf");
  recoTree->Draw("wc_x_corr>>h2(80,-20.5,19.5)","cef3_time_at_thresh[1]>0 && cef3_time_at_thresh[1]<50 && mcp_max_amplitude>200 && cef3_maxAmpl[2]>100 && TMath::Abs(cef3_time_at_frac50[1]-mcp_time_frac50)<11 && wc_x>-20 && (cef3_time_at_frac50[1]-mcp_time_frac50)>7 && (cef3_time_at_frac50[1]-mcp_time_frac50)>8.6","");
 c1->SaveAs("plots_timing_posOpt/x_upper.pdf");

  for (int i=0;i<10;++i){

    float lowerCut=-6.5+i*1.5;
    float upperCut=-6.5+(i+1)*1.5;

    TString lowerString,upperString,istring;

    istring+=i;
    lowerString+=lowerCut;
    upperString+=upperCut;

    

    recoTree->Draw("wc_y_corr>>h2(80,-20.5,19.5)","cef3_time_at_thresh[1]>0 && cef3_time_at_thresh[1]<50 && mcp_max_amplitude>200 && cef3_maxAmpl[2]>100 && TMath::Abs(cef3_time_at_frac50[1]-mcp_time_frac50)<11 && wc_x>-20 && (cef3_time_at_frac50[1]-mcp_time_frac50)>7 && (cef3_time_at_frac50[1]-mcp_time_frac50)<8.6 && wc_x_corr>"+lowerString+ " && wc_x_corr<"+upperString,"");
    c1->SaveAs("plots_timing_posOpt/wc_y_xbin"+istring+".pdf");

  }

}
