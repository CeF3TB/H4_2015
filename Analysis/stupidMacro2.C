{

  gStyle->SetOptStat(0);
  recoTree->Draw("cluster_pos_corr_hodoY1:cluster_pos_corr_hodoX1>>hnorm(18,-4.5,4.5,18,-4.5,4.5)"," wc_x>-20 && mcp_max_amplitude>200 ","");

  recoTree->Draw("cluster_pos_corr_hodoY1:cluster_pos_corr_hodoX1>>hmore(18,-4.5,4.5,18,-4.5,4.5)"," nino_LEtime-mcp_time_frac50>-10 && wc_x>-20 && (nino_LEtime-mcp_time_frac50)<-4 && mcp_max_amplitude>200 && (nino_LEtime-mcp_time_frac50)>-6 && (nino_LEtime-mcp_time_frac50)<-4 && cef3_maxAmpl[1]>700","");

  recoTree->Draw("cluster_pos_corr_hodoY1:cluster_pos_corr_hodoX1>>hless(18,-4.5,4.5,18,-4.5,4.5)"," nino_LEtime-mcp_time_frac50>-10 && wc_x>-20 && (nino_LEtime-mcp_time_frac50)<-4 && mcp_max_amplitude>200 && (nino_LEtime-mcp_time_frac50)<-6 && (nino_LEtime-mcp_time_frac50)>-10 && cef3_maxAmpl[1]>700","");


    


    hmore->Draw("colz");
    c1->SaveAs("plots_timing_posOpt/hmoreX_start.pdf");

    hless->Draw("colz");
    c1->SaveAs("plots_timing_posOpt/hlessX_start.pdf");


    TH2F* dummy = hnorm;
    hmore->Divide(dummy);
    hless->Divide(dummy);

    hmore->Draw("colz");
    c1->SaveAs("plots_timing_posOpt/hmore.pdf");

    hless->Draw("colz");
    c1->SaveAs("plots_timing_posOpt/hless.pdf");

    hnorm->Draw("colz");
    c1->SaveAs("plots_timing_posOpt/hnorm.pdf");




}
