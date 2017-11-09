#include <fstream>
#include <iostream>
#include <cstdio>
#include <sstream>
#include <cmath>
#include <iomanip>

#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TExec.h>
#include <TStyle.h>
#include <TMath.h>
#include <TROOT.h>
#include <THStack.h>
#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TF1.h>
#include <TClass.h>
#include <TError.h>

#include "DatacardMaker.h"
#include "AnalysisConfig.h"
#include "higgsUtils.h"
#include "Samples.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/RootFileReader.h"
#include "../../common/include/plotterUtils.h"

#include "HistoListReader.h"

#include <TList.h>
#include <TKey.h>
#include <TIterator.h>
#include <TObject.h>
#include <TObjString.h>
#include  <string>
#include "TSystem.h"
#include "TPRegexp.h"
#include "TClonesArray.h"
#include "TObjString.h"
#include <dirent.h>
#include <map>
#include <TObjArray.h>
#include <cstdlib>
#include <algorithm>

#include <string>



// Aids the find algorithm to access the first element of the pair object    
struct comp
{
  comp(std::string const& s) : _s(s) { }
  
  bool operator () (std::pair<std::string, std::string> const& p) {  
    return (p.first == _s);
  }

  std::string _s;
};


struct compare
{
  compare(std::string const& s, std::string const& t) : _s(s),_t(t) { }

  bool operator () (std::pair<std::string, std::string> const& p) {
    return ((p.first == _s) && (p.second == _t));
  }

  std::string _s;
  std::string _t;
};


DatacardMaker::DatacardMaker(const AnalysisConfig& analysisConfig,
                             const std::vector<std::string>& v_plot,
                             const std::vector<Channel::Channel>& v_channel,
                             const std::vector<Systematic::Systematic>& v_systematic,
                             const std::string& fileLists,
                             const std::string& configname) :
  signalModel_("ttH"),
  fileList_base_(fileLists),
  configname_(configname),
  outputBaseDir_("datacards"),
  outputFileName_("common/ttH_hbb_13TeV_dl.root"),
  addSystematicUncertainty_(true),
  addStatisticalUncertainty_(true),
  analysisConfig_(analysisConfig),
  inputFileLists_(NULL),
  fileReader_(RootFileReader::getInstance()),
  eventFileString_(""),
  datacardName_(""),
  datacard_(NULL),
  outputDirDatacard_(""),
  outputFile_(NULL),
  observableType_("BDT"),
  pruneBinByBin_(false),
  v_plot_(v_plot),
  v_channel_(v_channel),
  v_systematic_(v_systematic)
{
  // Intialize data members based on era (Switch on MC bin-by-bin statistical prunning by setting argument to true)
  if(analysisConfig_.general().era_ != Era::run1_8tev)
    initialization13TeV(pruneBinByBin_);

  // Setting the list of names of the input root files and processes (e.g. ttH & ttbar+XX) to be used
  const std::string histoListFile(tth::DATA_PATH_TTH() + "/" + configname_);
  std::cout << histoListFile << std::endl;

  const HistoListReader histoList(histoListFile.data());

  if(histoList.isZombie()){
    std::cerr << "Error in Histo! Cannot find HistoList with name: " << histoListFile << "\n...break\n" << std::endl;
    exit(12);      
  }
  
  // Loop over all histograms in histoList and print them  
  for(auto it = histoList.begin(); it != histoList.end(); ++it){
    
    // Access plot properties from histoList and check whether histogram name contains name pattern      
    const PlotProperties& plotProperties = it->second;

    std::cout << "checking " << plotProperties.name << std::endl;
    
    bool found = false;
    
    for(const auto& plot : v_plot_){
      
      if(plot.size() && plot[0] == '+'){
        
        if(plotProperties.name.CompareTo(&plot[1], TString::kIgnoreCase) == 0){
          found = true;                                                                                                                                     
          break;
        }  
      }
      else if(plotProperties.name.Contains(plot, TString::kIgnoreCase)){
        found = true;
        break;
      }
    }
    if(!found){ 
      std::cout << "... no histograms found, continue with next\n";
      continue;
    }
    
    std::cout << "finished checking " << plotProperties.name << "\n\n";
    
    processNames_.clear();
    
    TObjArray* strings = plotProperties.specialComment.Tokenize(" ");
    
    for(std::size_t i = 0; i < (std::size_t)strings->GetEntries(); ++i) {

      std::string sample = ((TObjString*)strings->At(i))->GetString().Data();
      
      if(std::find_if(convertSampleNames_.begin(), convertSampleNames_.end(), comp(sample)) !=  convertSampleNames_.end())
        processNames_.push_back(sample);
      else 
        continue;
    }
    fileNames_.push_back(plotProperties.name.Data());
  }

}


void DatacardMaker::initialization13TeV(bool pruneOption)
{
  // Set prunning option
  pruneBinByBin_ = pruneOption;

  // Event categories label conversion (Categories are based on jet/b-tag jet multiplicity)
  mapOfCategories_["cate0"] = "dl_3j2t";
  mapOfCategories_["cate1"] = "dl_3j3t";
  mapOfCategories_["cate2"] = "dl_ge4j2t";
  mapOfCategories_["cate3"] = "dl_ge4j3t";
  mapOfCategories_["cate4"] = "dl_ge4jge4t";
  mapOfCategories_["cate3Low"] = "dl_ge4j3t_low";
  mapOfCategories_["cate3High"] = "dl_ge4j3t_high";
  mapOfCategories_["cate4Low"] = "dl_ge4jge4t_low";
  mapOfCategories_["cate4High"] = "dl_ge4jge4t_high";

  // Convert internal process name to CMS ttH collaboration labeling
  convertSampleNames_["data"]         = "data_obs";
  convertSampleNames_["allmc"]        = "data_obs";
  convertSampleNames_["ttH"]          = "ttH";
  convertSampleNames_["notttH"]       = "nottH";
  convertSampleNames_["ttHnobb"]      = "ttH_nobb";
  convertSampleNames_["ttHbb"]         = "ttH_hbb";
  convertSampleNames_["ttHcc"]         = "ttH_hcc";
  convertSampleNames_["ttHww"]         = "ttH_hww";
  convertSampleNames_["ttHzz"]         = "ttH_hzz";
  convertSampleNames_["ttHtautau"]     = "ttH_htt";
  convertSampleNames_["ttHgluongluon"] = "ttH_hgluglu";
  convertSampleNames_["ttHgammagamma"] = "ttH_hgg";
  convertSampleNames_["ttHzgamma"]     = "ttH_hzg";
  convertSampleNames_["ttb"]      = "ttbarPlusB";
  convertSampleNames_["ttbb"]     = "ttbarPlusBBbar";
  convertSampleNames_["ttcc"]     = "ttbarPlusCCbar";
  convertSampleNames_["tt2b"]     = "ttbarPlus2B";
  convertSampleNames_["ttOther"]  = "ttbarOther";
  convertSampleNames_["ttZ"]      = "ttbarZ";
  convertSampleNames_["ttW"]      = "ttbarW";
  convertSampleNames_["singleTop"] = "singlet";
  convertSampleNames_["diboson"]   = "diboson";
  convertSampleNames_["wlnu"]    = "wjets";
  convertSampleNames_["w"]       = "wjets";
  convertSampleNames_["dy"]      = "zjets";
  convertSampleNames_["dyee"]    = "zjets";
  convertSampleNames_["dymumu"]    = "zjets";
  convertSampleNames_["dytautau"]  = "zjets";
  convertSampleNames_["ttGamma"]   = "ttGamma";
  convertSampleNames_["ww"]     = "ww";
  convertSampleNames_["wz"]     = "wz";
  convertSampleNames_["zz"]     = "zz";
  convertSampleNames_["www"]    = "ww";
  convertSampleNames_["wwz"]    = "wwz";
  convertSampleNames_["zzz"]    = "zzz";
  convertSampleNames_["qcd"]    = "qcd";
  convertSampleNames_["minorBkg"]   = "minorBkg";

  // Convert internal systematic uncertainty labeling to CMS ttH collaboration naming conventiion
  convertSystematicLabel_["BTAGDISCR_BPURITY_UP"]   = "CMS_btag_hfUp";//CMS_ttHbb_CSVHF
  convertSystematicLabel_["BTAGDISCR_BPURITY_DOWN"] = "CMS_btag_hfDown";
  convertSystematicLabel_["BTAGDISCR_BSTAT1_UP"]    = "CMS_btag_hfstats1Up";//CMS_ttHbb_CSVHFStats1
  convertSystematicLabel_["BTAGDISCR_BSTAT1_DOWN"]  = "CMS_btag_hfstats1Down";
  convertSystematicLabel_["BTAGDISCR_BSTAT2_UP"]    = "CMS_btag_hfstats2Up";//CMS_ttHbb_CSVHFStats2
  convertSystematicLabel_["BTAGDISCR_BSTAT2_DOWN"]  = "CMS_btag_hfstats2Down";
  convertSystematicLabel_["BTAGDISCR_LSTAT1_UP"]    = "CMS_btag_lfstats1Up";//CMS_ttHbb_CSVLFStats1
  convertSystematicLabel_["BTAGDISCR_LSTAT1_DOWN"]  = "CMS_btag_lfstats1Down";
  convertSystematicLabel_["BTAGDISCR_LSTAT2_UP"]    = "CMS_btag_lfstats2Up";//CMS_ttHbb_CSVLFStats2
  convertSystematicLabel_["BTAGDISCR_LSTAT2_DOWN"]  = "CMS_btag_lfstats2Down";
  convertSystematicLabel_["BTAGDISCR_LPURITY_UP"]   = "CMS_btag_lfUp";//CMS_ttHbb_CSVLF
  convertSystematicLabel_["BTAGDISCR_LPURITY_DOWN"] = "CMS_btag_lfDown";
  convertSystematicLabel_["BTAGDISCR_CERR1_UP"]     = "CMS_btag_cferr1Up";//CMS_ttHbb_CSVCErr1
  convertSystematicLabel_["BTAGDISCR_CERR1_DOWN"]   = "CMS_btag_cferr1Down";
  convertSystematicLabel_["BTAGDISCR_CERR2_UP"]     = "CMS_btag_cferr2Up";//CMS_ttHbb_CSVCErr2
  convertSystematicLabel_["BTAGDISCR_CERR2_DOWN"]   = "CMS_btag_cferr2Down";
  convertSystematicLabel_["Nominal"]                = "nominal";
  convertSystematicLabel_["JES_UP"]                 = "CMS_scale_jUp";
  convertSystematicLabel_["JES_DOWN"]               = "CMS_scale_jDown";
  convertSystematicLabel_["JER_UP"]   = "CMS_res_jUp";
  convertSystematicLabel_["JER_DOWN"] = "CMS_res_jDown";
  convertSystematicLabel_["LUMI_UP"]   = "lumi_13TeV_2016Up";
  convertSystematicLabel_["LUMI_DOWN"] = "lumi_13TeV_2016Down";

  convertSystematicLabel_["XSEC_TTBB_UP"]   = "bgnorm_ttbarPlusBBbarUp";
  convertSystematicLabel_["XSEC_TTBB_DOWN"] = "bgnorm_ttbarPlusBBbarDown";
  convertSystematicLabel_["XSEC_TTB_UP"]   = "bgnorm_ttbarPlusBUp";
  convertSystematicLabel_["XSEC_TTB_DOWN"] = "bgnorm_ttbarPlusBDown";
  convertSystematicLabel_["XSEC_TT2B_UP"]   = "bgnorm_ttbarPlus2BUp";
  convertSystematicLabel_["XSEC_TT2B_DOWN"] = "bgnorm_ttbarPlus2BDown";
  convertSystematicLabel_["XSEC_TTCC_UP"]   = "bgnorm_ttbarPlusCCbarUp";
  convertSystematicLabel_["XSEC_TTCC_DOWN"] = "bgnorm_ttbarPlusCCbarDown";
  convertSystematicLabel_["XSEC_TTOTHER_UP"]   = "bgnorm_ttbarPlusOtherUp";
  convertSystematicLabel_["XSEC_TTOTHER_DOWN"] = "bgnorm_ttbarPlusOtherDown";
  convertSystematicLabel_["XSEC_TTZ_UP"]   = "CMS_ttHbb_QCDscale_ttbarZUp";
  convertSystematicLabel_["XSEC_TTZ_DOWN"] = "CMS_ttHbb_QCDscale_ttbarZDown";
  convertSystematicLabel_["XSEC_TTW_UP"]   = "CMS_ttHbb_QCDscale_ttbarWUp";
  convertSystematicLabel_["XSEC_TTW_DOWN"] = "CMS_ttHbb_QCDscale_ttbarWDown";
  convertSystematicLabel_["XSEC_TTH_UP"]   = "QCDscale_ttHUp";
  convertSystematicLabel_["XSEC_TTH_DOWN"] = "QCDscale_ttHDown";
  convertSystematicLabel_["XSEC_TT_UP"]   = "QCDscale_ttbarUp";
  convertSystematicLabel_["XSEC_TT_DOWN"] = "QCDscale_ttbarDown";

  convertSystematicLabel_["NORMPDFGG_UP"]   = "pdf_ggUp";
  convertSystematicLabel_["NORMPDFGG_DOWN"] = "pdf_ggDown";
  convertSystematicLabel_["NORMPDFGQ_UP"]   = "pdf_qgUp";
  convertSystematicLabel_["NORMPDFGQ_DOWN"] = "pdf_qgDown";
  convertSystematicLabel_["NORMPDFQQ_UP"]   = "pdf_qqbarUp";
  convertSystematicLabel_["NORMPDFQQ_DOWN"] = "pdf_qqbarDown";
  convertSystematicLabel_["NORMPDFTTH_UP"]   = "pdf_Higgs_ttHUp";//pdf_gg_ttH
  convertSystematicLabel_["NORMPDFTTH_DOWN"] = "pdf_Higgs_ttHDown";//pdf_gg_ttH

  convertSystematicLabel_["XSEC_V_UP"]   = "QCDscale_VUp";
  convertSystematicLabel_["XSEC_V_DOWN"] = "QCDscale_VDown";
  convertSystematicLabel_["XSEC_VV_UP"]   = "QCDscale_VVUp";
  convertSystematicLabel_["XSEC_VV_DOWN"] = "QCDscale_VVDown";

  convertSystematicLabel_["XSEC_T_UP"] = "QCDscale_singletUp";
  convertSystematicLabel_["XSEC_T_DOWN"] = "QCDscale_singletDown";

  convertSystematicLabel_["SCALE_TTB_UP"]   = "CMS_ttHbb_Scale_ttbarPlusBUp";
  convertSystematicLabel_["SCALE_TTB_DOWN"] = "CMS_ttHbb_Scale_ttbarPlusBDown";
  convertSystematicLabel_["SCALE_TT2B_UP"]   = "CMS_ttHbb_Scale_ttbarPlus2BUp";
  convertSystematicLabel_["SCALE_TT2B_DOWN"] = "CMS_ttHbb_Scale_ttbarPlus2BDown";
  convertSystematicLabel_["SCALE_TTBB_UP"]   = "CMS_ttHbb_Scale_ttbarPlusBBbarUp";
  convertSystematicLabel_["SCALE_TTBB_DOWN"] = "CMS_ttHbb_Scale_ttbarPlusBBbarDown";
  convertSystematicLabel_["SCALE_TTCC_UP"]   = "CMS_ttHbb_Scale_ttbarPlusCCbarUp";
  convertSystematicLabel_["SCALE_TTCC_DOWN"] = "CMS_ttHbb_Scale_ttbarPlusCCbarDown";
  convertSystematicLabel_["SCALE_TTOTHER_UP"]   = "CMS_ttHbb_Scale_ttbarOtherUp"; 
  convertSystematicLabel_["SCALE_TTOTHER_DOWN"] = "CMS_ttHbb_Scale_ttbarOtherDown";

  convertSystematicLabel_["MESCALE_TTB_UP"]   = "CMS_ttHbb_Q2scale_ttbarPlusBUp";
  convertSystematicLabel_["MESCALE_TTB_DOWN"] = "CMS_ttHbb_Q2scale_ttbarPlusBDown";
  convertSystematicLabel_["MESCALE_TT2B_UP"]   = "CMS_ttHbb_Q2scale_ttbarPlus2BUp";
  convertSystematicLabel_["MESCALE_TT2B_DOWN"] = "CMS_ttHbb_Q2scale_ttbarPlus2BDown";
  convertSystematicLabel_["MESCALE_TTBB_UP"]   = "CMS_ttHbb_Q2scale_ttbarPlusBBbarUp";
  convertSystematicLabel_["MESCALE_TTBB_DOWN"] = "CMS_ttHbb_Q2scale_ttbarPlusBBbarDown";
  convertSystematicLabel_["MESCALE_TTCC_UP"]   = "CMS_ttHbb_Q2scale_ttbarPlusCCbarUp";
  convertSystematicLabel_["MESCALE_TTCC_DOWN"] = "CMS_ttHbb_Q2scale_ttbarPlusCCbarDown";
  convertSystematicLabel_["MESCALE_TTOTHER_UP"]   = "CMS_ttHbb_Q2scale_ttbarOtherUp";
  convertSystematicLabel_["MESCALE_TTOTHER_DOWN"] = "CMS_ttHbb_Q2scale_ttbarOtherDown";

  convertSystematicLabel_["MEFACSCALE_DOWN"] = "CMS_ttHbb_scaleMuFDown";//meFacScale
  convertSystematicLabel_["MEFACSCALE_UP"]   = "CMS_ttHbb_scaleMuFUp";//meFacScale

  convertSystematicLabel_["MERENSCALE_DOWN"] = "CMS_ttHbb_scaleMuRDown";//meRenScale
  convertSystematicLabel_["MERENSCALE_UP"]   = "CMS_ttHbb_scaleMuRUp";

  convertSystematicLabel_["PSFSRSCALE_UP"]   = "CMS_ttHbb_FSRUp";//psFSRScale
  convertSystematicLabel_["PSFSRSCALE_DOWN"] = "CMS_ttHbb_FSRDown";//psFSRScale  

  convertSystematicLabel_["PSISRSCALE_UP"]   = "CMS_ttHbb_ISRUp";//psISRScale
  convertSystematicLabel_["PSISRSCALE_DOWN"] = "CMS_ttHbb_ISRDown";//psISRScale

  convertSystematicLabel_["UETUNE_UP"]   = "CMS_ttHbb_UEUp";//uetune
  convertSystematicLabel_["UETUNE_DOWN"] = "CMS_ttHbb_UEDown";//uetune

  convertSystematicLabel_["UETUNE_TTBB_UP"]      = "CMS_ttHbb_UE_ttbarPlusBBbarUp";
  convertSystematicLabel_["UETUNE_TTBB_DOWN"]    = "CMS_ttHbb_UE_ttbarPlusBBbarDown";
  convertSystematicLabel_["UETUNE_TT2B_UP"]      = "CMS_ttHbb_UE_ttbarPlus2BUp";
  convertSystematicLabel_["UETUNE_TT2B_DOWN"]    = "CMS_ttHbb_UE_ttbarPlus2BDown";
  convertSystematicLabel_["UETUNE_TTB_UP"]       = "CMS_ttHbb_UE_ttbarPlusBUp";
  convertSystematicLabel_["UETUNE_TTB_DOWN"]     = "CMS_ttHbb_UE_ttbarPlusBDown";
  convertSystematicLabel_["UETUNE_TTCC_UP"]      = "CMS_ttHbb_UE_ttbarPlusCCbarUp";
  convertSystematicLabel_["UETUNE_TTCC_DOWN"]    = "CMS_ttHbb_UE_ttbarPlusCCbarDown";
  convertSystematicLabel_["UETUNE_TTOTHER_UP"]   = "CMS_ttHbb_UE_ttbarOtherUp";
  convertSystematicLabel_["UETUNE_TTOTHER_DOWN"] = "CMS_ttHbb_UE_ttbarOtherDown";


  convertSystematicLabel_["PSISRSCALE_TTBB_UP"]      = "CMS_ttHbb_ISR_ttbarPlusBBbarUp";
  convertSystematicLabel_["PSISRSCALE_TTBB_DOWN"]    = "CMS_ttHbb_ISR_ttbarPlusBBbarDown";
  convertSystematicLabel_["PSISRSCALE_TT2B_UP"]      = "CMS_ttHbb_ISR_ttbarPlus2BUp";
  convertSystematicLabel_["PSISRSCALE_TT2B_DOWN"]    = "CMS_ttHbb_ISR_ttbarPlus2BDown";
  convertSystematicLabel_["PSISRSCALE_TTB_UP"]       = "CMS_ttHbb_ISR_ttbarPlusBUp";
  convertSystematicLabel_["PSISRSCALE_TTB_DOWN"]     = "CMS_ttHbb_ISR_ttbarPlusBDown";
  convertSystematicLabel_["PSISRSCALE_TTCC_UP"]      = "CMS_ttHbb_ISR_ttbarPlusCCbarUp";
  convertSystematicLabel_["PSISRSCALE_TTCC_DOWN"]    = "CMS_ttHbb_ISR_ttbarPlusCCbarDown";
  convertSystematicLabel_["PSISRSCALE_TTOTHER_UP"]   = "CMS_ttHbb_ISR_ttbarOtherUp";
  convertSystematicLabel_["PSISRSCALE_TTOTHER_DOWN"] = "CMS_ttHbb_ISR_ttbarOtherDown";

  convertSystematicLabel_["PSFSRSCALE_TTBB_UP"]      = "CMS_ttHbb_FSR_ttbarPlusBBbarUp";
  convertSystematicLabel_["PSFSRSCALE_TTBB_DOWN"]    = "CMS_ttHbb_FSR_ttbarPlusBBbarDown";
  convertSystematicLabel_["PSFSRSCALE_TT2B_UP"]      = "CMS_ttHbb_FSR_ttbarPlus2BUp";
  convertSystematicLabel_["PSFSRSCALE_TT2B_DOWN"]    = "CMS_ttHbb_FSR_ttbarPlus2BDown";
  convertSystematicLabel_["PSFSRSCALE_TTB_UP"]       = "CMS_ttHbb_FSR_ttbarPlusBUp";
  convertSystematicLabel_["PSFSRSCALE_TTB_DOWN"]     = "CMS_ttHbb_FSR_ttbarPlusBDown";
  convertSystematicLabel_["PSFSRSCALE_TTCC_UP"]      = "CMS_ttHbb_FSR_ttbarPlusCCbarUp";
  convertSystematicLabel_["PSFSRSCALE_TTCC_DOWN"]    = "CMS_ttHbb_FSR_ttbarPlusCCbarDown";
  convertSystematicLabel_["PSFSRSCALE_TTOTHER_UP"]   = "CMS_ttHbb_FSR_ttbarOtherUp";
  convertSystematicLabel_["PSFSRSCALE_TTOTHER_DOWN"] = "CMS_ttHbb_FSR_ttbarOtherDown";

  convertSystematicLabel_["PSFSRSCALE_DOWN"] = "CMS_ttHbb_FSRDown";
  convertSystematicLabel_["PSISRSCALE_UP"]   = "CMS_ttHbb_ISRUp";
  convertSystematicLabel_["PSISRSCALE_DOWN"] = "CMS_ttHbb_ISRDown";
  
  convertSystematicLabel_["PSSCALE_TTB_UP"]       = "CMS_ttHbb_PSscale_ttbarPlusBUp";
  convertSystematicLabel_["PSSCALE_TTB_DOWN"]     = "CMS_ttHbb_PSscale_ttbarPlusBDown";
  convertSystematicLabel_["PSSCALE_TT2B_UP"]      = "CMS_ttHbb_PSscale_ttbarPlus2BUp";
  convertSystematicLabel_["PSSCALE_TT2B_DOWN"]    = "CMS_ttHbb_PSscale_ttbarPlus2BDown";
  convertSystematicLabel_["PSSCALE_TTBB_UP"]      = "CMS_ttHbb_PSscale_ttbarPlusBBbarUp";
  convertSystematicLabel_["PSSCALE_TTBB_DOWN"]    = "CMS_ttHbb_PSscale_ttbarPlusBBbarDown";
  convertSystematicLabel_["PSSCALE_TTCC_UP"]      = "CMS_ttHbb_PSscale_ttbarPlusCCbarUp";
  convertSystematicLabel_["PSSCALE_TTCC_DOWN"]    = "CMS_ttHbb_PSscale_ttbarPlusCCbarDown";
  convertSystematicLabel_["PSSCALE_TTOTHER_UP"]   = "CMS_ttHbb_PSscale_ttbarOtherUp";
  convertSystematicLabel_["PSSCALE_TTOTHER_DOWN"] = "CMS_ttHbb_PSscale_ttbarOtherDown";

  convertSystematicLabel_["MATCH_UP"]   = "CMS_ttHbb_HDAMPUp";
  convertSystematicLabel_["MATCH_DOWN"] = "CMS_ttHbb_HDAMPDown";

  convertSystematicLabel_["MATCH_TTBB_UP"]      = "CMS_ttHbb_HDAMP_ttbarPlusBBbarUp";
  convertSystematicLabel_["MATCH_TTBB_DOWN"]    = "CMS_ttHbb_HDAMP_ttbarPlusBBbarDown";
  convertSystematicLabel_["MATCH_TT2B_UP"]      = "CMS_ttHbb_HDAMP_ttbarPlus2BUp";
  convertSystematicLabel_["MATCH_TT2B_DOWN"]    = "CMS_ttHbb_HDAMP_ttbarPlus2BDown";
  convertSystematicLabel_["MATCH_TTB_UP"]       = "CMS_ttHbb_HDAMP_ttbarPlusBUp";
  convertSystematicLabel_["MATCH_TTB_DOWN"]     = "CMS_ttHbb_HDAMP_ttbarPlusBDown";
  convertSystematicLabel_["MATCH_TTCC_UP"]      = "CMS_ttHbb_HDAMP_ttbarPlusCCbarUp";
  convertSystematicLabel_["MATCH_TTCC_DOWN"]    = "CMS_ttHbb_HDAMP_ttbarPlusCCbarDown";
  convertSystematicLabel_["MATCH_TTOTHER_UP"]   = "CMS_ttHbb_HDAMP_ttbarOtherUp";
  convertSystematicLabel_["MATCH_TTOTHER_DOWN"] = "CMS_ttHbb_HDAMP_ttbarOtherDown";


  convertSystematicLabel_["PU_UP"]    = "CMS_ttHbb_PUUp";
  convertSystematicLabel_["PU_DOWN"]  = "CMS_ttHbb_PUDown";

  convertSystematicLabel_["SCALE_UP"]   = "CMS_ttHbb_scaleUp";
  convertSystematicLabel_["SCALE_DOWN"] = "CMS_ttHbb_scaleDown";

  convertSystematicLabel_["MESCALE_UP"]   = "CMS_ttHbb_Q2scaleUp";
  convertSystematicLabel_["MESCALE_DOWN"] = "CMS_ttHbb_Q2scaleDown";

  convertSystematicLabel_["PSSCALE_UP"]   = "CMS_ttHbb_psScaleUp";
  convertSystematicLabel_["PSSCALE_DOWN"] = "CMS_ttHbb_psScaleDown";

  convertSystematicLabel_["PDF_UP"]   = "CMS_ttHbb_PDFUp";//NNPDF
  convertSystematicLabel_["PDF_DOWN"] = "CMS_ttHbb_PDFDown";//NNPDF

  convertSystematicLabel_["TRIG_UP"]   = "CMS_ttHbb_effTrigger_dlUp";
  convertSystematicLabel_["TRIG_DOWN"] = "CMS_ttHbb_effTrigger_dlDown";

  convertSystematicLabel_["LEPT_UP"]   = "CMS_ttHbb_eff_leptonUp";
  convertSystematicLabel_["LEPT_DOWN"] = "CMS_ttHbb_eff_leptonDown";

  convertSystematicLabel_["JESAbsoluteStat_UP"]   = "CMS_scaleAbsoluteStat_jUp";
  convertSystematicLabel_["JESAbsoluteStat_DOWN"] = "CMS_scaleAbsoluteStat_jDown";
  convertSystematicLabel_["JESAbsoluteScale_UP"]   = "CMS_scaleAbsoluteScale_jUp";
  convertSystematicLabel_["JESAbsoluteScale_DOWN"] = "CMS_scaleAbsoluteScale_jDown";
  convertSystematicLabel_["JESAbsoluteFlavMap_UP"]   = "CMS_scaleAbsoluteFlavMap_jUp";
  convertSystematicLabel_["JESAbsoluteFlavMap_DOWN"] = "CMS_scaleAbsoluteFlavMap_jDown";
  convertSystematicLabel_["JESAbsoluteMPFBias_UP"]   = "CMS_scaleAbsoluteMPFBias_jUp";
  convertSystematicLabel_["JESAbsoluteMPFBias_DOWN"] = "CMS_scaleAbsoluteMPFBias_jDown";
  convertSystematicLabel_["JESFragmentation_UP"]   = "CMS_scaleFragmentation_jUp";
  convertSystematicLabel_["JESFragmentation_DOWN"] = "CMS_scaleFragmentation_jDown";
  convertSystematicLabel_["JESSinglePionECAL_UP"]   = "CMS_scaleSinglePionECAL_jUp";
  convertSystematicLabel_["JESSinglePionECAL_DOWN"] = "CMS_scaleSinglePionECAL_jDown";
  convertSystematicLabel_["JESSinglePionHCAL_UP"]   = "CMS_scaleSinglePionHCAL_jUp";
  convertSystematicLabel_["JESSinglePionHCAL_DOWN"] = "CMS_scaleSinglePionHCAL_jDown";
  convertSystematicLabel_["JESFlavorQCD_UP"]   = "CMS_scaleFlavorQCD_jUp";
  convertSystematicLabel_["JESFlavorQCD_DOWN"] = "CMS_scaleFlavorQCD_jDown";
  convertSystematicLabel_["JESRelativeJEREC1_UP"]   = "CMS_scaleRelativeJEREC1_jUp";
  convertSystematicLabel_["JESRelativeJEREC1_DOWN"] = "CMS_scaleRelativeJEREC1_jDown";
  convertSystematicLabel_["JESRelativeJEREC2_UP"]   = "CMS_scaleRelativeJEREC2_jUp";
  convertSystematicLabel_["JESRelativeJEREC2_DOWN"] = "CMS_scaleRelativeJEREC2_jDown";
  convertSystematicLabel_["JESRelativeJERHF_UP"]   = "CMS_scaleRelativeJERHF_jUp";
  convertSystematicLabel_["JESRelativeJERHF_DOWN"] = "CMS_scaleRelativeJERHF_jDown";
  convertSystematicLabel_["JESRelativePtBB_UP"]   = "CMS_scaleRelativePtBB_jUp";
  convertSystematicLabel_["JESRelativePtBB_DOWN"] = "CMS_scaleRelativePtBB_jDown";
  convertSystematicLabel_["JESRelativePtEC1_UP"]   = "CMS_scaleRelativePtEC1_jUp";
  convertSystematicLabel_["JESRelativePtEC1_DOWN"] = "CMS_scaleRelativePtEC1_jDown";
  convertSystematicLabel_["JESRelativePtEC2_UP"]   = "CMS_scaleRelativePtEC2_jUp";
  convertSystematicLabel_["JESRelativePtEC2_DOWN"] = "CMS_scaleRelativePtEC2_jDown";
  convertSystematicLabel_["JESRelativePtHF_UP"]   = "CMS_scaleRelativePtHF_jUp";
  convertSystematicLabel_["JESRelativePtHF_DOWN"] = "CMS_scaleRelativePtHF_jDown";
  convertSystematicLabel_["JESRelativeFSR_UP"]   = "CMS_scaleRelativeFSR_jUp";
  convertSystematicLabel_["JESRelativeFSR_DOWN"] = "CMS_scaleRelativeFSR_jDown";
  convertSystematicLabel_["JESRelativeStatFSR_UP"]   = "CMS_scaleRelativeStatFSR_jUp";
  convertSystematicLabel_["JESRelativeStatFSR_DOWN"] = "CMS_scaleRelativeStatFSR_jDown";
  convertSystematicLabel_["JESRelativeStatEC_UP"]   = "CMS_scaleRelativeStatEC_jUp";
  convertSystematicLabel_["JESRelativeStatEC_DOWN"] = "CMS_scaleRelativeStatEC_jDown";
  convertSystematicLabel_["JESRelativeStatHF_UP"]   = "CMS_scaleRelativeStatHF_jUp";
  convertSystematicLabel_["JESRelativeStatHF_DOWN"] = "CMS_scaleRelativeStatHF_jDown";
  convertSystematicLabel_["JESRelativeBal_UP"]   = "CMS_scaleRelativeBal_jUp";
  convertSystematicLabel_["JESRelativeBal_DOWN"] = "CMS_scaleRelativeBal_jDown";
  convertSystematicLabel_["JESPileUpDataMC_UP"]   = "CMS_scalePileUpDataMC_jUp";
  convertSystematicLabel_["JESPileUpDataMC_DOWN"] = "CMS_scalePileUpDataMC_jDown";
  convertSystematicLabel_["JESPileUpPtRef_UP"]   = "CMS_scalePileUpPtRef_jUp";
  convertSystematicLabel_["JESPileUpPtRef_DOWN"] = "CMS_scalePileUpPtRef_jDown";
  convertSystematicLabel_["JESPileUpPtEC1_UP"]   = "CMS_scalePileUpPtEC1_jUp";
  convertSystematicLabel_["JESPileUpPtEC1_DOWN"] = "CMS_scalePileUpPtEC1_jDown";
  convertSystematicLabel_["JESPileUpPtEC2_UP"]   = "CMS_scalePileUpPtEC2_jUp";
  convertSystematicLabel_["JESPileUpPtEC2_DOWN"] = "CMS_scalePileUpPtEC2_jDown";
  convertSystematicLabel_["JESPileUpPtHF_UP"]   = "CMS_scalePileUpPtHF_jUp";
  convertSystematicLabel_["JESPileUpPtHF_DOWN"] = "CMS_scalePileUpPtHF_jDown";  
  convertSystematicLabel_["JESPileUpPtBB_UP"]   = "CMS_scalePileUpPtBB_jUp";
  convertSystematicLabel_["JESPileUpPtBB_DOWN"] = "CMS_scalePileUpPtBB_jDown";
  convertSystematicLabel_["JESPileUpMuZero_UP"]   = "CMS_scalePileUpMuZero_jUp";
  convertSystematicLabel_["JESPileUpMuZero_DOWN"] = "CMS_scalePileUpMuZero_jDown";
  convertSystematicLabel_["JESPileUpEnvelope_UP"]   = "CMS_scalePileUpEnvelope_jUp";
  convertSystematicLabel_["JESPileUpEnvelope_DOWN"] = "CMS_scalePileUpEnvelope_jDown";
  convertSystematicLabel_["JESFlavorZJet_UP"]   = "CMS_scaleFlavorZJet_jUp";
  convertSystematicLabel_["JESFlavorZJet_DOWN"] = "CMS_scaleFlavorZJet_jDown";
  convertSystematicLabel_["JESFlavorPhotonJet_UP"]   = "CMS_scaleFlavorPhotonJet_jUp";
  convertSystematicLabel_["JESFlavorPhotonJet_DOWN"] = "CMS_scaleFlavorPhotonJet_jDown";
  convertSystematicLabel_["JESFlavorPureGluon_UP"]   = "CMS_scaleFlavorPureGluon_jUp";
  convertSystematicLabel_["JESFlavorPureGluon_DOWN"] = "CMS_scaleFlavorPureGluon_jDown";
  convertSystematicLabel_["JESFlavorPureQuark_UP"]   = "CMS_scaleFlavorPureQuark_jUp";
  convertSystematicLabel_["JESFlavorPureQuark_DOWN"] = "CMS_scaleFlavorPureQuark_jDown";
  convertSystematicLabel_["JESFlavorPureCharm_UP"]   = "CMS_scaleFlavorPureCharm_jUp";
  convertSystematicLabel_["JESFlavorPureCharm_DOWN"] = "CMS_scaleFlavorPureCharm_jDown";
  convertSystematicLabel_["JESFlavorPureBottom_UP"]   = "CMS_scaleFlavorPureBottom_jUp";
  convertSystematicLabel_["JESFlavorPureBottom_DOWN"] = "CMS_scaleFlavorPureBottom_jDown";
  convertSystematicLabel_["JESTimePtEta_UP"]   = "CMS_scaleTimePtEta_jUp";
  convertSystematicLabel_["JESTimePtEta_DOWN"] = "CMS_scaleTimePtEta_jDown";

  // Indicate if systematic uncertainty is of type "shape" or "lnN"
  listOfSystematicsByType_.push_back(std::make_pair("BTAGDISCR_BPURITY","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("BTAGDISCR_LPURITY","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("BTAGDISCR_BSTAT1","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("BTAGDISCR_BSTAT2","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("BTAGDISCR_LSTAT1","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("BTAGDISCR_LSTAT2","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("BTAGDISCR_CERR1","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("BTAGDISCR_CERR2","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JES","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JER","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("SCALE_TTB","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("SCALE_TTBB","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("SCALE_TT2B","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("SCALE_TTCC","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("SCALE_TTOTHER","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("SCALE","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("MESCALE_TTB","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("MESCALE_TTBB","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("MESCALE_TT2B","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("MESCALE_TTCC","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("MESCALE_TTOTHER","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("MESCALE","shape"));

  listOfSystematicsByType_.push_back(std::make_pair("MEFACSCALE","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("MERENSCALE","shape"));

  listOfSystematicsByType_.push_back(std::make_pair("PSFSRSCALE","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("PSFSRSCALE_TTBB","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("PSFSRSCALE_TT2B","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("PSFSRSCALE_TTB","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("PSFSRSCALE_TTCC","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("PSFSRSCALE_TTOTHER","lnN"));

  listOfSystematicsByType_.push_back(std::make_pair("PSISRSCALE","lnN")); 
  listOfSystematicsByType_.push_back(std::make_pair("PSISRSCALE_TTBB","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("PSISRSCALE_TT2B","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("PSISRSCALE_TTB","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("PSISRSCALE_TTCC","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("PSISRSCALE_TTOTHER","lnN"));

  listOfSystematicsByType_.push_back(std::make_pair("UETUNE","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("UETUNE_TTBB","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("UETUNE_TT2B","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("UETUNE_TTB","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("UETUNE_TTCC","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("UETUNE_TTOTHER","lnN"));

  listOfSystematicsByType_.push_back(std::make_pair("MATCH","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("MATCH_TTBB","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("MATCH_TT2B","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("MATCH_TTB","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("MATCH_TTCC","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("MATCH_TTOTHER","lnN"));

  /*
  listOfSystematicsByType_.push_back(std::make_pair("PSFSRSCALE","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("PSISRSCALE","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("UETUNE","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("MATCH","shape"));
  */
  listOfSystematicsByType_.push_back(std::make_pair("PSSCALE_TTB","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("PSSCALE_TTBB","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("PSSCALE_TT2B","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("PSSCALE_TTCC","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("PSSCALE_TTOTHER","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("PSSCALE","shape"));

  listOfSystematicsByType_.push_back(std::make_pair("PU","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("TRIG","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("LEPT","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("PDF","shape"));

  listOfSystematicsByType_.push_back(std::make_pair("JESAbsoluteStat","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESAbsoluteScale","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESAbsoluteFlavMap","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESAbsoluteMPFBias","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESFragmentation","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESSinglePionECAL","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESSinglePionHCAL","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESFlavorQCD","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESRelativeJEREC1","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESRelativeJEREC2","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESRelativeJERHF","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESRelativePtBB","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESRelativePtEC1","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESRelativePtEC2","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESRelativePtHF","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESRelativeFSR","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESRelativeStatFSR","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESRelativeStatEC","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESRelativeStatHF","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESRelativeBal","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESPileUpDataMC","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESPileUpPtRef","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESPileUpPtEC1","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESPileUpPtEC2","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESPileUpPtHF","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESPileUpMuZero","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESPileUpEnvelope","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESPileUpPtBB","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESFlavorZJet","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESFlavorPhotonJet","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESFlavorPureGluon","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESFlavorPureQuark","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESFlavorPureCharm","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESFlavorPureBottom","shape"));
  listOfSystematicsByType_.push_back(std::make_pair("JESTimePtEta","shape"));

  // Those of type lnN
  listOfSystematicsByType_.push_back(std::make_pair("LUMI","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("NORMPDFGG","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("NORMPDFGQ","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("NORMPDFQQ","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("NORMPDFTTH","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("XSEC_TTB","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("XSEC_TTBB","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("XSEC_TT2B","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("XSEC_TTCC","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("XSEC_TTOTHER","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("XSEC_TT","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("XSEC_TTH","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("XSEC_T","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("XSEC_V","lnN"));
  listOfSystematicsByType_.push_back(std::make_pair("XSEC_VV","lnN"));
  //listOfSystematicsByType_.push_back(std::make_pair("PDF","lnN"));
  //listOfSystematicsByType_.push_back(std::make_pair("LEPT","lnN"));
  //listOfSystematicsByType_.push_back(std::make_pair("TRIG","lnN"));

  // Assign uncertainty values for the corresponding uncertainty type for the given processes
  std::string lumiUncert = TString::Format("%f",1+analysisConfig_.general().luminosityUncertainty_).Data();

  valueOfSystematicBasedOnProcess_["LUMI"].push_back(std::make_pair("all",lumiUncert));
  valueOfSystematicBasedOnProcess_["JES"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JER"].push_back(std::make_pair("all","1.000000"));

  valueOfSystematicBasedOnProcess_["BTAGDISCR_BPURITY"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["BTAGDISCR_LPURITY"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["BTAGDISCR_BSTAT1"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["BTAGDISCR_BSTAT2"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["BTAGDISCR_LSTAT1"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["BTAGDISCR_LSTAT2"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["BTAGDISCR_CERR1"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["BTAGDISCR_CERR2"].push_back(std::make_pair("all","1.000000"));

  valueOfSystematicBasedOnProcess_["XSEC_TTB"].push_back(std::make_pair("ttb","1.50"));
  // tt+bb cross section uncertainty reduction to 35% (from 50%) based on this measurement https://arxiv.org/abs/1705.10141
  valueOfSystematicBasedOnProcess_["XSEC_TTBB"].push_back(std::make_pair("ttbb","1.35"));
  valueOfSystematicBasedOnProcess_["XSEC_TT2B"].push_back(std::make_pair("tt2b","1.50"));
  valueOfSystematicBasedOnProcess_["XSEC_TTCC"].push_back(std::make_pair("ttcc","1.50"));
  valueOfSystematicBasedOnProcess_["XSEC_TTOTHER"].push_back(std::make_pair("ttOther","1.030"));
  
  valueOfSystematicBasedOnProcess_["XSEC_TTH"].push_back(std::make_pair("ttH","0.908/1.058"));
  valueOfSystematicBasedOnProcess_["XSEC_TTH"].push_back(std::make_pair("ttHbb","0.908/1.058"));
  valueOfSystematicBasedOnProcess_["XSEC_TTH"].push_back(std::make_pair("ttHcc","0.908/1.058"));
  valueOfSystematicBasedOnProcess_["XSEC_TTH"].push_back(std::make_pair("ttHww","0.908/1.058"));
  valueOfSystematicBasedOnProcess_["XSEC_TTH"].push_back(std::make_pair("ttHzz","0.908/1.058"));
  valueOfSystematicBasedOnProcess_["XSEC_TTH"].push_back(std::make_pair("ttHtautau","0.908/1.058"));
  valueOfSystematicBasedOnProcess_["XSEC_TTH"].push_back(std::make_pair("ttHgluongluon","0.908/1.058"));
  valueOfSystematicBasedOnProcess_["XSEC_TTH"].push_back(std::make_pair("ttHgammagamma","0.908/1.058"));
  valueOfSystematicBasedOnProcess_["XSEC_TTH"].push_back(std::make_pair("ttHzgamma","0.908/1.058"));

  valueOfSystematicBasedOnProcess_["XSEC_TT"].push_back(std::make_pair("ttb","0.96/1.02")); 
  valueOfSystematicBasedOnProcess_["XSEC_TT"].push_back(std::make_pair("ttbb","0.96/1.02"));
  valueOfSystematicBasedOnProcess_["XSEC_TT"].push_back(std::make_pair("tt2b","0.96/1.02"));
  valueOfSystematicBasedOnProcess_["XSEC_TT"].push_back(std::make_pair("ttcc","0.96/1.02"));
  valueOfSystematicBasedOnProcess_["XSEC_TT"].push_back(std::make_pair("ttOther","0.96/1.02"));
  valueOfSystematicBasedOnProcess_["XSEC_TT"].push_back(std::make_pair("ttW","0.88/1.13")); 
  valueOfSystematicBasedOnProcess_["XSEC_TT"].push_back(std::make_pair("ttZ","0.88/1.10")); 
  
  valueOfSystematicBasedOnProcess_["NORMPDFTTH"].push_back(std::make_pair("ttH","1.036"));
  valueOfSystematicBasedOnProcess_["NORMPDFTTH"].push_back(std::make_pair("ttHbb","1.036"));
  valueOfSystematicBasedOnProcess_["NORMPDFTTH"].push_back(std::make_pair("ttHcc","1.036"));
  valueOfSystematicBasedOnProcess_["NORMPDFTTH"].push_back(std::make_pair("ttHww","1.036"));
  valueOfSystematicBasedOnProcess_["NORMPDFTTH"].push_back(std::make_pair("ttHzz","1.036"));
  valueOfSystematicBasedOnProcess_["NORMPDFTTH"].push_back(std::make_pair("ttHtautau","1.036"));
  valueOfSystematicBasedOnProcess_["NORMPDFTTH"].push_back(std::make_pair("ttHgluongluon","1.036"));
  valueOfSystematicBasedOnProcess_["NORMPDFTTH"].push_back(std::make_pair("ttHgammagamma","1.036"));
  valueOfSystematicBasedOnProcess_["NORMPDFTTH"].push_back(std::make_pair("ttHzgamma","1.036"));

  valueOfSystematicBasedOnProcess_["NORMPDFGG"].push_back(std::make_pair("ttb","1.04"));
  valueOfSystematicBasedOnProcess_["NORMPDFGG"].push_back(std::make_pair("ttbb","1.04"));
  valueOfSystematicBasedOnProcess_["NORMPDFGG"].push_back(std::make_pair("tt2b","1.04"));
  valueOfSystematicBasedOnProcess_["NORMPDFGG"].push_back(std::make_pair("ttcc","1.04"));
  valueOfSystematicBasedOnProcess_["NORMPDFGG"].push_back(std::make_pair("ttOther","1.04"));
  valueOfSystematicBasedOnProcess_["NORMPDFGG"].push_back(std::make_pair("ttZ","1.03")); 
  
  valueOfSystematicBasedOnProcess_["NORMPDFGQ"].push_back(std::make_pair("singleTop","1.03"));
    
  valueOfSystematicBasedOnProcess_["NORMPDFQQ"].push_back(std::make_pair("wlnu","1.04")); 
  valueOfSystematicBasedOnProcess_["NORMPDFQQ"].push_back(std::make_pair("dy","1.04")); 
  valueOfSystematicBasedOnProcess_["NORMPDFQQ"].push_back(std::make_pair("ttW","1.02")); 
  valueOfSystematicBasedOnProcess_["NORMPDFQQ"].push_back(std::make_pair("diboson","1.02"));  
      
  valueOfSystematicBasedOnProcess_["XSEC_T"].push_back(std::make_pair("singleTop","0.98/1.03")); 
  
  valueOfSystematicBasedOnProcess_["XSEC_V"].push_back(std::make_pair("wlnu","1.01"));
  valueOfSystematicBasedOnProcess_["XSEC_V"].push_back(std::make_pair("dy","1.01"));
  
  valueOfSystematicBasedOnProcess_["XSEC_VV"].push_back(std::make_pair("diboson","1.02"));
  
  valueOfSystematicBasedOnProcess_["SCALE"].push_back(std::make_pair("ttb","1.000000"));
  valueOfSystematicBasedOnProcess_["SCALE"].push_back(std::make_pair("ttbb","1.000000"));
  valueOfSystematicBasedOnProcess_["SCALE"].push_back(std::make_pair("tt2b","1.000000"));
  valueOfSystematicBasedOnProcess_["SCALE"].push_back(std::make_pair("ttcc","1.000000"));
  valueOfSystematicBasedOnProcess_["SCALE"].push_back(std::make_pair("ttOther","1.000000"));

  valueOfSystematicBasedOnProcess_["SCALE_TTB"].push_back(std::make_pair("ttb","1.000000"));
  valueOfSystematicBasedOnProcess_["SCALE_TTBB"].push_back(std::make_pair("ttbb","1.000000"));
  valueOfSystematicBasedOnProcess_["SCALE_TT2B"].push_back(std::make_pair("tt2b","1.000000"));
  valueOfSystematicBasedOnProcess_["SCALE_TTCC"].push_back(std::make_pair("ttcc","1.000000"));
  valueOfSystematicBasedOnProcess_["SCALE_TTOTHER"].push_back(std::make_pair("ttOther","1.000000"));

  valueOfSystematicBasedOnProcess_["MESCALE"].push_back(std::make_pair("ttb","1.000000"));
  valueOfSystematicBasedOnProcess_["MESCALE"].push_back(std::make_pair("ttbb","1.000000"));
  valueOfSystematicBasedOnProcess_["MESCALE"].push_back(std::make_pair("tt2b","1.000000"));
  valueOfSystematicBasedOnProcess_["MESCALE"].push_back(std::make_pair("ttcc","1.000000"));
  valueOfSystematicBasedOnProcess_["MESCALE"].push_back(std::make_pair("ttOther","1.000000"));

  valueOfSystematicBasedOnProcess_["MESCALE_TTB"].push_back(std::make_pair("ttb","1.000000"));
  valueOfSystematicBasedOnProcess_["MESCALE_TTBB"].push_back(std::make_pair("ttbb","1.000000"));
  valueOfSystematicBasedOnProcess_["MESCALE_TT2B"].push_back(std::make_pair("tt2b","1.000000"));
  valueOfSystematicBasedOnProcess_["MESCALE_TTCC"].push_back(std::make_pair("ttcc","1.000000"));
  valueOfSystematicBasedOnProcess_["MESCALE_TTOTHER"].push_back(std::make_pair("ttOther","1.000000"));

  valueOfSystematicBasedOnProcess_["MEFACSCALE"].push_back(std::make_pair("ttb","1.000000"));
  valueOfSystematicBasedOnProcess_["MEFACSCALE"].push_back(std::make_pair("ttbb","1.000000"));
  valueOfSystematicBasedOnProcess_["MEFACSCALE"].push_back(std::make_pair("tt2b","1.000000"));
  valueOfSystematicBasedOnProcess_["MEFACSCALE"].push_back(std::make_pair("ttcc","1.000000"));
  valueOfSystematicBasedOnProcess_["MEFACSCALE"].push_back(std::make_pair("ttOther","1.000000"));

  valueOfSystematicBasedOnProcess_["MERENSCALE"].push_back(std::make_pair("ttb","1.000000"));
  valueOfSystematicBasedOnProcess_["MERENSCALE"].push_back(std::make_pair("ttbb","1.000000"));
  valueOfSystematicBasedOnProcess_["MERENSCALE"].push_back(std::make_pair("tt2b","1.000000"));
  valueOfSystematicBasedOnProcess_["MERENSCALE"].push_back(std::make_pair("ttcc","1.000000"));
  valueOfSystematicBasedOnProcess_["MERENSCALE"].push_back(std::make_pair("ttOther","1.000000"));

  // NJET >=3
  //version 2
  if(false) {
  valueOfSystematicBasedOnProcess_["MATCH"].push_back(std::make_pair("ttb","0.96/1.04"));
  valueOfSystematicBasedOnProcess_["MATCH"].push_back(std::make_pair("ttbb","0.86/1.16"));
  valueOfSystematicBasedOnProcess_["MATCH"].push_back(std::make_pair("tt2b","0.94/1.06"));
  valueOfSystematicBasedOnProcess_["MATCH"].push_back(std::make_pair("ttcc","0.98/1.02"));
  valueOfSystematicBasedOnProcess_["MATCH"].push_back(std::make_pair("ttOther","1.02/0.984"));

  valueOfSystematicBasedOnProcess_["MATCH_TTBB"].push_back(std::make_pair("ttbb","0.856/1.056"));
  valueOfSystematicBasedOnProcess_["MATCH_TTB"].push_back(std::make_pair("ttb","0.978/1.049"));
  valueOfSystematicBasedOnProcess_["MATCH_TT2B"].push_back(std::make_pair("tt2b","0.959/1.042"));
  valueOfSystematicBasedOnProcess_["MATCH_TTCC"].push_back(std::make_pair("ttcc","0.987/1.013"));
  valueOfSystematicBasedOnProcess_["MATCH_TTOTHER"].push_back(std::make_pair("ttOther","1.006/0.994"));

  valueOfSystematicBasedOnProcess_["PSFSRSCALE"].push_back(std::make_pair("ttb","0.95/1.07"));
  valueOfSystematicBasedOnProcess_["PSFSRSCALE"].push_back(std::make_pair("ttbb","0.92/1.12"));
  valueOfSystematicBasedOnProcess_["PSFSRSCALE"].push_back(std::make_pair("tt2b","0.91/1.11"));
  valueOfSystematicBasedOnProcess_["PSFSRSCALE"].push_back(std::make_pair("ttcc","0.91/1.07"));
  valueOfSystematicBasedOnProcess_["PSFSRSCALE"].push_back(std::make_pair("ttOther","1.09/0.85")); 

  valueOfSystematicBasedOnProcess_["PSFSRSCALE_TTBB"].push_back(std::make_pair("ttbb","0.932/1.12"));
  valueOfSystematicBasedOnProcess_["PSFSRSCALE_TTB"].push_back(std::make_pair("ttb","0.967/1.078"));
  valueOfSystematicBasedOnProcess_["PSFSRSCALE_TT2B"].push_back(std::make_pair("tt2b","0.926/1.129"));
  valueOfSystematicBasedOnProcess_["PSFSRSCALE_TTCC"].push_back(std::make_pair("ttcc","0.929/1.096"));
  valueOfSystematicBasedOnProcess_["PSFSRSCALE_TTOTHER"].push_back(std::make_pair("ttOther","1.102/0.86"));

  valueOfSystematicBasedOnProcess_["PSISRSCALE"].push_back(std::make_pair("ttbb","1.1/0.9"));
  valueOfSystematicBasedOnProcess_["PSISRSCALE"].push_back(std::make_pair("ttb","1.04/0.95"));
  valueOfSystematicBasedOnProcess_["PSISRSCALE"].push_back(std::make_pair("tt2b","1.08/0.91"));
  valueOfSystematicBasedOnProcess_["PSISRSCALE"].push_back(std::make_pair("ttcc","1.02/0.97"));
  valueOfSystematicBasedOnProcess_["PSISRSCALE"].push_back(std::make_pair("ttOther","1.015/0.96"));
  
  valueOfSystematicBasedOnProcess_["PSISRSCALE_TTB"].push_back(std::make_pair("ttb","1.055/0.945"));
  valueOfSystematicBasedOnProcess_["PSISRSCALE_TTBB"].push_back(std::make_pair("ttbb","1.022/0.953"));
  valueOfSystematicBasedOnProcess_["PSISRSCALE_TT2B"].push_back(std::make_pair("tt2b","1.068/0.959"));
  valueOfSystematicBasedOnProcess_["PSISRSCALE_TTCC"].push_back(std::make_pair("ttcc","1.035/0.987"));
  valueOfSystematicBasedOnProcess_["PSISRSCALE_TTOTHER"].push_back(std::make_pair("ttOther","1.027/0.97"));

  valueOfSystematicBasedOnProcess_["UETUNE"].push_back(std::make_pair("ttb","1.04/0.96"));
  valueOfSystematicBasedOnProcess_["UETUNE"].push_back(std::make_pair("ttbb","1.10/0.90"));
  valueOfSystematicBasedOnProcess_["UETUNE"].push_back(std::make_pair("tt2b","1.08/0.92"));
  valueOfSystematicBasedOnProcess_["UETUNE"].push_back(std::make_pair("ttcc","1.02/0.98"));
  valueOfSystematicBasedOnProcess_["UETUNE"].push_back(std::make_pair("ttOther","1.01/0.99"));
  
  valueOfSystematicBasedOnProcess_["UETUNE_TTBB"].push_back(std::make_pair("ttbb", "0.943/1.055"));
  valueOfSystematicBasedOnProcess_["UETUNE_TTB"].push_back(std::make_pair("ttb","0.978/1.022"));
  valueOfSystematicBasedOnProcess_["UETUNE_TT2B"].push_back(std::make_pair("tt2b","1.063/0.959"));
  valueOfSystematicBasedOnProcess_["UETUNE_TTCC"].push_back(std::make_pair("ttcc","0.987/1.013"));
  valueOfSystematicBasedOnProcess_["UETUNE_TTOTHER"].push_back(std::make_pair("ttOther","1.009/0.996"));
  }

  // NJET >=4
  // version 2
  if(true) {
  valueOfSystematicBasedOnProcess_["MATCH"].push_back(std::make_pair("ttb","0.91/1.10"));
  valueOfSystematicBasedOnProcess_["MATCH"].push_back(std::make_pair("ttbb","0.92/1.08"));
  valueOfSystematicBasedOnProcess_["MATCH"].push_back(std::make_pair("tt2b","0.95/1.05"));
  valueOfSystematicBasedOnProcess_["MATCH"].push_back(std::make_pair("ttcc","0.94/1.06"));
  valueOfSystematicBasedOnProcess_["MATCH"].push_back(std::make_pair("ttOther","0.96/1.04"));

  valueOfSystematicBasedOnProcess_["MATCH_TTBB"].push_back(std::make_pair("ttbb","0.921/1.052"));
  valueOfSystematicBasedOnProcess_["MATCH_TTB"].push_back(std::make_pair("ttb","0.932/1.024"));
  valueOfSystematicBasedOnProcess_["MATCH_TT2B"].push_back(std::make_pair("tt2b","0.966/1.049"));
  valueOfSystematicBasedOnProcess_["MATCH_TTCC"].push_back(std::make_pair("ttcc","0.957/1.023"));
  valueOfSystematicBasedOnProcess_["MATCH_TTOTHER"].push_back(std::make_pair("ttOther","0.97/1.02"));

  valueOfSystematicBasedOnProcess_["PSFSRSCALE"].push_back(std::make_pair("ttb","0.94/1.04"));
  valueOfSystematicBasedOnProcess_["PSFSRSCALE"].push_back(std::make_pair("ttbb","0.95/1.05"));
  valueOfSystematicBasedOnProcess_["PSFSRSCALE"].push_back(std::make_pair("tt2b","0.91/1.1"));
  valueOfSystematicBasedOnProcess_["PSFSRSCALE"].push_back(std::make_pair("ttcc","0.92/1.07"));
  valueOfSystematicBasedOnProcess_["PSFSRSCALE"].push_back(std::make_pair("ttOther","1.12/0.81"));

  valueOfSystematicBasedOnProcess_["PSFSRSCALE_TTBB"].push_back(std::make_pair("ttbb","0.945/1.042"));
  valueOfSystematicBasedOnProcess_["PSFSRSCALE_TTB"].push_back(std::make_pair("ttb","0.947/1.055"));
  valueOfSystematicBasedOnProcess_["PSFSRSCALE_TT2B"].push_back(std::make_pair("tt2b","0.925/1.036"));
  valueOfSystematicBasedOnProcess_["PSFSRSCALE_TTCC"].push_back(std::make_pair("ttcc","0.931/1.085"));
  valueOfSystematicBasedOnProcess_["PSFSRSCALE_TTOTHER"].push_back(std::make_pair("ttOther","1.137/0.819"));

  valueOfSystematicBasedOnProcess_["PSISRSCALE"].push_back(std::make_pair("ttb","0.92/1.05"));
  valueOfSystematicBasedOnProcess_["PSISRSCALE"].push_back(std::make_pair("ttbb","0.95/1.06"));
  valueOfSystematicBasedOnProcess_["PSISRSCALE"].push_back(std::make_pair("tt2b","0.91/1.1"));
  valueOfSystematicBasedOnProcess_["PSISRSCALE"].push_back(std::make_pair("ttcc","0.93/1.04"));
  valueOfSystematicBasedOnProcess_["PSISRSCALE"].push_back(std::make_pair("ttOther","0.966/1.13"));

  valueOfSystematicBasedOnProcess_["PSISRSCALE_TTBB"].push_back(std::make_pair("ttbb","0.949/1.055"));
  valueOfSystematicBasedOnProcess_["PSISRSCALE_TTB"].push_back(std::make_pair("ttb","0.933/1.065"));
  valueOfSystematicBasedOnProcess_["PSISRSCALE_TT2B"].push_back(std::make_pair("tt2b","0.923/1.109"));
  valueOfSystematicBasedOnProcess_["PSISRSCALE_TTCC"].push_back(std::make_pair("ttcc","0.945/1.055"));
  valueOfSystematicBasedOnProcess_["PSISRSCALE_TTOTHER"].push_back(std::make_pair("ttOther","0.975/1.044"));

  valueOfSystematicBasedOnProcess_["UETUNE"].push_back(std::make_pair("ttb","1.02/0.98"));
  valueOfSystematicBasedOnProcess_["UETUNE"].push_back(std::make_pair("ttbb","1.03/0.97"));
  valueOfSystematicBasedOnProcess_["UETUNE"].push_back(std::make_pair("tt2b","1.06/0.94"));
  valueOfSystematicBasedOnProcess_["UETUNE"].push_back(std::make_pair("ttcc","1.04/0.96"));
  valueOfSystematicBasedOnProcess_["UETUNE"].push_back(std::make_pair("ttOther","1.01/0.99"));
  
  valueOfSystematicBasedOnProcess_["UETUNE_TTBB"].push_back(std::make_pair("ttbb","1.021/0.979"));
  valueOfSystematicBasedOnProcess_["UETUNE_TTB"].push_back(std::make_pair("ttb","0.984/1.016"));
  valueOfSystematicBasedOnProcess_["UETUNE_TT2B"].push_back(std::make_pair("tt2b","1.025/0.954"));
  valueOfSystematicBasedOnProcess_["UETUNE_TTCC"].push_back(std::make_pair("ttcc","0.981/1.018"));
  valueOfSystematicBasedOnProcess_["UETUNE_TTOTHER"].push_back(std::make_pair("ttOther","1.005/0.996"));
  }

  // Shape
  /*
  valueOfSystematicBasedOnProcess_["PSFSRSCALE"].push_back(std::make_pair("ttb","1.000000"));
  valueOfSystematicBasedOnProcess_["PSFSRSCALE"].push_back(std::make_pair("ttbb","1.000000"));
  valueOfSystematicBasedOnProcess_["PSFSRSCALE"].push_back(std::make_pair("tt2b","1.000000"));
  valueOfSystematicBasedOnProcess_["PSFSRSCALE"].push_back(std::make_pair("ttcc","1.000000"));
  valueOfSystematicBasedOnProcess_["PSFSRSCALE"].push_back(std::make_pair("ttOther","1.000000"));

  valueOfSystematicBasedOnProcess_["PSISRSCALE"].push_back(std::make_pair("ttb","1.000000"));
  valueOfSystematicBasedOnProcess_["PSISRSCALE"].push_back(std::make_pair("ttbb","1.000000"));
  valueOfSystematicBasedOnProcess_["PSISRSCALE"].push_back(std::make_pair("tt2b","1.000000"));
  valueOfSystematicBasedOnProcess_["PSISRSCALE"].push_back(std::make_pair("ttcc","1.000000"));
  valueOfSystematicBasedOnProcess_["PSISRSCALE"].push_back(std::make_pair("ttOther","1.000000"));

  valueOfSystematicBasedOnProcess_["UETUNE"].push_back(std::make_pair("ttb","1.000000"));
  valueOfSystematicBasedOnProcess_["UETUNE"].push_back(std::make_pair("ttbb","1.000000"));
  valueOfSystematicBasedOnProcess_["UETUNE"].push_back(std::make_pair("tt2b","1.000000"));
  valueOfSystematicBasedOnProcess_["UETUNE"].push_back(std::make_pair("ttcc","1.000000"));
  valueOfSystematicBasedOnProcess_["UETUNE"].push_back(std::make_pair("ttOther","1.000000"));

  valueOfSystematicBasedOnProcess_["MATCH"].push_back(std::make_pair("ttb","1.000000"));  
  valueOfSystematicBasedOnProcess_["MATCH"].push_back(std::make_pair("ttbb","1.000000"));
  valueOfSystematicBasedOnProcess_["MATCH"].push_back(std::make_pair("tt2b","1.000000"));
  valueOfSystematicBasedOnProcess_["MATCH"].push_back(std::make_pair("ttcc","1.000000")); 
  valueOfSystematicBasedOnProcess_["MATCH"].push_back(std::make_pair("ttOther","1.000000")); 
  */


  valueOfSystematicBasedOnProcess_["PSSCALE"].push_back(std::make_pair("ttb","1.000000"));
  valueOfSystematicBasedOnProcess_["PSSCALE"].push_back(std::make_pair("ttbb","1.000000"));
  valueOfSystematicBasedOnProcess_["PSSCALE"].push_back(std::make_pair("tt2b","1.000000"));
  valueOfSystematicBasedOnProcess_["PSSCALE"].push_back(std::make_pair("ttcc","1.000000"));
  valueOfSystematicBasedOnProcess_["PSSCALE"].push_back(std::make_pair("ttOther","1.000000"));

  valueOfSystematicBasedOnProcess_["PSSCALE_TTB"].push_back(std::make_pair("ttb","1.000000"));
  valueOfSystematicBasedOnProcess_["PSSCALE_TTBB"].push_back(std::make_pair("ttbb","1.000000"));
  valueOfSystematicBasedOnProcess_["PSSCALE_TT2B"].push_back(std::make_pair("tt2b","1.000000"));
  valueOfSystematicBasedOnProcess_["PSSCALE_TTCC"].push_back(std::make_pair("ttcc","1.000000"));
  valueOfSystematicBasedOnProcess_["PSSCALE_TTOTHER"].push_back(std::make_pair("ttOther","1.000000"));
  
  //valueOfSystematicBasedOnProcess_["PDF"].push_back(std::make_pair("ttb","1.036"));
  //valueOfSystematicBasedOnProcess_["PDF"].push_back(std::make_pair("ttbb","1.036"));
  //valueOfSystematicBasedOnProcess_["PDF"].push_back(std::make_pair("tt2b","1.036"));
  //valueOfSystematicBasedOnProcess_["PDF"].push_back(std::make_pair("ttcc","1.036"));
  //valueOfSystematicBasedOnProcess_["PDF"].push_back(std::make_pair("ttOther","1.036"));

  // PDF applied as shape: ttH signal
  valueOfSystematicBasedOnProcess_["PDF"].push_back(std::make_pair("ttHbb","1.000000"));
  valueOfSystematicBasedOnProcess_["PDF"].push_back(std::make_pair("ttHcc","1.000000"));
  valueOfSystematicBasedOnProcess_["PDF"].push_back(std::make_pair("ttHgluongluon","1.000000"));
  valueOfSystematicBasedOnProcess_["PDF"].push_back(std::make_pair("ttHgammagamma","1.000000"));
  valueOfSystematicBasedOnProcess_["PDF"].push_back(std::make_pair("ttHtautau","1.000000"));
  valueOfSystematicBasedOnProcess_["PDF"].push_back(std::make_pair("ttHww","1.000000"));
  valueOfSystematicBasedOnProcess_["PDF"].push_back(std::make_pair("ttHzz","1.000000"));
  valueOfSystematicBasedOnProcess_["PDF"].push_back(std::make_pair("ttHzgamma","1.000000"));

  // PDF applied as shape: ttJets background
  valueOfSystematicBasedOnProcess_["PDF"].push_back(std::make_pair("ttb","1.000000"));
  valueOfSystematicBasedOnProcess_["PDF"].push_back(std::make_pair("ttbb","1.000000"));
  valueOfSystematicBasedOnProcess_["PDF"].push_back(std::make_pair("tt2b","1.000000"));
  valueOfSystematicBasedOnProcess_["PDF"].push_back(std::make_pair("ttcc","1.000000"));
  valueOfSystematicBasedOnProcess_["PDF"].push_back(std::make_pair("ttOther","1.000000"));

  valueOfSystematicBasedOnProcess_["PU"].push_back(std::make_pair("all","1.000000"));

  valueOfSystematicBasedOnProcess_["TRIG"].push_back(std::make_pair("all","1.000000"));
  //valueOfSystematicBasedOnProcess_["TRIG"].push_back(std::make_pair("all","1.030000"));

  valueOfSystematicBasedOnProcess_["LEPT"].push_back(std::make_pair("all","1.000000"));
  //valueOfSystematicBasedOnProcess_["LEPT"].push_back(std::make_pair("all","1.040000"));

  valueOfSystematicBasedOnProcess_["JESAbsoluteStat"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESAbsoluteScale"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESAbsoluteFlavMap"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESAbsoluteMPFBias"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESFragmentation"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESSinglePionECAL"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESSinglePionHCAL"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESFlavorQCD"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESRelativeJEREC1"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESRelativeJEREC2"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESRelativeJERHF"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESRelativePtBB"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESRelativePtEC1"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESRelativePtEC2"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESRelativePtHF"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESRelativeFSR"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESRelativeStatFSR"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESRelativeStatEC"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESRelativeStatHF"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESRelativeBal"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESPileUpDataMC"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESPileUpPtRef"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESPileUpPtEC1"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESPileUpPtEC2"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESPileUpPtHF"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESPileUpPtBB"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESPileUpMuZero"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESPileUpEnvelope"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESFlavorZJet"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESFlavorPhotonJet"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESFlavorPureGluon"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESFlavorPureQuark"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESFlavorPureCharm"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESFlavorPureBottom"].push_back(std::make_pair("all","1.000000"));
  valueOfSystematicBasedOnProcess_["JESTimePtEta"].push_back(std::make_pair("all","1.000000"));

}


void DatacardMaker::writeDatacards()
{
  std::map<Channel::Channel, std::map<Systematic::Systematic, std::map<std::string,std::string > > > m_inputRootFileNames; 

  // Set up systematic variations
  std::vector<Systematic::Variation> v_variation;
  v_variation.push_back(Systematic::Variation::up);
  v_variation.push_back(Systematic::Variation::down);

  for(auto fileNamesConfig : fileNames_) {    

    for(Channel::Channel channel : v_channel_) { 

      for(Systematic::Systematic systematic : v_systematic_) {   
      
        // Constructing a neutral systematic for which two variations will be stored             
        Systematic::Systematic systematicToStore = Systematic::Systematic(systematic.type(), Systematic::undefinedVariation, systematic.variationNumber());
        TString inputFileListName;

        if(systematic.type() == Systematic::scale_ttb       || systematic.type() == Systematic::scale_ttbb  || systematic.type() == Systematic::scale_tt2b
           || systematic.type() == Systematic::scale_ttcc   || systematic.type() == Systematic::scale_ttother) {
          TString tmp = systematic.name().Contains("_UP") ? "SCALE_UP" : "SCALE_DOWN";
          inputFileListName = fileList_base_+"/"+"HistoFileList_"+tmp+"_"+Channel::convert(channel)+".txt";
        }
        else if(systematic.type() == Systematic::meScale_ttb     || systematic.type() == Systematic::meScale_ttbb  || systematic.type() == Systematic::meScale_tt2b
                || systematic.type() == Systematic::meScale_ttcc || systematic.type() == Systematic::meScale_ttother) {
          TString tmp = systematic.name().Contains("_UP") ? "MESCALE_UP" : "MESCALE_DOWN";
          inputFileListName = fileList_base_+"/"+"HistoFileList_"+tmp+"_"+Channel::convert(channel)+".txt";
        }
        else if(systematic.type() == Systematic::psScale_ttb     || systematic.type() == Systematic::psScale_ttbb  || systematic.type() == Systematic::psScale_tt2b || systematic.type() == Systematic::psScale_ttcc || systematic.type() == Systematic::psScale_ttother) {
          TString tmp = systematic.name().Contains("_UP") ? "PSSCALE_UP" : "PSSCALE_DOWN";
          inputFileListName = fileList_base_+"/"+"HistoFileList_"+tmp+"_"+Channel::convert(channel)+".txt";
        }
        else {
          inputFileListName = fileList_base_+"/"+"HistoFileList_"+systematic.name()+"_"+Channel::convert(channel)+".txt";
        }

        std::ifstream file(inputFileListName.Data());
      
        std::string lnN("lnN");
        TString lnNTypeUncertain = systematic.name().Contains("_UP") ? systematic.name().ReplaceAll("_UP","") : systematic.name().ReplaceAll("_DOWN","");

        // Check if systematic type is lnN (i.e. rate) if so skip
        if(std::find_if(listOfSystematicsByType_.begin(),listOfSystematicsByType_.end(),compare(lnNTypeUncertain.Data(),lnN)) != listOfSystematicsByType_.end()) {
          continue;
        }       

        if(!file) {
          std::cerr << "### File list not found: " << inputFileListName << " Breaking...\n\n";
          exit(1);
        }
        // Reading each line of the file corresponding to a separate histogram
        std::string line_;
        
        while(std::getline(file, line_)) {   
          std::string fileNameWithDirPath(line_);

          // Extracting the histogram name
          std::size_t found = std::string(fileNameWithDirPath).find_last_of("/\\");
          std::string fileName(std::string(fileNameWithDirPath).substr(found+1));

          if(TString(fileName).Contains(fileNamesConfig)) {
            m_inputRootFileNames[channel][systematic][fileName] = fileNameWithDirPath;
            inputFileLists_ = &m_inputRootFileNames;
          }
        }
      }
   
      TString labelString;
      TObjString labels;

      // Fill maps and meta information
      for(auto convertedSysLabel : convertSystematicLabel_) {
        labelString  += TString::Format("%s\t%s\n", (convertedSysLabel.first).c_str(), (convertedSysLabel.second).c_str());
      }

      TObjArray* token  = TString(fileNamesConfig).Tokenize("_");
      TString filename  = TString(fileNamesConfig);
      TString eventCategory  = ((TObjString*)token->At(token->GetLast()))->GetString();

      // Begin creating directory structure
      TString path("");
    
      // Create all subdirectories contained in output baseDir
      TObjArray* a_currentDir = TString(outputBaseDir_).Tokenize("/");
      for(Int_t iCurrentDir = 0; iCurrentDir < a_currentDir->GetEntriesFast(); ++iCurrentDir){
        const TString& currentDir = a_currentDir->At(iCurrentDir)->GetName();
        path.Append(currentDir);
        path.Append("/");
        gSystem->MakeDirectory(path);
      }

      // Create subdirectories for channel
      path.Append(Channel::convert(channel));
      path.Append("/");
      gSystem->MakeDirectory(path);

      outputDirDatacard_ = path;

      // Set datacard outfile name with partial directory path
      eventFileString_   = outputDirDatacard_+"ttH_hbb_13TeV_"+mapOfCategories_[eventCategory.Data()];
      datacardName_      = eventFileString_+".txt";

      // Create directory for output root file storage
      gSystem->MakeDirectory(TString(outputDirDatacard_+"common"));
      
      std::cout << "\nOpening file: " << datacardName_ << std::endl;

      // Update (i.e. append) to output file      
      if(!outputFile_) {
        outputFile_ = new TFile(TString(outputDirDatacard_+outputFileName_),"NEW");
      }
      else {
        outputFile_ = new TFile(TString(outputDirDatacard_+outputFileName_),"UPDATE"); 
      }

      // Write histograms to correct event category directory
      if(outputFile_->GetDirectory(TString(mapOfCategories_[eventCategory.Data()]+"_"+observableType_))) {
        outputFile_->cd(TString(mapOfCategories_[eventCategory.Data()]+"_"+observableType_));
      }                                                                                 
      else {
        outputFile_->mkdir(TString(mapOfCategories_[eventCategory.Data()]+"_"+observableType_), TString(mapOfCategories_[eventCategory.Data()]+"_"+observableType_));
        outputFile_->cd(TString(mapOfCategories_[eventCategory.Data()]+"_"+observableType_));
      }

      // Fill meta information on what mapping labeling used
      labels = labelString.Data();
      labels.Write("MetaInfoLabelConvert",TObject::kOverwrite);
      outputFile_->Close();
    }

    // Start writing datacard
    writeHeader(fileNamesConfig);
    extractYields(fileNamesConfig);

    if(!datacard_.is_open())
      datacard_.open(datacardName_, std::ios::out | std::ios::app);

    datacard_ << "#Source of uncertainty\t\t pdf\t\t";

    for(auto process : processNames_) {
      if(process != "allmc" && process != "data") {
        datacard_ << TString::Format("%-8s\t", convertSampleNames_[process].c_str());
      }
    }
    datacard_ << std::endl;

    writeYields(fileNamesConfig);
    if(addSystematicUncertainty_) writeSystematicUncertainties();
    if(addStatisticalUncertainty_) writeStatisticalUncertainties(fileNamesConfig);

    bool append_to_datacard_systematic_group_labels(false);

    if(!datacard_.is_open())
      datacard_.open(datacardName_, std::ios::out | std::ios::app);

    if(append_to_datacard_systematic_group_labels) {

      datacard_ << "\n"<< std::endl;
      datacard_ << "exp group = lumi_13TeV_2016 CMS_res_j CMS_ttHbb_effTrigger_dl CMS_scaleAbsoluteMPFBias_j CMS_scaleAbsoluteScale_j CMS_scaleAbsoluteStat_j CMS_scaleFlavorQCD_j CMS_scaleFragmentation_j CMS_scalePileUpDataMC_j CMS_scalePileUpPtBB_j CMS_scalePileUpPtEC1_j CMS_scalePileUpPtEC2_j CMS_scalePileUpPtHF_j CMS_scalePileUpPtRef_j CMS_scaleRelativeBal_j CMS_scaleRelativeFSR_j CMS_scaleRelativeJEREC1_j CMS_scaleRelativeJEREC2_j CMS_scaleRelativeJERHF_j CMS_scaleRelativePtBB_j CMS_scaleRelativePtEC1_j CMS_scaleRelativePtEC2_j CMS_scaleRelativePtHF_j CMS_scaleRelativeStatEC_j CMS_scaleRelativeStatFSR_j CMS_scaleRelativeStatHF_j CMS_scaleSinglePionECAL_j CMS_scaleSinglePionHCAL_j CMS_scaleTimePtEta_j CMS_btag_lf CMS_btag_hf CMS_btag_hfstats1 CMS_btag_hfstats2 CMS_btag_cferr1 CMS_btag_cferr2 CMS_btag_lfstats1 CMS_btag_lfstats2 CMS_ttHbb_PU" << std::endl;
      datacard_ << "syst group = QCDscale_V QCDscale_VV QCDscale_singlet QCDscale_ttH QCDscale_ttbar bgnorm_ttbarPlus2B bgnorm_ttbarPlusB bgnorm_ttbarPlusBBbar bgnorm_ttbarPlusCCbar pdf_gg pdf_qg pdf_qqbar pdf_Higgs_ttH lumi_13TeV_2016 CMS_res_j CMS_scaleAbsoluteMPFBias_j CMS_scaleAbsoluteScale_j CMS_scaleAbsoluteStat_j CMS_scaleFlavorQCD_j CMS_scaleFragmentation_j CMS_scalePileUpDataMC_j CMS_scalePileUpPtBB_j CMS_scalePileUpPtEC1_j CMS_scalePileUpPtEC2_j CMS_scalePileUpPtHF_j CMS_scalePileUpPtRef_j CMS_scaleRelativeBal_j CMS_scaleRelativeFSR_j CMS_scaleRelativeJEREC1_j CMS_scaleRelativeJEREC2_j CMS_scaleRelativeJERHF_j CMS_scaleRelativePtBB_j CMS_scaleRelativePtEC1_j CMS_scaleRelativePtEC2_j CMS_scaleRelativePtHF_j CMS_scaleRelativeStatEC_j CMS_scaleRelativeStatFSR_j CMS_scaleRelativeStatHF_j CMS_scaleSinglePionECAL_j CMS_scaleSinglePionHCAL_j CMS_scaleTimePtEta_j CMS_btag_lf CMS_btag_hf CMS_btag_hfstats1 CMS_btag_hfstats2 CMS_btag_cferr1 CMS_btag_cferr2 CMS_btag_lfstats1 CMS_btag_lfstats2 CMS_ttHbb_PU CMS_ttHbb_PDF CMS_ttHbb_scaleMuF CMS_ttHbb_scaleMuR CMS_ttHbb_UE_ttbarPlusBBbar CMS_ttHbb_UE_ttbarPlus2B CMS_ttHbb_UE_ttbarPlusB CMS_ttHbb_UE_ttbarPlusCCbar CMS_ttHbb_UE_ttbarOther CMS_ttHbb_ISR_ttbarPlusBBbar CMS_ttHbb_ISR_ttbarPlus2B CMS_ttHbb_ISR_ttbarPlusB CMS_ttHbb_ISR_ttbarPlusCCbar CMS_ttHbb_ISR_ttbarOther CMS_ttHbb_FSR_ttbarPlusBBbar CMS_ttHbb_FSR_ttbarPlus2B CMS_ttHbb_FSR_ttbarPlusB CMS_ttHbb_FSR_ttbarPlusCCbar CMS_ttHbb_FSR_ttbarOther CMS_ttHbb_HDAMP_ttbarPlusBBbar CMS_ttHbb_HDAMP_ttbarPlus2B CMS_ttHbb_HDAMP_ttbarPlusB CMS_ttHbb_HDAMP_ttbarPlusCCbar CMS_ttHbb_HDAMP_ttbarOther" << std::endl;
      datacard_ << "jes group = CMS_scaleAbsoluteMPFBias_j CMS_scaleAbsoluteScale_j CMS_scaleAbsoluteStat_j CMS_scaleFlavorQCD_j CMS_scaleFragmentation_j CMS_scalePileUpDataMC_j CMS_scalePileUpPtBB_j CMS_scalePileUpPtEC1_j CMS_scalePileUpPtEC2_j CMS_scalePileUpPtHF_j CMS_scalePileUpPtRef_j CMS_scaleRelativeBal_j CMS_scaleRelativeFSR_j CMS_scaleRelativeJEREC1_j CMS_scaleRelativeJEREC2_j CMS_scaleRelativeJERHF_j CMS_scaleRelativePtBB_j CMS_scaleRelativePtEC1_j CMS_scaleRelativePtEC2_j CMS_scaleRelativePtHF_j CMS_scaleRelativeStatEC_j CMS_scaleRelativeStatFSR_j CMS_scaleRelativeStatHF_j CMS_scaleSinglePionECAL_j CMS_scaleSinglePionHCAL_j CMS_scaleTimePtEta_j" << std::endl;
      datacard_<< "theory group = QCDscale_V QCDscale_VV QCDscale_singlet QCDscale_ttH QCDscale_ttbar bgnorm_ttbarPlus2B bgnorm_ttbarPlusB bgnorm_ttbarPlusBBbar bgnorm_ttbarPlusCCbar pdf_gg pdf_qg pdf_qqbar pdf_Higgs_ttH CMS_ttHbb_PDF CMS_ttHbb_scaleMuF CMS_ttHbb_scaleMuR CMS_ttHbb_UE_ttbarPlusBBbar CMS_ttHbb_UE_ttbarPlus2B CMS_ttHbb_UE_ttbarPlusB CMS_ttHbb_UE_ttbarPlusCCbar CMS_ttHbb_UE_ttbarOther CMS_ttHbb_ISR_ttbarPlusBBbar CMS_ttHbb_ISR_ttbarPlus2B CMS_ttHbb_ISR_ttbarPlusB CMS_ttHbb_ISR_ttbarPlusCCbar CMS_ttHbb_ISR_ttbarOther CMS_ttHbb_FSR_ttbarPlusBBbar CMS_ttHbb_FSR_ttbarPlus2B CMS_ttHbb_FSR_ttbarPlusB CMS_ttHbb_FSR_ttbarPlusCCbar CMS_ttHbb_FSR_ttbarOther CMS_ttHbb_HDAMP_ttbarPlusBBbar CMS_ttHbb_HDAMP_ttbarPlus2B CMS_ttHbb_HDAMP_ttbarPlusB CMS_ttHbb_HDAMP_ttbarPlusCCbar CMS_ttHbb_HDAMP_ttbarOther" << std::endl;
      datacard_<< "btag group = CMS_btag_lf CMS_btag_hf CMS_btag_hfstats1 CMS_btag_hfstats2 CMS_btag_cferr1 CMS_btag_cferr2 CMS_btag_lfstats1 CMS_btag_lfstats2" << std::endl;
      datacard_<< "bgnorm group = bgnorm_ttbarPlus2B bgnorm_ttbarPlusB bgnorm_ttbarPlusBBbar bgnorm_ttbarPlusCCbar" << std::endl;
      datacard_<< "pdf group = pdf_gg  pdf_qg pdf_qqbar pdf_Higgs_ttH" << std::endl;
      datacard_<< "QCDscale group = QCDscale_V QCDscale_VV QCDscale_singlet QCDscale_ttH QCDscale_ttbar" << std::endl;
      datacard_<< "misc group = CMS_ttHbb_UE_ttbarPlusBBbar CMS_ttHbb_UE_ttbarPlus2B CMS_ttHbb_UE_ttbarPlusB CMS_ttHbb_UE_ttbarPlusCCbar CMS_ttHbb_UE_ttbarOther CMS_ttHbb_ISR_ttbarPlusBBbar CMS_ttHbb_ISR_ttbarPlus2B CMS_ttHbb_ISR_ttbarPlusB CMS_ttHbb_ISR_ttbarPlusCCbar CMS_ttHbb_ISR_ttbarOther CMS_ttHbb_FSR_ttbarPlusBBbar CMS_ttHbb_FSR_ttbarPlus2B CMS_ttHbb_FSR_ttbarPlusB CMS_ttHbb_FSR_ttbarPlusCCbar CMS_ttHbb_FSR_ttbarOther CMS_ttHbb_HDAMP_ttbarPlusBBbar CMS_ttHbb_HDAMP_ttbarPlus2B CMS_ttHbb_HDAMP_ttbarPlusB CMS_ttHbb_HDAMP_ttbarPlusCCbar CMS_ttHbb_HDAMP_ttbarOther" << std::endl;
    }

    std::cout << "Closing file: " << datacardName_ << std::endl;
    datacard_.close();
  }

}


void DatacardMaker::writeHeader(const std::string& name)
{
  TObjArray* token  = TString(name).Tokenize("_");
  TString filename  = TString(name);
  TString eventCategory  = ((TObjString*)token->At(token->GetLast()))->GetString();

  if(!datacard_.is_open())
    datacard_.open(datacardName_, std::ios::out);

  datacard_ << "imax\t*\tnumber of categories" << std::endl;
  datacard_ << "jmax\t*\tnumber of samples minus one" << std::endl;
  datacard_ << "kmax\t*\tnumber of nuisance parameter" << std::endl;
  datacard_ << "----------------------------------------------------------------------------------------------------------------------" << std::endl;
  datacard_ << "\nshapes * * " << "common/ttH_hbb_13TeV_dl.root" << "\t$CHANNEL_"+observableType_+"/$PROCESS\t$CHANNEL_"+observableType_+"/$PROCESS_$SYSTEMATIC" << std::endl;
  //datacard_ << "\nshapes ttH$MASS_hbb * " << "ttH_hbb_13TeV_dl.root" << "\t$CHANNEL_"+observableType_+"/$PROCESS$MASS\t$CHANNEL_"+observableType_+"/$PROCESS$MASS_$SYSTEMATIC" << std::endl;
  datacard_ << "----------------------------------------------------------------------------------------------------------------------" << std::endl;

  datacard_.close();
}


void DatacardMaker::writeSystematicUncertainties() 
{ 
  if(!datacard_.is_open())
    datacard_.open(datacardName_, std::ios::out | std::ios::app);

  // Write out the sources of systematic/statistical uncertainty and their values to the datacard file 
  for(Systematic::Systematic systematic : v_systematic_) {

    // Constructing a neutral systematic for which two variations will be stored         
    Systematic::Systematic systematicToStore = Systematic::Systematic(systematic.type(), Systematic::undefinedVariation, systematic.variationNumber());

    std::string nameOfSystematic = systematicToStore.name().Data();
    std::string tmp0 = convertSystematicLabel_[nameOfSystematic+"_UP"];
    std::string tmp1;

    if(!tmp0.empty())
      tmp1 = tmp0.substr(0,tmp0.size()-2);
    //std::string tmp = TString(convertSystematicLabel_[nameOfSystematic+"_UP"]).ReplaceAll("Up","").Data();
    std::string tmp = tmp1;
    const char* systematType = tmp.c_str();    

    std::string shape("shape");
    std::string lnN("lnN");

    if(systematic.variation() == Systematic::down) continue;

    // Select if systematic uncertainty is of type shape and lnN
    if(std::find_if(listOfSystematicsByType_.begin(),listOfSystematicsByType_.end(),compare(nameOfSystematic,shape)) != listOfSystematicsByType_.end()) {

      if(systematic.type() == Systematic::btagDiscrBstat1     || systematic.type() == Systematic::btagDiscrBstat2 
         || systematic.type() == Systematic::btagDiscrLstat1  || systematic.type() == Systematic::btagDiscrLstat2
         || systematic.type() == Systematic::btagDiscrBpurity || systematic.type() == Systematic::btagDiscrLpurity
         || systematic.type() == Systematic::btagDiscrCerr1   || systematic.type() == Systematic::btagDiscrCerr2
         || systematic.type() == Systematic::jes || systematic.type() == Systematic::jer 
         || systematic.type() == Systematic::jesAbsoluteStat
         || systematic.type() == Systematic::jesAbsoluteScale
         || systematic.type() == Systematic::jesAbsoluteFlavMap
         || systematic.type() == Systematic::jesAbsoluteMPFBias
         || systematic.type() == Systematic::jesFragmentation
         || systematic.type() == Systematic::jesSinglePionECAL
         || systematic.type() == Systematic::jesSinglePionHCAL
         || systematic.type() == Systematic::jesFlavorQCD
         || systematic.type() == Systematic::jesTimePtEta
         //|| systematic.type() == Systematic::jesTimePt
         || systematic.type() == Systematic::jesRelativeJEREC1
         || systematic.type() == Systematic::jesRelativeJEREC2
         || systematic.type() == Systematic::jesRelativeJERHF
         || systematic.type() == Systematic::jesRelativePtBB
         || systematic.type() == Systematic::jesRelativePtEC1
         || systematic.type() == Systematic::jesRelativePtEC2
         || systematic.type() == Systematic::jesRelativePtHF
         || systematic.type() == Systematic::jesRelativeFSR
         || systematic.type() == Systematic::jesRelativeBal
         || systematic.type() == Systematic::jesRelativeStatFSR
         || systematic.type() == Systematic::jesRelativeStatEC
         || systematic.type() == Systematic::jesRelativeStatHF
         || systematic.type() == Systematic::jesPileUpDataMC
         || systematic.type() == Systematic::jesPileUpPtRef
         || systematic.type() == Systematic::jesPileUpPtBB
         || systematic.type() == Systematic::jesPileUpPtEC1
         || systematic.type() == Systematic::jesPileUpPtEC2
         || systematic.type() == Systematic::jesPileUpPtHF
         || systematic.type() == Systematic::jesPileUpMuZero
         || systematic.type() == Systematic::jesPileUpEnvelope
         || systematic.type() == Systematic::jesSubTotalPileUp
         || systematic.type() == Systematic::jesSubTotalRelative
         || systematic.type() == Systematic::jesSubTotalPt
         || systematic.type() == Systematic::jesSubTotalScale
         || systematic.type() == Systematic::jesSubTotalMC
         || systematic.type() == Systematic::jesSubTotalAbsolute
         || systematic.type() == Systematic::jesTotalNoFlavor
         || systematic.type() == Systematic::jesFlavorZJet
         || systematic.type() == Systematic::jesFlavorPhotonJet
         || systematic.type() == Systematic::jesFlavorPureGluon
         || systematic.type() == Systematic::jesFlavorPureQuark
         || systematic.type() == Systematic::jesFlavorPureCharm
         || systematic.type() == Systematic::jesFlavorPureBottom
         || systematic.type() == Systematic::jesCorrelationGroupMPFInSitu
         || systematic.type() == Systematic::jesCorrelationGroupIntercalibration
         || systematic.type() == Systematic::jesCorrelationGroupbJES
         || systematic.type() == Systematic::jesCorrelationGroupFlavor
         || systematic.type() == Systematic::jesCorrelationGroupUncorrelated
         || systematic.type() == Systematic::pdf
         || systematic.type() == Systematic::scale_ttb    || systematic.type() == Systematic::scale_ttbb      || systematic.type() == Systematic::scale_tt2b
         || systematic.type() == Systematic::scale_ttcc   || systematic.type() == Systematic::scale_ttother   || systematic.type() == Systematic::scale
         || systematic.type() == Systematic::meScale_ttb  || systematic.type() == Systematic::meScale_ttbb    || systematic.type() == Systematic::meScale_tt2b
         || systematic.type() == Systematic::meScale_ttcc || systematic.type() == Systematic::meScale_ttother || systematic.type() == Systematic::meScale 
         || systematic.type() == Systematic::psScale_ttb  || systematic.type() == Systematic::psScale_ttbb    || systematic.type() == Systematic::psScale_tt2b
         || systematic.type() == Systematic::psScale_ttcc || systematic.type() == Systematic::psScale_ttother || systematic.type() == Systematic::psScale
         || systematic.type() == Systematic::pu || systematic.type() == Systematic::trig || systematic.type() == Systematic::lept
         || systematic.type() == Systematic::meFacScale || systematic.type() == Systematic::meRenScale 
         || systematic.type() == Systematic::psFSRScale
         || systematic.type() == Systematic::psFSRScale_ttbb || systematic.type() == Systematic::psFSRScale_ttb || systematic.type() == Systematic::psFSRScale_tt2b || systematic.type() == Systematic::psFSRScale_ttcc || systematic.type() == Systematic::psFSRScale_ttother
         || systematic.type() == Systematic::psISRScale || systematic.type() == Systematic::ueTune || systematic.type() == Systematic::match
         ) {
        datacard_ << TString::Format("%-32s shape\t\t", systematType);

        for (auto p : processNames_) {

          if((p == "data") || (p == "allmc")) continue;

          // Uncertainty value for all processess or selected processes
          if(std::find_if(valueOfSystematicBasedOnProcess_[systematicToStore.name().Data()].begin(), valueOfSystematicBasedOnProcess_[systematicToStore.name().Data()].end(),comp("all")) != valueOfSystematicBasedOnProcess_[systematicToStore.name().Data()].end()) {
            std::string value = valueOfSystematicBasedOnProcess_[systematicToStore.name().Data()].at(0).second;

            datacard_ << TString::Format("%-8.11s\t", value.c_str());            
          }
          else if(std::find_if(valueOfSystematicBasedOnProcess_[systematicToStore.name().Data()].begin(), valueOfSystematicBasedOnProcess_[systematicToStore.name().Data()].end(), comp(p)) != valueOfSystematicBasedOnProcess_[systematicToStore.name().Data()].end()) {

            auto iter = std::find_if(valueOfSystematicBasedOnProcess_[systematicToStore.name().Data()].begin(), valueOfSystematicBasedOnProcess_[systematicToStore.name().Data()].end(), comp(p));
            auto index = std::distance(valueOfSystematicBasedOnProcess_[systematicToStore.name().Data()].begin(), iter);
            std::string value = valueOfSystematicBasedOnProcess_[systematicToStore.name().Data()].at(index).second;

            datacard_ << TString::Format("%-8.11s\t", value.c_str()); 
          }
          else {                                                                                                                          
            datacard_ << TString::Format("%-8.11s\t", "-");
          }
        }
        datacard_ << std::endl;
      }
    }
    else if(std::find_if(listOfSystematicsByType_.begin(),listOfSystematicsByType_.end(),compare(nameOfSystematic,lnN)) != listOfSystematicsByType_.end()) {

      if(systematic.type() == Systematic::xsec_ttb || systematic.type() == Systematic::xsec_ttbb || systematic.type() == Systematic::xsec_tt2b
         || systematic.type() == Systematic::xsec_ttcc || systematic.type() == Systematic::xsec_ttother
         || systematic.type() == Systematic::psFSRScale 
         || systematic.type() == Systematic::psFSRScale_ttbb
         || systematic.type() == Systematic::psFSRScale_tt2b
         || systematic.type() == Systematic::psFSRScale_ttb
         || systematic.type() == Systematic::psFSRScale_ttcc
         || systematic.type() == Systematic::psFSRScale_ttother
         || systematic.type() == Systematic::psISRScale 
         || systematic.type() == Systematic::psISRScale_ttbb
         || systematic.type() == Systematic::psISRScale_tt2b
         || systematic.type() == Systematic::psISRScale_ttb
         || systematic.type() == Systematic::psISRScale_ttcc
         || systematic.type() == Systematic::psISRScale_ttother
         || systematic.type() == Systematic::ueTune
         || systematic.type() == Systematic::ueTune_ttbb
         || systematic.type() == Systematic::ueTune_tt2b
         || systematic.type() == Systematic::ueTune_ttb
         || systematic.type() == Systematic::ueTune_ttcc
         || systematic.type() == Systematic::ueTune_ttother
         || systematic.type() == Systematic::match
         || systematic.type() == Systematic::match_ttbb
         || systematic.type() == Systematic::match_tt2b
         || systematic.type() == Systematic::match_ttb
         || systematic.type() == Systematic::match_ttcc
         || systematic.type() == Systematic::match_ttother
         ) {
        datacard_ << TString::Format("%-32s lnN\t\t", systematType); 
        
        for (auto p : processNames_) {

          if((p == "data") || (p == "allmc")) continue;

          // Cross section systematic uncertainties apply to only selective processes
          auto iter = std::find_if(valueOfSystematicBasedOnProcess_[systematicToStore.name().Data()].begin(), valueOfSystematicBasedOnProcess_[systematicToStore.name().Data()].end(), comp(p));

          if(iter != valueOfSystematicBasedOnProcess_[systematicToStore.name().Data()].end()) {

            auto index = std::distance(valueOfSystematicBasedOnProcess_[systematicToStore.name().Data()].begin(), iter);
            std::string value = valueOfSystematicBasedOnProcess_[systematicToStore.name().Data()].at(index).second;

            datacard_ << TString::Format("%-8.11s\t", value.c_str());
          }
          else datacard_ << TString::Format("%-8.11s\t", "-");
        }
        datacard_ << std::endl;
      } 
      else if (systematic.type() == Systematic::lumi || systematic.type() == Systematic::trig || systematic.type() ==  Systematic::xsec_ttH
               || systematic.type() == Systematic::normPdfGg  || systematic.type() == Systematic::normPdfGq || systematic.type() == Systematic::normPdfQq
               || systematic.type() == Systematic::normPdfTth || systematic.type() == Systematic::xsec_tt   || systematic.type() == Systematic::xsec_t
               || systematic.type() == Systematic::xsec_v || systematic.type() == Systematic::xsec_vv
               || systematic.type() == Systematic::xsec_ttZ || systematic.type() == Systematic::xsec_ttW || systematic.type() == Systematic::xsec_ttG
               //|| systematic.type() == Systematic::lept
               ) {
        datacard_ << TString::Format("%-32s lnN\t\t", systematType);

        for (auto p : processNames_) {

          if((p == "data") || (p == "allmc")) continue;

          // Uncertainty value for all processess or selected processes
          if(std::find_if(valueOfSystematicBasedOnProcess_[systematicToStore.name().Data()].begin(), valueOfSystematicBasedOnProcess_[systematicToStore.name().Data()].end(),comp("all")) != valueOfSystematicBasedOnProcess_[systematicToStore.name().Data()].end()) {

            auto iter = std::find_if(valueOfSystematicBasedOnProcess_[systematicToStore.name().Data()].begin(), valueOfSystematicBasedOnProcess_[systematicToStore.name().Data()].end(), comp("all"));
            auto index = std::distance(valueOfSystematicBasedOnProcess_[systematicToStore.name().Data()].begin(), iter);
            std::string value = valueOfSystematicBasedOnProcess_[systematicToStore.name().Data()].at(index).second;

            datacard_ << TString::Format("%-8.11s\t", value.c_str());
          }
          else if(std::find_if(valueOfSystematicBasedOnProcess_[systematicToStore.name().Data()].begin(), valueOfSystematicBasedOnProcess_[systematicToStore.name().Data()].end(), comp(p)) != valueOfSystematicBasedOnProcess_[systematicToStore.name().Data()].end()) {

            auto iter = std::find_if(valueOfSystematicBasedOnProcess_[systematicToStore.name().Data()].begin(), valueOfSystematicBasedOnProcess_[systematicToStore.name().Data()].end(), comp(p));
            auto index = std::distance(valueOfSystematicBasedOnProcess_[systematicToStore.name().Data()].begin(), iter);
            std::string value = valueOfSystematicBasedOnProcess_[systematicToStore.name().Data()].at(index).second;

            datacard_ << TString::Format("%-8.11s\t", value.c_str());               
          }
          else datacard_ << TString::Format("%-8.11s\t", "-");
        }
        datacard_ << std::endl;
      }
    }    
  }
  datacard_.close();
}


void DatacardMaker::writeStatisticalUncertainties(const std::string& name)
{
  // Force all histograms to use option Sumw2(), to switch histogram errors
  TH1::SetDefaultSumw2();

  TObjArray* token  = TString(name).Tokenize("_");
  TString filename  = TString(name);
  TString eventCategory  = ((TObjString*)token->At(token->GetLast()))->GetString();

  if(!datacard_.is_open())
    datacard_.open(datacardName_, std::ios::out | std::ios::app);

  std::map<TString,TH1*> mapOfNominalHisto;
  std::map<TString, int> shiftUpDown;

  shiftUpDown["Down"] = -1.0;
  shiftUpDown["Up"] = 1.0;

  for(auto channelCollection = inputFileLists_->begin(); channelCollection != inputFileLists_->end(); ++channelCollection) {

    double data_bin_error(-999.);

    double sig_bin_content(-999.);
    //double sig_bin_error(-999.);
    double bkg_bin_content(-999.);
    double bkg_bin_error(-999.);
    double sample_bin_content(-999.);
    double sample_bin_error(-999.);
    double other_frac(-999.);

    TH1* sig_hist(NULL);
    TH1* bkg_hist(NULL);
    TH1* data_hist(NULL);

    std::string obsProcess;

    for(TString processName : processNames_) {

      if (processName.Contains("data") || processName.Contains("allmc")) {

        obsProcess = processName;
        continue;
      }

      for(auto systematicCollection : (*channelCollection).second) {

        TString nominalFile;

        for(auto nameOfFile : systematicCollection.second) {
          nominalFile = TString(nameOfFile.second);
        }

        Systematic::Systematic systematic = systematicCollection.first;

        if(systematic.type() == Systematic::nominal) {

          data_hist =  addOrCreateHisto(data_hist, static_cast<TH1*>((TH1F*)fileReader_->GetClone<TH1>(nominalFile, TString(name)+"_"+TString(obsProcess), true, false)->Clone()));

          for (auto pro : processNames_) {

            if(TString(pro).Contains(signalModel_))
              sig_hist = addOrCreateHisto(sig_hist, static_cast<TH1*>(fileReader_->GetClone<TH1>(nominalFile, name+"_"+TString(pro), true, false)));
            
            if((pro == "data") || (pro == "allmc") || (pro == "signamlc") || TString(pro).Contains(signalModel_)) continue;

            bkg_hist = addOrCreateHisto(bkg_hist, (TH1*)fileReader_->GetClone<TH1>(nominalFile, TString(name)+"_"+TString(pro), true, false));
          }

          TH1* sample_hist(NULL);
          sample_hist = addOrCreateHisto(sample_hist, static_cast<TH1*>(fileReader_->GetClone<TH1>(nominalFile, name+"_"+processName, true, false)));

          int numBins = sample_hist->GetNbinsX();

          for(auto shift : shiftUpDown){

            // Set protection against down shift to zero for bin content
            if(shift.second < 0 && sample_hist->Integral() == 0.) 
              addOrCreateHisto(static_cast<TH1*>(mapOfNominalHisto[processName]), sample_hist);
            else
              mapOfNominalHisto[processName] = (TH1D*)sample_hist->Clone();

            for (int iBin = 0; iBin < numBins; ++iBin) {

              data_bin_error = data_hist->GetBinError(iBin+1);

              sig_bin_content = sig_hist->GetBinContent(iBin+1);
              //sig_bin_error   = sig_hist->GetBinError(iBin+1);

              bkg_bin_content = bkg_hist->GetBinContent(iBin+1);
              bkg_bin_error   = bkg_hist->GetBinError(iBin+1);

              sample_bin_content = mapOfNominalHisto[processName]->GetBinContent(iBin+1);
              sample_bin_error   = mapOfNominalHisto[processName]->GetBinError(iBin+1);

              other_frac = sqrt(bkg_bin_error*bkg_bin_error - sample_bin_error*sample_bin_error);             

              if(pruneBinByBin_) {

                // MC stat. prunining method
                /* NOTE: Information obtained from the following resources
                   https://twiki.cern.ch/twiki/bin/view/CMS/TTbarHbbRun2ReferenceAnalysisLimits
                   https://indico.cern.ch/event/373752/session/6/contribution/15/attachments/744533/1021297/Hbb_bbbStat_24022015.pdf
                   https://github.com/cms-ttH/ttH-Limits/blob/13TeV/python/datacard_ttbb13TeV.py
                 */
                if(bkg_bin_error < data_bin_error/5.
                   || sig_bin_content/bkg_bin_content < 0.01
                   || sample_bin_content < 0.01
                   || other_frac/bkg_bin_error > 0.95) {
                  continue;
                }
              }

              // Extract both histogram statistical error and bin content
              const float error   = mapOfNominalHisto[processName]->GetBinError(iBin+1);
              const float content = mapOfNominalHisto[processName]->GetBinContent(iBin+1);

              // Create shifted histogram
              TH1* histo = (TH1*)mapOfNominalHisto[processName]->Clone();
              const float shifted_content = std::max(content + shift.second*error, float(1e-8));
              histo->SetBinContent(iBin+1, shifted_content);

              const char* histo_Bin = convertSampleNames_[processName.Data()].c_str();
              TString histo_name  = TString::Format("CMS_ttH_%s_%s_13TeV_%sbin%d",histo_Bin, mapOfCategories_[eventCategory.Data()].c_str(), observableType_.c_str(), iBin+1);

              outputFile_ = new TFile(TString(outputDirDatacard_+outputFileName_),"UPDATE");
              outputFile_->cd(TString(mapOfCategories_[eventCategory.Data()]+"_"+observableType_));

              // Only need to include statistical uncertainties for both signal and background once in the datacard
              if(shift.first.Contains("Down")) {

                datacard_  << TString::Format("%-32s shape\t", histo_name.Data());

                for(auto process : processNames_) {

                  if(process == "data" || process == "allmc") continue;

                  // Default value of 1.0 assigned to MC stats shape systematics
                  if(convertSampleNames_[process] == histo_Bin)
                    datacard_  << TString::Format("%-8.11s\t", "1.000000");
                  else
                    datacard_  << TString::Format("%-8.11s\t", "-");
                }

                datacard_  << std::endl;
              }
              histo->Write(TString(histo_Bin)+"_"+histo_name+(shift.first));

              outputFile_->Write("",TObject::kOverwrite);
              outputFile_->Close();
            }
          }
        }
      }
    }
  }
  // End statistical uncertainty estimate
  datacard_.close();
}


void DatacardMaker::extractYields(const std::string& name)
{
  // Parse string into tokens
  TObjArray* token  = TString(name).Tokenize("_");
  TString filename  = TString(name);
  TString eventCategory = ((TObjString*)token->At(token->GetLast()))->GetString();

  if(!datacard_.is_open())
    datacard_.open(datacardName_, std::ios::out | std::ios::app);

  std::map<TString, TString>::iterator systematic; 
  std::map<TString,TH1D*> mapOfHistograms;

  TString observationRate;
  TString signalRate;

  std::map<TString, TString> processRates;

  TString process;

  // Begin iterating through the Plots directory for File System hierarchy
  if(filename.Contains(eventCategory)) {  

    // Loop over all systematics and channels
    for(auto fileCollection = inputFileLists_->begin(); fileCollection != inputFileLists_->end(); ++fileCollection) {

      for(auto uncertainty : (*fileCollection).second) {    

        Systematic::Systematic systematic = uncertainty.first;

        // Create input file handler
        TFile* sampleFile = new TFile((uncertainty.second.find(name+"_source.root")->second).c_str());

        bool debug_ = false;
        if(debug_)
          std::cout << sampleFile->GetName() << std::endl;

        // Check that file exists and is not corrupt
        if(sampleFile->IsZombie()) {
          std::cout<<"\n\tWe didn't find the "+TString(filename+"_source.root")+" input!!\n";
        }
        else {
        // Update (i.e. append) to output file
          if(!outputFile_) {
            outputFile_ = new TFile(TString(outputDirDatacard_+outputFileName_),"RECREATE");
          }
          else {
            outputFile_ = new TFile(TString(outputDirDatacard_+outputFileName_),"UPDATE");
          }
          
          // Create list from list of histogram from root file
          TList* list = sampleFile->GetListOfKeys();

          bool debug_ = false;
          if(debug_)list->Print();
          
          if(!list) {
            std::cout << TString::Format("<Error> No keys found in file\n");
            exit(1);
          }
          
          // Create iterator from list
          TIter next((TList*)list);
          TKey* key;

          // Iterate and extract histograms from file
          while((key=(TKey*)next())){
            TObjArray *subString;
            TString processType;

            if(filename == key->GetName()) continue; // ignore TCanvas object normally found in histogram root files
            subString = TPRegexp(filename+"_(\\w+)").MatchS(TString(key->GetName()));
            processType  = ((TObjString *)subString->At(1))->GetString();

            if(std::find(processNames_.begin(), processNames_.end(), processType)!= processNames_.end()){ 
              mapOfHistograms[processType] = (TH1D*)sampleFile->Get(key->GetName());
            }
          }

          // Check histogram exist
          for (std::map<TString,TH1D*>::iterator mvaHisto = mapOfHistograms.begin(); mvaHisto != mapOfHistograms.end(); ++mvaHisto) {
            if(!mapOfHistograms[mvaHisto->first])
              std::cout << "\n\tWe didn't find the "+mvaHisto->first+" histogram!!\n";
          }

          // Write histograms to file
          for(std::map<TString,TH1D*>::iterator mvaHisto = mapOfHistograms.begin(); mvaHisto != mapOfHistograms.end(); ++mvaHisto) {

            if(systematic.type() == Systematic::nominal) {

              if((mvaHisto->first == "data") || (mvaHisto->first == "allmc")) {
                observationRate = TString::Format("%f\t", mapOfHistograms[mvaHisto->first]->Integral());
              }
              else {
                processRates[mvaHisto->first] += TString::Format("%f\t", mapOfHistograms[mvaHisto->first]->Integral());
              }
            }
          }

          if(systematic.type() == Systematic::nominal) {

            datacard_ << TString::Format("\nbin\t\t%-s",mapOfCategories_[eventCategory.Data()].c_str()) << std::endl;
            datacard_ << TString::Format("observation\t%-s",observationRate.Data()) << std::endl;
            datacard_ << "----------------------------------------------------------------------------------------------------------------------" << std::endl;
            datacard_ << "\nbin\t";
            for(std::size_t i = 0; i < processNames_.size(); ++i) {
              if((processNames_[i] != "data") && (processNames_[i] != "allmc")) {
                datacard_  << TString::Format("%-8s\t", mapOfCategories_[eventCategory.Data()].c_str());
              }
            }
            datacard_ << "\nprocess\t";

            for(auto process : processNames_){
              if((process == "data") || (process == "allmc")) continue;
              datacard_  << TString::Format("%-8s\t", convertSampleNames_[process].c_str());
            }
            datacard_ << "\nprocess\t"; 

            int index = std::count_if(processNames_.begin(), processNames_.end(), [this](std::string word) {return (word.find(signalModel_) != std::string::npos);});
            for(std::size_t i = 1; i < processNames_.size(); ++i) {
              datacard_  << TString::Format("%-8d\t", (Int_t)(i-index));
            }

            datacard_ << "\nrate\t" ;

            for(auto process : processNames_) {
              if(!list->Contains(TString(process)))
                datacard_  << TString::Format("%s", processRates[process].Data());
              else
                datacard_  << TString::Format("%s", "-999.0");
            }

            datacard_ << "\n---------------------------------------------------------------------------------------------------------------------" << std::endl;
          } 

          outputFile_->Write("",TObject::kOverwrite);
          outputFile_->Close();
        } // End section where histograms are written

        sampleFile->Close();
        processRates.clear();
      }
    }
  }
  datacard_.close();
  return;
}


void DatacardMaker::writeYields(const std::string& name)
{
  // Parse string into tokens
  TObjArray* token  = TString(name).Tokenize("_");
  TString filename  = TString(name);
  TString eventCategory = ((TObjString*)token->At(token->GetLast()))->GetString();
  
  std::map<TString, TString>::iterator systematic; 
  std::map<TString,TH1D*> mapOfHistograms;

  TString observationRate;
  TString signalRate;

  TString process;

  // Begin iterating through the Plots directory for File System hierarchy
  if(filename.Contains(eventCategory)) {  

    // Loop over all systematics and channels
    for(auto fileCollection = inputFileLists_->begin(); fileCollection != inputFileLists_->end(); ++fileCollection) {

      for(auto uncertainty : (*fileCollection).second) {  
        
        Systematic::Systematic systematic = uncertainty.first;
        
        // Create input file handler
        TFile* sampleFile = new TFile((uncertainty.second.find(name+"_source.root")->second).c_str());

        bool debug_ = false;
        if(debug_)
          std::cout << sampleFile->GetName() << std::endl;

        // Check that file exists and is not corrupt
        if(sampleFile->IsZombie()) {
          std::cout<<"\n\tWe didn't find the "+TString(filename+"_source.root")+" input!!\n";
        }
        else {
          // Update (i.e. append) to output file
          if(!outputFile_) {
            outputFile_ = new TFile(TString(outputDirDatacard_+outputFileName_),"RECREATE");
          }
          else {
            outputFile_ = new TFile(TString(outputDirDatacard_+outputFileName_),"UPDATE");
          }

          // Create list from list of histogram from root file
          TList* list = sampleFile->GetListOfKeys();

          bool debug_ = false;
          if(debug_)list->Print();
          
          if(!list) {
            std::cout << TString::Format("<Error> No keys found in file\n");
            exit(1);
          }
          
          // Create iterator from list
          TIter next((TList*)list);
          TKey* key;
          
          // Iterate and extract histograms from file
          while((key=(TKey*)next())){

            TObjArray *subString;
            TString processType;
            
            if(filename == key->GetName()) continue; // ignore TCanvas object normally found in histogram root files
            subString = TPRegexp(filename+"_(\\w+)").MatchS(TString(key->GetName()));
            processType  = ((TObjString *)subString->At(1))->GetString();

            if(std::find(processNames_.begin(), processNames_.end(), processType)!= processNames_.end()) { 
              mapOfHistograms[processType] = (TH1D*)sampleFile->Get(key->GetName());
            }
          }

          // Check histogram exist
          for (std::map<TString,TH1D*>::iterator mvaHisto = mapOfHistograms.begin(); mvaHisto != mapOfHistograms.end(); ++mvaHisto) {
            if(!mapOfHistograms[mvaHisto->first])
              std::cout << "\n\tWe didn't find the "+mvaHisto->first+" histogram!!\n";
          }
          
          // Write histograms to correct event category directory
          if(outputFile_->GetDirectory(TString(mapOfCategories_[eventCategory.Data()]+"_"+observableType_))) {
              outputFile_->cd(TString(mapOfCategories_[eventCategory.Data()]+"_"+observableType_));
          }
          else {
            outputFile_->mkdir(TString(mapOfCategories_[eventCategory.Data()]+"_"+observableType_),TString(mapOfCategories_[eventCategory.Data()]+"_"+observableType_));
            outputFile_->cd(TString(mapOfCategories_[eventCategory.Data()]+"_"+observableType_));
          }

          // Write histograms to file
          for (std::map<TString,TH1D*>::iterator mvaHisto = mapOfHistograms.begin(); mvaHisto != mapOfHistograms.end(); ++mvaHisto) {
  
            if((mvaHisto->first == "data" || mvaHisto->first == "allmc") && systematic.type() == Systematic::nominal) {
  
              process =  convertSampleNames_[mvaHisto->first.Data()];
              mapOfHistograms[mvaHisto->first]->Write(process,TObject::kOverwrite);
            }
            else if (mvaHisto->first != "data" && mvaHisto->first != "allmc") {
  
              if(systematic.type() == Systematic::nominal) {
                process =  mvaHisto->first;
                mapOfHistograms[mvaHisto->first]->Write(TString(convertSampleNames_[process.Data()]),TObject::kOverwrite);
              }
              else {
                process =  convertSampleNames_[(mvaHisto->first).Data()]+"_"+convertSystematicLabel_[systematic.name().Data()]; 
                mapOfHistograms[mvaHisto->first]->Write(process,TObject::kOverwrite);
              }
            }
          }
          outputFile_->Write("",TObject::kOverwrite);
          outputFile_->Close();
        } // End section where histograms are written
        sampleFile->Close();
      }
    }
  }
  return;
}


TH1* DatacardMaker::addOrCreateHisto(TH1* base, const TH1* const add_histo) const
{
  TH1* tmp;
  if(!base) {
    tmp = static_cast<TH1*>(add_histo->Clone());
    for(int iBin = 0; iBin < tmp->GetNbinsX(); ++iBin) {
      const float content = std::max(tmp->GetBinContent(iBin+1), 1e-8);
      
      tmp->SetBinContent(iBin+1, content);
    }
    base = (TH1*)tmp->Clone();
  }
  else {
    TH1* tmp = static_cast<TH1*>(add_histo->Clone());
    for(int iBin = 0; iBin < tmp->GetNbinsX(); ++iBin) {
      const float content = std::max(tmp->GetBinContent(iBin+1), 1e-8);
      
      tmp->SetBinContent(iBin+1, content);
    }
    base->Add(tmp);
  }

  return base;
}

void DatacardMaker::setIncludeSystmeticUncertainties(bool useSys) {
  
  addSystematicUncertainty_ = useSys;
}

void DatacardMaker::setIncludeStatisticalUncertainties(bool useStat) {

  addStatisticalUncertainty_ = useStat;
}

