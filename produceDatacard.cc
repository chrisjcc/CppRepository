#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <set>
#include <string>
#include <algorithm>

#include <TString.h>

#include "AnalysisConfig.h"
#include "Samples.h"
#include "GlobalScaleFactors.h"
#include "EventYields.h"
#include "plotterHelpers.h"
#include "DatacardMaker.h"
#include "PlotterSystematic.h"
#include "HistoListReader.h"
#include "higgsUtils.h"
#include "../../common/include/utils.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/CommandLineParameters.h"




/// All systematics allowed for plotting
namespace Systematic{
  const std::vector<Type> allowedSystematics = {
    nominal, all, allAvailable,
    pu, lept, trig,
    
    jer, jes,

    jesAbsoluteStat, //0
    jesAbsoluteScale,//1
    jesAbsoluteFlavMap,//2
    jesAbsoluteMPFBias,//3
    jesFragmentation,   // "HighPtExtra,    //4 -----
    jesSinglePionECAL, //5
    jesSinglePionHCAL, //6
    jesFlavorQCD,        //7
    jesTimePtEta,               // "Time" 8 -----
    //jesTimePt,
    jesRelativeJEREC1,    //9
    jesRelativeJEREC2,    //10
    jesRelativeJERHF,   //11
    jesRelativePtBB,       //12
    jesRelativePtEC1,    //13
    jesRelativePtEC2,    //14
    jesRelativePtHF,    //15
    jesRelativeFSR,       //16
    jesRelativeStatFSR,   //new
    jesRelativeStatEC,    //2,    //17 -----
    jesRelativeStatHF,   //18
    jesRelativeBal,    //
    jesPileUpDataMC,       //19
    jesPileUpPtRef,
    jesPileUpPtBB,        //20
    jesPileUpPtEC1,        //21-----PileUpPtRef
    jesPileUpPtEC2, 
    jesPileUpPtHF,       //22
    //"PileUpBias,       //23-----
    jesPileUpMuZero,
    jesPileUpEnvelope,
    jesSubTotalPileUp,   //24
    jesSubTotalRelative,   //25
    jesSubTotalPt,       //26
    jesSubTotalScale,       //26
    jesSubTotalMC,       //27
    jesSubTotalAbsolute,       //26
    jesTotalNoFlavor,   //29
    jesFlavorZJet,       //30
    jesFlavorPhotonJet,   //31
    jesFlavorPureGluon,   //32
    jesFlavorPureQuark,   //33
    jesFlavorPureCharm,   //34
    jesFlavorPureBottom,//35
    jesCorrelationGroupMPFInSitu,//36
    jesCorrelationGroupIntercalibration,//37
    jesCorrelationGroupbJES,//38
    jesCorrelationGroupFlavor,//39
    jesCorrelationGroupUncorrelated,//40 

    btag, 
    btagPt, btagEta,
    btagLjet, 
    btagLjetPt, btagLjetEta,
    btagDiscrBstat1, btagDiscrBstat2,
    btagDiscrLstat1, btagDiscrLstat2,
    btagDiscrBpurity, btagDiscrLpurity,
    btagDiscrCerr1, btagDiscrCerr2,
    kin,
    
    lumi,
    xsec_ttbb, xsec_ttb, xsec_tt2b, xsec_ttcc, xsec_ttother, 
    xsec_ttZ, xsec_ttW, xsec_ttG,
    xsec_ttH, xsec_tt, xsec_t,
    xsec_v, xsec_vv,
    topPt,
    mass, 
    match, 
    match_ttbb, match_tt2b,  match_ttb,  match_ttcc,  match_ttother,
    powheg, powhegHerwig, mcatnlo, perugia11, perugia11NoCR,
    normPdfGg, normPdfGq, normPdfQq, normPdfTth,
    scale,     scale_ttbb,     scale_ttb, scale_tt2b,   scale_ttcc,   scale_ttother,
    meScale, meScale_ttbb, meScale_ttb, meScale_tt2b, meScale_ttcc, meScale_ttother,
    meFacScale, meRenScale,
    psScale, psScale_ttbb, psScale_ttb, psScale_tt2b, psScale_ttcc, psScale_ttother,
    psISRScale,
    psISRScale_ttbb, psISRScale_tt2b, psISRScale_ttb, psISRScale_ttcc, psISRScale_ttother,
    psFSRScale,
    psFSRScale_ttbb, psFSRScale_tt2b, psFSRScale_ttb, psFSRScale_ttcc, psFSRScale_ttother,
    ueTune,
    ueTune_ttbb, ueTune_tt2b, ueTune_ttb, ueTune_ttcc, ueTune_ttother,
    pdf,
  };
}



//int produceDatacard(int argc, char** argv){
int main(int argc, char** argv){
  
  // Get and check configuration parameters
  CLParameter<std::string> opt_config("t", "Name of histogram config file in data-directory, datacards_<tag>", false, 1, 100);
  CLParameter<std::string> opt_filelist("l", "Indicate which tag to use with FileLists_plot directory version (e.g. FileList_plot_<tag>)", false, 1, 1);
  CLParameter<std::string> opt_addStatUncertainty("stat", "Include statistical uncertianties in the datacards, default set to true", false, 1, 1);
  CLParameter<std::string> opt_addSysUncertainty("sys", "Include systematic  uncertianties in the datacards, default set to true", false, 1, 1);

  CLParameter<std::string> opt_plot("p", "Name (pattern) of plot; multiple patterns possible; use '+Name' to match name exactly", false, 1, 100);
  CLParameter<std::string> opt_channel("c", "Specify channel(s), valid: emu, ee, mumu, combined. Default: all channels", false, 1, 4,
				       common::makeStringCheck(Channel::convert(Channel::allowedChannelsPlotting)));
  CLParameter<std::string> opt_systematic("s", "Systematic variation - default is Nominal, use 'all' for all, use 'allAvailable' for all available", false, 1, 1000,
					  common::makeStringCheckBegin(Systematic::convertType(Systematic::allowedSystematics)));
  CLAnalyser::interpretGlobal(argc, argv);

  std::cout<<"\n"<<"--- Beginning setting up command line steering parameters\n";

  // Set FileLists_plot folder 
  std::string fileLists("FileLists_plot_systematic");
  if(opt_filelist.isSet() && (opt_filelist.getArguments())[0] != ""){
    fileLists += "_" + (opt_filelist.getArguments())[0];
  }
  else {
    fileLists += "_Plots_selectionRoot_mvaEventA";
  }
  //config file 
  std::string configname("HistoList");
  if(opt_config.isSet() && (opt_config.getArguments())[0] != ""){
    configname += "_"+ (opt_config.getArguments())[0];
  }
  else{
    configname += "_datacards";
  }
  
  // Set up plots
  std::vector<std::string> v_plot {""};
  
  if(opt_plot.isSet()){
    v_plot = opt_plot.getArguments();
    std::cout<< "Processing only histograms containing in name: ";
    
    for(auto plot : v_plot)std::cout<< plot << " ";
    std::cout << "\n\n";
  }

  // Set up channels
  std::vector<Channel::Channel> v_channel(Channel::allowedChannelsPlotting);
  
  if(opt_channel.isSet()) v_channel = Channel::convert(opt_channel.getArguments());
  std::cout << "Processing channels: ";
  
  for(auto channel : v_channel) std::cout << Channel::convert(channel) << " ";
  std::cout << "\n\n";
  
  // Set up systematics
  std::vector<Systematic::Systematic> v_systematic = Systematic::allowedSystematicsAnalysis(Systematic::allowedSystematics);
  if(opt_systematic.isSet()){
    
    if(opt_systematic[0] == Systematic::convertType(Systematic::all))
      ; // do nothing
    else if(opt_systematic[0] == Systematic::convertType(Systematic::allAvailable)) {
      
      // Adding systematics that do not require specific root files
      for(Systematic::Type type : Systematic::fileIndependentTypes) {
	if(std::find(Systematic::upDownTypes.begin(), Systematic::upDownTypes.end(), type) != Systematic::upDownTypes.end()) {
	  v_systematic.push_back(Systematic::Systematic(type, Systematic::up));
	  v_systematic.push_back(Systematic::Systematic(type, Systematic::down));
	} else if(std::find(Systematic::centralTypes.begin(), Systematic::centralTypes.end(), type) != Systematic::centralTypes.end()) {
	  v_systematic.push_back(Systematic::Systematic(type, Systematic::central));
	}
      }
      v_systematic = common::findSystematicsFromFilelists("FileLists_plot_systematic", v_channel, v_systematic);
    }
    else
      v_systematic = Systematic::setSystematics(opt_systematic.getArguments());
  }
  else{
    v_systematic.clear();
    v_systematic.push_back(Systematic::nominalSystematic());
  }
  std::cout << "Processing systematics: "; 
  
  for(auto systematic : v_systematic) std::cout << systematic.name() << " ";
  std::cout << "\n\n";


  // Start producing datacard
  std::cout << "--- Beginning with the datacard production\n\n";

  // Read analysis config from text file
  const AnalysisConfig analysisConfig;

  DatacardMaker datacard(analysisConfig, v_plot, v_channel, v_systematic, fileLists, configname);

  if(opt_addSysUncertainty.isSet()){
    bool param = (opt_addSysUncertainty.getArguments())[0] == "true" ? true : false;
    datacard.setIncludeSystmeticUncertainties(param);
  }
  if(opt_addStatUncertainty.isSet()){
    bool param = (opt_addStatUncertainty.getArguments())[0] == "true" ? true : false;
    datacard.setIncludeStatisticalUncertainties(param);
  }
  datacard.writeDatacards();

  std::cout << "\n=== Finishing with the datacard production\n\n";
}





