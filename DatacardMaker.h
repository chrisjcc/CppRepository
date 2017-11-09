#ifndef DatacardMaker_h
#define DatacardMaker_h

#include <vector>
#include <set>
#include <map>
#include <string>

#include <TString.h>
#include <TFile.h>

//class TLegend;
class RootFileReader;
class TH1;

#include "plotterHelpers.h"
#include "SamplesFwd.h"
#include "Sample.h"
#include "../../common/include/sampleHelpers.h"

#include <fstream>

class AnalysisConfig;





class DatacardMaker{
  
 public:
  
  /// Constructor
  DatacardMaker(const AnalysisConfig& analysisConfig,
		const std::vector<std::string>& v_plot,
		const std::vector<Channel::Channel>& v_channel,
		const std::vector<Systematic::Systematic>& v_systematic,
		const std::string& fileLists,
		const std::string& configname);

  /// Destructor
  ~DatacardMaker(){};
    
  /// Write the datacards for limit setting tool for all mva config
  void writeDatacards();
  
  /// Write datacard header
  void writeHeader(const std::string& name);

  /// Write histogram and body of datacard
  void writeSystematicUncertainties();

  /// Handles datacard header and signal, observed, and background yield printouts
  void extractYields(const std::string& name);
  
  /// Handles writing histogram to output root file 
  void writeYields(const std::string& name);

  void setIncludeSystmeticUncertainties(bool useSys);
  void setIncludeStatisticalUncertainties(bool useStat);

 private:
   
   /// Pair of a legend entry and the histogram for the corresponding sample
   typedef std::pair<TH1*, TH1*> HistoPair;
   typedef std::map<Systematic::Type, HistoPair> SystematicHistoMap;
   
   /// Process datacard writer
   void writeVariations(const SystematicHistoMap& histoCollection, const Channel::Channel channel, const std::string processName);

   /// Either creates or addes to the current histogram
   TH1* addOrCreateHisto(TH1* base, const TH1* const add_histo) const;

   /// Signal model (either SM or BSM)
   const char* signalModel_;

   /// Path to input files
   const TString fileList_base_;

   /// Path to input steering parameter file
   const std::string& configname_;

   /// Output folder name
   const char* outputBaseDir_;

   /// Output file name
   const char* outputFileName_;

   /// Include systematic uncertainty in datacard
   bool addSystematicUncertainty_;

   /// Include statistical uncertainty in datacard 
   bool addStatisticalUncertainty_;

   /// Reference to the analysis config                                                                                          
   const AnalysisConfig& analysisConfig_;

   /// Samples to be analysed
   std::map<Channel::Channel, std::map<Systematic::Systematic, std::map<std::string, std::string > > >* inputFileLists_;

   /// Write statistical uncertainties
   void writeStatisticalUncertainties(const std::string& name);

   /// Convert internal systematic uncertainty labeling to CMS ttH collaboration naming convention
   std::map<std::string, std::string> convertSampleNames_;

   /// Convert internal systematic uncertainty labeling to CMS ttH collaboration naming convention
   std::map<std::string, std::string> convertSystematicLabel_;

   /// Event categories based on the # of jets and # of b-tagged jets
   std::map<std::string, std::string> mapOfCategories_;

   /// File reader for accessing specific histogram from given file
   RootFileReader* fileReader_;

   /// Vector of process names obtained from steering parameter file
   std::vector<std::string> processNames_;
   
   /// Vector of root file names (i.e. corresponding to the BDT configs)
   std::vector<std::string> fileNames_;

   /// Complete directory path and input output file name
   std::string eventFileString_;
      
   /// File name of the datacard
   std::string datacardName_;
   
   /// Datacard output file handler
   std::ofstream datacard_;
   
   /// Datacard output directory
   std::string outputDirDatacard_;
   
   /// Datacard file hanlder
   TFile* outputFile_;
   
   /// Obserable type used in analysis (i.e. event classification)
   std::string observableType_;

   /// Set to true to apply MC bin-by-bin statistical shape uncertainty pruning
   bool pruneBinByBin_;
   
   /// Assign lnN type systematic value based on process type
   std::map<std::string,std::vector<std::pair<std::string, std::string>>> valueOfSystematicBasedOnProcess_;

   /// List of shape systematic uncertainties
   std::vector<std::string> listOfShapeSystematics_;

   /// List of lnN systematic uncertainties
   std::vector<std::string> listOflnNSystematics_;

   /// Alternative list of systematic uncertainty by type
   std::vector<std::pair<std::string,std::string>> listOfSystematicsByType_;

   /// Initialzer function for 13 TeV, to assign lnN systematic uncertainty values
   void initialization13TeV(bool pruneOption);

   /// Initial input from produceDatacards
   const std::vector<std::string>& v_plot_;
   const std::vector<Channel::Channel>& v_channel_;
   const std::vector<Systematic::Systematic>& v_systematic_;

  

};




#endif


