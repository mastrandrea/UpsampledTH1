#include <map>
#include <utility>
#include <iostream>
#include <vector>
#include <string>
#include <string_view>
#include "TH1F.h"
#include "TCanvas.h"
#include "TROOT.h"
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RResultPtr.hxx>


//////////////////////////////////////////////////
//
// INCREMENTAL version:
//     Each thread can periodically examines its buffer and "compact" the elements with the same genEvent in a temporary histo if the number of replicas is equal to the upsampling factor.
//     This implementation can shrink the memory footprint, in particular for situation with (quasi-)sequential identical genEvents and a small fraction of partially-selected genEvents.
//     (In the case of partially-selected genEvents, the buffer will be kept till the finalization of the task)
//
// OPEN ISSUES:
//   - if the value is missing, it is interpreted as present with a value of 0 (this might be a RDF "feature" - to be checked)
//       -> Fixed by adding a (logical) filter step to the RDF :)
//
//////////////////////////////////////////////////


// The template class inherits from ROOT::Detail::RDF::RActionImpl

template <typename TH1, typename ValT>
class UpsampledTH1 : public ROOT::Detail::RDF::RActionImpl<UpsampledTH1<TH1, ValT>> {

public:
  // Define the type of the histogram
  using Result_t         = TH1;

  using ValAndWeightT    = std::pair<ValT, float>;

  using UpsampledBufferT = std::map<long int, std::vector<ValAndWeightT>>;  // One vector of value_&_weight per event ID


private:
  // Number of slots (threads)
  int nSlots;

  // Upsampling factor - i.e. number of replicas per event (ID)
  int    UpsamplingFactor = 1;
  double UpsamplingWeight = 1.0;

  // Configuration parameters
  int IncrementalStep     = 0;
  int eventStep           = 0;
  int totElements         = 0;
  std::vector<int> processedEntries;

  // Histogram parameters
  std::string histo_name;
  int         nbin;
  double      xmin;
  double      xmax;

  // Map: 1 buffer per slot
  std::map<int, UpsampledBufferT>   fBuffers;              // One buffer per slot (vector of values for each eventID)


  // Shared pointer to the final histogram
  std::shared_ptr<TH1>              fFinalHisto;

  // Vector of shared pointer to the temporary histograms used to shrink the accumulate the incremental per-thread result
  std::vector<std::shared_ptr<TH1>> incrementalHistos;

  // Vector of shared pointer to the temporary histograms used to accumulate a single upsampled event
  std::vector<std::shared_ptr<TH1>> tempHistos;



public:

  // Constructor with parameters for histogram creation
  UpsampledTH1(std::string_view name, std::string_view title, int n_bin, double x_min, double x_max, int upsampling_factor = 1, int incremental_step = 0)
    : nSlots(ROOT::IsImplicitMTEnabled() ? ROOT::GetThreadPoolSize() : 1) {

    UpsamplingFactor = upsampling_factor;
    if (UpsamplingFactor > 1) UpsamplingWeight = (double)(1.0 / ((double)UpsamplingFactor));

    IncrementalStep = incremental_step;
    if (IncrementalStep > 0) eventStep = IncrementalStep * UpsamplingFactor;

    histo_name = std::string(name);
    nbin       = n_bin;
    xmin       = x_min;
    xmax       = x_max;

    fFinalHisto = std::make_shared<TH1>(histo_name.c_str(), std::string(title).c_str(), nbin, xmin, xmax);
    fFinalHisto->Sumw2();

    for (int i = 0; i < nSlots; i++) {
      fBuffers.emplace(i, UpsampledBufferT());

      std::string hN_slot      = histo_name + "_" + std::to_string(i);
      std::string hN_slot_temp = hN_slot + "_temp";

      std::shared_ptr<TH1> sH(std::make_shared<TH1>(hN_slot.c_str(), std::string(title).c_str(), nbin, xmin, xmax));
      sH->Sumw2();
      incrementalHistos.push_back(sH);

      std::shared_ptr<TH1> tH(std::make_shared<TH1>(hN_slot_temp.c_str(), std::string(title).c_str(), nbin, xmin, xmax));
      tH->Sumw2();
      tempHistos.push_back(tH);

      processedEntries.push_back(0);
    }


    std::cout << std::endl;
    std::cout << "[ UpsampledTH1 ] nSlots            = " << nSlots           << std::endl;
    std::cout << "[ UpsampledTH1 ] Upsampling Factor = " << UpsamplingFactor << std::endl;
    std::cout << "[ UpsampledTH1 ] Upsampling Weight = " << UpsamplingWeight << std::endl;
    std::cout << "[ UpsampledTH1 ] Incremental Step  = " << IncrementalStep  << std::endl;
    std::cout << "[ UpsampledTH1 ] Event Step        = " << eventStep        << std::endl;
    std::cout << std::endl;

    for (int i = 0; i < nSlots; i++) {
      std::cout << "[ UpsampledTH1 ] fBuffers[" << i << "].size() = " << fBuffers[i].size()    << std::endl;
      std::cout << std::endl;

      std::cout << "[ UpsampledTH1 ] tempHistos[" << i << "]->Print() :   nBins = " << tempHistos[i]->GetNbinsX() << std::endl;
      tempHistos[i]->Print();
      std::cout << std::endl;
    }
  }


  // Move constructor and deletion of copy constructor
  UpsampledTH1(UpsampledTH1 &&)      = default;
  UpsampledTH1(const UpsampledTH1 &) = delete;


  // Function to return the pointer to the final histogram
  std::shared_ptr<TH1> GetResultPtr() const { return fFinalHisto; }


  // Class initialization - only 1 per event loop
  void Initialize() {
    std::cout << "[ UpsampledTH1::Initialize ] fBuffers.size()    = " << fBuffers.size()    << std::endl;
    std::cout << "[ UpsampledTH1::Initialize ] fBuffers[0].size() = " << fBuffers[0].size() << std::endl;
    return;
  }


  // Task initialization - 1 per task (typically 2 tasks per thread)
  void InitTask(TTreeReader *, unsigned int slot) {
    std::cout << "[ UpsampledTH1::InitTask ] InitTask for slot = " << slot << std::endl;
    return;
  }


  // Exec function is called for each event and fills the buffer for the corresponding slot
  // Burrefed events are accumulated in the incremental histos if the processed entries (for the slot) is mod eventStep (= upsampling_factor * incremental_step)
  void Exec(unsigned int slot, unsigned long eventID, ValT value, float weight = 1.0) {

    processedEntries[slot]++;
    //if (eventID < 100) std::cout << "[ Exec  " << slot << " ] processed entries = " << processedEntries[slot] << "  /  eventId =  " << eventID << std::endl;

    UpsampledBufferT* buf = &fBuffers[slot];


    // Value buffering -----------------------

    if ( buf->find(eventID) == buf->end() ) {

      std::vector<ValAndWeightT> vv = {ValAndWeightT(value, weight)};
      buf->emplace(eventID, vv);

    } else {

      (*buf)[eventID].push_back(ValAndWeightT(value, weight));
    }


    // Incremental accumulation  ---------------------

    if ( (IncrementalStep > 0)  &&  ((processedEntries[slot]%eventStep) == 0) ) {

      // Loop over the bufferized events
      for (auto eV = (*buf).begin(); eV != (*buf).end(); ) {

	std::vector<ValAndWeightT> vVal = eV->second;

	// If an event is "complete", accumulate that and remove its entry from the buffer
	if (vVal.size() == UpsamplingFactor) {

	  getUpsampledEventDistribution(tempHistos[slot], vVal);

	  fillTH1withUpsampledEventDistribution(incrementalHistos[slot], tempHistos[slot]);

	  // erease element in buf - returned pointer to next location
	  eV = buf->erase(eV);

	// If the event is still "incomplete", simply move on the iterator to the next element 
	} else {

	  ++eV;
	}
      }
    }

    return;
  }


  void mergeBuffers() {

    std::cout << "[ UpsampledTH1::mergeBuffers ] mergeBuffers called " << std::endl;

    UpsampledBufferT* baseBuf = &fBuffers[0];

    for (int s = 1; s < nSlots; s++) {

      UpsampledBufferT* buf = &fBuffers[s];

      for ( const auto& eV : (*buf) ) {

	int eID                         = eV.first;
	std::vector<ValAndWeightT> vVal = eV.second;

	if ( baseBuf->find(eID) == baseBuf->end() ) {

	  baseBuf->emplace(eID, vVal);

	} else {

	  for (int i = 0; i < vVal.size(); i++) {  (*baseBuf)[eID].push_back(vVal[i]);  }
	}
      }
    }

    return;
  }



  // Fill temporary histo for 1 event ID
  void getUpsampledEventDistribution(std::shared_ptr<TH1> hTemp, std::vector<ValAndWeightT> vw) {

    hTemp->Reset();

    for (int i = 0; i < vw.size(); i++) { hTemp->Fill(vw[i].first, (vw[i].second * UpsamplingWeight)); }

    return;
  }


  // Fill the target histo with the *distribution* for 1 event ID
  void fillTH1withUpsampledEventDistribution(std::shared_ptr<TH1> hTarget, std::shared_ptr<TH1> hSource) {      

    for (int ib = 0; ib < (nbin+2); ib++) {        // This loop includes underflow and overflow bins

      double temp_v = hSource->GetBinContent(ib);

      if (temp_v > 0.0) {
	double temp_x = hSource->GetBinCenter(ib);
	hTarget->Fill(temp_x, temp_v);
      }
    }

    return;
  }




  void mergeUpsampledTH1() {

    std::cout << "[ UpsampledTH1::mergeUpsampledTH1 ] mergeUpsampleTH1 called " << std::endl;

    UpsampledBufferT* baseBuf = &fBuffers[0];

    totElements = 0;

    // Loop over all events to collect al the entries not yet "compacted" in the incremental histos
    // (at this stage all buffers are merged into baseBuf - for each event ID there is a vector of Val&Weight) 

    for (auto eV = (*baseBuf).begin(); eV != (*baseBuf).end(); ) {

      std::vector<ValAndWeightT> vVal = eV->second;

      // Use the slot 0

      getUpsampledEventDistribution(tempHistos[0], vVal);

      fillTH1withUpsampledEventDistribution(incrementalHistos[0], tempHistos[0]);

      // erease element in buf - returned pointer to next location
      eV = baseBuf->erase(eV);
    }

    // Sum all the incremental histograms in the final histo
    for (int s = 0; s < nSlots; s++) {

      fFinalHisto->Add( (incrementalHistos[s]).get() );

      totElements += processedEntries[s];
    }

    return;
  }




  void Finalize() {

    std::cout << "[ UpsampledTH1::Finalize ] Finalize called " << std::endl;

    std::cout << "[ UpsampledTH1::Finalize ]  fBuffers.size()    = " << fBuffers.size() << std::endl;

    for (int i = 0; i < nSlots; i++) {

      std::cout << "[ UpsampledTH1::Finalize ]  fBuffers[" << i << "].size() = " << fBuffers[i].size() << std::endl;
    }


    mergeBuffers();

    mergeUpsampledTH1();

    std::cout << std::endl << "[ UpsampledTH1::Finalize ]  TOT PROCESSED ENTRIES =  " << totElements << std::endl << std::endl;

    return;
  }

  // 
  std::string GetActionName() { return "UpsampledTH1"; }
};



// // RResultPtr<typename std::decay_t<UpsampledTH1<TH1F, float>>::Result_t> UpsampledTH1F(std::string_view hName,
// // auto UpsampledTH1F(ROOT::RDataFrame _rdf,
// auto UpsampledTH1F(ROOT::RDF::RInterface _rdf,
// 		   std::string_view hName,
// 		   std::string hTitle,
// 		   int hNbins, double hXmin, double hXmax,
// 		   int hNfold,
// 		   int hIstep,
// 		   std::string evIDVar,
// 		   std::string evValVar){

//   return _rdf.Book<long, float>(UpsampledTH1<TH1F, float>{"testTH1F", "Upsampled TH1F", 20, 0., 100., hNfold, hIstep},
// 			      {evIDVar, evValVar});
// }



//void UpsampledHisto() {
int main() {

  ROOT::EnableImplicitMT(4);

  ROOT::RDataFrame rdf{"Events", "test_oversampling.root"};

  auto dd = rdf.Filter("nJet > 0").Define("FirstJet_pt", "Jet_pt[0]");  // Adding the Filter step avoids the value for Jet_pt[0] = 0 when there is no jet :)


  int nFold = 4;
  int iStep = 10;

  std::string evIDVar( "genEventProgressiveNumber");
  std::string evValVar("FirstJet_pt");

  // Book the upsampled TH1F 

  // auto testTH1F = dd.Book<long, float>(UpsampledTH1<TH1F, float>{"testTH1F", "Upsampled TH1F", 20, 0., 100., nFold, iStep},
  // 				       {"genEventProgressiveNumber", "FirstJet_pt"});

  auto testTH1F = dd.Book<long, float>(UpsampledTH1<TH1F, float>{"testTH1F", "Upsampled TH1F", 20, 0., 100., nFold, iStep},
  				       {evIDVar, evValVar});

  // Simpler syntax - not working yet
  // auto testTH1F = UpsampledTH1F(dd, "testTH1F", "Upsampled TH1F", 20, 0., 100., nFold, iStep, evIDVar, evValVar);
  

  // Trigger the execution
  testTH1F->Print();


  // Print a few check info
  std::cout << " ------------- Integral  =  " << testTH1F->Integral()        << std::endl;
  std::cout << " ------------- Underflow =  " << testTH1F->GetBinContent(0)  << std::endl;
  std::cout << " ------------- Overflow  =  " << testTH1F->GetBinContent(testTH1F->GetNbinsX()+1) << std::endl;

  TCanvas* c1 = new TCanvas();

  testTH1F->Draw();

  c1->SaveAs("Up-test.pdf");


  std::cout << " ------------- result type :  " << testTH1F->Class_Name() << std::endl;

  return 0;
}
