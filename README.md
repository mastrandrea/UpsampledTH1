# UpsampledTH1

Implementation of TH1 for upsampled datasets (multiple reconstructions for each generated event) able to run with RDataFrame and compatible with the Root implementation of multithreading execution.

## Definitions

 - event ID          = unique identifier for each generated event
 - Upsampling factor = # of entries in the dataset reconstructed from the same generated event (same event ID)
 - Upsampling weight = event weight due to the upsampling (upsampling_weight = 1.0/upsampling_factor)
 - Incremental step  = step for event eventual accumulation in a thread
 - Event step        = number of events between 2 consecutive accumulations ( event_step = incremental_step * upsampling_factor)

If upsampling_factor is <= 1, the upsampling_weight is set to 1.
If the incremental_step is <= 0, the event_step is set to 0, and no incremental accumulation is attemped.

## Implementation

For each thread the following structures are instantiated:
  - a data buffer
      - the buffer is implemented as a  std::map<long int, std::vector<ValAndWeightT>> (the map key is the event ID).
      - the ValAndWeightT type is defined as  std::pair<ValT, float>  (where Val is for the variable used to fill the histogram, and the second value is the entry weight).
  - a temporary TH1 (used to build the values distribution for a single event ID)
  - an incremental TH1 (used to accumulate in a thread the distributions for each event ID - in case all the entries with the same event ID were processed by a single thread)

## Strategy

For each event the value and weight are stored in the buffer of the corresponding thread.
The weight includes the "upsampling weight".
If incremental_step is >= 1, every event_step entries an accumulation is attempted: all the vectors of ValAndWeight, in the buffer buffer of the current thread, have size = upsampling_factor, then the corresponding temporary histogram is reset, filled with the values and weigths in the vector. This distribution is then imported in the incremental histogram for the thread (** loop from n_bin-1 to n_bin+1, if bin content > 0, then the incremental histogram is filled with (bin_center, bin_content), with bin_content from the temporary histo **).
At the end of the loop over the whole dataset, all remaining buffers are merged and accumulated in the incremental histograms.
Finally the incremental histograms for all the threads are added (TH1->Add) in the final TH1.


## Constructor:

  UpsampledTH1(std::string_view name, std::string_view title, int n_bin, double x_min, double x_max, int upsampling_factor = 1, int incremental_step = 0)

 name  = TH1 name
 title = TH1 title
 n_bin = # bins
 x_min = minimum of the x-axis
 x_max = maximum of the x-axis
 upsampling_factor = # of entries in the dataset reconstructed from the same generated event (same event ID)
 incremental_step  = step for event eventual accumulation in a thread


## Example declaration for an RDataFrame

  auto testTH1F = rdf.Book<long, float>(UpsampledTH1<TH1F, float>{"testTH1F", "Upsampled TH1F", 20, 0., 100., nFold, iStep}, {evID_variable, evValue_variable});



# Testing recipe

```
git clone git@github.com:mastrandrea/UpsampledTH1.git
download test_oversampling.root file
g++ `root-config --cflags` `root-config --libs` UpsampledTH1.C -o test

./test

(or time ./test))

```
