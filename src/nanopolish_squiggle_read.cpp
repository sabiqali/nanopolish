//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_squiggle_read -- Class holding a squiggle (event)
// space nanopore read
//
#include <algorithm>
#include "nanopolish_common.h"
#include "nanopolish_squiggle_read.h"
#include "nanopolish_pore_model_set.h"
#include "nanopolish_methyltrain.h"
#include "nanopolish_extract.h"
#include "nanopolish_raw_loader.h"
#include "nanopolish_fast5_io.h"
#include "nanopolish_fast5_loader.h"

extern "C" {
#include "event_detection.h"
#include "scrappie_common.h"
}

#include <fast5.hpp>

//#define DEBUG_MODEL_SELECTION 1
//#define DEBUG_RECONSTRUCTION 1

// Track the number of skipped reads to warn the use at the end of the run
// Workaround for albacore issues.  Temporary, I hope
int g_total_reads = 0;
int g_unparseable_reads = 0;
int g_qc_fail_reads = 0;
int g_failed_calibration_reads = 0;
int g_failed_alignment_reads = 0;
int g_bad_fast5_file = 0;

const double MIN_CALIBRATION_VAR = 2.5;

void SquiggleScalings::set4(double _shift,
                            double _scale,
                            double _drift,
                            double _var)
{
    set6(_shift, _scale, _drift, _var, 1.0, 1.0);
}

void SquiggleScalings::set6(double _shift,
                            double _scale,
                            double _drift,
                            double _var,
                            double _scale_sd,
                            double _var_sd)
{
    // direct
    shift = _shift;
    scale = _scale;
    drift = _drift;
    var = _var;
    scale_sd = _scale_sd;
    var_sd = _var_sd;

    // derived
    log_var = log(var);
    scaled_var = var / scale;
    log_scaled_var = log(scaled_var);
}

//
SquiggleRead::SquiggleRead(const std::string& name, const ReadDB& read_db, const uint32_t flags)
{
    this->fast5_path = read_db.get_signal_path(name);
    g_total_reads += 1;
    if(this->fast5_path == "") {
        g_bad_fast5_file += 1;
        return;
    }

    std::string sequence = read_db.get_read_sequence(name);
    Fast5Data data = Fast5Loader::load_read(fast5_path, name);
    if(data.is_valid && !sequence.empty()) {
        init(sequence, data, flags);
    } else {
        fprintf(stderr, "[warning] fast5 file is unreadable and will be skipped: %s\n", fast5_path.c_str());
        g_bad_fast5_file += 1;
    }

    if(!this->events[0].empty()) {
        assert(this->base_model[0] != NULL);
    }
    free(data.rt.raw);
    data.rt.raw = NULL;
}

SquiggleRead::SquiggleRead(const ReadDB& read_db, const Fast5Data& data, const uint32_t flags)
{
    init(read_db.get_read_sequence(data.read_name), data, flags);
}

SquiggleRead::SquiggleRead(const std::string& sequence, const Fast5Data& data, const uint32_t flags)
{
    init(sequence, data, flags);
}

//
void SquiggleRead::init(const std::string& read_sequence, const Fast5Data& data, const uint32_t flags)
{
    this->nucleotide_type = SRNT_DNA;
    this->pore_type = PT_UNKNOWN;
    this->f_p = nullptr;

    this->events_per_base[0] = events_per_base[1] = 0.0f;
    this->base_model[0] = this->base_model[1] = NULL;
    g_total_reads += 1;

    this->read_name = data.read_name;
    this->read_sequence = read_sequence;

    // sometimes the basecaller will emit very short sequences, which causes problems
    // also there can be rare issues with the signal in the fast5 and we want to skip
    // such reads
    if(this->read_sequence.length() > 20 && data.is_valid && data.rt.n > 0) {
        load_from_raw(data, flags);
    } else {
        g_bad_fast5_file += 1;
    }

    if(!this->events[0].empty()) {
        assert(this->base_model[0] != NULL);
    }
}

SquiggleRead::~SquiggleRead()
{

}

// helper for get_closest_event_to
int SquiggleRead::get_next_event(int start, int stop, int stride, uint32_t strand) const
{
    while(start != stop) {

        int ei = base_to_event_map[start].indices[strand].start;
        if(ei != -1)
            return ei;
        start += stride;
    }
    return -1;
}

//
int SquiggleRead::get_closest_event_to(int k_idx, uint32_t strand) const
{
    int stop_before = std::max(0, k_idx - 1000);
    int stop_after = std::min(k_idx + 1000, (int)base_to_event_map.size() - 1);

    int event_before = get_next_event(k_idx, stop_before, -1, strand);
    int event_after = get_next_event(k_idx, stop_after, 1, strand);

    // TODO: better selection of "best" event to return
    if(event_before == -1)
        return event_after;
    return event_before;
}

//
/*void SquiggleRead::load_from_events(const uint32_t flags)
{
    assert(this->nucleotide_type != SRNT_RNA);

    assert(f_p->is_open());
    detect_pore_type();
    detect_basecall_group();
    assert(not basecall_group.empty());

    read_sequence = f_p->get_basecall_seq(read_type, basecall_group);

    // Load PoreModel for both strands
    std::vector<EventRangeForBase> event_maps_1d[NUM_STRANDS];
    std::string read_sequences_1d[NUM_STRANDS];

    for (size_t si = 0; si < 2; ++si) {

        // Do we want to load this strand?
        if(! (read_type == SRT_2D || read_type == si) ) {
            continue;
        }

        // Load the events for this strand
        auto f5_events = f_p->get_basecall_events(si, basecall_group);

        // copy events
        events[si].resize(f5_events.size());
        std::vector<double> p_model_states;

        for(size_t ei = 0; ei < f5_events.size(); ++ei) {
            auto const & f5_event = f5_events[ei];

            events[si][ei] = { static_cast<float>(f5_event.mean),
                               static_cast<float>(f5_event.stdv),
                               f5_event.start,
                               static_cast<float>(f5_event.length),
                               static_cast<float>(log(f5_event.stdv))
                             };
            assert(f5_event.p_model_state >= 0.0 && f5_event.p_model_state <= 1.0);
            p_model_states.push_back(f5_event.p_model_state);
        }

        // we need the 1D event map and sequence to calculate calibration parameters
        // these will be copied into the member fields later if this is a 1D read,
        // or discarded if this is a 2D read

        // NB we use event_group in this call rather than basecall_group as we want the 1D basecalls that match the events
        read_sequences_1d[si] = f_p->get_basecall_seq(si == 0 ? SRT_TEMPLATE : SRT_COMPLEMENT,
                                                      f_p->get_basecall_1d_group(basecall_group));
        event_maps_1d[si] = build_event_map_1d(read_sequences_1d[si], si, f5_events, 5);

        // Constructing the event map can fail due to an albacore bug
        // in this case, we have to set this strand to be invalid
        if(!event_maps_1d[si].empty()) {
            // run version-specific load
            if(pore_type == PT_R7) {
                _load_R7(si);
            } else {
                _load_R9(si, read_sequences_1d[si], event_maps_1d[si], p_model_states, flags);
            }
        } else {
            events[si].clear();
        }
    }

    // Build the map from k-mers of the read sequence to events
    if(read_type == SRT_2D) {
        if(pore_type == PT_R9) {
            build_event_map_2d_r9();
        } else {
            assert(pore_type == PT_R7);
            build_event_map_2d_r7();
        }
    } else {
        assert(read_type < NUM_STRANDS);
        this->base_to_event_map.swap(event_maps_1d[read_type]);
    }

    // Load raw samples if requested
    if(flags & SRF_LOAD_RAW_SAMPLES) {

        auto& sample_read_names = f_p->get_raw_samples_read_name_list();
        if(sample_read_names.empty()) {
            fprintf(stderr, "Error, no raw samples found\n");
            exit(EXIT_FAILURE);
        }

        // we assume the first raw sample read is the one we're after
        std::string sample_read_name = sample_read_names.front();

        samples = f_p->get_raw_samples(sample_read_name);
        sample_start_time = f_p->get_raw_samples_params(sample_read_name).start_time;

        // retrieve parameters
        auto channel_params = f_p->get_channel_id_params();
        sample_rate = channel_params.sampling_rate;
    }

    // Filter poor quality reads that have too many "stays"
    if(!events[0].empty() && events_per_base[0] > 5.0) {
        g_qc_fail_reads += 1;
        events[0].clear();
        events[1].clear();
    }
}*/

//
void SquiggleRead::load_from_raw(const Fast5Data& fast5_data, const uint32_t flags)
{

    // Try to detect whether this read is DNA or RNA
    // Fix issue 531: experiment_type in fast5 is "rna" for cDNA kit dcs108
    bool rna_experiment = fast5_data.experiment_type == "rna" || fast5_data.experiment_type == "internal_rna";
    this->nucleotide_type = rna_experiment && fast5_data.sequencing_kit != "sqk-dcs108" ? SRNT_RNA : SRNT_DNA;

    //
    this->read_type = SRT_TEMPLATE;
    this->pore_type = PT_R9;
    std::string strand_str = "template";
    size_t strand_idx = 0;

    // default to R9 parameters
    std::string alphabet = "nucleotide";
    std::string kit = "r9.4_450bps";
    size_t k = 6;
    const detector_param* ed_params = &event_detection_defaults;

    /*if(this->pore_type == PT_R10) {
        kit = "r10_450bps";
        k = 9;
        ed_params = &event_detection_r10;
    }*/

    if(this->nucleotide_type == SRNT_RNA) {
        assert(this->pore_type == PT_R9);
        kit = "r9.4_70bps";
        alphabet = "u_to_t_rna";
        k = 5;
        ed_params = &event_detection_rna;

        std::replace(this->read_sequence.begin(), this->read_sequence.end(), 'U', 'T');
    }


    // Set the base model for this read to either the nucleotide or U->T RNA model
    this->base_model[strand_idx] = PoreModelSet::get_model(kit, alphabet, strand_str, k);
    assert(this->base_model[strand_idx] != NULL);

    // Read the sample rate
    this->sample_rate = fast5_data.channel_params.sample_rate;
    this->channel_id = fast5_data.channel_params.channel_id;
    this->sample_start_time = fast5_data.start_time;

    // trim raw using scrappie's internal method
    // parameters taken directly from scrappie defaults
    int trim_start = 200;
    int trim_end = 10;
    int varseg_chunk = 100;
    float varseg_thresh = 0.0;
    trim_and_segment_raw(fast5_data.rt, trim_start, trim_end, varseg_chunk, varseg_thresh);
    event_table et = detect_events(fast5_data.rt, *ed_params);
    assert(et.n > 0);

    //
    this->scalings[strand_idx] = estimate_scalings_using_mom(this->read_sequence,
                                                             *this->base_model[strand_idx],
                                                             et);


    // copy events into nanopolish's format
    this->events[strand_idx].resize(et.n);
    double start_time = 0;
    for(size_t i = 0; i < et.n; ++i) {
        float length_in_seconds = et.event[i].length / this->sample_rate;
        this->events[strand_idx][i] = { et.event[i].mean, et.event[i].stdv, start_time, length_in_seconds, logf(et.event[i].stdv) };
        start_time += length_in_seconds;
    }

    if(flags & SRF_LOAD_RAW_SAMPLES) {
        this->sample_start_time = 0;
        this->samples.resize(fast5_data.rt.n);
        for(size_t i = 0; i < this->samples.size(); ++i) {
            assert(fast5_data.rt.start + i < fast5_data.rt.n);
            this->samples[i] = fast5_data.rt.raw[fast5_data.rt.start + i];
        }
    }

    // If sequencing RNA, reverse the events to be 5'->3'
    if(this->nucleotide_type == SRNT_RNA) {
        std::reverse(this->events[strand_idx].begin(), this->events[strand_idx].end());
    }

    // clean up event tables
    assert(et.event != NULL);
    free(et.event);

    // align events to the basecalled read
    std::vector<AlignedPair> event_alignment = adaptive_banded_simple_event_align(*this, *this->base_model[strand_idx], read_sequence);

    // transform alignment into the base-to-event map
    if(event_alignment.size() > 0) {

        // create base-to-event map
        size_t n_kmers = read_sequence.size() - this->get_model_k(strand_idx) + 1;
        this->base_to_event_map.clear();
        this->base_to_event_map.resize(n_kmers);

        size_t max_event = 0;
        size_t min_event = std::numeric_limits<size_t>::max();

        size_t prev_event_idx = -1;
        for(size_t i = 0; i < event_alignment.size(); ++i) {

            size_t k_idx = event_alignment[i].ref_pos;
            size_t event_idx = event_alignment[i].read_pos;
            IndexPair& elem = this->base_to_event_map[k_idx].indices[strand_idx];
            if(event_idx != prev_event_idx) {
                if(elem.start == -1) {
                    elem.start = event_idx;
                }
                elem.stop = event_idx;
            }

            max_event = std::max(max_event, event_idx);
            min_event = std::min(min_event, event_idx);
            prev_event_idx = event_idx;
        }

        events_per_base[strand_idx] = (double)(max_event - min_event) / n_kmers;

        // prepare data structures for the final calibration
        std::vector<EventAlignment> alignment =
            get_eventalignment_for_1d_basecalls(read_sequence, alphabet, this->base_to_event_map, this->base_model[strand_idx]->k, strand_idx, 0);

        // run recalibration to get the best set of scaling parameters and the residual
        // between the (scaled) event levels and the model.
        // internally this function will set shift/scale/etc of the pore model
        bool calibrated = recalibrate_model(*this, *this->base_model[strand_idx], strand_idx, alignment, true, false);

#ifdef DEBUG_MODEL_SELECTION
        fprintf(stderr, "[calibration] read: %s events: %zu"
                         " scale: %.2lf shift: %.2lf drift: %.5lf var: %.2lf\n",
                                read_name.substr(0, 6).c_str(), this->events[strand_idx].size(), this->scalings[strand_idx].scale,
                                this->scalings[strand_idx].shift, this->scalings[strand_idx].drift, this->scalings[strand_idx].var);
#endif

        // QC calibration
        if(!calibrated || this->scalings[strand_idx].var > MIN_CALIBRATION_VAR) {
            events[strand_idx].clear();
            g_failed_calibration_reads += 1;
        }
    } else {
        // Could not align, fail this read
        this->events[strand_idx].clear();
        this->events_per_base[strand_idx] = 0.0f;
        g_failed_alignment_reads += 1;
    }

    // Filter poor quality reads that have too many "stays"
    if(!this->events[strand_idx].empty() && this->events_per_base[strand_idx] > 5.0) {
        g_qc_fail_reads += 1;
        events[0].clear();
        events[1].clear();
    } 
}

// Return a vector of eventalignments for the events that made up the basecalls in the read
std::vector<EventAlignment> SquiggleRead::get_eventalignment_for_1d_basecalls(const std::string& read_sequence_1d,
                                                                              const std::string& alphabet_name,
                                                                              const std::vector<EventRangeForBase>& base_to_event_map_1d,
                                                                              const size_t k,
                                                                              const size_t strand_idx,
                                                                              const int shift_offset) const
{
    std::vector<EventAlignment> alignment;

    const Alphabet* alphabet = get_alphabet_by_name(alphabet_name);
    size_t n_kmers = read_sequence_1d.size() - k + 1;
    size_t prev_kmer_rank = -1;

    for(int ki = 0; ki < n_kmers; ++ki) {
        IndexPair event_range_for_kmer = base_to_event_map_1d[ki].indices[strand_idx];

        // skip kmers without events
        if(event_range_for_kmer.start == -1)
            continue;

        // skip k-mers that cannot be shifted to a valid position
        if(ki + shift_offset < 0 || ki + shift_offset >= n_kmers) {
            continue;
        }

        for(size_t event_idx = event_range_for_kmer.start;
            event_idx <= event_range_for_kmer.stop; event_idx++)
        {
            assert(event_idx < this->events[strand_idx].size());

            // since we use the 1D read seqence here we never have to reverse complement
            std::string kmer = read_sequence_1d.substr(ki + shift_offset, k);
            size_t kmer_rank = alphabet->kmer_rank(kmer.c_str(), k);

            EventAlignment ea;
            // ref data
            //ea.ref_name = "read";
            ea.read_idx = -1; // not needed
            ea.ref_kmer = kmer;
            ea.ref_position = ki;
            ea.strand_idx = strand_idx;
            ea.event_idx = event_idx;
            ea.rc = false;
            ea.model_kmer = kmer;
            ea.hmm_state = prev_kmer_rank != kmer_rank ? 'M' : 'E';
            alignment.push_back(ea);
            prev_kmer_rank = kmer_rank;
        }
    }

    return alignment;
}

size_t SquiggleRead::get_sample_index_at_time(size_t sample_time) const
{
    return sample_time - this->sample_start_time;
}

//
std::vector<float> SquiggleRead::get_scaled_samples_for_event(size_t strand_idx, size_t event_idx) const
{
    std::pair<size_t, size_t> sample_range = get_event_sample_idx(strand_idx, event_idx);

    std::vector<float> out;
    for(size_t i = sample_range.first; i < sample_range.second; ++i) {
        double curr_sample_time = (this->sample_start_time + i) / this->sample_rate;
        //fprintf(stderr, "event_start: %.5lf sample start: %.5lf curr: %.5lf rate: %.2lf\n", event_start_time, this->sample_start_time / this->sample_rate, curr_sample_time, this->sample_rate);
        double s = this->samples[i];
        // apply scaling corrections
        double scaled_s = s - this->scalings[strand_idx].shift;
        assert(curr_sample_time >= (this->sample_start_time / this->sample_rate));
        scaled_s -= (curr_sample_time - (this->sample_start_time / this->sample_rate)) * this->scalings[strand_idx].drift;
        scaled_s /= this->scalings[strand_idx].scale;
        out.push_back(scaled_s);
    }
    return out;
}

// return a pair of value corresponding to the start and end index of a given index on the signal
std::pair<size_t, size_t> SquiggleRead::get_event_sample_idx(size_t strand_idx, size_t event_idx) const
{
    double event_start_time = this->events[strand_idx][event_idx].start_time;
    double event_duration = this->events[strand_idx][event_idx].duration;

    size_t start_idx = this->get_sample_index_at_time(event_start_time * this->sample_rate);
    size_t end_idx = this->get_sample_index_at_time((event_start_time + event_duration) * this->sample_rate);

    return std::make_pair(start_idx, end_idx);
}

double median(const std::vector<double>& x)
{
    assert(!x.empty());

    //
    std::vector<double> y = x;
    std::sort(y.begin(), y.end());
    size_t n = y.size();
    return n % 2 == 1 ? y[n / 2] : (y[n / 2 - 1] + y[n / 2]) / 2.0f;
}


RateMetrics SquiggleRead::calculate_rate_metrics(size_t strand_idx) const
{
    RateMetrics rate_metrics;

    // get kmer stats
    size_t k = this->get_base_model(strand_idx)->k;
    size_t num_kmers = this->read_sequence.length() - k + 1;
    size_t num_skips = 0;

    // collect durations, collapsing by k-mer:
    std::vector<double> durations_per_kmer(num_kmers);
    std::vector<double> events_per_kmer(num_kmers);
    for (size_t i = 0; i < this->base_to_event_map.size(); ++i) {
        size_t start_idx = this->base_to_event_map[i].indices[strand_idx].start;
        size_t end_idx = this->base_to_event_map[i].indices[strand_idx].stop;
        // no events for this k-mer
        if (start_idx == -1) {
            num_skips += 1;
            continue;
        }
        assert(start_idx <= end_idx);
        for (size_t j = start_idx; j <= end_idx; ++j) {
            durations_per_kmer[i] += this->get_duration(j, strand_idx);
            events_per_kmer[i] += 1;
        }
    }

    //
    std::sort(durations_per_kmer.begin(), durations_per_kmer.end());
    assert(durations_per_kmer.size() > 0);
    rate_metrics.median_duration = median(durations_per_kmer);

    double stall_threshold = rate_metrics.median_duration * 10;

    double num_any_event = 0;
    double num_extra_event = 0;
    size_t num_stalls = 0;
    double sum_duration;
    double sum_events;
    for(size_t i = 0; i < num_kmers; ++i) {
        sum_duration += durations_per_kmer[i];
        sum_events += events_per_kmer[i];
        num_stalls += durations_per_kmer[i] > stall_threshold;
        num_extra_event += events_per_kmer[i] > 1 ? events_per_kmer[i] - 1 : 0;
        num_any_event += events_per_kmer[i] > 0;
    }
    double mean_duration = sum_duration / num_kmers;
    double mean_events = sum_events / num_kmers;

    // Calculate median duration over segments
    size_t segment_length = 100;
    std::vector<double> segment_mean_duration;
    std::vector<double> segment_mean_events;
    for(size_t i = 0; i < num_kmers; i += segment_length) {
        double segment_sum_duration = 0.0f;
        double segment_sum_events = 0.0f;

        for(size_t j = i; j < num_kmers && j < i + segment_length; ++j) {
            segment_sum_duration += durations_per_kmer[j];
            segment_sum_events += events_per_kmer[j];
        }

        segment_mean_duration.push_back(segment_sum_duration / segment_length);
        segment_mean_events.push_back(segment_sum_events / segment_length);
    }
    double segment_median_duration = median(segment_mean_duration);
    double segment_median_events = median(segment_mean_events);

    // this is our estimator of read rate, currently we use the median duration
    // per k-mer as its more robust to outliers caused by stalls
    double read_rate = 1.0 / rate_metrics.median_duration;
    rate_metrics.skip_frequency = (double)num_skips / num_kmers;
    rate_metrics.stall_frequency = (double)num_stalls / num_kmers;
    //this->extra_event_frequency = (double)num_extra_event / num_any_event;
    rate_metrics.extra_event_frequency = 1 - (1 / (num_extra_event / num_any_event + 1));
    rate_metrics.mean_speed = 1.0 / mean_duration;

    /*
       fprintf(stderr, "[readrate] %s med: [%.5lfs %.1lf] mean: [%.5lfs %.1lf] seg: [%.5lfs %.1f] "
       "ev_mean: %.2f seg med ev: %.2lf, stalls: %zu st f: %.2f skips: %zu skip f: %.4lf ee: %zu ee f: %.4f\n",
       this->read_name.substr(0, 6).c_str(), this->rate_metrics.median_duration, read_rate, mean_duration, 1.0 / mean_duration, segment_median_duration, 1.0 / segment_median_duration,
       mean_events, segment_median_events, num_stalls, (double)num_stalls / num_kmers, num_skips, this->rate_metrics.skip_frequency, num_extra_event, this->rate_metrics.extra_event_frequency);
    */

    /*
    fprintf(stderr, "[read-metadata]\treadid=%s\tmean_read_rate=%.1lf\tmean_events_per_base=%.2lf\tshift=%.2f\tscale=%.2lf\tvar=%.2lf\n",
            this->read_name.c_str(), 1.0 / mean_duration, mean_events, this->scalings[0].shift, this->scalings[0].scale, this->scalings[0].var);
    */
    return rate_metrics;
}

       
SNRMetrics SquiggleRead::calculate_snr_metrics(size_t strand_idx) const
{
    SNRMetrics out = { 0.0f, 0.0f };
    if(this->events[strand_idx].size() < 100) {
        return out;
    }

    // SNR estimates (from CW ONT)
    std::vector<double> event_means;
    std::vector<double> event_sds;

    for(size_t i = 0; i < this->events[strand_idx].size(); ++i) {
        event_means.push_back(this->events[strand_idx][i].mean);
        event_sds.push_back(this->events[strand_idx][i].stdv);
    }

    std::sort(event_means.begin(), event_means.end());
    std::sort(event_sds.begin(), event_sds.end());

    int idx10p = event_means.size() * 0.1;
    int idx90p = event_means.size() * 0.9;
    out.current_range = event_means[idx90p] - event_means[idx10p];
    
    int idx50p = event_sds.size() * 0.5;
    out.median_sd = event_sds[idx50p];
    return out;
}
