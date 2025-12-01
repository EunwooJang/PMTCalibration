void PMTCalib_Draw_Max_Amplitude_Histogram(const char* file = "dir_root/test.root",
                 int channel = 0,
                 double pedestal_start_ns = 0,
                 double pedestal_end_ns = 500,
                 double signal_peak_search_start_ns = 580,
                 double signal_peak_search_end_ns = 660){

    PMTCalib_Delete_Canvas();

    int pedestal_x_min = std::ceil(pedestal_start_ns / time_scale_ns);
    int pedestal_x_max = std::ceil(pedestal_end_ns / time_scale_ns);
    int signal_peak_search_x_min = std::ceil(signal_peak_search_start_ns / time_scale_ns);
    int signal_peak_search_x_max = std::ceil(signal_peak_search_end_ns / time_scale_ns);


    TFile *f = TFile::Open(file);
    TTree *tree = (TTree*)f->Get("AbsEvent");

    FChannelData *data = nullptr;
    tree->SetBranchAddress("FChannelData", &data);

    Long64_t nentries = tree->GetEntries();
    
    std::vector<double> amplitude_list;
    
    for (Long64_t i = 0; i < nentries; ++i) {

        tree->GetEntry(i);

        FChannel *ch = data->Get(channel);
        if (!ch) continue;

        int ndp = ch->GetNdp();
        const unsigned short *waveform = ch->GetWaveform();
        unsigned short pedestal = ch->GetPedestal();

        if (ndp <= 0 || !waveform) continue;

        float pedestal_avg_adc = 0;
        int pedestal_count = 0;

        if (pedestal_x_max > -1) {
            for (int j = pedestal_x_min; j <= pedestal_x_max && j < ndp; ++j) {
                pedestal_avg_adc += (float)(waveform[j]) - pedestal;
                pedestal_count++;
            }

            pedestal_avg_adc /= pedestal_count;
        }

        int peak_pos_x = signal_peak_search_x_min;
        float max_height = -1e9;

        for (int j = signal_peak_search_x_min; j <= signal_peak_search_x_max && j < ndp; ++j) {
            float height = (float)(waveform[j]) - (float)pedestal - pedestal_avg_adc;
            if (height > max_height) {
                max_height = height;
                peak_pos_x = j;
            }
        }

        amplitude_list.push_back(max_height);
    }

    gStyle->SetOptStat(10);
    gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);

    TCanvas *can = new TCanvas("can",
                               Form("Channel %d - Max Amplitude Histogram", channel),
                               1200, 600);

    int min_amplitude = std::floor(*std::min_element(amplitude_list.begin(), amplitude_list.end())) - 1;
    int max_amplitude = std::ceil(*std::max_element(amplitude_list.begin(), amplitude_list.end())) + 1;
    int range = max_amplitude - min_amplitude;

    TH1D *h = new TH1D("h",
                              Form("Channel %d - Max Amplitude Histogram ;Amplitude (mV);Counts", channel),
                              range, min_amplitude * amplitude_scale_mV, max_amplitude * amplitude_scale_mV);

    
    for (double amplitude : amplitude_list) {
        h->Fill(amplitude * amplitude_scale_mV);
    }

    h->SetLineWidth(2);
    h->SetLineColor(kBlue);
    h->Draw();

    can->SetLogy();
    can->Update();
}
void PMTCalib_Draw_TTS2(const char* file = "dir_root/test.root",
                 int channel = 0,
                 double pedestal_t_ns = 500,
                 double signal_peak_search_start_ns = 580,
                 double signal_peak_search_end_ns = 660,
                 double threshold_mV = 4,
                 double fraction = 0.5){
    PMTCalib_Delete_Canvas();
    int pedestal_x_max = std::ceil(pedestal_t_ns / time_scale_ns);
    int signal_peak_search_x_min = std::ceil(signal_peak_search_start_ns / time_scale_ns);
    int signal_peak_search_x_max = std::ceil(signal_peak_search_end_ns / time_scale_ns);
    TFile *f = TFile::Open(file);
    TTree *tree = (TTree*)f->Get("AbsEvent");
    FChannelData *data = nullptr;
    tree->SetBranchAddress("FChannelData", &data);
    Long64_t nentries = tree->GetEntries();
    
    std::vector<double> tts_list;
    
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);
        FChannel *ch = data->Get(channel);
        if (!ch) continue;
        int ndp = ch->GetNdp();
        const unsigned short *waveform = ch->GetWaveform();
        unsigned short pedestal = ch->GetPedestal();
        if (ndp <= 0 || !waveform) continue;
        float pedestal_avg_adc = 0;
        int pedestal_count = 0;
        if (pedestal_x_max > -1) {
            for (int j = 0; j <= pedestal_x_max && j < ndp; ++j) {
                pedestal_avg_adc += (float)(waveform[j]) - pedestal;
                pedestal_count++;
            }
            pedestal_avg_adc /= pedestal_count;
        }
        int peak_pos_x = signal_peak_search_x_min;
        float max_height = -1e9;
        for (int j = signal_peak_search_x_min; j <= signal_peak_search_x_max && j < ndp; ++j) {
            float height = (float)(waveform[j]) - (float)pedestal - pedestal_avg_adc;
            if (height > max_height) {
                max_height = height;
                peak_pos_x = j;
            }
        }
        
        double max_amplitude_mV = max_height * amplitude_scale_mV;
        if (max_amplitude_mV < threshold_mV) continue;
        
        double threshold_value = max_height * fraction;
        int timing_x = peak_pos_x;
        for (int j = peak_pos_x; j >= 0; --j) {
            float height = (float)(waveform[j]) - (float)pedestal - pedestal_avg_adc;
            if (height < threshold_value) {
                timing_x = j + 1;
                break;
            }
            timing_x = j;
        }
        
        double timing_ns = timing_x * time_scale_ns;
        tts_list.push_back(timing_ns);
    }
    
    gStyle->SetOptStat(10);
    gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);
    TCanvas *can = new TCanvas("can",
                               Form("Channel %d - TTS Histogram", channel),
                               1200, 600);
    double min_time = *std::min_element(tts_list.begin(), tts_list.end()) - 1;
    double max_time = *std::max_element(tts_list.begin(), tts_list.end()) + 1;
    int nbins = std::ceil((max_time - min_time) / 2.0);
    TH1D *h = new TH1D("h",
                       Form("Channel %d - TTS Histogram (%.0f%% of max, threshold > %.1f mV);Time (ns);Counts", 
                            channel, fraction * 100, threshold_mV),
                       nbins, min_time, max_time);
    
    for (double timing : tts_list) {
        h->Fill(timing);
    }
    h->SetLineWidth(2);
    h->SetLineColor(kBlue);
    h->Draw();
    can->SetLogy();
    can->Update();
}

void PMTCalib_Draw_TTS(const char* file = "dir_root/test.root",
                 int channel = 0,
                 double pedestal_t_ns = 500,
                 double signal_peak_search_start_ns = 580,
                 double signal_peak_search_end_ns = 660,
                 double threshold_mV = 4,
                 double fraction = 0.5){
    PMTCalib_Delete_Canvas();
    int pedestal_x_max = std::ceil(pedestal_t_ns / time_scale_ns);
    int signal_peak_search_x_min = std::ceil(signal_peak_search_start_ns / time_scale_ns);
    int signal_peak_search_x_max = std::ceil(signal_peak_search_end_ns / time_scale_ns);
    TFile *f = TFile::Open(file);
    TTree *tree = (TTree*)f->Get("AbsEvent");
    FChannelData *data = nullptr;
    tree->SetBranchAddress("FChannelData", &data);
    Long64_t nentries = tree->GetEntries();
    
    std::vector<double> tts_list;
    
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);
        FChannel *ch = data->Get(channel);
        if (!ch) continue;
        int ndp = ch->GetNdp();
        const unsigned short *waveform = ch->GetWaveform();
        unsigned short pedestal = ch->GetPedestal();
        if (ndp <= 0 || !waveform) continue;
        float pedestal_avg_adc = 0;
        int pedestal_count = 0;
        if (pedestal_x_max > -1) {
            for (int j = 0; j <= pedestal_x_max && j < ndp; ++j) {
                pedestal_avg_adc += (float)(waveform[j]) - pedestal;
                pedestal_count++;
            }
            pedestal_avg_adc /= pedestal_count;
        }
        int peak_pos_x = signal_peak_search_x_min;
        float max_height = -1e9;
        for (int j = signal_peak_search_x_min; j <= signal_peak_search_x_max && j < ndp; ++j) {
            float height = (float)(waveform[j]) - (float)pedestal - pedestal_avg_adc;
            if (height > max_height) {
                max_height = height;
                peak_pos_x = j;
            }
        }
        
        double max_amplitude_mV = max_height * amplitude_scale_mV;
        if (max_amplitude_mV < threshold_mV) continue;
        
        double threshold_value = max_height * fraction;
        int timing_x = peak_pos_x;
        for (int j = peak_pos_x; j >= 0; --j) {
            float height = (float)(waveform[j]) - (float)pedestal - pedestal_avg_adc;
            if (height < threshold_value) {
                timing_x = j + 1;
                break;
            }
            timing_x = j;
        }
        
        double timing_ns = timing_x * time_scale_ns;
        tts_list.push_back(timing_ns);
    }
    
    gStyle->SetOptStat(10);
    gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);
    TCanvas *can = new TCanvas("can",
                               Form("Channel %d - TTS Histogram", channel),
                               1200, 600);
    double min_time = *std::min_element(tts_list.begin(), tts_list.end()) - 1;
    double max_time = *std::max_element(tts_list.begin(), tts_list.end()) + 1;
    int nbins = std::ceil(max_time - min_time);
    TH1D *h = new TH1D("h",
                       Form("Channel %d - TTS Histogram (%.0f%% of max, threshold > %.1f mV);Time (ns);Counts", 
                            channel, fraction * 100, threshold_mV),
                       nbins, min_time, max_time);
    
    for (double timing : tts_list) {
        h->Fill(timing);
    }
    h->SetLineWidth(2);
    h->SetLineColor(kBlue);
    h->Draw();
    can->SetLogy();
    can->Update();
}
