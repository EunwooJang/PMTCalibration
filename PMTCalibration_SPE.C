


void PMTCalib_Calculate_SPE_Voltage_Height(int gain = 1e+7, 
                   int resistence_ohm = 50,
                   float pulse_width_ns = 20
                   )
{
    float e_charge = 1.602E-19;  // C
    float spe_charge = e_charge * gain; // C
    float mV_ns_area = spe_charge * resistence_ohm * 1E+12; // mV*ns (Area) or pC * Ohm
    float Pulse_Height = 2 * mV_ns_area / pulse_width_ns; // mV -> Area Shape as Triangle = Pulse Heigth * Pulse Width / 2

    cout << "Pulse Height (mV): " << Pulse_Height << endl;
}



void PMTCalib_Draw_Charge_Histogram(const char* file = "dir_root/test.root",
                 int channel = 0,
                 int pedestal_t_ns = 500.,
                 double signal_peak_search_start_ns = 580.,
                 double signal_peak_search_end_ns = 660.,
                 double signal_integral_start_ns = 4.,
                 double signal_integral_end_ns = 16.,
                 double threshold_mV = 8.){


    PMTCalib_Delete_Canvas();

    int pedestal_x_max = std::ceil(pedestal_t_ns / time_scale_ns);
    int signal_peak_search_x_min = std::ceil(signal_peak_search_start_ns / time_scale_ns);
    int signal_peak_search_x_max = std::ceil(signal_peak_search_end_ns / time_scale_ns);
    int signal_integral_x_min = std::ceil(signal_integral_start_ns / time_scale_ns);
    int signal_integral_x_max = std::ceil(signal_integral_end_ns / time_scale_ns);
    int threshold = std::ceil(threshold_mV / amplitude_scale_mV);


    TFile *f = TFile::Open(file);
    TTree *tree = (TTree*)f->Get("AbsEvent");

    FChannelData *data = nullptr;
    tree->SetBranchAddress("FChannelData", &data);

    Long64_t nentries = tree->GetEntries();
    
    std::vector<double> charge_list;
    int above_threshold_count = 0;
    
    double factor = amplitude_scale_mV * time_scale_ns / resistence_ohm;

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

            for (int j = 0; j < pedestal_x_max && j < ndp; ++j) {
                pedestal_avg_adc += (float)(waveform[j]) - pedestal;
                pedestal_count++;
            }

            if (pedestal_count > 0) {
                pedestal_avg_adc /= pedestal_count;
            }
        }

        int peak_pos_x = signal_peak_search_x_min;
        float max_height = -1e9;

        for (int j = signal_peak_search_x_min; j < signal_peak_search_x_max && j < ndp; ++j) {
            float height = (float)(waveform[j]) - (float)pedestal - pedestal_avg_adc;
            if (height > max_height) {
                max_height = height;
                peak_pos_x = j;
            }
        }

        if (max_height > threshold) {
            above_threshold_count++;
        }
        int integral_x_min = peak_pos_x - signal_integral_x_min;
        int integral_x_max = peak_pos_x + signal_integral_x_max;

        if (integral_x_min < 0) integral_x_min = 0;
        if (integral_x_max >= ndp) integral_x_max = ndp - 1;

        double charge = 0;
        for (int j = integral_x_min; j < integral_x_max; ++j) {
            float height = (float)(waveform[j]) - (float)pedestal - pedestal_avg_adc;
            charge += height;
        }

        charge *= factor;
        charge_list.push_back(charge);
    }


    gStyle->SetOptStat(10);

    TCanvas *can = new TCanvas("can",
                               Form("Channel %d - Charge Histogram", channel),
                               1200, 600);

    double min_charge = *std::min_element(charge_list.begin(), charge_list.end());
    double max_charge = *std::max_element(charge_list.begin(), charge_list.end());

    int min_val = (int)std::floor(min_charge / factor) - 1;
    int max_val = (int)std::ceil(max_charge / factor) + 1;
    int range = max_val - min_val;

    TH1D *h = new TH1D("h",
                          Form("Channel %d - Charge Histogram;Charge (pC);Counts", channel),
                          range, min_val * factor, max_val * factor); 
    for (double charge : charge_list) {
        h->Fill(charge);
    }

    h->SetLineWidth(2);
    h->SetLineColor(kBlue);
    h->Draw();
    can->SetLogy();
    can->Update();

    cout << "\n========================================" << endl;
    cout << "Threshold (mV): " << threshold_mV << endl;
    cout << "μ: " << (float)above_threshold_count / nentries << endl;
    cout << "========================================" << endl;

}
