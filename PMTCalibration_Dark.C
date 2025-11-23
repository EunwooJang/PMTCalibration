void PMTCalib_Create_DarkCounts_Root(const char* directory, double threshold_mV = 8.) {

    TSystemDirectory dir(directory, directory);
    TList *files = dir.GetListOfFiles();
    
    std::vector<std::string> input_files;
    std::regex pattern("RUN.*\\.[0-9]{5}$");
    
    TIter next(files);
    TSystemFile *file;
    while ((file = (TSystemFile*)next())) {
        std::string filename = file->GetName();
        if (!file->IsDirectory() && std::regex_match(filename, pattern)) {
            input_files.push_back(std::string(directory) + "/" + filename);
        }
    }
    
    std::sort(input_files.begin(), input_files.end());
    
    int threshold_adc = std::ceil(threshold_mV / amplitude_scale_mV);
    for (const auto& input_file : input_files) {
        
        size_t pos = input_file.rfind('.');
        if (pos == std::string::npos) continue;
        std::string file_id = input_file.substr(pos + 1);
        if (file_id.length() != 5) continue;
        
        std::string filename = input_file.substr(input_file.rfind('/') + 1);
        
        TChain *t = new TChain("AbsEvent");
        t->Add(input_file.c_str());
        
        FChannelData *data = new FChannelData();
        t->SetBranchAddress("FChannelData", &data);
        
        Long64_t nentries = t->GetEntries();
        
        std::cout << "Total entries: " << nentries << std::endl;
        
        std::vector<int> available_channels;
        t->GetEntry(0);
        for (int ch_id = 0; ch_id < 1000; ++ch_id) {
            FChannel *ch = data->Get(ch_id);
            if (ch && ch->GetNdp() > 0) {
                available_channels.push_back(ch_id);
            }
        }
        
        std::cout << "Available channels: ";
        for (int ch : available_channels) {
            std::cout << ch << " ";
        }
        std::cout << std::endl;
        std::cout << "Total channels: " << available_channels.size() << std::endl;
        
        std::string output_filename = std::string(directory) + "/Dark_Counts_th" + std::to_string(threshold_adc) + ".root." + file_id;
        TFile *outfile = new TFile(output_filename.c_str(), "RECREATE");
        TTree *outtree = new TTree("DarkCounts", "Dark Counts per Event");
        
        int event_id;
        std::vector<int> dark_counts(available_channels.size(), 0);
        
        outtree->Branch("event_id", &event_id, "event_id/I");
        
        for (size_t i = 0; i < available_channels.size(); ++i) {
            std::string branch_name = "ch_" + std::to_string(available_channels[i]);
            outtree->Branch(branch_name.c_str(), &dark_counts[i], (branch_name + "/I").c_str());
        }
        
        for (Long64_t i = 0; i < nentries; ++i) {
            
            t->GetEntry(i);
            event_id = i;
            
            for (size_t ch_idx = 0; ch_idx < available_channels.size(); ++ch_idx) {
                int ch_id = available_channels[ch_idx];
                FChannel *ch = data->Get(ch_id);
                
                dark_counts[ch_idx] = 0;
                
                if (!ch) continue;
                
                int ndp = ch->GetNdp();
                const unsigned short * waveform = ch->GetWaveform();
                unsigned short pedestal = ch->GetPedestal();
                
                if (ndp <= 0 || !waveform) continue;
                
                bool in_threshold_region = false;
                for (int j = 0; j < ndp; ++j) {
                    double bin_value = (double)(waveform[j]) - pedestal;
                    if (bin_value >= threshold_adc) {
                        if (!in_threshold_region) {
                            dark_counts[ch_idx]++;
                            in_threshold_region = true;
                        }
                    } else {
                        in_threshold_region = false;
                    }
                }
            }
            
            outtree->Fill();
        }
        
        outtree->Write();
        outfile->Close();
        delete outfile;
        delete t;
        
    }
}

void PMTCalib_Calculate_DarkRate(const char* directory,
                         int run_number = 0,
                         int channel = 0,
                         int start_event = 0,
                         int event_n = 100,
                         double threshold_mV = 8.)
{

    int threshold_adc = std::ceil(threshold_mV / amplitude_scale_mV);

    TString filename = Form("%s/Dark_Counts_th%d.root.%05d", directory, threshold_adc, run_number);
    
    TFile *tf = TFile::Open(filename);
    
    TTree *tree = (TTree*)tf->Get("DarkCounts");
    
    TString branch_name = Form("ch_%d", channel);
    
    int dark_count = 0;
    int event_id = 0;
    
    tree->SetBranchAddress("event_id", &event_id);
    tree->SetBranchAddress(branch_name, &dark_count);
    
    Long64_t total_entries = tree->GetEntries();
    Long64_t total_dark_count = 0;
    int count = 0;
    
    Long64_t end_event = (start_event + event_n > total_entries) ? total_entries : start_event + event_n;
    
    for (Long64_t i = start_event; i < end_event; ++i) {
        tree->GetEntry(i);
        total_dark_count += dark_count;
        count++;
    }
    
    cout << "---------------------------------" << endl;
    cout << "Run Number: " << run_number << endl;
    cout << "Channel: " << channel << endl;
    cout << "Event Range: " << start_event << " ~ " << end_event - 1 << endl;
    cout << "Threshold (mV): " << threshold_mV << endl;
    cout << "Total Dark Count: " << total_dark_count << endl;
    cout << "Dark Rate (Hz): " << total_dark_count / (2.0 * 16368 * event_n / 1e+9) << endl;
    cout << "Events Processed: " << count << endl;
    cout << "---------------------------------" << endl;
    
    tf->Close();
    delete tf;
}



void PMTCalib_Draw_DarkRate(const char* directory,
                    vector<int> run = {0,7,14},
                    int channel = 0,
                    vector<double> threshold_mV = {8., 4.},
                    int event_interval = 6000,
                    const char* name = "PMT"){


    PMTCalib_Delete_Canvas();

    TString filename = Form("%s/%s_Dark_Rate.png", directory, name);

    TCanvas *can = new TCanvas("can", "Dark Rate vs Time", 1200, 600);
    TMultiGraph *mg = new TMultiGraph();

    vector<int> colors = {kRed, kBlue, kMagenta, kCyan};

    double max_hour = 0;

    for (size_t t_idx = 0; t_idx < threshold_mV.size(); t_idx++) {
        vector<vector<pair<double, double>>> segments;
        vector<pair<double, double>> current_segment;

        int remainder_dark_count = 0;
        int remainder_count = 0;
        double current_hour = 0;
        int prev_file = -999;

        for (size_t run_idx = 0; run_idx < run.size(); run_idx++) {
            int current_file = run[run_idx];

            if (prev_file != -999 && current_file != prev_file + 1) {
                if (!current_segment.empty()) {
                    segments.push_back(current_segment);
                    current_segment.clear();
                }
                remainder_dark_count = 0;
                remainder_count = 0;
            }
            prev_file = current_file;

            int threshold_adc = std::ceil(threshold_mV[t_idx] / amplitude_scale_mV);

            TString filename = Form("%s/Dark_Counts_th%d.root.%05d", directory, threshold_adc, current_file);

            TFile *tf = TFile::Open(filename);

            TTree *tree = (TTree*)tf->Get("DarkCounts");

            TString branch_name = Form("ch_%d", channel);

            int dark_count = 0;
            int event_id = 0;

            tree->SetBranchAddress("event_id", &event_id);
            tree->SetBranchAddress(branch_name, &dark_count);

            Long64_t file_entries = tree->GetEntries();

            for (Long64_t i = 0; i < file_entries; i++) {
                tree->GetEntry(i);

                remainder_dark_count += dark_count;
                remainder_count++;

                if (remainder_count >= event_interval) {
                    double dark_rate = (double)remainder_dark_count / (time_scale_ns * DarkRate_SampleN * event_interval / 1e+9);
                    current_hour = current_file + (double)i / file_entries;

                    current_segment.push_back({current_hour, dark_rate});

                    remainder_dark_count = 0;
                    remainder_count = 0;
                }
            }

            current_hour = current_file + 1.0;
            tf->Close();
            delete tf;
        }

        if (remainder_count > 0) {
            double dark_rate = (double)remainder_dark_count / (time_scale_ns * DarkRate_SampleN * remainder_count / 1e+9);
            current_segment.push_back({current_hour, dark_rate});
        }

        if (!current_segment.empty()) {
            segments.push_back(current_segment);
        }

        max_hour = current_hour;

        for (auto& data_points : segments) {
            TGraph *gr = new TGraph();
            int idx = 0;

            for (size_t i = 0; i < data_points.size(); i++) {
                double x = data_points[i].first;
                double y = data_points[i].second;

                if (i > 0) {
                    double prev_x = data_points[i-1].first;
                    double prev_y = data_points[i-1].second;
                    gr->SetPoint(idx++, x, prev_y);
                }
                gr->SetPoint(idx++, x, y);
            }

            gr->SetLineColor(colors[t_idx % colors.size()]);
            gr->SetLineWidth(2);
            mg->Add(gr, "L");
        }
    }

    mg->Draw("A");

    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    for (size_t t_idx = 0; t_idx < threshold_mV.size(); t_idx++) {
        TGraph *dummy = new TGraph();
        dummy->SetLineColor(colors[t_idx % colors.size()]);
        dummy->SetLineWidth(2);
        leg->AddEntry(dummy, Form("Threshold (mV) = %.2f", threshold_mV[t_idx]), "L");
    }
    leg->Draw();

    int max_file = run[run.size()-1];
    mg->GetXaxis()->SetLimits(0, max_file + 1);
    mg->GetYaxis()->SetRangeUser(0, 50000);

    mg->GetXaxis()->SetTitle("Run Time (Hour)");
    mg->GetYaxis()->SetTitle("Dark Rate (Hz)");

    can->SaveAs(filename);
    can->Draw();

}
                    
