

void PMTCalib_Print_Single_Event_Info(const char* file, int event = 0) {
  
  TChain *t = new TChain("AbsEvent");
  t->Add(file);
  
  EventInfo *info = new EventInfo();
  FChannelData *fdata = new FChannelData();
  
  t->SetBranchAddress("EventInfo", &info);
  t->SetBranchAddress("FChannelData", &fdata);
  
  Long64_t nentries = t->GetEntries();
 
  t->GetEntry(event);
  
  std::cout << std::endl;
  std::cout << "======================================" << std::endl;
  std::cout << "      EVENT INFORMATION" << std::endl;
  std::cout << "======================================" << std::endl;
  std::cout << "Event Number:     " << info->GetEventNumber() << std::endl;
  std::cout << "Trigger Number:   " << info->GetTriggerNumber() << std::endl;
  std::cout << "Trigger Type:     " << info->GetTriggerType() << std::endl;
  std::cout << "Number of Hits:   " << info->GetNHit() << std::endl;
  std::cout << "Trigger Time:     " << info->GetTriggerTime() << std::endl;
  std::cout << "======================================" << std::endl;
  std::cout << std::endl;
  
  std::cout << "--- FADC CHANNEL DATA (Waveform) ---" << std::endl;
  int nfch = fdata->GetN();
  
  if (nfch == 0) {
    std::cout << "  No FADC channels found" << std::endl;
  } else {
    std::cout << "  Total FADC Channels: " << nfch << std::endl;
    std::cout << std::endl;
    
    for (int i = 0; i < nfch; i++) {
      FChannel *ch = fdata->Get(i);
      
      std::cout << "  Channel " << i << ":" << std::endl;
      std::cout << "    ID:           " << ch->GetID() << std::endl;
      std::cout << "    Trigger Bit:  " << ch->GetBit() << std::endl;
      std::cout << "    Pedestal:     " << ch->GetPedestal() << std::endl;
      std::cout << "    Data Points:  " << ch->GetNdp() << std::endl;
      
      const unsigned short *wave = ch->GetWaveform();
      double pedestal = ch->GetPedestal();
      if (wave && ch->GetNdp() > 0) {
        
          double sum = 0;
          for (int j = 0; j < ch->GetNdp(); j++) {
            sum += wave[j] - pedestal;
          }
          double mean = sum / ch->GetNdp();
          
          double sum_sq = 0;
          for (int j = 0; j < ch->GetNdp(); j++) {
            sum_sq += (wave[j] - pedestal - mean) * (wave[j] - pedestal - mean);
          }
          double std = std::sqrt(sum_sq / ch->GetNdp()); 
        
        std::cout << "    Waveform Statistics:" << std::endl;
        std::cout << "      Mean:  " << mean << std::endl;
        std::cout << "      Std:   " << std << std::endl;
        std::cout << "      First 10 samples: ";
        for (int j = 0; j < TMath::Min(10, ch->GetNdp()); j++) {
          std::cout << wave[j] - pedestal<< " ";
        }
        if (ch->GetNdp() > 10) std::cout << "...";
        std::cout << std::endl;
      }   
      std::cout << std::endl;
    }
  }
  delete t;
}



void PMTCalib_Draw_Amplitude_Histogram(const char* file, int channel = 0) {
  
  PMTCalib_Delete_Canvas();
  
  TChain *t = new TChain("AbsEvent");
  t->Add(file);
  
  FChannelData *data = nullptr;
  t->SetBranchAddress("FChannelData", &data);
  
  Long64_t nevt = t->GetEntries();
  
  t->GetEntry(0);
  FChannel *ch = data->Get(0);
  int ndp = ch->GetNdp();

  std::vector<int> adc_vals;
  adc_vals.reserve(nevt * ndp);
  
  int vmin = 0;
  int vmax = 0;
  
  for (Long64_t i = 0; i < nevt; ++i) {
    t->GetEntry(i);
    FChannel *ch = data->Get(channel);
    
    if (!ch) continue;
    
    const unsigned short *waveform = ch->GetWaveform();
    unsigned short pedestal = ch->GetPedestal();
    
    for (int j = 0; j < ndp; ++j) {
      int v = (int)waveform[j] - (int)pedestal;
      adc_vals.push_back(v);
      
        if (v < vmin) vmin = v;
        if (v > vmax) vmax = v;
    
    }
  }
  
  int bin_min = (int)std::floor(vmin - 1.);
  int bin_max = (int)std::ceil(vmax + 1.);
  int nbins = bin_max - bin_min + 1;
  
  TH1D *h = new TH1D("h",
                     Form("Channel %d - Amplitude Histogram;Voltage (mV);Counts", channel),
                     nbins,
                     bin_min * amplitude_scale_mV, 
                     (bin_max + 1) * amplitude_scale_mV);
  
  for (const int &v : adc_vals) {
    h->Fill(v * amplitude_scale_mV);
  }
  
  TCanvas *can = new TCanvas("can",
                             Form("Channel %d - Amplitude Histogram", channel),
                             1200, 
                             600);
 
  gStyle->SetOptStat(1110);
  gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);

  h->SetLineWidth(2);
  h->Draw();
  
  can->SetLogy();
  can->Modified();
  can->Update();
}



void PMTCalib_Draw_Pedestal_Histogram(const char* file, int channel = 0) {

  PMTCalib_Delete_Canvas();

  TChain *t = new TChain("AbsEvent");
  t->Add(file);

  FChannelData *data = nullptr;
  t->SetBranchAddress("FChannelData", &data);
  
  Long64_t nevt = t->GetEntries();
  
  std::vector<int> ped_vals;
  ped_vals.reserve(nevt);

  int vmin = 0;
  int vmax = 0;

  for (Long64_t i = 0; i < nevt; ++i) {
    t->GetEntry(i);

    FChannel *ch = data->Get(channel);
    if (!ch) continue;

    unsigned short pedestal = ch->GetPedestal();
    int v = (int)pedestal;

    ped_vals.push_back(v);

    if (v < vmin) vmin = v;
    if (v > vmax) vmax = v;
  }

  int bin_min = (int)std::floor(vmin - 1.);
  int bin_max = (int)std::ceil(vmax + 1.);
  int nbins = bin_max - bin_min + 1;

  TH1D *h = new TH1D("h",
                     Form("Channel %d - Pedestal Histogram;Voltage (mV);Counts", channel),
                     nbins,
                         bin_min * amplitude_scale_mV, (bin_max + 1) * amplitude_scale_mV);

  for (const int &v : ped_vals) {
    h->Fill(v * amplitude_scale_mV);
  }


  TCanvas *can = new TCanvas("can",
                             Form("Channel %d - Pedestal Histogram", channel),
                             1200, 
                             600);
  
  gStyle->SetOptStat(1110);
  gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);

  h->SetLineWidth(2);
  h->Draw();
  
  can->SetLogy();
  can->Modified();
  can->Update();
}



void PMTCalib_Draw_Single_Event(const char* file, int channel = 0, int event = 0) {

  PMTCalib_Delete_Canvas();

  TChain *t = new TChain("AbsEvent");
  t->Add(file);

  FChannelData *data = new FChannelData();
  t->SetBranchAddress("FChannelData", &data);

  t->GetEntry(event);
  FChannel *ch = data->Get(channel);

  int ndp = ch->GetNdp();
  const unsigned short *waveform = ch->GetWaveform();
  unsigned short pedestal = ch->GetPedestal();

  TH1D *h = new TH1D("h",
                     Form("Channel %d - Event %d;Time (ns);Voltage (mV)", channel, event),
                     ndp, 0, ndp * time_scale_ns);

  for (int j = 0; j < ndp; ++j) {
    int adc_value = (int)(waveform[j]) - (int)pedestal;
    h->SetBinContent(j + 1, adc_value * amplitude_scale_mV);
  }

  TCanvas *can = new TCanvas("can", "Event Viewer", 1200, 600);

  gStyle->SetOptStat(0);
  gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);

  h->SetLineWidth(2);
  h->Draw("HIST");

  can->Modified();
  can->Update();
}



void PMTCalib_Draw_Single_Event_Condition(const char* file,
                                          int channel = 0,
                                          int event_num = -1,
                                          int isAbove = 0,
                                          double threshold_mv = 2.) {

  PMTCalib_Delete_Canvas();
  
  TChain *t = new TChain("AbsEvent");
  t->Add(file);

  FChannelData *data = new FChannelData();
  t->SetBranchAddress("FChannelData", &data);

  Long64_t nentries = t->GetEntries();

  TH1D *h = nullptr;
  int found_event = -1;

  Long64_t start_event = (event_num >= 0) ? event_num : 0;

  for (Long64_t i = start_event; i < nentries; ++i) {

    t->GetEntry(i);
    FChannel *ch = data->Get(channel);

    int ndp = ch->GetNdp();
    const unsigned short *waveform = ch->GetWaveform();
    unsigned short pedestal = ch->GetPedestal();

    double max_val = -1e9;
    double min_val = 1e9;

    for (int j = 0; j < ndp; ++j) {
      int adc_value = (int)(waveform[j]) - (int)pedestal;
      if (adc_value > max_val) max_val = adc_value;
      if (adc_value < min_val) min_val = adc_value;
    }

    bool meet_threshold = false;
    if (isAbove) {
        meet_threshold = (max_val > threshold_mv / amplitude_scale_mV);
    } else {
        meet_threshold = (min_val < threshold_mv / amplitude_scale_mV);
    }

    if (meet_threshold) {
       h = new TH1D("h",
                    Form("Channel %d - Event %lld;Time (ns);Voltage (mV)", channel, i),
                    ndp, 0, ndp * time_scale_ns);

        for (int j = 0; j < ndp; ++j) {
          int adc_value = (int)(waveform[j]) - (int)pedestal;
          h->SetBinContent(j + 1, adc_value * amplitude_scale_mV);
        }

        found_event = i;
        break;
    }
  }
  
  TCanvas *can = new TCanvas("can", Form("Channel %d Event %d", channel, found_event), 1200, 600);

  gStyle->SetOptStat(0);
  gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);

  h->SetLineWidth(2);
  h->Draw("HIST");

  can->Modified();
  can->Update();
}



void PMTCalib_Draw_Events_Persistence(const char* file,
                                         int channel = 0,
                                         double min_mV = -10.,
                                         double max_mV = 60.) {


  PMTCalib_Delete_Canvas();

  TChain *t = new TChain("AbsEvent");
  t->Add(file);

  FChannelData *data = new FChannelData();
  t->SetBranchAddress("FChannelData", &data);

  Long64_t nevt = t->GetEntries();

  t->GetEntry(0);
  FChannel *ch = data->Get(0);
  int ndp = ch->GetNdp();
  
  double max_adc = std::floor(max_mV / amplitude_scale_mV);
  double min_adc = std::ceil(min_mV / amplitude_scale_mV);
  int nbins_adc = (int)(max_adc - min_adc) + 1;

  std::cout << "ADC range: [" << min_adc << ", " << max_adc << "]" << std::endl;
  std::cout << "Time bins: " << ndp << std::endl;

  TH2D *h = new TH2D("h",
                             Form("Channel %d Persistence;Time (ns);Voltage (mV)", channel),
                             ndp, 0, ndp * time_scale_ns,
                             nbins_adc, min_adc * amplitude_scale_mV, (max_adc + 1) * amplitude_scale_mV);

  int filled_count = 0;

  for (Long64_t i = 0; i < nevt; ++i) {
    t->GetEntry(i);

    FChannel *ch = data->Get(channel);

    const unsigned short *waveform = ch->GetWaveform();
    unsigned short pedestal = ch->GetPedestal();

    for (int j = 0; j < ndp; ++j) {
      int time_sample = j;
      int adc_value = (int)(waveform[j]) - (int)pedestal;
      h->Fill(time_sample * time_scale_ns, adc_value * amplitude_scale_mV);
    }

    filled_count++;

    if ((i + 1) % 10000 == 0) {
      std::cout << "  Processed " << (i + 1) << " / " << nevt << " events..." << std::endl;
    }
  }

  h->SetEntries(nevt);
  
  TCanvas *can = new TCanvas("can", Form("Channel %d - Persistence", channel), 1200, 600);
  
  gStyle->SetOptStat(10);
  gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);
  gStyle->SetPalette(kBird);

  h->Draw("COLZ");
  
  can->Modified();
  can->Update();
}

