R__LOAD_LIBRARY(build/libRawObjs.dylib)


TChain *g_t;
FChannelData *g_data;
Long64_t g_nentries;
int g_current_event;
int g_channel;
TCanvas *g_can;

void DrawCurrentEvent() {

    delete gROOT->FindObject("h");

    g_t->GetEntry(g_current_event);
    FChannel *ch = g_data->Get(g_channel);

    int ndp = ch->GetNdp();
    const unsigned short *waveform = ch->GetWaveform();
    unsigned short pedestal = ch->GetPedestal();

    TH1D *h = new TH1D("h",
                   Form("Channel %d - Event %d;Time (ns);Voltage (mV)", g_channel, g_current_event),
                    ndp, 0, ndp * time_scale_ns);

    for (int j = 0; j < ndp; ++j) {
        int adc_value = (int)(waveform[j]) - (int)pedestal;
        h->SetBinContent(j + 1, adc_value * amplitude_scale_mV);
    }

    h->SetLineWidth(2);
    
    g_can->cd();
    h->Draw("HIST");

    h->GetYaxis()->SetRangeUser(-10, 60);

    g_can->Modified();
    g_can->Update();
}

void KeyPressHandler() {
    int key = gPad->GetEvent();
    if (key == kKeyPress) {
        int keycode = gPad->GetEventX();
        if (keycode == 'a' || keycode == 'A') {
            g_current_event++;
            if (g_current_event >= g_nentries) {
                g_current_event = g_nentries - 1;
            }
            DrawCurrentEvent();
        }
    }
}

void PMTCalib_Draw_Events(const char* file, int channel = 0, int event = 0) {
    

    delete gROOT->FindObject("h");
    delete gROOT->FindObject("can");

    g_t = new TChain("AbsEvent");
    g_t->Add(file);
    g_data = new FChannelData();
    g_t->SetBranchAddress("FChannelData", &g_data);
    g_nentries = g_t->GetEntries();
    g_channel = channel;
    g_current_event = event;
    gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);

    g_can = new TCanvas("can", "Event Viewer", 1200, 600);
    
    DrawCurrentEvent();
    
    g_can->AddExec("ex", "KeyPressHandler()");
}
