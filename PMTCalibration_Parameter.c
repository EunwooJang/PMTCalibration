// Voltage - 2.5V pp & 12 bits resolution
double Vpp_mV = 2500.;
double bit_resolution = 12.;
double amplitude_scale_mV = Vpp_mV / pow(2, bit_resolution);
double resistence_ohm = 50.;

// Time - 500 Mhz -> 1 sample = 2ns 
double Sample_Rate_MHz = 500;
double time_scale_ns = 2.;

double DarkRate_SampleN = 16368;
double SPE_SampleN = 512; // FIXME
std::string Root_Directory = "dir_root"; // FIXME
