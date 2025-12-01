void PMTCalib_SPE_Fit_1(TH1D *h, TCanvas *can,
                    double fit_min = -2.0,
                    double fit_max = 12.0,

                    double mu = 0.061,

                    double w = 0.0175,
                    double Q_0 = 0.035,
                    double std_0 = 0.101,
                    double alpha = 0.65,
                    
                    int npe_max = 3,
                    double Q_1 = 2.92,
                    double std_1 = 0.72) {
    
    can->cd();
    
    // Helper function for error function
    auto erf_helper = [](double x) {
        return TMath::Erf(x);
    };
    
    // Helper function for sign
    auto sign = [](double x) {
        return (x >= 0) ? 1.0 : -1.0;
    };
    
    // Create TF1 for Pedestal component (n=0)
    // For n=0: (1-w)*G_0(x-Q_0) + w*I_{G_0⊗E}(x-Q_0)
    TF1 *f_pedestal = new TF1("f_pedestal", [](double *x, double *par) {
        double mu = par[0];
        double w = par[1];
        double Q_0 = par[2];
        double std_0 = par[3];
        double alpha = par[4];
        
        // Term 1: (1-w) * G_0(x-Q_0)
        double norm = 1.0 / (std_0 * TMath::Sqrt(2.0 * TMath::Pi()));
        double exponent = -TMath::Power(x[0] - Q_0, 2) / (2.0 * std_0 * std_0);
        double term1 = (1.0 - w) * norm * TMath::Exp(exponent);
        
        // Term 2: w * I_{G_0⊗E}(x-Q_0)
        // For n=0: I_{G_0⊗E} = α * exp[-α(x-Q_0)]
        double term2 = 0.0;
        if (x[0] >= Q_0) {
            term2 = w * alpha * TMath::Exp(-alpha * (x[0] - Q_0));
        }
        
        return TMath::Exp(-mu) * (term1 + term2);
    }, fit_min, fit_max, 5);
    f_pedestal->SetParameters(mu, w, Q_0, std_0, alpha);
    f_pedestal->SetNpx(1000);
    f_pedestal->SetLineColor(kBlack);
    f_pedestal->SetLineStyle(2);
    f_pedestal->SetLineWidth(2);
    
    // Create TF1 for I_{G_n⊗E}(x-Q_0) component for n >= 1
    // This implements Equation (8)
    auto IGnE_func = [](double x, double Q_0, int n, double Q_1, double std_1, double alpha) {
        if (n <= 0) return 0.0;
        
        double Q_n = Q_0 + n * Q_1;
        double sigma_n = std_1 * TMath::Sqrt(n);
        double sigma_n_sq = sigma_n * sigma_n;
        
        // I_{G_n⊗E}(x-Q_0) = (α/2) * exp[-α(x - Q_0 - α*σ_n²)]
        //                    × {erf[(Q_0 - Q_n - σ_n²*α)/(σ_n*√2)]
        //                       + sign(x - Q_n - σ_n²*α) * erf[|x - Q_n - σ_n²*α|/(σ_n*√2)]}
        
        double exp_term = (alpha / 2.0) * TMath::Exp(-alpha * (x - Q_0 - alpha * sigma_n_sq));
        
        double erf1_arg = (Q_0 - Q_n - sigma_n_sq * alpha) / (sigma_n * TMath::Sqrt(2.0));
        double erf1 = TMath::Erf(erf1_arg);
        
        double erf2_arg_raw = x - Q_n - sigma_n_sq * alpha;
        double erf2_sign = (erf2_arg_raw >= 0) ? 1.0 : -1.0;
        double erf2_arg = TMath::Abs(erf2_arg_raw) / (sigma_n * TMath::Sqrt(2.0));
        double erf2 = erf2_sign * TMath::Erf(erf2_arg);
        
        return exp_term * (erf1 + erf2);
    };
    
    // Create TF1 for S_total(x) = sum of [(1-w)*G_n + w*I_{G_n⊗E}] for n=1 to npe_max
    TF1 *f_stotal = new TF1("f_stotal", [IGnE_func](double *x, double *par) {
        double mu = par[0];
        double w = par[1];
        double Q_0 = par[2];
        double alpha = par[4];
        int npe_max = (int)par[5];
        double Q_1 = par[6];
        double std_1 = par[7];
        
        double sum = 0.0;
        for (int n = 1; n <= npe_max; n++) {
            // Poisson term: m1^n * e^(-m1) / n!
            double poisson = TMath::Power(mu, n) * TMath::Exp(-mu) / TMath::Gamma(n + 1);
            
            // G_n(x - Q_0)
            double Q_n = Q_0 + n * Q_1;
            double sigma_n = std_1 * TMath::Sqrt(n);
            double norm = 1.0 / (sigma_n * TMath::Sqrt(2.0 * TMath::Pi()));
            double G_n = norm * TMath::Exp(-TMath::Power(x[0] - Q_n, 2) / (2.0 * sigma_n * sigma_n));
            
            // I_{G_n⊗E}(x - Q_0) from Eq. (8)
            double IGnE = IGnE_func(x[0], Q_0, n, Q_1, std_1, alpha);
            
            // Sum: poisson * [(1-w)*G_n + w*I_{G_n⊗E}]
            sum += poisson * ((1.0 - w) * G_n + w * IGnE);
        }
        return sum;
    }, fit_min, fit_max, 8);
    f_stotal->SetParameters(mu, w, Q_0, std_0, alpha, npe_max, Q_1, std_1);
    f_stotal->SetNpx(1000);
    f_stotal->SetLineColor(kRed);
    f_stotal->SetLineStyle(2);
    f_stotal->SetLineWidth(2);
    
    // Create TF1 for Total T(x) = Pedestal (n=0) + S_total (n>=1)
    // This is the full Equation (7)
    TF1 *f_total = new TF1("f_total", [IGnE_func](double *x, double *par) {
        double mu = par[0];
        double w = par[1];
        double Q_0 = par[2];
        double std_0 = par[3];
        double alpha = par[4];
        int npe_max = (int)par[5];
        double Q_1 = par[6];
        double std_1 = par[7];
        
        // ============================================
        // n = 0 term (Pedestal)
        // ============================================
        double norm_0 = 1.0 / (std_0 * TMath::Sqrt(2.0 * TMath::Pi()));
        double G_0 = norm_0 * TMath::Exp(-TMath::Power(x[0] - Q_0, 2) / (2.0 * std_0 * std_0));
        
        double IGE_0 = 0.0;
        if (x[0] >= Q_0) {
            IGE_0 = alpha * TMath::Exp(-alpha * (x[0] - Q_0));
        }
        
        double n0_term = TMath::Exp(-mu) * ((1.0 - w) * G_0 + w * IGE_0);
        
        // ============================================
        // n >= 1 terms
        // ============================================
        double sum_n = 0.0;
        for (int n = 1; n <= npe_max; n++) {
            // Poisson
            double poisson = TMath::Power(mu, n) * TMath::Exp(-mu) / TMath::Gamma(n + 1);
            
            // G_n(x - Q_0)
            double Q_n = Q_0 + n * Q_1;
            double sigma_n = std_1 * TMath::Sqrt(n);
            double norm_n = 1.0 / (sigma_n * TMath::Sqrt(2.0 * TMath::Pi()));
            double G_n = norm_n * TMath::Exp(-TMath::Power(x[0] - Q_n, 2) / (2.0 * sigma_n * sigma_n));
            
            // I_{G_n⊗E}(x - Q_0)
            double IGnE = IGnE_func(x[0], Q_0, n, Q_1, std_1, alpha);
            
            sum_n += poisson * ((1.0 - w) * G_n + w * IGnE);
        }
        
        return n0_term + sum_n;
    }, fit_min, fit_max, 8);
    
    // Set initial parameters for fitting
    f_total->SetParameters(mu, w, Q_0, std_0, alpha, npe_max, Q_1, std_1);
    f_total->SetParNames("mu", "w", "Q_0", "std_0", "alpha", "npe_max", "Q_1", "std_1");
    
    // Fix npe_max since it's an integer constant
    f_total->FixParameter(5, npe_max);
    
    f_total->SetNpx(1000);
    f_total->SetLineColor(kMagenta);
    f_total->SetLineWidth(3);
    
    // Perform the fit
    h->Fit(f_total, "R");
    
    // Update component functions with fitted parameters
    f_pedestal->SetParameters(f_total->GetParameter(0), f_total->GetParameter(1), 
                              f_total->GetParameter(2), f_total->GetParameter(3),
                              f_total->GetParameter(4));
    
    f_stotal->SetParameters(f_total->GetParameter(0), f_total->GetParameter(1), 
                           f_total->GetParameter(2), f_total->GetParameter(3), 
                           f_total->GetParameter(4), npe_max, 
                           f_total->GetParameter(6), f_total->GetParameter(7));
    
    // Draw histogram and functions
    h->Draw("HIST");
    f_pedestal->Draw("SAME");
    f_stotal->Draw("SAME");
    f_total->Draw("SAME");
   
    // stat box 위치 가져오기
    TPaveStats *st = (TPaveStats*)h->GetListOfFunctions()->FindObject("stats");
    double x1 = st->GetX1NDC();
    double y2 = st->GetY2NDC();

    // Legend 생성
    TLegend *leg = new TLegend(0.5, y2-0.3, x1-0.05, y2-0.01);
    leg->SetBorderSize(0);
    leg->SetMargin(0.1);

    leg->AddEntry(h, "Data", "l");
    leg->AddEntry(f_pedestal, "Pedestal (n=0)", "l");
    leg->AddEntry(f_stotal, "S_{total}(x) (n#geq1)", "l");
    leg->AddEntry(f_total, "Total T(x)", "l");
    leg->Draw();
    
    can->Update();
    
    // Print fitted parameters
    std::cout << "Chi2/NDF   = " << f_total->GetChisquare() / (double) f_total->GetNDF() << std::endl;
}



void PMTCalib_SPE_Fit_2(TH1D *h, TCanvas *can,
                    double fit_min = -2.0,
                    double fit_max = 12.0,

                    double mu = 0.061,
                    double Q_0 = 0.035,
                    
                    double std_0 = 0.101,
                    
                    int npe_max = 3,
                    double Q_1 = 2.92,
                    double std_1 = 0.72) {

    can->cd();

    TF1 *f_pedestal = new TF1("f_pedestal", [](double *x, double *par) {
        double mu = par[0];
        double Q_0 = par[1];
        double std_0 = par[2];

        double norm = 1.0 / (std_0 * TMath::Sqrt(2.0 * TMath::Pi()));
        double exponent = -TMath::Power(x[0] - Q_0, 2) / (2.0 * std_0 * std_0);
        return TMath::Exp(-mu) * norm * TMath::Exp(exponent);
    }, fit_min, fit_max, 3);
    f_pedestal->SetParameters(mu, Q_0, std_0);
    f_pedestal->SetNpx(1000);
    f_pedestal->SetLineColor(kBlack);
    f_pedestal->SetLineStyle(2);
    f_pedestal->SetLineWidth(2);

    TF1 *f_stotal = new TF1("f_stotal", [](double *x, double *par) {
        double mu = par[0];
        double Q_0 = par[1];
        int npe_max = (int)par[3];
        double Q_1 = par[4];
        double std_1 = par[5];

        double sum = 0.0;
        for (int n = 1; n <= npe_max; n++) {
            if (n <= 0) continue;
            double m1_n = TMath::Power(mu, n) * TMath::Exp(-mu) / TMath::Gamma(n + 1);
            double norm = 1.0 / (std_1 * TMath::Sqrt(2.0 * TMath::Pi() * n));
            double arg = x[0] - n * Q_1 - Q_0;
            double exponent = -TMath::Power(arg, 2) / (2.0 * n * std_1 * std_1);
            sum += m1_n * norm * TMath::Exp(exponent);
        }
        return sum;
    }, fit_min, fit_max, 6);
    f_stotal->SetParameters(mu, Q_0, std_0, npe_max, Q_1, std_1);
    f_stotal->SetNpx(1000);
    f_stotal->SetLineColor(kRed);
    f_stotal->SetLineStyle(2);
    f_stotal->SetLineWidth(2);

    TF1 *f_total = new TF1("f_total", [](double *x, double *par) {
        double mu = par[0];
        double Q_0 = par[1];
        double std_0 = par[2];
        int npe_max = (int)par[3];
        double Q_1 = par[4];
        double std_1 = par[5];

        double norm_p = 1.0 / (std_0 * TMath::Sqrt(2.0 * TMath::Pi()));
        double exp_p = -TMath::Power(x[0] - Q_0, 2) / (2.0 * std_0 * std_0);
        double P = TMath::Exp(-mu) * norm_p * TMath::Exp(exp_p);

        double S_total = 0.0;
        for (int n = 1; n <= npe_max; n++) {
            if (n <= 0) continue;
            double m1_n = TMath::Power(mu, n) * TMath::Exp(-mu) / TMath::Gamma(n + 1);
            double norm = 1.0 / (std_1 * TMath::Sqrt(2.0 * TMath::Pi() * n));
            double arg = x[0] - n * Q_1 - Q_0;
            double exponent = -TMath::Power(arg, 2) / (2.0 * n * std_1 * std_1);
            S_total += m1_n * norm * TMath::Exp(exponent);
        }

        return P + S_total;
    }, fit_min, fit_max, 6);

    f_total->SetParameters(mu, Q_0, std_0, npe_max, Q_1, std_1);
    f_total->SetParNames("mu", "Q_0", "std_0", "npe_max", "Q_1", "std_1");
    f_total->FixParameter(3, npe_max);
    f_total->SetNpx(1000);
    f_total->SetLineColor(kMagenta);
    f_total->SetLineWidth(3);

    h->Fit(f_total, "R");

    f_pedestal->SetParameters(f_total->GetParameter(0), f_total->GetParameter(1),
                              f_total->GetParameter(2));

    f_stotal->SetParameters(f_total->GetParameter(0), f_total->GetParameter(1),
                           f_total->GetParameter(2), npe_max,
                           f_total->GetParameter(4), f_total->GetParameter(5));

    h->Draw("HIST");
    f_pedestal->Draw("SAME");
    f_stotal->Draw("SAME");
    f_total->Draw("SAME");

    TLegend *leg = new TLegend(0.6, 0.6, 0.88, 0.88);
    leg->AddEntry(h, "Data", "l");
    leg->AddEntry(f_pedestal, "Pedestal P(x)", "l");
    leg->AddEntry(f_stotal, "S_{total}(x)", "l");
    leg->AddEntry(f_total, "Total T(x)", "l");
    leg->Draw();

    can->Update();

    std::cout << "Chi2/NDF   = " << f_total->GetChisquare() / (double) f_total->GetNDF() << std::endl;
}



void PMTCalib_SPE_Fit_3(TH1D *h, TCanvas *can,
                    double fit_min = -2.0,
                    double fit_max = 12.0,


                    double mu = 0.061,
                    double w = 0.0175,
                    
                    double Q_0 = 0.035,
                    double std_0 = 0.101,
                    double alpha = 0.65,
                    
                    int npe_max = 3,
                    double Q_1 = 2.92,
                    double std_1 = 0.72) {
    
    can->cd();
    
    // Define S(x,n) calculation inline to avoid lambda capture issues
    
    // Create TF1 for Pedestal P(x)
    // P(x) = (1-w) * e^(-m1) / (s_pedestal * sqrt(2*pi)) * exp(-(x-Q_0)^2 / (2*s_pedestal^2))
    TF1 *f_pedestal = new TF1("f_pedestal", [](double *x, double *par) {
        double mu = par[0];
        double w = par[1];
        double Q_0 = par[2];
        double std_0 = par[3];
        
        double norm = 1.0 / (std_0 * TMath::Sqrt(2.0 * TMath::Pi()));
        double exponent = -TMath::Power(x[0] - Q_0, 2) / (2.0 * std_0 * std_0);
        return (1.0 - w) * TMath::Exp(-mu) * norm * TMath::Exp(exponent);
    }, fit_min, fit_max, 4);
    f_pedestal->SetParameters(mu, w, Q_0, std_0);
    f_pedestal->SetNpx(1000);
    f_pedestal->SetLineColor(kBlack);
    f_pedestal->SetLineStyle(2);
    f_pedestal->SetLineWidth(2);
    
    // Create TF1 for Exponential E(x)
    // E(x) = e^(-m1) * w * A * exp(-A*(x-Q_0)) for x >= Q_0, 0 otherwise
    TF1 *f_exponential = new TF1("f_exponential", [](double *x, double *par) {
        double mu = par[0];
        double w = par[1];
        double Q_0 = par[2];
        double alpha = par[4];
        
        if (x[0] >= Q_0) {
            return TMath::Exp(-mu) * w * alpha * TMath::Exp(-alpha * (x[0] - Q_0));
        } else {
            return 0.0;
        }
    }, fit_min, fit_max, 5);
    f_exponential->SetParameters(mu, w, Q_0, std_0, alpha);
    f_exponential->SetNpx(1000);
    f_exponential->SetLineColor(kBlack);
    f_exponential->SetLineStyle(10);
    f_exponential->SetLineWidth(2);
    
    // Create TF1 for S_total(x) = sum of S(x,n) for n=1 to 3
    TF1 *f_stotal = new TF1("f_stotal", [](double *x, double *par) {
        double mu = par[0];
         double w = par[1];
        double Q_0 = par[2];
        double alpha = par[4];
        int npe_max = (int)par[5];
        double Q_1 = par[6];
        double std_1 = par[7];
        
        double sum = 0.0;
        for (int n = 1; n <= npe_max; n++) {
            if (n <= 0) continue;
            // Poisson term: m1^n * e^(-m1) / n!
            double m1_n = TMath::Power(mu, n) * TMath::Exp(-mu) / TMath::Gamma(n + 1);
            // Normalization: 1/(s_spe * sqrt(2*pi*n))
            double norm = 1.0 / (std_1 * TMath::Sqrt(2.0 * TMath::Pi() * n));
            // Argument: (x - n*Q - Q_0 - w/A)^2
            double arg = x[0] - n * Q_1 - Q_0 - w / alpha;
            // Exponent: -arg^2 / (2*n*s_spe^2)
            double exponent = -TMath::Power(arg, 2) / (2.0 * n * std_1 * std_1);
            sum += m1_n * norm * TMath::Exp(exponent);
        }
        return sum;
    }, fit_min, fit_max, 8);
    f_stotal->SetParameters(mu, w, Q_0, std_0, alpha, npe_max, Q_1, std_1);
    f_stotal->SetNpx(1000);
    f_stotal->SetLineColor(kRed);
    f_stotal->SetLineStyle(2);
    f_stotal->SetLineWidth(2);
    
    // Create TF1 for Total T(x) = P(x) + E(x) + S_total(x)
    TF1 *f_total = new TF1("f_total", [](double *x, double *par) {
        double mu = par[0];
        double w = par[1];
        double Q_0 = par[2];
        double std_0 = par[3];
        double alpha = par[4];
        int npe_max = (int)par[5];
        double Q_1 = par[6];
        double std_1 = par[7];
        
        // Pedestal P(x)
        double norm_p = 1.0 / (std_0 * TMath::Sqrt(2.0 * TMath::Pi()));
        double exp_p = -TMath::Power(x[0] - Q_0, 2) / (2.0 * std_0 * std_0);
        double P = (1.0 - w) * TMath::Exp(-mu) * norm_p * TMath::Exp(exp_p);
        
        // Exponential E(x)
        double E = 0.0;
        if (x[0] >= Q_0) {
            E = TMath::Exp(-mu) * w * alpha * TMath::Exp(-alpha * (x[0] - Q_0));
        }
        
        // S_total(x)
        double S_total = 0.0;
        for (int n = 1; n <= npe_max; n++) {
            if (n <= 0) continue;
            // Poisson term: m1^n * e^(-m1) / n!
            double m1_n = TMath::Power(mu, n) * TMath::Exp(-mu) / TMath::Gamma(n + 1);
            // Normalization: 1/(s_spe * sqrt(2*pi*n))
            double norm = 1.0 / (std_1 * TMath::Sqrt(2.0 * TMath::Pi() * n));
            // Argument: (x - n*Q - Q_0 - w/A)^2
            double arg = x[0] - n * Q_1 - Q_0 - w / alpha;
            // Exponent: -arg^2 / (2*n*s_spe^2)
            double exponent = -TMath::Power(arg, 2) / (2.0 * n * std_1 * std_1);
            S_total += m1_n * norm * TMath::Exp(exponent);
        }
        
        return P + E + S_total;
    }, fit_min, fit_max, 8);
    
    // Set initial parameters for fitting
    f_total->SetParameters(mu, w, Q_0, std_0, alpha, npe_max, Q_1, std_1);
    f_total->SetParNames("mu", "w", "Q_0", "std_0", "alpha", "npe_max", "Q_1", "std_1");
    
    // Fix npe_max since it's an integer constant
    f_total->FixParameter(5, npe_max);
    
    f_total->SetNpx(1000);
    f_total->SetLineColor(kMagenta);
    f_total->SetLineWidth(3);
    
    // Perform the fit
    h->Fit(f_total, "R");  // "R" uses the range specified in TF1
    
    // Update component functions with fitted parameters
    f_pedestal->SetParameters(f_total->GetParameter(0), f_total->GetParameter(1), 
                              f_total->GetParameter(2), f_total->GetParameter(3));
    
    f_exponential->SetParameters(f_total->GetParameter(0), f_total->GetParameter(1), 
                                 f_total->GetParameter(2), f_total->GetParameter(3), 
                                 f_total->GetParameter(4));
    
    f_stotal->SetParameters(f_total->GetParameter(0), f_total->GetParameter(1), 
                           f_total->GetParameter(2), f_total->GetParameter(3), 
                           f_total->GetParameter(4), npe_max, 
                           f_total->GetParameter(6), f_total->GetParameter(7));
    
    // Draw histogram and functions
    h->Draw("HIST");
    f_pedestal->Draw("SAME");
    f_exponential->Draw("SAME");
    f_stotal->Draw("SAME");
    f_total->Draw("SAME");
    
    // Create legend
    TLegend *leg = new TLegend(0.6, 0.6, 0.88, 0.88);
    leg->AddEntry(h, "Data", "l");
    leg->AddEntry(f_pedestal, "Pedestal P(x)", "l");
    leg->AddEntry(f_exponential, "Exponential E(x)", "l");
    leg->AddEntry(f_stotal, "S_{total}(x)", "l");
    leg->AddEntry(f_total, "Total T(x)", "l");
    leg->Draw();
    
    can->Update();
    
    // Print fitted parameters
    std::cout << "Chi2/NDF   = " << f_total->GetChisquare() / (double) f_total->GetNDF() << std::endl;
}



void PMTCalib_SPE_Fit_4(TH1D *h, TCanvas *can,
                    double fit_min = -2.0,
                    double fit_max = 12.0,
                    double p0 = 4.0,   // pedestal amplitude
                    double p1 = 0.0,     // pedestal position
                    double p2 = 0.01,     // pedestal width^2
                    double p3 = 0.05,    // SPE amplitude
                    double p4 = 3.0,     // SPE position
                    double p5 = 0.5,     // SPE width^2
                    double p6 = 0.2) {   // inelastic fraction
    can->cd();
    
    // Pedestal peak function
    TF1 *f_pedestal = new TF1("f_pedestal", [](double *x, double *par) {
        double p0 = par[0];
        double p1 = par[1];
        double p2 = par[2];
        double arg = (x[0] - p1) / TMath::Sqrt(2.0 * p2);
        return p0 * TMath::Exp(-arg * arg);
    }, fit_min, fit_max, 3);
    
    f_pedestal->SetParameters(p0, p1, p2);
    f_pedestal->SetNpx(1000);
    f_pedestal->SetLineColor(kBlack);
    f_pedestal->SetLineStyle(2);
    f_pedestal->SetLineWidth(2);
    
    // 1 p.e peak function
    TF1 *f_stotal = new TF1("f_stotal", [](double *x, double *par) {
        double p3 = par[3];
        double p4 = par[4];
        double p5 = par[5];
        double arg = (x[0] - p4) / TMath::Sqrt(2.0 * p5);
        return p3 * TMath::Exp(-arg * arg);
    }, fit_min, fit_max, 6);
    
    f_stotal->SetParameters(p0, p1, p2, p3, p4, p5);
    f_stotal->SetNpx(1000);
    f_stotal->SetLineColor(kRed);
    f_stotal->SetLineStyle(2);
    f_stotal->SetLineWidth(2);
    
    // Inelastic scattering function
    TF1 *f_inelastic = new TF1("f_inelastic", [](double *x, double *par) {
        double p1 = par[1];
        double p2 = par[2];
        double p3 = par[3];
        double p4 = par[4];
        double p5 = par[5];
        double p6 = par[6];
        
        double arg1 = (x[0] - p1) / TMath::Sqrt(2.0 * p2);
        double arg2 = (x[0] - p4) / TMath::Sqrt(2.0 * p5);
        
        double erf1 = TMath::Erf(arg1);
        double erf2 = TMath::Erf(arg2);
        
        return (p6 * p3 / 2.0) * (erf1 - erf2);
    }, fit_min, fit_max, 7);
    
    f_inelastic->SetParameters(p0, p1, p2, p3, p4, p5, p6);
    f_inelastic->SetNpx(1000);
    f_inelastic->SetLineColor(kBlack);
    f_inelastic->SetLineStyle(10);
    f_inelastic->SetLineWidth(2);
    
    // Total function
    TF1 *f_total = new TF1("f_total", [](double *x, double *par) {
        double p0 = par[0];
        double p1 = par[1];
        double p2 = par[2];
        double p3 = par[3];
        double p4 = par[4];
        double p5 = par[5];
        double p6 = par[6];
        
        // Pedestal peak
        double arg_ped = (x[0] - p1) / TMath::Sqrt(2.0 * p2);
        double pedestal = p0 * TMath::Exp(-arg_ped * arg_ped);
        
        // 1 p.e peak
        double arg_spe = (x[0] - p4) / TMath::Sqrt(2.0 * p5);
        double spe = p3 * TMath::Exp(-arg_spe * arg_spe);
        
        // Inelastic scattering
        double arg1 = (x[0] - p1) / TMath::Sqrt(2.0 * p2);
        double arg2 = (x[0] - p4) / TMath::Sqrt(2.0 * p5);
        double erf1 = TMath::Erf(arg1);
        double erf2 = TMath::Erf(arg2);
        double inelastic = (p6 * p3 / 2.0) * (erf1 - erf2);
        
        return pedestal + spe + inelastic;
    }, fit_min, fit_max, 7);
    
    f_total->SetParameters(p0, p1, p2, p3, p4, p5, p6);
    f_total->SetParNames("p0", "p1", "p2", "p3", "p4", "p5", "p6");
    f_total->SetNpx(1000);
    f_total->SetLineColor(kMagenta);
    f_total->SetLineWidth(3);
    
    // Fit
    h->Fit(f_total, "R");
    
    // Update individual functions with fitted parameters
    f_pedestal->SetParameters(f_total->GetParameter(0), f_total->GetParameter(1),
                              f_total->GetParameter(2));
    f_stotal->SetParameters(f_total->GetParameter(0), f_total->GetParameter(1),
                           f_total->GetParameter(2), f_total->GetParameter(3),
                           f_total->GetParameter(4), f_total->GetParameter(5));
    f_inelastic->SetParameters(f_total->GetParameter(0), f_total->GetParameter(1),
                              f_total->GetParameter(2), f_total->GetParameter(3),
                              f_total->GetParameter(4), f_total->GetParameter(5),
                              f_total->GetParameter(6));
    
    // Draw
    h->Draw("HIST");
    f_pedestal->Draw("SAME");
    f_stotal->Draw("SAME");
    f_inelastic->Draw("SAME");
    f_total->Draw("SAME");
    
    TLegend *leg = new TLegend(0.6, 0.5, 0.88, 0.88);
    leg->AddEntry(h, "Data", "l");
    leg->AddEntry(f_pedestal, "Pedestal P(x)", "l");
    leg->AddEntry(f_inelastic, "Inelastic scattering", "l");
     leg->AddEntry(f_stotal, "1 p.e peak", "l");
    leg->AddEntry(f_total, "Total T(x)", "l");
    leg->Draw();
    
    can->Update();
    
    std::cout << "Chi2/NDF   = " << f_total->GetChisquare() / (double) f_total->GetNDF() << std::endl;
}
