// FILEPATH: symFIR.cpp
#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>

// Function to calculate frequency response H(e^jw) of the filter
void calculateFrequencyResponse(const std::vector<double>& h, std::vector<double>& w, std::vector<std::complex<double>>& H) {
    int N = h.size();
    w.resize(N);
    H.resize(N);

    for (int i = 0; i < N; i++) {
        w[i] = 2 * M_PI * i / N;
        H[i] = 0;

        for (int j = 0; j < N; j++) {
            H[i] += h[j] * std::polar(1.0, -w[i] * j);
        }
    }
}

// Function to calculate phase delay of the filter
void calculatePhaseDelay(const std::vector<std::complex<double>>& H, std::vector<double>& phase_delay) {
    int N = H.size();
    phase_delay.resize(N);

    for (int i = 0; i < N; i++) {
        phase_delay[i] = std::arg(H[i]);
    }

    std::transform(phase_delay.begin(), phase_delay.end(), phase_delay.begin(), [](double phase) {
        return -phase;
    });
}

// Function to calculate group delay of the filter
void calculateGroupDelay(const std::vector<double>& h, std::vector<double>& w_gd, std::vector<double>& group_delay) {
    int N = h.size();
    w_gd.resize(N);
    group_delay.resize(N);

    for (int i = 0; i < N; i++) {
        w_gd[i] = 2 * M_PI * i / N;
        group_delay[i] = 0;

        for (int j = 0; j < N; j++) {
            group_delay[i] += h[j] * (j - (N - 1) / 2) * std::sin(w_gd[i] * (j - (N - 1) / 2));
        }

        group_delay[i] /= -M_PI;
    }
}

int main() {
    std::vector<double> h = {1, 1};
    std::vector<double> w, phase_delay, w_gd, group_delay;
    std::vector<std::complex<double>> H;

    // Calculate frequency response H(e^jw) of the filter
    calculateFrequencyResponse(h, w, H);

    // Calculate phase delay of the filter
    calculatePhaseDelay(H, phase_delay);

    // Calculate group delay of the filter
    calculateGroupDelay(h, w_gd, group_delay);

    // Plotting Frequency Response
    // ...

    // Plotting Phase Delay
    // ...

    // Plotting Group Delay
    // ...

    return 0;
}
