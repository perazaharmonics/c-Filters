// FILEPATH: Untitled-1.cpp
#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>

// Function to calculate frequency response
void freqz(const std::vector<double>& h, std::vector<double>& w, std::vector<std::complex<double>>& H) {
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

// Function to calculate phase delay
void phaseDelay(const std::vector<std::complex<double>>& H, std::vector<double>& phaseDelay) {
    int N = H.size();
    phaseDelay.resize(N);

    for (int i = 0; i < N; i++) {
        phaseDelay[i] = std::arg(H[i]);
    }

    std::partial_sum(phaseDelay.begin(), phaseDelay.end(), phaseDelay.begin());
}

// Function to calculate group delay
void groupDelay(const std::vector<double>& h, std::vector<double>& w, std::vector<double>& groupDelay) {
    int N = h.size();
    w.resize(N);
    groupDelay.resize(N);

    for (int i = 0; i < N; i++) {
        w[i] = 2 * M_PI * i / N;
        groupDelay[i] = 0;
        for (int j = 0; j < N; j++) {
            groupDelay[i] += h[j] * (j - (N - 1) / 2) * std::sin(w[i] * (j - (N - 1) / 2));
        }
        groupDelay[i] /= -2 * M_PI;
    }
}

int main() {
    // Define the impulse response
    std::vector<double> h = {1, -1};

    // a) Determine the frequency response H(e^jw) of the filter.
    std::vector<double> w;
    std::vector<std::complex<double>> H;
    freqz(h, w, H);

    // b) Determine the phase delay of the filter.
    std::vector<double> phaseDelay;
    phaseDelay(H, phaseDelay);

    // c) Determine the group delay of the filter.
    std::vector<double> w_gd;
    std::vector<double> groupDelay;
    groupDelay(h, w_gd, groupDelay);

    // Print the results
    std::cout << "Frequency Response:" << std::endl;
    for (int i = 0; i < H.size(); i++) {
        std::cout << "w: " << w[i] << ", H: " << std::abs(H[i]) << std::endl;
    }

    std::cout << "Phase Delay:" << std::endl;
    for (int i = 0; i < phaseDelay.size(); i++) {
        std::cout << "w: " << w[i] << ", Phase Delay: " << phaseDelay[i] << std::endl;
    }

    std::cout << "Group Delay:" << std::endl;
    for (int i = 0; i < groupDelay.size(); i++) {
        std::cout << "w: " << w_gd[i] << ", Group Delay: " << groupDelay[i] << std::endl;
    }

    return 0;
}
