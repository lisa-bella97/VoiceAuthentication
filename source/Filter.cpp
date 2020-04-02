#include <cmath>
#include <algorithm>
#include <ctime>
#include "Filter.h"


std::vector<unsigned char> Filter::removePauses(const std::vector<unsigned char> &rawSamples, int sampleRate) {
    auto samples = convertFromRaw(rawSamples);

    // Определение окна для Фурье-преобразования
    size_t N = samples.size(); // Число отсчетов сигнала
    double T = 1.0 / sampleRate; // Период дискретизации
    double n = 0.02 / T; // Количество отсчетов в периоде стационарности
    auto s = (int) (log(n) / log(2)); // Степень двойки (log2(n) = log(n) / log(2))
    auto nInFrame = (size_t) pow(2, s); // Количество отсчетов в кадре
    size_t K = N / nInFrame; // Количество итераций Фурье-преобразования

    std::vector<std::pair<double, int>> framePowers(K);

    for (auto i = 0; i < K; ++i) {
        std::vector<double> frameSamples(nInFrame);
        for (auto j = 0; j < nInFrame; ++j) {
            frameSamples[j] = samples[i * nInFrame + j];
        }

        auto fftResult = FFTAnalysis(frameSamples);

        double power = 0.0;
        for (auto j = 0; j < nInFrame; ++j) {
            power += fftResult[j] * fftResult[j];
        }
        framePowers[i] = {0.5 * power, i};
    }

    auto clusters = kMeans(framePowers);
    std::vector<short> result;
    int prevIndex = 0;

    std::sort(clusters[2].begin(), clusters[2].end(), [](auto &left, auto &right) {
        return left.second < right.second;
    });

    for (const auto &el : clusters[2]) {
        for (auto i = prevIndex; i < el.second; ++i) {
            for (auto j = 0; j < nInFrame; ++j) {
                result.push_back(samples[i * nInFrame + j]);
            }
        }
        prevIndex = el.second + 1;
    }

    return convertToRaw(result);
}

std::vector<std::vector<double>>
Filter::getMFCC(const std::vector<unsigned char> &rawSamples, int sampleRate, int mfccNumber) {
    auto samples = convertFromRaw(rawSamples);
    int samplesSize = samples.size();

    int stride = 160; // шаг
    int length_frame = 400;
    int length_DFT = 512;
    int number_coefficients = 13;
    int number_filterbanks = 26;

    int number_feature_vectors;
    int nSamplesPerSec = sampleRate;

    double pi = 3.14159265358979323846;

    double **feature_vector;

    feature_vector = new double *[number_feature_vectors = (samplesSize - length_frame) / stride + 1];

    for (int i = 0; i < number_feature_vectors; i++) {
        feature_vector[i] = new double[3 * number_coefficients];
    }

    // MFCC
    for (int i = 0; i <= samplesSize - length_frame; i += stride) {
        auto *frame = new double[length_frame];
        // pre-emphasis
        for (int j = 0; j < length_frame; j++) {
            if (i + j < samplesSize) {
                frame[j] = samples[i + j] - 0.95 * samples[i + j - 1];
            } else {
                frame[j] = 0;
            }
        }

        // windowing
        for (int j = 0; j < length_frame; j++) {
            frame[j] *= 0.54 - 0.46 * cos(2 * pi * j / (length_frame - 1));
        }

        MFCC(length_frame, length_DFT, number_coefficients, number_filterbanks, nSamplesPerSec, frame,
             feature_vector[i / stride]);

        delete[] frame;
    }

    // deltas
    for (int i = 0; i < number_feature_vectors; i++) {
        int prev = (i == 0) ? (0) : (i - 1);
        int next = (i == number_feature_vectors - 1) ? (number_feature_vectors - 1) : (i + 1);

        for (int j = 0; j < number_coefficients; j++) {
            feature_vector[i][number_coefficients + j] = (feature_vector[next][j] - feature_vector[prev][j]) / 2;
        }
    }

    // delta-deltas
    for (int i = 0; i < number_feature_vectors; i++) {
        int prev = (i == 0) ? (0) : (i - 1);
        int next = (i == number_feature_vectors - 1) ? (number_feature_vectors - 1) : (i + 1);

        for (int j = number_coefficients; j < 2 * number_coefficients; j++) {
            feature_vector[i][number_coefficients + j] = (feature_vector[next][j] - feature_vector[prev][j]) / 2;
        }
    }

    std::vector<std::vector<double>> result;
    std::vector<double> avg(3 * number_coefficients);

    for (int i = 0; i < number_feature_vectors; i++) {
        result.emplace_back(3 * number_coefficients);
        for (int j = 0; j < 3 * number_coefficients; j++) {
            result[i][j] = feature_vector[i][j];
            avg[j] += feature_vector[i][j];
        }
    }

    for (int j = 0; j < 3 * number_coefficients; j++) {
        avg[j] /= 3 * number_coefficients;
    }

    for (int i = 0; i < number_feature_vectors; i++) {
        delete[] feature_vector[i];
    }
    delete[] feature_vector;

    return result;
}

std::vector<short> Filter::convertFromRaw(const std::vector<unsigned char> &samples) {
    const auto sampleSize = sizeof(short);
    const auto samplesSize = samples.size();
    std::vector<short> result(samplesSize / sampleSize);

    for (auto i = 0, j = 0; i < samplesSize; i += sampleSize) {
        short current = samples[i];
        for (auto k = i + 1; k < i + sampleSize; ++k) {
            current <<= 8;
            current |= samples[k];
        }
        result[j++] = current;
    }

    return result;
}

std::vector<unsigned char> Filter::convertToRaw(const std::vector<short> &samples) {
    const auto sampleSize = sizeof(short);
    const auto samplesSize = samples.size();
    std::vector<unsigned char> result(samplesSize * sampleSize);

    for (auto i = 0; i < samplesSize - 1; ++i) {
        short current = samples[i];
        result[i * sampleSize + sampleSize - 1] = current;
        for (int k = sampleSize - 2; k >= 0; --k) {
            current >>= 8;
            result[i * sampleSize + k] = current;
        }
    }

    return result;
}

// data - массив анализируемых данных, длина массива должна быть кратна степени 2.
std::vector<double> Filter::FFTAnalysis(const std::vector<double> &data) {
    int i, j, n, m, Mmax, Istp;
    double Tmpr, Tmpi, Wtmp, Theta;
    double Wpr, Wpi, Wr, Wi;
    n = data.size() * 2;
    std::vector<double> temp(n);

    for (i = 0; i < n; i += 2) {
        temp[i] = 0;
        temp[i + 1] = data[i / 2];
    }

    i = 1;
    j = 1;
    while (i < n) {
        if (j > i) {
            Tmpr = temp[i];
            temp[i] = temp[j];
            temp[j] = Tmpr;
            Tmpr = temp[i + 1];
            temp[i + 1] = temp[j + 1];
            temp[j + 1] = Tmpr;
        }
        i = i + 2;
        m = data.size();
        while ((m >= 2) && (j > m)) {
            j = j - m;
            m = m >> 1;
        }
        j = j + m;
    }

    Mmax = 2;
    const double TwoPi = 2 * M_PI;
    while (n > Mmax) {
        Theta = -TwoPi / Mmax;
        Wpi = sin(Theta);
        Wtmp = sin(Theta / 2);
        Wpr = Wtmp * Wtmp * 2;
        Istp = Mmax * 2;
        Wr = 1;
        Wi = 0;
        m = 1;

        while (m < Mmax) {
            i = m;
            m = m + 2;
            Tmpr = Wr;
            Tmpi = Wi;
            Wr = Wr - Tmpr * Wpr - Tmpi * Wpi;
            Wi = Wi + Tmpr * Wpi - Tmpi * Wpr;

            while (i < n) {
                j = i + Mmax;
                Tmpr = Wr * temp[j] - Wi * temp[j - 1];
                Tmpi = Wi * temp[j] + Wr * temp[j - 1];

                temp[j] = temp[i] - Tmpr;
                temp[j - 1] = temp[i - 1] - Tmpi;
                temp[i] = temp[i] + Tmpr;
                temp[i - 1] = temp[i - 1] + Tmpi;
                i = i + Istp;
            }
        }

        Mmax = Istp;
    }

    std::vector<double> result(data.size());

    for (i = 0; i < data.size(); i++) {
        j = i * 2;
        result[i] = 2 * sqrt(pow(temp[j], 2) + pow(temp[j + 1], 2)) / data.size();
    }

    return result;
}

std::vector<std::vector<std::pair<double, int>>> Filter::kMeans(const std::vector<std::pair<double, int>> &values) {
    std::vector<std::vector<std::pair<double, int>>> result = {std::vector<std::pair<double, int>>(),
                                                               std::vector<std::pair<double, int>>(),
                                                               std::vector<std::pair<double, int>>()};
    const auto cluster1 = &result[0];
    const auto cluster2 = &result[1];
    const auto cluster3 = &result[2];

    int c1 = static_cast<int>(std::max_element(values.begin(), values.end())->first);
    int c2 = static_cast<int>(getMeanValue(values));
    int c3 = static_cast<int>(std::min_element(values.begin(), values.end())->first);

    while (true) {
        for (auto v : values) {
            int p1 = static_cast<int>(abs(c1 - v.first));
            int p2 = static_cast<int>(abs(c2 - v.first));
            int p3 = static_cast<int>(abs(c3 - v.first));

            int min = std::min(p1, std::min(p2, p3));
            if (p1 == min)
                cluster1->push_back(v);
            else if (p2 == min)
                cluster2->push_back(v);
            else
                cluster3->push_back(v);
        }

        int last_c1 = !cluster1->empty() ? static_cast<int>(getMeanValue(*cluster1)) : c1;
        int last_c2 = !cluster2->empty() ? static_cast<int>(getMeanValue(*cluster2)) : c2;
        int last_c3 = !cluster3->empty() ? static_cast<int>(getMeanValue(*cluster3)) : c3;

        if (c1 == last_c1 && c2 == last_c2 && c3 == last_c3)
            break;

        c1 = last_c1;
        c2 = last_c2;
        c3 = last_c3;

        cluster1->clear();
        cluster2->clear();
        cluster3->clear();
    }

    return result;
}

double Filter::getMeanValue(const std::vector<std::pair<double, int>> &values) {
    double avg = 0.0;
    for (const auto &v : values) {
        avg += v.first / values.size();
    }
    return avg;
}

void Filter::DCT(int direction, int length, double X[]) {
    if (direction == 1 || direction == -1) {
        double pi = 3.14159265358979323846;

        auto *x = new double[length];

        for (int i = 0; i < length; i++) {
            x[i] = X[i];
        }
        for (int k = 0; k < length; k++) {
            double sum = 0;

            if (direction == 1) {
                for (int n = 0; n < length; n++) {
                    sum += ((k == 0) ? (sqrt(0.5)) : (1)) * x[n] * cos(pi * (n + 0.5) * k / length);
                }
            } else if (direction == -1) {
                for (int n = 0; n < length; n++) {
                    sum += ((n == 0) ? (sqrt(0.5)) : (1)) * x[n] * cos(pi * n * (k + 0.5) / length);
                }
            }
            X[k] = sum * sqrt(2.0 / length);
        }
        delete[] x;
    }
}

void Filter::FFT(int direction, int length, double Xr[], double Xi[]) {
    int log_length = (int) log2((double) length);

    if (direction != 1 && direction != -1) {
        return;
    }
    if (1 << log_length != length) {
        return;
    }

    double pi = 3.14159265358979323846;

    for (int i = 0, j = 0; i < length; i++, j = 0) {
        for (int k = 0; k < log_length; k++) {
            j = (j << 1) | (1 & (i >> k));
        }
        if (j < i) {
            double t;

            t = Xr[i];
            Xr[i] = Xr[j];
            Xr[j] = t;

            t = Xi[i];
            Xi[i] = Xi[j];
            Xi[j] = t;
        }
    }
    for (int i = 0; i < log_length; i++) {
        int L = (int) pow(2.0, i);

        for (int j = 0; j < length - 1; j += 2 * L) {
            for (int k = 0; k < L; k++) {
                double argument = direction * -pi * k / L;

                double xr = Xr[j + k + L] * cos(argument) - Xi[j + k + L] * sin(argument);
                double xi = Xr[j + k + L] * sin(argument) + Xi[j + k + L] * cos(argument);

                Xr[j + k + L] = Xr[j + k] - xr;
                Xi[j + k + L] = Xi[j + k] - xi;
                Xr[j + k] = Xr[j + k] + xr;
                Xi[j + k] = Xi[j + k] + xi;
            }
        }
    }
    if (direction == -1) {
        for (int k = 0; k < length; k++) {
            Xr[k] /= length;
            Xi[k] /= length;
        }
    }
}

double Filter::Mel_Scale(int direction, double x) {
    switch (direction) {
        case -1:
            return 700 * (exp(x / 1125.0) - 1);
        case 1:
            return 1125 * log(1 + x / 700.0);
    }
    return 0;
}

void Filter::MFCC(int length_frame, int length_DFT, int number_coefficients, int number_filterbanks, int sample_rate,
                  const double frame[], double feature_vector[]) {
    double max_Mels_frequency = Mel_Scale(1, sample_rate / 2);
    double min_Mels_frequency = Mel_Scale(1, 300);
    double interval = (max_Mels_frequency - min_Mels_frequency) / (number_filterbanks + 1);

    auto *filterbank = new double[number_filterbanks];
    auto *Xr = new double[length_DFT];
    auto *Xi = new double[length_DFT];

    for (int i = 0; i < number_filterbanks; i++) {
        filterbank[i] = 0;
    }
    for (int i = 0; i < length_DFT; i++) {
        Xr[i] = (i < length_frame) ? (frame[i]) : (0);
        Xi[i] = 0;
    }
    FFT(1, length_DFT, Xr, Xi);

    for (int i = 0; i < length_DFT / 2 + 1; i++) {
        double frequency = (sample_rate / 2) * i / (length_DFT / 2);
        double Mel_frequency = Mel_Scale(1, frequency);
        double power = (Xr[i] * Xr[i] + Xi[i] * Xi[i]) / length_frame;

        for (int j = 0; j < number_filterbanks; j++) {
            double frequency_boundary[] = {min_Mels_frequency + interval * (j + 0),
                                           min_Mels_frequency + interval * (j + 1),
                                           min_Mels_frequency + interval * (j + 2)};

            if (frequency_boundary[0] <= Mel_frequency && Mel_frequency <= frequency_boundary[1]) {
                double lower_frequency = Mel_Scale(-1, frequency_boundary[0]);
                double upper_frequency = Mel_Scale(-1, frequency_boundary[1]);

                filterbank[j] += power * (frequency - lower_frequency) / (upper_frequency - lower_frequency);
            } else if (frequency_boundary[1] <= Mel_frequency && Mel_frequency <= frequency_boundary[2]) {
                double lower_frequency = Mel_Scale(-1, frequency_boundary[1]);
                double upper_frequency = Mel_Scale(-1, frequency_boundary[2]);

                filterbank[j] += power * (upper_frequency - frequency) / (upper_frequency - lower_frequency);
            }
        }
    }

    for (int i = 0; i < number_filterbanks; i++) {
        filterbank[i] = log(filterbank[i]);
    }
    DCT(1, number_filterbanks, filterbank);

    for (int i = 0; i < number_coefficients; i++) {
        feature_vector[i] = filterbank[i];
    }

    delete[] filterbank;
    delete[] Xr;
    delete[] Xi;
}
