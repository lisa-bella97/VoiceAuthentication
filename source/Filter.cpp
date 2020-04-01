#include <cmath>
#include <fft.h>
#include "Filter.h"


std::vector<unsigned char> Filter::normalizeVolume(const std::vector<unsigned char> &samples, int sampleSize) {
    auto samplesSize = samples.size();
    std::vector<unsigned char> result(samplesSize);

    if (samplesSize > sampleSize - 1) {
        for (auto i = 0; i < sampleSize; ++i) {
            result[i] = samples[i];
        }
    }

    for (auto i = sampleSize; i < samplesSize; i += sampleSize) {
        // TODO: сейчас это реализовано только для sampleSize = 2
        int current = samples[i];
        current <<= 8;
        current |= samples[i + 1];

        int previous = samples[i - 2];
        previous <<= 8;
        previous |= samples[i - 1];

        //auto res = short(current - 0.95*previous);

        int res = int(current - int(previous));

        result[i + 1] = res;
        res >>= 8;
        result[i] = res;
    }

    return result;
}

std::vector<unsigned char> Filter::removePauses(const std::vector<unsigned char> &rawSamples, int sampleSize,
                                                int sampleRate) {
    std::vector<unsigned char> result;
    auto samples = convert(rawSamples);

    // Определение окна для Фурье-преобразования
    size_t N = samples.size(); // Число отсчетов сигнала
    double T = 1.0 / sampleRate; // Период дискретизации
    double n = 0.02 / T; // Количество отсчетов в периоде стационарности
    auto s = (int) (log(n) / log(2)); // Степень двойки (log2(n) = log(n) / log(2))
    auto nInFrame = (size_t) pow(2, s); // Количество отсчетов в кадре
    size_t K = N / nInFrame; // Количество итераций Фурье-преобразования

    std::vector<double> framePowers(K);

    for (auto i = 0; i < K; ++i) {
        std::vector<double> frameSamples(nInFrame);
        for (auto j = 0; j < nInFrame; ++j) {
            frameSamples[j] = samples[i * nInFrame + j];
        }

        std::vector<double> fftResult(nInFrame);
        FFTAnalysis(frameSamples.data(), fftResult.data(), nInFrame, nInFrame);

        std::vector<ShortComplex> frameSamples2(nInFrame);
        for (auto j = 0; j < nInFrame; ++j) {
            frameSamples2[j].re = samples[i * nInFrame + j];
            frameSamples2[j].im = 0.0;
        }

        fft(frameSamples2.data(), s, false);

        double power = 0.0;
        for (auto j = 0; j < nInFrame; ++j) {
            power += fftResult[j] * fftResult[j];
        }
        framePowers[i] = 0.5 * power;
    }

    return result;
}

std::vector<short> Filter::convert(const std::vector<unsigned char> &samples) {
    const auto sampleSize = sizeof(unsigned short);
    const auto samplesSize = samples.size();
    std::vector<short> result(samplesSize / sampleSize);

    for (auto i = 0, j = 0; i < samplesSize; i += sampleSize) {
        short current = samples[i + sampleSize - 1];
        for (auto k = i + 1; k < i + sampleSize; ++k) {
            current <<= 8;
            current |= samples[k];
        }

        /*unsigned short current = samples[i + 1];
        current <<= 8;
        current |= samples[i];*/
        result[j++] = current;
    }

    return result;
}

// AVal - массив анализируемых данных, Nvl - длина массива должна быть кратна степени 2.
// FTvl - массив полученных значений, Nft - длина массива должна быть равна Nvl.
void Filter::FFTAnalysis(const double *AVal, double *FTvl, int Nvl, int Nft) {
    int i, j, n, m, Mmax, Istp;
    double Tmpr, Tmpi, Wtmp, Theta;
    double Wpr, Wpi, Wr, Wi;
    double *Tmvl;

    n = Nvl * 2;
    Tmvl = new double[n];

    for (i = 0; i < n; i += 2) {
        Tmvl[i] = 0;
        Tmvl[i + 1] = AVal[i / 2];
    }

    i = 1;
    j = 1;
    while (i < n) {
        if (j > i) {
            Tmpr = Tmvl[i];
            Tmvl[i] = Tmvl[j];
            Tmvl[j] = Tmpr;
            Tmpr = Tmvl[i + 1];
            Tmvl[i + 1] = Tmvl[j + 1];
            Tmvl[j + 1] = Tmpr;
        }
        i = i + 2;
        m = Nvl;
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
                Tmpr = Wr * Tmvl[j] - Wi * Tmvl[j - 1];
                Tmpi = Wi * Tmvl[j] + Wr * Tmvl[j - 1];

                Tmvl[j] = Tmvl[i] - Tmpr;
                Tmvl[j - 1] = Tmvl[i - 1] - Tmpi;
                Tmvl[i] = Tmvl[i] + Tmpr;
                Tmvl[i - 1] = Tmvl[i - 1] + Tmpi;
                i = i + Istp;
            }
        }

        Mmax = Istp;
    }

    for (i = 0; i < Nft; i++) {
        j = i * 2;
        FTvl[i] = 2 * sqrt(pow(Tmvl[j], 2) + pow(Tmvl[j + 1], 2)) / Nvl;
    }

    delete[]Tmvl;
}
