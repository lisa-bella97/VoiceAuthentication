#pragma once

#include <vector>


class Filter {
public:
    static std::vector<unsigned char> removePauses(const std::vector<unsigned char> &rawSamples, int sampleRate);

    static std::vector<std::vector<double>>
    getMFCC(const std::vector<unsigned char> &rawSamples, int sampleRate, int mfccNumber);

private:
    static std::vector<short> convertFromRaw(const std::vector<unsigned char> &rawSamples);

    static std::vector<unsigned char> convertToRaw(const std::vector<short> &samples);

    /**
     * Вычисление модуля спектра действительного массива чисел на основе реализации быстрого преобразования Фурье
     * @param data - массив анализируемых данных, длина массива должна быть кратна степени 2
     */
    static std::vector<double> FFTAnalysis(const std::vector<double> &data);

    static std::vector<std::vector<std::pair<double, int>>> kMeans(const std::vector<std::pair<double, int>> &values);

    static double getMeanValue(const std::vector<std::pair<double, int>> &values);

    static void DCT(int direction, int length, double X[]);

    static void FFT(int direction, int length, double Xr[], double Xi[]);

    static double melScale(int direction, double x);

    static void MFCC(int length_frame, int length_DFT, int number_coefficients, int number_filterbanks, int sample_rate,
                     const double frame[], double feature_vector[]);
};
