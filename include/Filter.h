#pragma once

#include <vector>


class Filter {
public:
    static std::vector<unsigned char> normalizeVolume(const std::vector<unsigned char> &samples, int sampleSize);

    static std::vector<unsigned char> removePauses(const std::vector<unsigned char> &rawSamples, int sampleSize,
                                                   int sampleRate);

private:
    static std::vector<short> convert(const std::vector<unsigned char> &samples);

    /**
     * Вычисление модуля спектра действительного массива чисел на основе реализации быстрого преобразования Фурье
     * @param AVal - массив анализируемых данных
     * @param FTvl - длина массива; должна быть кратна степени 2
     * @param Nvl - массив полученных значений
     * @param Nft - длина массива; должна быть равна Nvl
     */
    static void FFTAnalysis(const double *AVal, double *FTvl, int Nvl, int Nft);
};
