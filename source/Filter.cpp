#include <cmath>
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

std::vector<unsigned char> Filter::removePauses(const std::vector<unsigned char> &samples, int sampleSize,
                                                int sampleRate) {
    std::vector<unsigned char> result;

    // Определение окна для Фурье-преобразования
    int N = samples.size() / sampleSize; // Число отсчетов сигнала
    double T = 1.0 / sampleRate; // Период дискретизации
    double n = 0.02 / T; // Количество отсчетов в периоде стационарности
    double s = log(n) / log(2); // Степень двойки (log2(n) = log(n) / log(2))
    int nInFrame = (int)pow(2, s); // Количество отсчетов в кадре
    int K = N / nInFrame; // Количество итераций Фурье-преобразования

    // Вычисление Фурье-преобразования



    return result;
}

unsigned long long Filter::getAmplitude(const std::vector<unsigned char> &samples, int sampleSize) {
    auto samplesSize = samples.size();
    unsigned long long max = 0;

    for (auto i = 0; i < samplesSize; i += sampleSize) {
        // TODO: сейчас это реализовано только для sampleSize = 2
        unsigned long long current = samples[i];
        current <<= 8;
        current |= samples[i + 1];

        if (current > max)
            max = current;
    }

    return max;
}
