#pragma once

#include <vector>


class Filter {
public:
    static std::vector<unsigned char> normalizeVolume(const std::vector<unsigned char> &samples, int sampleSize);

    static std::vector<unsigned char> removePauses(const std::vector<unsigned char> &samples, int sampleSize,
                                                   int sampleRate);

private:
    static unsigned long long getAmplitude(const std::vector<unsigned char> &samples, int sampleSize);

};
