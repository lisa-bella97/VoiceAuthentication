#include <OpenALRecorder.h>

#include <cstring>
#include <iostream>
#include <Filter.h>


int main(int argc, char **argv) {
    std::string deviceName;
    int channels = 1;
    int bits = 16;
    int sampleRate = 9600;
    int time = 4;
    std::string fileName;

    const std::string helpText = "Usage:\n" + std::string(argv[0]) + "\n"
                                                                     "[--device <device name>]\n"
                                                                     "[--channels <1|2>]\n"
                                                                     "[--bits <8|16|32>]\n"
                                                                     "[--rate <record rate in Hz>]\n"
                                                                     "[--time <record time in seconds>]\n"
                                                                     "[--file <file name (*.wav)>]\n\n"
                                                                     "Default options:\n"
                                                                     "--device default_device\n"
                                                                     "--channels 1\n"
                                                                     "--bits 16\n"
                                                                     "--rate 9600\n"
                                                                     "--time 5\n"
                                                                     "--file no_file\n";

    for (auto i = 1; i < argc; i += 2) {
        if (std::strcmp(argv[i], "--device") == 0)
            deviceName = argv[i + 1];
        else if (std::strcmp(argv[i], "--channels") == 0)
            channels = std::atoi(argv[i + 1]);
        else if (std::strcmp(argv[i], "--bits") == 0)
            bits = std::atoi(argv[i + 1]);
        else if (std::strcmp(argv[i], "--rate") == 0)
            sampleRate = std::atoi(argv[i + 1]);
        else if (std::strcmp(argv[i], "--time") == 0)
            time = std::atoi(argv[i + 1]);
        else if (std::strcmp(argv[i], "--file") == 0)
            fileName = std::string(argv[i + 1]);
        else if (std::strcmp(argv[i], "--help") == 0 || std::strcmp(argv[i], "-h") == 0) {
            std::cout << helpText;
            return 0;
        } else {
            std::cout << "Wrong args." << std::endl << helpText;
            return -1;
        }
    }

    try {
        OpenALRecorder recorder(deviceName, channels, bits, sampleRate);
        std::cout << "Opened device: " << recorder.getCaptureDeviceName() << std::endl;

        if (fileName.empty()) {
            auto result = recorder.record(time + 1);
            std::cout << "Recording is completed." << std::endl;
            recorder.writeToFile(result, "sample.wav");
            result = Filter::removePauses(result, sampleRate);
            recorder.writeToFile(result, "pausesDeleted.wav");
            auto res = Filter::getMFCC(result, sampleRate, 13);
            for (auto i = 0; i < res.size(); ++i) {
                for (auto j = 0; j < res[0].size(); ++j) {
                    std::cout << res[i][j] << ' ';
                }
                std::cout << std::endl;
            }
        } else {
            recorder.record(time, fileName);
            std::cout << "Recording is completed to file: " << fileName << std::endl;
        }
    } catch (const std::exception &e) {
        std::cout << "An exception occurred: " << e.what() << std::endl;
        return -2;
    }

    return 0;
}
