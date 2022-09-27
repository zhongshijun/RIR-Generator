#include "rir_generator_core.h"
#include <iostream>
#include <algorithm>
#include <cstdio>
#include <string>
#include <vector>
#include <fstream>
#include <ctime>
#include <istream>
#include <ostream>

using namespace std;

struct RIR_Params {
    double soundVelocity;
    double sampleRate;

    int micNum;
    double *receiverPosition;
    double *sourcePosition;
    double *roomDimensions;  // room w x h x d 房间长宽高

    double *orientation; // 麦克风整列角度
    
    double rt60;  // reverberation times(s)
    //int n;        // number of samples

    bool hp_filter;
    char microphoneType; 
    int reflectionOrder;

    int roomDimension;
};


int RirGenegatorFunction(double *rir_out, RIR_Params &prhs)
{
    /*if (nrhs == 0)
    {
        printf("--------------------------------------------------------------------\n"
            "| Room Impulse Response Generator                                  |\n"
            "|                                                                  |\n"
            "| Computes the response of an acoustic source to one or more       |\n"
            "| microphones in a reverberant room using the image method [1,2].  |\n"
            "|                                                                  |\n"
            "| Author    : dr.ir. Emanuel Habets (e.habets@ieee.org)            |\n"
            "|                                                                  |\n"
            "| Version   : 2.2.20201022                                         |\n"
            "|                                                                  |\n"
            "| Copyright (C) 2003-2020 E.A.P. Habets                            |\n"
            "|                                                                  |\n"
            "| [1] J.B. Allen and D.A. Berkley,                                 |\n"
            "|     Image method for efficiently simulating small-room acoustics,|\n"
            "|     Journal Acoustic Society of America,                         |\n"
            "|     65(4), April 1979, p 943.                                    |\n"
            "|                                                                  |\n"
            "| [2] P.M. Peterson,                                               |\n"
            "|     Simulating the response of multiple microphones to a single  |\n"
            "|     acoustic source in a reverberant room, Journal Acoustic      |\n"
            "|     Society of America, 80(5), November 1986.                    |\n"
            "--------------------------------------------------------------------\n\n"
            "function [h, beta_hat] = rir_generator(c, fs, r, s, L, beta, nsample,\n"
            " mtype, order, dim, orientation, hp_filter);\n\n"
            "Input parameters:\n"
            " c           : sound velocity in m/s.\n"
            " fs          : sampling frequency in Hz.\n"
            " r           : M x 3 array specifying the (x,y,z) coordinates of the\n"
            "               receiver(s) in m.\n"
            " s           : 1 x 3 vector specifying the (x,y,z) coordinates of the\n"
            "               source in m.\n"
            " L           : 1 x 3 vector specifying the room dimensions (x,y,z) in m.\n"
            " beta        : 1 x 6 vector specifying the reflection coefficients\n"
            "               [beta_x1 beta_x2 beta_y1 beta_y2 beta_z1 beta_z2] or\n"
            "               beta = reverberation time (T_60) in seconds.\n"
            " nsample     : number of samples to calculate, default is T_60*fs.\n"
            " mtype       : [omnidirectional, subcardioid, cardioid, hypercardioid,\n"
            "               bidirectional], default is omnidirectional.\n"
            " order       : reflection order, default is -1, i.e. maximum order.\n"
            " dim         : room dimension (2 or 3), default is 3.\n"
            " orientation : direction in which the microphones are pointed, specified using\n"
            "               azimuth and elevation angles (in radians), default is [0 0].\n"
            " hp_filter   : use 'false' to disable high-pass filter, the high-pass filter\n"
            "               is enabled by default.\n\n"
            "Output parameters:\n"
            " h           : M x nsample matrix containing the calculated room impulse\n"
            "               response(s).\n"
            " beta_hat    : In case a reverberation time is specified as an input parameter\n"
            "               the corresponding reflection coefficient is returned.\n\n");
        return;
    } else
    {
        printf("Room Impulse Response Generator (Version 2.2.20201022) by Emanuel Habets\n"
            "Copyright (C) 2003-2020 E.A.P. Habets\n");
    }*/


    // Load parameters
    double          c = prhs.soundVelocity;
    double          fs = prhs.sampleRate;
    int             nMicrophones = prhs.micNum;
    double*         rr = prhs.receiverPosition;
    double*         ss = prhs.sourcePosition;
    double*         LL = prhs.roomDimensions;
    double          beta_input = prhs.rt60;
    double          beta[6];
    int             nSamples;
    char            microphone_type;
    int             nOrder;
    int             nDimension;
    double          microphone_angle[2];
    int             isHighPassFilter;
    double          reverberation_time = 0;

    // Reflection coefficients or reverberation time?
    double V = LL[0] * LL[1] * LL[2];
    double S = 2 * (LL[0] * LL[2] + LL[1] * LL[2] + LL[0] * LL[1]);
    reverberation_time = beta_input;
    if (reverberation_time != 0) {
        double alfa = 24 * V*log(10.0) / (c*S*reverberation_time);
        if (alfa > 1)
            printf("Error: The reflection coefficients cannot be calculated using the current "
                "room parameters, i.e. room size and reverberation time.\n           Please "
                "specify the reflection coefficients or change the room parameters.");
        for (int i = 0; i < 6; i++)
            beta[i] = sqrt(1 - alfa);
    } else
    {
        for (int i = 0; i < 6; i++)
            beta[i] = 0;
    }

    // High-pass filter (optional)
    if (prhs.hp_filter == false)
    {
        isHighPassFilter = (int)prhs.hp_filter;
    } else
    {
        isHighPassFilter = 1;
    }

    // 3D Microphone orientation (optional)
    if (nMicrophones > 1)
    {
        double* orientation = prhs.orientation;
        if (nMicrophones == 2)
        {
            microphone_angle[0] = orientation[0];
            microphone_angle[1] = 0;
        } else
        {
            microphone_angle[0] = orientation[0];
            microphone_angle[1] = orientation[1];
        }
    } else
    {
        microphone_angle[0] = 0;
        microphone_angle[1] = 0;
    }

    // Room Dimension (optional)
    if (prhs.roomDimension > 0)
    {
        nDimension = (int)prhs.roomDimension;
        if (nDimension != 2 && nDimension != 3)
            printf("Invalid input arguments!");

        if (nDimension == 2)
        {
            beta[4] = 0;
            beta[5] = 0;
        }
    } else
    {
        nDimension = 3;
    }

    // Reflection order (optional)
    nOrder = prhs.reflectionOrder;

    // Type of microphone (optional)
    microphone_type = prhs.microphoneType;

    // Number of samples (optional)
    if (reverberation_time < 0.128)
        reverberation_time = 0.128;
    nSamples = (int)(reverberation_time * fs) + 4800;
    //nSamples = nSamples > 48000 ? 48000 : nSamples;
    //nSamples = 48000; // 1s

    printf("c: %f, fs: %f,micNum: %d, microphone_type: %c, nOrder: %d, isHighPassFilter: %d\n", c, fs, nMicrophones, microphone_type, nOrder, isHighPassFilter);
    computeRIR(rir_out, c, fs, rr, nMicrophones, nSamples, ss, LL, beta, microphone_type, nOrder, microphone_angle, isHighPassFilter);

    return nSamples;
}

double rir_out[500000];

int main(int argc, char* argv[])
{
    printf("Generator rir start!\n");

    string outRirPath = "D:\\HMScloud\\RIR-Generator-master\\RIR-Generator-master\\rir\\";

    int rirNum = 20000;
    // 房间大小
    int roomHeightMax = 4; // m
    int roomWidthMax = 5;
    int roomLengthMax = 5;
    int RT60Max = 2; // 0.35s
    for (int i = 0; i < rirNum; ++i) {
        int W = rand() % roomWidthMax + 1;
        int L = rand() % roomLengthMax + 1;
        int H = rand() % roomHeightMax + 1;
        H = max(H, 2); // 房间高度最低2m

        int rW = rand() % W;
        int rL = rand() % L;
        int rH = rand() % H;

        int sW = rand() % W;
        int sL = rand() % L;
        int sH = rand() % H;

        double roomDimensions[5] = { 5, 4, 4 };
        double receiverPosition[5] = { 2, 1.5, 2 };
        double sourcePosition[5] = { 2, 3.5, 2 };

        roomDimensions[0] = W + (rand() % 10 / 10.f);
        roomDimensions[1] = L + (rand() % 10 / 10.f);
        roomDimensions[2] = H + (rand() % 10 / 10.f);

        sourcePosition[0] = sW + (rand() % 10 / 10.f);
        sourcePosition[1] = sL + (rand() % 10 / 10.f);
        sourcePosition[2] =sH + (rand() % 10 / 10.f);

        receiverPosition[0] = rW + (rand() % 10 / 10.f);
        receiverPosition[1] = rL + (rand() % 10 / 10.f);
        receiverPosition[2] = rH + (rand() % 10 / 10.f);

        sourcePosition[2] = min(sourcePosition[2], 2.0);
        sourcePosition[2] = max(sourcePosition[2], 0.3);

        receiverPosition[2] = min(receiverPosition[2], 2.0);
        receiverPosition[2] = max(receiverPosition[2], 0.3);

        RIR_Params rirParams;

        rirParams.soundVelocity = 340;
        rirParams.sampleRate = 48000;
        rirParams.rt60 = (double)(rand() % RT60Max / 10.0 + rand() % 5 / 100.0 + 0.15); // MAX RT60=> 0.1~0.4s
        //rirParams.rt60 = max(0.15, rirParams.rt60);

        rirParams.receiverPosition = receiverPosition;
        rirParams.sourcePosition = sourcePosition;
        rirParams.roomDimensions = roomDimensions;
        rirParams.roomDimension = 3;

        rirParams.hp_filter = 1;
        rirParams.micNum = 1;

        rirParams.reflectionOrder = -1;
        rirParams.microphoneType = 'o';

        int nSamples = RirGenegatorFunction(rir_out, rirParams);

        int rt60 = (int)(rirParams.rt60 * 1000); // ms
        //stringstream ss;
        //ss << rt60;
        string rirName = "RT60_" + std::to_string(rt60) + "ms_" + std::to_string(i) + ".pcm";
        string outFilePath = outRirPath + rirName;
        ofstream pcmFile;

        std::cout << "path=>" << outFilePath << std::endl;

        pcmFile.open(outFilePath, std::ifstream::out | std::ifstream::binary);

        short* out_pcm = new short[nSamples];
        for (int i = 0; i < nSamples; ++i) {
            int tmp = (int)(32768 * rir_out[i]);
            tmp = tmp > 32768 ? 32760 : tmp;
            tmp = tmp < -32760 ? 32760 : tmp;

            out_pcm[i] = (short)(tmp);


            //printf("%hd, ", out_pcm[i]);
            //if ((i + 1) % 10 == 0) {
            //    printf("\n");
            //}
        }
        //printf("\n");

        pcmFile.write((char *)out_pcm, nSamples * sizeof(short));

        delete[] out_pcm;
        pcmFile.close();

    }

    //double receiverPosition[5] = {2, 1.5, 2};
    //double sourcePosition[5] = {2, 3.5, 2};
    //double roomDimensions[5] = {5, 4, 6};

    //RIR_Params rirParams;
    //
    //rirParams.soundVelocity = 340; 
    //rirParams.sampleRate = 48000;
    //rirParams.rt60 = 0.5; // s

    //rirParams.receiverPosition = receiverPosition;
    //rirParams.sourcePosition = sourcePosition;
    //rirParams.roomDimensions = roomDimensions;
    //rirParams.roomDimension = 3;

    //rirParams.hp_filter = 1;
    //rirParams.micNum = 1;

    //rirParams.reflectionOrder = -1;
    //rirParams.microphoneType = 'o';

    //double rir_out[50000];
    //int nSamples = RirGenegatorFunction(rir_out, rirParams);

    //std::cout << "nSample: " << nSamples << std::endl;
    //
    ///*for (int i = 0; i < nSamples; ++i) {
    //    printf("%f, ", rir_out[i]);
    //    if ((i + 1) % 10 == 0) {
    //        printf("\n");
    //    }
    //}
    //std::cout << std::endl;*/


    //ofstream pcmFile;
    //pcmFile.open("D:\\HMScloud\\RIR-Generator-master\\RIR-Generator-master\\rir_out.pcm", std::ifstream::out | std::ifstream::binary);

    //short* out_pcm = new short[nSamples];
    //for (int i = 0; i < nSamples; ++i) {
    //    int tmp = (int)(32768 * rir_out[i]);
    //    tmp = tmp > 32768 ? 32760 : tmp;
    //    tmp = tmp < -32760 ? 32760 : tmp;

    //    out_pcm[i] = (short)(tmp);


    //    //printf("%hd, ", out_pcm[i]);
    //    //if ((i + 1) % 10 == 0) {
    //    //    printf("\n");
    //    //}
    //}
    ////printf("\n");

    //pcmFile.write((char *)out_pcm, nSamples * sizeof(short));

    //pcmFile.close();


    printf("Generator rir END!\n");

    return 0;
}
