#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <math.h>
#include <stdlib.h>

#include "catalina.h"
//#include "alg1_psyho.h"
//#include "alg2_res1.h"
//#include "alg3_pfr.h"
//#include "alg4_nofto.h"
//#include "alg5_kushal1.h"

using namespace std;

struct SDetection
{
    double RA[4];
    double DEC[4];
    string imgID;
};


int main(int argc, char *argv[])
{
    char buf[65536];

    AsteroidDetector my;
    for (int i=0;i<100;i++)
    {
        int W,H,N;
        cin >> W >> H;
        vector<int> image[4];
        vector<string> header[4];
        vector<double> wcs[4];

        for (int f=0;f<4;f++)
        {
            image[f].resize(W*H);
            for (int i=0;i<W*H;i++)
                cin >> image[f][i];
            cin >> N;
            cerr << "N=" << N << endl;
            cin.getline(buf, 65536);
            for (int i=0;i<N;i++)
            {
                cin.getline(buf, 65536);
                header[f].push_back(buf);
            }
            wcs[f].resize(8);
            for (int i=0;i<8;i++)
            {
                cin >> wcs[f][i];
            }
        }
        vector<string> detections;
        cin >> N;
        cerr << "N=" << N << endl;
        cin.getline(buf, 65536);
        for (int i=0;i<N;i++)
        {
            cin.getline(buf, 65536);
            detections.push_back(buf);
        }
      
        int v = my.trainingData(W, H, image[0], header[0], wcs[0], image[1], header[1], wcs[1], image[2], header[2], wcs[2], image[3], header[3], wcs[3], detections);
        cout << v << endl;
        cout.flush();
        if (v==1) break;
    }

    for (int i=0;i<20;i++)
    {
        int W,H,N;
        string ID;
        cin.getline(buf, 65536);
        ID = buf;
        cin >> W >> H;
        vector<int> image[4];
        vector<string> header[4];
        vector<double> wcs[4];

        cerr << "S : W=" << W << " H=" << H << endl;
        for (int f=0;f<4;f++)
        {
            image[f].resize(W*H);
            for (int i=0;i<W*H;i++)
                cin >> image[f][i];
            cin >> N;
            cerr << "S : N=" << N << endl;
            cin.getline(buf, 65536);
            for (int i=0;i<N;i++)
            {
                cin.getline(buf, 65536);
                header[f].push_back(buf);
            }
            wcs[f].resize(8);
            for (int i=0;i<8;i++)
            {
                cin >> wcs[f][i];
            }
        }
        cin.getline(buf, 65536);
        int v = my.testingData(ID, W, H, image[0], header[0], wcs[0],
                                         image[1], header[1], wcs[1],
                                         image[2], header[2], wcs[2],
                                         image[3], header[3], wcs[3]);
        cout << v << endl;
        cout.flush();
    }

    vector<string> res = my.getAnswer();
    cout << res.size() << endl;
    for (int i=0; i < res.size(); i++)
        cout << res[i] << endl;
    cout.flush();


    return 0;
}
