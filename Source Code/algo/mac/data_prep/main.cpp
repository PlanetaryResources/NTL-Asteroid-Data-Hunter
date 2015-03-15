#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <math.h>
#include <stdlib.h>
#include <map>

using namespace std;

vector<string> tokenize(const string& str2, char delimiters)
{
    vector<string> tokens;

    int len = (int)str2.length();
    stringstream ss;
    for (int i=0;i<len;i++)
    {
        if (str2[i]!=delimiters)
        {
            ss << str2[i];
        } else
        {
            if (i+1<len && str2[i+1]!=delimiters)
            {
                ss << str2[i];
            }
        }
    }
    string str = ss.str();
  //  cerr << str << endl;
    len = str.length();

    int pos = 0;
    int i = 0;
    int quote = 0;
    for (;;)
    {
        while (pos<len && (str[pos]!=delimiters || (str[pos]==delimiters && quote!=0)))
        {
            if (str[pos]=='"')
                quote = 1-quote;
            pos++;
        }
        if (pos<len && str[pos]==delimiters)
        {
            if (str[i]=='"')
            {
                tokens.push_back(str.substr(i+1, pos-i-2));
            }
            else
            {
                tokens.push_back(str.substr(i, pos-i));
            }
            i = pos+1;
        }
        else
        {
            if (str[i]=='"')
            {
                tokens.push_back(str.substr(i+1, pos-i-2));
            }
            else
            {
                tokens.push_back(str.substr(i, pos-i));
            }
            break;
        }
        pos++;
    }
    return tokens;
}

void checkb()
{
    char buf[65536];
    vector<string> list;
    ifstream fi("../dt/mpcd.txt");
    while (!fi.eof())
    {
        fi.getline(buf, 65536);
        string s = buf;
        if (s.length()==0) break;
        list.push_back(buf);
    }
    fi.close();
    int cnt1 = 0;
    int cnt2 = 0;
    int cnt3 = 0;
    int tot = 0;
    for (int i=0;i<list.size();i++)
    {
        stringstream ss;
        ss << "../dt/" << list[i];
        fi.open(ss.str().c_str());
        while (!fi.eof())
        {
            fi.getline(buf, 65536);
            string s = buf;
            if (s.length()==0) break;
            vector<string> tok = tokenize(s, ' ');
            if (tok.size()<10) continue;
            float deg1 = atof(tok[7].c_str());
            float deg2 = atof(tok[8].c_str());
            float deg3 = atof(tok[9].c_str());

            if (deg1<0)
            {
                cerr << deg1 << " " << deg2 << " " << deg3 << " : ";
                float v =  deg1 + deg2 / 60.0 + deg3 / 3600.0;
                cerr << v << endl;
                cnt1++;
            }
            if (deg2<0) cnt2++;
            if (deg3<0) cnt3++;
            tot++;

        }
        fi.close();
    }
    cerr << cnt1 << " " << cnt2 << " " << cnt3 << endl;
    cerr << "tot " << tot << endl;
}

int main(int argc, char *argv[])
{
//    checkb();
 //   return 0;

    srand(123579L);
    char buf[65536];
    vector<string> list;
    map<string, int> allcheck;
    ifstream fi("list.txt");
    while (!fi.eof())
    {
        fi.getline(buf, 65536);
        string s = buf;
        if (s.length()==0) break;
        list.push_back(buf);
        allcheck[s] = 1;
    }
    fi.close();

    map<string, string> basedir;

    map<string, int> alldt;

    vector<string> ephm;
    map<string, int> mephm;
    for (int i=0;i<list.size();i++)
    {
        string sub = list[i].substr(list[i].length()-4, 4);
        if (sub=="ephm")
        {
            string id = list[i].substr(list[i].length()-4-23, 17);
            cerr << id << endl;
            string b = list[i].substr(0, list[i].length()-4-23);
            cerr << b << endl;
            basedir[id] = b;
            mephm[id] = 1;
            ephm.push_back(id);
          //  cerr << list[i] << endl;
        }
    }

    vector<string> neos;
    for (int i=0;i<list.size();i++)
    {
        string sub = list[i].substr(list[i].length()-4, 4);
        if (sub=="neos")
        {
            string id = list[i].substr(list[i].length()-4-23, 17);
            cerr << id << endl;
            string b = list[i].substr(0, list[i].length()-4-23);
            cerr << b << endl;
            if (mephm.count(id)==0) continue;
            neos.push_back(list[i]);
           // cerr << list[i] << endl;
            basedir[id] = b;
            alldt[id] = 1;
        }
    }

    while (alldt.size()<300)
    {
        int i = rand()%ephm.size();
        alldt[ephm[i]] = 1;
    }
    cerr << "alldt=" << alldt.size() << endl;
    cerr << "ephm=" << ephm.size() << endl;
    cerr << "neos=" << neos.size() << endl;
    cerr << "total=" << list.size() << endl;

    vector<string> valldt;
    for (map<string, int>::iterator itr=alldt.begin();itr!=alldt.end();++itr)
    {
        valldt.push_back(itr->first);
    }
    int used[300] = {0};
    for (int i=0;i<100;i++)
    {
        int i1;
        do
        {
            i1 = rand()%300;
        } while (used[i1]);
        used[i1] = 1;
    }
    // training set
    vector<string> trainID;
    vector<string> testID;
    for (int i=0;i<300;i++)
    {
        if (used[i])
            trainID.push_back(valldt[i]);
        else
            testID.push_back(valldt[i]);
    }
    {
        ofstream fo("traindata.txt");
        for (int i=0;i<trainID.size();i++)
            fo << trainID[i] << endl;
        fo.close();
    }
    {
        ofstream fo("testdata.txt");
        for (int i=0;i<testID.size();i++)
            fo << testID[i] << endl;
        fo.close();
    }

    // download script
    {
        cerr << "TRAINING SCRIPT\n";
        ofstream fo("wget_traindata.sh");
        for (int i=0;i<trainID.size();i++)
        {
            string ID = trainID[i];
            string bdir = basedir[ID];
            for (int f=1;f<=4;f++)
            {
                stringstream ss;
                ss << bdir << ID << "_000" << f << ".arch.H";
                if (allcheck.count(ss.str())==0) cerr << "!!!" << ss.str() << endl;
                fo << "wget " << ss.str() << endl;
            }
            {
                stringstream ss;
                ss << bdir << ID << "_0001.ephm";
                if (allcheck.count(ss.str())==0) {}//cerr << "!!!" << ss.str() << endl;
                else fo << "wget " << ss.str() << endl;
            }
            {
                stringstream ss;
                ss << bdir << ID << "_0001.neos";
                if (allcheck.count(ss.str())==0) {}//cerr << "!!!" << ss.str() << endl;
                else fo << "wget " << ss.str() << endl;
            }
            {
                stringstream ss;
                ss << bdir << ID << "_0001.mpcd";
                if (allcheck.count(ss.str())==0) {}//cerr << "!!!" << ss.str() << endl;
                else fo << "wget " << ss.str() << endl;
            }
            {
                stringstream ss;
                ss << bdir << ID << "_0001.dets";
                if (allcheck.count(ss.str())==0) {}//cerr << "!!!" << ss.str() << endl;
                else fo << "wget " << ss.str() << endl;
            }
        }
        fo.close();

    }
    {
        cerr << "TEST SCRIPT\n";
        ofstream fo("wget_testdata.sh");
        for (int i=0;i<testID.size();i++)
        {
            string ID = testID[i];
            string bdir = basedir[ID];
            for (int f=1;f<=4;f++)
            {
                stringstream ss;
                ss << bdir << ID << "_000" << f << ".arch.H";
                if (allcheck.count(ss.str())==0) cerr << "!!!" << ss.str() << endl;
                fo << "wget " << ss.str() << endl;
            }
            {
                stringstream ss;
                ss << bdir << ID << "_0001.ephm";
                if (allcheck.count(ss.str())==0) {}//cerr << "!!!" << ss.str() << endl;
                else fo << "wget " << ss.str() << endl;
            }
            {
                stringstream ss;
                ss << bdir << ID << "_0001.neos";
                if (allcheck.count(ss.str())==0) {}//cerr << "!!!" << ss.str() << endl;
                else fo << "wget " << ss.str() << endl;
            }
            {
                stringstream ss;
                ss << bdir << ID << "_0001.mpcd";
                if (allcheck.count(ss.str())==0) {}//cerr << "!!!" << ss.str() << endl;
                else fo << "wget " << ss.str() << endl;
            }
            {
                stringstream ss;
                ss << bdir << ID << "_0001.dets";
                if (allcheck.count(ss.str())==0) {}//cerr << "!!!" << ss.str() << endl;
                else fo << "wget " << ss.str() << endl;
            }
        }
        fo.close();

    }


    return 0;
}
