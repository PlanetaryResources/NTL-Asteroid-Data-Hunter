#ifdef LOCAL
#include <fstream>
#endif
#include <ctime>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <set>
#define SQDIST(x1___,y1___,x2___,y2___)(((x1___)-(x2___))*((x1___)-(x2___))+((y1___)-(y2___))*((y1___)-(y2___)))
#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()
#define GCOST(x,n) ( 0.5*pow( 4.0*(x)*(1 - (x)),1.0/(n) ) )

#ifdef LOCAL
#include <time.h>
#else
#include <sys/time.h>
#endif

using namespace std;

#ifdef LOCAL
ofstream subor;
ofstream vysledok;
ofstream ficury;
#endif

const int PV=3;
const int CLUSTERS[PV]={2,2,2};
double ROOT=3.0;
const int FEATURES=32;
const int FEATURETRY=6;
const int MINNODE=1;
const int MAXLEVEL=50;
#ifdef LOCAL
const int TREES=50;
const double deadline=36000.0;
#else
const int TREES=200;
const double deadline=36000.0;
#endif
const double EPS = 1e-10;
const int MAXSAMPLESIZE=30000;

const bool BLURING=true;
const bool MEANING=false;
const int MEANGRID=128;
const int MEDIANGRID=64;
const int MAXOVERMEDIAN=100;        // max limit over median when calculating mean
const int SHAPEOVERMEANLIMIT=200;   // how much over mean pixel must be to be in shape
const int ADEPTGRID=50; //70;  // was 50             // grid size for adept lists
const double MINDIST=1e-3;          // minimum RD-distance of adept1 and adept4 to be candidate
const double MAXDIST=0.03; //0.06; // was 0.03;          // maximum RD-distance of adept1 and adept4 to be candidate
const double MAXXY=2;				// +- range of x, y when searching adept2 and adept3 
const double MAXRD=0.001;           // maximum RD-difference from the exact point
const double STARLIMIT=0.0005; //0.0007;      // maximum RD-distance to decide whether the object is moving
const double NORMRATIOLIMIT=0.4;    // maximum allowed ratio of the max and min norm among 3 segments
const double ANGLEDIFLIMIT=0.8;     // maximum allowed difference among angles of 3 segments

const int NEARLIMIT=8;

vector <double> VALUE[PV];
vector <double> REPRE[PV];
vector <int> featureScore(FEATURES,0);
vector <int> electionCount[PV];
vector <int> counterTree(PV,0);


const double MINDISTSQ=MINDIST*MINDIST;
const double MAXDISTSQ=MAXDIST*MAXDIST;
const double MAXRDSQ=MAXRD*MAXRD;
const double STARLIMITSQ=STARLIMIT*STARLIMIT;

const double PI=3.141592653589793238;
const double MATCH_DISTANCE = 0.001;


const int MARGIN=5;

double startTime;
int number=0;
int testNumber=0;
int W=4110;
int H=4096;
int WH=W*H;

int totalDet=0;
int totalNEO=0;
int foundDet=0;
int foundNEO=0;
int expectedNEOs=0;
int expectedTwos=0;
int totalTwo=0;
int foundTwo=0;

int WM=ceil(static_cast<double>(W-MARGIN)/MEDIANGRID);
int HM=ceil(static_cast<double>(H)/MEDIANGRID);
int WHM=WM*HM;

int WA=ceil(static_cast<double>(W)/ADEPTGRID);
int HA=ceil(static_cast<double>(H)/ADEPTGRID);
int WHA=WA*HA;

		 
double getTime() {
#ifdef LOCAL
    return (0.0+clock())/CLK_TCK;
#else
    timeval t;
    gettimeofday(&t,NULL);
    return 1e-6*t.tv_usec + t.tv_sec;
#endif
}

// using random generator from some solution of OctaveClassifier
unsigned long long nowRand = 1;
void seedBig(unsigned long long seed){
	nowRand = seed;
}
unsigned long long randBig(){
	nowRand = ((nowRand * 6364136223846793005ULL + 1442695040888963407ULL) >> 1);
	return nowRand;
}

vector <string> split(const string& text){   // split string by space
	vector <string> vys;
	stringstream ss(text);
    string word;
    while(getline(ss,word,' ')){
        vys.push_back(word);
    }
    return vys;
}
vector <string> splitBy(const string& text,char by){   // split string by by
	vector <string> vys;
	stringstream ss(text);
    string word;
    while(getline(ss,word,by)){
        vys.push_back(word);
    }
    return vys;
}
string trim(const string& str,const string& whitespace=" \t"){
    size_t strBegin = str.find_first_not_of(whitespace);
    if (strBegin == string::npos)
        return ""; // no content

    size_t strEnd = str.find_last_not_of(whitespace);
    size_t strRange = strEnd - strBegin + 1;
    return str.substr(strBegin, strRange);
}

double totalLS(const double& x1,const double& x2,const double& x3,const double& x4,const double& y1,const double& y2,const double& y3,const double& y4,double &slope){
	double qx=(x1+x2+x3+x4)/4;
	double qy=(y1+y2+y3+y4)/4;
	double M11=(x1-qx)*(x1-qx)+(x2-qx)*(x2-qx)+(x3-qx)*(x3-qx)+(x4-qx)*(x4-qx);
	double M22=(y1-qy)*(y1-qy)+(y2-qy)*(y2-qy)+(y3-qy)*(y3-qy)+(y4-qy)*(y4-qy);
	double M12=(x1-qx)*(y1-qy)+(x2-qx)*(y2-qy)+(x3-qx)*(y3-qy)+(x4-qx)*(y4-qy);
	double D=(M11-M22)*(M11-M22)+4*M12*M12;
	double lambda=(M11+M22-sqrt(D))/2;
	if(abs(M12)<1e-20) slope=PI/2;
	else slope=atan((M22-lambda)/M12);
	double a=lambda-M22;
	double b=M12;
	if(a==0 && b==0){
		a=1;b=1;
	}
	double c=a*qx+b*qy;
	return ((a*x1+b*y1-c)*(a*x1+b*y1-c)+(a*x2+b*y2-c)*(a*x2+b*y2-c)+(a*x3+b*y3-c)*(a*x3+b*y3-c)+(a*x4+b*y4-c)*(a*x4+b*y4-c))/(a*a+b*b);
}

vector <int> selectFeatures(){
	vector <int> result;
	result.reserve(FEATURETRY);
	set <int> temp;
	pair<set<int>::iterator,bool> ret;
	while(result.size()!=FEATURETRY){
		int f=rand()%FEATURES;
		while(f==7 || f==10 || f==18) f=rand()%FEATURES;
		ret=temp.insert(f);
		if(ret.second) result.push_back(f);
	}
	return result;
}

class EdgePoint{
  public:
	int x;
	int y;
	EdgePoint(int xx,int yy){
		x=xx;y=yy;
	}
	EdgePoint(){
		x=0;y=0;
	}
};
class Shape{
  public:
  	vector <EdgePoint> point;
  	//int scaleX;
  	//int scaleY;
  	//int maxData;
  	Shape(EdgePoint e){
	  	point.clear();point.push_back(e);
	}
  	Shape(){
	    point.clear();
	}
};
class ShapeInfo{
  public:
  	int scaleX;
  	int scaleY;
  	int scaleN;
  	ShapeInfo(){
  	}
  	ShapeInfo(int x,int y,int n){
  		scaleX=x;scaleY=y;scaleN=n;
  	}
};
vector <ShapeInfo> shapeDetection(const vector <int>& image,const vector <int>& mean,int limit,vector <int>& getShape){
	
	vector <int> isEdge(WH,0);
	vector <int> edgeStack;
	edgeStack.reserve(WH/10);
	
	for(int i=0;i<WH;i++){
		if(image[i]>mean[i]+limit){
			isEdge[i]=1;
			edgeStack.push_back(i);
		}
	}
	
	int shapeId=1;
	//vector <int> shapeMinX;
	//vector <int> shapeMaxX;
	//vector <int> shapeMinY;
	//vector <int> shapeMaxY;
	//vector <int> shapeMaxData;
	
	vector <Shape> shape;
	
	while(!edgeStack.empty()){
		int pixel=edgeStack.back();
		edgeStack.pop_back();
		int x=pixel%W;
		int y=pixel/W;
		if(isEdge[pixel]>1){
			continue;
		}
		if(isEdge[pixel]==1){
			shapeId++;
			//shapeMinX.push_back(x);shapeMaxX.push_back(x);shapeMinY.push_back(y);shapeMaxY.push_back(y);shapeMaxData.push_back(image[pixel]);
			isEdge[pixel]=shapeId;
			shape.push_back(Shape(EdgePoint(x,y)));
		}
		else{
			if(isEdge[pixel]<-1){
				isEdge[pixel]=-isEdge[pixel];
				shape[shapeId-2].point.push_back(EdgePoint(x,y));
				//shapeMinX[shapeId-2]=min(x,shapeMinX[shapeId-2]);
				//shapeMaxX[shapeId-2]=max(x,shapeMaxX[shapeId-2]);
				//shapeMinY[shapeId-2]=min(y,shapeMinY[shapeId-2]);
				//shapeMaxY[shapeId-2]=max(y,shapeMaxY[shapeId-2]);
				//shapeMaxData[shapeId-2]=max(image[pixel],shapeMaxData[shapeId-2]);
			}
		}
		vector <int> testing;
		testing.clear();
		if(x>0 && y>0){testing.push_back((y-1)*W+x-1);}
		if(x>0 && y<H-1){testing.push_back((y+1)*W+x-1);}
		if(x<W-1 && y>0){testing.push_back((y-1)*W+x+1);}
		if(x<W-1 && y<H-1){testing.push_back((y+1)*W+x+1);}
		if(x>0){testing.push_back(y*W+x-1);}
		if(x<W-1){testing.push_back(y*W+x+1);}
		if(y>0){testing.push_back((y-1)*W+x);}
		if(y<H-1){testing.push_back((y+1)*W+x);}
		
		while(!testing.empty()){
			int testPoint=testing.back();
			testing.pop_back();
			if(isEdge[testPoint]==-1 || isEdge[testPoint]==1){
				isEdge[testPoint]=-shapeId;
				edgeStack.push_back(testPoint);
			}
		}
	}
		
	// first update shape properties
	vector <ShapeInfo> result;
	getShape=vector <int> (WH,-1);
	for(int i=0;i<shape.size();i++){
		int maxX=-1;
		int maxY=-1;
		int minX=W;
		int minY=H;
		for(int j=0;j<shape[i].point.size();j++){
			if(shape[i].point[j].x<minX) minX=shape[i].point[j].x;
			if(shape[i].point[j].x>maxX) maxX=shape[i].point[j].x;
			if(shape[i].point[j].y<minY) minY=shape[i].point[j].y;
			if(shape[i].point[j].y>maxY) maxY=shape[i].point[j].y;
			getShape[shape[i].point[j].x+shape[i].point[j].y*W]=i;
		}
		result.push_back(ShapeInfo(maxX-minX,maxY-minY,shape[i].point.size()));
		//shape[i].scaleX=shapeMaxX[i]-shapeMinX[i];
		//shape[i].scaleY=shapeMaxY[i]-shapeMinY[i];
		//shape[i].maxData=shapeMaxData[i];
		//if(true) result.push_back(shape[i]);
	}
	return result;
}

class Wcs{
  public:
    double crpix1;
    double crpix2;
    double crval1;
    double crval2;
    double cd11;
    double cd12;
    double cd21;
    double cd22;
    double invcd11;
    double invcd12;
    double invcd21;
    double invcd22;
    double D2R;
    double R2D;
	Wcs(){
		D2R = PI/180.0;
		R2D = 180.0/PI;
	}
    void update_invcd(){
        double inv_det = cd11*cd22 - cd12*cd21;
        inv_det = 1.0 / inv_det;
        invcd11 = cd22 * inv_det;
        invcd12 = -cd12 * inv_det;
        invcd21 = -cd21 * inv_det;
        invcd22 = cd11 * inv_det;
    }
    void init(const vector <double>& wcs){
        crpix1 = wcs[0];
        crpix2 = wcs[1];
        crval1 = wcs[2];
        crval2 = wcs[3];
        cd11 = wcs[4];
        cd12 = wcs[5];
        cd21 = wcs[6];
        cd22 = wcs[7];
        update_invcd();
    }
    vector <double> sphs2x(const double& lng, const double& lat){
        double coslat, coslat3, coslat4, coslng, dlng, dphi, sinlat, sinlat3, sinlat4, sinlng, x, y, z;
        double eul[5];
        vector <double> phiTheta(2);
        eul[0] = crval1;
        eul[1] = 90.0 - crval2;
        eul[2] = 180.0;
        eul[3] = cos( D2R * eul[1] );
        eul[4] = sin( D2R * eul[1] );
        /* Do lng dependency. */
        dlng = lng - eul[0];
        phiTheta[0] = dlng;
        /* Do lat dependency. */
        sinlat = sin(lat*D2R);
        coslat = cos(lat*D2R);
        coslat3 = coslat*eul[3];
        coslat4 = coslat*eul[4];
        sinlat3 = sinlat*eul[3];
        sinlat4 = sinlat*eul[4];
        dlng = phiTheta[0];
        sinlng = sin(dlng*D2R);
        coslng = cos(dlng*D2R);
        /* Compute the native longitude. */
        x = sinlat4 - coslat3*coslng;
        y = -coslat*sinlng;
        dphi = R2D*atan2(y, x);
        phiTheta[0] = (eul[2] + dphi);
		while(phiTheta[0]<0) phiTheta[0]+=360.0;
		while(phiTheta[0]>360) phiTheta[0]-=360.0;
        /* Normalize the native longitude. */
        if (phiTheta[0] > 180.0) {
            phiTheta[0] -= 360.0;
        } else if (phiTheta[0] < -180.0) {
            phiTheta[0] += 360.0;
        }
        /* Compute the native latitude. */
        z = sinlat3 + coslat4*coslng;
        if (fabs(z) > 0.99) {
            if (z < 0.0)
                phiTheta[1] = -fabs(R2D * acos(sqrt(x*x+y*y)));
            else
                phiTheta[1]= fabs(R2D * acos(sqrt(x*x+y*y)));
        } else {
            phiTheta[1] = R2D * asin(z);
        }
        return phiTheta;
    }
	vector <double> sphx2s(const double& phi, const double& theta) {
            double cosphi, costhe, costhe3, costhe4, dlng, dphi, sinphi, sinthe, sinthe3, sinthe4, x, y, z;
            vector <double> lngLat(2);
            double eul[5];
            eul[0] = crval1;
            eul[1] = 90.0 - crval2;
            eul[2] = 180.0;
            eul[3] = cos( D2R * eul[1] );
            eul[4] = sin( D2R * eul[1] );
            /* Do phi dependency. */
            dphi = phi - eul[2];
            /* Do theta dependency. */
            sinthe = sin( theta * D2R );
            costhe = cos( theta * D2R );
            costhe3 = costhe*eul[3];
            costhe4 = costhe*eul[4];
            sinthe3 = sinthe*eul[3];
            sinthe4 = sinthe*eul[4];
            sinphi = sin( dphi * D2R );
            cosphi = cos( dphi * D2R );
            /* Compute the celestial longitude. */
            x = sinthe4 - costhe3*cosphi;
            y = -costhe*sinphi;
            dlng = R2D * atan2( y, x);
            lngLat[0] = eul[0] + dlng;
            /* Normalize the celestial longitude. */
            if (eul[0] >= 0.0) {
                if (lngLat[0] < 0.0)
                    lngLat[0] += 360.0;
            } else {
                if (lngLat[0] > 0.0)
                    lngLat[0] -= 360.0;
            }
            if (lngLat[0] > 360.0) {
                lngLat[0] -= 360.0;
            } else if (lngLat[0] < -360.0) {
                lngLat[0] += 360.0;
            }
            /* Compute the celestial latitude. */
            z = sinthe3 + costhe4*cosphi;
            if (fabs(z) > 0.99) {
            /* Use an alternative formula for greater accuracy. */
                if (z<0.0)
                    lngLat[1] = -fabs( R2D * acos(sqrt(x*x+y*y)) );
                else
                    lngLat[1] = fabs( R2D * acos(sqrt(x*x+y*y)) );
            } else {
                lngLat[1] = R2D * asin(z);
            }
            return lngLat;
    }
    vector <double> tans2x(const double& phi, const double& theta) {
            vector <double> xy(2);
            double cotan_theta = 1.0 / tan( D2R * theta );
            xy[0] = R2D * sin( D2R * phi ) * cotan_theta ;
            xy[1] = -R2D * cos( D2R * phi ) * cotan_theta;
            return xy;
    }
    vector <double> tanx2s(const double& x, const double& y) {
            vector <double> phiTheta(2);
            phiTheta[0] = R2D * atan2(x, -y);
            phiTheta[1] = R2D * atan2(R2D, sqrt(x*x + y*y) );
            return phiTheta;
    }
	void convertRADEC2XY(const vector <double>& wcs, const double& RA, const double& DEC, double& X,double& Y){
            init(wcs);
            vector <double> phiTheta = sphs2x(RA, DEC);
            vector <double> XXYY = tans2x(phiTheta[0], phiTheta[1]);
            X = invcd11*XXYY[0] + invcd12*XXYY[1] + crpix1;
            Y = invcd21*XXYY[0] + invcd22*XXYY[1] + crpix2;
    }
    void convertXY2RADEC(const vector <double>& wcs, const double& x, const double& y,double& RA,double& DEC){
            init(wcs);
            double dx = x-crpix1;
            double dy = y-crpix2;
            double xx = cd11*dx + cd12*dy;
            double yy = cd21*dx + cd22*dy;
            vector <double> phiTheta = tanx2s(xx, yy);
            vector <double> RD = sphx2s(phiTheta[0], phiTheta[1]);
            RA=RD[0];
            DEC=RD[1];
    }
};

/*void convertRADEC2XY(const vector <double>& WCS, const double& RA, const double& DEC,double &x,double &y){
    double dR = RA - WCS[2];
    double dD = DEC - WCS[3];
    double dY = (dR * WCS[6] - dD * WCS[4]) / (WCS[6]*WCS[5] - WCS[4]*WCS[7]);
    double dX = (dR - dY * WCS[5]) / WCS[4];
    x = dX + WCS[0];
    y = dY + WCS[1];
}
void convertXY2RADEC(const vector <double>& WCS, const double& X, const double& Y,double &RA,double & DEC){
    double dX = X - WCS[0];
    double dY = Y - WCS[1];
    RA = dX * WCS[4] + dY * WCS[5] + WCS[2];
    DEC = dX * WCS[6] + dY * WCS[7] + WCS[3];
}*/
Wcs wc;


vector <char> peaks(const vector <int>& image){
	vector <char> result(WH);
	for(int y=0;y<H;y++) for(int x=0;x<W;x++){
		int i=W*y+x;
		if(x==0 || y==0 || x==W-1 || y==H-1
			|| image[i+1]>image[i] || image[i-1]>image[i] || image[i+W]>image[i] || image[i-W]>image[i]){
			result[i]=0;
			continue;		
		}
		if(image[i+1+W]>image[i] || image[i-1+W]>image[i] || image[i+1-W]>image[i] || image[i-1-W]>image[i]){
			result[i]=1;
			continue;
		}
		if(x==1 || y==1 || x==W-2 || y==H-2
			|| image[i+2]>image[i+1] || image[i-2]>image[i-1] || image[i+2*W]>image[i+W] || image[i-2*W]>image[i-W]){
			result[i]=2;
			continue;		
		}
		if(    (image[i+W+2]>image[i+W+1] && image[i+W+2]>image[i+1]) || (image[i+2*W+1]>image[i+W+1] && image[i+2*W+1]>image[i+W])
			|| (image[i-W+2]>image[i-W+1] && image[i-W+2]>image[i+1]) || (image[i-2*W+1]>image[i-W+1] && image[i-2*W+1]>image[i-W])
			|| (image[i+W-2]>image[i+W-1] && image[i+W-2]>image[i-1]) || (image[i+2*W-1]>image[i+W-1] && image[i+2*W-1]>image[i+W])
			|| (image[i-W-2]>image[i-W-1] && image[i-W-2]>image[i-1]) || (image[i-2*W-1]>image[i-W-1] && image[i-2*W-1]>image[i-W])
			){
			result[i]=3;
			continue;		
		}
		if(image[i+2*W+2]>image[i+W+1] || image[i+2*W-2]>image[i+W-1] || image[i-2*W+2]>image[i-W+1] || image[i-2*W-2]>image[i-W-1]){
			result[i]=4;
			continue;
		}
		if(x==2 || y==2 || x==W-3 || y==H-3
			|| image[i+3]>image[i+2] || image[i-3]>image[i-2] || image[i+3*W]>image[i+2*W] || image[i-3*W]>image[i-2*W]){
			result[i]=5;
			continue;		
		}
		if(    (image[i+W+3]>image[i+W+2] && image[i+W+3]>image[i+2]) || (image[i+3*W+1]>image[i+2*W+1] && image[i+3*W+1]>image[i+2*W])
			|| (image[i-W+3]>image[i-W+2] && image[i-W+3]>image[i+2]) || (image[i-3*W+1]>image[i-2*W+1] && image[i-3*W+1]>image[i-2*W])
			|| (image[i+W-3]>image[i+W-2] && image[i+W-3]>image[i-2]) || (image[i+3*W-1]>image[i+2*W-1] && image[i+3*W-1]>image[i+2*W])
			|| (image[i-W-3]>image[i-W-2] && image[i-W-3]>image[i-2]) || (image[i-3*W-1]>image[i-2*W-1] && image[i-3*W-1]>image[i-2*W])
			){
			result[i]=6;
			continue;		
		}
		if(   (image[i+2*W+3]>image[i+2*W+2] && image[i+2*W+3]>image[i+W+2]) || (image[i+3*W+2]>image[i+2*W+2] && image[i+3*W+2]>image[i+2*W+1])
		   || (image[i-2*W+3]>image[i-2*W+2] && image[i-2*W+3]>image[i-W+2]) || (image[i-3*W+2]>image[i-2*W+2] && image[i-3*W+2]>image[i-2*W+1])
	       || (image[i+2*W-3]>image[i+2*W-2] && image[i+2*W-3]>image[i+W-2]) || (image[i+3*W-2]>image[i+2*W-2] && image[i+3*W-2]>image[i+2*W-1])
		   || (image[i-2*W-3]>image[i-2*W-2] && image[i-2*W-3]>image[i-W-2]) || (image[i-3*W-2]>image[i-2*W-2] && image[i-3*W-2]>image[i-2*W-1])
			){
			result[i]=7;
			continue;		
		}
		if(x==3 || y==3 || x==W-4 || y==H-4
			|| image[i+4]>image[i+3] || image[i-4]>image[i-3] || image[i+4*W]>image[i+3*W] || image[i-4*W]>image[i-3*W]){
			result[i]=8;
			continue;		
		}
		if(    (image[i+W+4]>image[i+W+3] && image[i+W+4]>image[i+3]) || (image[i+4*W+1]>image[i+3*W+1] && image[i+4*W+1]>image[i+3*W])
			|| (image[i-W+4]>image[i-W+3] && image[i-W+4]>image[i+3]) || (image[i-4*W+1]>image[i-3*W+1] && image[i-4*W+1]>image[i-3*W])
			|| (image[i+W-4]>image[i+W-3] && image[i+W-4]>image[i-3]) || (image[i+4*W-1]>image[i+3*W-1] && image[i+4*W-1]>image[i+3*W])
			|| (image[i-W-4]>image[i-W-3] && image[i-W-4]>image[i-3]) || (image[i-4*W-1]>image[i-3*W-1] && image[i-4*W-1]>image[i-3*W])
			){
			result[i]=9;
			continue;		
		}
		if(image[i+3*W+3]>image[i+2*W+2] || image[i+3*W-3]>image[i+2*W-2] || image[i-3*W+3]>image[i-2*W+2] || image[i-3*W-3]>image[i-2*W-2]){
			result[i]=10;
			continue;
		}
		result[i]=11;
		/*int i=1;
		bool go=true;
		while(go && i<10){
			if(x-i>=0 && x+i<W && y-i>=0 && y+i<H
				&& image[W*(y+i-1)+x]>=image[W*(y+i)+x] && image[W*(y-i+1)+x]>=image[W*(y-i)+x]
			 	&& image[W*(y)+x+i-1]>=image[W*(y)+x+i] && image[W*(y)+x-i+1]>=image[W*(y)+x-i]
				&& image[W*(y+i-1)+x+i-1]>=image[W*(y+i)+x+i] && image[W*(y-i+1)+x+i-1]>=image[W*(y-i)+x+i]
				&& image[W*(y+i-1)+x-i+1]>=image[W*(y+i)+x-i] && image[W*(y-i+1)+x-i+1]>=image[W*(y-i)+x-i]) i++;
			else go=false;
		}
		result[W*y+x]=i-1;*/
	}
	vector <int> counter(12,0);
	for(int i=0;i<WH;i++) counter[result[i]]++;
	for(int i=0;i<12;i++){
		cerr << i << ": " << counter[i] << ",   ";
		if(i==5) cerr << "\n";
	}
	cerr << "\n";	
	return result;
}

vector <int> blur(const vector <int>& image){
	vector <int> result(WH);
	for(int x=0;x<W;x++){
		result[x]=image[x];
		result[(H-1)*W+x]=image[(H-1)*W+x];
	}
	for(int y=1;y<H-1;y++){
		result[y*W]=image[y*W];
		result[y*W+W-1]=image[y*W+W-1];
	}
	for(int y=1;y<H-1;y++) for(int x=1;x<W-1;x++){
		result[y*W+x]=(2*image[y*W+x]+image[y*W+x+1]+image[y*W+x-1]+image[(y+1)*W+x]+image[(y-1)*W+x])/6;
	}
	return result;
}

vector <int> medians(const vector <int>& image){
	vector <int> median(WM*HM);
	vector <int> count(WM*HM);
	for(int y=0;y<HM;y++) for(int x=0;x<WM;x++){
		vector <int> data;
		data.reserve(MEDIANGRID*MEDIANGRID);
		for(int j=y*MEDIANGRID;j<(y+1)*MEDIANGRID && j<H;j++) for(int i=max(MARGIN,x*MEDIANGRID);i<(x+1)*MEDIANGRID && i<W-MARGIN;i++){
			data.push_back(image[i+W*j]);
		}
		nth_element(data.begin(), data.begin()+data.size()/2, data.end());
		count[y*WM+x]=data.size();
		median[y*WM+x]=data[data.size()/2];
	}
	vector <int> result(WM*HM);
	for(int y=0;y<HM;y++) for(int x=0;x<WM;x++){
		int counter=count[y*WM+x];
		result[y*WM+x]=count[y*WM+x]*median[y*WM+x];
		if(x>0){
			counter+=count[y*WM+x-1];
			result[y*WM+x]+=count[y*WM+x-1]*median[y*WM+x-1];
			if(y>0){
				counter+=count[(y-1)*WM+x-1];
				result[y*WM+x]+=count[(y-1)*WM+x-1]*median[(y-1)*WM+x-1];
			}
			if(y<HM-1){
				counter+=count[(y+1)*WM+x-1];
				result[y*WM+x]+=count[(y+1)*WM+x-1]*median[(y+1)*WM+x-1];
			}
		}
		if(x<WM-1){
			counter+=count[y*WM+x+1];
			result[y*WM+x]+=count[y*WM+x+1]*median[y*WM+x+1];
			if(y>0){
				counter+=count[(y-1)*WM+x+1];
				result[y*WM+x]+=count[(y-1)*WM+x+1]*median[(y-1)*WM+x+1];
			}
			if(y<HM-1){
				counter+=count[(y+1)*WM+x+1];
				result[y*WM+x]+=count[(y+1)*WM+x+1]*median[(y+1)*WM+x+1];
			}
		}
		if(y>0){
			counter+=count[(y-1)*WM+x];
			result[y*WM+x]+=count[(y-1)*WM+x]*median[(y-1)*WM+x];
		}
		if(y<HM-1){
			counter+=count[(y+1)*WM+x];
			result[y*WM+x]+=count[(y+1)*WM+x]*median[(y+1)*WM+x];
		}
		result[y*WM+x]/=counter;
	}
	return result;
}

vector <int> means(const vector <int>& image,int limit){
	int WR=W-2*MARGIN-MEANGRID+1;
	int HR=H-MEANGRID+1;
	vector <int> sum(H*WR);
	vector <int> lcount(H*WR);
	for(int i=0;i<H;i++){
		sum[i*WR]=0;
		lcount[i*WR]=0;
		for(int j=MARGIN;j<MARGIN+MEANGRID;j++){
			if(image[i*W+j]<=limit){
				sum[i*WR]+=image[i*W+j];
				lcount[i*WR]++;
			}
		}
		for(int j=1;j<WR;j++){
			sum[i*WR+j]=sum[i*WR+j-1];
			lcount[i*WR+j]=lcount[i*WR+j-1];
			if(image[i*W+MARGIN+j-1]<=limit){
				sum[i*WR+j]-=image[i*W+MARGIN+j-1];
				lcount[i*WR+j]--;
			}
			if(image[i*W+MARGIN+MEANGRID+j-1]<=limit){
				sum[i*WR+j]+=image[i*W+MARGIN+MEANGRID+j-1];
				lcount[i*WR+j]++;
			}
		}
	}
	vector <int> result(HR*WR);
	vector <int> count(HR*WR);
	for(int j=0;j<WR;j++){
		result[j]=0;
		count[j]=0;
		for(int i=0;i<MEANGRID;i++){
			result[j]+=sum[i*WR+j];
			count[j]+=lcount[i*WR+j];
		}
		for(int i=1;i<HR;i++){
			result[i*WR+j]=result[(i-1)*WR+j]-sum[(i-1)*WR+j]+sum[(i-1+MEANGRID)*WR+j];
			count[i*WR+j]=count[(i-1)*WR+j]-lcount[(i-1)*WR+j]+lcount[(i-1+MEANGRID)*WR+j];
		}
	}
	vector <int> mean(WH);
	for(int y=0;y<H;y++) for(int x=0;x<W;x++){
		int i=x;
		int j=y;
		if(x<MARGIN+MEANGRID/2) i=MARGIN+MEANGRID/2;
		if(x>=W-MARGIN-MEANGRID/2) i=W-MARGIN-MEANGRID/2-1;
		if(y<MEANGRID/2) j=MEANGRID/2;
		if(y>=H-MEANGRID/2) j=H-MEANGRID/2-1;
		int index=i-(MARGIN+MEANGRID/2)+(j-MEANGRID/2)*WR;
		if(count[index]==0) mean[x+y*W]=limit;
		else mean[x+y*W]=result[index]/count[index];
	}
	return mean;
}

vector <int> medians2(const vector <int>& image,const vector<int>& median){
	vector <int> result(WH);
	for(int y=0;y<H;y++) for(int x=0;x<W;x++){
		result[y*W+x]=median[WM*(y/MEDIANGRID)+x/MEDIANGRID];
	}
	return result;
}
int medianAt(int x,int y,const vector<int>& median){
	return median[WM*(y/MEDIANGRID)+x/MEDIANGRID];
}

vector <int> means2(const vector <int>& image,const vector <int>& median,int limit){
	if(MEANING){
		int WR=W-2*MARGIN-MEANGRID+1;
		int HR=H-MEANGRID+1;
		vector <int> sum(H*WR);
		vector <int> lcount(H*WR);
		for(int i=0;i<H;i++){
			sum[i*WR]=0;
			lcount[i*WR]=0;
			for(int j=MARGIN;j<MARGIN+MEANGRID;j++){
				if(image[i*W+j]<=median[WM*(i/MEDIANGRID)+j/MEDIANGRID]+limit){
					sum[i*WR]+=image[i*W+j];
					lcount[i*WR]++;
				}
			}
			for(int j=1;j<WR;j++){
				sum[i*WR+j]=sum[i*WR+j-1];
				lcount[i*WR+j]=lcount[i*WR+j-1];
				if(image[i*W+MARGIN+j-1]<=median[WM*(i/MEDIANGRID)+(MARGIN+j-1)/MEDIANGRID]+limit){
					sum[i*WR+j]-=image[i*W+MARGIN+j-1];
					lcount[i*WR+j]--;
				}
				if(image[i*W+MARGIN+MEANGRID+j-1]<=median[WM*(i/MEDIANGRID)+(MARGIN+MEANGRID+j-1)/MEDIANGRID]+limit){
					sum[i*WR+j]+=image[i*W+MARGIN+MEANGRID+j-1];
					lcount[i*WR+j]++;
				}
			}
		}
		vector <int> result(HR*WR);
		vector <int> count(HR*WR);
		for(int j=0;j<WR;j++){
			result[j]=0;
			count[j]=0;
			for(int i=0;i<MEANGRID;i++){
				result[j]+=sum[i*WR+j];
				count[j]+=lcount[i*WR+j];
			}
			for(int i=1;i<HR;i++){
				result[i*WR+j]=result[(i-1)*WR+j]-sum[(i-1)*WR+j]+sum[(i-1+MEANGRID)*WR+j];
				count[i*WR+j]=count[(i-1)*WR+j]-lcount[(i-1)*WR+j]+lcount[(i-1+MEANGRID)*WR+j];
			}
		}
		vector <int> mean(WH);
		for(int y=0;y<H;y++) for(int x=0;x<W;x++){
			int i=x;
			int j=y;
			if(x<MARGIN+MEANGRID/2) i=MARGIN+MEANGRID/2;
			if(x>=W-MARGIN-MEANGRID/2) i=W-MARGIN-MEANGRID/2-1;
			if(y<MEANGRID/2) j=MEANGRID/2;
			if(y>=H-MEANGRID/2) j=H-MEANGRID/2-1;
			int index=i-(MARGIN+MEANGRID/2)+(j-MEANGRID/2)*WR;
			if(count[index]==0) mean[x+y*W]=limit;
			else mean[x+y*W]=result[index]/count[index];
		}
		return mean;
	}
	else{
		vector <int> mean(WH);
		for(int y=0;y<H;y++) for(int x=0;x<W;x++) mean[x+y*W]=medianAt(x,y,median);
		return mean;
	}
}


void getRange(const vector <int>& imageData,int& dmin,int& dmax,int& median){
	vector <int> histo(65536,0);
	for(int i=0;i<WH;i++){
		int x = i%W;
		if(x<MARGIN || x>=W-MARGIN) continue;
		int p = max(0,min(65535,imageData[i]));
		++histo[p];
	}
	int cumsum=0;
	int WHm=(W-2*MARGIN)*H/2;
	for(int i=0;i<65536;i++){
		cumsum+= histo[i];
		if(cumsum>WHm){
			median=i;
			break;
		}
	}
	for(int i=1;i<65536;++i) histo[i]+=histo[i-1];
	int center = 0;
	int cmax = 0;
	int spread = 1000;
	for(int i=spread+1; i<65536; i++) {
		int c = histo[i]-histo[i-spread-1];
		if(c>cmax) {
			cmax = c;
			center = i-spread/2;
		}
	}
	dmin = max(0,center-spread/2);
	dmax = min(65536, center+spread/2);
}

void visualizeArea(const vector <int>& imageData, int dmin, int dmax, int x0,int y0,int d,const string& fileName){
	unsigned char *data;
	data=new unsigned char [d*d*3];
    for (int y=0;y<d;y++) for (int x=0;x<d;x++)
    {
        int i=x0+x+(y0+y)*W;
        int j=x+y*d;
		if(x0+x<0 || x0+x>=W || y0+y<0 || y0+y>=H){
        	data[3*j]=0;
        	data[3*j+1]=0;
        	data[3*j+2]=0;
        	continue;
        }
		int ival = max(0, min( 255, (255*(imageData[i]-dmin))/(dmax-dmin) ) );
        if (ival<0) ival = 0;
        if (ival>255) ival = 255;
    	data[3*j]=ival;
       	data[3*j+1]=ival;
       	data[3*j+2]=ival;
    }
    for(int x=0;x<16;x++){
    	int j=d/2-8+x+(d/2-8)*d;
    	data[3*j]=255;
       	data[3*j+1]=0;
       	data[3*j+2]=0;
    	j=d/2-8+x+(d/2+8)*d;
    	data[3*j]=255;
       	data[3*j+1]=0;
       	data[3*j+2]=0;
    	j=d/2-8+(d/2-8+x)*d;
    	data[3*j]=255;
       	data[3*j+1]=0;
       	data[3*j+2]=0;
    	j=d/2+8+(d/2-8+x)*d;
    	data[3*j]=255;
       	data[3*j+1]=0;
       	data[3*j+2]=0;
    }
    #ifdef LOCAL
	stbi_write_bmp(fileName.c_str(),d,d,3,data);
	#endif
	delete[] data;
}
void visualizeAreaShapes(const vector <int>& imageData,const vector <int>& mean, int limit, int x0,int y0,int d,const string& fileName){
	unsigned char *data;
	data=new unsigned char [d*d*3];
    for (int y=0;y<d;y++) for (int x=0;x<d;x++)
    {
        int i=x0+x+(y0+y)*W;
        int j=x+y*d;
		if(x0+x<0 || x0+x>=W || y0+y<0 || y0+y>=H){
        	data[3*j]=0;
        	data[3*j+1]=0;
        	data[3*j+2]=0;
        	continue;
        }
		int ival = imageData[i]>mean[i]+limit ? 0 : (imageData[i]>mean[i]+limit/3 ? 196 : 255);
    	data[3*j]=ival;
       	data[3*j+1]=ival;
       	data[3*j+2]=ival;
    }
    for(int x=0;x<16;x++){
    	int j=d/2-8+x+(d/2-8)*d;
    	data[3*j]=255;
       	data[3*j+1]=0;
       	data[3*j+2]=0;
    	j=d/2-8+x+(d/2+8)*d;
    	data[3*j]=255;
       	data[3*j+1]=0;
       	data[3*j+2]=0;
    	j=d/2-8+(d/2-8+x)*d;
    	data[3*j]=255;
       	data[3*j+1]=0;
       	data[3*j+2]=0;
    	j=d/2+8+(d/2-8+x)*d;
    	data[3*j]=255;
       	data[3*j+1]=0;
       	data[3*j+2]=0;
    }
    #ifdef LOCAL
    stbi_write_bmp(fileName.c_str(),d,d,3,data);
	#endif
	delete[] data;
}
void visualizeAreaPeaks(const vector <char>& peak, int x0,int y0,int d,const string& fileName){
	unsigned char *data;
	data=new unsigned char [d*d*3];
	for(int i=0;i<d*d*3;i++) data[i]=255;
    for (int y=0;y<d;y++) for (int x=0;x<d;x++)
    {
        int i=x0+x+(y0+y)*W;
        int j=x+y*d;
		if(x0+x<0 || x0+x>=W || y0+y<0 || y0+y>=H){
        	data[3*j]=0;
        	data[3*j+1]=0;
        	data[3*j+2]=0;
        	continue;
        }
		if(peak[i]>=5) for(int p=1;p<=peak[i]-3;p++){
			if(x+p<d){
				data[3*(j+p)]=0;
    	   		data[3*(j+p)+1]=0;
       			data[3*(j+p)+2]=0;
			}
			if(x-p>=0){
				data[3*(j-p)]=0;
    	   		data[3*(j-p)+1]=0;
       			data[3*(j-p)+2]=0;
			}
			if(y+p<d){
				data[3*(j+p*d)]=0;
    	   		data[3*(j+p*d)+1]=0;
       			data[3*(j+p*d)+2]=0;
			}
			if(y-p>=0){
				data[3*(j-p*d)]=0;
    	   		data[3*(j-p*d)+1]=0;
       			data[3*(j-p*d)+2]=0;
			}
		}
    }
    for(int x=0;x<16;x++){
    	int j=d/2-8+x+(d/2-8)*d;
    	data[3*j]=255;
       	data[3*j+1]=0;
       	data[3*j+2]=0;
    	j=d/2-8+x+(d/2+8)*d;
    	data[3*j]=255;
       	data[3*j+1]=0;
       	data[3*j+2]=0;
    	j=d/2-8+(d/2-8+x)*d;
    	data[3*j]=255;
       	data[3*j+1]=0;
       	data[3*j+2]=0;
    	j=d/2+8+(d/2-8+x)*d;
    	data[3*j]=255;
       	data[3*j+1]=0;
       	data[3*j+2]=0;
    }
    #ifdef LOCAL
    stbi_write_bmp(fileName.c_str(),d,d,3,data);
	#endif
	delete[] data;
}


int maxInRadius(const vector <int>& imageData,int x,int y){
	int SEARCHRADIUS=15;
	int maxData=0;
	for(int i=max(0,x-SEARCHRADIUS);i<=x+SEARCHRADIUS && i<W;i++) for(int j=max(0,y-SEARCHRADIUS);j<=y+SEARCHRADIUS && j<H;j++){
		if(imageData[i+W*j]>maxData) maxData=imageData[i+W*j];
	}
	return maxData;
}

class Detection{
  public:
  	int id;
  	vector <double> ra;
  	vector <double> dec;
  	vector <int> x;
  	vector <int> y;
  	vector <double> magnitude;
  	bool isNEO;
  	bool isTwo;
  	Detection(){
  	}
  	Detection(const vector <string>& det,int position){
  		vector <string> value;
  		ra.resize(4,0.0);
  		dec.resize(4,0.0);
  		x.resize(4,-1);
  		y.resize(4,-1);
  		magnitude.resize(4,0.0);
  		isNEO=false;
  		isTwo=false;
		for(int i=0;i<4;i++){
  			int j=position+i;
  			vector <string> value=split(det[j]);
  			if(i==0){
  				id=atoi(value[0].c_str());
  				if(value[7]=="0") isNEO=false; else isNEO=true;
  			}
  			if(atoi(value[0].c_str())==id){
	  			int f=atoi(value[1].c_str())-1;
				ra[f]=atof(value[2].c_str());
				dec[f]=atof(value[3].c_str());
				x[f]=atoi(value[4].c_str());
				y[f]=atoi(value[5].c_str());
				magnitude[f]=atof(value[6].c_str());
  			}
  		}
  	}
  	void visualize(int number,const vector <int>& imageData1,int dmin1,int dmax1,const vector <int>& imageData2,int dmin2,int dmax2,const vector <int>& imageData3,int dmin3,int dmax3,const vector <int>& imageData4,int dmin4,int dmax4){
  		int d=128;
		visualizeArea(imageData1,dmin1,dmax1,x[0]-d/2,y[0]-d/2,d,"vis_"+SSTR(number)+"_"+SSTR(id)+"_1.bmp");
		visualizeArea(imageData2,dmin2,dmax2,x[1]-d/2,y[1]-d/2,d,"vis_"+SSTR(number)+"_"+SSTR(id)+"_2.bmp");
		visualizeArea(imageData3,dmin3,dmax3,x[2]-d/2,y[2]-d/2,d,"vis_"+SSTR(number)+"_"+SSTR(id)+"_3.bmp");
		visualizeArea(imageData4,dmin4,dmax4,x[3]-d/2,y[3]-d/2,d,"vis_"+SSTR(number)+"_"+SSTR(id)+"_4.bmp");
  	}
  	void visualizeShapes(int number,int d,const vector <int>& limit,const vector <int>& imageData1,const vector <int>& imageData2,const vector <int>& imageData3,const vector <int>& imageData4,const vector <int>& mean1,const vector <int>& mean2,const vector <int>& mean3,const vector <int>& mean4){
  		//int d=2000;
		visualizeAreaShapes(imageData1,mean1,limit[0],x[0]-d/2,y[0]-d/2,d,"viss_"+SSTR(number)+"_"+SSTR(id)+"_1.bmp");
		visualizeAreaShapes(imageData2,mean2,limit[1],x[1]-d/2,y[1]-d/2,d,"viss_"+SSTR(number)+"_"+SSTR(id)+"_2.bmp");
		visualizeAreaShapes(imageData3,mean3,limit[2],x[2]-d/2,y[2]-d/2,d,"viss_"+SSTR(number)+"_"+SSTR(id)+"_3.bmp");
		visualizeAreaShapes(imageData4,mean4,limit[3],x[3]-d/2,y[3]-d/2,d,"viss_"+SSTR(number)+"_"+SSTR(id)+"_4.bmp");
  	}
  	void visualizePeaks(int number,int d,const vector <char>& peaks1,const vector <char>& peaks2,const vector <char>& peaks3,const vector <char>& peaks4){
		visualizeAreaPeaks(peaks1,x[0]-d/2,y[0]-d/2,d,"viss_"+SSTR(number)+"_"+SSTR(id)+"_1p.bmp");
		visualizeAreaPeaks(peaks2,x[1]-d/2,y[1]-d/2,d,"viss_"+SSTR(number)+"_"+SSTR(id)+"_2p.bmp");
		visualizeAreaPeaks(peaks3,x[2]-d/2,y[2]-d/2,d,"viss_"+SSTR(number)+"_"+SSTR(id)+"_3p.bmp");
		visualizeAreaPeaks(peaks4,x[3]-d/2,y[3]-d/2,d,"viss_"+SSTR(number)+"_"+SSTR(id)+"_4p.bmp");
  	}
  	void getProperties(int number,const vector <int>& median,const vector <int>& imageData1,const vector <int>& imageData2,const vector <int>& imageData3,const vector <int>& imageData4){
  		vector <int> maxData(4);
		maxData[0]=maxInRadius(imageData1,x[0],y[0]);
  		maxData[1]=maxInRadius(imageData2,x[1],y[1]);
  		maxData[2]=maxInRadius(imageData3,x[2],y[2]);
  		maxData[3]=maxInRadius(imageData4,x[3],y[3]);
  		#ifdef LOCAL
		subor << number << ", " << id;
		for(int i=0;i<4;i++) subor << ", " << maxData[i];
		for(int i=0;i<4;i++) subor << ", " << maxData[i]-median[i];
		for(int i=0;i<4;i++) subor << ", " << static_cast<double>(maxData[i])/median[i];
		subor << "\n";
  		#endif
	}
};

class Adept{
  public:
  	double x;
  	double y;
  	double ra;
  	double dec;
  	char peak;
  	int lum;
	int overmed;
	int overmean;
	int scaleX;
	int scaleY;
	int scaleN;
	bool isStar;
	Adept(){
	}
	Adept(const double& xx,const double& yy,const double& rra,const double& ddec,char ppeak,int llum,int omed,int omean,int scx,int scy,int scn){
		x=xx;y=yy;ra=rra;dec=ddec;peak=ppeak;lum=llum;overmed=omed;overmean=omean;scaleX=scx;scaleY=scy;scaleN=scn;isStar=false;
	}
};
bool checkCandidate(const Adept& ad1,const Adept& ad2,const Adept& ad3,const Adept& ad4,double & angleDif){
	double norm1=sqrt(SQDIST(ad1.ra,ad1.dec,ad2.ra,ad2.dec));
	double norm2=sqrt(SQDIST(ad2.ra,ad2.dec,ad3.ra,ad3.dec));
	double norm3=sqrt(SQDIST(ad3.ra,ad3.dec,ad4.ra,ad4.dec));
	double minNorm=min(norm1,min(norm2,norm3));
	double maxNorm=max(norm1,max(norm2,norm3));
	if(minNorm<maxNorm*NORMRATIOLIMIT) return false;
	double angle1=atan2(ad2.ra-ad1.ra,ad2.dec-ad1.dec);
	double angle2=atan2(ad3.ra-ad2.ra,ad3.dec-ad2.dec);
	double angle3=atan2(ad4.ra-ad3.ra,ad4.dec-ad3.dec);
	double difa1=fabs(angle2-angle1); if(difa1>PI) difa1=2*PI-difa1;
	double difa2=fabs(angle3-angle2); if(difa2>PI) difa2=2*PI-difa2;
	double difa3=fabs(angle1-angle3); if(difa3>PI) difa3=2*PI-difa3;
	angleDif=max(difa1,max(difa2,difa3));
	if(difa1>ANGLEDIFLIMIT || difa2>ANGLEDIFLIMIT || difa3>ANGLEDIFLIMIT) return false;
	return true;	
}
double span(const double& v1,const double& v2,const double& v3,const double& v4){
	vector <double> val(4);
	val[0]=v1;val[1]=v2;val[2]=v3;val[3]=v4;
	sort(val.begin(),val.end());
	//mean=(val[1]+val[2])*0.5;
	return val[3]-val[0];
}
class Candidate{
  public:
	vector <double> ra;
	vector <double> dec;
	vector <double> x;
	vector <double> y;
	vector <int> data;
	bool isDetection;
	bool isNEO;
	bool isTwo;
	char starCount;
	double score;
	string imageID;
	vector <double> feature;
	vector <double> result;
	int temp;
	Candidate(){
		ra.resize(4);
		dec.resize(4);
	}
	Candidate(const Adept& ad1,const Adept& ad2,const Adept& ad3,const Adept& ad4,const string& imid,const double& duration,const double& season,const double& zd,const double& angleDif){
		result.resize(PV,0.0);
		feature.resize(FEATURES,0.0);
		imageID=imid;
		starCount=0;
		if(ad1.isStar) starCount++;
		if(ad2.isStar) starCount++;
		if(ad3.isStar) starCount++;
		if(ad4.isStar) starCount++;
		ra.resize(4);
		dec.resize(4);
		x.resize(4);
		y.resize(4);
		ra[0]=ad1.ra;
		ra[1]=ad2.ra;
		ra[2]=ad3.ra;
		ra[3]=ad4.ra;
		dec[0]=ad1.dec;
		dec[1]=ad2.dec;
		dec[2]=ad3.dec;
		dec[3]=ad4.dec;
		x[0]=ad1.x;
		x[1]=ad2.x;
		x[2]=ad3.x;
		x[3]=ad4.x;
		y[0]=ad1.y;
		y[1]=ad2.y;
		y[2]=ad3.y;
		y[3]=ad4.y;
		data.resize(4*8);
		data[0]=ad1.peak;data[1]=ad2.peak;data[2]=ad3.peak;data[3]=ad4.peak;
		data[4]=ad1.lum;data[5]=ad2.lum;data[6]=ad3.lum;data[7]=ad4.lum;
		data[8]=ad1.overmed;data[9]=ad2.overmed;data[10]=ad3.overmed;data[11]=ad4.overmed;
		data[12]=ad1.overmean;data[13]=ad2.overmean;data[14]=ad3.overmean;data[15]=ad4.overmean;
		data[16]=ad1.scaleX;data[17]=ad2.scaleX;data[18]=ad3.scaleX;data[19]=ad4.scaleX;
		data[20]=ad1.scaleY;data[21]=ad2.scaleY;data[22]=ad3.scaleY;data[23]=ad4.scaleY;
		data[24]=ad1.scaleN;data[25]=ad2.scaleN;data[26]=ad3.scaleN;data[27]=ad4.scaleN;
		
		score=1;
		for(int f=8;f<12;f++) if(data[f]<46 || data[f]>8813) score*=0.5;
		for(int f=8;f<12;f++) if(data[f]<62 || data[f]>1352) score*=0.5;
		for(int f=24;f<28;f++) if(data[f]>50) score*=0.5;
		for(int f=24;f<28;f++) if(data[f]>18) score*=0.5;
		for(int f=16;f<20;f++) if((data[f]+1)*2<(data[f+4]+1) || (data[f]+1)>2*(data[f+4]+1)) score*=0.5;
		
		isDetection=false;
		
		feature[0]=sqrt(SQDIST(ra[0],dec[0],ra[3],dec[3]));
		feature[1]=(ra[0]+ra[1]+ra[2]+ra[3])/4;
		feature[2]=(dec[0]+dec[1]+dec[2]+dec[3])/4;
		feature[3]=(x[0]+x[1]+x[2]+x[3])/4;
		feature[4]=(y[0]+y[1]+y[2]+y[3])/4;
		feature[5]=static_cast<double>(data[0]+data[1]+data[2]+data[3])/4;
		feature[6]=static_cast<double>(data[8]+data[9]+data[10]+data[11])/4;
		feature[7]=static_cast<double>(data[16]+data[17]+data[18]+data[19])/4;
		feature[8]=static_cast<double>(data[20]+data[21]+data[22]+data[23])/4;
		feature[9]=static_cast<double>(data[24]+data[25]+data[26]+data[27])/4;
		feature[10]=(static_cast<double>(data[16]+1)/(data[20]+1)+static_cast<double>(data[17]+1)/(data[21]+1)+static_cast<double>(data[18]+1)/(data[22]+1)+static_cast<double>(data[19]+1)/(data[23]+1))/4;
		feature[11]=totalLS(ra[0],ra[1],ra[2],ra[3],dec[0],dec[1],dec[2],dec[3],feature[12]);
		feature[13]=feature[11]/SQDIST(ra[0],dec[0],ra[3],dec[3]);
		feature[14]=atan2(dec[3]-dec[0],ra[3]-ra[0]);
		feature[15]=atan2(dec[3]-dec[0],ra[0]-ra[3]);
		feature[16]=span(data[0],data[1],data[2],data[3]);
		feature[17]=span(data[8],data[9],data[10],data[11]);
		feature[18]=span(data[16],data[17],data[18],data[19]);
		feature[19]=span(data[20],data[21],data[22],data[23]);
		feature[20]=span(data[24],data[25],data[26],data[27]);
		feature[21]=duration;
		feature[22]=feature[0]/duration;
		feature[23]=season;
		feature[24]=zd;
		feature[25]=min(min(sqrt(SQDIST(ra[0],dec[0],ra[1],dec[1])),sqrt(SQDIST(ra[1],dec[1],ra[2],dec[2]))),sqrt(SQDIST(ra[2],dec[2],ra[3],dec[3])))/feature[0];
		feature[26]=max(max(sqrt(SQDIST(ra[0],dec[0],ra[1],dec[1])),sqrt(SQDIST(ra[1],dec[1],ra[2],dec[2]))),sqrt(SQDIST(ra[2],dec[2],ra[3],dec[3])))/feature[0];
		feature[27]=sqrt(SQDIST(ra[1],dec[1],2*ra[0]/3+ra[3]/3,2*dec[0]/3+dec[3]/3));
		feature[28]=sqrt(SQDIST(ra[2],dec[2],2*ra[3]/3+ra[0]/3,2*dec[3]/3+dec[0]/3));
		feature[29]=feature[27]/feature[0];
		feature[30]=feature[28]/feature[0];
		feature[31]=angleDif;
	}
	int checkDetection(const vector <Detection>& det){
		isDetection=false;
		isNEO=false;
		isTwo=false;
		for(int i=0;i<det.size();i++){
			double sum = 0;
            for(int f=0;f<4;f++){
                sum += (ra[f] - det[i].ra[f]) * (ra[f] - det[i].ra[f]);
                sum += (dec[f] - det[i].dec[f]) * (dec[f] - det[i].dec[f]);
            }
            if(sum<MATCH_DISTANCE){
            	double difra=ra[3]-ra[0];
            	double difdec=dec[3]-dec[0];
            	double difdetra=det[i].ra[3]-det[i].ra[0];
            	double difdetdec=det[i].dec[3]-det[i].dec[0];
            	double ratioabs=(difra*difra+difdec*difdec)/(difdetra*difdetra+difdetdec*difdetdec);
            	double angledif=fabs(atan2(difra,difdec)-atan2(difdetra,difdetdec));
            	if(angledif>PI) angledif=2*PI-angledif;
				if(ratioabs<1.2 && ratioabs>0.8 && angledif<PI/18){
	            	isDetection=true;
	            	isNEO=det[i].isNEO;
	            	isTwo=det[i].isTwo;
	            	result[0]=1.0;
	            	result[1]=isNEO?1.0:0.0;
	            	result[2]=isTwo?1.0:0.0;
	            	//cerr << i << " ; ";
	            	return i;
	            }
            }
		}
		result[0]=0.0;
		result[1]=0.0;
		return -1;
	}
  	void visualize(int number,int id,const vector <int>& imageData1,int dmin1,int dmax1,const vector <int>& imageData2,int dmin2,int dmax2,const vector <int>& imageData3,int dmin3,int dmax3,const vector <int>& imageData4,int dmin4,int dmax4){
  		int d=128;
		visualizeArea(imageData1,dmin1,dmax1,x[0]-d/2,y[0]-d/2,d,"vis_"+SSTR(number)+"_"+SSTR(id)+"_1.bmp");
		visualizeArea(imageData2,dmin2,dmax2,x[1]-d/2,y[1]-d/2,d,"vis_"+SSTR(number)+"_"+SSTR(id)+"_2.bmp");
		visualizeArea(imageData3,dmin3,dmax3,x[2]-d/2,y[2]-d/2,d,"vis_"+SSTR(number)+"_"+SSTR(id)+"_3.bmp");
		visualizeArea(imageData4,dmin4,dmax4,x[3]-d/2,y[3]-d/2,d,"vis_"+SSTR(number)+"_"+SSTR(id)+"_4.bmp");
  	}
 	bool operator < (const Candidate& str) const{
        return result[0]>str.result[0];
		//return score>str.score;
    }
};
bool sortByNEO(const Candidate& c1,const Candidate& c2){
	return c1.result[1]>c2.result[1] || (c1.result[1]==c2.result[1] && c1.result[0]>c2.result[0]);
}
bool sortByTwo(const Candidate& c1,const Candidate& c2){
	return c1.result[2]>c2.result[2] || (c1.result[2]==c2.result[2] && c1.result[0]>c2.result[0]);
}
vector <Candidate> cand;
vector <Candidate> candTest;
bool comp0(int i,int j){return (cand[i].feature[0]<cand[j].feature[0]);}
bool comp1(int i,int j){return (cand[i].feature[1]<cand[j].feature[1]);}
bool comp2(int i,int j){return (cand[i].feature[2]<cand[j].feature[2]);}
bool comp3(int i,int j){return (cand[i].feature[3]<cand[j].feature[3]);}
bool comp4(int i,int j){return (cand[i].feature[4]<cand[j].feature[4]);}
bool comp5(int i,int j){return (cand[i].feature[5]<cand[j].feature[5]);}
bool comp6(int i,int j){return (cand[i].feature[6]<cand[j].feature[6]);}
bool comp7(int i,int j){return (cand[i].feature[7]<cand[j].feature[7]);}
bool comp8(int i,int j){return (cand[i].feature[8]<cand[j].feature[8]);}
bool comp9(int i,int j){return (cand[i].feature[9]<cand[j].feature[9]);}
bool comp10(int i,int j){return (cand[i].feature[10]<cand[j].feature[10]);}
bool comp11(int i,int j){return (cand[i].feature[11]<cand[j].feature[11]);}
bool comp12(int i,int j){return (cand[i].feature[12]<cand[j].feature[12]);}
bool comp13(int i,int j){return (cand[i].feature[13]<cand[j].feature[13]);}
bool comp14(int i,int j){return (cand[i].feature[14]<cand[j].feature[14]);}
bool comp15(int i,int j){return (cand[i].feature[15]<cand[j].feature[15]);}
bool comp16(int i,int j){return (cand[i].feature[16]<cand[j].feature[16]);}
bool comp17(int i,int j){return (cand[i].feature[17]<cand[j].feature[17]);}
bool comp18(int i,int j){return (cand[i].feature[18]<cand[j].feature[18]);}
bool comp19(int i,int j){return (cand[i].feature[19]<cand[j].feature[19]);}
bool comp20(int i,int j){return (cand[i].feature[20]<cand[j].feature[20]);}
bool comp21(int i,int j){return (cand[i].feature[21]<cand[j].feature[21]);}
bool comp22(int i,int j){return (cand[i].feature[22]<cand[j].feature[22]);}
bool comp23(int i,int j){return (cand[i].feature[23]<cand[j].feature[23]);}
bool comp24(int i,int j){return (cand[i].feature[24]<cand[j].feature[24]);}
bool comp25(int i,int j){return (cand[i].feature[25]<cand[j].feature[25]);}
bool comp26(int i,int j){return (cand[i].feature[26]<cand[j].feature[26]);}
bool comp27(int i,int j){return (cand[i].feature[27]<cand[j].feature[27]);}
bool comp28(int i,int j){return (cand[i].feature[28]<cand[j].feature[28]);}
bool comp29(int i,int j){return (cand[i].feature[29]<cand[j].feature[29]);}
bool comp30(int i,int j){return (cand[i].feature[30]<cand[j].feature[30]);}
bool comp31(int i,int j){return (cand[i].feature[31]<cand[j].feature[31]);}
bool comp32(int i,int j){return (cand[i].feature[32]<cand[j].feature[32]);}
bool comp33(int i,int j){return (cand[i].feature[33]<cand[j].feature[33]);}
bool comp34(int i,int j){return (cand[i].feature[34]<cand[j].feature[34]);}
bool comp35(int i,int j){return (cand[i].feature[35]<cand[j].feature[35]);}
bool comp36(int i,int j){return (cand[i].feature[36]<cand[j].feature[36]);}
bool comp37(int i,int j){return (cand[i].feature[37]<cand[j].feature[37]);}
bool comp38(int i,int j){return (cand[i].feature[38]<cand[j].feature[38]);}
bool comp39(int i,int j){return (cand[i].feature[39]<cand[j].feature[39]);}
bool comp40(int i,int j){return (cand[i].feature[40]<cand[j].feature[40]);}
bool comp41(int i,int j){return (cand[i].feature[41]<cand[j].feature[41]);}
bool comp42(int i,int j){return (cand[i].feature[42]<cand[j].feature[42]);}
bool comp43(int i,int j){return (cand[i].feature[43]<cand[j].feature[43]);}
bool comp44(int i,int j){return (cand[i].feature[44]<cand[j].feature[44]);}
bool comp45(int i,int j){return (cand[i].feature[45]<cand[j].feature[45]);}
bool comp46(int i,int j){return (cand[i].feature[46]<cand[j].feature[46]);}
bool comp47(int i,int j){return (cand[i].feature[47]<cand[j].feature[47]);}
bool comp48(int i,int j){return (cand[i].feature[48]<cand[j].feature[48]);}
bool comp49(int i,int j){return (cand[i].feature[49]<cand[j].feature[49]);}
bool comp50(int i,int j){return (cand[i].feature[50]<cand[j].feature[50]);}
bool comp51(int i,int j){return (cand[i].feature[51]<cand[j].feature[51]);}
bool comp52(int i,int j){return (cand[i].feature[52]<cand[j].feature[52]);}
bool comp53(int i,int j){return (cand[i].feature[53]<cand[j].feature[53]);}
bool comp54(int i,int j){return (cand[i].feature[54]<cand[j].feature[54]);}
bool comp55(int i,int j){return (cand[i].feature[55]<cand[j].feature[55]);}
bool comp56(int i,int j){return (cand[i].feature[56]<cand[j].feature[56]);}
bool comp57(int i,int j){return (cand[i].feature[57]<cand[j].feature[57]);}
bool comp58(int i,int j){return (cand[i].feature[58]<cand[j].feature[58]);}
bool comp59(int i,int j){return (cand[i].feature[59]<cand[j].feature[59]);}
bool comp60(int i,int j){return (cand[i].feature[60]<cand[j].feature[60]);}
bool comp61(int i,int j){return (cand[i].feature[61]<cand[j].feature[61]);}
bool comp62(int i,int j){return (cand[i].feature[62]<cand[j].feature[62]);}
bool comp63(int i,int j){return (cand[i].feature[63]<cand[j].feature[63]);}
bool comp64(int i,int j){return (cand[i].feature[64]<cand[j].feature[64]);}
bool comp65(int i,int j){return (cand[i].feature[65]<cand[j].feature[65]);}
bool comp66(int i,int j){return (cand[i].feature[66]<cand[j].feature[66]);}
bool comp67(int i,int j){return (cand[i].feature[67]<cand[j].feature[67]);}
bool comp68(int i,int j){return (cand[i].feature[68]<cand[j].feature[68]);}
bool comp69(int i,int j){return (cand[i].feature[69]<cand[j].feature[69]);}
bool comp70(int i,int j){return (cand[i].feature[70]<cand[j].feature[70]);}
bool comp71(int i,int j){return (cand[i].feature[71]<cand[j].feature[71]);}
bool comp72(int i,int j){return (cand[i].feature[72]<cand[j].feature[72]);}
bool comp73(int i,int j){return (cand[i].feature[73]<cand[j].feature[73]);}
bool comp74(int i,int j){return (cand[i].feature[74]<cand[j].feature[74]);}
bool comp75(int i,int j){return (cand[i].feature[75]<cand[j].feature[75]);}
bool comp76(int i,int j){return (cand[i].feature[76]<cand[j].feature[76]);}
bool comp77(int i,int j){return (cand[i].feature[77]<cand[j].feature[77]);}
bool comp78(int i,int j){return (cand[i].feature[78]<cand[j].feature[78]);}
bool comp79(int i,int j){return (cand[i].feature[79]<cand[j].feature[79]);}
class Sample{
  public:
  	vector <int> id;
  	void sortBy(int feature){
  		switch(feature){
  			case 0: sort(id.begin(),id.end(),comp0);break;
  			case 1: sort(id.begin(),id.end(),comp1);break;
  			case 2: sort(id.begin(),id.end(),comp2);break;
  			case 3: sort(id.begin(),id.end(),comp3);break;
  			case 4: sort(id.begin(),id.end(),comp4);break;
  			case 5: sort(id.begin(),id.end(),comp5);break;
  			case 6: sort(id.begin(),id.end(),comp6);break;
  			case 7: sort(id.begin(),id.end(),comp7);break;
  			case 8: sort(id.begin(),id.end(),comp8);break;
  			case 9: sort(id.begin(),id.end(),comp9);break;
  			case 10: sort(id.begin(),id.end(),comp10);break;
  			case 11: sort(id.begin(),id.end(),comp11);break;
  			case 12: sort(id.begin(),id.end(),comp12);break;
  			case 13: sort(id.begin(),id.end(),comp13);break;
  			case 14: sort(id.begin(),id.end(),comp14);break;
  			case 15: sort(id.begin(),id.end(),comp15);break;
  			case 16: sort(id.begin(),id.end(),comp16);break;
  			case 17: sort(id.begin(),id.end(),comp17);break;
  			case 18: sort(id.begin(),id.end(),comp18);break;
  			case 19: sort(id.begin(),id.end(),comp19);break;
  			case 20: sort(id.begin(),id.end(),comp20);break;
  			case 21: sort(id.begin(),id.end(),comp21);break;
  			case 22: sort(id.begin(),id.end(),comp22);break;
  			case 23: sort(id.begin(),id.end(),comp23);break;
  			case 24: sort(id.begin(),id.end(),comp24);break;
  			case 25: sort(id.begin(),id.end(),comp25);break;
  			case 26: sort(id.begin(),id.end(),comp26);break;
  			case 27: sort(id.begin(),id.end(),comp27);break;
  			case 28: sort(id.begin(),id.end(),comp28);break;
  			case 29: sort(id.begin(),id.end(),comp29);break;
  			case 30: sort(id.begin(),id.end(),comp30);break;
  			case 31: sort(id.begin(),id.end(),comp31);break;
  			case 32: sort(id.begin(),id.end(),comp32);break;
  			case 33: sort(id.begin(),id.end(),comp33);break;
  			case 34: sort(id.begin(),id.end(),comp34);break;
  			case 35: sort(id.begin(),id.end(),comp35);break;
  			case 36: sort(id.begin(),id.end(),comp36);break;
  			case 37: sort(id.begin(),id.end(),comp37);break;
  			case 38: sort(id.begin(),id.end(),comp38);break;
  			case 39: sort(id.begin(),id.end(),comp39);break;
  			case 40: sort(id.begin(),id.end(),comp40);break;
  			case 41: sort(id.begin(),id.end(),comp41);break;
  			case 42: sort(id.begin(),id.end(),comp42);break;
  			case 43: sort(id.begin(),id.end(),comp43);break;
  			case 44: sort(id.begin(),id.end(),comp44);break;
  			case 45: sort(id.begin(),id.end(),comp45);break;
  			case 46: sort(id.begin(),id.end(),comp46);break;
  			case 47: sort(id.begin(),id.end(),comp47);break;
  			case 48: sort(id.begin(),id.end(),comp48);break;
  			case 49: sort(id.begin(),id.end(),comp49);break;
  			case 50: sort(id.begin(),id.end(),comp50);break;
  			case 51: sort(id.begin(),id.end(),comp51);break;
  			case 52: sort(id.begin(),id.end(),comp52);break;
  			case 53: sort(id.begin(),id.end(),comp53);break;
  			case 54: sort(id.begin(),id.end(),comp54);break;
  			case 55: sort(id.begin(),id.end(),comp55);break;
  			case 56: sort(id.begin(),id.end(),comp56);break;
  			case 57: sort(id.begin(),id.end(),comp57);break;
  			case 58: sort(id.begin(),id.end(),comp58);break;
  			case 59: sort(id.begin(),id.end(),comp59);break;
  			case 60: sort(id.begin(),id.end(),comp60);break;
  			case 61: sort(id.begin(),id.end(),comp61);break;
  			case 62: sort(id.begin(),id.end(),comp62);break;
  			case 63: sort(id.begin(),id.end(),comp63);break;
  			case 64: sort(id.begin(),id.end(),comp64);break;
  			case 65: sort(id.begin(),id.end(),comp65);break;
  			case 66: sort(id.begin(),id.end(),comp66);break;
  			case 67: sort(id.begin(),id.end(),comp67);break;
  			case 68: sort(id.begin(),id.end(),comp68);break;
  			case 69: sort(id.begin(),id.end(),comp69);break;
  			case 70: sort(id.begin(),id.end(),comp70);break;
  			case 71: sort(id.begin(),id.end(),comp71);break;
  			case 72: sort(id.begin(),id.end(),comp72);break;
  			case 73: sort(id.begin(),id.end(),comp73);break;
  			case 74: sort(id.begin(),id.end(),comp74);break;
  			case 75: sort(id.begin(),id.end(),comp75);break;
  			case 76: sort(id.begin(),id.end(),comp76);break;
  			case 77: sort(id.begin(),id.end(),comp77);break;
  			case 78: sort(id.begin(),id.end(),comp78);break;
  			case 79: sort(id.begin(),id.end(),comp79);break;
  		}
  	}
  	Sample(){}
  	Sample(int size){
  		id.clear();
  		id.reserve(size);
  	}
};
class Node{
  public:
	int left;
	int right;
	int feature;
	double value;
	int level;
	vector <int> counts;
	int total;
	Node(int by){
		left=-1;
		right=-1;
		feature=-1;
		value=0.0;
		level=-1;
		counts.resize(CLUSTERS[by],0);
		total=0;
	}
	Node(int lev,vector <int> cou,int tot){
		left=-1;
		right=-1;
		feature=-1;
		value=0.0;
		level=lev;
		counts=cou;
		total=tot;
	}
};
class Tree{
  public:
  	vector<Node> node;
    int clusterAtNode(int nodeId,const Candidate& ad,int by){
    	if(node[nodeId].left==-1){
    		int cluster=0;
    		int cumCount=node[nodeId].counts[0];
			while(cluster<CLUSTERS[by]-1 && 2*cumCount<node[nodeId].total){
				cluster++;
				cumCount+=node[nodeId].counts[cluster];
			}
			return cluster;
		}
		if(ad.feature[node[nodeId].feature]<node[nodeId].value) return clusterAtNode(node[nodeId].left,ad,by);
		return clusterAtNode(node[nodeId].right,ad,by);
    }
	int assignCluster(const Candidate& ad,int by){
		return clusterAtNode(0,ad,by);
    }
	void divideNode(int nodeIndex,Sample& sample,int by){
		int n=sample.id.size();
		int nonzero=0;
		for(int i=0;i<CLUSTERS[by];i++){
			if(node[nodeIndex].counts[i]>0) nonzero++;
			if(nonzero>1) break;
		}
		if(node[nodeIndex].level<MAXLEVEL-1 && nonzero>1 && node[nodeIndex].total>MINNODE){
	    	vector <int> feaID=selectFeatures();
	    	double minCost=1000000.0; 
	    	int bestF=-1;
	    	double bestValue=0.0;
	    	int bestI=0;
	    	vector <int> bestC1(CLUSTERS[by],0);
	    	int bestTotalL=0;
			for(int f=0;f<FEATURETRY;f++){
				sample.sortBy(feaID[f]);
	    		vector <int> c1(CLUSTERS[by],0);
	    		int totalL=0;
	    		for(int i=0;i<n-1;i++){
	    			for(int cl=0;cl<CLUSTERS[by];cl++){
	    				if(cand[sample.id[i]].result[by]<VALUE[by][cl+1]){
	    					c1[cl]++;
	    					break;
	    				}
	    			}
	    			totalL++;
					if(cand[sample.id[i+1]].feature[feaID[f]]-cand[sample.id[i]].feature[feaID[f]]>EPS){
	    				// Gini 
						//double cost=1.0/n*(d1*r1/static_cast<double>(d1+r1)+(node[nodeIndex].detections-d1)*(node[nodeIndex].rejections-r1)/static_cast<double>(n-d1-r1));
	    			    double costL=0.0;
						double costR=0.0;
						for(int cl=0;cl<CLUSTERS[by];cl++){
							costL+=GCOST(c1[cl]/static_cast<double>(totalL),ROOT);
							costR+=GCOST((node[nodeIndex].counts[cl]-c1[cl])/static_cast<double>(n-totalL),ROOT);
						}
						double cost=(totalL*costL+(n-totalL)*costR)/n;
						if(cost<minCost){
	    			    	minCost=cost;
	    			    	bestF=feaID[f];
	    			    	bestValue=(cand[sample.id[i+1]].feature[feaID[f]]+cand[sample.id[i]].feature[feaID[f]])/2.0;
							bestI=i;
	    			    	bestC1=c1;
	    			    	bestTotalL=totalL;
	    			    }
					}
	    		}
	    	}
	    	if(bestF>=0){
	    		featureScore[bestF]+=n;
		    	Sample sampleLeft(bestI+1);
		    	Sample sampleRight(n-bestI-1);
				for(int i=0;i<n;i++){
		    		if(cand[sample.id[i]].feature[bestF]<bestValue){
						sampleLeft.id.push_back(sample.id[i]);
					}
		    		else sampleRight.id.push_back(sample.id[i]);
		    	}
		        node[nodeIndex].feature=bestF;
		    	node[nodeIndex].value=bestValue;
		    	node.push_back(Node(node[nodeIndex].level+1,bestC1,bestTotalL));
		    	node[nodeIndex].left=node.size()-1;
		    	vector<int> c2(CLUSTERS[by],0);
		    	for(int i=0;i<CLUSTERS[by];i++){
		    		c2[i]=node[nodeIndex].counts[i]-bestC1[i];	
		    	}
		    	node.push_back(Node(node[nodeIndex].level+1,c2,node[nodeIndex].total-bestTotalL));
		    	node[nodeIndex].right=node.size()-1;
			    divideNode(node[nodeIndex].left,sampleLeft,by);
				divideNode(node[nodeIndex].right,sampleRight,by);
			}
		}
	}
	Tree(){
  	}
};
vector <Tree> randomTree[PV];
Tree buildTree(int by){
	int n=cand.size();
	int ns=min(n,MAXSAMPLESIZE);
	Sample sample;
	sample.id.resize(ns);
	Tree tree;
	tree.node.resize(1,Node(by));
	tree.node[0].level=0;
	for(int i=0;i<ns;i++){
		sample.id[i]=randBig()%n;
		for(int cl=0;cl<CLUSTERS[by];cl++){
			if(cand[sample.id[i]].result[by]<VALUE[by][cl+1]){
				tree.node[0].counts[cl]++;
				break;
			}
		}
	}
	tree.node[0].total=ns;
	tree.divideNode(0,sample,by);
	return tree;
}

vector <Adept> calculateAdepts(const vector <int>& image,const vector <int>& mean,const vector <int>& median,const vector <ShapeInfo>& shape,const vector <int>& getShape,const vector <double>& wcs,vector < vector <int> >& adeptList,vector <int>& getAdept){
	adeptList.clear();
	adeptList.resize(WHA);
	vector <Adept> result;
	getAdept=vector <int> (W*H,-1);
	for(int y=0;y<H;y++) for(int x=0;x<W;x++){
		int i=W*y+x;
		char peak;
		if(x==0 || y==0 || x==W-1 || y==H-1
			|| image[i+1]>image[i] || image[i-1]>image[i] || image[i+W]>image[i] || image[i-W]>image[i]){
			peak=0;
		}
		else if(image[i+1+W]>image[i] || image[i-1+W]>image[i] || image[i+1-W]>image[i] || image[i-1-W]>image[i]){
			peak=1;
		}
		else if(x==1 || y==1 || x==W-2 || y==H-2
			|| image[i+2]>image[i+1] || image[i-2]>image[i-1] || image[i+2*W]>image[i+W] || image[i-2*W]>image[i-W]){
			peak=2;
		}
		else if(    (image[i+W+2]>image[i+W+1] && image[i+W+2]>image[i+1]) || (image[i+2*W+1]>image[i+W+1] && image[i+2*W+1]>image[i+W])
			|| (image[i-W+2]>image[i-W+1] && image[i-W+2]>image[i+1]) || (image[i-2*W+1]>image[i-W+1] && image[i-2*W+1]>image[i-W])
			|| (image[i+W-2]>image[i+W-1] && image[i+W-2]>image[i-1]) || (image[i+2*W-1]>image[i+W-1] && image[i+2*W-1]>image[i+W])
			|| (image[i-W-2]>image[i-W-1] && image[i-W-2]>image[i-1]) || (image[i-2*W-1]>image[i-W-1] && image[i-2*W-1]>image[i-W])
			){
			peak=3;
		}
		else if(image[i+2*W+2]>image[i+W+1] || image[i+2*W-2]>image[i+W-1] || image[i-2*W+2]>image[i-W+1] || image[i-2*W-2]>image[i-W-1]){
			peak=4;
		}
		else if(x==2 || y==2 || x==W-3 || y==H-3
			|| image[i+3]>image[i+2] || image[i-3]>image[i-2] || image[i+3*W]>image[i+2*W] || image[i-3*W]>image[i-2*W]){
			peak=5;
		}
		else if(    (image[i+W+3]>image[i+W+2] && image[i+W+3]>image[i+2]) || (image[i+3*W+1]>image[i+2*W+1] && image[i+3*W+1]>image[i+2*W])
			|| (image[i-W+3]>image[i-W+2] && image[i-W+3]>image[i+2]) || (image[i-3*W+1]>image[i-2*W+1] && image[i-3*W+1]>image[i-2*W])
			|| (image[i+W-3]>image[i+W-2] && image[i+W-3]>image[i-2]) || (image[i+3*W-1]>image[i+2*W-1] && image[i+3*W-1]>image[i+2*W])
			|| (image[i-W-3]>image[i-W-2] && image[i-W-3]>image[i-2]) || (image[i-3*W-1]>image[i-2*W-1] && image[i-3*W-1]>image[i-2*W])
			){
			peak=6;
		}
		else if(   (image[i+2*W+3]>image[i+2*W+2] && image[i+2*W+3]>image[i+W+2]) || (image[i+3*W+2]>image[i+2*W+2] && image[i+3*W+2]>image[i+2*W+1])
		   || (image[i-2*W+3]>image[i-2*W+2] && image[i-2*W+3]>image[i-W+2]) || (image[i-3*W+2]>image[i-2*W+2] && image[i-3*W+2]>image[i-2*W+1])
	       || (image[i+2*W-3]>image[i+2*W-2] && image[i+2*W-3]>image[i+W-2]) || (image[i+3*W-2]>image[i+2*W-2] && image[i+3*W-2]>image[i+2*W-1])
		   || (image[i-2*W-3]>image[i-2*W-2] && image[i-2*W-3]>image[i-W-2]) || (image[i-3*W-2]>image[i-2*W-2] && image[i-3*W-2]>image[i-2*W-1])
			){
			peak=7;
		}
		else if(x==3 || y==3 || x==W-4 || y==H-4
			|| image[i+4]>image[i+3] || image[i-4]>image[i-3] || image[i+4*W]>image[i+3*W] || image[i-4*W]>image[i-3*W]){
			peak=8;
		}
		else if(    (image[i+W+4]>image[i+W+3] && image[i+W+4]>image[i+3]) || (image[i+4*W+1]>image[i+3*W+1] && image[i+4*W+1]>image[i+3*W])
			|| (image[i-W+4]>image[i-W+3] && image[i-W+4]>image[i+3]) || (image[i-4*W+1]>image[i-3*W+1] && image[i-4*W+1]>image[i-3*W])
			|| (image[i+W-4]>image[i+W-3] && image[i+W-4]>image[i-3]) || (image[i+4*W-1]>image[i+3*W-1] && image[i+4*W-1]>image[i+3*W])
			|| (image[i-W-4]>image[i-W-3] && image[i-W-4]>image[i-3]) || (image[i-4*W-1]>image[i-3*W-1] && image[i-4*W-1]>image[i-3*W])
			){
			peak=9;
		}
		else if(image[i+3*W+3]>image[i+2*W+2] || image[i+3*W-3]>image[i+2*W-2] || image[i-3*W+3]>image[i-2*W+2] || image[i-3*W-3]>image[i-2*W-2]){
			peak=10;
		}
		else peak=11;
		
		// HERE DECIDE WHETHER IT IS ADEPT
		// available info: peak, mean, median, shape
		if(peak>=3 && x>=MARGIN && x<W-MARGIN){
			int scaleX=0;
			int scaleY=0;
			int scaleN=0;
			if(getShape[i]>=0){
				scaleX=shape[getShape[i]].scaleX;
				scaleY=shape[getShape[i]].scaleY;
				scaleN=shape[getShape[i]].scaleN;
			}
			bool take=true;
			int overmed=image[i]-medianAt(x,y,median);
			//if(overmed<-67 || overmed>24875) take=false;
			//else
			if(peak<=4 && (overmed<46 || overmed>8813)) take=false;
			if(overmed<-67 || overmed>24875) take=false;
			//else if(overmed<62 || overmed>1352) score*=0.1;
			//if(scaleN>4642) take=false;
			//else if(scaleN>50) score*=0.02;
			//else if(scaleN>18) score*=0.1;
			
			double ratioXY=(scaleX+1.0)/(scaleY+1.0);
			if(ratioXY<0.1621 || ratioXY>3) take=false;
			//else if(ratioXY<0.5 || ratioXY>2) score*=0.02;
			//else if(ratioXY<0.8333 || ratioXY>1.5) score*=0.1;
			
			if(scaleX>77) take=false;
			if(scaleY>474) take=false;
			if(scaleN>4642) take=false;
			
			if(take){
				double exx=x+0.7*(-2*image[i-2]-image[i-1]+image[i+1]+2*image[i+2])/(-2*image[i-2]+image[i-1]+2*image[i]+image[i+1]-2*image[i+2]);
				double exy=y+0.7*(-2*image[i-2*W]-image[i-W]+image[i+W]+2*image[i+2*W])/(-2*image[i-2*W]+image[i-W]+2*image[i]+image[i+W]-2*image[i+2*W]);
				double ra,dec;
				wc.convertXY2RADEC(wcs,exx,exy,ra,dec);
				result.push_back(Adept(exx,exy,ra,dec,peak,image[i],overmed,image[i]-mean[i],scaleX,scaleY,scaleN));
				adeptList[x/ADEPTGRID+y/ADEPTGRID*WA].push_back(result.size()-1);
				getAdept[i]=result.size()-1;
			}
		}		
	}
	return result;
}

vector <string> getHeader(const vector <string>& header,const vector <string>& item){
	vector <string> result(item.size(),"");
	for(int i=0;i<header.size();i++){
		vector <string> word=splitBy(header[i],'=');
		if(word.size()>=2){
			string title=trim(word[0]);
			for(int j=0;j<item.size();j++) if(item[j]==title){
				vector <string> value=splitBy(word[1],'/');
				if(value.size()>=1){
					result[j]=trim(value[0]," '\t");
				}
				break;
			}
		}
	}
	return result;
}

vector <Candidate> getCandidates(int width, int height, const vector <int>& imageData_1, const vector <string>& header_1, const vector <double>& wcs_1, const vector <int>& imageData_2, const vector <string>& header_2, const vector <double>& wcs_2, const vector <int>& imageData_3, const vector <string>& header_3, const vector <double>& wcs_3, const vector <int>& imageData_4, const vector <string>& header_4, const vector <double>& wcs_4,const string& imageID){
	W=width;
	H=height;
	WH=W*H;
	WM=ceil(static_cast<double>(W-MARGIN)/MEDIANGRID);
	HM=ceil(static_cast<double>(H)/MEDIANGRID);
	WHM=WM*HM;
	WA=ceil(static_cast<double>(W)/ADEPTGRID);
	HA=ceil(static_cast<double>(H)/ADEPTGRID);
	WHA=WA*HA;
			
	vector <double> wcs[4];
	wcs[0]=wcs_1;
	wcs[1]=wcs_2;
	wcs[2]=wcs_3;
	wcs[3]=wcs_4;
	
	vector <Adept> adept[4];
	vector <vector <int> > adeptList[4];
	vector <int> getAdept[4];

	double st=getTime();
	for(int i=0;i<4;i++){
		vector <int> image;
		switch(i){
			case 0: image=imageData_1;break;
			case 1: image=imageData_2;break;
			case 2: image=imageData_3;break;
			case 3: image=imageData_4;break;
		}
		if(BLURING) image=blur(image);
		vector <int> median=medians(image);
		vector <int> mean=means2(image,median,MAXOVERMEDIAN);
		vector <int> getShape;
		vector <ShapeInfo> shape=shapeDetection(image,mean,SHAPEOVERMEANLIMIT,getShape);
		adept[i]=calculateAdepts(image,mean,median,shape,getShape,wcs[i],adeptList[i],getAdept[i]);
		cerr << getTime()-st << " sec.\n";st=getTime();
		cerr << adept[i].size() << "\n";
	}
	
	int starCount=0;
	for(int i=0;i<adept[0].size();i++){
		vector <double> xf(4);
		vector <double> yf(4);
		vector <bool> is(4,false);
		vector <int> id(4,-1);
		for(int f=1;f<4;f++){
			wc.convertRADEC2XY(wcs[f],adept[0][i].ra,adept[0][i].dec,xf[f],yf[f]);
			int xx=round(xf[f]);int yy=round(yf[f]);
			if(xx<1) xx=1; if(yy<1) yy=1; if(xx>W-2) xx=W-2; if(yy>H-2) yy=H-2;
			for(int x=xx-1;x<=xx+1;x++) for(int y=yy-1;y<=yy+1;y++){
				if(getAdept[f][x+W*y]>=0){
					int idf=getAdept[f][x+W*y];
					double radif=adept[0][i].ra-adept[f][idf].ra;
					double decdif=adept[0][i].dec-adept[f][idf].dec;
					if(radif*radif+decdif*decdif<STARLIMITSQ){
						is[f]=true;
						id[f]=idf;
					}
				}
			}
		}
		if(is[1] && is[2] && is[3]){
			adept[0][i].isStar=true;
			for(int f=1;f<4;f++) adept[f][id[f]].isStar=true;
			starCount++;
		}
	}
	cerr << getTime()-st << " sec.\n";st=getTime();
	cerr << starCount << " stars.\n";
	
	vector <string> item(11);
	item[0]="MJD";
	item[1]="DATE-OBS";
	item[2]="TIME-OBS";
	item[3]="ST";
	item[4]="ZD";
	item[5]="HA";
	item[6]="EXPTIME";
	item[7]="AIRMASS";
	item[8]="TELFOCUS";
	item[9]="TELTEMP";
	item[10]="CAMTEMP";
		
	vector <string> info[4];
	info[0]=getHeader(header_1,item);
	info[1]=getHeader(header_2,item);
	info[2]=getHeader(header_3,item);
	info[3]=getHeader(header_4,item);
	
	vector<double> mjd(4);
	for(int f=0;f<4;f++){
		mjd[f]=atof(info[f][0].c_str());
	}
	double duration=mjd[3]-mjd[0];
	double season=mjd[0]/365.25-floor(mjd[0]/365.25);
	double zd=(atof(info[0][4].c_str())+atof(info[3][4].c_str()))/2;
	
	vector <Candidate> candidate;
	
	for(int i=0;i<adept[0].size();i++) if(!adept[0][i].isStar){
		double x4,y4;
		//int sc1=adept[0][i].isStar?1:0;
		wc.convertRADEC2XY(wcs[3],adept[0][i].ra,adept[0][i].dec,x4,y4);
		int x=floor(x4/ADEPTGRID);
		int y=floor(y4/ADEPTGRID);
		for(int dx=-1;dx<=1;dx++) for(int dy=-1;dy<=1;dy++) if(x+dx>=0 && y+dy>=0 && x+dx<WA && y+dy<HA){
			int a=x+dx+WA*(y+dy);
			for(int j=0;j<adeptList[3][a].size();j++) if(!adept[3][adeptList[3][a][j]].isStar){
				//int sc2=sc1+(adept[3][adeptList[3][a][j]].isStar?1:0);
				//if(sc2>=2) continue;
				double radif=adept[3][adeptList[3][a][j]].ra-adept[0][i].ra;
				double decdif=adept[3][adeptList[3][a][j]].dec-adept[0][i].dec;
				double dist=radif*radif+decdif*decdif;
				if(dist>MINDISTSQ && dist<MAXDISTSQ){
					bool is2=false;
					bool is3=false;
					double ra2=adept[0][i].ra+radif/3;
					double dec2=adept[0][i].dec+decdif/3;
					double x2,y2;
					wc.convertRADEC2XY(wcs[1],ra2,dec2,x2,y2);
					int x2i=x2;
					int y2i=y2;
					int bestAdept2=-1;
					double minDist2=1e10;
					for(int xx2=x2i-MAXXY;xx2<=x2i+MAXXY;xx2++) for(int yy2=y2i-MAXXY;yy2<=y2i+MAXXY;yy2++) if(xx2>=0 && xx2<W && yy2>=0 && yy2<H){
						if(getAdept[1][xx2+W*yy2]>=0){
							int ai=getAdept[1][xx2+W*yy2];
							if(adept[1][ai].isStar) continue;
							//if(sc2+(adept[1][ai].isStar?1:0)>=2) continue;
							double dist2=(adept[1][ai].ra-ra2)*(adept[1][ai].ra-ra2)+(adept[1][ai].dec-dec2)*(adept[1][ai].dec-dec2);
							if(dist2<minDist2 && dist2<MAXRDSQ){
								minDist2=dist2;
								bestAdept2=ai;
								is2=true;
							}
						} 
					}
					double ra3=adept[0][i].ra+2*radif/3;
					double dec3=adept[0][i].dec+2*decdif/3;
					double x3,y3;
					wc.convertRADEC2XY(wcs[2],ra3,dec3,x3,y3);
					int x3i=x3;
					int y3i=y3;
					int bestAdept3=-1;
					double minDist3=1e10;
					for(int xx3=x3i-MAXXY;xx3<=x3i+MAXXY;xx3++) for(int yy3=y3i-MAXXY;yy3<=y3i+MAXXY;yy3++) if(xx3>=0 && xx3<W && yy3>=0 && yy3<H){
						if(getAdept[2][xx3+W*yy3]>=0){
							int ai=getAdept[2][xx3+W*yy3];
							if(adept[2][ai].isStar) continue;
							//if(sc2+(adept[2][ai].isStar?1:0)>=2) continue;
							double dist3=(adept[2][ai].ra-ra3)*(adept[2][ai].ra-ra3)+(adept[2][ai].dec-dec3)*(adept[2][ai].dec-dec3);
							if(dist3<minDist3 && dist3<MAXRDSQ){
								minDist3=dist3;
								bestAdept3=ai;
								is3=true;
							}
						} 
					}
					if(is2 && is3){
						// && ((adept[0][i].isStar?1:0)+(adept[1][bestAdept2].isStar?1:0)+(adept[2][bestAdept3].isStar?1:0)+(adept[3][adeptList[3][a][j]].isStar?1:0)<2)){
						double angleDif=0.0;
						if(checkCandidate(adept[0][i],adept[1][bestAdept2],adept[2][bestAdept3],adept[3][adeptList[3][a][j]],angleDif))
							candidate.push_back(Candidate(adept[0][i],adept[1][bestAdept2],adept[2][bestAdept3],adept[3][adeptList[3][a][j]],imageID,duration,season,zd,angleDif));
					}
				}
			}
		}
	}
	cerr << getTime()-st << " sec. pairing\n";st=getTime();
	cerr << candidate.size() << " candidates.\n";
		
	return candidate;
}

class AsteroidDetector{
  public:
	int trainingData(int width, int height, vector <int> imageData_1, vector <string> header_1, vector <double> wcs_1, vector <int> imageData_2, vector <string> header_2, vector <double> wcs_2, vector <int> imageData_3, vector <string> header_3, vector <double> wcs_3, vector <int> imageData_4, vector <string> header_4, vector <double> wcs_4, vector <string> detections){
		
		int neos=0;
		
		vector <Detection> detection(detections.size()/4);
		for(int i=0;i<detections.size()/4;i++){
			detection[i]=Detection(detections,4*i);
			if(detection[i].isNEO) neos++;
			//vysledok << number << "," << i;
			//for(int f=0;f<4;f++) vysledok << "," << detection[i].ra[f] << "," << detection[i].dec[f];
			//vysledok << "," << detection[i].isNEO << "\n";
		}
		
		if(detection.size()>1) for(int i=0;i<detection.size()-1;i++) for(int j=i+1;j<detection.size();j++){
			if(abs(detection[i].x[0]-detection[j].x[0])<NEARLIMIT && abs(detection[i].y[0]-detection[j].y[0])<NEARLIMIT
			   && abs(detection[i].x[1]-detection[j].x[1])<NEARLIMIT && abs(detection[i].y[1]-detection[j].y[1])<NEARLIMIT
			   && abs(detection[i].x[2]-detection[j].x[2])<NEARLIMIT && abs(detection[i].y[2]-detection[j].y[2])<NEARLIMIT
			   && abs(detection[i].x[3]-detection[j].x[3])<NEARLIMIT && abs(detection[i].y[3]-detection[j].y[3])<NEARLIMIT){
				detection[i].isTwo=true;
				detection[j].isTwo=true;
				totalTwo++;
			}
		}

		//ficury << number << "," << detection.size() << "," << neos;
		//for(int i=0;i<11;i++) for(int f=0;f<4;f++) ficury << "," << info[f][i];
		//ficury << "\n";
		
		//number++; if(number==100) return 2; else return 0;
				
		if(number==0){
			startTime=getTime();
			
			for(int j=0;j<PV;j++){
				VALUE[j].resize(CLUSTERS[j]+1);
				REPRE[j].resize(CLUSTERS[j],0);
				VALUE[j][0]=-1;
				VALUE[j][CLUSTERS[j]]=1.5;
				for(int i=1;i<CLUSTERS[j];i++){
					VALUE[j][i]=0.5;
				}
			}
			vector <int> clusterSize[PV];
			for(int j=0;j<PV;j++){
				clusterSize[j].resize(CLUSTERS[j],0);
			}
			for(int j=0;j<PV;j++){
				for(int k=0;k<CLUSTERS[j];k++){
					REPRE[j][k]=k;
				}
			}
		}
		/*int dmin1,dmin2,dmin3,dmin4,dmax1,dmax2,dmax3,dmax4;
		vector <int> median(4);
		getRange(imageData_1,dmin1,dmax1,median[0]);
		getRange(imageData_2,dmin2,dmax2,median[1]);
		getRange(imageData_3,dmin3,dmax3,median[2]);
		getRange(imageData_4,dmin4,dmax4,median[3]);
		*/
		
		vector <Candidate> candidate=getCandidates(width,height,imageData_1,header_1,wcs_1,imageData_2,header_2,wcs_2,imageData_3,header_3,wcs_3,imageData_4,header_4,wcs_4,"");
				
		int counter=0;
		int counterNEO=0;
		int counterTwo=0;
		for(int i=0;i<candidate.size();i++){
			int di=candidate[i].checkDetection(detection);
			cand.push_back(candidate[i]);
			
			if(candidate[i].isDetection){
				counter++;
			}
			if(candidate[i].isNEO){
				counterNEO++;
			}
			if(candidate[i].isTwo){
				counterTwo++;
			}
			#ifdef LOCAL
			subor << number << "," << candidate[i].isDetection << "," << di;
			for(int f=0;f<4;f++) subor << "," << candidate[i].ra[f] << "," << candidate[i].dec[f];
			for(int f=0;f<4;f++) subor << "," << candidate[i].x[f] << "," << candidate[i].y[f];
			for(int j=0;j<32;j++) subor << "," << candidate[i].data[j];
			subor << "," << candidate[i].score;
			subor << "," << candidate[i].isNEO;
			subor << "," << candidate[i].isTwo;
			subor << "\n";
			//if(i<10) candidate[i].visualize(number,i,imageData_1,dmin1,dmax1,imageData_2,dmin2,dmax2,imageData_3,dmin3,dmax3,imageData_4,dmin4,dmax4);
			#endif
		}
		foundDet+=counter;
		foundNEO+=counterNEO;
		totalDet+=detection.size();
		totalNEO+=neos;
		foundTwo+=counterTwo;

		cerr << counter << " of " << detection.size() << ".\n";
		
		
		number++;
		return 0;
		
	}
	int testingData(string imageID, int width, int height, vector <int> imageData_1, vector <string> header_1, vector <double> wcs_1, vector <int> imageData_2, vector <string> header_2, vector <double> wcs_2, vector <int> imageData_3, vector <string> header_3, vector <double> wcs_3, vector <int> imageData_4, vector <string> header_4, vector <double> wcs_4){
		if(testNumber==0){
			cerr << "Number of candidates: " << cand.size() << "\n";
			cerr << "Number of detections: " << totalDet << "\n";
			cerr << "Number of double det: " << totalTwo << "\n";
			cerr << "Number of NEOs      : " << totalNEO << "\n";
			cerr << "Number of matched detections: " << foundDet << "\n";
			cerr << "Number of matched NEOs      : " << foundNEO << "\n";
			cerr << "Number of matched double det: " << foundTwo << "\n";
			expectedNEOs=round(static_cast<double>(foundNEO)/5);
	  		expectedTwos=round(static_cast<double>(foundTwo)/5);
	  		cerr << "Number of expected NEOs     : " << expectedNEOs << "\n";
	  		#ifdef LOCAL
	  		vysledok << "Number of candidates: " << cand.size() << "\n";
			vysledok << "Number of detections: " << totalDet << "\n";
			vysledok << "Number of double det: " << totalTwo << "\n";
			vysledok << "Number of NEOs      : " << totalNEO << "\n";
			vysledok << "Number of matched detections: " << foundDet << "\n";
			vysledok << "Number of matched NEOs      : " << foundNEO << "\n";
			vysledok << "Number of matched double det: " << foundTwo << "\n";
			vysledok << "Number of expected NEOs     : " << expectedNEOs << "\n";
			#endif
			for(int j=0;j<PV;j++){
	    		randomTree[j].resize(TREES);
				for(int i=0;i<TREES;i++){
					if(getTime()-startTime>deadline*(j+1)) break;
	  				randomTree[j][i]=buildTree(j);
	  				counterTree[j]++;
					//#ifdef LOCAL
					if(i%10==0) cerr << getTime()-startTime << " sec. " << i+1 << " trees completed...\n";
					//#endif
				}
				for(int i=0;i<FEATURES;i++) cerr << i << ", " << featureScore[i] << "\n";
				#ifdef LOCAL
				for(int i=0;i<FEATURES;i++) vysledok << i << ", " << featureScore[i] << "\n";
				#endif
	    	}

		}
		vector <Candidate> candidate=getCandidates(width,height,imageData_1,header_1,wcs_1,imageData_2,header_2,wcs_2,imageData_3,header_3,wcs_3,imageData_4,header_4,wcs_4,imageID);
		for(int i=0;i<candidate.size();i++){
			#ifdef LOCAL
			candidate[i].temp=testNumber;
			#endif
			candTest.push_back(candidate[i]);
		}	
		testNumber++;
		return 0;
	}
	vector <string> getAnswer(){
		vector <string> ret;
		
		for(int j=0;j<PV;j++){
			electionCount[j].resize(CLUSTERS[j],0);
	    	for(int i=0;i<candTest.size();i++){
				vector <int> election(CLUSTERS[j],0);
	    		for(int t=0;t<counterTree[j];t++) election[randomTree[j][t].assignCluster(candTest[i],j)]++;
    			candTest[i].result[j]=0.0;
    			for(int cl=0;cl<CLUSTERS[j];cl++){
    				candTest[i].result[j]+=election[cl]*REPRE[j][cl];
    				electionCount[j][cl]+=election[cl];
    			}
    			candTest[i].result[j]/=counterTree[j];
    		}
    	}
		cerr << "Number of test candidates: " << candTest.size() << "\n";
		expectedNEOs=round(static_cast<double>(foundNEO)*candTest.size()/cand.size());
  		expectedTwos=round(static_cast<double>(foundTwo)*candTest.size()/cand.size());
  		cerr << "Number of expected NEOs     : " << expectedNEOs << "\n";
  		cerr << "Number of expected twos     : " << expectedTwos << "\n";
  		#ifdef LOCAL
  		vysledok << "Number of test candidates: " << candTest.size() << "\n";
		vysledok << "Number of expected NEOs     : " << expectedNEOs << "\n";
  		vysledok << "Number of expected twos     : " << expectedTwos << "\n";
		#endif
		sort(candTest.begin(),candTest.end(),sortByTwo);
		for(int i=0;i<candTest.size();i++){
			if(i<expectedTwos){
				candTest[i].isTwo=true;
			}
			else{
				candTest[i].isTwo=false;
			}
		}
		sort(candTest.begin(),candTest.end(),sortByNEO);
		double minNEO=1.0;
		for(int i=0;i<candTest.size();i++){
			if(i<expectedNEOs){
				candTest[i].isNEO=true;
				minNEO=min(candTest[i].result[0],minNEO);	
			}
			else{
				if(candTest[i].result[0]<=minNEO){
					candTest[i].isNEO=true;
					minNEO=min(candTest[i].result[0],minNEO);	
				}
				else candTest[i].isNEO=false;
			}
		}
		
		stable_sort(candTest.begin(),candTest.end());
		int ansCount=0;
		for(int i=0;i<candTest.size();i++){
			string neo=" 0";
			if(candTest[i].isNEO) neo=" 1";
			if(ansCount>=100000) break;
			ret.push_back(candTest[i].imageID+" "+SSTR(candTest[i].ra[0])+" "+SSTR(candTest[i].dec[0])+" "+SSTR(candTest[i].ra[1])+" "+SSTR(candTest[i].dec[1])+" "+SSTR(candTest[i].ra[2])+" "+SSTR(candTest[i].dec[2])+" "+SSTR(candTest[i].ra[3])+" "+SSTR(candTest[i].dec[3])+neo);
			ansCount++;
			#ifdef LOCAL
			vysledok << SSTR(candTest[i].temp)+","+SSTR(candTest[i].ra[0])+","+SSTR(candTest[i].dec[0])+","+SSTR(candTest[i].ra[1])+","+SSTR(candTest[i].dec[1])+","+SSTR(candTest[i].ra[2])+","+SSTR(candTest[i].dec[2])
			          +","+SSTR(candTest[i].ra[3])+","+SSTR(candTest[i].dec[3]) << "," << candTest[i].result[0] << "," << candTest[i].result[1] << "," << candTest[i].result[2] << "\n";
		    #endif
			if(candTest[i].result[2]>=0.03999){
				if(ansCount>=100000) break;
				ret.push_back(candTest[i].imageID+" "+SSTR(candTest[i].ra[0])+" "+SSTR(candTest[i].dec[0])+" "+SSTR(candTest[i].ra[1])+" "+SSTR(candTest[i].dec[1])+" "+SSTR(candTest[i].ra[2])+" "+SSTR(candTest[i].dec[2])+" "+SSTR(candTest[i].ra[3])+" "+SSTR(candTest[i].dec[3])+neo);
				ansCount++;
				#ifdef LOCAL
				//vysledok << SSTR(candTest[i].temp)+","+SSTR(candTest[i].ra[0])+","+SSTR(candTest[i].dec[0])+","+SSTR(candTest[i].ra[1])+","+SSTR(candTest[i].dec[1])+","+SSTR(candTest[i].ra[2])+","+SSTR(candTest[i].dec[2])
				//          +","+SSTR(candTest[i].ra[3])+","+SSTR(candTest[i].dec[3]) << "," << candTest[i].result[0] << "," << candTest[i].result[1] << "\n";
			    #endif
			}
		}
		
		
		return ret;
	}
};

