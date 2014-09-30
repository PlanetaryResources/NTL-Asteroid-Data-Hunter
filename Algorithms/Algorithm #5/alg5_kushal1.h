#define _CRT_SECURE_NO_WARNINGS
#ifdef OFFLINE_DEBUGGING
#define REFERENCE &
#else
#define REFERENCE
#endif
# define M_PI 3.14159265358979323846

#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <queue>
#include <map>
#include <cmath>
#include <assert.h>
#include <algorithm>
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock

using namespace std;

template <class T>
class Vector{
public:
	vector<T> dat;
	int N;
	Vector():N(2){
		dat = vector<T>(2,0);
	}
	Vector(int n):N(n){
		dat = vector<T>(N,0);
	}

	Vector(int n, T* data):N(n){
		dat = vector<T>(data,data+n);
	}

	Vector(T a, T b):N(2){
		dat.push_back(a);
		dat.push_back(b);
	}

	Vector(T a, T b, T c):N(3){
		dat.push_back(a);
		dat.push_back(b);
		dat.push_back(c);
	}

	template <class R>
	operator Vector<R>(){
		Vector<R> result(N);
		for(int i=0;i<N;i++){
			result.dat[i] = dat[i];
		}
		return result;
	}

	Vector& operator+=(const Vector &sv){
		assert(sv.N == N);
		for(int i=0;i<N;i++)
			dat[i] += sv.dat[i];
		return *this;
	}

	Vector& operator*=(const Vector &sv){
		assert(sv.N == N);
		for(int i=0;i<N;i++)
			dat[i] *= sv.dat[i];
		return *this;
	}

	Vector& operator*=(const T sv){
		for(int i=0;i<N;i++)
			dat[i] *= sv;
		return *this;
	}
	Vector& operator/=(const T sv){
		for(int i=0;i<N;i++)
			dat[i] /= sv;
		return *this;
	}

	Vector& operator-=(const Vector &sv){
		assert(sv.N == N);
		for(int i=0;i<2;i++)
			dat[i] -= sv.dat[i];
		return *this;
	}

	const Vector operator+(const Vector &sv) const{
		Vector result(*this);
		result += sv;
		return result;
	}
	
	const Vector operator-(const Vector &sv) const{
		Vector result(*this);
		result -= sv;
		return result;
	}
	
	const Vector operator*(const Vector &sv) const{
		Vector result(*this);
		result *= sv;
		return result;
	}

	const Vector operator*(const T sv) const{
		Vector result(*this);
		result *= sv;
		return result;
	}

	const Vector operator/(const T sv) const{
		Vector result(*this);
		result /= sv;
		return result;
	}

	bool operator==(const Vector &sv) const{
		for(int i=0;i<3;i++){
			if(dat[i] != sv.dat[i])
				return false;
		}
		return true;
	}
	
	bool operator!=(const Vector &sv) const{
		return !(sv == *this);
	}

	inline T& operator[](unsigned int i){ return dat[i];}

	double dot(const Vector &sv){
		T sum = 0;
		for(int i=0;i<N;i++){
			sum += dat[i]*sv.dat[i];
		}
		return sum;
	}

	Vector<T> homoVector() const{
		Vector result(*this);
		assert(N == 2);
		result.dat.push_back(1);
		result.N++;
		return result;
	}

	Vector<T> unHomoVector() const{
		Vector result(*this);
		assert(N == 3);
		result.dat.pop_back();
		result.N--;
		return result;
	}

	double normSq() const{
		double sum = 0;
		for(int i=0;i<N;i++){
			sum+=dat[i]*dat[i];
		}
		return sum;
	}

	double norm() const{
		return sqrt(normSq());
	}
};
template <class T>
class Matrix{
public:
	int R;
	int C;
	vector<T> mat;

	Matrix(int r, int c){
		R = r;
		C = c;
		mat.resize(R*C);
		setIdentity();
	}

	Matrix(int r, int c, T* data){
		R = r;
		C = c;
		mat = vector<T>(data, data+r*c);
	}

	void setIdentity(){
		for(int i=0;i<R;i++){
			for(int j=0;j<C;j++)
				mat[i*C+j] = i==j;
		}
	}
	
	Matrix& operator+=(const Matrix &sv){
		assert(sv.R == R && sv.C == C);
		for(int i=0;i<R*C;i++)
			mat[i] += sv.mat[i];
		return *this;
	}

	Matrix& operator*=(const T sv){
		for(int i=0;i<R*C;i++)
			mat[i] *= sv;
		return *this;
	}

	Matrix& operator/=(const T sv){
		for(int i=0;i<R*C;i++)
			mat[i] /= sv;
		return *this;
	}

	Matrix& operator-=(const Matrix &sv){
		assert(sv.R == R && sv.C == C);
		for(int i=0;i<R*C;i++)
			mat[i] -= sv.mat[i];
		return *this;
	}

	Matrix& operator*=(const Matrix &sv){
		assert(C == sv.R);
		vector<T> oldData = mat;
		
		mat = vector<T>(R*sv.C, 0);

		for(int i=0;i<R;i++){
			for(int j=0;j<sv.C;j++){
				for(int k=0;k<C;k++)
					mat[i*sv.C + j] += oldData[i*C + k]*sv.mat[k*sv.C + j];
			}
		}

		C = sv.C;

		return *this;
	}

	const Matrix operator+(const Matrix &sv) const{
		Matrix result(*this);
		result += sv;
		return result;
	}
	
	const Matrix operator-(const Matrix &sv) const{
		Matrix result(*this);
		result -= sv;
		return result;
	}
	
	const Matrix operator*(const Matrix &sv) const{
		Matrix result(*this);
		result *= sv;
		return result;
	}

	const Matrix operator*(const T sv) const{
		Matrix result(*this);
		result *= sv;
		return result;
	}

	const Matrix operator/(const T sv) const{
		Matrix result(*this);
		result /= sv;
		return result;
	}


	const Vector<T> operator*(const Vector<T> &sv) const{
		assert(sv.N == C);
		Vector<T> result = Vector<T>(R);
		for(int i=0;i<R;i++){
			for(int j=0;j<C;j++){
				result.dat[i] += mat[i*C + j]*sv.dat[j];
			}
		}
		return result;
	}

	bool operator==(const Matrix &sv) const{
		for(int i=0;i<R*C;i++){
			if(mat[i] != sv.mat[i])
				return false;
		}
		return true;
	}
	
	bool operator!=(const Matrix &sv) const{
		return !(sv == *this);
	}

	inline T* operator[](unsigned int i) const{ return (T*)mat.data()+i*C;}
	
	const Matrix inverse() const{
		assert(R==C && R>=2);

		Matrix result(*this);


		//Following code was adapted from users.erols.com/mdinolfo/matrix.htm
		for (int i=1; i < R; i++)
			result.mat[i] /= result.mat[0]; // normalize row 0
		for (int i=1; i < R; i++)  { 
			for (int j=i; j < R; j++)  { // do a column of L
				T sum = 0.0;
				for (int k = 0; k < i; k++)  
					sum += result.mat[j*R+k] * result.mat[k*R+i];
				result.mat[j*R+i] -= sum;
			}
			if (i == R-1)
				continue;
			for (int j=i+1; j < R; j++)  {  // do a row of U
				T sum = 0.0;
				for (int k = 0; k < i; k++)
					sum += result.mat[i*R+k]*result.mat[k*R+j];
				result.mat[i*R+j] = 
					(result.mat[i*R+j]-sum) / result.mat[i*R+i];
			}
		}
		for ( int i = 0; i < R; i++ ){  // invert L
			for ( int j = i; j < R; j++ )  {
				T x = 1.0;
				if ( i != j ) {
					x = 0.0;
					for ( int k = i; k < j; k++ ) 
						x -= result.mat[j*R+k]*result.mat[k*R+i];
					}
				result.mat[j*R+i] = x / result.mat[j*R+j];
			}
		}
		for ( int i = 0; i < R; i++ ){   // invert U
			for ( int j = i; j < R; j++ )  {
				if ( i == j ) continue;
				T sum = 0.0;
				for ( int k = i; k < j; k++ )
					sum += result.mat[k*R+j]*( (i==k) ? 1.0 : result.mat[i*R+k] );
				result.mat[i*R+j] = -sum;
			}
		}
		for ( int i = 0; i < R; i++ ){   // final inversion
			for ( int j = 0; j < R; j++ )  {
				T sum = 0.0;
				for ( int k = ((i>j)?i:j); k < R; k++ )  
					sum += ((j==k)?1.0:result.mat[j*R+k])*result.mat[k*R+i];
				result.mat[j*R+i] = sum;
			}
		}
		return result;
	}

	const Vector<T> solveAxEqualsB(Vector<T> b) const{
		assert(R==C && R>=2);
		assert(b.N == R);

		Matrix m(*this);
		Vector<T> x(b);

		vector<int> order(R,-1);
		vector<bool> used(R, false);
		for(int r=0;r<R;r++){
			T maxVal;
			int maxLoc = -1;
			for(int rr=0;rr<R;rr++){
				if(!used[rr] && (maxLoc == -1 || fabs(m.mat[rr*C + r]) > maxVal)){
					maxVal = fabs(m.mat[rr*C + r]);
					maxLoc = rr;
				}
			}
			assert(maxLoc != -1);
			order[r] = maxLoc;
			used[maxLoc] = true;

			for(int rr=0;rr<R;rr++){
				if(rr != maxLoc){
					T scaleFac = m.mat[rr*C + r] / m.mat[maxLoc*C + r];
					for(int c = r; c< C;c++){
						m.mat[rr*C + c] -= scaleFac*m.mat[maxLoc*C + c];
					}
					x.dat[rr] -= scaleFac*x.dat[maxLoc];
				}
			}
		}
		Vector<T> answer(b);
		for(int i=0;i<R;i++){
			answer.dat[i] = x[order[i]] / m.mat[order[i]*C + i];
		}
		return answer;
	}
		

	const Matrix<T> transpose() const{
		Matrix<T> result(C,R);
		for ( int i = 0; i < R; i++ ){   // final inversion
			for ( int j = 0; j < C; j++ )  {
				result.mat[j*R+i] = mat[i*C + j];
			}
		}
		return result;
	}
	T trace() const{
		T result = 0;
		for(int i=0; i<min(R,C);i++){
			result += mat[i*C + i];
		}
		return result;
	}
	T determinent() const{
		assert(R == C && R>= 2);
		if(R == 2){
			return mat[0]*mat[3] - mat[1]*mat[2];
		}
		int sign = 1;
		T det = 0;
		for(int i=0;i<R;i++){
			Matrix<T> subMat(R-1,C-1);
			int subR = 0;
			for(int r=0;r<R;r++){
				if(r != i){
					for(int c=1;c<C;c++){
						subMat.mat[subR*(C-1) + c-1] = mat[r*C + c];
					}
					subR++;
				}
			}
			det += sign*mat[i*C]*subMat.determinent();
			sign *= -1;
		}
		return det;
	}
};
class Transformation : private Matrix<double>{
public:
	Transformation():Matrix<double>(3,3){
	}

	operator Matrix<double>(){
		return (Matrix<double>)(*this);
	}

	Transformation(double m[2][2]):Matrix<double>(3,3){
		setIdentity();
		for(int i=0;i<2;i++){
			for(int j=0;j<2;j++){
					mat[i*C + j] = m[i][j];
			}
		}
	}

	Transformation(Matrix<double> m):Matrix<double>(m){
		assert(m.R == 3 && m.C == 3);
	}

	static const Transformation idenity(){
		Transformation t;
		return t;
	}

	static const Transformation translation(Vector<double> translation){
		assert(translation.N == 2);
		Transformation t;
		for(int i=0;i<2;i++)
			t.mat[i*t.C+2] = translation.dat[i];
		return t;
	}

	static const Transformation rotation(double theta){
		Transformation t;
		t.mat[0] = cos(theta);
		t.mat[1] = -sin(theta);
		t.mat[3] = sin(theta);
		t.mat[4] = cos(theta);
		return t;
	}

	static const Transformation dialation(Vector<double> scaleFac){
		assert(scaleFac.N >= 2);
		Transformation t;
		t.setIdentity();
		for(int i=0;i<2;i++)
			t.mat[i*t.C+i]*=scaleFac[i];
		return t;
	}
	const Transformation operator*(const Transformation &sv) const{
		return (Transformation)((Matrix)(*this)*(Matrix)sv);
	}

	const Vector<double> operator*(const Vector<double> &sv) const{
		assert(sv.N == 2);
		Vector<double> result = ((Matrix)*this) * sv.homoVector();
		assert(abs(result[2] - 1) < 0.001);
		return result.unHomoVector();
	}

	double* operator[](unsigned int i) const{ return (double*)mat.data()+i*3;}

	const Transformation inverse() const{
		return (Transformation)(((Matrix)(*this)).inverse());
	}
	
	Matrix<double> getMatrix() const{
		return (Matrix<double>) (*this);
	}
};
template <class T>
std::ostream &operator<<(std::ostream &os, Vector<T> const &v) { 
    os << "(";
	for(int i=0;i<v.N-1; i++)
		os << v.dat[i] << ", ";
	return os << v.dat[v.N-1] << ")";
}
template <class T>
std::ostream &operator<<(std::ostream &os, Matrix<T> const &m) { 
	for(int i=0;i<m.R;i++){
		for(int j=0;j<m.C;j++){
			if(j==0){
				if(i==0)
					os << "[ ";
				else
					os << "  ";
			}

			os << m[i][j];

			if(j == m.C - 1){
				if(i== m.R - 1)
					os << " ]";
				else
					os << endl;
			}
			else
				os << '\t';
			
			
		}
	}
    return os;
}

// element accesses are i, j
// vector accesses are <x, y>
template <class T>
class Image{
public:
	vector<T> p;
	const int w, h;
	const static int maxVal = 65535;
	Image(int _w, int _h, T initVal): w(_w), h(_h),p(_w*_h,initVal){
	}

	Image(int _w, int _h, vector<T>& data): w(_w), h(_h), p(data){
		assert(data.size() == _w*_h);
	}

	Image& operator+=(const Image &sv){
		assert(sv.w == w && sv.h == h);
		for(int i=0;i<h;i++){
			for(int j=0;j<w;j++)
				p[i*w + j] += sv.p[i*w + j];
		}
		return *this;
	}

	Image& operator*=(const Image &sv){
		assert(sv.w == w && sv.h == h);
		for(int i=0;i<h;i++){
			for(int j=0;j<w;j++)
				p[i*w + j] *= sv.p[i*w + j];
		}
		return *this;
	}

	Image& operator/=(const Image &sv){
		assert(sv.w == w && sv.h == h);
		for(int i=0;i<h;i++){
			for(int j=0;j<w;j++)
				p[i*w + j] /= sv.p[i*w + j];
		}
		return *this;
	}

	Image& operator-=(const Image &sv){
		assert(sv.w == w && sv.h == h);
		for(int i=0;i<h;i++){
			for(int j=0;j<w;j++)
				p[i*w + j] -= sv.p[i*w + j];
		}
		return *this;
	}

	void unApply(Transformation t, Image& result){
		for(int i=0;i<h;i++){
			for(int j=0;j<w;j++)
				result.p[i*w + j] = bilinear(t*Vector<double>(j+0.5,i+0.5));
		}
		//Get new image bounds
	}

	inline T* operator[](unsigned int i){ return (T*)&(p[i*w]); }
	inline T& operator[](Vector<int> v){ return p[v.dat[1]*w + v.dat[0]]; }
	const T safeRead(Vector<int> v){
		if(v.dat[0] < 0 || v.dat[0] >= w || v.dat[1] < 0 || v.dat[1] >= h){
			return 0;
		}
		else{
			return (*this)[v];
		}
	}

	const double bilinear(Vector<double> v){
		v *= 2;
		if(v.dat[0] < 1 || v.dat[1] < 1 || v.dat[0] > 2*w - 1 || v.dat[1] > 2*h - 1){
			return 0;
		}
		int x = (v.dat[0] - 1)/2;
		int y = (v.dat[1] - 1)/2;

		if(x == w - 1)
			x--;
		if(y == h - 1)
			y--;

		double x_cor = x + 0.5;
		double y_cor = y + 0.5;

		v*=0.5;

		//Bilinear Interpolation (en.wikipedia.org/wiki/Bilinear_interpolation)
		return p[y*w + x]*(x_cor+1-v.dat[0])*(y_cor+1-v.dat[1])
			+p[y*w + x+1]*(v.dat[0]-x_cor)*(y_cor+1-v.dat[1])
			+p[(y+1)*w + x]*(x_cor+1-v.dat[0])*(v.dat[1]-y_cor)
			+p[(y+1)*w + x+1]*(v.dat[0]-x_cor)*(v.dat[1]-y_cor);
	}
	void saveAs(int ind, string id){
		if(w > 1000){//quick fix for strange bug
			std::ostringstream filename;
			filename << ind << ' ' << id << ' ' << w << 'x' << h << ".raw";
			ofstream fout(filename.str(),ios_base::binary);
			int pixelVal;
			for(int i=0;i<w*h;i++){
				pixelVal = max(0,min((int)(p[i]),maxVal));
				fout.write((char*)&(pixelVal),2);
			}
			fout.close();
		}
	}
	void threshold(T val){
		for(int i=0;i<w*h;i++){
			if(p[i] > val)
				p[i] = maxVal;
			else
				p[i] = 0;
		}
	}
};
class Detection{
public:
	Vector<int> gridBase;
	Vector<double> pos;
	Vector<double> raDec;
	double time;
	double constOffset;
	double scaleFac;
	double meanSqError;
	const static double var;
	bool converged;
	Detection() {
		
	}
	Detection(Vector<int> aproxLoc, double aproxOffset, double aproxScaleFac, Image<int>& img): raDec(0,0), pos(((Vector<double>)aproxLoc) + Vector<double>(0.5,0.5)), gridBase(aproxLoc - Vector<int>(2,2)) {
		scaleFac = aproxScaleFac;
		constOffset = aproxOffset;
		converged = fineTune(img);
		meanSqError = calcError(img)/25.0;
	}
	bool fineTune(Image<int>& img){
		Vector<double> step(4);
		for(int i=0;i<10;i++){
			Matrix<double> J = jacobian();
			Matrix<double> JT = J.transpose();
			step = (JT*J).solveAxEqualsB(JT*residuals(img));
			constOffset -= step[0];
			scaleFac -= step[1];
			pos[0] -= step[2];
			pos[1] -= step[3];
		}
		return step.normSq() < 0.1;
	}
	Matrix<double> jacobian(){
		Matrix<double> result(25,4);
		for(int i=0;i<25;i++){
			double diffVect_i = i%5 + 0.5 + gridBase.dat[1] - pos[1];
			double diffVect_j = i/5 + 0.5 + gridBase.dat[0] - pos[0];
			result[i][0] = 1;
			result[i][1] = rawGaussian(diffVect_j,diffVect_i);
			double tmp = scaleFac * result[i][1] / var;
			result[i][2] = tmp*diffVect_j;
			result[i][3] = tmp*diffVect_i;
		}
		return result;
	}
	inline double calcError(Image<int>& img){
		return residuals(img).normSq();
	}
	Vector<double> residuals(Image<int>& img){
		Vector<double> residuals(25);
		int count = 0;
		for(int j=gridBase.dat[0];j<5+gridBase.dat[0];j++){
			for(int i=gridBase.dat[1];i<5+gridBase.dat[1];i++){
				residuals[count] = gaussian(j+0.5 - pos[0], i+0.5 - pos[1]) - img[i][j];
				count++;
			}
		}
		return residuals;
	}
	inline double gaussian(double x,double y){
		return constOffset + scaleFac*rawGaussian(x, y);
	}
	inline double rawGaussian(double x,double y){
		return (1.0/(2*M_PI*var))*exp(-(x*x + y*y)/(2*var));
	}
};
const double Detection::var = 1.05*1.05;
class DetectionSet{
public:
	Detection dets[4];
	int size;
	string imgId;

	//Quality factors
	double error;//displacements are in the same direction
	double avgLen; //Displacement in pixels/day
	double scaleFac[4];//Brightness
	double constOffset[4];//Background brightness

	DetectionSet(){
		size = false;
	}

	void addDet(Detection d){
		assert(size < 4);
		dets[size] = d;
		size++;
	}
	
	bool operator< (const DetectionSet &other) const{
		return badness() < other.badness();
	}

	double badness() const{
		return error;// + (avgLen/144.0) + ((scaleFac[0] + scaleFac[1] + scaleFac[2] + scaleFac[3]) / 60000.0);
	}
};
string trim(string s){
	//Coppied from codereview.stackexchange.com/questions/40124/trim-white-space-from-string
	const string whitespace = "\t\f\v  ";
	int start = s.find_first_not_of(whitespace);
	int end = s.find_last_not_of(whitespace);
	s.erase(0,start);
	s.erase((end - start) + 1);
	return s;
}
class Header{
public:
	map<string,string> dict;
	Header(vector<string> headers){
		for(int i = 0; i<headers.size();i++){
			int equalsPos = headers[i].find('=');
    		if(equalsPos != string::npos){
				int commentPos = headers[i].find('/');
				if(commentPos == -1){
					commentPos = headers[i].length();
				}
				dict[trim(headers[i].substr(0,equalsPos))] = trim(headers[i].substr(equalsPos+1,commentPos-equalsPos-1));
			}
		}
	}
	const string operator[](string key){ return dict[key]; }
};
class WCSTransformation{
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
        const double D2R;
        const double R2D;

        WCSTransformation(vector<double> wcs):
		D2R(M_PI/180.0), R2D(180.0/M_PI)
		{
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

        void update_invcd() {
            double inv_det = cd11*cd22 - cd12*cd21;
            inv_det = 1.0 / inv_det;
            invcd11 = cd22 * inv_det;
            invcd12 = -cd12 * inv_det;
            invcd21 = -cd21 * inv_det;
            invcd22 = cd11 * inv_det;
        }


        vector<double> sphs2x(double lng, double lat) {
            double coslat, coslat3, coslat4, coslng, dlng, dphi, sinlat, sinlat3, sinlat4, sinlng, x, y, z;
            double eul[5];
			vector<double> phiTheta(2);
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
            phiTheta[0] = fmod((eul[2] + dphi),360.0);
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

        vector<double> sphx2s(double phi, double theta) {
            double cosphi, costhe, costhe3, costhe4, dlng, dphi, sinphi, sinthe, sinthe3, sinthe4, x, y, z;
            vector<double> lngLat(2);
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
            if (abs(z) > 0.99) {
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


        vector<double> tans2x(double phi, double theta) {
            vector<double> xy(2);
            double cotan_theta = 1.0 / tan( D2R * theta );
            xy[0] = R2D * sin( D2R * phi ) * cotan_theta ;
            xy[1] = -R2D * cos( D2R * phi ) * cotan_theta;
            return xy;
        }

        vector<double> tanx2s(double x, double y) {
            vector<double> phiTheta(2);
            phiTheta[0] = R2D * atan2(x, -y);
            phiTheta[1] = R2D * atan2(R2D, sqrt(x*x + y*y) );
            return phiTheta;
        }

        Vector<double> convertRADEC2XY( Vector<double> raDec) {
            Vector<double> XY(2);
            vector<double> phiTheta = sphs2x(raDec[0], raDec[1]);
            vector<double> XXYY = tans2x(phiTheta[0], phiTheta[1]);
            XY[0] = invcd11*XXYY[0] + invcd12*XXYY[1] + crpix1;
            XY[1] = invcd21*XXYY[0] + invcd22*XXYY[1] + crpix2;
            return XY;
        }

        Vector<double> convertXY2RADEC( Vector<double> XY) {
            double dx = XY[0]-crpix1;
            double dy = XY[1]-crpix2;
            double xx = cd11*dx + cd12*dy;
            double yy = cd21*dx + cd22*dy;
            vector<double> phiTheta = tanx2s(xx, yy);
            vector<double> RD = sphx2s(phiTheta[0], phiTheta[1]);
			Vector<double> raDec(RD[0],RD[1]);
            return raDec;
        }
};


class AsteroidDetector{
	vector<DetectionSet> allAsteroids;


public:

	int trainingData(int width, int height, vector<int> REFERENCE imageData_1, vector<string> header_1, vector<double> wcs_1, vector<int> REFERENCE imageData_2, vector<string> header_2, vector<double> wcs_2, vector<int> REFERENCE imageData_3, vector<string> header_3, vector<double> wcs_3, vector<int> REFERENCE imageData_4, vector<string> header_4, vector<double> wcs_4, vector<string> detections){
		cerr << "called trainingData()" << endl;
		
		#ifdef TRAIN
		Header headers[4] = {Header(header_1),Header(header_2),Header(header_3),Header(header_4)};
		
		for(int i=0;i<detections.size();i++){
			fout << headers[i%4]["MJD"] << " " << detections[i] << endl;
		}
		return 0;
		#endif

		allAsteroids.clear();
		return 1;
	}

	int testingData(string imageID, int width, int height, vector<int> imageData_1, vector<string> header_1, vector<double> wcs_1, vector<int> imageData_2, vector<string> header_2, vector<double> wcs_2, vector<int> imageData_3, vector<string> header_3, vector<double> wcs_3, vector<int> imageData_4, vector<string> header_4, vector<double> wcs_4){
		// Setup my classes
		Image<int> images[4] = {Image<int>(width,height,imageData_1),Image<int>(width,height,imageData_2),Image<int>(width,height,imageData_3),Image<int>(width,height,imageData_4)};
		clear(imageData_1);
		clear(imageData_2);
		clear(imageData_3);
		clear(imageData_4);
		Header headers[4] = {Header(header_1),Header(header_2),Header(header_3),Header(header_4)};

		// Compute the perliminary tranformations		
		WCSTransformation WCSTransforms[4] = {WCSTransformation(wcs_1),WCSTransformation(wcs_2),WCSTransformation(wcs_3),WCSTransformation(wcs_4)};		
        
		// Find potential asteroids and align them
		vector<Detection> detections[4];
		Image<Detection*> detectionsImg[4] = {Image<Detection*>(width,height,NULL),Image<Detection*>(width,height,NULL),Image<Detection*>(width,height,NULL),Image<Detection*>(width,height,NULL)};
		for(int i=0;i<4;i++){
			//fillInStarCenters();
			findSources(images[i],detections[i]);
			for(int d=0;d<detections[i].size();d++)
				detections[i][d].pos = WCSTransforms[3].convertRADEC2XY(WCSTransforms[i].convertXY2RADEC(detections[i][d].pos));
			recomputeDetectionsImage(detections[i],detectionsImg[i]);
			//detectionsImg[i].saveAs(i,"detections");
			//images[i].saveAs(i,"orig");
		}
		

		int nStars = removeStars(detectionsImg);

		// Find asteroids by looking for movements
		vector<DetectionSet> asteroids;
		cerr << "Ratio = " << ( nStars / (double) detections[0].size() ) << endl;
		if( ( nStars / (double) detections[0].size() ) > 0.15){
			findAsteroids(detectionsImg, headers, asteroids);
			cerr << "found " << asteroids.size() << " asteroids" << endl;
			for(int i=0;i<asteroids.size();i++){
				asteroids[i].imgId = imageID;
				for(int j=0;j<4;j++){
					//asteroids[i].dets[j].raDec = XYToRaDec[3]*asteroids[i].dets[j].pos;
					asteroids[i].dets[j].raDec = WCSTransforms[3].convertXY2RADEC(asteroids[i].dets[j].pos);
					asteroids[i].dets[j].time = atof(headers[j]["MJD"].c_str());
				}
			}
			allAsteroids.insert(allAsteroids.end(),asteroids.begin(), asteroids.end());
		}
		else
			cerr << "Bad wcs data" << endl;
		return 1;
	}
	
	vector<string> getAnswer(){
		cerr << "called getAnswer()" << endl;
		vector<string> answer;
		
		sort(allAsteroids.begin(),allAsteroids.end());

		for(int i=0;i<2*allAsteroids.size();i++){
			//cerr << allDetections[i].badness() << endl;

			Vector<double> vel = (allAsteroids[i/2].dets[3].raDec -allAsteroids[i/2].dets[0].raDec)/(allAsteroids[i/2].dets[3].time -allAsteroids[i/2].dets[0].time);

			if(vel[1] < 0.25 && vel[1] > -0.6 && vel[0] < 0.5 && vel[0] > -0.5){

				std::ostringstream line;
				line << allAsteroids[i/2].imgId;
				for(int img=0;img<4;img++){
					for(int j=0;j<2;j++){
						 line << ' ' << allAsteroids[i/2].dets[img].raDec[j] + (i%2)*0.002;
					}
				}
			
				line << ' ' << (vel[0] < -0.4)?1:0;
				//cerr << "entery was: " << line.str() << endl;
				answer.push_back(line.str());
			}
		}
		return answer;
	}
private:
	void clear(vector<int>& v){
		vector<int> tmp;
		swap(v, tmp);
	}

	void transformDetectionsImage(vector<Detection>& detections, Image<Detection*>& detectionsImg, Transformation transformation){
		for(int i=0;i<detectionsImg.h;i++){
			for(int j=0;j<detectionsImg.w;j++){
				detectionsImg[i][j] = NULL;
			}
		}
		for(int j=0;j<detections.size();j++){
			detections[j].pos = transformation*detections[j].pos;
			if(detections[j].pos[0] >=0 && detections[j].pos[0] < detectionsImg.w && detections[j].pos[1] >=0 && detections[j].pos[1] < detectionsImg.h){
				detectionsImg[(Vector<int>)(detections[j].pos)] = &detections[j];
			}
		}
	}

	void recomputeDetectionsImage(vector<Detection>& detections, Image<Detection*>& detectionsImg){
		for(int i=0;i<detectionsImg.h;i++){
			for(int j=0;j<detectionsImg.w;j++){
				detectionsImg[i][j] = NULL;
			}
		}
		for(int j=0;j<detections.size();j++){
			if(detections[j].pos[0] >=0 && detections[j].pos[0] < detectionsImg.w && detections[j].pos[1] >=0 && detections[j].pos[1] < detectionsImg.h){
				detectionsImg[(Vector<int>)(detections[j].pos)] = &detections[j];
			}
		}
	}

	Transformation pixelLevelAlign(Image<Detection*>& base, Image<Detection*>& img){
		Vector<double> bestDisp;
		double bestDispMatches = -1;

		Image<int> counts(41,41,0);
		Vector<int> imgOffset(20,20);
		
		for(int pos_i=0; pos_i<base.h; pos_i++){
			for(int pos_j=0; pos_j<base.w; pos_j++){
				if(base[pos_i][pos_j]){
					for(int disp_i = -16; disp_i <= 16; disp_i++){
						for(int disp_j = -16; disp_j <= 16; disp_j++){
							int i = pos_i - disp_i;
							int j = pos_j - disp_j;
							if(i >= 0 && j>=0 && i<img.h && j<img.w && img[i][j])
								counts[disp_i+imgOffset[0]][disp_j+imgOffset[1]]++;
						}
					}
				}
			}
		}
		
		Vector<int> maxLoc;
		int maxVal = -1;
		for(int i=0;i<counts.h;i++){
			for(int j=0;j<counts.w;j++){
				if(counts[i][j] > maxVal){
					maxVal = counts[i][j];
					maxLoc = Vector<int>(j,i) - imgOffset;
				}
			}
		}
		return Transformation::translation(maxLoc);
	}

	double median(vector<double>& data){
		sort(data.begin(),data.end());
		return data[data.size()/2];
	}

	double approxMedian(vector<double>& data, int downsampleFactor){
		vector<int> newData(data.size()/downsampleFactor);
		for(int i=0;i<data.size()/downsampleFactor;i++){
			newData[i] = data[i*downsampleFactor];
		}
		sort(newData.begin(),newData.end());
		return newData[newData.size()/2];
	}

	Transformation subPixelLevelAlign(Image<Detection*>& base, Image<Detection*>& img){
		vector<double> xDisps;
		vector<double> yDisps;
		for(int i=0;i<base.h;i++){
			for(int j=0;j<base.w; j++){
				if(base[i][j]){
					double bestDist = 1e9;
					Vector<double> bestDisplacement;
					bool found = false;
					for(int di=-2;di<=2;di++){
						for(int dj=-2;dj<=2;dj++){
							if(img.safeRead(Vector<int>(j+dj,i+di))){
								Vector<double> displacement = base[i][j]->pos - img[i+di][j+dj]->pos;
								if( displacement.normSq() < bestDist){
									bestDist = displacement.normSq();
									bestDisplacement = displacement;
									found = true;
								}
							}
						}
					}
					if(found){
						xDisps.push_back(bestDisplacement[0]);
						yDisps.push_back(bestDisplacement[1]);
					}
				}
			}
		}
		Transformation result = Transformation::translation(Vector<double>(median(xDisps), median(yDisps)));
		return result;
	}

	int removeStars(Image<Detection*> dets[4]){
		int nStars = 0;
		for(int i=0;i<dets[0].h;i++){
			for(int j=0;j<dets[0].w; j++){
				if(dets[3][i][j]){
					Vector<int> objects[4] = {Vector<int>(-1,-1), Vector<int>(-1,-1), Vector<int>(-1,-1), Vector<int>(j,i)};
					bool star = true;
					for(int d=0;d<3;d++){
						double dist = 1e9;
						for(int di=-1;di<=1;di++){
							for(int dj=-1;dj<=1; dj++){
								if(dets[d].safeRead(Vector<int>(j+dj,i+di))){
									if( (dets[3][objects[3]]->pos - dets[d][i+di][j+dj]->pos).normSq() < dist){
										dist = (dets[3][objects[3]]->pos - dets[d][i+di][j+dj]->pos).normSq();
										objects[d] = Vector<int>(j+dj, i+di);
									}
								}
							}
						}
						//cerr << dist << endl;
						if(dist > 1){
							star = false;
						}
					}
					if(star){
						nStars++;
						for(int d=0;d<4;d++){
							dets[d][objects[d]] = 0;
						}
					}
				}
			}
		}
		#ifdef OFFLINE_DEBUGGING
		for(int i=0;i<4;i++){
			//dets[i].saveAs(i, "detections (no stars)");
		}
		#endif
		return nStars;
	}


	void findAsteroids(Image<Detection*> dets[4], Header headers[4], vector<DetectionSet>& outputDetections){
		static double times[4] = {atof(headers[0]["MJD"].c_str()),atof(headers[1]["MJD"].c_str()),atof(headers[2]["MJD"].c_str()),atof(headers[3]["MJD"].c_str())};
		static double timeDiffs[3] = {times[1]-times[0],times[2]-times[1],times[3]-times[2]};
	
		static Vector<double> asteroids[4];
		static Vector<int> inds[4];

		for(inds[0][0] = 0; inds[0][0] < dets[0].w; inds[0][0]++){
			for(inds[0][1] = 0; inds[0][1] < dets[0].h; inds[0][1]++){
				if(inds[0][0] == 1807 && inds[0][1] == 3402){
					//break point here
					//cerr << "found" << endl;
				}
				if(dets[0].safeRead(inds[0])){// #1 Asteroid
					asteroids[0] = ((Detection*)dets[0][inds[0]])->pos;
					for(inds[1][0] = inds[0][0] - 5; inds[1][0] <= inds[0][0] + 5; inds[1][0]++){
						for(inds[1][1] = inds[0][1] - 5; inds[1][1] <= inds[0][1] + 5;inds[1][1]++){
							if(dets[1].safeRead(inds[1])){// #2 Asteroid
								asteroids[1] = ((Detection*)dets[1][inds[1]])->pos;
								Vector<double> asteroid2SearchCenter = asteroids[1] + (asteroids[1] - asteroids[0])*timeDiffs[1]/timeDiffs[0];
								for(inds[2][0]=asteroid2SearchCenter[0] - 1; inds[2][0] <= asteroid2SearchCenter[0] + 1; inds[2][0]++){
									for(inds[2][1]=asteroid2SearchCenter[1] - 1; inds[2][1] <= asteroid2SearchCenter[1] + 1; inds[2][1]++){
										if(dets[2].safeRead(inds[2])){// #3 Asteroid
											asteroids[2] = ((Detection*)dets[2][inds[2]])->pos;
											Vector<double> asteroid3SearchCenter = asteroids[2] + (asteroids[2] - asteroids[0])*timeDiffs[2]/(timeDiffs[0]+timeDiffs[1]);
											for(inds[3][0] = asteroid3SearchCenter[0] - 1; inds[3][0] <= asteroid3SearchCenter[0] + 1; inds[3][0]++){
												for(inds[3][1] = asteroid3SearchCenter[1] - 1; inds[3][1] <= asteroid3SearchCenter[1] + 1; inds[3][1]++){
													if(dets[3].safeRead(inds[3])){// #4 Asteroid
														asteroids[3] = ((Detection*)dets[3][inds[3]])->pos;
														
														vector< Vector<double> > normalizedMovements;
														double avgLen = 0;
														for(int i=0;i<3;i++){
															normalizedMovements.push_back((asteroids[i+1] - asteroids[i])/timeDiffs[i]);
															avgLen += normalizedMovements[i].norm();
														}
														avgLen /= 3;
														for(int i=0;i<3;i++){
															normalizedMovements[i] /= avgLen;
														}
														double errors = 0;
														for(int i=0;i<3;i++){
															errors += (normalizedMovements[i] - normalizedMovements[(i+1)%3]).norm();
														}
														
														//cerr << "ERROR: " << errors <<',' << avgLen << endl;

														if(errors < 2.5 && avgLen > 70){//TODO better filter
															//cerr << inds[0][0] << " " << inds[0][1] << endl << inds[3][0] << " " << inds[3][1] << endl;
															DetectionSet result;
															for(int i=0;i<4;i++){
																result.addDet(*(Detection*)dets[i][inds[i]]);
																result.scaleFac[i] = result.dets[i].scaleFac;
																result.constOffset[i] = result.dets[i].constOffset;
															}
															result.error = errors;
															result.avgLen = avgLen;
															outputDetections.push_back(result);
														}
													}
												}
											}
										}
									}
								}
								
							}
						}
					}
				}
			}
		}
	}

	void findSources(Image<int>& img, vector<Detection>& detections){
		// Consider finding a better point spread function (it would improve detections, and help a whole lot with the gauss newton fitting)

		double gaussian[5][5];

		Image<double> scaleFacImage(img.w, img.h,0);
		Image<double> offsetImage(img.w, img.h,0);
		Image<double> potential(img.w, img.h,0);
		Image<double> mask(img.w, img.h,0);

		//This is to efficiently compute the best offset and scaling for the gaussian psf we are fitting
		//It was derived using calculus.
		double sumOfG = 0;
		double sumOfG2 = 0;
		double sumOfVG = 0;
		double sumOfV = 0;
		double sumOf1 = 25;
		
		for(int offsetX = -2; offsetX<=2; offsetX++){
			for(int offsetY = -2; offsetY<=2; offsetY++){
				gaussian[offsetX+2][offsetY+2] = (1.0/(2*M_PI*Detection::var))*exp(-(offsetX*offsetX+offsetY*offsetY)/(2*Detection::var));
				sumOfG += gaussian[offsetX+2][offsetY+2];
				sumOfG2 += gaussian[offsetX+2][offsetY+2] * gaussian[offsetX+2][offsetY+2];
			}
		}
		for(int i=2;i<img.h-2;i++){
			for(int j=2;j<img.w-2;j++){
				sumOfV = 0;
				sumOfVG = 0;
				for(int offsetX = -2; offsetX<=2; offsetX++){
					for(int offsetY = -2; offsetY<=2; offsetY++){
						double v = img[i+offsetY][j+offsetX];
						sumOfV += v;
						sumOfVG += v * gaussian[offsetX+2][offsetY+2];
					}
				}

				//This comes from solving a system of two equations
				//Old code explains what I did, the new code below is optimized
					//double mat[2][2] = {{sumOf1, sumOfG},{sumOfG,sumOfG2}};
					//Vector<double> result = Matrix<double>(2,2,(double*)mat).solveAxEqualsB(Vector<double>(sumOfV, sumOfVG));
					//double constOffset = result[0];
					//double scaleFac = result[1];
					//offsetImage.p[i*img.w + j] = constOffset;
					//scaleFacImage.p[i*img.w + j] = scaleFac;

				double det = sumOf1*sumOfG2 - sumOfG*sumOfG;				
				offsetImage.p[i*img.w + j] = (sumOfG2*sumOfV - sumOfG*sumOfVG)/det;
				scaleFacImage.p[i*img.w + j] = (-sumOfG*sumOfV+sumOf1*sumOfVG)/det;
				potential.p[i*img.w + j] = scaleFacImage.p[i*img.w + j]/sqrt(offsetImage.p[i*img.w + j]);
			}
		}

		double threshold = 0.9*(2*M_PI*1.05*1.05)*70.0/53.3;//70 was the original cutoff, 53.3 accounts for the dividing by the sqrt of brightness, the pi factor takes care of the scalling factor that I changed in the gaussian. The leftmost factor is a correction to try to improve this threshold amount
		double maskThreshold = 0.7*threshold;
		Vector<int> dirs[4] = {Vector<int>(-1,0),Vector<int>(1,0),Vector<int>(0,-1),Vector<int>(0,1)};
		for(int i=2;i<img.h-2;i++){
			for(int j=2;j<img.w-2;j++){
				if(potential[i][j] > maskThreshold){
					vector< Vector<int> > points;
					int minX = j, maxX = j, minY = i, maxY = i;
					Vector<int> startPix(j,i);
					potential[startPix] = 0;
					queue< Vector<int> > q;
					q.push(startPix);
					points.push_back(startPix);
					while(q.size()){
						for(int i=0;i<4;i++){
							Vector<int> newPix = q.front() + dirs[i];
							if(newPix[0] >= 2 && newPix[0] < img.w-2 && newPix[1] >= 2 && newPix[1] < img.h-2){
								if(potential[newPix] > maskThreshold){
									
									minX = min(minX, newPix[0]);
									maxX = max(maxX, newPix[0]);
									minY = min(minY, newPix[1]);
									maxY = max(maxY, newPix[1]);

									potential[newPix] = 0;
									
									points.push_back(newPix);
									q.push(newPix);
								}
							}
						}
						q.pop();
					}
					if(abs((maxY - minY) - (maxX - minX))  > 3 || points.size() > 16){
						for(int i=0;i<points.size();i++){
							mask[points[i]] = Image<int>::maxVal;
						}
					}
				}
			}
		}

		//Image<int> allPoints(img.w, img.h, 0);
		for(int i=2;i<img.h-2;i++){
			for(int j=2;j<img.w-2;j++){
				if(scaleFacImage[i][j]/sqrt(offsetImage[i][j]) > threshold){
					int count = 1;
					int minX = j, maxX = j, minY = i, maxY = i;
					Vector<int> bestPix(j,i);
					int bestPixVal = scaleFacImage[bestPix];
					scaleFacImage[bestPix] = 0;
					queue< Vector<int> > q;
					q.push(bestPix);
					while(q.size()){
						for(int i=0;i<4;i++){
							Vector<int> newPix = q.front() + dirs[i];
							if(newPix[0] >= 2 && newPix[0] < img.w-2 && newPix[1] >= 2 && newPix[1] < img.h-2){
								if(scaleFacImage[newPix]/sqrt(offsetImage[newPix]) > threshold){
									count++;
									minX = min(minX, newPix[0]);
									maxX = max(maxX, newPix[0]);
									minY = min(minY, newPix[1]);
									maxY = max(maxY, newPix[1]);
									if(scaleFacImage[newPix] > bestPixVal){
										bestPixVal = scaleFacImage[newPix];
										bestPix = newPix;
									}
									scaleFacImage[newPix] = 0;
									q.push(newPix);
								}
							}
						}
						q.pop();
					}
					if(abs((maxY - minY) - (maxX - minX))  <= 2){
						if(count <= 16){
							if(!mask[bestPix]){
								Detection provDet = Detection(bestPix,offsetImage[bestPix],bestPixVal,img);
								//allPoints[bestPix] = Image<int>::maxVal;
								if(provDet.converged && provDet.scaleFac/sqrt(provDet.constOffset) > threshold /*&& noise is low*/)
									detections.push_back(provDet);
							}
						}
					}
				}
			}
		}
		
		#ifdef OFFLINE_DEBUGGING
		//static int count = 0;
		//allPoints.saveAs(count, "All Points");
		//count++;
		//count %= 4;
		#endif
	}
};


#ifdef OFFLINE_DEBUGGING
int main(){
	AsteroidDetector a;
	int W, H, N;
	char buffer[100000];
	vector<int> imageData_[4];
	vector<string> header_[4];
	vector<double> wcs_[4];

#ifdef TRAIN    
	for (int i=0; i < 100; i++)
    {
        
		gets(buffer);
		W = atoi(buffer);
		gets(buffer);
		H = atoi(buffer);
        for (int f=0; f < 4; f++)
        {
			imageData_[f].resize(W*H);
            for (int j=0; j < W*H; j++){
				gets(buffer);
				imageData_[f][j] = atoi(buffer);
			}
			
            gets(buffer);
			N = atoi(buffer);
			header_[f].resize(N);
			cerr << N << endl;
			for (int j=0; j < N; j++){
				gets(buffer);
				header_[f][j] = string(buffer);
			}

			wcs_[f].resize(8);
            for (int j=0; j < 8; j++){
				gets(buffer);
				wcs_[f][j] = atof(buffer);
			}
        }
        gets(buffer);
		N = atoi(buffer);
		vector<string> detections(N);
        for (int j=0; j < N; j++){
			gets(buffer);
			detections[j] = string(buffer);
		}
		int result = a.trainingData(W, H, imageData_[0], header_[0], wcs_[0], imageData_[1], header_[1], wcs_[1], imageData_[2], header_[2], wcs_[2], imageData_[3], header_[3], wcs_[3], detections);
        cout << result << endl;
        if (result==1) break;
    }
#endif
    for (int i=1; i <= 20; i++)
    {
		gets(buffer);
		string imageID(buffer);
		gets(buffer);
		W = atoi(buffer);
		gets(buffer);
		H = atoi(buffer);

        for (int f=0; f < 4; f++)
        {
			imageData_[f].resize(W*H);
			for (int j=0; j < W*H; j++){
				gets(buffer);
				imageData_[f][j] = atoi(buffer);
			}

			gets(buffer);
			N = atoi(buffer);
			header_[f].resize(N);
			for (int j=0; j < N; j++){
				gets(buffer);
				header_[f][j] = string(buffer);
			}

			wcs_[f].resize(8);
            for (int j=0; j < 8; j++){
				gets(buffer);
				wcs_[f][j] = atof(buffer);
			}
        }
		cerr << "Doing #" << i << endl;
        int result = 0;
		//if(i==1)
			result = a.testingData(imageID, W, H, imageData_[0], header_[0], wcs_[0], imageData_[1], header_[1], wcs_[1], imageData_[2], header_[2], wcs_[2], imageData_[3], header_[3], wcs_[3]);
        cout << result << endl;
    }
    vector<string> results = a.getAnswer();
	cout << results.size() << endl;
    for (int i=0; i<results.size(); i++)
        cout << results[i] << endl;
	return 0;
}
#endif

