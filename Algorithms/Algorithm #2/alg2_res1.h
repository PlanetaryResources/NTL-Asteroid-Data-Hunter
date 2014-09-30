//rks: My file combiner strips comments except ones that start with //rks:
#pragma once
#include <string>
#include <vector>
#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <unordered_map>
#include <exception>
using namespace std;
#include <limits>
#include <tuple>
#include <exception>
#include <stdexcept>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-overflow"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#include <string>
#include <cstdio>
#include <cstring>
#include <cstdarg>
#define MAX_LOG_MESSAGE_LENGTH 800
#define MAX_TOTAL_LOG_MESSAGE_LENGTH 4096
#define VA_COPY_BUF_A char buf1[MAX_TOTAL_LOG_MESSAGE_LENGTH]; va_list vl; va_start(vl, format); vsprintf(buf1, format, vl); va_end(vl);
class Logd
{
private:
	string name;
	void FormatLogMessageA(const char* levelName, string& name, const char* str, char* buf)
	{
		sprintf(buf, "Logd: ");
		strncat(buf, levelName,255);
		strncat(buf, ": ",255);
		strncat(buf, name.c_str(),255);
		strncat(buf, ": ",255);
		strncat(buf, str, MAX_LOG_MESSAGE_LENGTH);
		strncat(buf, "\n",256);
	}
	void outSub(const char* levelName, const char* str)
	{
		char buf[MAX_TOTAL_LOG_MESSAGE_LENGTH];
		FormatLogMessageA(levelName, name, str, buf);
		fprintf(stderr, buf);
		fflush(stderr);
	}
public:
	Logd()
	{
		name = "Main";
	}
	int printf(const char* format, ...)
	{
		VA_COPY_BUF_A;
		fprintf(stderr, buf1);
		return 1;
	}
	int debug(const char* format, ...)
	{
		VA_COPY_BUF_A;
		outSub("Debug", buf1);
		return 1;
	}
	int info(const char* format, ...)
	{
		VA_COPY_BUF_A;
		outSub("Info", buf1);
		return 1;
	}
	int warn(const char* format, ...)
	{
		VA_COPY_BUF_A;
		outSub("Warn", buf1);
		return 1;
	}
	int error(const char* format, ...)
	{
		VA_COPY_BUF_A;
		outSub("Error", buf1);
		return 1;
	}
};
extern Logd logd;
Logd logd;
string resultFileTag = "unset";
#define PRODUCTION
#ifndef PRODUCTION
#define LOCAL
#define USE_IPP
#define USE_IPP_REGISTER
#endif
#ifndef LOCAL
#define ErrorExit(x,...) { char rthstr[512]; sprintf(rthstr, x, ##__VA_ARGS__); logd.error(rthstr); }
#endif
#define DO_CROP_INPUT_IMAGE true
#define REGISTER_N_DOWN_SAMPLE 3
#define DO_SUBPIXEL_SHIFT true
#define DO_IMAGE_REGISTER false
#define DO_IMAGE_AFFINE_ALIGN true
#define IMAGE_AFFINE_ALIGN_ADD_SCALE 0.00033f
#define SBO_BLOB_MIN_AREA 150
#define CONTOUR_DS_SEG_SIZE 32
#define NSI_NSTDDEV_MIN (2.0f)
#define NSI_NSTDDEV_MAX (12.0f)
#define NSI_DIFF_FLOOR 5
#define NSI_DIFF_MIN_NEIGHBOR_COUNT 3
#define NSI_DIFF_MIN_NEIGHBOR_MEAN 20
#define NSI_DIFF_MIN_NEIGHBOR_MAX 20
#define NSI_MEDIAN_COMPENSATION_SCALE 1.4f
#define DO_8_CONNECTED_NSI_BLOBS true
#define MIN_NSI_BLOB_PREFILTER_ROUNDNESS_SCORE 0.2f
#define MAX_NSI_BLOB_DIAMETER_PIX 16
#define NSI_BLOB_INDEX_GRID_DIM_PIX 32
#define MIN_NSI_BLOB_VELOCITY_PIXELS_PER_IMG 1
#define MAX_NSI_BLOB_VELOCITY_PIXELS_PER_IMG 8
#define BLOB_SIM_SCORE_MIN_THRESHOLD 0.01f
#define BLOB_DIST_MAX_THRESHOLD_PIX 2.5f
#define MIN_N_BLOBS_PER_CA (2.5f)
#define REMOVE_DET_TRUTH_DIST_THRESHOLD 16
#define CA_SBO_RADIUS_FACTOR	1.2f
#define CA_PROX_SEARCH_RADIUS_PIX	64
#ifdef LOCAL_JAVA_RUN
#define TIME_LIMIT_MIN_TRAIN_ACCUM 9999950.0f
#define TIME_LIMIT_MIN_TRAIN_RF 999998.0f
#define TIME_LIMIT_MIN_ALL 9999975.0f
#else
#define TIME_LIMIT_MIN_TRAIN_ACCUM 45.0f
#define TIME_LIMIT_MIN_TRAIN_RF 10.0f
#define TIME_LIMIT_MIN_ALL 74.0f
#endif
#ifndef PRODUCTION
#define MAX_CAS_PER_IMAGE_SET_TO_USE_RF_TO_TOP_N 100000
#define MAX_CAS_PER_IMAGE_SET_TRAIN 500
#define MAX_CAS_PER_IMAGE_SET_TEST 2500
#else
#define MAX_CAS_PER_IMAGE_SET_TO_USE_RF_TO_TOP_N 100000
#define MAX_CAS_PER_IMAGE_SET_TRAIN 2000
#define MAX_CAS_PER_IMAGE_SET_TEST 20000
#endif
#define MIN_DETECTIONS_FOR_TRAINING	1
#define DO_PREFILTER_CLASSIFIER_IN_TEST false
#define RF_N_TREES 500
#define EXCLUDE_ZERO_DET_SETS_FROM_TRAIN true
#define DO_NEO_CLASSIFIER true
#define DO_NEO_BLOCKERS false
#define N_NEOS_TO_SELECT_PER_IMAGE_SET (0.5f)
#define NEO_SPEED_THRESHOLD_PIX (99996.0f)
#define MAX_N_SPEED_NEOS_TO_SELECT_PER_IMAGE_SET (0.5f)
#define BASE_INDEX_FOR_DUP_CAS 1000000
#define DO_OVERLAP_CLASSIFIER true
#define P_OVERLAP_THRESHOLD 0.075f
#define MAX_OVERLAP_DUPS_TO_ADD	10000
#define DO_SLOW_MOVERS true
#define BASE_BLOB_INDEX_FOR_SLOW_MOVERS 2000000
#define SLOW_MOVER_MIN_TOTAL_MOVE_DIST_PIX 2.0f
#define SLOW_MOVERS_MIN_NEIGHBOR_PIX 4
#define SLOW_MOVERS_MIN_NEIGHBOR_PIX_MEAN 25
#define SLOW_MOVERS_MIN_NEIGHBOR_PIX_MAX 35
#define MAX_ANSWERS 100000
#define DET_TRUTH_ALIGN_THRESHOLD_PIX 8
#define SCORING_DIST_THRESHOLD_DEG 0.001
#define DEBUG_COORDS_TOLERANCE_PIX 4
#define CX 31
#define CY 31
#define IDIM 64
#define ILEN 4096
#define IBYTES 8192
#define BG_FRACTION (0.65f) 
#define FG_FRACTION (0.01f) 
#define UNDEF_FEATURE (-999.0f)
#define MIN_FEAT_CLUSTER_COUNT (2)
#define DO_CC_CLUSTER true
#define CENTER_TOL 2
#define CENTER_DESCENT_TOL_SHIFT 2
#define MEAN_PX_PER_DEGREE (1440)
#define CENTER_MASK_RADIUS 9
#define CSV_FEAT(x) (x == x ? (x == std::numeric_limits<float>::infinity() ? UNDEF_FEATURE : x) : UNDEF_FEATURE)
const int RegionIntensityLevels[] = {0, 100, 1000, 1500, 2200, 5000, 22000, 32767, 45000, 50000, 60000, 65000, 65535, 65536};
#ifdef USE_MT
#include <thread>
#include <mutex>
#endif
#include <string>
#include <vector>
#include <algorithm>
#include <cassert>
#include <limits>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cassert>
#include <cfloat>
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (a > b ? b : a)
#define SIGN(a) (a >= 0 ? 1 : -1)
#define SMAX(a,b) ((a) > 0 ? MAX(a,b) : (-(a) > (b) ? a : -(b)))
#define SMIN(a,b) ((a) > 0 ? MIN(a,b) : (-(a) < (b) ? a : -(b)))
#define PI_DOUBLE (3.14159265358979323846)
#define PI_FLOAT (3.14159265358979f)
#define HALF_PI_FLOAT (1.570796370506f)
#define QUARTER_PI_FLOAT (0.78539818525f)
#define TWO_PI_FLOAT (6.28318530717959f)
#define SQRT_TWO_PI_FLOAT (2.506628274631f)
#define E_FLOAT (2.71828182845904f)
#define CLIP(x,min,max) MIN(MAX(x, min), max)
#define CLOSE_ENOUGH(x,y,e) ((x >= y-e) && (x <= y+e))
#define roundd(dbl) ((int)floor(dbl + 0.5))
#define roundf(x) ((int)floor(x + 0.5f))
#ifndef byte
typedef unsigned char byte;
#endif
struct LeastSquaresFitResult
{
	float a; 
	float b; 
	float stddev; 
	float r2; 
	LeastSquaresFitResult()
	{
		a = b = stddev = r2 = FLT_MAX;
	}
	std::string ToString()
	{
		char tmp[1024];
		if (a == FLT_MAX)
		{
			sprintf(tmp, "unset");
		}
		else
		{
			sprintf(tmp, "a: %.3f, b: %.3f, stddev: %.3f, r2: %.3f", a, b, stddev, r2);
		}
		return (string)tmp;
	}
};
float degreesToRadiansF(float d)
{
	return d * PI_FLOAT / 180.0f;
}
float radiansToDegreesF(float d)
{
	return d * 180.0f / PI_FLOAT;
}
bool checkIsNotNaN(float v)
{
	return (v == v);
}
bool checkIsNormalFloat(float v)
{
	if (v != v) return false; 
	if (v == std::numeric_limits<float>::infinity()) return false;
	if (v == -std::numeric_limits<float>::infinity()) return false;
	return true;
}
float rangeScore(float lowEnd, float highEnd, float score)
{
	if (highEnd <= lowEnd) ErrorExit("rangeScore: inverted bounds");
	float dx = highEnd - lowEnd;
	float dv = score - lowEnd;
	float ranged = dv / dx;
	ranged = CLIP(ranged, 0.0f, 1.0f);
	return ranged;
}
void rotatePoints(float* x, float* y, int n, float angleRadians)
{
	float s = sin(angleRadians);
	float c = cos(angleRadians);
	for (int i = 0; i < n; i++)
	{
		x[i] = x[i] * c - y[i] * s;
		y[i] = x[i] * s + y[i] * c;
	}
}
void leastSquaresFitWeighted(float* x, float* y, float* weights, int n, bool doRotateForR2, LeastSquaresFitResult& result, float* dys)
{
	assert(n > 1);
	result.a = result.b = result.stddev = result.r2 = 0;
	double sumX = 0;
	double sumY = 0;
	double sumWeight = 0;
	float minX = FLT_MAX, maxX = -FLT_MAX, minY = FLT_MAX, maxY = -FLT_MAX;
	for (int i = 0; i < n; i++)
	{
		sumWeight += weights[i];
		sumX += x[i] * weights[i];
		sumY += y[i] * weights[i];
		minX = MIN(minX, x[i]);
		maxX = MAX(maxX, x[i]);
		minY = MIN(minY, y[i]);
		maxY = MAX(maxY, y[i]);
	}
	float dx = maxX - minX;
	float dy = maxY - minY;
	if (doRotateForR2 && (dy > 2 * dx))
	{
		float* xr = new float[n];
		float* yr = new float[n];
		memcpy(xr, x, sizeof(float) * n);
		memcpy(yr, y, sizeof(float) * n);
		dx = x[n-1] - x[0];
		dy = y[n-1] - y[0];
		float aa = PI_FLOAT / 4 - atan2(dy, dx);
		float angle = aa;
		rotatePoints(xr, yr, n, angle);
		LeastSquaresFitResult rotatedResult;
		leastSquaresFitWeighted(xr, yr, weights, n, false, rotatedResult, dys);
		result.stddev = rotatedResult.stddev;
		result.r2 = rotatedResult.r2;
		delete [] xr;
		delete [] yr;
		return;
	}
	if (sumWeight <= 0)
	{
		return;
	}
	double meanX = sumX / sumWeight;
	double ss = 0;
	double s2 = 0;
	double meanY = sumY / sumWeight;
	int nNonZeroWeight = 0;
	for (int i = 0; i < n; i++)
	{
		double xv = x[i] - meanX;
		ss += xv * xv * weights[i];
		double v2 = xv * y[i] * weights[i];
		s2 += v2;
		assert(weights[i] >= 0);
		if (weights[i] > 0.000001f)
		{
			nNonZeroWeight++;
		}
	}
	if (nNonZeroWeight == 1)
	{
		result.b = 0;
		result.a = sumY / sumWeight;
		result.stddev = 0;
		result.r2 = 1; 
		return;
	}
	result.b = s2 / ss;
	result.a = (sumY - sumX * result.b) / sumWeight;
	if (doRotateForR2 && (abs(result.b) < 0.5f)) 
	{
		LeastSquaresFitResult rotatedResult;
		float* xr = new float[n];
		float* yr = new float[n];
		memcpy(xr, x, sizeof(float) * n);
		memcpy(yr, y, sizeof(float) * n);
		float aa = PI_FLOAT / 6 - atan(abs(result.b));
		float angle = (result.b < 0 ? -aa : aa);
		rotatePoints(xr, yr, n, angle); 
		leastSquaresFitWeighted(xr, yr, weights, n, false, rotatedResult, dys);
		result.stddev = rotatedResult.stddev;
		result.r2 = rotatedResult.r2;
		delete [] xr;
		delete [] yr;
		return;
	}
	double chi2 = 0;
	double yhv = 0, yov = 0;
	for (int i = 0; i < n; i++)
	{
		double ycomp = result.a + result.b * x[i];
		double t = (ycomp - meanY) * weights[i];
		yhv += t * t;
		double tyov = (y[i] - meanY) * weights[i];
		yov += tyov * tyov;
		double v = (y[i] - ycomp);
		chi2 += v * v * weights[i];
		if (dys != NULL)
		{
			dys[i] = v;
		}
	}
	chi2 /= sumWeight;
	result.stddev = sqrt(chi2);
	if (yov > 0)
	{
		result.r2 = yhv / yov;
	}
	else
	{
		result.r2 = 1;
	}
}
float weightedMean(float* v, float* weights, int n)
{
	float sum = 0;
	float sumWeight = 0;
	for (int i = 0; i < n; i++)
	{
		sum += v[i] * weights[i];
		sumWeight += weights[i];
	}
	return sum / sumWeight;
}
float geometricMean(float x, float y)
{
	float prod = x * y;
	if (prod <= 0) return 0;
	return sqrtf(prod);
}
float geometricMean(float x, float y, float z)
{
	float prod = x * y * z;
	if (prod <= 0) return 0;
	return pow(prod, 0.333333333f);
}
float geometricMean(vector<float> v)
{
	int n = v.size();
	float prod = 1;
	for (int i = 0; i < n; i++) prod *= v[i];
	if (prod <= 0) return 0;
	return pow(prod, 1.0f/n);
}
float geometricMean(vector<float> v, vector<float> weights)
{
	int n = v.size();
	float prod = 1;
	float sumWeight = 0;
	for (int i = 0; i < n; i++)
	{
		prod *= pow(v[i], weights[i]);
		sumWeight += weights[i];
	}
	if ((prod <= 0) || (sumWeight <= 0)) return 0;
	return pow(prod, 1.0f / sumWeight);
}
float pyth(int x, int y)
{
	float prod = (float)(x * x + y * y);
	return sqrtf(prod);
}
float pyth(float x, float y)
{
	float prod = x * x + y * y;
	return sqrtf(prod);
}
double pyth(double x, double y)
{
	double prod = x * x + y * y;
	return sqrt(prod);
}
float evalGaussianPdf(float mean, float stddev, float val)
{
	float s2 = stddev * stddev;
	float exp1 = -1 / (2 * s2);
	float denom = 1 / (stddev * SQRT_TWO_PI_FLOAT);
	float b = val - mean;
	float y = denom * pow(E_FLOAT, exp1 * b * b);
	return y;
}
double evalGaussianPdf(double mean, double stddev, double val)
{
	double s2 = stddev * stddev;
	double exp1 = -1 / (2 * s2);
	double denom = 1 / (stddev * SQRT_TWO_PI_FLOAT);
	double b = val - mean;
	double y = denom * pow(E_FLOAT, exp1 * b * b);
	return y;
}
float sumv(vector<float>& v)
{
	float sum = 0;
	int n = v.size();
	for (int i = 0; i < n; i++) sum += v[i];
	return sum;
}
float meanv(vector<float>& v)
{
	return sumv(v) / v.size();
}
float maxv(vector<float>& v)
{
	float max = -FLT_MAX;
	int n = v.size();
	for (int i = 0; i < n; i++) max = MAX(max, v[i]);
	return max;
}
float minv(vector<float>& v)
{
	float m = FLT_MAX;
	int n = v.size();
	for (int i = 0; i < n; i++) m = MIN(m, v[i]);
	return m;
}
float medianv(vector<float>& v)
{
	int n = v.size();
	if (n == 1) return v[0];
	if (n == 2) return (v[0] + v[1]) / 2;
	float m = FLT_MAX;
	vector<float> dup = v;
	if ((n % 2) == 0)
	{
		std::nth_element(dup.begin(), dup.begin() + n/2 - 1, dup.end());
		float a = dup[n/2 - 1];
		std::nth_element(dup.begin(), dup.begin() + n/2, dup.end());
		float b = dup[n/2];
		m = (a + b)/2;
	}
	else
	{
		std::nth_element(dup.begin(), dup.begin() + n/2, dup.end());
		m = dup[n/2];
	}
	return m;
}
void setv(vector<float>& v, float val)
{
	int n = v.size();
	for (int i = 0; i < n; i++) v[i] = val;
}
vector<float> absv(vector<float>& v)
{
	int n = v.size();
	vector<float> a(n);
	for (int i = 0; i < n; i++) a[i] = abs(v[i]);
	return a;
}
vector<float> exceptv(vector<float>& v, float val)
{
	int n = v.size();
	vector<float> a;
	for (int i = 0; i < n; i++) if (v[i] != val) a.push_back(v[i]);
	return a;
}
bool anygoodv(vector<float>& v)
{
	int n = v.size();
	for (int i = 0; i < n; i++) if (checkIsNormalFloat(v[i])) return true;
	return false;
}
float pctDiff(float a, float b)
{
	float dv = abs(a - b);
	float m = (a + b) / 2;
	if (m == 0)
	{
		if (dv == 0) return 0;
		else return dv; 
	}
	return abs(dv / m);
}
class FArrayStats
{
public:
	float mean;
	float stddev; 
	float meandev;
	float max;
	float min;
	int count;
	float* data;
	bool isDataSorted;
	void Init()
	{
		count = 0;
		max = -3.402823466e+38F;
		min = +3.402823466e+38F;
		mean = 0.0f;
		stddev = 0.0f;
		meandev = 0.0f;
		this->data = NULL;
		isDataSorted = false;
	}
	FArrayStats()
	{
		Init();
	}
	FArrayStats(const FArrayStats& other)
	{
		count = other.count;
		max = other.max;
		min = other.min;
		mean = other.mean;
		stddev = other.stddev;
		meandev = other.meandev;
		this->data = NULL;
		isDataSorted = false;
		if (other.data != NULL)
		{
			this->data = new float[count];
			memcpy(this->data, other.data, count * sizeof(float));
			isDataSorted = other.isDataSorted;
		}
	}
	FArrayStats& operator=(const FArrayStats& other)
	{
		count = other.count;
		max = other.max;
		min = other.min;
		mean = other.mean;
		stddev = other.stddev;
		meandev = other.meandev;
		isDataSorted = false;
		if (this->data != NULL)
		{
			delete [] this->data;
			this->data = NULL;
		}
		if (other.data != NULL)
		{
			this->data = new float[count];
			memcpy(this->data, other.data, count * sizeof(float));
			isDataSorted = other.isDataSorted;
		}
		return *this;
	}
	FArrayStats(std::vector<float>& v, bool includePercentiles)
	{
		Init();
		compute(v, includePercentiles);
	}
	~FArrayStats()
	{
		if (this->data != NULL)
		{
			delete [] this->data;
			this->data = NULL;
		}
	}
	bool compute(int* array, int count)
	{
		return compute(array, count, false);
	}
	bool compute(vector<int>& v, bool includePercentiles)
	{
		return compute(v.data(), v.size(), includePercentiles);
	}
	bool compute(int* array, int count, bool includePercentiles)
	{
		float* fa = new float[count];
		for (int i = 0; i < count; i++) { fa[i] = array[i]; }
		bool ret = compute(fa, count, includePercentiles);
		delete [] fa;
		return ret;
	}
	bool compute(float* a, int count)
	{
		return compute(a, count, false);
	}
	bool compute(vector<float>& v, bool includePercentiles)
	{
		return compute(v.data(), (int)v.size(), includePercentiles);
	}
	bool compute(float* array, int count, bool includePercentiles)
	{
		float sum = 0.0f;
		float sqrSum = 0.0f;
		float val;
		int nGood = 0;
		for (int i = 0; i < count; i++)
		{
			val = array[i];
			if (checkIsNormalFloat(val))
			{
				if (val > max) max = val;
				if (val < min) min = val;
				sum += val;
				sqrSum += val * val;
				nGood++;
			}
		}
		this->count = nGood;
		mean = sum / nGood;
		float variance = (sqrSum / nGood) - mean * mean;
		if (variance <= 0)
		{
			stddev = 0;
		}
		else
		{
			stddev = sqrt(variance);
		}
		float sumDev = 0.0f;
		for (int i = 0; i < count; i++)
		{
			if (checkIsNormalFloat(array[i]))
			{
				sumDev += abs(array[i] - mean);
			}
		}
		meandev = sumDev / nGood;
		if (includePercentiles)
		{
			if (this->data != NULL)
			{
				delete [] this->data;
			}
			this->data = new float[nGood];
			int j = 0;
			for (int i = 0; i < count; i++)
			{
				if (checkIsNormalFloat(array[i]))
				{
					this->data[j++] = array[i];
				}
			}
		}
		return true;
	}
	std::string toString()
	{
		char tmp[1000];
		sprintf(tmp, "n: %d, min: %.3f, max: %.3f, mean: %.3f, stddev: %.3f", count, min, max, mean, stddev);
		return (std::string)tmp;
	}
	void sortData()
	{
		assert(this->data != NULL);
		std::sort(this->data, this->data + count);
		this->isDataSorted = true;
	}
	float getPercentile(float percentile)
	{
		assert(this->data != NULL);
		assert(percentile >= 0.0f);
		assert(percentile <= 1.0f);
		if (count == 0)
		{
			ErrorExit("FArrayStats: No values, can't get percentile.");
		}
		if (this->data != NULL)
		{
			int n = (int)(percentile * (count - 1));
			return getNthValue(n);
		}
		else
		{
			ErrorExit("FArrayStats doesn't have any data to getPercentile of.");
			return getNthValue(0);
		}
	}
	float getNthValue(int n)
	{
		if (n >= count)
		{
			ErrorExit("n is gte count");
			n = count - 1;
		}
		if (!this->isDataSorted)
		{
			nth_element(data, data + n, data + count);
		}
		return data[n];
	}
};
#include <vector>
#include <string>
#include <cstdio>
#include <sstream>
#include <iostream>
#include <functional> 
#include <cctype>
using std::string;
using std::vector;
class StringUtils
{
public:
	static std::vector<std::string>& split(const std::string &s, char delim, std::vector<std::string> &elems, bool trim = true) 
	{
		std::stringstream ss(s);
		std::string item;
		while (std::getline(ss, item, delim)) 
		{
			elems.push_back(TrimLeft(item));
		}
		return elems;
	}
	static void splitPreAlloc(const std::string &s, char delim, std::vector<std::string> &elems, bool trim = true)
	{
		std::stringstream ss(s);
		std::string item;
		int i = 0;
		while (std::getline(ss, item, delim)) 
		{
			elems[i++] = TrimLeft(item);
		}
	}
	static std::vector<std::string> split(const std::string &s, char delim, bool trim = true)
	{
		std::vector<std::string> elems;
		return split(s, delim, elems, trim);
	}
	static bool splitString(string& inStr, char splitChar, vector<string>& resultList, bool trim = true)
	{
		size_t nextPos, pos = 0;
		size_t len;
		if (inStr.size() <= 0)
		{
			return false;;
		}
		while (pos < inStr.size())
		{
			nextPos = inStr.find(splitChar, pos);
			if (nextPos == (size_t)(-1))
			{
				nextPos = inStr.size();
			}
			if (nextPos != pos)
			{
				len = nextPos - pos;
				string subString(inStr, pos, len);
				resultList.push_back(TrimLeft(subString));
			}
			pos = nextPos + 1;
		}
		return true;
	};
	static std::string& TrimLeft(std::string &s)
	{
		s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
		return s;
	}
	static string& TrimRight(string &s)
	{
		s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
		return s;
	}
};
#include <cmath>
#include <time.h>
string toStringBits(byte b)
{
	string s;
	for (int i = 0; i < 8; i++) s += (((b >> (7-i)) & 1) ? "1" : "0");
	return s;
}
string toStringBits(unsigned int b)
{
	string s;
	for (int i = 0; i < 32; i++) s += (((b >> (31-i)) & 1) ? "1" : "0");
	return s;
}
string toStringBits(unsigned short b)
{
	string s;
	for (int i = 0; i < 16; i++) s += (((b >> (15-i)) & 1) ? "1" : "0");
	return s;
}
bool checkIsNormalFeatureValue(float v)
{
	if (!checkIsNormalFloat(v)) return false; 
	if (v < -998) return false;
	return true;
}
float getTimeSeconds()
{
	return (0.0f + clock()) / CLOCKS_PER_SEC;
}
int singleCluster(vector<float>& values, int nMin, byte initialSelectionMask, byte& selectedMask, float& mean, float& stddev, vector<float>& weights, bool debug = false)
{
	if (values.size() > 8)
	{
		ErrorExit("This can't cluster more than 8 numbers.");
	}
	int n = values.size();
	float sum = 0, s2 = 0;
	weights.resize(n);
	for (int i = 0; i < n; i++) weights[i] = 1;
	selectedMask = (initialSelectionMask & ((1 << n) - 1));
	if (debug) logd.printf("singleCluster: Initial mask: %s\n", toStringBits(selectedMask).c_str());
	int ns = 0; 
	for (int i = 0; i < n; i++) 
	{
		if ((selectedMask >> i) & 1)
		{
			if (values[i] < -998)
			{
				assert(false && "Probable magic value coming into singleCluster");
				values[i] = 0;
			}
			sum += values[i];
			s2 += values[i] * values[i];
			ns++;
		}
	}
	if (ns == 0)
	{
		ErrorExit("Input selection mask is 0 so no sense calling this.");
	}
	mean = sum / ns;
	stddev = sqrt(((ns * s2) - sum * sum) / (ns * (ns - 1)));
	if (debug)
	{
		logd.printf("Start\n");
		logd.printf("  mean: %f\n", mean);
		logd.printf("  stddev: %f\n", stddev);
	}
	if (stddev < 0.0001f) return ns;
	int nIter = 4; 
	for (int i = 0; i < n; i++)
	{
		if ((selectedMask >> i) & 1)
		{
			weights[i] = 1;
		}
		else
		{
			weights[i] = 0;
		}
	}
	for (int iter = 0; iter < nIter; iter++)
	{
		int minIndex = -1;
		float minWeight = 1;
		for (int i = 0; i < n; i++) 
		{
			if ((selectedMask >> i) & 1)
			{
				float nstd = abs(values[i] - mean) / stddev;
				weights[i] = 1.0f - rangeScore(0.5f, 2.75f, nstd); 
				if (weights[i] < minWeight) { minWeight = weights[i]; minIndex = i; }
			}
		}
		ns = 0;
		for (int i = 0; i < n; i++) ns += ((selectedMask >> i) & 1);
		if ((ns > nMin) && (minWeight < 0.66f))
		{
			if (debug) logd.printf("drop from %d by dropping index %d\n", ns, minIndex);
			selectedMask &= ~(1 << minIndex);
			weights[minIndex] = 0;
			ns--;
		}
		sum = s2 = 0;
		float sumWeight = 0;
		for (int i = 0; i < n; i++) 
		{
			if ((selectedMask >> i) & 1)
			{
				sum += weights[i] * values[i];
				s2 += weights[i] * values[i] * values[i];
				sumWeight += weights[i];
			}
		}
		mean = sum / sumWeight;
		assert(sumWeight >= 1 && "sumWeight less than 1");
		float t1 = (sumWeight * s2) - sum * sum;
		if (t1 > 0)
		{
			stddev = sqrt(t1 / (sumWeight * (sumWeight - 1)));
		}
		else
		{
			stddev = 0;
		}
		if (stddev < 0.0001f) return ns;
		if (mean != mean)
		{
			assert(false && "Bad mean");
		}
		if (stddev != stddev)
		{
			assert(false && "Bad stddev");
		}
		if (debug)
		{
			logd.printf("Iter %d\n", iter);
			logd.printf("  ns: %d\n", ns);
			logd.printf("  weighted mean: %f\n", mean);
			logd.printf("  weighted stddev: %f\n", stddev);
			logd.printf("  mask: "); for (int i = 0; i < n; i++) logd.printf("%d", ((selectedMask >> i) & 1)); logd.printf("\n");
			for (int i = 0; i < n; i++)
			{
				logd.printf("    [%d]: %f %f\n", i, values[i], weights[i]);
			}
		}
	}
	for (int i = 0; i < n; i++) 
	{
		if (weights[i] < 0.5f) selectedMask &= ~(1 << i);
	}
	ns = 0;
	for (int i = 0; i < n; i++) ns += ((selectedMask >> i) & 1);
	return ns;
}
int singleCluster(int nMin, byte initialSelectionMask, byte& selectedMask, float& mean, float& stddev, vector<float>& weights, float v0, float v1, float v2, float v3, bool debug = false)
{
	vector<float> values(4);
	values[0] = v0;
	values[1] = v1;
	values[2] = v2;
	values[3] = v3;
	return singleCluster(values, nMin, initialSelectionMask, selectedMask, mean, stddev, weights, debug);
}
int singleCluster(int nMin, byte initialSelectionMask, byte& selectedMask, float& mean, float& stddev, vector<float>& weights, float v0, float v1, float v2, bool debug = false)
{
	vector<float> values(3);
	values[0] = v0;
	values[1] = v1;
	values[2] = v2;
	return singleCluster(values, nMin, initialSelectionMask, selectedMask, mean, stddev, weights, debug);
}
void normMax(vector<float>& values)
{
	float max = -FLT_MAX;
	int n = values.size();
	for (int i = 0; i < n; i++) max = MAX(max, values[i]);
	for (int i = 0; i < n; i++) values[i] = values[i] / max;
}
void testSingleCluster()
{
	byte initialMask = 0xFF;
	byte mask;
	float mean, stddev;
    int ns;
	vector<float> weights;
    ns = singleCluster(2, initialMask, mask, mean, stddev, weights, 0.324032992f, 0.324647695f, 0.324000239f, true);
	logd.printf("n selected: %d, mean: %f, stddev: %f\n", ns, mean, stddev);
}
void RandomizeIndices(int n, int nTrain, int nTest, vector<int>& trainIndices, vector<int>& testIndices)
{
	for (int i = 0; i < nTrain; i++)
	{
		int idx = rand() % n;
		trainIndices.push_back(idx);
	}
	for (int i = 0; i < nTest; i++)
	{
		int idx = rand() % n;
		testIndices.push_back(idx);
	}
}
void FoldIndices(int n, int fold, int nFolds, vector<int>& trainIndices, vector<int>& testIndices)
{
	if (nFolds == 1)
	{
		for (int i = 0; i < n; i++)
		{
			trainIndices.push_back(i);
			testIndices.push_back(i);
		}
	}
	else
	{
		int npf = n / nFolds;
		for (int i = 0; i < n; i++)
		{
			bool isTest = (i >= npf * fold) && (i < (fold + 1) * npf);
			if (!isTest) trainIndices.push_back(i);
			else testIndices.push_back(i);
		}
	}
}
template<typename T>
inline std::vector<T> erase_indices(const std::vector<T>& data, std::vector<size_t>& indicesToDelete)
{
	if (indicesToDelete.empty())
		return data;
	std::vector<T> ret;
	ret.reserve(data.size() - indicesToDelete.size());
	std::sort(indicesToDelete.begin(), indicesToDelete.end());
	auto itBlockBegin = data.begin();
	for (std::vector<size_t>::const_iterator it = indicesToDelete.begin(); it != indicesToDelete.end(); ++it)
	{
		auto itBlockEnd = data.begin() + *it;
		if (itBlockBegin != itBlockEnd)
		{
			std::copy(itBlockBegin, itBlockEnd, std::back_inserter(ret));
		}
		itBlockBegin = itBlockEnd + 1;
	}
	if (itBlockBegin != data.end())
	{
		std::copy(itBlockBegin, data.end(), std::back_inserter(ret));
	}
	return ret;
}
float NanIfBadOrLteZero(float x)
{
	if (checkIsNormalFeatureValue(x) && (x > 0)) return x;
	else return std::numeric_limits<float>::quiet_NaN();
}
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <limits>
#include <vector>
#include <queue>
#include <memory>
#include <stdexcept>
#include <vector>
#include <queue>
#include <memory>
#include <cassert>
class IRoi
{
public:
	int x0;
	int y0;
	int x1; 
	int y1; 
	int width() { return(x1 - x0); } 
	int height() { return(y1 - y0); } 
	int getWidth() { return(x1 - x0); }
	int getHeight() { return(y1 - y0); }
	int centerX() { return(x0 + width()/2); }
	int centerY() { return(y0 + height()/2); }
	IRoi()
	{	
		x0 = y0 = x1 = y1 = 0;	
	}
	IRoi(int x0, int y0, int x1, int y1)
	{
		this->set(x0, y0, x1, y1);
	}
	void set(int x0, int y0, int x1, int y1)
	{
		assert(x1 >= x0);
		assert(y1 >= y0);
		this->x0 = x0;
		this->y0 = y0;
		this->x1 = x1;
		this->y1 = y1;
	}
	IRoi(const IRoi& other)
	{
		this->set(other.x0, other.y0, other.x1, other.y1);
	}
	IRoi& operator=(const IRoi& other)
	{
		this->set(other.x0, other.y0, other.x1, other.y1);
		return(*this);
	}
	string ToString()
	{
		char tmp[256];
		sprintf(tmp, "(%d,%d to %d,%d)", x0, y0, x1, y1);
		return (string)tmp;
	}
	void intersect(int ox0, int oy0, int ox1, int oy1)
	{
		x0 = MAX(x0, ox0);
		y0 = MAX(y0, oy0);
		x1 = MIN(x1, ox1);
		y1 = MIN(y1, oy1);
	}
	void intersect(IRoi& other)
	{
		x0 = MAX(x0, other.x0);
		y0 = MAX(y0, other.y0);
		x1 = MIN(x1, other.x1);
		y1 = MIN(y1, other.y1);
	}
	void intersect(int ox0, int oy0, int ox1, int oy1, IRoi& parallelRoi)
	{
		int newX0 = MAX(x0, ox0);
		parallelRoi.x0 += (newX0 - x0);
		x0 = newX0;
		int newY0 = MAX(y0, oy0);
		parallelRoi.y0 += (newY0 - y0);
		y0 = newY0;
		int newX1 = MIN(x1, ox1);
		parallelRoi.x1 += (newX1 - x1);
		x1 = newX1;
		int newY1 = MIN(y1, oy1);
		parallelRoi.y1 += (newY1 - y1);
		y1 = newY1;
	}
	bool checkIsInImage(int width, int height)
	{
		return ((this->x0 >= 0) && (this->x1 <= width) && (this->x0 < this->x1)
			&& (this->y0 >= 0) && (this->y1 <= height) && (this->y0 < this->y1));
	}
	void offset(int dx, int dy)
	{
		this->x0 += dx;
		this->x1 += dx;
		this->y0 += dy;
		this->y1 += dy;
	}
	void offsetCenter(int dx, int dy, int minX, int minY, int maxX, int maxY)
	{
		offset(dx, dy);
		shrinkToFit(minX, minY, maxX, maxY);
	}
	void shrinkToFit(int minX, int minY, int maxX, int maxY)
	{
		int xShrink = MAX(minX - this->x0, this->x1 - maxX);
		int yShrink = MAX(minY - this->y0, this->y1 - maxY);
		xShrink = MAX(0, xShrink);
		yShrink = MAX(0, yShrink);
		this->x0 += xShrink;
		this->x1 -= xShrink;
		this->y0 += yShrink;
		this->y1 -= yShrink;
	}
	bool checkOverlap(IRoi& other)
	{
		return !((this->y0 >= other.y1)
			|| (this->y1 <= other.y0)
			|| (this->x0 >= other.x1)
			|| (this->x1 <= other.x0));
	}
	int distTaxi(IRoi other)
	{
		if (checkOverlap(other)) return 0;
		int m = 0;
		if (x0 >= other.x1) m = MAX(m, abs(x0 - other.x1));
		else if (x1 <= other.x0) m = MAX(m, abs(x1 - other.x0));
		if (y0 >= other.y1) m = MAX(m, abs(y0 - other.y1));
		else if (y1 <= other.y0) m = MAX(m, abs(y1 - other.y0));
		return m;
	}
	float dist(IRoi other)
	{
		if (checkOverlap(other)) return 0;
		if (!(this->y0 >= other.y1 || this->y1 <= other.y0))
		{
			return distTaxi(other);
		}
		if (!(this->x0 >= other.x1 || this->x1 <= other.x0))
		{
			return distTaxi(other);
		}
		if (x0 >= other.x1)
		{
			if (y0 >= other.y1)
			{
				return pyth((float)x0 - other.x1, (float)y0 - other.y1);
			}
			else
			{
				return pyth((float)x0 - other.x1, (float)y1 - other.y0);
			}
		}
		else if (x1 <= other.x0)
		{
			if (y0 >= other.y1)
			{
				return pyth((float)x1 - other.x0, (float)y0 - other.y1);
			}
			else
			{
				return pyth((float)x1 - other.x0, (float)y1 - other.y0);
			}
		}
		else
		{
			ErrorExit("Oops, logic bug in dist.");
			return FLT_MAX;
		}
	}
};
class Img
{
protected:
	int width;
	int height;
	int depthInBits;
	int pitch; 
	int channels; 
	int bufferSize; 
	unsigned char* data;
	bool isDataOwned; 
public:
	int getWidth() { return(width); }
	int getHeight() { return(height); }
	int getDepthInBits() { return(depthInBits); }
	int getPitch() { return pitch; } 
	int getBufferSize() { return bufferSize; } 
	unsigned char* getData() { return(data); }
	Img()
	{
		init();
	}
	~Img()
	{
		deallocate();
	}
	Img(int width, int height, int depthInBits=8)
	{
		init();
		assertDepthInBits(depthInBits);
		this->width = width;
		this->height = height;
		this->depthInBits = depthInBits;
		this->pitch = computePitch(width);
		assertPitch(this->pitch);
		allocate();
	}
	Img(int width, int height, int depthInBits, int pitch)
	{
		init();
		assertDepthInBits(depthInBits);
		this->width = width;
		this->height = height;
		this->depthInBits = depthInBits;
		assertPitch(pitch);
		this->pitch = pitch;
		if (!allocate()) return;
	}
	Img(int width, int height, int depthInBits, int pitch, const unsigned char* data)
	{
		init();
		assertDepthInBits(depthInBits);
		this->width = width;
		this->height = height;
		this->depthInBits = depthInBits;
		this->pitch = pitch;
		assertPitch(pitch);
		if (!allocate()) return;
		memcpy(this->data, data, bufferSize);
	}
	Img(const Img& image)
	{
		init();
		copyImg(image);
	}
	void init()
	{
		data = NULL;
		bufferSize = 0;
		width = 0;
		height = 0;
		pitch = 0;
		depthInBits = 0;
		channels = 1;
		isDataOwned = true; 
	}
	bool allocate()
	{
		assert((depthInBits % 8) == 0);
		deallocate();
		bufferSize = height * pitch;
		if (bufferSize == 0)
		{
			return true;
		}
		data = (unsigned char*)malloc(bufferSize);
		isDataOwned = true;
		if (data == NULL) 
		{
			logd.error("Img: mem alloc failed on size %d", bufferSize);
			ErrorExit("Img: mem alloc failed");
			return false;
		}
		else 
		{
			return true;
		}
	}
	private:
	bool deallocate()
	{
		if ((data != NULL) && isDataOwned)
		{
			free(data);
		}
		data = NULL;
		bufferSize = 0;
		return true;
	}
	public:
	void destroy()
	{
		deallocate();
		init();
	}
	bool setByCopy(int width, int height, int depthInBits, int pitch, const unsigned char* buffer)
	{
		int newPitch = staticComputePitch(width, depthInBits);
		resize(width, height, depthInBits, newPitch);
		if (this->pitch == pitch)
		{
			memcpy(this->data, buffer, height * pitch);
		}
		else
		{
			int n = this->width * this->depthInBits >> 3;
			for (int y = 0; y < height; y++)
			{
				memcpy(this->data + y * this->pitch, buffer + y * pitch, n);
			}
		}
		return true;
	}
	bool setByPtr(int width, int height, int depthInBits, int pitch, unsigned char* buffer)
	{
		deallocate();
		init();
		this->width = width;
		this->height = height;
		this->pitch = pitch;
		this->depthInBits = depthInBits;
		this->data = buffer;
		this->bufferSize = height * pitch;
		isDataOwned = false;
		return true;
	}
	unsigned char* getData(int x, int y)
	{
		int db = this->depthInBits >> 3;
		return data + y * pitch + x * db;
	}
	uint16_t* getData16u(int x, int y)
	{
		assert(this->depthInBits == 16);
		int pp = pitch >> 1;
		return (uint16_t*)data + y * pp + x;
	}
	bool SetRoi(IRoi& thisRoi, Img& other, IRoi& otherRoi)
	{
		assert(this->depthInBits == other.depthInBits);
		bool roisSameSize = ((thisRoi.width() == otherRoi.width()) && (thisRoi.height() == otherRoi.height()));
		bool thisRoiFits = thisRoi.checkIsInImage(this->width, this->height);
		bool otherRoiFits = otherRoi.checkIsInImage(other.getWidth(), other.getHeight());
		if (!(roisSameSize && thisRoiFits && otherRoiFits))
		{
			assert(roisSameSize && thisRoiFits && otherRoiFits);
			return false;
		}
		int thisIndex, otherIndex;
		int dx = otherRoi.x0 - thisRoi.x0;
		int dy = otherRoi.y0 - thisRoi.y0;
		int otherPitch = other.getPitch();
		unsigned char* ps = this->data;
		unsigned char* pd = other.data;
		int depthBytes = (this->depthInBits >> 3);
		int len = (thisRoi.x1 - thisRoi.x0) * depthBytes;
		for (int y = thisRoi.y0; y < thisRoi.y1; y++)
		{
			thisIndex = thisRoi.x0 * depthBytes + y * pitch;
			otherIndex = (thisRoi.x0 + dx) * depthBytes + (y + dy) * otherPitch;
			memcpy(&ps[thisIndex], &pd[otherIndex], len);
		}
		return true;
	}
	bool SetRoiAllowOffSrc(IRoi& dstRoi, Img& src, IRoi& srcRoi)
	{
		assert(this->depthInBits == src.depthInBits);
		bool roisSameSize = ((dstRoi.width() == srcRoi.width()) && (dstRoi.height() == srcRoi.height()));
		bool dstRoiFits = dstRoi.checkIsInImage(this->width, this->height);
		if (!(roisSameSize && dstRoiFits))
		{
			assert(roisSameSize && dstRoiFits);
			return false;
		}
		int srcPitch = src.getPitch();
		unsigned char* pd = this->data;
		unsigned char* ps = src.data;
		int depthBytes = (this->depthInBits >> 3);
		IRoi clippedSrcRoi(srcRoi);
		IRoi clippedDstRoi(dstRoi);
		clippedSrcRoi.intersect(0, 0, src.getWidth(), src.getHeight(), clippedDstRoi); 
		int nBytes = (clippedSrcRoi.x1 - clippedSrcRoi.x0) * depthBytes;
		int dstIndex, srcIndex;
		int nRows = clippedSrcRoi.y1 - clippedSrcRoi.y0;
		for (int i = 0; i < nRows; i++)
		{
			dstIndex = clippedDstRoi.x0 * depthBytes + (i + clippedDstRoi.y0) * pitch; 
			srcIndex = clippedSrcRoi.x0 * depthBytes + (i + clippedSrcRoi.y0) * srcPitch; 
			memcpy(&pd[dstIndex], &ps[srcIndex], nBytes);
		}
		return true;
	}
	int computePitch(int width)
	{
		assertDepthInBits(depthInBits);
		assert(width >= 0);
		int align = 32;
		int s = (width - 1) / align;
		int w = (s + 1) * align;
		return w * (depthInBits >> 3);
	}
	static int staticComputePitch(int width, int depthInBits)
	{
		assert(width >= 0);
		assert(depthInBits >= 8);
		int align = 32;
		int db = depthInBits >> 3;
		int s = (width - 1) / align;
		int w = (s + 1) * align;
		return w * db;
	}
	void assertPitch(int pitch)
	{
		assert(pitch >= width * (depthInBits >> 3));
		assert(pitch % 16 == 0);
	}
	void assertDepthInBits(int depthInBits)
	{
		assert(depthInBits > 1);
		assert((depthInBits % 8) == 0);
	}
	bool resize(int width, int height)
	{
		if (depthInBits == 0)
		{
			depthInBits = 8;
		}
		assertDepthInBits(depthInBits);
		int newPitch = computePitch(width);
		return resize(width, height, depthInBits, newPitch);
	}
	bool resize(int width, int height, int depthInBits)
	{
		if ((this->width == width)
			&& (this->height == height)
			&& (this->depthInBits == depthInBits)
			&& (this->channels == 1))
		{
			if (data != NULL) return true; 
			else return allocate();
		}
		this->depthInBits = depthInBits;
		assertDepthInBits(depthInBits);
		int newPitch = computePitch(width);
		assert(newPitch >= width * depthInBits >> 3);
		return resize(width, height, depthInBits, newPitch);
	}
	bool checkSameSize(Img& other)
	{
		if ((this->width == other.width)
			&& (this->height == other.height)
			&& (this->depthInBits == other.depthInBits)
			&& (this->pitch == other.pitch)
			&& (this->channels == other.channels))
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	bool resize(int width, int height, int depthInBits, int pitch)
	{
		return resize(width, height, depthInBits, 1, pitch);
	}
	bool resize(int width, int height, int depthInBits, int channels, int pitch)
	{
		int minPitch = width * depthInBits / 8;
		assert(pitch >= minPitch);
		if (pitch < minPitch)
		{
			ErrorExit("resize: Specified pitch is less than the minimum possible.");
		}
		if ((this->width == width)
			&& (this->height == height)
			&& (this->depthInBits == depthInBits)
			&& (this->pitch == pitch)
			&& (this->channels == channels))
		{
			if (data != NULL) return true; 
			else return allocate();
		}
		if (data != NULL) deallocate();
		assert(width >= 0);
		assert(height >= 0);
		this->width = width;
		this->height = height;
		this->depthInBits = depthInBits;
		this->pitch = pitch;
		this->channels = channels;
		return allocate();
	}
	bool copyImg(const Img& image)
	{
		if (!resize(image.width, image.height, image.depthInBits, image.pitch))
		{
			logd.error("copyImg: Unable to resize.");
			return false;
		}
		memcpy(data, image.data, bufferSize);
		return true;
	}
	int getPixel(int i)
	{
		return getPixel(i % width, i / width);
	}
	bool checkIsInImage(int i)
	{
		return ((i >= 0) && (i < width * height));
	}
	int getPixel(int x, int y) 
	{ 
		assert(x >= 0);
		assert(y >= 0);
		assert(x < width);
		assert(y < height);
		int db = this->depthInBits >> 3;
		if (db == 1)
		{
			return((int)data[x + y * pitch]); 
		}
		else if (db == 2)
		{
			short v = (*(short*)&data[x * db + y * pitch]); 
			return v;
		}
		else if (db == 3)
		{
			int v = (*(int*)&data[x * db + y * pitch]); 
			v = v & 0x00FFFFFF;
			return v;
		}
		else if (db == 4)
		{
			return(*(int*)&data[x * db + y * pitch]); 
		}
		else
		{
			ErrorExit("Unimplemented depth in bits.");
			return 0;
		}
	}
	int getPixel16u(int i) 
	{ 
		return getPixel16u(i % width, i / width);
	}
	int getPixel16u(int x, int y) 
	{ 
		assert(this->depthInBits == 16);
		unsigned short* pd = (unsigned short*)data;
		return pd[x + y * (pitch >> 1)]; 
	}
	int getPixel32s(int x, int y) 
	{ 
		assert(this->depthInBits == 32);
		assert(x >= 0);
		assert(y >= 0);
		assert(x < width);
		assert(y < height);
		int* pd = (int*)data;
		return pd[x + y * (pitch >> 2)]; 
	}
	float getPixel32f(int x, int y) 
	{ 
		assert(this->depthInBits == 32);
		assert(x >= 0);
		assert(y >= 0);
		assert(x < width);
		assert(y < height);
		float* pd = (float*)data;
		return pd[x + y * (pitch >> 2)]; 
	}
	void setPixel(int i, int val) 
	{ 
		assert(this->depthInBits == 8);
		assert((i >= 0) && (i < width * pitch));
		data[i] = val; 
	}
	void setPixel(int x, int y, int val) 
	{ 
		assert(this->depthInBits == 8);
		assert((x >= 0) && (y >= 0) && (x < width) && (y < height));
		data[x + y * pitch] = val; 
	}
	void setRow(int y, int x0, int x1, int val)
	{
		assert(this->depthInBits == 8);
		assert((x0 >= 0) && (y >= 0) && (x0 < width) && (y < height));
		assert((x1 >= 0) && (x1 < width));
		memset(&data[x0 + y * pitch], val, (x1 - x0 + 1));
	}
	void setPixel32s(int x, int y, int val)
	{
		assert(this->depthInBits == 32);
		int* pd = (int*)data;
		pd[x + y * (pitch >> 2)] = val; 
	}
	void setPixel32f(int x, int y, float val)
	{
		assert(this->depthInBits == 32);
		float* pd = (float*)data;
		pd[x + y * (pitch >> 2)] = val; 
	}
	void setPixel16u(int x, int y, int val) 
	{ 
		assert(this->depthInBits == 16);
		unsigned short* pd = (unsigned short*)data;
		pd[x + y * (pitch >> 1)] = (unsigned short)val; 
	}
	void setPixel16u(int i, int val)
	{
		setPixel16u(i % width, i / width, val);
	}
	bool changePixels(IRoi roi, int currentValue, int newValue)
	{
		int index;
		for (int y = roi.y0; y < roi.y1; y++)
		{
			for (int x = roi.x0; x < roi.x1; x++)
			{
				index = x + y * pitch;
				if (data[index] == currentValue) data[index] = newValue;
			}
		}
		return true;
	}
	bool clear(unsigned char value = 0)
	{
		if (bufferSize <= 0) return true;
		memset(data, value, bufferSize);
		return true;
	}
	void Shift(int dx, int dy, int bg)
	{
		unsigned char* newBuffer = (unsigned char*)malloc(bufferSize);
		if (newBuffer == NULL) 
		{
			logd.error("Shift: mem alloc failed on size %d", bufferSize);
			ErrorExit("Shift: mem alloc failed");
		}
		if (getDepthInBits() == 16)
		{
			CopyAndShiftBytes16((uint16_t*)this->data, this->width, this->height, this->pitch, (uint16_t*)newBuffer, dx, dy, (uint16_t)bg);
		}
		else if (getDepthInBits() == 8)
		{
			CopyAndShiftBytes8(this->data, this->width, this->height, this->pitch, newBuffer, dx, dy, (unsigned char)bg);
		}
		else
		{
			ErrorExit("Depth not handled.");
		}
		unsigned char* oldBuffer = this->data;
		this->data = newBuffer;
		free(oldBuffer);
	}
	void CopyAndShiftBytes16(const uint16_t* srcp, int width, int height, int pitch, uint16_t* dstp, int dx, int dy, uint16_t bg)
	{
		int pp = pitch / 2; 
		int sy0, nRows, dy0;
		if (dy > 0)
		{
			sy0 = 0;
			dy0 = dy;
			nRows = height - dy;
		}
		else
		{
			sy0 = -dy;
			dy0 = 0;
			nRows = height + dy;
		}
		int sx0, dx0, nColumns;
		if (dx > 0)
		{
			sx0 = 0;
			dx0 = dx;
			nColumns = width - dx;
		}
		else
		{
			sx0 = -dx;
			dx0 = 0;
			nColumns = width + dx;
		}
		const uint16_t* sp = (srcp + sx0 + sy0 * pp);
		uint16_t* dp = (dstp + dx0 + dy0 * pp);
		for (int h = 0; h < nRows; h++)
		{
			memcpy(dp, sp, nColumns*2);
			sp += pp;
			dp += pp;
		}
		if (dy > 0)
		{
			dp = dstp;
			for (int y = 0; y < dy; y++)
			{
				for (int x = 0; x < width; x++)
				{
					dp[x] = bg;
				}
				dp += pp;
			}
		}
		else if (dy < 0)
		{
			dp = (dstp + (height + dy) * pp);
			for (int y = height+dy; y < height; y++)
			{
				for (int x = 0; x < width; x++)
				{
					dp[x] = bg;
				}
				dp += pp;
			}
		}
		if (dx > 0)
		{
			dp = dstp;
			for (int y = 0; y < height; y++)
			{
				for (int x = 0; x < dx; x++)
				{
					dp[x] = bg;
				}
				dp += pp;
			}
		}
		else if (dx < 0)
		{
			dp = dstp;
			for (int y = 0; y < height; y++)
			{
				for (int x = width + dx; x < width; x++)
				{
					dp[x] = bg;
				}
				dp += pp;
			}
		}
	}
	void CopyAndShiftBytes8(const unsigned char* srcp, int width, int height, int pitch,  unsigned char* dstp, int dx, int dy, unsigned char bg)
	{
		int sy0, nRows, dy0;
		if (dy > 0)
		{
			sy0 = 0;
			dy0 = dy;
			nRows = height - dy;
		}
		else
		{
			sy0 = -dy;
			dy0 = 0;
			nRows = height + dy;
		}
		int sx0, dx0, nColumns;
		if (dx > 0)
		{
			sx0 = 0;
			dx0 = dx;
			nColumns = width - dx;
		}
		else
		{
			sx0 = -dx;
			dx0 = 0;
			nColumns = width + dx;
		}
		const byte* sp = (srcp + sx0 + sy0 * pitch);
		byte* dp = (dstp + dx0 + dy0 * pitch);
		for (int h = 0; h < nRows; h++) 
		{
			memcpy(dp, sp, nColumns);
			sp += pitch;
			dp += pitch;
		}
		if (dy > 0)
		{
			dp = dstp;
			for (int y = 0; y < dy; y++)
			{
				for (int x = 0; x < width; x++)
				{
					dp[x] = bg;
				}
				dp += pitch;
			}
		}
		else if (dy < 0)
		{
			dp = (dstp + (height + dy) * pitch);
			for (int y = height + dy; y < height; y++)
			{
				for (int x = 0; x < width; x++)
				{
					dp[x] = bg;
				}
				dp += pitch;
			}
		}
		if (dx > 0)
		{
			dp = dstp;
			for (int y = 0; y < height; y++)
			{
				for (int x = 0; x < dx; x++)
				{
					dp[x] = bg;
				}
				dp += pitch;
			}
		}
		else if (dx < 0)
		{
			dp = dstp;
			for (int y = 0; y < height; y++)
			{
				for (int x = width + dx; x < width; x++)
				{
					dp[x] = bg;
				}
				dp += pitch;
			}
		}
	}
	Img& operator=(const Img& other)
	{
		if (!resize(other.width, other.height, other.depthInBits, other.pitch)) return(*this);
		memcpy(data, other.data, bufferSize);
		return(*this);
	}
};
#include <cstring>
class VecUtil
{
public:
	template <class T>
	static string ToString(const char* prefix, const char* fmt, std::vector<T>& values)
	{
		string s = prefix;
		char tmp[4096];
		for (size_t i = 0; i < values.size(); i++)
		{
			sprintf(tmp, fmt, values[i]);
			s += (string)tmp;
		}
		return s;
	}
};
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-parameter"
template <>
string VecUtil::ToString<bool>(const char* prefix, const char* fmt, std::vector<bool>& values)
{
	string s = prefix;
	char tmp[512];
	fmt = NULL; 
	for (size_t i = 0; i < values.size(); i++)
	{
		sprintf(tmp, "%c, ", (values[i] ? 'T' : 'F'));
		s += (string)tmp;
	}
	return s;
}
#pragma GCC diagnostic pop
class Stats
{
public:
	int n;
	uint16_t min;
	uint16_t max;
	float mean;
	float stddev;
	float meandev;
	Stats()
	{
		n = 0;
		min = max = 0;
		mean = stddev = meandev = 0;
	}
	Stats(vector<uint16_t>& values)
	{
		this->n = values.size();
		ComputeStats(values, min, max, mean, stddev, meandev);
	}
	void Compute(vector<uint16_t>& values)
	{
		this->n = values.size();
		ComputeStats(values, min, max, mean, stddev, meandev);
	}
	std::string ToString()
	{
		char tmp[512];
		sprintf(tmp, "N: %d, min: %d, mean: %.2f, max: %d, stddev: %.3f, meandev: %.3f", n, min, mean, max, stddev, meandev);
		return (string)tmp;
	}
	static std::string ToCsvHeaderString(string prefix = "")
	{
		char tmp[512];
		sprintf(tmp, "%sN, %sMin, %sMean, %sMax, %sStddev, %sMeandev", prefix.c_str(), prefix.c_str(), prefix.c_str(), prefix.c_str(), prefix.c_str(), prefix.c_str());
		return (string)tmp;
	}
	std::string ToCsvString()
	{
		char tmp[512];
		sprintf(tmp, "%d, %d, %f, %d, %f, %f", n, min, mean, max, stddev, meandev);
		return (string)tmp;
	}
	static void ComputeStats(vector<uint16_t>& values, uint16_t& min, uint16_t& max, float& mean, float& stddev, float& meandev)
	{
		int n = values.size();
		double sum = 0, s2 = 0; 
		min = 65535;
		max = 0;
		for (int i = 0; i < n; i++) 
		{
			min = MIN(min, values[i]);
			max = MAX(max, values[i]);
			sum += values[i];
			s2 += values[i] * values[i];
		}
		mean = (float)(sum / n);
		double t1 = ((n * s2) - sum * sum);
		if (t1 > 0)
		{
			stddev = sqrt(t1 / ((double)n * (n - 1))); 
		}
		else
		{
			stddev = 0;
		}
		sum = 0;
		for (int i = 0; i < n; i++) sum += abs(values[i] - mean);
		meandev = (float)(sum / n);
	}
	static void ComputeStats(vector<float>& values, float& min, float& max, float& mean, float& stddev, float& meandev)
	{
		int n = values.size();
		double sum = 0, s2 = 0;
		min = FLT_MAX;
		max = -FLT_MAX;
		for (int i = 0; i < n; i++) 
		{
			if (values[i] < -998)
			{
				ErrorExit("Probable magic value coming into singleCluster");
			}
			min = MIN(min, values[i]);
			max = MAX(max, values[i]);
			sum += values[i];
			s2 += values[i] * values[i];
		}
		mean = (float)(sum / n);
		double t1 = (double)((n * s2) - sum * sum);
		if (t1 > 0)
		{
			stddev = (float)(sqrt(t1 / ((double)n * (n - 1))));
		}
		else
		{
			stddev = 0;
		}
		sum = 0;
		for (int i = 0; i < n; i++) sum += abs(values[i] - mean);
		meandev = (float)(sum / n);
	}
	static void ComputeStats(vector<float>& values, vector<float>& weights, float& min, float& max, float& mean, float& stddev, float& meandev)
	{
		int n = values.size();
		double sum = 0, s2 = 0;
		double sumWeight = 0;
		min = FLT_MAX;
		max = -FLT_MAX;
		for (int i = 0; i < n; i++) 
		{
			if (weights[i] > 0)
			{
				if (values[i] < -998)
				{
					assert(false && "Probable magic value coming into singleCluster");
					values[i] = 0;
				}
				min = MIN(min, values[i]);
				max = MAX(max, values[i]);
			}
			sumWeight += weights[i];
			sum += weights[i] * values[i];
			s2 += weights[i] * values[i] * values[i];
		}
		mean = (float)(sum / sumWeight);
		double t1 = ((sumWeight * s2) - sum * sum);
		if (t1 > 0)
		{
			stddev = (float)(sqrt(t1 / (sumWeight * (sumWeight - 1))));
		}
		else
		{
			stddev = 0;
		}
		sum = 0;
		for (int i = 0; i < n; i++) sum += abs(weights[i] * (values[i] - mean));
		meandev = (float)(sum / sumWeight);
	}
	static void testComputeWeighted()
	{
		vector<float> values;
		vector<float> weights;
		int n = 5;
		for (int i = 0; i < n; i++) 
		{
			values.push_back(i);
			weights.push_back(i);
		}
		float min, max, mean, stddev, meandev;
		Stats::ComputeStats(values, weights, min, max, mean, stddev, meandev);
		logd.printf("Values: %s\n", VecUtil::ToString("", "%.3f, ", values).c_str());
		logd.printf("Weights: %s\n", VecUtil::ToString("", "%.3f, ", weights).c_str());
		logd.printf("Stats: min: %.3f, max: %.3f, mean: %.3f, stddev: %.3f, meandev: %.3f\n", min, max, mean, stddev, meandev);
	}
};
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
class FileUtil
{
public:
	static bool WriteFile(string fileName, vector<string>& strings, bool addNewlines)
	{
		FILE* out = fopen(fileName.c_str(), "w");
		if (out == NULL)
		{
			cout << "Error: Unable to open and write file." << endl;
			return(false);
		}
		for (int i = 0; i < strings.size(); i++)
		{
			fprintf(out, strings[i].c_str());
			if (addNewlines) fprintf(out, "\n");
		}
		fclose(out);
		return true;
	}
	static bool GetFileSize(const char* pathname, int64_t* size)
	{
		ifstream in(pathname, ios::in | ios::binary);
		if (!in)
		{
			logd.printf("Error: FileUtil: GetFileSize: Unable to open file '%s'\n", pathname);
			return false;;
		}
		in.seekg(0, ios::end);
		*size = (int64_t)in.tellg();
		in.close();
		return true;
	}
	static bool SnarfFile(const char* pathname, int64_t length, unsigned char* buffer, bool binaryMode = true)
	{
		int64_t realFileSize;
		if (!GetFileSize(pathname, &realFileSize))
		{
			logd.printf("Error: Unable to get file size for file '%s'.\n", pathname);
			return false;;
		}
		if (realFileSize > length)
		{
			logd.printf("Length argument %I64d is less than real file size %I64d.", length, realFileSize);
			return false;;
		}
		char mode[3];
		if (binaryMode)	strcpy(mode, "rb");
		else strcpy(mode, "r");
		FILE* fptr = fopen(pathname, mode);
		if (fptr == NULL)
		{
			logd.printf("Error: Unable to open and read file '%s'.\n", pathname);
			return false;;
		}
		int readCount = fread(buffer, 1, realFileSize, fptr);
		fclose(fptr);
		if (readCount < realFileSize)
		{
			logd.printf("Error: Unexpected end of file at %d", readCount);
			return false;;
		}
		return true;
	}
	static bool SnarfFile(string path, vector<string>& strings)
	{
		return SnarfFile(path.c_str(), strings);
	}
	static bool SnarfFile(const char* pathname, vector<string>& strings)
	{
		int64_t realFileSize;
		if (!GetFileSize(pathname, &realFileSize))
		{
			char tmp[512];
			sprintf(tmp, "Error: Unable to get file size: %s", pathname);
			ErrorExit(tmp);
		}
		if (realFileSize < 0)
		{
			char tmp[512];
			sprintf(tmp, "Error: Unable to get file size: %s", pathname);
			ErrorExit(tmp);
		}
		unsigned char* buffer = new unsigned char[(uint64_t)realFileSize + 1];
		buffer[realFileSize] = 0;
		if (buffer == NULL)
		{
			char tmp[512];
			sprintf(tmp, "Error: Unable to allocate buffer for whole file of size: %lld", realFileSize);
			ErrorExit(tmp);
		}
		if (!SnarfFile(pathname, realFileSize, buffer))
		{
			char tmp[512];
			sprintf(tmp, "Error: Unable to snarf file");
			logd.printf("%s\n", tmp);
			delete [] buffer;
			ErrorExit(tmp);
		}
		char c;
		char line[512];
		int subLineIndex = 0;
		int index = 0;
		while((c = buffer[index++]) != 0)
		{
			if (c == 10 || c == 13)  
			{
				if (c == 13) index++; 
				line[subLineIndex] = (char)0;
				strings.push_back(string(line));
				subLineIndex = 0;
			}
			else
			{
				line[subLineIndex++] = c;
			}
		}
		if (subLineIndex != 0)
		{
			line[subLineIndex] = (char)0;
			strings.push_back(string(line));
		}
		delete [] buffer;
		return true;
	}
};
#include <vector>
#include <algorithm>
template <class T>
class sort_indices
{
	private:
		std::vector<T>& pv;
	public:
		sort_indices(vector<T>& parr) : pv(parr) {}
		bool operator()(int i, int j) { return pv[i] < pv[j]; }
};
class SortUtil
{
public:
	template <class T>
	static void Sort(std::vector<T>& v)
	{
		std::sort(v.begin(), v.end());
	}
	template <class T>
	static void SortIndices(std::vector<T>& v, std::vector<int>& indices)
	{
		int n = v.size();
		indices.resize(n);
		for (int i = 0; i < n; i++) indices[i] = i;
		std::sort(indices.begin(), indices.begin() + n, sort_indices<T>(v));
	}
	template <class T1, class T2>
	static void ParallelSort(std::vector<T1>& v1, std::vector<T2>& v2)
	{
		assert(v1.size() == v2.size());
		vector<int> indices;
		SortIndices(v1, indices);
		int n = v1.size();
		vector<T1> tmp1(n);
		vector<T2> tmp2(n);
		for (int i = 0; i < n; i++)
		{
			tmp1[i] = v1[indices[i]];
			tmp2[i] = v2[indices[i]];
		}
		v1 = tmp1;
		v2 = tmp2;
	}
};
#include <tuple>
#include <vector>
class PointIndex
{
	int origSideLength; 
	int sideLength; 
	int cellSideLength; 
	int dim; 
	vector<vector<int>> gridCellsKeys;
	vector<vector<pair<int, int>>> gridCellsPoints;
public:
	PointIndex()
	{
		origSideLength = sideLength = cellSideLength = dim = -1;
	}
	PointIndex(int sideLength, int gridCellSideLength)
	{
		Init(sideLength, gridCellSideLength);
	}
	void Init(int sideLength, int gridCellSideLength)
	{
		if (sideLength <= 0) ErrorExit("sideLength must be gt 0");
		if (sideLength <= gridCellSideLength) ErrorExit("sideLength must be gt gridCellSideLength");
		this->dim = (sideLength - 1) / gridCellSideLength + 1;
		this->origSideLength = sideLength;
		this->sideLength = this->dim * gridCellSideLength;
		this->cellSideLength = gridCellSideLength;
		gridCellsKeys.clear();
		gridCellsPoints.clear();
		int nCells = dim * dim;
		gridCellsKeys.resize(nCells);
		gridCellsPoints.resize(nCells);
	}
	int GetGridDim()
	{
		return dim;
	}
	void Add(pair<int, int> pt, int key)
	{
		if (pt.first < 0) ErrorExit("This requires points be gte 0,0");
		if (pt.second < 0) ErrorExit("This requires points be gte 0,0");
		if (pt.first >= sideLength) ErrorExit("This requires points be lt sideLength");
		if (pt.second >= sideLength) ErrorExit("This requires points be lt sideLength");
		int idx = GetCellIndexFromPoint(pt);
		gridCellsKeys[idx].push_back(key);
		gridCellsPoints[idx].push_back(pt);
	}
	void AddRound(float x, float y, int key)
	{
		pair<int, int> ipt = pair<int, int>(roundf(x), roundf(y));
		Add(ipt, key);
	}
	void Add(vector<pair<int, int>>& points)
	{
		int n = points.size();
		for (int i = 0; i < n; i++)
		{
			Add(points[i], i);
		}
	}
	void AddRound(vector<pair<float, float>>& points)
	{
		int n = points.size();
		for (int i = 0; i < n; i++)
		{
			pair<int, int> ipt = pair<int, int>(roundf(points[i].first), roundf(points[i].second));
			Add(ipt, i);
		}
	}
	vector<int> GetBlobsInGridCell(int xCellIndex, int yCellIndex)
	{
		int idx = GetCellIndexFromCellCoords(xCellIndex, yCellIndex);
		return gridCellsKeys[idx];
	}
	vector<int> FindNearby(float x, float y, float distThreshold, bool debug = false)
	{
		return FindNearby(roundf(x), roundf(y), distThreshold, false, debug);
	}
	vector<int> FindNearby(pair<float, float> pt, float distThreshold, bool debug = false)
	{
		return FindNearby(roundf(pt.first), roundf(pt.second), distThreshold, false, debug);
	}
	vector<int> FindNearby(pair<int, int> pt, float distThreshold, bool debug = false)
	{
		return FindNearby(pt.first, pt.second, distThreshold, false, debug);
	}
	vector<int> FindNearby(int xc, int yc, float distThreshold, bool doSort, bool debug)
	{
		float distThresholdSq = distThreshold * distThreshold;
		vector<int> v;
		vector<float> distances;
		int x0 = xc - (int)ceil(distThreshold);
		int y0 = yc - (int)ceil(distThreshold);
		int x1 = xc + (int)ceil(distThreshold);
		int y1 = yc + (int)ceil(distThreshold);
		x0 = MAX(0, x0);
		y0 = MAX(0, y0);
		x1 = MIN(origSideLength - 1, x1); 
		y1 = MIN(origSideLength - 1, y1);
		if (debug) logd.printf("FindNearby: User-space bounds: (%d, %d) to (%d, %d)\n", x0, y0, x1, y1);
		int xg0 = x0 / cellSideLength;
		int yg0 = y0 / cellSideLength;
		int xg1 = x1 / cellSideLength;
		int yg1 = y1 / cellSideLength;
		if (debug) logd.printf("FindNearby: Grid-space bounds: (%d, %d) to (%d, %d)\n", xg0, yg0, xg1, yg1);
		for (int y = yg0; y <= yg1; y++)
		{
			for (int x = xg0; x <= xg1; x++)
			{
				int idx = GetCellIndexFromCellCoords(x, y);
				if (debug) logd.printf("FindNearby: Search grid cell: %d,%d  (idx: %d)\n", x, y, idx);
				int n = gridCellsKeys[idx].size();
				for (int i = 0; i < n; i++)
				{
					pair<int, int>& p1 = gridCellsPoints[idx][i];
					float dx = p1.first - xc;
					float dy = p1.second - yc;
					float distSq = dx*dx + dy*dy;
					if (distSq <= distThresholdSq)
					{
						if (debug) logd.printf("    Close enough to pt: (%d,%d), dist: %.1f\n", p1.first, p1.second, sqrtf(distSq));
						v.push_back(gridCellsKeys[idx][i]);
						distances.push_back(distSq);
					}
					else
					{
						if (debug) logd.printf("    Too far to pt: (%d,%d), dist: %.1f\n", p1.first, p1.second, sqrtf(distSq));
					}
				}
			}
		}
		if (doSort)
		{
			SortUtil::ParallelSort(distances, v);
		}
		return v;
	}
	void Display(int mode = 0)
	{
		logd.printf("PointIndex: %dx%d length, %dx%d cells (each %dx%d)\n", sideLength, sideLength, dim, dim, cellSideLength, cellSideLength);
		for (int y = 0; y < dim; y++)
		{
			for (int x = 0; x < dim; x++)
			{
				pair<int, int> pt(x, y);
				int idx = GetCellIndexFromCellCoords(pt);
				if (mode == 0)
				{
					logd.printf("%4d, ", gridCellsKeys[idx].size());
				}
				else
				{
					if (gridCellsKeys[idx].size() > 0)
					{
					}
					else
					{
						logd.printf("(--),  ");
					}
				}
			}
			logd.printf("\n");
		}
	}
	pair<int, int> GetCellCoordsFromUserSpaceCoords(int x, int y)
	{
		if ((x < 0) || (x >= sideLength)) ErrorExit("GetCellCoordsFromUserSpaceCoords: x %d is out of range", x);
		if ((y < 0) || (y >= sideLength)) ErrorExit("GetCellCoordsFromUserSpaceCoords: y %d is out of range", y);
		int idx = GetCellIndexFromPoint(pair<int, int>(x, y));
		if ((idx < 0) || (idx >= gridCellsKeys.size())) ErrorExit("GetCellCoordsFromUserSpaceCoords: idx %d is out of range", idx);
		return GetCellCoords(idx);
	}
private:
	inline int GetCellIndexFromPoint(int x, int y)
	{
		int xc = (int)(x / cellSideLength);
		int yc = (int)(y / cellSideLength);
		return xc + dim * yc;
	}
	inline int GetCellIndexFromPoint(pair<int, int> pt)
	{
		int xc = (int)(pt.first / cellSideLength);
		int yc = (int)(pt.second / cellSideLength);
		return xc + dim * yc;
	}
	inline int GetCellIndexFromCellCoords(pair<int, int> pt)
	{
		return pt.first + dim * pt.second;
	}
	inline int GetCellIndexFromCellCoords(int x, int y)
	{
		return x + dim * y;
	}
	inline pair<int, int> GetCellCoords(int idx)
	{
		if ((idx < 0) || (idx >= gridCellsKeys.size())) ErrorExit("GetCellCoords: idx %d is out of range", idx);
		pair<int, int> pt;
		pt.first = idx % dim;
		pt.second = idx / dim;
		return pt;
	}
	vector<int> GetCellNeighbors(int cellIndex)
	{
		vector<int> v;
		int offsets[8] = {-dim - 1, -dim, -dim + 1, -1, 1, dim - 1, dim, dim + 1};
		int nCells = dim * dim;
		for (int i = 0; i < 8; i++)
		{
			int idx = cellIndex + offsets[i];
			if ((idx >= 0) && (idx < nCells)) v.push_back(idx);
		}
		return v;
	}
	vector<int> GetKeysInCell(int x, int y)
	{
		int idx = GetCellIndexFromCellCoords(x, y);
		return gridCellsKeys[idx];
	}
public:
};
class Transform2d
{
	float m[6];
public:
	Transform2d();
	void Clear();
	void SetScale(float xZoom, float yZoom);
	void SetTranslation(float dx, float dy);
	void SetRotation(float radians);
	Transform2d Mult(Transform2d& other);
	void Concat(Transform2d& other);
	void Transform(float& x, float& y);
	pair<float,float> TransformPt(float x, float y);
	string ToString()
	{
		char tmp[2048];
		sprintf(tmp, "%.3f:%.3f:%.3f_%.3f:%.3f:%.3f", m[0], m[1], m[2], m[3], m[4], m[5]);
		return (string)tmp;
	}
	void Display();
	static bool TestSelf();
};
Transform2d::Transform2d()
{
	for (int i = 0; i < 6; i++) m[i] = 0;
	m[0] = m[4] = 1;
}
void Transform2d::Clear()
{
	for (int i = 0; i < 6; i++) m[i] = 0;
	m[0] = m[4] = 1;
}
void Transform2d::SetTranslation(float dx, float dy)
{
	Clear();
	m[2] = dx;
	m[5] = dy;
}
void Transform2d::SetScale(float xZoom, float yZoom)
{
	Clear();
	m[0] = xZoom;
	m[4] = yZoom;
}
void Transform2d::SetRotation(float radians)
{
	Clear();
	m[0] = m[4] = cos(radians);
	m[3] = sin(radians);
	m[1] = -m[3];
}
void Transform2d::Transform(float& x, float& y)
{
	float xp = m[0] * x + m[1] * y + m[2];
	y = m[3] * x + m[4] * y + m[5];
	x = xp;
}
pair<float,float> Transform2d::TransformPt(float x, float y)
{
	pair<float, float> pt;
	float xp = m[0] * x + m[1] * y + m[2];
	pt.second = m[3] * x + m[4] * y + m[5];
	pt.first = xp;
	return pt;
}
void Transform2d::Concat(Transform2d& other)
{
	*this = other.Mult(*this);
}
void Transform2d::Display()
{
	printf("%.3f, %.3f, %.3f\n", m[0], m[1], m[2]);
	printf("%.3f, %.3f, %.3f\n", m[3], m[4], m[5]);
	printf("%.3f, %.3f, %.3f\n", 0.0f, 0.0f, 1.0f);
}
Transform2d Transform2d::Mult(Transform2d& other)
{
	Transform2d out;
	out.m[0] = m[0] * other.m[0] + m[1] * other.m[3];
	out.m[1] = m[0] * other.m[1] + m[1] * other.m[4];
	out.m[2] = m[0] * other.m[2] + m[1] * other.m[5] + m[2];
	out.m[3] = m[3] * other.m[0] + m[4] * other.m[3];
	out.m[4] = m[3] * other.m[1] + m[4] * other.m[4];
	out.m[5] = m[3] * other.m[2] + m[4] * other.m[5] + m[5];
	return out;
}
bool Transform2d::TestSelf()
{
	return true;
}
#include <string>
#include <vector>
class Feature
{
public:
	string name;
	float value;
	bool isSelected;
	int csvIndex; 
	Feature()
	{
		name = "";
		value = 0;
		isSelected = true;
		csvIndex = -1;
	}
	Feature(string& n, float v)
	{
		name = n;
		value = v;
		isSelected = true;
		csvIndex = -1;
	}
	Feature(const char* n, float v)
	{
		name = n;
		value = v;
		isSelected = true;
		csvIndex = -1;
	}
	std::string ToString()
	{
		char tmp[128];
		sprintf(tmp, "name: %s, val: %.3f, sel: %c", name.c_str(), value, (isSelected ? 'Y' : 'N'));
		return (string)tmp;
	}
	static vector<float> GetSelectedValues(vector<Feature>& v, vector<Feature>& selectionPrototype)
	{
		int nf = 0;
		for (int i = 0; i < v.size(); i++)
		{
			if (selectionPrototype[i].isSelected) nf++;
		}
		vector<float> values(nf);
		int k = 0;
		for (int j = 0; j < v.size(); j++)
		{
			if (selectionPrototype[j].isSelected) values[k++] = v[j].value;
		}
		return values;
	}
	static vector<float> GetAllValues(vector<Feature>& v)
	{
		int nf = v.size();
		vector<float> values(nf);
		for (int j = 0; j < v.size(); j++)
		{
			values[j] = v[j].value;
		}
		return values;
	}
	static vector<string> GetSelectedFeatureNames(vector<Feature>& v)
	{
		int nf = 0;
		for (int i = 0; i < v.size(); i++)
		{
			if (v[i].isSelected) nf++;
		}
		vector<string> names(nf);
		int k = 0;
		for (int j = 0; j < v.size(); j++)
		{
			if (v[j].isSelected) names[k++] = v[j].name;
		}
		return names;
	}
};
class Wcs
{
private:
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
	double D2R = PI_DOUBLE / 180.0;
	double R2D = 180.0 / PI_DOUBLE;
public:
	void init(vector<double>& wcs)
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
	void update_invcd()
	{
		double inv_det = cd11*cd22 - cd12*cd21;
		inv_det = 1.0 / inv_det;
		invcd11 = cd22 * inv_det;
		invcd12 = -cd12 * inv_det;
		invcd21 = -cd21 * inv_det;
		invcd22 = cd11 * inv_det;
	}
	vector<double> sphs2x(double lng, double lat)
	{
		double coslat, coslat3, coslat4, coslng, dlng, dphi, sinlat, sinlat3, sinlat4, sinlng, x, y, z;
		vector<double> eul(5);
		vector<double> phiTheta(2);
		eul[0] = crval1;
		eul[1] = 90.0 - crval2;
		eul[2] = 180.0;
		eul[3] = cos(D2R * eul[1]);
		eul[4] = sin(D2R * eul[1]);
		dlng = lng - eul[0];
		phiTheta[0] = dlng;
		sinlat = sin(lat*D2R);
		coslat = cos(lat*D2R);
		coslat3 = coslat*eul[3];
		coslat4 = coslat*eul[4];
		sinlat3 = sinlat*eul[3];
		sinlat4 = sinlat*eul[4];
		dlng = phiTheta[0];
		sinlng = sin(dlng*D2R);
		coslng = cos(dlng*D2R);
		x = sinlat4 - coslat3*coslng;
		y = -coslat*sinlng;
		dphi = R2D*atan2(y, x);
		phiTheta[0] = fmod(eul[2] + dphi, 360.0);
		if (phiTheta[0] > 180.0)
		{
			phiTheta[0] -= 360.0;
		}
		else if (phiTheta[0] < -180.0)
		{
			phiTheta[0] += 360.0;
		}
		z = sinlat3 + coslat4*coslng;
		if (fabs(z) > 0.99)
		{
			if (z < 0.0)
				phiTheta[1] = -fabs(R2D * acos(sqrt(x*x + y*y)));
			else
				phiTheta[1] = fabs(R2D * acos(sqrt(x*x + y*y)));
		}
		else
		{
			phiTheta[1] = R2D * asin(z);
		}
		return phiTheta;
	}
	vector<double> sphx2s(double phi, double theta)
	{
		double cosphi, costhe, costhe3, costhe4, dlng, dphi, sinphi, sinthe, sinthe3, sinthe4, x, y, z;
		vector<double> lngLat(2);
		vector<double> eul(5);
		eul[0] = crval1;
		eul[1] = 90.0 - crval2;
		eul[2] = 180.0;
		eul[3] = cos(D2R * eul[1]);
		eul[4] = sin(D2R * eul[1]);
		dphi = phi - eul[2];
		sinthe = sin(theta * D2R);
		costhe = cos(theta * D2R);
		costhe3 = costhe*eul[3];
		costhe4 = costhe*eul[4];
		sinthe3 = sinthe*eul[3];
		sinthe4 = sinthe*eul[4];
		sinphi = sin(dphi * D2R);
		cosphi = cos(dphi * D2R);
		x = sinthe4 - costhe3*cosphi;
		y = -costhe*sinphi;
		dlng = R2D * atan2(y, x);
		lngLat[0] = eul[0] + dlng;
		if (eul[0] >= 0.0)
		{
			if (lngLat[0] < 0.0)
				lngLat[0] += 360.0;
		}
		else
		{
			if (lngLat[0] > 0.0)
				lngLat[0] -= 360.0;
		}
		if (lngLat[0] > 360.0)
		{
			lngLat[0] -= 360.0;
		}
		else if (lngLat[0] < -360.0)
		{
			lngLat[0] += 360.0;
		}
		z = sinthe3 + costhe4*cosphi;
		if (fabs(z) > 0.99)
		{
			if (z<0.0)
				lngLat[1] = -fabs(R2D * acos(sqrt(x*x + y*y)));
			else
				lngLat[1] = fabs(R2D * acos(sqrt(x*x + y*y)));
		}
		else
		{
			lngLat[1] = R2D * asin(z);
		}
		return lngLat;
	}
	vector<double> tans2x(double phi, double theta)
	{
		vector<double> xy(2);
		double cotan_theta = 1.0 / tan(D2R * theta);
		xy[0] = R2D * sin(D2R * phi) * cotan_theta;
		xy[1] = -R2D * cos(D2R * phi) * cotan_theta;
		return xy;
	}
	vector<double> tanx2s(double x, double y)
	{
		vector<double> phiTheta(2);
		phiTheta[0] = R2D * atan2(x, -y);
		phiTheta[1] = R2D * atan2(R2D, sqrt(x*x + y*y));
		return phiTheta;
	}
	vector<double> convertRADEC2XY(double RA, double DEC)
	{
		vector<double> XY(2);
		vector<double> phiTheta = sphs2x(RA, DEC);
		vector<double> XXYY = tans2x(phiTheta[0], phiTheta[1]);
		XY[0] = invcd11*XXYY[0] + invcd12*XXYY[1] + crpix1;
		XY[1] = invcd21*XXYY[0] + invcd22*XXYY[1] + crpix2;
		return XY;
	}
	vector<double> convertXY2RADEC(double x, double y)
	{
		double dx = x - crpix1;
		double dy = y - crpix2;
		double xx = cd11*dx + cd12*dy;
		double yy = cd21*dx + cd22*dy;
		vector<double> phiTheta = tanx2s(xx, yy);
		vector<double> RD = sphx2s(phiTheta[0], phiTheta[1]);
		return RD;
	}
};
struct DebugSpec
{
	int debugImageSetIndex;
	int debugDetNumber;
	int debugFrameIndex;
	int caIndex;
	int maxToProcess;
	bool debug; 
	bool debugc; 
	bool info; 
	string casFeatureCsv; 
	DebugSpec()
	{
		debugImageSetIndex = debugDetNumber = caIndex = maxToProcess = -1;
		debug = debugc = false;
		info = true;
	}
	pair<int, int> debugBlobCoords; 
	bool CheckDebugCA(int imageSetIndex, int caIndex, int detNum = -1)
	{
		return CheckDebugImageSet(imageSetIndex)
			&& ((caIndex == this->caIndex) || ((detNum >= 0) && (detNum == this->debugDetNumber)));
	}
	bool CheckDebugDet(int imageSetIndex, int detNum)
	{
		return CheckDebugImageSet(imageSetIndex)
			&& (detNum >= 0) && (detNum == this->debugDetNumber);
	}
	bool CheckDebugImageSet(int imageSetIndex)
	{
		return (debugImageSetIndex == -2) || (debugImageSetIndex == imageSetIndex);
	}
	bool CheckDebugImageFrame(int imageSetIndex, int frameIndex)
	{
		return CheckDebugImageSet(imageSetIndex) && ((debugFrameIndex == -2) || (debugFrameIndex == frameIndex));
	}
};
class Img8Blob
{
public:
	int xCentroid; 
	int yCentroid; 
	int cx; 
	int cy; 
	byte label;
	int area;
	IRoi roi;
	string ToString()
	{
		char tmp[256];
		sprintf(tmp, "label: %d, centroid: (%d,%d), contact: (%d,%d), area: %d", label, xCentroid, yCentroid, cx, cy, area);
		return (string)tmp;
	}
};
class ImgBlob
{
public:
	float xCentroid; 
	float yCentroid; 
	int cx; 
	int cy; 
	uint16_t label;
	int area;
	IRoi roi;
	ImgBlob() : roi()
	{
		xCentroid = yCentroid = 0;
		cx = cy = area = 0;
		label = 0;
	}
	float Dist(pair<float, float> pt)
	{
		return sqrtf(DistSq(pt));
	}
	float DistInt(pair<int, int> pt)
	{
		return sqrtf(DistIntSq(pt));
	}
	float DistSq(ImgBlob& other)
	{
		pair<float, float> pt;
		pt.first = other.xCentroid;
		pt.second = other.yCentroid;
		return DistSq(pt);
	}
	float DistSq(pair<float, float> pt)
	{
		float dx = pt.first - xCentroid;
		float dy = pt.second - yCentroid;
		return dx * dx + dy * dy;
	}
	float DistIntSq(pair<int, int> pt)
	{
		float dx = (float)pt.first - xCentroid;
		float dy = (float)pt.second - yCentroid;
		return dx * dx + dy * dy;
	}
	static void AddAllToPointIndex(vector<ImgBlob>& blobs, PointIndex& pi)
	{
		int n = blobs.size();
		for (int i = 0; i < n; i++)
		{
			pi.AddRound(blobs[i].xCentroid, blobs[i].yCentroid, i);
		}
	}
	bool CheckDebug(DebugSpec dspec)
	{
		return (abs(xCentroid - dspec.debugBlobCoords.first) <= DEBUG_COORDS_TOLERANCE_PIX)
			&& (abs(yCentroid - dspec.debugBlobCoords.second) <= DEBUG_COORDS_TOLERANCE_PIX);
	}
	string ToString()
	{
		char tmp[256];
		sprintf(tmp, "label: %d, ctd: (%.1f,%.1f), area: %d", label, xCentroid, yCentroid, area);
		return (string)tmp;
	}
};
#include <string>
#include <vector>
#include <string>
#include <vector>
class RaDec
{
public:
	double RaDegrees; 
	double DecDegrees; 
	RaDec()
	{
		RaDegrees = DecDegrees = 0;
	}
	RaDec(double ra, double dec)
	{
		RaDegrees = ra;
		DecDegrees = dec;
	}
	double Dist(RaDec other)
	{
		double dx = other.RaDegrees - this->RaDegrees;
		double dy = other.DecDegrees - this->DecDegrees;
		return pyth(dx, dy);
	}
	double DistPixelsApprox(RaDec other)
	{
		return Dist(other) / MEAN_PX_PER_DEGREE;
	}
	std::string ToString()
	{
		char tmp[512];
		sprintf(tmp, "%f/%f", RaDegrees, DecDegrees);
		return (string)tmp;
	}
};
#include <string>
class JulianTime
{
public:
	double days;
	JulianTime()
	{
		days = 0;
	}
	std::string ToString()
	{
		char tmp[512];
		sprintf(tmp, "%.7f", days);
		return (string)tmp;
	}
	static bool TryParse(std::string& s, JulianTime& t)
	{
		vector<string> words;
		size_t period = s.find('.', 0);
		if ((period > 0) && (period != s.npos))
		{
			t.days = atof(s.c_str());
			return true;
		}
		else
		{
			return false;
		}
	}
	void ConvertModifiedToNormal()
	{
		days += 2400000.5;
	}
};
#include <string>
#include <vector>
class ImgMoment
{
private:
	float m00, m01, m10, m11, m02, m20;
public:
	float xCentroid, yCentroid; 
	float u00, u01, u10, u11, u02, u20;
	float mTheta; 
	float mElong;
	float wEllipse,lEllipse;
	float mRadius; 
	ImgMoment()
	{
		m00 = m01 = m10 = m11 = m02 = m20 = 0;
		u00 = u01 = u10 = u11 = u02 = u20 = 0;
		xCentroid = yCentroid = 0;
		mTheta = mElong = 0;
		wEllipse = lEllipse = mRadius = 0;
	}
	bool CheckAllGood()
	{
		return mElong > 0;
	}
	float ComputeRoundness()
	{
		if (mElong <= 0)
		{
			return 0;
		}
		else
		{
			return 1.0f / mElong;
		}
		float musum = u20 + u02;
		float mudiff = u20 - u02;
		float round = sqrtf(mudiff * mudiff + 4.0 * u11 * u11) / musum;
		return round;
	}
	static float angleDeltaDegrees(float a, float b)
	{
		float da = abs(a - b);
		if (da > 90) da = 180 - da;
		return da;
	}
	void Compute(Img& image, IRoi& roi)
	{
		for (int y = roi.y0; y < roi.y1; y++) 
		{
			float y2 = y * y;
			for (int x = roi.x0; x < roi.x1; x++) 
			{
				byte v = image.getPixel(x, y);
				m00 += v;
				m10 += x * v;
				m01 += y * v;
				m11 += x * y * v;
				m20 += x * x * v;
				m02 += y2 * v;
			}
		}
		if (m00 > 0)
		{
			xCentroid = m10/m00;
			yCentroid = m01/m00;
			u00 = m00;
			u01 = u10 = 0;
			u11 = m11 - xCentroid * m01; 
			u20 = m20 - xCentroid * m10; 
			u02 = m02 - yCentroid * m01; 
			float up20 = u20/u00;
			float up02 = u02/u00;
			float up11 = u11/u00;
			float v1 = 2 * up11;
			float v2 = up20 - up02;
			if (v2 != 0)
			{
				mTheta = radiansToDegreesF(0.5f * atan2(v1, v2));
			}
			float a = m20 / m00 - xCentroid * xCentroid;
			float b = 2 * (m11 / m00 - xCentroid * yCentroid);
			float c = m02 / m00 - yCentroid * yCentroid;
			float b2 = b * b;
			float ac2 = (a-c)*(a-c);
			float we2 = 6 * (a + c - sqrtf(b2 + ac2));
			float le2 = 6 * (a + c + sqrtf(b2 + ac2));
			if (we2 > numeric_limits<float>::epsilon())
			{
				wEllipse = sqrtf(we2);
				if (le2 > numeric_limits<float>::epsilon())
				{
					lEllipse = sqrtf(le2);
					float xe = MAX(wEllipse, lEllipse);
					float ne = MIN(wEllipse, lEllipse);
					if (ne > numeric_limits<float>::epsilon()) 
					{
						mElong = xe / ne;
					}
				}
			}
			mRadius = (wEllipse + lEllipse) / 2;
		}
	}
	float ComputeSimScore(ImgMoment& other, float& sumWeight, bool debug)
	{
		sumWeight = 0;
		if (!CheckAllGood() || !other.CheckAllGood())
		{
			sumWeight = 0;
			return 1;
		}
		float delta, deltaPct;
		float totalWeightSum = 0;
		if (debug) logd.printf("ComputeSimScore -----------\n");
		float radiusMean = (mRadius + other.mRadius) / 2;
		float wEllipseMean = (wEllipse + other.wEllipse) / 2;
		delta = abs(wEllipse - other.wEllipse);
		deltaPct = delta / wEllipseMean;
		float wEllipseScore = 1.0f - rangeScore(0.1f, 0.66f, deltaPct);
		float wEllipseWeight = 1.0f;
		totalWeightSum += 1;
		if (debug) logd.printf("wEllipse1: %.2f, wEllipse2: %.2f, delta: %.2f, wEllipseMean: %.2f, deltaPct: %.2f, score: %.2f, weight: %.2f\n", wEllipse, other.wEllipse, delta, wEllipseMean, deltaPct, wEllipseScore, wEllipseWeight);
		float lEllipseMean = (lEllipse + other.lEllipse) / 2;
		delta = abs(lEllipse - other.lEllipse);
		deltaPct = delta / lEllipseMean;
		float lEllipseScore = 1.0f - rangeScore(0.1f, 0.66f, deltaPct);
		float lEllipseWeight = 1.0f;
		totalWeightSum += 1;
		if (debug) logd.printf("lEllipse1: %.2f, lEllipse2: %.2f, delta: %.2f, lEllipseMean: %.2f, deltaPct: %.2f, score: %.2f, weight: %.2f\n", lEllipse, other.lEllipse, delta, lEllipseMean, deltaPct, lEllipseScore, lEllipseWeight);
		float elongMean = (mElong + other.mElong) / 2;
		delta = abs(mElong - other.mElong);
		deltaPct = delta / elongMean;
		float elongScore = 1.0f - rangeScore(0.1f, 0.5f, deltaPct);
		float elongThetaWeight = rangeScore(2.0f, 5.0f, radiusMean); 
		totalWeightSum += 1;
		if (debug) logd.printf("elong1: %.2f, elong2: %.2f, delta: %.2f, elongMean: %.2f, deltaPct: %.2f, score: %.2f, weight: %.2f\n", mElong, other.mElong, delta, elongMean, deltaPct, elongScore, elongThetaWeight);
		delta = angleDeltaDegrees(mTheta, other.mTheta);
		deltaPct = delta / 90;
		float thetaScore = 1.0f - rangeScore(0.1f, 0.4f, deltaPct);
		float thetaWeight = elongThetaWeight * rangeScore(1.2f, 1.6f, elongMean); 
		totalWeightSum += 1;
		if (debug) logd.printf("theta1: %.2f, theta2: %.2f, delta: %.2f, deltaPct: %.2f, weight: %.2f, score: %.2f, weight: %.2f\n", mTheta, other.mTheta, delta, deltaPct, thetaScore, thetaWeight);
		float sum = wEllipseScore * wEllipseWeight + lEllipseScore * lEllipseWeight + elongScore * elongThetaWeight + thetaScore * thetaWeight;
		sumWeight = wEllipseWeight + lEllipseWeight + elongThetaWeight + thetaWeight;
		float wm = sum / sumWeight;
		sumWeight /= totalWeightSum; 
		if (debug) logd.printf("sum: %.2f, sumWeight: %.2f, score: %.2f\n", sum, sumWeight, wm);
		return wm;
	}
	std::string ToString()
	{
		char tmp[4096];
		sprintf(tmp, "sum: %.2f, mean: (%.1f, %.1f), theta: %.2f, elong: %.2f, ell: (%.1f, %.1f), radius: %.1f", m00, xCentroid, yCentroid, mTheta, mElong, wEllipse, lEllipse, mRadius);
		return (string)tmp;
	}
	static std::string ToCsvHeaderString(string baseName)
	{
		char tmp[4096];
		sprintf(tmp, "%sSum, %sXMean, %sYMean, %sTheta, %sElong, %sWEllipse, %sLEllipse, %sRadius", 
			baseName.c_str(), baseName.c_str(), baseName.c_str(), baseName.c_str(), baseName.c_str(), baseName.c_str(), baseName.c_str(), baseName.c_str());
		return (string)tmp;
	}
	std::string ToCsvString()
	{
		char tmp[4096];
		sprintf(tmp, "%f, %f, %f, %f, %f, %f, %f, %f", m00, xCentroid, yCentroid, mTheta, mElong, wEllipse, lEllipse, mRadius);
		return (string)tmp;
	}
};
#include <vector>
#include <queue>
#include <emmintrin.h>
#define SortSwap(a,b,c) {c=b;b=max(a,b);a=min(a,c);}
class Img8Util
{
public:
	static void CheckImageArgsNotInPlace(Img& src, Img& dst)
	{
		if (src.getWidth() != dst.getWidth())
		{
			ErrorExit("CheckImageArgsWidth: Different width.");
		}
		if (src.getHeight() != dst.getHeight())
		{
			ErrorExit("CheckImageArgsHeight: Different height.");
		}
		if (src.getDepthInBits() != dst.getDepthInBits())
		{
			ErrorExit("CheckImageArgsDepthInBits: Different.");
		}
		if (src.getData() == dst.getData())
		{
			ErrorExit("CheckImageArgsPitch: Not supposed to be same buffer.");
		}
	}
	static void CheckImageArgsPitch(Img& src, Img& dst)
	{
		if (src.getPitch() != dst.getPitch())
		{
			ErrorExit("CheckImageArgsPitch: Different pitch.");
		}
	}
	static void CheckImageArgsDepth8(Img& src, Img& dst)
	{
		if (src.getDepthInBits() != 8)
		{
			ErrorExit("CheckImageArgsPitch: Src not depth 8.");
		}
		if (dst.getDepthInBits() != 8)
		{
			ErrorExit("CheckImageArgsPitch: Dst not depth 8.");
		}
	}
	static void Clear(Img& src, byte value = 0)
	{
		int h = src.getHeight();
		int p = src.getPitch();
		memset(src.getData(), value, h * p);
	}
	static void ClearRoi(Img& src, int x0, int y0, int x1, int y1, byte value = 0)
	{
		IRoi roi(x0, y0, x1, y1);
		SetRoi(src, roi, value);
	}
	static void SetRoi(Img& src, IRoi& roi, byte value)
	{
		for (int y = roi.y0; y < roi.y1; y++)
		{
			for (int x = roi.x0; x < roi.x1; x++)
			{
				src.setPixel(x, y, value);
			}
		}
	}
	static void SetTile(Img& src, IRoi roi, Img& tile)
	{
		tile.setByPtr(roi.width(), roi.height(), src.getDepthInBits(), src.getPitch(), src.getData(roi.x0, roi.y0));
	}
	static void Copy(Img& src, Img& dst)
	{
		int w = src.getWidth();
		int h = src.getHeight();
		int spp = src.getPitch();
		int dpp = dst.getPitch();
		uint8_t* ps = (uint8_t*)src.getData();
		uint8_t* pd = (uint8_t*)dst.getData();
		for (int y = 0; y < h; y++)
		{
			memcpy(pd, ps, w);
			ps += spp;
			pd += dpp;
		}
	}
	static void CopyRoi(Img& src, IRoi srcRoi, Img& dst)
	{
		if (dst.getWidth() < srcRoi.getWidth()) ErrorExit("Dest image is too small for roi.");
		if (dst.getHeight() < srcRoi.getHeight()) ErrorExit("Dest image is too small for roi.");
		IRoi dstRoi(0, 0, srcRoi.getWidth(), srcRoi.getHeight());
		Clear(dst);
		dst.SetRoiAllowOffSrc(dstRoi, src, srcRoi);
	}
	static void And(Img& src1, Img& src2, Img& dst)
	{
		assert("Does not handle pitch bytes.");
		int h = src1.getHeight();
		int p = src1.getPitch();
		byte* ps1 = src1.getData();
		byte* ps2 = src2.getData();
		byte* pd = dst.getData();
		for (int i = 0; i < h*p; i++) pd[i] = ps1[i] & ps2[i];
	}
	static void Stats(Img& src, byte& min, byte& max, int& sum, int& countNonZero)
	{
		int w = src.getWidth();
		int h = src.getHeight();
		int p = src.getPitch();
		byte* ps = src.getData();
		min = 255;
		max = 0;
		sum = countNonZero = 0;
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				byte v = ps[x];
				max = MAX(max, v);
				min = MIN(min, v);
				sum += v;
				countNonZero += (v > 0);
			}
			ps += p;
		}
	}
	static void Max(Img& src1, Img& src2, Img& dst)
	{
		int h = src1.getHeight();
		int p = src1.getPitch();
		byte* ps1 = src1.getData();
		byte* ps2 = src2.getData();
		byte* pd = dst.getData();
		assert(p == src2.getPitch());
		assert(p == dst.getPitch());
		for (int i = 0; i < h*p; i++) pd[i] = MAX(ps1[i], ps2[i]);
	}
	static void Set(Img& src, byte value)
	{
		assert("Does not handle pitch bytes.");
		int h = src.getHeight();
		int p = src.getPitch();
		byte* ps = src.getData();
		for (int i = 0; i < h*p; i++) ps[i] = value;
	}
	static void Not(Img& src, Img& dst)
	{
		CheckImageArgsDepth8(src, dst);
		int w = src.getWidth();
		int h = src.getHeight();
		int spp = src.getPitch();
		int dpp = dst.getPitch();
		uint8_t* ps = (uint8_t*)src.getData();
		uint8_t* pd = (uint8_t*)dst.getData();
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				pd[x] = ~ps[x];
			}
			ps += spp;
			pd += dpp;
		}
	}
	static void Multiply(Img& src, float v, Img& dst)
	{
		CheckImageArgsDepth8(src, dst);
		int w = src.getWidth();
		int h = src.getHeight();
		int spp = src.getPitch();
		int dpp = dst.getPitch();
		uint8_t* ps = (uint8_t*)src.getData();
		uint8_t* pd = (uint8_t*)dst.getData();
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				pd[x] = (byte)(MIN(255, ps[x] * v + 0.5f));
			}
			ps += spp;
			pd += dpp;
		}
	}
	static void MinAccum(Img& src1, Img& dst)
	{
		assert("Does not handle pitch bytes.");
		int h = src1.getHeight();
		int p = src1.getPitch();
		byte* ps1 = src1.getData();
		byte* pd = dst.getData();
		for (int i = 0; i < h*p; i++) pd[i] = MIN(pd[i], ps1[i]);
	}
	static int Sum(Img& src)
	{
		assert("Does not handle pitch bytes.");
		int h = src.getHeight();
		int p = src.getPitch();
		int sum = 0;
		byte* ps = src.getData();
		for (int i = 0; i < h*p; i++) sum += ps[i];
		return sum;
	}
	static int SumRoi(Img& src, IRoi roi)
	{
		int w = roi.getHeight();
		int h = roi.getHeight();
		int p = src.getPitch();
		int sum = 0;
		byte* ps = src.getData(roi.x0, roi.y0);
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				sum += ps[x];
			}
			ps += p;
		}
		return sum;
	}
	static void AndNot(Img& src1, Img& src2, Img& dst)
	{
		int w = src1.getWidth();
		int h = src1.getHeight();
		int sp1 = src1.getPitch();
		int sp2 = src2.getPitch();
		int dp = dst.getPitch();
		byte* ps1 = src1.getData();
		byte* ps2 = src2.getData();
		byte* pd = dst.getData();
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++) pd[x] = ps1[x] & ~ps2[x];
			ps1 += sp1;
			ps2 += sp2;
			pd += dp;
		}
	}
	static int CountNonZero(Img& src)
	{
		IRoi roi(0, 0, src.getWidth(), src.getHeight());
		return CountNonZero(src, roi);
	}
	static int CountNonZero(Img& src, IRoi& roi)
	{
		int n = 0;
		for (int y = roi.y0; y < roi.y1; y++)
		{
			for (int x = roi.x0; x < roi.x1; x++)
			{
				byte v1 = src.getPixel(x, y);
				if (v1 > 0) n++;
			}
		}
		return n;
	}
	static void BitShift(Img& src, int shift, Img& dst)
	{
		CheckImageArgsDepth8(src, dst);
		int w = src.getWidth();
		int h = src.getHeight();
		int spp = src.getPitch();
		int dpp = dst.getPitch();
		uint8_t* ps = (uint8_t*)src.getData();
		uint8_t* pd = (uint8_t*)dst.getData();
		if (shift > 0)
		{
			for (int y = 0; y < h; y++)
			{
				for (int x = 0; x < w; x++)
				{
					pd[x] = (ps[x] << shift);
				}
				ps += spp;
				pd += dpp;
			}
		}
		else if (shift < 0)
		{
			shift = -shift;
			for (int y = 0; y < h; y++)
			{
				for (int x = 0; x < w; x++)
				{
					pd[x] = (ps[x] >> shift);
				}
				ps += spp;
				pd += dpp;
			}
		}
	}
	static void BitShiftMasked(Img& src1, Img& src2, int shift, Img& dst)
	{
		ErrorExit("Not upated for size.");
		byte* ps1 = src1.getData();
		byte* ps2 = src2.getData();
		byte* pd = dst.getData();
		if (shift > 0)
		{
			for (int i = 0; i < ILEN; i++) 
			{
				if (ps2[i]) pd[i] = (ps1[i] << shift);
				else pd[i] = ps1[i];
			}
		}
		else if (shift < 0)
		{
			shift = -shift;
			for (int i = 0; i < ILEN; i++) 
			{
				if (ps2[i]) pd[i] = (ps1[i] >> shift);
				else pd[i] = ps1[i];
			}
		}
	}
	static void Flood(Img& src, bool is8Connected, int cx, int cy, byte val, Img& dst, bool clearDest = true)
	{
		int sp = src.getPitch();
		int offsets4[4] = { -1, -sp, sp, 1 };
		int offsets8[8] = { -sp-1, -sp, -sp+1, -1, 1, sp-1,sp,sp+1};
		int nn = 4;
		int* offsets = offsets4;
		if (is8Connected)
		{
			nn = 8;
			offsets = offsets8;
		}
		if (clearDest)
		{
			Img8Util::Clear(dst);
		}
		if (src.getPixel(cx, cy) > 0)
		{
			queue<int> qi;
			qi.push(cx + cy * sp);
			dst.setPixel(cx, cy, val);
			while (qi.size() > 0)
			{
				int ci = qi.front();
				qi.pop();
				for (int i = 0; i < nn; i++)
				{
					int idx = ci + offsets[i];
					if (src.checkIsInImage(idx) && (src.getPixel(idx) > 0) && (dst.getPixel(idx) == 0))
					{
						dst.setPixel(idx, val);
						qi.push(idx);
					}
				}
			}
		}
	}
	static bool FindCenterMax(Img& src, int radius, int& cx, int& cy)
	{
		byte max = 0;
		cx = CX;
		cy = CY;
		bool ret = false;
		for (int y = CY - radius; y <= CY + radius; y++)
		{
			for (int x = CX - radius; x <= CX + radius; x++)
			{
				byte v = src.getPixel(x, y);
				if (v > max)
				{
					max = v;
					cx = x;
					cy = y;
					ret = true;
				}
			}
		}
		return ret;
	}
	static void FloodDescend(Img& src, bool is8Connected, int cx, int cy, byte val, int descentTolerancePower, Img& dst)
	{
		assert(false && "Not updated for image sizes.");
		int offsets4[4] = {-1, -IDIM, IDIM, 1};
		int offsets8[8] = { -IDIM-1, -IDIM, -IDIM+1, -1, 1, IDIM-1,IDIM,IDIM+1};
		const int descentToleranceMinDist = 2; 
		int nn = 4;
		int* offsets = offsets4;
		if (is8Connected)
		{
			nn = 8;
			offsets = offsets8;
		}
		Img8Util::Clear(dst);
		if (src.getPixel(cx, cy) > 0)
		{
			queue<int> qi;
			qi.push(cx + cy * IDIM);
			dst.setPixel(cx, cy, val);
			while (qi.size() > 0)
			{
				int ci = qi.front();
				qi.pop();
				byte vi = src.getPixel(ci);
				int x = ci % IDIM;
				int y = ci / IDIM;
				int cdist = roundf(pyth((float)x - cx, (float)y - cy));
				int dtpCdist = 0;
				if (cdist > descentToleranceMinDist) dtpCdist = (cdist - descentToleranceMinDist);
				int toleranceShift = dtpCdist + descentTolerancePower;
				for (int i = 0; i < nn; i++)
				{
					int idx = ci + offsets[i];
					byte nv = src.getPixel(idx);
					if ((nv > 0) && (dst.getPixel(idx) == 0))
					{
						if (nv <= (vi + (vi >> toleranceShift)))
						{
							dst.setPixel(idx, val);
							qi.push(idx);
						}
					}
				}
			}
		}
	}
	static void CopyBorders(Img& src, int borderWidth, Img& dst)
	{
		int w = src.getWidth();
		int h = src.getHeight();
		int spp = src.getPitch();
		int dpp = dst.getPitch();
		CheckImageArgsDepth8(src, dst);
		byte* ps = (byte*)src.getData();
		byte* pd = (byte*)dst.getData();
		for (int y = 0; y < h; y++)
		{
			for (int i = 0; i < borderWidth; i++)
			{
				pd[i] = ps[i];
			}
			for (int i = w - borderWidth; i < w; i++)
			{
				pd[i] = ps[i];
			}
			ps += spp;
			pd += dpp;
		}
		ps = (byte*)src.getData();
		pd = (byte*)dst.getData();
		for (int y = 0; y < borderWidth; y++)
		{
			memcpy(&pd[y*dpp], &ps[y*spp], w);
		}
		for (int y = h - borderWidth; y < h; y++)
		{
			memcpy(&pd[y*dpp], &ps[y*spp], w);
		}
	}
	static void Mean3x3(Img& src, Img& dst, bool includeBorders = true)
	{
		CheckImageArgsDepth8(src, dst);
		int w = src.getWidth();
		int h = src.getHeight();
		int spp = src.getPitch();
		int dpp = dst.getPitch();
		byte* ps = (byte*)src.getData(0, 1);
		byte* pd = (byte*)dst.getData(0, 1);
		const int off[9] = {-spp - 1, -spp, -spp + 1, -1, 0, 1, spp - 1, spp, spp + 1};
		for (int y = 1; y < h - 1; y++)
		{
			for (int x = 1; x < w - 1; x++)
			{
				int sum = 0;
				for (int j = 0; j < 9; j++)
				{
					byte v = ps[x + off[j]];
					sum += v;
				}
				pd[x] = (byte)(sum / 9);
			}
			ps += spp;
			pd += dpp;
		}
		if (!includeBorders)
		{
			CopyBorders(src, 1, dst);
		}
		else
		{
			ps = (byte*)src.getData();
			pd = (byte*)dst.getData();
			const int boff[6] = {-1, 0, 1, spp - 1, spp, spp + 1};
			const int bloff[4] = {0, 1, spp, spp + 1};
			int x = 0;
			int sum = 0;
			for (int j = 0; j < 4; j++)
			{
				byte v = ps[x + bloff[j]];
				sum += v;
			}
			pd[x] = (byte)(sum / 4);
			for (int x = 1; x < w - 1; x++)
			{
				sum = 0;
				for (int j = 0; j < 6; j++)
				{
					byte v = ps[x + boff[j]];
					sum += v;
				}
				pd[x] = (byte)(sum / 6);
			}
			const int broff[4] = {-1, 0, spp - 1, spp};
			x = w - 1;
			sum = 0;
			for (int j = 0; j < 4; j++)
			{
				byte v = ps[x + broff[j]];
				sum += v;
			}
			pd[x] = (byte)(sum / 4);
			ps = (byte*)src.getData(0, h - 1);
			pd = (byte*)dst.getData(0, h - 1);
			const int toff[9] = {-spp - 1, -spp, -spp + 1, -1, 0, 1};
			const int tloff[4] = {0, 1, -spp, -spp + 1};
			x = 0;
			sum = 0;
			for (int j = 0; j < 4; j++)
			{
				byte v = ps[x + tloff[j]];
				sum += v;
			}
			pd[x] = (byte)(sum / 4);
			for (int x = 1; x < w - 1; x++)
			{
				sum = 0;
				for (int j = 0; j < 6; j++)
				{
					byte v = ps[x + toff[j]];
					sum += v;
				}
				pd[x] = (byte)(sum / 6);
			}
			const int troff[4] = {-1, 0, -spp - 1, -spp};
			x = w - 1;
			sum = 0;
			for (int j = 0; j < 4; j++)
			{
				byte v = ps[x + troff[j]];
				sum += v;
			}
			pd[x] = (byte)(sum / 4);
			ps = (byte*)src.getData(0, 1);
			pd = (byte*)dst.getData(0, 1);
			const int loff[6] = {-spp, -spp + 1, 0, 1, spp, spp + 1};
			for (int y = 1; y < h - 1; y++)
			{
				sum = 0;
				for (int j = 0; j < 6; j++)
				{
					byte v = ps[loff[j]];
					sum += v;
				}
				pd[0] = (byte)(sum / 6);
				ps += spp;
				pd += dpp;
			}
			ps = (byte*)src.getData(w - 1, 1);
			pd = (byte*)dst.getData(w - 1, 1);
			const int roff[6] = {-spp - 1, -spp, -1, 0, spp - 1, spp};
			for (int y = 1; y < h - 1; y++)
			{
				sum = 0;
				for (int j = 0; j < 6; j++)
				{
					byte v = ps[roff[j]];
					sum += v;
				}
				pd[0] = (byte)(sum / 6);
				ps += spp;
				pd += dpp;
			}
		}
	}
	static void Mean3x3(Img& src, bool includeBorders = true)
	{
		Img tmp(src.getWidth(), src.getHeight(), 8);
		Mean3x3(src, tmp, includeBorders);
		Copy(tmp, src);
	}
	static void Mean3x3GoodValues(Img& src, Img& dst, bool includeBorders = true)
	{
		CheckImageArgsDepth8(src, dst);
		int w = src.getWidth();
		int h = src.getHeight();
		int spp = src.getPitch();
		int dpp = dst.getPitch();
		byte* ps = (byte*)src.getData(0, 1);
		byte* pd = (byte*)dst.getData(0, 1);
		const int off[9] = {-spp - 1, -spp, -spp + 1, -1, 0, 1, spp - 1, spp, spp + 1};
		for (int y = 1; y < h - 1; y++)
		{
			for (int x = 1; x < w - 1; x++)
			{
				int n = 0;
				int sum = 0;
				for (int j = 0; j < 9; j++)
				{
					byte v = ps[x + off[j]];
					if (v > 0)
					{
						n++;
						sum += v;
					}
				}
				if (n > 0)
				{
					pd[x] = (byte)(sum / n);
				}
			}
			ps += spp;
			pd += dpp;
		}
		if (!includeBorders)
		{
			CopyBorders(src, 1, dst);
		}
		else
		{
			ps = (byte*)src.getData();
			pd = (byte*)dst.getData();
			const int boff[6] = {-1, 0, 1, spp - 1, spp, spp + 1};
			if (ps[0] > 0)
			{
				const int bloff[4] = {0, 1, spp, spp + 1};
				int x = 0;
				int n = 0;
				int sum = 0;
				for (int j = 0; j < 4; j++)
				{
					byte v = ps[x + bloff[j]];
					if (v > 0) { n++; sum += v; }
				}
				if (n > 0)
				{
					pd[x] = (byte)(sum / n);
				}
				else
				{
					pd[x] = 0;
				}
			}
			else
			{
				pd[0] = 0;
			}
			for (int x = 1; x < w - 1; x++)
			{
				int n = 0;
				int sum = 0;
				if (ps[x] > 0)
				{
					for (int j = 0; j < 6; j++)
					{
						byte v = ps[x + boff[j]];
						if (v > 0) { n++; sum += v; }
					}
					if (n > 0)
					{
						pd[x] = (byte)(sum / n);
					}
					else
					{
						pd[x] = 0;
					}
				}
				else
				{
					pd[x] = 0;
				}
			}
			if (ps[w - 1] > 0)
			{
				const int broff[4] = {-1, 0, spp - 1, spp};
				int x = w - 1;
				int n = 0;
				int sum = 0;
				for (int j = 0; j < 4; j++)
				{
					byte v = ps[x + broff[j]];
					if (v > 0) { n++; sum += v; }
				}
				if (n > 0)
				{
					pd[x] = (byte)(sum / n);
				}
				else
				{
					pd[x] = 0;
				}
			}
			else
			{
				pd[w - 1] = 0;
			}
			ps = (byte*)src.getData(0, h - 1);
			pd = (byte*)dst.getData(0, h - 1);
			const int toff[9] = {-spp - 1, -spp, -spp + 1, -1, 0, 1};
			if (ps[0] > 0)
			{
				const int tloff[4] = {0, 1, -spp, -spp + 1};
				int x = 0;
				int n = 0;
				int sum = 0;
				for (int j = 0; j < 4; j++)
				{
					byte v = ps[x + tloff[j]];
					if (v > 0) { n++; sum += v; }
				}
				if (n > 0)
				{
					pd[x] = (byte)(sum / n);
				}
				else
				{
					pd[x] = 0;
				}
			}
			else
			{
				pd[0] = 0;
			}
			for (int x = 1; x < w - 1; x++)
			{
				int n = 0;
				int sum = 0;
				if (ps[x] > 0)
				{
					for (int j = 0; j < 6; j++)
					{
						byte v = ps[x + toff[j]];
						if (v > 0) { n++; sum += v; }
					}
					if (n > 0)
					{
						pd[x] = (byte)(sum / n);
					}
					else
					{
						pd[x] = 0;
					}
				}
				else
				{
					pd[x] = 0;
				}
			}
			if (ps[w - 1] > 0)
			{
				const int troff[4] = {-1, 0, -spp - 1, -spp};
				int x = w - 1;
				int n = 0;
				int sum = 0;
				for (int j = 0; j < 4; j++)
				{
					byte v = ps[x + troff[j]];
					if (v > 0) { n++; sum += v; }
				}
				if (n > 0)
				{
					pd[x] = (byte)(sum / n);
				}
				else
				{
					pd[x] = 0;
				}
			}
			else
			{
				pd[w - 1] = 0;
			}
			ps = (byte*)src.getData(0, 1);
			pd = (byte*)dst.getData(0, 1);
			const int loff[6] = {-spp, -spp + 1, 0, 1, spp, spp + 1};
			for (int y = 1; y < h - 1; y++)
			{
				int n = 0;
				int sum = 0;
				if (ps[0] > 0)
				{
					for (int j = 0; j < 6; j++)
					{
						byte v = ps[loff[j]];
						if (v > 0) { n++; sum += v; }
					}
					if (n > 0)
					{
						pd[0] = (byte)(sum / n);
					}
					else
					{
						pd[0] = 0;
					}
				}
				else
				{
					pd[0] = 0;
				}
				ps += spp;
				pd += dpp;
			}
			ps = (byte*)src.getData(w - 1, 1);
			pd = (byte*)dst.getData(w - 1, 1);
			const int roff[6] = {-spp - 1, -spp, -1, 0, spp - 1, spp};
			for (int y = 1; y < h - 1; y++)
			{
				int n = 0;
				int sum = 0;
				if (ps[0] > 0)
				{
					for (int j = 0; j < 6; j++)
					{
						byte v = ps[roff[j]];
						if (v > 0) { n++; sum += v; }
					}
					if (n > 0)
					{
						pd[0] = (byte)(sum / n);
					}
					else
					{
						pd[0] = 0;
					}
				}
				else
				{
					pd[0] = 0;
				}
				ps += spp;
				pd += dpp;
			}
		}
	}
	static void Binomial3x3(Img& image, Img& dst)
	{
		ErrorExit("Not upated for size.");
		assert(image.getData() != dst.getData());
		byte* ps = (byte*)image.getData();
		byte* pd = (byte*)dst.getData();
		const int off[8] = { -IDIM-1, -IDIM, -IDIM+1, -1, 1, IDIM-1,IDIM,IDIM+1};
		memset(pd, 0, ILEN);
		for (int y = 1; y < IDIM - 1; y++)
		{
			int si = y * IDIM;
			for (int x = 1; x < IDIM - 1; x++)
			{
				int n = 2;
				int sum = ps[si + x] * 2;
				for (int j = 0; j < 8; j++)
				{
					byte v = ps[si + x + off[j]];
					n++; 
					sum += v;
				}
				if (n > 0)
				{
					pd[si + x] = (byte)(sum / n);
				}
			}
		}
	}
	static void LocalMax(Img& image, Img& dst)
	{
		ErrorExit("Not upated for size.");
		assert(image.getData() != dst.getData());
		byte* ps = (byte*)image.getData();
		byte* pd = (byte*)dst.getData();
		const int off[8] = { -IDIM-1, -IDIM, -IDIM+1, -1, 1, IDIM-1,IDIM,IDIM+1};
		memset(pd, 0, ILEN);
		for (int y = 1; y < IDIM - 1; y++)
		{
			int si = y * IDIM;
			for (int x = 1; x < IDIM - 1; x++)
			{
				byte vc = ps[si + x];
				byte max = vc;
				for (int j = 0; j < 8; j++)
				{
					byte v = ps[si + x + off[j]];
					max = MAX(max, v);
				}
				if ((max > 0) && (vc >= max))
				{
					pd[si + x] = (byte)255;
				}
			}
		}
	}
	static void Label(Img& image, Img& dst, bool is8Connected)
	{
		CheckImageArgsDepth8(image, dst);
		CheckImageArgsNotInPlace(image, dst);
		int w = image.getWidth();
		int h = image.getHeight();
		int sp = image.getPitch();
		byte* ps = (byte*)image.getData();
		byte* pd = (byte*)dst.getData();
		Clear(dst);
		int label = 1;
		for (int y = 1; y < h - 1; y++)
		{
			int si = y * sp;
			for (int x = 1; x < w - 1; x++)
			{
				byte vs = ps[si + x];
				byte vd = pd[si + x];
				if ((vs > 0) && (vd == 0))
				{
					Img8Util::Flood(image, is8Connected, x, y, label++, dst, false);
					if (label > 255) ErrorExit("The input image has more than 255 blobs.");
				}
			}
		}
	}
	static void FindBlobs(Img& image, vector<Img8Blob>& blobs)
	{
		int w = image.getWidth();
		int h = image.getHeight();
		int sp = image.getPitch();
		byte* ps = (byte*)image.getData();
		vector<float> xsums(255);
		vector<float> ysums(255);
		vector<int> counts(255);
		vector<int> xs(255);
		vector<int> ys(255);
		for (int y = 1; y < h - 1; y++)
		{
			int si = y * sp;
			for (int x = 1; x < w - 1; x++)
			{
				byte vs = ps[si + x];
				if (vs > 0)
				{
					counts[vs]++;
					xsums[vs] += x;
					ysums[vs] += y;
					xs[vs] = x;
					ys[vs] = y;
				}
			}
		}
		int n = 0;
		for (int i = 1; i < 255; i++) if (counts[i] > 0) n++;
		blobs.resize(n);
		int j = 0;
		for (int i = 1; i < 255; i++)
		{
			if (counts[i] > 0)
			{
				Img8Blob blob;
				blob.area = counts[i];
				blob.label = i;
				blob.xCentroid = roundf(xsums[i] / counts[i]);
				blob.yCentroid = roundf(ysums[i] / counts[i]);
				if (image.getPixel(roundf(blob.xCentroid), roundf(blob.yCentroid)) == i)
				{
					blob.cx = roundf(blob.xCentroid);
					blob.cy = roundf(blob.yCentroid);
				}
				else
				{
					blob.cx = xs[i];
					blob.cy = ys[i];
				}
				blobs[j++] = blob;
			}
		}
	}
	static void DrawBlobCenters(vector<Img8Blob>& blobs, Img& dst)
	{
		ErrorExit("Not upated for size.");
		byte* pd = (byte*)dst.getData();
		memset(pd, 0, ILEN);
		int n = blobs.size();
		for (int i = 0; i < n; i++)
		{
			dst.setPixel(blobs[i].cx, blobs[i].cy, i);
		}
	}
	static void IncrementMasked(Img& mask, int increment, Img& dst)
	{
		ErrorExit("Not upated for size.");
		byte* ps1 = mask.getData();
		byte* pd = dst.getData();
		for (int i = 0; i < ILEN; i++) 
		{
			if (ps1[i] > 0) pd[i] += increment;
		}
	}
	static void Threshold(Img& image, byte threshold, bool isLte, Img& dst)
	{
		int h = image.getHeight();
		int p = image.getPitch();
		assert(p == dst.getPitch());
		byte* ps = (byte*)image.getData();
		byte* pd = (byte*)dst.getData();
		int np = h * p;
		if (isLte)
		{
			for (int i = 0; i < np; i++) 
			{
				pd[i] = (ps[i] <= threshold ? 255 : 0);
			}
		}
		else
		{
			for (int i = 0; i < np; i++) 
			{
				pd[i] = (ps[i] >= threshold ? 255 : 0);
			}
		}
	}
	static float NormCrossCorr(Img& src1, IRoi& roi, Img& src2)
	{
		float sum12 = 0;
		float sum1 = 0, sum2 = 0;
		float sqSum1 = 0, sqSum2 = 0;
		int n = 0;
		for (int y = roi.y0; y < roi.y1; y++)
		{
			for (int x = roi.x0; x < roi.x1; x++)
			{
				byte v1 = src1.getPixel(x, y);
				byte v2 = src2.getPixel(x, y);
				sum1 += v1;
				sum2 += v2;
				sum12 += v1 * v2;
				sqSum1 += v1 * v1;
				sqSum2 += v2 * v2;
				n++;
			}
		}
		float numerator = n * sum12 - sum1 * sum2;
		float denom = sqrt(n * sqSum1 - sum1 * sum1) * sqrt(n * sqSum2 - sum2 * sum2);
		if (abs(denom < 0.000001f))
		{
			return 0.0f;
		}
		else
		{
			return numerator / denom;
		}
	}
	static void Shift(Img& src, int dx, int dy, Img& dst)
	{
		Img8Util::Copy(src, dst);
		dst.Shift(dx,dy,0);
	}
	static void Max(unsigned char* p1, int n, unsigned char* p2, unsigned char* dst)
	{
		for (int i = 0; i < n; i++)
		{
			dst[i] = MAX(p1[i], p2[i]);
		}
	}
	static void Min(unsigned char* p1, int n, unsigned char* p2, unsigned char* dst)
	{
		for (int i = 0; i < n; i++)
		{
			dst[i] = MIN(p1[i], p2[i]);
		}
	}
	static void Dilate3x3(Img& src, Img& dst)
	{
		int w = src.getWidth();
		int h = src.getHeight();
		Img tmp(w,h,8);
		Dilate1x3(src, tmp);
		Dilate3x1(tmp, dst);
	}
	static void Erode3x3(Img& src, Img& dst)
	{
		int w = src.getWidth();
		int h = src.getHeight();
		Img tmp(w, h, 8);
		Erode1x3(src, tmp);
		Erode3x1(tmp, dst);
	}
	static void Dilate1x3(Img& src, Img& dst)
	{
		int width = src.getWidth();
		int height = src.getHeight();
		Img tmp(width, height);
		Copy(src, tmp);
		for (int y = 0; y < height - 1; y++)
		{
			Max(src.getData(0,y), width, tmp.getData(0,y+1), dst.getData(0,y));
		}
		memcpy(dst.getData(0,height-1), src.getData(0, height-1), width); 
		for (int y = 0; y < height - 1; y++)
		{
			Max(src.getData(0, y), width, dst.getData(0, y + 1), tmp.getData(0, y + 1));
		}
		for (int y = 1; y < height; y++)
		{
			memcpy(dst.getData(0, y), tmp.getData(0, y), width);
		}
	}
	static void Erode1x3(Img& src, Img& dst)
	{
		assert(src.getDepthInBits() == 8);
		assert(dst.getDepthInBits() == 8);
		int width = src.getWidth();
		int height = src.getHeight();
		Img tmp(width, height);
		Copy(src, tmp);
		for (int y = 0; y < height - 1; y++)
		{
			Min(src.getData(0, y), width, tmp.getData(0, y + 1), dst.getData(0, y));
		}
		memcpy(dst.getData(0, height - 1), src.getData(0, height - 1), width); 
		for (int y = 0; y < height - 1; y++)
		{
			Min(src.getData(0, y), width, dst.getData(0, y + 1), tmp.getData(0, y + 1));
		}
		for (int y = 1; y < height; y++)
		{
			memcpy(dst.getData(0, y), tmp.getData(0, y), width);
		}
	}
	static void Dilate3x1(Img& src, Img& dst)
	{
		int width = src.getWidth();
		int height = src.getHeight();
		int srcPitch = src.getPitch();
		int dstPitch = dst.getPitch();
		unsigned char* ps = src.getData();
		unsigned char* pd = dst.getData();
		for (int y = 0; y < height; y++)
		{
			pd[0] = MAX(ps[0], ps[1]);
			for (int x = 1; x < width - 1; x++)
			{
				pd[x] = MAX(MAX(ps[x - 1], ps[x]), ps[x + 1]);
			}
			pd[width - 1] = MAX(ps[width - 1], ps[width - 2]);
			ps += srcPitch;
			pd += dstPitch;
		}
	}
	static void Erode3x1(Img& src, Img& dst)
	{
		int width = src.getWidth();
		int height = src.getHeight();
		int srcPitch = src.getPitch();
		int dstPitch = dst.getPitch();
		unsigned char* ps = src.getData();
		unsigned char* pd = dst.getData();
		for (int y = 0; y < height; y++)
		{
			pd[0] = MIN(ps[0], ps[1]);
			for (int x = 1; x < width - 1; x++)
			{
				pd[x] = MIN(MIN(ps[x - 1], ps[x]), ps[x + 1]);
			}
			pd[width - 1] = MIN(ps[width - 1], ps[width - 2]);
			ps += srcPitch;
			pd += dstPitch;
		}
	}
	static void Median4Masked(Img& i1, Img& i2, Img& i3, Img& i4, Img& m1, Img& m2, Img& m3, Img& m4, Img& dst, Img& dstMask)
	{
		int nGood = 2;
		int threshold = 255 * nGood;
		CheckImageArgsNotInPlace(i1, i2);
		CheckImageArgsNotInPlace(i1, i3);
		CheckImageArgsNotInPlace(i1, i4);
		CheckImageArgsNotInPlace(i1, dst);
		CheckImageArgsNotInPlace(i1, dstMask);
		CheckImageArgsDepth8(i1, i2);
		CheckImageArgsDepth8(i3, i4);
		CheckImageArgsPitch(i1, i2);
		CheckImageArgsPitch(i1, i3);
		CheckImageArgsPitch(i1, i4);
		CheckImageArgsPitch(m1, m2);
		CheckImageArgsPitch(m1, m3);
		CheckImageArgsPitch(m1, m4);
		int w = i1.getWidth();
		int h = i1.getHeight();
		int sp = i1.getPitch(); 
		int smp = m1.getPitch(); 
		int dp = dst.getPitch();
		int dmp = dstMask.getPitch();
		unsigned char* ps1 = i1.getData();
		unsigned char* ps2 = i2.getData();
		unsigned char* ps3 = i3.getData();
		unsigned char* ps4 = i4.getData();
		unsigned char* pm1 = m1.getData();
		unsigned char* pm2 = m2.getData();
		unsigned char* pm3 = m3.getData();
		unsigned char* pm4 = m4.getData();
		unsigned char* pd = dst.getData();
		unsigned char* pmd = dstMask.getData();
		byte tmp;
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				byte v1 = ps1[x] & ~pm1[x];
				byte v2 = ps2[x] & ~pm2[x];
				byte v3 = ps3[x] & ~pm3[x];
				byte v4 = ps4[x] & ~pm4[x];
				SortSwap(v1, v2, tmp);
				SortSwap(v1, v3, tmp);
				SortSwap(v2, v3, tmp);
				SortSwap(v1, v4, tmp);
				SortSwap(v2, v4, tmp);
				int maskCount = pm1[x] + pm2[x] + pm3[x] + pm4[x];
				pmd[x] = (byte)(maskCount > threshold ? 0xFF : 0); 
				if (maskCount > threshold)
				{
					pd[x] = 0;
				}
				else if (maskCount == threshold)
				{
					SortSwap(v3, v4, tmp);
					pd[x] = v3;
				}
				else
				{
					pd[x] = v2;
				}
			}
			ps1 += sp;
			ps2 += sp;
			ps3 += sp;
			ps4 += sp;
			pm1 += smp;
			pm2 += smp;
			pm3 += smp;
			pm4 += smp;
			pd += dp;
			pmd += dmp;
		}
	}
	static void Downsample(Img& src, Img& dst, bool doMean2x2Filter = false)
	{
		int width = src.getWidth();
		int height = src.getHeight();
		assert(width / 2 == dst.getWidth());
		assert(height / 2 == dst.getHeight());
		assert(src.getDepthInBits() == dst.getDepthInBits());
		assert(src.getDepthInBits() == 8);
		assert(src.getData() != dst.getData());
		int snp = src.getPitch();
		int dnp = dst.getPitch();
		byte* __restrict ps = src.getData();
		byte* __restrict pd = dst.getData();
		if (!doMean2x2Filter)
		{
			for (int y = 0; y < height; y += 2)
			{
				int idx = 0;
				for (int i = 0; i < width; i += 2)
				{
					pd[idx++] = ps[i];
				}
				ps += snp;
				pd += dnp;
			}
		}
		else
		{
			for (int y = 0; y < height; y += 2)
			{
				int idx = 0;
				for (int i = 0; i < width; i += 2)
				{
					int val = ps[i] + ps[i + 1] + ps[i + snp] + ps[i + snp + 1];
					pd[idx++] = val >> 2;
				}
				ps += snp * 2;
				pd += dnp;
			}
		}
	}
	static void Subtract(Img& src1, Img& src2, Img& dst)
	{
		int w = src1.getWidth();
		int h = src1.getHeight();
		int sp1 = src1.getPitch();
		int sp2 = src2.getPitch();
		int dp = dst.getPitch();
		byte* ps1 = src1.getData();
		byte* ps2 = src2.getData();
		byte* pd = dst.getData();
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				pd[x] = MAX(0, ps1[x] - ps2[x]);
			}
			ps1 += sp1;
			ps2 += sp2;
			pd += dp;
		}
	}
	static void Add(Img& src1, Img& src2, Img& dst)
	{
		int w = src1.getWidth();
		int h = src1.getHeight();
		int sp1 = src1.getPitch();
		int sp2 = src2.getPitch();
		int dp = dst.getPitch();
		byte* ps1 = src1.getData();
		byte* ps2 = src2.getData();
		byte* pd = dst.getData();
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				pd[x] = MIN(255, ps1[x] + ps2[x]);
			}
			ps1 += sp1;
			ps2 += sp2;
			pd += dp;
		}
	}
	static void AddNonZero(Img& src1, Img& src2, Img& dst)
	{
		int w = src1.getWidth();
		int h = src1.getHeight();
		int sp1 = src1.getPitch();
		int sp2 = src2.getPitch();
		int dp = dst.getPitch();
		byte* ps1 = src1.getData();
		byte* ps2 = src2.getData();
		byte* pd = dst.getData();
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				if (ps1[x] > 0)
				{
					pd[x] = MIN(255, ps1[x] + ps2[x]);
				}
				else
				{
					pd[x] = 0;
				}
			}
			ps1 += sp1;
			ps2 += sp2;
			pd += dp;
		}
	}
	static void ZeroCol(Img& src, int colIndex)
	{
		int h = src.getHeight();
		int spp = src.getPitch();
		byte* ps = (byte*)src.getData();
		for (int y = 0; y < h; y++)
		{
			ps[y * spp + colIndex] = 0;
		}
	}
	static void ZeroRow(Img& src, int rowIndex)
	{
		int sp = src.getPitch();
		byte* ps = (byte*)src.getData(0, rowIndex);
		memset(ps, 0, sp);
	}
	static void CountNeighbors8(Img& src, Img& dst, int nNonZeroThreshold, int meanThreshold, int minNeighborMax)
	{
		CheckImageArgsDepth8(src, dst);
		int w = src.getWidth();
		int h = src.getHeight();
		int srcPitch = src.getPitch();
		int dstPitch = dst.getPitch();
		dst.setPixel(0, 0, 0);
		dst.setPixel(0, h-1, 0);
		dst.setPixel(w-1, 0, 0);
		dst.setPixel(w-1, h-1, 0);
		Img8Util::ZeroRow(dst, 0);
		Img8Util::ZeroRow(dst, h - 1);
		Img8Util::ZeroCol(dst, 0);
		Img8Util::ZeroCol(dst, w - 1);
		unsigned char* ps = src.getData(0, 1); 
		unsigned char* pd = dst.getData(0, 1);
		for (int y = 1; y < h - 1; y++)
		{
			for (int x = 1; x < w - 1; x++)
			{
				byte count = 1;
				count += (ps[x - 1 - srcPitch] > 0);
				count += (ps[x - srcPitch] > 0);
				count += (ps[x + 1 - srcPitch] > 0);
				count += (ps[x - 1] > 0);
				count += (ps[x + 1] > 0);
				count += (ps[x - 1 + srcPitch] > 0);
				count += (ps[x + srcPitch] > 0);
				count += (ps[x + 1 + srcPitch] > 0);
				int sum = ps[x];
				sum += (ps[x - 1 - srcPitch]);
				sum += (ps[x - srcPitch]);
				sum += (ps[x + 1 - srcPitch]);
				sum += (ps[x - 1]);
				sum += (ps[x + 1]);
				sum += (ps[x - 1 + srcPitch]);
				sum += (ps[x + srcPitch]);
				sum += (ps[x + 1 + srcPitch]);
				int max = ps[x];
				max = MAX(max, ps[x - 1 - srcPitch]);
				max = MAX(max, ps[x - srcPitch]);
				max = MAX(max, ps[x + 1 - srcPitch]);
				max = MAX(max, ps[x - 1]);
				max = MAX(max, ps[x + 1]);
				max = MAX(max, ps[x - 1 + srcPitch]);
				max = MAX(max, ps[x + srcPitch]);
				max = MAX(max, ps[x + 1 + srcPitch]);
				byte mean = (byte)(sum / count);
				pd[x] = (ps[x] > 0 ? 255 : 0) 
					& (max >= minNeighborMax ? 255 : 0)
					& (count >= nNonZeroThreshold ? 255 : 0)
					& (mean >= meanThreshold ? 255 : 0);
			}
			ps += srcPitch;
			pd += dstPitch;
		}
	}
	static void CountNeighbors4(Img& src, Img& dst, int nNonZeroThreshold)
	{
		CheckImageArgsDepth8(src, dst);
		int width = src.getWidth();
		int height = src.getHeight();
		int srcPitch = src.getPitch();
		int dstPitch = dst.getPitch();
		unsigned char* ps = src.getData();
		unsigned char* pd = dst.getData();
		ps += srcPitch;
		pd += dstPitch;
		for (int y = 1; y < height - 1; y++)
		{
			for (int x = 1; x < width - 1; x++)
			{
				byte sum = 1 + 
					(ps[x - srcPitch] > 0) +
					(ps[x - 1] > 0) +
					(ps[x + 1] > 0) +
					(ps[x + srcPitch] > 0);
				pd[x] = (ps[x] > 0 ? 255 : 0) & (sum >= nNonZeroThreshold ? 255 : 0);
			}
			ps += srcPitch;
			pd += dstPitch;
		}
	}
	static void DropSmallBlobs(Img& src, Img& dst, int minNeighborCount, int minPixelMean, int minNeighborMax)
	{
		CheckImageArgsDepth8(src, dst);
		int w = src.getWidth();
		int h = src.getHeight();
		Img i1(w, h, 8);
		Img i2(w, h, 8);
		Img8Util::CountNeighbors8(src, i2, minNeighborCount, minPixelMean, minNeighborMax);
		Img8Util::Dilate3x3(i2, i1); 
		Img8Util::And(src, i1, dst);
	}
	static void ApplyMask(Img& src, Img& mask, bool isBadMask = true)
	{
		CheckImageArgsDepth8(src, mask);
		int w = src.getWidth();
		int h = src.getHeight();
		int spp = src.getPitch();
		int mpp = mask.getPitch();
		uint8_t* ps = (uint8_t*)src.getData();
		uint8_t* pm = (uint8_t*)mask.getData();
		if (isBadMask)
		{
			for (int y = 0; y < h; y++)
			{
				for (int x = 0; x < w; x++)
				{
					if (pm[x] > 0) ps[x] = 0;
				}
				ps += spp;
				pm += mpp;
			}
		}
		else
		{
			for (int y = 0; y < h; y++)
			{
				for (int x = 0; x < w; x++)
				{
					if (pm[x] == 0) ps[x] = 0;
				}
				ps += spp;
				pm += mpp;
			}
		}
	}
	static void DrawRoi(Img& img, IRoi roi, byte value)
	{
		DrawBox(img, roi.x0, roi.y0, roi.x1, roi.y1, value);
	}
	static void DrawBox(Img& img, int x0, int y0, int x1, int y1, byte value)
	{
		int w = img.getWidth();
		int h = img.getHeight();
		int pb = img.getPitch(); 
		byte* data = (byte*)img.getData();
		x0 = CLIP(x0, 0, w - 1);
		x1 = CLIP(x1, 0, w);
		y0 = CLIP(y0, 0, h - 1);
		y1 = CLIP(y1, 0, h);
		for (int i = x0 + y0 * pb; i < x1 + y0 * pb; i++)
		{
			data[i] = value;
		}
		for (int y = y0 + 1; y < y1; y++)
		{
			data[x0 + y * pb] = value;
			data[x1 - 1 + y * pb] = value;
		}
		for (int i = x0 + (y1 - 1) * pb; i < x1 + (y1 - 1) * pb; i++)
		{
			data[i] = value;
		}
	}
	static void FillUnderMask(Img& src, Img& posMask, bool interpolate = false)
	{
		CheckImageArgsDepth8(src, posMask);
		int w = src.getWidth();
		int h = src.getHeight();
		int spp = src.getPitch();
		int mpp = posMask.getPitch();
		uint8_t* ps = (uint8_t*)src.getData();
		uint8_t* pm = (uint8_t*)posMask.getData();
		if (!interpolate)
		{
			for (int y = 0; y < h; y++)
			{
				int lgx = -1; 
				for (int x = 0; x < w; x++)
				{
					if (pm[x] > 0)
					{
						if (lgx > 0)
						{
							ps[x] = ps[lgx];
						}
					}
					else
					{
						if (lgx < 0)
						{
							for (int j = x - 1; j >= 0; j--)
							{
								ps[j] = ps[x];
							}
						}
						lgx = x;
					}
				}
				ps += spp;
				pm += mpp;
			}
		}
		else
		{
			for (int y = 0; y < h; y++)
			{
				int lgx = -1; 
				for (int x = 0; x < w; x++)
				{
					if (pm[x] > 0)
					{
					}
					else
					{
						if (lgx < 0)
						{
							for (int j = x - 1; j >= 0; j--)
							{
								ps[j] = ps[x];
							}
						}
						else
						{
							if ((lgx >= 0) && (x > lgx + 1))
							{
								float len = x - lgx;
								int dv = ps[x] - ps[lgx];
								for (int i = lgx + 1; i < x; i++)
								{
									float frac = (i - lgx) / len;
									uint16_t val = (uint16_t)(ps[lgx] + roundf(frac * dv));
									ps[i] = val;
								}
							}
						}
						lgx = x;
					}
				}
				if ((lgx < w - 1) && (lgx >= 0))
				{
					for (int i = lgx + 1; i < w; i++)
					{
						ps[i] = ps[lgx];
					}
				}
				ps += spp;
				pm += mpp;
			}
		}
	}
	static void ShiftSmartRowRight(byte* ps, int n, float dx, byte* pd, bool debug)
	{
		assert(dx >= 0); 
		assert(dx <= 0.5f); 
		memcpy(pd, ps, n);
		if (dx < 0.005f)
		{
			return;
		}
		int x = 1; 
		while (x < n)
		{
			while ((ps[x] == 0) && (x < n))
			{
				x++;
			}
			if (x < n)
			{
				int runStart = x;
				byte maxValue = 0;
				int peakStart = -1; 
				int peakEnd = -1; 
				while ((ps[x] != 0) && (x < n))
				{
					if (ps[x] > maxValue)
					{
						peakStart = x;
						maxValue = ps[x];
						peakEnd = x;
					}
					else if (ps[x] == maxValue)
					{
						peakEnd = x;
					}
					x++;
				}
				peakEnd++; 
				int runEnd = x; 
				if (debug) logd.printf("Run from %d to %d\n", runStart, runEnd);
				if (debug) logd.printf("Peak: %d from %d to %d\n", maxValue, peakStart, peakEnd);
				int i = runStart;
				for (; i < peakStart; i++)
				{
					pd[i] = (byte)(dx * ps[i - 1] + (1 - dx) * ps[i] + 0.5f);
				}
				for (; i < peakEnd; i++)
				{
					pd[i] = ps[i];
				}
				for (; i < runEnd; i++)
				{
					pd[i] = (byte)(dx * ps[i - 1] + (1 - dx) * ps[i] + 0.5f);
				}
			}
		}
	}
	static void Transpose(Img& src, Img& dst)
	{
		int w = src.getWidth();
		int h = src.getHeight();
		int sp = src.getPitch();
		int dp = src.getPitch();
		byte* ps = src.getData();
		byte* pd = dst.getData();
		for (int y = 0; y < h; y++)
		{
			pd = dst.getData(y,0);
			for (int x = 0; x < w; x++)
			{
				pd[x * dp] = ps[x];
			}
			ps += sp;
		}
	}
	static void ShiftInt(Img& src, pair<int, int> offset, int value)
	{
		src.Shift(offset.first, offset.second, value);
	}
	static void ShiftSmart(Img& src, float dx, float dy, Img& dst)
	{
		if ((dx <= -1.0f) || (dx >= 1.0f) || (dy <= -1.0f) || (dy >= 1.0f))
		{
			pair<int, int> offsetInt;
			offsetInt.first = (int)roundf(dx);
			offsetInt.second = (int)roundf(dy);
			Img ti1(src);
			ShiftInt(ti1, offsetInt, 0);
			ShiftSmart(ti1, dx - offsetInt.first, dy - offsetInt.second, dst);
		}
		else
		{
			int w = src.getWidth();
			int h = src.getHeight();
			int sp = src.getPitch();
			int dp = src.getPitch();
			byte* ps = src.getData();
			byte* pd = dst.getData();
			byte* buf1 = new byte[w];
			byte* buf2 = new byte[w];
			if (CLOSE_ENOUGH(dx, 0.0f, 0.005f))
			{
				Copy(src, dst);
			}
			else if (dx >= 0)
			{
				if (dx <= 0.5f)
				{
					for (int y = 0; y < h; y++)
					{
						ShiftSmartRowRight(ps, w, dx, pd, false);
						ps += sp;
						pd += dp;
					}
				}
				else 
				{
					for (int y = 0; y < h; y++)
					{
						buf1[0] = 0;
						for (int i = 0; i < w - 1; i++) buf1[i + 1] = ps[i];
						for (int i = 0; i < w; i++) buf2[i] = buf1[w - i - 1];
						ShiftSmartRowRight(buf2, w, 1.0f - dx, buf1, false);
						for (int i = 0; i < w; i++) pd[w - i - 1] = buf1[i];
						ps += sp;
						pd += dp;
					}
				}
			}
			else 
			{
				if (dx >= -0.5f)
				{
					for (int y = 0; y < h; y++)
					{
						for (int i = 0; i < w; i++) buf2[i] = ps[w - i - 1];
						ShiftSmartRowRight(buf2, w, abs(dx), buf1, false);
						for (int i = 0; i < w; i++) pd[w - i - 1] = buf1[i];
						ps += sp;
						pd += dp;
					}
				}
				else 
				{
					for (int y = 0; y < h; y++)
					{
						buf1[w-1] = 0;
						for (int i = 0; i < w - 1; i++) buf1[i] = ps[i + 1];
						ShiftSmartRowRight(buf1, w, 1.0f + dx, pd, false);
						ps += sp;
						pd += dp;
					}
				}
			}
			delete[] buf1;
			delete[] buf2;
			if (CLOSE_ENOUGH(dy, 0.0f, 0.005f))
			{
			}
			else
			{
				Img trans1(w, h, 8);
				Img trans2(w, h, 8);
				Transpose(dst, trans1);
				ShiftSmart(trans1, dy, 0.0f, trans2);
				Transpose(trans2, dst);
			}
		}
	}
private:
	static void upSampleSubSse2(Img& src, IRoi& srcRoi, Img& dst)
	{
		assert(srcRoi.x0 == 0);
		assert(srcRoi.y0 == 0);
		assert(dst.getWidth() >= srcRoi.getWidth() * 2);
		assert(dst.getHeight() >= srcRoi.getHeight() * 2);
		int w = srcRoi.getWidth();
		int h = srcRoi.getHeight();
		int srcPitch = src.getPitch();
		int dstPitch = dst.getPitch();
		int nb = w / 16;
		int nTail = w - (nb * 16); 
		int srcExtra = srcPitch - w;
		int dstExtra = dstPitch - w * 2;
		byte* srcp = (byte*)src.getData(srcRoi.x0, srcRoi.y0);
		byte* dstp1 = (byte*)dst.getData();
		byte* dstp2 = dstp1 + dstPitch;
		for (int y = 0; y < h; y++)
		{
			for (int i = 0; i < nb; i++)
			{
				__m128i a = _mm_loadu_si128((const __m128i*)srcp);
				__m128i t1 = _mm_unpacklo_epi8(a, a); 
				_mm_storeu_si128((__m128i*)dstp1, t1); 
				_mm_storeu_si128((__m128i*)dstp1, t1); 
				_mm_storeu_si128((__m128i*)dstp2, t1);
				dstp1 += 16;
				dstp2 += 16;
				t1 = _mm_unpackhi_epi8(a, a); 
				_mm_storeu_si128((__m128i*)dstp1, t1); 
				_mm_storeu_si128((__m128i*)dstp2, t1);
				dstp1 += 16;
				dstp2 += 16;
				srcp += 16;
			}
			for (int i = 0; i < nTail; i++)
			{
				byte val = *srcp;
				dstp1[0] = val;
				dstp1[1] = val;
				dstp2[0] = val;
				dstp2[1] = val;
				srcp++;
				dstp1 += 2;
				dstp2 += 2;
			}
			srcp += srcExtra;
			dstp1 += dstExtra + dstPitch;
			dstp2 += dstExtra + dstPitch;
		}
	}
public:
	static void Upsample(Img& src, Img& dst, bool doInterpolate = false)
	{
		assert(!doInterpolate && "interpolation not implemented in upsample yet");
		int width = src.getWidth();
		int height = src.getHeight();
		IRoi roi(0, 0, width, height);
		upSampleSubSse2(src, roi, dst);
		int dw = dst.getWidth();
		int dh = dst.getHeight();
		if (dw > width * 2)
		{
			for (int y = 0; y < dh; y++)
			{
				byte val = dst.getPixel(width * 2 - 1, y);
				for (int x = width * 2; x < dw; x++)
				{
					dst.setPixel(x, y, val);
				}
			}
		}
		if (dh > height * 2)
		{
			for (int y = height * 2; y < dh; y++)
			{
				for (int x = 0; x < dw; x++)
				{
					byte val = dst.getPixel(x, height * 2 - 1);
					dst.setPixel(x, y, val);
				}
			}
		}
	}
	static void Upsample(Img& src, int nTimes, Img& dst, bool doFilter = false)
	{
		int w = src.getWidth();
		int h = src.getHeight();
		if (dst.getWidth() < (w << nTimes))
		{
			ErrorExit("Dest image isn't wide enough for the upsample.");
		}
		if (dst.getHeight() < (h << nTimes))
		{
			ErrorExit("Dest image isn't high enough for the upsample.");
		}
		Img tmp;
		if (nTimes > 1)
		{
			tmp.resize(w << (nTimes - 1), h << (nTimes - 1), 8);
		}
		bool inDst = false;
		IRoi srcRoi(0, 0, w, h);
		for (int i = 0; i < nTimes; i++)
		{
			if (i == 0)
			{
				if (nTimes % 2 == 0)
				{
					upSampleSubSse2(src, srcRoi, tmp);
					if (doFilter) Mean3x3(tmp, true);
					inDst = false;
				}
				else
				{
					upSampleSubSse2(src, srcRoi, dst);
					if (doFilter) Mean3x3(dst, true);
					inDst = true;
				}
			}
			else
			{
				if (!inDst)
				{
					upSampleSubSse2(tmp, srcRoi, dst);
					if (doFilter) Mean3x3(dst, true);
					inDst = true;
				}
				else
				{
					upSampleSubSse2(dst, srcRoi, tmp);
					if (doFilter) Mean3x3(tmp, true);
					inDst = false;
				}
			}
			srcRoi.x1 *= 2;
			srcRoi.y1 *= 2;
		}
	}
	static uint64_t computeProduct(byte* __restrict psrc1, int src1Pitch, byte* __restrict psrc2, int src2Pitch, pair<int, int> src2RoiSize)
	{
		int w = src2RoiSize.first;
		int h = src2RoiSize.second;
		uint64_t sum = 0;
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				sum += psrc1[x] * psrc2[x];
			}
			psrc1 += src1Pitch;
			psrc2 += src2Pitch;
		}
		return sum;
	}
	static void myCrossCorrValid(byte* psrc1, int src1Pitch, pair<int, int> src1RoiSize,
		byte* psrc2, int src2Pitch, pair<int, int> src2RoiSize,
		pair<int, int>& maxOffset)
	{
		int xp = src1RoiSize.first - src2RoiSize.first + 1;
		int yp = src1RoiSize.second - src2RoiSize.second + 1;
		uint64_t maxProduct = 0;
		for (int yi = 0; yi < yp; yi++)
		{
			for (int xi = 0; xi < xp; xi++)
			{
#ifndef USE_SSE41
				uint64_t prod = computeProduct(psrc1 + yi * src1Pitch + xi, src1Pitch, psrc2, src2Pitch, src2RoiSize);
#endif
				if (prod > maxProduct)
				{
					maxProduct = prod;
					maxOffset.first = xi;
					maxOffset.second = yi;
				}
			}
		}
	}
	static void MyFindMaxCrossCorrelation(Img& src1, IRoi& src1Roi, Img& src2, IRoi& src2Roi, int& dx, int& dy, float& matchScore)
	{
		assert(src1Roi.getWidth() >= src2Roi.getWidth());
		assert(src1Roi.getHeight() >= src2Roi.getHeight());
		assert((src1Roi.getHeight() > src2Roi.getHeight()) || (src1Roi.getWidth() >= src2Roi.getWidth()));
		dx = dy = -1;
		matchScore = -1;
		pair<int, int> src1RoiSize = {src1Roi.width(), src1Roi.height()};
		pair<int, int> src2RoiSize = {src2Roi.width(), src2Roi.height()};
		int src1Pitch = src1.getPitch();
		unsigned char* psrc1 = src1.getData(src1Roi.x0, src1Roi.y0);
		int src2Pitch = src2.getPitch();
		unsigned char* psrc2 = src2.getData(src2Roi.x0, src2Roi.y0);
		pair<int, int> bestMatchOffset;
		myCrossCorrValid(psrc1, src1Pitch, src1RoiSize, psrc2, src2Pitch, src2RoiSize, bestMatchOffset);
		dx = bestMatchOffset.first;
		dy = bestMatchOffset.second;
	}
	static bool MyFindMaxCrossCorrelation(Img& src1, Img& src2, int offsetX, int offsetY, int searchRadius, int& dx, int& dy, float& matchScore)
	{
		dx = dy = 0;
		matchScore = 0;
		int width = src1.getWidth();
		int height = src1.getHeight();
		assert(src2.getWidth() == width);
		assert(src2.getHeight() == height);
		IRoi roi2(0, 0, width, height);
		roi2.offsetCenter(offsetX, offsetY, 0, 0, width, height); 
		roi2.shrinkToFit(searchRadius, searchRadius, width - searchRadius, height - searchRadius);
		int insetX = (width - roi2.width()) / 2 - searchRadius;
		int insetY = (height - roi2.height()) / 2 - searchRadius;
		IRoi roi1(insetX, insetY, width - insetX, height - insetY);
		int odx, ody;
		MyFindMaxCrossCorrelation(src1, roi1, src2, roi2, odx, ody, matchScore);
		dx = odx - searchRadius;
		dy = ody - searchRadius;
		dx -= offsetX;
		dy -= offsetY;
		dx = -dx;
		dy = -dy;
		return true;
	}
	static bool MyFindMaxCrossCorrelationSubPixel(Img& src1, IRoi& src1IRoi, Img& src2, IRoi& src2IRoi, int nUpsamples, float& dx, float& dy, float& matchScore)
	{
		assert(nUpsamples <= 2); 
		int width = src1.getWidth();
		int height = src1.getHeight();
		Img src1Up, src2Up;
		int upWidth = width << nUpsamples;
		int upHeight = height << nUpsamples;
		src1Up.resize(upWidth, upHeight, 8);
		src2Up.resize(upWidth, upHeight, 8);
		Upsample(src1, nUpsamples, src1Up, true); 
		Upsample(src2, nUpsamples, src2Up, true); 
		IRoi src1UpIRoi(src1IRoi.x0 << nUpsamples, src1IRoi.y0 << nUpsamples, src1IRoi.x1 << nUpsamples, src1IRoi.y1 << nUpsamples);
		IRoi src2UpIRoi(src2IRoi.x0 << nUpsamples, src2IRoi.y0 << nUpsamples, src2IRoi.x1 << nUpsamples, src2IRoi.y1 << nUpsamples);
		int idx, idy;
		MyFindMaxCrossCorrelation(src1Up, src1UpIRoi, src2Up, src2UpIRoi, idx, idy, matchScore);
		dx = (float)idx / (1 << nUpsamples);
		dy = (float)idy / (1 << nUpsamples);
		return true;
	}
	static bool MyFindMaxCrossCorrelationSubPixel(Img& src1, Img& src2, int nUpsamples, int offsetX, int offsetY, int searchRadius, float& dx, float& dy, float& matchScore)
	{
		assert(src1.getDepthInBits() == 8);
		assert(src2.getDepthInBits() == 8);
		dx = dy = 0;
		matchScore = 0;
		int width = src1.getWidth();
		int height = src1.getHeight();
		assert(src2.getWidth() == width);
		assert(src2.getHeight() == height);
		IRoi roi2(0, 0, width, height);
		roi2.offsetCenter(offsetX, offsetY, 0, 0, width, height); 
		roi2.shrinkToFit(searchRadius, searchRadius, width - searchRadius, height - searchRadius);
		int insetX = (width - roi2.width()) / 2 - searchRadius;
		int insetY = (height - roi2.height()) / 2 - searchRadius;
		IRoi roi1(insetX, insetY, width - insetX, height - insetY);
		float odx, ody;
		MyFindMaxCrossCorrelationSubPixel(src1, roi1, src2, roi2, nUpsamples, odx, ody, matchScore);
		dx = odx - searchRadius;
		dy = ody - searchRadius;
		dx -= offsetX;
		dy -= offsetY;
		dx = -dx;
		dy = -dy;
		return true;
	}
	public:
	static void MeanFilter(Img& i1, int dim, Img& i2)
	{
		CheckImageArgsNotInPlace(i1, i2);
#ifndef USE_IPP
		ErrorExit("MeanFilter: No IPP.");
#endif
	}
#ifdef USE_IPP_REGISTER
#endif
#ifndef USE_IPP_REGISTER
#endif
	public:
};
class DetectionRecord
{
public:
	int uniqueId; 
	int detectionNumber; 
	int frameIndex; 
	int sourceExtractorId; 
	JulianTime julianTime;
	RaDec raDec;
	bool isNeo; 
	float x; 
	float y; 
	float xReg; 
	float yReg; 
	bool hasOverlap;  
	int overlapDetNum;  
	int dupCount; 
	float minDistNonDup; 
	float brightnessMag; 
	float gaussianFwhm; 
	float elongRatio; 
	float theta; 
	float lineFitRmsError;  
	float deltaMu; 
	bool isReject; 
	bool isComputed; 
	bool isCenterBlob; 
	bool isOffImage; 
	float caVis; 
	ImgMoment nsim;
	vector<int> refs; 
	float dtDays; 
	float vra; 
	float vdec; 
	float vTheta; 
	float getVelocityCombined() 
	{ 
		if ((vra == UNDEF_FEATURE) || (vdec == UNDEF_FEATURE)) return UNDEF_FEATURE;
		return sqrtf(vra * vra + vdec*vdec); 
	}
	float DistReg(pair<int, int>& pt)
	{
		float dx = pt.first - xReg;
		float dy = pt.second - yReg;
		return sqrtf(dx*dx + dy*dy);
	}
	Img nsi; 
	DetectionRecord()
	{
		isComputed = isCenterBlob = isNeo = false;
		isReject = false;
		uniqueId = detectionNumber = frameIndex = sourceExtractorId = -1;
		x = y = xReg = yReg = brightnessMag = gaussianFwhm = elongRatio = theta = lineFitRmsError = deltaMu = -1;
		isOffImage = hasOverlap = false;
		overlapDetNum = -1;
		dtDays = vra = vdec = caVis = 0;
		nsi.resize(IDIM, IDIM, 8);
		Img8Util::Clear(nsi);
		dupCount = 0;
		minDistNonDup = 0;
	}
	bool CheckFeatureValues()
	{
		bool ret = true;
		if ((caVis != caVis) || (caVis == std::numeric_limits<float>::infinity()))
		{ 
			logd.warn("Bad caVis in det %d", this->detectionNumber); 
			ret = false; 
		}
		return ret;
	}
	float ComputePGood()
	{
		if (!nsim.CheckAllGood()) return 0; 
		if (nsim.mRadius < 1.25f) return 0; 
		if (nsim.mRadius < 13) return 0; 
		return 1.0f;
	}
	void GetVPixels(float& dx, float& dy, float& dist)
	{
		double pxPerDegree = MEAN_PX_PER_DEGREE; 
		dx = pxPerDegree * -vra * dtDays;
		dy = pxPerDegree * vdec * dtDays;
		double v = getVelocityCombined();
		dist = v * pxPerDegree * dtDays;
	}
	float GetExpectedElong()
	{
		return 2.2f * getVelocityCombined() / 6.9f;
	}
	std::string ToString()
	{
		char tmp[512];
		sprintf(tmp, "uniqueId: %d, num: %d, frame: %d, isNeo: %c, isReject: %c, ra: %f, dec: %f, pt: (%.1f, %.1f), ptReg: (%.1f, %.1f), t: %s", uniqueId, detectionNumber, frameIndex, 
			isNeo ? 'Y' : 'N',
			isReject ? 'Y' : 'N',
			raDec.RaDegrees, raDec.DecDegrees, x, y, xReg, yReg, julianTime.ToString().c_str());
		return (string)tmp;
	}
	static bool TryParse(std::string& s, DetectionRecord& rec)
	{
		vector<string> words;
		StringUtils::splitString(s, ' ', words);
		if (words.size() == 8)
		{
			int i = 0;
			rec.detectionNumber = atoi(words[i++].c_str());
			rec.frameIndex = atoi(words[i++].c_str()) - 1; 
			rec.raDec.RaDegrees = atof(words[i++].c_str()); 
			rec.raDec.DecDegrees = atof(words[i++].c_str()); 
			rec.x = atof(words[i++].c_str());
			rec.y = atof(words[i++].c_str());
			rec.brightnessMag = atof(words[i++].c_str());
			rec.isNeo = (atoi(words[i++].c_str()) > 0);
			return true;
		}
		else
		{
			logd.printf("DetectionRecord: TryParse: Skip: %s\n", s.c_str());
			return false;
		}
	}
	static void WriteTruthCsv(std::string& truthCsvPath, vector<DetectionRecord>& dets)
	{
		FILE* csvp = fopen(truthCsvPath.c_str(), "w");
		fprintf(csvp, "dn,isReject\n");
		for (int i = 0; i < dets.size(); i++)
		{
			int dn = dets[i].detectionNumber;
			if (dets[i].frameIndex == 4)
			{
				fprintf(csvp, "%d,%d\n", dn, dets[i].isReject);
			}
		}
		fclose(csvp);
	}
};
class CFeatureStats
{
public:
	float minimum;
	float maximum;
	float mean;
	float sigma;
	float sigmaPct;
	int count; 
	byte clusterSelectionMask;
	byte clusterSelectedCount;
	vector<float> clusterSelectionWeights;
	CFeatureStats()
	{
		minimum = maximum = mean = sigma = sigmaPct = FLT_MAX;
		clusterSelectionMask = clusterSelectedCount = 0;
	}
	static std::string ToCsvHeaderString(string baseName)
	{
		char tmp[4096];
		sprintf(tmp, "%sMin, %sMax, %sMean, %sSigmaPct", baseName.c_str(), baseName.c_str(), baseName.c_str(), baseName.c_str());
		return (string)tmp;
	}
	std::string ToString()
	{
		char tmp[4096];
		sprintf(tmp, "count: %d, selCount: %d, min: %f, max: %f, mean: %f, sigma: %f, sigmaPct: %f, mask: %s", 
			count, clusterSelectedCount, minimum, maximum, mean, sigma, sigmaPct, toStringBits(clusterSelectionMask).c_str());
		return (string)tmp;
	}
	std::string ToCsvString()
	{
		char tmp[4096];
		sprintf(tmp, "%f, %f, %f, %f", 
			(minimum != FLT_MAX ? minimum : UNDEF_FEATURE), 
			(maximum != FLT_MAX ? maximum : UNDEF_FEATURE), 
			(mean != FLT_MAX ? mean : UNDEF_FEATURE), 
			(sigmaPct != FLT_MAX ? sigmaPct : UNDEF_FEATURE)
			);
		return (string)tmp;
	}
	void UpdateSigmaPct()
	{
		if (mean != 0)
		{
			sigmaPct = abs(sigma / mean);
		}
		else
		{
			sigmaPct = 0;
		}
	}
	void ComputeClustered(vector<float>& values, bool skipZeros, bool debug)
	{
		int n = values.size();
		count = n;
		byte initialMask = ((1 << n) - 1);
		for (int i = 0; i < n; i++) 
		{
			if ((values[i] != UNDEF_FEATURE) && (!skipZeros || (values[i] != 0)))
			{
				initialMask = initialMask | (1 << i);
			}
			else
			{
				initialMask = initialMask & ~(1 << i);
			}
		}
		if (initialMask == 0)
		{
			minimum = 0;
			maximum = 0;
			mean = UNDEF_FEATURE;
			sigma = UNDEF_FEATURE;
			sigmaPct = UNDEF_FEATURE;
		}
		else
		{
			float cmean, cstddev;
			clusterSelectedCount = singleCluster(values, MIN_FEAT_CLUSTER_COUNT, initialMask, clusterSelectionMask, cmean, cstddev, clusterSelectionWeights, false);
			minimum = FLT_MAX;
			maximum = -FLT_MAX;
			for (int i = 0; i < n; i++)
			{
				if ((clusterSelectionMask >> i) & 1)
				{
					minimum = MIN(minimum, values[i]);
					maximum = MAX(maximum, values[i]);
				}
			}
			mean = cmean;
			sigma = cstddev;
			UpdateSigmaPct();
		}
	}
	void ComputeClustered(float a, float b, float c, float d, bool skipZeros, bool debug)
	{
		vector<float> values(4);
		values[0] = a;
		values[1] = b;
		values[2] = c;
		values[3] = d;
		ComputeClustered(values, skipZeros, debug);
	}
	void ComputeClustered(float a, float b, float c, bool skipZeros, bool debug)
	{
		vector<float> values(3);
		values[0] = a;
		values[1] = b;
		values[2] = c;
		ComputeClustered(values, skipZeros, debug);
	}
	void ComputeWeighted(vector<float>& values, vector<float>& weights, bool skipZeros, bool debug)
	{
		count = values.size();
		clusterSelectionMask = 0;
		vector<float> finalWeights = weights;
		for (int i = 0; i < count; i++)
		{
			if ((values[i] != UNDEF_FEATURE) && (!skipZeros || (values[i] != 0)))
			{
				clusterSelectionMask |= (1 << i);
			}
			else
			{
				finalWeights[i] = 0;
				clusterSelectionMask &= ~(1 << i);
			}
		}
		clusterSelectedCount = 0;
		for (int i = 0; i < count; i++)
		{
			if (finalWeights[i] == 0) clusterSelectionMask &= ~(1 << i);
			else clusterSelectedCount++;
		}
		float meandev;
		Stats::ComputeStats(values, finalWeights, minimum, maximum, mean, sigma, meandev);
		UpdateSigmaPct();
	}
	void ComputeWeighted(vector<float>& weights, float a, float b, float c, float d, bool skipZeros, bool debug)
	{
		vector<float> values(4);
		values[0] = a;
		values[1] = b;
		values[2] = c;
		values[3] = d;
		ComputeWeighted(values, weights, skipZeros, debug);
	}
	void ComputeWeighted(vector<float>& weights, float a, float b, float c, bool skipZeros, bool debug)
	{
		vector<float> values(3);
		values[0] = a;
		values[1] = b;
		values[2] = c;
		ComputeWeighted(values, weights, skipZeros, debug);
	}
	static void Test()
	{
		CFeatureStats s;
		vector<float> weights;
		weights.push_back(1.0f);
		weights.push_back(1.0f);
		weights.push_back(1.0f);
		s.ComputeWeighted(weights, 3.2f, 3.1f, 5.0f, true, true);
		logd.printf("s: %s\n", s.ToString().c_str());
		s.ComputeWeighted(weights, 3.2f, 0.0f, 5.0f, true, true);
		logd.printf("s: %s\n", s.ToString().c_str());
		weights[2] = 0;
		s.ComputeWeighted(weights, 3.2f, 3.1f, 5.0f, true, true);
		logd.printf("s: %s\n", s.ToString().c_str());
	}
};
static Img ts1, ts2, ts3, ts4; 
static Img tb1, tb2, tb3, tb4; 
#include <vector>
#include <queue>
#include <emmintrin.h>
#include <string>
#include <vector>
#include <cstdint>
class Hist
{
private:
	bool isStatsComputed;
	float meanBinIndex;
	float variance; 
	double weightSum; 
public:
	vector<int> bins;
	vector<int> levels; 
	int minBinIndex; 
	int maxBinIndex; 
	uint64_t sampleCount; 
	Hist()
	{
		Init();
	}
private:
	void Init()
	{
		minBinIndex = maxBinIndex = 0;
		sampleCount = 0;
		isStatsComputed = false;
		meanBinIndex = weightSum = variance = 0;
	}
	void ComputeStats()
	{
		if (!isStatsComputed)
		{
			int n = bins.size();
			assert(n > 0);
			weightSum = 0;
			sampleCount = 0;
			for (int i = 0; i < n; i++)
			{
				sampleCount += (uint64_t)bins[i];
				weightSum += (uint64_t)i * bins[i];
			}
			if (sampleCount > 0)
			{
				this->meanBinIndex = weightSum / sampleCount;
			}
			else
			{
				this->meanBinIndex = 0;
				this->variance = 0;
			}
			isStatsComputed = true;
			float meanVal = GetMean();
			double varSum = 0;
			for (int i = 0; i < n; i++)
			{
				if (bins[i] > 0)
				{
					double level = (levels[i] + levels[i + 1]) / 2;
					double dv = level - meanVal;
					varSum += bins[i] * dv * dv;
				}
			}
			variance = varSum / sampleCount;
		}
	}
public:
	void Resize(int nBins)
	{
		if (nBins != this->bins.size())
		{
			this->bins.resize(nBins);
			this->levels.resize(nBins + 1);
		}
		Init();
	}
	void Clear()
	{
		int n = this->bins.size();
		for (int i = 0; i < n; i++)
		{
			this->bins[i] = 0;
		}
		n = this->levels.size();
		for (int i = 0; i < n; i++)
		{
			this->levels[i] = 0;
		}
		Init();
	}
	void Copy(const Hist& other)
	{
		Init();
		this->bins = other.bins;
		this->levels = other.levels;
		this->maxBinIndex = other.maxBinIndex;
		this->minBinIndex = other.minBinIndex;
		this->sampleCount = other.sampleCount;
		this->isStatsComputed = other.isStatsComputed;
		this->meanBinIndex = other.meanBinIndex;
		this->weightSum = other.weightSum;
		this->variance = other.variance;
	}
	void Recompute()
	{
		sampleCount = 0;
		int n = this->bins.size();
		minBinIndex = n;
		maxBinIndex = -1;
		for (int i = 0; i < n; i++)
		{
			if (bins[i] > 0)
			{
				if (minBinIndex == n)
				{
					minBinIndex = i;
				}
				maxBinIndex = i;
				sampleCount += bins[i];
			}
		}
	}
	void SetData(vector<int>& bins, vector<int>& levels)
	{
		assert(levels.size() == bins.size() + 1);
		this->bins = bins;
		this->levels = levels;
		isStatsComputed = false;
		Recompute();
		ComputeStats();
	}
	void SetData(int* bins, int nBins, int* levels)
	{
		this->bins.clear();
		this->bins.reserve(nBins);
		this->levels.clear();
		this->levels.reserve(nBins + 1);
		for (int i = 0; i < nBins; i++)
		{
			this->bins.push_back(bins[i]);
			this->levels.push_back(levels[i]);
		}
		this->levels.push_back(levels[nBins]);
		isStatsComputed = false;
		Recompute();
		ComputeStats();
	}
	void ClearBinValues()
	{
		isStatsComputed = false;
		SetBinValues(0, bins.size() - 1, 0);
	}
	void SetBinValues(int startIndex, int endIndex, int newValue)
	{
		assert(startIndex >= 0);
		assert(startIndex < endIndex);
		assert(endIndex <= bins.size());
		for (int i = startIndex; i < endIndex; i++)
		{
			bins[i] = newValue;
		}
		isStatsComputed = false;
		Recompute();
	}
	void SetBinValue(int binIndex, int newValue)
	{
		isStatsComputed = false;
		int n = bins.size();
		assert(binIndex >= 0);
		assert(binIndex < n);
		int oldValue = bins[binIndex];
		sampleCount -= oldValue;
		sampleCount += newValue;
		bins[binIndex] = newValue;
		if (newValue == 0)
		{
			if (minBinIndex == binIndex)
			{
				for (int i = minBinIndex; i < n; i++)
				{
					if (bins[i] > 0)
					{
						minBinIndex = i;
						break;
					}
				}
			}
			else if (maxBinIndex == binIndex)
			{
				for (int i = maxBinIndex; i >= 0; i--)
				{
					if (bins[i] > 0)
					{
						maxBinIndex = i;
						break;
					}
				}
			}
		}
	}
	float BinIndexToLevel(float binIndex)
	{
		int i = (int)binIndex;
		float f = binIndex - i;
		float v = (1 - f) * levels[i] + f * levels[i + 1];
		return v;
	}
	int GetBinIndexWithMaxCount()
	{
		int maxCount = 0;
		int maxIndex = 0;
		int n = bins.size();
		for (int i = 0; i < n; i++)
		{
			if (bins[i] > maxCount)
			{
				maxCount = bins[i];
				maxIndex = i;
			}
		}
		return maxIndex;
	}
	float GetMinLevel()
	{
		return BinIndexToLevel(minBinIndex);
	}
	float GetMaxLevel()
	{
		return BinIndexToLevel(maxBinIndex);
	}
	float GetMeanBinIndex()
	{
		ComputeStats();
		return meanBinIndex;
	}
	float GetMean()
	{
		ComputeStats();
		return BinIndexToLevel(meanBinIndex);
	}
	float GetVariance()
	{
		ComputeStats();
		return variance;
	}
	float GetStddev()
	{
		return sqrt(GetVariance());
	}
	float ComputeWeightedSum()
	{
		ComputeStats();
		return (float)weightSum;
	}
	int GetPercentileBinIndex(float pct)
	{
		int targetCount = roundf(pct * sampleCount);
		int n = bins.size();
		if (targetCount >= sampleCount) return maxBinIndex;
		int result = -1;
		int runningSum = 0;
		for (int i = 0; i < n; i++)
		{
			runningSum += bins[i];
			if (runningSum > targetCount)
			{
				result = i;
				break;
			}
		}
		return result;
	}
	int GetPercentileBinLevel(float pct)
	{
		int targetCount = roundf(pct * sampleCount);
		int n = bins.size();
		if (targetCount >= sampleCount) return levels[maxBinIndex];
		int result = -1;
		int runningSum = 0;
		for (int i = 0; i < n; i++)
		{
			runningSum += bins[i];
			if (runningSum > targetCount)
			{
				result = levels[i];
				break;
			}
		}
		return result;
	}
	void SaveToCsv(const string& path)
	{
		logd.debug("Hist: Save to csv: %s", path.c_str());
	}
	string ToString()
	{
		char tmp[256];
		float var = GetVariance();
		sprintf(tmp, "nBins: %d, min: %.1f (bin %u), mean: %.1f (bin %.1f), max: %.1f (bin %d), sampleCount: %lld, var/stddev: %.1f/%.1f", 
			this->bins.size(), GetMinLevel(), minBinIndex, GetMean(), meanBinIndex, GetMaxLevel(), maxBinIndex, sampleCount, var, sqrtf(var));
		return(string(tmp));
	}
	void Display(const string& prefix)
	{
		logd.printf("%s hist = %s\n", prefix.c_str(), ToString().c_str());
		if (bins.size() == 0)
		{
			logd.printf("%s hist = NULL\n", prefix.c_str());
		}
		else
		{
			logd.printf("%s bins:\n", prefix.c_str());
			for (int i = 0; i < bins.size(); i++)
			{
				logd.printf("%s    %d: %d to %d: %d\n", prefix.c_str(), i, levels[i], levels[i + 1], bins[i]);
			}
		}
	}
};
class ImgUtil
{
public:
	static void CheckImageArgsSize(Img& src, Img& dst)
	{
		if (src.getWidth() != dst.getWidth())
		{
			ErrorExit("CheckImageArgsWidth: Different width.");
		}
		if (src.getHeight() != dst.getHeight())
		{
			ErrorExit("CheckImageArgsHeight: Different height.");
		}
	}
	static void CheckImageArgs16InPlace(Img& src, Img& dst)
	{
		CheckImageArgsSize(src, dst);
		if (src.getDepthInBits() != dst.getDepthInBits())
		{
			ErrorExit("CheckImageArgsDepthInBits: Different.");
		}
		if (src.getDepthInBits() != 16)
		{
			ErrorExit("CheckImageArgsPitch: Supposed to be 16 bit.");
		}
	}
	static void CheckImageArgsNotInPlace(Img& src, Img& dst)
	{
		CheckImageArgsSize(src, dst);
		if (src.getDepthInBits() != dst.getDepthInBits())
		{
			ErrorExit("CheckImageArgsDepthInBits: Different.");
		}
		if (src.getData() == dst.getData())
		{
			ErrorExit("CheckImageArgsPitch: Not supposed to be same buffer.");
		}
	}
	static void CheckImageArgs16NotInPlace(Img& src, Img& dst)
	{
		CheckImageArgs16InPlace(src, dst);
		if (src.getData() == dst.getData())
		{
			ErrorExit("CheckImageArgsPitch: Not supposed to be same buffer.");
		}
	}
	static void CheckImageArgsMask(Img& src, Img& mask)
	{
		CheckImageArgsSize(src, mask);
		if (src.getDepthInBits() != 16)
		{
			ErrorExit("CheckImageArgsPitch: Supposed to be 16 bit.");
		}
		if (mask.getDepthInBits() != 8)
		{
			ErrorExit("CheckImageArgsPitch: Supposed to be 8 bit.");
		}
		if (src.getData() == mask.getData())
		{
			ErrorExit("CheckImageArgsPitch: Not supposed to be same buffer.");
		}
	}
	static void Clear(Img& src)
	{
		memset(src.getData(), 0, src.getBufferSize());
	}
	static void Copy(Img& src, Img& dst)
	{
		if (src.getDepthInBits() != 16) ErrorExit("Source depth must be 16");
		if (dst.getDepthInBits() != 16) ErrorExit("Dest depth must be 16");
		int w = src.getWidth();
		int h = src.getHeight();
		int spp = src.getPitch() / 2;
		int dpp = dst.getPitch() / 2;
		uint16_t* ps = (uint16_t*)src.getData();
		uint16_t* pd = (uint16_t*)dst.getData();
		for (int y = 0; y < h; y++)
		{
			memcpy(pd, ps, w*2);
			ps += spp;
			pd += dpp;
		}
	}
	static void SetFromInt(int width, int height, vector<int>& rawData, Img& dst)
	{
		if (rawData.size() != width * height)
		{
			logd.error("SetFromInt: Input raw data is not correct size. It is %d instead of %d.", rawData.size(), width * height);
			return;
		}
		dst.resize(width, height, 16);
		for (int y = 0; y < height; y++)
		{
			for (int x = 0; x < width; x++)
			{
				dst.setPixel16u(x, y, rawData[y * width + x]);
			}
		}
	}
	static void ConvertToInt(Img& src, vector<int>& v)
	{
		int w = src.getWidth();
		int h = src.getHeight();
		int spp = src.getPitch() / 2;
		v.resize(w*h);
		uint16_t* ps = (uint16_t*)src.getData();
		int i = 0;
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				v[i++] = (int)ps[x];
			}
			ps += spp;
		}
	}
	static void CopyRoi(Img& src, IRoi srcRoi, Img& dst)
	{
		if (dst.getWidth() < srcRoi.getWidth()) ErrorExit("Dest image is too small for roi.");
		if (dst.getHeight() < srcRoi.getHeight()) ErrorExit("Dest image is too small for roi.");
		IRoi dstRoi(0, 0, srcRoi.getWidth(), srcRoi.getHeight());
		Clear(dst);
		dst.SetRoiAllowOffSrc(dstRoi, src, srcRoi);
	}
	static void SetRoiByPtr(Img& src, IRoi roi, Img& tile)
	{
		SetTile(src, roi, tile);
	}
	static void SetTile(Img& src, IRoi roi, Img& tile)
	{
		tile.setByPtr(roi.width(), roi.height(), src.getDepthInBits(), src.getPitch(), src.getData(roi.x0, roi.y0));
	}
	static int GetGoodValues(Img& src, vector<uint16_t>& values)
	{
		ErrorExit("ILEN");
		assert(false); 
		int n = 0;
		values.resize(ILEN);
		for (int i = 0; i < ILEN; i++) 
		{
			uint16_t v = src.getPixel16u(i);
			if ((v > 0) && (v < 65535))
			{
				values[n++] = v;
			}
		}
		values.resize(n);
		return n;
	}
	static void GetGoodValuesMasked(Img& src, Img& mask, bool isGoodMask, vector<uint16_t>& values)
	{
		CheckImageArgsMask(src, mask);
		int w = src.getWidth();
		int h = src.getHeight();
		int spp = src.getPitch() / 2;
		int mpp = mask.getPitch();
		CheckImageArgsMask(src, mask);
		uint16_t* ps = (uint16_t*)src.getData();
		uint8_t* pm = (uint8_t*)mask.getData();
		values.reserve(w * h);
		if (isGoodMask)
		{
			for (int y = 0; y < h; y++)
			{
				for (int x = 0; x < w; x++)
				{
					if (pm[x]) values.push_back(ps[x]);
				}
				ps += spp;
				pm += mpp;
			}
		}
		else
		{
			for (int y = 0; y < h; y++)
			{
				for (int x = 0; x < w; x++)
				{
					if (pm[x] == 0) values.push_back(ps[x]);
				}
				ps += spp;
				pm += mpp;
			}
		}
		values.shrink_to_fit();
	}
	static int GetGoodValuesRoi(Img& src, IRoi& roi, vector<uint16_t>& values)
	{
		values.resize(roi.getWidth() * roi.getHeight());
		for (int y = roi.y0; y < roi.y1; y++) 
		{
			for (int x = roi.x0; x < roi.x1; x++) 
			{
				uint16_t v = src.getPixel16u(x, y);
				if ((v > 0) && (v < 65535))
				{
					values.push_back(v);
				}
			}
		}
		values.shrink_to_fit();
		return values.size();
	}
	static int GetGoodValuesCenter(Img& src, int centerDiameter, vector<uint16_t>& values)
	{
		ErrorExit("Not updated per IDIM");
		int centerRadius = centerDiameter / 2;
		int cv1 = IDIM/2 - 1 - centerRadius;
		int cv2 = IDIM/2 + centerRadius;
		IRoi roi(cv1, cv1, cv2, cv2);
		int ngc = GetGoodValuesRoi(src, roi, values);
		return ngc;
	}
	static float MeanRoiGoodValues(Img& src, IRoi& roi, int& n)
	{
		float sum = 0;
		vector<uint16_t> values;
		n = GetGoodValuesRoi(src, roi, values); 
		for (int i = 0; i < n; i++) sum += values[i];
		return sum / n;
	}
	static void ZeroMagic(Img& src)
	{
		ErrorExit("ILEN");
		for (int i = 0; i < ILEN; i++) 
		{
			uint16_t v = src.getPixel16u(i);
			src.setPixel16u(i, (v == 65535 ? 0 : v));
		}
	}
	static void ZeroCol(Img& src, int colIndex)
	{
		int h = src.getHeight();
		int spp = src.getPitch() / 2;
		uint16_t* ps = (uint16_t*)src.getData();
		for (int y = 0; y < h; y++)
		{
			ps[y * spp + colIndex] = 0;
		}
	}
	static void ZeroRow(Img& src, int rowIndex)
	{
		int sp = src.getPitch();
		uint16_t* ps = (uint16_t*)src.getData(0,rowIndex);
		memset(ps, 0, sp);
	}
	static void ZeroEdges(Img& src, int borderWidth)
	{
		int w = src.getWidth();
		int h = src.getHeight();
		for (int i = 0; i < borderWidth; i++)
		{
			ZeroCol(src, i);
			ZeroCol(src, w - 1 - i);
			ZeroRow(src, i);
			ZeroRow(src, h - 1 - i);
		}
	}
	static void Mask(Img& src, Img& mask)
	{
		ErrorExit("ILEN");
		for (int i = 0; i < ILEN; i++)
		{
			uint16_t v = src.getPixel16u(i);
			byte m = mask.getPixel(i);
			src.setPixel16u(i, (m > 0 ? v : 0));
		}
	}
	static void MaskNot(Img& src, Img& mask)
	{
		Mask(src, mask, true);
	}
	static void Mask(Img& src, Img& mask, bool doNot)
	{
		int w = src.getWidth();
		int h = src.getHeight();
		int spp = src.getPitch() / 2;
		int mpp = mask.getPitch();
		CheckImageArgsMask(src, mask);
		uint16_t* ps = (uint16_t*)src.getData();
		uint8_t* pm = (uint8_t*)mask.getData();
		if (doNot)
		{
			for (int y = 0; y < h; y++)
			{
				for (int x = 0; x < w; x++)
				{
					uint16_t m = (pm[x] | (pm[x] << 8));
					ps[x] = ps[x] & ~m;
				}
				ps += spp;
				pm += mpp;
			}
		}
		else
		{
			for (int y = 0; y < h; y++)
			{
				for (int x = 0; x < w; x++)
				{
					uint16_t m = (pm[x] | (pm[x] << 8));
					ps[x] = ps[x] & m;
				}
				ps += spp;
				pm += mpp;
			}
		}
	}
	static void MaskSaturate(Img& src, Img& mask)
	{
		int w = src.getWidth();
		int h = src.getHeight();
		int spp = src.getPitch() / 2;
		int mpp = mask.getPitch();
		CheckImageArgsMask(src, mask);
		uint16_t* ps = (uint16_t*)src.getData();
		uint8_t* pm = (uint8_t*)mask.getData();
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				uint16_t m = (pm[x] | (pm[x] << 8));
				ps[x] = ps[x] | m;
			}
			ps += spp;
			pm += mpp;
		}
	}
	static void Mean2x2(Img& src, Img& dst)
	{
		ErrorExit("Not updated per IDIM");
		assert(src.getWidth() == IDIM);
		assert(src.getData() != dst.getData());
		uint16_t* ps = (uint16_t*)src.getData();
		uint16_t* pd = (uint16_t*)dst.getData();
		memset(pd, 0, IBYTES);
		for (int y = 0; y < IDIM - 1; y++)
		{
			int i = y * IDIM;
			for (int x = 0; x < IDIM - 1; x++)
			{
				int sum = ps[i] + ps[i + 1] + ps[i + IDIM] + ps[i + IDIM + 1];
				pd[i] = (uint16_t)(sum >> 2);
				i++;
			}
		}
	}
	static void Mean2x2GoodValues(Img& src, Img& dst)
	{
		ErrorExit("Not updated per IDIM");
		assert(src.getWidth() == IDIM);
		assert(src.getData() != dst.getData());
		uint16_t* ps = (uint16_t*)src.getData();
		uint16_t* pd = (uint16_t*)dst.getData();
		memset(pd, 0, IBYTES);
		for (int y = 0; y < IDIM - 1; y++)
		{
			int i = y * IDIM;
			for (int x = 0; x < IDIM - 1; x++)
			{
				uint16_t a = ps[i];
				uint16_t b = ps[i + 1];
				uint16_t c = ps[i + IDIM];
				uint16_t d = ps[i + IDIM + 1];
				int n = 0;
				int sum = 0;
				if (a > 0) { n++; sum += a; }
				if (b > 0) { n++; sum += b; }
				if (c > 0) { n++; sum += c; }
				if (d > 0) { n++; sum += d; }
				if ((n > 0) && (a > 0))
				{
					pd[i] = (uint16_t)(sum / n);
				}
				i++;
			}
		}
	}
	static void Mean2x2(Img& src, Img& mask, Img& dst)
	{
		int w = src.getWidth();
		int h = src.getHeight();
		int spp = src.getPitch() / 2;
		int dpp = dst.getPitch() / 2;
		int mpp = mask.getPitch();
		CheckImageArgs16NotInPlace(src, dst);
		CheckImageArgsMask(src, mask);
		uint16_t* ps = (uint16_t*)src.getData();
		uint8_t* pm = (uint8_t*)mask.getData();
		uint16_t* pd = (uint16_t*)dst.getData();
		for (int y = 0; y < h - 1; y++)
		{
			for (int x = 0; x < w - 1; x++)
			{
				uint16_t a = ps[x];
				uint16_t b = ps[x + 1];
				uint16_t c = ps[x + spp];
				uint16_t d = ps[x + spp + 1];
				int n = 0;
				int sum = 0;
				if (pm[x] == 0) { n++; sum += a; }
				if (pm[x + 1] == 0) { n++; sum += b; }
				if (pm[x + spp] == 0) { n++; sum += c; }
				if (pm[x + spp + 1] == 0) { n++; sum += d; }
				if ((n > 0) && (a > 0))
				{
					pd[x] = (uint16_t)(sum / n);
				}
				else
				{
					pd[x] = 0;
				}
			}
			ps += spp;
			pm += mpp;
			pd += dpp;
		}
	}
	static void CopyBorders(Img& src, int borderWidth, Img& dst)
	{
		int w = src.getWidth();
		int h = src.getHeight();
		int spp = src.getPitch() / 2;
		int dpp = dst.getPitch() / 2;
		CheckImageArgs16NotInPlace(src, dst);
		uint16_t* ps = (uint16_t*)src.getData();
		uint16_t* pd = (uint16_t*)dst.getData();
		for (int y = 0; y < h; y++)
		{
			for (int i = 0; i < borderWidth; i++)
			{
				pd[i] = ps[i];
			}
			for (int i = w - borderWidth; i < w; i++)
			{
				pd[i] = ps[i];
			}
			ps += spp;
			pd += dpp;
		}
		ps = (uint16_t*)src.getData();
		pd = (uint16_t*)dst.getData();
		for (int y = 0; y < borderWidth; y++)
		{
			memcpy(&pd[y*dpp], &ps[y*spp], w * 2);
		}
		for (int y = h - borderWidth; y < h; y++)
		{
			memcpy(&pd[y*dpp], &ps[y*spp], w * 2);
		}
	}
	static void Mean3x3GoodValues(Img& src, Img& dst, bool includeBorders = true)
	{
		CheckImageArgs16NotInPlace(src, dst);
		int w = src.getWidth();
		int h = src.getHeight();
		int spp = src.getPitch() / 2;
		int dpp = dst.getPitch() / 2;
		uint16_t* ps = (uint16_t*)src.getData();
		uint16_t* pd = (uint16_t*)dst.getData();
		const int off[9] = { -spp-1, -spp, -spp+1, -1, 0, 1, spp-1,spp,spp+1};
		ps += spp;
		pd += dpp;
		for (int y = 1; y < h - 1; y++)
		{
			for (int x = 1; x < w - 1; x++)
			{
				int n = 0;
				int sum = 0;
				if (ps[x] > 0)
				{
					for (int j = 0; j < 9; j++)
					{
						uint16_t v = ps[x + off[j]];
						if (v > 0) { n++; sum += v; }
					}
					if (n > 0)
					{
						pd[x] = (uint16_t)(sum / n);
					}
					else
					{
						pd[x] = 0;
					}
				}
				else
				{
					pd[x] = 0;
				}
			}
			ps += spp;
			pd += dpp;
		}
		if (!includeBorders)
		{
			CopyBorders(src, 1, dst);
		}
		else
		{
			ps = (uint16_t*)src.getData();
			pd = (uint16_t*)dst.getData();
			const int boff[6] = {-1, 0, 1, spp - 1, spp, spp + 1};
			if (ps[0] > 0)
			{
				const int bloff[4] = {0, 1, spp, spp + 1};
				int x = 0;
				int n = 0;
				int sum = 0;
				for (int j = 0; j < 4; j++)
				{
					uint16_t v = ps[x + bloff[j]];
					if (v > 0) { n++; sum += v; }
				}
				if (n > 0)
				{
					pd[x] = (uint16_t)(sum / n);
				}
				else
				{
					pd[x] = 0;
				}
			}
			else
			{
				pd[0] = 0;
			}
			for (int x = 1; x < w - 1; x++)
			{
				int n = 0;
				int sum = 0;
				if (ps[x] > 0)
				{
					for (int j = 0; j < 6; j++)
					{
						uint16_t v = ps[x + boff[j]];
						if (v > 0) { n++; sum += v; }
					}
					if (n > 0)
					{
						pd[x] = (uint16_t)(sum / n);
					}
					else
					{
						pd[x] = 0;
					}
				}
				else
				{
					pd[x] = 0;
				}
			}
			if (ps[w-1] > 0)
			{
				const int broff[4] = {-1, 0, spp - 1, spp};
				int x = w-1;
				int n = 0;
				int sum = 0;
				for (int j = 0; j < 4; j++)
				{
					uint16_t v = ps[x + broff[j]];
					if (v > 0) { n++; sum += v; }
				}
				if (n > 0)
				{
					pd[x] = (uint16_t)(sum / n);
				}
				else
				{
					pd[x] = 0;
				}
			}
			else
			{
				pd[w-1] = 0;
			}
			ps = (uint16_t*)src.getData(0, h - 1);
			pd = (uint16_t*)dst.getData(0, h - 1);
			const int toff[9] = {-spp - 1, -spp, -spp + 1, -1, 0, 1};
			if (ps[0] > 0)
			{
				const int tloff[4] = {0, 1, -spp, -spp + 1};
				int x = 0;
				int n = 0;
				int sum = 0;
				for (int j = 0; j < 4; j++)
				{
					uint16_t v = ps[x + tloff[j]];
					if (v > 0) { n++; sum += v; }
				}
				if (n > 0)
				{
					pd[x] = (uint16_t)(sum / n);
				}
				else
				{
					pd[x] = 0;
				}
			}
			else
			{
				pd[0] = 0;
			}
			for (int x = 1; x < w - 1; x++)
			{
				int n = 0;
				int sum = 0;
				if (ps[x] > 0)
				{
					for (int j = 0; j < 6; j++)
					{
						uint16_t v = ps[x + toff[j]];
						if (v > 0) { n++; sum += v; }
					}
					if (n > 0)
					{
						pd[x] = (uint16_t)(sum / n);
					}
					else
					{
						pd[x] = 0;
					}
				}
				else
				{
					pd[x] = 0;
				}
			}
			if (ps[w - 1] > 0)
			{
				const int troff[4] = {-1, 0, -spp - 1, -spp};
				int x = w - 1;
				int n = 0;
				int sum = 0;
				for (int j = 0; j < 4; j++)
				{
					uint16_t v = ps[x + troff[j]];
					if (v > 0) { n++; sum += v; }
				}
				if (n > 0)
				{
					pd[x] = (uint16_t)(sum / n);
				}
				else
				{
					pd[x] = 0;
				}
			}
			else
			{
				pd[w-1] = 0;
			}
			ps = (uint16_t*)src.getData(0, 1);
			pd = (uint16_t*)dst.getData(0, 1);
			const int loff[6] = {-spp, -spp + 1, 0, 1, spp, spp + 1};
			for (int y = 1; y < h - 1; y++)
			{
				int n = 0;
				int sum = 0;
				if (ps[0] > 0)
				{
					for (int j = 0; j < 6; j++)
					{
						uint16_t v = ps[loff[j]];
						if (v > 0) { n++; sum += v; }
					}
					if (n > 0)
					{
						pd[0] = (uint16_t)(sum / n);
					}
					else
					{
						pd[0] = 0;
					}
				}
				else
				{
					pd[0] = 0;
				}
				ps += spp;
				pd += dpp;
			}
			ps = (uint16_t*)src.getData(w-1, 1);
			pd = (uint16_t*)dst.getData(w-1, 1);
			const int roff[6] = {-spp - 1, -spp, -1, 0, spp - 1, spp};
			for (int y = 1; y < h - 1; y++)
			{
				int n = 0;
				int sum = 0;
				if (ps[0] > 0)
				{
					for (int j = 0; j < 6; j++)
					{
						uint16_t v = ps[roff[j]];
						if (v > 0) { n++; sum += v; }
					}
					if (n > 0)
					{
						pd[0] = (uint16_t)(sum / n);
					}
					else
					{
						pd[0] = 0;
					}
				}
				else
				{
					pd[0] = 0;
				}
				ps += spp;
				pd += dpp;
			}
		}
	}
	static void NStddevOver8bit(Img& src, float mean, float stddev, float minRange, float maxRange, Img& dst)
	{
		CheckImageArgsMask(src, dst);
		int imean = (int)ceilf(mean);
		int w = src.getWidth();
		int h = src.getHeight();
		int spp = src.getPitch() / 2;
		int dpp = dst.getPitch();
		uint16_t* ps = (uint16_t*)src.getData();
		uint8_t* pd = (uint8_t*)dst.getData();
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				if (ps[x] > imean)
				{
					float d = (ps[x] - mean) / stddev;
					float val = rangeScore(minRange, maxRange, d); 
					pd[x] = (byte)CLIP(roundf(val * 255), 0, 255);
				}
				else
				{
					pd[x] = 0;
				}
			}
			ps += spp;
			pd += dpp;
		}
	}
	static void NStddevOver8bit(Img& src, Img& meanImage, float stddev, float minRange, float maxRange, Img& dst)
	{
		CheckImageArgsMask(src, dst);
		int w = src.getWidth();
		int h = src.getHeight();
		int spp = src.getPitch() / 2;
		int mpp = meanImage.getPitch() / 2;
		int dpp = dst.getPitch();
		uint16_t* ps = (uint16_t*)src.getData();
		uint16_t* pm = (uint16_t*)meanImage.getData();
		uint8_t* pd = (uint8_t*)dst.getData();
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				uint16_t m = pm[x];
				if (ps[x] > m)
				{
					float d = (ps[x] - m) / stddev;
					float val = rangeScore(minRange, maxRange, d); 
					pd[x] = (byte)CLIP(roundf(val * 255), 0, 255);
				}
				else
				{
					pd[x] = 0;
				}
			}
			ps += spp;
			pm += mpp;
			pd += dpp;
		}
	}
	static void Equals(Img& src, uint16_t value, Img& dst)
	{
		CheckImageArgsMask(src, dst);
		int w = src.getWidth();
		int h = src.getHeight();
		int spp = src.getPitch() / 2;
		int dpp = dst.getPitch();
		uint16_t* ps = (uint16_t*)src.getData();
		uint8_t* pd = (uint8_t*)dst.getData();
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				if (ps[x] == value)
				{
					pd[x] = (byte)255;
				}
				else
				{
					pd[x] = 0;
				}
			}
			ps += spp;
			pd += dpp;
		}
	}
	static void ThresholdZero(Img& src, uint16_t threshold, bool isLte, Img& dst)
	{
		CheckImageArgs16InPlace(src, dst);
		uint16_t* __restrict ps = (uint16_t*)src.getData();
		uint16_t* __restrict pd = (uint16_t*)dst.getData();
		int n = src.getHeight() * src.getPitch() / 2;
		if (isLte)
		{
			for (int i = 0; i < n; i++) 
			{
				pd[i] = (ps[i] <= threshold ? (uint16_t)0 : ps[i]);
			}
		}
		else
		{
			for (int i = 0; i < n; i++) 
			{
				pd[i] = (ps[i] >= threshold ? (uint16_t)0 : ps[i]);
			}
		}
	}
	static void ThresholdZeroSse2(Img& src, uint16_t threshold, bool isLte, Img& dst)
	{
		CheckImageArgs16NotInPlace(src, dst);
		uint16_t* __restrict ps = (uint16_t*)src.getData();
		uint16_t* __restrict pd = (uint16_t*)dst.getData();
		int np = src.getHeight() * src.getPitch() >> 4; 
		short offset = (short)(1 << 15);
		short signedOffset = (short)threshold + offset;
		__m128i t = _mm_set1_epi16(signedOffset);
		__m128i voffset = _mm_set1_epi16(offset);
		int idx = 0;
		if (isLte)
		{
			for (int i = 0; i < np; i++)
			{
				__m128i v1 = _mm_loadu_si128((__m128i*)(ps + idx));
				__m128i v2 = _mm_add_epi16(v1, voffset);
				__m128i mask = _mm_cmpgt_epi16(v2, t);
				__m128i res = _mm_and_si128(v1, mask);
				_mm_storeu_si128((__m128i*)(pd + idx), res);
				idx += 8;
			}
		}
		else
		{
			for (int i = 0; i < np; i++)
			{
				__m128i v1 = _mm_loadu_si128((__m128i*)(ps + idx)); 
				__m128i v2 = _mm_add_epi16(v1, voffset);
				__m128i mask = _mm_cmplt_epi16(v2, t);
				__m128i res = _mm_and_si128(v1, mask);
				_mm_storeu_si128((__m128i*)(pd + idx), res);
				idx += 8;
			}
		}
	}
	static void Threshold2(Img& src, uint16_t threshold, bool isLte, Img& dst)
	{
		int w = src.getWidth();
		int h = src.getHeight();
		int spp = src.getPitch() / 2;
		int dpp = dst.getPitch();
		uint16_t* ps = (uint16_t*)src.getData();
		uint8_t* pd = (uint8_t*)dst.getData();
		if (isLte)
		{
			for (int y = 0; y < h; y++)
			{
				for (int x = 0; x < w; x++)
				{
					pd[x] = (byte)((ps[x] <= threshold) ? 255 : 0);
				}
				ps += spp;
				pd += dpp;
			}
		}
		else
		{
			for (int y = 0; y < h; y++)
			{
				for (int x = 0; x < w; x++)
				{
					pd[x] = (byte)((ps[x] >= threshold) ? 255 : 0);
				}
				ps += spp;
				pd += dpp;
			}
		}
	}
	static void ThresholdRange(Img& src, int low, int high, bool inRange, Img& dst)
	{
		int w = src.getWidth();
		int h = src.getHeight();
		int p = src.getPitch() / 2;
		int dp = dst.getPitch();
		uint16_t* ps = (uint16_t*)src.getData();
		uint8_t* pd = (uint8_t*)dst.getData();
		if (inRange)
		{
			for (int y = 0; y < h; y++)
			{
				for (int x = 0; x < w; x++)
				{
					pd[x] = (byte)(((ps[x] >= low) && (ps[x] < high)) ? 255 : 0);
				}
				pd += dp;
				ps += p;
			}
		}
		else
		{
			for (int y = 0; y < h; y++)
			{
				for (int x = 0; x < w; x++)
				{
					pd[x] = (byte)(((ps[x] < low) || (ps[x] >= high)) ? 255 : 0);
				}
				pd += dp;
				ps += p;
			}
		}
	}
	static void ThresholdRangeSse2(Img& src, int low, int high, bool inRange, Img& dst)
	{
		int h = src.getHeight();
		int spp = src.getPitch() / 2;
		int dpp = dst.getPitch();
		uint16_t* ps = (uint16_t*)src.getData();
		uint8_t* pd = (uint8_t*)dst.getData();
		int np = (spp >> 4); 
		short offset = (short)(1 << 15);
		short sLow = (short)low + offset;
		short sHigh = (short)high + offset;
		__m128i mlow = _mm_set1_epi16(sLow);
		__m128i mhigh = _mm_set1_epi16(sHigh);
		__m128i voffset = _mm_set1_epi16(offset);
		if (inRange)
		{
			for (int y = 0; y < h; y++)
			{
				for (int i = 0; i < np; i++) 
				{
					__m128i v1 = _mm_load_si128((__m128i*)(ps + (i << 4)));
					__m128i v2 = _mm_load_si128((__m128i*)(ps + (i << 4) + 8));
					v1 = _mm_add_epi16(v1, voffset);
					v2 = _mm_add_epi16(v2, voffset);
					__m128i lt1 = _mm_cmplt_epi16(v1, mhigh);
					__m128i gt1 = _mm_cmpgt_epi16(v1, mlow);
					__m128i res1 = _mm_and_si128(gt1, lt1);
					__m128i lt2 = _mm_cmplt_epi16(v2, mhigh);
					__m128i gt2 = _mm_cmpgt_epi16(v2, mlow);
					__m128i res2 = _mm_and_si128(gt2, lt2);
					__m128i res = _mm_packs_epi16(res1, res2);
					_mm_storeu_si128((__m128i*)(pd + (i << 4)), res);
				}
				pd += dpp;
				ps += spp;
			}
		}
		else
		{
			for (int y = 0; y < h; y++)
			{
				for (int i = 0; i < np; i++) 
				{
					__m128i v1 = _mm_load_si128((__m128i*)(ps + (i << 4)));
					__m128i v2 = _mm_load_si128((__m128i*)(ps + (i << 4) + 8));
					v1 = _mm_add_epi16(v1, voffset);
					v2 = _mm_add_epi16(v2, voffset);
					__m128i gt1 = _mm_cmpgt_epi16(v1, mhigh);
					__m128i lt1 = _mm_cmplt_epi16(v1, mlow);
					__m128i res1 = _mm_or_si128(gt1, lt1);
					__m128i gt2 = _mm_cmpgt_epi16(v2, mhigh);
					__m128i lt2 = _mm_cmplt_epi16(v2, mlow);
					__m128i res2 = _mm_or_si128(gt2, lt2);
					__m128i res = _mm_packs_epi16(res1, res2);
					_mm_storeu_si128((__m128i*)(pd + (i << 4)), res);
				}
				pd += dpp;
				ps += spp;
			}
		}
	}
	static void ColProfile(Img& src, vector<int>& counts, vector<uint16_t>& means, bool debug = false)
	{
		int w = src.getWidth();
		int h = src.getHeight();
		int spp = src.getPitch() / 2;
		int nc = w;
		vector<int> sums(nc);
		counts.resize(nc);
		means.resize(nc);
		uint16_t* ps = (uint16_t*)src.getData();
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				uint16_t v = ps[x];
				if (v > 0)
				{
					counts[x]++;
					sums[x] += v;
				}
			}
			ps += spp;
		}
		if (debug) logd.printf("Col profile: ");
		for (int i = 0; i < nc; i++) 
		{
			if (counts[i] > 0)
			{
				means[i] = (uint16_t)(sums[i] / counts[i]);
			}
			if (debug) logd.printf("%d,", means[i]);
		}
		if (debug) logd.printf("\n");
	}
	static void RowProfile(Img& src, vector<int>& counts, vector<uint16_t>& means, bool debug = false)
	{
		int w = src.getWidth();
		int h = src.getHeight();
		int spp = src.getPitch() / 2;
		int nc = h;
		vector<int> sums(nc);
		counts.resize(nc);
		means.resize(nc);
		uint16_t* ps = (uint16_t*)src.getData();
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				uint16_t v = ps[x];
				if (v > 0)
				{
					counts[y]++;
					sums[y] += v;
				}
			}
			ps += spp;
		}
		if (debug) logd.printf("Row profile: ");
		for (int i = 0; i < nc; i++)
		{
			if (counts[i] > 0)
			{
				means[i] = (uint16_t)(sums[i] / counts[i]);
			}
			if (debug) logd.printf("%d,", means[i]);
		}
		if (debug) logd.printf("\n");
	}
	static bool ZeroColsByProfile(Img& src, float lowNstddevThreshold, float highNstddevThreshold, bool debug)
	{
		int w = src.getWidth();
		int debugColIndex = -1; 
		vector<int> cpCounts;
		vector<uint16_t> cpMeans;
		ColProfile(src, cpCounts, cpMeans, false);
		vector<float> goodMeansf;
		for (int i = 0; i < cpMeans.size(); i++) if (cpCounts[i] > 0) goodMeansf.push_back((float)cpMeans[i]);
		vector<int> ci; 
		float lowThreshold = 0, highThreshold = 0;
		if (goodMeansf.size() > w/2)
		{
			float min, max, mean, stddev, meandev;
			Stats::ComputeStats(goodMeansf, min, max, mean, stddev, meandev);
			int n2 = goodMeansf.size();
			vector<float> gm2;
			gm2.reserve(n2);
			float t1 = mean - 3 * stddev;
			float t2 = mean + 3 * stddev;
			for (int i = 0; i < goodMeansf.size(); i++) 
				if ((goodMeansf[i] > t1) && (goodMeansf[i] < t2)) 
					gm2.push_back(goodMeansf[i]);
			gm2.shrink_to_fit();
			Stats::ComputeStats(gm2, min, max, mean, stddev, meandev);
			lowThreshold = roundf(mean - stddev * lowNstddevThreshold);
			highThreshold = roundf(mean + stddev * highNstddevThreshold);
			if (debug) logd.printf("ZeroColsByProfile: thresholds %.1f and %.1f from min/mean/max: %.1f / %.1f / %.1f, stddev: %.1f, meandev: %.1f, nstddev: %.1f/%.1f\n", 
				lowThreshold, highThreshold, min, mean, max, stddev, meandev, lowNstddevThreshold, highNstddevThreshold);
			for (int i = 0; i < w; i++)
			{
				if (debug && (debugColIndex == i))
				{
					logd.printf("ZeroColsByProfile: Col %d: count: %d, mean: %d\n", i, cpCounts[i], cpMeans[i]);
				}
				if ((cpCounts[i] > 0) && ((cpMeans[i] < lowThreshold) || (cpMeans[i] > highThreshold)))
				{
					ci.push_back(i);
				}
			}
			if (ci.size() > w / 32)
			{
				logd.warn("ZeroColsByProfile: Selected %d cols to zero based on thresholds %.1f to %.1f from mean %f, stddev: %f, nstddev: %f/%f",
					ci.size(), lowThreshold, highThreshold, mean, stddev, lowNstddevThreshold, highNstddevThreshold);
				logd.warn("ZeroColsByProfile: Col profile: ");
				for (int i = 0; i < MIN(cpMeans.size(), 200); i++)
				{
					logd.printf("%d, ", cpMeans[i]);
				}
				logd.printf("\n");
				return false;
			}
			for (int i = 0; i < ci.size(); i++)
			{
				ZeroCol(src, ci[i]);
				if (debug) logd.printf("ZeroColsByProfile: Zero col %d on mean %d vs thresholds %.1f to %.1f\n", ci[i], cpMeans[ci[i]], lowThreshold, highThreshold);
			}
			return true;
		}
		else
		{
			logd.warn("ZeroColsByProfile: Only %d good cols, not enough to make decisions on.", goodMeansf.size());
			return false;
		}
	}
	static bool ZeroRowsByProfile(Img& src, float lowNstddevThreshold, float highNstddevThreshold, bool debug)
	{
		int h = src.getHeight();
		int debugRowIndex = -1; 
		vector<int> cpCounts;
		vector<uint16_t> cpMeans;
		RowProfile(src, cpCounts, cpMeans, false);
		vector<float> goodMeansf;
		for (int i = 0; i < cpMeans.size(); i++) if (cpCounts[i] > 0) goodMeansf.push_back((float)cpMeans[i]);
		vector<int> ci; 
		float lowThreshold = 0, highThreshold = 0;
		if (goodMeansf.size() > h / 2)
		{
			float min, max, mean, stddev, meandev;
			Stats::ComputeStats(goodMeansf, min, max, mean, stddev, meandev);
			int n2 = goodMeansf.size();
			vector<float> gm2;
			gm2.reserve(n2);
			float t1 = mean - 3 * stddev;
			float t2 = mean + 3 * stddev;
			for (int i = 0; i < goodMeansf.size(); i++)
				if ((goodMeansf[i] > t1) && (goodMeansf[i] < t2))
					gm2.push_back(goodMeansf[i]);
			gm2.shrink_to_fit();
			Stats::ComputeStats(gm2, min, max, mean, stddev, meandev);
			lowThreshold = roundf(mean - stddev * lowNstddevThreshold);
			highThreshold = roundf(mean + stddev * highNstddevThreshold);
			if (debug) logd.printf("ZeroRowsByProfile: thresholds %.1f and %.1f from min/mean/max: %.1f / %.1f / %.1f, stddev: %.1f, meandev: %.1f, nstddev: %.1f/%.1f\n",
				lowThreshold, highThreshold, min, mean, max, stddev, meandev, lowNstddevThreshold, highNstddevThreshold);
			for (int i = 0; i < h; i++)
			{
				if (debug && (debugRowIndex == i))
				{
					logd.printf("ZeroRowsByProfile: Row %d: count: %d, mean: %d\n", i, cpCounts[i], cpMeans[i]);
				}
				if ((cpCounts[i] > 0) && ((cpMeans[i] < lowThreshold) || (cpMeans[i] > highThreshold)))
				{
					ci.push_back(i);
				}
			}
			if (ci.size() > h / 32)
			{
				logd.warn("ZeroRowsByProfile: Selected %d rows to zero based on thresholds %.1f to %.1f from mean %f, stddev: %f, nstddev: %f/%f",
					ci.size(), lowThreshold, highThreshold, mean, stddev, lowNstddevThreshold, highNstddevThreshold);
				logd.warn("ZeroRowsByProfile: Row profile: ");
				for (int i = 0; i < MIN(cpMeans.size(), 200); i++)
				{
					logd.printf("%d, ", cpMeans[i]);
				}
				logd.printf("\n");
				return false;
			}
			for (int i = 0; i < ci.size(); i++)
			{
				ZeroRow(src, ci[i]);
				if (debug) logd.printf("ZeroRowsByProfile: Zero row %d on mean %d vs thresholds %.1f to %.1f\n", ci[i], cpMeans[ci[i]], lowThreshold, highThreshold);
			}
			return true;
		}
		else
		{
			logd.warn("ZeroRowsByProfile: Only %d good rows, not enough to make decisions on.", goodMeansf.size());
			return false;
		}
	}
	static void Convert16sTo16u(Img& src, Img& dst)
	{
		CheckImageArgs16InPlace(src, dst);
		int16_t* ps = (int16_t*)src.getData();
		uint16_t* pd = (uint16_t*)dst.getData();
		int h = src.getHeight();
		int pp = src.getPitch() / 2;
		int n = pp * h; 
		for (int i = 0; i < n; i++)
		{
			short v = ps[i] - 32768;
			pd[i] = (unsigned short)v;
		}
	}
	static void HistCompute(Img& src, vector<int>& bins)
	{
		int nb = bins.size();
		int shift = (int)round(log2(nb));
		if (bins.size() != (1 << shift)) ErrorExit("Bins size must be power of 2, not %d (shift %d).", bins.size(), shift);
		shift = 16 - shift;
		for (int i = 0; i < nb; i++) bins[i] = 0;
		uint16_t* ps = (uint16_t*)src.getData();
		int w = src.getWidth();
		int h = src.getHeight();
		int p = src.getPitch() / 2;
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				bins[ps[x] >> shift]++;
			}
			ps += p;
		}
	}
	static void HistCompute(Img& src, int nBins, Hist& outputHist)
	{
		vector<int> bins(nBins);
		int shift = (int)round(log2(nBins));
		if (bins.size() != (1 << shift)) ErrorExit("Bins size must be power of 2, not %d (shift %d).", bins.size(), shift);
		shift = 16 - shift;
		for (int i = 0; i < nBins; i++) bins[i] = 0;
		vector<int> levels(nBins + 1);
		levels[0] = 0;
		for (int i = 1; i <= nBins; i++) levels[i] = i << shift;
		uint16_t* ps = (uint16_t*)src.getData();
		int w = src.getWidth();
		int h = src.getHeight();
		int p = src.getPitch() / 2;
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				bins[ps[x] >> shift]++;
			}
			ps += p;
		}
		outputHist.SetData(bins, levels);
	}
	static void HistCompute(Img& src, Img& mask, int nBins, Hist& outputHist, bool isGoodMask = false, bool debug = false)
	{
		vector<int> bins(nBins);
		int shift = (int)round(log2(nBins));
		if (bins.size() != (1 << shift)) ErrorExit("Bins size must be power of 2, not %d (shift %d).", bins.size(), shift);
		shift = 16 - shift;
		for (int i = 0; i < nBins; i++) bins[i] = 0;
		vector<int> levels(nBins + 1);
		levels[0] = 0;
		for (int i = 1; i <= nBins; i++) levels[i] = i << shift;
		uint16_t* ps = (uint16_t*)src.getData();
		uint8_t* pm = (uint8_t*)mask.getData();
		int w = src.getWidth();
		int h = src.getHeight();
		int p = src.getPitch() / 2;
		int mpp = mask.getPitch();
		assert(w == mask.getWidth());
		assert(h == mask.getHeight());
		assert(src.getDepthInBits() == 16);
		assert(mask.getDepthInBits() == 8);
		if (isGoodMask)
		{
			for (int y = 0; y < h; y++)
			{
				for (int x = 0; x < w; x++)
				{
					int idx = ps[x] >> shift;
					if (debug && (idx == 0) && (1 & pm[x])) logd.printf("HistCompute: Index 0 at %d,%d\n", x, y);
					bins[idx] = bins[idx] + (1 & pm[x]);
				}
				ps += p;
				pm += mpp;
			}
		}
		else
		{
			for (int y = 0; y < h; y++)
			{
				for (int x = 0; x < w; x++)
				{
					int idx = ps[x] >> shift;
					if (debug && (idx == 0) && (1 & ~pm[x])) logd.printf("HistCompute: Index 0 at %d,%d\n", x, y);
					bins[idx] = bins[idx] + (1 & ~pm[x]);
				}
				ps += p;
				pm += mpp;
			}
		}
		outputHist.SetData(bins, levels);
	}
	static void HistCustomLevels(Img& src, vector<int>& levels, vector<int>& bins)
	{
		int nb = bins.size();
		for (int i = 0; i < nb; i++) bins[i] = 0;
		assert(levels.size() == nb + 1);
		uint16_t* ps = (uint16_t*)src.getData();
		int w = src.getWidth();
		int h = src.getHeight();
		int pp = src.getPitch() / 2;
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				uint16_t v = ps[x];
				for (int b = 1; b < nb + 1; b++)
				{
					if (v < levels[b])
					{
						bins[b - 1]++;
						break;
					}
				}
			}
			ps += pp;
		}
	}
	static void Downsample(Img& src, Img& dst, bool doMean2x2Filter = false)
	{
		int width = src.getWidth();
		int height = src.getHeight();
		assert(width / 2 <= dst.getWidth());
		assert(height / 2 <= dst.getHeight());
		assert(src.getDepthInBits() == dst.getDepthInBits());
		assert(src.getDepthInBits() == 16);
		assert(src.getData() != dst.getData());
		int snp = src.getPitch() >> 1; 
		int sStride = snp * 2;
		int dnp = dst.getPitch() >> 1;
		uint16_t* __restrict ps = (uint16_t*)src.getData();
		uint16_t* __restrict pd = (uint16_t*)dst.getData();
		if (!doMean2x2Filter)
		{
			for (int y = 0; y < height; y += 2)
			{
				int idx = 0;
				for (int i = 0; i < width; i += 2)
				{
					pd[idx++] = ps[i];
				}
				ps += sStride;
				pd += dnp;
			}
		}
		else
		{
			for (int y = 0; y < height; y += 2)
			{
				int idx = 0;
				for (int i = 0; i < width; i += 2)
				{
					int val = ps[i] + ps[i + 1] + ps[i + snp] + ps[i + snp + 1];
					pd[idx++] = val >> 2;
				}
				ps += sStride;
				pd += dnp;
			}
		}
	}
	static void DownsampleSse2(Img& src, Img& dst, bool doMean2x2Filter = false)
	{
		int width = src.getWidth();
		int height = src.getHeight();
		assert(width / 2 <= dst.getWidth());
		assert(height / 2 <= dst.getHeight());
		assert(src.getDepthInBits() == dst.getDepthInBits());
		assert(src.getDepthInBits() == 16);
		assert(src.getData() != dst.getData());
		int sp = src.getPitch() >> 1; 
		int sStride = sp * 2; 
		int dnp = dst.getPitch() >> 1; 
		__m128i mask = _mm_set1_epi32(0x0000FFFF);
		uint16_t* __restrict ps = (uint16_t*)src.getData();
		uint16_t* __restrict pd = (uint16_t*)dst.getData();
		if (!doMean2x2Filter)
		{
			for (int y = 0; y < height; y += 2)
			{
				int i;
				for (i = 0; i < width - 15; i += 16)
				{
					__m128i v1 = _mm_load_si128((__m128i*)(ps + i));
					__m128i v2 = _mm_load_si128((__m128i*)(ps + i + 8));
					v1 = _mm_and_si128(v1, mask);
					v2 = _mm_and_si128(v2, mask);
					__m128i res = _mm_packs_epi32(v1, v2); 
					_mm_storeu_si128((__m128i*)(pd + (i >> 1)), res);
				}
				for (; i < width; i += 2)
				{
					pd[i >> 1] = ps[i];
				}
				ps += sStride;
				pd += dnp;
			}
		}
		else
		{
			for (int y = 0; y < height - 1; y += 2)
			{
				for (int i = 0; i < width - 15; i += 16)
				{
					__m128i v1 = _mm_load_si128((__m128i*)(ps + i));
					__m128i v2 = _mm_load_si128((__m128i*)(ps + i + 8));
					__m128i v3 = _mm_load_si128((__m128i*)(ps + i + sp));
					__m128i v4 = _mm_load_si128((__m128i*)(ps + i + 8 + sp));
					__m128i v1s = _mm_avg_epu16(v1, v3);
					__m128i v2s = _mm_avg_epu16(v2, v4);
					v1 = _mm_srli_epi32(v1s, 16);
					v1s = _mm_avg_epu16(v1, v1s);
					v1s = _mm_and_si128(v1s, mask);
					v2 = _mm_srli_epi32(v2s, 16);
					v2s = _mm_avg_epu16(v2, v2s);
					v2s = _mm_and_si128(v2s, mask);
					__m128i res = _mm_packs_epi32(v1s, v2s);
					_mm_storeu_si128((__m128i*)(pd + (i >> 1)), res);
				}
				ps += sStride;
				pd += dnp;
			}
		}
	}
	private:
	static void upSampleSubSse2(Img& src, IRoi& srcRoi, Img& dst)
	{
		assert(srcRoi.x0 == 0);
		assert(srcRoi.y0 == 0);
		assert(dst.getWidth() >= srcRoi.getWidth() * 2);
		assert(dst.getHeight() >= srcRoi.getHeight() * 2);
		int w = srcRoi.getWidth();
		int h = srcRoi.getHeight();
		int srcPitch = src.getPitch() >> 1;
		int dstPitch = dst.getPitch() >> 1;
		int nb = w / 8;
		int nTail = w - (nb * 8); 
		int srcExtra = srcPitch - w;
		int dstExtra = dstPitch - w * 2;
		uint16_t* srcp = (uint16_t*)src.getData(srcRoi.x0, srcRoi.y0);
		uint16_t* dstp1 = (uint16_t*)dst.getData();
		uint16_t* dstp2 = dstp1 + dstPitch;
		for (int y = 0; y < h; y++)
		{
			for (int i = 0; i < nb; i++)
			{
				__m128i a = _mm_loadu_si128((const __m128i*)srcp);
				__m128i t1 = _mm_unpacklo_epi16(a, a); 
				_mm_storeu_si128((__m128i*)dstp1, t1); 
				_mm_storeu_si128((__m128i*)dstp1, t1); 
				_mm_storeu_si128((__m128i*)dstp2, t1);
				dstp1 += 8;
				dstp2 += 8;
				t1 = _mm_unpackhi_epi16(a, a); 
				_mm_storeu_si128((__m128i*)dstp1, t1); 
				_mm_storeu_si128((__m128i*)dstp2, t1);
				dstp1 += 8;
				dstp2 += 8;
				srcp += 8;
			}
			for (int i = 0; i < nTail; i++)
			{
				uint16_t val = *srcp;
				dstp1[0] = val;
				dstp1[1] = val;
				dstp2[0] = val;
				dstp2[1] = val;
				srcp++;
				dstp1 += 2;
				dstp2 += 2;
			}
			srcp += srcExtra;
			dstp1 += dstExtra + dstPitch;
			dstp2 += dstExtra + dstPitch;
		}
	}
	public:
	static void Upsample(Img& src, Img& dst, bool doInterpolate = false)
	{
		assert(!doInterpolate && "interpolation not implemented in upsample yet");
		int width = src.getWidth();
		int height = src.getHeight();
		IRoi roi(0, 0, width, height);
		upSampleSubSse2(src, roi, dst);
		int dw = dst.getWidth();
		int dh = dst.getHeight();
		if (dw > width * 2)
		{
			for (int y = 0; y < dh; y++)
			{
				uint16_t val = dst.getPixel16u(width * 2 - 1, y);
				for (int x = width * 2; x < dw; x++)
				{
					dst.setPixel16u(x, y, val);
				}
			}
		}
		if (dh > height * 2)
		{
			for (int y = height * 2; y < dh; y++)
			{
				for (int x = 0; x < dw; x++)
				{
					uint16_t val = dst.getPixel16u(x, height * 2 - 1);
					dst.setPixel16u(x, y, val);
				}
			}
		}
	}
	static Stats ComputeStats(Img& src, Img& mask, bool isGoodMask = false)
	{
		int w = src.getWidth();
		int h = src.getHeight();
		int spp = src.getPitch() / 2;
		int mpp = mask.getPitch();
		uint16_t* ps = (uint16_t*)src.getData();
		uint8_t* pm = (uint8_t*)mask.getData();
		assert(w == mask.getWidth());
		assert(h == mask.getHeight());
		assert(src.getDepthInBits() == 16);
		assert(mask.getDepthInBits() == 8);
		uint16_t mn = 65535;
		uint16_t mx = 0;
		double sum = 0;
		double ssum = 0;
		int n = 0;
		if (!isGoodMask)
		{
			for (int y = 0; y < h; y++)
			{
				for (int x = 0; x < w; x++)
				{
					uint16_t m = (pm[x] | (pm[x] << 8)); 
					uint16_t v = ~m & ps[x];
					sum += v;
					ssum += ((double)v * v);
					n += (1 & ~m);
					mn = min(mn, (uint16_t)(v | m));
					mx = max(mx, (uint16_t)(v & ~m));
				}
				ps += spp;
				pm += mpp;
			}
		}
		else
		{
			for (int y = 0; y < h; y++)
			{
				for (int x = 0; x < w; x++)
				{
					uint16_t m = (pm[x] | (pm[x] << 8)); 
					uint16_t v = m & ps[x];
					sum += v;
					ssum += ((double)v * v);
					n += (1 & m);
					mn = min(mn, (uint16_t)(v | ~m));
					mx = max(mx, (uint16_t)(v & m));
				}
				ps += spp;
				pm += mpp;
			}
		}
		Stats stats;
		stats.max = mx;
		stats.min = mn;
		stats.mean = (float)(sum / n);
		stats.n = n;
		double t1 = ((n * ssum) - sum * sum);
		if (t1 > 0)
		{
			stats.stddev = (float)sqrt(t1 / ((double)n * (n - 1)));
		}
		else
		{
			stats.stddev = 0;
		}
		ps = (uint16_t*)src.getData();
		pm = (uint8_t*)mask.getData();
		sum = 0;
		if (!isGoodMask)
		{
			for (int y = 0; y < h; y++)
			{
				for (int x = 0; x < w; x++)
				{
					uint16_t dv = (uint16_t)abs(ps[x] - stats.mean);
					uint16_t m = (pm[x] | (pm[x] << 8)); 
					dv = ~m & dv;
					sum += dv;
				}
				ps += spp;
				pm += mpp;
			}
		}
		else
		{
			for (int y = 0; y < h; y++)
			{
				for (int x = 0; x < w; x++)
				{
					uint16_t dv = (uint16_t)abs(ps[x] - stats.mean);
					uint16_t m = (pm[x] | (pm[x] << 8)); 
					dv = m & dv;
					sum += dv;
				}
				ps += spp;
				pm += mpp;
			}
		}
		stats.meandev = sum / n;
		return stats;
	}
	static void ComputeStats(Img& src, Img& mask, bool isGoodMask, Stats& stats)
	{
		vector<uint16_t> values;
		GetGoodValuesMasked(src, mask, isGoodMask, values);
		stats.Compute(values);
	}
	static void Unwrap(Img& src, Stats stats, float nmeandev)
	{
		int w = src.getWidth();
		int h = src.getHeight();
		int spp = src.getPitch() / 2;
		uint16_t* ps = (uint16_t*)src.getData();
		int lhx = -1; 
		uint16_t mean = (uint16_t)stats.mean;
		uint16_t highThreshold = (uint16_t)(stats.mean + nmeandev * stats.meandev);
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				uint16_t v = ps[x];
				if (v > highThreshold)
				{
					lhx = x;
				}
				else if ((v < mean) && (lhx >= 0))
				{
					ps[x] = ps[lhx];
				}
				else
				{
					lhx = -1;
				}
			}
			ps += spp;
		}
	}
	static void MaskWrappedValues(Img& src, Stats stats, float nmeandev, Img& mask)
	{
		int w = src.getWidth();
		int h = src.getHeight();
		int spp = src.getPitch() / 2;
		int mpp = mask.getPitch();
		uint16_t* ps = (uint16_t*)src.getData();
		uint8_t* pm = (uint8_t*)mask.getData();
		assert(w == mask.getWidth());
		assert(h == mask.getHeight());
		assert(src.getDepthInBits() == 16);
		assert(mask.getDepthInBits() == 8);
		uint16_t mean = (uint16_t)stats.mean;
		uint16_t highThreshold = (uint16_t)(stats.mean + nmeandev * stats.meandev);
		const int MaxRunLength = 10; 
		int lhx = -1; 
		int nlow = 0;
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				uint16_t v = ps[x];
				if (v > highThreshold)
				{
					if ((nlow > 0) && (nlow < MaxRunLength))
					{
						for (int xx = x - nlow; xx < x; xx++)
						{
							pm[xx] = (uint8_t)0xFF;
						}
						nlow = 0;
					}
					else
					{
						lhx = x;
					}
				}
				else if ((v < mean) && (lhx >= 0))
				{
					nlow++;
				}
				else
				{
					lhx = -1;
					nlow = 0;
				}
			}
			ps += spp;
			pm += mpp;
		}
	}
	static void ApplyMask(Img& src, Img& mask, bool isBadMask = true)
	{
		CheckImageArgsMask(src, mask);
		int w = src.getWidth();
		int h = src.getHeight();
		int spp = src.getPitch() / 2;
		int mpp = mask.getPitch();
		uint16_t* ps = (uint16_t*)src.getData();
		uint8_t* pm = (uint8_t*)mask.getData();
		assert(w == mask.getWidth());
		assert(h == mask.getHeight());
		assert(src.getDepthInBits() == 16);
		assert(mask.getDepthInBits() == 8);
		if (isBadMask)
		{
			for (int y = 0; y < h; y++)
			{
				for (int x = 0; x < w; x++)
				{
					if (pm[x] > 0) ps[x] = 0;
				}
				ps += spp;
				pm += mpp;
			}
		}
		else
		{
			for (int y = 0; y < h; y++)
			{
				for (int x = 0; x < w; x++)
				{
					if (pm[x] == 0) ps[x] = 0;
				}
				ps += spp;
				pm += mpp;
			}
		}
	}
	static void FillUnderMask(Img& src, Img& posMask, bool interpolate = false)
	{
		int w = src.getWidth();
		int h = src.getHeight();
		int spp = src.getPitch() / 2;
		int mpp = posMask.getPitch();
		uint16_t* ps = (uint16_t*)src.getData();
		uint8_t* pm = (uint8_t*)posMask.getData();
		assert(w == posMask.getWidth());
		assert(h == posMask.getHeight());
		assert(src.getDepthInBits() == 16);
		assert(posMask.getDepthInBits() == 8);
		if (!interpolate)
		{
			for (int y = 0; y < h; y++)
			{
				int lgx = -1; 
				for (int x = 0; x < w; x++)
				{
					if (pm[x] > 0)
					{
						if (lgx > 0)
						{
							ps[x] = ps[lgx];
						}
					}
					else
					{
						if (lgx < 0)
						{
							for (int j = x - 1; j >= 0; j--)
							{
								ps[j] = ps[x];
							}
						}
						lgx = x;
					}
				}
				ps += spp;
				pm += mpp;
			}
		}
		else
		{
			for (int y = 0; y < h; y++)
			{
				int lgx = -1; 
				for (int x = 0; x < w; x++)
				{
					if (pm[x] > 0)
					{
					}
					else
					{
						if (lgx < 0)
						{
							for (int j = x - 1; j >= 0; j--)
							{
								ps[j] = ps[x];
							}
						}
						else
						{
							if ((lgx >= 0) && (x > lgx + 1))
							{
								float len = x - lgx;
								int dv = ps[x] - ps[lgx];
								for (int i = lgx + 1; i < x; i++)
								{
									float frac = (i - lgx) / len;
									uint16_t val = (uint16_t)(ps[lgx] + roundf(frac * dv));
									ps[i] = val;
								}
							}
						}
						lgx = x;
					}
				}
				if ((lgx < w - 1) && (lgx >= 0))
				{
					for (int i = lgx + 1; i < w; i++)
					{
						ps[i] = ps[lgx];
					}
				}
				ps += spp;
				pm += mpp;
			}
		}
	}
	static void ShiftInt(Img& src, pair<int, int> offset, int value)
	{
		src.Shift(offset.first, offset.second, value);
	}
	static void DrawBox(Img& img, int x0, int y0, int x1, int y1, uint16_t value)
	{
		int w = img.getWidth();
		int h = img.getHeight();
		int pb = img.getPitch(); 
		int pp = pb >> 1; 
		uint16_t* data = (uint16_t*)img.getData();
		x0 = CLIP(x0, 0, w - 1);
		x1 = CLIP(x1, 0, w);
		y0 = CLIP(y0, 0, h - 1);
		y1 = CLIP(y1, 0, h);
		for (int i = x0 + y0 * pp; i < x1 + y0 * pp; i++)
		{
			data[i] = value;
		}
		for (int y = y0 + 1; y < y1; y++)
		{
			data[x0 + y * pp] = value;
			data[x1 - 1 + y * pp] = value;
		}
		for (int i = x0 + (y1 - 1) * pp; i < x1 + (y1 - 1) * pp; i++)
		{
			data[i] = value;
		}
	}
	static void Flood(Img& src, bool is8Connected, int cx, int cy, uint16_t val, Img& dst, bool clearDest = true)
	{
		CheckImageArgsSize(src, dst);
		assert(src.getDepthInBits() == 8);
		assert(dst.getDepthInBits() == 16);
		int w = src.getWidth();
		int offsets4[4] = {-1, -w, w, 1};
		int offsets8[8] = {-w - 1, -w, -w + 1, -1, 1, w - 1, w, w + 1};
		int nn = 4;
		int* offsets = offsets4;
		if (is8Connected)
		{
			nn = 8;
			offsets = offsets8;
		}
		if (clearDest)
		{
			Clear(dst);
		}
		if (src.getPixel(cx, cy) > 0)
		{
			queue<int> qi;
			qi.push(cx + cy * w);
			dst.setPixel16u(cx, cy, val);
			while (qi.size() > 0)
			{
				int ci = qi.front();
				qi.pop();
				for (int i = 0; i < nn; i++)
				{
					int idx = ci + offsets[i];
					if (src.checkIsInImage(idx) && (src.getPixel(idx) > 0) && (dst.getPixel16u(idx) == 0))
					{
						dst.setPixel16u(idx, val);
						qi.push(idx);
					}
				}
			}
		}
	}
	static bool Label(Img& src, Img& dst, bool is8Connected)
	{
		CheckImageArgsSize(src, dst);
		assert(src.getDepthInBits() == 8);
		assert(dst.getDepthInBits() == 16);
		int w = src.getWidth();
		int h = src.getHeight();
		Clear(dst);
		int label = 1;
		for (int y = 1; y < h - 1; y++)
		{
			for (int x = 1; x < w - 1; x++)
			{
				byte vs = src.getPixel(x, y);
				uint16_t vd = dst.getPixel16u(x, y);
				if ((vs > 0) && (vd == 0))
				{
					Flood(src, is8Connected, x, y, label++, dst, false);
					if (label > 65535)
					{
						return false;
					}
				}
			}
		}
		return true;
	}
	static void FindBlobs(Img& image, vector<ImgBlob>& blobs)
	{
		assert(image.getDepthInBits() == 16);
		int w = image.getWidth();
		int h = image.getHeight();
		int spp = image.getPitch() >> 1;
		uint16_t* ps = (uint16_t*)image.getData();
		vector<float> xsums(65535);
		vector<float> ysums(65535);
		vector<int> counts(65535);
		vector<int> xmin(65535);
		vector<int> xmax(65535);
		vector<int> ymin(65535);
		vector<int> ymax(65535);
		for (int i = 0; i < 65535; i++)
		{
			xmin[i] = 999999;
			xmax[i] = 0;
			ymin[i] = 999999;
			ymax[i] = 0;
		}
		vector<int> xs(65535);
		vector<int> ys(65535);
		for (int y = 1; y < h - 1; y++)
		{
			int si = y * spp;
			for (int x = 1; x < w - 1; x++)
			{
				uint16_t vs = ps[si + x];
				if (vs > 0)
				{
					counts[vs]++;
					xsums[vs] += x;
					ysums[vs] += y;
					xs[vs] = x;
					ys[vs] = y;
					xmin[vs] = MIN(xmin[vs], x);
					xmax[vs] = MAX(xmax[vs], x);
					ymin[vs] = MIN(ymin[vs], y);
					ymax[vs] = MAX(ymax[vs], y);
				}
			}
		}
		int n = 0;
		for (int i = 1; i < 65535; i++) if (counts[i] > 0) n++;
		blobs.resize(n);
		int j = 0;
		for (int i = 1; i < 65535; i++)
		{
			if (counts[i] > 0)
			{
				ImgBlob blob;
				blob.area = counts[i];
				blob.label = i;
				blob.xCentroid = xsums[i] / counts[i];
				blob.yCentroid = ysums[i] / counts[i];
				if (image.getPixel(roundf(blob.xCentroid), roundf(blob.yCentroid)) == i)
				{
					blob.cx = roundf(blob.xCentroid);
					blob.cy = roundf(blob.yCentroid);
				}
				else
				{
					blob.cx = xs[i];
					blob.cy = ys[i];
				}
				blob.roi.x0 = xmin[i];
				blob.roi.x1 = xmax[i] + 1;
				blob.roi.y0 = ymin[i];
				blob.roi.y1 = ymax[i] + 1;
				blobs[j++] = blob;
			}
		}
	}
	static void CopyBlob(Img& src, Img& label, IRoi& roi, uint16_t labelValue, Img& dst, bool sameLocationInDst = false)
	{
		CheckImageArgsSize(src, label);
		assert(roi.getWidth() <= dst.getWidth()); 
		assert(roi.getHeight() <= dst.getHeight()); 
		assert(src.getDepthInBits() == 8);
		assert(label.getDepthInBits() == 16);
		assert(dst.getDepthInBits() == 8);
		int w = roi.getWidth();
		int h = roi.getHeight();
		int sp = src.getPitch();
		int lp = label.getPitch() >> 1;
		int dp = dst.getPitch();
		byte* ps = (byte*)src.getData(roi.x0, roi.y0);
		uint16_t* pl = label.getData16u(roi.x0, roi.y0);
		byte* pd;
		if (sameLocationInDst)
		{
			pd = (byte*)dst.getData(roi.x0, roi.y0);
		}
		else
		{
			pd = (byte*)dst.getData();
		}
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				pd[x] = (pl[x] == labelValue ? ps[x] : pd[x]);
			}
			ps += sp;
			pl += lp;
			pd += dp;
		}
	}
	static void CopyBlobs(Img& src, Img& label, vector<ImgBlob>& blobs, Img& dst)
	{
		int n = blobs.size();
		for (int i = 0; i < n; i++)
		{
			CopyBlob(src, label, blobs[i].roi, blobs[i].label, dst, true);
		}
	}
	static void MeanMasked(Img& src, Img& goodMask, int dim, Img& dst)
	{
		int w = src.getWidth();
		int h = src.getHeight();
		int spp = src.getPitch() / 2;
		int dpp = dst.getPitch() / 2;
		int mpp = goodMask.getPitch();
		if (spp != mpp)
		{
			ErrorExit("This doesn't work with different pitch");
		}
		int rad = dim / 2;
		CheckImageArgs16NotInPlace(src, dst);
		CheckImageArgsMask(src, goodMask);
		uint16_t* ps = (uint16_t*)src.getData();
		uint8_t* pm = (uint8_t*)goodMask.getData();
		uint16_t* pd = (uint16_t*)dst.getData();
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				int sum = 0;
				int n = 0;
				for (int ym = y - rad; ym <= y + rad; ym++)
				{
					for (int xm = x - rad; xm <= x + rad; xm++)
					{
						if ((xm >= 0) && (xm < w) && (ym >= 0) && (ym < h))
						{
							int idx = ym * spp + xm;
							if (pm[idx]) 
							{ 
								n++; 
								sum += ps[idx]; 
							}
						}
					}
				}
				if ((n > 0) && (pm[y * spp + x] > 0))
				{
					pd[x] = (uint16_t)(sum / n);
				}
				else
				{
					pd[x] = 0;
				}
			}
			pd += dpp;
		}
	}
	static void Mean(Img& src, int dim, Img& dst)
	{
		int w = src.getWidth();
		int h = src.getHeight();
		Img goodMask(w, h, 8);
		goodMask.clear(255);
		MeanMasked(src, goodMask, dim, dst);
	}
	static void WrapMin(Img& src, int segmentSize, Img& dst)
	{
		int debugYIndex = -1;
		CheckImageArgs16NotInPlace(src, dst);
		int w = src.getWidth();
		int h = src.getHeight();
		Clear(dst);
		int* lastPixelIndices = new int[h];
		for (int y = 0; y < h; y++)
		{
			lastPixelIndices[y] = -1;
			for (int x = w - 1; x >= 0; x--)
			{
				if (src.getPixel16u(x, y) > 0)
				{
					lastPixelIndices[y] = x;
					break;
				}
			}
			if ((y == debugYIndex)) logd.printf("y=%d: lastPixelIndex: %d\n", y, lastPixelIndices[y]);
		}
		for (int y = 0; y < h; y++)
		{
			int anchorIndex = 0; 
			uint16_t anchorVal = 0;
			float lastSlope = 0;
			while (anchorIndex < w)
			{
				for (int x = anchorIndex; x < w; x++)
				{
					anchorVal = src.getPixel16u(x, y);
					if (anchorVal > 0)
					{
						anchorIndex = x;
						break;
					}
				}
				float lowSlope = FLT_MAX;
				int lowIndex = anchorIndex;
				for (int x = anchorIndex + 1; x < MIN(anchorIndex + segmentSize, w); x++)
				{
					uint16_t vs = src.getPixel16u(x, y);
					if (vs > 0)
					{
						float slope = (float)(vs - anchorVal) / (x - anchorIndex);
						if (slope < lowSlope)
						{
							lowSlope = slope;
							lowIndex = x;
						}
					}
				}
				if (lowIndex == anchorIndex)
				{
					dst.setPixel16u(anchorIndex, y, anchorVal);
					anchorIndex++;
				}
				else
				{
					if ((y == debugYIndex)) logd.printf("y=%d: segment: %d to %d, lowSlope: %.3f, lastSlope: %.3f\n", y, anchorIndex, lowIndex, lowSlope, lastSlope);
					if ((lowIndex == lastPixelIndices[y]) && (lowIndex - anchorIndex < 3))
					{
						lowSlope = MIN(lowSlope, lastSlope);
					}
					for (int x = anchorIndex; x <= lowIndex; x++)
					{
						uint16_t val = anchorVal + (uint16_t)roundf(lowSlope * (x - anchorIndex));
						dst.setPixel16u(x, y, val);
					}
					lastSlope = lowSlope;
					if (lowIndex == lastPixelIndices[y])
					{
						anchorIndex = w;
					}
					else
					{
						anchorIndex = lowIndex;
					}
				}
			}
		}
		delete[] lastPixelIndices;
	}
	static void MirrorHorizontal(Img& src, Img& dst)
	{
		CheckImageArgs16NotInPlace(src, dst);
		int w = src.getWidth();
		int h = src.getHeight();
		int sp = src.getPitch() >> 1;
		int dp = dst.getPitch() >> 1;
		uint16_t* ps = (uint16_t*)src.getData();
		uint16_t* pd = (uint16_t*)dst.getData();
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				pd[w - x - 1] = ps[x];
			}
			ps += sp;
			pd += dp;
		}
	}
	static void Min(Img& src1, Img& src2, Img& dst)
	{
		CheckImageArgs16NotInPlace(src1, dst);
		CheckImageArgs16NotInPlace(src2, dst);
		int w = src1.getWidth();
		int h = src1.getHeight();
		int sp1 = src1.getPitch() >> 1;
		int sp2 = src2.getPitch() >> 1;
		int dp = dst.getPitch() >> 1;
		uint16_t* ps1 = (uint16_t*)src1.getData();
		uint16_t* ps2 = (uint16_t*)src2.getData();
		uint16_t* pd = (uint16_t*)dst.getData();
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				pd[x] = MIN(ps1[x], ps2[x]);
			}
			ps1 += sp1;
			ps2 += sp2;
			pd += dp;
		}
	}
	static int EqualsMask(Img& src1, Img& src2, Img& dst)
	{
		CheckImageArgsMask(src1, dst);
		CheckImageArgsMask(src2, dst);
		int w = src1.getWidth();
		int h = src1.getHeight();
		int sp1 = src1.getPitch() >> 1;
		int sp2 = src2.getPitch() >> 1;
		int dp = dst.getPitch();
		uint16_t* ps1 = (uint16_t*)src1.getData();
		uint16_t* ps2 = (uint16_t*)src2.getData();
		byte* pd = (byte*)dst.getData();
		int dcount = 0;
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				if (ps1[x] == ps2[x])
				{
					pd[x] = 0;
				}
				else
				{
					pd[x] = 255;
					dcount++;
				}
			}
			ps1 += sp1;
			ps2 += sp2;
			pd += dp;
		}
		return dcount;
	}
	static void WrapMinBiDir(Img& src, int segmentSize, Img& dst)
	{
		int w = src.getWidth();
		int h = src.getHeight();
		Img tmp1(w, h, 16);
		Img tmp2(w, h, 16);
		WrapMin(src, segmentSize, tmp1);
		MirrorHorizontal(src, tmp2);
		WrapMin(tmp2, segmentSize, dst);
		MirrorHorizontal(dst, tmp2);
		Min(tmp1, tmp2, dst);
	}
	static void WrapMinAllDir(Img& src, int segmentSize, Img& dst)
	{
		int w = src.getWidth();
		int h = src.getHeight();
		Img tmp1(w, h, 16);
		WrapMinBiDir(src, segmentSize, tmp1);
		Img tmp2(w, h, 16);
		Img tmp3(w, h, 16);
		Transpose(tmp1, tmp2);
		WrapMinBiDir(tmp2, segmentSize, tmp3);
		Transpose(tmp3, tmp2);
		Min(tmp1, tmp2, dst);
	}
	static void ClearRoi(Img& image, int x0, int y0, int x1, int y1)
	{
		int sp = image.getPitch() >> 1;
		uint16_t* ps = (uint16_t*)image.getData(0,y0);
		for (int y = y0; y < y1; y++)
		{
			for (int x = x0; x < x1; x++)
			{
				ps[x] = 0;
			}
			ps += sp;
		}
	}
	static const int FuzzMaskCount = 14;
	static void FuzzMask(Img& src, int index, Img& dst)
	{
		int w = src.getWidth();
		int h = src.getHeight();
		Img mask(w, h, 8);
		CreateFuzzMask(index, true, mask);
		ApplyMask(src, mask, true);
	}
	static void CreateFuzzMask(int index, bool smallRegionsNonZero, Img& mask)
	{
		int w = mask.getWidth();
		int h = mask.getHeight();
		byte bgValue = 0;
		byte maskValue = 255;
		if (!smallRegionsNonZero)
		{
			bgValue = 255;
			maskValue = 0;
		}
		mask.clear(bgValue);
		switch (index)
		{
		case 0: 
			break;
		case 1: 
			mask.clear(maskValue);
			break;
		case 2: 
			Img8Util::ClearRoi(mask, 0, 0, 11, 17, maskValue);
			break;
		case 3: 
			Img8Util::ClearRoi(mask, 0, h - 21, 13, h, maskValue);
			break;
		case 4: 
			Img8Util::ClearRoi(mask, w - 32, 0, w, 17, maskValue);
			break;
		case 5: 
			Img8Util::ClearRoi(mask, w - 32, h - 21, w, h, maskValue);
			break;
		case 6: 
			Img8Util::ClearRoi(mask, 0, 0, 1, h, maskValue);
			break;
		case 7: 
			Img8Util::ClearRoi(mask, w - 1, 0, w, h, maskValue);
			break;
		case 8: 
			Img8Util::ClearRoi(mask, 0, 0, w, 1, maskValue);
			break;
		case 9: 
			Img8Util::ClearRoi(mask, 0, h - 1, w, h, maskValue);
			break;
		case 10: 
			Img8Util::ClearRoi(mask, 0, 0, 1, h, maskValue);
			Img8Util::ClearRoi(mask, w - 1, 0, w, h, maskValue);
			Img8Util::ClearRoi(mask, 0, 0, w, 1, maskValue);
			Img8Util::ClearRoi(mask, 0, h - 1, w, h, maskValue);
			break;
		case 11: 
			Img8Util::ClearRoi(mask, 11, 31, 53, 51, maskValue);
			break;
		case 12: 
			Img8Util::ClearRoi(mask, 11, 31, 12, 51, maskValue);
			break;
		case 13: 
			Img8Util::ClearRoi(mask, 11, 31, 53, 32, maskValue);
			break;
		default:
			ErrorExit("Undefined fuzz mask index.");
		}
	}
	static void FillZeros(Img& src, Img& dst)
	{
		CheckImageArgs16NotInPlace(src, dst);
		int w = src.getWidth();
		int h = src.getHeight();
		int sp = src.getPitch() >> 1;
		int dp = dst.getPitch() >> 1;
		uint16_t* ps = (uint16_t*)src.getData();
		uint16_t* pd = (uint16_t*)dst.getData();
		for (int y = 0; y < h; y++)
		{
			int firstNonZero;
			for (firstNonZero = 0; firstNonZero < w; firstNonZero++)
			{
				if (ps[firstNonZero] > 0) break;
			}
			if (firstNonZero == w)
			{
				memset(pd, 0, dp * 2);
			}
			else
			{
				for (int x = firstNonZero; x >= 0; x--)
				{
					pd[x] = ps[firstNonZero];
				}
				uint16_t lastVal = ps[firstNonZero];
				for (int x = firstNonZero + 1; x < w; x++)
				{
					if (ps[x] > 0)
					{
						pd[x] = ps[x];
						lastVal = ps[x];
					}
					else
					{
						pd[x] = lastVal;
					}
				}
			}
			ps += sp;
			pd += dp;
		}
		ps = (uint16_t*)src.getData();
		pd = (uint16_t*)dst.getData();
		for (int y = 0; y < h; y++)
		{
			if (pd[0] == 0)
			{
				if (y == 0)
				{
					memcpy(pd, pd + dp, dp * 2);
				}
				else
				{
					memcpy(pd, pd - dp, dp * 2);
				}
			}
			ps += sp;
			pd += dp;
		}
	}
	static void Transpose(Img& src, Img& dst)
	{
		CheckImageArgs16NotInPlace(src, dst);
		int w = src.getWidth();
		int h = src.getHeight();
		if (w != h) ErrorExit("This only handles w == h");
		int sp = src.getPitch() >> 1;
		int dp = dst.getPitch() >> 1;
		uint16_t* ps = (uint16_t*)src.getData();
		uint16_t* pd = (uint16_t*)dst.getData();
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				pd[y + x * dp] = ps[x];
			}
			ps += sp;
		}
	}
	static void Subtract(Img& src1, Img& src2, Img& dst)
	{
		CheckImageArgs16NotInPlace(src1, dst);
		CheckImageArgs16NotInPlace(src2, dst);
		int w = src1.getWidth();
		int h = src1.getHeight();
		int sp1 = src1.getPitch() >> 1;
		int sp2 = src2.getPitch() >> 1;
		int dp = dst.getPitch() >> 1;
		uint16_t* ps1 = (uint16_t*)src1.getData();
		uint16_t* ps2 = (uint16_t*)src2.getData();
		uint16_t* pd = (uint16_t*)dst.getData();
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				int v = (int)ps1[x] - (int)ps2[x];
				pd[x] = MAX(0, v);
			}
			ps1 += sp1;
			ps2 += sp2;
			pd += dp;
		}
	}
	static void Add(Img& src1, Img& src2, Img& dst)
	{
		CheckImageArgs16NotInPlace(src1, dst);
		CheckImageArgs16NotInPlace(src2, dst);
		int w = src1.getWidth();
		int h = src1.getHeight();
		int sp1 = src1.getPitch() >> 1;
		int sp2 = src2.getPitch() >> 1;
		int dp = dst.getPitch() >> 1;
		uint16_t* ps1 = (uint16_t*)src1.getData();
		uint16_t* ps2 = (uint16_t*)src2.getData();
		uint16_t* pd = (uint16_t*)dst.getData();
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				int v = (int)ps1[x] + (int)ps2[x];
				pd[x] = MIN(65535, v);
			}
			ps1 += sp1;
			ps2 += sp2;
			pd += dp;
		}
	}
	static void Add(Img& src, int v, Img& dst)
	{
		CheckImageArgs16NotInPlace(src, dst);
		int w = src.getWidth();
		int h = src.getHeight();
		int sp = src.getPitch() >> 1;
		int dp = dst.getPitch() >> 1;
		uint16_t* ps = (uint16_t*)src.getData();
		uint16_t* pd = (uint16_t*)dst.getData();
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				uint16_t dv = 0;
				if (ps[x] > 0)
				{
					int nv = ps[x] + v;
					if (nv > 65535) dv = 65535;
					else if (nv < 0) dv = 0;
					else dv = nv;
				}
				pd[x] = dv;
			}
			ps += sp;
			pd += dp;
		}
	}
	static void Median3x3GoodValues(Img& src, Img& dst)
	{
		CheckImageArgs16NotInPlace(src, dst);
		int w = src.getWidth();
		int h = src.getHeight();
		int spp = src.getPitch() / 2;
		int dpp = dst.getPitch() / 2;
		uint16_t* psBegin = (uint16_t*)src.getData();
		uint16_t* psEnd = psBegin + h * spp;
		uint16_t* pd = (uint16_t*)dst.getData();
		const int off[9] = {-spp - 1, -spp, -spp + 1, -1, 0, 1, spp - 1, spp, spp + 1};
		ZeroEdges(dst, 1);
		uint16_t values[9];
		uint16_t* ps = psBegin;
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				int n = 0;
				if (ps[x] > 0)
				{
					for (int j = 0; j < 9; j++)
					{
						int idx = x + off[j];
						if ((ps + idx >= psBegin) && (ps + idx < psEnd))
						{
							uint16_t v = ps[x + off[j]];
							if (v > 0) values[n++] = v;
						}
					}
					if (n == 0)
					{
						pd[x] = 0;
					}
					else if (n > 2)
					{
						int mi = n / 2;
						std::nth_element(values, values + mi, values + n);
						pd[x] = values[mi];
					}
					else
					{
						pd[x] = values[0];
					}
				}
				else
				{
					pd[x] = 0;
				}
			}
			ps += spp;
			pd += dpp;
		}
	}
	static void DistanceTransform(Img& src, Img& dst)
	{
	}
	static void EraseBlob(Img& img, Img& label, ImgBlob& blob)
	{
		IRoi& roi = blob.roi;
		int w = roi.getWidth();
		int h = roi.getHeight();
		byte* ps = img.getData(roi.x0, roi.y0);
		int sp = img.getPitch();
		uint16_t* pl = label.getData16u(roi.x0, roi.y0);
		int lp = label.getPitch() / 2;
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				if (pl[x] == blob.label)
				{
					ps[x] = 0;
				}
			}
			ps += sp;
			pl += lp;
		}
	}
	static void CopyExceptBlobs(Img& src, Img& label, vector<ImgBlob>& droppedBlobs, Img& dst)
	{
		Img8Util::Copy(src, dst);
		int n = droppedBlobs.size();
		for (int i = 0; i < n; i++)
		{
			EraseBlob(dst, label, droppedBlobs[i]);
		}
	}
	static void Affine(Img& src, vector<pair<float, float>> corners, Img& dst)
	{
		CheckImageArgs16NotInPlace(src, dst);
#ifndef USE_IPP
		ErrorExit("Not implemented without IPP yet.");
#endif
	}
	static Transform2d CreateForwardTransform(int w, int h, pair<float, float>& translate, double& scale, double& rotate)
	{
		Transform2d toOrigin;
		toOrigin.SetTranslation(-w / 2, -h / 2);
		Transform2d trans;
		trans.SetTranslation(w / 2 + translate.first, h / 2 + translate.second);
		Transform2d zoom;
		zoom.SetScale((float)scale, (float)scale);
		Transform2d rot;
		rot.SetRotation((float)rotate);
		Transform2d xform;
		xform.Concat(toOrigin);
		xform.Concat(zoom);
		xform.Concat(rot);
		xform.Concat(trans);
		return xform;
	}
	static Transform2d CreateReverseTransform(int w, int h, pair<float, float>& translate, double& scale, double& rotate)
	{
		Transform2d toOrigin;
		toOrigin.SetTranslation(w / 2, h / 2);
		Transform2d trans;
		trans.SetTranslation(-w / 2 - translate.first, -h / 2 - translate.second);
		Transform2d zoom;
		zoom.SetScale(1.0f / (float)scale, 1.0f / (float)scale);
		Transform2d rot;
		rot.SetRotation(-(float)rotate);
		Transform2d xform;
		xform.Concat(trans);
		xform.Concat(rot);
		xform.Concat(zoom);
		xform.Concat(toOrigin);
		return xform;
	}
	static void Affine(Img& src, pair<float, float>& translate, double& scale, double& rotate, Img& dst)
	{
		CheckImageArgs16NotInPlace(src, dst);
		int w = src.getWidth();
		int h = src.getHeight();
		Transform2d revXform = CreateReverseTransform(w, h, translate, scale, rotate);
		Affine(src, revXform, dst);
	}
	static void Affine(Img& src, Transform2d revXform, Img& dst)
	{
		CheckImageArgs16NotInPlace(src, dst);
		int w = src.getWidth();
		int h = src.getHeight();
		int spp = src.getPitch() / 2;
		int dpp = dst.getPitch() / 2;
		uint16_t* ps = (uint16_t*)src.getData();
		uint16_t* pd = (uint16_t*)dst.getData();
		float sx, sy;
		int xi, yi;
		for (int y = 0; y < h; y++)
		{
			for (int x = 0; x < w; x++)
			{
				sx = (float)x;
				sy = (float)y;
				revXform.Transform(sx, sy);
				xi = roundf(sx);
				yi = roundf(sy);
				if ((xi >= 0) && (xi < w) && (yi >= 0) && (yi < h))
				{
					pd[x] = ps[yi * spp + xi];
				}
				else
				{
					pd[x] = 0;
				}
			}
			pd += dpp;
		}
	}
};
struct NsiBlob
{
	ImgBlob imgBlob;
	bool isGood; 
	int sum; 
	int max; 
	ImgMoment moments;
	float strengthScore; 
	float roundScore; 
	float salienceScore; 
	float bsmScore; 
	NsiBlob(): imgBlob()
	{
		isGood = true;
		sum = max = 0;
		strengthScore = roundScore = salienceScore = bsmScore = 0;
	}
	bool CheckDebug(DebugSpec dspec)
	{
		return imgBlob.CheckDebug(dspec);
	}
	float Dist(pair<float, float> pt)
	{
		return imgBlob.Dist(pt);
	}
	float DistInt(pair<int, int> pt)
	{
		return imgBlob.DistInt(pt);
	}
	float DistSq(NsiBlob& other)
	{
		return imgBlob.DistSq(other.imgBlob);
	}
	string ToString()
	{
		char tmp[1024];
		sprintf(tmp, "%s, sum: %d, max: %d, strength: %.1f, round: %.2f, salience: %.2f, bsmScore: %.2f", 
			imgBlob.ToString().c_str(), sum, max, strengthScore, roundScore, salienceScore, bsmScore);
		return (string)tmp;
	}
	string ToLongString()
	{
		char tmp[1024];
		sprintf(tmp, "%s, sum: %d, max: %d, %s, strength: %.1f, round: %.2f, salience: %.2f, bsmScore: %.2f", 
			imgBlob.ToString().c_str(), sum, max, moments.ToString().c_str(), strengthScore, roundScore, salienceScore, bsmScore);
		return (string)tmp;
	}
	float CompSim(NsiBlob& other, bool debug = false)
	{
		if (sum == 0) ErrorExit("Cannot compare blobs where sum is not set.");
		if (max == 0) ErrorExit("Cannot compare blobs where max is not set.");
		if (imgBlob.area == 0) ErrorExit("Cannot compare blobs where imgBlob.area is not set.");
		if (!isGood) ErrorExit("Cannot compare blobs where blob is not good.");
		if (isGood != other.isGood) return 0;
		if ((sum == 0) != (other.sum == 0)) return 0;
		if ((max == 0) != (other.max == 0)) return 0;
		if ((imgBlob.area == 0) != (other.imgBlob.area == 0)) return 0;
		float sumNoiseFloor = 100;
		float maxNoiseFloor = 50;
		float areaNoiseFloor = 4;
		float sumSimScore = MAX(0.0f, 1.0f - pctDiff(sumNoiseFloor + sum, sumNoiseFloor + other.sum)); 
		float maxSimScore = MAX(0.0f, 1.0f - pctDiff(max + maxNoiseFloor, other.max + maxNoiseFloor)); 
		float imgBlobAreaSimScore = MAX(0.0f, 1.0f - pctDiff(imgBlob.area + areaNoiseFloor, other.imgBlob.area + areaNoiseFloor)); 
		float momentSumWeight;
		float momentSimScore = moments.ComputeSimScore(other.moments, momentSumWeight, false);
		float BlobSimScoreFloor = 0.001f;
		sumSimScore = MAX(sumSimScore, BlobSimScoreFloor);
		maxSimScore = MAX(maxSimScore, BlobSimScoreFloor);
		imgBlobAreaSimScore = MAX(imgBlobAreaSimScore, BlobSimScoreFloor);
		momentSimScore = MAX(momentSimScore, BlobSimScoreFloor);
		vector<float> values;
		vector<float> weights;
		values.push_back(sumSimScore);
		weights.push_back(1.0f);
		values.push_back(maxSimScore);
		weights.push_back(1.0f);
		values.push_back(imgBlobAreaSimScore);
		weights.push_back(1.0f);
		values.push_back(momentSimScore);
		weights.push_back(momentSumWeight);
		float score = geometricMean(values);
		if (debug) logd.printf("CompSim: Score: %.3f, SumSS: %.3f, MaxSS: %.3f, AreaSS: %.3f, MomentSS: %.3f\n", score, sumSimScore, maxSimScore, imgBlobAreaSimScore, momentSimScore);
		return score;
	}
	void UpdateScores()
	{
		if (imgBlob.area > 0)
		{
			strengthScore = sum * max * imgBlob.area;
			int roiMaxSide = MAX(imgBlob.roi.getWidth(), imgBlob.roi.getHeight());
			if ((imgBlob.area == 1)
				|| ((imgBlob.area <= 4) && (roiMaxSide == 2))
				)
			{
				roundScore = 1;
			}
			else
			{
				roundScore = moments.ComputeRoundness();
				float rawAreaRoundScore = (float)imgBlob.area / (roiMaxSide * roiMaxSide) * 4 / PI_FLOAT;
				roundScore = (rawAreaRoundScore + roundScore) / 2;
				roundScore = rangeScore(0.2f, 0.95f, roundScore); 
				roundScore = MAX(0.05f, roundScore);
			}
		}
	}
};
#include <string>
#include <vector>
class Detection
{
public:
	string imageSetId;
	int detectionNumber; 
	bool isNeoTruth; 
	vector<DetectionRecord> recs;
	bool hasOverlap; 
	int overlapDetNum; 
	int rank; 
	bool isFound; 
	float velocity; 
	float brightnessMagMean;
	bool isNeoFound; 
	Detection()
	{
		imageSetId = "";
		rank = detectionNumber = -1;
		isNeoTruth = false;
		isFound = isNeoFound = hasOverlap = false;
		overlapDetNum = -1;
		velocity = brightnessMagMean = 0;
	}
	void UpdateDerivedFields()
	{
		float sum = 0;
		for (int i = 1; i < recs.size(); i++)
		{
			float dx = recs[i].xReg - recs[i - 1].xReg;
			float dy = recs[i].yReg - recs[i - 1].yReg;
			sum += pyth(dx, dy);
		}
		float bsum = 0;
		for (int i = 0; i < recs.size(); i++)
		{
			bsum += recs[i].brightnessMag;
		}
		this->velocity = sum / 4;
		this->brightnessMagMean = bsum / 4;
	}
	std::string ToString()
	{
		char tmp[512];
		sprintf(tmp, "imageSetId: %s, detNum: %d, isNeo: %c, fr0 reg coords: %.0f, %.0f", imageSetId.c_str(), detectionNumber, isNeoTruth ? 'Y' : 'N', 
			(recs.size() > 1 ? recs[0].xReg : -1), (recs.size() > 1 ? recs[0].yReg : -1));
		return (string)tmp;
	}
	static vector<Detection> Build(string imageSetId, vector<DetectionRecord>& drs)
	{
		vector<Detection> dets;
		int n = drs.size();
		int lastDetectionNumber = -1;
		for (int i = 0; i < n;)
		{
			Detection det;
			while (drs[i].detectionNumber == lastDetectionNumber)
			{
				assert(false && "Out of sync in build, too many with same det num.");
				i++;
			}
			det.imageSetId = imageSetId;
			det.detectionNumber = drs[i].detectionNumber;
			lastDetectionNumber = drs[i].detectionNumber;
			det.isNeoTruth = drs[i].isNeo;
			det.hasOverlap = drs[i].hasOverlap;
			det.overlapDetNum = drs[i].overlapDetNum;
			det.recs.push_back(drs[i]);
			i++;
			for (int j = 0; j < 3; j++) 
			{
				if (drs[i].detectionNumber == det.detectionNumber)
				{
					det.recs.push_back(drs[i]);
					i++;
				}
				else
				{
					assert(false && "Out of sync in build, not enough with same det num.");
					break;
				}
			}
			if (det.recs.size() == 4)
			{
				det.UpdateDerivedFields();
				dets.push_back(det);
			}
		}
		return dets;
	}
	static std::string ToCsvHeaderString()
	{
		const char* fmt = "imageSetId, detNum, isNeoTruth, hasOverlap, overlapDetNum, isFound, rank, isNeoFound, vPix, brightnessMagMean, xReg, yReg";
		char tmp[4096];
		sprintf(tmp, fmt);
		return (string)tmp;
	}
	std::string ToCsvString()
	{
		char tmp[4096];
		sprintf(tmp, "%s, %d, %c, %c, %d, %c, %d, %c, %.1f, %.1f, %.1f, %.1f",
			imageSetId.c_str(), 
			detectionNumber, 
			(isNeoTruth ? 'Y' : 'N'),
			(hasOverlap ? 'Y' : 'N'),
			overlapDetNum,
			(isFound ? 'Y' : 'N'),
			rank,
			(isNeoFound ? 'Y' : 'N'),
			velocity,
			brightnessMagMean,
			(recs.size() > 1 ? recs[0].xReg : -1), 
			(recs.size() > 1 ? recs[0].yReg : -1)
			);
		return (string)tmp;
	}
	static void WriteCsv(std::string& csvPath, vector<Detection>& dets)
	{
		FILE* csvp = fopen(csvPath.c_str(), "w");
		fprintf(csvp, "%s\n", ToCsvHeaderString().c_str());
		AppendToCsv(csvp, dets);
		fclose(csvp);
	}
	static void AppendToCsv(FILE* fp, vector<Detection>& dets)
	{
		for (int i = 0; i < dets.size(); i++)
		{
			fprintf(fp, "%s\n", dets[i].ToCsvString().c_str());
		}
	}
	void ParseCsvWords(vector<string>& words)
	{
		int i = 0;
		imageSetId = words[i++];
		detectionNumber = atoi(words[i++].c_str());
		isNeoTruth = (words[i++] == "Y");
		hasOverlap = (words[i++] == "Y");
		if (words.size() >= 12)
		{
			overlapDetNum = atoi(words[i++].c_str());
		}
		isFound = (words[i++] == "Y");
		rank = atoi(words[i++].c_str());
		isNeoFound = (words[i++] == "Y");
		velocity = atof(words[i++].c_str());
		brightnessMagMean = atof(words[i++].c_str());
		DetectionRecord dr;
		dr.xReg = atof(words[i++].c_str());
		dr.yReg = atof(words[i++].c_str());
		recs.push_back(dr);
	}
	static vector<Detection> LoadResultCsv(string path)
	{
		vector<Detection> dets;
		return dets;
	}
};
#define NOMINMAX
int debugImageSerialNumber = 0;
void SaveDebugTiff(const char* baseName, Img& image)
{
	if (image.getWidth() == 0) return;
}
void SaveDebugBmp(const char* baseName, Img& image)
{
}
void SaveToImageCatalog(int candidateIndex, int indexOf4, DetectionRecord& det, Img& image)
{
}
string BuildDebugImageBaseName(const char* tag, int candidateIndex, DetectionRecord& det)
{
	char baseName[512];
	sprintf(baseName, "ca-%04d_%01d_det-%05d_%c_%s", candidateIndex, det.frameIndex, det.uniqueId, (det.isReject ? 'r' : 'n'), tag);
	return (string)baseName;
}
void SaveDebugImage(const char* tag, int candidateIndex, DetectionRecord& det, Img& image)
{
}
void SaveDebugImage(int imageSetIndex, int frameIndex, const char* tag, DetectionRecord& det, Img& image)
{
}
void SaveDebugImageRawName(const char* baseName, Img& image)
{
}
void SaveDebugImage(const char* baseName, Img& image)
{
}
void SaveTiff(Img& image, string path)
{
}
void LoadTiff(string path, Img& image)
{
}
void SaveHist(const char* tag, int caIndex, DetectionRecord& det, Img& image)
{
}
void LoadRawImageFile(string& path, vector<int>& imageData)
{
	FILE* fp = fopen(path.c_str(), "rb");
	fseek (fp, 0, SEEK_END);
	size_t lSize = ftell(fp);
	rewind (fp);
	int nImages = (lSize - 27) / IBYTES;
	unsigned short buf[27];
	int read = fread(buf, 1, 27, fp); 
	if (read != 27)
	{
		ErrorExit("Bad read.");
	}
	int restBytes = lSize - 27;
	int nPixels = IDIM * IDIM * nImages;
	assert(restBytes == nPixels * 2);
	imageData.resize(nPixels);
	vector<uint16_t> sbuffer(nPixels);
	read = fread((byte*)sbuffer.data(), 1, restBytes, fp);
	if (read != restBytes)
	{
		if (feof(fp)) logd.printf("Early EOF on read. Tried for %d bytes and got %d\n", restBytes, read);
		if (ferror(fp)) logd.printf("Error flag set on stream\n");
		ErrorExit("Bad read.");
	}
	for (int i = 0; i < nPixels; i++) imageData[i] = (int)sbuffer[i];
	fclose(fp);
}
int LoadRawImageFile(string& path, vector<Img>& images)
{
	FILE* fp = fopen(path.c_str(), "rb");
	fseek (fp, 0, SEEK_END);
	size_t lSize = ftell(fp);
	rewind (fp);
	int nImages = (lSize - 27) / IBYTES;
	unsigned short buf[27];
	int read = fread(buf, 1, 27, fp); 
	if (read != 27)
	{
		ErrorExit("Bad read.");
	}
	Img image(IDIM, IDIM, 16);
	for (int i = 0; i < nImages; i++)
	{
		read = fread(image.getData(), 1, IBYTES, fp);
		if (read != IBYTES)
		{
			if (feof(fp)) logd.printf("Early EOF on read %d\n", i);
			if (ferror(fp)) logd.printf("Error flag set on stream\n");
			ErrorExit("Bad read.");
		}
		images.push_back(image);
	}
	fclose(fp);
	return nImages;
}
void SaveDebugBlobRoi(const char* tag, int frameIndex, int blobIndex, Img& img, ImgBlob& blob, Img* pLabel = NULL)
{
	int dim = 64;
	Img i1(dim, dim, img.getDepthInBits());
	int x0 = roundf(blob.xCentroid - dim / 2);
	int y0 = roundf(blob.yCentroid - dim / 2);
	IRoi roi(x0, y0, x0 + dim, y0 + dim);
	if (img.getDepthInBits() == 16)
	{
		ImgUtil::CopyRoi(img, roi, i1);
	}
	else
	{
		Img8Util::CopyRoi(img, roi, i1);
	}
	if (pLabel != NULL)
	{
		Img label(dim, dim, 16);
		Img labelMask(dim, dim, 8);
		ImgUtil::CopyRoi(*pLabel, roi, label);
		ImgUtil::Equals(label, blob.label, labelMask);
		if (img.getDepthInBits() == 16)
		{
			ImgUtil::ApplyMask(i1, labelMask, false);
		}
		else
		{
			Img8Util::ApplyMask(i1, labelMask, false);
		}
	}
	char tmp[512];
	sprintf(tmp, "%s_%d-%03d", tag, frameIndex, blobIndex);
	SaveDebugImage(tmp, i1);
}
void DrawText8(Img& img, int x0, int y0, const char* text, int ptSize, bool useBlack)
{
	if (img.getDepthInBits() != 8) ErrorExit("Needs to be 8 bit");
}
void RenderColoredRois(Img& img, const char* baseName, int radius, vector<pair<int, int>>& blue, vector<int>& blueLabels,
	vector<pair<int, int>>& magenta, vector<int>& magentaLabels,
	vector<pair<int, int>>& green, vector<int>& greenLabels,
	vector<pair<int, int>>& red, vector<int>& redLabels)
{
}
void RenderColoredRois(Img& img, const char* baseName, int radius, vector<tuple<int, int, uint32_t, string>>& points)
{
}
#include <string>
class FitsImage2
{
public:
	int bzero;
	int width;
	int height;
	int getWidth()
	{
		return width;
	}
	int getHeight()
	{
		return height;
	}
	void UpdateDimsFromOrigImage()
	{
		this->width = origImage.getWidth();
		this->height = origImage.getHeight();
	}
	vector<string> header; 
	vector<double> wcsDoubles; 
	double centerRa; 
	double centerDec;  
	Img origImage; 
	pair<int, int> cropOffset;
	Img fixedImage; 
	Img mask; 
	Img sboMask; 
	Img bgMask; 
	Img fgMask; 
	Img dipMask; 
	Img nsi; 
	Img nsiDiff; 
	vector<Img*> dsPyramid; 
	Wcs wcsConverter;
	FitsImage2();
	~FitsImage2();
	RaDec PixelToWorld(float x, float y);
	void WorldToPixel(RaDec& rd, float& x, float& y);
	RaDec GetImageOriginRaDec();
	void FinalizeLoad(vector<string>& header, vector<double>& wcsDoubles)
	{
		if (origImage.getWidth() == 0) ErrorExit("FitsImage2: Image must already be set.");
		this->header = header;
		UpdateDimsFromOrigImage();
		this->wcsDoubles = wcsDoubles; 
		wcsConverter.init(wcsDoubles);
	}
	void CropImage()
	{
		int dim = 4096;
		int margin = (origImage.getWidth() - dim) / 2;
		Img cropped(dim, dim, 16);
		if ((origImage.getWidth() == 4110) && (origImage.getHeight() == dim))
		{
			cropped.setByCopy(dim, dim, 16, origImage.getPitch(), origImage.getData() + margin * 2);
			cropOffset.first = 8;
			cropOffset.second = 0;
			this->origImage = cropped;
		}
		else
		{
		}
		this->width = origImage.getWidth();
		this->height = origImage.getHeight();
	}
	void DeleteDsPyramid()
	{
		for (int i = 0; i < dsPyramid.size(); i++)
		{
			delete dsPyramid[i];
			dsPyramid[i] = NULL;
		}
		dsPyramid.clear();
	}
	void GetDegreesPerPixel(double& raDegreesPerPixel, double& decDegreesPerPixel)
	{
		RaDec rd1 = PixelToWorld(0, 0);
		RaDec rd2 = PixelToWorld(1, 0);
		RaDec rd3 = PixelToWorld(0, 1);
		raDegreesPerPixel = rd2.RaDegrees - rd1.RaDegrees;
		decDegreesPerPixel = rd3.DecDegrees - rd1.DecDegrees;
	}
};
FitsImage2::FitsImage2()
{
	bzero = 0;
	centerRa = centerDec = 0;
	cropOffset.first = cropOffset.second = 0;
	wcsDoubles.resize(8);
}
FitsImage2::~FitsImage2()
{
}
RaDec FitsImage2::GetImageOriginRaDec()
{
	return PixelToWorld(0, 0);
}
RaDec FitsImage2::PixelToWorld(float x, float y)
{
	x = x + cropOffset.first;
	y = y + cropOffset.second;
	vector<double> world = wcsConverter.convertXY2RADEC(x, y);
	RaDec rd(world[0], world[1]);
	return rd;
}
void FitsImage2::WorldToPixel(RaDec& rd, float& x, float& y)
{
	vector<double> pix = wcsConverter.convertRADEC2XY(rd.RaDegrees, rd.DecDegrees);
	x = (float)(pix[0]) - cropOffset.first;
	y = (float)(pix[1]) - cropOffset.second;
}
class NsiBlobSet
{
public:
	Img labelImage;
	vector<NsiBlob> nsiBlobs;
	PointIndex pointIndex;
	void Create(Img& img, Img& filteredImg, DebugSpec dspec)
	{
		int w = img.getWidth();
		int h = img.getHeight();
		labelImage.resize(w, h, 16);
		ImgUtil::Clear(labelImage);
		if (!ImgUtil::Label(img, labelImage, DO_8_CONNECTED_NSI_BLOBS))
		{
			logd.warn("The input image has more than 65535 blobs. Bailing on additional blobs.");
		}
		vector<ImgBlob> blobs;
		ImgUtil::FindBlobs(labelImage, blobs);
		nsiBlobs.clear();
		vector<ImgBlob> droppedBlobs = MeasureNsiBlobs(img, labelImage, blobs, nsiBlobs, dspec.debug);
		filteredImg.resize(w, h, 8);
		if (droppedBlobs.size() > 0)
		{
			if (dspec.debug) logd.debug("NsiBlobSet: Create: Filtered out %d blobs.", droppedBlobs.size());
			ImgUtil::CopyExceptBlobs(img, labelImage, droppedBlobs, filteredImg);
		}
		else
		{
			if (dspec.debug) logd.debug("NsiBlobSet: Create: No blobs filtered out.");
			Img8Util::Copy(img, filteredImg);
		}
		RebuildPointIndex();
		ComputeSalienceScores(dspec);
	}
	int size()
	{
		return nsiBlobs.size();
	}
	void CullNsiBlobs(Img& nsiDiffImage, Img& filteredImage, DebugSpec dspec)
	{
		int cullRadius = 32;
		int cullMaxBlobCount = 10; 
		int w = nsiDiffImage.getWidth();
		int h = nsiDiffImage.getHeight();
		int nb = nsiBlobs.size();
		vector<NsiBlob> nsiKeepers;
		vector<ImgBlob> keepers; 
		for (int i = 0; i < nsiBlobs.size(); i++)
		{
			pair<int, int> pt = pair<int, int>(roundf(nsiBlobs[i].imgBlob.xCentroid), roundf(nsiBlobs[i].imgBlob.yCentroid));
			vector<int> neighbors = pointIndex.FindNearby(pt, cullRadius);
			bool myScoreIsHighest = true;
			float myScore = nsiBlobs[i].roundScore;
			if (neighbors.size() > cullMaxBlobCount)
			{
				for (int j = 0; j < neighbors.size(); j++)
				{
					if (nsiBlobs[neighbors[j]].roundScore > myScore)
					{
						myScoreIsHighest = false;
						break;
					}
				}
			}
			if (myScoreIsHighest)
			{
				nsiKeepers.push_back(nsiBlobs[i]);
				keepers.push_back(nsiBlobs[i].imgBlob);
			}
		}
		Img8Util::Clear(filteredImage);
		ImgUtil::CopyBlobs(nsiDiffImage, labelImage, keepers, filteredImage);
		if (dspec.debug)
		{
			logd.printf("Culled %d of %d blobs.\n", nb - keepers.size(), nb);
		}
		Create(filteredImage, nsiDiffImage, dspec);
		Img8Util::Copy(nsiDiffImage, filteredImage);
	}
	float ComputeBsmScore(NsiBlob& blob, int blobIndex, float meanMean, DebugSpec dspec)
	{
		float thisMean = 0;
		if (blob.imgBlob.area > 0)
		{
			thisMean = (float)blob.sum / blob.imgBlob.area;
			if (thisMean <= 0)
			{
				ErrorExit("Blob with mean lte 0");
			}
		}
		float meanMeanScore = 1;
		if (meanMean > 0) 
		{
			meanMeanScore = rangeScore(0.0f, 2.0f, thisMean / meanMean);
		}
		float bsmScore = meanMeanScore;
		if (blob.CheckDebug(dspec))
		{
			logd.debug("ComputeSalienceScore: Blob[%d]: label: %d, pt: (%.0f, %.0f): meanMean: %.2f, thisMean: %.2f, bsmScore: %.3f",
				blobIndex, blob.imgBlob.label, blob.imgBlob.xCentroid, blob.imgBlob.yCentroid,
				meanMean, thisMean, bsmScore);
		}
		return bsmScore;
	}
	float ComputeSalienceScore(NsiBlob& blob, int blobIndex, int totalBlobSum, DebugSpec dspec)
	{
		float ratio = (float)blob.sum / totalBlobSum;
		float salienceScore = rangeScore(0.0f, 1.0f, ratio);
		if (blob.CheckDebug(dspec))
		{
			logd.debug("ComputeSalienceScore: Blob[%d]: label: %d, pt: (%.0f, %.0f): sum: %d, tileSum: %d, salienceScore: %.2f",
				blobIndex, blob.imgBlob.label, blob.imgBlob.xCentroid, blob.imgBlob.yCentroid,
				blob.sum, totalBlobSum, salienceScore);
		}
		return salienceScore;
	}
	void ComputeSalienceScores(DebugSpec dspec)
	{
		int gd = pointIndex.GetGridDim();
		for (int y = 0; y < gd; y++)
		{
			for (int x = 0; x < gd; x++)
			{
				vector<int> bi = pointIndex.GetBlobsInGridCell(x, y);
				int totalBlobSum = 0;
				float meanMean = 0;
				ComputeGridCellStats(bi, totalBlobSum, meanMean);
				for (int i = 0; i < bi.size(); i++)
				{
					NsiBlob& blob = nsiBlobs[bi[i]];
					assert(blob.sum > 0);
					assert(blob.sum <= totalBlobSum);
					if (totalBlobSum > 0)
					{
						blob.salienceScore = ComputeSalienceScore(blob, i, totalBlobSum, dspec);
						blob.bsmScore = ComputeBsmScore(blob, i, meanMean, dspec);
					}
					else
					{
					}
				}
			}
		}
	}
	void ComputeSingleSalienceScore(NsiBlob& blob, float& salienceScore, float& bsmScore, DebugSpec dspec)
	{
		salienceScore = 0;
		bsmScore = 0;
		pair<int, int> cell = pointIndex.GetCellCoordsFromUserSpaceCoords(roundf(blob.imgBlob.xCentroid), roundf(blob.imgBlob.yCentroid));
		vector<int> bi = pointIndex.GetBlobsInGridCell(cell.first, cell.second);
		int totalBlobSum = 0;
		float meanMean = 0;
		ComputeGridCellStats(bi, totalBlobSum, meanMean);
		totalBlobSum += blob.sum;
		assert(blob.sum > 0);
		assert(blob.sum <= totalBlobSum);
		if (totalBlobSum > 0)
		{
			salienceScore = ComputeSalienceScore(blob, -1, totalBlobSum, dspec);
			bsmScore = ComputeBsmScore(blob, -1, meanMean, dspec);
		}
	}
private:
	void ComputeGridCellStats(vector<int>& blobIndices, int& totalBlobSum, float& meanMean)
	{
		totalBlobSum = 0;
		meanMean = 0;
		int n = blobIndices.size();
		int meanSum = 0;
		for (int i = 0; i < n; i++)
		{
			totalBlobSum += nsiBlobs[blobIndices[i]].sum;
			if (nsiBlobs[blobIndices[i]].imgBlob.area > 0)
			{
				meanSum += nsiBlobs[blobIndices[i]].sum / nsiBlobs[blobIndices[i]].imgBlob.area;
			}
		}
		if (n > 0)
		{
			meanMean = (float)meanSum / n;
		}
	}
	void RebuildPointIndex()
	{
		int w = labelImage.getWidth();
		int h = labelImage.getHeight();
		pointIndex.Init(MAX(w, h), NSI_BLOB_INDEX_GRID_DIM_PIX); 
		int n = nsiBlobs.size();
		for (int i = 0; i < n; i++)
		{
			pointIndex.AddRound(nsiBlobs[i].imgBlob.xCentroid, nsiBlobs[i].imgBlob.yCentroid, i);
		}
	}
public:
	static vector<ImgBlob> MeasureNsiBlobs(Img& img, Img& label, vector<ImgBlob>& blobs, vector<NsiBlob>& nblobs, bool debug)
	{
		int nb = blobs.size();
		vector<ImgBlob> droppers;
		nblobs.reserve(nb);
		for (int i = 0; i < nb; i++)
		{
			if ((blobs[i].roi.getWidth() < MAX_NSI_BLOB_DIAMETER_PIX) && (blobs[i].roi.getHeight() < MAX_NSI_BLOB_DIAMETER_PIX))
			{
				NsiBlob nb;
				nb.imgBlob = blobs[i];
				nb.isGood = true;
				nblobs.push_back(nb);
			}
			else
			{
				droppers.push_back(blobs[i]);
			}
		}
		nblobs.shrink_to_fit();
		nb = nblobs.size();
		Img tile(MAX_NSI_BLOB_DIAMETER_PIX, MAX_NSI_BLOB_DIAMETER_PIX, 8);
		Img i1(MAX_NSI_BLOB_DIAMETER_PIX, MAX_NSI_BLOB_DIAMETER_PIX, 8);
		Img i2(MAX_NSI_BLOB_DIAMETER_PIX, MAX_NSI_BLOB_DIAMETER_PIX, 8);
		int w = img.getWidth();
		int h = img.getHeight();
		for (int i = 0; i < nb; i++)
		{
			IRoi roi = nblobs[i].imgBlob.roi;
			roi.intersect(0, 0, w, h); 
			Img8Util::Clear(tile);
			ImgUtil::CopyBlob(img, label, roi, nblobs[i].imgBlob.label, tile); 
			byte min, max;
			int sum, countNonZero;
			Img8Util::Stats(tile, min, max, sum, countNonZero);
			nblobs[i].max = max;
			nblobs[i].sum = sum;
			IRoi iro(0, 0, roi.getWidth(), roi.getHeight());
			nblobs[i].moments.Compute(tile, iro);
			nblobs[i].UpdateScores(); 
			nblobs[i].imgBlob.xCentroid = roi.x0 + nblobs[i].moments.xCentroid;
			nblobs[i].imgBlob.yCentroid = roi.y0 + nblobs[i].moments.yCentroid;
		}
		nblobs = PreFilterOnRoundness(nblobs, droppers);
		return droppers;
	}
	static vector<NsiBlob> PreFilterOnRoundness(vector<NsiBlob>& blobs, vector<ImgBlob>& droppers)
	{
		vector<NsiBlob> keepers;
		int n = blobs.size();
		for (int i = 0; i < n; i++)
		{
			if ((blobs[i].imgBlob.area < 6) || (blobs[i].roundScore > MIN_NSI_BLOB_PREFILTER_ROUNDNESS_SCORE))
			{
				keepers.push_back(blobs[i]);
			}
			else
			{
				droppers.push_back(blobs[i].imgBlob);
			}
		}
		return keepers;
	}
};
#define LOG_FEATURE(x) if (checkIsNormalFloat(x)) x = log(x)
#define LOG_FEATURE_OFFSET(x,y) if (checkIsNormalFloat(x)) x = log(x+y)
class CAFeat
{
public:
	int index; 
	string imageSetId; 
	float simpleScore; 
	float sboProxScore; 
	float caProxScore; 
	int nCAsThisImageSet; 
	float salienceScore; 
	float caBsmScore; 
	float lineFitRmsError; 
	CFeatureStats blobAreaStats;
	CFeatureStats blobSumStats;
	CFeatureStats blobMaxStats;
	float velocity;
	float brightness;
	float roundness;
	int nBlobs;
	float nfc;
	pair<int, int> fr0RegCoords;
	int corrDetNumber; 
	int classIndex; 
	bool isNeoComp; 
	float pClass1;
	float pNeo;
	float pOverlap;
private:
	void Init()
	{
		index = -1;
		imageSetId = "unset";
		corrDetNumber = -1;
		classIndex = -1;
		isNeoComp = false;
		simpleScore = 0;
		sboProxScore = 0;
		caProxScore = 0;
		nCAsThisImageSet = 0;
		salienceScore = 0;
		caBsmScore = 0;
		lineFitRmsError = 0;
		velocity = brightness = roundness = 0;
		nBlobs = 0;
		nfc = 0;
		pClass1 = pNeo = pOverlap = -1;
	}
public:
	CAFeat()
	{
		Init();
	}
	void SetDetectionNumber(int corrDetNumber)
	{
		this->corrDetNumber = corrDetNumber;
		if (corrDetNumber < 0)
		{
			this->classIndex = 0;
		}
		else
		{
			this->classIndex = 1;
		}
	}
	std::string ToString()
	{
		char tmp[1024];
		sprintf(tmp, "%s, %d: nBlobs: %d, nfc: %.2f, ss: %.2f, sbops: %.2f, caPs: %.2f, salience: %.2f, caBsmScore: %.2f, pClass1: %.2f, pNeo: %.2f, pOverlap: %.2f, det: %d",
			imageSetId.c_str(), index, nBlobs, nfc, simpleScore, sboProxScore, caProxScore, salienceScore, caBsmScore,
			pClass1, pNeo, pOverlap, corrDetNumber);
		return (string)tmp;
	}
	void ParseCsvWords(vector<string>& words)
	{
		int i = 0;
		imageSetId = words[i++];
		index = atoi(words[i++].c_str());
		corrDetNumber = atoi(words[i++].c_str());
		classIndex = (corrDetNumber < 0 ? 0 : 1);
		lineFitRmsError = atof(words[i++].c_str());
		simpleScore = atof(words[i++].c_str());
		if (words.size() >= 7)
		{
			pClass1 = atof(words[i++].c_str());
		}
		if (words.size() >= 8)
		{
			pNeo = atof(words[i++].c_str());
		}
		if (words.size() >= 9)
		{
			pOverlap = atof(words[i++].c_str());
		}
		if (words.size() >= 11)
		{
			fr0RegCoords.first = atoi(words[i++].c_str());
			fr0RegCoords.second = atoi(words[i++].c_str());
		}
		isNeoComp = (words[i++] == "Y");
	}
	static vector<CAFeat> LoadResultCsv(string path)
	{
		vector<CAFeat> cas;
		return cas;
	}
	static const int CAF_N_FEATURES = 17; 
	static const int CAF_N_CSV_COLS = CAF_N_FEATURES + 3; 
	void ParseFeatureCsvWords(vector<string>& words)
	{
		int i = 0;
		index = atoi(words[i++].c_str());
		if (words.size() == CAF_N_CSV_COLS)
		{
			imageSetId = words[i++];
		}
		else if (words.size() != CAF_N_CSV_COLS - 1)
		{
			ErrorExit("Wrong word count in csv.");
		}
		nBlobs = atoi(words[i++].c_str());
		nfc = atof(words[i++].c_str());
		blobAreaStats.mean = atof(words[i++].c_str());
		blobMaxStats.mean = atof(words[i++].c_str());
		blobSumStats.mean = atof(words[i++].c_str());
		blobAreaStats.sigmaPct = atof(words[i++].c_str());
		blobMaxStats.sigmaPct = atof(words[i++].c_str());
		blobSumStats.sigmaPct = atof(words[i++].c_str());
		simpleScore = atof(words[i++].c_str());
		sboProxScore = atof(words[i++].c_str());
		caProxScore = atof(words[i++].c_str());
		salienceScore = atof(words[i++].c_str());
		caBsmScore = atof(words[i++].c_str());
		lineFitRmsError = atof(words[i++].c_str());
		velocity = atof(words[i++].c_str());
		brightness = atof(words[i++].c_str());
		roundness = atof(words[i++].c_str());
		classIndex = atoi(words[i++].c_str());
		corrDetNumber = (classIndex == 1 ? 9999999 : -1);
	}
	static void LoadFeatureCsv(std::string& path, vector<CAFeat>& cs)
	{
		vector<string> lines;
		logd.printf("Loading csv file...\n");
		if (!FileUtil::SnarfFile(path.c_str(), lines))
		{
			ErrorExit("Failed to load CSV.");
		}
		int n = lines.size();
		logd.printf("Processing %d lines.\n", n);
		for (int i = 1; i < n; i++)
		{
			vector<string> words = StringUtils::split(lines[i], ',');
			CAFeat c;
			c.ParseFeatureCsvWords(words);
			cs.push_back(c);
		}
	}
	vector<Feature> GetFeatures()
	{
		vector<Feature> fs;
		int nBlobs = this->nBlobs;
		float nfc = this->nfc;
		float areaMean = blobAreaStats.mean;
		float maxMean = blobMaxStats.mean;
		float sumMean = blobSumStats.mean;
		float areaSigmaPct = blobAreaStats.sigmaPct;
		float maxSigmaPct = blobMaxStats.sigmaPct;
		float sumSigmaPct = blobSumStats.sigmaPct;
		float simpleScore = this->simpleScore;
		float sboProxScore = this->sboProxScore;
		float caProxScore = this->caProxScore;
		float salienceScore = this->salienceScore;
		float caBsmScore = this->caBsmScore;
		float lineFitRmsError = this->lineFitRmsError;
		float velocity = this->velocity;
		float brightness = this->brightness;
		float roundness = this->roundness;
		nfc = NanIfBadOrLteZero(nfc);
		simpleScore = NanIfBadOrLteZero(simpleScore);
		salienceScore = NanIfBadOrLteZero(salienceScore);
		caBsmScore = NanIfBadOrLteZero(caBsmScore);
		velocity = NanIfBadOrLteZero(velocity);
		roundness = NanIfBadOrLteZero(roundness);
		fs.push_back(Feature("nBlobs", nBlobs));
		fs.push_back(Feature("nfc", nfc));
		fs.push_back(Feature("areaMean", areaMean));
		fs.push_back(Feature("maxMean", maxMean));
		fs.push_back(Feature("sumMean", sumMean));
		fs.push_back(Feature("areaSigmaPct", areaSigmaPct));
		fs.push_back(Feature("maxSigmaPct", maxSigmaPct));
		fs.push_back(Feature("sumSigmaPct", sumSigmaPct));
		fs.push_back(Feature("simpleScore", simpleScore));
		fs.push_back(Feature("sboProxScore", sboProxScore));
		fs.push_back(Feature("caProxScore", caProxScore));
		fs.push_back(Feature("salienceScore", salienceScore));
		fs.push_back(Feature("caBsmScore", caBsmScore));
		fs.push_back(Feature("lineFitRmsError", lineFitRmsError));
		fs.push_back(Feature("velocity", velocity));
		fs.push_back(Feature("brightness", brightness));
		fs.push_back(Feature("roundness", roundness));
		if (fs.size() != CAF_N_FEATURES)
		{
			ErrorExit("Number of features is wrong, see CAF_N_FEATURES");
		}
		for (int i = 0; i < fs.size(); i++)
		{
			fs[i].isSelected = true;
		}
		return fs;
	}
	static void GetSamples(vector<CAFeat>& cas, vector<vector<float>>& samples, vector<int>& truthClasses)
	{
		int n = cas.size();
		samples.clear();
		samples.reserve(n);
		truthClasses.clear();
		truthClasses.reserve(n);
		for (int i = 0; i < n; i++)
		{
			vector<Feature> features = cas[i].GetFeatures();
			samples.push_back(Feature::GetAllValues(features));
			truthClasses.push_back(cas[i].classIndex);
		}
	}
	struct CAFClass1GreaterThan
	{
		inline bool operator() (const CAFeat& c1, const CAFeat& c2)
		{
			if ((c1.pClass1 < 0) || (c2.pClass1 < 0)) ErrorExit("Comparing a CAFeat pClass1 that is not set.");
			return (c1.pClass1 > c2.pClass1);
		}
	};
	static void SortOnPClass1(vector<CAFeat>& v)
	{
		std::sort(v.begin(), v.end(), CAFClass1GreaterThan());
	}
	struct CFSimpleScoreGreaterThan
	{
		inline bool operator() (const CAFeat& c1, const CAFeat& c2)
		{
			return (c1.simpleScore > c2.simpleScore);
		}
	};
	static void Sort(vector<CAFeat>& v)
	{
		std::sort(v.begin(), v.end(), CFSimpleScoreGreaterThan());
	}
	struct CAFNeoGreaterThan
	{
		inline bool operator() (const CAFeat& c1, const CAFeat& c2)
		{
			if ((c1.pNeo < 0) || (c2.pNeo < 0)) ErrorExit("Comparing a CAFeat pNeo that is not set.");
			return (c1.pNeo > c2.pNeo);
		}
	};
	static void SortOnPNeo(vector<CAFeat>& v)
	{
		std::sort(v.begin(), v.end(), CAFNeoGreaterThan());
	}
};
#define LOG_FEATURE(x) if (checkIsNormalFloat(x)) x = log(x)
#define LOG_FEATURE_OFFSET(x,y) if (checkIsNormalFloat(x)) x = log(x+y)
class CA
{
#pragma region Fields
	public:
		int index; 
		string imageSetId; 
		vector<NsiBlob> blobs; 
		vector<float> emBlobScores; 
		vector<float> simScores; 
		float simpleScore; 
		float sboProxScore; 
		float caProxScore; 
		int nCAsThisImageSet; 
		float salienceScore; 
		float caBsmScore; 
		int corrDetNumber; 
		bool isNeo; 
		LeastSquaresFitResult xLsr;
		LeastSquaresFitResult yLsr;
		float lineFitRmsError; 
		vector<int> blobIndices;
		float pClass1;
		float pNeo;
		float pOverlap;
		float pFinal; 
	private:
		pair<float, float> vec; 
		vector<pair<float, float>> coords; 
public:
		vector<RaDec> radecs;  
public: 
		bool featureFieldsSet;
		CFeatureStats blobAreaStats;
		CFeatureStats blobSumStats;
		CFeatureStats blobMaxStats;
		float velocity;
		float brightness;
		float roundness;
		int nBlobs;
		float nfc;
		int classIndex; 
#pragma endregion
#pragma region Private
	private:
		void init()
		{
			xLsr.a = FLT_MAX;
			yLsr.a = FLT_MAX;
			index = -1;
			corrDetNumber = -1;
			isNeo = false;
			imageSetId = "unset";
			coords.resize(4);
			blobs.resize(4);
			emBlobScores.resize(4);
			radecs.resize(4);
			for (int i = 0; i < 4; i++)
			{
				blobs[i].imgBlob.label = 0; 
				coords[i].first = -1;
				coords[i].second = -1;
			}
			simScores.resize(3);
			for (int i = 0; i < 3; i++) simScores[i] = -1;
			vec.first = 0;
			vec.second = 0;
			simpleScore = 0;
			sboProxScore = 0;
			caProxScore = 0;
			nCAsThisImageSet = 0;
			salienceScore = 0;
			caBsmScore = 0;
			lineFitRmsError = 0;
			featureFieldsSet = false;
			velocity = brightness = roundness = 0;
			nBlobs = 0;
			nfc = 0;
			classIndex = -1;
			pClass1 = pNeo = pOverlap = -1;
			pFinal = -1;
		}
		void SetBlobSub(int frameIndex, NsiBlob& b)
		{
			if (blobs[frameIndex].imgBlob.label > 0)
			{
				ErrorExit("Blob already set in CA, can't overwrite with this method.");
			}
			blobs[frameIndex] = b;
			coords[frameIndex].first = b.imgBlob.xCentroid;
			coords[frameIndex].second = b.imgBlob.yCentroid;
			xLsr.a = FLT_MAX;
		}
		void SetTwoBlobsSub(int frameIndexA, NsiBlob& a, int frameIndexB, NsiBlob& b)
		{
			SetBlobSub(frameIndexA, a);
			SetBlobSub(frameIndexB, b);
			int di = abs(frameIndexB - frameIndexA);
			if (di == 1)
			{
				int j = MIN(frameIndexA, frameIndexB);
				simScores[j] = a.CompSim(b);
			}
			UpdatePerBlobChanges();
		}
		void FindFirstAndSecondSetBlobs(int& i1, int& i2)
		{
			i1 = i2 = -1;
			for (int i = 0; i < 4; i++)
			{
				if (CheckIsBlobSet(i) && (i1 < 0))
				{
					i1 = i;
				}
				else if (CheckIsBlobSet(i) && (i2 < 0))
				{
					i2 = i;
				}
			}
		}
		vector<pair<float, float>> GetPresentCentroids()
		{
			vector<pair<float, float>> list;
			for (int i = 0; i < 4; i++)
			{
				if (CheckIsBlobSet(i))
				{
					list.push_back(pair<float, float>(blobs[i].imgBlob.xCentroid, blobs[i].imgBlob.yCentroid));
				}
			}
			return list;
		}
		void DoLineFit()
		{
			if (xLsr.a == FLT_MAX)
			{
				vector<pair<float, float>> centroids = GetPresentCentroids();
				int n = centroids.size();
				if (n < 2)
				{
					ErrorExit("Need at least 2 points for this.");
					return;
				}
				float times[4];
				int j = 0;
				for (int i = 0; i < 4; i++)
				{
					if (CheckIsBlobSet(i))
					{
						times[j++] = i;
					}
				}
				float xs[4], ys[4], weights[4];
				for (int i = 0; i < n; i++)
				{
					xs[i] = centroids[i].first;
					ys[i] = centroids[i].second;
					weights[i] = 1; 
				}
				leastSquaresFitWeighted(times, xs, weights, n, false, xLsr, NULL);
				leastSquaresFitWeighted(times, ys, weights, n, false, yLsr, NULL);
				lineFitRmsError = geometricMean(xLsr.stddev, yLsr.stddev);
			}
		}
		void UpdatePerBlobChanges()
		{
			if (xLsr.a == FLT_MAX)
			{
				DoLineFit();
				for (int i = 0; i < 4; i++)
				{
					coords[i].first = xLsr.a + xLsr.b * i;
					coords[i].second = yLsr.a + yLsr.b * i;
				}
				vec.first = xLsr.b;
				vec.second = yLsr.b;
			}
			UpdateScores();
			UpdateFeatureFields();
		}
#pragma endregion
#pragma region Main
	public:
		CA()
		{
			init();
		}
		CA(int frameIndexA, NsiBlob& a, int frameIndexB, NsiBlob& b)
		{
			init();
			SetTwoBlobsSub(frameIndexA, a, frameIndexB, b);
		}
		bool CheckDebug(DebugSpec dspec)
		{
			if (this->index == dspec.caIndex) return true;
			for (int i = 0; i < 4; i++)
			{
				if (CheckIsBlobSet(i) && (blobs[i].DistInt(dspec.debugBlobCoords) <= DEBUG_COORDS_TOLERANCE_PIX))
				{
					return true;
				}
			}
			return false;
		}
		void InitBlobIndices()
		{
			blobIndices.resize(4);
			for (int i = 0; i < 4; i++) blobIndices[i] = -1;
		}
		bool CheckBlobIndicesEqual(CA& other)
		{
			for (int i = 0; i < 4; i++)
			{
				if ((blobIndices[i] != -1) && (other.blobIndices[i] != -1) && (blobIndices[i] != other.blobIndices[i])) return false;
			}
			return true;
		}
		void SetBlob(int frameIndex, NsiBlob& b)
		{
			SetBlobSub(frameIndex, b);
			UpdatePerBlobChanges();
		}
		void SetDetectionNumber(int corrDetNumber)
		{
			this->corrDetNumber = corrDetNumber;
			if (corrDetNumber < 0)
			{
				this->classIndex = 0;
			}
			else
			{
				this->classIndex = 1;
			}
		}
		int GetDetectionNumber() 
		{
			return this->corrDetNumber;
		}
		bool CheckIsBlobSet(int frameIndex)
		{
			return blobs[frameIndex].imgBlob.label > 0;
		}
		bool CheckIsSetBlobAdjacent(int frameIndex)
		{
			if (CheckIsBlobSet(frameIndex)) return true;
			if ((frameIndex > 0) && CheckIsBlobSet(frameIndex - 1)) return true;
			if ((frameIndex < 3) && CheckIsBlobSet(frameIndex + 1)) return true;
			return false;
		}
		int GetClassIndex()
		{
			if ((corrDetNumber >= 0) != (classIndex == 1))
			{
				ErrorExit("GetClassIndex: corrDetNumber and classIndex don't agree.");
			}
			return classIndex;
		}
		int GetBlobCount()
		{
			int sum = 0;
			for (int i = 0; i < blobs.size(); i++) sum += CheckIsBlobSet(i);
			return sum;
		}
		float GetFrameCountFloat()
		{
			float sum = 0;
			for (int i = 0; i < 4; i++)
			{
				if (CheckIsBlobSet(i))
				{
					sum++;
				}
				else
				{
					sum += emBlobScores[i];
				}
			}
			return sum;
		}
		int GetFrameCount()
		{
			int sum = 0;
			for (int i = 0; i < blobs.size(); i++)
			{
				if (CheckIsBlobSet(i) || (emBlobScores[i] > 0.5f))
				{
					sum++;
				}
			}
			return sum;
		}
		vector<pair<float, float>> GetFittedCentroids()
		{
			UpdatePerBlobChanges();
			return coords;
		}
		pair<float, float> GetVector()
		{
			if (xLsr.a == FLT_MAX)
			{
				ErrorExit("Cannot call GetVector before doing line fit.");
			}
			return vec;
		}
		vector<RaDec> GetRaDecs() const
		{
			return radecs;
		}
		void SetRaDec(int i, RaDec rd)
		{
			radecs[i] = rd;
		}
		void UpdateSimpleScore()
		{
			float blobSimScore = ComputeBlobSimScore();
			float meanRoundness = ComputeBlobRoundnessScore();
			this->simpleScore = meanRoundness * blobSimScore * GetBlobCount();
		}
		pair<float, float> GetFittedCoords(int frameIndex)
		{
			UpdatePerBlobChanges();
			return coords[frameIndex];
		}
		void SetTwoBlobs(int frameIndexA, NsiBlob& a, int frameIndexB, NsiBlob& b)
		{
			SetTwoBlobsSub(frameIndexA, a, frameIndexB, b);
		}
		float ComputeSimScore(CA& other)
		{
			if ((coords[0].first < 0) || (coords[0].second < 0)) ErrorExit("Cannot compare if coords[0] is not set.");
			if ((other.coords[0].first < 0) || (other.coords[0].second < 0)) ErrorExit("Cannot compare if coords[0] is not set.");
			float dx = coords[0].first - other.coords[0].first;
			float dy = coords[0].second - other.coords[0].second;
			float dp = sqrtf(dx*dx + dy*dy);
			float distScore = 1.0f - rangeScore(2.0f, 5.0f, dp);
			float dvx = vec.first - other.vec.first;
			float dvy = vec.second - other.vec.second;
			float dv = sqrtf(dvx*dvx + dvy*dvy);
			float dvScore = 1.0f - rangeScore(1.0f, 5.0f, dv);
			return geometricMean(distScore, dvScore);
		}
		static void AddAllToPointIndex(vector<CA>& cas, int frameIndex, PointIndex& index)
		{
			int n = cas.size();
			for (int i = 0; i < n; i++)
			{
				if (cas[i].CheckIsBlobSet(frameIndex))
				{
					index.AddRound(cas[i].blobs[frameIndex].imgBlob.xCentroid, cas[i].blobs[frameIndex].imgBlob.yCentroid, i);
				}
			}
		}
		float ComputeBlobSimScore()
		{
			vector<float> scores;
			int i1, i2;
			FindFirstAndSecondSetBlobs(i1, i2);
			if ((i1 < 0) || (i2 < 0)) ErrorExit("Needs to be at least 2 blobs set.");
			int lastSetBlob = i1;
			for (int i = i1 + 1; i < 4; i++)
			{
				if (CheckIsBlobSet(i))
				{
					scores.push_back(blobs[i].CompSim(blobs[lastSetBlob]));
					lastSetBlob = i;
				}
			}
			return geometricMean(scores);
		}
		void UpdateScores()
		{
			int n = 0;
			float ss = 0;
			float bsmSum = 0;
			for (int i = 0; i < 4; i++)
			{
				if (CheckIsBlobSet(i))
				{
					n++;
					ss += blobs[i].salienceScore;
					bsmSum += blobs[i].bsmScore;
				}
			}
			if (n > 0)
			{
				this->salienceScore = ss / n;
				this->caBsmScore = bsmSum / n;
			}
			else
			{
				this->salienceScore = 0;
			}
		}
		void ApplyPostProcessingScores(bool debug)
		{
			float origSimpleScore = simpleScore;
			float bright = ComputeBlobBrightnessScore();
			float sboFarAwayScore = (1.0f - sboProxScore);
			float sboDamper = CLIP(bright + sboFarAwayScore, 0.0f, 1.0f);
			simpleScore *= sboDamper;
			if (debug) logd.debug("CA::ApplyPostProcessingScores: CA: %d, Orig simpleScore: %.3f, sboProxScore: %.3f, bright: %.3f, damper: %.3f, new simpleScore: %.3f",
				index, origSimpleScore, sboProxScore, bright, sboDamper, simpleScore);
			origSimpleScore = simpleScore;
			float caRangedProxScore = rangeScore(1.0f, 15.0f, caProxScore); 
			float caDamper = rangeScore(-0.5f, 1.0f, 1.0f - caRangedProxScore); 
			simpleScore *= caDamper;
			if (debug) logd.debug("CA::ApplyPostProcessingScores: CA: %d, Orig simpleScore: %.3f, caProxScore: %.3f, caRangedProxScore: %.3f, damper: %.3f, new simpleScore: %.3f",
				index, origSimpleScore, caProxScore, caRangedProxScore, caDamper, simpleScore);
			origSimpleScore = simpleScore;
			simpleScore *= salienceScore;
			if (debug) logd.debug("CA::ApplyPostProcessingScores: CA: %d, Orig simpleScore: %.3f, salience: %.3f, new simpleScore: %.3f",
				index, origSimpleScore, salienceScore, simpleScore);
		}
		float ComputeBlobSimScore(NsiBlob& b)
		{
			vector<float> scores;
			for (int i = 0; i < 4; i++)
			{
				if (CheckIsBlobSet(i))
				{
					scores.push_back(blobs[i].CompSim(b));
				}
			}
			return geometricMean(scores);
		}
		static CA CreateDummy()
		{
			NsiBlob b1;
			b1.imgBlob.xCentroid = 11;
			b1.imgBlob.yCentroid = 20;
			b1.sum = 100;
			b1.max = 200;
			b1.imgBlob.area = 300;
			b1.imgBlob.label = 1;
			NsiBlob b2;
			b2.imgBlob.xCentroid = 10;
			b2.imgBlob.yCentroid = 15;
			b2.sum = 100;
			b2.max = 200;
			b2.imgBlob.area = 300;
			b2.imgBlob.label = 2;
			CA ca(0, b1, 1, b2);
			return ca;
		}
#pragma endregion
#pragma region Standard Methods
		std::string ToString()
		{
			char tmp[1024];
			sprintf(tmp, "%d: b/fcount: %d/%d, ss: %.2f, sbops: %.2f, caPs: %.2f, salience: %.2f, caBsm: %.2f, start: (%.0f, %.0f), pClass1: %.2f, pNeo: %.2f, pOverlap: %.2f, det: %d",
				index, GetBlobCount(), GetFrameCount(), simpleScore, sboProxScore, caProxScore, salienceScore, caBsmScore,
				coords[0].first, coords[0].second,
				pClass1, pNeo, pOverlap, corrDetNumber);
			return (string)tmp;
		}
		void Display()
		{
			logd.printf("    CA: %s\n", ToString().c_str());
			logd.printf("    Blobs:\n");
			for (int i = 0; i < blobs.size(); i++)
			{
				if (CheckIsBlobSet(i))
				{
					logd.printf("        blob[%d]: %s\n", i, blobs[i].ToString().c_str());
				}
				else
				{
					logd.printf("        blob[%d]: emb score: %.3f\n", i, emBlobScores[i]);
				}
			}
			logd.printf("    Coords:\n");
			for (int i = 0; i < 4; i++)
			{
				logd.printf("        coords[%d]: %.1f, %.1f\n", i, coords[i].first, coords[i].second);
			}
			logd.printf("    LSR X: %s\n", xLsr.ToString().c_str());
			logd.printf("    LSR Y: %s\n", yLsr.ToString().c_str());
		}
#pragma endregion
#pragma region Results Csv
		static std::string ToCsvHeaderString()
		{
			const char* fmt = "imageSetId, caIndex, detNum, lineFitRmsError, simpleScore, pClass1, pNeo, pOverlap, isNeoComp, xReg0, yReg0";
			char tmp[4096];
			sprintf(tmp, fmt);
			return (string)tmp;
		}
		void AppendToCsv(FILE* fp)
		{
			pair<float, float> coords = pair<float, float>(0.0f, 0.0f);
			if (GetBlobCount() > 1)
			{
				coords = GetFittedCoords(0);
			}
			fprintf(fp, "%s, %d, %d, %.3f, %.3f, %.3f, %.3f, %.3f, %c, %d, %d\n", 
				imageSetId.c_str(), index, corrDetNumber, lineFitRmsError, simpleScore, pClass1, pNeo, pOverlap, (isNeo ? 'Y' : 'N'), roundf(coords.first), roundf(coords.second));
		}
		static void WriteCsv(std::string& path, vector<CA>& cas)
		{
			if (cas.size() == 0) return;
			FILE* csvp = fopen(path.c_str(), "w");
			fprintf(csvp, "%s\n", cas[0].ToCsvHeaderString().c_str());
			for (int i = 0; i < cas.size(); i++)
			{
				cas[i].AppendToCsv(csvp);
			}
			fclose(csvp);
		}
		void ParseCsvWords(vector<string>& words)
		{
			int i = 0;
			imageSetId = words[i++];
			index = atoi(words[i++].c_str());
			corrDetNumber = atoi(words[i++].c_str());
			classIndex = (corrDetNumber < 0 ? 0 : 1);
			float unused = atof(words[i++].c_str()); 
			simpleScore = atof(words[i++].c_str());
			if (words.size() >= 7)
			{
				pClass1 = atof(words[i++].c_str());
			}
			if (words.size() >= 8)
			{
				pNeo = atof(words[i++].c_str());
			}
			if (words.size() >= 9)
			{
				pOverlap = atof(words[i++].c_str());
			}
			if (words.size() >= 11)
			{
				i += 2;
			}
			isNeo = (words[i++] == "Y");
		}
		static vector<CA> LoadResultCsv(string path)
		{
			vector<CA> cas;
			return cas;
		}
#pragma endregion
#pragma region Compute Features
		void ComputePrimaryBlobStats(bool debug = false)
		{
			int n = 0;
			vector<float> blobAreas;
			vector<float> blobSums;
			vector<float> blobMaxes;
			for (int i = 0; i < 4; i++)
			{
				if (CheckIsBlobSet(i))
				{
					n++;
					blobAreas.push_back(blobs[i].imgBlob.area);
					blobSums.push_back(blobs[i].sum);
					blobMaxes.push_back(blobs[i].max);
				}
			}
			blobAreaStats.ComputeClustered(blobAreas, false, debug);
			blobSumStats.ComputeClustered(blobSums, false, debug);
			blobMaxStats.ComputeClustered(blobMaxes, false, debug);
		}
		float ComputeMeanBlobArea()
		{
			float sum = 0;
			int n = 0;
			for (int i = 0; i < 4; i++)
			{
				if (CheckIsBlobSet(i))
				{
					n++;
					sum += blobs[i].imgBlob.area;
				}
			}
			if (n > 0)
			{
				return sum / n;
			}
			else
			{
				return 0.0f;
			}
		}
		int GetMeanBlobSum()
		{
			int sum = 0;
			int n = 0;
			for (int i = 0; i < 4; i++)
			{
				if (CheckIsBlobSet(i))
				{
					n++;
					sum += blobs[i].sum;
				}
			}
			if (n > 0)
			{
				return sum / n;
			}
			else
			{
				return 0;
			}
		}
		float ComputeBlobRoundnessScore()
		{
			vector<float> scores;
			for (int i = 0; i < 4; i++)
			{
				if (CheckIsBlobSet(i))
				{
					scores.push_back(blobs[i].roundScore);
				}
			}
			return meanv(scores);
		}
		float ComputeBlobBrightnessScore()
		{
			vector<float> areas;
			vector<float> maxes;
			for (int i = 0; i < 4; i++)
			{
				if (CheckIsBlobSet(i))
				{
					areas.push_back(blobs[i].imgBlob.area);
					maxes.push_back(blobs[i].max);
				}
			}
			if (areas.size() < 3) return 0.0f;
			float areaScore = rangeScore(0.0f, 5.0f, meanv(areas));
			float maxScore = rangeScore(0.0f, 100.0f, meanv(maxes));
			return geometricMean(areaScore, maxScore);
		}
		float ComputeVelocityPix()
		{
			pair<float, float> dp = GetVector();
			return pyth(dp.first, dp.second);
		}
		float GetFrameCountScore()
		{
			float sum = GetBlobCount();
			for (int i = 0; i < 4; i++)
			{
				sum += emBlobScores[i];
			}
			return sum;
		}
		void UpdateFeatureFields()
		{
			ComputePrimaryBlobStats();
			this->velocity = ComputeVelocityPix();
			this->brightness = ComputeBlobBrightnessScore();
			this->roundness = ComputeBlobRoundnessScore();
			this->nBlobs = GetBlobCount();
			this->nfc = GetFrameCountScore();
			featureFieldsSet = true;
		}
#pragma endregion
#pragma region Dead
#pragma endregion
#pragma region Csv : Feature Csv
		std::string ToFeatureCsvHeaderString(bool withImageSetId = true)
		{
			vector<Feature> features = GetFeatures();
			char tmp[4096];
			sprintf(tmp, "index, ");
			if (withImageSetId)
			{
				strcat(tmp, "imageSet, "); 
			}
			for (int i = 0; i < features.size(); i++)
			{
				strcat(tmp, features[i].name.c_str());
				strcat(tmp, ", ");
			}
			strcat(tmp, "class");
			return (string)tmp;
		}
		std::string ToFeatureCsvString(bool withImageSetId = true)
		{
			char tmp[4096];
			char tmp1[1024];
			vector<Feature> features = GetFeatures();
			sprintf(tmp, "%d, ", index);
			if (withImageSetId)
			{
				sprintf(tmp1, "%s, ", imageSetId.c_str()); 
				strcat(tmp, tmp1);
			}
			for (int i = 0; i < features.size(); i++)
			{
				sprintf(tmp1, "%f, ", features[i].value);
				strcat(tmp, tmp1);
			}
			sprintf(tmp1, "%d", classIndex);
			strcat(tmp, tmp1);
			return (string)tmp;
		}
		static void WriteFeatureCsv(std::string& path, vector<CA>& cas, bool withImageSetId = true)
		{
			if (cas.size() == 0) return;
			FILE* csvp = fopen(path.c_str(), "w");
			fprintf(csvp, "%s\n", cas[0].ToFeatureCsvHeaderString(withImageSetId).c_str());
			for (int i = 0; i < cas.size(); i++)
			{
				fprintf(csvp, "%s\n", cas[i].ToFeatureCsvString(withImageSetId).c_str());
			}
			fclose(csvp);
		}
		void AppendToFeatureCsv(FILE* fp)
		{
			fprintf(fp, "%s\n", this->ToFeatureCsvString().c_str());
		}
		void ParseFeatureCsvWords(vector<string>& words)
		{
			int i = 0;
			index = atoi(words[i++].c_str());
			if (words.size() == CAN_N_CSV_COLS)
			{
				imageSetId = words[i++];
			}
			else if (words.size() != CAN_N_CSV_COLS - 1)
			{
				ErrorExit("Wrong word count in csv.");
			}
			nBlobs = atoi(words[i++].c_str());
			nfc = atof(words[i++].c_str());
			blobAreaStats.mean = atof(words[i++].c_str());
			blobMaxStats.mean = atof(words[i++].c_str());
			blobSumStats.mean = atof(words[i++].c_str());
			blobAreaStats.sigmaPct = atof(words[i++].c_str());
			blobMaxStats.sigmaPct = atof(words[i++].c_str());
			blobSumStats.sigmaPct = atof(words[i++].c_str());
			simpleScore = atof(words[i++].c_str());
			sboProxScore = atof(words[i++].c_str());
			caProxScore = atof(words[i++].c_str());
			salienceScore = atof(words[i++].c_str());
			caBsmScore = atof(words[i++].c_str());
			lineFitRmsError = atof(words[i++].c_str());
			velocity = atof(words[i++].c_str());
			brightness = atof(words[i++].c_str());
			roundness = atof(words[i++].c_str());
			classIndex = atoi(words[i++].c_str());
			corrDetNumber = (classIndex == 1 ? 9999999 : -1);
			featureFieldsSet = true;
		}
		static void LoadFeatureCsv(std::string& path, vector<CA>& cs)
		{
			vector<string> lines;
			logd.printf("Loading csv file...\n");
			if (!FileUtil::SnarfFile(path.c_str(), lines))
			{
				ErrorExit("Failed to load CSV.");
			}
			int n = lines.size();
			logd.printf("Processing %d lines.\n", n);
			for (int i = 1; i < n; i++)
			{
				vector<string> words = StringUtils::split(lines[i], ',');
				CA c;
				c.ParseFeatureCsvWords(words);
				cs.push_back(c);
			}
		}
#pragma endregion
#pragma region Features
		static const int CAN_N_FEATURES = 17; 
		static const int CAN_N_CSV_COLS = CAN_N_FEATURES + 3; 
		vector<Feature> GetFeatures()
		{
			vector<Feature> fs;
			if (!featureFieldsSet)
			{
				ErrorExit("Feature fields not set.");
				return fs;
			}
			int nBlobs = this->nBlobs;
			float nfc = this->nfc;
			float areaMean = blobAreaStats.mean;
			float maxMean = blobMaxStats.mean;
			float sumMean = blobSumStats.mean;
			float areaSigmaPct = blobAreaStats.sigmaPct;
			float maxSigmaPct = blobMaxStats.sigmaPct;
			float sumSigmaPct = blobSumStats.sigmaPct;
			float simpleScore = this->simpleScore;
			float sboProxScore = this->sboProxScore;
			float caProxScore = this->caProxScore;
			float salienceScore = this->salienceScore;
			float caBsmScore = this->caBsmScore;
			float lineFitRmsError = this->lineFitRmsError;
			float velocity = this->velocity;
			float brightness = this->brightness;
			float roundness = this->roundness;
			fs.push_back(Feature("nBlobs", nBlobs));
			fs.push_back(Feature("nfc", nfc));
			fs.push_back(Feature("areaMean", areaMean));
			fs.push_back(Feature("maxMean", maxMean));
			fs.push_back(Feature("sumMean", sumMean));
			fs.push_back(Feature("areaSigmaPct", areaSigmaPct));
			fs.push_back(Feature("maxSigmaPct", maxSigmaPct));
			fs.push_back(Feature("sumSigmaPct", sumSigmaPct));
			fs.push_back(Feature("simpleScore", simpleScore));
			fs.push_back(Feature("sboProxScore", sboProxScore));
			fs.push_back(Feature("caProxScore", caProxScore));
			fs.push_back(Feature("salienceScore", salienceScore));
			fs.push_back(Feature("caBsmScore", caBsmScore));
			fs.push_back(Feature("lineFitRmsError", lineFitRmsError));
			fs.push_back(Feature("velocity", velocity));
			fs.push_back(Feature("brightness", brightness));
			fs.push_back(Feature("roundness", roundness));
			if (fs.size() != CAN_N_FEATURES)
			{
				ErrorExit("Number of features is wrong, see CAN_N_FEATURES");
			}
			for (int i = 0; i < fs.size(); i++)
			{
				fs[i].isSelected = true;
			}
			return fs;
		}
#pragma endregion
};
class CAList
{
private:
	vector<unordered_map<int, vector<int>>> refMap;
	void AddRefMapEntry(CA& ca, int casIndex, bool debug)
	{
		for (int frameIndex = 0; frameIndex < 4; frameIndex++)
		{
			int blobIdx = ca.blobIndices[frameIndex];
			if (blobIdx >= 0)
			{
				refMap[frameIndex][blobIdx].push_back(casIndex);
				if (debug) logd.debug("CAList: Add: Fr %d, blobIdx %d -> CA %d at index %d", frameIndex, blobIdx, ca.index, casIndex);
			}
		}
	}
	bool FindRefMapEntry(CA& ca, bool debug)
	{
		bool foundASetBlobIndex = false;
		for (int frameIndex = 0; frameIndex < 4; frameIndex++)
		{
			int blobIdx = ca.blobIndices[frameIndex];
			if (blobIdx >= 0)
			{
				foundASetBlobIndex = true;
				if (refMap[frameIndex].find(blobIdx) != refMap[frameIndex].end())
				{
					for (int i = 0; i < refMap[frameIndex][blobIdx].size(); i++) 
					{
						int otherCasIndex = refMap[frameIndex][blobIdx][i];
						if (ca.CheckBlobIndicesEqual(cas[otherCasIndex]))
						{
							if (debug) logd.debug("CAList: FindRefMapEntry: CA %d has match %d already in ref map referring to same blob indices.", ca.index, cas[otherCasIndex].index);
							return true;
						}
					}
				}
			}
		}
		if (!foundASetBlobIndex)
		{
			ErrorExit("Don't call this with a CA that has no blobs.");
			return true; 
		}
		if (debug) logd.debug("CAList: FindRefMapEntry: CA %d has no match already in ref map.", ca.index);
		return false;
	}
public:
	vector<CA> cas; 
	CAList()
	{
		refMap.resize(4);
	}
	int size()
	{
		return cas.size();
	}
	vector<int> Find(CA& ca)
	{
		vector<int> list;
		int n = cas.size();
		for (int i = 0; i < n; i++)
		{
			if (cas[i].ComputeSimScore(ca) > 0.5f)
			{
				list.push_back(i);
			}
		}
		return list;
	}
	int FindIndexByIndex(int caIdx)
	{
		for (int i = 0; i < cas.size(); i++)
		{
			if (cas[i].index == caIdx) return i;
		}
		return -1;
	}
	void Append(vector<CA>& v)
	{
		cas.insert(cas.end(), v.begin(), v.end());
	}
	void AccumTopNPClass1(vector<CA>& v, int nMax)
	{
		Append(v);
		SortOnPClass1();
		if (cas.size() > nMax)
		{
			cas.resize(nMax);
		}
	}
	void AddSimple(CA& ca)
	{
		cas.push_back(ca);
	}
	bool Add(CA& ca, DebugSpec dspec)
	{
		bool debug = ca.CheckDebug(dspec);
		if (debug)
		{
			logd.debug("CAList: Add: Adding CA %d with blob indices %d, %d, %d, %d and emscores %.1f, %.1f, %.1f, %.1f", 
				ca.index, 
				ca.blobIndices[0], ca.blobIndices[1], ca.blobIndices[2], ca.blobIndices[3],
				ca.emBlobScores[0], ca.emBlobScores[1], ca.emBlobScores[2], ca.emBlobScores[3]
				);
		}
		if (FindRefMapEntry(ca, debug))
		{
			return false;
		}
		int thisIndex = cas.size();
		cas.push_back(ca);
		AddRefMapEntry(ca, thisIndex, debug);
		return true;
	}
	void DropWeakest(int nToKeep, bool sortOnPClass1, DebugSpec dspec)
	{
		if (sortOnPClass1)
		{
			SortOnPClass1();
		}
		else
		{
			Sort();
		}
		if (size() > nToKeep)
		{
			if (dspec.debugDetNumber >= 0)
			{
				for (int i = nToKeep; i < size(); i++)
				{
					if (cas[i].corrDetNumber == dspec.debugDetNumber)
					{
						logd.warn("DropWeakest: Dropping CA with det number %d at index %d (over limit %d) on %s %.3f",
							cas[i].corrDetNumber, i, nToKeep, (sortOnPClass1 ? "pClass1" : "simpleScore"),
							(sortOnPClass1 ? cas[i].pClass1 : cas[i].simpleScore));
					}
				}
			}
			cas.resize(nToKeep);
		}
	}
	struct CASimpleScoreGreaterThan
	{
		inline bool operator() (const CA& c1, const CA& c2)
		{
			return (c1.simpleScore > c2.simpleScore);
		}
	};
	static void Sort(vector<CA>& v)
	{
		std::sort(v.begin(), v.end(), CASimpleScoreGreaterThan());
	}
	void Sort()
	{
		Sort(this->cas);
	}
	struct CAPClass1GreaterThan
	{
		inline bool operator() (const CA& c1, const CA& c2)
		{
			if ((c1.pClass1 < 0) || (c2.pClass1 < 0)) ErrorExit("Comparing a CA pClass1 that is not set.");
			return (c1.pClass1 > c2.pClass1);
		}
	};
	static void SortOnPClass1(vector<CA>& v)
	{
		std::sort(v.begin(), v.end(), CAPClass1GreaterThan());
	}
	void SortOnPClass1()
	{
		SortOnPClass1(this->cas);
	}
	struct CANeoGreaterThan
	{
		inline bool operator() (const CA& c1, const CA& c2)
		{
			if ((c1.pNeo < 0) || (c2.pNeo < 0)) ErrorExit("Comparing a CA pNeo that is not set.");
			return (c1.pNeo > c2.pNeo);
		}
	};
	static void SortOnPNeo(vector<CA>& v)
	{
		std::sort(v.begin(), v.end(), CANeoGreaterThan());
	}
	void SortOnPNeo()
	{
		SortOnPNeo(this->cas);
	}
	struct CAOverlapGreaterThan
	{
		inline bool operator() (const CA& c1, const CA& c2)
		{
			if ((c1.pOverlap < 0) || (c2.pOverlap < 0)) ErrorExit("Comparing a CA pOverlap that is not set.");
			return (c1.pOverlap > c2.pOverlap);
		}
	};
	static void SortOnPOverlap(vector<CA>& v)
	{
		std::sort(v.begin(), v.end(), CAOverlapGreaterThan());
	}
	void SortOnPOverlap()
	{
		SortOnPOverlap(this->cas);
	}
	struct CAVelocityGreaterThan
	{
		inline bool operator() (const CA& c1, const CA& c2)
		{
			if ((c1.velocity < 0) || (c2.velocity < 0)) ErrorExit("Comparing a CA velocity that is not set.");
			return (c1.velocity > c2.velocity);
		}
	};
	static void SortOnVelocityDescending(vector<CA>& v)
	{
		std::sort(v.begin(), v.end(), CAVelocityGreaterThan());
	}
	void ApplyPostProcessingScores(DebugSpec dspec)
	{
		for (int i = 0; i < cas.size(); i++)
		{
			bool debug = cas[i].CheckDebug(dspec);
			cas[i].ApplyPostProcessingScores(debug);
		}
		Sort();
	}
	void UpdateSimpleScores()
	{
		int n = cas.size();
		for (int i = 0; i < n; i++)
		{
			cas[i].UpdateSimpleScore();
		}
	}
	void UpdateNCAsCount()
	{
		int n = cas.size();
		for (int i = 0; i < n; i++)
		{
			cas[i].nCAsThisImageSet = n;
		}
	}
	void Display()
	{
		logd.printf("CAList: Size: %d\n", cas.size());
		int n = cas.size();
		n = MIN(10, n);
		for (int i = 0; i < n; i++)
		{
			logd.printf("    %d: %s\n", i, cas[i].ToString().c_str());
		}
	}
	void DisplayCA(int caIndex)
	{
		int n = cas.size();
		for (int i = 0; i < n; i++)
		{
			if (cas[i].index == caIndex)
			{
				cas[i].Display();
				break;
			}
		}
	}
	std::string ToString()
	{
		char tmp[512];
		sprintf(tmp, "Count: %d", cas.size());
		return (string)tmp;
	}
	void GetSamples(vector<vector<float>>& samples, vector<int>& truthClasses)
	{
		GetSamples(this->cas, samples, truthClasses);
	}
	static void GetSamples(vector<CA>& cas, vector<vector<float>>& samples, vector<int>& truthClasses)
	{
		int n = cas.size();
		samples.reserve(n);
		truthClasses.reserve(n);
		for (int i = 0; i < n; i++)
		{
			vector<Feature> features = cas[i].GetFeatures();
			samples.push_back(Feature::GetAllValues(features));
			truthClasses.push_back(cas[i].GetClassIndex());
		}
	}
	static int FindCAMatchedWithDet(vector<CA>& cas, int detNumber)
	{
		int n = cas.size();
		for (int i = 0; i < n; i++)
		{
			if (cas[i].corrDetNumber == detNumber) return i;
		}
		return -1;
	}
	struct CAPFinalGreaterThan
	{
		inline bool operator() (const CA& c1, const CA& c2)
		{
			if ((c1.pFinal < 0) || (c2.pFinal < 0)) ErrorExit("Comparing a CA pFinal that is not set.");
			return (c1.pFinal > c2.pFinal);
		}
	};
	static void SortOnPFinal(vector<CA>& v)
	{
		std::sort(v.begin(), v.end(), CAPFinalGreaterThan());
	}
	static void DoFinalSort(vector<CA>& cas)
	{
		int n = cas.size();
		for (int i = 0; i < n; i++)
		{
			cas[i].pFinal = cas[i].pClass1;
			if (cas[i].isNeo)
			{
				cas[i].pFinal = cas[i].pClass1 + cas[i].pNeo / 4;
			}
		}
		SortOnPFinal(cas);
	}
};
#define _USE_MATH_DEFINES
#include <cmath>
static int privateUniqueId = 0;
class ImageSet
{
public:
	int index;
	string imageSetId;
	vector<FitsImage2> images; 
	vector<NsiBlobSet> blobSets; 
	vector<string> detStrings; 
	vector<DetectionRecord> origDetRecords; 
	vector<DetectionRecord> detRecords; 
	vector<Detection> detectionsWithOverlaps; 
	vector<Detection> detectionsWithoutOverlaps; 
	vector<pair<float, float>> registerShifts; 
	vector<Transform2d> fwdXforms; 
	vector<Transform2d> revXforms; 
	Img nsiCM;
private:
	void ClearDerivedDets()
	{
		detRecords.clear();
		detectionsWithOverlaps.clear();
		detectionsWithoutOverlaps.clear();
	}
public:
	ImageSet()
	{
		this->imageSetId = "";
		images.resize(4);
		registerShifts.resize(4);
		blobSets.resize(4);
	}
	int getWidth()
	{
		if (images.size() > 0)
		{
			return images[0].getWidth();
		}
		else
		{
			return 0;
		}
	}
	int getHeight()
	{
		if (images.size() > 0)
		{
			return images[0].getHeight();
		}
		else
		{
			return 0;
		}
	}
	void UpdateDimsFromImages()
	{
		for (int i = 0; i < images.size(); i++)
		{
			images[i].UpdateDimsFromOrigImage();
		}
	}
	void ParseDets()
	{
		origDetRecords.resize(detStrings.size());
		for (int i = 0; i < detStrings.size(); i++)
		{
			DetectionRecord rec;
			if (!rec.TryParse(detStrings[i], rec))
			{
				ErrorExit("DetectionRecord parse failed.");
			}
			rec.uniqueId = privateUniqueId++;
			origDetRecords[i] = rec;
		}
		ClearDerivedDets();
		RebuildDownstreamDets();
	}
	void RebuildDownstreamDets()
	{
		detectionsWithOverlaps = Detection::Build(imageSetId, origDetRecords);
		unordered_map<int, vector<int>> keeperMap; 
		DebugSpec dspec;
		RemoveOverlappingDets(keeperMap, dspec); 
		detectionsWithoutOverlaps = Detection::Build(imageSetId, detRecords);
	}
	int GetOverlapDetsCount()
	{
		return (origDetRecords.size() - detRecords.size()) / 4;
	}
	void Display()
	{
		logd.printf("ImageSet: Index: %d, ID: %s, images: %d (%dx%d), origDets: %d, Dets: %d\n", index, imageSetId.c_str(), images.size(),
			images.size() > 0 ? images[0].getWidth() : 0,
			images.size() > 0 ? images[0].getHeight() : 0,
			origDetRecords.size() / 4, detRecords.size() / 4);
		for (int i = 0; i < images.size(); i++)
		{
			int cx = images[i].getWidth() / 2 - 1;
			int cy = images[i].getHeight() / 2 - 1;
			logd.printf("    %d: Origin: %s, Center: %s, RegisterShift: (%.3f, %.3f)\n", i, images[i].PixelToWorld(0, 0).ToString().c_str(), 
				images[i].PixelToWorld(cx, cy).ToString().c_str(), registerShifts[i].first, registerShifts[i].second);
		}
	}
	void SaveOrigImages()
	{
		char tmp[64];
		for (int i = 0; i < images.size(); i++)
		{
			sprintf(tmp, "imageSet_%03d_%d", index, i);
			SaveDebugTiff(tmp, images[i].origImage);
		}
	}
	void CheckOrigDetsXY(DebugSpec dspec)
	{
	}
	void CropImages()
	{
		for (int i = 0; i < images.size(); i++)
		{
			images[i].CropImage();
		}
	}
	void AdjustDetsForCrop()
	{
		int n = origDetRecords.size();
		for (int i = 0; i < n; i++)
		{
			int ox, oy;
			ox = images[origDetRecords[i].frameIndex].cropOffset.first;
			oy = images[origDetRecords[i].frameIndex].cropOffset.second;
			origDetRecords[i].x -= ox;
			origDetRecords[i].y -= oy;
		}
		RebuildDownstreamDets();
	}
	void UpdateDetsForRegisterShift()
	{
		int n = origDetRecords.size();
		for (int i = 0; i < n; i++)
		{
			int ox, oy;
			ox = roundf(registerShifts[origDetRecords[i].frameIndex].first);
			oy = roundf(registerShifts[origDetRecords[i].frameIndex].second);
			origDetRecords[i].xReg = origDetRecords[i].x + ox;
			origDetRecords[i].yReg = origDetRecords[i].y + oy;
		}
		RebuildDownstreamDets();
	}
	void RecomputeXYInOrigDetRecords(bool doEmitWarnings = true)
	{
		int nd = origDetRecords.size();
		for (int i = 0; i < 4; i++)
		{
			for (int j = i; j < nd; j += 4)
			{
				float ox, oy;
				ConvertWorldToPixel(origDetRecords[j].raDec, i, false, ox, oy);
				float dx = ox - origDetRecords[j].x;
				float dy = oy - origDetRecords[j].y;
				float dist = sqrtf(dx*dx + dy*dy);
				if ((j < 10) && (dist > 3))
				{
					if (doEmitWarnings) logd.printf("Warn: RecomputeXYInOrigDetRecords: Det: %d, wcslib off by %.1f pixels compared to orig det truth pixel points.\n", origDetRecords[j].detectionNumber, dist);
				}
				int cx, cy;
				cx = images[i].cropOffset.first;
				cy = images[i].cropOffset.second;
				origDetRecords[j].x = ox - cx;
				origDetRecords[j].y = oy - cy;
			}
		}
		UpdateDetsForRegisterShift();
	}
	void SaveRoiImage(const char* tag, IRoi roi, Img& image, int dim, int frameIndex)
	{
		if (image.getDepthInBits() == 16)
		{
			Img i1(dim, dim, 16);
			ImgUtil::CopyRoi(image, roi, i1);
			SaveDebugImage(tag, i1);
		}
		else
		{
			Img i2(dim, dim, 8);
			ImgUtil::CopyRoi(image, roi, i2);
			SaveDebugImage(tag, i2);
		}
	}
	void SaveDetRoiImage(const char* tag, int origDetRecordIndex, Img& image, bool useRegCoords)
	{
		char tmp[256];
		int dim = 128;
		int frameIndex = origDetRecords[origDetRecordIndex].frameIndex;
		int x0, y0;
		if (useRegCoords)
		{
			x0 = roundf(origDetRecords[origDetRecordIndex].xReg - dim / 2);
			y0 = roundf(origDetRecords[origDetRecordIndex].yReg - dim / 2);
		}
		else
		{
			x0 = roundf(origDetRecords[origDetRecordIndex].x - dim / 2);
			y0 = roundf(origDetRecords[origDetRecordIndex].y - dim / 2);
		}
		IRoi roi(x0, y0, x0 + dim, y0 + dim);
		if (image.getDepthInBits() == 16)
		{
			Img i1(dim, dim, 16);
			ImgUtil::CopyRoi(image, roi, i1);
			SaveDebugImage(this->index, frameIndex, tag, origDetRecords[origDetRecordIndex], i1);
		}
		else
		{
			Img i2(dim, dim, 8);
			ImgUtil::CopyRoi(image, roi, i2);
			SaveDebugImage(this->index, frameIndex, tag, origDetRecords[origDetRecordIndex], i2);
		}
	}
	void SaveDetRoiImage(const char* tag, int origDetRecordIndex, int whichImageToSave = 0)
	{
		char tmp[256];
		int dim = 128;
		Img i1(dim, dim, 16);
		Img i2(dim, dim, 8);
		int frameIndex = origDetRecords[origDetRecordIndex].frameIndex;
		int x0Reg = roundf(origDetRecords[origDetRecordIndex].xReg - dim / 2);
		int y0Reg = roundf(origDetRecords[origDetRecordIndex].yReg - dim / 2);
		IRoi regRoi(x0Reg, y0Reg, x0Reg + dim, y0Reg + dim);
		if (whichImageToSave == 0) 
		{
			sprintf(tmp, "trueDetOrig_%s", tag);
			SaveDetRoiImage(tmp, origDetRecordIndex, images[frameIndex].origImage, false);
		}
		else if (whichImageToSave == 1) 
		{
			sprintf(tmp, "trueDetNsiDiff_%s", tag);
			SaveDetRoiImage(tmp, origDetRecordIndex, images[frameIndex].nsiDiff, true);
		}
		else if (whichImageToSave == 2) 
		{
			sprintf(tmp, "trueDetNsi_%s", tag);
			SaveDetRoiImage(tmp, origDetRecordIndex, images[frameIndex].nsi, true);
		}
		else if (whichImageToSave == 3) 
		{
			sprintf(tmp, "trueDetFixedImg_%s", tag);
			SaveDetRoiImage(tmp, origDetRecordIndex, images[frameIndex].fixedImage, true);
		}
	}
	void SaveDetRoiImages(const char* tag, int detectionNumber, int whichImageToSave = 0)
	{
		int n = origDetRecords.size();
		for (int i = 0; i < n; i++)
		{
			if ((detectionNumber == -2) || (origDetRecords[i].detectionNumber == detectionNumber))
			{
				SaveDetRoiImage(tag, i, whichImageToSave);
			}
		}
	}
	void SaveFnNeoImages(vector<CA>& cas)
	{
		int n = origDetRecords.size();
		int whichImageToSave = 2;
		for (int i = 0; i < n; i++)
		{
			int dn = origDetRecords[i].detectionNumber;
			int frameIndex = origDetRecords[i].frameIndex;
			bool isMatchedWithCAThatSurvivedCulling = (CAList::FindCAMatchedWithDet(cas, dn) >= 0);
			if (origDetRecords[i].isNeo)
			{
				if (!isMatchedWithCAThatSurvivedCulling)
				{
					char tag[512];
					sprintf(tag, "unmatchedNeo_det%03d_fr%d", dn, frameIndex);
					SaveDetRoiImage(tag, i, whichImageToSave);
					logd.printf("**** Unmatched neo: %s det%03d fr%d is (%d, %d)\n", imageSetId.c_str(), dn, frameIndex, roundf(origDetRecords[i].xReg), roundf(origDetRecords[i].yReg));
				}
				else
				{
					logd.printf("**** Matched neo: %s det%03d fr%d is (%d, %d)\n", imageSetId.c_str(), dn, frameIndex, roundf(origDetRecords[i].xReg), roundf(origDetRecords[i].yReg));
				}
			}
		}
	}
	void SaveDetVizImage(int frameIndex = -1, int detNumber = -1, bool drawNum = false)
	{
		char tmp[64];
		int radius = 8;
		Img viz(images[0].getWidth(), images[1].getHeight(), 16);
		Img text(images[0].getWidth(), images[1].getHeight(), 8);
		Img8Util::Clear(text);
		int textRand = 6;
		int textSize = 12;
		int textMargin = textRand;
		if (images[0].origImage.getWidth() > 0)
		{
			for (int i = 0; i < images.size(); i++)
			{
				if ((frameIndex < 0) || (frameIndex == i))
				{
					ImgUtil::Copy(images[i].origImage, viz);
					Img8Util::Clear(text);
					for (int j = 0; j < origDetRecords.size(); j++)
					{
						if ((origDetRecords[j].frameIndex == i) && ((detNumber < 0) || (detNumber == origDetRecords[j].detectionNumber)))
						{
							int x = roundf(origDetRecords[j].x); 
							int y = roundf(origDetRecords[j].y); 
							ImgUtil::DrawBox(viz, x - radius, y - radius, x + radius, y + radius, 65535);
							sprintf(tmp, "%03d", origDetRecords[j].detectionNumber);
							if (drawNum)
							{
								int tx = x - roundf(textSize * 1.5) + (rand() % textRand);
								int ty = y - radius - textSize - textMargin + (rand() % textRand);
								DrawText8(text, tx, ty, tmp, 12, false);
							}
						}
					}
					ImgUtil::MaskSaturate(viz, text);
					sprintf(tmp, "imageSet_%03d_%d_detviz", index, i);
					SaveDebugImage(tmp, viz);
				}
			}
		}
		for (int i = 0; i < images.size(); i++)
		{
			if ((frameIndex < 0) || (frameIndex == i))
			{
				vector<tuple<int, int, uint32_t, string>> nsiRois;
				for (int j = 0; j < origDetRecords.size(); j++)
				{
					if ((origDetRecords[j].frameIndex == i) && ((detNumber < 0) || (detNumber == origDetRecords[j].detectionNumber)))
					{
						sprintf(tmp, "%03d", origDetRecords[j].detectionNumber);
						nsiRois.push_back(tuple<int, int, uint32_t, string>(roundf(origDetRecords[j].xReg), roundf(origDetRecords[j].yReg), 0xFFFFFFFF, (string)tmp));
					}
				}
				sprintf(tmp, "imageSet_%03d_%d_detvizNsi", index, i);
				Img8Util::Copy(images[i].nsiDiff, text);
				RenderColoredRois(text, tmp, 8, nsiRois);
			}
		}
	}
	int CountNeosInOrigDetRecords()
	{
		int n = origDetRecords.size();
		int count = 0;
		for (int i = 0; i < n; i+=4)
		{
			if (origDetRecords[i].isNeo) count++;
		}
		return count;
	}
	vector<int> FindAllOverlappingDetsLeadIndices(int thisLeadIndex, vector<bool>& dups, bool debug)
	{
		float distSqThreshold = (REMOVE_DET_TRUTH_DIST_THRESHOLD * REMOVE_DET_TRUTH_DIST_THRESHOLD);
		int n = origDetRecords.size();
		vector<int> otherDetOverlapLeadIndices;
		for (int i = thisLeadIndex + 4; i < n; i += 4)
		{
			float maxDist = 0;
			for (int k = 0; k < 4; k++)
			{
				float dx = origDetRecords[i + k].x - origDetRecords[thisLeadIndex + k].x;
				float dy = origDetRecords[i + k].y - origDetRecords[thisLeadIndex + k].y;
				float distSq = dx*dx + dy*dy;
				if (distSq > maxDist)
				{
					maxDist = distSq;
				}
				origDetRecords[i + k].dupCount = 0;
			}
			if (maxDist < distSqThreshold)
			{
				if (!dups[i])
				{
					if (debug) logd.printf("    Det %d overlaps current %d with dist %.1f\n", origDetRecords[i].detectionNumber, 
						origDetRecords[thisLeadIndex].detectionNumber, sqrtf(maxDist));
					otherDetOverlapLeadIndices.push_back(i);
				}
				else
				{
					if (debug) logd.printf("    Det %d overlaps current %d with dist %.1f, but it is already marked as a dup.\n", origDetRecords[i].detectionNumber,
						origDetRecords[thisLeadIndex].detectionNumber, sqrtf(maxDist));
				}
			}
		}
		return otherDetOverlapLeadIndices;
	}
	void RemoveOverlappingDets(unordered_map<int, vector<int>>& keeperMap, DebugSpec dspec, bool debugAll = false)
	{
		int n = origDetRecords.size();
		float distSqThreshold = (REMOVE_DET_TRUTH_DIST_THRESHOLD * REMOVE_DET_TRUTH_DIST_THRESHOLD);
		vector<size_t> overlapLeadIndices; 
		vector<bool> dups(n); 
		for (int i = 0; i < n; i++) dups[i] = false;
		for (int i = 0; i < n; i += 4)
		{
			bool debugThis = dspec.CheckDebugDet(index, origDetRecords[i].detectionNumber);
			if (dups[i])
			{
				if (debugAll || debugThis)
					logd.printf("For det %d it is already marked as a dup, continuing\n", origDetRecords[i].detectionNumber);
				continue;
			}
			vector<int> otherDetOverlapLeadIndices = FindAllOverlappingDetsLeadIndices(i, dups, debugAll || debugThis);
			int no = otherDetOverlapLeadIndices.size();
			if (debugAll || debugThis) logd.printf("For det %d there are %d overlaps\n", origDetRecords[i].detectionNumber, no);
			if (no > 0)
			{
				vector<int> dupLeadIndicesToMark;
				if (origDetRecords[i].isNeo)
				{
					if (debugAll || debugThis) logd.printf("    For det %d: This is a neo\n", origDetRecords[i].detectionNumber);
					dupLeadIndicesToMark = otherDetOverlapLeadIndices;
					origDetRecords[i].hasOverlap = true;
					origDetRecords[i].overlapDetNum = origDetRecords[otherDetOverlapLeadIndices[0]].detectionNumber;
				}
				else
				{
					if (debugAll || debugThis) logd.printf("    For det %d: This is not a neo\n", origDetRecords[i].detectionNumber);
					int neoIndex = -1;
					for (int j = 0; j < otherDetOverlapLeadIndices.size(); j++)
					{
						if (origDetRecords[otherDetOverlapLeadIndices[j]].isNeo)
						{
							neoIndex = otherDetOverlapLeadIndices[j];
						}
					}
					if (neoIndex < 0)
					{
						if (debugAll || debugThis) logd.printf("    For det %d: No other neo's\n", origDetRecords[i].detectionNumber);
						dupLeadIndicesToMark = otherDetOverlapLeadIndices;
						origDetRecords[i].hasOverlap = true;
						origDetRecords[i].overlapDetNum = origDetRecords[otherDetOverlapLeadIndices[0]].detectionNumber;
					}
					else
					{
						if (debugAll || debugThis) logd.printf("    For det %d: This is not a neo, but %d is a neo\n", origDetRecords[i].detectionNumber, origDetRecords[neoIndex].detectionNumber);
						dupLeadIndicesToMark.push_back(i);
						origDetRecords[neoIndex].hasOverlap = true;
						origDetRecords[neoIndex].overlapDetNum = origDetRecords[i].detectionNumber;
						for (int j = 0; j < otherDetOverlapLeadIndices.size(); j++)
						{
							int idx = otherDetOverlapLeadIndices[j];
							if (idx != neoIndex)
							{
								dupLeadIndicesToMark.push_back(idx);
							}
						}
					}
				}
				for (int j = 0; j < dupLeadIndicesToMark.size(); j++)
				{
					int idx = dupLeadIndicesToMark[j];
					if (debugAll || debugThis) logd.printf("    Det %d chosen as an overlap when evaluating %d\n", origDetRecords[idx].detectionNumber, origDetRecords[i].detectionNumber);
					overlapLeadIndices.push_back(idx);
					dups[idx] = true;
					keeperMap[i].push_back(idx);
				}
			}
		}
		vector<size_t> overlapIndices; 
		for (int i = 0; i < overlapLeadIndices.size(); i++)
		{
			if (debugAll) logd.printf("    %d chosen as an overlap, removing it.\n", origDetRecords[overlapLeadIndices[i]].detectionNumber);
			for (int j = 0; j < 4; j++)
			{
				overlapIndices.push_back(overlapLeadIndices[i] + j);
			}
		}
		int nto = overlapIndices.size();
		std::sort(overlapIndices.begin(), overlapIndices.end());
		auto it = std::unique(overlapIndices.begin(), overlapIndices.end());
		overlapIndices.resize(std::distance(overlapIndices.begin(), it));
		if (overlapIndices.size() != nto) ErrorExit("Got dup overlap indices.");
		detRecords = erase_indices(origDetRecords, overlapIndices);
		if (debugAll) logd.printf("Remove overlapping dets took it from %d to %d det sets.\n", origDetRecords.size() / 4, detRecords.size() / 4);
	}
	RaDec PixelToWorld(float x, float y, int frameIndex, bool isRegistrationShifted)
	{
		if (isRegistrationShifted)
		{
			x -= registerShifts[frameIndex].first;
			y -= registerShifts[frameIndex].second;
		}
#ifndef USE_WCSLIB
		RaDec rd = images[frameIndex].PixelToWorld(x, y);
#endif
		return rd;
	}
	pair<float, float> ConvertWorldToPixel(RaDec rd, int frameIndex, bool doRegistrationShift)
	{
#ifndef USE_WCSLIB
		float x, y;
		images[frameIndex].WorldToPixel(rd, x, y);
#endif
		if (doRegistrationShift)
		{
			x += registerShifts[frameIndex].first;
			y += registerShifts[frameIndex].second;
		}
		return pair<float, float>((float)x, (float)y);
	}
	void ConvertWorldToPixel(RaDec rd, int frameIndex, bool doRegistrationShift, float& x, float& y)
	{
		pair<float, float> pt = ConvertWorldToPixel(rd, frameIndex, doRegistrationShift);
		x = pt.first;
		y = pt.second;
	}
	vector<RaDec> ConvertPointsToRaDec(vector<pair<float, float>>& points, bool isRegistrationShifted)
	{
		vector<RaDec> radecs;
		if (points.size() != 4) ErrorExit("Expect 4 points to find matching det set.");
		for (int i = 0; i < 4; i++)
		{
			float px = points[i].first;
			float py = points[i].second;
			RaDec rd = PixelToWorld(px, py, i, isRegistrationShifted);
			radecs.push_back(rd);
		}
		return radecs;
	}
	struct FindMdsAscending
	{
		inline bool operator() (const pair<int, float>& c1, const pair<int, float>& c2)
		{
			return (c1.second < c2.second);
		}
	};
	vector<pair<int, float>> FindAllMatchingDetSets(bool useNonOverlapSet, vector<pair<float, float>>& points, bool debug)
	{
		if (points.size() != 4) ErrorExit("Expect 4 points to find matching det set.");
		vector<RaDec> radecs = ConvertPointsToRaDec(points, true);
		vector<DetectionRecord>& dets = origDetRecords;
		if (useNonOverlapSet)
		{
			dets = detRecords;
		}
		vector<pair<int, float>> results;
		int n = dets.size();
		float sqThreshold = DET_TRUTH_ALIGN_THRESHOLD_PIX * DET_TRUTH_ALIGN_THRESHOLD_PIX;
		for (int i = 0; i < n; i += 4)
		{
			int matchCount = 0;
			float meanDistSq = 0;
			for (int j = 0; j < 4; j++)
			{
				float dx, dy;
				dx = dets[i + j].xReg - points[j].first;
				dy = dets[i + j].yReg - points[j].second;
				float distSq = dx*dx + dy*dy;
				if (distSq < sqThreshold)
				{
					if (debug) logd.printf("    FindAllMatchingDetSets: Dist from dr %d with %s (%.1f, %.1f) to %s (%.1f, %.1f) is %.5f is under threshold %.5f.\n",
						dets[i + j].detectionNumber,
						dets[i + j].raDec.ToString().c_str(), 
						dets[i + j].xReg, dets[i + j].yReg,
						radecs[j].ToString().c_str(), points[j].first, points[j].second,
						sqrtf(distSq), sqrtf(sqThreshold));
					matchCount++;
					meanDistSq += distSq;
				}
				else
				{
				}
			}
			if (matchCount == 4)
			{
				meanDistSq /= 4;
				results.push_back(pair<int, float>(i, meanDistSq));
			}
		}
		std::sort(results.begin(), results.end(), FindMdsAscending());
		return results;
	}
	int FindMatchingDetSet(bool useNonOverlapSet, vector<pair<float, float>>& points, vector<bool>& matched, int caIndex, bool debug)
	{
		vector<pair<int, float>> sets = FindAllMatchingDetSets(useNonOverlapSet, points, debug);
		if (sets.size() == 0)
		{
			if (debug) logd.printf("FindMatchingDetSet: Found no matching det sets for CA %d: (%.1f, %.1f)\n", caIndex, points[0].first, points[1].second);
			return -1; 
		}
		vector<DetectionRecord>& dets = origDetRecords;
		if (useNonOverlapSet)
		{
			dets = detRecords;
		}
		for (int i = 0; i < sets.size(); i++)
		{
			if (!matched[sets[i].first])
			{
				matched[sets[i].first] = true;
				if (debug) logd.printf("FindMatchingDetSet: Found %d matching det sets (%d, ...) for CA %d: (%.1f, %.1f) and %d was not already matched.\n",
					sets.size(), dets[sets[0].first].detectionNumber, caIndex, points[0].first, points[1].second, dets[sets[i].first].detectionNumber);
				return sets[i].first;
			}
		}
		if (debug) logd.printf("FindMatchingDetSet: Found %d matching det sets (%d, ...) for CA %d: (%.1f, %.1f) but all were already matched.\n",
			sets.size(), dets[sets[0].first].detectionNumber, caIndex, points[0].first, points[1].second);
		return -1;
	}
	vector<pair<float,float>> GetDetsRegisteredPositionDeltas()
	{
		int n = origDetRecords.size();
		vector<pair<float, float>> v;
		for (int i = 0; i < n; i += 4)
		{
			for (int j = 1; j < 4; j++)
			{
				float dx = origDetRecords[i + j].xReg - origDetRecords[i].xReg;
				float dy = origDetRecords[i + j].yReg - origDetRecords[i].yReg;
				v.push_back(pair<float,float>(dx,dy));
			}
		}
		return v;
	}
	void ComputeCARaDecs(CAList& cal, DebugSpec dspec)
	{
		int n = cal.cas.size();
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				pair<float, float> pt = cal.cas[i].GetFittedCoords(j); 
				RaDec rd = PixelToWorld(pt.first, pt.second, j, true);
				cal.cas[i].SetRaDec(j, rd);
			}
		}
	}
	static int GetDetIndexByDetNumber(vector<DetectionRecord>& dets, int detectionNumber, int frameIndex)
	{
		int n = dets.size();
		for (int i = 0; i < n; i++)
		{
			if ((dets[i].detectionNumber == detectionNumber) && (dets[i].frameIndex == frameIndex)) return i;
		}
		return -1;
	}
	void LoadWcsHeaders()
	{
#ifndef USE_WCSLIB
#endif
	}
	bool CheckPixelScale()
	{
			return true;
	}
	bool CheckCoordsConversion()
	{
		float tolerance = 0.000001f;
		int width = images[0].getWidth();
		int height = images[0].getHeight();
		float xs[4] = {0.0f, 1.0f, 100.0f, (float)width - 1};
		float ys[4] = {0.0f, 1.0f, 100.0f, (float)height - 1};
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < sizeof(xs) / sizeof(float); j++)
			{
				float x = xs[j];
				float y = ys[j];
				RaDec rd = images[i].PixelToWorld(x, y);
				float rtx, rty;
				images[i].WorldToPixel(rd, rtx, rty);
				float dx = abs(x - rtx);
				float dy = abs(y - rty);
				if ((dx > tolerance) || (dy > tolerance))
				{
					logd.error("Pixel coordinates (%f, %f) round-tripped to (%f, %f) which is off by (%f, %f) which is over tolerance %f",
						x, y, rtx, rty, dx, dy, tolerance);
					return false;
				}
			}
		}
		return true;
	}
	pair<float, float> ComputeImageShift(int frameIndex, int refFrameIndex, int x = 0, int y = 0)
	{
		RaDec rd1Origin = images[frameIndex].PixelToWorld(x, y);
		float x2, y2;
		images[refFrameIndex].WorldToPixel(rd1Origin, x2, y2);
		float dx = x2 - x;
		float dy = y2 - y;
		return pair<float, float>(dx, dy);
	}
	pair<float, float> ComputeLocationInOtherImage(int frameIndex, int refFrameIndex, int x = 0, int y = 0)
	{
		RaDec rd1Origin = images[frameIndex].PixelToWorld(x, y);
		float x2, y2;
		images[refFrameIndex].WorldToPixel(rd1Origin, x2, y2);
		return pair<float, float>(x2, y2);
	}
	vector<pair<float, float>> ComputeRegisteredImageCorners(int frameIndex, int refFrameIndex)
	{
		vector<pair<float, float>> corners;
		int w = images[frameIndex].getWidth();
		int h = images[frameIndex].getHeight();
		corners.push_back(ComputeLocationInOtherImage(frameIndex, refFrameIndex, 0, 0));
		corners.push_back(ComputeLocationInOtherImage(frameIndex, refFrameIndex, w - 1, 0));
		corners.push_back(ComputeLocationInOtherImage(frameIndex, refFrameIndex, w - 1, h - 1));
		corners.push_back(ComputeLocationInOtherImage(frameIndex, refFrameIndex, 0, h - 1));
		return corners;
	}
	void ComputeImageAffine(int frameIndex, int refFrameIndex, pair<float, float>& translate, double& scale, double& rotate)
	{
		int xc = images[frameIndex].getWidth() / 2;
		int yc = images[frameIndex].getHeight() / 2;
		int xtr = images[frameIndex].getWidth() - 1;
		int ytr = images[frameIndex].getHeight() - 1;
		pair<float, float> dCenter = ComputeImageShift(frameIndex, refFrameIndex, xc, yc);
		pair<float, float> dTopRight = ComputeImageShift(frameIndex, refFrameIndex, xtr, ytr);
		translate = dCenter;
		double dx = (dTopRight.first + xtr) - (dCenter.first + xc);
		double dy = (dTopRight.second + ytr) - (dCenter.second + yc);
		double dist = pyth(dx, dy);
		double normDist = pyth(xtr - xc, ytr - yc);
		scale = normDist / dist;
		rotate = atan2(dy, dx) - QUARTER_PI_FLOAT;
	}
	void CheckCentRaOffsets(int index1, int index2, DebugSpec dspec)
	{
		pair<float, float> pt = ComputeImageShift(index1, index2);
		if (registerShifts[index1].first != 0.0f)
		{
			float distx = registerShifts[index1].first - pt.first;
			float disty = registerShifts[index1].second - pt.second;
			float dist = pyth(distx, disty);
			if (dspec.debug) logd.printf("RaDec offset from %d to %d: dx, dy = %.2f, %.2f. Dist to register shift: %.2f\n", index1, index2, pt.first, pt.second, dist);
		}
		RaDec i1Origin = images[index1].PixelToWorld(0, 0);
		RaDec i1Top = images[index1].PixelToWorld(0, 4096);
		RaDec i1Right = images[index1].PixelToWorld(4096, 0);
		RaDec i2Origin = images[index2].PixelToWorld(0, 0);
		RaDec i2Top = images[index2].PixelToWorld(0, 4096);
		RaDec i2Right = images[index2].PixelToWorld(4096, 0);
		float dVert1 = i1Top.RaDegrees - i1Origin.RaDegrees;
		float dHorz1 = i1Right.DecDegrees - i1Origin.DecDegrees;
		float dVert2 = i2Top.RaDegrees - i2Origin.RaDegrees;
		float dHorz2 = i2Right.DecDegrees - i2Origin.DecDegrees;
		double radpp, decdpp;
		images[index1].GetDegreesPerPixel(radpp, decdpp);
		float dhp = abs(dVert1 - dVert2) / (float)radpp;
		float dvp = abs(dHorz1 - dHorz2) / (float)decdpp;
		if (dspec.debug) logd.printf("Horizontal and vertical scales are off from fr%d to fr%d by: %f, %f\n", index1, index2, dhp, dvp);
	}
	void CheckCentRaOffsets(DebugSpec dspec)
	{
		CheckCentRaOffsets(0, 1, dspec);
		CheckCentRaOffsets(2, 1, dspec);
		CheckCentRaOffsets(3, 1, dspec);
	}
	vector<pair<float, float>> ApplyAdditionalScaleToCorners(vector<pair<float, float>> corners)
	{
		vector<pair<float, float>> sc(4);
		int cx = getWidth() / 2;
		int cy = getHeight() / 2;
		float scale = IMAGE_AFFINE_ALIGN_ADD_SCALE;
		for (int i = 0; i < 4; i++)
		{
			float dx = (corners[i].first - cx) * scale;
			float dy = (corners[i].second - cy) * scale;
			sc[i].first = corners[i].first + dx;
			sc[i].second = corners[i].second + dy;
		}
		return sc;
	}
	void AlignOrigImagesAffineIpp(int frameIndex, int refFrameIndex, Img& tmpImage, DebugSpec dspec)
	{
	}
	void BuildImageXforms(DebugSpec dspec)
	{
		int w = getWidth();
		int h = getHeight();
		fwdXforms.resize(4);
		revXforms.resize(4);
		for (int frameIndex = 0; frameIndex < 4; frameIndex++)
		{
			if (frameIndex == 1)
			{
				Transform2d identity;
				fwdXforms[frameIndex] = identity;
				revXforms[frameIndex] = identity;
			}
			else
			{
				double rotate, scale;
				pair<float, float> translate;
				ComputeImageAffine(frameIndex, 1, translate, scale, rotate);
				fwdXforms[frameIndex] = ImgUtil::CreateForwardTransform(w, h, translate, scale, rotate);
				revXforms[frameIndex] = ImgUtil::CreateReverseTransform(w, h, translate, scale, rotate);
			}
			if (dspec.CheckDebugImageSet(this->index)) logd.debug("Frame %d revXform: %s", frameIndex, revXforms[frameIndex].ToString().c_str());
		}
	}
	pair<float, float> TransformPointPreToPostRegister(pair<float, float> pt, int frameIndex)
	{
		return fwdXforms[frameIndex].TransformPt(pt.first, pt.second);
	}
	pair<float, float> TransformPointPostToPreRegister(pair<float, float> pt, int frameIndex)
	{
		return revXforms[frameIndex].TransformPt(pt.first, pt.second);
	}
	void AlignOrigImagesAffineNonIpp(int frameIndex, int refFrameIndex, Img& tmpImage, DebugSpec dspec)
	{
		ImgUtil::Affine(images[frameIndex].origImage, revXforms[frameIndex], tmpImage);
		ImgUtil::Copy(tmpImage, images[frameIndex].origImage);
		int w = getWidth();
		int h = getHeight();
		registerShifts[frameIndex] = ComputeImageShift(frameIndex, refFrameIndex, w / 2, h / 2);
		if (dspec.CheckDebugImageSet(this->index))
		{
			logd.debug("Frame %d center shift is (%.1f, %.1f)", frameIndex, registerShifts[frameIndex].first, registerShifts[frameIndex].second);
		}
	}
	void AlignOrigImagesAffine(DebugSpec dspec)
	{
		int w = getWidth();
		int h = getHeight();
		Img i1(w, h, 16);
		BuildImageXforms(dspec);
		AlignOrigImagesAffineNonIpp(0, 1, i1, dspec);
		AlignOrigImagesAffineNonIpp(2, 1, i1, dspec);
		AlignOrigImagesAffineNonIpp(3, 1, i1, dspec);
	}
};
#include <string>
class ImageSetResult
{
public:
	string imageSetId;
	int imageSetIndex;
	int nOverlaps;
	int nTP;
	int nFP;
	int nFN;
	int nTPNeo;
	int nFPNeo;
	int nFNNeo;
	int nMatchedNeo;
	float AP;
	ImageSetResult()
	{
		imageSetId = "";
		imageSetIndex = -1;
		AP = -1;
		nOverlaps = nTP = nFP = nFN = nTPNeo = nFPNeo = nFNNeo = nMatchedNeo = 0;
	}
	std::string ToShortString()
	{
		char tmp[512];
		sprintf(tmp, "nDets: %d, nTP: %d, nFP: %d, nFN: %d, nNeo: %d, nTPNeo: %d, nFPNeo: %d, nFNNeo: %d, AP: %.0f",
			nTP + nFN, nTP, nFP, nFN, nTPNeo + nFNNeo, nTPNeo, nFPNeo, nFNNeo, AP);
		return (string)tmp;
	}
	std::string ToString()
	{
		char tmp[512];
		sprintf(tmp, "id: %s, idx: %d, nDets: %d, nOverlaps: %d, nTP: %d, nFP: %d, nFN: %d, nNeo: %d, nTPNeo: %d, nFPNeo: %d, nFNNeo: %d, AP: %.0f", 
			imageSetId.c_str(), imageSetIndex, nTP + nFN, nOverlaps, nTP, nFP, nFN, nTPNeo + nFNNeo, nTPNeo, nFPNeo, nFNNeo, AP);
		return (string)tmp;
	}
	static std::string ToCsvHeaderString()
	{
		const char* fmt = "imageSetId, imageSetIndex, nDets, nOverlaps, nTP, nFP, nFN, nNeo, nTPNeo, nFPNeo, nFNNeo, nMNeo, AP";
		char tmp[4096];
		sprintf(tmp, fmt);
		return (string)tmp;
	}
	void AppendToCsv(FILE* fp)
	{
		fprintf(fp, "%s, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %.0f\n", 
			imageSetId.c_str(), imageSetIndex, nTP + nFN, nOverlaps, nTP, nFP, nFN, nTPNeo + nFNNeo, nTPNeo, nFPNeo, nFNNeo, nMatchedNeo, AP);
	}
};
struct SlowMover
{
	vector<int> blobIndices;
	vector<pair<float, float>> centroids;
	float totalDist;
	int corrCAIndex; 
	SlowMover()
	{
		blobIndices.resize(4);
		centroids.resize(4);
		totalDist = -1;
	}
	string ToString()
	{
		char tmp[512];
		sprintf(tmp, "blobIndices: (%d, %d, %d, %d), start: (%d, %d)", blobIndices[0], blobIndices[1], blobIndices[2], blobIndices[3],
			roundf(centroids[0].first), roundf(centroids[0].second));
		return (string)tmp;
	}
};
void FindSlowMovers(vector<NsiBlobSet>& bsets, vector<SlowMover>& sms, DebugSpec dspec)
{
	int n0 = bsets[0].size();
	for (int i = 0; i < n0; i++)
	{
		NsiBlob& f0Blob = bsets[0].nsiBlobs[i];
		bool debugThis = f0Blob.CheckDebug(dspec);
		float x = bsets[0].nsiBlobs[i].imgBlob.xCentroid;
		float y = bsets[0].nsiBlobs[i].imgBlob.yCentroid;
		SlowMover sm;
		sm.blobIndices[0] = i;
		sm.centroids[0] = pair<float, float>(f0Blob.imgBlob.xCentroid, f0Blob.imgBlob.yCentroid);
		float blobRadius = sqrtf(f0Blob.imgBlob.area) / 2;
		float maxDist = blobRadius; 
		if (debugThis) logd.debug("FindSlowMovers: Searching on blob %d at (%.1f, %.1f) with radius %.1f",
			i, f0Blob.imgBlob.xCentroid, f0Blob.imgBlob.yCentroid, blobRadius);
		for (int frameIndex = 1; frameIndex < 4; frameIndex++)
		{
			vector<int> nis = bsets[frameIndex].pointIndex.FindNearby(x, y, maxDist, true, false);
			if (nis.size() > 0)
			{
				bool foundOne = false;
				for (int j = 0; j < nis.size(); j++)
				{
					NsiBlob& otherBlob = bsets[frameIndex].nsiBlobs[nis[j]];
					float ss = otherBlob.CompSim(f0Blob, false);
					if (ss > BLOB_SIM_SCORE_MIN_THRESHOLD)
					{
						foundOne = true;
						sm.blobIndices[frameIndex] = nis[j];
						sm.centroids[frameIndex] = pair<float, float>(otherBlob.imgBlob.xCentroid, otherBlob.imgBlob.yCentroid);
						if (debugThis) logd.debug("FindSlowMovers: fr%d: For blob %d: found neighbor in frame %d %d at (%.1f, %.1f) with simScore %.3f",
							frameIndex, i, frameIndex, nis[j], otherBlob.imgBlob.xCentroid, otherBlob.imgBlob.yCentroid, ss);
						if (frameIndex == 3)
						{
							sm.totalDist = pyth(f0Blob.imgBlob.xCentroid - otherBlob.imgBlob.xCentroid, f0Blob.imgBlob.xCentroid - otherBlob.imgBlob.xCentroid);
							if (sm.totalDist > SLOW_MOVER_MIN_TOTAL_MOVE_DIST_PIX)
							{
								sms.push_back(sm);
								if (debugThis) logd.debug("FindSlowMovers: fr%d: Added slow mover %d based on blob %d with totalDist: %.3f",
									frameIndex, sms.size() - 1, i, sm.totalDist);
							}
						}
						x = otherBlob.imgBlob.xCentroid;
						y = otherBlob.imgBlob.yCentroid;
					}
					else
					{
						if (debugThis) logd.debug("FindSlowMovers: fr%d: For blob %d: found neighbor in frame %d %d at (%.1f, %.1f) but bad simScore %.3f",
							frameIndex, i, frameIndex, nis[0], otherBlob.imgBlob.xCentroid, otherBlob.imgBlob.yCentroid, ss);
					}
				}
				if (!foundOne)
				{
					if (debugThis) logd.debug("FindSlowMovers: fr%d: For blob %d: at (%.1f, %.1f) with radius %.1f, found multiple neighbors in range, but none similar enough.",
						frameIndex, i, f0Blob.imgBlob.xCentroid, f0Blob.imgBlob.yCentroid, blobRadius);
					break;
				}
			}
			else
			{
				if (debugThis) logd.debug("FindSlowMovers: fr%d: For blob %d at (%.1f, %.1f) with radius %.1f, from ref pt (%.1f, %.1f), no neighbors in range.",
					frameIndex, i, f0Blob.imgBlob.xCentroid, f0Blob.imgBlob.yCentroid, blobRadius, x, y);
				break;
			}
		}
	}
	if (dspec.debug) logd.debug("FindSlowMovers: Found %d candidate slow movers.", sms.size());
}
void AddSlowMoversToCAList(string imageSetId, vector<SlowMover>& sms, vector<NsiBlobSet>& bsets, CAList& cal, DebugSpec dspec)
{
	const int BlobIndexOffset = BASE_BLOB_INDEX_FOR_SLOW_MOVERS;
	for (int i = 0; i < sms.size(); i++)
	{
		int fr0BlobIndex = sms[i].blobIndices[0];
		int fr1BlobIndex = sms[i].blobIndices[1];
		bool debugThis = bsets[0].nsiBlobs[fr0BlobIndex].CheckDebug(dspec);
		CA ca(0, bsets[0].nsiBlobs[fr0BlobIndex], 1, bsets[1].nsiBlobs[fr1BlobIndex]);
		ca.imageSetId = imageSetId;
		ca.index = cal.size();
		sms[i].corrCAIndex = ca.index;
		ca.InitBlobIndices(); 
		ca.blobIndices[0] = fr0BlobIndex + BlobIndexOffset;
		ca.blobIndices[1] = fr1BlobIndex + BlobIndexOffset;
		for (int j = 2; j < 4; j++)
		{
			int blobIndex = sms[i].blobIndices[j];
			ca.SetBlob(j, bsets[j].nsiBlobs[blobIndex]); 
			ca.blobIndices[j] = blobIndex + BlobIndexOffset;
		}
		if (debugThis || (dspec.caIndex == ca.index)) logd.debug("AddSlowMoversToCAList: Adding CA: %s", ca.ToString().c_str());
		if (!cal.Add(ca, dspec))
		{
			logd.error("AddSlowMoversToCAList: Error trying to add: %s", ca.ToString().c_str());
			ErrorExit("AddSlowMoversToCAList: Somehow had a blob indices collision on adding this.");
		}
	}
}
void SlowMoversFilterNsiAndCreateBlobSet(Img& nsi, Img& tmp1, Img& tmp2, NsiBlobSet& bset, DebugSpec dspec)
{
	Img8Util::DropSmallBlobs(nsi, tmp1, SLOW_MOVERS_MIN_NEIGHBOR_PIX, SLOW_MOVERS_MIN_NEIGHBOR_PIX_MEAN, SLOW_MOVERS_MIN_NEIGHBOR_PIX_MAX);
	bset.Create(tmp1, tmp2, dspec);
}
void SaveSlowMoverImages(ImageSet& iset, vector<SlowMover>& sms)
{
	int n = sms.size();
	int dim = 128;
	char tag[512];
	for (int frameIndex = 0; frameIndex < 4; frameIndex++)
	{
		vector<tuple<int, int, uint32_t, string>> nsiRois;
		for (int i = 0; i < n; i++)
		{
			int x0, y0;
			x0 = roundf(sms[i].centroids[frameIndex].first);
			y0 = roundf(sms[i].centroids[frameIndex].second);
			sprintf(tag, "%03d", sms[i].corrCAIndex);
			nsiRois.push_back(tuple<int, int, uint32_t, string>(roundf(x0), roundf(y0), 0xFFFFFFFF, (string)tag));
		}
		Img text(iset.getWidth(), iset.getHeight(), 8);
		Img8Util::Clear(text);
		sprintf(tag, "slowMovers_fr%d", frameIndex);
		Img8Util::Copy(iset.images[frameIndex].nsi, text);
		RenderColoredRois(text, tag, 8, nsiRois);
	}
}
void FindSlowMoversTop(ImageSet& iset, CAList& cal, DebugSpec dspec)
{
	bool debug = dspec.CheckDebugImageSet(iset.index);
	int w = iset.getWidth();
	int h = iset.getHeight();
	Img i1(w, h, 8), i2(w, h, 8);
	vector<NsiBlobSet> bsets(4);
	if (debug) SaveDebugTiff("f0", iset.images[0].nsi);
	SlowMoversFilterNsiAndCreateBlobSet(iset.images[0].nsi, i1, i2, bsets[0], dspec);
	if (debug) SaveDebugTiff("filtered-0", i2);
	if (debug) SaveDebugTiff("f1", iset.images[1].nsi);
	SlowMoversFilterNsiAndCreateBlobSet(iset.images[1].nsi, i1, i2, bsets[1], dspec);
	if (debug) SaveDebugTiff("filtered-1", i2);
	if (debug) SaveDebugTiff("f2", iset.images[2].nsi);
	SlowMoversFilterNsiAndCreateBlobSet(iset.images[2].nsi, i1, i2, bsets[2], dspec);
	if (debug) SaveDebugTiff("filtered-2", i2);
	if (debug) SaveDebugTiff("f3", iset.images[3].nsi);
	SlowMoversFilterNsiAndCreateBlobSet(iset.images[3].nsi, i1, i2, bsets[3], dspec);
	if (debug) SaveDebugTiff("filtered-3", i2);
	i1.destroy();
	i2.destroy();
	vector<SlowMover> sms;
	FindSlowMovers(bsets, sms, dspec);
	int nToFix = MIN(sms.size(), 2000);
	for (int i = 0; i < nToFix; i++)
	{
		bool debugThis = bsets[0].nsiBlobs[sms[i].blobIndices[0]].CheckDebug(dspec);
		if (debugThis) logd.printf("    After FindSlowMovers: %d: %s\n", i, sms[i].ToString().c_str());
		for (int frameIndex = 0; frameIndex < 4; frameIndex++)
		{
			int bi = sms[i].blobIndices[frameIndex];
			NsiBlob& blob = bsets[frameIndex].nsiBlobs[bi];
			float oldSalienceScore = blob.salienceScore;
			float oldCaBsmScore = blob.bsmScore;
			float newSalienceScore, newCaBsmScore;
			iset.blobSets[frameIndex].ComputeSingleSalienceScore(blob, newSalienceScore, newCaBsmScore, dspec);
			if (debugThis) logd.debug("Updating salienceScore of blob (%.0f, %.0f) from %.3f to %.3f and bsm from %.3f to %.3f.",
				blob.imgBlob.xCentroid, blob.imgBlob.yCentroid, oldSalienceScore, newSalienceScore, oldCaBsmScore, newCaBsmScore);
			blob.salienceScore = newSalienceScore;
			blob.bsmScore = newCaBsmScore;
		}
	}
	AddSlowMoversToCAList(iset.imageSetId, sms, bsets, cal, dspec);
	if (dspec.info) logd.info("FindSlowMoversTop: Added %d slow mover CA's to the cal.", sms.size());
	if (dspec.CheckDebugImageSet(iset.index))
	{
		SaveSlowMoverImages(iset, sms);
	}
}
void FixupAndMask(Img& image, Img& mask, bool debug = false)
{
	int w = image.getWidth();
	int h = image.getHeight();
	if (!DO_CROP_INPUT_IMAGE)
	{
		if (debug) logd.printf("ZeroMargin...\n");
		int margin = 8;
		for (int y = 0; y < margin; y++) ImgUtil::ZeroCol(image, y);
		for (int i = 0; i < margin; i++) ImgUtil::ZeroCol(image, w - 1 - i);
		if (debug) SaveDebugTiff("zeroMargin", image);
	}
	if (debug) logd.printf("ZeroColsByProfile...\n");
	ImgUtil::ZeroColsByProfile(image, 3.0f, 1000.0f, debug); 
	if (debug) SaveDebugTiff("zeroCols", image);
	if (debug) logd.printf("ZeroRowsByProfile...\n");
	ImgUtil::ZeroRowsByProfile(image, 3.0f, 10.0f, debug); 
	if (debug) SaveDebugTiff("zeroRows", image);
	if (debug) logd.printf("ThresholdRangeSse2...\n");
	ImgUtil::ThresholdRangeSse2(image, 1000, 60000, false, mask); 
	if (debug) SaveDebugTiff("extrema mask", mask);
	if (debug) logd.printf("ComputeStats...\n");
	Stats stats1 = ImgUtil::ComputeStats(image, mask);
	if (debug) logd.printf("Stats1: %s\n", stats1.ToString().c_str());
	if (debug) logd.printf("MaskWrappedValues...\n");
	ImgUtil::MaskWrappedValues(image, stats1, 10, mask); 
	if (debug) SaveDebugTiff("wrapMask", mask);
	if (debug) logd.printf("Dilate...\n");
	tb1.resize(w, h, 8);
	Img8Util::Dilate3x3(mask, tb1);
	Img8Util::Dilate3x3(tb1, mask);
	if (debug) SaveDebugTiff("dilatedMask", mask);
	if (debug) logd.printf("Apply...\n");
	ImgUtil::ApplyMask(image, mask);
	if (debug) SaveDebugTiff("maskedOrig", image);
}
void CreateContourImage(Img& mf, Img& contour, int segmentSize, bool debug)
{
	int w = mf.getWidth();
	int h = mf.getHeight();
	int smallest = 8; 
	ts1.resize(w, h, 16);
	Img i1(w / 2, h / 2, 16); 
	Img i2(w / 4, h / 4, 16); 
	Img tmp1(w / smallest, h / smallest, 16);
	Img tmp2(w / smallest, h / smallest, 16);
	bool doFilter = false; 
	ImgUtil::DownsampleSse2(mf, i1, doFilter); 
	ImgUtil::DownsampleSse2(i1, i2, doFilter);
	ImgUtil::DownsampleSse2(i2, tmp1, doFilter);
	if (debug) SaveDebugTiff("contour-ds", tmp1);
	ImgUtil::FillZeros(tmp1, tmp2);
	if (debug) SaveDebugTiff("contour-fillzero", tmp2);
	ImgUtil::Median3x3GoodValues(tmp2, tmp1);
	ImgUtil::Median3x3GoodValues(tmp1, tmp2);
	ImgUtil::Median3x3GoodValues(tmp2, tmp1);
	if (debug) SaveDebugTiff("contour-median", tmp1);
	ImgUtil::Mean(tmp1, 17, tmp2);
	if (debug) SaveDebugTiff("contour-mean", tmp2);
	ImgUtil::Mean3x3GoodValues(tmp2, tmp1, true);
	ImgUtil::Mean3x3GoodValues(tmp1, tmp2, true);
	if (debug) SaveDebugTiff("contour-mean", tmp2);
	ImgUtil::WrapMinAllDir(tmp2, segmentSize, tmp1);
	if (debug) SaveDebugTiff("contour-wrapMin", tmp1);
	ImgUtil::Mean(tmp1, 7, tmp2);
	if (debug) SaveDebugTiff("fillAndMean", tmp2);
	ImgUtil::Upsample(tmp2, i2);
	ImgUtil::Upsample(i2, i1);
	ImgUtil::Upsample(i1, ts1);
	ImgUtil::Add(ts1, -1000, contour);
}
void CreateContourImage2(Img& fixedImage, Img& bgMask, Img& contour, bool debug)
{
	int w = fixedImage.getWidth();
	int h = fixedImage.getHeight();
	ImgUtil::MeanMasked(fixedImage, bgMask, 17, contour);
	tb1.resize(w, h, 8);
	Img8Util::Not(bgMask, tb1);
	ImgUtil::FillUnderMask(contour, tb1);
}
void CreateFgBgMasks(Img& image, Img& mf, Img& mask, Img& fgMask, Img& bgMask, Img& dipMask, bool debug)
{
	int w = image.getWidth();
	int h = image.getHeight();
	tb1.resize(w, h, 8);
	Hist hist1;
	ImgUtil::HistCompute(mf, mask, 4096, hist1); 
	if (debug) logd.printf("Hist: %s\n", hist1.ToString().c_str());
	uint16_t dipthresh = (uint16_t)hist1.GetPercentileBinLevel(0.05f);
	uint16_t bgthresh = (uint16_t)hist1.GetPercentileBinLevel(0.85f);
	if (debug) logd.printf("dipthresh: %d, bgthresh: %d\n", dipthresh, bgthresh);
	dipMask.resize(w, h, 8);
	ImgUtil::Threshold2(mf, dipthresh, true, dipMask);
	Img8Util::Erode3x3(dipMask, tb1);
	Img8Util::Dilate3x3(tb1, dipMask);
	Img8Util::AndNot(dipMask, mask, tb1);
	Img8Util::Copy(tb1, dipMask);
	if (debug) SaveDebugTiff("prelimDipMask", dipMask);
	bgMask.resize(w, h, 8);
	ImgUtil::Threshold2(mf, bgthresh, true, bgMask);
	Img8Util::AndNot(bgMask, dipMask, tb1); 
	Img8Util::AndNot(tb1, mask, bgMask);
	Img8Util::Erode3x3(bgMask, tb1);
	Img8Util::Dilate3x3(tb1, bgMask);
	if (debug) SaveDebugTiff("prelimBgMask", bgMask);
	Stats bgStats = ImgUtil::ComputeStats(image, bgMask, true);
	uint16_t fgthresh = bgStats.mean + 2 * bgStats.stddev;
	dipthresh = (uint16_t)roundf(bgStats.mean - 2 * bgStats.stddev);
	if (debug) logd.printf("bgStats: %s, dipthresh: %d, fgthresh: %d\n", bgStats.ToString().c_str(), dipthresh, fgthresh);
	ImgUtil::Threshold2(mf, dipthresh, true, dipMask);
	Img8Util::AndNot(dipMask, mask, tb1);
	Img8Util::Copy(tb1, dipMask);
	if (debug) SaveDebugTiff("dipMask", dipMask);
	fgMask.resize(w, h, 8);
	ImgUtil::Threshold2(mf, fgthresh, false, fgMask);
	if (debug) SaveDebugTiff("fgThreshold", fgMask);
	Img8Util::Copy(fgMask, tb1);
	Img8Util::AndNot(tb1, mask, fgMask);
	tb2.resize(w, h, 8);
	Img8Util::Dilate3x3(fgMask, tb2);
	Img8Util::Not(tb2, tb1);
	Img8Util::AndNot(tb1, mask, tb2);
	Img8Util::AndNot(tb2, dipMask, bgMask);
}
void CreateSboMask(Img& fgMask, Img& badMask, Img& sboMask, bool debug)
{
	int w = fgMask.getWidth();
	int h = fgMask.getHeight();
	tb1.resize(w, h, 8);
	tb2.resize(w, h, 8);
	ts1.resize(w, h, 16); 
	Img8Util::Erode3x3(fgMask, tb1);
	ImgUtil::Label(tb1, ts1, true);
	vector<ImgBlob> blobs;
	ImgUtil::FindBlobs(ts1, blobs);
	int n = blobs.size();
	vector<ImgBlob> keepers;
	for (int i = 0; i < n; i++)
	{
		if (blobs[i].area > SBO_BLOB_MIN_AREA)
		{
			keepers.push_back(blobs[i]);
		}
	}
	Img8Util::Clear(tb1);
	ImgUtil::CopyBlobs(fgMask, ts1, keepers, tb1); 
	Img8Util::Dilate3x3(tb1, tb2); 
	Img8Util::Copy(tb2, sboMask);
	Img8Util::FillUnderMask(sboMask, badMask);
}
void CreateNsiImageSub(Img& image, Img& bgMask, Img& nsi, bool debug = false)
{
	if (debug) SaveDebugTiff("imageSub", image);
	if (debug) SaveDebugTiff("bgMaskSub", bgMask);
	Stats bgStats = ImgUtil::ComputeStats(image, bgMask, true);
	if (debug) logd.printf("CreateNsiImageSub: ComputeStats on bg: %s\n", bgStats.ToString().c_str());
	float mean = bgStats.mean;
	float stddev = bgStats.meandev;
	if (debug) logd.printf("CreateNsiImage: mean: %f, stddev: %f\n", mean, stddev);
	ImgUtil::NStddevOver8bit(image, mean, stddev, NSI_NSTDDEV_MIN, NSI_NSTDDEV_MAX, nsi); 
}
void CreateNsiImageTop(FitsImage2* p, bool debug)
{
	int w = p->getWidth();
	int h = p->getHeight();
	if (debug) logd.printf("CreateNsiImageTop: FixupAndMask...\n");
	p->fixedImage.resize(w, h, 16);
	p->mask.resize(w, h, 8);
	ImgUtil::Copy(p->origImage, p->fixedImage);
	if (debug) SaveDebugTiff("orig", p->origImage);
	FixupAndMask(p->fixedImage, p->mask, false);
	if (debug) SaveDebugTiff("fixed", p->fixedImage);
	if (debug) SaveDebugTiff("mask", p->mask);
	p->origImage.destroy();
	Img mf3(w, h, 16);
	ImgUtil::Mean3x3GoodValues(p->fixedImage, mf3); 
	if (debug) SaveDebugTiff("mean3x3masked", mf3);
	Img contour(w, h, 16); 
	if (debug) logd.printf("CreateNsiImageTop: CreateContourImage...\n");
	CreateContourImage(mf3, contour, CONTOUR_DS_SEG_SIZE, false);
	if (debug) SaveDebugTiff("mean3x3", mf3);
	if (debug) SaveDebugTiff("contour", contour);
	ts3.resize(w, h, 16);
	Img decont(w, h, 16);
	ImgUtil::Subtract(p->fixedImage, contour, ts3);
	if (debug) SaveDebugTiff("de-contour", ts3);
	ImgUtil::Copy(ts3, p->fixedImage); 
	if (debug) logd.printf("CreateNsiImageTop: Mean3x3GoodValues...\n");
	ImgUtil::Mean3x3GoodValues(ts3, mf3); 
	if (debug) SaveDebugTiff("mean3x3masked", mf3);
	if (debug) logd.printf("BG mask tiles...\n");
	CreateFgBgMasks(p->fixedImage, mf3, p->mask, p->fgMask, p->bgMask, p->dipMask, debug);
	p->sboMask.resize(w, h, 8);
	CreateSboMask(p->fgMask, p->mask, p->sboMask, debug);
	if (debug) SaveDebugTiff("sboMask", p->sboMask);
	ts1.destroy(); 
	if (debug)
	{
		SaveDebugTiff("fgMask", p->fgMask);
		SaveDebugTiff("bgMask", p->bgMask);
		SaveDebugTiff("dipMask", p->dipMask);
		ts3.copyImg(p->fixedImage);
		ImgUtil::ApplyMask(ts3, p->bgMask, false);
		SaveDebugTiff("bg", ts3);
		ImgUtil::Copy(p->fixedImage, ts3);
		ImgUtil::ApplyMask(ts3, p->dipMask, false);
		SaveDebugTiff("dip", ts3);
		ImgUtil::Copy(p->fixedImage, ts3);
		ImgUtil::ApplyMask(ts3, p->fgMask, false);
		SaveDebugTiff("fg", ts3);
	}
	ts3.destroy();
	p->fgMask.destroy();
	p->nsi.resize(w, h, 8);
	CreateNsiImageSub(p->fixedImage, p->bgMask, p->nsi, debug);
	p->bgMask.destroy();
	p->fixedImage.destroy();
}
void CreateDSPyramid(FitsImage2* p, int nDown, bool debug)
{
	if (p->dsPyramid.size() > 0) return; 
	int w = p->getWidth();
	int h = p->getHeight();
	ImgUtil::FillUnderMask(p->fixedImage, p->mask);
	if (debug) SaveDebugTiff("fillUnderMask", p->fixedImage);
	for (int i = 1; i <= nDown; i++)
	{
		Img* down = new Img();
		down->resize(w >> i, h >> i, 8);
		if (i == 1)
		{
			Img8Util::Downsample(p->nsi, *down, true);
		}
		else
		{
			Img8Util::Downsample(*(p->dsPyramid[i - 2]), *down, true);
		}
		p->dsPyramid.push_back(down);
		if (debug)
		{
			char tmp[64];
			sprintf(tmp, "ds_%d", i);
			SaveDebugTiff(tmp, *down);
		}
	}
}
void ShiftFitsImageInt(FitsImage2* pi, pair<int, int> offset)
{
	int maskBgValue = 0;
	if (pi->fixedImage.getWidth() > 0) ImgUtil::ShiftInt(pi->fixedImage, offset, 0);
	if (pi->mask.getWidth() > 0) ImgUtil::ShiftInt(pi->mask, offset, 255); 
	if (pi->fgMask.getWidth() > 0) ImgUtil::ShiftInt(pi->fgMask, offset, maskBgValue);
	if (pi->bgMask.getWidth() > 0) ImgUtil::ShiftInt(pi->bgMask, offset, maskBgValue);
	if (pi->dipMask.getWidth() > 0) ImgUtil::ShiftInt(pi->dipMask, offset, maskBgValue);
	if (pi->nsi.getWidth() > 0) ImgUtil::ShiftInt(pi->nsi, offset, 0);
	if (pi->sboMask.getWidth() > 0) ImgUtil::ShiftInt(pi->sboMask, offset, 0);
	pi->DeleteDsPyramid();
}
void ShiftFitsImageFloat(FitsImage2* pi, pair<float, float> offset)
{
	int maskBgValue = 0;
	pair<int, int> offsetInt;
	offsetInt.first = (int)roundf(offset.first);
	offsetInt.second = (int)roundf(offset.second);
	if ((offsetInt.first != 0) || (offsetInt.second != 0))
	{
		if (pi->fixedImage.getWidth() > 0) ImgUtil::ShiftInt(pi->fixedImage, offsetInt, 0);
		if (pi->mask.getWidth() > 0) Img8Util::ShiftInt(pi->mask, offsetInt, 255); 
		if (pi->fgMask.getWidth() > 0) Img8Util::ShiftInt(pi->fgMask, offsetInt, maskBgValue);
		if (pi->bgMask.getWidth() > 0) Img8Util::ShiftInt(pi->bgMask, offsetInt, maskBgValue);
		if (pi->sboMask.getWidth() > 0) Img8Util::ShiftInt(pi->sboMask, offsetInt, 0);
	}
	Img8Util::ShiftSmart(pi->nsi, offset.first, offset.second, tb1);
	Img8Util::Copy(tb1, pi->nsi);
	pi->DeleteDsPyramid();
}
void ShiftFitsImageAllAtOnce(ImageSet& iset, int frameIndex, pair<float, float> offset)
{
	if (DO_SUBPIXEL_SHIFT)
	{
		pair<int, int> offsetInt;
		offsetInt.first = (int)floor(offset.first);
		offsetInt.second = (int)floor(offset.second);
		pair<float, float> offsetFloat;
		offsetFloat.first = offset.first - offsetInt.first;
		offsetFloat.second = offset.second - offsetInt.second;
		ShiftFitsImageInt(&iset.images[frameIndex], offsetInt);
		ShiftFitsImageFloat(&iset.images[frameIndex], offsetFloat);
	}
	else
	{
		pair<int, int> offsetInt;
		offsetInt.first = (int)roundf(offset.first);
		offsetInt.second = (int)roundf(offset.second);
		ShiftFitsImageInt(&iset.images[frameIndex], offsetInt);
		offset.first = offsetInt.first;
		offset.second = offsetInt.second;
	}
	iset.registerShifts[frameIndex] = offset;
}
void SaveFitsImage(FitsImage2* p)
{
	SaveDebugTiff("fixedImage", p->fixedImage);
	SaveDebugTiff("mask", p->mask);
	SaveDebugTiff("dipMask", p->dipMask);
	SaveDebugTiff("nsi", p->nsi);
}
void SaveFitsImages(ImageSet& iset, bool justNsi = false)
{
	if (!justNsi)
	{
		for (int i = 0; i < iset.images.size(); i++)
		{
			SaveDebugTiff("fixedImage", iset.images[i].fixedImage);
		}
		for (int i = 0; i < iset.images.size(); i++)
		{
			SaveDebugTiff("mask", iset.images[i].mask);
		}
		for (int i = 0; i < iset.images.size(); i++)
		{
			SaveDebugTiff("dipMask", iset.images[i].dipMask);
		}
		for (int i = 0; i < iset.images.size(); i++)
		{
			SaveDebugTiff("sboMask", iset.images[i].sboMask);
		}
		for (int i = 0; i < iset.images.size(); i++)
		{
			SaveDebugTiff("fgMask", iset.images[i].fgMask);
		}
		for (int i = 0; i < iset.images.size(); i++)
		{
			SaveDebugTiff("bgMask", iset.images[i].bgMask);
		}
	}
	for (int i = 0; i < iset.images.size(); i++)
	{
		if (iset.images[i].nsi.getWidth() > 0) SaveDebugTiff("nsi", iset.images[i].nsi);
	}
	for (int i = 0; i < iset.images.size(); i++)
	{
		if (iset.images[i].nsiDiff.getWidth() > 0) SaveDebugTiff("nsiDiff", iset.images[i].nsiDiff);
	}
}
bool RegisterAndShift(ImageSet& iset, int frameIndex1, int frameIndex2, bool debug)
{
	FitsImage2* pImg1 = &iset.images[frameIndex1];
	FitsImage2* pImg2 = &iset.images[frameIndex2];
	if (!DO_IMAGE_REGISTER)
	{
		pair<float, float> foffset = iset.ComputeImageShift(frameIndex1, frameIndex2);
		if (debug) logd.printf("RaDec shift %d to %d: %f, %f\n", frameIndex1, frameIndex2, foffset.first, foffset.second);
		foffset.first = -foffset.first;
		foffset.second = -foffset.second;
		ShiftFitsImageFloat(pImg2, foffset);
		iset.registerShifts[frameIndex2] = foffset;
	}
	else
	{
		ErrorExit("Commented out");
	}
	return true;
}
void RegisterImagesTop(ImageSet& imageSet, bool shortcutRegister, DebugSpec dspec)
{
	if (!DO_IMAGE_AFFINE_ALIGN)
	{
		ErrorExit("Commented out");
	}
}
void EnactDipMasks(ImageSet& iset, bool debug)
{
	int n = iset.images.size();
	int w = iset.images[0].dipMask.getWidth();
	int h = iset.images[0].dipMask.getHeight();
	tb1.resize(w, h, 8);
	Img combinedDipMask(iset.images[0].dipMask);
	for (int i = 1; i < n; i++)
	{
		Img8Util::Max(combinedDipMask, iset.images[i].dipMask, tb1);
		Img8Util::Copy(tb1, combinedDipMask);
		iset.images[i].dipMask.destroy();
	}
	if (debug) SaveDebugTiff("combinedDipMask", combinedDipMask);
	for (int i = 0; i < n; i++)
	{
		Img8Util::Max(combinedDipMask, iset.images[i].mask, tb1);
		Img8Util::Copy(tb1, iset.images[i].mask);
	}
}
void FilterNsiDiffImage(Img& img, bool debug)
{
	int w = img.getWidth();
	int h = img.getHeight();
	tb1.resize(w, h, 8);
	tb2.resize(w, h, 8);
	Img8Util::Threshold(img, NSI_DIFF_FLOOR, false, tb2);
	Img8Util::And(img, tb2, tb1);
	if (debug) SaveDebugTiff("FilterNsiDiffImage-thresholded", tb1);
	Img8Util::DropSmallBlobs(tb1, tb2, NSI_DIFF_MIN_NEIGHBOR_COUNT, NSI_DIFF_MIN_NEIGHBOR_MEAN, NSI_DIFF_MIN_NEIGHBOR_MAX);
	if (debug) SaveDebugTiff("FilterNsiDiffImage-dropSmallBlobs", tb2);
	Img8Util::Copy(tb2, img);
}
void CreateNsiDiffImages(ImageSet& iset, DebugSpec dspec)
{
	bool debug = dspec.CheckDebugImageSet(iset.index);
	int n = iset.images.size();
	int w = iset.images[0].getWidth();
	int h = iset.images[0].getHeight();
	Img nsic(w, h, 8);
	Img nsicMask(w, h, 8);
	Img tmp(w, h, 8);
	Img8Util::Median4Masked(
		iset.images[0].nsi,
		iset.images[1].nsi,
		iset.images[2].nsi,
		iset.images[3].nsi,
		iset.images[0].mask,
		iset.images[1].mask,
		iset.images[2].mask,
		iset.images[3].mask,
		nsic, nsicMask);
	Img8Util::Dilate3x3(nsic, tmp);
	Img8Util::Multiply(tmp, NSI_MEDIAN_COMPENSATION_SCALE, nsic);
	if (debug) SaveDebugTiff("createNsiDiff-nsiCombinedMask", nsicMask);
	if (debug) SaveDebugTiff("createNsiDiff-nsiCombined", nsic);
	iset.nsiCM.resize(w, h, 8);
	Img8Util::Add(nsicMask, nsic, iset.nsiCM);
	if (debug) SaveDebugTiff("createNsiDiff-nsiCM", iset.nsiCM);
	for (int i = 0; i < n; i++)
	{
		if (debug) SaveDebugTiff("createNsiDiff-pre-diff", iset.images[i].nsi);
		iset.images[i].nsiDiff.resize(w, h, 8);
		Img8Util::Subtract(iset.images[i].nsi, nsic, tmp);
		if (debug) SaveDebugTiff("createNsiDiff-diff", tmp);
		Img8Util::AndNot(tmp, nsicMask, iset.images[i].nsiDiff);
		if (debug) SaveDebugTiff("createNsiDiff-masked", iset.images[i].nsiDiff);
		FilterNsiDiffImage(iset.images[i].nsiDiff, debug); 
		Img8Util::Copy(iset.images[i].nsiDiff, tmp);
		Img8Util::AndNot(tmp, iset.images[i].sboMask, iset.images[i].nsiDiff);
		if (debug) SaveDebugTiff("createNsiDiff-notSbo", iset.images[i].nsiDiff);
	}
}
float GetExplainedMissingBlobScore(ImageSet& iset, CA& ca, int frameIndex, bool debug)
{
	if (ca.CheckIsBlobSet(frameIndex)) ErrorExit("Already a blob set, no need to search for corresponding blob.");
	pair<float, float> epoint = ca.GetFittedCoords(frameIndex); 
	float meanBlobArea = ca.ComputeMeanBlobArea();
	float dim = MAX(1.0f, sqrt(meanBlobArea));
	float rad = dim/2.0f;
	int x0 = roundf(epoint.first - rad);
	int y0 = roundf(epoint.second - rad);
	int idim = roundf(dim);
	IRoi roi(x0, y0, x0 + idim, y0 + idim);
	if (!roi.checkIsInImage(iset.getWidth(), iset.getHeight()))
	{
		if (debug) logd.debug("GetExplainedMissingBlobScore: Roi off image, so score 1. CA %d, epoint: (%.1f, %.1f) in fr%d, roi: %s", 
			ca.index, epoint.first, epoint.second, frameIndex, roi.ToString().c_str());
		return 1;
	}
	IRoi unshiftedRoi = roi;
	unshiftedRoi.offset(-roundf(iset.registerShifts[frameIndex].first), -roundf(iset.registerShifts[frameIndex].second));
	if (!unshiftedRoi.checkIsInImage(iset.getWidth(), iset.getHeight()))
	{
		if (debug) logd.debug("GetExplainedMissingBlobScore: Shifted roi off image, so score 1. CA %d, epoint: (%.1f, %.1f) in fr%d, roi: %s, shiftedRoi: %s",
			ca.index, epoint.first, epoint.second, frameIndex, roi.ToString().c_str(), unshiftedRoi.ToString().c_str());
		return 1;
	}
	int sumNsicm = Img8Util::SumRoi(iset.nsiCM, roi);
	int meanBlobSum = ca.GetMeanBlobSum();
	float score = (float)sumNsicm / meanBlobSum;
	if (debug) logd.debug("GetExplainedMissingBlobScore: CA %d, epoint: (%.1f, %.1f) in fr%d, roi: %s, sumNsicm: %d, meanBlobArea: %.1f, meanBlobSum: %d, score: %.3f",
		ca.index, epoint.first, epoint.second, frameIndex, roi.ToString().c_str(), sumNsicm, meanBlobArea, meanBlobSum, score);
	score = CLIP(score, 0.0f, 1.0f);
	return score;
}
int GetCorrespondingBlob(ImageSet& iset, CA& ca, int frameIndex, bool debug)
{
	if (ca.CheckIsBlobSet(frameIndex)) ErrorExit("Already a blob set, no need to search for corresponding blob.");
	float maxDist = (float)(BLOB_DIST_MAX_THRESHOLD_PIX);
	if (!ca.CheckIsSetBlobAdjacent(frameIndex))
	{
		maxDist *= 1.5f; 
	}
	pair<float, float> epoint = ca.GetFittedCoords(frameIndex); 
	if ((epoint.first < 0) || (epoint.second < 0) || (epoint.first > iset.getWidth() - 1) || (epoint.second > iset.getHeight() - 1))
	{
		if (debug) logd.printf("    GetCorrespondingBlob: Extrapolated coords for are off image, so no nearby blob.\n");
		return -1;
	}
	vector<int> ci = iset.blobSets[frameIndex].pointIndex.FindNearby(epoint, maxDist, false);
	if (debug) logd.printf("    GetCorrespondingBlob: Extrapolated coords to frame %d are (%.1f, %.1f) and there are %d nearby blobs.\n", 
		frameIndex, epoint.first, epoint.second, ci.size());
	if (ci.size() == 0)
	{
		if (debug) logd.printf("    GetCorrespondingBlob: No nearby blobs within dist %f.\n", maxDist);
		return -1;
	}
	float bestMatchScore = -1;
	int bestIndex = -1;
	for (int i = 0; i < ci.size(); i++)
	{
		NsiBlob& otherBlob = iset.blobSets[frameIndex].nsiBlobs[ci[i]];
		float dist = otherBlob.Dist(epoint);
		float ss = ca.ComputeBlobSimScore(otherBlob);
		float minDistForRangeScore = 1.5f; 
		float maxDistForRangeScore = MAX(maxDist, 4.0f); 
		float distScore = 1.0f - rangeScore(minDistForRangeScore, maxDistForRangeScore, dist);
		float matchScore = geometricMean(distScore, ss);
		if (debug)
		{
			logd.printf("    GetCorrespondingBlob: To blob index %d is dist: %.1f, distScore: %.2f, matchScore: %.3f\n", ci[i], dist, distScore, matchScore);
			logd.printf("    GetCorrespondingBlob: which is: %s\n", otherBlob.ToString().c_str());
		}
		if (matchScore > bestMatchScore)
		{
			bestMatchScore = matchScore;
			bestIndex = i;
		}
	}
	if (debug) logd.printf("    GetCorrespondingBlob: bestMatchScore: %.3f, to blob index %d\n", bestMatchScore, ci[bestIndex]);
	if (bestMatchScore > 0.5f)
	{
		NsiBlob& otherBlob = iset.blobSets[frameIndex].nsiBlobs[ci[bestIndex]];
		if (debug)
		{
			logd.printf("    GetCorrespondingBlob: Found blob %d-%d with best match score %.3f:\n", frameIndex, ci[bestIndex], bestMatchScore);
		}
		return ci[bestIndex];
	}
	else
	{
		if (debug) logd.printf("    GetCorrespondingBlob: None found.\n");
		return -1;
	}
}
void AnalyzeImagePairBlobSets(ImageSet& iset, int thisFrameIndex, int otherFrameIndex, CAList& cal, int minDistPixels, int maxDistPixels, DebugSpec dspec)
{
	if (dspec.debug) logd.debug("AnalyzeImagePairBlobSets: Frame pair %d - %d", thisFrameIndex, otherFrameIndex);
	float distMinSqThreshold = minDistPixels * minDistPixels;
	float distMaxSqThreshold = maxDistPixels * maxDistPixels;
	int caIndex = 0;
	int n = iset.blobSets[thisFrameIndex].size();
	for (int i = 0; i < n; i++)
	{
		NsiBlob& thisBlob = iset.blobSets[thisFrameIndex].nsiBlobs[i];
		float x = thisBlob.imgBlob.xCentroid;
		float y = thisBlob.imgBlob.yCentroid;
		bool debug = thisBlob.imgBlob.CheckDebug(dspec);
		vector<int> otherFrameCloseBlobIndices = iset.blobSets[otherFrameIndex].pointIndex.FindNearby(x, y, maxDistPixels);
		int no = otherFrameCloseBlobIndices.size();
		if (debug)
		{
			logd.printf("AnalyzeImagePairBlobSets: -----------\n");
			logd.printf("AnalyzeImagePairBlobSets: This blob: %d-%d (blob index) has %d neighbors in range:\n", thisFrameIndex, i, no);
			logd.printf("    %s\n", thisBlob.ToString().c_str());
		}
		for (int io = 0; io < no; io++) 
		{
			int j = otherFrameCloseBlobIndices[io];
			NsiBlob& otherBlob = iset.blobSets[otherFrameIndex].nsiBlobs[j];
			debug = (thisBlob.imgBlob.CheckDebug(dspec) || otherBlob.imgBlob.CheckDebug(dspec));
			float dq = thisBlob.DistSq(otherBlob);
			if ((dq >= distMinSqThreshold) && (dq < distMaxSqThreshold))
			{
				float ss = otherBlob.CompSim(thisBlob, false);
				if (ss > BLOB_SIM_SCORE_MIN_THRESHOLD)
				{
					float ox = otherBlob.imgBlob.xCentroid;
					float oy = otherBlob.imgBlob.yCentroid;
					float dx = x - ox;
					float dy = y - oy;
					if (debug)
					{
						logd.printf("AnalyzeImagePairBlobSets: This blob: %d-%d (blob index)\n", thisFrameIndex, i);
						logd.printf("AnalyzeImagePairBlobSets: Possible other blob: %d-%d dist: %d, simScore: %.3f, offset: (%.1f, %.1f)\n",
							otherFrameIndex, j, (int)sqrt(dq), ss, dx, dy);
						logd.printf("    %s\n", otherBlob.ToString().c_str());
					}
					CA ca(thisFrameIndex, thisBlob, otherFrameIndex, otherBlob);
					ca.imageSetId = iset.imageSetId;
					ca.index = cal.size();
					ca.InitBlobIndices(); 
					ca.blobIndices[thisFrameIndex] = i;
					ca.blobIndices[otherFrameIndex] = j;
					for (int k = 0; k < 6; k++)
					{
						int sign = ((k + 1) % 2) * 2 - 1;
						int ffi = (thisFrameIndex + (k + 2) / 2 * sign);
						if ((ffi < 0) || (ffi > 3) || (ffi == thisFrameIndex) || (ffi == otherFrameIndex)) continue;
						int cbi = GetCorrespondingBlob(iset, ca, ffi, debug);
						if (cbi >= 0)
						{
							ca.SetBlob(ffi, iset.blobSets[ffi].nsiBlobs[cbi]); 
							if (debug) logd.printf("AnalyzeImagePairBlobSets: Found frame %d blob at index %d\n", ffi, cbi);
							ca.blobIndices[ffi] = cbi;
						}
						else
						{
							if (debug) logd.printf("AnalyzeImagePairBlobSets: No corresponding blob found in frame %d.\n", ffi);
							ca.emBlobScores[ffi] = GetExplainedMissingBlobScore(iset, ca, ffi, debug);
							if (ca.emBlobScores[ffi] > 0.5f)
							{
								if (debug) logd.printf("AnalyzeImagePairBlobSets: Found explained-missing blob in frame %d with score %.2f.\n", ffi, ca.emBlobScores[ffi]);
							}
							else
							{
								if (debug) logd.printf("AnalyzeImagePairBlobSets: No explained-missing blob either in frame %d.\n", ffi);
							}
						}
					}
					if (ca.GetFrameCountFloat() >= MIN_N_BLOBS_PER_CA)
					{
						if (cal.Add(ca, dspec))
						{
							if (debug) logd.printf("AnalyzeImagePairBlobSets: Added a new CA: %s\n", ca.ToString().c_str());
							if (debug) ca.Display();
						}
					}
					else
					{
						if (debug) logd.printf("AnalyzeImagePairBlobSets: Too few blobs to keep this CA: %s\n", ca.ToString().c_str());
					}
				}
				else
				{
					if (debug)
					{
						logd.printf("AnalyzeImagePairBlobSets: thisBlob:  %s\n", thisBlob.ToString().c_str());
						logd.printf("AnalyzeImagePairBlobSets: otherBlob: %s\n", otherBlob.ToString().c_str());
						logd.printf("AnalyzeImagePairBlobSets: SimScore: %.3f, so other blob not similar enough\n", ss);
					}
				}
			}
		}
	}
}
void SaveCAImages(CA& ca, ImageSet& iset)
{
	for (int i = 0; i < 4; i++)
	{
		char tmp[512];
		sprintf(tmp, "ca-%d", ca.index);
		SaveDebugBlobRoi(tmp, i, 0, iset.images[i].nsiDiff, ca.blobs[i].imgBlob);
		SaveDebugBlobRoi(tmp, i, 0, iset.images[i].nsiDiff, ca.blobs[i].imgBlob, &iset.blobSets[i].labelImage);
	}
}
void CullAndFindMoreCas(ImageSet& iset, CAList& cal, int dropSmallBlobsMinPixelMean, int dropSmallBlobsMinPixelMax, int minDistPix, int maxDistPix, DebugSpec dspec)
{
	int w = iset.getWidth();
	int h = iset.getHeight();
	Img filteredImg(w, h, 8);
	for (int frameIndex = 0; frameIndex < 4; frameIndex++)
	{
		int oldBlobCount = iset.blobSets[frameIndex].size();
		if (dspec.debug) SaveDebugTiff("preCull-nsiDiff", iset.images[frameIndex].nsiDiff);
		Img8Util::DropSmallBlobs(iset.images[frameIndex].nsiDiff, filteredImg, 4, dropSmallBlobsMinPixelMean, dropSmallBlobsMinPixelMax);
		Img8Util::Copy(filteredImg, iset.images[frameIndex].nsiDiff);
		DebugSpec subspec = dspec;
		subspec.debug = false;
		iset.blobSets[frameIndex].Create(iset.images[frameIndex].nsiDiff, filteredImg, subspec); 
		iset.blobSets[frameIndex].CullNsiBlobs(filteredImg, iset.images[frameIndex].nsiDiff, subspec);
		if (dspec.debug) SaveDebugTiff("postCull-nsiDiff", iset.images[frameIndex].nsiDiff);
		if (dspec.debug) logd.debug("CullAndFindMoreCas: Culled from %d to %d blobs.", oldBlobCount, iset.blobSets[frameIndex].size());
		for (int i = 0; i < iset.blobSets[frameIndex].size(); i++)
		{
			NsiBlob& nb = iset.blobSets[frameIndex].nsiBlobs[i];
			if (nb.CheckDebug(dspec))
			{
				logd.debug("CullAndFindMoreCas: Blob %d fr %d at (%.0f, %.0f) survived cull", i, frameIndex, nb.imgBlob.xCentroid, nb.imgBlob.yCentroid);
			}
		}
	}
	int oldSize = cal.size();
	AnalyzeImagePairBlobSets(iset, 1, 0, cal, minDistPix, maxDistPix, dspec);
	if (dspec.debug) logd.debug("CullAndFindMoreCas: Frame pair 0-1 add %d CA's (for total %d)", cal.size() - oldSize, cal.size());
	oldSize = cal.size();
	AnalyzeImagePairBlobSets(iset, 2, 1, cal, minDistPix, maxDistPix, dspec);
	if (dspec.debug) logd.debug("CullAndFindMoreCas: Frame pair 1-2 add %d CA's (for total %d)", cal.size() - oldSize, cal.size());
	oldSize = cal.size();
	AnalyzeImagePairBlobSets(iset, 3, 2, cal, minDistPix, maxDistPix, dspec);
	if (dspec.debug) logd.debug("CullAndFindMoreCas: Frame pair 1-2 add %d CA's (for total %d)", cal.size() - oldSize, cal.size());
}
void FindCandidatesSimple(ImageSet& iset, CAList& cal, DebugSpec dspec)
{
	AnalyzeImagePairBlobSets(iset, 1, 0, cal, MIN_NSI_BLOB_VELOCITY_PIXELS_PER_IMG, MAX_NSI_BLOB_VELOCITY_PIXELS_PER_IMG, dspec);
	if (dspec.debug) logd.debug("FindCandidatesSimple: After frame pair 0-1 there are %d CA's.", cal.size());
	AnalyzeImagePairBlobSets(iset, 2, 1, cal, MIN_NSI_BLOB_VELOCITY_PIXELS_PER_IMG, MAX_NSI_BLOB_VELOCITY_PIXELS_PER_IMG, dspec);
	if (dspec.debug) logd.debug("FindCandidatesSimple: After frame pair 1-2 there are %d CA's.", cal.size());
	AnalyzeImagePairBlobSets(iset, 3, 2, cal, MIN_NSI_BLOB_VELOCITY_PIXELS_PER_IMG, MAX_NSI_BLOB_VELOCITY_PIXELS_PER_IMG, dspec);
	if (dspec.debug) logd.debug("FindCandidatesSimple: After frame pair 2-3 there are %d CA's.", cal.size());
	CullAndFindMoreCas(iset, cal, 30, 40, MAX_NSI_BLOB_VELOCITY_PIXELS_PER_IMG, 32, dspec);
	CullAndFindMoreCas(iset, cal, 50, 60, 32, 96, dspec);
	if (DO_SLOW_MOVERS)
	{
		FindSlowMoversTop(iset, cal, dspec);
	}
	cal.UpdateSimpleScores();
	cal.Sort();
	for (int i = 0; i < cal.cas.size(); i++)
	{
		bool debug = cal.cas[i].CheckDebug(dspec);
		if (debug)
		{
			logd.printf("Debug CA: %s\n", cal.cas[i].ToString().c_str());
			SaveCAImages(cal.cas[i], iset);
		}
	}
}
void FindBlobsTop(ImageSet& iset, DebugSpec dspec)
{
	for (int i = 0; i < 4; i++)
	{
		bool debug = dspec.CheckDebugImageFrame(iset.index, i);
		int w = iset.getWidth();
		int h = iset.getHeight();
		Img filteredImg(w,h,8);
		iset.blobSets[i].Create(iset.images[i].nsiDiff, filteredImg, dspec);
	}
}
void RenderComparison(bool useNonOverlapSet, CAList& cal, ImageSet& iset, vector<bool>& matchedDets, bool useNsiDiff, DebugSpec dspec)
{
	cal.Sort(); 
	vector<DetectionRecord>& dets = iset.origDetRecords;
	if (useNonOverlapSet)
	{
		dets = iset.detRecords;
	}
	for (int frameIndex = 0; frameIndex < 4; frameIndex++)
	{
		vector<pair<int, int>> blue, magenta, green, red;
		vector<int> blueLabels, magentaLabels, greenLabels, redLabels;
		int nd = dets.size();
		for (int i = 0; i < nd; i += 4)
		{
			if (matchedDets[i])
			{
			}
			else
			{
				magenta.push_back(pair<int, int>(roundf(dets[i + frameIndex].xReg), roundf(dets[i + frameIndex].yReg)));
				magentaLabels.push_back(dets[i].detectionNumber);
			}
		}
		int nc = cal.cas.size();
		for (int i = 0; i < nc; i++)
		{
			int detNumber = cal.cas[i].GetDetectionNumber();
			if (detNumber >= 0)
			{
				if (cal.cas[i].CheckIsBlobSet(frameIndex))
				{
					green.push_back(pair<int, int>(roundf(cal.cas[i].blobs[frameIndex].imgBlob.xCentroid), roundf(cal.cas[i].blobs[frameIndex].imgBlob.yCentroid)));
					greenLabels.push_back(cal.cas[i].index);
				}
				else if (cal.cas[i].emBlobScores[frameIndex] > 0.5f)
				{
					pair<float, float> fc = cal.cas[i].GetFittedCoords(frameIndex); 
					pair<int, int> fci = pair<int, int>(roundf(fc.first), roundf(fc.second));
					green.push_back(fci);
					greenLabels.push_back(cal.cas[i].index);
					red.push_back(fci);
					redLabels.push_back(detNumber);
				}
				int detIndex = ImageSet::GetDetIndexByDetNumber(dets, detNumber, frameIndex);
				blue.push_back(pair<int, int>(roundf(dets[detIndex].xReg), roundf(dets[detIndex].yReg)));
				blueLabels.push_back(detNumber);
			}
			else
			{
				if ((dspec.caIndex == cal.cas[i].index) || (redLabels.size() < 500))
				{
					pair<float, float> fc = cal.cas[i].GetFittedCoords(frameIndex); 
					pair<int, int> fci = pair<int, int>(roundf(fc.first), roundf(fc.second));
					red.push_back(fci);
					redLabels.push_back(cal.cas[i].index);
				}
			}
		}
		char tmp[512];
		if (!useNsiDiff && (iset.images[frameIndex].nsi.getWidth() > 0))
		{
			sprintf(tmp, "compare-nsi_%d", frameIndex);
			RenderColoredRois(iset.images[frameIndex].nsi, tmp, 16, blue, blueLabels, magenta, magentaLabels, green, greenLabels, red, redLabels);
		}
		else
		{
			sprintf(tmp, "compare-nsiDiff_%d", frameIndex);
			RenderColoredRois(iset.images[frameIndex].nsiDiff, tmp, 16, blue, blueLabels, magenta, magentaLabels, green, greenLabels, red, redLabels);
		}
	}
}
float CompareToTruth(CAList& cal, ImageSet& iset, bool useNonOverlapSet, DebugSpec dspec)
{
	bool debugThisSet = dspec.CheckDebugImageSet(iset.index);
	vector<DetectionRecord>& dets = iset.origDetRecords;
	if (useNonOverlapSet)
	{
		dets = iset.detRecords;
	}
	int nd = dets.size();
	if (nd == 0)
	{
		ErrorExit("CompareToTruth: No dets.");
		return 0.0f;
	}
	vector<bool> matched(nd); 
	for (int i = 0; i < nd; i++) matched[i] = false;
	int nc = cal.cas.size();
	vector<int> fpList;
	for (int i = 0; i < nc; i++) 
	{
		bool debugThis = (cal.cas[i].index == dspec.caIndex);
		vector<pair<float, float>> points = cal.cas[i].GetFittedCentroids();
		int detIndex = iset.FindMatchingDetSet(useNonOverlapSet, points, matched, cal.cas[i].index, debugThis);
		if (detIndex < 0)
		{
			cal.cas[i].SetDetectionNumber(-1);
			fpList.push_back(i);
			if (debugThis) logd.printf("CompareToTruth: FP: This CA was not found in det's: %s\n", cal.cas[i].ToString().c_str());
		}
		else
		{
			cal.cas[i].SetDetectionNumber(dets[detIndex].detectionNumber);
			if (debugThis || dspec.CheckDebugDet(iset.index, dets[detIndex].detectionNumber)) 
				logd.printf("CompareToTruth: TP: CA %d was matched with det number: %d\n", cal.cas[i].index, cal.cas[i].GetDetectionNumber());
		}
	}
	if (debugThisSet)
	{
		logd.printf("Compare found %d FP's.\n", fpList.size());
		logd.printf("Some of the FP's: ");
		for (int i = 0; i < MIN(10, fpList.size()); i++) logd.printf("%d, ", cal.cas[fpList[i]].index);
		logd.printf("\n");
	}
	vector<int> fnList; 
	for (int i = 0; i < nd; i += 4) 
	{
		float dist = dets[i].DistReg(dspec.debugBlobCoords);
		bool debugThis = (dist <= DEBUG_COORDS_TOLERANCE_PIX);
		if (!matched[i])
		{
			fnList.push_back(i);
			if (debugThis) logd.printf("FN: This det was not found in CA's: %s\n", dets[i].ToString().c_str());
		}
	}
	if (debugThisSet) logd.printf("Compare found %d FN's (det sets).\n", fnList.size());
	if (debugThisSet)
	{
		RenderComparison(useNonOverlapSet, cal, iset, matched, false, dspec);
		RenderComparison(useNonOverlapSet, cal, iset, matched, true, dspec);
	}
	int nds = nd / 4;
	int fn = fnList.size();
	float sensitivity = (nds - fn) / (float)nds;
	return sensitivity;
}
ImageSetResult ComputeAP(vector<CA>& cas, vector<Detection>& dets, DebugSpec dspec, bool debug)
{
	ImageSetResult imageSetResult;
	int nCA = cas.size();
	int nDets = dets.size();
	if (debug) logd.printf("ComputeAP: %d CA's and %d dets.\n", nCA, nDets);
	int numOfNEOs = 0;
	for (int i = 0; i < nDets; i++) if (dets[i].isNeoTruth) numOfNEOs++;
	if (debug) logd.printf("ComputeAP: There are %d true NEO's.\n", numOfNEOs);
	bool* matched = new bool[nDets];
	for (int i = 0; i < nDets; i++) matched[i] = false;
	double score = 0;
	double detected = 0;
	int neo_count = 0;
	double neo_detected = 0;
	int tp, fp, fn, tpNeo, fpNeo, fnNeo;
	tp = fn = fp = tpNeo = fpNeo = fnNeo = 0;
	for (int i = 0; i < nCA; i++)
	{
		CA& ca = cas[i];
		bool debugThisCA = (dspec.caIndex == ca.index);
		if (debugThisCA) logd.printf("   ComputeAP: Starting CA %d\n", ca.index);
		const vector<RaDec>& caRadecs = ca.GetRaDecs();
		ca.SetDetectionNumber(-1); 
		if (ca.isNeo) neo_count++;
		bool matchedThisCA = false;
		for (int j = 0; j < nDets; j++)
		{
			if (ca.imageSetId == dets[j].imageSetId)
			{
				bool debugThisDet = (dspec.debugDetNumber == dets[j].detectionNumber);
				Detection& det = dets[j];
				double sum = 0;
				for (int f = 0; f < 4; f++)
				{
					float dx = caRadecs[f].RaDegrees - det.recs[f].raDec.RaDegrees;
					float dy = caRadecs[f].DecDegrees - det.recs[f].raDec.DecDegrees;
					sum += dx * dx + dy * dy;
					if (debugThisDet && debugThisCA)
					{
						logd.printf("       ComputeAP: CA pt degrees: (%.6f, %.6f), det pt: (%.6f, %.6f)\n", caRadecs[f].RaDegrees, caRadecs[f].DecDegrees,
							det.recs[f].raDec.RaDegrees, det.recs[f].raDec.DecDegrees);
						logd.printf("       ComputeAP: dx/dy degrees: (%.6f, %.6f), dx/dy pixels: %.1f, %.1f   -> dqs: %.5f\n", dx, dy, dx * MEAN_PX_PER_DEGREE, dy * MEAN_PX_PER_DEGREE, dx*dx + dy*dy);
					}
				}
				if (debugThisDet && debugThisCA)
				{
					logd.printf("    ComputeAP: Dist sum from CA %d to det %d is %.4f vs threshold %.4f\n", ca.index, det.detectionNumber, sum, SCORING_DIST_THRESHOLD_DEG);
				}
				if (sum < SCORING_DIST_THRESHOLD_DEG)
				{
					if (!matched[j])
					{
						if (debug || debugThisCA || debugThisDet) logd.printf("ComputeAP: TP CA: %s: CA %d matched det %d\n", ca.imageSetId.c_str(), ca.index, j);
						matchedThisCA = true;
						ca.SetDetectionNumber(dets[j].detectionNumber);
						matched[j] = true;
						detected += 1.0;
						score += (1000000.0 / nDets) * (detected / (i + 1));;
						if (det.isNeoTruth)
						{
							imageSetResult.nMatchedNeo++;
							if (ca.isNeo)
							{
								neo_detected += 1.0;
								score += (100000.0 / numOfNEOs) * (neo_detected / neo_count);
								if (debug || debugThisCA) logd.printf("ComputeAP: Neo bonus: CA %d to det %d\n", ca.index, j);
								tpNeo++;
							}
							else
							{
								fnNeo++;
							}
						}
						else if (ca.isNeo)
						{
							fpNeo++;
						}
						det.isFound = true;
						det.isNeoFound = ca.isNeo;
						det.rank = i;
						break;
					}
					else
					{
						if (debugThisDet && debugThisCA)
						{
							logd.printf("ComputeAP: For CA %d the det %d is already matched despite being a close neighbor.\n", ca.index, dets[j].detectionNumber, sum);
						}
					}
				}
				else if (sum < SCORING_DIST_THRESHOLD_DEG * 10)
				{
					if (debugThisCA || debugThisDet) logd.printf("ComputeAP: Too far from CA %d to det %d. Sum is: %.4f\n", ca.index, j, sum);
				}
			}
		}
		if (!matchedThisCA)
		{
			fp++;
			if (ca.isNeo)
			{
				if (debugThisCA) logd.printf("ComputeAP: CA: %d is FP and FP Neo\n", ca.index);
				fpNeo++;
			}
			else
			{
				if (debugThisCA) logd.printf("ComputeAP: CA: %d is FP\n", ca.index);
			}
		}
		else
		{
			tp++;
		}
	}
	for (int i = 0; i < nDets; i++)
	{
		if (!matched[i])
		{
			fn++;
			if (dspec.debugDetNumber == dets[i].detectionNumber)
			{
				logd.printf("ComputeAP: Det was not matched: %s\n", dets[i].ToString().c_str());
			}
		}
	}
	delete[] matched;
	imageSetResult.AP = (float)score;
	imageSetResult.nFN = fn;
	imageSetResult.nFP = fp;
	imageSetResult.nTP = tp;
	imageSetResult.nFNNeo = numOfNEOs - tpNeo;
	imageSetResult.nFPNeo = fpNeo;
	imageSetResult.nTPNeo = tpNeo;
	return imageSetResult;
}
ImageSetResult ComputeAPWithoutRaDec(vector<CAFeat>& cas, vector<Detection>& dets, DebugSpec dspec, bool debug)
{
	ImageSetResult imageSetResult;
	int nCA = cas.size();
	int nDets = dets.size();
	if (debug) logd.printf("ComputeAP: %d CA's and %d dets.\n", nCA, nDets);
	int numOfNEOs = 0;
	for (int i = 0; i < nDets; i++) if (dets[i].isNeoTruth) numOfNEOs++;
	if (debug) logd.printf("ComputeAP: There are %d true NEO's.\n", numOfNEOs);
	bool* matched = new bool[nDets];
	for (int i = 0; i < nDets; i++) matched[i] = false;
	double score = 0;
	double detected = 0;
	int neo_count = 0;
	double neo_detected = 0;
	int tp, fp, fn, tpNeo, fpNeo, fnNeo;
	tp = fn = fp = tpNeo = fpNeo = fnNeo = 0;
	for (int i = 0; i < nCA; i++)
	{
		CAFeat& ca = cas[i];
		bool debugThisCA = (dspec.caIndex == ca.index);
		if (debugThisCA) logd.printf("   ComputeAP: Starting CA %d\n", ca.index);
		if (ca.isNeoComp) neo_count++;
		bool matchedThisCA = false;
		for (int j = 0; j < nDets; j++)
		{
			if (ca.imageSetId == dets[j].imageSetId)
			{
				bool debugThisDet = (dspec.debugDetNumber == dets[j].detectionNumber);
				if (!matched[j])
				{
					Detection& det = dets[j];
					if (ca.corrDetNumber == det.detectionNumber)
					{
						matchedThisCA = true;
						ca.SetDetectionNumber(dets[j].detectionNumber);
						matched[j] = true;
						detected += 1.0;
						float add = (1000000.0 / nDets) * (detected / (i + 1));
						score += add;
						if (debug || debugThisCA || debugThisDet)
							logd.printf("ComputeAP: TP CA: %s: CA %d matched det %d at index %d with add %f to total %f\n", ca.imageSetId.c_str(), ca.index, j, i, add, score);
						if (det.isNeoTruth)
						{
							if (ca.isNeoComp)
							{
								neo_detected += 1.0;
								score += (100000.0 / numOfNEOs) * (neo_detected / neo_count);
								if (dspec.info || debug || debugThisCA || debugThisDet) logd.printf("ComputeAP: Neo bonus: CA %d to det %d\n", ca.index, j);
								tpNeo++;
							}
							else
							{
								fnNeo++;
							}
						}
						else if (ca.isNeoComp)
						{
							fpNeo++;
						}
						det.isFound = true;
						det.isNeoFound = ca.isNeoComp;
						det.rank = i;
						break;
					}
				}
				else
				{
					if (debugThisDet && debugThisCA)
					{
						logd.printf("ComputeAP: For CA %d the det %d is already matched.\n", ca.index, dets[j].detectionNumber);
					}
				}
			}
		}
		if (!matchedThisCA)
		{
			fp++;
			if (ca.isNeoComp)
			{
				if (debugThisCA) logd.printf("ComputeAP: CA: %d is FP and FP Neo\n", ca.index);
				fpNeo++;
			}
			else
			{
				if (debugThisCA) logd.printf("ComputeAP: CA: %d is FP\n", ca.index);
			}
		}
		else
		{
			tp++;
		}
	}
	for (int i = 0; i < nDets; i++)
	{
		if (!matched[i])
		{
			fn++;
			if (dspec.debugDetNumber == dets[i].detectionNumber)
			{
				logd.printf("ComputeAP: Det was not matched: %s\n", dets[i].ToString().c_str());
			}
		}
		else
		{
			if (dspec.debugDetNumber == dets[i].detectionNumber)
			{
				logd.printf("ComputeAP: Det was matched: %s\n", dets[i].ToString().c_str());
			}
		}
	}
	delete[] matched;
	imageSetResult.AP = (float)score;
	imageSetResult.nFN = fn;
	imageSetResult.nFP = fp;
	imageSetResult.nTP = tp;
	imageSetResult.nFNNeo = numOfNEOs - tpNeo;
	imageSetResult.nFPNeo = fpNeo;
	imageSetResult.nTPNeo = tpNeo;
	return imageSetResult;
}
ImageSetResult ComputeAP(CAList& cal, ImageSet& iset, DebugSpec dspec, bool useOverlapSet)
{
	bool debug = false; 
	ImageSetResult isr;
	if (useOverlapSet)
	{
		isr = ComputeAP(cal.cas, iset.detectionsWithOverlaps, dspec, debug);
		isr.nOverlaps = iset.GetOverlapDetsCount();
	}
	else
	{
		isr = ComputeAP(cal.cas, iset.detectionsWithoutOverlaps, dspec, debug);
	}
	isr.imageSetId = iset.imageSetId;
	isr.imageSetIndex = iset.index;
	return isr;
}
string CAToAnswer(CA& ca)
{
	char tmp[512];
	string s = ca.imageSetId;
	vector<RaDec> radecs = ca.GetRaDecs();
	for (int j = 0; j < 4; j++)
	{
		sprintf(tmp, " %f %f", radecs[j].RaDegrees, radecs[j].DecDegrees);
		s += (string)tmp;
	}
	if (ca.isNeo)
	{
		s += " 1";
	}
	else
	{
		s += " 0";
	}
	return s;
}
vector<string> BuildAnswers(CAList& cal, DebugSpec dspec)
{
	vector<string> a;
	int n = cal.size();
	for (int i = 0; i < n; i++)
	{
		string s = CAToAnswer(cal.cas[i]);
		a.push_back(s);
	}
	if (dspec.debug) logd.info("Built %d answers.", a.size());
	return a;
}
CA ParseAnswer(string s)
{
	vector<string> words;
	StringUtils::splitString(s, ' ', words);
	if (words.size() != 10)
	{
		ErrorExit("Wrong number of words in answer.");
	}
	int i = 0;
	CA ca;
	ca.imageSetId = words[i++];
	for (int j = 0; j < 4; j++)
	{
		RaDec rd;
		rd.RaDegrees = atof(words[i++].c_str());
		rd.DecDegrees = atof(words[i++].c_str());
		ca.SetRaDec(j, rd);
	}
	ca.isNeo = (atoi(words[i++].c_str()) != 0);
	return ca;
}
void ParseAnswers(vector<string> answers, vector<CA>& allCAs)
{
	for (int i = 0; i < answers.size(); i++)
	{
		CA ca = ParseAnswer(answers[i]);
		ca.index = i;
		allCAs.push_back(ca);
	}
}
void TestParseAnswer()
{
	string s = "IS_000 27.115137 3.734035 27.115327 3.734889 27.115024 3.734315 27.115446 3.734419 0";
	CA ca = ParseAnswer(s);
	logd.printf("Parsed CA: %s\n", CAToAnswer(ca).c_str());
}
float ComputeSboProxScore(int caIndex, Img& sboMask, Img& sboLabel, ImgBlob& sboBlob, NsiBlob& casBlob, bool debug)
{
	float roiDist = sboBlob.roi.dist(casBlob.imgBlob.roi);
	int sboRadius = roundf(sqrtf(sboBlob.area / PI_FLOAT)); 
	int casBlobRadius = roundf(sqrtf(casBlob.imgBlob.area / PI_FLOAT)); 
	float radialDist = pyth(sboBlob.xCentroid - casBlob.imgBlob.xCentroid, sboBlob.yCentroid - casBlob.imgBlob.yCentroid);
	radialDist = radialDist - sboRadius - casBlobRadius; 
	float dist = MAX(radialDist, roiDist);
	float score;
	sboRadius = (sboRadius * CA_SBO_RADIUS_FACTOR); 
	if (dist > sboRadius)
	{
		score = 0; 
	}
	else if (dist < 2)
	{
		score = 1; 
	}
	else
	{
		score = (sboRadius - dist) / sboRadius;
	}
	if (debug) logd.debug("ComputeSboProxScore: CA: %d, sbo: (%.0f, %.0f), sboRadius: %d, casBlobRadius: %d, roiDist: %.1f, radialDist: %.1f, dist: %.1f, score: %.3f",
		caIndex, sboBlob.xCentroid, sboBlob.yCentroid, sboRadius, casBlobRadius, roiDist, radialDist, dist, score);
	return score;
}
vector<float> ComputeSboProxScores(int frameIndex, Img& sboMask, vector<CA>& cas, DebugSpec dspec)
{
	int nca = cas.size();
	int w = sboMask.getWidth();
	int h = sboMask.getHeight();
	int nEighth = (w * h) / 8;
	Img label(w, h, 16);
	ImgUtil::Label(sboMask, label, true);
	vector<ImgBlob> blobs;
	ImgUtil::FindBlobs(label, blobs);
	vector<float> scores(nca);
	if (blobs.size() > 0)
	{
		int maxArea = 0;
		for (int i = 0; i < blobs.size(); i++)
		{
			maxArea = MAX(maxArea, blobs[i].area);
		}
		int maxRadius = (int)(sqrtf(maxArea / PI_FLOAT) * 1.5f); 
		int gridSize = MIN(maxRadius, 128);
		PointIndex index(MAX(w, h), gridSize);
		ImgBlob::AddAllToPointIndex(blobs, index);
		int searchRadius = MIN(512, maxRadius * 2); 
		for (int i = 0; i < nca; i++)
		{
			bool debug = cas[i].CheckDebug(dspec);
			if (cas[i].CheckIsBlobSet(frameIndex))
			{
				if (debug) logd.debug("ComputeSboProxScores: Frame %d: Max sbo radius %d on max area %d", frameIndex, maxRadius, maxArea);
				vector<int> neighbors = index.FindNearby(cas[i].blobs[frameIndex].imgBlob.xCentroid, cas[i].blobs[frameIndex].imgBlob.yCentroid, searchRadius, false);
				if (debug) logd.debug("ComputeSboProxScores: CA %d: Found %d neighbor sbo's on search radius %d",
					cas[i].index, neighbors.size(), searchRadius);
				float proxScore = 0;
				if (neighbors.size() > 0)
				{
					for (int j = 0; j < neighbors.size(); j++)
					{
						if (blobs[neighbors[j]].area > nEighth)
						{
							if (debug) logd.debug("ComputeSboProxScores: CA %d: Neighbor has huge area %d, probably junk, ignoring.",
								cas[i].index, blobs[neighbors[j]].area);
						}
						else
						{
							proxScore = MAX(proxScore, ComputeSboProxScore(cas[i].index, sboMask, label, blobs[neighbors[j]], cas[i].blobs[frameIndex], false));
						}
					}
				}
				scores[i] = proxScore;
			}
			else
			{
				if (debug) logd.debug("ComputeSboProxScores: CA %d: Frame %d CA blob not set.", cas[i].index, frameIndex);
				scores[i] = -1;
			}
		}
	}
	else
	{
		ErrorExit("No blobs found.");
	}
	return scores;
}
void SetSboProxScores(ImageSet& imageSet, CAList& cal, DebugSpec dspec)
{
	int nca = cal.size();
	if (nca == 0)
	{
		logd.info("SetSboProxScores: No CA's to set scores on.");
		return;
	}
	vector<vector<float>> framesScores(4); 
	for (int fi = 0; fi < 4; fi++)
	{
		int nnz = Img8Util::CountNonZero(imageSet.images[fi].sboMask);
		int w = imageSet.images[fi].sboMask.getWidth();
		int h = imageSet.images[fi].sboMask.getHeight();
		if (nnz == 0)
		{
			logd.info("SetSboProxScores: No SBO's to use for scores.");
			framesScores[fi].resize(nca); 
		}
		else if (nnz > (w * h / 4))
		{
			logd.warn("SetSboProxScores: More than 25 percent of image %s fr %d is sbo masked, bailing on sbo scoring this frame.", imageSet.imageSetId.c_str(), fi);
			framesScores[fi].resize(nca); 
		}
		else
		{
			framesScores[fi] = ComputeSboProxScores(fi, imageSet.images[fi].sboMask, cal.cas, dspec);
		}
	}
	vector<float> scores(4);
	for (int i = 0; i < nca; i++)
	{
		bool debug = cal.cas[i].CheckDebug(dspec);
		for (int j = 0; j < 4; j++)
		{
			scores[j] = framesScores[j][i];
		}
		vector<float> goodScores = exceptv(scores, -1);
		cal.cas[i].sboProxScore = meanv(goodScores);
		if (debug) logd.debug("SetSboProxScores: CA %d: scores: %.3f, %.3f, %.3f, %.3f  -> sboProxScore: %.3f",
			cal.cas[i].index, scores[0], scores[1], scores[2], scores[3], cal.cas[i].sboProxScore);
	}
}
vector<float> ComputeCAProxScores(int frameIndex, int width, int height, vector<CA>& cas, DebugSpec dspec)
{
	int nca = cas.size();
	vector<float> scores(nca);
	int searchRadius = CA_PROX_SEARCH_RADIUS_PIX;
	int gridSize = CA_PROX_SEARCH_RADIUS_PIX; 
	PointIndex index(MAX(width, height), gridSize);
	CA::AddAllToPointIndex(cas, frameIndex, index);
	for (int i = 0; i < nca; i++)
	{
		bool debug = cas[i].CheckDebug(dspec);
		if (cas[i].CheckIsBlobSet(frameIndex))
		{
			vector<int> neighbors = index.FindNearby(cas[i].blobs[frameIndex].imgBlob.xCentroid, cas[i].blobs[frameIndex].imgBlob.yCentroid, searchRadius, false);
			float proxScore = 0;
			if (neighbors.size() == 0)
			{
				ErrorExit("ComputeCAProxScores: CA: %d, got 0 neighbors when self is also in index.", cas[i].index);
			}
			else
			{
				proxScore = neighbors.size() - 1; 
			}
			scores[i] = proxScore;
		}
		else
		{
			if (debug) logd.debug("ComputeCAProxScores: CA %d: Frame %d CA blob not set.", cas[i].index, frameIndex);
			scores[i] = -1;
		}
	}
	return scores;
}
void SetCAProxScores(ImageSet& imageSet, CAList& cal, DebugSpec dspec)
{
	int nca = cal.size();
	if (nca == 0)
	{
		logd.info("SetCAProxScores: No CA's to set scores on.");
		return;
	}
	vector<vector<float>> framesScores(4);
	for (int frameIndex = 0; frameIndex < 4; frameIndex++)
	{
		if (nca < 2)
		{
			logd.info("SetCAProxScores: less than 2 CA's to use, leave CA prox score 0.");
			framesScores[frameIndex].resize(nca); 
		}
		else
		{
			int w = imageSet.getWidth();
			int h = imageSet.getHeight();
			framesScores[frameIndex] = ComputeCAProxScores(frameIndex, w, h, cal.cas, dspec);
		}
	}
	vector<float> scores(4);
	for (int i = 0; i < nca; i++) 
	{
		bool debug = cal.cas[i].CheckDebug(dspec);
		for (int frameIndex = 0; frameIndex < 4; frameIndex++) 
		{
			scores[frameIndex] = framesScores[frameIndex][i];
		}
		vector<float> goodScores = exceptv(scores, -1);
		cal.cas[i].caProxScore = meanv(goodScores);
		if (debug) logd.debug("SetCAProxScores: CA %d: scores: %.3f, %.3f, %.3f, %.3f  -> caProxScore: %.3f",
			cal.cas[i].index, scores[0], scores[1], scores[2], scores[3], cal.cas[i].caProxScore);
	}
}
void PostProcessCandidates(ImageSet& imageSet, CAList& cal, DebugSpec dspec)
{
	SetSboProxScores(imageSet, cal, dspec);
	SetCAProxScores(imageSet, cal, dspec);
	cal.UpdateNCAsCount(); 
	cal.ApplyPostProcessingScores(dspec);
}
#include <random>
#include <cstdint>
class RNG
{
private:
	std::mt19937 mt;
	unsigned int randSub()
	{
		return mt();
	}
public:
	RNG(int seed = 1)
	{
		init(seed);
	}
	void init(int seed = 1)
	{
		mt.seed(seed);
	}
	int next()
	{
		return randSub();
	}
	int next(int x)
	{
		return randSub() % x;
	}
	int next(int a, int b)
	{
		return a + (randSub() % (b - a));
	}
	double nextDouble()
	{
		return (rand() + 0.5) * (1.0 / 4294967296.0);
	}
};
#include <vector>
#include <thread>
#include <mutex>
#include <algorithm>
#include <set>
//rks: Based on Psyho's MM1 classifier
#pragma region RandomForestParams
struct RandomForestParams
{
	int randSeed;
	float timeLimitSeconds;
	int nTrees;
	int maxDepth;
	int maxNodeSize;
	float bootstrapRatio;
	int nRandFeat;
	int nRandPos;
	RandomForestParams()
	{
		randSeed = 1;
		timeLimitSeconds = 0.0f;
		nTrees = 200;
		maxDepth = 25;
		maxNodeSize = 1;
		bootstrapRatio = 1.2f;
		nRandFeat = 2;
		nRandPos = 2;
	}
	void Display()
	{
		logd.printf("        randSeed: %d\n", randSeed);
		logd.printf("        nTrees: %d\n", nTrees);
		logd.printf("        maxDepth: %d\n", maxDepth);
		logd.printf("        maxNodeSize: %d\n", maxNodeSize);
		logd.printf("        bootstrapRatio: %.3f\n", bootstrapRatio);
		logd.printf("        timeLimitSeconds: %.1f\n", timeLimitSeconds);
		logd.printf("        nRandFeat: %d\n", nRandFeat);
		logd.printf("        nRandPos: %d\n", nRandPos);
	}
};
#pragma endregion
#pragma region RandomForestTreeNode
	struct RandomForestTreeNode
	{
		int level; 
		int featureIndex; 
		float splitVal; 
		float result; 
		int left;
		int right;
		float error; 
		int padding; 
		RandomForestTreeNode()
		{
			level = -1;
			featureIndex = -1;
			splitVal = 0;
			result = 0;
			left = -1;
			right = -1;
			error = 0;
		}
	};
#pragma endregion
#pragma region RandomForestDecisionTree
	struct RandomForestDecisionTree
	{
		vector<RandomForestTreeNode> nodes;
		static float Gini(int n0, int n1)
		{
			int n = n0 + n1;
			float p0 = ((float)n0 / n);
			float p1 = ((float)n1 / n);
			float sum = p0 * (1.0f - p0) + p1 * (1.0f - p1);
			return sum;
		}
		static float Gini(int n0Left, int n1Left, int n0Right, int n1Right)
		{
			float gLeft = Gini(n0Left, n1Left);
			float gRight = Gini(n0Right, n1Right);
			int nleft = n0Left + n1Left;
			int nright = n0Right + n1Right;
			float g = (gLeft * nleft + gRight * nright) / (nleft + nright);
			return g;
		}
		RandomForestDecisionTree(vector<vector<float>> &samples, vector<int> &truthClasses, RandomForestParams &params, int seed)
		{
			RNG rng(seed);
			int nSamples = samples.size();
			int nf = samples[0].size();
			int nSamplesChosen = (int)(nSamples * params.bootstrapRatio);
			vector<int> chosenSampleIndices;
			chosenSampleIndices.reserve(nSamplesChosen);
			for (int i = 0; i < nSamplesChosen; i++)
			{
				chosenSampleIndices.push_back(rng.next(nSamples));
			}
			RandomForestTreeNode root;
			root.level = 0;
			root.left = 0;
			root.right = nSamplesChosen;
			nodes.push_back(root);
			vector<int> nodeIndicesStack; 
			nodeIndicesStack.push_back(0);
			while (nodeIndicesStack.size())
			{
				int currentNode = nodeIndicesStack.back(); 
				nodeIndicesStack.pop_back();
				int leftIndex = nodes[currentNode].left;
				int rightIndex = nodes[currentNode].right;
				int nodeSize = rightIndex - leftIndex;
				int truthSumThisNode = 0;
				for (int i = leftIndex; i < rightIndex; i++)
				{
					truthSumThisNode += truthClasses[chosenSampleIndices[i]];
				}
				bool allEqual = ((truthSumThisNode == nodeSize) || (truthSumThisNode == 0));
				if (allEqual || (nodeSize <= params.maxNodeSize) || (nodes[currentNode].level >= params.maxDepth))
				{
					nodes[currentNode].result = (float)truthSumThisNode / nodeSize;
					continue; 
				}
				int bestFeature = -1;
				int bestNSmallerThanSplitValue = 0;
				int bestNGreaterThanSplitValue = 0;
				float bestSplitValue = 0;
				float bestError = FLT_MAX;
				for (int i = 0; (i < params.nRandFeat) && (bestError > 0); i++)
				{
					int featureIndex = rng.next(samples[0].size());
					float minValue, maxValue;
					minValue = FLT_MAX;
					maxValue = -FLT_MAX;
					for (int j = leftIndex; j < rightIndex; j++)
					{
						minValue = MIN(minValue, samples[chosenSampleIndices[j]][featureIndex]);
						maxValue = MAX(maxValue, samples[chosenSampleIndices[j]][featureIndex]);
					}
					if (minValue == maxValue) continue; 
					int nPosTries = MIN(params.nRandPos + nodes[currentNode].level / 2, MAX(params.nRandPos, nodeSize / 2));
					for (int j = 0; (j < nPosTries) && (bestError > 0); j++)
					{
						int selOffset = rng.next(nodeSize);
						float splitValue = samples[chosenSampleIndices[leftIndex + selOffset]][featureIndex];
						float sumLeft = 0, sumRight = 0;
						int nSmallerThanSplitValue = 0, nGreaterThanSplitValue = 0;
						int countLeft = 0, countRight = 0;
						for (int k = leftIndex; k < rightIndex; k++)
						{
							int p = chosenSampleIndices[k];
							if (samples[p][featureIndex] < splitValue)
							{
								sumLeft += truthClasses[p];
								countLeft++;
								nSmallerThanSplitValue++;
							}
							else
							{
								sumRight += truthClasses[p];
								countRight++;
								nGreaterThanSplitValue++;
							}
						}
						if ((nSmallerThanSplitValue == 0) || (nGreaterThanSplitValue == 0)) continue;
						int n0Left = countLeft - sumLeft;
						int n1Left = sumLeft;
						int n0Right = countRight - sumRight;
						int n1Right = sumRight;
						float error = Gini(n0Left, n1Left, n0Right, n1Right);
						if (error < bestError)
						{
							bestError = error;
							bestFeature = featureIndex;
							bestSplitValue = splitValue;
							bestNSmallerThanSplitValue = nSmallerThanSplitValue;
							bestNGreaterThanSplitValue = nGreaterThanSplitValue;
						}
					} 
				} 
				if ((bestNSmallerThanSplitValue == 0) || (bestNGreaterThanSplitValue == 0))
				{
					nodes[currentNode].result = (float)truthSumThisNode / nodeSize;
					continue;
				}
				nodes[currentNode].featureIndex = bestFeature;
				float nextLowerValue = -FLT_MAX;
				for (int i = leftIndex; i < rightIndex; i++)
				{
					float featureVal = samples[chosenSampleIndices[i]][bestFeature];
					if (featureVal < bestSplitValue)
					{
						nextLowerValue = max(nextLowerValue, featureVal);
					}
				}
				nodes[currentNode].splitVal = (bestSplitValue + nextLowerValue) / 2.0;
				if (nodes[currentNode].splitVal <= nextLowerValue)
				{
					nodes[currentNode].splitVal = bestSplitValue;
				}
				RandomForestTreeNode leftChild, rightChild;
				int newLevel = nodes[currentNode].level + 1;
				leftChild.level = newLevel;
				leftChild.error = bestError;
				rightChild.level = newLevel;
				rightChild.error = bestError;
				nodes[currentNode].left = nodes.size();
				nodes[currentNode].right = nodes.size() + 1;
				int featIdx = nodes[currentNode].featureIndex;
				float splitVal = nodes[currentNode].splitVal;
				auto mi = std::partition(chosenSampleIndices.begin() + leftIndex, chosenSampleIndices.begin() + rightIndex,
					[&](int xi) { return samples[xi][featIdx] < splitVal; });
				int middleIndex = (mi - chosenSampleIndices.begin());
				leftChild.left = leftIndex;
				leftChild.right = middleIndex;
				rightChild.left = middleIndex;
				rightChild.right = rightIndex;
				nodeIndicesStack.push_back(nodes.size());
				nodeIndicesStack.push_back(nodes.size() + 1);
				nodes.push_back(leftChild);
				nodes.push_back(rightChild);
			}
		}
		float Eval(vector<float> &samples)
		{
			int ni = 0;
			while (true)
			{
				if (nodes[ni].featureIndex < 0) return nodes[ni].result;
				ni = (samples[nodes[ni].featureIndex] < nodes[ni].splitVal ? nodes[ni].left : nodes[ni].right);
			}
		}
	};
#pragma endregion
#pragma region RandomForest
	class RandomForest
	{
	private:
		RNG rng;
		vector<RandomForestDecisionTree> trees;
		RandomForestParams params;
	public:
		void SetParams(RandomForestParams p)
		{
			params = p;
		}
		void Train(vector<vector<float>>& samples, vector<int>& truthClasses, bool debug = false)
		{
			if (trees.size() > 0)
			{
				ErrorExit("This can only be called once.");
			}
			else if (samples.size() != truthClasses.size())
			{
				ErrorExit("Different counts of samples and truthClasses.");
			}
			else if (samples.size() == 0)
			{
				ErrorExit("No samples to train on.");
			}
			else
			{
				for (int i = 0; i < truthClasses.size(); i++)
				{
					if ((truthClasses[i] < 0) || (truthClasses[i] > 1))
					{
						ErrorExit("Bad truth class value %d at index %d", truthClasses[i], i);
					}
				}
				int nt = params.nTrees;
				float startTime = getTimeSeconds();
				rng.init(params.randSeed);
				for (int i = 0; i < nt; i++)
				{
					if ((params.timeLimitSeconds > 0) && ((getTimeSeconds() - startTime) > params.timeLimitSeconds))
					{
						break;
					}
					trees.push_back(RandomForestDecisionTree(samples, truthClasses, params, rng.next()));
				}
			}
		}
		void Eval(vector<vector<float>>& samples, vector<float>& pClass1, bool debug = false)
		{
			int n = samples.size();
			if (n > 0)
			{
				int nf = samples[0].size();
				int nt = trees.size();
				vector<float> values(nf);
				for (int s = 0; s < n; s++)
				{
					float sum = 0;
					for (int i = 0; i < nt; i++) sum += trees[i].Eval(samples[s]);
					if (sum < 0)
					{
						ErrorExit("Eval: Got sum of tree's responses lt 0 for sample %d", s);
					}
					float score = sum / nt;
					pClass1[s] = score;
				}
			}
		}
		void Display()
		{
			int nt = trees.size();
			int nn = 0;
			for (int i = 0; i < nt; i++) nn += trees[i].nodes.size();
			logd.printf("    treeCount: %d\n", nt);
			logd.printf("    nodeCount: %d\n", nn);
			params.Display();
		}
	};
#pragma endregion
RandomForestParams CreateRFParams()
{
	RandomForestParams p;
	p.nTrees = RF_N_TREES;
	p.randSeed = 1;
	p.timeLimitSeconds = TIME_LIMIT_MIN_TRAIN_RF * 60.0f;
	return p;
}
void TrainCAClassifier(vector<vector<float>>& samples, vector<int>& truthClasses, RandomForest& rf, DebugSpec dspec)
{
	if (samples.size() != truthClasses.size())
	{
		ErrorExit("TrainCAClassifier: There are %d samples and %d truth classes.", samples.size(), truthClasses.size());
		return;
	}
	int ns = samples.size();
	for (int i = 0; i < ns; i++)
	{
		if (truthClasses[i] < 0)
		{
			ErrorExit("Got a truth class lt 0 at index %d", i);
		}
	}
	float st = getTimeSeconds();
	rf.Train(samples, truthClasses, dspec.debugc);
	if (dspec.info) logd.info("TrainCAClassifier: Done training classifier with %d CA's in %.1f seconds", samples.size(), getTimeSeconds() - st);
	if (dspec.debugc)
	{
	}
}
void TrainCAClassifier(vector<CA>& cas, RandomForest& rf, DebugSpec dspec)
{
	vector<vector<float>> samples;
	vector<int> truthClasses;
	CAList::GetSamples(cas, samples, truthClasses); 
	TrainCAClassifier(samples, truthClasses, rf, dspec);
}
int GetNeoTruthClasses(vector<CA>& cas, vector<Detection>& dets, vector<int>& neoTruthClasses, DebugSpec dspec)
{
	int n = cas.size();
	neoTruthClasses.reserve(n);
	int neoFoundCount = 0;
	vector<bool> matchedNeo(dets.size());
	for (int i = 0; i < n; i++)
	{
		int ntc = 0;
		if (cas[i].corrDetNumber >= 0)
		{
			for (int j = 0; j < dets.size(); j++)
			{
				if ((cas[i].corrDetNumber == dets[j].detectionNumber) && (cas[i].imageSetId == dets[j].imageSetId))
				{
					if (dets[j].isNeoTruth)
					{
						ntc = 1;
						if (dspec.debug) logd.info("GetNeoTruthClasses: Found a true neo: %s: CA %d, det %d", cas[i].imageSetId.c_str(), cas[i].index, dets[j].detectionNumber);
						neoFoundCount++;
						matchedNeo[j] = true;
					}
				}
			}
		}
		neoTruthClasses.push_back(ntc);
	}
	int neoTruthCount = 0;
	for (int j = 0; j < dets.size(); j++)
	{
		if (dets[j].isNeoTruth)
		{
			neoTruthCount++;
			if (!matchedNeo[j])
			{
				if (dspec.debug) logd.info("GetNeoTruthClasses: No CA for NEO %s det %d", dets[j].imageSetId.c_str(), dets[j].detectionNumber);
			}
		}
	}
	if (dspec.info) logd.info("GetNeoTruthClasses: CA's are matched with %d of %d NEO's in the det list.", neoFoundCount, neoTruthCount);
	return neoFoundCount;
}
void TrainRFNeoClassifier(vector<CA>& cas, vector<Detection>& dets, RandomForest& rfNeo, DebugSpec dspec)
{
	vector<int> neoTruthClasses;
	int nTruthNeo = GetNeoTruthClasses(cas, dets, neoTruthClasses, dspec);
	if (nTruthNeo > 0)
	{
		if (dspec.info) logd.info("TrainRFNeoClassifier: %d NEO's to train on.", nTruthNeo);
		vector<vector<float>> samples;
		vector<int> unusedTruthClasses;
		CAList::GetSamples(cas, samples, unusedTruthClasses);
		TrainCAClassifier(samples, neoTruthClasses, rfNeo, dspec);
	}
	else
	{
		if (dspec.info) logd.info("TrainRFNeoClassifier: No NEO's to train on.");
	}
}
int GetOverlapTruthClasses(vector<CA>& cas, vector<Detection>& dets, vector<int>& overlapTruthClasses, DebugSpec dspec)
{
	int n = cas.size();
	overlapTruthClasses.resize(n);
	int nOverlaps = 0;
	for (int i = 0; i < n; i++)
	{
		overlapTruthClasses[i] = 0;
		if (cas[i].corrDetNumber >= 0)
		{
			bool foundDet = false;
			for (int j = 0; j < dets.size(); j++)
			{
				if (cas[i].corrDetNumber == dets[j].detectionNumber)
				{
					if (dets[j].hasOverlap)
					{
						overlapTruthClasses[i] = 1;
						nOverlaps++;
						if (dspec.caIndex == cas[i].index)
							logd.printf("GetOverlaptruthClasses: CA %d has overlap like det %d\n", cas[i].index, dets[j].detectionNumber);
					}
					foundDet = true;
				}
			}
			if (!foundDet)
			{
				logd.warn("GetOverlaptruthClasses: Did not find det %d for CA %d", cas[i].corrDetNumber, cas[i].index);
			}
		}
	}
	return nOverlaps;
}
void TrainRFOverlapClassifier(vector<CA>& cas, vector<Detection>& dets, RandomForest& rfOverlap, DebugSpec dspec)
{
	vector<int> overlapTruthClasses;
	int nTruthOverlap = GetOverlapTruthClasses(cas, dets, overlapTruthClasses, dspec);
	if (nTruthOverlap > 0)
	{
		if (dspec.info) logd.info("TrainRFOverlapClassifier: %d Overlap's to train on.", nTruthOverlap);
		vector<vector<float>> samples;
		vector<int> unusedTruthClasses;
		CAList::GetSamples(cas, samples, unusedTruthClasses);
		TrainCAClassifier(samples, overlapTruthClasses, rfOverlap, dspec);
	}
	else
	{
		if (dspec.info) logd.info("TrainRFOverlapClassifier: No Overlap's to train on.");
	}
}
void ClassifyCAs(vector<vector<float>>& samples, RandomForest& rf, vector<float>& pClass1, DebugSpec dspec)
{
	int ns = samples.size();
	float st = getTimeSeconds();
	pClass1.resize(ns);
	rf.Eval(samples, pClass1, dspec.debugc);
	if (dspec.debug) logd.info("ClassifySamples: Classified %d CA's in %.1f seconds", ns, getTimeSeconds() - st);
}
void ClassifyCAs(vector<CA>& cas, RandomForest& rf, DebugSpec dspec)
{
	vector<vector<float>> samples;
	vector<int> truthClasses;
	CAList::GetSamples(cas, samples, truthClasses); 
	vector<float> pClass1(cas.size());
	ClassifyCAs(samples, rf, pClass1, dspec);
	for (int i = 0; i < cas.size(); i++)
	{
		if (pClass1[i] < 0)
		{
			ErrorExit("ClassifyCAs: Got a lt 0 result in pClass1. Index %d", i);
			cas[i].pClass1 = 0.0f;
		}
		else
		{
			cas[i].pClass1 = pClass1[i];
		}
	}
}
void ClassifyCAsRFNeo(int imageSetIndex, vector<CA>& cas, RandomForest& rfNeo, DebugSpec dspec)
{
	vector<vector<float>> samples;
	vector<int> truthClasses;
	CAList::GetSamples(cas, samples, truthClasses); 
	vector<float> pNeo(cas.size());
	ClassifyCAs(samples, rfNeo, pNeo, dspec);
	for (int i = 0; i < cas.size(); i++)
	{
		if (pNeo[i] < 0)
		{
			ErrorExit("ClassifyCAs: Got a lt 0 result in pNeo. Index %d", i);
			cas[i].pNeo = 0.0f;
		}
		else
		{
			cas[i].pNeo = pNeo[i];
			if (dspec.CheckDebugCA(imageSetIndex, cas[i].index)) logd.debug("ClassifyCAsRFNeo: CA %d gets pNeo %.3f", cas[i].index, pNeo[i]);
		}
	}
}
void ClassifySetIsNeo(vector<CA>& cas, int testImageSetCount, DebugSpec dspec)
{
	CAList::SortOnPNeo(cas);
	int n = cas.size();
	int nNeosToKeep = MAX(1, roundf(N_NEOS_TO_SELECT_PER_IMAGE_SET * testImageSetCount));
	int i = 0;
	for (; i < MIN(nNeosToKeep, n); i++)
	{
		cas[i].isNeo = true;
		if (dspec.info) logd.info("ClassifyCAsRFNeo: %s CA %d is classified NEO on pNeo %.3f", cas[i].imageSetId.c_str(), cas[i].index, cas[i].pNeo);
	}
	for (; i < n; i++)
	{
		if ((dspec.caIndex == cas[i].index) || ((dspec.debugDetNumber >= 0) && (dspec.debugDetNumber == cas[i].corrDetNumber)))
			logd.debug("ClassifyCAsRFNeo: CA %d (%s) is classified NOT NEO on pNeo %.3f", cas[i].index, cas[i].imageSetId.c_str(), cas[i].pNeo);
	}
	CAList::SortOnVelocityDescending(cas);
	int nSpeedNeos = MAX(0, roundf(MAX_N_SPEED_NEOS_TO_SELECT_PER_IMAGE_SET * testImageSetCount));
	for (int i = 0; i < MIN(n, nSpeedNeos); i++)
	{
		float v = cas[i].ComputeVelocityPix();
		if (v > NEO_SPEED_THRESHOLD_PIX)
		{
			cas[i].isNeo = true;
			if (dspec.info) logd.info("ClassifyCAsRFNeo: %s CA %d is classified NEO on speed %.1f pix", cas[i].imageSetId.c_str(), cas[i].index, v);
		}
		else
		{
			break;
		}
	}
	CAList::SortOnPClass1(cas);
}
void AddDupCas(vector<CA>& cas, vector<int>& indicesToDup, DebugSpec dspec)
{
	std::sort(indicesToDup.begin(), indicesToDup.end());
	int n = indicesToDup.size();
	for (int i = n - 1; i >= 0; i--)
	{
		int idx = indicesToDup[i];
		CA dup = cas[idx];
		dup.index = BASE_INDEX_FOR_DUP_CAS + idx;
		if (dup.isNeo)
		{
			dup.pClass1 -= 0.000001f; 
			dup.pClass1 = MAX(0.0f, dup.pClass1);
			dup.isNeo = false;
		}
		if (dspec.caIndex == cas[idx].index)
		{
			logd.debug("AddDupCas: Added dup CA at index %d for CA %d. New CA index is %d", idx, cas[idx].index, dup.index);
		}
		cas.insert(cas.begin() + idx, dup);
	}
	if (cas.size() > MAX_ANSWERS)
	{
		cas.resize(MAX_ANSWERS);
	}
}
void ClassifyCAsRFOverlap(int imageSetIndex, vector<CA>& cas, RandomForest& rfOverlap, DebugSpec dspec)
{
	vector<vector<float>> samples;
	vector<int> truthClasses;
	CAList::GetSamples(cas, samples, truthClasses); 
	vector<float> pOverlap(cas.size());
	ClassifyCAs(samples, rfOverlap, pOverlap, dspec);
	for (int i = 0; i < cas.size(); i++)
	{
		if (pOverlap[i] < 0)
		{
			ErrorExit("ClassifyCAs: Got a lt 0 result in pOverlap. Index %d", i);
			cas[i].pOverlap = 0.0f;
		}
		else
		{
			cas[i].pOverlap = pOverlap[i];
			if (dspec.CheckDebugCA(imageSetIndex, cas[i].index)) logd.debug("ClassifyCAsRFNeo: CA %d gets pOverlap %.3f", cas[i].index, pOverlap[i]);
		}
	}
}
void AddNeoBlockers(vector<CA>& cas, DebugSpec dspec)
{
	int n = cas.size();
	for (int i = n - 1; i >= 0; i--)
	{
		if (cas[i].isNeo)
		{
			CA blocker = cas[i];
			blocker.isNeo = false;
			blocker.index = BASE_INDEX_FOR_DUP_CAS + n - i;
			cas.insert(cas.begin() + i, blocker);
			if (dspec.caIndex == cas[i].index)
			{
				logd.debug("AddNeoBlockers: Added blocker CA at index %d for CA %d", i, cas[i].index);
			}
		}
	}
	if (cas.size() > MAX_ANSWERS)
	{
		cas.resize(MAX_ANSWERS);
	}
}
void AddDupsForOverlaps(int nImageSets, vector<CA>& cas, DebugSpec dspec)
{
	CAList::SortOnPOverlap(cas);
	int nToDup = 0; 
	for (int i = 0; i < cas.size(); i++)
	{
		if (cas[i].pOverlap > P_OVERLAP_THRESHOLD) nToDup++;
	}
	nToDup = MIN(MAX_OVERLAP_DUPS_TO_ADD, nToDup);
	if (dspec.debug) logd.debug("AddDupsForOverlaps: Adding %d dups to get to pOverlap 0.5", nToDup);
	vector<int> sortedOverlapIndices(nToDup);
	for (int i = 0; i < nToDup; i++) sortedOverlapIndices[i] = i;
	if (dspec.debug)
	{
		for (int i = 0; i < sortedOverlapIndices.size(); i++)
		{
			int idx = sortedOverlapIndices[i];
			if (dspec.debug) logd.printf("    Adding dup for CA %d which has pOverlap %.3f\n", cas[idx].index, cas[idx].pOverlap);
		}
	}
	AddDupCas(cas, sortedOverlapIndices, dspec); 
	if (dspec.info) logd.info("AddDupsForOverlaps: Added %d dups on pOverlap threshold %.3f, now there are %d CA's", nToDup, P_OVERLAP_THRESHOLD, cas.size());
}
ImageSetResult ProcessImageSet(ImageSet& imageSet, vector<CA>& cas, RandomForest* prefilterRf, int maxCasToKeep, DebugSpec dspec)
{
	if (dspec.debug) logd.printf("------------------------------------------------\n");
	if (dspec.info) logd.printf("ProcessImageSet: Starting set %s\n", imageSet.imageSetId.c_str());
	ImageSetResult imageSetResult;
	imageSetResult.imageSetId = imageSet.imageSetId;
	imageSetResult.imageSetIndex = imageSet.index;
	bool debugThisSet = dspec.CheckDebugImageSet(imageSet.index);
	bool shortcutRegister = false;
	float startSeconds = getTimeSeconds();
	if (imageSet.images.size() != 4)
	{
		ErrorExit("Didn't get 4 images.");
	}
	for (int i = 0; i < 4; i++)
	{
		if ((imageSet.images[i].getWidth() <= 0) || (imageSet.images[i].getHeight() <= 0))
		{
			logd.error("ProcessImageSet: Zero size input image.");
			return imageSetResult;
		}
		if ((dspec.debugImageSetIndex == imageSet.index) && (dspec.debugFrameIndex < 0))
		{
			SaveDebugTiff("loadedOrig", imageSet.images[i].origImage);
		}
	}
	if (DO_CROP_INPUT_IMAGE)
	{
		imageSet.CropImages(); 
		imageSet.AdjustDetsForCrop();
	}
	if (DO_IMAGE_AFFINE_ALIGN)
	{
		if (dspec.debug) logd.printf("------------------------\n");
		if (dspec.debug) logd.printf("Affine align...\n");
		imageSet.AlignOrigImagesAffine(dspec);
		if (dspec.CheckDebugImageSet(imageSet.index))
		{
			for (int i = 0; i < 4; i++)
			{
				SaveDebugTiff("affineOrig", imageSet.images[i].origImage);
			}
		}
		imageSet.UpdateDetsForRegisterShift();
	}
	if (dspec.debugDetNumber >= 0)
	{
		int idx = imageSet.GetDetIndexByDetNumber(imageSet.origDetRecords, dspec.debugDetNumber, 0);
		if (idx >= 0) logd.printf("ProcessImageSet: Original debug det: %s\n", imageSet.origDetRecords[idx].ToString().c_str());
	}
	if (dspec.debug) logd.printf("------------------------\n");
	if (dspec.debug) logd.printf("Creating nsi images...\n");
	for (int i = 0; i < 4; i++)
	{
		if (dspec.debug) logd.printf("Process image %d\n", i);
		bool debug = dspec.CheckDebugImageFrame(imageSet.index, i);
		CreateNsiImageTop(&imageSet.images[i], debug);
	}
	if (dspec.debug) logd.printf("Creating dip mask...\n");
	EnactDipMasks(imageSet, false);
	if (!DO_IMAGE_AFFINE_ALIGN)
	{
		if (dspec.debug) logd.printf("------------------------\n");
		if (dspec.debug) logd.printf("Registration...\n");
		RegisterImagesTop(imageSet, shortcutRegister, dspec);
		if (dspec.debug)
		{
			for (int i = 0; i < 4; i++)
			{
				logd.printf("    Register: image %d: %.1f, %.1f\n", i, imageSet.registerShifts[i].first, imageSet.registerShifts[i].second);
			}
		}
	}
	if (dspec.debug) logd.printf("------------------------\n");
	if (dspec.debug) logd.printf("Diff NSI images...\n");
	CreateNsiDiffImages(imageSet, dspec);
	if (dspec.debug) logd.printf("------------------------\n");
	if (dspec.debug) logd.printf("Find blobs...\n");
	FindBlobsTop(imageSet, dspec);
	if (dspec.debug) logd.printf("Blob counts per frame: %d, %d, %d, %d\n", 
		imageSet.blobSets[0].size(), imageSet.blobSets[1].size(), imageSet.blobSets[2].size(), imageSet.blobSets[3].size());
	if (dspec.debug) logd.printf("------------------------\n");
	if (dspec.debug) logd.printf("Find candidates...\n");
	CAList cal;
	FindCandidatesSimple(imageSet, cal, dspec);
	imageSet.ComputeCARaDecs(cal, dspec); 
	if (dspec.debug) logd.printf("------------------------\n");
	if (dspec.debug) logd.printf("Post-process candidates...\n");
	PostProcessCandidates(imageSet, cal, dspec);
	if (imageSet.origDetRecords.size() > 0)
	{
		if (dspec.debug) logd.printf("------------------------\n");
		if (dspec.debug) logd.printf("Compare to truth...\n");
		bool useNonOverlapDets = true;
		float sensitivity = CompareToTruth(cal, imageSet, useNonOverlapDets, dspec);
		if (dspec.debug) logd.printf("Sensitivity: %.3f\n", sensitivity);
	}
	else
	{
		if (dspec.debug) logd.printf("No dets to compare to truth.\n");
		for (int i = 0; i < cal.size(); i++)
		{
			cal.cas[i].SetDetectionNumber(-1);
		}
	}
	if (dspec.debug) logd.debug("Found %d CA's:", cal.cas.size());
	bool doPrefilter = (DO_PREFILTER_CLASSIFIER_IN_TEST && (prefilterRf != NULL));
	if (doPrefilter)
	{
		int nBeforeDrop = cal.size();
		cal.DropWeakest(MAX_CAS_PER_IMAGE_SET_TO_USE_RF_TO_TOP_N, false, dspec);
		float st = getTimeSeconds();
		ClassifyCAs(cal.cas, *prefilterRf, dspec);
		if (dspec.debug) logd.debug("Done prefilter RF on %d CA's (down from %d) in %.0f seconds", cal.cas.size(), nBeforeDrop, getTimeSeconds() - st);
	}
	bool usePClass1 = doPrefilter;
	int nBeforeDw = cal.size();
	cal.DropWeakest(maxCasToKeep, usePClass1, dspec); 
	if (dspec.info) logd.debug("Drop from %d to %d CA's on %s", nBeforeDw, cal.cas.size(), (usePClass1 ? "pClass1" : "simpleScore"));
	cas.insert(cas.end(), cal.cas.begin(), cal.cas.end());
	float durationSec = getTimeSeconds() - startSeconds;
	if (dspec.debug) logd.debug("Duration for image set %s: %.1f seconds", imageSet.imageSetId.c_str(), durationSec);
	return imageSetResult;
}
bool CheckImageDims(int w, int h)
{
	if (((w == 4110) && (h == 4096))
		|| ((w == 4096) && (h == 4096))
		|| ((w == h) && (w >= 2048) && ((w % 16) == 0)) 
		)
	{
		return true;
	}
	else
	{
		logd.error("Unsupported image dims: %dx%d", w, h);
		return false;
	}
}
bool CheckHeader(vector<string>& header)
{
	if (header.size() < 10)
	{
		logd.error("CheckHeader: Header has too few lines to be a valid fits header.");
		return false;
	}
	return true;
}
bool LoadAndCheckImageSet(int& width, int& height,
	vector<int>& imageData_1, vector<string>& header_1, vector<double>& wcs_1,
	vector<int>& imageData_2, vector<string>& header_2, vector<double>& wcs_2,
	vector<int>& imageData_3, vector<string>& header_3, vector<double>& wcs_3,
	vector<int>& imageData_4, vector<string>& header_4, vector<double>& wcs_4,
	ImageSet& imageSet)
{
	if (!CheckImageDims(width, height))
	{
		return false;
	}
	ImgUtil::SetFromInt(width, height, imageData_1, imageSet.images[0].origImage);
	ImgUtil::SetFromInt(width, height, imageData_2, imageSet.images[1].origImage);
	ImgUtil::SetFromInt(width, height, imageData_3, imageSet.images[2].origImage);
	ImgUtil::SetFromInt(width, height, imageData_4, imageSet.images[3].origImage);
	for (int i = 0; i < 4; i++)
	{
		imageSet.UpdateDimsFromImages(); 
		if ((imageSet.images[i].getWidth() <= 0) || (imageSet.images[i].getHeight() <= 0))
		{
			logd.error("An image was not loaded. Bailing on this image set.");
			return false;
		}
	}
	imageSet.images[0].FinalizeLoad(header_1, wcs_1);
	imageSet.images[1].FinalizeLoad(header_2, wcs_2);
	imageSet.images[2].FinalizeLoad(header_3, wcs_3);
	imageSet.images[3].FinalizeLoad(header_4, wcs_4);
	if (!CheckHeader(header_1)) return 0;
	if (!CheckHeader(header_2)) return 0;
	if (!CheckHeader(header_3)) return 0;
	if (!CheckHeader(header_4)) return 0;
	imageSet.LoadWcsHeaders();
	if (!imageSet.CheckCoordsConversion()) return false;
	return true;
}
class AsteroidDetector
{
private:
	CAList trainCAs;
	CAList testCAs;
	vector<Detection> trainDets;
	unordered_map<string, bool> seenImageSetIds; 
	RandomForest rf;
	RandomForest rfNeo;
	RandomForest rfOverlap;
	float startTrainingAccumTimeSeconds;
	float startTestingTimeSeconds;
	int testImageSetIndex;
	int trainingSetIndex;
	DebugSpec debugSpec;
public:
	AsteroidDetector()
	{
		testImageSetIndex = 0;
		trainingSetIndex = 0;
		debugSpec.debug = true;
		debugSpec.info = true;
		debugSpec.debugBlobCoords = pair<int, int>(-100, -100);
		startTestingTimeSeconds = 0;
		startTrainingAccumTimeSeconds = 0;
#ifdef PRODUCTION
		debugSpec.debug = false;
#endif
	}
	void SetDebugSpec(DebugSpec spec)
	{
		debugSpec = spec;
	}
public:
	void ResetCAIndices(vector<CA>& cas)
	{
		float degTol = 0.00001f;
		for (int i = 0; i < cas.size(); i++)
		{
			cas[i].index = -1;
			for (int j = 0; j < testCAs.size(); j++)
			{
				for (int frameIndex = 0; frameIndex < 4; frameIndex++)
				{
					if (CLOSE_ENOUGH(cas[i].radecs[frameIndex].RaDegrees, testCAs.cas[j].radecs[frameIndex].RaDegrees, degTol)
						&& CLOSE_ENOUGH(cas[i].radecs[frameIndex].DecDegrees, testCAs.cas[j].radecs[frameIndex].DecDegrees, degTol))
					{
						cas[i].index = testCAs.cas[j].index;
						break;
					}
				}
			}
			if (cas[i].index < 0)
			{
				ErrorExit("Failed to find match for CA new index %d: %s", i, cas[i].ToString().c_str());
			}
		}
	}
	int trainingData(int& width, int& height,
		vector<int>& imageData_1, vector<string>& header_1, vector<double>& wcs_1,
		vector<int>& imageData_2, vector<string>& header_2, vector<double>& wcs_2,
		vector<int>& imageData_3, vector<string>& header_3, vector<double>& wcs_3,
		vector<int>& imageData_4, vector<string>& header_4, vector<double>& wcs_4,
		vector<string>& detections)
	{
		if (startTrainingAccumTimeSeconds == 0)
		{
			startTrainingAccumTimeSeconds = getTimeSeconds();
		}
		if (debugSpec.casFeatureCsv.length() > 0)
		{
			CA::LoadFeatureCsv(debugSpec.casFeatureCsv, trainCAs.cas);
			if (trainCAs.size() > 0)
			{
				logd.info("Loaded training csv, shortcutting trainingData.");
				trainingSetIndex++;
				return 1;
			}
			else
			{
				ErrorExit("Failed to load shortcut training csv.");
			}
		}
		float minSoFar = (getTimeSeconds() - startTrainingAccumTimeSeconds) / 60;
		if (minSoFar >= TIME_LIMIT_MIN_TRAIN_ACCUM)
		{
			if (debugSpec.info) logd.info("trainingData: Set index %d: Timeout at %.1f minutes with timeout %.1f", trainingSetIndex, minSoFar, TIME_LIMIT_MIN_TRAIN_ACCUM);
			trainingSetIndex++;
			return 1;
		}
		else
		{
			if (debugSpec.info) logd.info("trainingData: Set index %d: Training accum time so far: %.1f minutes (with timeout %.1f)", trainingSetIndex, minSoFar, TIME_LIMIT_MIN_TRAIN_ACCUM);
		}
		if (EXCLUDE_ZERO_DET_SETS_FROM_TRAIN && (detections.size() < MIN_DETECTIONS_FOR_TRAINING * 4))
		{
			if (debugSpec.info) logd.info("trainingData: Skipping low-detection set %d (has %d dets)", trainingSetIndex, detections.size());
			trainingSetIndex++;
			return 0;
		}
		ImageSet imageSet;
		imageSet.index = trainingSetIndex;
		char tmp[64];
		sprintf(tmp, "IS_%03d", trainingSetIndex);
		imageSet.imageSetId = tmp;
		if (LoadAndCheckImageSet(width, height,
			imageData_1, header_1, wcs_1,
			imageData_2, header_2, wcs_2,
			imageData_3, header_3, wcs_3,
			imageData_4, header_4, wcs_4,
			imageSet))
		{
			imageSet.detStrings = detections;
			imageSet.ParseDets();
			vector<CA> cas;
			ProcessImageSet(imageSet, cas, NULL, MAX_CAS_PER_IMAGE_SET_TRAIN, debugSpec);
			if (debugSpec.info) logd.info("trainingData: Done image set index: %d, ID: %s, found %d CA's", trainingSetIndex, imageSet.imageSetId.c_str(), cas.size());
			trainCAs.Append(cas);
			if (DO_NEO_CLASSIFIER)
			{
				trainDets.insert(trainDets.end(), imageSet.detectionsWithoutOverlaps.begin(), imageSet.detectionsWithoutOverlaps.end());
			}
		}
		else
		{
		}
		trainingSetIndex++;
		return 0; 
	}
	int testingData(string& imageID, int& width, int& height,
		vector<int>& imageData_1, vector<string>& header_1, vector<double>& wcs_1,
		vector<int>& imageData_2, vector<string>& header_2, vector<double>& wcs_2,
		vector<int>& imageData_3, vector<string>& header_3, vector<double>& wcs_3,
		vector<int>& imageData_4, vector<string>& header_4, vector<double>& wcs_4)
	{
		if (debugSpec.info) logd.printf("----------------------------------------------------------------------------------\n");
		StringUtils::TrimRight(imageID);
		if (debugSpec.info) logd.info("testingData: Image set index: %d, ID: %s, image size %d, %d", testImageSetIndex, imageID.c_str(), width, height);
		if (startTestingTimeSeconds == 0)
		{
			startTestingTimeSeconds = getTimeSeconds();
		}
		float minSoFar = (getTimeSeconds() - startTrainingAccumTimeSeconds) / 60;
		float testMinSoFar = (getTimeSeconds() - startTestingTimeSeconds) / 60;
		if (minSoFar >= TIME_LIMIT_MIN_ALL)
		{
			if (debugSpec.info) logd.info("testingData: Timeout at %.1f minutes (test %.1f minutes) with timeout %.1f", minSoFar, testMinSoFar, TIME_LIMIT_MIN_ALL);
			testImageSetIndex++;
			return 1;
		}
		else
		{
			if (debugSpec.info) logd.info("testingData: Time so far: %.1f minutes (test %.1f minutes) with timeout %.1f", minSoFar, testMinSoFar, TIME_LIMIT_MIN_ALL);
		}
		if (testImageSetIndex == 0)
		{
			if (debugSpec.info) logd.info("testingData: First test image set.");
			rf.SetParams(CreateRFParams());
			TrainCAClassifier(trainCAs.cas, rf, debugSpec);
			if (DO_NEO_CLASSIFIER)
			{
				if (debugSpec.info) logd.info("testingData: Training NEO classifier on %d CA's and %d dets.", trainCAs.size(), trainDets.size());
				rfNeo.SetParams(CreateRFParams());
				TrainRFNeoClassifier(trainCAs.cas, trainDets, rfNeo, debugSpec);
			}
			if (DO_OVERLAP_CLASSIFIER)
			{
				if (debugSpec.info) logd.info("testingData: Training Overlap classifier on %d CA's and %d dets.", trainCAs.size(), trainDets.size());
				rfOverlap.SetParams(CreateRFParams());
				TrainRFOverlapClassifier(trainCAs.cas, trainDets, rfOverlap, debugSpec);
			}
		}
		if (seenImageSetIds.find(imageID) != seenImageSetIds.end())
		{
			logd.error("testingData: Already seen image set ID '%s'", imageID.c_str());
			testImageSetIndex++;
			return 0;
		}
		seenImageSetIds[imageID] = true;
		ImageSet imageSet;
		imageSet.index = testImageSetIndex;
		imageSet.imageSetId = imageID;
		if (LoadAndCheckImageSet(width, height,
			imageData_1, header_1, wcs_1,
			imageData_2, header_2, wcs_2,
			imageData_3, header_3, wcs_3,
			imageData_4, header_4, wcs_4,
			imageSet))
		{
			vector<CA> cas;
			ProcessImageSet(imageSet, cas, &rf, MAX_CAS_PER_IMAGE_SET_TEST, debugSpec);
			if (debugSpec.info) logd.info("testingData: Done image set index: %d, ID: %s, kept %d CA's", testImageSetIndex, imageID.c_str(), cas.size());
			if (!DO_PREFILTER_CLASSIFIER_IN_TEST)
			{
				ClassifyCAs(cas, rf, debugSpec);
			}
			if (DO_NEO_CLASSIFIER)
			{
				if (debugSpec.info) logd.info("testingData: Applying NEO classifier.");
				ClassifyCAsRFNeo(imageSet.index, cas, rfNeo, debugSpec);
			}
			if (DO_OVERLAP_CLASSIFIER)
			{
				if (debugSpec.info) logd.info("testingData: Applying Overlap classifier.");
				ClassifyCAsRFOverlap(imageSet.index, cas, rfOverlap, debugSpec);
			}
			testCAs.AccumTopNPClass1(cas, MAX_ANSWERS);
		}
		else
		{
		}
		testImageSetIndex++;
		return 0;
	}
	vector<string> getAnswer()
	{
		if (debugSpec.info) logd.info("getAnswer: Have %d CA's", testCAs.size());
		if (DO_NEO_CLASSIFIER)
		{
			ClassifySetIsNeo(testCAs.cas, testImageSetIndex, debugSpec);
			if (DO_NEO_BLOCKERS)
			{
				AddNeoBlockers(testCAs.cas, debugSpec);
			}
		}
		if (DO_OVERLAP_CLASSIFIER)
		{
			if (debugSpec.info) logd.info("getAnswer: Applying Overlap classifier.");
			AddDupsForOverlaps(testImageSetIndex, testCAs.cas, debugSpec);
		}
		CAList::DoFinalSort(testCAs.cas);
		if (testCAs.cas.size() > MAX_ANSWERS)
		{
			testCAs.cas.resize(MAX_ANSWERS);
		}
		vector<string> answers = BuildAnswers(testCAs, debugSpec);
		for (int i = 0; i < MIN(5, answers.size()); i++)
		{
			if (debugSpec.info) logd.printf("    Answer[%d]: %s\n", i, answers[i].c_str());
		}
		return answers;
	}
};

