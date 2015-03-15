bool CFG_FEATURES_USE_POSITION = false; //use position (RA/DEC)
bool CFG_FEATURES_USE_SPEED = true; //use absolute speed
bool CFG_FEATURES_USE_DIRECTION = true; //use direction of the speed (angle + speed in each direction)

#include <sys/time.h>
#include <pthread.h>
#include "boost/archive/text_oarchive.hpp"
#include "boost/archive/text_iarchive.hpp"
#include "boost/serialization/vector.hpp"
#include <map>
using namespace std;

/** Psyho Code **/

#define INLINE   inline __attribute__ ((always_inline))
#define NOINLINE __attribute__ ((noinline))

#define FOR(i,a,b)  for(int i=(a);i<(b);++i)
#define REP(i,a)    FOR(i,0,a)
#define ZERO(m)     memset(m,0,sizeof(m))
#define ALL(x)      x.begin(),x.end()
#define PB          push_back
#define SZ          size()
#define LL          long long
#define ULL         unsigned long long
#define LD          long double
#define MP          make_pair
#define X           first
#define Y           second
#define VC          vector
#define PII         pair <int, int>
#define PDD         pair <double, double>
#define PFF         pair <float, float>
#define VI          VC < int >
#define VVI         VC < VI >
#define VVVI        VC < VVI >
#define VVVVI       VC < VVVI >
#define VPII        VC < PII >
#define VD          VC < double >
#define VVD         VC < VD >
#define VVVD        VC < VVD >
#define VF          VC < float >
#define VVF         VC < VF >
#define VVVF        VC < VVF >
#define VS          VC < string >
#define VVS         VC < VS >
#define DB(a)       cerr << #a << ": " << (a) << endl;
#define DBM(a)       cerr << (a) << endl;

static double t0 = 0;
double getTime() {
	timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + 1e-6 * tv.tv_usec;
}

VS splt(string s, char c = ' ') {VS all; int p = 0, np; while (np = s.find(c, p), np >= 0) {if (np != p) all.PB(s.substr(p, np - p)); p = np + 1;} if (p < s.SZ) all.PB(s.substr(p)); return all;}

struct RNG {
    unsigned int MT[624];
    int index;
	
	RNG(int seed = 1) {
		init(seed);
	}
    
    void init(int seed = 1) {
        MT[0] = seed;
        FOR(i, 1, 624) MT[i] = (1812433253UL * (MT[i-1] ^ (MT[i-1] >> 30)) + i);
        index = 0;
    }
    
    void generate() {
        const unsigned int MULT[] = {0, 2567483615UL};
        REP(i, 227) {
            unsigned int y = (MT[i] & 0x8000000UL) + (MT[i+1] & 0x7FFFFFFFUL);
            MT[i] = MT[i+397] ^ (y >> 1);
            MT[i] ^= MULT[y&1];
        }
        FOR(i, 227, 623) {
            unsigned int y = (MT[i] & 0x8000000UL) + (MT[i+1] & 0x7FFFFFFFUL);
            MT[i] = MT[i-227] ^ (y >> 1);
            MT[i] ^= MULT[y&1];
        }
        unsigned int y = (MT[623] & 0x8000000UL) + (MT[0] & 0x7FFFFFFFUL);
        MT[623] = MT[623-227] ^ (y >> 1);
        MT[623] ^= MULT[y&1];
    }
    
    unsigned int rand() {
        if (index == 0) {
            generate();
        }
        
        unsigned int y = MT[index];
        y ^= y >> 11;
        y ^= y << 7  & 2636928640UL;
        y ^= y << 15 & 4022730752UL;
        y ^= y >> 18;
        index = index == 623 ? 0 : index + 1;
        return y;
    }
    
    INLINE int next() {
        return rand();
    }
    
    INLINE int next(int x) {
        return rand() % x;
    }
    
    INLINE int next(int a, int b) {
        return a + (rand() % (b - a));
    }
    
    INLINE double nextDouble() {
        return (rand() + 0.5) * (1.0 / 4294967296.0);
    }
};

static RNG rng;

struct TreeNode {
    friend class boost::serialization::access;
    // Serialize the class
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
    	ar & level;
    	ar & feature;
    	ar & value;
    	ar & result;
    	ar & left;
    	ar & right;
    }

	int level;
	int feature;
	float value;
	float result;
	int left;
	int right;
	
	TreeNode() {
		level = -1;
		feature = -1;
		value = 0;
		result = 0;
		left = -1;
		right = -1;
	}
};

struct RandomForestConfig {
    friend class boost::serialization::access;
    // Serialize the class
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
    	ar & randomFeatures;
    	ar & randomPositions;
    	ar & featuresIgnored;
    	ar & groupFeature;
    	ar & maxLevel;
    	ar & maxNodeSize;
    	ar & threadsNo;
    	ar & bagSize;
    	ar & timeLimit;
    	ar & useBootstrapping;
    	ar & computeImportances;
    	ar & computeOOB;
    }

	VI randomFeatures = {5};
	VI randomPositions = {2};
	int featuresIgnored = 0;
	int groupFeature = -1; //NOT IMPLEMENTED
	int maxLevel = 100;
	int maxNodeSize = 1;
	int threadsNo = 1;
	double bagSize = 1.0;
	double timeLimit = 0;
	bool useBootstrapping = true;
	bool computeImportances = false; 
	bool computeOOB = false; //NOT IMPLEMENTED
};

class DecisionTree {
private:
    friend class boost::serialization::access;
    // Serialize the class
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & nodes;
        ar & importances;
    }

	public:
	VC < TreeNode > nodes;
	VD importances;
	
	DecisionTree() { }
	
	template < class T > DecisionTree(VC < VC < T > > &features, VC < T > &results, RandomForestConfig &config, int seed) {
		RNG r(seed);
		
		if (config.computeImportances) {
			importances = VD(features[0].SZ);
		}
	
		VI chosenGroups(features.SZ);
		REP(i, (int)(features.SZ * config.bagSize)) chosenGroups[r.next(features.SZ)]++;
		
		int bagSize = 0;
		REP(i, features.SZ) if (chosenGroups[i]) bagSize++;
		
		VI bag(bagSize);
		VI weight(features.SZ);
		
		int pos = 0;
		
		REP(i, features.SZ) {
			weight[i] = config.useBootstrapping ? chosenGroups[i] : min(1, chosenGroups[i]);
			if (chosenGroups[i]) bag[pos++] = i;
		}
		
		TreeNode root;
		root.level = 0;
		root.left = 0;
		root.right = pos;
		nodes.PB(root);
		
		VI stack;
		stack.PB(0);
		
		while (stack.SZ) {
			bool equal = true;
			
			int curNode = stack[stack.SZ - 1];
			stack.pop_back();
			
			int bLeft = nodes[curNode].left;
			int bRight = nodes[curNode].right;
			int bSize = bRight - bLeft;
			
			int totalWeight = 0; 
			T totalSum = 0;
			FOR(i, bLeft, bRight) {
				totalSum += results[bag[i]] * weight[bag[i]];
				totalWeight += weight[bag[i]];
			}
			
			assert(bSize > 0);
			
			FOR(i, bLeft + 1, bRight) if (results[bag[i]] != results[bag[i - 1]]) {
				equal = false;
				break;
			}
			
			if (equal || bSize <= config.maxNodeSize || nodes[curNode].level >= config.maxLevel) {
				if (totalWeight < 10) {
					nodes[curNode].result = totalSum / totalWeight;
				} else {
					T mn = +1e10;
					T mx = -1e10;
					FOR(i, bLeft + 1, bRight) {
						mn = min(mn, results[bag[i]]);
						mx = max(mx, results[bag[i]]);
					}
					nodes[curNode].result = (totalSum - mn - mx) / (totalWeight - 2);
				}
				continue;
			}
			
			int bestFeature = -1;
			int bestLeft = 0;
			int bestRight = 0;
			T bestValue = 0;
			T bestMSE = 1e99;
			
			const int randomFeatures = config.randomFeatures[min(nodes[curNode].level, (int)config.randomFeatures.SZ - 1)];
			REP(i, randomFeatures) {
			
				int featureID = config.featuresIgnored + r.next(features[0].SZ - config.featuresIgnored);
				
				T vlo, vhi;
				vlo = vhi = features[bag[bLeft]][featureID];
				FOR(j, bLeft + 1, bRight) {
					vlo = min(vlo, features[bag[j]][featureID]);
					vhi = max(vhi, features[bag[j]][featureID]);
				}
				if (vlo == vhi) continue;
				
				const int randomPositions = config.randomPositions[min(nodes[curNode].level, (int)config.randomPositions.SZ - 1)];
				REP(j, randomPositions) {
					T splitValue = features[bag[bLeft + r.next(bSize)]][featureID];
					if (splitValue == vlo) {
						j--;
						continue;
					}
					
					T sumLeft = 0;
					int totalLeft = 0;
					int weightLeft = 0;
					FOR(k, bLeft, bRight) {
						int p = bag[k];
						if (features[p][featureID] < splitValue) {
							sumLeft += results[p] * weight[p];
							weightLeft += weight[p];
							totalLeft++;
						}
					}
					
					T sumRight = totalSum - sumLeft;
					int weightRight = totalWeight - weightLeft;
					int totalRight = bSize - totalLeft;
					
					if (totalLeft == 0 || totalRight == 0)
						continue;
					
					T meanLeft = sumLeft / weightLeft;
					T meanRight = sumRight / weightRight;
					T mse = 0;
					
					mse = (1 - meanLeft) * (1 - meanLeft) * (1 - meanLeft) * sumLeft + 
						meanLeft * meanLeft * meanLeft * (weightLeft - sumLeft) +
						(1 - meanRight) * (1 - meanRight) * (1 - meanRight) * sumRight +
						meanRight * meanRight * meanRight * (weightRight - sumRight);
					
					if (mse < bestMSE) {
						bestMSE = mse;
						bestValue = splitValue;
						bestFeature = featureID;
						bestLeft = totalLeft;
						bestRight = totalRight;
						if (mse == 0) goto outer;
					}
				}
			}
			outer: 
			
			if (bestLeft == 0 || bestRight == 0) {
				if (totalWeight < 10) {
					nodes[curNode].result = totalSum / totalWeight;
				} else {
					T mn = +1e10;
					T mx = -1e10;
					FOR(i, bLeft + 1, bRight) {
						mn = min(mn, results[bag[i]]);
						mx = max(mx, results[bag[i]]);
					}
					nodes[curNode].result = (totalSum - mn - mx) / (totalWeight - 2);
				}
				continue;
			}
			
			if (config.computeImportances) {
				importances[bestFeature] += 1. * (bRight - bLeft) / features.SZ;
			}
			
			T mean = totalSum / totalWeight;
			
			T nextValue = -1e99;
			FOR(i, bLeft, bRight) if (features[bag[i]][bestFeature] < bestValue) nextValue = max(nextValue, features[bag[i]][bestFeature]);
			
			TreeNode left;
			TreeNode right;
			
			left.level = right.level = nodes[curNode].level + 1;
			nodes[curNode].feature = bestFeature;
			nodes[curNode].value = (bestValue + nextValue) / 2.0;
			if (!(nodes[curNode].value > nextValue)) nodes[curNode].value = bestValue;
			nodes[curNode].left = nodes.SZ;
			nodes[curNode].right = nodes.SZ + 1;
			
			int bMiddle = bRight;
			FOR(i, bLeft, bMiddle) {
				if (features[bag[i]][nodes[curNode].feature] >= nodes[curNode].value) {
					swap(bag[i], bag[--bMiddle]);
					i--;
					continue;
				}
			}
			
			assert(bestLeft == bMiddle - bLeft);
			assert(bestRight == bRight - bMiddle);
			
			left.left = bLeft;
			left.right = bMiddle;
			right.left = bMiddle;
			right.right = bRight;
			
			stack.PB(nodes.SZ);
			stack.PB(nodes.SZ + 1);
			
			nodes.PB(left);
			nodes.PB(right);
			
		}
		
		nodes.shrink_to_fit();
	}
	
	template < class T > double estimate(VC < T > &features) {
		TreeNode *pNode = &nodes[0];
		while (true) {
			if (pNode->feature < 0) return pNode->result;
			pNode = &nodes[features[pNode->feature] < pNode->value ? pNode->left : pNode->right];
		}
	}
	
};

RNG gRNG(1);

class RandomForest {
private:
    friend class boost::serialization::access;
    // Serialize the class
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & trees;
        ar & importances;
        ar & config;
    }

	public:
	
	VC < DecisionTree > trees;
	VD importances;
	RandomForestConfig config;

	void clear() {
		trees.clear();
		trees.shrink_to_fit();
	}
	
	template < class T > void train(VC < VC < T > > &features, VC < T > &results, RandomForestConfig &_config, int treesNo) {
		double startTime = getTime();
		config = _config;
		
		if (config.threadsNo == 1) {
			REP(i, treesNo) {	
				if (config.timeLimit && getTime() - startTime > config.timeLimit) break;
				trees.PB(DecisionTree(features, results, config, gRNG.next()));
			}
		} 		
		if (config.computeImportances) {
			importances = VD(features[0].SZ);
			for (DecisionTree tree : trees)
				REP(i, importances.SZ)
					importances[i] += tree.importances[i];
			double sum = 0;
			REP(i, importances.SZ) sum += importances[i];
			REP(i, importances.SZ) importances[i] /= sum;
		}
	}
	
	template <class T> double estimate(VC<T> &features) {
		assert(trees.SZ);
	
		double sum = 0;
		REP(i, trees.SZ) sum += trees[i].estimate(features);
		return sum / trees.SZ;
	}
	
	template <class T> VD estimate(VC<VC<T>> &features) {
		assert(trees.SZ);
	
		int samplesNo = features.SZ;
	
		VD rv(samplesNo);
		if (config.threadsNo == 1) {
			REP(j, samplesNo) {
				REP(i, trees.SZ) rv[j] += trees[i].estimate(features[j]);
				rv[j] /= trees.SZ;
			}
		} 		
        return rv;
	}
	
	LL countTotalNodes() {
		LL rv = 0;
		REP(i, trees.SZ) rv += trees[i].nodes.SZ;
		return rv;
	}
	
};



/** Pfr Code **/

class FitsHeader {
protected:
	map<string, double> _double_value;
public:
	vector<string> _raw_header;
	double double_value(std::string key) const {
		return _double_value.at(key);
	}
	void add_line(std::string raw) {
		_raw_header.push_back(raw);
		istringstream iss(raw);
		string untrimmed_key;
		std::getline(iss, untrimmed_key, '=');
		istringstream key_iss(untrimmed_key);
		string key;
		double value;
		key_iss >> key;
		iss >> value;
		if (!iss.fail())
			_double_value[key] = value;
	}
};

namespace catalina {
	const int ADC_DISCONTINUITY_X = 2055;
	const int ADC_DISCONTINUITY_RECOVER_WIDTH = 3;
};

const double PI = 3.1415926535897932384626433832795;
namespace config {
	const int FRAMES = 4;
	const int MAX_ALLOWED_DETECTIONS = 100000; 
	const int FrameWidth = 4110;
	const int FrameHeight = 4096;
	const int W = FrameWidth;
	const int H = FrameHeight;
	const int MAX_FRAME_DISPLACEMENT = 64;
	const int Margin = 128;
	const int MarginWidth = FrameWidth + 2*Margin;
	const int MarginHeight = FrameHeight + 2*Margin;
	const int NotEigenStride = MarginWidth;
	const int UNGLOW_BLOCK_SIZE = 128;
	const int UNGLOW_MIN_EDGE = UNGLOW_BLOCK_SIZE;
	const int blur_n = 5;
	const float blur_radius = 0.95;
	const float STANDARD_QUAL_FACTOR = 50;
};

namespace astro {
	const double DEG_TO_RAD = PI/180;
	const double RAD_TO_DEG = 180/PI;
};

using namespace config;
template<class T>
struct AnyVec2 {
	T value[2];
	AnyVec2() {}
	AnyVec2(T x, T y) {
		value[0] = x;
		value[1] = y;
	}
	T& operator[](int i) { return value[i]; }
	T operator[](int i) const { return value[i]; }
};

template<class S, class T>
S& operator<<(S& stream, AnyVec2<T> v) {
	return stream << v[0] << " " << v[1];
}

typedef AnyVec2<double> Vec2;
typedef AnyVec2<int> Vec2i;

struct Mat2d {
	typedef double T;
	T a,b,c,d;
	
	Mat2d() {}
	
	Mat2d(T a, T b, T c, T d): a(a), b(b), c(c), d(d) {}
	
	static Mat2d Zero() {
		return Mat2d(0,0,0,0);
	}
	
	T operator()(int i, int j) const {
		if (i==0)
			return j==0 ? a : b;
		else
			return j==0 ? c : d;
	}
	
	T at(int i, int j) const {
		return operator()(i,j);
	}
	
	Vec2 row(int i) {
		return Vec2(at(i,0), at(i,1));
	}
	
	Vec2 col(int j) {
		return Vec2(at(0,j), at(1,j));
	}
	
	Mat2d inverse() const {
		T det = a*d - b*c;
		return Mat2d(d/det, -b/det,	-c/det, a/det);
	}
	
	Mat2d transpose() const {
		return Mat2d(a,c,b,d);
	}
	
	Vec2 operator*(Vec2 x) const {
		return Vec2(a*x[0] + b*x[1], c*x[0] + d*x[1]);
	}
	
	Mat2d operator*(Mat2d m) const {
		Vec2 col0 = (*this) * m.col(0);
		Vec2 col1 = (*this) * m.col(1);
		return Mat2d(col0[0], col1[0], col0[1], col1[1]);
	}
	
	Mat2d& operator+=(Mat2d y) {
		a += y.a;
		b += y.b;
		c += y.c;
		d += y.d;
		return *this;
	}
};

template<class S>
S& operator<<(S& stream, Mat2d m) {
	return stream << m.row(0) << "\n" << m.row(1);
}

Vec2 operator+(Vec2 a, Vec2 b) { return Vec2(a[0]+b[0], a[1]+b[1]); }
Vec2 operator-(Vec2 a, Vec2 b) { return Vec2(a[0]-b[0], a[1]-b[1]); }
Vec2 operator*(Vec2 a, Vec2 b) { return Vec2(a[0]*b[0], a[1]*b[1]); }
Vec2 operator/(Vec2 a, Vec2 b) { return Vec2(a[0]/b[0], a[1]/b[1]); }
Vec2 operator*(double b, Vec2 a) { return Vec2(a[0]*b, a[1]*b); }
Vec2 operator*(Vec2 a, double b) { return Vec2(a[0]*b, a[1]*b); }
Vec2 operator/(Vec2 a, double b) { return Vec2(a[0]/b, a[1]/b); }
Vec2 operator-(Vec2 a) { return Vec2(-a[0], -a[1]); }
typedef Vec2 Point;
typedef Vec2 DarkPoint;
Vec2 operator+(Vec2 a, float b) { return Vec2(a[0]+b, a[1]+b); }
Vec2 operator-(Vec2 a, float b) { return Vec2(a[0]-b, a[1]-b); }
Vec2i operator+(Vec2i a, int b) { return Vec2i(a[0]+b, a[1]+b); }
Vec2i operator-(Vec2i a, int b) { return Vec2i(a[0]-b, a[1]-b); }
float norm2(Vec2 a) { return a[0]*a[0] + a[1]*a[1]; }

float sqr(float x) { return x*x; }
double sqr(double x) { return x*x; }

static inline void assert_invariant(bool condition) {
	if (!condition)
		throw std::logic_error("assert_invariant failed");
}

static inline void range_check(bool condition) {
	if (!condition)
		throw std::out_of_range("range_check failed");
}

struct PlaneToSphereMapping {
	Vec2 ref_pix;
	Vec2 ref_val;
	Mat2d linear_pix_to_val;
	Vec2 operator()(Vec2 yx) const {
		Vec2 yx_pix = yx - ref_pix;
		Vec2 xy_val = -astro::DEG_TO_RAD * (linear_pix_to_val * yx_pix);
		double x = xy_val[0], y = xy_val[1];
		double a0 = astro::DEG_TO_RAD * ref_val[0];
		double d0 = astro::DEG_TO_RAD * ref_val[1];
		double D = atan(sqrt(x*x+y*y));
		double B = atan2(-x, y);
		double XX = sin(d0)*sin(D)*cos(B) + cos(d0)*cos(D);
		double YY = sin(D)*sin(B);
		double a = a0 + atan2(YY, XX);
		double d = asin(sin(d0)*cos(D) - cos(d0)*sin(D)*cos(B));
		return astro::RAD_TO_DEG * Vec2(a,d);
	}
};

struct SphereToPlaneMapping {
	Vec2 ref_val;
	Vec2 ref_pix;
	Mat2d linear_val_to_pix;
	Vec2 operator()(Vec2 rd) const {
		double a0 = astro::DEG_TO_RAD * ref_val[0];
		double d0 = astro::DEG_TO_RAD * ref_val[1];
		double a = astro::DEG_TO_RAD * rd[0];
		double d = astro::DEG_TO_RAD * rd[1];
		double A = cos(d) * cos(a-a0);
		double F = astro::RAD_TO_DEG / (sin(d0)*sin(d) + cos(d0)*A);
		double x_val = F * cos(d) * sin(a-a0);
		double y_val = F * (cos(d0)*sin(d) - A * sin(d0));
		Vec2 xy_val = { x_val, y_val };
		Vec2 yx_pix = (linear_val_to_pix * xy_val);
		return yx_pix + ref_pix;
	}
};

struct FullImage;
struct AffineMap {
	Mat2d linear_matrix;
	Vec2 offset;
	Vec2 operator()(Vec2 yx) const {
		return linear_matrix * yx + offset;
	}
	AffineMap inverse() const {
		Mat2d inv = linear_matrix.inverse();
		return AffineMap { inv, inv * -offset };
	}
	void transform(FullImage src, FullImage dst) const;
};

Mat2d outer(Vec2 x, Vec2 y) {
	Mat2d m(x[0]*y[0], x[0]*y[1], x[1]*y[0], x[1]*y[1]);
	return m;
}

static inline AffineMap operator*(const SphereToPlaneMapping& g, const PlaneToSphereMapping& f) {
	Vec2 center = { 2055, 2048 };
	Vec2 im_center = g(f(center));
	Mat2d left = Mat2d::Zero();
	Mat2d right = Mat2d::Zero();
	const int w = 4;
	for (int i=-w; i<=w; i++) {
		for (int j=-w; j<=w; j++) {
			Vec2 a = { i/(double)w * center[0], j/(double)w * center[1] };
			Vec2 b = g(f(center + a)) - im_center;
			left += outer(a, a);
			right += outer(a, b);
		}
	}
	Mat2d linear_matrix = (left.inverse() * right).transpose();
	Vec2 offset = im_center - linear_matrix * center;
	return AffineMap { linear_matrix, offset };
}

struct FullImageBuffer {
	float raw_data[MarginHeight * MarginWidth];
};

struct FullImage {
	float* origin_ptr;
	shared_ptr<FullImageBuffer> buf;
	
	FullImage(): buf(make_shared<FullImageBuffer>()) {
		origin_ptr = &buf->raw_data[Margin*NotEigenStride + Margin];
	}
	
	float& unsafe_at(int y, int x) const {
		return origin_ptr[y*NotEigenStride + x];
	}
	
	float& at(int y, int x) const {
		range_check(-Margin <= y && y < H+Margin);
		range_check(-Margin <= x && x < W+Margin);
		return unsafe_at(y, x);
	}
	
	float atf(float y, float x) {
		int iy = (int)y;
		int ix = (int)x;
		float wy = y - iy;
		float wx = x - ix;
		return at(iy + 1, ix + 1) * wy * wx + at(iy, ix + 1) * (1 - wy) * wx + at(iy + 1, ix) * wy * (1 - wx) + at(iy, ix) * (1 - wy) * (1 - wx);
	}
	
	float& at(Vec2i yx) const {
		return at(yx[0], yx[1]);
	}
	
	template<class T>
	void fill_with(const T& callable) {
		for (int i=0; i<H; i++)
			for (int j=0; j<W; j++)
				unsafe_at(i,j) = callable(i,j);
	}
};

struct FullImageSet {
	FullImage image[FRAMES];
	FullImage& operator[](int f) {
		return image[f];
	}
	
	template<class T>
	void fill_with(const T& callable) {
		for (int f=0; f<FRAMES; f++)
			image[f].fill_with([f,callable](int i, int j) -> float { return callable(f,i,j); });
	}
};

struct Series {
	std::string id;
	PlaneToSphereMapping darkToSphere[FRAMES];
	SphereToPlaneMapping sphereToDark[FRAMES];
	AffineMap darkToStar[FRAMES];
	AffineMap starToDark[FRAMES];
	FitsHeader frame_head[FRAMES];
	float frame_bias[FRAMES];
	double time_delta;
	
	void setup() {
		time_delta = (frame_head[FRAMES-1].double_value("MJD") - frame_head[0].double_value("MJD")) / (FRAMES-1);
		for (int f=0; f<FRAMES; f++) {
			const FitsHeader& head = frame_head[f];
			Vec2 ref_pix = Vec2(head.double_value("CRPIX2"), head.double_value("CRPIX1"));
			Vec2 ref_val = Vec2(head.double_value("CRVAL1"), head.double_value("CRVAL2"));
			Mat2d pix_to_val(
					head.double_value("CD1_2"), head.double_value("CD1_1"),
					head.double_value("CD2_2"), head.double_value("CD2_1"));
			darkToSphere[f] = PlaneToSphereMapping { ref_pix, ref_val, pix_to_val };
			sphereToDark[f] = SphereToPlaneMapping { ref_val, ref_pix, pix_to_val.inverse() };
			darkToStar[f] = sphereToDark[0] * darkToSphere[f];
			starToDark[f] = darkToStar[f].inverse();
		}
	}
};

struct Detection {
	Series* series;
	Point yx0;
	Vec2 dyx;
	float score;
	float gain_Kobj = NAN;
	float gain_Vobj = NAN;
	float gain_histogram[FRAMES] = { NAN, NAN, NAN, NAN };
	bool neo_candidate = false;
	bool is_duplicate = false;
	
	Detection(Series* series, Point yx0, Vec2 dyx, float score) : series(series), yx0(yx0), dyx(dyx), score(score) {}
	
	bool operator<(const Detection& o) const {
		return score > o.score;
	}
	
	Point yx(int f) const {
		return yx0 + dyx*f;
	}
};

struct DiskSet {
	vector<Point> points;
	float radius2;
	
	DiskSet(float radius): radius2(radius*radius) {}
	
	bool includes(Point p) const {
		for (const Point& q : points) {
			Vec2 d = p-q;
			if (d[0]*d[0] + d[1]*d[1] <= radius2)
				return true;
		}
		return false;
	}
	
	void add(Point p) {
		points.push_back(p);
	}
};

template<class T>
void resize_if_larger(vector<T>& v, int max_size) {
	if ((int)v.size() > max_size) {
		v.erase(v.begin() + max_size, v.end());
	}
}

struct BrightAnchor {
	Series* series;
	int fr;
	Point yxr;
	float score;
	
	bool operator<(const BrightAnchor& o) const {
		return score > o.score;
	}
};

vector<Detection> detectBright(Series* series, FullImageSet Z, FullImageSet &Z2) {
	for (int f=0; f<FRAMES; f++)
		for (int y=0; y<H; y++)
			for (int x=0; x<W; x++) {
				float z = Z[f].at(y,x);
				z = max(0.f, z);
				Z2[f].at(y,x) = sqr(z) / (2 * sqr(STANDARD_QUAL_FACTOR));
			}
			
	const int anchor_step = 9; //TODO: splitting to 9x9 blocks far from good
	const float anchor_thresh = sqr(1.5);
	const int max_frame_anchors = 3000; //TODO: increase
	const int max_anchors = 4*max_frame_anchors; 
	const int SPEED_RANGE = 30; //TODO: increase max speed range
	const float TRACK_MERIT = .075; //TODO: tweak (decrease -> more detections)
	const int max_tracks = 3000; //TODO: increase 
	const int min_track_distance = 12; //TODO: tweak
	const float MEAN_SPEED_THRESHOLD_SCORE = 10; //it's median, d'oh
	const float delta_sigma = .23;
	const float speed_sigma = 70;
	
	double minimumDistanceFromEdge = 0;
	for (int f=0; f<FRAMES; f++) {
		Point p = series->starToDark[f](Point { H/2, W/2 }) - Point{H/2,W/2};
		minimumDistanceFromEdge = max(minimumDistanceFromEdge, max(abs(p[0]), abs(p[1])));
	}
	
	
	vector<BrightAnchor> anchors;
	for (int f=0; f<FRAMES; f++) {
		const int D = anchor_step;
		for (int bi=0; bi<H/D; bi++) {
			for (int bj=0; bj<W/D; bj++) {
				float best_val = 0;
				Vec2i best_pos = Vec2i{0,0};
				for (int i=0; i<D; i++)
					for (int j=0; j<D; j++) {
						float val = Z2[f].at(bi*D+i, bj*D+j);
						Vec2i pos = Vec2i(bi*D+i, bj*D+j);
						if (val >= best_val) {
							best_val = val;
							best_pos = Vec2i(bi*D+i, bj*D+j);
						}
					}
					
				//TODO: consult if we should remove asteroids close to border??
				if (best_pos[0]<15 || best_pos[1]<38)
					continue;
				if (best_pos[0]>H-15 || best_pos[1]>W-25)
					continue;
					
				
				// if (best_pos[0]<minimumDistanceFromEdge || best_pos[1]<minimumDistanceFromEdge)
					// continue;
				// if (best_pos[0]>H-minimumDistanceFromEdge || best_pos[1]>W-minimumDistanceFromEdge)
					// continue;
					
				if (best_val >= anchor_thresh)
					anchors.push_back(BrightAnchor { series, f, series->darkToStar[f](Vec2(best_pos[0], best_pos[1])), best_val });
			}
		}
	}
	
	sort(anchors.begin(), anchors.end());
	if ((int)anchors.size() > max_anchors) {
		resize_if_larger(anchors, max_anchors);
	}
	
	vector<Detection> allTracks;
	for (BrightAnchor& anchor : anchors) {
		int fr = anchor.fr;
		Point yxr = anchor.yxr;
		Vec2 offset_f[FRAMES];
		for (int f=0; f<FRAMES; f++)
			offset_f[f] = series->starToDark[f](yxr);
		float best_merit = 0;
		Vec2 best_dyx = Vec2(0,0);
		float best_x[FRAMES];
		for (int vi=-SPEED_RANGE; vi<=SPEED_RANGE; vi++) {
			for (int vj=-SPEED_RANGE; vj<=SPEED_RANGE; vj++) {
				Vec2 dyx = { vi/5.0, vj/5.0 }; //TODO: increase speed resolution
				float x[FRAMES];
				for (int f=0; f<FRAMES; f++)
					x[f] = Z2[f].at(
							(int)(offset_f[f][0] + dyx[0] * (f-fr) + .5),
							(int)(offset_f[f][1] + dyx[1] * (f-fr) + .5)
							);
				std::sort(x, x+FRAMES);
				float merit = 4*x[0] + 2*x[1] + x[2]; //TODO: tweak? probably pointless
				if (merit > best_merit) { 
					best_merit = merit;
					best_dyx = dyx;
					for (int f=0; f<FRAMES; f++)
						best_x[f] = x[f];
				}
			}
		}
		
		if (best_merit > TRACK_MERIT) {
			Detection d = Detection { series, yxr - best_dyx * fr, best_dyx, best_merit };
			std::copy(&best_x[0], &best_x[FRAMES], d.gain_histogram);
			allTracks.push_back(d);
		}
	}
	
	//TODO: removes "duplicates" - can be improved if it's going to be the bottleneck, since it's naive O(n^2)
	sort(allTracks.begin(), allTracks.end());
	cerr << (int)allTracks.size() << " raw tracks" << endl;
	vector<Detection> tracks;
	DiskSet diskSet(min_track_distance);
	for (Detection& d : allTracks) {
		if (diskSet.includes(d.yx0))
			continue;
		diskSet.add(d.yx0);
		tracks.push_back(d);
		if ((int)tracks.size() >= max_tracks)
			break;
	}
	
	return tracks;
}


vector<string> finish(vector<Detection>& detections) {
	std::sort(detections.begin(), detections.end());
	resize_if_larger(detections, MAX_ALLOWED_DETECTIONS);
	
	//TODO: remove duplications
	vector<Detection> extra_detections;
	for (Detection& d : detections) {
		Detection d2 = d;
		float s = d.score;
		float M = 1.2e9;
		float s2 = M/sqr(6000*(1 - .95 * max(1e-4f, 1 - 1/sqrt(s/M)/6000) * max(1e-4f, 1 - 1/sqrt(s/M)/4800)));
		s2 = min(s2, s - 20);
		d2.score = s2;
		d2.is_duplicate = true;
		extra_detections.push_back(d2);
	}
	for (Detection& d : extra_detections) {
		detections.push_back(d);
	}
	
	std::sort(detections.begin(), detections.end());
	resize_if_larger(detections, MAX_ALLOWED_DETECTIONS);
	vector<string> lines;
	for (int i=0; i<(int)detections.size(); i++) {
		Detection& d = detections[i];
		stringstream ss;
		ss << d.series->id << ' ';
		for (int f=0; f<FRAMES; f++) {
			Vec2 ra_dec = d.series->darkToSphere[f](d.series->starToDark[f](d.yx(f)));
			ss << setprecision(6) << fixed;
			ss << ra_dec[0] << ' ' << ra_dec[1] << ' ';
		}
		ss << (int)d.neo_candidate;
		lines.push_back(ss.str());
	}
	return lines;
}

template<class T>
static inline T lerp(T t, T a, T b) {
	return a + t * (b-a);
}

template<class T>
static inline T cubic_delta(T t, T bdelta, T b, T c, T cdelta) {
	return lerp(t*t*(3-2*t), b, c) + 0.5*(t - t*t) * lerp(t, bdelta, -cdelta);
}

template<class T>
static inline T cubic(T t, T a, T b, T c, T d) {
	return cubic_delta(t, c-a, b, c, d-b);
}

float arrayinterp1d_cubic_nearest(const float* ptr, int size, int stride, float fi) {
	int i = (int)fi;
	float t = fi - i;
	float y;
	if (i < 1) {
		if (fi < 0) {
			y = ((ptr[stride*(0)]));
		} else {
			y = cubic_delta(t, (2*((ptr[stride*(1)]) - (ptr[stride*(0)]))), (ptr[stride*(0)]), (ptr[stride*(1)]), (size==2 ? (2*((ptr[stride*(size-1)]) - (ptr[stride*(size-2)]))) : (ptr[stride*(2)]) - (ptr[stride*(0)])));
		}
	} else if (i >= size-2) {
		if (i >= size-1) {
			y = ((ptr[stride*(size-1)]));
		} else {
			y = cubic_delta(t, (ptr[stride*(i+1)]) - (ptr[stride*(i-1)]), (ptr[stride*(i)]), (ptr[stride*(i+1)]), (2*((ptr[stride*(size-1)]) - (ptr[stride*(size-2)]))));
		}
	} else {
		y = cubic(t, (ptr[stride*(i-1)]), (ptr[stride*(i)]), (ptr[stride*(i+1)]), (ptr[stride*(i+2)]));
	}
	return y;
}

float arrayinterp1d_cubic_nan(const float* ptr, int size, int stride, float fi) {
	int i = (int)fi;
	float t = fi - i;
	float y;
	if (i < 1) {
		if (fi < 0) {
			y = NAN;
		} else {
			y = cubic_delta(t, (2*((ptr[stride*(1)]) - (ptr[stride*(0)]))), (ptr[stride*(0)]), (ptr[stride*(1)]), (size==2 ? (2*((ptr[stride*(size-1)]) - (ptr[stride*(size-2)]))) : (ptr[stride*(2)]) - (ptr[stride*(0)])));
		}
	} else if (i >= size-2) {
		if (i >= size-1) {
			y = NAN;
		} else {
			y = cubic_delta(t, (ptr[stride*(i+1)]) - (ptr[stride*(i-1)]), (ptr[stride*(i)]), (ptr[stride*(i+1)]), (2*((ptr[stride*(size-1)]) - (ptr[stride*(size-2)]))));
		}
	} else {
		y = cubic(t, (ptr[stride*(i-1)]), (ptr[stride*(i)]), (ptr[stride*(i+1)]), (ptr[stride*(i+2)]));
	}
	return y;
}

static void partial_transform(float* src, float* dst, int nx, int ny, int sx, int sy, float s, float a, float b) {
	float idx[max(W,H)];
	for (int x=0; x<nx; x++)
		idx[x] = x*s + b;
	float offset = 0;
	// cerr << s << "," << a << "," << b << endl;
	for (int y=0; y<ny; y++) {
		for (int x=0; x<nx; x++) {
			dst[y*sy + x*sx] = arrayinterp1d_cubic_nan(&src[y*sy], nx, sx, idx[x] + offset);
		}
		offset += a;
	}
}

void AffineMap::transform(FullImage img, FullImage tmp) const {
	AffineMap inv = this->inverse();
	double xx = inv.linear_matrix(0,0);
	double xy = inv.linear_matrix(0,1);
	double x0 = inv.offset[0];
	double yx = inv.linear_matrix(1,0);
	double yy = inv.linear_matrix(1,1);
	double y0 = inv.offset[1];
	if (sqr(xx-1) + sqr(yy-1) + sqr(xy) + sqr(yx) + sqr(x0) + sqr(y0) < 1e-16)
		return;
	double k = xy / yy;
	partial_transform(&img.at(0,0), &tmp.at(0,0), H, W, NotEigenStride, 1, xx - k*yx, k, x0 - k *y0);
	partial_transform(&tmp.at(0,0), &img.at(0,0), W, H, 1, NotEigenStride, yy, yx, y0);
}

void clean_catalina_discontinuity(FullImage& img) {
	int x_disc = catalina::ADC_DISCONTINUITY_X;
	int w_disc = catalina::ADC_DISCONTINUITY_RECOVER_WIDTH;
	float delta[H];
	for (int i=0; i<H; i++) {
		float s = 0;
		for (int j=x_disc-w_disc; j<x_disc; j++)
			s -= img.at(i,j);
		for (int j=x_disc; j<x_disc+w_disc; j++)
			s += img.at(i,j);
		delta[i] = s/w_disc;
	}
	std::nth_element(&delta[0], &delta[H/2], &delta[H]);
	float offset = roundf(delta[H/2]);
	// cerr << "ADC offset: " << offset << endl;
	for (int i=0; i<H; i++)
		for (int j=x_disc; j<W; j++)
			img.at(i,j) -= offset;
}

void computeUnglow(Series* series, FullImageSet frames) {
	for (int f=0; f<FRAMES; f++) {
		FullImage& data = frames[f];
		clean_catalina_discontinuity(data);
		const int S = config::UNGLOW_BLOCK_SIZE;
		const int E = config::UNGLOW_MIN_EDGE;
		const int w = 1;
		// Point center = series->starToDark[f](Point { H/2, W/2 }); //Source of the bug, the idea to align the blocks between the frames doesn't make much sense anyway
		Point center = Point{H/2,W/2}; 
		
		vector<int> Y, X;
		int i_mid = 1 + (H/2 - E)/S;
		int ni_full = 2*i_mid + 1;
		int ni = ni_full-w;
		int j_mid = 1 + (W/2 - E)/S;
		int nj_full = 2*j_mid + 1;
		int nj = nj_full-w;
		for (int i=0; i<ni_full; i++)
			Y.push_back(i==0 ? 0 : i<ni_full-1 ? center[0] + S*(i-i_mid) : H);
		for (int j=0; j<nj_full; j++)
			X.push_back(j==0 ? 0 : j<nj_full-1 ? center[1] + S*(j-j_mid) : W);
		
		// cerr << "about to median" << endl;
		vector<float> glow(ni*nj);
		const int M = UNGLOW_MIN_EDGE*2 + UNGLOW_BLOCK_SIZE;
		float glow_input[M*M];
		for (int i=0; i<ni; i++)
			for (int j=0; j<nj; j++) {
				int height = Y[i+w] - X[i];
				int width = X[j+w] - X[j];
				range_check(0 < height && height <= M && 0 < width && width <= M);
				for (int y=Y[i]; y<Y[i+w]; y++)
					for (int x=X[j]; x<X[j+w]; x++)
						glow_input[(y-Y[i])*width + (x-X[j])] = data.at(y,x);
				int p = width*height / 2;
				std::nth_element(&glow_input[0], &glow_input[p], &glow_input[width*height]);
				glow[i*nj + j] = glow_input[p];
			}
			
		// cerr << "about to expand" << endl;
		vector<float> glow_half_expanded(ni*W);
		for (int i=0; i<ni; i++) {
			for (int x=0; x<W; x++) {
				glow_half_expanded[i*W + x] =
					arrayinterp1d_cubic_nearest(&glow[i*nj + 0], nj, 1, (x-center[1])/S + j_mid - 0.5f * w);
			}
		}
		
		// cerr << "about to expand" << endl;
		for (int y=0; y<H; y++) {
			for (int x=0; x<W; x++) {
				data.at(y,x) -=
					arrayinterp1d_cubic_nearest(&glow_half_expanded[0*W + x], ni, W, (y-center[0])/S + i_mid - 0.5f * w);
			}
		}
		
		std::nth_element(&glow[0], &glow[ni*nj/2], &glow[ni*nj]);
		series->frame_bias[f] = glow[ni*nj/2];
	}
}

void computeLayers(Series* series, FullImageSet unglow, FullImage star, FullImage dark, FullImageSet transformed, FullImageSet mov) {
	const float SATURATION = 15000;
	for (int f=0; f<FRAMES; f++) {
		float frame_bias = series->frame_bias[f];
		// cerr << "frame_bias = " << frame_bias << endl;
		for (int i=0; i<H; i++) {
			for (int j=0; j<W; j++) {
				float val = unglow[f].at(i,j);
				float biased_val = val + frame_bias;
				bool invalid = biased_val<1000 || val>SATURATION || j==1190 || j<7 || j>W-7;
				unglow[f].at(i,j) = invalid ? NAN : val;
			}
		}
	}
	FullImage tmp;
	for (int f=0; f<FRAMES; f++) {
		for (int i=0; i<H; i++) {
			for (int j=0; j<W; j++) {
				transformed[f].at(i,j) = unglow[f].at(i,j);
			}
		}
		series->darkToStar[f].transform(transformed[f], tmp);
	}
	for (int i=0; i<H; i++) {
		for (int j=0; j<W; j++) {
			float buf[FRAMES];
			int p = 0;
			for (int f=0; f<FRAMES; f++) {
				float x = transformed[f].at(i,j);
				if (x==x)
					buf[p++] = x;
			}
			float r;
			if (p==0)
				r = 0;
			else {
				int k = (p-1)/2;
				std::nth_element(&buf[0], &buf[k], &buf[p]);
				if ((p&1)==0) {
					float x = buf[k];
					float y = *std::min_element(&buf[k+1], &buf[p]);
					r = x*y < 0 ? 0 : std::abs(x)<std::abs(y) ? x : y;
				} else {
					r = buf[k];
				}
			}
			star.at(i,j) = r;
		}
	}
	for (int f=0; f<FRAMES; f++) {
		for (int i=0; i<H; i++) {
			for (int j=0; j<W; j++) {
				transformed[f].at(i,j) = star.at(i,j);
			}
		}
		series->starToDark[f].transform(transformed[f], tmp);
	}
	for (int i=0; i<H; i++) {
		for (int j=0; j<W; j++) {
			float buf[FRAMES];
			int p = 0;
			for (int f=0; f<FRAMES; f++) {
				float x = unglow[f].at(i,j) - transformed[f].at(i,j);
				if (x==x)
					buf[p++] = x;
			}
			float r;
			if (p==0)
				r = 0;
			else {
				int k = (p-1)/2;
				std::nth_element(&buf[0], &buf[k], &buf[p]);
				if ((p&1)==0) {
					float x = buf[k];
					float y = *std::min_element(&buf[k+1], &buf[p]);
					r = x*y < 0 ? 0 : std::abs(x)<std::abs(y) ? x : y;
				} else {
					r = buf[k];
				}
			}
			dark.at(i,j) = r;
			for (int f=0; f<FRAMES; f++) {
				mov[f].at(i,j) = unglow[f].at(i,j) - transformed[f].at(i,j) - r;
			}
		}
	}
}

template<int n>
struct SeparableKernel {
	float weight[n];
	float& operator()(int i) { return weight[n/2+i]; }
};

template<int n>
SeparableKernel<n> make_gaussian(float sigma) {
	SeparableKernel<n> r;
	for (int i=-n/2; i<=n/2; i++)
		r(i) = expf(-i*i/(2*sigma*sigma));
	float s = 0;
	for (int i=-n/2; i<=n/2; i++)
		s += r(i);
	for (int i=-n/2; i<=n/2; i++)
		r(i) /= s;
	return r;
}

template<int n>
void convolve(FullImage src, FullImage dst, FullImage tmp, SeparableKernel<n> ker) {
	for (int i=0; i<H; i++) for (int j=0; j<W; j++) {
		float s = 0;
		for (int k=-n/2; k<=n/2; k++)
			s += src.unsafe_at(i,j-k) * ker(k);
		tmp.at(i,j) = s;
	}
	for (int i=0; i<H; i++) for (int j=0; j<W; j++) {
		float s = 0;
		for (int k=-n/2; k<=n/2; k++)
			s += tmp.unsafe_at(i-k,j) * ker(k);
		dst.at(i,j) = s;
	}
}

template<int n>
void box_filter(FullImage img, FullImage tmp) {
	for (int i=0; i<H; i++) for (int j=0; j<W; j++) {
		float s = 0;
		for (int k=-n/2; k<=n/2; k++)
			s += img.unsafe_at(i,j-k);
		tmp.at(i,j) = s;
	}
	for (int i=0; i<H; i++) for (int j=0; j<W; j++) {
		float s = 0;
		for (int k=-n/2; k<=n/2; k++)
			s += tmp.unsafe_at(i-k,j);
		img.at(i,j) = s * (1.0f / (n*n));
	}
}

template<int n>
void max_filter(FullImage img, FullImage tmp) {
	for (int i=0; i<H; i++) for (int j=0; j<W; j++) {
		float s = -9e99;
		for (int k=-n/2; k<=n/2; k++)
			s = max(s, img.unsafe_at(i,j-k));
		tmp.at(i,j) = s;
	}
	for (int i=0; i<H; i++) for (int j=0; j<W; j++) {
		float s = -9e99;
		for (int k=-n/2; k<=n/2; k++)
			s = max(s, tmp.unsafe_at(i-k,j));
		img.at(i,j) = s;
	}
}

SeparableKernel<blur_n> gauss_kernel = make_gaussian<blur_n>(blur_radius);
float fast_invsqrt(float x) {
	union {
		float f;
		uint32_t i;
	} u;
	u.f = x;
	u.i = 0x5f3759df - u.i/2;
	float y = u.f;
	return y * (1.5f - 0.5f * x * y * y);
}

void computeFiltered(Series* series, FullImageSet tstar, FullImage dark, FullImageSet mov, FullImageSet z) {
	FullImage dark_blur;
	FullImage dark_lowpass;
	FullImage w;
	FullImage tmp;
	FullImage frame_lowpass;
	float sigma_base = 29;
	float filtered_sigma_base = sigma_base / (blur_radius * sqrt(4*PI));
	float dark_uncertainty_for_lowpass = .2;
	float dark_uncertainty = .4;
	float star_uncertainty = .055;
	const int B = 9;
	for (int i=0; i<H; i++) for (int j=0; j<W; j++) {
		w.at(i,j) = fast_invsqrt(sqr(sigma_base) + sqr(dark.at(i,j)));
		dark_lowpass.at(i,j) = dark.at(i,j) * w.at(i,j);
	}
	box_filter<B>(dark_lowpass, tmp);
	box_filter<B>(w, tmp);
	for (int i=0; i<H; i++) for (int j=0; j<W; j++) {
		float low = dark_lowpass.at(i,j) / w.at(i,j);
		dark.at(i,j) -= low;
		dark_lowpass.at(i,j) = low;
	}
	convolve<blur_n>(dark, dark_blur, tmp, gauss_kernel);
	for (int f=0; f<FRAMES; f++) {
		for (int i=0; i<H; i++) {
			mov[f].at(i,1190) = 0.5f*(mov[f].at(i,1189) + mov[f].at(i,1191));
		}
		for (int i=0; i<H; i++) for (int j=0; j<W; j++) {
			if (mov[f].at(i,j)!=mov[f].at(i,j)) {
				w.at(i,j) = 0;
				frame_lowpass.at(i,j) = 0;
			} else {
				w.at(i,j) = fast_invsqrt(sqr(sigma_base) + sqr(dark_uncertainty_for_lowpass) * sqr(dark_blur.at(i,j)) + sqr(mov[f].at(i,j)));
				frame_lowpass.at(i,j) = mov[f].at(i,j) * w.at(i,j);
			}
		}
		box_filter<B>(frame_lowpass, tmp);
		box_filter<B>(w, tmp);
		for (int i=0; i<H; i++) for (int j=0; j<W; j++) {
			float low = frame_lowpass.at(i,j) / w.at(i,j);
			mov[f].at(i,j) -= low;
			frame_lowpass.at(i,j) = low;
		}
		convolve<blur_n>(mov[f], mov[f], tmp, gauss_kernel);
		max_filter<3>(tstar[f], tmp);
		convolve<blur_n>(tstar[f], tstar[f], tmp, gauss_kernel);
		for (int i=0; i<H; i++) for (int j=0; j<W; j++) {
			float sigma2 = sqr(filtered_sigma_base) + sqr(star_uncertainty) * sqr(tstar[f].at(i,j)) + sqr(dark_uncertainty) * sqr(dark.at(i,j));
			z[f].at(i,j) = mov[f].at(i,j) * fast_invsqrt(sigma2);
			tmp.at(i,j) = z[f].at(i,j) >= -3.5;
		}
		box_filter<5>(tmp, w);
		for (int i=0; i<H; i++) for (int j=0; j<W; j++) {
			z[f].at(i,j) = tmp.at(i,j)!=1.f ? 0 : z[f].at(i,j);
			z[f].at(i,j) *= STANDARD_QUAL_FACTOR;
		}
	}
}

VF extractFeatures(Detection &detection, Series *series, FullImageSet &Z2) {
	//Reserved Features (not used in actual prediction)
	VF rv;
	rv.PB(0);
	rv.PB(0);
	rv.PB(0);
	
	REP(i, 4) {
		Vec2 ra_dec = series->darkToSphere[i](series->starToDark[i](detection.yx(i)));
		rv.PB(ra_dec[0]);
		rv.PB(ra_dec[1]);
	}
	
	
	//Actual Features
	float v1[4] = {0};
	float v3[4] = {0};
	
	REP(i, 4) {
		Vec2 offset = series->starToDark[i](detection.yx0);
		v1[i] = Z2[i].atf(offset[0] + detection.dyx[0] * i, offset[1] + detection.dyx[1] * i);
		FOR(dx, -1, 2) FOR(dy, -1, 2) 
			v3[i] += Z2[i].atf(offset[0] + detection.dyx[0] * i + dx, offset[1] + detection.dyx[1] * i + dy);
	}
	
	sort(v1, v1 + 4);
	sort(v3, v3 + 4);
	
	rv.PB(v1[0]);
	rv.PB(v1[1]);
	rv.PB(v1[2]);
	rv.PB(v1[3]);
	
	rv.PB(v1[3] - v1[0]);
	rv.PB(v1[3] - v1[2]);
	rv.PB(v1[2] - v1[1]);
	rv.PB(v1[1] - v1[0]);
	rv.PB(v1[0] * 4 + v1[1] * 2 + v1[2]);
	rv.PB(v1[0] * 10 + v1[1] * 5 + v1[2] * 3 + v1[3]);
	
	rv.PB(v3[0]);
	rv.PB(v3[1]);
	rv.PB(v3[2]);
	rv.PB(v3[3]);
	
	rv.PB(v3[3] - v3[0]);
	rv.PB(v3[3] - v3[2]);
	rv.PB(v3[2] - v3[1]);
	rv.PB(v3[1] - v3[0]);
	rv.PB(v3[0] * 4 + v3[1] * 2 + v3[2]);
	rv.PB(v3[0] * 10 + v3[1] * 5 + v3[2] * 3 + v3[3]);
	
	
	if (CFG_FEATURES_USE_POSITION) {
		double sx = 0, sy = 0;
		REP(i, 4) {
			Vec2 ra_dec = series->darkToSphere[i](series->starToDark[i](detection.yx(i)));
			sx += ra_dec[0];
			sy += ra_dec[1];
		}
		rv.PB(sx);
		rv.PB(sy);
	}
	
	if (CFG_FEATURES_USE_SPEED) {
		rv.PB(norm2(detection.dyx));
		rv.PB(norm2(detection.dyx) / sqr(series->time_delta));
	}
	
	if (CFG_FEATURES_USE_DIRECTION) {
		rv.PB(atan2(detection.dyx[0], detection.dyx[1]));
		rv.PB(detection.dyx[0] / sqr(series->time_delta));
		rv.PB(detection.dyx[1] / sqr(series->time_delta));
	}
	
	return rv;
}

VVF extractAllFeatures(Series* series, FullImageSet frames) {
	double startTime = getTime();
	FullImageSet z;
	computeUnglow(series, frames); 
	// cerr << "Time Passed Unglow: " << (getTime() - startTime) << endl;
	FullImage star, dark;
	FullImageSet mov;
	FullImageSet tstar;
	computeLayers(series, frames, star, dark, tstar, mov);
	// cerr << "Time Passed Layers: " << (getTime() - startTime) << endl;
	computeFiltered(series, tstar, dark, mov, z);
	// cerr << "Time Passed Filtered: " << (getTime() - startTime) << endl;
	FullImageSet Z2;
	vector<Detection> detections = detectBright(series, z, Z2);
	// cerr << "Time Passed Bright: " << (getTime() - startTime) << endl;
	
	VVF rv;
	for (auto &detection : detections) rv.PB(extractFeatures(detection, series, Z2));
	return rv;
}

class InputSession {
protected:
	bool isTrainingMode;
	int numberOfSeriesRead;
public:
	const int expectedNumberOfSeries;
	InputSession(int expectedNumberOfSeries): isTrainingMode(true), numberOfSeriesRead(0), expectedNumberOfSeries(expectedNumberOfSeries) {}
	bool finished() {
		return numberOfSeriesRead == expectedNumberOfSeries;
	}
	Series* readSeries(FullImageSet& frame_data) {
		Series* series = new Series(); //TODO: memory leak?
		int height, width;
		int dummy;
		if (!isTrainingMode)
			cin >> series->id;
		else
			series->id = "<invalid>";
		cin >> width >> height;
		if (cin.eof()) {
			cerr << "premature end of file" << endl;
			std::exit(1);
		}
		for (int f=0; f<FRAMES; f++) {
			for (int i=0; i<H; i++) for (int j=0; j<W; j++) {
				int val;
				cin >> val;
				frame_data[f].at(i,j) = val;
			}
			int headerSize;
			cin >> headerSize;
			std::string line;
			std::getline(cin, line);
			for (int i=0; i<headerSize; i++) {
				std::getline(cin, line);
				series->frame_head[f].add_line(line);
			}
			for (int i=0; i<8; i++)
				std::getline(cin, line);
		}
		if (isTrainingMode) {
			int nDetections;
			cin >> nDetections;
			std::string line;
			std::getline(cin, line);
			for (int i=0; i<nDetections; i++)
				std::getline(cin, line);
		}
		if (cin.fail())
			cerr << "WARNING: failure in cin!" << endl;
		if (isTrainingMode) {
			cout << 1 << endl;
			isTrainingMode = false;
			delete series;
			return readSeries(frame_data);
		} else {
			cout << "dummy" << endl;
			numberOfSeriesRead++;
			series->setup();
			return series;
		}
	}
};

static bool initialized = false;
static vector<Detection> detections;

VVF trainFeatures;
VVF testFeatures;
VF trainResults;

struct AsteroidDetector {
	

	AsteroidDetector() {
		t0 = getTime();
	}

	VC<VC<PDD>> convertDetections(VS &v) {
		VC<VC<PDD>> rv;
		VC<PDD> current;
		for (string s : v) {
			VS vs = splt(s, ' ');
			double ra = atof(vs[2].c_str());
			double dec = atof(vs[3].c_str());
			current.PB(MP(ra, dec));
			if (current.SZ == 4) {
				rv.PB(current);
				current.clear();
			}
		}
		return rv;
	}

	double resultsDistance(VC<PDD> &v1, VC<PDD> &v2) {
		if (v1[0].X > v2[0].X + 0.1 || v1[0].X < v2[0].X - 0.1 || v1[0].Y > v2[0].Y + 0.1 || v1[0].Y < v2[0].Y - 0.1) return 1e9;

		double rv = 0;
		REP(i, 4) rv += (v1[i].X - v2[i].X) * (v1[i].X - v2[i].X) + (v1[i].Y - v2[i].Y) * (v1[i].Y - v2[i].Y);
		return rv;
	}

	double resultsDistance(VF &v1, VF &v2) {
		if (v1[3] > v2[3] + 0.1 || v1[3] < v2[3] - 0.1 || v1[4] > v2[4] + 0.1 || v1[4] < v2[4] - 0.1) return 1e9;
		
		double rv = 0;
		FOR(i, 3, 11) rv += (v1[i] - v2[i]) * (v1[i] - v2[i]);
		return rv;
	}

	double resultsDistance(VF &v1, VC<PDD> &v2) {
		if (v1[3] > v2[0].X + 0.1 || v1[3] < v2[0].X - 0.1 || v1[4] > v2[0].Y + 0.1 || v1[4] < v2[0].Y - 0.1) return 1e9;
		
		double rv = 0;
		REP(i, 4) rv += (v1[i*2+3] - v2[i].X) * (v1[i*2+3] - v2[i].X) + (v1[i*2+4] - v2[i].Y) * (v1[i*2+4] - v2[i].Y);
		return rv;
	}


	int trainCall = 0;
	int testCall = 0;

	int trainingData(int width, int height,
			vector<int> data0, vector<string> head0, vector<double> wcs0,
			vector<int> data1, vector<string> head1, vector<double> wcs1,
			vector<int> data2, vector<string> head2, vector<double> wcs2,
			vector<int> data3, vector<string> head3, vector<double> wcs3,
			vector<string> detections) {
			
		DB(trainCall);
		cerr << "AsteroidDetector: time = " << getTime() - t0 << endl;
		FullImageSet frames;
		Series* series = new Series(); //TODO: memory leak?
		series->id = "training";
		assert(width == W);
		assert(height == H);
		vector<int>* pData[FRAMES] = { &data0, &data1, &data2, &data3 };
		vector<string>* pHead[FRAMES] = { &head0, &head1, &head2, &head3 };
		for (int f=0; f<FRAMES; f++) {
			for (int i=0; i<H; i++) for (int j=0; j<W; j++) {
				frames[f].at(i,j) = (*pData[f])[i*W + j];
			}
			for (string line : *pHead[f]) {
				series->frame_head[f].add_line(line);
			}
		}
		series->setup();
		// cerr << "frontend: Got series " << series->id << endl;
		
		const double TRAIN_CORRECT_DISTANCE = 0.0001;
		VC<VC<PDD>> correct = convertDetections(detections);
		
		VVF features = extractAllFeatures(series, frames);
		for (VF &v : features) {
			REP(i, correct.SZ) if (resultsDistance(v, correct[i]) < TRAIN_CORRECT_DISTANCE) {
				v[1] = 1;
				if (splt(detections[i * 4], ' ')[7] == "1") v[2] = 1;
			}
			trainFeatures.PB(v);
		}
		
		trainCall++;
		
		return 0;
		// return trainCall >= 80 ? 1 : 0;
		// return 1;
	}
	
	VS testIDs;
	int testingData(string id, int width, int height,
			vector<int> data0, vector<string> head0, vector<double> wcs0,
			vector<int> data1, vector<string> head1, vector<double> wcs1,
			vector<int> data2, vector<string> head2, vector<double> wcs2,
			vector<int> data3, vector<string> head3, vector<double> wcs3) {
			
		DB(testCall);
		
		cerr << "AsteroidDetector: time = " << getTime() - t0 << endl;
		testIDs.PB(id);
		FullImageSet frames;
		Series* series = new Series(); //TODO: memory leak?
		series->id = id;
		assert(width == W);
		assert(height == H);
		vector<int>* pData[FRAMES] = { &data0, &data1, &data2, &data3 };
		vector<string>* pHead[FRAMES] = { &head0, &head1, &head2, &head3 };
		for (int f=0; f<FRAMES; f++) {
			for (int i=0; i<H; i++) for (int j=0; j<W; j++) {
				frames[f].at(i,j) = (*pData[f])[i*W + j];
			}
			for (string line : *pHead[f]) {
				series->frame_head[f].add_line(line);
			}
		}
		series->setup();
		
		VVF features = extractAllFeatures(series, frames);
		for (VF &v : features) {
			v[0] += testCall;
			testFeatures.PB(v);
		}
		
		testCall++;
		
		return 0;
	}
	
	string RFFolder;
	void saveRF() {
		if (RFFolder == "") {
			// not neccessary to save
			return;
		}

		VF trainResultsDetection;
		VF trainResultsNeo;

		for (VF v : trainFeatures) {
			trainResultsDetection.PB(v[1]);
		}

		DB(trainFeatures.SZ);

		double xtime = getTime();
		double rftrainTimer = 0;
		double neotrainTimer = 0;
		RandomForestConfig cfg;
		RandomForest RFDetection;
		RandomForest RFNeo;

		cfg.featuresIgnored = 11;
		cfg.bagSize = 1.5;
		cfg.randomFeatures = {1, 2, 3, 4, 6, 8};
		cfg.randomPositions = {4};
		cfg.maxNodeSize = 15;

		const int TREES_NO = 1000;

		RFDetection.train(trainFeatures, trainResultsDetection, cfg, TREES_NO);
		rftrainTimer = getTime() - xtime;
		DB(RFDetection.countTotalNodes());
		// Serialize
		ofstream ofs(RFFolder + "TrainedRF");
		boost::archive::text_oarchive oa(ofs);
		oa << RFDetection;
		ofs.close();

		/*xtime = getTime();
		cfg.maxNodeSize = 50;
		RFNeo.train(trainFeatures, trainResultsNeo, cfg, TREESNEO_NO);
		neotrainTimer = getTime() - xtime;
		DB(RFNeo.countTotalNodes());

		// Serialize
		ofstream neoOfs(RFFolder + "TrainedNeoRF");
		boost::archive::text_oarchive neoOa(neoOfs);
		neoOa << RFNeo;
		neoOfs.close();
		*/
		DB(rftrainTimer);
		DB(neotrainTimer);
	}

	int loadRF(RandomForest& RFDetection) {
		if (RFFolder.size() == 0) {
			return -1;
		}
		ifstream ifs(RFFolder + "TrainedRF");
		if (ifs.good()) {
			// Trained RF is found, deserialize
			boost::archive::text_iarchive ia(ifs);
			ia >> RFDetection;
			DB(RFDetection.countTotalNodes());
		} else {
			DBM("Failed to load the TrainedRF");
			return -1;
		}
		ifs.close();

		/*ifstream neoIfs(RFFolder + "TrainedNeoRF");
		if (neoIfs.good()) {
			// Trained RF is found, deserialize
			boost::archive::text_iarchive ia(neoIfs);
			ia >> RFNeo;
			DB(RFNeo.countTotalNodes());
		} else {
			DBM("Failed to load the TrainedNeoRF");
			return -1;
		}
		neoIfs.close();*/
		DBM("Success to load the RF data");
		return 0;
	}

	vector<string> getAnswer() {
		cerr << "AsteroidDetector: time = " << getTime() - t0 << endl;
		double xtime;
		RandomForestConfig cfg;
		cfg.featuresIgnored = 11;
		cfg.bagSize = 1.5;
		cfg.randomFeatures = {1, 2, 3, 4, 6, 8};
		cfg.randomPositions = {4};
		cfg.maxNodeSize = 10;
		
		VF trainResults;
		for (VF v : trainFeatures) trainResults.PB(v[1]);
		
		int goodSamples = 0;
		for (VF v : trainFeatures) goodSamples += v[1] == 1;
		DB(goodSamples);		
		DB(trainFeatures.SZ);
		DB(testFeatures.SZ);
		
		VC<pair<float,int>> vp;
		
		const int TREES_NO = 1000;
		
		RandomForest RF;
		int result = loadRF(RF);
		if (result != 0) {
			DBM("Failed to load the RF files. we use the train data on the fly");
			xtime = getTime();
			RF.train(trainFeatures, trainResults, cfg, TREES_NO);
		}

		DB(RF.countTotalNodes());

		REP(i, testFeatures.SZ) {
			float v = RF.estimate(testFeatures[i]);
			vp.PB(MP(v, i));
		}

		sort(vp.rbegin(), vp.rend());
		
		DB(vp.SZ);
		
		//TODO: move duplication removal here
		// VC<unordered_map<int,VI>> choice(100);
		
		VS rv;
		REP(i, vp.SZ) {
			if (rv.SZ >= MAX_ALLOWED_DETECTIONS) break;
			int sampleID = vp[i].Y;
			float score = vp[i].X;
			int testID = ((int)testFeatures[sampleID][0]);
			
			string s = testIDs[testID];
			char str[100];
			FOR(j, 3, 11) {
				
				sprintf(str, " %.10f", testFeatures[sampleID][j]);
				s += str;
			}
			s += " 0";
			
			sprintf(str, " %.10f", score);
			s += str;
			
			sprintf(str, " %.10f", (testFeatures[sampleID][21] + testFeatures[sampleID][22] + testFeatures[sampleID][23] + testFeatures[sampleID][24]) / 4);
			s += str;
			
			rv.PB(s);
		}
		
		DB(rv.SZ);
		
		cerr << "frontend: Finished, time = " << getTime()-t0 << endl;
		return rv;
	}
	
	vector<string> getAnswerNEO() {
		cerr << "AsteroidDetector: time = " << getTime() - t0 << endl;
		
		RandomForestConfig cfg;
		cfg.featuresIgnored = 11;
		cfg.bagSize = 1.5;
		cfg.randomFeatures = {1, 2, 3, 4, 6, 8};
		cfg.randomPositions = {4};
		cfg.maxNodeSize = 10;
		
		VF trainResults;
		for (VF v : trainFeatures) trainResults.PB(v[1]);
		
		int goodSamples = 0;
		for (VF v : trainFeatures) goodSamples += v[2] == 1; //This is the only place where the functions getAnswer() & getAnswerNEO() differ
		DB(goodSamples);		
		DB(trainFeatures.SZ);
		DB(testFeatures.SZ);
		
		VC<pair<float,int>> vp;
		
		const int TREES_NO = 1000;
		
		if (trainFeatures.SZ) {
			RandomForest RF;
			RF.train(trainFeatures, trainResults, cfg, TREES_NO);
			DB(RF.countTotalNodes());
			
			REP(i, testFeatures.SZ) {
				float v = RF.estimate(testFeatures[i]);
				vp.PB(MP(v, i));
			}
		} else {
			//that's how original pfr solution worked
			REP(i, testFeatures.SZ) {
				float v = testFeatures[i][10] * 52 + testFeatures[i][11] * 8 + testFeatures[i][12];
				vp.PB(MP(v, i));
			}
		}
		sort(vp.rbegin(), vp.rend());
		
		DB(vp.SZ);
		
		//TODO: move duplication removal here
		// VC<unordered_map<int,VI>> choice(100);
		
		VS rv;
		REP(i, vp.SZ) {
			if (rv.SZ >= MAX_ALLOWED_DETECTIONS) break;
			int sampleID = vp[i].Y;
			int testID = ((int)testFeatures[sampleID][0]);
			
			string s = testIDs[testID];
			FOR(j, 3, 11) {
				char str[100];
				sprintf(str, " %.10f", testFeatures[sampleID][j]);
				s += str;
			}
			s += " 0";
			rv.PB(s);
		}
		
		DB(rv.SZ);
		
		cerr << "frontend: Finished, time = " << getTime()-t0 << endl;
		return rv;
	}
	
};

