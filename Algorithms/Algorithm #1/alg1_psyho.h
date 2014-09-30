#include <bits/stdc++.h>
#include <sys/time.h>

const double MIN_DISTANCE = 0.00025;
const int MAX_RESULTS = 100000;
const double targetSamples = 2000; //2000
const double CORRECT_DISTANCE = 0.001;
const double TRAIN_CORRECT_DISTANCE = 0.0001;
const double MIN_OBJECT_DISTANCE = 0.000005;
const double MAX_OBJECT_MOVE = 0.007;
const double MAX_OBJECT_TOTAL_MOVE = 0.0125;
const int MAX_OBJECT_SIZE = 250;
const int TREES_NO = 450; //400
const int TREESNEO_NO = 500; //400
const int NEOS_USED = 6;

const int BACK_VALUE = ((1<<16)-1)/6;

using namespace std;

#define INLINE   inline __attribute__ ((always_inline))
#define NOINLINE __attribute__ ((noinline))

#define ALIGNED __attribute__ ((aligned(16)))

#define likely(x)   __builtin_expect(!!(x),1)
#define unlikely(x) __builtin_expect(!!(x),0)

#define SSELOAD(a)     _mm_load_si128((__m128i*)&a)
#define SSESTORE(a, b) _mm_store_si128((__m128i*)&a, b)

#define FOR(i,a,b)  for(int i=(a);i<(b);++i)
#define REP(i,a)    FOR(i,0,a)
#define ZERO(m)     memset(m,0,sizeof(m))
#define ALL(x)      x.begin(),x.end()
#define PB          push_back
#define S           size()
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

template<class T> void print(VC < T > v) {cerr << "[";if (v.S) cerr << v[0];FOR(i, 1, v.S) cerr << ", " << v[i];cerr << "]" << endl;}
template<class T> string i2s(T x) {ostringstream o; o << x; return o.str();}
VS splt(string s, char c = ' ') {VS all; int p = 0, np; while (np = s.find(c, p), np >= 0) {if (np != p) all.PB(s.substr(p, np - p)); p = np + 1;} if (p < s.S) all.PB(s.substr(p)); return all;}

double getTime() {
	timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + tv.tv_usec * 1e-6;
}

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
	public:
	VC < TreeNode > nodes;
	VD importances;
	
	DecisionTree() { }
	
	template < class T > DecisionTree(VC < VC < T > > &features, VC < T > &results, RandomForestConfig &config, int seed) {
		RNG r(seed);
		
		if (config.computeImportances) {
			importances = VD(features[0].S);
		}
	
		VI chosenGroups(features.S);
		REP(i, (int)(features.S * config.bagSize)) chosenGroups[r.next(features.S)]++;
		
		int bagSize = 0;
		REP(i, features.S) if (chosenGroups[i]) bagSize++;
		
		VI bag(bagSize);
		VI weight(features.S);
		
		int pos = 0;
		
		REP(i, features.S) {
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
		
		while (stack.S) {
			bool equal = true;
			
			int curNode = stack[stack.S - 1];
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
			
			const int randomFeatures = config.randomFeatures[min(nodes[curNode].level, (int)config.randomFeatures.S - 1)];
			REP(i, randomFeatures) {
			
				int featureID = config.featuresIgnored + r.next(features[0].S - config.featuresIgnored);
				
				T vlo, vhi;
				vlo = vhi = features[bag[bLeft]][featureID];
				FOR(j, bLeft + 1, bRight) {
					vlo = min(vlo, features[bag[j]][featureID]);
					vhi = max(vhi, features[bag[j]][featureID]);
				}
				if (vlo == vhi) continue;
				
				const int randomPositions = config.randomPositions[min(nodes[curNode].level, (int)config.randomPositions.S - 1)];
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
				importances[bestFeature] += 1. * (bRight - bLeft) / features.S;
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
			nodes[curNode].left = nodes.S;
			nodes[curNode].right = nodes.S + 1;
			
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
			
			stack.PB(nodes.S);
			stack.PB(nodes.S + 1);
			
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
	
	// template < class T > void addSample(VC < T > &features, T result, double weight = 1.0) {
		// TreeNode *pNode = &nodes[0];
		// while (true) {
			// if (pNode->feature < 0) {
				// pNode->result = (pNode->result * pNode->weight + result * weight) / (pNode->weight + weight);
				// pNode->weight += weight;
				// return;
			// }
			// pNode = &nodes[features[pNode->feature] < pNode->value ? pNode->left : pNode->right];
		// }
	// }
};

RNG gRNG(1);

class RandomForest {
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
		} else {
			thread *threads = new thread[config.threadsNo];
			mutex mutex;
			REP(i, config.threadsNo) 
				threads[i] = thread([&] {
					while (true) {
						auto tree = DecisionTree(features, results, config, gRNG.next());
						mutex.lock();
						if (trees.S < treesNo)
							trees.PB(tree);
						bool done = trees.S >= treesNo || config.timeLimit && getTime() - startTime > config.timeLimit;
						mutex.unlock();
						if (done) break;
					}
				});
			REP(i, config.threadsNo) threads[i].join();
			delete[] threads;
		}
		
		if (config.computeImportances) {
			importances = VD(features[0].S);
			for (DecisionTree tree : trees)
				REP(i, importances.S)
					importances[i] += tree.importances[i];
			double sum = 0;
			REP(i, importances.S) sum += importances[i];
			REP(i, importances.S) importances[i] /= sum;
		}
	}
	
	template <class T> double estimate(VC<T> &features) {
		assert(trees.S);
	
		double sum = 0;
		REP(i, trees.S) sum += trees[i].estimate(features);
		return sum / trees.S;
	}
	
	template <class T> VD estimate(VC<VC<T>> &features) {
		assert(trees.S);
	
		int samplesNo = features.S;
	
		VD rv(samplesNo);
		if (config.threadsNo == 1) {
			REP(j, samplesNo) {
				REP(i, trees.S) rv[j] += trees[i].estimate(features[j]);
				rv[j] /= trees.S;
			}
		} else {
			thread *threads = new thread[config.threadsNo];
			mutex mutex;
			int order = 0;
			REP(i, config.threadsNo) 
				threads[i] = thread([&] {
					mutex.lock();
					int offset = order++;
					mutex.unlock();
					for (int j = offset; j < samplesNo; j += config.threadsNo) {
						REP(k, trees.S) rv[j] += trees[k].estimate(features[j]);
						rv[j] /= trees.S;
					}
				});
			REP(i, config.threadsNo) threads[i].join();
			delete[] threads;
		}
		return rv;
	}
	
	LL countTotalNodes() {
		LL rv = 0;
		REP(i, trees.S) rv += trees[i].nodes.S;
		return rv;
	}
	
	// template < class T > void addSample(VC < T > &features, T result) {
		// REP(i, trees.S) {
			// int weight = gRNG.poisson(1);
			// if (weight > 0)	trees[i].addSample(features, result, weight);
		// }
	// }
	
};


//WCS stuff
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
double PI = 2 * acos(0);
double D2R = PI/180.0;
double R2D = 180.0/PI;

void init(VD &wcs) {
	crpix1 = wcs[0];
	crpix2 = wcs[1];
	crval1 = wcs[2];
	crval2 = wcs[3];
	cd11 = wcs[4];
	cd12 = wcs[5];
	cd21 = wcs[6];
	cd22 = wcs[7];
	double inv_det = cd11*cd22 - cd12*cd21;
	inv_det = 1.0 / inv_det;
	invcd11 = cd22 * inv_det;
	invcd12 = -cd12 * inv_det;
	invcd21 = -cd21 * inv_det;
	invcd22 = cd11 * inv_det;
}

PDD sphs2x(double lng, double lat) {
	double coslat, coslat3, coslat4, coslng, dlng, dphi, sinlat, sinlat3, sinlat4, sinlng, x, y, z;
	VD eul(5);
	PDD phiTheta;
	eul[0] = crval1;
	eul[1] = 90.0 - crval2;
	eul[2] = 180.0;
	eul[3] = cos(D2R * eul[1]);
	eul[4] = sin(D2R * eul[1]);
	/* Do lng dependency. */
	dlng = lng - eul[0];
	phiTheta.X = dlng;
	/* Do lat dependency. */
	sinlat = sin(lat*D2R);
	coslat = cos(lat*D2R);
	coslat3 = coslat*eul[3];
	coslat4 = coslat*eul[4];
	sinlat3 = sinlat*eul[3];
	sinlat4 = sinlat*eul[4];
	dlng = phiTheta.X;
	sinlng = sin(dlng*D2R);
	coslng = cos(dlng*D2R);
	/* Compute the native longitude. */
	x = sinlat4 - coslat3*coslng;
	y = -coslat*sinlng;
	dphi = R2D*atan2(y, x);
	phiTheta.X = eul[2] + dphi;
	if (phiTheta.X < 0) phiTheta.X += 360;
	if (phiTheta.X >= 360) phiTheta.X -= 360;
	
	/* Normalize the native longitude. */
	if (phiTheta.X > 180.0) {
		phiTheta.X -= 360.0;
	} else if (phiTheta.X < -180.0) {
		phiTheta.X += 360.0;
	}
	/* Compute the native latitude. */
	z = sinlat3 + coslat4*coslng;
	if (abs(z) > 0.99) {
		if (z < 0.0)
			phiTheta.Y = -abs(R2D * acos(sqrt(x*x+y*y)));
		else
			phiTheta.Y = abs(R2D * acos(sqrt(x*x+y*y)));
	} else {
		phiTheta.Y = R2D * asin(z);
	}
	return phiTheta;
}

PDD sphx2s(double phi, double theta) {
	double cosphi, costhe, costhe3, costhe4, dlng, dphi, sinphi, sinthe, sinthe3, sinthe4, x, y, z;
	PDD lngLat;
	VD eul(5);
	eul[0] = crval1;
	eul[1] = 90.0 - crval2;
	eul[2] = 180.0;
	eul[3] = cos(D2R * eul[1]);
	eul[4] = sin(D2R * eul[1]);
	/* Do phi dependency. */
	dphi = phi - eul[2];
	/* Do theta dependency. */
	sinthe = sin(theta * D2R);
	costhe = cos(theta * D2R);
	costhe3 = costhe*eul[3];
	costhe4 = costhe*eul[4];
	sinthe3 = sinthe*eul[3];
	sinthe4 = sinthe*eul[4];
	sinphi = sin(dphi * D2R);
	cosphi = cos(dphi * D2R);
	/* Compute the celestial longitude. */
	x = sinthe4 - costhe3*cosphi;
	y = -costhe*sinphi;
	dlng = R2D * atan2( y, x);
	lngLat.X = eul[0] + dlng;
	/* Normalize the celestial longitude. */
	if (eul[0] >= 0.0) {
		if (lngLat.X < 0.0)
			lngLat.X += 360.0;
	} else {
		if (lngLat.X > 0.0)
			lngLat.X -= 360.0;
	}
	if (lngLat.X > 360.0) {
		lngLat.X -= 360.0;
	} else if (lngLat.X < -360.0) {
		lngLat.X += 360.0;
	}
	/* Compute the celestial latitude. */
	z = sinthe3 + costhe4*cosphi;
	if (abs(z) > 0.99) {
	/* Use an alternative formula for greater accuracy. */
		if (z<0.0)
			lngLat.Y = -abs(R2D * acos(sqrt(x*x+y*y)));
		else
			lngLat.Y = abs(R2D * acos(sqrt(x*x+y*y)));
	} else {
		lngLat.Y = R2D * asin(z);
	}
	return lngLat;
}


PDD tans2x(double phi, double theta) {
	PDD xy;
	double cotan_theta = 1.0 / tan(D2R * theta);
	xy.X = R2D * sin(D2R * phi) * cotan_theta ;
	xy.Y = -R2D * cos(D2R * phi) * cotan_theta;
	return xy;
}

PDD tanx2s(double x, double y) {
	PDD phiTheta;
	phiTheta.X = R2D * atan2(x, -y);
	phiTheta.Y = R2D * atan2(R2D, sqrt(x*x + y*y) );
	return phiTheta;
}

PDD convertRADEC2XY(VD &wcs, double RA, double DEC) {
	init(wcs);
	PDD phiTheta = sphs2x(RA, DEC);
	PDD XXYY = tans2x(phiTheta.X, phiTheta.Y);
	PDD XY;
	XY.X = invcd11*XXYY.X + invcd12*XXYY.Y + crpix1;
	XY.Y = invcd21*XXYY.X + invcd22*XXYY.Y + crpix2;
	return XY;
}

PDD convertXY2RADEC(VD &wcs, double x, double y) {
	init(wcs);
	double dx = x-crpix1;
	double dy = y-crpix2;
	double xx = cd11*dx + cd12*dy;
	double yy = cd21*dx + cd22*dy;
	PDD phiTheta = tanx2s(xx, yy);
	PDD RD = sphx2s(phiTheta.X, phiTheta.Y);
	return RD;
}

struct WCSData {
	static const int GRID_POINTS = 500;
	bool valid;
	double minValueX;
	double maxValueX;
	double minValueY;
	double maxValueY;
	PDD grid[GRID_POINTS+1][GRID_POINTS+1];
	VD wcs;
	bool isXY;
};

bool wrapRA;
bool wrapDEC;
WCSData wcsRADEC[4];
WCSData wcsXY[4];

const double BAD_COORD = -1e9;
const double MAX_COORD_DIFF = 10;

void createWCSData(VVD &wcs, int width, int height) {
	REP(i, 4) {
		wcsRADEC[i].minValueX = 1e9;
		wcsRADEC[i].minValueY = 1e9;
		wcsRADEC[i].maxValueX = -1e9;
		wcsRADEC[i].maxValueY = -1e9;

		wcsXY[i].minValueX = 0;
		wcsXY[i].minValueY = 0;
		wcsXY[i].maxValueX = width;
		wcsXY[i].maxValueY = height;
		REP(j, WCSData::GRID_POINTS+1) REP(k, WCSData::GRID_POINTS+1) {
			double x = wcsXY[i].minValueX + (wcsXY[i].maxValueX - wcsXY[i].minValueX) / WCSData::GRID_POINTS * j;
			double y = wcsXY[i].minValueY + (wcsXY[i].maxValueY - wcsXY[i].minValueY) / WCSData::GRID_POINTS * k;
			PDD p = convertXY2RADEC(wcs[i], x, y);
			wcsXY[i].grid[k][j] = p;
			wcsRADEC[i].minValueX = min(wcsRADEC[i].minValueX, p.X);
			wcsRADEC[i].maxValueX = max(wcsRADEC[i].maxValueX, p.X);
			wcsRADEC[i].minValueY = min(wcsRADEC[i].minValueY, p.Y);
			wcsRADEC[i].maxValueY = max(wcsRADEC[i].maxValueY, p.Y);
		}
		
		bool valid = true;
		if (wcsRADEC[i].maxValueX - wcsRADEC[i].minValueX > MAX_COORD_DIFF) valid = false;
		if (wcsRADEC[i].maxValueY - wcsRADEC[i].minValueY > MAX_COORD_DIFF) valid = false;
		
		wcsXY[i].valid = wcsRADEC[i].valid = valid;
		
		REP(j, WCSData::GRID_POINTS+1) REP(k, WCSData::GRID_POINTS+1) {
			double x = wcsRADEC[i].minValueX + (wcsRADEC[i].maxValueX - wcsRADEC[i].minValueX) / WCSData::GRID_POINTS * j;
			double y = wcsRADEC[i].minValueY + (wcsRADEC[i].maxValueY - wcsRADEC[i].minValueY) / WCSData::GRID_POINTS * k;
			PDD p = convertRADEC2XY(wcs[i], x, y);
			wcsRADEC[i].grid[k][j] = p;
		}
		
		wcsXY[i].wcs = wcs[i];
		wcsXY[i].isXY = true;
		wcsRADEC[i].wcs = wcs[i];
		wcsRADEC[i].isXY = false;
		
	}
}

PDD convertCoord(WCSData &wcs, double x, double y) {
	if (!wcs.valid)	return wcs.isXY ? convertXY2RADEC(wcs.wcs, x, y) : convertRADEC2XY(wcs.wcs, x, y);

	x = (x - wcs.minValueX) / (wcs.maxValueX - wcs.minValueX) * WCSData::GRID_POINTS;
	y = (y - wcs.minValueY) / (wcs.maxValueY - wcs.minValueY) * WCSData::GRID_POINTS;
	if (x < 0 || x >= WCSData::GRID_POINTS || y < 0 || y >= WCSData::GRID_POINTS) return MP(BAD_COORD, BAD_COORD);
	int ix = (int)x;
	int iy = (int)y;
	double rx = wcs.grid[iy][ix].X + (wcs.grid[iy][ix+1].X - wcs.grid[iy][ix].X) * (x - ix) + (wcs.grid[iy+1][ix].X - wcs.grid[iy][ix].X) * (y - iy);
	double ry = wcs.grid[iy][ix].Y + (wcs.grid[iy][ix+1].Y - wcs.grid[iy][ix].Y) * (x - ix) + (wcs.grid[iy+1][ix].Y - wcs.grid[iy][ix].Y) * (y - iy);
	return MP(rx, ry);
}

	
template <class T> double getValue(VC<T> &img, int w, double x, double y) {
	int x0 = (int)x;
	int y0 = (int)y;
	if (x0 < 0 || x0 >= w - 1 || y0 < 0 || y0 >= img.S / w - 1) return -1;
	double wx = 1 - (x - x0);
	double wy = 1 - (y - y0);
	double v = 0;
	v += img[y0 * w + x0] * wx * wy;
	v += img[y0 * w + x0 + 1] * (1 - wx) * wy;
	v += img[(y0 + 1) * w + x0] * wx * (1 - wy);
	v += img[(y0 + 1) * w + x0 + 1] * (1 - wx) * (1 - wy);
	return v;
}

template <class T> double getValue3x3(VC<T> &img, int w, double x, double y) {
	double sumW = 0;
	double sum = 0;
	FOR(dx, -1, 2) FOR(dy, -1, 2) {
		double v = getValue(img, w, x + dx, y + dy);
		if (v < 0) continue;
		double w = 1 * (dx ? 0.66 : 1) * (dy ? 0.66 : 1);
		sum += w * v;
		sumW += w;
	}
	return sumW == 0 ? -1 : sum / sumW;
}

template <class T> double getValue5x5(VC<T> &img, int w, double x, double y) {
	double sumW = 0;
	double sum = 0;
	FOR(dx, -2, 3) FOR(dy, -2, 3) {
		double v = getValue(img, w, x + dx, y + dy);
		if (v < 0) continue;
		double w = 1 * (abs(dx) == 2 ? 0.4 : dx ? 0.66 : 1) * (abs(dy) == 2 ? 0.4 : dy ? 0.66 : 1);
		sum += w * v;
		sumW += w;
	}
	return sumW == 0 ? -1 : sum / sumW;
}

/*
template <class T> VC<T> medianFilter(VC<T> &img, int w, int size) {
	T data[500];
	int h = img.S / w;
	VC<T> rv(img.S);
	REP(y, h) REP(x, w) {
		int no = 0;
		FOR(dy, -size, +size + 1) FOR(dx, -size, +size + 1) {
			int ny = y + dy;
			int nx = x + dx;
			if (ny >= 0 && ny < h && nx >= 0 && nx < w)
				data[no++] = img[ny * w + nx];
		}
		nth_element(data, data + no / 2, data + no);
		rv[y * w + x] = data[no/2];
	}
	return rv;
}
*/

template <class T> VC<T> medianFilterFast(VC<T> &img, int w) {
	int h = img.S / w;
	
	VC<T> rv(img.S);
	
	unsigned short data[500];
	
	VC<VC<T>> avg(h / 3 + 1, VC<T>(w / 3 + 1));
	REP(y, h/3 + 1) REP(x, w/3 + 1) {
		int sum = 0;
		int sumW = 0;
		FOR(yy, y*3, y*3+3) FOR(xx, x*3, x*3+3) {
			if (yy >= h || xx >= w) continue;
			sumW++;
			sum += img[yy*w+xx];
		}
		avg[y][x] = max(0, min((1<<16)-1, sumW == 0 ? 0 : sum / sumW));
	}
	
	REP(y, h) REP(x, w) {
		if (x % 3 || y % 3) {
			rv[y*w+x] = rv[(y-y%3)*w+x-x%3];
			continue;
		}
		int no = 0;
		int yy = y / 3;
		int xx = x / 3;
		FOR(dy, -4, 5) FOR(dx, -4, 5) {
			int ny = yy + dy;
			int nx = xx + dx;
			if (ny >= 0 && ny < avg.S && nx >= 0 && nx < avg[0].S)
				data[no++] = avg[ny][nx];
		}
		nth_element(data, data + no / 2, data + no);
		rv[y * w + x] = data[no/2];
	}
	return rv;
}

template <class T> void stretchFilter(VC<T> &img, VC<T> &back, int w) {
	int MAX = (1 << 16) - 1;
	
	int h = img.S / w;

	int spread = 925;
	
	REP(y, h) REP(x, w) {
		int med = back[x+y*w];
		int dmin = med - spread / 6;
		int dmax = dmin + spread;
		int ival = img[x+y*w];
		ival = max(0, min(MAX, (MAX * (ival - dmin)) / (dmax - dmin)));
		img[x+y*w] = ival;
	}
}

template <class T> void clearObject(VC<T> &img, int w, int x, int y) {
	const int h = img.S / w;
	const int dx[] = {-1,1,0,0};
	const int dy[] = {0,0,-1,1};
	queue<int> q;
	q.push(x), q.push(y);
	// assert(img[y*w+x] != 0);
	while (!q.empty()) {
		int px = q.front(); q.pop();
		int py = q.front(); q.pop();
		if (px < 0 || px >= w || py < 0 || py >= h) continue;
		if (img[py*w+px] == 0) continue;
		img[py*w+px] = 0;
		REP(d, 4) q.push(px+dx[d]), q.push(py+dy[d]);
	}
}

template <class T> VC<T> findObjects(VC<T> &img, int w, double tolerance = 0.048) {
	int h = img.S / w;
	VC<T> rv(img.S);
	REP(y, h - 1) REP(x, w - 1) {
		int p = y*w+x;
		int sumA = (int)img[p] + (int)img[p+1] + (int)img[p+w] + (int)img[p+w+1];
		bool full = false;
		full |= img[p] == (1<<16) - 1;
		full |= img[p+1] == (1<<16) - 1;
		full |= img[p+w] == (1<<16) - 1;
		full |= img[p+w+1] == (1<<16) - 1;
		if (full || (sumA - 4 * BACK_VALUE) / 4 > (1 << 16) * tolerance) {
			rv[p] = img[p];
			rv[p+1] = img[p+1];
			rv[p+w] = img[p+w];
			rv[p+w+1] = img[p+w+1];
		}
		
		// if (x == 0 || y == 0) continue;
		// int sum2A = (int)img[p-w] + (int)img[p-1] + (int)img[p] + (int)img[p+1] + (int)img[p+w];
		// if ((sum2A - 5 * BACK_VALUE) / 5 > (1 << 16) * tolerance) {
			// rv[p-w] = img[p-w];
			// rv[p-1] = img[p-1];
			// rv[p] = img[p];
			// rv[p+1] = img[p+1];
			// rv[p+w] = img[p+w];
		// }
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

template <class T> VC<pair<VC<PDD>, VF>> findMovingObjects(VC<VC<T>> &img, VC<VC<T>> &orig, VVD &wcs, int w, bool removeRepeats, double tolerance = 1.12, double movingThreshold = 0.7) {
	int h = img[0].S / w;
	VC<VC<T>> rv = img;
	
	//clear stationary objects
	VC<PDD> objects[4];
	VF objectsSize[4];
	VF objectsBrightness[4];
	VF objectsRelBrightness[4];
	REP(i, 4) {
		VI vs(img[0].S, 0);
		REP(y, h) REP(x, w) {
			int p = y*w+x;
			if (vs[p] || rv[i][p] == 0) continue;
			queue<int> q;
			q.push(x), q.push(y);
			int total = 0;
			VD sum(4, 0);
			bool error = false;
			double sumRA = 0;
			double sumDEC = 0;
			double sumW = 0;
			const int dx[] = {-1,1,0,0};
			const int dy[] = {0,0,-1,1};
			while (!q.empty()) { 
				int px = q.front(); q.pop();
				int py = q.front(); q.pop();
				if (px < 0 || px >= w || py < 0 || py >= h) continue;
				if (vs[py*w+px] || img[i][py*w+px] == 0) continue;
				vs[py*w+px] = 1;
				// PDD radeg = convertXY2RADEC(wcs[i], px, py); 
				PDD radeg = convertCoord(wcsXY[i], px, py); 
				REP(j, 4) {
					// PDD xy = convertRADEC2XY(wcs[j], radeg.X, radeg.Y);
					PDD xy = convertCoord(wcsRADEC[j], radeg.X, radeg.Y);
					double v = getValue(orig[j], w, xy.X, xy.Y);
					if (v == -1) {
						error = true;
						break;
					}
					sum[j] += max(v - BACK_VALUE, 0.0);
				}
				if (error) break;
				total++;
				double weight = img[i][py*w+px];
				sumRA += radeg.X * weight;
				sumDEC += radeg.Y * weight;
				sumW += weight;
				REP(d, 4) q.push(px + dx[d]), q.push(py + dy[d]);
			}
			
			bool valid = true;
			if (total > MAX_OBJECT_SIZE) valid = false;
			if (error) valid = false;
			
			double maxSum = 0;
			REP(j, 4) if (i != j) maxSum = max(maxSum, sum[j]);
			valid &= maxSum < sum[i] * tolerance;
			
			if (!valid) {
				clearObject(rv[i], w, x, y);
				continue;
			}
			
			/*
			PDD xy = convertCoord(wcsRADEC[i], sumRA / sumW, sumDEC / sumW);
			
			if (xy.X <= -1e8) continue;
			
			double bv = getValue3x3(orig[i], w, xy.X, xy.Y);
			double bx = xy.X, by = xy.Y;
			FOR(xx, -2, 3) FOR(yy, -2, 3) {
				double nx = xy.X + xx * 0.15;
				double ny = xy.Y + yy * 0.15;
				double av = getValue3x3(orig[i], w, nx, ny);
				if (av > bv + 1e-3) {
					bv = av;
					bx = nx;
					by = ny;
				}
			}
			
			PDD radec = convertCoord(wcsXY[i], bx, by);
			if (radec.X <= -1e8) continue;
			
			objects[i].PB(radec);
			*/
			
			objects[i].PB(MP(sumRA / sumW, sumDEC / sumW));
			objectsSize[i].PB(total);
			objectsBrightness[i].PB(sum[i]);
			objectsRelBrightness[i].PB(sum[i] / maxSum);
		}
	}
	
	double minRA = 1e9;
	double maxRA = -1e9;
	double minDEC = 1e9;
	double maxDEC = -1e9;
	REP(i, 4) for (PDD p : objects[i]) {
		minRA = min(minRA, p.X);
		maxRA = max(maxRA, p.X);
		minDEC = min(minDEC, p.Y);
		maxDEC = max(maxDEC, p.Y);
	}
	
	int bucketsRA =  (int)((maxRA - minRA)   / MAX_OBJECT_TOTAL_MOVE) + 1;
	int bucketsDEC = (int)((maxDEC - minDEC) / MAX_OBJECT_TOTAL_MOVE) + 1;
	VVVVI buckets(4, VVVI(bucketsRA));
	VVVI bucketsRes(bucketsRA);
	REP(i, 4) REP(j, objects[i].S) {
		PDD &p = objects[i][j];
		if (buckets[i][(p.X - minRA) / MAX_OBJECT_TOTAL_MOVE].S == 0) buckets[i][(p.X - minRA) / MAX_OBJECT_TOTAL_MOVE] = VVI(bucketsDEC);
		buckets[i][(p.X - minRA) / MAX_OBJECT_TOTAL_MOVE][(p.Y - minDEC) / MAX_OBJECT_TOTAL_MOVE].PB(j);
	}
	
	VC<int> objectsUsed[4];
	REP(i, 4) objectsUsed[i] = VI(objects[i].S, 0);
	
	int movObjectsTested = 0;
	int movObjectsCorrect = 0;
	
	VC<pair<VC<PDD>, VF>> results;
	VI resultsRepeat;
	
	for (int dd = 3; dd >= 1; dd--) REP(i, 4) REP(x, bucketsRA) if (buckets[i][x].S) REP(y, bucketsDEC) for (int p : buckets[i][x][y]) 
		for (int j = max(0, i - dd); j <= i - dd; j++) FOR(x2, max(0, x - 1), min(bucketsRA, x + 2)) if (buckets[j][x2].S) FOR(y2, max(0, y - 1), min(bucketsDEC, y + 2)) for (int p2 : buckets[j][x2][y2]) {
			double x1 = objects[i][p].X;
			double y1 = objects[i][p].Y;
			double x2 = objects[j][p2].X;
			double y2 = objects[j][p2].Y;
			double dist = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
			if (dist > min(MAX_OBJECT_TOTAL_MOVE * MAX_OBJECT_TOTAL_MOVE, MAX_OBJECT_MOVE * MAX_OBJECT_MOVE * (i - j) * (i - j))) continue;
			movObjectsTested++;
			
			double dx = (x2 - x1) / (j - i);
			double dy = (y2 - y1) / (j - i);
			double x0 = x1 - dx * i;
			double y0 = y1 - dy * i;
			double v[4];
			REP(k, 4) {
				PDD p = convertCoord(wcsRADEC[k], x0 + k * dx, y0 + k * dy);
				v[k] = getValue(rv[k], w, p.X, p.Y);
			}
			
			double mn = 1e9, mn2 = 1e9;
			double ok = min(v[i], v[j]);
			REP(k, 4) {
				if (v[k] < mn) {
					mn2 = mn;
					mn = v[k];
				} else if (v[k] < mn2) {
					mn2 = v[k];
				}
			}
			
			bool correct = mn2 > ok * movingThreshold && mn >= 0;
			movObjectsCorrect += correct;
			if (correct) {
				objectsUsed[i][p] = objectsUsed[j][p2] = 1;
				VC<PDD> vp(4);
				REP(k, 4) {
					vp[k] = MP(x0 + k * dx, y0 + k * dy);
				}
				
				int bx = (vp[0].X - minRA)  / MAX_OBJECT_TOTAL_MOVE;
				int by = (vp[0].Y - minDEC) / MAX_OBJECT_TOTAL_MOVE;
				
				if (bx < 0 || bx >= bucketsRA || by < 0 || by >= bucketsDEC) continue;
				
				bool repeats = false;
				FOR(x2, max(0, bx - 1), min(bucketsRA, bx + 2)) if (bucketsRes[x2].S) FOR(y2, max(0, by - 1), min(bucketsDEC, by + 2)) for (int p : bucketsRes[x2][y2])
					if (!resultsRepeat[p] && resultsDistance(results[p].X, vp) < MIN_OBJECT_DISTANCE) {
						repeats = true;
						break;
					}
				if (removeRepeats && repeats) continue;
				
				VF bf = {repeats ? 1.0f : 0.0f, (float)(ok / mn2), objectsSize[i][p] + objectsSize[j][p2], objectsBrightness[i][p] + objectsBrightness[j][p2], objectsRelBrightness[i][p] + objectsRelBrightness[j][p2]};
				results.PB(MP(vp, bf));
				resultsRepeat.PB(repeats);
				if (bucketsRes[bx].S == 0) bucketsRes[bx] = VVI(bucketsDEC);
				bucketsRes[bx][by].PB(results.S - 1);
				
			}
		}
	
	return results;
}

template <class T> VC<VC<T>> alignImages(VC<VC<T>> &img, VVD &wcs, int w, int margin) {
	int h = img[0].S / w;
	VC<VC<T>> rv(4, VC<T>((w - 2 * margin) * (h - 2 * margin)));
	int nw = w - 2 * margin;
	REP(y, h - 2 * margin) REP(x, w - 2 * margin) {
		rv[0][y * nw + x] = img[0][(y + margin) * w + (x + margin)];
		// PDD rd = convertXY2RADEC(wcs[0], x + margin, y + margin);
		PDD rd = convertCoord(wcsXY[0], x + margin, y + margin);
		FOR(i, 1, 4) {
			// PDD p = convertRADEC2XY(wcs[i], rd.X, rd.Y);
			PDD p = convertCoord(wcsRADEC[i], rd.X, rd.Y);
			assert(p.X >= 0 && p.X < w - 1);
			assert(p.Y >= 0 && p.Y < h - 1);
			int x0 = (int)p.X;
			int y0 = (int)p.Y;
			double wx = 1 - (p.X - x0);
			double wy = 1 - (p.Y - y0);
			double v = 0;
			v += img[i][y0 * w + x0] * wx * wy;
			v += img[i][y0 * w + x0 + 1] * (1 - wx) * wy;
			v += img[i][(y0 + 1) * w + x0] * wx * (1 - wy);
			v += img[i][(y0 + 1) * w + x0 + 1] * (1 - wx) * (1 - wy);
			rv[i][y * nw + x] = (int)(v + 0.5);
		}
	}
	return rv;
}

double average(VD &v) {
	double x = 0; REP(i, v.S) x += v[i]; return x / v.S;
}

double stddev(VD &v) {
	double avg = average(v);
	double x = 0; 
	REP(i, v.S) x += (v[i] - avg) * (v[i] - avg); 
	return sqrt(x / v.S);
}

template <class T> VF extractFeatures(VC<PDD> &v, VF bf, int w, VVD &wcs, VC<VC<T>> &obj, VC<VC<T>> &orig, VC<LD> &dates, bool original) {
	VF rv;
	
	//reserved
	rv.PB(original ? 100 : bf[0] ? 200 : 0);
	rv.PB(0);
	rv.PB(0);
	
	//add solution
	REP(i, v.S) {
		rv.PB(v[i].X);
		rv.PB(v[i].Y);
	}
	
	//bonus features
	FOR(i, 1, bf.S) rv.PB(bf[i]);
	
	//speed
	rv.PB(sqrt((v[0].X - v[1].X) * (v[0].X - v[1].X) + (v[0].Y - v[1].Y) * (v[0].Y - v[1].Y)));
	rv.PB(v[0].X - v[1].X);
	rv.PB(v[0].Y - v[1].Y);
	rv.PB(atan2(v[0].Y - v[1].Y, v[0].X - v[1].X));	
	
	//speed / time
	if (dates[0] > 0 && dates[3] > 0) {
		rv.PB(sqrt(((v[0].X - v[1].X) * (v[0].X - v[1].X) + (v[0].Y - v[1].Y) * (v[0].Y - v[1].Y))) / (dates[3] - dates[0]));
		rv.PB((v[0].X - v[1].X) / (dates[3] - dates[0]));
		rv.PB((v[0].Y - v[1].Y) / (dates[3] - dates[0]));
	} else {
		rv.PB(0);
		rv.PB(0);
		rv.PB(0);
	}
	
	//time between frames
	rv.PB((dates[3] - dates[0]) / 3);
	
	//date
	rv.PB(dates[0]);
	
	LL yearLength = 366 * 24 * 60 * 60;
	rv.PB(dates[0] - (int)(dates[0] / yearLength) * yearLength);
	
	
	
	//coordinates
	rv.PB(v[0].X);
	rv.PB(v[0].Y);
	
	//brightness at center
	// VD objBrightness(4);
	// REP(i, 4) {
		// PDD p = convertRADEC2XY(wcs[i], v[i].X, v[i].Y);
		// objBrightness[i] = getValue(obj[i], w, p.X, p.Y);
	// }
	// sort(ALL(objBrightness));
	// rv.PB(objBrightness[0]);
	// rv.PB(objBrightness[1]);
	// rv.PB(objBrightness[2]);
	// rv.PB(objBrightness[3]);
	
	//brightness at center
	VD objBrightness(4);
	VD objBrightness2(4);
	VD objBrightness3(4);
	VD maxBrightness(4);
	REP(i, 4) {
		PDD p = convertCoord(wcsRADEC[i], v[i].X, v[i].Y);
		objBrightness[i] = getValue(orig[i], w, p.X, p.Y);
		objBrightness2[i] = getValue3x3(orig[i], w, p.X, p.Y);
		maxBrightness[i] = getValue3x3(obj[i], w, p.X, p.Y);
		// objBrightness3[i] = getValue5x5(orig[i], w, p.X, p.Y);
	}
	
	
	VD objBrightnessDiff(4);
	VD objBrightnessDiff2(4);
	VD objBrightnessDiff3(4);
	REP(i, 4) {
		// objBrightnessDiff[i] = 1e9;
		// REP(j, 4) if (i != j) {
			// PDD p = convertCoord(wcsRADEC[j], v[i].X, v[i].Y);
			// objBrightnessDiff[i] = min(objBrightnessDiff[i], getValue(orig[j], w, p.X, p.Y));
		// }
		// objBrightnessDiff[i] = objBrightness[i] - objBrightnessDiff[i]; 
		
		objBrightnessDiff2[i] = 1e9;
		REP(j, 4) if (i != j) {
			PDD p = convertCoord(wcsRADEC[j], v[i].X, v[i].Y);
			objBrightnessDiff2[i] = min(objBrightnessDiff2[i], getValue3x3(orig[j], w, p.X, p.Y));
		}
		
		objBrightnessDiff[i] = objBrightness2[i] - objBrightnessDiff2[i]; 
		
		// objBrightnessDiff3[i] = 1e9;
		// REP(j, 4) if (i != j) {
			// PDD p = convertCoord(wcsRADEC[j], v[i].X, v[i].Y);
			// objBrightnessDiff3[i] = min(objBrightnessDiff3[i], getValue5x5(orig[j], w, p.X, p.Y));
		// }
	}
	
	// REP(i, 4) {
		// objBrightnessDiff2[i] = -1e9;
		// REP(j, 4) if (i != j) {
			// PDD p = convertCoord(wcsRADEC[j], v[i].X, v[i].Y);
			// objBrightnessDiff2[i] = max(objBrightnessDiff2[i], getValue(orig[j], w, p.X, p.Y));
		// }
		// objBrightnessDiff2[i] = objBrightness[i] - objBrightnessDiff2[i]; 
	// }
	
	// sort(ALL(objBrightness));
	// rv.PB(objBrightness[0]);
	// rv.PB(objBrightness[1]);
	// rv.PB(objBrightness[2]);
	// rv.PB(objBrightness[3]);
	
	sort(ALL(maxBrightness));
	rv.PB(maxBrightness[0]);
	rv.PB(maxBrightness[1]);
	rv.PB(maxBrightness[2]);
	rv.PB(maxBrightness[3]);
	
	sort(ALL(objBrightness2));
	rv.PB(objBrightness2[0]);
	rv.PB(objBrightness2[1]);
	rv.PB(objBrightness2[2]);
	rv.PB(objBrightness2[3]);
	
	sort(ALL(objBrightnessDiff));
	rv.PB(objBrightnessDiff[0]);
	rv.PB(objBrightnessDiff[1]);
	rv.PB(objBrightnessDiff[2]);
	rv.PB(objBrightnessDiff[3]);
	
	sort(ALL(objBrightnessDiff2));
	rv.PB(objBrightnessDiff2[0]);
	rv.PB(objBrightnessDiff2[1]);
	rv.PB(objBrightnessDiff2[2]);
	rv.PB(objBrightnessDiff2[3]);
	
	// sort(ALL(objBrightness3));
	// rv.PB(objBrightness3[0]);
	// rv.PB(objBrightness3[1]);
	// rv.PB(objBrightness3[2]);
	// rv.PB(objBrightness3[3]);
	
	// sort(ALL(objBrightnessDiff3));
	// rv.PB(objBrightnessDiff3[0]);
	// rv.PB(objBrightnessDiff3[1]);
	// rv.PB(objBrightnessDiff3[2]);
	// rv.PB(objBrightnessDiff3[3]);
	
	
	//
	//day?
	//noise?
	
	return rv;
}

VC<VC<PDD>> convertDetections(VS &v) {
	VC<VC<PDD>> rv;
	VC<PDD> current;
	for (string s : v) {
		VS vs = splt(s, ' ');
		double ra = atof(vs[2].c_str());
		double dec = atof(vs[3].c_str());
		current.PB(MP(ra, dec));
		if (current.S == 4) {
			rv.PB(current);
			current.clear();
		}
	}
	return rv;
}

LD convertDate(string s) {
	const int monthLen[] = {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
	int monthSum[12];
	ZERO(monthSum);
	FOR(i, 1, 12) monthSum[i] = monthSum[i-1] + monthLen[i-1];
	
	VS v = splt(s, 'T');
	
	if (v.S != 2) return 0;
	
	VS v0 = splt(v[0], '-');
	VS v1 = splt(v[1], ':');
	if (v0.S != 3) return 0;
	if (v1.S != 3) return 0;
	
	int year = atoi(v0[0].c_str());
	int month = atoi(v0[1].c_str());
	int day = atoi(v0[2].c_str());
	int hour = atoi(v1[0].c_str());
	int minute = atoi(v1[1].c_str());
	int seconds = atoi(v1[2].c_str());
	if (month <= 0 || month > 12) return 0;
	
	int dayOfYear = monthSum[month - 1] + day;
	return (LD)seconds + 60.0 * (minute + 60.0 * (hour + 24.0 * ((LD)dayOfYear + 366.0 * year)));
}

string extractDate(VS &header) {
	for (string s : header) {
		if (s.find("DATE    =") == string::npos) continue;
		return splt(s, ' ')[2];
	}
	return "";
}


VVF xFeatures;

template <class T> void solve(VC<VC<T>> &data, VVD &wcs, int width, VS detections, VVS headers, bool train) {
	double startTime;
	double totalTimer = 0;
	double stretchFilterTimer = 0;
	double correctOffsetTimer = 0;
	double medianFilterTimer = 0;
	double findObjectsTimer = 0;
	double findMovingObjectsTimer = 0;
	
	double entryTime = getTime();

	assert(sizeof(T) == 2);
	
	createWCSData(wcs, width, data[0].S / width);
	
	VC<VC<unsigned short>> objects(4);
	
	REP(i, 4) {
		startTime = getTime();
		VC<unsigned short> back = medianFilterFast(data[i], width);
		medianFilterTimer += getTime() - startTime;
		
		startTime = getTime();
		stretchFilter(data[i], back, width);
		stretchFilterTimer += getTime() - startTime;
	}
	
	REP(i, 4) {
		startTime = getTime();
		objects[i] = findObjects(data[i], width);
		findObjectsTimer += getTime() - startTime;
	}
	
	startTime = getTime();
	auto xDetected = findMovingObjects(objects, data, wcs, width, train);
	findMovingObjectsTimer += getTime() - startTime;
	
	VC<LD> dates;
	REP(i, 4) dates.PB(convertDate(extractDate(headers[i])));
	
	xFeatures.clear();
	for (auto v : xDetected) xFeatures.PB(extractFeatures(v.X, v.Y, width, wcs, objects, data, dates, false));
	
	// VS neodetections;
	// for (string s : detections) if (s.S && s[s.S-1] == '1') neodetections.PB(s);
	// VC<VC<PDD>> correct = convertDetections(neodetections);
	// for (auto v : correct) xFeatures.PB(extractFeatures(v, {0, 0, 0, 0, 0}, width, wcs, objects, data, dates, true));
	
	// VVF oldFeatures = xFeatures;
	// xFeatures.clear();
	// FOR(i, -3, 4) {
		// for (auto v : oldFeatures) {
			// auto u = v;
			// REP(j, 4) u[j*2+3] += i * 0.02;
			// u.PB(i);
			// xFeatures.PB(u);
		// }
	// }
	
	// xFeatures.clear();
	// FOR(i, -3, 4) {
		// for (auto v : xDetected) {
			// auto u = v;
			// REP(j, 4) u[j].X += i * 0.02;
			// VF vd = extractFeatures(u, width, wcs, objects, data, back);
			// vd.PB(i);
			// xFeatures.PB(vd);
		// }
	// }
	
	totalTimer = getTime() - entryTime;
	
	// DB(stretchFilterTimer);
	// DB(correctOffsetTimer);
	// DB(medianFilterTimer);
	// DB(findObjectsTimer);
	// DB(findMovingObjectsTimer);
	// DB(totalTimer);
}

class AsteroidDetector {public:

VS testIDs;
VVF testFeatures;
VVF trainFeatures;

VC<unsigned short> preprocess(VI &v) {
	VC<unsigned short> rv(v.S);
	int MAX_VALUE = (1<<16)-1;
	REP(i, v.S) rv[i] = (unsigned short)max(0, min(MAX_VALUE, v[i]));
	return rv;
}

double trainTimer = 0;
double testTimer = 0;

int testCall = 0;
int trainCall = 0;

int trainingData(int width, int height, VI imageData_1, VS header_1, VD wcs_1, VI imageData_2, VS header_2, VD wcs_2, VI imageData_3, VS header_3, VD wcs_3, VI imageData_4, VS header_4, VD wcs_4, VS detections) {
	double xtime = getTime();

	VC<VC<unsigned short>> data(4);
	data[0] = preprocess(imageData_1);
	data[1] = preprocess(imageData_2);
	data[2] = preprocess(imageData_3);
	data[3] = preprocess(imageData_4);
	
	VVD wcs(4);
	wcs[0] = wcs_1;
	wcs[1] = wcs_2;
	wcs[2] = wcs_3;
	wcs[3] = wcs_4;
	
	VVS headers(4);
	headers[0] = header_1;
	headers[1] = header_2;
	headers[2] = header_3;
	headers[3] = header_4;
	
	solve(data, wcs, width, detections, headers, true);
	
	VC<VC<PDD>> correct = convertDetections(detections);
	
	VI sampleGood(xFeatures.S);
	VI sampleNeo(xFeatures.S);
	
	int total = xFeatures.S;
	int totalBad = 0;
	
	REP(i, xFeatures.S) {
		REP(j, correct.S) if (resultsDistance(xFeatures[i], correct[j]) < TRAIN_CORRECT_DISTANCE) {
			sampleGood[i] = true;
			sampleNeo[i] |= splt(detections[j * 4], ' ')[7] == "1";
		}
		xFeatures[i][1] = sampleGood[i];
		xFeatures[i][2] = sampleNeo[i];
		totalBad += !sampleGood[i];
	}
	
	REP(i, xFeatures.S) {
		if (sampleGood[i] || rng.nextDouble() < targetSamples / totalBad)
			trainFeatures.PB(xFeatures[i]);
	}
	
	trainCall++;
	
	trainTimer += getTime() - xtime;
	
	return 0;
}

int testingData(string imageID, int width, int height, VI imageData_1, VS header_1, VD wcs_1, VI imageData_2, VS header_2, VD wcs_2, VI imageData_3, VS header_3, VD wcs_3, VI imageData_4, VS header_4, VD wcs_4) {
	double xtime = getTime();

	testIDs.PB(imageID);
	
	VC<VC<unsigned short>> data(4);
	data[0] = preprocess(imageData_1);
	data[1] = preprocess(imageData_2);
	data[2] = preprocess(imageData_3);
	data[3] = preprocess(imageData_4);
	
	VVD wcs(4);
	wcs[0] = wcs_1;
	wcs[1] = wcs_2;
	wcs[2] = wcs_3;
	wcs[3] = wcs_4;
	
	VVS headers(4);
	headers[0] = header_1;
	headers[1] = header_2;
	headers[2] = header_3;
	headers[3] = header_4;
	
	solve(data, wcs, width, VS(), headers, false);
	
	REP(i, xFeatures.S) {
		xFeatures[i][0] += testCall;
		testFeatures.PB(xFeatures[i]);
	}
	
	testCall++;
	
	testTimer += getTime() - xtime;
	
	return 0;
}

VS getAnswer() {
	double xtime;

	RandomForestConfig cfg;
	cfg.featuresIgnored = 11;
	cfg.bagSize = 1.5;
	cfg.randomFeatures = {1, 2, 3, 4, 6, 8};
	cfg.randomPositions = {4};
	cfg.maxNodeSize = 15;
	
	
	VF trainResultsDetection;
	VF trainResultsNeo;
	
	for (VF v : trainFeatures) {
		trainResultsDetection.PB(v[1]);
		trainResultsNeo.PB(v[2]);
	}
	
	DB(trainFeatures.S);
	DB(testFeatures.S);
	
	xtime = getTime();
	RandomForest RFDetection;
	RFDetection.train(trainFeatures, trainResultsDetection, cfg, TREES_NO);
	double rftrainTimer = getTime() - xtime;
	
	DB(RFDetection.countTotalNodes());
	
	xtime = getTime();
	VC < pair < float, int> > vp;
	REP(i, testFeatures.S) {
		float v0 = RFDetection.estimate(testFeatures[i]);
		vp.PB(MP(v0, i));
	}
	double rfestimateTimer = getTime() - xtime;
	
	sort(vp.rbegin(), vp.rend());
	VS rv;
	
	VC<unordered_map<int,VI>> choice(100);
	
	VI samplesUsed;
	
	xtime = getTime();
	REP(i, vp.S) {
		if (rv.S >= MAX_RESULTS) break;
		int sampleID = vp[i].Y;
		int testID = ((int)testFeatures[sampleID][0]) % 100;
		int px = (int)(testFeatures[sampleID][3] * 10);
		int py = (int)(testFeatures[sampleID][4] * 10);
		
		int repeats = 0;
		FOR(gx, px - 1, px + 2) FOR(gy, py - 1, py + 2) if (choice[testID].count(gx*10000+gy)) for (int x : choice[testID][gx*10000+gy]) 
			if (resultsDistance(testFeatures[sampleID], testFeatures[x]) < MIN_DISTANCE) {
				repeats++;
				if (repeats >= 2) break;
			}
		if (repeats >= 2) continue;
		choice[testID][px*10000+py].PB(sampleID);
		samplesUsed.PB(sampleID);
		string s = testIDs[testID];
		FOR(j, 3, 11) {
			char str[100];
			sprintf(str, " %.10f", testFeatures[sampleID][j]);
			s += str;
		}
		rv.PB(s);
	}
	double constructTimer = getTime() - xtime;
	
	
	xtime = getTime();
	cfg.maxNodeSize = 50;
	RandomForest RFNeo;
	RFNeo.train(trainFeatures, trainResultsNeo, cfg, TREESNEO_NO);
	double neotrainTimer = getTime() - xtime;
	
	DB(RFNeo.countTotalNodes());
	
	xtime = getTime();
	VC < pair < float, int> > vn;
	REP(i, samplesUsed.S) {
		float v0 = RFNeo.estimate(testFeatures[samplesUsed[i]]);
		vn.PB(MP(v0, samplesUsed[i]));
	}
	double neoestimateTimer = getTime() - xtime;
	
	sort(vn.rbegin(), vn.rend());
	
	unordered_set<int> neos;
	REP(i, min(NEOS_USED, (int)samplesUsed.S)) neos.insert(vn[i].Y);
	REP(i, rv.S) rv[i] += neos.count(samplesUsed[i]) ? " 1" : " 0";
	
	DB(trainTimer);
	DB(testTimer);
	DB(rftrainTimer);
	DB(rfestimateTimer);
	DB(neotrainTimer);
	DB(neoestimateTimer);
	DB(constructTimer);
	
	return rv;
}

};
