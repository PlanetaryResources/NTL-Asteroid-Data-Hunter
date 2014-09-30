#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>
#include <sstream>
#include <map>
#include <algorithm>
#include <memory>
#include <sys/time.h>
using namespace std;
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
}
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
 return Mat2d(d/det, -b/det,
      -c/det, a/det);
    }
    Mat2d transpose() const {
 return Mat2d(a,c,b,d);
    }
    Vec2 operator*(Vec2 x) const {
 return Vec2(a*x[0] + b*x[1],
      c*x[0] + d*x[1]);
    }
    Mat2d operator*(Mat2d m) const {
 Vec2 col0 = (*this) * m.col(0);
 Vec2 col1 = (*this) * m.col(1);
 return Mat2d(col0[0], col1[0],
       col0[1], col1[1]);
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
    Mat2d m(x[0]*y[0], x[0]*y[1],
     x[1]*y[0], x[1]*y[1]);
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
vector<Detection> detectBright(Series* series, FullImageSet Z) {
    FullImageSet Z2;
    for (int f=0; f<FRAMES; f++)
        for (int y=0; y<H; y++)
            for (int x=0; x<W; x++) {
                float z = Z[f].at(y,x);
                z = max(0.f, z);
                Z2[f].at(y,x) = sqr(z) / (2 * sqr(STANDARD_QUAL_FACTOR));
            }
    const int anchor_step = 9;
    const float anchor_thresh = sqr(1.5);
    const int max_frame_anchors = 800*2;
    const int max_anchors = 4*max_frame_anchors;
    const int SPEED_RANGE = 7*2;
    const float TRACK_MERIT = .15;
    const int max_tracks = 600*2;
    const int min_track_distance = 12;
    const float MEAN_SPEED_THRESHOLD_SCORE = 10;
    const float delta_sigma = .23;
    const float speed_sigma = 70;
    vector<BrightAnchor> anchors;
    for (int f=0; f<FRAMES; f++) {
        vector<BrightAnchor> frame_anchors;
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
                if (best_pos[0]<15 || best_pos[1]<38)
                    continue;
                if (best_pos[0]>H-15 || best_pos[1]>W-25)
                    continue;
                if (best_val >= anchor_thresh)
                    frame_anchors.push_back(BrightAnchor { series, f, series->darkToStar[f](Vec2(best_pos[0], best_pos[1])), best_val });
            }
        }
        sort(frame_anchors.begin(), frame_anchors.end());
        if ((int)frame_anchors.size() > max_frame_anchors) {
            resize_if_larger(frame_anchors, max_frame_anchors);
        }
        for (BrightAnchor& anchor : frame_anchors)
            anchors.push_back(anchor);
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
                Vec2 dyx = { vi/3.0, vj/3.0 };
                float x[FRAMES];
                for (int f=0; f<FRAMES; f++)
                    x[f] = Z2[f].at(
                            (int)(offset_f[f][0] + dyx[0] * (f-fr) + .5),
                            (int)(offset_f[f][1] + dyx[1] * (f-fr) + .5)
                            );
                std::sort(x, x+FRAMES);
                float merit = 4*x[0] + 2*x[1] + x[2];
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
    vector<float> track_dy;
    vector<float> track_dx;
    for (Detection& d : tracks) {
        if (d.score >= MEAN_SPEED_THRESHOLD_SCORE) {
            track_dy.push_back(d.dyx[0]);
            track_dx.push_back(d.dyx[1]);
        }
    }
    int p = track_dy.size() > 0 ? (int)track_dy.size()/2 : -1;
    if (p >= 0) {
        std::nth_element(track_dy.begin(), track_dy.begin() + p, track_dy.end());
        std::nth_element(track_dx.begin(), track_dx.begin() + p, track_dx.end());
    }
    cerr << (int)tracks.size() << " final tracks" << endl;
    Vec2 mean_dyx = { p>=0 ? track_dy[p] : 0, p>=0 ? track_dx[p] : 0 };
    cerr << mean_dyx << " mean speed" << endl;
    for (Detection& d : tracks) {
        float g_base = 52*d.gain_histogram[0] + 8*d.gain_histogram[1] + d.gain_histogram[2];
        float g_series = 0;
        float delta2 = norm2(d.dyx - mean_dyx);
        float g_delta = -delta2 / (2 * sqr(delta_sigma));
        float g_speed = -delta2 / (2 * sqr(speed_sigma) * sqr(series->time_delta));
        float y = d.yx0[0];
        float x = d.yx0[1];
        float g_position = 0;
        float g_total = g_base + g_series + g_delta + g_speed + g_position;
        bool skip = false;
        if (y > 4096 - 50 || y < 10 ||
            x > 4096 - 40 || x < 10)
            skip = true;
        if (skip) {
     if (y<=0 || y>=4096 || x<=0 || x>=4110)
  d.score = -9999;
     else
  d.score = g_total/5 - 6;
        } else {
            d.score = g_total;
        }
    }
    sort(tracks.begin(), tracks.end());
    return tracks;
}
vector<string> finish(vector<Detection>& detections) {
    std::sort(detections.begin(), detections.end());
    resize_if_larger(detections, MAX_ALLOWED_DETECTIONS);
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
    cerr << s << "," << a << "," << b << endl;
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
    cerr << "ADC offset: " << offset << endl;
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
        Point center = series->starToDark[f](Point { H/2, W/2 });
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
        cerr << "about to median" << endl;
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
        cerr << "about to expand" << endl;
        vector<float> glow_half_expanded(ni*W);
        for (int i=0; i<ni; i++) {
            for (int x=0; x<W; x++) {
                glow_half_expanded[i*W + x] =
                    arrayinterp1d_cubic_nearest(&glow[i*nj + 0], nj, 1, (x-center[1])/S + j_mid - 0.5f * w);
            }
        }
        cerr << "about to expand" << endl;
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
        cerr << "frame_bias = " << frame_bias << endl;
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
void computeEverything(Series* series, FullImageSet frames, FullImageSet z) {
    cerr << "computeEverything started" << endl;
    computeUnglow(series, frames);
    FullImage star, dark;
    FullImageSet mov;
    FullImageSet tstar;
    computeLayers(series, frames, star, dark, tstar, mov);
    computeFiltered(series, tstar, dark, mov, z);
    cerr << "computeEverything done" << endl;
}
vector<Detection> detectEverything(Series* series, FullImageSet frames) {
    FullImageSet z;
    computeEverything(series, frames, z);
    return detectBright(series, z);
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
 Series* series = new Series();
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
static double t0;
double getTime() {
    timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + 1e-6 * tv.tv_usec;
}
static bool initialized = false;
static vector<Detection> detections;
struct AsteroidDetector {
    static int trainingData(int, int, vector <int>, vector <string>, vector <double>, vector <int>, vector <string>, vector <double>, vector <int>, vector <string>, vector <double>, vector <int>, vector <string>, vector <double>, vector <string>) {
 t0 = getTime();
 return 1;
    }
    static int testingData(string id, int width, int height,
     vector<int> data0, vector<string> head0, vector<double> wcs0,
     vector<int> data1, vector<string> head1, vector<double> wcs1,
     vector<int> data2, vector<string> head2, vector<double> wcs2,
     vector<int> data3, vector<string> head3, vector<double> wcs3) {
 cerr << "AsteroidDetector: time = " << getTime() - t0 << endl;
 FullImageSet frames;
 Series* series = new Series();
 series->id = id;
 if (width != W)
     cerr << "Invalid width" << endl;
 if (height != H)
     cerr << "Invalid height" << endl;
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
 cerr << "frontend: Got series " << series->id << endl;
 for (Detection& d : detectEverything(series, frames))
     detections.push_back(d);
 return 0;
    }
    static vector<string> getAnswer() {
 cerr << "AsteroidDetector: time = " << getTime() - t0 << endl;
 vector<string> r = finish(detections);
 cerr << "frontend: Finished, time = " << getTime()-t0 << endl;
 return r;
    }
};

