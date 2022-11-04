/* -----------------------------------------------------------------
 * File:    morphing.cpp
 * Created: 2015-09-25
 * -----------------------------------------------------------------
 *
 *
 *
 * ---------------------------------------------------------------*/

#include "morphing.h"
#include <cassert>

using namespace std;

Vec2f operator+(const Vec2f &a, const Vec2f &b) {
  // --------- HANDOUT  PS05 ------------------------------
  // Return the vector sum of a an b
  return Vec2f(a.x + b.x, a.y + b.y); 
}

Vec2f operator-(const Vec2f &a, const Vec2f &b) {
  // --------- HANDOUT  PS05 ------------------------------
  // Return a-b
  return Vec2f(a.x - b.x, a.y - b.y); 
}

Vec2f operator*(const Vec2f &a, float f) {
  // --------- HANDOUT  PS05 ------------------------------
  // Return a*f
  return Vec2f(a.x * f, a.y * f);
}

Vec2f operator/(const Vec2f &a, float f) {
  // --------- HANDOUT  PS05 ------------------------------
  // Return a/f
  return Vec2f(a.x / f, a.y / f); 
}

float dot(const Vec2f &a, const Vec2f &b) {
  // --------- HANDOUT  PS05 ------------------------------
  // Return the dot product of a and b.
  return a.x * b.x + a.y * b.y; 
}

float length(const Vec2f &a) {
  // --------- HANDOUT  PS05 ------------------------------
  // Return the length of a.
  return sqrt(pow(a.x, 2) + pow(a.y, 2));
}

Vec2f perpendicular(const Vec2f &a) {
  // --------- HANDOUT  PS05 ------------------------------
  // Return a vector that is perpendicular to a.
  // Either direction is fine.
  return Vec2f(-a.y, a.x);
}

// The Segment constructor takes in 2 points P(x1,y1) and Q(x2,y2) corresponding
// to the ends of a segment and initialize the local reference frame e1,e2.
Segment::Segment(Vec2f P_, Vec2f Q_) : P(P_), Q(Q_) {
  // // --------- HANDOUT  PS05 ------------------------------
  // // The initializer list above ": P(P_), Q(Q_)" already copies P_
  // // and Q_, so you don't have to do it in the body of the constructor.
  // You should:
  // * Initialize the local frame e1,e2 (see header file)
  Vec2f PQ = Q - P; 
  lPQ = length(PQ); 
  e1 = PQ / lPQ; 
  e2 = perpendicular(PQ) / lPQ;   
}

Vec2f Segment::XtoUV(Vec2f X) const {
  // --------- HANDOUT  PS05 ------------------------------
  // Compute the u,v coordinates of a point X with
  // respect to the local frame of the segment.
  // e2 ^
  //    |
  // v  +  * X
  //    | /
  //    |/
  //    *--+------>-----*
  //    P  u     e1     Q
  //                    u=1
  //
  // * Be careful with the different normalization for u and v 
  float u = dot(X - P, Q - P)/ pow(lPQ, 2); 
  float v =  dot(X - P, perpendicular(Q - P))/lPQ; 
  return Vec2f(u, v); 
}

Vec2f Segment::UVtoX(Vec2f uv) const {
  // --------- HANDOUT  PS05 ------------------------------
  // compute the (x, y) position of a point given by the (u,v)
  // location relative to this segment.
  // * Be careful with the different normalization for u and v
  Vec2f preX = e1 * (lPQ * uv.x); 
  Vec2f X = (preX + (e2 * uv.y)) + P;
  return X;
}

float Segment::distance(Vec2f X) const {
  // --------- HANDOUT  PS05 ------------------------------
  // Implement distance from a point X(x,y) to the segment. Remember the 3
  // cases from class.
  Vec2f PQ = Q - P; 
  float squared_PQ = pow(length(PQ), 2); 
  if(squared_PQ == 0){
    return length(Q - X); 
  }

  float projection_p = dot(X - P, PQ)/squared_PQ; 

  if(projection_p < 0.0){
    return length(X - P); 
  }else if (projection_p > 1.0){
    return length(X - Q); 
  }

  Vec2f projection = P + (PQ * projection_p); 
  return length(projection - X);
}

Image warpBy1(const Image &im, const Segment &segBefore,
              const Segment &segAfter) {
  // --------- HANDOUT  PS05 ------------------------------
  // Warp an entire image according to a pair of segments.
  Image output(im.width(), im.height(), im.channels()); 
  for(int x = 0; x < im.width(); x++){
    for(int y = 0; y < im.height(); y++){
      Vec2f uv = segAfter.XtoUV(Vec2f(x, y)); 
      Vec2f X = segBefore.UVtoX(uv);
      for(int z = 0; z < im.channels(); z++){
        output(x,y,z) = interpolateLin(im, X.x, X.y, z, true);
      }
    }
  }
  return output;
}

float Segment::weight(Vec2f X, float a, float b, float p) const {
  // --------- HANDOUT  PS05 ------------------------------
  // compute the weight of a segment to a point X(x,y) given the weight
  // parameters a,b, and p (see paper for details).
  float numerator = pow(lPQ, p); 
  float denominator = a + distance(X); 
  return pow(numerator/denominator,b); 
}

Image warp(const Image &im, const vector<Segment> &src_segs,
           const vector<Segment> &dst_segs, float a, float b, float p) {
  // --------- HANDOUT  PS05 ------------------------------
  // Warp an image according to a vector of before and after segments using
  // segment weighting
  Image output(im.width(), im.height(), im.channels()); 
  for(int x = 0; x < im.width(); x++){
    for(int y = 0; y < im.height(); y++){
      for(int z = 0; z < im.channels(); z++){
        Vec2f DSUM = Vec2f(0,0); 
        Vec2f X = Vec2f(x,y); 
        float weightsum = 0; 
        for(int i = 0; i < dst_segs.size(); i++){
          Segment dst = dst_segs.at(i); 
          Segment src = src_segs.at(i); 

          Vec2f uv = dst.XtoUV(X); 
          Vec2f Xprime_i = src.UVtoX(uv); 

          Vec2f dis = Xprime_i - X;
          float weight = dst.weight(X, a, b, p); 
          DSUM = DSUM + (dis * weight); 
          weightsum += weight; 
        }
        Vec2f X_prime = X + (DSUM / weightsum);
        output(x,y,z) = interpolateLin(im, X_prime.x, X_prime.y, z, true);
      }
    }
  }
  return output;
}


vector<Image> morph(const Image &im_before, const Image &im_after,
                    const vector<Segment> &segs_before,
                    const vector<Segment> &segs_after, int N, float a, float b,
                    float p) {
  // --------- HANDOUT  PS05 ------------------------------
  // return a vector of N+2 images: the two inputs plus N images that morphs
  // between im_before and im_after for the corresponding segments. im_before
  // should be the first image, im_after the last.
  vector<Image> output; 
  output.push_back(im_before); 
  vector<float> interpolateFactors; 
  for(int n = 0; n < N; n++){
    if(n == 0){
      interpolateFactors.push_back(1.0/(N+1));
    }else{ 
      interpolateFactors.push_back(interpolateFactors.at(n - 1) + 1.0 / (N + 1));
    }
  }
  for(int n = 0; n < interpolateFactors.size(); n++){
    vector<Segment> segments; 
    for(int seg = 0; seg < segs_before.size(); seg++){
      Vec2f P1 = segs_before.at(seg).getP();
      Vec2f Q1 = segs_before.at(seg).getQ(); 
      Vec2f P2 = segs_after.at(seg).getP(); 
      Vec2f Q2 = segs_after.at(seg).getQ(); 

      Vec2f P = P1 + (P2 - P1) * interpolateFactors.at(n); 
      Vec2f Q = Q1 + (Q2 - Q1) * interpolateFactors.at(n); 

      Segment PQ = Segment(P, Q); 
      segments.push_back(PQ); 
    }
    Image before = warp(im_before, segs_before, segments, a,b,p); 
    Image after = warp(im_after, segs_after, segments, a, b, p); 
    Image morphed = before*(1 - interpolateFactors.at(n)) + after * interpolateFactors.at(n); 
    output.push_back(morphed); 
  }
  output.push_back(im_after); 
  return output;
}
