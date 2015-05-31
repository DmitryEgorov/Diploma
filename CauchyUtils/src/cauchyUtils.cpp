///////////////////////////////
// C++ Function File

#include "cauchyUtils.h"
#include "R.h" // R memory io
#include "Rmath.h" // R math functions
#include <algorithm>
#include <iostream>
#include <assert.h>
#include <time.h>

using namespace std;

static unsigned int I1=time(NULL)/*1234*/, I2=time(NULL)/*5678*/;  
static double start_time;
 
void set_seed(unsigned int i1, unsigned int i2) { 
  I1 = i1; I2 = i2; 
} 
 
void get_seed(unsigned int *i1, unsigned int *i2) { 
  *i1 = I1; *i2 = I2; 
} 
 
double unif_rand(void) { 
  I1= 36969*(I1 & 0177777) + (I1>>16); 
  I2= 18000*(I2 & 0177777) + (I2>>16); 
  return ((I1 << 16)^(I2 & 0177777)) * 2.328306437080797e-10; /* in [0,1) */ 
} 

double mrcauchy(double location, double scale) { 
 	return location + scale * tan(M_PI * unif_rand()); 
} 

double sqr (double a) {
  return a * a;
}

void init_timer () {
  start_time = clock() * 1.0 / CLOCKS_PER_SEC;
}

ans_t generate_ans (vector<double> &V) {
  ans_t res = mp(mp(0.0, 0.0), clock() * 1.0 / CLOCKS_PER_SEC - start_time);
  assert(V.size()%20 == 0);
  sort(V.begin(),V.end());
  double s = 0;
  for (int i = 0; i < (int)V.size(); i++)
    s += V[i];
  res.first = mp(sqrt(s / V.size()), sqrt(V[V.size() * 19 / 20]));
  return res;
}



//methods:
//0 - median;
//1 - mean;
//2 - truncated mean;
double cauchyEstimateLocation (double *V, int n, int method) {
  double *Vcopy = (double *)malloc(n * sizeof(double));
  memcpy(Vcopy, V, n * sizeof(double));
  sort(Vcopy, Vcopy + n);
  double res = 0.0;
  switch (method)
  {
    case 0:
    {
      if (n & 1)
        res = Vcopy[n >> 1];
      else
        res = (Vcopy[(n >> 1)] + Vcopy[(n >> 1) - 1]) / 2.0;
      break;
    }
    case 1:
    {
      for (int i = 0; i < n; i++)
        res += V[i];
      res /= n;
      break;
    }
    case 2:
    {
      int c = (int)(0.38 * n + 0.5);
      if ((c << 1) >= n)
        c = ((n - 1) >> 1);
      for (int i = c; i + c < n; i++)
        res += Vcopy[i];
      res /= (n - (c << 1));
      break;
    }
    default:
      assert(0 && "Unknown method for estimating location");
  }
  free(Vcopy);
  return res;
}

double cauchyLikelyhoodDiff(double *V, int n, double loc, double sc, double log_sc, int f, int x) {
  double res = 0;
  double sc2 = sc * sc;
  if (f == 0) {
    if (x == 0) {
      for (int i = 0; i < n; i++)
        res += (-sc2+sqr(loc-V[i]))/sqr(sc2+sqr(loc-V[i]));
    } else {
      for (int i = 0; i < n; i++)
        res -= (2*sc2*(-loc+V[i]))/sqr(sc2+sqr(loc-V[i]));
    }
  } else {
    if (x == 0) {
      for (int i = 0; i < n; i++)
        res += (2*sc2*(-loc+V[i]))/sqr(sc2+sqr(loc-V[i]));
    } else {
      for (int i = 0; i < n; i++)
        res += (2*sc2*sqr(-loc+V[i]))/sqr(sc2+sqr(loc-V[i]));
    }
  }
  return res;
}

double cauchyLikelyhoodEval(double *V, int n, double loc, double sc, int f) {
  double res = 0;
  sc*=sc;
  if (f == 0)
    for (int i = 0; i < n; i++)
      res += (V[i] - loc) / (sc + sqr(V[i] - loc));
  else {
    for (int i = 0; i < n; i++)
      res += 1.0 / (sc + sqr(V[i] - loc));
    res *= sc, res -= n / 2.0;
  }
  return res;
}

pair<double, double> cauchyOptimizeLikelyhoodScale (double *V, int n, double location) {
  double l = 1e100, r = 0.0;
  for (int i = 0; i < n; i++)
    l = min(l, fabs(location - V[i])), r = max(r, fabs(location - V[i]));
  double log_sc = (log(l) + log(r)) / 2.0;
  double sc = exp(log_sc);
  int it = 0;
  while (it < 30) {
    it ++;
    assert(it < 30);
    double val = cauchyLikelyhoodEval(V, n, location, sc, 1);
    double df = cauchyLikelyhoodDiff(V, n, location, sc, log_sc, 1, 1);
    if (fabs(df) < eps)
      break;
    val = -val / df;
    if (fabs(val) < eps)
      break;
    if (fabs(val) > 10.0)
      val = -log_sc;
    log_sc += val;
    sc = exp(log_sc);
  }
  return mp(sc, log_sc);
  //slow binary search version for debug
  /*for (int it = 0; it < 60; it++) {
    double m = (l + r) / 2.0;
    double m2 = sqr(m);
    double val = -n / 2.0;
    for (int j = 0; j < n; j++)
      val += m2 / (m2 + sqr(location - V[j]));
    (val <= 0.0) ? (l = m) : (r = m);
  }
  double t = (l + r) / 2.0;
  return mp(t, log(t));*/
}

double cauchyOptimizeLikelyhoodLocation (double *V, int n, double scale) {
  double l = 1e100, r = -1e100;
  for (int i = 0; i < n; i++)
    l = min(l, V[i]), r = max(r, V[i]);
  for (int it = 0; it < 60; it++) {
    double m = (l + r) / 2.0;
    double val = 0.0;
    for (int j = 0; j < n; j++)
      val += (V[j] - m) / (sqr(scale) + sqr(m - V[j]));
    (val >= 0.0) ? (l = m) : (r = m);
  }
  return (l + r) / 2.0;
}

//methods:
//0 - (Q_3 - Q_1) / 2;
//1 - Hodges-Lehmann slow
//2 - maximum likelyhood
//3 - Hodges-Lehmann fast
double cauchyEstimateScale (double *V, int n, int method, double location) {
  double *Vcopy = (double *)malloc(n * sizeof(double));
  memcpy(Vcopy, V, n * sizeof(double));
  double res = 0.0;
  switch (method)
  {
    case 0:
    {
      sort(Vcopy, Vcopy + n);
      res = Vcopy[(n >> 2) + (n >> 1)] - Vcopy[(n >> 2)];
      break;
    }
    case 1:
    {
      if (n > 1000)
        break;
      for (int i = 0; i < n; i++)
        Vcopy[i] -= location;
      int e = (n * (n + 1) / 2);
      double *Tmp = (double *)malloc(e * sizeof(double));
      e = 0;
      for (int i = 0; i < n; i++)
        for (int j = i; j < n; j++)
          Tmp[e++] = fabs(V[i] * V[j]);
      sort(Tmp, Tmp + e);
      (e & 1) ? (res = Tmp[e >> 1]) : (res = (Tmp[e >> 1] + Tmp[(e >> 1) - 1]) / 2.0);
      res = sqrt(res);
      free(Tmp);
      break;
    }
    case 2:
    {
      res = cauchyOptimizeLikelyhoodScale(V, n, location).first;
      break;
    }
    case 3:
    {
      for (int i = 0; i < n; i++)
        Vcopy[i] -= location, Vcopy[i] = fabs(Vcopy[i]);
      sort(Vcopy, Vcopy + n);
      double l = Vcopy[0] * Vcopy[0], r = Vcopy[n-1] * Vcopy[n-1];
      int it = 0;
      while (sqrt(r) - sqrt(l) > 1e-12 && it < 50) {
        it++;
        assert(it < 50);
        double m = (l + r) / 2;
        long long cnt = 0;
        int j = n - 1;
        for (int i = 0; i <= j; i++) {
          while (j >= 0 && Vcopy[i] * Vcopy[j] > m)
            j--;
          cnt += j - i + 1;
        }
        (cnt <= n * 1ll * (n + 1)/ 4 + 1) ? (l = m) : (r = m);    
      }
      res = (l + r) / 2;
      res = sqrt(res);
      break;
    }
                                                       
    default:
      assert(0 && "Unknown method for estimating scale");
  }
  free(Vcopy);
  return res;
}

//slow version of MLE for debug
/*void cauchyEstimateMaxLikelyhood (double *V, int n, double *location, double *scale) {
  int it = 0;
  *location = cauchyEstimateLocation(V, n, 2);
  *scale = cauchyEstimateScale(V, n, 2, *location);
  while (it < 10) {
    it++;
    assert(it < 10);
    double loc = *location, sc = *scale;
    sc = cauchyOptimizeLikelyhoodScale(V, n, loc).first;
    loc = cauchyOptimizeLikelyhoodLocation(V, n, sc);
    if (fabs(loc - *location) + fabs(sc - *scale) < eps)
      break;
    *location = loc, *scale = sc;
  }
}*/

void cauchyEstimateMaxLikelyhood (double *V, int n, double *location, double *scale) {
  int it = 0;
  double loc = cauchyEstimateLocation(V, n, 2);
  pair<double, double> t = cauchyOptimizeLikelyhoodScale(V, n, *location);
  double sc = t.first, log_sc = t.second;
  double dloc, dlog_sc, a, b, c, d, v1, v2, det;
  while (it < 10) {
    it++;
    assert(it < 10);
    v1 = -cauchyLikelyhoodEval(V, n, loc, sc, 0);
    v2 = -cauchyLikelyhoodEval(V, n, loc, sc, 1);
    a = cauchyLikelyhoodDiff(V, n, loc, sc, log_sc, 0, 0);
    b = cauchyLikelyhoodDiff(V, n, loc, sc, log_sc, 0, 1);
    c = cauchyLikelyhoodDiff(V, n, loc, sc, log_sc, 1, 0);
    d = cauchyLikelyhoodDiff(V, n, loc, sc, log_sc, 1, 1);            
    det = a * d - b * c;
    if (fabs(det) < eps)
      break;
    dloc = (d * v1 - b * v2) / det;
    dlog_sc = (- c * v1 + a * v2) / det;
    if (fabs(dloc) + fabs(dlog_sc) < eps)
      break;
    if (fabs(dloc) > 1.0)
      dloc /= fabs(dloc);
     if (fabs(dlog_sc) > 0.5)
      dlog_sc /= 2.0 * fabs(dlog_sc);
    loc += dloc, log_sc += dlog_sc, sc = exp(log_sc);
  }
  t = cauchyOptimizeLikelyhoodScale(V, n, loc);
  *scale = t.first, *location = loc;
}

void cauchyEstimateNagy (double *V, int n, double *location, double *scale) {
  int it = 0;
  double u, e0, e1, loc, sc;
  *location = 0.0;
  *scale = 1.0;
  while (it < 50) {
    it++;
    assert(it < 50);
    e0 = 0.0, e1 = 0.0;
    for (int i = 0; i < n; i++) {
      u = (V[i] - *location) / (*scale);
      e0 += 1.0 / (1.0 + u * u);
      e1 += u / (1.0 + u * u);
    }
    e0 /= n, e1 /= n;
    loc = *scale * e1 / e0, sc = sqrt(1.0 / e0 - 1.0);
    if (fabs(loc) + fabs(*scale * (sc - 1)) < eps)
      break;
    *location += loc, *scale *= sc;
  }
}
