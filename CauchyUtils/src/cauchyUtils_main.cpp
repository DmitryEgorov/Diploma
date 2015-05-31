///////////////////////////////
// XDemo_main.cc
// C interaction file

#include "CauchyUtils.h"
#include "R.h" // R functions
#include "Rmath.h" // R math
#include <iostream>
#include <utility>
#include <vector>

using namespace std;

// Functions Passed to C++ from R must be passed in C extern format
// All variables are passed to C by reference (pointers);
// All output of functions is "void" (adjustments made via reference change)
extern "C" { 

  void testExact (double *location, double *scale, int *size, int *iterations_num) {
    Rprintf("\n\n\ntestExact location=%.2lf, scale=%.2lf, n=%d, it=%d\n\n", *location, *scale, *size, *iterations_num);
    int n = *size, it = *iterations_num;
    double *V = (double *)malloc(n * it * sizeof(double));
    for (int i = 0; i < n * it; i++)
      V[i] = mrcauchy(*location, *scale);
    for (int i = 0; i < LOCATION_METHODS; i++) {
      if (i == 1)
        continue;
      for (int j = 0; j < SCALE_METHODS; j++) {
        if (j == 1 && n > 100)
          continue;
        vector<double> R;
        R.resize(it);
        init_timer();
        for (int k = 0; k < it; k++) {
          double estLocation = cauchyEstimateLocation(V + n * k, n, i);
          double estScale = cauchyEstimateScale(V + n * k, n, j, estLocation);
          R[k] = sqr(estLocation - *location) + sqr(log(estScale) - log(*scale));
        }
        ans_t ans = generate_ans(R);
        Rprintf("location_t=%d, scale_t=%d, ans={%.5lf/%.5lf/%.3lfs} &\n", i, j, ans.first.first, ans.first.second, ans.second);
      }    
    }
    vector<double> R;
    R.resize(it);
    ans_t ans;

    init_timer();
    for (int k = 0; k < it; k++) {                                         
      double estLocation, estScale;
      cauchyEstimateNagy(V + n * k, n, &estLocation, &estScale);
      R[k] = sqr(estLocation - *location) + sqr(log(estScale) - log(*scale));
    }
    ans = generate_ans(R);
    Rprintf("Nagy, ans={%.5lf/%.5lf/%.5lfs}\n", ans.first.first, ans.first.second, ans.second);

    init_timer();
    for (int k = 0; k < it; k++) {                                         
      double estLocation, estScale;
      cauchyEstimateMaxLikelyhood(V + n * k, n, &estLocation, &estScale);
      R[k] = sqr(estLocation - *location) + sqr(log(estScale) - log(*scale));
    }
    ans = generate_ans(R);
    Rprintf("MLE, ans={%.5lf/%.5lf/%.5lfs}\n", ans.first.first, ans.first.second, ans.second);

    free(V);
  }

}
