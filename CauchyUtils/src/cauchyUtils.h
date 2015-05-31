#pragma once

#include<iostream>
#include<vector>
#include<utility>

using namespace std;

#define LOCATION_METHODS 3
#define SCALE_METHODS 4
#define eps 1e-12

#define mp make_pair
#define pb push_back

typedef pair <pair<double, double>, double> ans_t;

double mrcauchy (double location, double scale);
double sqr (double a);
void init_timer();
ans_t generate_ans (vector<double> &V);


double cauchyEstimateLocation (double *V, int n, int method);
double cauchyEstimateScale (double *V, int n, int method, double location);
void cauchyEstimateMaxLikelyhood (double *V, int n, double *location, double *scale);
void cauchyEstimateNagy (double *V, int n, double *location, double *scale);
