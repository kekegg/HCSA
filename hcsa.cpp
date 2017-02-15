#include <map>
#include <cmath>
#include <ctime>
#include <string>
#include <iostream>
#include <glpk.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics.h>
#include "func.h"

using namespace std;

typedef std::pair<double, int> Pair;

struct CmpPair
{
    bool operator()(const Pair& a, const Pair& b)
    { return a.first > b.first; }
};

void sortingPermutation(
    const std::vector<double>& values,
    std::vector<int>& permutation)
{
    std::vector<Pair> pairs;
    for (int i = 0; i < (int)values.size(); i++)
        pairs.push_back(Pair(values[i], i));

    std::sort(pairs.begin(), pairs.end(), CmpPair());

    typedef std::vector<Pair>::const_iterator I;
    for (I p = pairs.begin(); p != pairs.end(); ++p)
        permutation.push_back(p->second);
}

int hcsa(glp_prob* lp,map<int,double> lbs,map<int,double> ubs,vector<double> * dt, vector<double> * tar, fstream* out, string m, double mxx)
{
	glp_set_obj_coef(lp,1,0.);
	int N = glp_get_num_cols(lp);
	if(!(N>0)) {
		cout<<"Blank LP model, exit"<<endl;
		exit(0);
	}
	vector<double> tar_copy;
	tar_copy.push_back(0.0);
	for(int k=1;k<=ROW;k++) {
		tar_copy.push_back((*tar)[k]);
	}
	srand48(time(NULL));
	std::vector<int> permutation, perm;
	sortingPermutation((*tar), perm);
	vector<double>::iterator mx_iter = std::max_element(tar->begin(),tar->end());
	double mx = (*tar)[distance(tar->begin(),mx_iter)];
	if(mxx>0) {
		mx = mxx;
		cout<<"Using input maximum enzyme level : "<<mxx<<endl;
	} else {
		cout<<"Maximum enzyme level : "<<mx<<" at :"<<distance(tar->begin(),mx_iter)<<endl;
	}
	int dim = perm.size();
	for(int i=0;i<dim+1;i++) {
		permutation.push_back(0);
	}
	for(int i=0;i<dim;i++) {
		permutation[perm[i]] = i+1;
	}

	///////////////////////////////
	//Fake Objective Reaction
	glp_add_cols(lp,1);
	N = glp_get_num_cols(lp);
	glp_set_col_name(lp,N,"M");
	glp_set_col_bnds(lp,N,GLP_LO,0,0);

	srand(time(0));
	map<int,int> col_row;
	std::string s;
	std::stringstream ss;
	int ind[3]={0,0,N};
	double val[3]={0.,1.,0.};
	double *coef = new double[N];
	double *coefup = new double[N];
	double *coef_pre = new double[N];
	double *coefup_pre = new double[N];
	double *down = new double[N];
	double *up = new double[N];
	double *ran = new double[N];
	double *freq = new double[N];  //Frequence of acceptted changing
	double *freq_t = new double[N];  //Frequence of touched changing
	int *state= new int[N];
	double modf = 1.;

	for(int i=1;i<N;i++) {
		double range = ubs[i] - lbs[i];
		ran[i] = range; //Initialization
		freq[i] = 0;
		freq_t[i] = 0;
		double flag = 0;
		double r=0, r2=0;
		r = 1-permutation[i]/(dim+1.0);
		r2 = 1-r;

		if(range > 1e-10 && (*tar)[i]>=0 && (glp_get_col_type(lp,i) == GLP_DB || glp_get_col_type(lp,i) == GLP_LO) && lbs[i] >= -1e-12) {
			flag = 1;
			coef[i] = r;//(*tar)[i]/mx/1.0001;
			coef_pre[i] = r;
			coefup[i] = r2;
			coefup_pre[i] = r2;
			state[i] = 1;
		} else if(range > 1e-10 && (*tar)[i]>0 && (glp_get_col_type(lp,i) == GLP_DB || glp_get_col_type(lp,i) == GLP_UP) && ubs[i] <= 1e-12) {
			//negative flux bounds, flip coefficient
			flag = 1;
			coef[i] = r;//1-(*tar)[i]/mx/1.0001;
			coef_pre[i] = r;
			coefup[i] = r2;
			coefup_pre[i] = r2;
			state[i] = -1;
		} else if(range<1e-10) {
			flag = 0;
			coef[i] = -2;
			state[i] = 0;
		} else {
			flag = 0;
			coef[i] = -2;
			state[i] = -2;
		}
	}
	double obj_val = 0;// obj_pre=0;
	vector<double> rt;
	vector<int> rtind;
	map<int,int> lr;
	double ret;
	//(*out)<<"T step : "<<Tstep<<endl;
	for(int i=1;i<N;i++) {
		if(state[i] ==1 || state[i] == -1) {
			double range = ubs[i]-lbs[i];
			int r_n = glp_get_num_rows(lp);
			cout<<"Add row "<<r_n<<endl;
			glp_add_rows(lp,2);
			ss.str("");
			ss<<r_n+1;
			s = ss.str();
			lr[i] = r_n+1;
			glp_set_row_name(lp,r_n+1,s.c_str());
			col_row[i] = r_n+1; //Store the row number
			ind[1] = i;
			val[2] = -1*coef[i]*state[i]*range/modf;
			glp_set_mat_row(lp,r_n+1,2,ind,val);  //ind: column numbers; val: coeffiecients
			if(state[i] == 1) {
				glp_set_row_bnds(lp,r_n+1,GLP_LO,lbs[i],0);
			} else if(state[i] == -1) {
				glp_set_row_bnds(lp,r_n+1,GLP_UP,0,ubs[i]);
			}
			ss<<"_2";
			s = ss.str();
			glp_set_row_name(lp,r_n+2,s.c_str());
			ind[1] = i;
			val[2] = coefup[i]*state[i]*range/modf;
			glp_set_mat_row(lp,r_n+2,2,ind,val);
			if(state[i] == 1) {
				glp_set_row_bnds(lp,r_n+2,GLP_UP,0.,ubs[i]);
			} else if(state[i] == -1) {
				glp_set_row_bnds(lp,r_n+2,GLP_LO,lbs[i],0.);
			}
		}
	}
	glp_set_col_bnds(lp,N,GLP_LO,0,0);
	glp_set_obj_dir(lp,GLP_MAX);
	glp_set_obj_coef(lp,N,1);
	//Run the LP
	ret = glp_simplex(lp,NULL);
	//value of V
	obj_val = glp_get_obj_val(lp);
	//Get the boundaries

	if(ret == 0 && obj_val > 0) {
		for(int i=1;i<N;i++) {
			(*dt)[i] = glp_get_col_prim(lp,i);
			if(state[i] == 1) {
				down[i] = obj_val*coef[i]*ran[i]/modf+lbs[i]; 
				up[i] = ubs[i]-obj_val*coefup[i]*ran[i]/modf;
			} else if(state[i] == -1) {
				up[i] = obj_val*coef[i]*state[i]*ran[i]/modf+ubs[i]; 
				down[i] = lbs[i]-state[i]*obj_val*coefup[i]*ran[i]/modf;
			}
		}
		cout<<"Obj val: "<<"\t"<<obj_val<<endl;
		if(BIO) {
			cout<<"Growth rate: "<<"\t"<<glp_get_col_prim(lp,BIO)<<endl;
		}
		(*out)<<"UP"<<"\t"<<"Down\t"<<"Range\t"<<"Ratio\t"<<"DownCut\t"<<"UP\t"<<"flux\ttar\tstate\tfreq\tfreq_touched\tnRange\tnTar"<<endl;
		for(int i=1;i<N;i++) {
			double r=coef[i];//r2=coefup[i];
			(*out)<<ubs[i]<<"\t"<<lbs[i]<<"\t"<<ubs[i]-lbs[i]<<"\t"<<r<<
				"\t"<<down[i]<<"\t"<<up[i]<<"\t"<<(*dt)[i]<<"\t"<<(*tar)[i]<<"\t"<<
				state[i]<<"\t"<<freq[i]<<"\t"<<freq_t[i]<<"\t"<<ran[i]<<"\t"<<tar_copy[i]<<"\n";
		}
	} else {
		cout<<"obj_val: "<<obj_val<<"; NO FEASIBLE SOLUTION!!!!!!!"<<endl;
		for(int i=1;i<N;i++) {
			(*out)<<state[i]<<endl;
		}
	}

	glp_write_lp(lp,NULL,"model_final.lp");

	delete []coef;
	return 0;
}

