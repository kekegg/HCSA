/***********************************************************************
//      AUTHOR:  Zhengwei Xie (Mr), <xiezhengwei@hsc.pku.edu.cn>
//      COMPANY:  Department of Pharmacology, School of Basic Medical
//				Sciences, Peking University, 
//				38 Xueyuan Lu, Haidian District, 
//				Beijing, 100191, China
//      VERSION:  1.0
//      PAPER :  
//		"Genome-scale fluxes were predicted under the guide of enzyme 
//		abundance using a novel Hy-per-Cube Shrink Algorithm"
//		Zhengwei Xie, Tianyu Zhang and Qi Ouyang
***********************************************************************/

#include <map>
#include <cmath>
#include <ctime>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <glpk.h>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <iomanip>

#include "func.h"

using namespace std;
using namespace boost;
#define INFI 1e20
#define ERR 1e-20
#define UP 1
#define LOW 0

int BIO=1; 
int ROW = -1;
double LOWCUT = 1e-1;
double MaxFlux=0.0;

map<int,int> tp;
void tp_init() {
    tp[0] = GLP_FX;
    tp[1] = GLP_DB;
    tp[2] = GLP_LO;
    tp[3] = GLP_UP;
}

int
main (int argc, char **argv)
{	
	BIO = 0;
	string fname="";
	string confile="";
	string boundsfile="";
	string target_file="";
	string fvafile="";
	string flux_na="";
	string trace_na="";
	string method ="ratio";
	string help="";
	double mxx = 0;
	tp_init();
	for(int i=1; i<argc; i++) {
		if(strcmp(argv[i],"-i") == 0) {
			i++;
			cout<<"Model file : "<<argv[i]<<endl;
			fname += argv[i];
		} else if(strcmp(argv[i],"-c") == 0) {
			i++;
			cout<<"Constraint file (optional) : "<<argv[i]<<endl;
			confile += argv[i];
		} else if(strcmp(argv[i],"-l") == 0) {
			i++;
			cout<<"Inputted boudary file, optional: "<<argv[i]<<endl;
			boundsfile+= argv[i];
		} else if(strcmp(argv[i],"-b") == 0) {
			i++;
			cout<<"Enzyme abundance file (EV) : "<<argv[i]<<endl;
			target_file = argv[i];
		} else if(strcmp(argv[i],"-o") == 0) {
			i++;
			cout<<"Fluxes output file (FV) : "<<argv[i]<<endl;
			flux_na += argv[i];
		} else if(strcmp(argv[i],"--row") == 0) {
			i++;
			cout<<"Number of row in EV file : "<<argv[i]<<endl;
			ROW = atoi(argv[i]);
		}  else if(strcmp(argv[i],"--max") == 0) {
			i++;
			cout<<"Maximum value in EV (optional) : "<<argv[i]<<endl;
			mxx = atof(argv[i]);
		} else if(strcmp(argv[i],"--bio") == 0) {
			i++;
			cout<<"Column index of biomass flux: "<<argv[i]<<endl;
			BIO = atoi(argv[i]);
		} else if(strcmp(argv[i],"--fva") == 0) {
			i++;
			cout<<"Perform FVA and output the file in (optional) : "<<argv[i]<<endl;
			fvafile += argv[i];
		} else if(strcmp(argv[i],"-t") == 0) {
			i++;
			cout<<"Print trace and print to (optional) : "<<argv[i]<<endl;
			trace_na += argv[i];
		}  else if(strcmp(argv[i],"-h") == 0 || strcmp(argv[i],"--help") == 0 ) {
			help += argv[i];
		}  else {
			cout<<"Unrecognized option : "<<argv[i]<<endl;
			exit(0);
		}
	}
	if(fname.length() < 1 || help.length() > 0) {
		cout<<"		Usage : ./hcsa -i toy.mps --row 3 -b EV.txt -o flux.txt -t trace.txt\n";
		cout<<"		-i model.mps\n";
		cout<<"		-c constraint.txt: Extra Constraint file (optional)\n";
		cout<<"		-l boundary.txt: FVA file, Inputted boudary file, optional\n";
		cout<<"		-b EV.txt: Enzyme abundance file (EV)\n";
		cout<<"		-o flux.txt: Fluxes output file (FV) \n";
		cout<<"		--row N: Number of row in EV file\n";
		cout<<"		--max val: Maximum value in EV (optional) \n";
		cout<<"		--bio N: Column index of biomass flux\n";
		cout<<"		--fva variable.txt: Perform FVA and output the file in (optional) \n";
		cout<<"		-t tracefile.txt: print trace and print to (optional) \n";
		cout<<"		-h or --help : Print this help.\n";
		exit(-1);
	}
	if(ROW <= 0) {
		cout<<"Please specify the row number of target file\n";
		cout<<"     --row N \n";
		exit(-2);
	}
	//////////////////////////////////////////////////////////////
	//Loading input files
	glp_prob *lp=NULL;
	lp=glp_create_prob();
	glp_read_mps(lp,GLP_MPS_FILE,NULL,fname.c_str());

	cout<<"Successfully read model from "<<fname<<endl;
	if(strlen(confile.c_str()) > 0) {
		cout<<"Constraint file is "<<confile<<endl;
		get_constraints(lp,confile,tp);
	}
	fstream out2;
	if(flux_na.length() > 0) {
		out2.open(flux_na.c_str(),ios::out);
		cout<<"Output file is :"<<flux_na<<endl;
	} else {
		out2.open("internal_points",ios::out);
		cout<<"Output file is : internal_points. "<<endl;
	}	
	fstream out3;
	if(trace_na.length() > 0) {
		out3.open(trace_na.c_str(),ios::out);
	}    
	//////////////////////////////////////////////////////////////
	//FVA
	int N = glp_get_num_cols(lp);
	map<int,double> lbs;
	map<int,double> ubs;
	if(boundsfile.length()>0) {
		cout<<"Got bounds from file : "<<boundsfile<<endl;
		readbounds(&lbs,&ubs,boundsfile);
	} else {
		for(int i=1;i<=N;i++) {
			glp_set_obj_coef(lp,i,0.);
		}
		cout<<"Generating column bounds . . .\n";
		for(int i=1;i<=N;i++) {
			glp_set_obj_coef(lp,i,1);
			glp_set_obj_dir(lp,GLP_MAX);
			glp_simplex(lp,NULL);
			ubs[i] = glp_get_col_prim(lp,i);
			glp_set_obj_dir(lp,GLP_MIN);
			glp_simplex(lp,NULL);
			lbs[i] = glp_get_col_prim(lp,i);
			glp_set_obj_coef(lp,i,0);
		}
	}

	//////////////////////////////////////////////////////////////
	//
	cout<<"Initializing directions . . .\n";
	srand(time(0));

	//////////////////////////////////////////////////////////////
	vector<double> target;// = new double[row+1];
	for(int i=1;i<=N+1;i++) {
		target.push_back(0.0);
	} 	
	read_target(target_file,&target,ROW);

	cout<<"Running HCSA ... "<<endl;
	vector<double> flux;
	for(int i=1;i<=N+1;i++) {
		flux.push_back(0.0);
	} 
	hcsa(lp,lbs,ubs,&flux,&target,&out3,method,mxx);
	cout<<"Printing out the flux distribution.\n";
	for(int i=1;i<=N;i++) {
		out2<<flux[i]<<"\n";
	}

	//////////////////////////////////////////////////////////////
	return 0;
}


