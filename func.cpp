#include <map>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <glpk.h>
#include <algorithm>
#include <boost/algorithm/string.hpp>

#include "func.h"

using namespace std;
using namespace boost;

void readbounds(map<int,double> *lbs,map<int,double> *ubs,string boundsfile) {
    ifstream infile(boundsfile.c_str(),ios::in);
    char str[1024];
    string words = "";
	if (infile.is_open()) {
		while(!infile.eof()) {
			words = "";
			infile.getline(str,1024);
			if(strlen(str) < 1) {
				continue;
			}

			words += str;
			vector<string> x;
			split(x,words,is_any_of(" \t"));
			/*for(vector<string>::iterator viter=x.begin();viter!=x.end();viter++) {
			  cout<<(*viter)<<"\n";
			  }*/
			int ifir = atoi(x[0].c_str());
			double fsec = atof(x[1].c_str());
			double fthd = atof(x[2].c_str());
			//cout<<ifir<<"\t"<<fsec<<endl;
			if(ifir>0) {
				(*ubs)[ifir] = fsec;
				(*lbs)[ifir] = fthd;
			}
		}
	}
	else {
		cout<<"Can't open boundsfile from : "<<boundsfile<<endl;
		exit(-1);
	}
}
void read_mut(vector<int> *ind,vector<double> *mut,string mut_inputfile) {
    ifstream infile(mut_inputfile.c_str(),ios::in);
    char str[1024];
    string words = "";
	if (infile.is_open()) {
		while(!infile.eof()) {
			words = "";
			infile.getline(str,1024);
			if(strlen(str) < 1) {
				continue;
			}

			words += str;
			vector<string> x;
			split(x,words,is_any_of(" \t"));
			/*for(vector<string>::iterator viter=x.begin();viter!=x.end();viter++) {
			  cout<<(*viter)<<"\n";
			  }*/
			int ifir = atoi(x[0].c_str());
			double fsec = atof(x[1].c_str());
			//cout<<ifir<<"\t"<<fsec<<endl;
			if(ifir>0) {
				(*ind).push_back(ifir);
				(*mut).push_back(fsec);
			}
		}
	}
	else {
		cout<<"Can't open boundsfile from : "<<mut_inputfile<<endl;
		exit(-1);
	}
}

void get_constraints(glp_prob *lp,string confile,map<int,int> tp) {
	ifstream in(confile.c_str(),ios::in);
	//cout<<"Constraint file is "<<confile<<endl;
	char str[1024];
	string words = "";
	if (in.is_open()) {
		while(!in.eof()) {
			words = "";
			in.getline(str,1024);
			//cout<<strlen(str)<<endl;
			if(strlen(str) < 1) {
				continue;
			}
			words += str;
			vector<string> x;
			split(x,words,is_any_of(" \t"));
			for(vector<string>::iterator viter=x.begin();viter!=x.end();viter++) {
				cout<<(*viter)<<"\t";
			}
			cout<<"\n";
			int ind = atoi(x[0].c_str());
			int type = atoi(x[1].c_str());
			double lb = atof(x[2].c_str());
			double ub = atof(x[3].c_str());
			glp_set_col_bnds(lp,ind,tp[type],lb,ub);
		}
	} else {
		cout<<"Can't open constraint file from : "<<confile<<endl;
		exit(-1);
	}
}

void read_target(string confile,vector<double> *dt,int col) {
    ifstream in(confile.c_str(),ios::in);
    //cout<<"Constraint file is "<<confile<<endl;
    char str[10240];
    string words = "";
    int h_col = 0;
    cout<<"row is "<<col<<"\n";
	if (in.is_open()) {
		while(!in.eof()) {
			words = "";
			in.getline(str,10240);
			//cout<<strlen(str)<<endl;
			if(strlen(str) < 1) {
				continue;
			}
			h_col++;
			if(h_col>col) {
				cout<<"Reading target file error!!\n"<<endl;
			}
			(*dt)[h_col] = atof(str);
		}
		cout<<h_col<<" columns have been loaded for targets. "<<"\n";
	} else {
		cout<<"Can't open target file from : "<<confile<<endl;
		exit(-1);
	}
}

