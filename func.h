#include <map>
#include <cmath>
#include <string>

using namespace std;
extern int FIX;
extern int SAM;
extern int BIO;
extern int ROW;
extern int LOG;
extern int SCH;
extern int COST_TYPE;
extern int ADAPT;
extern double RANGE;
extern double COSTCOEF;
extern double BIO_MAG;
extern double LOWCUT;
extern double MaxFlux;
extern double MUT;
extern double DEL;
extern string MUTFILE;
extern string MUTFILEOBJ;
extern string DELFILE;
extern string DELFILEOBJ;
extern string VARIFILE;
extern string MUT_INPUTFILE;

void readbounds(map<int,double> *lbs,map<int,double> *ubs,string bounds_file);
void read_mut(vector<int> *,vector<double> *,string mutfile);
void get_constraints(glp_prob *lp,string confile,map<int,int>);
void read_target(string confile,vector<double> *dt,int row);
int hcsa(glp_prob* lp,map<int,double> lbs,map<int,double> ubs,vector<double> * dt, vector<double> *tar,fstream* out,string m,double mxx);
int efva(glp_prob* lp,vector<double> * dt, vector<double> *tar,string m);
int test_flux(glp_prob *lp,vector<double> flux,vector<double> * dt, vector<double> *tar,string file);
