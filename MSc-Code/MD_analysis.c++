using namespace std;

#include<iostream>
#include<stdlib.h>
#include<cmath>
#include<sstream>
#include<fstream>
#include<string>
#include<iomanip>

const double PI  = 3.141592653589793238462;

long int k=0;
long int Kmax;
int Kstart;
const int Nmax=100;
const int Nbeadsmax=100;

int N,Nbeads;
int DeltaT=1;
double dt=0.01; 
double Lx,Ly,Lz;

double Rg[Nmax];
double Rg2[Nmax];
double Rg4[Nmax];
double Rgave;
double Rg2ave;
double Rg4ave;

int Ktot;
int Kfinal;
int Ktotspec;
int Kfinalspec;

int meetings[100][100][100][100] = {0};
int Kinst[250] = {0};

//Observables
//t=0
double COMRing0[Nmax][3];
double position0[Nmax][Nbeadsmax][3];
//t>0
double COMRing[Nmax][3];
double COMRingWrapped[Nmax][3];
double position[Nmax][Nbeadsmax][3];
double positionWrapped[Nmax][Nbeadsmax][3];
double COM[3];
//
double MSDRing[Nmax][3];
double DISRing[Nmax][3];
double MSD[3];
double DIS[3];
double Rfinal=0;
double R2final=0;
//

void computeRg(); 
void computeMSD(); 
void computeK();
void computeKspec();

double gaussLnkNumb(double c1[][3], double c2[][3]);
double min(double , double);
double Compdist(double, double, double, double,double,double);
double reduceInt(int v[], int s);

int main(int argc, char* argv[]){

cout << "Write argv[1]: datain; argv[2] = dataout; argv[3]=Tmax; argv[4]=append"<< endl;
cout << "Num of files to convert:";
Kmax=int(atoi(argv[3]));
cout << "Number of Polymers:";
N=atoi(argv[4]);
cout << N <<endl;
cout<<"Number of Beads in each polymer:";
Nbeads=atoi(argv[5]);
cout << Nbeads<<endl; 
cout<<"Starting time:";
Kstart=atoi(argv[6]);
cout << Kstart<<endl;


for(int n=0;n<Nmax;n++)for(int i=0;i<3;i++) COMRing0[n][i]= COMRingWrapped[n][i]=0;

stringstream writeFile;
writeFile <<"Details2_"<<argv[2];
ofstream write(writeFile.str().c_str());
write << "#t/1e3\tMSD\tR_g^2\tDIS"<< endl;

stringstream writeFileMSD;
writeFileMSD <<"MSD2_"<<argv[2];
ofstream writeMSD(writeFileMSD.str().c_str());
writeMSD << "#t/1e3\tMSD"<< endl;

stringstream writeFileK;
writeFileK <<"Details2_K_"<<argv[2];
ofstream writeK(writeFileK.str().c_str());
writeK << "#t/1e3\tK\tKinst"<< endl;

stringstream writeFileKspec;
writeFileKspec <<"Details2_K_spec"<<argv[2];
ofstream writeKspec(writeFileKspec.str().c_str());
writeKspec << "#t/1e3\tKspec"<< endl;



///////////////////////////////////////////////////
//MAIN LOOP OVER TIME!!!
///////////////////////////////////////////////////

//Kmax-= Kstart;
for(k=0;k<Kmax;k++){

int blob = k;
ifstream indata;
stringstream readFile;
readFile.str("");
readFile.clear();
long long int sum = int((blob+Kstart));
if(sum==0) readFile << argv[1] <<sum;
if(sum>0)readFile << argv[1] <<sum ;
indata.open(readFile.str().c_str());
cout << readFile.str().c_str()<<endl;

for(int n=0;n<N;n++)for(int i=0;i<3;i++)COMRing[n][i]=COM[i]=0;

long int id,type,mol;
double num1,num2,num3;
double x,y,z;
string dummy;
long double time;
long int Ntot;
double l1,l2;

//read 10 lines
for(int i=0;i<10;i++){
if(i==1) {
indata >> time;
time = k*DeltaT*dt;
cout << "time " << time <<endl;
}
if(i==3) {
indata >> Ntot;
cout << "Ntot " << Ntot<<endl;
}
if(i==5) {
indata >> l1 >> l2;
cout << "L " << l1<< " " << l2 <<endl;
Lx = l2-l1;
cout << "Lx " << Lx <<endl;
}
if(i==6) {
indata >> l1 >> l2;
cout << "L " << l1<< " " << l2 <<endl;
Ly = l2-l1;
cout << "Ly " << Ly <<endl;
}
if(i==7) {
indata >> l1 >> l2;
cout << "L " << l1<< " " << l2 <<endl;
Lz = l2-l1;
cout << "Lz " << Lz<<endl;
}

else getline(indata,dummy);
cout << dummy<<endl;
}
 
//////////////////////////
//READ FILES
////////////////////////////
for(int n=0; n<Ntot; n++){
	indata >> id>> type>> x>>y>>z>>num1>>num2>>num3;
	//indata >> id>> mol>>type>>x>>y>>z>>num1>>num2>>num3;
	//cout << id << " " << x <<endl; cin.get();

	int Nring = floor((double)(id-1)/(1.0*Nbeads));
	//int Nring = floor((double)(type-2));
	//int Nring = floor((double)(type-1));
	
	if(Nring<N){
	position[Nring][(id-1)%Nbeads][0] = (x + Lx*num1);
	position[Nring][(id-1)%Nbeads][1] = (y + Ly*num2);
	position[Nring][(id-1)%Nbeads][2] = (z + Lz*num3);

	positionWrapped[Nring][(id-1)%Nbeads][0] = x;
	positionWrapped[Nring][(id-1)%Nbeads][1] = y;
	positionWrapped[Nring][(id-1)%Nbeads][2] = z;

	if(k==0){
	position0[Nring][(id-1)%Nbeads][0] = position[Nring][(id-1)%Nbeads][0];
	position0[Nring][(id-1)%Nbeads][1] = position[Nring][(id-1)%Nbeads][1];
	position0[Nring][(id-1)%Nbeads][2] = position[Nring][(id-1)%Nbeads][2];

	COMRing0[Nring][0] += (position0[Nring][(id-1)%Nbeads][0])/Nbeads;
	COMRing0[Nring][1] += (position0[Nring][(id-1)%Nbeads][1])/Nbeads;
	COMRing0[Nring][2] += (position0[Nring][(id-1)%Nbeads][2])/Nbeads;
	}
	
	COMRing[Nring][0] += (x + Lx*num1)/Nbeads;
	COMRing[Nring][1] += (y + Ly*num2)/Nbeads;
	COMRing[Nring][2] +=(z + Lz*num3)/Nbeads;

	COM[0] += (x + Lx*num1)/(double(N*Nbeads));
	COM[1] += (y + Ly*num2)/(double(N*Nbeads));	
	COM[2] += (z + Lz*num3)/(double(N*Nbeads));
	}	
}

//////////////
//CENTRE wrt initial COM
//////////////
//resetCoord();
//
//computeImages();
//
//RG/////////////////////////
Rgave = Rg2ave = Rg4ave = 0.0;
computeRg();
//////////////////////////////

//MSD and DIS/////////////////////
MSD[0] = MSD[1] = MSD[2] =0.0;
DIS[0] = DIS[1] = DIS[2] =0.0;
computeMSD();
///////////////////////////

///K calculation////////////////
Ktot = 0;
computeK();

Ktotspec = 0;
computeKspec();
////////////////////////////////

write << time << " " <<(MSD[0] + MSD[1] + MSD[2]) << " " << Rg2ave << " " << (DIS[0]+DIS[1]+DIS[2]) << endl; 

writeMSD << time << " ";
for(int np=0;np<N;np++){
writeMSD <<  MSDRing[np][0]+MSDRing[np][1]+MSDRing[np][2] << " ";
}
writeMSD<<endl;

writeK << time << " " << Ktot << " " << Kinst[k] << endl;

writeKspec << time << " " << Ktotspec << endl;

////////////////////////////


Rfinal += Rg2ave;
R2final += Rg2ave*Rg2ave;


Kfinal += Ktot;

Kfinalspec += Ktotspec;

}//closes time (loop over k)

Rfinal=Rfinal*1.0/Kmax;
R2final=R2final*1.0/Kmax;

cout << "<R^2_g> \t\t" << Rfinal << " +/- " << sqrt(R2final-Rfinal*Rfinal) <<endl;

cout << "K_final \t\t" << Kfinal <<endl;

cout << "K_final_spec \t\t" << Kfinalspec <<endl;
//int tcomp =100;

return 0 ;
}


//////////////////////////////////////
//FUNCTIONS
////////////////////////////////////////
double reduceInt(int v[], int s){
double r=0.;
	for(int n=0;n<s;n++) r+=v[n];
return r;	
}

void computeRg(){
for(int n=0; n<N; n++){
Rg[n]=0;
	for(int m =0;m<Nbeads; m++){	
		Rg[n] = Rg[n] + ((position[n][m][0] - COMRing[n][0])*(position[n][m][0] - COMRing[n][0]) + (position[n][m][1] - COMRing[n][1])*(position[n][m][1] - COMRing[n][1]) + (position[n][m][2] - COMRing[n][2])*(position[n][m][2] - COMRing[n][2]));		
	//cout << n << " " << m <<" "<< Rg[n] <<endl;//cin.get();	
	}
Rg[n]=Rg[n]*1.0/Nbeads;
Rg2ave += Rg[n];	
Rg4ave += pow(Rg[n],2.0);	
//cout << "Rgave " << Rg2ave <<endl;//cin.get();	
}
Rg2ave=Rg2ave*1.0/N;
Rg4ave=Rg4ave*1.0/N;
}
	
void computeMSD(){
	
for(int n=0 ; n<N; n++){
MSDRing[n][0]= pow((COMRing[n][0] - COMRing0[n][0]),2.0);
MSDRing[n][1]= pow((COMRing[n][1] - COMRing0[n][1]),2.0);
MSDRing[n][2]= pow((COMRing[n][2] - COMRing0[n][2]),2.0);

for(int m=0; m<Nbeads; m++){
DISRing[n][0] += pow((position[n][m][0] - position0[n][m][0]),2.0)/Nbeads;
DISRing[n][1] += pow((position[n][m][1] - position0[n][m][1]),2.0)/Nbeads;
DISRing[n][2] += pow((position[n][m][2] - position0[n][m][2]),2.0)/Nbeads;
}

MSD[0] += MSDRing[n][0]/(1.0*(N));
MSD[1] += MSDRing[n][1]/(1.0*(N));
MSD[2] += MSDRing[n][2]/(1.0*(N));

DIS[0] += DISRing[n][0]/(1.0*(N));
DIS[1] += DISRing[n][1]/(1.0*(N));
DIS[2] += DISRing[n][2]/(1.0*(N));
}	
}	

void computeK(){
Ktot = 0;
extern long int Kmax;
extern int meetings[100][100][100][100];
extern int Kinst[250];
for(int n=0; n<N; n++){
        for(int m=0;m<Nbeads; m++){
		for(int j=n+1; j<N; j++){
			for(int l=0; l<Nbeads; l++){
                		if((pow((position[j][l][0]-position[n][m][0]),2)+pow((position[j][l][1]-position[n][m][1]),2)+pow((position[j][l][2]-position[n][m][2]),2))<1.5){
				Ktot++;
				if(meetings[n][m][j][l]==0){
				meetings[n][m][j][l] = 1;
				Kinst[k]++;
				}
				}
        //cout << n << " " << m <<" "<< j << " " << l <<" " Ktot <<endl;//cin.get();       
        }
}
}
}
}
        

void computeKspec(){
Ktotspec = 0;

for(int n=0; n<N; n++){
        for(int m=0;m<Nbeads; m+=(Nbeads-1)){
                for(int j=n+1; j<N; j++){
                        for(int l=0; l<Nbeads; l+=(Nbeads-1)){
                                if((pow((position[j][l][0]-position[n][m][0]),2)+pow((position[j][l][1]-position[n][m][1]),2)+pow((position[j][l][2]-position[n][m][2]),2))<7){
                                Ktotspec++;
                                }
        //cout << n << " " << m <<" "<< j << " " << l <<" " Ktot_spec <<endl;//cin.get();       
        }
}
}
}
}
        



double min(double a , double b){
	if(a<b) return a;
	else return b;
}
double Compdist(double cx,double cy, double cz, double c1x, double c1y, double c1z){
	double d=0.0;
	d = sqrt((cx-c1x)*(cx-c1x) + (cy-c1y)*(cy-c1y) + (cz-c1z)*(cz-c1z));
	return d;
}


