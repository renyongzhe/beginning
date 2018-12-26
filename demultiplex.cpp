#include <iostream>
#include <string>
#include <string.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <zlib.h>
#define BUF_SIZE 300
using namespace std;


// split filename
vector<string> split(const string &str,const string &pattern);
//get prefix of filename
string getPrefix(string name);
bool seqCompare(string s1,string s2);

//vector & fq_read(string f);
void demultiplex(char* r1,char* o1,char* index1,char* index2);

int main(int argc,char* argv[])
{
    //demultiplex(string r1,string r2,string p5,string p7);
    //newr1=getPrefix(argv[1]);
    //newr2=getPrefix(argv[2]);
    if (argc==5){
	demultiplex(argv[1],argv[2],argv[3],argv[4]);
    }
    else{
	cout<<"\nUsage: " << argv[0]<<" test.fq.gz out.fq.gz index1 index2.\nPaired-end fastq file should demultiplex one by one.\n"<<endl;
	exit(1);
    }
    return 0;
}

vector<string> split(const string &str,const string &pattern)
{
    char *strc = new char[strlen(str.c_str())+1];
    strcpy(strc,str.c_str());
    vector<string> res;
    char *tmpStr = strtok(strc,pattern.c_str());
    while (tmpStr != NULL)
    {
	res.push_back(string(tmpStr));
	tmpStr = strtok(NULL,pattern.c_str());
    }
    delete[] strc;
    return res;
}

string getPrefix(string name)
{
    vector<string> nameVec;
    string newName;
    nameVec=split(name,".");
    return nameVec[0];
}

void demultiplex(char *r1,char *o1,char *index1,char* index2)
{
    gzFile fq1,out1;
    out1 = gzopen(o1,"wb");
    char seq[BUF_SIZE];
    unsigned int lineNum=0;
    fq1 = gzopen(r1,"rb");
    bool isWrite = false;
    while(gzgets(fq1,seq,BUF_SIZE)){
	lineNum++;
	if(lineNum%4 == 1)
	{
	    if(strstr(seq,index1) && strstr(seq,index2))
	    {	
		//cout << lineNum<<seq<<endl;
		isWrite = true;
		gzwrite(out1,seq,strlen(seq));
	    }else{
		isWrite = false;
	    }
	}
	else
	{
	    if(isWrite == true)
	    {
		//cout << lineNum << seq <<endl;
		gzwrite(out1,seq,strlen(seq));
	    }
	}
    }
    gzclose(fq1);
    gzclose(out1);
}
