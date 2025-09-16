/***************************************************************************************************************************************************
*
*                                 This program will read segment and output condense files
*                                 
****************************************************************************************************************************************************/
#include <time.h>
#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <errno.h>
#include <sys/ipc.h> 
#include <sys/shm.h> 
#include <sys/sysinfo.h>
#include <unistd.h>
#include <sched.h>
#include <sys/syscall.h>
#include <omp.h>
#include <stdint.h>

//#include "readinteger.h"
//#include "readreal.h"
//#include "readnegativereal.h"

#define NBMAXMEOISIS 10
#define NMAX 664 
#define NCHR 23
#define NSNPPERCHR 410000
#define MAXPOP 173200
#define MAXLENGHTNAME 100
#define MAXNBID 173200

#define MAXFAM 80000 
#define NBGROUPPERPERSON NMAX+1
#define MAXSEGMENTPAIR 200 
#define MAXSUMLENPAIR 35000
#define MAXCLOSERELAT 100000
#define MAXCLOSERELATTEMP 500000
#define MAXID 6026457
#define MINPIHAT 0
#define MAXPIHAT 1
#define NBINMATRIX1 100
#define NBINMATRIX2 200
#define NBINMATRIX3 500
#define NBINMATRIX4 1000
#define NBINMATRIX5 2000
#define NBINMATRIX6 5000
#define NBINMATRIX7 10000
#define nbbyteforsegment 7
#define MAXSEGINDIVTEMP 2000000
#define MAXSEGINDIV 120000
#define MAXNBTRIO 978
#define NBCLOSER 4351
#define MAXFILE 6
#define NBINDIVEA 20
#define NBTOURNAMENT NBINDIVEA/4
#define MAXGEN 50
#define MAXVARIABLE 20
#define MAXNBPAIR  152252 
#define MAXNBPO 3719

#define X1SIZE 18
#define Y1SIZE 25
const int XSTART[10] = {648, 71, 136, 191, 263,327,390,448,519,576}; 
#define Y1START 30


//unsigned char snpcall[(int64_t ) MAXPOP][NSNPPERCHR];//
unsigned char * snpcall;

char nameindiv[MAXPOP][MAXLENGHTNAME];
	
#define degree0to1 0.9
#define degree1to2 0.375
#define degree2to3 (0.25+0.125)/2
#define degree3to4 (0.125+0.06125)/2
#define degree4to5 (0.06125+0.06125/2)/2
#define degree5to6 (0.06125/2+0.06125/4)/2 
#define degree6to7 (0.06125/4+0.06125/8)/2
#define degree7to8 (0.06125/8+0.06125/16)/2
#define degree8to9 (0.06125/16+0.06125/32)/2
#define degree9 (0.06125/32+0.06125/64)/2
int IDjob=0;

unsigned int focaltotlenghtpos=0;
unsigned int focaltotseg=0;



int printdetail;
int IDtorank[MAXID];
//int ranktoIDiffam[MAXID]={0};
//int IDtorankinsnpcall[MAXID]={0};

int nbsnpperchr[23];
unsigned long long sizechr[23];
	
//int nblistrelat[MAXPOP];
//int listrelat[MAXPOP][13];



typedef struct
{	int deb;
	int fin;
	int hap;
} typesegments;


typesegments tabseg[23][10];


typedef struct
{	int chi;
	unsigned char haptohap;
} resultphaing;

char namefamily[MAXFAM][MAXLENGHTNAME];

typedef struct
{	int ID;
	int off;
	int p1;
	int p2;
} strucfam;

strucfam famstrt[MAXFAM];

unsigned char transmit[MAXFAM][NSNPPERCHR]; // 0 to 15 first p1 16 to 255 p2

typedef struct
{	int ID;
	unsigned char coef;
} indivclose;
	
//char tabIDtext[MAXNBID][100];
 
typedef struct 
{	int IDoffspring;
	int IDp1;
	float score;
} structPO;

//structPO tabPO[MAXNBPO];

//int isparent[MAXPOP]; //-1 not offspring;>-1 trio number 
//int isoffspring[MAXPOP]; //-1 not offspring;>-1 trio number 



//int seuil[23][23];

unsigned char * genomes;

uint64_t * tempload;


/*
typedef struct 
{	unsigned long long pos;
//	int cm;
} structSNP;

structSNP map[NSNPPERCHR];*/


//snpinteger genomeinteger[MAXPOP][5230];//NSNPPERCHR
//5227
//int listrelativestart[MAXPOP];
//int listrelativeend[MAXPOP];
//int listrelativeP1[MAXPOP];
//int listrelativeP2[MAXPOP];
	
//structSNP map[NCHR][NSNPPERCHR];

void delay(int number_of_seconds)  
{ 
    // Converting time into milli_seconds 
    int milli_seconds = 1000 * number_of_seconds; 
  
    // Stroing start time 
    clock_t start_time = clock(); 
  
    // looping till required time is not acheived 
    while (clock() < start_time + milli_seconds); 
} 

int GetCPUCount()
{ 	cpu_set_t cs;
	CPU_ZERO(&cs);
	sched_getaffinity(0, sizeof(cs), &cs);
	int count = 0;
	for (int i = 0; i < 20; i++)
	{  if (CPU_ISSET(i, &cs))
		count++;
	};
	return count;
}

//
//  main program
//	

int main(int argc, char *argv[])
{	//return 0;
/*	int count=0;
	for(int IDrun=0;IDrun<255;IDrun++)// 23
	{	for(int IDrun1=0;IDrun1<255;IDrun1++)// 23
		{	for(int IDrun2=0;IDrun2<255;IDrun2++)// 23
			{	for(int IDrun3=0;IDrun3<255;IDrun3++)// 23
				{	if (IDrun==IDrun1) {} else if (IDrun==IDrun2) {} else	if (IDrun2==IDrun1) count++;
					//if IDrun3==IDrun1;
				};
			};
		};	
	};
	 
	clock_t end = clock();
	float elapsed_secs = (float)(end - begin0);
    printf("\nTotal time 4:%f cpu click so %f seconds",elapsed_secs,elapsed_secs/CLOCKS_PER_SEC );
	return 0;*/
	//printf("This program takes six possible arguments: \n -min_snp N where N+1 is the minimal number of SNPs in a segment, \n -hap_file /path/hapfile.hap where /path/hapfile.hap");
	//printf("is a hap file,\n -IDs /path/IDs.txt where /path/IDs.txt is a file with a list of IDs that was used to generated the hap file,\n -pair_file /path/pair.txt where");
	//printf("/path/pair.txt is a file of pairs: ID1 ID2 (corresponding to IDs.txt) for which segments are outputed,\n -result /path/result.txt where /path/result.txt is the file ");  
	//printf("to store the result\n");
	char result[500]; 
	char famfile[500]; 
	int input;
	char hapfile[500];
	char listID[500];
	int off;
	int parent1;
	int parent2;
	char flagPO[500]; 
	//	readID(pathID);
	//readtrio();
	int chr;
	printf("Reading parameters\n");
	
	for(input=1;input<argc;input++)
	{/*	if( strncmp(argv[input], "-off", strlen("-off")) == 0 && input < argc-1) 
		{	off = atoi(argv[++input]);
		
		} 
		else if( strncmp(argv[input], "-parent1", strlen("-parent1")) == 0 && input < argc-1) 
		{	parent1 = atoi(argv[++input]);
		} 
		else if( strncmp(argv[input], "-parent2", strlen("-parent2")) == 0 && input < argc-1) 
		{	parent2 = atoi(argv[++input]);
		} */
		if( strncmp(argv[input], "-fam_file", strlen("-fam_file")) == 0 && input < argc-1) strcpy(famfile,argv[++input]);
		else if( strncmp(argv[input], "-hap_file", strlen("-hap_file")) == 0 && input < argc-1) strcpy(hapfile,argv[++input]);
		else if( strncmp(argv[input], "-list_ID", strlen("-list_ID")) == 0 && input < argc-1) strcpy(listID,argv[++input]);
		else if( strncmp(argv[input], "-output", strlen("-result")) == 0 && input < argc-1) strcpy(result,argv[++input]);
		else if( strncmp(argv[input], "-CHR", strlen("-CHR")) == 0 && input < argc-1) 	chr = atoi(argv[++input]); 
		else if( strncmp(argv[input], "-OPright", strlen("-OPright")) == 0 && input < argc-1) strcpy(flagPO,argv[++input]);	 
	};
	int  chrtemp1=chr;
	//printf("offspring is %d parent1 is %d and parent2 is %d\n",off,parent1,parent2);
	/*int64_t * scorechr;
	scorechr = (int64_t*) calloc ((unsigned long long) 23*2*2, sizeof(int64_t));
	for(int chrtemp1=1;chrtemp1<23;chrtemp1++)// if (ID<10)
	{	 int64_t *p = (scorechr+chrtemp1*4); 
		(*(p+0*2+0  ))=7+chrtemp1; 
		(*(scorechr+chrtemp1*4+0*2+1  ))=13+chrtemp1; 
		(*(scorechr+chrtemp1*4+1*2+0  ))=51+chrtemp1; 
		(*(scorechr+chrtemp1*4+1*2+1  ))=69+chrtemp1; 
	};
	for(int chrtemp1=1;chrtemp1<23;chrtemp1++)// if (ID<10)
	{	printf("%d %d %d %d\n",(*(scorechr+chrtemp1*4+0*2+0  )),(*(scorechr+chrtemp1*4+0*2+1  )),(*(scorechr+chrtemp1*4+1*2+0  )),(*(scorechr+chrtemp1*4+1*2+1  )));
	};
	return 0;*/
	
	/*seuil[22][21]=790;seuil[21][22]=790;
	seuil[20][21]=794;	seuil[21][20]=794;
	seuil[19][22]=786;	seuil[22][19]=786;
	seuil[20][22]=788;	seuil[22][20]=788;
	seuil[20][19]=788;	seuil[19][20]=788;
	seuil[21][19]=790;	seuil[19][21]=790;
	seuil[22][18]=790;	seuil[18][22]=790;
	seuil[21][18]=790;	seuil[18][21]=790;
	seuil[19][18]=790;	seuil[18][19]=790;
	seuil[20][18]=790;	seuil[18][20]=790;
	for(int chrtemp1=chrtodo1;chrtemp1<chr2+1;chrtemp1++)// if (ID<10)
	{	for(int chrtemp2=chrtodo1;chrtemp2<chr2+1;chrtemp2++) 
		{	seuil[chrtemp1][chrtemp2]=790;
		};
	};*/
	
	FILE * MAFfile;
	char first;

	printf("Reading ID file\n");
	
	
/*	if ((MAFfile = fopen(listID, "r")) == NULL) 
	{	if (printdetail) printf("file triofile is not found\n");		
		return (1);
	};
	do {first=getc(MAFfile);} while (first!='\n' && first!=EOF); 	
	do {first=getc(MAFfile);} while (first!='\n' && first!=EOF); 	
	int rank=0;
	do
	{	first=getc(MAFfile);
		if (first!=EOF)
		{	int ID=0;
			int place=0;
			do 
			{	//ID=ID*10+first-48;
				first=getc(MAFfile);
				tabIDtext[rank][place]=first;
				place++;
			} while (first!=' ' && first!=9);
			tabIDtext[rank][place]='\0';
			
		//	IDtorank[ID]=rank;
			
			printf("ID %d rank %d\n",ID, rank);
			rank++;
		};
		do {first=getc(MAFfile);} while (first!='\n' && first!=EOF); 	
	}
	while (first!=EOF); 
	fclose(MAFfile);*/
	
	/*if ((MAFfile = fopen(listID, "r")) == NULL) 
	{	if (printdetail) printf("file triofile is not found\n");		
		return (1);
	};
	do {first=getc(MAFfile);} while (first!='\n' && first!=EOF); 	
	do {first=getc(MAFfile);} while (first!='\n' && first!=EOF); 	
	int rank=0;
	first=getc(MAFfile);
	do
	{	int numname=0;
		do 
		{	nameindiv[rank][numname]=first;
			numname++;
			first=getc(MAFfile);
			printf("%d %c\n",numname,first);
		} while  (first!=EOF && first!=' ' && first!=9  && first!='\n' && first!='\0' ); 
		nameindiv[rank][numname]='\0';
		
		printf("%d %s\n",rank,nameindiv[rank]);	
		rank++;
		do {first=getc(MAFfile);} while (first!='\n' && first!=EOF); 	
		first=getc(MAFfile);
	}
	while (first!=EOF);
	fclose(MAFfile);
	//nameindiv
	
	printf("Reading pedigree file\n");
	*/
	
	
	
	/*
	if ((MAFfile = fopen(famfile, "r")) == NULL) 
	{	if (printdetail) printf("file triofile is not found\n");		
		return (1);
	};
	int nbfam=0;
	int nbline=0;
	do {first=getc(MAFfile);} while (first!='\n' && first!=EOF); 	
	
	
	do
	{	first=getc(MAFfile);
		printf("%c\n",first);
		nbline++;
		int IDfam=0;
		if (first!=EOF)
		{	char IDtext[100];
			IDtext[0]='\0';
			int place=0;
			int ID=0;
			do 
			{	printf("%d %c\n",place,first);
				namefamily[nbfam][place]=first; //ID=ID*10+first-48;
				place++;
				first=getc(MAFfile);
			} while (first!=' ' && first!='\t');
			namefamily[nbfam][place]='\0';
			
			famstrt[nbfam].ID=nbfam;
			printf("Family n %d ID %s\n",nbfam,namefamily[nbfam]);
			IDfam++;
			
			
			
			first=getc(MAFfile);
			
			IDtext[0]='\0';
			place=0;
			do 
			{	//ID=ID*10+first-48;
				IDtext[place]=first;
				place++;
				first=getc(MAFfile);
			} while (first!=' ' && first!='\t');
			IDtext[place]='\0';
			int rank=0;
			while (strcmp(IDtext,nameindiv[rank])!=0 && rank<MAXNBID) 
			{	//printf("%s %s\n",IDtext,nameindiv[rank]);
				rank++;
			}
			if (rank==MAXNBID) 
			{	printf("ID %s not found in ID file\n",IDtext);
				exit(0);
			}
			printf("ID of offspring %s found as %s is in rank %d\n",IDtext,nameindiv[rank],rank);
			famstrt[nbfam].off=rank;
		
			//do {first=getc(MAFfile);} while (first!='\t' && first!=EOF); 	
			first=getc(MAFfile);
			int atleastaparent=0;
			
			if (first=='N') 
			{	atleastaparent=1;
				famstrt[nbfam].p1=-1;
				getc(MAFfile);getc(MAFfile);
			} else 
			{	int rank=0;
				IDtext[0]='\0';
				place=0;
				do 
				{	//ID=ID*10+first-48;
					printf("%d %c\n",place,first);
					IDtext[place]=first;
					place++;
					first=getc(MAFfile);
				
				} while (first!=' ' && first!='\t' && first!='\n' && first!=EOF);
				IDtext[place]='\0';
				rank=0;
				while (strcmp(IDtext,nameindiv[rank])!=0 && rank<MAXNBID) 
				{	//printf("%s %s\n",IDtext,nameindiv[rank]);
					rank++;
				}
				if (rank==MAXNBID) 
				{	printf("ID %s of dad not found in ID file\n",IDtext);
					exit(0);
				}

				printf("%s %s\n",IDtext,nameindiv[rank]);
				famstrt[nbfam].p1=rank;
			};
		//	do {first=getc(MAFfile);} while (first!='\t' && first!=EOF); 	
			first=getc(MAFfile);
			if (first!='N' || atleastaparent==0)
			{	if (first!='N')
				{	int rank=0;
					IDtext[0]='\0';
					place=0;
					do 
					{	//ID=ID*10+first-48;
						printf("%d %c\n",place,first);
						IDtext[place]=first;
						place++;
						first=getc(MAFfile);
					} while (first!=' ' && first!='\t' && first!='\n' && first!=EOF);
					if (rank==MAXNBID) 
					{	printf("ID %s of mum not found in ID file\n",IDtext);
						exit(0);
					}

					rank=0;
					IDtext[place]='\0';
					while (strcmp(IDtext,nameindiv[rank])!=0 && rank<MAXNBID)
					{//printf("%s %s\n",IDtext,nameindiv[rank]);
						rank++;
					}
					printf("%s %s\n",IDtext,nameindiv[rank]);
					famstrt[nbfam].p2=rank;
					if (rank==MAXNBID) 
					{	printf("ID %s not found in ID file\n",IDtext);
						exit(0);
					}
				} else 
				{	famstrt[nbfam].p2=-1;
					getc(MAFfile);
					printf("%c\n",first);
					getc(MAFfile);
					printf("%c\n",first);
				};
				printf("fam %d ID %s IDoff %d IDp1 %d IDp2 %d\n",nbfam,namefamily[nbfam],famstrt[nbfam].off,famstrt[nbfam].p1,famstrt[nbfam].p2); 
			//	ranktoIDiffam[IDtorank[famstrt[nbfam].off]]=famstrt[nbfam].off;
			//	if (famstrt[nbfam].p1>-1) ranktoIDiffam[IDtorank[famstrt[nbfam].p1]]=famstrt[nbfam].p1;
			//	if (famstrt[nbfam].p2>-1) ranktoIDiffam[IDtorank[famstrt[nbfam].p2]]=famstrt[nbfam].p2;
			//	famstrt[nbfam].ID=nbline;
				nbfam++;
				if (MAXFAM<=nbfam) 
				{	printf("Too many families %d %d \n",MAXFAM,nbfam);
					exit(0);
				}
				//do {first=getc(MAFfile);} while (first!='\n' && first!=EOF); 	
			}	else			
			{	famstrt[nbfam].p2=-1;
				getc(MAFfile);
				printf("%c\n",first);
				getc(MAFfile);
				printf("%c\n",first);
			}
		}
		//if (first!=EOF) first=getc(MAFfile);
	} while (first!=EOF);
	
	fclose(MAFfile);
	for (int fam=0;fam<nbfam;fam++)//nbfam
	{	printf("%d %d %d \n",famstrt[fam].off,famstrt[fam].p1,famstrt[fam].p2);
	}	
	//exit(0);
	//nbfam--;
	
	
	*/
	
	
	
	/*
	if ((MAFfile = fopen(hapfile, "r")) == NULL) 
	{	return (1);
	};
	 
	int nbsnpinfile=0;
	int nbindivinfile=0;
	char namesnp[NSNPPERCHR][100];
	printf("Read Hap file\n");
	char allele1[NSNPPERCHR];
	char allele2[NSNPPERCHR];
	for (int fam=0;fam<nbfam;fam++)//nbfam
	{	printf("%d %d %d \n",famstrt[fam].off,famstrt[fam].p1,famstrt[fam].p2);
	}	
	snpcall = (unsigned char*) calloc ((unsigned long long) MAXPOP*NSNPPERCHR, sizeof(char));
	for(int snp=0;snp<NSNPPERCHR;snp++)		//22
	{	printf("Reading SNP %d\n",snp);
		nbsnpinfile++;
		int nb0=0;
		int nb1=0;
		int nb2=0;
		int nb3=0;
		
		do {first=getc(MAFfile);} while (first!=32 && first!=EOF); 	
		int pos=0;
		do 
		{	first=getc(MAFfile);
			namesnp[snp][pos]=first;
			pos++;
		} while (first!=32 && first!='\t' && first!=EOF); 	
		namesnp[snp][pos]='\0';
		printf("%s\n",namesnp[snp]);
		
		do {first=getc(MAFfile);} while (first!=32 && first!='\t' && first!=EOF); 	//229133-223882
		allele1[snp]=getc(MAFfile);
		getc(MAFfile);
		allele2[snp]=getc(MAFfile);
		getc(MAFfile);
		printf("%c %c \n",allele1[snp],allele2[snp]);
	//	do {first=getc(MAFfile);} while (first!=32 && first!=EOF); 	
	//	do {first=getc(MAFfile);} while (first!=32 && first!=EOF); 	
		printf("%d %d %d \n",famstrt[133].off,famstrt[133].p1,famstrt[133].p2);
		
		nbindivinfile=0;
		if (first!=EOF) 
		{	int nbnb2=0;
			int nbnb3=0;
			int nbnb0=0;
			int nbnb1=0;
	
			for(int place=0;place<MAXID;place++)		//22
			{	if (snp==0) printf("place %d \n",place);
				char c1=getc(MAFfile);	
			//	printf("%c\n",c1);
			//	if (snp<2)
				//printf("a%c ",c1);
				if (c1!=EOF && c1!='\n')
				{	getc(MAFfile);
					char c2=getc(MAFfile);
					//if (snp<2)
					//printf("b%c ",c2);
					char c3=getc(MAFfile);
					
			//		if (ranktoIDiffam[place]) 
					{	snpcall[(int64_t ) nbindivinfile+(int64_t )MAXPOP*snp]=4;	
						if (c1=='1' && c2=='1') 
						{	snpcall[(int64_t ) nbindivinfile+(int64_t )MAXPOP*snp]=3;
							nb2++;
							nbnb2++;
						};
						if (c1=='0' && c2=='1') 
						{	 snpcall[(int64_t ) nbindivinfile+(int64_t )MAXPOP*snp]=1;
							nb1++;
							nbnb1++;
						};	
						if (c1=='1' && c2=='0') 
						{	snpcall[(int64_t ) nbindivinfile+(int64_t )MAXPOP*snp]=2;
							nb3++;
							nbnb3++;
						};	
						if (c1=='0' && c2=='0') 
						{	 snpcall[(int64_t ) nbindivinfile+(int64_t )MAXPOP*snp]=0;
							nb0++;
							nbnb0++;
						};
					//	IDtorankinsnpcall[(int64_t ) ranktoIDiffam[place]]=nbindivinfile;
					//	if (snp==0) printf("rank in snp call %d ID %d \n",IDtorankinsnpcall[(int64_t ) ranktoIDiffam[place]],ranktoIDiffam[place] );
						nbindivinfile++;
					};
					if (c3!=' ') 
					{	printf("%d\n",place);
						place=MAXID;
					};		
				} else
				{	printf("%d\n",place);
					place=MAXID;
				};		
			};
			printf("%d %d %d %d %d\n",nbnb0,nbnb1,nbnb2,nbnb3,nbindivinfile);
		} else snp=NSNPPERCHR;
		// while (first!='\n'  && first!=EOF) {first=getc(MAFfile);};
		//printf("%d %d %d %d %d \n",snp,nb0,nb1,nb2,nb3);
		
	};
	fclose(MAFfile);
	printf("there is %d SNPs in file.Reading complete !",nbsnpinfile);
	nbsnpinfile--;
	*/
	
	
	
	
	
	//readrelat();
/*	if (printdetail) printf("read trio done \n");
	if (printdetail) printf("read ID done \n");
	//genomes = (unsigned char*) calloc ((unsigned long long) (1+nbsnpperchr[chrtodo1]/4)*MAXPOP, sizeof(char));
	printf("%d\n",sizeof(genomes));
	clock_t step1 = clock(); 
	float elapsed_secs1 = (float)(step1-begin0);
    printf("\nTotal time 1:%f cpu click so %f seconds\n",elapsed_secs1,elapsed_secs1/CLOCKS_PER_SEC );
	int nbbadphase[MAXNBPO];
	for(int offspring=0;offspring<MAXNBPO;offspring++)//
	{	nbbadphase[offspring]=0;
	};
	int nbphaseerror[NSNPPERCHR][23];
	int nbSNPerror[NSNPPERCHR][23];
	for(int chrtemp1=1;chrtemp1<22+1;chrtemp1++)		//22
	{	for(int snp=0;snp<nbsnpperchr[chrtemp1];snp++)		
		{	nbphaseerror[snp][chrtemp1]=0;
			nbSNPerror[snp][chrtemp1]=0;
		};
	};
	int nbsitePO=0;
	int nbsitetrio=0;
	int categoryphaseeroor[2][3]={0};
	int categoryphaseeroorindiv[MAXNBPO][3];
	int nbsnpineachphase[MAXNBPO][23][2];
	 for(int offspring=0;offspring<MAXNBTRIO;offspring++)//MAXNBTRIO
	{	categoryphaseeroorindiv[offspring][0]=0;
		categoryphaseeroorindiv[offspring][1]=0;
		categoryphaseeroorindiv[offspring][2]=0;
		
		for(int chrtemp1=1;chrtemp1<22+1;chrtemp1++)		//22
		{	nbsnpineachphase[offspring][chrtemp1][0]=0;
			nbsnpineachphase[offspring][chrtemp1][1]=0;
		};
	}*/

		//printf("%d\n",version);
	int IDoffspring=101104;
	int IDp1=300319;
	int IDp2=300322;
	int IDavenc=161722;

	
	int countunphase=0;
	int phaseparent1=1;
	int lastphase=1;
	int dephasefound=0;
	int MEerrorfind[NSNPPERCHR][MAXFAM];
	for (int fam=0;fam<MAXFAM;fam++)//MAXFAM
	{	for(int snp=0;snp<nbsnpinfile;snp++)		
		{	MEerrorfind[snp][fam]=0;
		}
	}
	printf("Start analysing");
	
	//for (int fam=0;fam<nbfam;fam++)//nbfam
	{	int off=IDoffspring;
		int parent1=IDp1;
		int parent2=IDp2;
		//printf("%d %d %d\n",ID,IDp1loop,IDp2loop);
		


		for(int seg=0;seg<10;seg++)		
		{	for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
			{	tabseg[	chrtemp1][seg].deb=0		
				tabseg[	chrtemp1][seg].fin=0		
				tabseg[	chrtemp1][seg].hap=0		
			}
	
		}
		
		phaseparent1=1;
		int transmittedfromparent1=-1;
		int transmittedfromparent2=-1;
		int lasthetsnp1=-1;
		int lasthetsnp2=-1;
		int lasthetsnp1mayo=-1;
		int lasthetsnp2mayo=-1;
		printf("%d %d %d\n",off,parent1,parent2);
		for(int snp=0;snp<nbsnpinfile;snp++)		
		{	unsigned char marker=snpcall[(int64_t ) off+(int64_t )MAXPOP*snp];//(*((genomes+(unsigned long long) ID*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)>>(((snp%4)*2)))&3;
			unsigned char markerp1=(famstrt[fam].p1>-1)?snpcall[(int64_t ) parent1+(int64_t )MAXPOP*snp]:snpcall[(int64_t ) off+(int64_t )MAXPOP*snp];//(*((genomes+(unsigned long long) IDp1loop*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)>>(((snp%4)*2)))&3;
			unsigned char markerp2=(famstrt[fam].p2>-1)?snpcall[(int64_t ) parent2+(int64_t )MAXPOP*snp]:snpcall[(int64_t ) off+(int64_t )MAXPOP*snp];//(*((genomes+(unsigned long long) IDp2loop*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)>>(((snp%4)*2)))&3;
		//	if (off==5) printf("%d %d %d %d\n",snp,marker,markerp1,markerp2);
			if (transmittedfromparent1==-1)
			{	if ((phaseparent1==1 && (marker>>1)==(markerp1>>1) ) || (phaseparent1==2 && (marker&1)==(markerp1>>1) ) )
				{	transmittedfromparent1=1;
				} else if  ((phaseparent1==1 && (marker>>1)==(markerp1&1)) || (phaseparent1==2 && (marker&1)==(markerp1&1)) )
				{	transmittedfromparent1=2;
				} else 
				{	transmittedfromparent1=5;
				};

			} else 
			{	if ( markerp1==1 || markerp1==2 )
				{	if ((phaseparent1==1 && (marker>>1)==(markerp1>>1)) || (phaseparent1==2 && (marker&1)==(markerp1>>1)) )
					{	if (transmittedfromparent1!=1)
						{	//	if (lasthetsnp1mayo!=lasthetsnp1) printf("%d %d %d %d %d\n",tabtrio[offspring].IDoffspring,tabtrio[offspring].IDp1,chrtemp1,snp,lasthetsnp1);
							transmittedfromparent1=1;
							lasthetsnp1mayo=snp;
						};
					} else if  ((phaseparent1==1 && (marker>>1)==(markerp1&1)) || (phaseparent1==2 && (marker&1)==(markerp1&1)) )
					{	if (transmittedfromparent1!=2)
						{//	if (lasthetsnp1mayo!=lasthetsnp1) printf("%d %d %d %d %d\n",tabtrio[offspring].IDoffspring,tabtrio[offspring].IDp1,chrtemp1,snp,lasthetsnp1);
							transmittedfromparent1=2;
							lasthetsnp1mayo=snp;
						};
					};
				} else 			
				{	if ((phaseparent1==1 && (marker>>1)!=(markerp1>>1)) || (phaseparent1==2 && (marker&1)!=(markerp1>>1)) )
					{	printf("ME snp %d %d\n",snp,transmittedfromparent1);
						if (transmittedfromparent1==1) transmittedfromparent1=3;
						if (transmittedfromparent1==2) transmittedfromparent1=4;
					}		
				};
			};
			lasthetsnp1=snp;
			int Merror=0;
			
			if (transmittedfromparent1==-1)
			{	transmit[fam][snp]=0b00000011;
			} else if (transmittedfromparent1==1)
			{	transmit[fam][snp]=0b00000001;
			} else  if (transmittedfromparent1==2)
			{	transmit[fam][snp]=0b00000000;
			} else if (transmittedfromparent1==3)
			{	transmit[fam][snp]=0b00000101;
				transmittedfromparent1=transmittedfromparent1-2;
			} else if (transmittedfromparent1==4)
			{	transmit[fam][snp]=0b00000100;
				transmittedfromparent1=transmittedfromparent1-2;
			} else if (transmittedfromparent1==5)
			{	transmit[fam][snp]=0b00001000;
				transmittedfromparent1=-1;
			};
			
			if (transmittedfromparent2==-1)
			{	if ((phaseparent1==1 && (marker&1)==(markerp2>>1)) || (phaseparent1==2 && (marker>>1)==(markerp2>>1)) )
				{	transmittedfromparent2=1;
				} else if ((phaseparent1==1 && (marker&1)==(markerp2&1)) || (phaseparent1==2 && (marker>>1)==(markerp2&1)) )
				{	transmittedfromparent2=2;
				} else 
				{	transmittedfromparent2=5;
				};
			} else 
			{	if (markerp2==1 || markerp2==2 )
				{	if ((phaseparent1==1 && (marker&1)==(markerp2>>1)) || (phaseparent1==2 && (marker>>1)==(markerp2>>1)) )
					{	if (transmittedfromparent2!=1)
						{	//	if (lasthetsnp2mayo!=lasthetsnp2) printf("%d %d %d %d %d\n",tabtrio[offspring].IDoffspring,tabtrio[offspring].IDp2,chrtemp1,snp,lasthetsnp2);
							transmittedfromparent2=1;
							lasthetsnp2mayo=snp;
						};
					} else if ((phaseparent1==1 && (marker&1)==(markerp2&1)) || (phaseparent1==2 && (marker>>1)==(markerp2&1)) )
					{	if (transmittedfromparent2!=2)
						{	//	if (lasthetsnp2mayo!=lasthetsnp2) printf("%d %d %d %d %d\n",tabtrio[offspring].IDoffspring,tabtrio[offspring].IDp2,chrtemp1,snp,lasthetsnp2);
							transmittedfromparent2=2;
							lasthetsnp2mayo=snp;
						};
					};
				}	else 
				{	if ((phaseparent1==2 && (marker>>1)!=(markerp2>>1)) || (phaseparent1==1 && (marker&1)!=(markerp2>>1)) )
					{	printf("ME snp %d %d\n",snp,transmittedfromparent2);
						if (transmittedfromparent2==1) transmittedfromparent2=3;
						if (transmittedfromparent2==2) transmittedfromparent2=4;
					}	
				};
			};
			lasthetsnp2=snp;
			if (off==5) printf("%d\n",transmit[fam][snp]);
			if (transmittedfromparent2==-1)
			{	transmit[fam][snp]=0b00100000 | transmit[fam][snp];
			} else if (transmittedfromparent2==1)
			{	transmit[fam][snp]=0b00010000 | transmit[fam][snp];
			} else if (transmittedfromparent2==2)
			{	transmit[fam][snp]=0b00000000 | transmit[fam][snp];
			} else if (transmittedfromparent2==3)
			{	transmit[fam][snp]=0b01010000  | transmit[fam][snp];
				transmittedfromparent2=transmittedfromparent2-2;  
			} else if (transmittedfromparent2==4)
			{	transmit[fam][snp]=0b01000000  | transmit[fam][snp];
				transmittedfromparent2=transmittedfromparent2-2;
			} else if (transmittedfromparent2==5)
			{	transmit[fam][snp]=0b10000000  | transmit[fam][snp];
				transmittedfromparent2=-1;
			};
		//	if ( (marker==0 && markerp1==3  ) || (marker==3 && markerp1==0 ) || ( marker==2 && markerp1==0) || (marker==3 && markerp2==3) )))
		//	{	transmit[fam][snp]=transmit[fam][snp] & 0b11110100;
		//	
		//	};
		//	if ( (marker==0 && (markerp1==3 || markerp2==3) ) || (marker==3 && (markerp1==0 || markerp2==0)) || ((marker==1 || marker==2) && ((markerp1==0 && markerp2==0) || (markerp1==3 && markerp2==3) )))
		//	{	transmit[fam][snp]=transmit[fam][snp] & 0b01000100;
		//	
		//	};
			if (transmit[fam][snp]==4)  printf("%d %d %d\n",snp ,fam,transmit[fam][snp]);
			if (off==16) printf("%d %d %d %d %d %d\n",fam,snp,transmit[fam][snp],marker,markerp1,markerp2);
		//	printf("%d %d %d %d %d %d %d\n",chrtemp1,snp,marker,markerp1,markerp2,transmittedfromparent1,transmittedfromparent2);
		};
	};
		printf("Analysing done");

	//exit(1);
	
	
	FILE * output;
	FILE * output2;
	char resultfam[100];	
	strcpy(resultfam,result);
	
	strcat(resultfam,"Wmissing");
	if ((output2 = fopen(resultfam, "w")) == NULL) 
	{	printf("file %s can not be opened\n",resultfam);		
		return (0);
	};
	fprintf(output2,"CHR SNP A1 A2 ");	
	 for (int fam=0;fam<nbfam;fam++)//MAXFAM
	{	int off=famstrt[fam].off;
		int parent1=(famstrt[fam].p1>-1)?famstrt[fam].p1:-1;
		int parent2=(famstrt[fam].p2>-1)?famstrt[fam].p2:-1;
		if (famstrt[fam].p1>-1 && famstrt[fam].p2>-1)
		{	fprintf(output2,"%s-%s.%s-%s.T %s-%s.%s-%s.NT %s-%s.%s-%s.T %s-%s.%s-%s.NT ",
						namefamily[fam],nameindiv[famstrt[fam].off],namefamily[fam],nameindiv[famstrt[fam].p1],namefamily[fam],nameindiv[famstrt[fam].off],namefamily[fam],nameindiv[famstrt[fam].p1],
						namefamily[fam],nameindiv[famstrt[fam].off],namefamily[fam],nameindiv[famstrt[fam].p2],namefamily[fam],nameindiv[famstrt[fam].off],namefamily[fam],nameindiv[famstrt[fam].p2]);
		};
		if (famstrt[fam].p1==-1)
		{	fprintf(output2,"%s-%s.%s-NA.T %s-%s.%s-NA.NT %s-%s.%s-%s.T %s-%s.%s-%s.NT ",
						namefamily[fam],nameindiv[famstrt[fam].off],namefamily[fam],namefamily[fam],nameindiv[famstrt[fam].off],namefamily[fam],
						namefamily[fam],nameindiv[famstrt[fam].off],namefamily[fam],nameindiv[famstrt[fam].p2],namefamily[fam],nameindiv[famstrt[fam].off],namefamily[fam],nameindiv[famstrt[fam].p2]);
		} else if (famstrt[fam].p2==-1)
		{	fprintf(output2,"%s-%s.%s-%s.T %s-%s.%s-%s.NT %s-%s.%s-NA.T %s-%s.%s-NA.NT ",
						namefamily[fam],nameindiv[famstrt[fam].off],namefamily[fam],nameindiv[famstrt[fam].p1],namefamily[fam],nameindiv[famstrt[fam].off],namefamily[fam],nameindiv[famstrt[fam].p1],
						namefamily[fam],nameindiv[famstrt[fam].off],namefamily[fam],namefamily[fam],nameindiv[famstrt[fam].off],namefamily[fam]);
		}
	}
	fprintf(output2,"\n");
		
	for(int snp=0;snp<nbsnpinfile;snp++)		
	{	fprintf(output2,"%d %s %c %c ",chr,namesnp[snp],allele1[snp],allele2[snp]);

		for (int fam=0;fam<nbfam;fam++)//MAXFAM
		{			//printf("place %d ",place);
			int off=famstrt[fam].off;
			int parent1=(famstrt[fam].p1>-1)?famstrt[fam].p1:-1;
			int parent2=(famstrt[fam].p2>-1)?famstrt[fam].p2:-1;
		
		//	if (off==5) printf("%d %d %d %d %d %d \n",snp,transmit[fam][snp],transmit[fam][snp]&0b00010000,transmit[fam][snp]&0b00100000,transmit[fam][snp]&0b01000000,transmit[fam][snp]&0b10000000); 
				 
			 if (off==16) printf("SNP %d %d %d %d %d\n",fam,snp,transmit[fam][snp],snpcall[(int64_t ) parent1+(int64_t )MAXPOP*snp],snpcall[(int64_t ) off+(int64_t )MAXPOP*snp]);
		//	if (transmit[fam][snp]==4) printf("%d %d %d\n",snp,fam ,transmit[fam][snp],famstrt[fam].p1,famstrt[fam].p2);
			if (famstrt[fam].p1>-1 && famstrt[fam].p2>-1)
			{	if (transmit[fam][snp]&0b00001000)
				{	fprintf(output2,"NA NA "); 			
				} else if (transmit[fam][snp]&0b00000010)
				{	if (phaseparent1==1) fprintf(output2,"%d %d ",snpcall[(int64_t ) off+(int64_t )MAXPOP*snp]>>1,snpcall[(int64_t ) parent1+(int64_t )MAXPOP*snp]&1); 
					else fprintf(output2,"%d %d ",snpcall[(int64_t ) off+(int64_t )MAXPOP*snp]&1,snpcall[(int64_t ) parent1+(int64_t )MAXPOP*snp]&1); 		
			//
					
				} else if ( transmit[fam][snp]&0b00000100 )
				{	fprintf(output2,"NA NA ");
				}	else 			
				{	if (transmit[fam][snp]&0b00000001) fprintf(output2,"%d %d ",snpcall[(int64_t ) parent1+(int64_t )MAXPOP*snp]>>1,snpcall[(int64_t ) parent1+(int64_t )MAXPOP*snp]&1); 
					else fprintf(output2,"%d %d ",snpcall[(int64_t ) parent1+(int64_t )MAXPOP*snp]&1,snpcall[(int64_t ) parent1+(int64_t )MAXPOP*snp]>>1); 
				};
				if (transmit[fam][snp]&0b10000000)
				{	fprintf(output2,"NA NA "); 			
				} else if (transmit[fam][snp]&0b00100000)
				{	if (phaseparent1==2) fprintf(output2,"%d %d ",snpcall[(int64_t ) off+(int64_t )MAXPOP*snp]&1,snpcall[(int64_t ) parent2+(int64_t )MAXPOP*snp]&1); 
					else fprintf(output2,"%d %d ",snpcall[(int64_t ) off+(int64_t )MAXPOP*snp]>>1,snpcall[(int64_t ) parent2+(int64_t )MAXPOP*snp]&1); 			
				} else if (transmit[fam][snp]&0b01000000)
				{	fprintf(output2,"NA NA ");
				} else
				{	if (transmit[fam][snp]&0b00010000) fprintf(output2,"%d %d ",snpcall[(int64_t ) parent2+(int64_t )MAXPOP*snp]>>1,snpcall[(int64_t ) parent2+(int64_t )MAXPOP*snp]&1); 
					else fprintf(output2,"%d %d ",snpcall[(int64_t ) parent2+(int64_t )MAXPOP*snp]&1,snpcall[(int64_t ) parent2+(int64_t )MAXPOP*snp]>>1); 			
				}
			} else if (famstrt[fam].p1==-1)
			{	if (phaseparent1==1) fprintf(output2,"%d NA ",snpcall[(int64_t ) off+(int64_t )MAXPOP*snp]>>1); 
				else fprintf(output2,"%d NA ",snpcall[(int64_t ) off+(int64_t )MAXPOP*snp]&1); 			
				if (transmit[fam][snp]&0b10000000)
				{	fprintf(output2,"NA NA "); 			
				} else if (transmit[fam][snp]&0b00100000)
				{	if (phaseparent1==2) fprintf(output2,"%d %d ",snpcall[(int64_t ) off+(int64_t )MAXPOP*snp]&1,snpcall[(int64_t ) parent2+(int64_t )MAXPOP*snp]&1); 
					else fprintf(output2,"%d %d ",snpcall[(int64_t ) off+(int64_t )MAXPOP*snp]>>1,snpcall[(int64_t ) parent2+(int64_t )MAXPOP*snp]&1); 	
				//		if (off==5) printf(	"%d %d ",snpcall[(int64_t ) off+(int64_t )MAXPOP*snp]>>1,snpcall[(int64_t ) parent2+(int64_t )MAXPOP*snp]&1); 	
				} else if (transmit[fam][snp]&0b01000000)
				{	fprintf(output2,"NA NA ");
				//	if (off==5) printf(	"%d %d ",-9,-9); 	
				} else
				{	if (transmit[fam][snp]&0b00010000) fprintf(output2,"%d %d ",snpcall[(int64_t ) parent2+(int64_t )MAXPOP*snp]>>1,snpcall[(int64_t ) parent2+(int64_t )MAXPOP*snp]&1); 
					else fprintf(output2,"%d %d ",snpcall[(int64_t ) parent2+(int64_t )MAXPOP*snp]&1,snpcall[(int64_t ) parent2+(int64_t )MAXPOP*snp]>>1); 	
				//	if (off==5) printf("%d %d ",snpcall[(int64_t ) parent2+(int64_t )MAXPOP*snp]&1,snpcall[(int64_t ) parent2+(int64_t )MAXPOP*snp]>>1); 	 			
				}
			} else if (famstrt[fam].p2==-1)
			{	if (transmit[fam][snp]&0b00001000)
				{	fprintf(output2,"NA NA "); 			
				} else if (transmit[fam][snp]&0b00000010)
				{	if (phaseparent1==1) fprintf(output2,"%d %d ",snpcall[(int64_t ) off+(int64_t )MAXPOP*snp]>>1,snpcall[(int64_t ) parent1+(int64_t )MAXPOP*snp]&1); 
					else fprintf(output2,"%d %d ",snpcall[(int64_t ) off+(int64_t )MAXPOP*snp]&1,snpcall[(int64_t ) parent1+(int64_t )MAXPOP*snp]&1); 			
				} else if ( transmit[fam][snp]&0b00000100 )
				{	fprintf(output2,"NA NA ");
				}	else 			
				{	if (transmit[fam][snp]&0b00000001) fprintf(output2,"%d %d ",snpcall[(int64_t ) parent1+(int64_t )MAXPOP*snp]>>1,snpcall[(int64_t ) parent1+(int64_t )MAXPOP*snp]&1); 
					else fprintf(output2,"%d %d ",snpcall[(int64_t ) parent1+(int64_t )MAXPOP*snp]&1,snpcall[(int64_t ) parent1+(int64_t )MAXPOP*snp]>>1); 
				};
				if (phaseparent1==2) fprintf(output2,"%d NA ",snpcall[(int64_t ) off+(int64_t )MAXPOP*snp]>>1); 
				else fprintf(output2,"%d NA ",snpcall[(int64_t ) off+(int64_t )MAXPOP*snp]&1); 			
				
			};			
				
				
		};
	
		fprintf(output2,"\n");
	};
	fclose(output2);

	for (int fam=0;fam<nbfam;fam++)//MAXFAM
	{	int off=famstrt[fam].off;
		int parent1=(famstrt[fam].p1>-1)?famstrt[fam].p1:-1;
		int parent2=(famstrt[fam].p2>-1)?famstrt[fam].p2:-1;
		//printf("%d %d %d\n",ID,IDp1loop,IDp2loop);
		
		for(int snp=0;snp<nbsnpinfile;snp++)		
		{//	printf("%d %d %d %d \n",snp,snpcall[(int64_t ) off+MAXPOP*snp],snpcall[(int64_t ) parent1+MAXPOP*snp],snpcall[(int64_t ) parent2+MAXPOP*snp]);
			//	printf("%d %d %d %d %d %d\n",dephasefound,anotherdephasefound,atleastonephasefind,categoryphaseeroor[0][0],categoryphaseeroor[0][1],categoryphaseeroor[0][2]);
			unsigned char marker=snpcall[(int64_t ) off+(int64_t )MAXPOP*snp];//(*((genomes+(unsigned long long) ID*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)>>(((snp%4)*2)))&3;
			unsigned char markerp1=(famstrt[fam].p1>-1)?snpcall[(int64_t ) parent1+(int64_t )MAXPOP*snp]:snpcall[(int64_t ) off+(int64_t )MAXPOP*snp];//(*((genomes+(unsigned long long) IDp1loop*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)>>(((snp%4)*2)))&3;
			unsigned char markerp2=(famstrt[fam].p2>-1)?snpcall[(int64_t ) parent2+(int64_t )MAXPOP*snp]:snpcall[(int64_t ) off+(int64_t )MAXPOP*snp];//(*((genomes+(unsigned long long) IDp2loop*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)>>(((snp%4)*2)))&3;
		//	printf("%d %d %d\n",marker,markerp1,markerp2);
			int snperrofind=0;
		/*	if (marker==4 || markerp2==4 || markerp1==4) snperrofind=1;
			else 
			{	if (marker==0 || marker==3)
				{	if (markerp1==3-marker) 
					{	snperrofind=1;
					} else if (markerp2==3-marker) 
					{	snperrofind=1;
					};
				} else  
				{	if (markerp1==3 && markerp2==3) 
					{	snperrofind=1;
					};
					if (markerp1==0 && markerp2==0) 
					{	snperrofind=1;
					};
				};
			};	
			if (snperrofind==0)
			{	int phase=-1;
				if (marker!=0 && marker!=3)
				{	if (markerp1==0) phase=marker;
					else if (markerp1==3) phase=3-marker;	
					else if (markerp2==0) phase=3-marker;
					else if (markerp2==3) phase=marker;	
				};
				//	 printf("snp %d: %d %d %d\n",snp,phaseparent1temp,phaseparent1,phase);
				if (phaseparent1==-1) //initiate phase
				{	if (phase>-1) 
					{	phaseparent1=phase;
						lastphase=phase;
					};
				} else 
				{	if (phase==-1) 
					{	if (marker!=0 && marker!=3 && phaseparent1!=lastphase)
						{	snpcall[(int64_t ) off+(int64_t )MAXPOP*snp]=3-snpcall[(int64_t ) off+(int64_t )MAXPOP*snp];
						}
					} else 
					{	lastphase=phase;
						if (phaseparent1!=phase) //phase detected not previous phase
						{	//	printf("%d %d %d %d %d %d %d\n ",offspring,snp,marker,allele1run,allele2run,phaseparent1,phase);
							snpcall[(int64_t ) off+(int64_t )MAXPOP*snp]=3-snpcall[(int64_t ) off+(int64_t )MAXPOP*snp];
					//		printf("%d %d %d %d %d %d %d\n ",offspring,snp,marker,allele1run,allele2run,phaseparent1,phase);
						};		
					};
				//printf("%d s %d se %d 1 %d 0 %d P %d r %d w %d\n",chrtemp1, snp,nbsnperror,nbone,nbzero,phaseparent1,right[chrtemp1],wrong[chrtemp1]);
				};	
			} else 
			{	MEerrorfind[snp][fam]=1;
			};*/
			phaseparent1=1;
			/*for(int offspring=0;offspring<MAXNBPO;offspring++)//
			{	int countunphase=0;
				int ID=UKBID[tabPO[offspring].IDoffspring]-1;
				int IDp1loop=UKBID[tabPO[offspring].IDp1]-1;
			//	printf("%d %d\n",ID,IDp1loop);
				int phaseparent1=-1;
				int phaseparent1temp=-1;
				int phasewell=1;
				int lastphase=-1;
				int dephasefound=0;
				int anotherdephasefound=0;
				int atleastonephasefind=0;
				for(int snp=0;snp<nbsnpperchr[chrtemp1];snp++)		
				{	unsigned char marker=(*((genomes+(unsigned long long) ID*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)>>(((snp%4)*2)))&3;
					unsigned char markerp1=(*((genomes+(unsigned long long) IDp1loop*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)>>(((snp%4)*2)))&3;
					unsigned char snperrofind=0;
					if (marker==0 || marker==3)
					{	if (markerp1==3-marker) 
						{	snperrofind=1;
							nbSNPerror[snp][chrtemp1]++;
						};
					};
					if (snperrofind==0)
					{	int phase=-1;
						if (marker!=0 && marker!=3)
						{	nbsitePO++;
							if (markerp1==0) phase=marker;
							else if (markerp1==3) phase=3-marker;	
						};
					//	 printf("snp %d: %d %d %d\n",snp,phaseparent1temp,phaseparent1,phase);
						if (phaseparent1==-1) //initiate phase
						{	if (phase>-1) 
							{	phaseparent1=phase;
								lastphase=phase;
							};
						} else 
						{	if (phase==-1) 
							{	if (marker!=0 && marker!=3 && dephasefound)
								{	dephasefound--;
									if (dephasefound==0)
									{	if (atleastonephasefind)
										{	if (anotherdephasefound)	
											{	categoryphaseeroor[1][1]++;
												dephasefound=0;
												atleastonephasefind=0;
												anotherdephasefound=0;
											} else 
											{	categoryphaseeroor[1][0]++;
												dephasefound=0;
												atleastonephasefind=0;
												anotherdephasefound=0;
											}	
										} else 
										{	categoryphaseeroor[1][2]++;
											dephasefound=0;
											atleastonephasefind=0;
											anotherdephasefound=0;	
										};
									};
								}
								if (marker!=0 && marker!=3 && phaseparent1!=lastphase)
								{	unsigned char allele1run=marker&1;
									unsigned char allele2run=marker>>1;
									*( genomes+(unsigned long long) (((((unsigned long long) ID)*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)))=
										(*( genomes+(unsigned long long) (((((unsigned long long) ID)*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4))) & 
											(~(3<<((snp%4)*2))))
												|	((2-allele2run*2+1-allele1run)<<((snp%4)*2));	
								}
							} else 
							{	if (lastphase!=phase) 
								{	nbphaseerror[snp][chrtemp1]++;
									if (dephasefound) 
									{	anotherdephasefound=1;
										atleastonephasefind=1;
									} else 
									{	atleastonephasefind=0; 
										anotherdephasefound=0;
									}
									dephasefound=2;
								} else 
								{	if (dephasefound)
									{	atleastonephasefind=1;
										dephasefound--;
										if (dephasefound==0)
										{	if (atleastonephasefind)
											{	if (anotherdephasefound)	
												{	categoryphaseeroor[1][1]++;
													dephasefound=0;
													atleastonephasefind=0;
													anotherdephasefound=0;
												} else 
												{	categoryphaseeroor[1][0]++;
													dephasefound=0;
													atleastonephasefind=0;
													anotherdephasefound=0;
												}	
											} else 
											{	categoryphaseeroor[1][2]++;
												dephasefound=0;
												atleastonephasefind=0;
												anotherdephasefound=0;	
											};
										};
									};
								};
								lastphase=phase;
								if (phaseparent1!=phase) //phase detected not previous phase
								{	unsigned char allele1run=marker&1;
									unsigned char allele2run=marker>>1;
								//	printf("%d %d %d %d %d %d %d\n ",offspring,snp,marker,allele1run,allele2run,phaseparent1,phase);
									*( genomes+(unsigned long long) (((((unsigned long long) ID)*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)))=
										(*( genomes+(unsigned long long) (((((unsigned long long) ID)*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4))) & 
											(~(3<<((snp%4)*2))))
												|	((2-allele2run*2+1-allele1run)<<((snp%4)*2));	
									 marker=(*((genomes+(unsigned long long) ID*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)>>(((snp%4)*2)))&3;
							//		printf("%d %d %d %d %d %d %d\n ",offspring,snp,marker,allele1run,allele2run,phaseparent1,phase);
								};		
							};
						};
					};
					
					
					//printf("%d s %d se %d 1 %d 0 %d P %d r %d w %d\n",chrtemp1, snp,nbsnperror,nbone,nbzero,phaseparent1,right[chrtemp1],wrong[chrtemp1]);
				};
			//	printf("%d %d\n",offspring,countunphase);
				
			};*/
			
		//	writegenome(chrtemp1,run,1000);
			
		}
		int transmittedfromparent1=-1;
		int transmittedfromparent2=-1;
		int lasthetsnp1=-1;
		int lasthetsnp2=-1;
		int lasthetsnp1mayo=-1;
		int lasthetsnp2mayo=-1;
		
		for(int snp=0;snp<nbsnpinfile;snp++)		
		{	unsigned char marker=snpcall[(int64_t ) off+(int64_t )MAXPOP*snp];//(*((genomes+(unsigned long long) ID*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)>>(((snp%4)*2)))&3;
			unsigned char markerp1=(famstrt[fam].p1>-1)?snpcall[(int64_t ) parent1+(int64_t )MAXPOP*snp]:snpcall[(int64_t ) off+(int64_t )MAXPOP*snp];//(*((genomes+(unsigned long long) IDp1loop*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)>>(((snp%4)*2)))&3;
			unsigned char markerp2=(famstrt[fam].p2>-1)?snpcall[(int64_t ) parent2+(int64_t )MAXPOP*snp]:snpcall[(int64_t ) off+(int64_t )MAXPOP*snp];//(*((genomes+(unsigned long long) IDp2loop*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)>>(((snp%4)*2)))&3;
			if (transmittedfromparent1==-1)
			{	if ((phaseparent1==1 && (marker>>1)==(markerp1>>1) ) || (phaseparent1==2 && (marker&1)==(markerp1>>1) ) )
				{	transmittedfromparent1=1;
				} else if  ((phaseparent1==1 && (marker>>1)==(markerp1&1)) || (phaseparent1==2 && (marker&1)==(markerp1&1)) )
				{	transmittedfromparent1=2;
				} else 
				{	transmittedfromparent1=5;
				};

			} else 
			{	if ( markerp1==1 || markerp1==2 )
				{	if ((phaseparent1==1 && (marker>>1)==(markerp1>>1)) || (phaseparent1==2 && (marker&1)==(markerp1>>1)) )
					{	if (transmittedfromparent1!=1)
						{	//	if (lasthetsnp1mayo!=lasthetsnp1) printf("%d %d %d %d %d\n",tabtrio[offspring].IDoffspring,tabtrio[offspring].IDp1,chrtemp1,snp,lasthetsnp1);
							transmittedfromparent1=1;
							lasthetsnp1mayo=snp;
						};
					} else if  ((phaseparent1==1 && (marker>>1)==(markerp1&1)) || (phaseparent1==2 && (marker&1)==(markerp1&1)) )
					{	if (transmittedfromparent1!=2)
						{//	if (lasthetsnp1mayo!=lasthetsnp1) printf("%d %d %d %d %d\n",tabtrio[offspring].IDoffspring,tabtrio[offspring].IDp1,chrtemp1,snp,lasthetsnp1);
							transmittedfromparent1=2;
							lasthetsnp1mayo=snp;
						};
					};
				} else 			
				{	if ((phaseparent1==1 && (marker>>1)!=(markerp1>>1)) || (phaseparent1==2 && (marker&1)!=(markerp1>>1)) )
					{	printf("ME snp %d %d\n",snp,transmittedfromparent1);
						if (transmittedfromparent1==1) transmittedfromparent1=3;
						if (transmittedfromparent1==2) transmittedfromparent1=4;
					}		
				};
			};
			lasthetsnp1=snp;
			int Merror=0;
			
			if (transmittedfromparent1==-1)
			{	transmit[fam][snp]=0b00000011;
			} else if (transmittedfromparent1==1)
			{	transmit[fam][snp]=0b00000001;
			} else  if (transmittedfromparent1==2)
			{	transmit[fam][snp]=0b00000000;
			} else if (transmittedfromparent1==3)
			{	transmit[fam][snp]=0b00000101;
				transmittedfromparent1=transmittedfromparent1-2;
			} else if (transmittedfromparent1==4)
			{	transmit[fam][snp]=0b00000100;
				transmittedfromparent1=transmittedfromparent1-2;
			} else if (transmittedfromparent1==5)
			{	transmit[fam][snp]=0b00001000;
				transmittedfromparent1=-1;
			};
			
			if (transmittedfromparent2==-1)
			{	if ((phaseparent1==1 && (marker&1)==(markerp2>>1)) || (phaseparent1==2 && (marker>>1)==(markerp2>>1)) )
				{	transmittedfromparent2=1;
				} else if ((phaseparent1==1 && (marker&1)==(markerp2&1)) || (phaseparent1==2 && (marker>>1)==(markerp2&1)) )
				{	transmittedfromparent2=2;
				} else 
				{	transmittedfromparent2=5;
				};
			} else 
			{	if (markerp2==1 || markerp2==2 )
				{	if ((phaseparent1==1 && (marker&1)==(markerp2>>1)) || (phaseparent1==2 && (marker>>1)==(markerp2>>1)) )
					{	if (transmittedfromparent2!=1)
						{	//	if (lasthetsnp2mayo!=lasthetsnp2) printf("%d %d %d %d %d\n",tabtrio[offspring].IDoffspring,tabtrio[offspring].IDp2,chrtemp1,snp,lasthetsnp2);
							transmittedfromparent2=1;
							lasthetsnp2mayo=snp;
						};
					} else if ((phaseparent1==1 && (marker&1)==(markerp2&1)) || (phaseparent1==2 && (marker>>1)==(markerp2&1)) )
					{	if (transmittedfromparent2!=2)
						{	//	if (lasthetsnp2mayo!=lasthetsnp2) printf("%d %d %d %d %d\n",tabtrio[offspring].IDoffspring,tabtrio[offspring].IDp2,chrtemp1,snp,lasthetsnp2);
							transmittedfromparent2=2;
							lasthetsnp2mayo=snp;
						};
					};
				}	else 
				{	if ((phaseparent1==2 && (marker>>1)!=(markerp2>>1)) || (phaseparent1==1 && (marker&1)!=(markerp2>>1)) )
					{	printf("ME snp %d %d\n",snp,transmittedfromparent2);
						if (transmittedfromparent2==1) transmittedfromparent2=3;
						if (transmittedfromparent2==2) transmittedfromparent2=4;
					}	
				};
			};
			lasthetsnp2=snp;
			if (off==5) printf("%d\n",transmit[fam][snp]);
			if (transmittedfromparent2==-1)
			{	transmit[fam][snp]=0b00100000 | transmit[fam][snp];
			} else if (transmittedfromparent2==1)
			{	transmit[fam][snp]=0b00010000 | transmit[fam][snp];
			} else if (transmittedfromparent2==2) 
			{	transmit[fam][snp]=0b00000000 | transmit[fam][snp];
			} else if (transmittedfromparent2==3)
			{	transmit[fam][snp]=0b01010000  | transmit[fam][snp];
				transmittedfromparent2=transmittedfromparent2-2;  
			} else if (transmittedfromparent2==4)
			{	transmit[fam][snp]=0b01000000  | transmit[fam][snp];
				transmittedfromparent2=transmittedfromparent2-2;
			} else if (transmittedfromparent2==5)
			{	transmit[fam][snp]=0b10000000  | transmit[fam][snp];
				transmittedfromparent2=-1;
			};
		//	if ( (marker==0 && (markerp1==3 || markerp2==3) ) || (marker==3 && (markerp1==0 || markerp2==0)) || ((marker==1 || marker==2) && ((markerp1==0 && markerp2==0) || (markerp1==3 && markerp2==3) )))
		//	{	transmit[fam][snp]=transmit[fam][snp] & 0b01000100;
		//	
		//	};
			if (transmit[fam][snp]==4)  printf("%d %d %d\n",snp ,fam,transmit[fam][snp]);
		//	printf("%d %d %d %d %d %d %d\n",chrtemp1,snp,marker,markerp1,markerp2,transmittedfromparent1,transmittedfromparent2);
		};
	};
	
	strcpy(resultfam,result);
	if ((output = fopen(resultfam, "w")) == NULL) 
	{	printf("file %s can not be opened\n",resultfam);		
		return (0);
	};
	fprintf(output,"CHR SNP A1 A2 ");
	 for (int fam=0;fam<nbfam;fam++)//MAXFAM
	{	int off=famstrt[fam].off;
		int parent1=(famstrt[fam].p1>-1)?famstrt[fam].p1:-1;
		int parent2=(famstrt[fam].p2>-1)?famstrt[fam].p2:-1;
		if (famstrt[fam].p1>-1 && famstrt[fam].p2>-1)
		{	fprintf(output,"%s-%s.%s-%s.T %s-%s.%s-%s.NT %s-%s.%s-%s.T %s-%s.%s-%s.NT ",
						namefamily[fam],nameindiv[famstrt[fam].off],namefamily[fam],nameindiv[famstrt[fam].p1],namefamily[fam],nameindiv[famstrt[fam].off],namefamily[fam],nameindiv[famstrt[fam].p1],
						namefamily[fam],nameindiv[famstrt[fam].off],namefamily[fam],nameindiv[famstrt[fam].p2],namefamily[fam],nameindiv[famstrt[fam].off],namefamily[fam],nameindiv[famstrt[fam].p2]);
		};
		if (famstrt[fam].p1==-1)
		{	fprintf(output,"%s-%s.%s-NA.T %s-%s.%s-NA.NT %s-%s.%s-%s.T %s-%s.%s-%s.NT ",
						namefamily[fam],nameindiv[famstrt[fam].off],namefamily[fam],namefamily[fam],nameindiv[famstrt[fam].off],namefamily[fam],
						namefamily[fam],nameindiv[famstrt[fam].off],namefamily[fam],nameindiv[famstrt[fam].p2],namefamily[fam],nameindiv[famstrt[fam].off],namefamily[fam],nameindiv[famstrt[fam].p2]);
		} else if (famstrt[fam].p2==-1)
		{	fprintf(output,"%s-%s.%s-%s.T %s-%s.%s-%s.NT %s-%s.%s-NA.T %s-%s.%s-NA.NT ",
						namefamily[fam],nameindiv[famstrt[fam].off],namefamily[fam],nameindiv[famstrt[fam].p1],namefamily[fam],nameindiv[famstrt[fam].off],namefamily[fam],nameindiv[famstrt[fam].p1],
						namefamily[fam],nameindiv[famstrt[fam].off],namefamily[fam],namefamily[fam],nameindiv[famstrt[fam].off],namefamily[fam]);
		}
	}
	fprintf(output,"\n");
	
	for(int snp=0;snp<nbsnpinfile;snp++)		
	{	fprintf(output,"%d %s %c %c ",chr,namesnp[snp],allele1[snp],allele2[snp]);
	
		for (int fam=0;fam<nbfam;fam++)//MAXFAM
		{			//printf("place %d ",place);
			
			int off=famstrt[fam].off;
			int parent1=(famstrt[fam].p1>-1)?famstrt[fam].p1:-1;
			int parent2=(famstrt[fam].p2>-1)?famstrt[fam].p2:-1;
		 if (off==3) printf("SNP %d %d %d  %d %d %d %d %d\n",fam,snp,off,parent1,parent2,transmit[fam][snp],snpcall[(int64_t ) parent1+(int64_t )MAXPOP*snp],snpcall[(int64_t ) off+(int64_t )MAXPOP*snp]);
		
			if (famstrt[fam].p1>-1)
			{	if (transmit[fam][snp]&0b00000010)
				{	if (phaseparent1==1) fprintf(output,"%d %d ",snpcall[(int64_t ) off+(int64_t )MAXPOP*snp]>>1,snpcall[(int64_t ) parent1+(int64_t )MAXPOP*snp]&1); 
					else fprintf(output,"%d %d ",snpcall[(int64_t ) off+(int64_t )MAXPOP*snp]&1,snpcall[(int64_t ) parent1+(int64_t )MAXPOP*snp]&1); 			
				} else if ((transmit[fam][snp]&0b00000100) || MEerrorfind[snp][fam]==1)
				{	if (phaseparent1==1) 
					{	if (flagPO[0]=='O')
						{	if ((snpcall[(int64_t ) off+(int64_t )MAXPOP*snp]&1)==0 )
							{	if (snpcall[(int64_t ) parent1+(int64_t )MAXPOP*snp]==0) fprintf(output,"%d %d ",0,0); 
								else fprintf(output,"%d %d ",0,1); 
							} else 
							{	if (snpcall[(int64_t ) parent1+(int64_t )MAXPOP*snp]==3) fprintf(output,"%d %d ",1,1); 
								else fprintf(output,"%d %d ",1,0); 
							}								
						} else 
						{	if (snpcall[(int64_t ) parent1+(int64_t )MAXPOP*snp]==0)
							{	fprintf(output,"%d %d ",0,0); 
								
							} else 
							{	
								fprintf(output,"%d %d ",1,1); 
							}								
						}
					}
					else 
					{	if (flagPO[0]=='O')
						{	if ((snpcall[(int64_t ) off+(int64_t )MAXPOP*snp]>>1)==0 )
							{	if (snpcall[(int64_t ) parent1+(int64_t )MAXPOP*snp]==0)	fprintf(output,"%d %d ",0,0); 
								else fprintf(output,"%d %d ",0,1); 
							} else 
							{	if (snpcall[(int64_t ) parent1+(int64_t )MAXPOP*snp]==3)	fprintf(output,"%d %d ",1,1); 
								else fprintf(output,"%d %d ",1,0); 
							}								
						} else 
						{	if (snpcall[(int64_t ) parent1+(int64_t )MAXPOP*snp]==0)
							{	fprintf(output,"%d %d ",0,0); 	
							} else 
							{	fprintf(output,"%d %d ",1,1); 
							}								
						}
					};
					
				} else 
				{	if (transmit[fam][snp]&0b00000001) fprintf(output,"%d %d ",snpcall[(int64_t ) parent1+(int64_t )MAXPOP*snp]>>1,snpcall[(int64_t ) parent1+(int64_t )MAXPOP*snp]&1); 
					else fprintf(output,"%d %d ",snpcall[(int64_t ) parent1+(int64_t )MAXPOP*snp]&1,snpcall[(int64_t ) parent1+(int64_t )MAXPOP*snp]>>1); 
				};
			} else if (famstrt[fam].p1==-1)
			{	if (phaseparent1==1) fprintf(output,"%d NA ",snpcall[(int64_t ) off+(int64_t )MAXPOP*snp]>>1); 
				else fprintf(output,"%d NA ",snpcall[(int64_t ) off+(int64_t )MAXPOP*snp]&1); 			
			};
				
			if (famstrt[fam].p2>-1)
			{	if (transmit[fam][snp]&0b00100000)
				{	if (phaseparent1==2) fprintf(output,"%d %d ",snpcall[(int64_t ) off+(int64_t )MAXPOP*snp]>>1,snpcall[(int64_t ) parent2+(int64_t )MAXPOP*snp]&1); 
					else fprintf(output,"%d %d ",snpcall[(int64_t ) off+(int64_t )MAXPOP*snp]&1,snpcall[(int64_t ) parent2+(int64_t )MAXPOP*snp]&1); 			
				} else if ((transmit[fam][snp]&0b01000000) || MEerrorfind[snp][fam]==1)
				{	if (phaseparent1==2) 
					{	if (flagPO[0]=='O')
						{	if ((snpcall[(int64_t ) off+(int64_t )MAXPOP*snp]&1)==0 )
							{	if (snpcall[(int64_t ) parent2+(int64_t )MAXPOP*snp]==0)	fprintf(output,"%d %d ",0,0); 
								else fprintf(output,"%d %d ",0,1); 
							} else 
							{	if (snpcall[(int64_t ) parent2+(int64_t )MAXPOP*snp]==3)	fprintf(output,"%d %d ",1,1); 
								else fprintf(output,"%d %d ",1,0); 
							};							
						} else 
						{	if (snpcall[(int64_t ) parent2+(int64_t )MAXPOP*snp]==0)
							{	fprintf(output,"%d %d ",0,0); 	
							} else 
							{	fprintf(output,"%d %d ",1,1); 
							}		
						}
					}
					else 
					{	if (flagPO[0]=='O')
						{	if ((snpcall[(int64_t ) off+(int64_t )MAXPOP*snp]>>1)==0 )
							{	if (snpcall[(int64_t ) parent2+(int64_t )MAXPOP*snp]==0)	fprintf(output,"%d %d ",0,0); 
								else fprintf(output,"%d %d ",0,1); 
							} else 
							{	if (snpcall[(int64_t ) parent2+(int64_t )MAXPOP*snp]==3)	fprintf(output,"%d %d ",1,1); 
								else fprintf(output,"%d %d ",1,0); 
							}	
						} else 
						{	if (snpcall[(int64_t ) parent2+(int64_t )MAXPOP*snp]==0)
							{	fprintf(output,"%d %d ",0,0); 	
							} else 
							{	fprintf(output,"%d %d ",1,1); 
							}		
						}						
					};
					
				} else 
				{	if (transmit[fam][snp]&0b00010000)
					fprintf(output,"%d %d ",snpcall[(int64_t ) parent2+(int64_t )MAXPOP*snp]>>1,snpcall[(int64_t ) parent2+(int64_t )MAXPOP*snp]&1); 
					else fprintf(output,"%d %d ",snpcall[(int64_t ) parent2+(int64_t )MAXPOP*snp]&1,snpcall[(int64_t ) parent2+(int64_t )MAXPOP*snp]>>1); 
				}
			} else if (famstrt[fam].p2==-1)	
			{	if (phaseparent1==2) fprintf(output,"%d NA ",snpcall[(int64_t ) off+(int64_t )MAXPOP*snp]>>1); 
				else fprintf(output,"%d NA ",snpcall[(int64_t ) off+(int64_t )MAXPOP*snp]&1); 			
				
			};
			if (transmit[fam][snp]==4) printf("%d %d %d %d %d\n",snp,fam ,transmit[fam][snp],famstrt[fam].p1,famstrt[fam].p2);
		};
		fprintf(output,"\n");
	};
	fclose(output);
	printf("executable run normally");
	/*
	if ((MAFfile = fopen(hapfile, "r")) == NULL) 
	{	if (printdetail) printf("file triofile is not found\n");		
		return (1);
	};
	FILE * output;
	if ((output = fopen(result, "w")) == NULL) 
	{	if (printdetail) printf("file triofile is not found\n");		
		return (1);
	};
	do {first=getc(MAFfile);fprintf(output,"%c",first);} while (first!='\n' ); 	
	do {first=getc(MAFfile);fprintf(output,"%c",first);} while (first!='\n'); 	
	do {first=getc(MAFfile);fprintf(output,"%c",first);} while (first!='\n'); 	
	for(int snp=0;snp<nbsnpinfile;snp++)		
	{	printf("%d\n",snp);
		do {first=getc(MAFfile);fprintf(output,"%c",first);} while (first!=9 && first!=EOF); 	
		do {first=getc(MAFfile);fprintf(output,"%c",first);} while (first!=9 && first!=EOF); 	
		do {first=getc(MAFfile);fprintf(output,"%c",first);} while (first!=9 && first!=EOF); 	//229133-223882

		do {first=getc(MAFfile);fprintf(output,"%c",first);} while (first!=9 && first!=EOF); 	
		do {first=getc(MAFfile);fprintf(output,"%c",first);} while (first!=9 && first!=EOF); 	
		do {first=getc(MAFfile);fprintf(output,"%c",first);} while (first!=9 && first!=EOF);
		
		do {first=getc(MAFfile);fprintf(output,"%c",first);} while (first!=9 && first!=EOF);
		do {first=getc(MAFfile);fprintf(output,"%c",first);} while (first!=9 && first!=EOF);
		do {first=getc(MAFfile);fprintf(output,"%c",first);} while (first!=9 && first!=EOF);	
		
		if (first!=EOF) 	//22
		{	//printf("place %d ",place); 
			fprintf(output,"%d",snpcall[(int64_t ) off][snp]>>1);
		//	if (snp<2) printf("a%c ",c1);
			fprintf(output,"/");
			fprintf(output,"%d ",snpcall[(int64_t ) off][snp]%2);
			fprintf(output,"%d",snpcall[(int64_t ) parent1][snp]>>1);
		//	if (snp<2) printf("a%c ",c1);
			fprintf(output,"/");
			fprintf(output,"%d ",snpcall[(int64_t ) parent1][snp]%2);
			fprintf(output,"%d",snpcall[(int64_t ) parent2][snp]>>1);
		//	if (snp<2) printf("a%c ",c1);
			fprintf(output,"/");
			fprintf(output,"%d\n",snpcall[(int64_t ) parent2][snp]%2);			
			do {first=getc(MAFfile);} while (first!='\n' && first!=EOF); 	
		};
	};
	fclose(MAFfile);
	fclose(output);
	*/
		/*FILE * errors;
	if ((errors = fopen("/pl/active/KellerLab/Emmanuel/gameticphasing/errors.txt", "w")) == NULL) 
	{	if (printdetail) printf("file triofile is not found\n");		
		return (1);
	};
	for(int  chrtemp1=1;chrtemp1<22+1;chrtemp1++)		
	{	for(int snp=0;snp<nbsnpperchr[chrtemp1];snp++)		
		{	fprintf(errors,"%d %d %d %d\n",chrtemp1,snp,nbSNPerror[snp][chrtemp1],nbphaseerror[snp][chrtemp1]);
		};
	};
	fclose(	errors);*/
//	compareressultindivallcombin( 1,chrtodo1,  chr2);
//	compareressultindivallcombin( 2,chrtodo1,  chr2);

	//compareressult(chrtodo1,chr2);
//	compareressultalea(chrtodo1,chr2);
//	compareressult3(chrtodo1,chr2,20);
//	compareressultn(chrtodo1,chr2);
	return 0;
}
