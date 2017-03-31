#ifndef MAIN_H
#define MAIN_H

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
using namespace std;
#include "sqlite3.h"
#include <Bpp/Seq/Alphabet.all> /* this includes all alphabets in one shot */
#include <Bpp/Seq/Container.all> /* this includes all containers */
#include <Bpp/Seq/AlphabetIndex/DefaultNucleotideScore.h>
#include <Bpp/Seq/Io.all> /* this includes all sequence readers and writers */
#include <Bpp/Phyl/Model.all> /* this includes all models */

#include <Bpp/Phyl/Distance.all>

#include <Bpp/Numeric/Prob.all> /* include all probability distributions */
#include <Bpp/Numeric.all>
#include <Bpp/Numeric/Matrix.all>

using namespace bpp;

typedef struct {
    char model[10];
    char command[10];
    char inFile[100];
    char outFile[100];
    char dataBase[100];
    char roughBarcode[100];
    char mdBarcode[100];
}InData;//input paramaters

typedef struct {
	char model[10];
    double md;
    double theta1;
    double theta2;
	char dataBase[100];
    char roughBarcode[100];
    char mdBarcode[100];
    char taxoName[100];
    char roughResult[3][100];
	char roughFamily[100];
	char roughGenus[100];
	char speciesName[100];
}mdResult;//compute MD

typedef struct {
	char model[10];
    double theta1;
    char genusName[100];
	char speciesName[100];
    char mdBarcode[100];
}theta1Result;//compute theta1

typedef struct {
	char model[10];
    double theta2;
    char taxoName[100];
    char genusName[100];
	char speciesName[100];
    char mdBarcode[100];
}theta2Result;//compute theta2

typedef struct {
    char host[20];
    char userName[20];
    char passWord[20];
    int port;
}databaseManage;//db connection paramaters

extern ofstream logFile;//log file
extern ofstream outFile;//ID result file

extern sqlite3 *conn;//database connection
extern char *dbErrMsg;

//query sequence name repeative check
extern string checkQuerySeqName(char * filename);

//string compare ignore case
extern int strcmpIgnoreCase(char *str1, char *str2);
//find string ignore case
extern int findIgnoreCase(char *str1, char *str2);
//readFileIntoString
extern string readFileIntoString(char * filename);
//split string
extern void split( char **arr, char *str, const char *del);

//init get input option
extern void init(int argc, char **argv, InData *optIn);

//estimate pairwise distance
extern double pairwiseDistance(Sequence &seq1,Sequence &seq2,string model);

//compute MF
extern double calMF(double x, double theta1, double theta2);

//main function
//ID
extern void ID(InData *optIn);
//Theta1
extern void Theta1(InData *optIn);
//Theta2
extern void Theta2(InData *optIn);

//roughGenus
extern int roughGenus(Sequence *seqQuery, mdResult *mdr);

//getMD
extern int getMD(Sequence *seqQuery, mdResult *mdr);

//getTheta1Sequences
extern void getTheta1function(theta1Result *t1r);
extern string getTheta1AfterMD(mdResult *mdr);
//get theta1
extern void getTheta1(VectorSequenceContainer *sequences, theta1Result *t1r);

//update theta1
extern void updateTheta1(theta1Result *t1r);

extern void updateSingletonTheta1(theta1Result *t1r);

//getTeta2Sequences
extern int getTheta2function(theta2Result *t2r);
extern string getTheta2AfterMD(mdResult *mdr);

//get min distance taxo
extern int getMinTaxo(VectorSequenceContainer *sequences, int i, theta2Result *t2r);

extern void getTaxoSequences(VectorSequenceContainer *sequences, const char *query, theta2Result *t2r);

extern void getTheta2(VectorSequenceContainer *sequencesA, VectorSequenceContainer *sequencesB, theta2Result *t2r);
extern void updateTheta2(theta2Result *t2r);

extern void updateSingletonTheta2(theta2Result *t2r);
#endif // MAIN_H
