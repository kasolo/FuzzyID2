#include "main.h"

ofstream logFile;
ofstream outFile;
sqlite3 *conn;
char *dbErrMsg;

int main(int argc, char *argv[])
{
    InData *optIn;
    optIn = (InData *)malloc(sizeof(InData));

    init(argc, argv, optIn);
    logFile << "option check pass" << endl;

    if(strcmp(optIn->command, "ID")==0)
    {
        logFile << "Function ID start" << endl;
        cout<<"ID st"<<endl;
        ID(optIn);
        logFile << "Function ID end" << endl;
    }
    else if(strcmp(optIn->command, "Theta1")==0)
    {
        logFile << "Function Theta1 start" << endl;
        Theta1(optIn);
        logFile << "Function Theta1 end" << endl;
    }
    else if(strcmp(optIn->command, "Theta2")==0)
    {
        logFile << "Function Theta2 start" << endl;
        Theta2(optIn);
        logFile << "Function Theta2 end" << endl;
    }
    //close files
	logFile.close();
	outFile.close();
    return 0;
}

