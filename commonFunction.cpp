#include "main.h"

//init get input option
void init(int argc, char **argv, InData *optIn)
{
    char tmpArg[10];
    char tmpOpt[10];

    //get options
    for(int i=1;i<argc;i=i+2)
    {
        strcpy(tmpArg,argv[i]);

        if(strcmp(tmpArg, "-h")==0)
        {
            //read help file
            string help = readFileIntoString("HELP");
            cout << help << endl;

            exit(0);
        }
        else if(strcmp(tmpArg, "-m")==0)
        {
            strcpy(tmpOpt, argv[i+1]);
            if(strcmp(tmpOpt, "K2P")==0 || strcmp(tmpOpt, "JC69")==0 || strcmp(tmpOpt, "GTR")==0)
            {
                strcpy(optIn->model, tmpOpt);
                cout << "Model:" << optIn->model << endl;
            }
            else
            {
                cout << "Error: Model " << tmpOpt << " does not exist, please use JC69, K2P or GTR models"<< endl;
                exit(EXIT_FAILURE);
            }
        }
        else if(strcmp(tmpArg, "-in")==0)
        {
            strcpy(optIn->inFile , argv[i+1]);
            cout << "Query file:" <<optIn->inFile << endl;
        }
        else if(strcmp(tmpArg, "-out")==0)
        {
            strcpy(optIn->outFile , argv[i+1]);
            cout << "Output file:" <<optIn->outFile << endl;
        }
        else if(strcmp(tmpArg, "-c")==0)
        {
            strcpy(tmpOpt, argv[i+1]);
            if(strcmp(tmpOpt, "ID")==0 || strcmp(tmpOpt, "Theta1")==0 || strcmp(tmpOpt, "Theta2")==0)
            {
                strcpy(optIn->command, tmpOpt);
                cout << "Command:" <<optIn->command << endl;
            }
            else
            {
                cout << "Error: command " << tmpOpt << " does not correct, please use ID, Theta1 or Theta2 commands" << endl;
                exit(EXIT_FAILURE);
            }
        }
        else if(strcmp(tmpArg, "-d")==0)
        {
            strcpy(optIn->dataBase , argv[i+1]);
            cout << "Database name:" <<optIn->dataBase << endl;
        }
        else if(strcmp(tmpArg, "-rb")==0)
        {
            string tmp=argv[i+1];
            transform(tmp.begin(), tmp.end(), tmp.begin(), towupper);  
            sprintf(optIn->roughBarcode,"%s",tmp.c_str());
            //strcpy(optIn->roughBarcode , argv[i+1]);
            cout << "Rough barcode:" <<optIn->roughBarcode << endl;
        }
        else if(strcmp(tmpArg, "-mb")==0)
        {
            string tmp=argv[i+1];
            transform(tmp.begin(), tmp.end(), tmp.begin(), towupper);  
            sprintf(optIn->mdBarcode,"%s",tmp.c_str());
            //strcpy(optIn->mdBarcode , argv[i+1]);
            cout << "MD barcode:" <<optIn->mdBarcode << endl;
        }
        else
        {
            cout << "Error, illegal option:" << tmpArg << " Please check the option name" << endl;
            exit(EXIT_FAILURE);
        }
    }//for

    //check options
    if(strcmp(optIn->command, "")==0)
    {
        cout << "Error: please specify a command for this run with -c.(Theta1, Theta2 or ID)" << endl;
        exit(EXIT_FAILURE);
    }
    //if the command is ID, you must specify inFile,outFile,model,database
    else if(strcmp(optIn->command, "ID")==0)
    {
        if(strcmp(optIn->inFile, "")==0)
        {
            cout << "Error: please specify an input query file for this run with -in.(query sequences in fasta file)" << endl;
            exit(EXIT_FAILURE);
        }
        else
        {
            char *myArray[3];
            char tmpFile[100];
            strcpy(tmpFile , optIn->inFile);
            memset(myArray, 0x0, sizeof(myArray));
            split(myArray, tmpFile, "+");
            for(int i=0;i<3;i++){
                if(myArray[i]!='\0'){
                    if(findIgnoreCase(myArray[i],optIn->roughBarcode)!=std::string::npos){
                        fstream _file;
                        _file.open(myArray[i], std::ios_base::in);
                        if(!_file)
                        {
                            cout << "Error: the file " << myArray[i] << " you want to open for reading does not exist." << endl;
                            exit(EXIT_FAILURE);
                        }
                        else
                        {
                            //query sequence name repeative check
                            string p=checkQuerySeqName(myArray[i]);
                            if(strcmp(p.c_str(), "")!=0)
                            {
                                cout<< "Query sequence repeat:"<< p << endl;
                                exit(EXIT_FAILURE);
                            }
                        }
                    }
                    if(findIgnoreCase(myArray[i],optIn->mdBarcode)!=std::string::npos){
                        fstream _file;
                        _file.open(myArray[i], std::ios_base::in);
                        if(!_file)
                        {
                            cout << "Error: the file " << myArray[i] << " you want to open for reading does not exist." << endl;
                            exit(EXIT_FAILURE);
                        }
                        else
                        {
                            //query sequence name repeative check
                            string p=checkQuerySeqName(myArray[i]);
                            if(strcmp(p.c_str(), "")!=0)
                            {
                                cout<< "Query sequence repeat:"<< p << endl;
                                exit(EXIT_FAILURE);
                            }
                        }
                    }
                }
            }
        }
        if(strcmp(optIn->model, "")==0)
        {
            cout << "Error: please specify a model for this run with -m.(JC69, K2P or GTR)" << endl;
            exit(EXIT_FAILURE);
        }
        if(strcmp(optIn->outFile, "")==0)
        {
            cout << "Error: please specify an output file for this run with -out." << endl;
            exit(EXIT_FAILURE);
        }
        if(strcmp(optIn->dataBase, "")==0)
        {
            cout << "Error: please specify an database for this run with -d." << endl;
            exit(EXIT_FAILURE);
        }
    }
    //if the command is Theta1 or Theta2, you must specify the database
    else if(strcmp(optIn->command, "Theta1")==0 || strcmp(optIn->command, "Theta2")==0)
    {
        if(strcmp(optIn->dataBase, "")==0)
        {
            cout << "Error: please specify an database for this run with -d." << endl;
            exit(EXIT_FAILURE);
        }
        if(strcmp(optIn->model, "")==0)
        {
            cout << "Error: please specify a model for this run with -m.(JC69, K2P or GTR)" << endl;
            exit(EXIT_FAILURE);
        }
        if(strcmp(optIn->mdBarcode, "")==0)
        {
            cout << "Error: please specify the marker for species identification." << endl;
            exit(EXIT_FAILURE);
        }
    }
    //open log file
    logFile.open("./log");
    //open result file
    outFile.open(optIn->outFile);

    //open database
	int rc;
	char dbName[100];
    char hmmDBName[100];
	sprintf(dbName,"./DB/%s",optIn->dataBase);
	rc = sqlite3_open(dbName, &conn); //open db
	if(rc)
	{
		cout << "Can't open database:" << sqlite3_errmsg(conn) << endl;
	    sqlite3_close(conn);
	    exit(EXIT_FAILURE);
	}
	//check HmmDB existence
    if(strcmp(optIn->command, "ID")==0){
        sprintf(hmmDBName,"./HmmDB/%s_%s",optIn->dataBase,optIn->roughBarcode);
        std::ifstream f(hmmDBName);
        if(f){
            f.close();
        }
        else{
            cout << "Can't find HmmDB:" << optIn->dataBase<<"_"<<optIn->roughBarcode << endl;
            exit(EXIT_FAILURE);
        }
    }
    logFile << "db connect ok" << endl;
}

int strcmpIgnoreCase(char *str1, char *str2)
{
    string s1=str1;
    string s2=str2;
    transform(s1.begin(), s1.end(), s1.begin(), towupper);  
    transform(s2.begin(), s2.end(), s2.begin(), towupper); 
    return s1==s2;
}
int findIgnoreCase(char *str1, char *str2)
{
    string s1=str1;
    string s2=str2;
    transform(s1.begin(), s1.end(), s1.begin(), towupper);  
    transform(s2.begin(), s2.end(), s2.begin(), towupper); 
    return s1.find(s2);
}

void split( char **arr, char *str, const char *del)//split string
{ 
	char *s =NULL;  
	s=strtok(str,del); 
	while(s != NULL) 
	{  
		*arr++ = s; 
		s = strtok(NULL,del); 
	}
}

//从文件读入到string里
string readFileIntoString(char * filename)
{
	ifstream ifile(filename);
	//将文件读入到ostringstream对象buf中
	ostringstream buf;
	char ch;
	while(buf&&ifile.get(ch))
		buf.put(ch);
	//返回与流对象buf关联的字符串
	return buf.str();
}

//trim
string& LeftTrim(string &str)
{
    string::iterator iter=find_if(str.begin(),str.end(),not1(ptr_fun<int>(::isspace)));
    str.erase(str.begin(),iter);
    return str;
}
string& RightTrim(string &str)
{
    string::reverse_iterator rev_iter=find_if(str.rbegin(),str.rend(),not1(ptr_fun<int>(::isspace)));
    str.erase(rev_iter.base(),str.end());
    return str;
}
string& trim(string &str)
{
     return LeftTrim(RightTrim(str));
}

//query sequence name repeative check
string checkQuerySeqName(char * filename)
{
    string a;
    string seqNames;
    string::size_type ps;
	ifstream infile;
    infile.open(filename,ios::in);
    while(!infile.eof())
    {
        getline(infile, a, '\n');//读取一行，以换行符结束，存入 a中
        if(a[0]=='>'){
            ps=seqNames.find(a);
            if (ps!=std::string::npos)
            {
                return a;
            }
            else{
                seqNames=seqNames+a;
            }
        }
    }
    return "";
}

//estimate pairwise distance
double pairwiseDistance(const Sequence &seq1,const Sequence &seq2, string model)
{
    SubstitutionModel *evolutionModel;
    if(model == "K2P")
    {
        evolutionModel = new K80(&AlphabetTools::DNA_ALPHABET, 2.5);
    }
    else if(model == "JC69")
    {
        evolutionModel = new JCnuc(&AlphabetTools::DNA_ALPHABET);
    }
    else if(model == "GTR")
    {
        evolutionModel = new GTR(&AlphabetTools::DNA_ALPHABET);
    }
    //add the two sequences to a container
    VectorSequenceContainer* vsc = new VectorSequenceContainer(&AlphabetTools::DNA_ALPHABET);
    vsc->addSequence(seq1);
    vsc->addSequence(seq2);

    //alignment
    SiteContainer *sites = SiteContainerTools::alignNW(vsc->getSequence(0),vsc->getSequence(1),DefaultNucleotideScore(&AlphabetTools::DNA_ALPHABET), -5);
    //delete gaps
    SiteContainer * completeSites = SiteContainerTools::getSitesWithoutGaps(* sites);

    //estimate distance
    DiscreteDistribution *rateDist = new ConstantDistribution(1.);
    
    DistanceEstimation distanceMethod(evolutionModel, rateDist, completeSites);
    DistanceMatrix *distances = distanceMethod.getMatrix();

    return distances->row(1)[0];
}

//calculate MF
double calMF(double x, double theta1, double theta2)
{
    if(x <= theta1)
    {
        return 1.0;
    }
    else if(theta1<=x && x<=(theta1+theta2)/2)
    {
        return 1-2*((x-theta1)/(theta2-theta1))*((x-theta1)/(theta2-theta1));
    }
    else if((theta1+theta2)/2<=x && x<=theta2)
    {
       return  2*((x-theta2)/(theta2-theta1))*((x-theta2)/(theta2-theta1));
    }
    else if(x>=theta2)
    {
        return 0;
    }
    else{
        return 0;
    }
}
void getRoughResult(string &str, mdResult *mdr)
{
	string::size_type ps,pe;
    char *myArrayTaxo[3];
    char tmpChar[100];

	ps = str.find("\n");
	string tmpStr=str.substr(0,ps);
    string resultRemain=str.substr(ps+2);

	trim(tmpStr);
	ps=tmpStr.find_last_of("  ");
	pe=tmpStr.find("\n");
	tmpStr=tmpStr.substr(ps+1,pe-ps);
    memset(myArrayTaxo, 0x0, sizeof(myArrayTaxo));
    sprintf(tmpChar,"%s",tmpStr.c_str());
    split(myArrayTaxo, tmpChar, "_");
    sprintf(mdr->roughResult[0],"%s_%s",myArrayTaxo[1],myArrayTaxo[2]);

    ps = resultRemain.find("\n");
    tmpStr=resultRemain.substr(0,ps);
    resultRemain=resultRemain.substr(ps+2);
	trim(tmpStr);
	ps=tmpStr.find_last_of("  ");
	pe=tmpStr.find("\n");
	tmpStr=tmpStr.substr(ps+1,pe-ps);
    memset(myArrayTaxo, 0x0, sizeof(myArrayTaxo));
    sprintf(tmpChar,"%s",tmpStr.c_str());
    split(myArrayTaxo, tmpChar, "_");
    sprintf(mdr->roughResult[1],"%s_%s",myArrayTaxo[1],myArrayTaxo[2]);

    ps = resultRemain.find("\n");
    tmpStr=resultRemain.substr(0,ps);
    resultRemain=resultRemain.substr(ps+2);
	trim(tmpStr);
	ps=tmpStr.find_last_of("  ");
	pe=tmpStr.find("\n");
	tmpStr=tmpStr.substr(ps+1,pe-ps);
    memset(myArrayTaxo, 0x0, sizeof(myArrayTaxo));
    sprintf(tmpChar,"%s",tmpStr.c_str());
    split(myArrayTaxo, tmpChar, "_");
    sprintf(mdr->roughResult[2],"%s_%s",myArrayTaxo[1],myArrayTaxo[2]);
}
int roughGenus(Sequence *seqQuery, mdResult *mdr)
{
	int ret;
	char *tmpifile = "./tmp/query.fasta";
	char *tmpofile = "./tmp/result.fasta";
	char *hmmpath = "./bin";
	ofstream tmpQuery;
	tmpQuery.open(tmpifile);
	tmpQuery << ">query" << endl;

	VectorSequenceContainer *Item_sequences = new VectorSequenceContainer (&AlphabetTools::DNA_ALPHABET);
	Item_sequences->addSequence(*seqQuery);

	cout<<"query sequence:"<<Item_sequences->getSequence(0).getName()<<endl;
	tmpQuery << Item_sequences->getSequence(0).toString() << endl;
	tmpQuery.close();
	char hmmcommand[200];
    #ifdef __x86_64__
        sprintf(hmmcommand,"%s/hmmscan_x86_64 ./HmmDB/%s_%s %s > %s",hmmpath,mdr->dataBase,mdr->roughBarcode,tmpifile,tmpofile);
    #elif __i386__
        sprintf(hmmcommand,"%s/hmmscan_i386 ./HmmDB/%s_%s %s > %s",hmmpath,mdr->dataBase,mdr->roughBarcode,tmpifile,tmpofile);
    #endif
	//sprintf(hmmcommand,"%s/hmmscan ./HmmDB/%s %s > %s",hmmpath,mdr->dataBase,tmpifile,tmpofile);
	ret = system(hmmcommand);
	if (ret)
	{
		logFile<<"query sequence:"<<seqQuery->getName()<<" get Rough genus error!"<<endl;
		cout<<"query sequence:"<<seqQuery->getName()<<" get Rough genus error!"<<endl;
		exit(EXIT_FAILURE);
	}
	string result = readFileIntoString(tmpofile);
	string::size_type ps,pe;
    ps=result.find("No hits detected that satisfy reporting thresholds");
    if (ps!=std::string::npos)
    {
        return 1;
    }
	ps = result.find("-----------");
	pe = result.find("Domain annotation");
	string resultStr = result.substr(ps+20,pe-ps);
    //cout<<"resultStr:"<<resultStr<<endl;
    getRoughResult(resultStr, mdr);
    return 0;
}

//getMD
int getMD(Sequence *seqQuery, mdResult *mdr)
{
    char sqlStr[1000];
	char sequenceName[500];
    char **selectResult; //2d array
    int nrow = 0, ncolumn = 0;
    int i=0;

	//get Item selection dataset to identify to species level,get out every species
	VectorSequenceContainer *Item_sequences = new VectorSequenceContainer (&AlphabetTools::DNA_ALPHABET);
	sprintf(sqlStr, "select species_name from Item where family_name='%s' and genus_name='%s' and marker='%s' group by species_name",mdr->roughFamily, mdr->roughGenus,mdr->mdBarcode);
	logFile << "  +get species_name SQL:" << sqlStr << endl;
	sqlite3_get_table(conn, sqlStr, &selectResult, &nrow, &ncolumn, &dbErrMsg);
	if(nrow<1){
		return 2;//have not get data
	}
	for(i=1;i<nrow+1;i++){
	   	int index=ncolumn*i;
        char **itemResult; //2d array
        int nr = 0, nc = 0;
		//from individual species get 5 sequences if it has
		sprintf(mdr->speciesName,"%s", selectResult[index+0]);
        sprintf(sqlStr, "select sequenceID,family_name,genus_name,species_name,theta1_species,theta2_species,nucleotides from Item where family_name='%s' and genus_name='%s' and species_name='%s' and marker='%s' order by random() limit 20",mdr->roughFamily,mdr->roughGenus, mdr->speciesName,mdr->mdBarcode);
		logFile << "  +get sequences SQL:" << sqlStr << endl;
	    sqlite3_get_table(conn, sqlStr, &itemResult, &nr, &nc, &dbErrMsg);
	    for(int row=1;row<nr+1;row++){
    		int index=nc*row;
	    	sprintf(sequenceName,"%s_%s_%s_%s_%s_%s",itemResult[index+0],itemResult[index+1],itemResult[index+2],itemResult[index+3],itemResult[index+4],itemResult[index+5]);
            Sequence *seq = new BasicSequence(sequenceName, itemResult[index+6], &AlphabetTools::DNA_ALPHABET);
            Item_sequences->addSequence(*seq);
	    }
	    sqlite3_free_table(itemResult);
	}
	sqlite3_free_table(selectResult);

    //get MD
    double minDistance = 1.0;
    double tmpDistance = 1.0;
    int tmpPoint = 0;
    int num = Item_sequences->getNumberOfSequences();
    for(i=0;i<num;i++)
    {
        tmpDistance = pairwiseDistance(*seqQuery,Item_sequences->getSequence(i), mdr->model);
        if(minDistance >= tmpDistance)
        {
            minDistance = tmpDistance;
            tmpPoint = i;
        }
    }
    mdr->md = minDistance;
    //get MD sequence name
    char *myArray[8];
    sprintf(sequenceName,"%s",Item_sequences->getSequence(tmpPoint).getName().c_str());
    memset(myArray, 0x0, sizeof(myArray));
    split(myArray, sequenceName, "_");
    sprintf(mdr->taxoName,"%s_%s_%s",myArray[1],myArray[2],myArray[3]);

    sprintf(mdr->speciesName,"%s", myArray[3]);
    string theta1=myArray[4];
    string theta2=myArray[5];
    cout<<"theta1:"<<theta1<<endl;
    cout<<"theta2:"<<theta1<<endl;
    if(strcmp(myArray[4], "2.0")==0)//theta1==2, means parameter theta1 this species has not bean calculated
    {
        theta1=getTheta1AfterMD(mdr);
    }
    if(strcmp(myArray[5], "2.0")==0)//theta2==2, means parameter theta2 this species has not bean calculated
    {
        theta2=getTheta2AfterMD(mdr);
    }
    mdr->theta1 = atof(theta1.c_str());
    mdr->theta2 = atof(theta2.c_str());
	delete Item_sequences;
    logFile<<"    result:"<<mdr->taxoName<<endl;
}

//get theta1 sequences
void getTheta1function(theta1Result *t1r)
{
    //get genus
    char sqlStr[1000];
    char **genusResult; //2d array
    int nrow_genus = 0, ncolumn_genus = 0;

	sprintf(sqlStr, "select genus_name from Item where marker='%s' group by genus_name",t1r->mdBarcode);
	logFile << "  +get genus_name SQL:" << sqlStr << endl;
	sqlite3_get_table(conn, sqlStr, &genusResult, &nrow_genus, &ncolumn_genus, &dbErrMsg);
	for(int i=1;i<nrow_genus+1;i++){
		int index1=i*ncolumn_genus;
		char **speciesResult; //2d array
		int nrow_species = 0, ncolumn_species = 0;

		sprintf(t1r->genusName,"%s", genusResult[index1+0]);

		sprintf(sqlStr, "select species_name from Item where genus_name='%s' and marker='%s' group by species_name",t1r->genusName,t1r->mdBarcode);
		logFile << "    +get species_name SQL:" << sqlStr << endl;
		sqlite3_get_table(conn, sqlStr, &speciesResult, &nrow_species, &ncolumn_species, &dbErrMsg);
        if (nrow_species>11)
        {
            nrow_species=11;
        }
		for(int j=1;j<nrow_species+1;j++){
			int index2=j*ncolumn_species;
			char **sequenceResult; //2d array
			int nrow_sequence = 0, ncolumn_sequence = 0;

			sprintf(t1r->speciesName,"%s", speciesResult[index2+0]);

			//get the sequences existing in db
			VectorSequenceContainer *DB_sequences = new VectorSequenceContainer (&AlphabetTools::DNA_ALPHABET);
			sprintf(sqlStr, "select sequenceID,species_name,nucleotides from Item where genus_name='%s' and species_name='%s' and marker='%s' order by random() limit 20",t1r->genusName,t1r->speciesName,t1r->mdBarcode);
			logFile << "    +get sequences dataset SQL:" << sqlStr << endl;
			sqlite3_get_table(conn, sqlStr, &sequenceResult, &nrow_sequence, &ncolumn_sequence, &dbErrMsg);
			for(int k=1;k<nrow_sequence+1;k++){
				int index3=k*ncolumn_sequence;
				char name[100];
				sprintf(name,"%s_%s",sequenceResult[index3+0],sequenceResult[index3+1]);
				Sequence *seq = new BasicSequence(name, sequenceResult[index3+2], &AlphabetTools::DNA_ALPHABET);
				DB_sequences->addSequence(*seq);
				delete seq;
			}
			sqlite3_free_table(sequenceResult);
			if(DB_sequences->getNumberOfSequences()>0){
				getTheta1(DB_sequences, t1r);
				updateTheta1(t1r);
			}
			delete DB_sequences;
		}
		sqlite3_free_table(speciesResult);
	}
	sqlite3_free_table(genusResult);
	//update singleton theta1 with average value in the database
	updateSingletonTheta1(t1r);
}

//get theta1
void getTheta1(VectorSequenceContainer *sequences, theta1Result *t1r)
{
    double maxDistance1 = 0;
    double maxDistance2 = 0;
    double tmpDistance = 0;

    int num = sequences->getNumberOfSequences();
    //logFile << "    +fetch " << num << " samples for genus:"<<t1r->genusName<<" and species:" << t1r->speciesName << endl;
	if(num == 1){
		t1r->theta1 = 2;
		return;
	}
    for(int i = 0; i<num; i++)
    {
        for(int j = i+1; j<num; j++)
        {
            tmpDistance = pairwiseDistance(sequences->getSequence(i),sequences->getSequence(j),t1r->model);
            if(maxDistance1<tmpDistance){
                maxDistance2=maxDistance1;
                maxDistance1=tmpDistance;
            }else if(maxDistance2<tmpDistance){
                maxDistance2=tmpDistance;
            }else{
                continue;
            }
        }
    }
    cout<<"max D1:"<<maxDistance1<<endl;
    cout<<"max D2:"<<maxDistance2<<endl;
    if(maxDistance2!=0){
        t1r->theta1 = maxDistance2;
    }else{
        t1r->theta1 = maxDistance1;
    }
    cout<<"t1r-theta1:"<<t1r->theta1<<endl;
}

//update theta1
void updateTheta1(theta1Result *t1r)
{
	char sqlStr[1000];
	sprintf(sqlStr, "update Item set theta1_species=%.8f where genus_name='%s' and species_name='%s' and marker='%s'",t1r->theta1,t1r->genusName,t1r->speciesName,t1r->mdBarcode);
	logFile << "    +update theta1:" << sqlStr << endl;
	int nRes = sqlite3_exec(conn , sqlStr, 0, 0, &dbErrMsg);
	if (nRes != SQLITE_OK){
		cout<<"update Item theta1 fail: "<<dbErrMsg<<endl;
		exit(EXIT_FAILURE);
	}
}

//update singleton theta1 with average value in the database
void updateSingletonTheta1(theta1Result *t1r)
{
	char sqlStr[1000];
	char **familyResult; //2d array
    int nrow_family = 0, ncolumn_family = 0;
	double theta1;

    sprintf(sqlStr, "select family_name from Item where marker='%s' group by family_name",t1r->mdBarcode);
    sqlite3_get_table(conn, sqlStr, &familyResult, &nrow_family, &ncolumn_family, &dbErrMsg);
	for(int i=1;i<nrow_family+1;i++){
        int index1=i*ncolumn_family;

        char **speciesResult; //2d array
        int nrow_species = 0, ncolumn_species = 0;

        sprintf(sqlStr, "select count(*),sum(theta1_species) from Item where family_name='%s' and marker='%s' and theta1_species<>2.0",familyResult[index1+0],t1r->mdBarcode);
        sqlite3_get_table(conn, sqlStr, &speciesResult, &nrow_species, &ncolumn_species, &dbErrMsg);
	    if(atoi(speciesResult[2])==0){
		    continue;
	    }
        theta1=atof(speciesResult[3])/atof(speciesResult[2]);
        //cout<<"Family:"<<familyResult[index1+0]<<", sum(theta1_species):"<<speciesResult[3]<<", count(*):"<<speciesResult[2]<<", theta1:"<<theta1<<endl;
        sprintf(sqlStr, "update Item set theta1_species=%.8f where family_name='%s' and marker='%s' and theta1_species=2.0",theta1,familyResult[index1+0],t1r->mdBarcode);
        //logFile<<"  +update singleton theta1:" << sqlStr << endl;
        int nRes = sqlite3_exec(conn , sqlStr, 0, 0, &dbErrMsg);
        if (nRes != SQLITE_OK){
            cout<<"updateSingletonTheta1 fail: "<<dbErrMsg<<endl;
            exit(EXIT_FAILURE);
        }
        sqlite3_free_table(speciesResult);
    }
    sqlite3_free_table(familyResult);

	char **Result; //2d array
    int nrow = 0, ncolumn = 0;

	sprintf(sqlStr, "select count(*),sum(theta1_species) from Item where theta1_species<>2.0 and marker='%s'",t1r->mdBarcode);
	sqlite3_get_table(conn, sqlStr, &Result, &nrow, &ncolumn, &dbErrMsg);
    if(atoi(Result[2])==0){
        theta1=1.0;
    }else{
	    theta1=atof(Result[3])/atof(Result[2]);
    }
	//cout<<"sum(theta1_species):"<<Result[3]<<", count(*):"<<Result[2]<<", theta1:"<<theta1<<endl;
	sprintf(sqlStr, "update Item set theta1_species=%.8f where theta1_species=0.2 and marker='%s'",theta1,t1r->mdBarcode);
	//logFile<<"  +update singleton theta1:" << sqlStr << endl;
	int nRes = sqlite3_exec(conn , sqlStr, 0, 0, &dbErrMsg);
	if (nRes != SQLITE_OK){
		cout<<"updateSingletonTheta1 fail: "<<dbErrMsg<<endl;
		exit(EXIT_FAILURE);
	}
    sqlite3_free_table(Result);
}

//getTeta2Sequences
int getTheta2function(theta2Result *t2r)
{

    char sqlStr[1000];
    char **genusResult; //2d array
    int nrow_genus = 0, ncolumn_genus = 0;

    //get the genus name existing in db
    sprintf(sqlStr, "select genus_name from Item where marker='%s' group by genus_name",t2r->mdBarcode);
    //logFile << "  +get genus SQL:" << sqlStr << endl;
    sqlite3_get_table(conn, sqlStr, &genusResult, &nrow_genus, &ncolumn_genus, &dbErrMsg);
	for(int i=1;i<nrow_genus+1;i++){
		int index1=i*ncolumn_genus;
		char **speciesResult; //2d array
		int nrow_species = 0, ncolumn_species = 0;
		sprintf(t2r->genusName,"%s", genusResult[index1+0]);
		VectorSequenceContainer *DB_sequences = new VectorSequenceContainer (&AlphabetTools::DNA_ALPHABET);

		sprintf(sqlStr, "select species_name from Item where genus_name='%s' and marker='%s' group by species_name", t2r->genusName,t2r->mdBarcode);
		//logFile << "  +get species SQL:" << sqlStr << endl;
		sqlite3_get_table(conn, sqlStr, &speciesResult, &nrow_species, &ncolumn_species, &dbErrMsg);
        if (nrow_species>11)
        {
            sqlite3_free_table(speciesResult);
            delete DB_sequences;
            continue;
        }
		for(int j=1;j<nrow_species+1;j++){
			int index2=ncolumn_species*j;
			char **sequenceResult; //2d array
			int nrow_sequence = 0, ncolumn_sequence = 0;

			sprintf(t2r->speciesName,"%s", speciesResult[index2+0]);

			sprintf(sqlStr, "select nucleotides from Item where genus_name='%s' and species_name='%s' and marker='%s' order by random() limit 1", t2r->genusName,t2r->speciesName,t2r->mdBarcode);
			//logFile << "  +get theta2 sequences SQL:" << sqlStr << endl;
			sqlite3_get_table(conn, sqlStr, &sequenceResult, &nrow_sequence, &ncolumn_sequence, &dbErrMsg);
			for(int k=1;k<nrow_sequence+1;k++){
				int index3=ncolumn_sequence*k;
				Sequence *seq = new BasicSequence(t2r->speciesName, sequenceResult[index3+0], &AlphabetTools::DNA_ALPHABET);
				DB_sequences->addSequence(*seq);
				delete seq;
			}
			sqlite3_free_table(sequenceResult);
		}
		sqlite3_free_table(speciesResult);
		//the genus has only one species
		if(DB_sequences->getNumberOfSequences()<=1){
			sprintf(t2r->taxoName,"%s",DB_sequences->getName(0).c_str());
			t2r->theta2 = 2.0;
			updateTheta2(t2r);
		}
		else{
			//calculate theta2 and update db
			for(int n=0; n<DB_sequences->getNumberOfSequences(); n++)
			{
				int minTaxo = getMinTaxo(DB_sequences, n, t2r);
				VectorSequenceContainer *sequencesA = new VectorSequenceContainer (&AlphabetTools::DNA_ALPHABET);
				VectorSequenceContainer *sequencesB = new VectorSequenceContainer (&AlphabetTools::DNA_ALPHABET);
				getTaxoSequences(sequencesA, DB_sequences->getName(n).c_str(), t2r);
				getTaxoSequences(sequencesB, DB_sequences->getName(minTaxo).c_str(), t2r);
				getTheta2(sequencesA, sequencesB, t2r);
				updateTheta2(t2r);
				delete sequencesA;
				delete sequencesB;
			}
		}
		delete DB_sequences;
	}
	sqlite3_free_table(genusResult);

	//update singleton theta2 with average value in the database
	updateSingletonTheta2(t2r);
}

//get min distance taxo
int getMinTaxo(VectorSequenceContainer *sequences, int query, theta2Result *t2r)
{
    double minDistance = 2;
    double tmpDistance = 0;
    int point = 0;
    for(int i=0; i<sequences->getNumberOfSequences(); i++)
    {
        if(i != query)
        {
            tmpDistance = pairwiseDistance(sequences->getSequence(query),sequences->getSequence(i),t2r->model);
            if(minDistance > tmpDistance && tmpDistance >= 0.000001)
            {
                minDistance = tmpDistance;
                point = i;
            }
        }
    }
    sprintf(t2r->taxoName,"%s",sequences->getName(query).c_str());
    t2r->theta2 = minDistance;
    return point;
}

void getTaxoSequences(VectorSequenceContainer *sequences, const char *query, theta2Result *t2r)
{
    char sqlStr[1000];
    char **selectResult; //2d array
    int nrow = 0, ncolumn = 0;

    //get random 3 sequences for every species
    sprintf(sqlStr, "select sequenceID, nucleotides from Item where genus_name='%s' and species_name='%s' and marker='%s' order by random() limit 5", t2r->genusName,query,t2r->mdBarcode);
    //logFile << "  +get sequences SQL:" << sqlStr << endl;
    sqlite3_get_table(conn, sqlStr, &selectResult, &nrow, &ncolumn, &dbErrMsg);
	for(int i=1;i<nrow+1;i++){
		int index=i*ncolumn;
		char name[100];
		sprintf(name,"%s",selectResult[index+0]);
		Sequence *seq = new BasicSequence(name, selectResult[index+1], &AlphabetTools::DNA_ALPHABET);
		sequences->addSequence(*seq);
		delete seq;
	}
    sqlite3_free_table(selectResult);
}

void getTheta2(VectorSequenceContainer *sequencesA, VectorSequenceContainer *sequencesB, theta2Result *t2r)
{
    double tmpDistance = 0;
    double minDistance1 = 2;
    double minDistance2 = 2;

    for(int i=0; i<sequencesA->getNumberOfSequences(); i++)
    {
        for(int j=0; j<sequencesB->getNumberOfSequences(); j++)
        {
            tmpDistance = pairwiseDistance(sequencesA->getSequence(i),sequencesB->getSequence(j),t2r->model);
            if(minDistance1 > tmpDistance && tmpDistance >= 0.00001)
            {
                minDistance2 = minDistance1;
                minDistance1 = tmpDistance;
            }else if(minDistance2 > tmpDistance && tmpDistance >= 0.00001){
                minDistance2 = tmpDistance;
            }else{
                continue;
            }
        }
    }
    cout<<"min D1:"<<minDistance1<<endl;
    cout<<"min D2:"<<minDistance2<<endl;
    if(minDistance2!=2){
        t2r->theta2 = minDistance2;
    }else{
        t2r->theta2 = minDistance1;
    }
    cout<<"t2r-theta2:"<<t2r->theta2<<endl;
}

void updateTheta2(theta2Result *t2r)
{
	char sqlStr[1000];
//	if(t2r->theta2<0.05){
//		t2r->theta2 = 2.0;
//	}
	sprintf(sqlStr, "update Item set theta2_species=%.8f where genus_name='%s' and species_name='%s' and marker='%s'",t2r->theta2,t2r->genusName,t2r->taxoName,t2r->mdBarcode);
	logFile << "    +update theta2:" << sqlStr << endl;
	int nRes = sqlite3_exec(conn , sqlStr, 0, 0, &dbErrMsg);
	if (nRes != SQLITE_OK){
		cout<<"update Item theta2 fail: "<<dbErrMsg<<endl;
		exit(EXIT_FAILURE);
	}
}

//update singleton theta2 with average value in the database
void updateSingletonTheta2(theta2Result *t2r)
{
	char sqlStr[1000];
	char **familyResult; //2d array
    int nrow_family = 0, ncolumn_family = 0;
	double theta2;

    sprintf(sqlStr, "select family_name from Item where marker='%s' group by family_name",t2r->mdBarcode);
    sqlite3_get_table(conn, sqlStr, &familyResult, &nrow_family, &ncolumn_family, &dbErrMsg);
	for(int i=1;i<nrow_family+1;i++){
        int index1=i*ncolumn_family;

        char **speciesResult; //2d array
        int nrow_species = 0, ncolumn_species = 0;

	    sprintf(sqlStr, "select count(*),sum(theta2_species) from Item where family_name='%s' and marker='%s' and theta2_species<>2.0",familyResult[index1+0],t2r->mdBarcode);
        //logFile<<"sql:"<<sqlStr<<endl;
	    sqlite3_get_table(conn, sqlStr, &speciesResult, &nrow_species, &ncolumn_species, &dbErrMsg);
	    if(atoi(speciesResult[2])==0){
		    continue;
	    }
	    theta2=atof(speciesResult[3])/atof(speciesResult[2]);
	    //cout<<"Family:"<<familyResult[index1+0]<<", sum(theta2_species):"<<speciesResult[3]<<", count(*):"<<speciesResult[2]<<", theta2:"<<theta2<<endl;
	    sprintf(sqlStr, "update Item set theta2_species=%.8f where family_name='%s' and marker='%s' and theta2_species=2.0",theta2,familyResult[index1+0],t2r->mdBarcode);
        //logFile<<"update:"<<sqlStr<<endl;
	    //logFile<<"  +update singleton theta2:" << sqlStr << endl;
	    int nRes = sqlite3_exec(conn , sqlStr, 0, 0, &dbErrMsg);
	    if (nRes != SQLITE_OK){
		    cout<<"updateSingletonTheta2 fail: "<<dbErrMsg<<endl;
		    exit(EXIT_FAILURE);
        }
        sqlite3_free_table(speciesResult);
	}
    sqlite3_free_table(familyResult);

	char **Result; //2d array
    int nrow = 0, ncolumn = 0;

	sprintf(sqlStr, "select count(*),sum(theta2_species) from Item where marker='%s' theta2_species<>2.0",t2r->mdBarcode);
	sqlite3_get_table(conn, sqlStr, &Result, &nrow, &ncolumn, &dbErrMsg);
    if(atoi(Result[2])==0){
        theta2=1.0;
    }else{
        theta2=atof(Result[3])/atof(Result[2]);
    }
	//cout<<"sum(theta2_species):"<<Result[3]<<", count(*):"<<Result[2]<<", theta2:"<<theta2<<endl;
	sprintf(sqlStr, "update Item set theta2_species=%.8f where marker='%s' theta2_species=2.0",theta2,t2r->mdBarcode);
	//logFile<<"  +update singleton theta2:" << sqlStr << endl;
	int nRes = sqlite3_exec(conn , sqlStr, 0, 0, &dbErrMsg);
	if (nRes != SQLITE_OK){
		cout<<"updateSingletonTheta2 fail: "<<dbErrMsg<<endl;
		exit(EXIT_FAILURE);
	}
    sqlite3_free_table(Result);
}

string getTheta1AfterMD(mdResult *mdr)
{
    theta1Result *t1r;
    t1r = (theta1Result *)malloc(sizeof(theta1Result));
    char sqlStr[1000];
    char **sequenceResult; //2d array
    int nrow_sequence = 0, ncolumn_sequence = 0;
    VectorSequenceContainer *DB_sequences = new VectorSequenceContainer (&AlphabetTools::DNA_ALPHABET);

    sprintf(t1r->model,"%s",mdr->model);
    sprintf(t1r->mdBarcode,"%s",mdr->mdBarcode);
    sprintf(t1r->genusName,"%s",mdr->roughGenus);
    sprintf(t1r->speciesName,"%s",mdr->speciesName);

    sprintf(sqlStr, "select sequenceID,species_name,marker,nucleotides from Item where genus_name='%s' and species_name='%s' and marker='%s' order by random() limit 20",t1r->genusName,t1r->speciesName,t1r->mdBarcode);
    sqlite3_get_table(conn, sqlStr, &sequenceResult, &nrow_sequence, &ncolumn_sequence, &dbErrMsg);
    for(int k=1;k<nrow_sequence+1;k++){
        int index3=k*ncolumn_sequence;
        char name[100];
        sprintf(name,"%s_%s",sequenceResult[index3+0],sequenceResult[index3+1]);
        Sequence *seq = new BasicSequence(name, sequenceResult[index3+3], &AlphabetTools::DNA_ALPHABET);
        DB_sequences->addSequence(*seq);
        delete seq;
    }
    sqlite3_free_table(sequenceResult);
    if(DB_sequences->getNumberOfSequences()<2){
        delete DB_sequences;
        t1r->theta1 = 3.0;//This species has only one sequence, singleton
        updateTheta1(t1r);
        char theta1[20];
        sprintf(theta1,"%.8f",t1r->theta1);
        return theta1;
    }else{
        getTheta1(DB_sequences, t1r);
        updateTheta1(t1r);

        delete DB_sequences;
        char theta1[20];
        sprintf(theta1,"%.8f",t1r->theta1);
        return theta1;
    }
}

string getTheta2AfterMD(mdResult *mdr)
{
    theta2Result *t2r;
    t2r = (theta2Result *)malloc(sizeof(theta2Result));
    char sqlStr[1000];
    char **speciesResult; //2d array
    int nrow_species = 0, ncolumn_species = 0;
    VectorSequenceContainer *DB_sequences = new VectorSequenceContainer (&AlphabetTools::DNA_ALPHABET);
    Sequence *seqMD;

    sprintf(t2r->model,"%s",mdr->model);
    sprintf(t2r->mdBarcode,"%s",mdr->mdBarcode);
    sprintf(t2r->genusName,"%s",mdr->roughGenus);
    sprintf(t2r->taxoName,"%s",mdr->speciesName);

    sprintf(sqlStr, "select species_name from Item where genus_name='%s' and marker='%s' group by species_name", t2r->genusName,t2r->mdBarcode);
    //logFile << "  +get species SQL:" << sqlStr << endl;
    sqlite3_get_table(conn, sqlStr, &speciesResult, &nrow_species, &ncolumn_species, &dbErrMsg);
    for(int j=1;j<nrow_species+1;j++){
        int index=ncolumn_species*j;
        char **sequenceResult; //2d array
        int nrow_sequence = 0, ncolumn_sequence = 0;

        sprintf(sqlStr, "select nucleotides from Item where genus_name='%s' and species_name='%s' and marker='%s' order by random() limit 1", t2r->genusName,speciesResult[index+0],t2r->mdBarcode);
        //logFile << "  +get theta2 sequences SQL:" << sqlStr << endl;
        sqlite3_get_table(conn, sqlStr, &sequenceResult, &nrow_sequence, &ncolumn_sequence, &dbErrMsg);
        for(int k=1;k<nrow_sequence+1;k++){
            int index1=ncolumn_sequence*k;
            Sequence *seq = new BasicSequence(speciesResult[index+0], sequenceResult[index1+0], &AlphabetTools::DNA_ALPHABET);
            if(strcmp(speciesResult[index+0],t2r->taxoName)==0)
            {
                seqMD=seq->clone();
            }else{
                DB_sequences->addSequence(*seq);
            }
            delete seq;
        }
        sqlite3_free_table(sequenceResult);
    }
    sqlite3_free_table(speciesResult);

    double minDistance = 2;
    double tmpDistance = 0;
    if(DB_sequences->getNumberOfSequences()<1){
        delete DB_sequences;
        t2r->theta2 = 4.0;//This genus has only one species
        updateTheta2(t2r);
        char theta2[20];
        sprintf(theta2,"%.8f",t2r->theta2);
        return theta2;
    }
    else{
        int minPoint = 0;
        for(int n=0; n<DB_sequences->getNumberOfSequences(); n++){
            tmpDistance = pairwiseDistance(*seqMD,DB_sequences->getSequence(n),t2r->model);
            if(minDistance > tmpDistance && tmpDistance >= 0.000001)
            {
                minDistance = tmpDistance;
                minPoint = n;
            }
        }
        VectorSequenceContainer *sequencesA = new VectorSequenceContainer (&AlphabetTools::DNA_ALPHABET);
        VectorSequenceContainer *sequencesB = new VectorSequenceContainer (&AlphabetTools::DNA_ALPHABET);
        getTaxoSequences(sequencesA, seqMD->getName().c_str(), t2r);
        getTaxoSequences(sequencesB, DB_sequences->getName(minPoint).c_str(), t2r);
        getTheta2(sequencesA, sequencesB, t2r);
        updateTheta2(t2r);
        delete DB_sequences;
        delete sequencesA;
        delete sequencesB;

        char theta2[20];
        sprintf(theta2,"%.8f",t2r->theta2);
        return theta2;
    }
}
