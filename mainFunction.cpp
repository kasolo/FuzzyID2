#include "main.h"

void ID(InData *optIn)
{
    Fasta fasReader;
    VectorSequenceContainer *roughQuerySeqs;
    VectorSequenceContainer *MDQuerySeqs;

    //read the query sequences
    char *myArray[8];
    char inFileName[100];
    strcpy(inFileName , optIn->inFile);
    memset(myArray, 0x0, sizeof(myArray));
    cout<<"in file:"<<inFileName<<endl;
    split(myArray, inFileName, "+");
    for(int i=0;i<6;i++)
    {
        if(myArray[i]!='\0'){
            string inFile=myArray[i];
            if(findIgnoreCase(myArray[i],optIn->roughBarcode)!=std::string::npos){
                try {
                    roughQuerySeqs = fasReader.readSequences(inFile, &AlphabetTools::DNA_ALPHABET);
                }catch(exception& e){
                    cout<<"Please check sequence DNA alphabet! Exit..."<<endl;
                    cout<<e.what()<<endl;
                    exit(EXIT_FAILURE);
                }
                if(roughQuerySeqs->getNumberOfSequences()==0){
                    cout<<"Rough marker query seqFile:"<< inFile <<" has 0 sequence! Exit..."<<endl;
                    exit(EXIT_FAILURE);
                }
                //cout<<"rough seq ok!"<<endl;
            }
            if(findIgnoreCase(myArray[i],optIn->mdBarcode)!=std::string::npos){
                try {
                    MDQuerySeqs = fasReader.readSequences(inFile, &AlphabetTools::DNA_ALPHABET);
                }catch(exception& e){
                    cout<<"Please check sequence DNA alphabet! Exit..."<<endl;
                    cout<<e.what()<<endl;
                    exit(EXIT_FAILURE);
                }
                if(MDQuerySeqs->getNumberOfSequences()==0){
                    cout<<"MD marker query seqFile:"<< inFile <<" has 0 sequence! Exit..."<<endl;
                    exit(EXIT_FAILURE);
                }
                //cout<<"md seq ok!"<<endl;
            }
        }
    }

    int querySequenceNum = roughQuerySeqs->getNumberOfSequences();
    logFile << "The rough marker query file has " << querySequenceNum << " sequences!" << endl;

    for(int i=0;i<querySequenceNum;i++)
    {
        int flag=0;
        mdResult *mdr;
        mdr = (mdResult *)malloc(sizeof(mdResult));

        outFile << roughQuerySeqs->getSequence(i).getName();
        logFile << roughQuerySeqs->getSequence(i).getName() << endl;
        if(roughQuerySeqs->getSequence(i).toString().length()==0)
        {
            outFile << "\t No sequence!" << endl;
            continue;
        }
        Sequence *seqRoughQuery = new BasicSequence(roughQuerySeqs->getName(i),roughQuerySeqs->getSequence(i).toString(),&AlphabetTools::DNA_ALPHABET);

        //calculate order,family,genus,species MD and get the theta1,theta2

        sprintf(mdr->model,"%s",optIn->model);
        sprintf(mdr->dataBase,"%s",optIn->dataBase);
        sprintf(mdr->roughBarcode,"%s",optIn->roughBarcode);
        sprintf(mdr->mdBarcode,"%s",optIn->mdBarcode);

        flag=roughGenus(seqRoughQuery, mdr);
        if (flag==1)
        {
            outFile << "\t The query sequence is too far from the reference dataset!" << endl;
            continue;
        }

        string roughResultList="";
        for(int j=0;j<3;j++){
            //cout<<"roughResultList:"<<roughResultList<<endl;
            //cout<<"roughResult:"<<mdr->roughResult[j]<<endl;
            if(roughResultList.find(mdr->roughResult[j])!=std::string::npos){
                continue;
            }else{
                string roughResultStr=mdr->roughResult[j];
                roughResultList=roughResultList+roughResultStr;
            }

            //get rough genusName
            char *myArrayTaxo[3];
            char tmpChar[100];
            memset(myArrayTaxo, 0x0, sizeof(myArrayTaxo));
            sprintf(tmpChar,"%s",mdr->roughResult[j]);
            split(myArrayTaxo, tmpChar, "_");
            sprintf(mdr->roughFamily,"%s",myArrayTaxo[0]);
            sprintf(mdr->roughGenus,"%s",myArrayTaxo[1]);
            //cout<<"roughFamily"<<mdr->roughFamily<<endl;
            //cout<<"roughGenus"<<mdr->roughGenus<<endl;
            if(strcmp(mdr->roughFamily, "(null)")==0 or strcmp(mdr->roughGenus, "(null)")==0){
                continue;
            }
            logFile<<"    Rough Family:"<<mdr->roughFamily<<endl;
            logFile<<"    Rough genus:"<<mdr->roughGenus<<endl;

            outFile << "\t" <<"Rough Family:" << mdr->roughFamily;
            outFile << "\t" <<"Rough Genus:" << mdr->roughGenus;

            Sequence *seqMDQuery = new BasicSequence(roughQuerySeqs->getName(i),MDQuerySeqs->getSequence(roughQuerySeqs->getName(i)).toString(),&AlphabetTools::DNA_ALPHABET);
            //get MD,theta1,theta2
            //logFile<<"get MD start!"<<endl;
            flag = getMD(seqMDQuery, mdr);
            if(flag == 2)//doesn't found the MD sequence
            {
                outFile << "\t" << mdr->roughGenus << " No appropriate species data in the reference dataset!" << endl;
                continue;
            }

            //get FMF score
            outFile << "\t" << mdr->taxoName;
            outFile << "\t" <<"FMF:" << calMF(mdr->md,mdr->theta1,mdr->theta2);
            outFile << "\t" <<"MD:" << mdr->md;
            outFile << "\t" <<"Theta1:" << mdr->theta1;
            outFile << "\t" <<"Theta2:" << mdr->theta2 << endl;
        }
    }//for,loop query sequences
}

//Theta1
void Theta1(InData *optIn)
{
    theta1Result *t1r;
    t1r = (theta1Result *)malloc(sizeof(theta1Result));

    logFile << "getTheta1function start " << endl;
    //get the sequences existing in db
	sprintf(t1r->model,"%s",optIn->model);
	sprintf(t1r->mdBarcode,"%s",optIn->mdBarcode);
    getTheta1function(t1r);
}

//Theta2
void Theta2(InData *optIn)
{
    theta2Result *t2r;
    t2r = (theta2Result *)malloc(sizeof(theta2Result));
    logFile << "getTheta2function start " << endl;

    //get the sequences existing in db
	sprintf(t2r->model,"%s",optIn->model);
	sprintf(t2r->mdBarcode,"%s",optIn->mdBarcode);
    int flag = getTheta2function(t2r);
}
