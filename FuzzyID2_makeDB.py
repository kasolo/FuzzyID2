#!/usr/bin/env python
#encoding: utf-8
"""
approach:FuzzyID2_makeDB.py
version:1.0
author:zhiyong shi
email:zhiyong-shi@163.com
descripe:FuzzyID2_makeDB.py is the init approach of FuzzyID2 pipline. This approach converts bold tsv format to fasta format and save sequences in DB, and then form rough scan database.
date:2017.2.6

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
The following restriction to the GNU GPL applies: contact the author of this program (http://willpearse.github.com/phyloGenerator) for a suitable citation when publishing work that uses this code.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import time
import sys
import re
import sqlite3
import os
import platform

#codeCheck##########################################################################################
#check sequence DNA alphabet
def codeCheck(sequence):
    DNA="atgcumrwsykvhdbnATGCUMRWSYKVHDBN-"
    for i in range(0,len(sequence)):
        if sequence[i] not in DNA:
            return 1
    return 0

#seqInOne#############################################################################################
#convert multiple sequence line in one for fasta file
def seqInOne(refFile):
    if ".tsv" in refFile:
        fileName=refFile.split(".")[0]+".fasta"
    else:
        fileName=refFile
    #open files
    try:
        fi=open(fileName,"r")
        fo=open("./tmp/"+fileName.split("/")[-1].split('.')[0]+".one.fasta","w")
    except IOError as e:
        print("file open error:"+str(e))
        end = input("press any key to exit.")
        exit(1)

    flag=0
    for line in fi:
        if line.startswith('>'):
            fo.write("\n"+line)
            flag=1
        else:
            if flag==1:
                fo.write(line.strip())
    fo.write("\n")
    fi.close()
    fo.close()

    try:
        fi=open("./tmp/"+fileName.split("/")[-1].split('.')[0]+".one.fasta","r")
        fo=open("./tmp/"+fileName.split("/")[-1],"w")
    except IOError as e:
        print("file open error:"+str(e))
        end = input("press any key to exit.")
        exit(1)
    a=fi.readlines()
    b=''.join(a[1:])
    fo.write(b)
    fi.close()
    fo.close()

    command = "rm -f ./tmp/"+fileName.split("/")[-1].split('.')[0]+".one.fasta"
    os.system(command)

#tsv2fas############################################################################################
#get sequence from tsv file
def tsv2fas(refFile,paraList):
    """
    Reference dataset tsv file can be downloaded from BOLD system.
    Function tsv2fas extract target fields from tsv file and save them as fasta format. Fasta format: >SequenceID_Family_Genus_Species
    Out fasta file example:
    >GBCH3233-09_Araneidae_Verrucosa_Verrucosa arenata
    ACTTTATATTTGATTTTTGGAGCTTGGGCTGCTATAGTAGGAACTGCAATAAGAGT...
    """
    flog=open("log.txt","a")
    flog.write("#"*50+"\n")
    flog.write("#function tsv2fas start.\n")
    #open files
    if ".tsv" in refFile:
        fileName=refFile
    else:
        print("Reference dataset file format error, should be tsv format! Exiting...")
        sys.exit()
    try:
        fi=open(fileName,"r") # Opens file for reading
        fo=open("./tmp/"+fileName.split("/")[-1].split('.')[0]+".fasta","w")
        fe=open("./tmp/"+fileName.split("/")[-1].split('.')[0]+".info","w")
    except IOError as e:
        print("file open error:"+str(e))
        end = input("press any key to exit.")
        sys.exit()

    fe.write("****************************************************\n")
    fe.write("tsv2fas function error:")
    fe.write("DNA code check error: please check the DNA sequence!\n")
    fe.write("Data error: please check whether marker is not refMarker, DNA length is less than 500bp, Genus name is not given!\n")
    fe.write("****************************************************\n\n")

    #read file per line
    for line in fi:
        taxo={}
        item=line.split('\t')
        taxo['SequenceID']=item[0].strip()
        taxo['Family']=item[14].strip()
        taxo['Genus']=item[18].strip()
        taxo['Species']=item[20].strip()
        taxo['Marker']=item[40].strip()
        taxo['DNA']=item[42].strip().replace("-","").replace("N","")
        if len(taxo['DNA'])>int(paraList["Max_length_DNA"]):
            taxo['DNA']=taxo['DNA'][0:int(paraList["Max_length_DNA"])]

        if codeCheck(taxo['DNA'])==1:
            fe.write("DNA alphabet check error! processid:"+taxo['SequenceID']+"\n")
            continue

        if re.match('^[0-9A-Z-]+$',taxo['SequenceID']):#code check
            if taxo['Marker'].upper() in paraList["refMarker"].upper() and len(taxo['DNA'])>int(paraList["Min_length_DNA"]) and len(taxo['Family'])>0 and len(taxo['Genus'])>0:#sequence quality filter
                #write data to file
                fo.write(">"+taxo['SequenceID']+"_"+taxo['Family']+"_"+taxo['Genus']+"_"+taxo['Species']+"\n")
                fo.write(taxo['DNA']+"\n")
            else:
                fe.write("Data error! processid:"+taxo['SequenceID']+"\n")
    #close file
    fi.close()
    fo.close()
    fe.close()
    flog.close()

#deleteRepeatSeq#############################################################################################
#deal with Haplotype sequences
def deleteRepeatSeq(refFile):
    """
    Function deleteRepeatSeq delete haplotype sequences and reserve two.
    """
    flog=open("log.txt","a")
    flog.write("#"*50+"\n")
    flog.write("#function deleteRepeatSeq start.\n")
    #open files
    if ".tsv" in refFile:
        fileName=refFile.split("/")[-1].split(".")[0]+".fasta"
    else:
        fileName=refFile.split("/")[-1]
    try:
        fi=open("./tmp/"+fileName,"r")
        fo=open("./tmp/"+fileName.split('.')[0]+".out.fasta","w")
        fs=open("./tmp/"+fileName.split('.')[0]+".sum.txt","w")
    except IOError as e:
        print("file open error:"+str(e))
        end = input("press any key to exit.")
        sys.exit()

    dictDNA={}#Dictionary with sequence as key and comment as value
    flag=0
    key=""
    value=""
    for line in fi:
        if len(line.strip())==0:
            continue
        if line.startswith('>'):
            value=line.strip()
            flag=0
        else:
            key=line.strip()
            flag=1
        if flag==1:
            if key in dictDNA:#add new Haplotype into Dict 
                tmp=dictDNA[key]+"|"+value
                dictDNA[key]=tmp
            else:
                dictDNA[key]=value
    for key in dictDNA:
        fo.write(dictDNA[key].split('|')[0]+"\n")
        fo.write(key+"\n")
        fs.write(dictDNA[key]+"\n")

    fi.close()
    fo.close()
    fs.close()

    command = "mv -f ./tmp/"+fileName.split('.')[0]+".out.fasta ./tmp/"+fileName
    os.system(command)
    flog.close()

#fas2db#############################################################################################
#save data in database
def fas2db(refFile,paraList):
    """
    Function fas2db extract fields from fasta file and save them in database. Fasta format: >SequenceID_Family_Genus_Species
    fasta file example:
    >GBCH3233-09_Araneidae_Verrucosa_Verrucosa arenata
    ACTTTATATTTGATTTTTGGAGCTTGGGCTGCTATAGTAGGAACTGCAATAAGAGT...
    """
    flog=open("log.txt","a")
    flog.write("#"*50+"\n")
    flog.write("#function fas2db start.\n")
    #open file and database
    if ".tsv" in refFile:
        fileName=refFile.split("/")[-1].split(".")[0]+".fasta"
    else:
        fileName=refFile.split("/")[-1]

    if os.path.exists("./tmp/"+fileName.split('.')[0]+".error"):
        fe=open("./tmp/"+fileName.split('.')[0]+".error","a")
    else:
        fe=open("./tmp/"+fileName.split('.')[0]+".error","w")
    try:
        fi=open("./tmp/"+fileName,"r") # Opens file for reading
        conn = sqlite3.connect("./DB/"+paraList["refDBName"])#connect database
        cur = conn.cursor()#create cursor
    except IOError as e:
        print("file open error:"+str(e))
        end = input("press any key to exit.")
        sys.exit()
        
    #conn.isolation_level = None
    fe.write("****************************************************\n")
    fe.write("fas2db function error:")
    fe.write("DNA code check error: please check the DNA sequence!\n")
    fe.write("****************************************************\n\n")

    ###########################################################
    #create Item table
    command=("CREATE TABLE if not exists Item("
        "recordID integer primary key autoincrement,"
        "sequenceID varchar(50) NOT NULL,"
        "nucleotides varchar(1200) NOT NULL,"
        "marker varchar(50) NOT NULL,"
        "theta1_species double(10,8) DEFAULT 2,"
        "theta2_species double(10,8) DEFAULT 2,"
        "family_name varchar(50) NOT NULL,"
        "genus_name varchar(50) NOT NULL,"
        "species_name varchar(50) NOT NULL)")
    cur.execute(command)
    #create Rough table
    command=("CREATE TABLE if not exists Rough("
        "recordID integer primary key autoincrement,"
        "sequenceID varchar(50) NOT NULL,"
        "nucleotides varchar(1200) NOT NULL,"
        "marker varchar(50) NOT NULL,"
        "family_name varchar(50) NOT NULL,"
        "genus_name varchar(50) NOT NULL,"
        "species_name varchar(50) NOT NULL)")
    cur.execute(command)
    conn.commit()
    ###############################################################
    #read fasta file and save data in database

    refMarker=""
    for marker in paraList["refMarker"].split("|"):
        if marker in refFile:
            refMarker=marker.strip().upper()
            break

    cmdItem=""
    taxo={}
    flag=0
    sql_lines=[]
    for line in fi:
        if line.startswith('>'):
            item=line.split('>')[1].split('_')
            taxo['SequenceID']=item[0].strip()
            taxo['Family']=item[1].strip()
            taxo['Genus']=item[2].strip()
            taxo['Species']=item[3].strip()

            if len(taxo['Genus'])==0 or len(taxo['Species'])==0:
                flag=0
                continue
            #insert data in Item table
            cmdItem=("insert into Item("
                "SequenceID,nucleotides,marker,family_name,genus_name,species_name)"
                "values('"+taxo['SequenceID']+"','nucleotides_replace','"+refMarker+"','"+taxo['Family']+"','"+taxo['Genus']+"','"+taxo['Species']+"');")
            flag=1
        else:
            if flag==1:
                #add DNA check
                sequence=line.strip().replace("-","").replace("N","")
                if len(sequence)>int(paraList["Max_length_DNA"]):
                    sequence=sequence[0:int(paraList["Max_length_DNA"])]
                if len(sequence)<int(paraList["Min_length_DNA"]):
                    cmdItem=""
                    sequence=""
                    flag=0
                    continue

                if codeCheck(sequence)==1:
                    fe.write("DNA code check error! SequenceID:"+taxo['SequenceID']+"\n")
                    cmdItem=""
                    sequence=""
                    flag=0
                    continue
                cmdItem=cmdItem.replace("nucleotides_replace",sequence)
                sql_lines.append(cmdItem)
                cmdItem=""
                sequence=""
                flag=0
    flog.write('\r\n'.join(sql_lines))
    for sql in sql_lines:
        cur.execute(sql)
    conn.commit()
    flog.write("Create index on Item by genus_name and species_name.")
    sql="CREATE INDEX Item_genus on Item(genus_name)"
    cur.execute(sql)
    sql="CREATE INDEX Item_species on Item(species_name)"
    cur.execute(sql)
    conn.commit()

    #close file and db
    fi.close()
    cur.close()
    conn.close()
    flog.close()

#makeHmm###########################################################################################
#Fetch data from Rough table and convert to hmm data format
def makeHmm(paraList):
    """
    Function makeHmm fetch data from Rough table and convert to hmm data format.
    """
    flog=open("log.txt","a")
    flog.write("#"*50+"\n")
    flog.write("#function makeHmm start.\n")
    #connect database
    try:
        conn = sqlite3.connect("./DB/"+paraList["refDBName"])
        cur = conn.cursor()
    except IOError as e:
        print("database open error:"+str(e))
        end = input("press any key to exit.")
        sys.exit()
    #conn.isolation_level = None

    #################################################################
    #Fetch data randomly from Item table and insert into Rough table#
    #cur.execute("BEGIN TRANSACTION")
    flog.write("Fetch data randomly from Item table and insert into Rough table.\n")
    sql_lines=[]
    genusSql="select genus_name from Item where marker='"+paraList["refRoughMarker"].upper()+"' group by genus_name"
    cur.execute(genusSql)
    resGenus = cur.fetchall()
    for genus_name in resGenus:
        speciesSql="select species_name from Item where genus_name='"+genus_name[0]+"' and marker='"+paraList["refRoughMarker"].upper()+"' group by species_name"
        cur.execute(speciesSql)
        resSpecies = cur.fetchall()
        for species_name in resSpecies:
            countSql="select count(*) from Item where species_name='"+species_name[0]+"' and marker='"+paraList["refRoughMarker"].upper()+"'"
            cur.execute(countSql)
            resSeq=cur.fetchone()
            seqNum=int(resSeq[0])
            if seqNum<int(paraList["roughMinNum"]):
                roughNum=seqNum
            elif seqNum*float(paraList["roughPercent"])<int(paraList["roughMinNum"]):
                roughNum=int(paraList["roughMinNum"])
            elif seqNum*float(paraList["roughPercent"])>int(paraList["roughMaxNum"]):
                roughNum=int(paraList["roughMaxNum"])
            else:
                roughNum=int(seqNum*float(paraList["roughPercent"]))
            insertRough=("insert into Rough"
                     "(sequenceID,nucleotides,marker,family_name,genus_name,species_name) "
                     "select sequenceID,nucleotides,marker,family_name,genus_name,species_name "
                     "from Item where marker='"+paraList["refRoughMarker"].upper()+"' and genus_name='"+str(genus_name[0])+"' and species_name='"+str(species_name[0])+"' order by RANDOM() limit "+str(roughNum))
            #print(insertRough)
            sql_lines.append(insertRough)
    flog.write('\r\n'.join(sql_lines))
    for sql in sql_lines:
        cur.execute(sql)
    conn.commit()
    #cur.execute("COMMIT")

    #################################################################
    #Make HmmDB using rough table data#
    flog.write("Make HmmDB using rough table data.\n")
    sqlstr = "select sequenceID,family_name,genus_name,nucleotides from Rough where marker='"+paraList["refRoughMarker"].upper()+"' order by sequenceID"
    head = "# STOCKHOLM 1.0\r\n"

    hmmbuild="hmmbuild"
    hmmpress="hmmpress"
    if platform.uname().processor=="x86_64":
        hmmbuild=hmmbuild+"_x86_64"
        hmmpress=hmmpress+"_x86_64"
    elif platform.uname().processor=="i386":
        hmmbuild=hmmbuild+"_i386"
        hmmpress=hmmpress+"_i386"

    n = cur.execute(sqlstr)
    for row in cur.fetchall():
        fo=open("./HmmDB/"+str(row[0])+"_"+str(row[1])+"_"+str(row[2])+".fasta","w")
        fo.write(head)
        fo.write(str(row[0])+"_"+str(row[1])+"_"+str(row[2])+" "+str(row[3])+"\r\n"+"//\r\n")
        fo.close()
        hmmpath = "./bin/"
        ifile = "./HmmDB/" + str(row[0])+"_"+str(row[1])+"_"+str(row[2]) + ".fasta"
        ofile = "./HmmDB/" + str(row[0])+"_"+str(row[1])+"_"+str(row[2]) + ".hmm"
        hmmcommand = hmmpath + hmmbuild + " " + ofile + " " + ifile
        os.system(hmmcommand)

    #hmmcommand = "cat " + "./HmmDB/*.hmm > " + "../HmmDB/"+paraList["refDBName"]
    hmmcommand = "find ./HmmDB -type f -name \"*.hmm\"|xargs cat > ./HmmDB/"+paraList["refDBName"]+"_"+paraList["refRoughMarker"].upper()
    os.system(hmmcommand)
    hmmcommand = hmmpath + hmmpress + " " + "./HmmDB/"+paraList["refDBName"]+"_"+paraList["refRoughMarker"].upper()
    os.system(hmmcommand)
    #hmmcommand="rm -fr ./HmmDB/*.fasta"
    hmmcommand = "find ./HmmDB -type f -name \"*.fasta\"|xargs rm -f"
    os.system(hmmcommand)
    #hmmcommand="rm -fr ./HmmDB/*.hmm"
    hmmcommand = "find ./HmmDB -type f -name \"*.hmm\"|xargs rm -f"
    os.system(hmmcommand)

    #close db
    cur.close()
    conn.close()
    flog.close()

################################################################
def getPara(paraFileName):
    """
    get parameters from config file
    """
    flog=open("log.txt","a")
    flog.write("#"*50+"\n")
    flog.write("#function getPara start.\n")
    try:
        paraFile=open(paraFileName,"r")
    except IOError as e:
        print("\nERROR: config file [config.txt] not found. Exiting...")
        sys.exit()
    para={}
    for line in paraFile:
        if  len(line.strip())==0:
            continue
        if line.strip().startswith('#'):
            continue
        pKey=line.split('=')[0].strip()
        pValue=line.split('=')[1].split('#')[0].strip()
        para[pKey]=pValue
        flog.write("parameter "+pKey+":"+pValue+"\n")
    paraFile.close()
    flog.write("#"*30+"\n")
    flog.close()

    ############################
    #parameter check
    for refFile in para["refFileName"].split("|"):
        if os.path.exists(refFile):
            pass
        else:
            print("Input reference dataset file:"+refFile+" dose not exist! Exiting...")
            sys.exit()
        if ".tsv" in refFile or ".fas" in refFile or ".fasta" in refFile:
            pass
        else:
            print("Please check file Format of reference dataset file:"+refFile+", should be tsv or fasta! Exiting...")
            sys.exit()

    for refMarker in para["refMarker"].split("|"):
        if refMarker.upper() in para["refFileName"].upper():
            pass
        else:
            print("refMarker:"+refMarker+" should be contained in refFileName! Exiting...")
            sys.exit()

    for refMDMarker in para["refMDMarker"].split("|"):
        if refMDMarker.upper() in para["refFileName"].upper():
            pass
        else:
            print("refMDMarker:"+refMDMarker+" should be contained in refFileName! Exiting...")
            sys.exit()

    if int(para["roughMaxNum"])>0 and int(para["roughMaxNum"])<=20:
        pass
    else:
        print("Please specify parameter roughMaxNum, int 1-20! Exiting...")
        sys.exit()

    if int(para["roughMinNum"])>0 and int(para["roughMinNum"])<=5:
        pass
    else:
        print("Please specify parameter roughMixNum, int 1-20! Exiting...")
        sys.exit()

    if float(para["roughPercent"])>0.1 and float(para["roughPercent"])<0.51:
        pass
    else:
        print("Please specify parameter roughPercent, float 0.1-0.5! Exiting...")
        sys.exit()

    if para["deleteRepeatSeq"]=="on" or para["deleteRepeatSeq"]=="off":
        pass
    else:
        print("Please swich deleteRepeatSeq, on or off! Exiting...")
        sys.exit()

    if para["estimateFuzzyPara"]=="on" or para["estimateFuzzyPara"]=="off":
        pass
    else:
        print("Please swich estimateFuzzyPara, on or off! Exiting...")
        sys.exit()

    return para

#main#############################################################################################
if __name__ == '__main__':
    #main
    flog=open("log.txt","w")
    flog.write("FuzzyID2_makeDB.py: make reference database of FuzzyID2.\n")
    flog.write("Any problems, please email me.[zhiyong-shi@163.com]\n")
    flog.close()
    print("FuzzyID2_makeDB.py: make reference database of FuzzyID2.")
    print("Any problems, please email me.[zhiyong-shi@163.com]")
    print("*"*60)

    #process start time
    timeb=time.strftime('%H%M%S', time.localtime(time.time()))
    print("Process start at "+timeb)
    print("-"*60)

    #input
    paraFile=input("Please input config file name:")
    paraList=getPara(paraFile)#parameters
    print("parameters:")
    print(paraList)
    print("-"*60)
    flog=open("log.txt","a")
    flog.write("The file format of reference dataset is "+paraList["refFileFormat"]+"\n")
    flog.close()

    for refFile in paraList["refFileName"].split("|"):
        if ".tsv" in refFile:
            print("Processing bold_tsv2fas start.........")
            tsv2fas(refFile,paraList)
        seqInOne(refFile)
        print("Processing fas2db start.........")
        if paraList["deleteRepeatSeq"]=="on":
            deleteRepeatSeq(refFile)
        fas2db(refFile,paraList)
        if paraList["refRoughMarker"].upper() in refFile.upper():
            print("Processing makeHmm start........")
            makeHmm(paraList)

    if paraList["estimateFuzzyPara"]=="on":
        for MDMarker in paraList["refMDMarker"].split("|"):
            command="./FuzzyID2 -c Theta1 -m K2P -d "+paraList["refDBName"]+" -mb "+MDMarker.upper()
            os.system(command)
            command="./FuzzyID2 -c Theta2 -m K2P -d "+paraList["refDBName"]+" -mb "+MDMarker.upper()
            os.system(command)

    #process end time
    print("-"*60)
    timee=time.strftime('%H%M%S', time.localtime(time.time()))
    print("Process end at "+timee)
    end = input("press any key to exit.")
