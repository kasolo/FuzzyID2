# FuzzyID2
FuzzyID2 is an approach for species identification using DNA barcoding tech.

FuzzyID2 installation:

    X86_64 compile command:
        g++ main.cpp mainFunction.cpp commonFunction.cpp -o FuzzyID2_X86_64 -lsqlite3 -lbpp-core -lbpp-seq -lbpp-phyl -Wl,-rpath ./lib64 -L ./lib64 -I ./include
        ln -s FuzzyID2_X86_64 FuzzyID2

    i686 compile command:
        g++ main.cpp mainFunction.cpp commonFunction.cpp -o FuzzyID2_i686 -lsqlite3 -lbpp-core -lbpp-seq -lbpp-phyl -Wl,-rpath ./lib -L ./lib -I ./include
        ln -s FuzzyID2_i686 FuzzyID2

FuzzyID2 execute commands:

    1.make reference database
        python3 FuzzyID2_makeDB.py
    2.fuzzy formula paramaters estimation
        ./FuzzyID2 -c Theta1 -m K2P -d AAL -mb COI
        ./FuzzyID2 -c Theta2 -m K2P -d AAL -mb COI
    3.query sequence identification
        ./FuzzyID2 -c ID -in ./test/AAL_COI_query.fas -m K2P -out ./test/AAL-out.txt -d AAL -rb COI -mb COI

For more information, please read the FuzzyID2-guide.pdf
