#include <cstring>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;


string get_line_fasta(ifstream* partToExamine)
{
        string line,justTheSequence;
        getline(*partToExamine,line);
        char caracter=partToExamine->peek();
        while(caracter!='>' and caracter!=EOF) //avoid line with '>' and End Of File
        {
                getline(*partToExamine,line);
                justTheSequence+=line;
                caracter=partToExamine->peek();
        }
        return justTheSequence;
}


vector<string> parse_this_fasta(string fastaFile)
{
        vector<string> genomeVector;
        ifstream theRead(fastaFile);
        while(not theRead.eof()) //put sequences string in genome vector while it's not End of File
        {
                genomeVector.pushback(get_line_fasta(&theRead));
        }
        return genomeVector;
}
