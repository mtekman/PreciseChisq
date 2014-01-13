#include <boost/tokenizer.hpp>
#include <fstream>
#include <string>

#include "ChiSquare.h"

using namespace boost;

//inline expands the code inside the scope of where it is called. Almost a macro function
inline double convertToDouble(std::string const& s){
  std::istringstream i(s); double x;
  if (!(i >> x)){
    cerr << "Could not convert to double: " << s << endl;
    exit(-1);
  }
  return x;
}

inline int convertToInt(std::string const& s){
  std::istringstream i(s); int x;
  if (!(i >> x)){
    cerr << "Could not convert to int: " << s << endl;
    exit(-1);
  }
  return x;
}

void checkArg(char * arg, string valid_arg_start){
    string setarg = string(arg).substr(0,valid_arg_start.size());
    if (setarg!=valid_arg_start){
        cerr << "Please give " << valid_arg_start << "N argument!" << endl;
        exit(-1);
    }
}


int main(int argc, char ** argv)
{
    //Args
    string setnum_id("--#sets=");
    string valset_id("--#valsper=");
    string sigfig_id("--#sf=");

    if (argc<4) {
        cerr << endl;
        cerr << "Takes in a text file where each row contains several different sets and calculates the chi-squared and p-value across all sets for that row.";
        cerr << endl;
        cerr << "\nEach row of the text file is delimited by spaces in the format:";
        cerr << "\n\trow1:\tA11 A12 B11 B12 C11 C12 D11 D12 ... etc";
        cerr << "\n\trow2:\tA21 A22 B21 B22 C21 C22 D21 D22 ... etc" << endl;
        cerr << "\tetc..\n";
        cerr << endl;
        cerr << "Set number MUST be specified to indicate how many values there are per set (in the above example there are 2).\n" << endl;
        cerr << "  usage:  " << argv[0] << " <text.input> "<< setnum_id << "N " << valset_id << "N [" << sigfig_id << "N]" <<  endl;
        cerr << endl;
        exit(-1);
    }

    //File arg
    ifstream file(argv[1]);

    //////////Flags
    //Reqs
    checkArg(argv[2],setnum_id);
    checkArg(argv[3],valset_id);

    int setnum = convertToInt(string(argv[2]).substr(setnum_id.size()));
    int valset = convertToInt(string(argv[3]).substr(valset_id.size()));

    //Opts
    int sigfig=3; //default
    if (argc==5){
        checkArg(argv[4],sigfig_id);
        sigfig = convertToInt(string(argv[4]).substr(sigfig_id.size()));
    }

    /////////Begin
    chi_squared_distribution<long double> dist(setnum-1);  // Make chisq once

    //Read in file, process lines
    typedef tokenizer<char_separator<char> > tokenizer;

    string line;
    while(getline(file, line))
    {
        vector<vector<int> > row_data;
        tokenizer tok(line);

        for (tokenizer::iterator beg=tok.begin(); beg!=tok.end(); false) {
            vector<int> val_in_set;
            for (int r=0; r < valset; r++){
                int val = convertToInt(*beg);

                val_in_set.push_back(val);

                if ((++beg).at_end())break;
            }
            row_data.push_back( val_in_set );
        }
        ChiSquare(row_data, dist, valset, sigfig);
    }
    return 0;
}

