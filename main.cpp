#include <boost/tokenizer.hpp>
#include <fstream>
#include <string>

#include "ChiSquareALL.h"

using namespace boost;

//inline expands the code inside the scope of where it is called. Almost a macro function
inline Precision convertToPrecision(std::string const& s){
  std::istringstream i(s); Precision x;
  if (!(i >> x)){
    cerr << "Could not convert to precision format: " << s << endl;
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

inline Pos convertToPosition(std::string const& s){
  std::istringstream i(s);  Pos x;
  if (!(i >> x)){
    cerr << "Could not convert to position: " << s << endl;
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
    string setnum_id("--sets=");
    string sigfig_id("--sf=");

    if (argc<4) {
        cerr << endl;
        cerr << "Takes in a text file where each row contains several different sets and calculates the chi-squared and p-value across all sets for that row.";
        cerr << endl;
        cerr << "The first column contains position, and will not be processed." << endl;
        cerr << "\nEach row of the text file is delimited by spaces in the format:";
        cerr << "\n\tpos1\tA11 A12 B11 B12 ...";
        cerr << "\n\tpos2\tA21 A22 B21 B22 ..." << endl;
        cerr << "\n\tpos3\tA31 A32 B31 B32 ..." << endl;
        cerr << "\tetc..\n";
        cerr << endl;
        cerr << "In the above example there are 2 sets (A and B) with 4 values per set\n" << endl;
        cerr << "  usage:  " << argv[0] << " <text.input> "<< setnum_id << "N " << sigfig_id << "N" <<  endl;
        cerr << "\ndefault output is to 5 sigfig" << endl;
        cerr << endl;
        exit(-1);
    }

    //File arg
    ifstream file(argv[1]); //

    // Flags
    //Reqs
    checkArg(argv[2],setnum_id);
    checkArg(argv[3],sigfig_id);

    int setnum = convertToInt(string(argv[2]).substr(setnum_id.size()));
    int sigfig = convertToInt(string(argv[3]).substr(sigfig_id.size()));


    //Read in file, process lines
    typedef tokenizer<char_separator<char> > tokenizer;

    ChiSquareALL ch(setnum, sigfig);

    vector<Data> &observeds = ch.observeds;
    Positions    &positions = ch.positions;

    string line;
    while(getline(file, line))
    {
        tokenizer tok(line);

        for (tokenizer::iterator beg=tok.begin(); beg!=tok.end(); false)
        {
            Data obs_line;

            Pos position = convertToPosition(*beg);
            positions.push_back(position);

            while(!((++beg).at_end())){
                int val = convertToPrecision(*beg);
                obs_line.push_back(val);
            }

            observeds.push_back(obs_line);
        }
    }   
    ch.process();

    return 0;
}

