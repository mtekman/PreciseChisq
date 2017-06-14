#include "ChiSquare.h"

typedef long unsigned int Pos;
typedef vector<Pos> Positions;



class ChiSquareALL {

public:
    Positions    positions; // N x 1 -- populated from main
    vector<Data> observeds; // N x M -- populated from main
    vector<Data> expecteds; // N X M
    vector<Data> chisqvals; // N x M
    vector<Data> totals;    // N x 2  [[total chisq, pval]]
    int sigfig;
    int numsets;

    chi_squared_distribution<Precision > *dist;


    ChiSquareALL(int num_sets, int sf)
    {
        dist = new chi_squared_distribution<Precision>(num_sets-1);

        sigfig  = sf;
        numsets = num_sets;
    }

    void initialise_tables()
    {
        int N = observeds.size(),
            M = observeds[0].size();

        for (int i=0; i < N; i++){

            expecteds.push_back( Data(M,0) );
            chisqvals.push_back( Data(M,0) );
            totals.push_back(    Data(2,0) );
        }
    }


    void process()
    {
        initialise_tables();

        if (observeds.size() == 0){
            cerr << "Observeds have not yet been populated" << endl;
            exit(-1);
        }

        for (int i=0; i < observeds.size(); i++)
        {
            Data &observed_row = observeds[i],
                 &expected_row = expecteds[i],
                 &chisquar_row = chisqvals[i],
                 &total_row    = totals[i];

            ChiSquare(numsets, dist,
                      observed_row, expected_row, chisquar_row, total_row);
        }

        print_tables();
    }


    void print_tables(){
        int rows = observeds.size(),
            cols = observeds[0].size();

        cout << setprecision(sigfig) << flush;

        for (int i=0; i < rows; i++)
        {
            cout << positions[i] << '\t';

            int j;
            for (j=0; j < cols; j++) cout << observeds[i][j] << '\t';
            cout << "|\t";
            for (j=0; j < cols; j++) cout << expecteds[i][j] << '\t';
            cout << "|\t";
            for (j=0; j < cols; j++) cout << chisqvals[i][j] << '\t';
            cout << "|\t";
            for (j=0; j <    2; j++) cout <<    totals[i][j] << '\t';
            cout << endl;
        }
    }

};
