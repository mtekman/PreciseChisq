#ifndef CHISQUARE_H
#define CHISQUARE_H

#include <boost/math/distributions/chi_squared.hpp>
#include <iostream>
#include <vector>

using namespace std;
using namespace boost::math;

static bool printed_yates = false;

class ChiSquare {
private:

    inline long double single_chi(double &obs, double &col_total, long double &set_fract, int num_vals_per_set){
        long double\
                expected = col_total*set_fract,\
                interm_chisquar = obs-expected;

        //Yates correction
        if (num_vals_per_set==2)
            interm_chisquar += interm_chisquar>0?-0.5:0.5;
            if (!printed_yates){
                cerr << "Nb: Applied Yates correction\n" << endl;
                printed_yates=true;
            }

//        cerr << "obs=" << obs << ", set_fract=" << set_fract << endl;
//        cerr << "exp=" << expected << ",  interm=" << interm_chisquar << endl;

        if (expected==0) return -1;

        return (interm_chisquar*interm_chisquar)/expected;
    }



public:
    //Each row contains S sets of data each containing 2 numbers
    ChiSquare(
            vector<vector<int> > &row_data,
            chi_squared_distribution<long double> &dist,
            int &num_vals_per_set, int &sigfig)
    {
        int num_sets = row_data.size();

//        cerr << "num_sets=" << num_sets << endl;

        //1. Get totals
        int total=0;
        vector<long double> set_totals(num_sets,0);
        vector<double> column_totals(num_vals_per_set,0); // {A1+B1+C1, A2+B2+C2, A3+...}

        for (int s=0; s< num_sets; s++){
            int s_total=0;

            for (int c=0; c< num_vals_per_set; c++){
                double val = row_data[s][c];

                column_totals[c] += val;
                s_total += val;
            }
            set_totals[s] = s_total;
            total += s_total;

//            cerr << "s_total=" << s_total << endl;
        }

        //2. Compute fractions
        vector<long double> set_fracts(num_sets,0);
        for (int s=0; s< num_sets; s++)  set_fracts[s] = set_totals[s]/total;

        //3. Compute individual chisquares, add to running total;
        cout << setprecision(sigfig) << flush;

        long double total_chisq=0;
        for (int s=0; s< num_sets; s++){
            for (int c=0; c < num_vals_per_set; c++){
                double val = row_data[s][c];
                cout << val << ' ' << flush;

                total_chisq += single_chi(val, column_totals[c], set_fracts[s], num_vals_per_set);
            }
        }
        cout << "  chisq=" << total_chisq
             << "  p-val=" << (total_chisq>0?cdf(complement(dist,total_chisq)):-1)
             << endl;
//        exit(0);

    }
};

#endif // CHISQUARE_H
