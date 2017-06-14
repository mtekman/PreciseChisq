#ifndef CHISQUARE_H
#define CHISQUARE_H

#include <boost/math/distributions/chi_squared.hpp>
#include <iostream>
#include <vector>

using namespace std;
using namespace boost::math;

typedef long double       Precision;  // good enough..?
typedef vector<Precision> Data;       //


static bool printed_yates = false;


class ChiSquare {
private:
    int num_sets;
public:
    void calcExpecteds(Data &observed, Data &expected)
    {
        int obs_size = observed.size();
        Data set_totals(num_sets, 0);
        Data col_totals(num_sets, 0);

        Precision grand_total = 0;
        int vals_per_set = obs_size / num_sets;


        // Make set and column totals
        /*   e.g. for 2 set observeds: 667  3   1004  1
         *        set_totals = [ 667 + 3   , 1004 + 1]
         *        col_totals = [ 667 + 1004,    3 + 1]
         */
        for (int i=0; i < obs_size; i++)
        {
            int set = i / vals_per_set;
            int col = i % vals_per_set;

            set_totals[set] += observed[i];
            col_totals[col] += observed[i];
            grand_total     += observed[i];
        }

        // Generate expecteds
        for (int i=0; i < obs_size; i++)
        {
            int set = i / vals_per_set,
                col = i % vals_per_set;

            Precision &set_total = set_totals[set],
                        &col_total = col_totals[col];

            Precision expect = set_total * (col_total / grand_total );
            expected[i] = expect;
        }
    }


    void calcChisquareds(Data &observed, Data &expected, Data &chisq)
    {
        for (int i=0; i < observed.size(); i++)
        {
               Precision &obs = observed[i],
                           &exp = expected[i];

               Precision interm_chisquar = obs-exp;
               if(interm_chisquar<0){interm_chisquar *= -1;}

               if (num_sets == 2){
                   interm_chisquar -= 0.5;
                   if (!printed_yates){
                       cerr << "Applied Yates correction\n" << endl;
                       printed_yates=true;
                   }
               }

               Precision chi = (interm_chisquar*interm_chisquar)/exp;
               chisq[i] = chi;
        }
    }


    void grandChiAndPval(Data &chisq, Data &result, chi_squared_distribution<Precision> *dist)
    {
        long double total = 0;
        for (int i=0; i < chisq.size(); i ++){
            total += chisq[i];
        }

        long double pval =(total>0?cdf(complement(*dist,total)):-1);
        result[0] = total;
        result[1] = pval;
    }


    // Deals with individual rows/sets
    ChiSquare(int &numsets,
              chi_squared_distribution<Precision> *dist,
              Data &observed, // size N
              Data &expected, // size N
              Data &chisq,    // size N
              Data &result)   // size 2 (total chisq, pval)
    {
        num_sets = numsets;
        calcExpecteds(observed, expected);
        calcChisquareds(observed, expected, chisq);

        grandChiAndPval(chisq, result, dist);
    }
};

#endif // CHISQUARE_H
