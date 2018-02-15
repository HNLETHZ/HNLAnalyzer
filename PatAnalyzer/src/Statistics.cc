#include "SUSYAnalyzer/PatAnalyzer/interface/Statistics.h"
#include "TMath.h"
using namespace std;

namespace stats {
  double LowBound(int n, int N)
  {
    if (N <= 0 ) return 0.;
    double eff = N ? (double)n/(double)N : 0.;
    double low = N ? TMath::Sqrt( eff*(1.-eff)*(double)N ) / (double)N : 0.;
    

    
    if( eff < 0.5 && n <= 20 )
      {
	switch (n) {
	case 0 :
	  low=0.00;
	  break;
	case 1 :
	  low=eff - 0.37/(double)N;
	  break ;
      case 2 :
	low=eff -0.74/(double)N;
	break;
      case 3 :
	low=eff -1.10/(double)N;
	break;
      case 4 :
	low=eff -2.34/(double)N;
	break;
      case 5 :
	low=eff -2.75/(double)N;
	break ;
      case 6 :
	low=eff -3.82/(double)N;
	break ;
      case 7 :
	low=eff -4.25/(double)N;
	break ;
      case 8 :
	low=eff -5.3/(double)N;
	break ;
      case 9 :
	low=eff -6.33/(double)N;
	break;
      case 10 :
	low=eff -6.78/(double)N;
	break;
      case 11 :
	low=eff -7.81/(double)N;
	break;
      case 12 :
	low=eff -8.83/(double)N;
	break;
      case 13 :
	low=eff -9.28/(double)N;
	break;
      case 14 :
	low=eff -10.30/(double)N;
	break;
      case 15 :
	low=eff -11.32/(double)N;
	break;
      case 16 :
	low=eff -12.33/(double)N;
	break;
      case 17 :
	low=eff -12.79/(double)N;
	break;
      case 18 :
	low=eff -13.81/(double)N;
	break;
      case 19 :
	low=eff -14.82/(double)N;
	break;
      case 20 :
	low=eff -15.83/(double)N;
	break;
      default :
	break;
      }
    }
  /*
  else if (N-n <= 20 ) 
    {
      switch (N-n) {
      case 0 :
	low=1.29/(double)N;
	break;
      case 1 :
	low=2.75/(double)N;
	break ;
      case 2 :
	low=4.25/(double)N;
	break;
      case 3 :
	low=5.30/(double)N;
	break;
      case 4 :
	low=6.78/(double)N;
	break;
      case 5 :
	low=7.81/(double)N;
	break ;
      case 6 :
	low=9.28/(double)N;
	break ;
      case 7 :
	low=10.30/(double)N;
	break ;
      case 8 :
	low=11.32/(double)N;
	break ;
      case 9 :
	low=12.79/(double)N;
	break;
      case 10 :
	low=13.81/(double)N;
	break;
      case 11 :
	low=14.82/(double)N;
	break;
      case 12 :
	low=16.29/(double)N;
	break;
      case 13 :
	low=17.30/(double)N;
	break;
      case 14 :
	low=18.32/(double)N;
	break;
      case 15 :
	low=19.32/(double)N;
	break;
      case 16 :
	low=20.80/(double)N;
	break;
      case 17 :
	low=21.81/(double)N;
	break;
      case 18 :
	low=22.82/(double)N;
	break;
      case 19 :
	low=23.82/(double)N;
	break;
      case 20 :
	low=25.30/(double)N;
	break;
      default :
	break;
	}
	}
  */
  return low;

}


double HighBound(int n, int N)
{
  if (N <= 0 ) return 0.;
  double eff = N ? (double)n/(double)N : 0.;
  double high = N ? TMath::Sqrt( eff*(1.-eff)* (double)N ) / (double)N : 0. ;
  
  if( eff < 0.5 && n <= 20 )
    {
      switch (n) {
      case 0 :
	high=1.29/(double)N - eff;
	break;
      case 1 :
	high=2.75/(double)N - eff;
	break ;
      case 2 :
	high=4.25/(double)N - eff;
	break;
      case 3 :
	high=5.30/(double)N - eff;
	break;
      case 4 :
	high=6.78/(double)N - eff;
	break;
      case 5 :
	high=7.81/(double)N - eff;
	break ;
      case 6 :
	high=9.28/(double)N - eff;
	break ;
      case 7 :
	high=10.30/(double)N - eff;
	break ;
      case 8 :
	high=11.32/(double)N - eff;
	break ;
      case 9 :
	high=12.79/(double)N - eff;
	break;
      case 10 :
	high=13.81/(double)N - eff;
	break;
      case 11 :
	high=14.82/(double)N - eff;
	break;
      case 12 :
	high=16.29/(double)N - eff;
	break;
      case 13 :
	high=17.30/(double)N - eff;
	break;
      case 14 :
	high=18.32/(double)N - eff;
	break;
      case 15 :
	high=19.32/(double)N - eff;
	break;
      case 16 :
	high=20.80/(double)N - eff;
	break;
      case 17 :
	high=21.81/(double)N - eff;
	break;
      case 18 :
	high=22.82/(double)N - eff;
	break;
      case 19 :
	high=23.82/(double)N - eff;
	break;
      case 20 :
	high=25.30/(double)N - eff;
	break;
      default :
	break;
      }
    }
  /*
  else if ( N -n <= 20) 
    {
      switch(N-n){
      case 0 :
	high=0.00;
	break;
      case 1 :
	high=0.37/(double)N;
	break ;
      case 2 :
	high=0.74/(double)N;
	break;
      case 3 :
	high=1.10/(double)N;
	break;
      case 4 :
	high=2.34/(double)N;
	break;
      case 5 :
	high=2.75/(double)N;
	break ;
      case 6 :
	high=3.82/(double)N;
	break ;
      case 7 :
	high=4.25/(double)N;
	break ;
      case 8 :
	high=5.3/(double)N;
	break ;
      case 9 :
	high=6.33/(double)N;
	break;
      case 10 :
	high=6.78/(double)N;
	break;
      case 11 :
	high=7.81/(double)N;
	break;
      case 12 :
	high=8.83/(double)N;
	break;
      case 13 :
	high=9.28/(double)N;
	break;
      case 14 :
	high=10.30/(double)N;
	break;
      case 15 :
	high=11.32/(double)N;
	break;
      case 16 :
	high=12.33/(double)N;
	break;
      case 17 :
	high=12.79/(double)N;
	break;
      case 18 :
	high=13.81/(double)N;
	break;
      case 19 :
	high=14.82/(double)N;
	break;
      case 20 :
	high=15.83/(double)N;
	break;
      default :
	break;
      }
    }
  */

  return high;
}




double LowBound(int n)
{
  if( n < 0 ) return 0;
  double low = TMath::Sqrt(n);
  if( n <= 20 ) 
    {
      switch (n) {
      case 0 :
	low=0.00;
	break;
      case 1:
	low= n - 0.37;
	break ;
      case 2 :
	low=n -0.74;
	break;
      case 3 :
	low=n -1.10;
	break;
      case 4 :
	low=n -2.34;
	break;
      case 5 :
	low=n -2.75;
	break ;
      case 6 :
	low=n -3.82;
	break ;
      case 7 :
	low=n -4.25;
	break ;
      case 8 :
	low=n -5.3;
	break ;
      case 9 :
	low=n -6.33;
	break;
      case 10 :
	low=n -6.78;
	break;
      case 11 :
	low=n -7.81;
	break;
      case 12 :
	low=n -8.83;
	break;
      case 13 :
	low=n -9.28;
	break;
      case 14 :
	low=n -10.30;
	break;
      case 15 :
	low=n -11.32;
	break;
      case 16 :
	low=n -12.33;
	break;
      case 17 :
	low=n -12.79;
	break;
      case 18 :
	low=n -13.81;
	break;
      case 19 :
	low=n -14.82;
	break;
      case 20 :
	low=n -15.83;
	break;
      default :
	break;
      }
    }
  else 
    low = TMath::Sqrt(n);

  return low;
}


double HighBound(int n)
{
  if (n < 0 ) return 0.;
  double high = TMath::Sqrt(n);
  if( n <= 20 ) 
    {
      switch (n) {
      case 0 :
	high=1.29 -n;
	break;
      case 1 :
	high=2.75 -n;
	break ;
      case 2 :
	high=4.25 -n;
	break;
      case 3 :
	high=5.30 -n;
	break;
      case 4 :
	high=6.78 -n;
	break;
      case 5 :
	high=7.81 -n;
	break ;
      case 6 :
	high=9.28 -n;
	break ;
      case 7 :
	high=10.30 -n;
	break ;
      case 8 :
	high=11.32 -n;
	break ;
      case 9 :
	high=12.79 -n;
	break;
      case 10 :
	high=13.81 -n;
	break;
      case 11 :
	high=14.82 -n;
	break;
      case 12 :
	high=16.29 -n;
	break;
      case 13 :
	high=17.30 -n;
	break;
      case 14 :
	high=18.32 -n;
	break;
      case 15 :
	high=19.32 -n;
	break;
      case 16 :
	high=20.80 -n;
	break;
      case 17 :
	high=21.81 -n;
	break;
      case 18 :
	high=22.82 -n;
	break;
      case 19 :
	high=23.82 -n;
	break;
      case 20 :
	high=25.30 -n;
	break;
      default :
	break;
      }
    }
  else
    high = TMath::Sqrt(n);

  return high;
}


}
