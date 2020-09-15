#ifndef ALGS_CPP
#define ALGS_CPP

#include "mygraph.cpp"
#include <set>
#include <map>
#include <string>
#include <vector>
#include <fstream>

using namespace std;
using namespace mygraph;

enum Algs {IG=0, TG, RG, SG, BLITS, FIG, FRG, ATG, LATG, ANM};

uniform_real_distribution< double > unidist(0, 1);

resultsHandler allResults;
vector< double > alpha; 

void init_alpha( tinyGraph& g ) {
   alpha.assign( g.n, 0.0 );
   mt19937 gen( 0 ); //same sequence each time
   
   for (node_id u = 0; u < g.n; ++u) {
      alpha[u] = unidist( gen );
   }
}

struct Args {
   Algs alg;
   string graphFileName;
   string outputFileName = "";
   size_t k;
   tinyGraph g;
   double tElapsed;
   double wallTime;
   Logger logg;
   bool steal = true;
   double epsi = 0.1;
   double delta = 0.1;
   size_t N = 1;
   bool fast = false;
   bool reportRounds = false;
};

vector< bool > emptyBoolVector;
vector< bool > emptySetVector;
vector< double > emptyDoubleVector;
vector< size_t > emptySize_tVector;

class MyPair {
public:
   node_id u;
   double  gain; //may be negative

   MyPair() {}
   MyPair( node_id a,
	   int64_t g ) {
      u = a;
      gain = g;
   }

   MyPair( const MyPair& rhs ) {
      u = rhs.u;
      gain = rhs.gain;
   }

   void operator=( const MyPair& rhs ) {
      u = rhs.u;
      gain = rhs.gain;
   }
};



struct gainLT {
   bool operator() (const MyPair& p1, const MyPair& p2) {
      return p1.gain < p2.gain;
   }
} gainLTobj;

struct revgainLT {
   bool operator() (const MyPair& p1, const MyPair& p2) {
      return (p1.gain > p2.gain);
   }
} revgainLTobj;

#ifndef REVMAX
signed long marge( size_t& nEvals, tinyGraph& g, node_id u, vector<bool>& set) {
   
   if (set[u])
      return 0;
   
   ++nEvals;
   
   signed long m;
   double mx = 2 * g.getWeightedDegreeMinusSet(u,set);
   double my = g.getWeightedDegree( u );

   m = (mx - my);
      
   return m;
}

size_t compute_valSet( size_t& nEvals, tinyGraph& g, vector<bool>& set ) {
   ++nEvals;
   size_t val = 0;
   for (node_id u = 0 ; u < g.n; ++u) {
      vector< tinyEdge >& neis = g.adjList[u].neis;
      for (size_t j = 0; j < neis.size(); ++j) {
	 node_id v = neis[j].target;
	 if ( ( set[u] && !set[v] ) || (!set[u] && set[v]) ) 
	    val += neis[j].weight;
      }
   }

   return val / 2;
}
#else
double compute_valSet( size_t& nEvals, tinyGraph& g, vector<bool>& set,
		       vector< bool >& cov = emptySetVector ) {
   if (alpha.size() == 0) {
      init_alpha( g );
   }
   
   ++nEvals;
   cov.assign( g.n, false );
   double val = 0;
   
   for (node_id u = 0 ; u < g.n; ++u) {
      if (!set[u]) {
	 vector< tinyEdge >& neis = g.adjList[u].neis;
	 double valU = 0.0;
	 for (size_t j = 0; j < neis.size(); ++j) {
	    node_id v = neis[j].target;
	    if (set[v]) {
	       valU += neis[j].weight;
	    }
	 }
	 valU = pow( valU, alpha[u] );
	 val += valU;
      }
   }

   return val;
}

double compute_valSet( size_t& nEvals, tinyGraph& g, vector<node_id>& sset ) {
   if (alpha.size() == 0) {
      init_alpha( g );
   }
   vector< bool > set(g.n, false);
   for (size_t i = 0; i < sset.size(); ++i) {
      set[ sset[i] ] = true;
   }
   
   ++nEvals;

   double val = 0;
   for (node_id u = 0 ; u < g.n; ++u) {
      if (!set[u]) {
	 vector< tinyEdge >& neis = g.adjList[u].neis;
	 double valU = 0.0;
	 for (size_t j = 0; j < neis.size(); ++j) {
	    node_id v = neis[j].target;
	    if (set[v]) {
	       valU += neis[j].weight;
	    }
	 }
	 valU = pow( valU, alpha[u] );
	 val += valU;
      }
   }

   return val;
}

double marge( size_t& nEvals, tinyGraph& g, node_id x, vector<bool>& set,
		   vector< bool >& cov = emptySetVector ) {
   if (alpha.size() == 0) {
      init_alpha( g );
   }
   if (set[x])
      return 0;

   double loss = 0.0;
   
   vector< tinyEdge >& neis = g.adjList[x].neis;
   double valX = 0.0;
   for (size_t j = 0; j < neis.size(); ++j) {
      node_id v = neis[j].target;
      if (set[v]) {
	 valX += neis[j].weight;
      }
   }

   valX = pow( valX, alpha[x] );

   loss = valX;

   double gain = 0.0;
   for (size_t j = 0; j < neis.size(); ++j) {
      node_id v = neis[j].target;
      vector< tinyEdge >& neisV = g.adjList[ v ].neis;
      double valV = 0.0;
      double valVwithX = 0.0;
      for (size_t k = 0; k < neisV.size(); ++k) {
	 node_id w = neisV[k].target;
	 if (w != x) {
	    if (set[w]) {
	       valV += neisV[k].weight;
	       valVwithX += neisV[k].weight;
	    }
	 } else {
	    valVwithX += neisV[k].weight;
	 }
      }

      gain += pow( valVwithX, alpha[v] ) - pow( valV, alpha[v] );
   }

   ++nEvals;
   return gain - loss; 
}
#endif


void reportResults( size_t nEvals, double obj, size_t rounds = 0) {
   allResults.add( "obj", obj );
   allResults.add( "nEvals", nEvals );
   allResults.add( "rounds", rounds );
}

void filter( size_t& nEvals,
	     tinyGraph& g,
	     vector< bool >& A,
	     vector< bool >& S,
	     double tau,
	     vector< size_t >& idsA ) {
  vector< size_t > newIdsA;
  for (size_t i = 0; i < idsA.size(); ++i) {
    size_t x = idsA[i];
    if (S[x]) {
      // filter x
      A[x] = false;
    } else {
      if (marge( nEvals, g, x, S ) < tau) {
	//filter x
	A[x] = false;
      } else {
	//keep x
	newIdsA.push_back( x );
      }
    }
  }

  idsA.swap( newIdsA );
}

void random_set( tinyGraph& g, vector< bool >& C, vector< size_t >& A ) {
  C.assign( g.n, false );
  uniform_int_distribution<size_t> dist(0, A.size() - 1);
  size_t size_of_c = dist( gen );
  double prob = 1.0 / size_of_c;
      
  for (size_t i = 0; i < A.size(); ++i) {
    if (unidist(gen) < prob) {
      C[ A[i] ] = true;
    }
  }
}

void random_set( tinyGraph& g, vector< size_t >& C, vector< size_t >& A ) {
  C.clear();
  uniform_int_distribution<size_t> dist(0, A.size() - 1);
  size_t size_of_c = dist( gen );
  double prob = 1.0 / size_of_c;
      
  for (size_t i = 0; i < A.size(); ++i) {
    if (unidist(gen) < prob) {
      C.push_back( A[i] );
    }
  }
}

void unc_max( size_t& nEvals, tinyGraph& g, vector< size_t >& A, double epsi, double delta, vector< size_t >& result ) {
  result.clear();
  size_t ell = (log( 1 /delta ) )/ (log( 1 + (4.0/3)*epsi )) + 1;
  vector< bool > tmp;
  vector< bool > max;
  double tmpVal = 0;
  double maxVal = 0;
  for (size_t i = 0; i < ell; ++i) {
    random_set( g, tmp, A );
    tmpVal = compute_valSet( nEvals, g, tmp );
    if (tmpVal >= maxVal) {
      max = tmp;
      maxVal = tmpVal;
    }
  }

  for (size_t i = 0; i < g.n; ++i) {
    if (max[i])
      result.push_back(i);
    
  }
}

void sampleUt( vector<size_t>& R, vector<size_t>& A, size_t t ) {
  //g.logg << "Starting sampleUX..." << endL;
  R.clear();
  uniform_int_distribution<size_t> dist(0, A.size() - 1);
  vector< bool > alreadySampled( A.size(), false );

  for (size_t i = 0; i < t; ++i) {
    size_t pos;
    do {
      pos = dist( gen );
    } while ((alreadySampled[pos]) );

    alreadySampled[pos] = true;
    R.push_back( A[pos] );
  }
}

unsigned samp_Dt( size_t& nEvals,
		  tinyGraph& g,
		  vector< size_t >& idsA,
		  vector< size_t >& idsS,
		  size_t t,
		  double tau ) {
  vector< size_t > idsT;
  vector< bool > T( g.n, false );
  uniform_int_distribution<size_t> dist(0, idsA.size() - 1);
  while (idsT.size() < t - 1) {
    size_t pos;
    pos = dist( gen );
    if ( !T[ idsA[ pos ] ] ) {
      T[ idsA[ pos ] ] = true;
      idsT.push_back( idsA[ pos ] );
    }
  }

  size_t x = 0;
  do {
    size_t pos = dist( gen );
    x = idsA[ pos ];
  } while ( T[ x ] );

  vector< bool > ScupT(g.n, false);
      
  for (size_t i = 0; i < idsS.size(); ++i) {
    ScupT[idsS[i]] = true;
  }
  for (size_t i = 0; i < idsT.size(); ++i) {
    ScupT[idsT[i]] = true;
  }

  if ( marge( nEvals, g, x, ScupT ) >= tau ) {

    return 1;
  }

  return 0;
}

double reduced_mean_dagum( size_t& nEvals,
			   tinyGraph& g,
			   vector< size_t >& idsA,
			   vector< size_t >& idsS,
			   size_t t,
			   double tau,
			   double epsi,
			   double delta,
			   bool fast = false) {
  g.logg << DEBUG << "RMD: Chernoff requires: " << 16*( log( 2 / delta )/ epsi / epsi + 1 ) << endL;
  g.logg << DEBUG << "t: " << t << endL;
											       
  
  double lambda = exp(1) - 2.0;
  double upsilon = 4.0 * lambda * log(2.0 / delta) / pow(epsi, 2);
  double upsilon_1;
  if (!fast)
    upsilon_1 = 1.0 + (1.0 + epsi) * upsilon;
  else
    upsilon_1 = 100;
  
  size_t sum = 0;
  size_t nSamps = 0;
  size_t goal = ceil(upsilon_1);
											  
  while(sum < goal) {
    nSamps += 1;
    sum += samp_Dt( nEvals, g, idsA, idsS, t, tau );

    if (static_cast<double>( goal ) / nSamps < (1.0 - 1.5*epsi)) {
      break;
    }
  }
  
    g.logg << DEBUG << "Dagum required: " << nSamps << endL;
						       g.logg << DEBUG << "Dagum returning: " << static_cast<double>( goal ) / nSamps << endL;
  return static_cast<double>( goal ) / nSamps;
}

double reduced_mean_chernoff( size_t& nEvals,
			      tinyGraph& g,
			      vector< size_t >& idsA,
			      vector< size_t >& idsS,
			      size_t t,
			      double tau,
			      double epsi,
			      double delta,
			      bool fast = false) {

    g.logg << DEBUG << "Chernoff requires: " << 16*( log( 2 / delta )/ epsi / epsi + 1 ) << endL;
  size_t ell = 0;
  if (fast)
    ell = 100;
  else
    ell = 16*( log( 2 / delta )/ epsi / epsi + 1);

  double ubar = 0;
  for (size_t l = 0; l < ell; ++l) {
    ubar += samp_Dt(nEvals, g, idsA, idsS, t, tau );
	  
  }
  ubar /= ell;

										  
  return ubar;
}
			      




void report_rounds( vector< bool >& S,
		    vector< double >& valRounds,
		    tinyGraph& g ) {
   // size_t tmp = 0;
   // if (valRounds.size() == 0) {
   //    valRounds.push_back( compute_valSet( tmp, g, S ) );
   //    return;
   // }
   
   // size_t tmpVal = valRounds[ valRounds.size() - 1 ];
   // size_t tmpVal2 = compute_valSet( tmp, g, S );
   // if (tmpVal2 > tmpVal) {
   //    valRounds.push_back( tmpVal2 );
   // } else {
   //    valRounds.push_back( tmpVal );
   // }
   valRounds.push_back( 0 );
}

void threshold_sample( size_t& nEvals,
		       tinyGraph& g,
		       vector< bool >& S,
		       size_t k,
		       double tau,
		       double epsi,
		       double delta,
		       vector< bool >& exclude,
		       bool bexcl,
		       vector< size_t >& idsS,
		       bool fast,
		       bool returnA = false,
		       vector< bool >& retA = emptyBoolVector,
		       vector< size_t >& retidsA = emptySize_tVector,
		       bool reportRounds = false,
		       vector< double >& valRounds = emptyDoubleVector) {
  double hatepsi;
  if (!fast) {
    hatepsi = epsi / 3;
  } else {
    hatepsi = epsi;
  }
  g.logg << DEBUG << "hatepsi = " << hatepsi << endL;
  size_t r = log( 2*g.n / delta ) + 1;
  size_t m = log( k ) / hatepsi + 1;
  double hatdelta = delta / (2*r*(m+1));
  size_t adaRounds = 0;

  size_t sizeS = 0;
  vector< bool > A( g.n, true );
  vector< size_t > idsA;
  //vector< size_t > idsS;

  vector< size_t > idsT;
  for (size_t i = 0; i < g.n; ++i) {
    if (bexcl) {
      if (!exclude[i])
	idsA.push_back(i);
      else
	A[i] = false;
    } else
      idsA.push_back( i );
  }

  filter( nEvals, g, A, S, tau, idsA );
  //nothing added during filter
  if (reportRounds) {
     report_rounds( S, valRounds, g );
     ++adaRounds;
  }

  // if (idsA.size() <= log( g.n )) {
  //   can run in sequential mode while maintaining log(n)-
  //   adaptivity
  //   for (size_t j = 0; j < idsA.size(); ++j) {
  //     size_t& x = idsA[j];
  //     if (!S[x]) {
  // 	 signed long margeTmp = marge(nEvals, g, x, S);
  // 	if (margeTmp >= static_cast<signed long>( tau )) {
  // 	  S[x] = true;
  // 	  ++sizeS;
  // 	  idsS.push_back( x );
	  
  // 	} else {

  // 	}

  // 	if (reportRounds) {
  // 	   report_rounds( S, valRounds, g );
  // 	}
	
  // 	if (sizeS >= k) {
  // 	  if (returnA) {
  // 	    retA.swap( A );
  // 	    retidsA.swap( idsA );
  // 	  }
  // 	  return;
  // 	}
  //     }
  //   }

  //   if (returnA) {
  //     retA.swap( A );
  //     retidsA.swap( idsA );
  //   }
    
  //   return;
  // }
  
  for (unsigned j = 0; j < r; ++j) {

    g.logg << DEBUG << "Post filter size of A: " << idsA.size() << endL;
    //	 g.logg << DEBUG << "(1 - hatepsi): " << 1 - hatepsi << endL;
    //  g.logg << DEBUG << "log n = " << log( g.n ) << endL;
						   
    if (idsA.size() == 0) {
      return;
    } else {
      if (idsA.size() == 1) {

	S[ idsA[0] ] = true;
	idsS.push_back( idsA[0] );
	++sizeS;
	
	return;
      }
    }
    double tmpT = 1;//= (1+hatepsi)^i
    size_t t = 0;
    for (unsigned i = 0; i < m; ++i) {
	    
      t = idsA.size();
      if (t > static_cast<size_t>(tmpT))
	t = static_cast<size_t>(tmpT);
      double oldTmpT = tmpT;
      if (fast) {
	 tmpT = tmpT * (1 + hatepsi);//	 tmpT = tmpT * 2;
      } else {
	 tmpT = tmpT * (1 + hatepsi);
      }
      if (tmpT < oldTmpT + 1)
	tmpT = oldTmpT + 1;
      
      double ubar = 0;

      if (t > k - sizeS)
	  break;
      
      if ( t > 1) {
	
	ubar = reduced_mean_dagum( nEvals,
				   g,
				   idsA,
				   idsS,
				   t,
				   tau,
				   hatepsi,
				   hatdelta,
				   fast );

	
	if (ubar <= 1 - 1.5*hatepsi)
	  break;

	if (t == idsA.size()) {
	  //all possible elements are good in expectation
	  //can break
	  break;
	}

	
      } 
    }

    if ( t > k - sizeS )
      t = k - sizeS;

	 
    sampleUt( idsT, idsA, t );
    g.logg << DEBUG << "Adding random set of size: "
	   << t << " " << idsT.size() << endL;
    for (size_t i = 0; i < t; ++i){
      if (!S[ idsT[i] ]) {
	S[ idsT[i] ] = true;
	++sizeS;
	idsS.push_back( idsT[i] );
      }
	    
    }

    if (reportRounds) {
       report_rounds( S, valRounds, g );
       ++adaRounds;
    }

    if (sizeS >= k) {
      if (returnA) {
	retA.swap( A );
	retidsA.swap( idsA );
      }
      return;
    }

    if (j < r - 1) {
       filter( nEvals, g, A, S, tau, idsA );
       if (reportRounds) {
	  report_rounds( S, valRounds, g );
	  ++adaRounds;
       }
    }
  }
  
  if (returnA) {
    retA.swap( A );
    retidsA.swap( idsA );
  }
}


class Ig {
   size_t k;
   tinyGraph& g;
   bool steal;
   size_t nEvals = 0;
public:
   Ig( Args& args ) : g( args.g ) {
      k = args.k;
      steal = args.steal;
   }

   double leastBenefit( node_id u, vector<bool>& set) {
      set[u] = false;
      double m = marge( nEvals, g, u, set );
      set[u] = true;
      return m;
   }
   
 
   void run() {
      vector<bool> A( g.n, false );
      vector<bool> B( g.n, false );
      vector<bool> C( g.n, false );
      vector<bool> D( g.n, false );
      vector<bool> E( g.n, false );

      double valA = 0;
      double valB = 0;
      double valD;
      double valE;
      
      node_id maxSingle;
      for (size_t i = 0; i < k; ++i) {
	 size_t maxAid = 0;
	 double maxMargeA = 0;
	 for (node_id u = 0; u < g.n; ++u) {
	    if (!( B[u] || A[u] )) {
	       if (marge( nEvals, g, u, A ) >= maxMargeA) {
		  maxMargeA = marge( nEvals, g,u, A);
		  maxAid = u;
	       }
	    }
	 }

	 if (maxMargeA > 0) {
	    A[maxAid] = true;
	 
	    valA += maxMargeA;
	 }

	 if (i == 0)
	    maxSingle = maxAid;
	 
	 size_t maxBid = 0;
	 double maxMargeB = 0;
	 for (node_id u = 0; u < g.n; ++u) {
	    if (!( B[u] || A[u] )) {
	       if (marge( nEvals, g, u, B ) >= maxMargeB) {
		  maxMargeB = marge( nEvals, g, u, B );
		  maxBid = u;
	       }
	    }
	 }
	 
	 if (maxMargeB > 0) {
	    B[ maxBid ] = true;
	    valB += maxMargeB;
	 }
      }

      g.logg << "IG: First interlacing complete." << endL;

      //Begin second interlacing
      g.logg << "IG: Adding maxSingle to D,E: " << maxSingle << endL;
      valD = marge( nEvals, g, maxSingle, D );
      valE = valD;
      D[ maxSingle ] = true;
      E[ maxSingle ] = true;
      for (size_t i = 0; i < k - 1; ++i) {
	 size_t maxDid = 0;
	 long maxMargeD = 0;
	 for (node_id u = 0; u < g.n; ++u) {
	    if (!( D[u] || E[u] )) {
	       if (marge( nEvals, g, u, D ) > maxMargeD) {
		  maxMargeD = marge( nEvals, g,u, D);
		  maxDid = u;
	       }
	    }
	 }

	 if (maxMargeD > 0) {
	    D[ maxDid ] = true;
	 
	    valD += maxMargeD;
	 }
	 
	 size_t maxEid = 0;
	 double maxMargeE = 0;
	 for (node_id u = 0; u < g.n; ++u) {
	    if (!( E[u] || D[u])) {
	       if (marge( nEvals, g, u, E ) > maxMargeE) {
		  maxMargeE = marge( nEvals, g, u, E );
		  maxEid = u;
	       }
	    }
	 }
	 if (maxMargeE > 0) {
	    E[ maxEid ] = true;
	    valE += maxMargeE;
	 }

	 	 
      }

      g.logg << "IG: Second interlacing complete." << endL;

      valA = compute_valSet( nEvals,  g, A );
      valB = compute_valSet( nEvals,  g, B );
      valD = compute_valSet( nEvals,  g, D );
      valE = compute_valSet( nEvals,  g, E );
      vector <double> vC;
      vC.push_back( valA );
      vC.push_back( valB );
      vC.push_back( valD );
      vC.push_back( valE );

      double valC = 0;
      size_t jj = 0;
      for (size_t i = 0; i < vC.size(); ++i) {
	 if (vC[i] > valC) {
	    valC = vC[i];
	    jj = i;
	 }
      }
      
      if (jj == 0) 
	 C = A;
      if (jj == 1)
	 C = B;
      if (jj == 2)
	 C = D;
      if (jj == 3)
	 C = E;
      g.logg << "C: " << compute_valSet( nEvals,  g, C ) << endL;

      //steal      
      if (steal) {
	 vector< MyPair > possibleGain;
	 vector< MyPair > Cbenefits;
	 MyPair tmp;
	 for (size_t i = 0; i < g.n; ++i) {
	    if (A[i] || B[i] || D[i] || E[i]) {
	       tmp.u = i;
	       tmp.gain = marge( nEvals, g, i, C );

	       
	       possibleGain.push_back( tmp );
	    }

	    if (C[i]) {
	       tmp.u = i;
	       tmp.gain = leastBenefit(i,C);

	       Cbenefits.push_back ( tmp );
	    }
	 }

	 //g.logg << "IG: Sorting..." << endL;
	 std::sort( Cbenefits.begin(), Cbenefits.end(), gainLT() );
	 std::sort( possibleGain.begin(), possibleGain.end(), revgainLT() );



	 //Attempt to replace elements of C
	 size_t nStolen = 0;
	 for (size_t i = 0; i < Cbenefits.size(); ++i) {

	    if ( Cbenefits[ i ].gain < possibleGain[ i ].gain ) {

	       if ( C[ possibleGain[i].u ] ) {
		  C[ possibleGain[i].u ] = false;
	       } else {
		  ++nStolen;
		  C[ Cbenefits[i].u ] = false;
		  C[ possibleGain[i].u ] = true;
	       }
	    }
	 }
	 g.logg << "IG: Stealing complete: " << nStolen << " stolen." << endL;
      }

      g.logg << "C: " << compute_valSet( nEvals,  g, C ) << endL;
      g.logg << "Evals: " << nEvals << endL;

      reportResults( nEvals, compute_valSet(nEvals, g, C) );
   }
};

class Fig {
   size_t k;
   tinyGraph& g;
   bool steal;
   size_t nEvals = 0;
   double epsi;
   double stopGain;
   bool reportRounds = false;
public:
   Fig( Args& args ) : g( args.g ) {
      k = args.k;
      steal = args.steal;
      epsi = args.epsi;
      reportRounds = args.reportRounds;
   }

   double leastBenefit( node_id u, vector<bool>& set, vector< double >& valRounds ) {
      set[u] = false;
      double m = marge( nEvals, g, u, set );
      set[u] = true;

      if (reportRounds) {
	 report_rounds( set, valRounds, g );
      }
      return m;
   }

   bool swap( node_id u, node_id v, vector<bool>& set) {
      double init = compute_valSet( nEvals, g, set );
      set[u] = false;
      set[v] = true;
      double m = compute_valSet( nEvals, g, set );
      if (m > init) {
	 return true;
      }
      set[u] = true;
      set[v] = false;
      return false;
   }

   size_t sizeSet( vector<bool>& S ) {
      size_t ssize = 0;
      for (size_t i = 0; i < g.n; ++i) {
	 if (S[i])
	    ++ssize;
      }
      return ssize;
   }
   
   void add( vector<bool>& S, vector<bool>& T, node_id& j, double& tau, vector< double >& valRounds ) {
      if (sizeSet( S ) == k) {
	 j = 0;
	 tau = stopGain;
	 return;
      }

      while ( tau > stopGain ) {
	 for (node_id x = j; x < g.n; ++x) {
	    if (!T[x]) {
	       if (reportRounds) {
		  report_rounds( S, valRounds, g );
	       }
	       if (marge(nEvals, g, x, S) >= ( tau )) {
		  S[ x ] = true;
		  j = x;
		  return;
	       }

	       
	    }
	 }
	 tau = ( 1 - epsi ) * tau;
	 j = 0;
      }
      j = 0;
      return;
   }
  
   void run() {
      vector<bool> A( g.n, false );
      vector<bool> B( g.n, false );
      vector<bool> C( g.n, false );
      vector<bool> D( g.n, false );
      vector<bool> E( g.n, false );

      g.logg << "FIG: epsi = " << epsi << ", k = " << k << endL;
      
      //Get max singleton
      g.logg << "FIG: Determining max singleton..." << endL;
      double M = 0;
      node_id a0;
      for (size_t x = 0; x < g.n; ++x) {
	 if ( marge( nEvals, g, x, A ) > (M) ) {
	    a0 = x;
	    M = marge( nEvals, g, x, A );
	 }
      }

      vector< double > valRounds;
      
      g.logg << "FIG: M = " << M << endL;
      g.logg << "FIG: Stopping condition: " << stopGain << endL;
      
      g.logg << "FIG: Starting first interlacing..." << endL;
      double tauA = M;
      double tauB = M;
      stopGain = epsi * M / g.n;
      node_id a = 0;
      node_id b = 0;
      while ( tauA > stopGain || tauB > stopGain) {
	 //g.logg << "FIG: tauA = " << tauA << ", tauB = " << tauB << endL;
	 add( A, B, a, tauA, valRounds );
	 add( B, A, b, tauB, valRounds );
      }

      //g.logg << "FIG: First interlacing complete." << endL;
      g.logg << "FIG: Starting second interlacing..." << endL;
      double tauD = M;
      double tauE = M;
      node_id d = 0;
      node_id e = 0;
      D[ a0 ] = true;
      E[ a0 ] = true;
      while ( tauD > stopGain || tauE > stopGain) {
	 //g.logg << "FIG: tauD = " << tauD << ", tauE = " << tauE << endL;
	 add( D, E, d, tauD, valRounds );
	 add( E, D, e, tauE, valRounds );
      }
      
      double valA = compute_valSet( nEvals,  g, A );
      double valB = compute_valSet( nEvals,  g, B );
      double valD = compute_valSet( nEvals,  g, D );
      double valE = compute_valSet( nEvals,  g, E );
      vector <double> vC;
      vC.push_back( valA );
      vC.push_back( valB );
      vC.push_back( valD );
      vC.push_back( valE );

      double valC = 0;
      size_t jj = 0;
      for (size_t i = 0; i < vC.size(); ++i) {
	 if (vC[i] > valC) {
	    valC = vC[i];
	    jj = i;
	 }
      }
      
      if (jj == 0) 
	 C = A;
      if (jj == 1)
	 C = B;
      if (jj == 2)
	 C = D;
      if (jj == 3)
	 C = E;
      g.logg << "FIG: f(C) = " << compute_valSet( nEvals,  g, C ) << endL;

      //steal      
      if (steal) {
	 vector< MyPair > possibleGain;
	 vector< MyPair > Cbenefits;
	 MyPair tmp;
	 for (size_t i = 0; i < g.n; ++i) {
	    if (A[i] || B[i] || D[i] || E[i]) {
	       tmp.u = i;
	       tmp.gain = marge( nEvals, g, i, C );

	       
	       possibleGain.push_back( tmp );
	    }
	    
	    if (reportRounds) {
	       report_rounds( C, valRounds, g );
	    }
	    
	    if (C[i]) {
	       tmp.u = i;
	       tmp.gain = leastBenefit(i,C, valRounds);
	       
	       Cbenefits.push_back ( tmp );
	    }
	 }

	 std::sort( Cbenefits.begin(), Cbenefits.end(), gainLT() );
	 std::sort( possibleGain.begin(), possibleGain.end(), revgainLT() );



	 //Attempt to replace elements of C
	 size_t nStolen = 0;
	 for (size_t i = 0; i < Cbenefits.size(); ++i) {

	    if ( Cbenefits[ i ].gain < possibleGain[ i ].gain ) {
	       
	       if ( C[ possibleGain[i].u ] ) {
		  C[ possibleGain[i].u ] = false;
	       } else {
		  if (this->swap( Cbenefits[i].u,
				  possibleGain[i].u,
				  C )) {
		     ++nStolen;
		     //C[ Cbenefits[i].u ] = false;
		     //C[ possibleGain[i].u ] = true;
		  }

		  if (reportRounds) {
		     report_rounds( C, valRounds, g );
		  }
	       }
	    }
	 }
	 g.logg << "FIG: Stealing complete: " << nStolen << " stolen." << endL;
	 g.logg << "FIG: f(C) = " << compute_valSet( nEvals,  g, C ) << endL;
      }

      g.logg << "FIG: # evals = " << nEvals << endL;
      reportResults( nEvals, compute_valSet(nEvals, g, C), valRounds.size() );

      if (reportRounds) {
	 for (size_t j = 0; j < valRounds.size(); ++j) {
	    //allResults.add( to_string( j ), valRounds[ j ] );
	 }
      }
      
   }
};


class Rg {
   size_t k;
   tinyGraph& g;
   size_t nEvals = 0;
   bool reportRounds = false;
public:
   Rg( Args& args ) : g( args.g ) {
      k = args.k;
      reportRounds = args.reportRounds;
   }

   double leastBenefit( node_id u, vector<bool>& set ) {
      ++nEvals;
      set[u] = false;
      double m = marge( nEvals, g, u, set );
      set[u] = true;
      return m;
   }
   
   void run() {
      vector<bool> A( g.n, false );
      vector< MyPair > margeGains;
      MyPair tmp;

      vector< double > valRounds;
      for (size_t i = 0; i < k; ++i) {
	 if (reportRounds) {
	    report_rounds( A, valRounds, g );
	 }
	 
	 margeGains.clear();
	 for (node_id u = 0; u < g.n; ++u) {
	    if (!( A[u] )) {
	       tmp.gain = marge( nEvals, g,u, A);
	       tmp.u = u;
	       margeGains.push_back( tmp );
	    }
	 }

	 std::sort( margeGains.begin(), margeGains.end(), revgainLT() );
	 uniform_int_distribution< size_t > dist(0, k - 1);
	 size_t rand = dist( gen );
	 node_id u = margeGains[ rand ].u;
	 A[u] = true;
      }

      g.logg << "A: " << compute_valSet( nEvals,  g, A ) << endL;
      g.logg << "Evals: " << nEvals << endL;

      reportResults( nEvals, compute_valSet(nEvals, g, A), k );

      if (reportRounds) {
	 // for (size_t j = 0; j < valRounds.size(); ++j) {
	 //    allResults.add( to_string( j ), valRounds[ j ] );
	 // }
      }
   }
};

class Frg {
   Args& myArgs;
   size_t k;
   tinyGraph& g;
   size_t nEvals = 0;
   double epsi;
   double w;
   double W;
public:
   Frg( Args& args ) : myArgs( args ), g( args.g ) {
      k = args.k;
      epsi = args.epsi;
   }

   double leastBenefit( node_id u, vector<bool>& set ) {
      ++nEvals;
      set[u] = false;
      double m = marge( nEvals, g, u, set );
      set[u] = true;
      return m;
   }

   void fillM( vector< node_id >& M, vector< bool >& S ) {

      while (w > epsi*W / k) {
	 for (node_id x = 0; x < g.n; ++x) {
	    if (marge( nEvals, g, x, S ) > (1 - epsi)*w) {
	       M.push_back( x );
	       if ( M.size() >= k )
		  return;
	    }
	 }
	 
	 w = (1 - epsi)*w;
      }
   }

   void run() {
      runRandom();
   }

   void randomSampling( double p, size_t s, vector<bool>& A ) {
      size_t rho = p * g.n + 1;
      vector< bool > M;
      vector< MyPair > margeGains;
      MyPair tmp;
      for (size_t i = 0; i < k; ++i) {
	 sampleUnifSize( M, rho );

	 margeGains.clear();
	 for (node_id u = 0; u < g.n; ++u) {
	    if ( M[u] ) {
	       tmp.gain = marge( nEvals, g,u, A);
	       tmp.u = u;
	       margeGains.push_back( tmp );
	    }
	 }

	 std::sort( margeGains.begin(), margeGains.end(), revgainLT() );
	 uniform_int_distribution< size_t > dist(0, s - 1);
	 size_t rand = dist( gen );
	 node_id u = margeGains[ rand ].u;
	 if ( marge( nEvals, g, u, A ) >= 0.0 ) {
	    A[u] = true;
	 }
      }
   }

   void sampleUnifSize( vector<bool>& R, size_t Size ) {

      uniform_int_distribution<size_t> dist(0, g.n-1);
      R.assign(g.n, false);

      for (size_t i = 0; i < Size; ++i) {
	 size_t pos;
	 do {
	    pos = dist( gen );
	 } while ( R[ pos ] );
	 R[ pos ] = true;
      }
    
   }

   void runRandom() {
      double p = 8.0 / (k * epsi * epsi) * log( 2 / (epsi) );
      g.logg << "FastRandom: p = " << p << endL;
      if (p >= 1.0) {
	 //run RandomGreedy
	 g.logg << "FastRandom: Running RandomGreedy..." << endL;
	 Rg rg( myArgs );
	 rg.run();
	 
      } else {
	 g.logg << "FastRandom: Running RandomSampling..." << endL;
	 vector<bool> S(g.n, false );
	 size_t rho = p * g.n + 1;
	 size_t s = rho * k / g.n;

	 randomSampling( p, s, S );

	 g.logg << "S: " << compute_valSet( nEvals,  g, S ) << endL;
	 g.logg << "Evals: " << nEvals << endL;

	 reportResults( nEvals, compute_valSet(nEvals, g, S), k );
      }
   }   
   
   void runSimple() {
      vector<node_id> M;
      vector<node_id> newM;
      vector<bool> S( g.n, false );

      //Get max singleton
      g.logg << "FRG: Determining max singleton..." << endL;
      W = 0;


      for (size_t x = 0; x < g.n; ++x) {
	 if ( marge( nEvals, g, x, S ) > W ) {
	    W = marge( nEvals, g, x, S );
	 }
      }
      w = W;
      
      for (size_t i = 0; i < k; ++i) {
	 g.logg << "FRG: iteration " << i << endL;
	 fillM( M, S );
	 uniform_int_distribution< size_t > dist(0, k - 1 );
	 size_t rand = dist( gen );
	 if (rand >= M.size())
	    continue;
	 
	 node_id u = M[ rand ];
	 S[u] = true;
	 //M.erase( M.begin() + u );
	 for (size_t j = 0; j < M.size(); ++j) {
	    if ( marge( nEvals, g, M[j], S ) <= (1 - epsi)*w ) {
	       //
	    } else {
	       newM.push_back( M[j] );
	    }
	 }
	 M = newM;
      }

      g.logg << "S: " << compute_valSet( nEvals,  g, S ) << endL;
      g.logg << "Evals: " << nEvals << endL;

      reportResults( nEvals, compute_valSet(nEvals, g, S) );
   }

   void runImproved() {
      vector<node_id> M;
      vector<bool> S( g.n, false );

      //Get max singleton
      g.logg << "FRG: Determining max singleton..." << endL;
      W = 0;


      for (size_t x = 0; x < g.n; ++x) {
	 if ( marge( nEvals, g, x, S ) > static_cast<signed long>(W) ) {
	    W = marge( nEvals, g, x, S );
	 }
      }

      //w = W;

      fillM( M, S );
      
      for (size_t i = 0; i < k; ++i) {
	 uniform_int_distribution< size_t > dist(0, M.size() - 1 );
	 size_t rand = dist( gen );
	 node_id u = M[ rand ];
	 if ( marge( nEvals, g, u, S ) > (1 - epsi)*w ) {
	    S[u] = true;
	 } else {
	    vector< node_id > newM;
	    for (size_t j = 0; j < M.size(); ++j) {
	       if ( marge( nEvals, g, M[j], S ) <= (1 - epsi)*w ) {
		  //
	       } else {
		  newM.push_back( M[j] );
	       }
	    }
	    M = newM;
	    size_t sizeOld = M.size();
	    fillM( M, S );
	    size_t sizeInc = M.size() - sizeOld;
	    if (sizeInc > 0) {
	       uniform_int_distribution< size_t > dist(0, sizeInc - 1 );
	       size_t rand = dist( gen );
	       node_id u = M[ sizeOld + rand ];
	       S[u] = true;
	    }
	 }
      }

      g.logg << "S: " << compute_valSet( nEvals,  g, S ) << endL;
      g.logg << "Evals: " << nEvals << endL;

      reportResults( nEvals, compute_valSet(nEvals, g, S) );
   }
};

class Sg {
   size_t k;
   tinyGraph& g;
   size_t nEvals = 0;
   bool reportRounds = false;
public:
   Sg( Args& args ) : g( args.g ) {
      k = args.k;
      reportRounds = args.reportRounds;
   }

   long leastBenefit( node_id u, vector<bool>& set ) {
      set[u] = false;
      long m = marge( nEvals, g, u, set );
      set[u] = true;
      return m;
   }
   
   
   void run() {
      vector<bool> A( g.n, false );

      double maxGain;
      node_id maxIdx;
      MyPair tmp;

      vector< double > valRounds;

      for (size_t i = 0; i < k; ++i) {
	 if (reportRounds) {
	      report_rounds( A, valRounds, g );
	   }
	 
	 maxGain = 0;
	 for (node_id u = 0; u < g.n; ++u) {
	    
	    if (marge( nEvals, g,u,A) > maxGain) {
	       maxIdx = u;
	       maxGain = marge( nEvals, g,u,A);
	    }
	 }

	 if (maxGain > 0) {
	    A[maxIdx] = true;
	 } else {
	    break;
	 }
      }

      g.logg << "A: " << compute_valSet( nEvals,  g, A ) << endL;
      g.logg << "Evals: " << nEvals << endL;

      reportResults( nEvals, compute_valSet(nEvals, g, A) );

      if (reportRounds) {
	 for (size_t j = 0; j < valRounds.size(); ++j) {
	    allResults.add( to_string( j ), valRounds[ j ] );
	 }
      }
   }
};

class Anm {
public:
   size_t k;
   tinyGraph& g;
   double epsi;
   double delta;
   size_t r;
   double OPT; //guess for opt
   size_t nEvals = 0;
   bool reportRounds = false;
   
  bool fast = false;
  
   Anm( Args& args ) : g( args.g ) {
      k = args.k;
      //OPT = 130180;
      OPT = g.m;
      epsi = args.epsi;
      delta = args.delta;
      fast = args.fast;
      reportRounds = args.reportRounds;
      g.logg << "AdaptiveNonmonotoneMax initialized:" << endL;
      g.logg << "epsi=" << epsi << endL;
      g.logg << "delta=" << delta << endL;
      g.logg << "k=" << k << endL;
      if (fast) {
	 g.logg << WARN << "Fast mode enabled. Theoretical guarantees will not hold!" << endL << INFO;
	
      }
   }

   size_t get_size_set( vector< bool >& S ) {
      size_t sizeSol = 0;
      for (size_t i = 0; i < g.n; ++i)
	 if (S[i])
	    ++sizeSol;
      return sizeSol;
   }
   
  void run() {
      g.logg << "ANM starting run..." << endL;
      vector<bool> sol( g.n, false );
      double solVal = 0;

      vector< MyPair > margeGains;
      MyPair tmp;

      margeGains.clear();
      for (node_id u = 0; u < g.n; ++u) {
	 tmp.gain = marge( nEvals, g,u, sol);
	 tmp.u = u;
	 margeGains.push_back( tmp );
      }

      double topKgains = 0;
      std::sort( margeGains.begin(), margeGains.end(), revgainLT() );
      for (size_t i = 0; i < k; ++i) {
	 topKgains += margeGains[ i ].gain;
      }

      double M = margeGains[0].gain; //topKgains / k;

      size_t r = 0;
      if (!fast) {
	epsi = epsi / 6.0;
	r = log(  k ) / epsi + 1;
	delta = delta / (2*(r + 1));
      } else {
	r = log(  k ) / epsi + 1;
      }

      double c_1 = 1.0 / 7;
      double c_3 = 3;
      
      g.logg << "topKgains = " << topKgains << endL;
      g.logg << "M = " << M << endL;


      double tau_i = M / (1 - epsi) *c_1;

      g.logg << "r = " << r << endL;
      double solRound = 0;
      vector< double > valueTmp;      
      vector< vector< double > > valueAllRounds( r + 1, valueTmp );
     
      for (unsigned i = 0; i <= r; ++i) {
	 vector< double >& valRounds = valueAllRounds[i];
	 valRounds.clear();
	 valRounds.push_back( 0 );
	 
	vector<bool> A( g.n, false );
	vector<bool> S( g.n, false );
	vector<size_t> idsA;
	vector<size_t> idsS;

	threshold_sample( nEvals, g, S, k, tau_i, epsi, delta, A, false, idsS, fast, true, A, idsA, reportRounds, valRounds );
	double tempVal = compute_valSet( nEvals, g, S );
	if ( tempVal >= solVal ) {
	    solVal = tempVal;
	    sol = S;
	    solRound = i;
	 }
	
	tau_i = tau_i * (1 - epsi);

	size_t sizeA = idsA.size();

	if (sizeA < c_3 * k) {
	  vector< size_t > idsU;
	  random_set(g, idsU, idsA );
	  //unc_max( nEvals, g, idsA, epsi, delta, idsU );
	  if (reportRounds)
	     report_rounds( sol, valRounds, g );
	  if (idsU.size() > k) {
	    vector< size_t > idsD;
	    sampleUt( idsD, idsU, k );
	    idsU = idsD;
	  }

	  random_shuffle ( idsU.begin(), idsU.end() );

	  vector <bool> empty( g.n, false );
	  vector <bool> prefix( g.n, false );
	  double tempVal = 0;
	  
	  for (size_t j = 0; j < idsU.size(); ++j) {
	    prefix[ idsU[j] ] = true;

	    tempVal = compute_valSet( nEvals, g, prefix );

	    if ( tempVal >= solVal ) {
	      solVal = tempVal;
	      sol = prefix;
	      solRound = i;
	    }
	  }

	  if (reportRounds) {
	     report_rounds( sol, valRounds, g );
	  }

	}
      }

      g.logg << INFO << "ANM: solVal=" << solVal << endL;
      g.logg << INFO << "ANM: queries=" << nEvals << endL;
      g.logg << INFO << "ANM: solSize=" << get_size_set( sol ) << endL;


      size_t rounds = 0;
      if (reportRounds) {
	 g.logg << INFO << "ANM: Adaptive rounds=" << valueAllRounds[ solRound ].size() << endL;
	 rounds = valueAllRounds[ solRound ].size();
	 for (size_t j = 0; j < valueAllRounds[ solRound ].size(); ++j) {
	    allResults.add( to_string( j ), valueAllRounds[ solRound ][ j ] );
	 }
      }

      reportResults( nEvals, solVal, rounds );
   }
};

class Atg {
public:
   size_t k;
   tinyGraph& g;
   double epsi;
   double delta;
   size_t r;
   double OPT; //guess for opt
   size_t nSamps = 30;
   size_t nEvals = 0;
  bool fast = false;
   bool reportRounds = false;
   
   Atg( Args& args ) : g( args.g ) {
      k = args.k;
      //OPT = 130180;
      OPT = g.m;
      epsi = args.epsi;
      delta = args.delta;
      fast = args.fast;
      reportRounds = args.reportRounds;
      g.logg << "ATG initialized:" << endL;
      g.logg << "epsi=" << epsi << endL;
      g.logg << "delta=" << delta << endL;
      g.logg << "k=" << k << endL;
      if (fast) {
	 g.logg << WARN << "Fast mode enabled. Theoretical guarantees will not hold!" << endL << INFO;
	 
      } else {
	delta = delta / 2;
      }
   }

   size_t get_size_set( vector< bool >& S ) {
      size_t sizeSol = 0;
      for (size_t i = 0; i < g.n; ++i)
	 if (S[i])
	    ++sizeSol;
      return sizeSol;
   }
   
  void run() {
      g.logg << "ATG: starting run..." << endL;
      vector<bool> A( g.n, false );
      vector<size_t> idsA;
      vector<size_t> idsB;
      vector<bool> B( g.n, false );
      vector<bool> C( g.n, false );
      vector<bool> sol( g.n, false );
      size_t solVal = 0;
      size_t solRound = 0;

      //g.logg.set_level( DEBUG );
      
      // //Get max singleton
      // g.logg << "ATG: Determining max singleton..." << endL;
      // size_t M = 0;
      // //node_id a0;

      // for (size_t x = 0; x < g.n; ++x) {
      // 	 if ( marge( nEvals, g, x, A ) > static_cast<signed long>(M) ) {
      // 	    //a0 = x;
      // 	    M = marge( nEvals, g, x, A );
      // 	 }
      // }

      vector< MyPair > margeGains;
      MyPair tmp;

      margeGains.clear();
      for (node_id u = 0; u < g.n; ++u) {
	 tmp.gain = marge( nEvals, g,u, A);
	 tmp.u = u;
	 margeGains.push_back( tmp );
      }

      double topKgains = 0;
      std::sort( margeGains.begin(), margeGains.end(), revgainLT() );
      for (size_t i = 0; i < k; ++i) {
	 topKgains += margeGains[ i ].gain;
      }

      double M = topKgains / k;

      g.logg << "ATG: topKgains = " << topKgains << endL;
      g.logg << "ATG: M = " << M << endL;

      size_t m = log( 1.0 / (6 * k) ) / log( 1 - epsi );
      //      size_t m = 1;

      double tau_i = M / (1 - epsi);

      g.logg << "ATG: m = " << m << endL;

      vector< double > valueTmp;      
      vector< vector< double > > valueAllRounds( m + 1, valueTmp );

      for (unsigned i = 0; i <= m; ++i) {
	 vector< double_t >& valueThisRound = valueAllRounds[i];
	 valueThisRound.clear();
	 valueThisRound.push_back( 0 ); //Finding topKgains was first adaptive round
	 
	 g.logg << DEBUG << "ATG: iteration i = " << i << endL;

	 tau_i = tau_i * (1 - epsi);
	 g.logg << "ATG: iteration tau_i = " << tau_i << endL;

	 g.logg << "ATG: current solVal = " << solVal << endL;
	 g.logg << "ATG: Current size of solution = " << get_size_set( sol ) << endL;

	g.logg << "ATG: Stopping condition: " << solVal * (1 - epsi) / (6*k) << endL << INFO;
	
	if (tau_i < solVal * (1 - epsi) / (6*k)) {
	  //solVal is a lower bound on OPT and
	  //telling us we can stop now.
	  break;
	}

	A.assign( g.n, false );
	idsA.clear();
	threshold_sample( nEvals, g, A, k, tau_i, epsi, delta, A, false, idsA, fast, false, emptyBoolVector, emptySize_tVector, reportRounds, valueThisRound );

	B.assign( g.n, false );
	idsB.clear();
	threshold_sample( nEvals, g, B, k, tau_i, epsi, delta, A, true , idsB, fast, false, emptyBoolVector, emptySize_tVector, reportRounds, valueThisRound );
	
	random_set( g,C, idsA );

	double tempVal = compute_valSet( nEvals, g, A );

	 g.logg << DEBUG << "f(A)=" << tempVal << endL;
	 if ( tempVal >= solVal ) {
	    solVal = tempVal;
	    sol = A;

	    solRound = i;
	 }

	 tempVal = compute_valSet( nEvals, g, B );
	 g.logg << DEBUG << "f(B)=" << tempVal << endL;
	 if ( tempVal >= solVal ) {
	    solVal = tempVal;
	    sol = B;

	    solRound = i;
	 }

	 tempVal = compute_valSet( nEvals, g, C );
	 g.logg << DEBUG << "f(C)=" << tempVal << endL;
	 if ( tempVal >= solVal ) {
	    solVal = tempVal;
	    sol = C;
	    solRound = i;
	 }

	 if (get_size_set( sol ) == k ) {
	    //We can quit now.
	    break;
	 }

	 
	 // if (!fast) {
	 //    if (solVal > topKgains / 6.0) {
	 //       break;
	 //    }
	 // }
	 
      }

      g.logg << INFO << "ATG: solVal=" << solVal << endL;
      g.logg << INFO << "ATG: queries=" << nEvals << endL;
      g.logg << INFO << "ATG: solSize=" << get_size_set( sol ) << endL;

      size_t rounds = 0;
      


      if (reportRounds) {
	 g.logg << INFO << "ATG: Adaptive rounds=" << valueAllRounds[ solRound ].size() << endL;
	 rounds = valueAllRounds[ solRound ].size();
	 for (size_t j = 0; j < valueAllRounds[ solRound ].size(); ++j) {
	    allResults.add( to_string( j ), valueAllRounds[ solRound ][ j ] );
	 }
      }

      reportResults( nEvals, solVal, rounds );
  }
};

class Latg {
public:
   size_t k;
   tinyGraph& g;
   double epsi;
   double delta;
   size_t nEvals = 0;
  bool fast = false;
  size_t c = 0;
   bool reportRounds = false;
   Latg( Args& args ) : g( args.g ) {
     k = args.k;
      epsi = args.epsi;
      delta = args.delta;
      fast = args.fast;
      reportRounds = args.reportRounds;
      g.logg << "LATG initialized:" << endL;
      g.logg << "epsi=" << epsi << endL;
      g.logg << "delta=" << delta << endL;
      g.logg << "k=" << k << endL;
      if (fast) {
	g.logg << WARN << "Fast mode enabled. Theoretical guarantees will not hold!" << endL;

	c = 1 / epsi;
      } else {
	c = 8 / epsi;
	epsi = 0.63 * epsi / 8;
	delta = delta / 2;
      }
   }


   size_t get_size_set( vector< bool >& S ) {
      size_t sizeSol = 0;
      for (size_t i = 0; i < g.n; ++i)
	 if (S[i])
	    ++sizeSol;
      return sizeSol;
   }
   
   void run() {
     g.logg << INFO << "LATG is starting run..." << endL;
      vector<bool> A( g.n, false );
      vector<size_t> idsA;
      vector<size_t> idsB;
      vector<bool> B( g.n, false );
      vector<bool> C( g.n, false );
      vector<bool> sol( g.n, false );
      double solVal = 0;

      //Get max singleton
      g.logg << "Determining max singleton M..." << endL;
      double M = 0;
      //node_id a0;

      for (size_t x = 0; x < g.n; ++x) {
      	 if ( marge( nEvals, g, x, A ) > (M) ) {
      	    M = marge( nEvals, g, x, A );
      	 }
      }

      g.logg << "M = " << M << endL;
      double d = 1.0;
      if (fast)
	 d = 1.0;
      
      size_t m = log( 1.0 / (c * k) ) / log( pow(1 - epsi, d) );

      double tau_i = M / (1 - epsi);

      g.logg << DEBUG << "m = " << m << endL << INFO;

      if (!fast) {
	 delta = delta / m;
      }
      
      vector< double > valueThisRound;

      valueThisRound.push_back( 0 ); //Finding topKgains was first adaptive round
      for (unsigned i = 0; i <= m; ++i) {
	 tau_i = tau_i * pow(1 - epsi, d);
	 
	threshold_sample( nEvals, g, A, k - idsA.size(),
			  tau_i,
			  epsi, delta,
			  A,
			  false,
			  idsA, fast,
			  false,
			  emptyBoolVector,
			  emptySize_tVector,
			  reportRounds,
			  valueThisRound );

	if (idsA.size() == k)
	  break;
      }

      tau_i = M / (1 - epsi);

      for (unsigned i = 0; i <= m; ++i) {

	 tau_i = tau_i * pow(1 - epsi, d);
	 
	threshold_sample( nEvals, g, B, k - idsB.size(),
			  tau_i,
			  epsi, delta,
			  A,
			  true,
			  idsB, fast,
			  false,
			  emptyBoolVector,
			  emptySize_tVector,
			  reportRounds,
			  valueThisRound );

	if (idsB.size() == k)
	  break;
      }

      
      random_set( g, C, idsA );

      double tempVal = compute_valSet( nEvals, g, A );

      g.logg << DEBUG << "f(A)=" << tempVal << endL;
      if ( tempVal >= solVal ) {
	solVal = tempVal;
	sol = A;

      }

      tempVal = compute_valSet( nEvals, g, B );
      g.logg << DEBUG << "f(B)=" << tempVal << endL;
      if ( tempVal >= solVal ) {
	solVal = tempVal;
	sol = B;

      }

      tempVal = compute_valSet( nEvals, g, C );
      g.logg << DEBUG << "f(C)=" << tempVal << endL;
      if ( tempVal >= solVal ) {
	solVal = tempVal;
	sol = C;

      }

      g.logg << INFO << "LATG: solVal=" << solVal << endL;
      g.logg << INFO << "LATG: queries=" << nEvals << endL;
      g.logg << INFO << "LATG: solSize=" << get_size_set( sol ) << endL;

      size_t rounds = 0;
      

      if (reportRounds) {

	 g.logg << INFO << "LATG: Adaptive rounds=" << valueThisRound.size() << endL;
	 rounds = valueThisRound.size();
	 //for (size_t j = 0; j < valueThisRound.size(); ++j) {
	 //allResults.add( to_string( j ), valueThisRound[ j ] );
	 //}
      } else {
	 
      }

      reportResults( nEvals, solVal, rounds );
   }
};



class Blits {
   size_t k;
   tinyGraph& g;
   double epsi;
   size_t r;
   double OPT; //guess for opt
   size_t nSamps = 100;
   size_t nEvals = 0;
public:
   Blits( Args& args ) : g( args.g ) {
      k = args.k;
      //OPT = 130180;
      OPT = g.m;
      epsi = args.epsi;
      r = 10; //20/ epsi * log( g.n ) / log ( 1 + epsi /2 );
      g.logg << "Blits initialized, r = " << r << endL;
   }

   size_t sizeSet( vector<bool>& S ) {
      size_t ssize = 0;
      for (size_t i = 0; i < g.n; ++i) {
	 if (S[i])
	    ++ssize;
      }
      return ssize;
   }
   
   void sampleUX( vector<bool>& R, vector<bool>& X ) {
      //g.logg << "Starting sampleUX..." << endL;
      uniform_int_distribution<size_t> dist(0, g.n-1);
      R.assign(g.n, false);
      for (size_t i = 0; i < k / r; ++i) {
	 size_t pos;
	 do {
	    pos = dist( gen );
	 } while (!X[ pos ] || R[ pos ] );
	 R[ pos ] = true;
      }
      //g.logg << "Finished sampleUX." << endL;
   }
   
   double Delta( node_id a, vector<bool>& S, vector<bool>& X ) {
      double avg = 0.0;
      for (size_t i = 0; i < nSamps; ++i) {
	 vector<bool> R;
	 if (sizeSet(X) >= k/r)
	    sampleUX( R, X );
	 else
	    R = X;
	 R[a] = false;
	 vector<bool> RcupS;
	 setunion( RcupS, R, S );
	 //RcupS[a] = false;
	 avg += marge( nEvals, g, a, RcupS );
      }
      
      return avg / nSamps;
   }

   double exMarge( vector<bool>& S, vector<bool>& Xplus, vector<bool>& X ) {
      double sum = 0.0;
      for (size_t i = 0; i < nSamps; ++i) {
	 vector<bool> R;
	 sampleUX( R, X );
	 vector<bool> base;
	 vector<bool> larger;
	 setintersection( base, R, Xplus );
	 int64_t valBase = compute_valSet( nEvals,  g, S );
	 setunion( larger, S, base );
	 int64_t valLarger = compute_valSet( nEvals,  g, larger );
	 sum += (valLarger - valBase);
      }
      
      return sum / nSamps;
   }
   
   bool sieve( vector<bool>& res, vector<bool>& A, size_t i, size_t& rounds ) {
      g.logg << "Starting sieve, iteration " << i << endL;
      vector<bool> X( g.n, true );
      res.assign( g.n, false );
      size_t sizeX = g.n;
      double base = 1.0 - 1 / static_cast<double>(r);
      double t = (0.5 - epsi / 4)*( pow( base, i - 1)*(1 - epsi / 2.0) * OPT - compute_valSet( nEvals,  g, A ));
      //if (t < 0)
      //return false;
      g.logg << "t = " << t << endL;
      size_t sizeXprior = 0;
      while (sizeX > k) {
	 g.logg << "Size of X: " << sizeX << endL;
	 vector<bool> Xplus(g.n, false);
	 for (size_t a = 0; a < g.n; ++a) {
	    if (Delta(a, A, X) >= 0.0)
	       Xplus[a] = true;
	 }
	 ++rounds;
	 
	 if (sizeXprior == sizeX) {
	    g.logg << WARN <<"SizeX is not decreasing..." << endL;
	    g.logg << INFO;

	    vector<bool> R;
	    sampleUX( R, X );
	    setintersection( res, R, Xplus );
	    return true;

	    //	    return false;
	 }

	 double exMar = exMarge( A, Xplus, X );
	 ++rounds;
	 g.logg << INFO <<"exMarge: " << exMar << endL;
	 g.logg << INFO <<"t/r: " << t/r << endL;
	 if ( exMar >= t / r ) {
	    vector<bool> R;
	    sampleUX( R, X );
	    setintersection( res, R, Xplus );
	    return true;
	 }

	 Xplus.assign(g.n, false);
	 sizeXprior = sizeX;
	 sizeX = 0;
	 for (size_t a = 0; a < g.n; ++a) {
	    if (Delta(a, A, X) >= (1 + epsi / 4)*t/k ) {
	       Xplus[a] = true;
	       ++sizeX;
	    }
	 }
	 X = Xplus;
	 ++rounds;
      }
      
      g.logg << "Size of X: " << sizeX << endL;
      
      vector<bool> Xplus(g.n, false);
      for (size_t a = 0; a < g.n; ++a) {
	 if (Delta(a, A, X) >= 0.0)
	    Xplus[a] = true;
      }
      vector<bool> R;
      ++rounds;
      if (sizeSet(X) >= k / r)
	 sampleUX( R, X );
      else
	 R = X;

      setintersection( res, R, Xplus );
      return true;
   }
   
   double leastBenefit( node_id u, vector<bool>& set ) {
      set[u] = false;
      double m = marge( nEvals, g, u, set );
      set[u] = true;
      return m;
   }
   
      void setunion( vector<bool>& res, vector<bool>& set1, vector<bool>& set2) {
      res.assign(g.n , false);
      for (node_id u = 0; u < g.n; ++u) {
	 if (set1 [u] || set2[u] ) {
	    res[u] = true;
	 }
      }
   }

   void setintersection( vector<bool>& res, vector<bool>& set1, vector<bool>& set2) {
      res.assign(g.n , false);
      for (node_id u = 0; u < g.n; ++u) {
	 if (set1 [u] && set2[u] ) {
	    res[u] = true;
	 }
      }
   }
   
   void run() {
      if (r > k) {
	 g.logg << ERROR << "r > k" << endL;
	 g.logg << "Exiting Blits." << endL;
	 g.logg << INFO;
	 return;
      }

      vector<bool> S(g.n,false);
      
      vector< MyPair > margeGains;
      MyPair tmp;

      margeGains.clear();
      for (node_id u = 0; u < g.n; ++u) {
	 tmp.gain = marge( nEvals, g,u, S);
	 tmp.u = u;
	 margeGains.push_back( tmp );
      }

      double topKgains = 0;
      std::sort( margeGains.begin(), margeGains.end(), revgainLT() );
      for (size_t i = 0; i < k; ++i) {
	 topKgains += margeGains[ i ].gain;
      }

      OPT = topKgains; //upper bound on OPT
      
      vector< vector<bool > > vSols;
      vector <size_t> roundsByGuess;
      size_t rounds;
      while ( OPT > compute_valSet( nEvals,  g, S ) ) {
	 g.logg << "Starting Blits, with OPT estimate: " << OPT << endL;
	 rounds = 1; //Guess OPT in parallel
	 S.assign(g.n, false);
	 for (size_t i = 1; i <= r; ++i) {
	    vector<bool> step;
	    if (!sieve( step, S, i, rounds ))
	       break;
	    vector<bool> Splus;
	    setunion ( Splus, S, step );
	    S = Splus;
	 }

	 vSols.push_back( S );
	 OPT = OPT*(1 - epsi);
	 roundsByGuess.push_back( rounds );
      }
      
      vector<bool> Smax( g.n, false );
      double valS = 0;
      for (size_t i = 0; i < vSols.size(); ++i) {
	 if (compute_valSet( nEvals,  g, vSols[i] ) > valS) {
	    valS = compute_valSet( nEvals,  g, vSols[i] );
	    Smax = vSols[i];
	    rounds = roundsByGuess[i];
	 }
      }

      g.logg << "S: " << valS << endL;
      g.logg << "Evals: " << nEvals << endL;
      g.logg << "Rounds: " << rounds << endL;

      reportResults( nEvals, valS, rounds );
   }
   
};

class Tg {
   size_t k;
   tinyGraph& g;
   size_t nEvals = 0;
   bool reportRounds = false;
public:
   Tg( Args& args ) : g( args.g ) {
      k = args.k;
      reportRounds = args.reportRounds;
   }

   double leastBenefit( node_id u, vector<bool>& set ) {
      set[u] = false;
      double m = marge( nEvals, g, u, set );
      set[u] = true;
      return m;
   }
   

   void uncMax( vector<bool>& set, vector< double >& valRounds ) {
      vector<bool> all(set);
      vector<bool> none( g.n, false );
      
      for (size_t u = 0; u < g.n; ++u) {
	 if (set[u]) {
	    if (reportRounds) {
	       report_rounds( none, valRounds, g );
	    }
	    double margeA = marge( nEvals, g, u, none ); //adding u to A
	    vector<bool> allMinus ( all );
	    allMinus[u] = false;
	    double margeB = static_cast<double>(compute_valSet(nEvals, g, allMinus)) - compute_valSet(nEvals,g, all );
	    if (margeA >= margeB) {
	       none[u] = true;
	    } else {
	       all[u] = false;
	    }
	 }
      }

      set = all;
   }
   
   void run() {
      vector<bool> A( g.n, false );

      double maxGain;
      node_id maxIdx;
      MyPair tmp;

      vector< double > valRounds;
      
      for (size_t i = 0; i < k; ++i) {
	 if (reportRounds) {
	    report_rounds( A, valRounds, g );
	 }
	 maxGain = -1.0 * g.n * g.n;
	 for (node_id u = 0; u < g.n; ++u) {

	    if (marge( nEvals, g,u,A) > maxGain) {

	       maxIdx = u;
	       maxGain = marge( nEvals, g,u,A);
	    }
	 }

	 A[maxIdx] = true;
      }

      vector<bool> B( g.n, false );
      for (size_t i = 0; i < k; ++i) {
	 if (reportRounds) {
	    report_rounds( A, valRounds, g );
	 }
	 maxGain = -1.0 * g.n * g.n;
	 for (node_id u = 0; u < g.n; ++u) {
	    if (! A[ u ] ) {
	       if ( marge( nEvals, g,u,B) > maxGain ) {
		  maxIdx = u;
		  maxGain = marge( nEvals, g,u,B);
	       }
	    }
	 }

	 B[maxIdx] = true;
      }

      vector<bool> Amax = A;
      vector<bool> Bmax = B;

      uncMax( Amax, valRounds );
      uncMax( Bmax, valRounds );

      vector< double > vals;
      vals.push_back( compute_valSet( nEvals,  g, A ) );
      vals.push_back( compute_valSet( nEvals,  g, Amax ) );
      vals.push_back( compute_valSet( nEvals,  g, B ) );
      vals.push_back( compute_valSet( nEvals,  g, Bmax ) );
      size_t maxVal = 0;
      g.logg << "Vals: ";
      for (size_t i = 0; i < vals.size(); ++i) {
	 g.logg << vals[i] << " ";
	 if (vals[i] > maxVal) {
	    maxVal = vals[i];
	 }
      }
      
      g.logg << endL;
      g.logg << "Evals: " << nEvals << endL;

      reportResults( nEvals, maxVal, valRounds.size() );

      if (reportRounds) {
	 for (size_t j = 0; j < valRounds.size(); ++j) {
	    //allResults.add( to_string( j ), valRounds[ j ] );
	 }
      }
   }
};



#endif
