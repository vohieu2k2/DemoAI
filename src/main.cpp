#include "mygraph.cpp"
#include "algs.cpp"

#include <iostream>
#include <string>
#include <unistd.h>
#include <chrono>

using namespace mygraph;
using namespace std;

void print_help() {
  cout << "Options: " << endl;
  cout << "-G <graph filename in binary format>" << endl
       << "-k <cardinality constraint>" << endl
       << "-A [run ATG]" << endl
       << "-L [run LATG]" << endl
       << "-I [run InterlaceGreedy]" << endl
       << "-F [run FastInterlaceGreedy]" << endl
       << "-T [run Gupta et al. (2010)]" << endl
       << "-R [run RandomizedGreedy]" << endl
       << "-Q [run FastRandomizedGreedy]" << endl
       << "-B [run Blits]" << endl
       << "-l [turn off stealing behavior of (Fast)InterlaceGreedy]" << endl
       << "-o <outputFileName>" << endl
       << "-N <repetitions>" << endl
       << "-e <epsilon (default 0.1)>" << endl
       << "-d <delta (default 0.1)>" << endl
       << "-v [verbose]>" << endl
       << "-f [fast mode (for ATG)]>" << endl;
}

void parseArgs( int argc, char** argv, Args& arg ) {
  int c;
  extern char *optarg;

  if (argc == 1) {
    print_help();
    exit( 2 );
  }

  string sarg;
  
  while ((c = getopt( argc, argv, ":G:k:IQTRlSBALFN:o:e:d:vf") ) != -1) {
    switch(c) {
    case 'f':
      arg.fast = true;
      break;
    case 'v':
      arg.g.logg.set_level( DEBUG );
      break;
    case 'e':
       sarg.assign( optarg );
       arg.epsi = stod( sarg );
       break;
    case 'd':
       sarg.assign( optarg );
       arg.delta = stod( sarg );
       break;
    case 'o':
       sarg.assign( optarg );
       arg.outputFileName = sarg;
       break;
    case 'G':
      //graph specification
      arg.graphFileName.assign( optarg );
      break;
    case 'k':
       sarg.assign( optarg );
       arg.k = stoi( sarg );
       break;
    case 'N':
       sarg.assign( optarg );
       arg.N = stoi( sarg );
       break;
    case 'l':
       arg.steal = false;
       break;
    case 'I':
       arg.alg = IG;
       break;
    case 'F':
       arg.alg = FIG;
       break;
    case 'Q':
       arg.alg = FRG;
       break;
    case 'S':
       arg.alg = SG;
       break;
    case 'B':
       arg.alg = BLITS;
       break;
    case 'T':
       arg.alg = TG;
       break;
    case 'R':
       arg.alg = RG;
       break;
    case 'A':
       arg.alg = ATG;
       break;
    case 'L':
       arg.alg = LATG;
       break;
    case '?':
      print_help();
      exit( 2 );
      break;
    }
  }
}

void readGraph( Args& args ) {
   args.logg( INFO, "Reading graph from file: " + args.graphFileName + "..." );
   args.g.read_bin( args.graphFileName );
   args.logg( INFO, "Input finished.");
   args.logg << "Nodes: " << args.g.n << ", edges: " << args.g.m << endL;
}

void runAlg( Args& args ) {
   size_t N = args.N;
   allResults.init( "obj" );
   allResults.init( "nEvals" );
   allResults.init( "k" );
   allResults.add( "k", args.k );
   
   for (size_t i = 0; i < N; ++i) {
      args.g.logg << "runAlg: Repetition = " << i << endL;
      clock_t t_start = clock();
      std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
      switch (args.alg) {
      case IG:
	 {
	    args.logg(INFO, "Starting InterlaceGreedy..." );
	    Ig ig( args );
	    ig.run();
	 }
	 break;
      case FRG:
	 {
	    args.logg(INFO, "Starting FastRandomizedGreedy..." );
	    Frg frg( args );
	    frg.run();
	 }
	 break;
      case FIG:
	 {
	    args.logg(INFO, "Starting FastInterlaceGreedy..." );
	    Fig fig( args );
	    fig.run();
	 }
	 break;
      case BLITS:
	 {
	    args.logg(INFO, "Starting Blits..." );
	    Blits blits( args );
	    blits.run();
	 }
	 break;
      case SG:
	 {
	    args.logg(INFO, "Starting StandardGreedy..." );
	    Sg sg( args );
	    sg.run();
	 }
	 break;

      case TG:
	 {
	    args.logg(INFO, "Starting TripleGreedy..." );
	    Tg tg( args );
	    tg.run();
	 }
	 break;

      case ATG:
	 {
	    args.logg(INFO, "Starting TripleGreedy..." );
	    Atg atg( args );
	    atg.run();
	 }
	 break;

      case LATG:
	 {
	    args.logg(INFO, "Starting LATG..." );
	    Latg latg( args );
	    latg.run();
	 }
	 break;

      case RG:
	 {
	    args.logg(INFO, "Starting RandomGreedy..." );
	    Rg rg( args );
	    rg.run();
	 }
	 break;

      }
   
      
      args.tElapsed = elapsedTime( t_start );
      std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
      size_t WallTimeMillis = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
      args.wallTime = WallTimeMillis / 1000.0;
   }

 
}


void outputResults( Args& args ) {
   if (args.outputFileName != "") {
      args.g.logg << "Writing output to file: " << args.outputFileName << endL;
      ofstream of( args.outputFileName.c_str(), ofstream::out | ofstream::app );
      allResults.print( of );
   } else {
      //allResults.print( cout );
   }
}

int main(int argc, char** argv) {
  Args args;
  parseArgs( argc, argv, args );
  readGraph( args );
  runAlg( args );
  outputResults( args );
}
