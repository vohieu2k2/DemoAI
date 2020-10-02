import matplotlib.pyplot as plt
import numpy as np
import sys

which=0; # 0 = solValue

scaleDiv=1;
Nodes=875713;
if (len(sys.argv) == 1):
    outFName='val-google.pdf';
    colObj = 5;
    colObjStd = 6;
    scaleDiv=1000000;
    normalize=True
else:
    normalize=False
    if (sys.argv[1][0] == 'q'):
        which = 1; #1 = queries
        outFName='query-google.pdf';
        colObj = 7;
        colObjStd = 8;
    else:
        which = 2; #2 = rounds
        outFName='rounds-google.pdf';
        colObj=9;
        colObjStd = 10;

colX = 4;

fnames = [ 'resultsGoogle-T.txt',
           'resultsGoogle-E.txt',
           'resultsGoogle-Q.txt',
           'resultsGoogle-M.txt',
           'resultsGoogle-A.txt',
           'resultsGoogle-L.txt' ];
algnames = [ 'Gupta et al. 2010',
             'Ene et al. 2020',
             'FRG (Buchbinder et al. 2015)',
             'ANM (Fahrbach et al. 2019',
             'AST (ours)',
             'ATG (ours)' ];

colors = [ 'm',
           'b',
           'r',
           'g',
           'c',
           'k'];

markers = [ 's',
            'v',
            'o',
            '<',
            '>',
            'd' ];

nalgs = 6;
X = [];
Obj = [];
ObjStd = [];


for i in range( 0, nalgs ):
    fname = fnames[ i ];
    print ("Reading from file", fname);
    with open(fname) as f:
        lines = f.readlines();
        XS = [float(line.split()[ colX ]) for line in lines if line[0] != '#']
        SObj = [float(line.split()[ colObj ])/scaleDiv for line in lines if line[0] != '#']
        SObjStd = [float(line.split()[ colObjStd ])/scaleDiv for line in lines if line[0] != '#']
        X.append( XS );
        Obj.append (SObj);
        ObjStd.append( SObjStd );

if normalize:
    for i in range( 1, nalgs ):
        for j in range( 0, len( Obj[ i ] ) ):
            if j < len( Obj[ 0 ] ):
                Obj[i][j] = Obj[i][j] / Obj[0][j];
                ObjStd[i][j] = ObjStd[i][j] / Obj[0][j];
            else:
                Obj[i].pop();
                X[i].pop();
                ObjStd[i].pop();
        
plt.gcf().clear();
plt.rcParams['pdf.fonttype'] = 42;
plt.rcParams.update({'font.size': 20});
plt.rcParams.update({'font.weight': "bold"});
plt.rcParams["axes.labelweight"] = "bold";
#plt.xscale('log');


plt.xlabel( 'k' );
#plt.ticklabel_format(axis='both', style='sci' );
if (which == 1):
    plt.ylabel( 'Total Queries' );
    plt.yscale('log');
else:
    if (which == 2):
        plt.ylabel( 'Adaptive Rounds' );
        plt.yscale('log');
    else:
        plt.ylabel( 'Value / IteratedGreedy' );


markSize=15

if normalize:
    algmin = 1;
else:
    algmin = 0;
for i in range(algmin,nalgs):
    if i % 2 == 0:
        mi=round(len(X[i]) / 2) - 1;
    else:
        mi=round(len(X[i]) / 3) - 1;
    plt.plot( X[i], Obj[i], '-', marker=markers[i],  label=algnames[i],ms = markSize,color = colors[i],markevery=mi);
    BObj = np.asarray( Obj[i] );
    BObjStd = np.asarray( ObjStd[i] );
    plt.fill_between( X[i], BObj - BObjStd, BObj + BObjStd,
                      alpha=0.5, edgecolor=colors[i], facecolor=colors[i]);

#plt.errorbar( X, Obj, yerr=BObjStd, fmt='-');


#plt.xlim( 0.05, 1.05 );
plt.gca().grid(which='major', axis='both', linestyle='--')

#plt.grid(color='grey', linestyle='--' );

#plt.legend(loc='best', numpoints=1,prop={'size':10},framealpha=1.0);
plt.savefig( outFName, bbox_inches='tight' );

#plt.savefig( 'WithLegend.png', bbox_inches='tight',dpi=1000 );
