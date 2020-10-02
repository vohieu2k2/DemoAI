import matplotlib.pyplot as plt
import numpy as np
import sys

which=0; # 0 = solValue

if (len(sys.argv) == 1):
    outFName='val-ba.pdf';
    colObj = 5;
    colObjStd = 6;
    normalize=True;
else:
    normalize=False;
    if (sys.argv[1][0] == 'q'):
        which = 1; #1 = queries
        outFName='query-ba.pdf';
        colObj = 7;
        colObjStd = 8;
    else:
        which = 2; #2 = rounds
        outFName='rounds-ba.pdf';
        colObj=9;
        colObjStd = 10;

colX = 4;

Enefnames = [ #'resultsBA-E-10000.txt',
    #'resultsBA-E-1000.txt',
              'resultsBA-E-100.txt',
              'resultsBA-E-10.txt',
              'resultsBA-E-1.txt' ];

EneSamples =[ #10000,
    100, 10, 1 ];

algnames = [ 'Gupta et al.',
             'FIG',
             'FRG',
             'Blits',
             'ANM',
             'ATG',
             'LATG' ];

colors = [ 'm',
           'b',
           'r',
           'y',
           'g',
           'c',
           'k'];

markers = [ 's',
            'v',
            'o',
            'p',
            '<',
            '>',
            'd' ];

nalgs = 7;
X = [];
Obj = [];
ObjStd = [];

scaleDiv=1;

# Read ENE results

print("Reading results for Ene et al. (2020):");
for i in range( 0, len(Enefnames ) ):
    fname = Enefnames[ i ];
    print ("Reading from file", fname);
    with open(fname) as f:
        lines = f.readlines();
        SObj = [float(line.split()[ colObj ])/scaleDiv for line in lines if line[0] != '#']
        SObjStd = [float(line.split()[ colObjStd ])/scaleDiv for line in lines if line[0] != '#' ]
        SObj = SObj[0];
        SObjStd = SObjStd[0];
        Obj.append (SObj);
        ObjStd.append( SObjStd );

X = EneSamples;
print("X: ", X)
print("Obj: ", Obj)


print("Reading results for Gupta et al. (2010):");
fname="resultsBA-T.txt"
with open(fname) as f:
    lines = f.readlines();
    GObj = [float(line.split()[ colObj ])/scaleDiv for line in lines if line[0] != '#']
    GObj = GObj[0];

if which == 0:
    scaleDiv = GObj;
else:
    scaleDiv = 1;

print("Reading results for Ene & Nguyen (2020) [exact]:");
fname="resultsBA-E-0.txt"
with open(fname) as f:
    lines = f.readlines();
    E2Obj = [float(line.split()[ colObj ])/scaleDiv for line in lines if line[0] != '#']

    E2Obj = [E2Obj[0] for x in X];
    E2ObjStd = [float(line.split()[ colObjStd ])/scaleDiv for line in lines if line[0] != '#']
    E2ObjStd = [E2ObjStd[0] for x in X];
    
print("Reading results for ATG:");
fname="resultsBA-L.txt"
with open(fname) as f:
    lines = f.readlines();
    LObj = [float(line.split()[ colObj ])/scaleDiv for line in lines if line[0] != '#']

    LObj = [LObj[0] for x in X];
    LObjStd = [float(line.split()[ colObjStd ])/scaleDiv for line in lines if line[0] != '#']
    LObjStd = [LObjStd[0] for x in X];

print("Reading results for AST:");
fname="resultsBA-A.txt"
with open(fname) as f:
    lines = f.readlines();
    AObj = [float(line.split()[ colObj ])/scaleDiv for line in lines if line[0] != '#']
    AObjStd = [float(line.split()[ colObjStd ])/scaleDiv for line in lines if line[0] != '#']
    AObj = [AObj[0]  for x in X];
    AObjStd = [AObjStd[0] for x in X];


if normalize:
    for j in range( 0, len( Obj ) ):
        Obj[j] = Obj[j] / GObj;
        ObjStd[j] = ObjStd[j] / GObj;

        
plt.gcf().clear();
plt.rcParams['pdf.fonttype'] = 42;
plt.rcParams.update({'font.size': 20});
plt.rcParams.update({'font.weight': "bold"});
plt.rcParams["axes.labelweight"] = "bold";
#plt.xscale('log');


plt.xlabel( 'Number of Samples' );
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

plt.xscale('log');
markSize=15

if normalize:
    algmin = 1;
else:
    algmin = 0;

mi=round(len(X) / 2) - 1;
plt.plot( X, Obj, '-', marker='^',  label='Ene & Nguyên 2020 [sampling]',ms = markSize,color = 'r',markevery=mi);
BObj = np.asarray( Obj );
BObjStd = np.asarray( ObjStd );
plt.fill_between( X, BObj - BObjStd, BObj + BObjStd, alpha=0.5, edgecolor='r', facecolor='r');
plt.xticks( X );

if (which >= 0):
    plt.plot( X, E2Obj, '-', marker='v',  label='Ene & Nguyên 2020 [exact]',ms = markSize,color = 'b',markevery=mi);
    BObj = np.asarray( E2Obj );
    BObjStd = np.asarray( E2ObjStd );
    plt.fill_between( X, BObj - BObjStd, BObj + BObjStd, alpha=0.5, edgecolor='b', facecolor='b');
    plt.xticks( X );

plt.plot( X, AObj, '-', marker='>',  label='AdaptiveSimpleThreshold',ms = markSize,color = 'c',markevery=mi);
BObj = np.asarray( AObj );
BObjStd = np.asarray( AObjStd );
plt.fill_between( X, BObj - BObjStd, BObj + BObjStd, alpha=0.5, edgecolor='c', facecolor='c');
plt.xticks( X );    

plt.plot( X, LObj, '-', marker='d',  label='AdaptiveThresholdGreedy',ms = markSize,color = 'k',markevery=mi);
BObj = np.asarray( LObj );
BObjStd = np.asarray( LObjStd );
plt.fill_between( X, BObj - BObjStd, BObj + BObjStd, alpha=0.5, edgecolor='r', facecolor='r');
plt.xticks( X );

# for i in range(algmin,nalgs):
#     if i % 2 == 0:
#         mi=round(len(X[i]) / 2) - 1;
#     else:
#         mi=round(len(X[i]) / 3) - 1;
#     plt.plot( X[i], Obj[i], '-', marker=markers[i],  label=algnames[i],ms = markSize,color = colors[i],markevery=mi);
#     BObj = np.asarray( Obj[i] );
#     BObjStd = np.asarray( ObjStd[i] );
#     plt.fill_between( X[i], BObj - BObjStd, BObj + BObjStd,
#                       alpha=0.5, edgecolor=colors[i], facecolor=colors[i]);

#plt.errorbar( X, Obj, yerr=BObjStd, fmt='-');

if (which==0):
    plt.ylim(0.0,1.05);


#plt.xlim( 0.05, 1.05 );
plt.gca().grid(which='major', axis='both', linestyle='--')

#plt.grid(color='grey', linestyle='--' );
plt.legend(loc='best', numpoints=1,prop={'size':15});
plt.savefig( outFName, bbox_inches='tight' );

