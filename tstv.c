#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
//#include "norm.h"  // include normal variate
#define MAXN 101 //251        //max genome length
#define MAXSTEPS 500    //max number of steps in adaptive walk
#define MAXMUTES 251   //max number of mutations in DFE
#define MAXWALKS 10    //max number of adaptive walks 
#define MAXSPECS 10    //max number of different ts:tv values for output files
#define MAXREPL 101     // Kishore: Number of repitions for bias average calculations, must be greater than numrepl
#define MAXDFES 6 // Kishore: number of short walks, greater than ndfes
#define NREVERSE 101 // Kishore: number of steps in the 'short' walk, including the first
#define MAXNUMSHIFTS 10 // Kishore: maximum number of shifted DFEs
#define MAXFTS 10 // Kishore: maximum number of fts we test during the long walk. 
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)


// give sequence seq, make a random transition at position pos to form newseq
void transition(int seq[], int n, int newseq[]);
// give sequence seq, make a random transversion at position pos to form newseq
void transversion(int seq[], int n, int newseq[]);
// given some fixed AT->AT and GC->GC fraction "ffixed",
// and GC->AT fraction "fts"
// make a mutation to seq to form newseq
void mutation(int seq[],int newseq[],double ffixed, double fts);
//double * dfe(double currentBias, int ancDFE[MAXN], int neighs[MAXN], int wis[4][4][MAXN]);
double computeVariance(double arr[], int n);
double compute_autoc(double arr[MAXREPL], int N, int lag);


float ran1(long *idum);
float gamdev(int ia, long *idum);
float gasdev(long *idum);
void nrerror(char error_text[]);

int N=100;  // genome length, global variable

void main() 
{

int gc_flag = 0;  // if this flag = 0, run ts:tv case; if 1, run at:gc case
int nsteps;  // nsteps is number of accepted steps in adaptive walk
int maxwalklength = 1e3;  // maxwalklength is how many possible new mutations
                           // are generated and tested during the adaptive walk.
                           // These are only accepted if s>0 and if a random
                           // number is less than 2s (mimicking fixation).
			   // Thus "walklengths" are much longer than "steps".
double ftstep=.8;  // increment for fraction of transitions
double ffixed = .1; // fixed fraction of GC-GC or AT-AT mutations in AT:GC case
int nwalks = 1;//60;//500; // number of adaptive walks
//int nreverse = 200; // Kishore: number of steps to go for reversal of bias
double fts, ftsmu;  // Ts fraction in wt and mutant
double fwt[MAXWALKS];  // final wt fitness at the end of the walk, for each walk
double dfe[MAXWALKS][MAXMUTES];  // s values (DFE) for all mutants for each walk 
double dfe_ftsmu[MAXDFES][MAXMUTES][MAXREPL];  // Output the DFE at different moments in the short walk
//double dfe_shifts[MAXMUTES][MAXNUMSHIFTS][MAXDFES][MAXREPL]; // KISHORE: num mutations x number of shifted dfes x step in short walk x number of replicates
double dfe_shift_data[MAXNUMSHIFTS][6][MAXREPL]; // output matrix has data pts: currentBias, shiftedBias, sbarben,sbardet,fben,fdet
double dfe_shift_data_summary[MAXNUMSHIFTS][7]; // output matrix has data pts: shiftedBias, mean(sbarben), mean(sbardel), mean(fben), mean(fdel), var(fben)
double biasEvolution[NREVERSE][MAXREPL];  // Kishore: matrix to track all bias iterations
double fbenDFEs[MAXDFES][MAXREPL]; // Kishore: Record fben for each DFE across replicates
int nstepsall[MAXWALKS];  // number of steps taken in each walk
int i,j,k, ift, iwalk, ispec, itry, pos, n;  // temporary counters 
double oneoverNe = .000001;  // deleterious mutations can also fix with
  // probability exp(2s)*oneoverNe i.e. exp(2s)/Ne for Ne = 1e5 here
  // This happened so rarely results were indistinguishable whether this was
  // allowed or commented out.
double wis[4][4][MAXN];  // w_i values for each locus, given base at that locus
                        //  and base at the epistatically coupled locus
int nmutes = 250;   // number of fitnesses to use to calculate DFE
int neighs[MAXN];   // locus of epistatic neighbour for each locus in sequence
int anc[MAXN], tmpseq[MAXN], ancEvolve[MAXN], ancDFE[MAXN];  // ancestor sequence, temporary sequence
double w0, wseq[MAXSTEPS];  // initial fitness, fitness at each step in the long walk;
double wseqShort[NREVERSE]; // fitness at each step in the short walk
double winit;
char filename[100];  // output file name
FILE *fpout;   // output file pointer
FILE *fpout2; // KISHORE: file pointer for tstv bias calculation (biases as they evolve)
FILE *dfeBias; // Kishore: File pointer for DFE at start of 'short' walk
FILE *dfeShift; // Kishore: File pointer for s values at each shifted DFE
FILE *unifStat; // Kishore: File pointer to look at uniformly changing biases
FILE *unifAutoCor; // Kishore: File pointer to look at the autocorrelations
FILE *fbenFile; // Kishore: File pointer to look at the fben across short walk
double fben[MAXWALKS], fdel[MAXWALKS];  // fractions of ben and del mutns in wt DFE
double sbarben[MAXWALKS],sbardel[MAXWALKS]; // mean s_ben and s_del in wt DFE
double fbenmu[MAXWALKS][MAXSPECS], fdelmu[MAXWALKS][MAXSPECS];  // as above for mutants
double sbarbenmu[MAXWALKS][MAXSPECS],sbardelmu[MAXWALKS][MAXSPECS]; // ditto
double var[MAXREPL];
double autocor[NREVERSE];
double uniformFit[NREVERSE][5][MAXREPL];  // Record bias, wtmp, binary classifier if accepted step (1)
// note that for mutant variables we keep track for mutants of "MAXSPECS" possible
// ts fractions
double sumpos, sumneg;  // housekeeping, counting, arithmetic
int npos, nneg, imute, ilength, istep;
double sums[20], dfemu;
double wtmp, s;
double currentBias, oldBias; // Kishore: Double to represent the current tstv bias, and older bias
char filenamets[100];  // Kishore: output file name for ts bias evolution
char dfefilenamebias[100]; //Kishore: output file name for DFEs
long idum=(1); // Kishore: Seed for the Normal 
double sigma = 0.01;  // increase to minimic the uniform dist. 
int dfeNum; // counter for dfe number
int numreverse = 100; // MUST BE LESS than NREVERSE
int ishort;
int numrepl = 100; // MUST be less than MAXREPL
double shiftedbias;
double firstUnif = 100; // Number of uniform variates I look at
int noReverseFlag = 0;  // 0: Reversals of bias allowed, 1: no reversal allowed in short walk
int uniformFlag = 1; // 0: Evolve the bias with a normal distribution (better), 1: Evolve the bias with a uniform distribution 
int autocorFlag = 0; // 0: Do not compute autocorrelation, 1: Compute autocorrelation
int fbenFlag = 1; // 0: Do not look at fbens for DFEs in short walk, 1: DO look at fbens
double s_temp; // Kishore: Anonther s value placeholder
// Edit the line below (several examples shown) to set the walklengths required.
// Use only one walklengths and nlengths=1 for waterfall plots, for example.
// Note that this code is extremely inefficient but runs a completely new independent
// set of walks for each length.  So if you have walklengths = {1e4, 2e4} for example,
// the code will generate the required number of walks to length 1e4, then start over
// and generate NEW walks from time 0 to time 2e4, rather than continue the walks used
// for 1e4.  This could be easily changed.

int walklengths[1] = {5e5}, nlengths = 1;
// int walklengths[6] = {1e4, 2e4, 5e4, 1e5, 2e5, 5e5}, nlengths = 6;
// int walklengths[12] = {3e4, 4e4, 6e4, 7e4, 8e4, 9e4, 12e4, 14e4, 16e4, 18e4, 3e5, 4e5}, nlengths = 12;


// The line below is for making DFEs at specific points in the reversal of mutational bias experiment
// Each element in the array dfenumbers represents a different time point to sample and assess the DFE shifts

//int dfenumbers[5] = {0, 50, 100, 150, 175}, ndfes = 5;  
//int dfenumbers[5] = {0, 25, 50, 100, 150}, ndfes = 5;  
int dfenumbers[5] = {0, 10, 25, 50, 75}, ndfes = 5;  
//int dfenumbers[1] = {10}, ndfes = 1;  
//int dfenumbers[4] = {0, 2, 5, 7}, ndfes = 4;  

// The line below is for setting the transition biases that we shift to at certain steps in the short walk (steps above)
// nshifts = 1 if do not want to shift at all, and just want to test out various DFEs
double shftBias[5] = {0,0.1,0.3,0.5,0.7}, nshifts = 1; // make sure to change to 5 nshifts

srand(6);  // initialize random number generator

// Checks to make sure there is no memory overflow
if (N >= MAXN) fprintf(stdout, "Increase the size of MAXN\n");
if (numreverse >= NREVERSE) fprintf(stdout, "Please increase numreverse, it is >= NREVERSE right now\n");
if (numrepl >= MAXREPL) fprintf(stdout, "numrepl >= MAXREPL. Fix this \n");
if (nwalks >= MAXWALKS) fprintf(stdout,"iwalk is greater than MAXWALKS. Will cause memory overflows.\n");
if (nmutes >= MAXMUTES) fprintf(stdout, "nmutes >=  MAXMUTES");
if (firstUnif >= MAXREPL) fprintf(stdout, "firstUnif >= MAXREPL");

if (gc_flag)  fprintf(stdout,"NOTE: ts = bias GC -> AT (0,1), tv = bias AC->GC (2,3), ffixed = bias AT->AT or GC->GC fixed to %f\n",ffixed); 

for (ilength=0;ilength<nlengths;ilength++) {  // loop over walk lengths
   maxwalklength = walklengths[ilength];  // how many adaptive steps will be accepted

    // uncomment only one of the two lines below, depending on whether we want
    // to investigate a range of wt fts (waterfall plot), or a single fixed wt fts
  
// for (fts=ftstep;fts<=1-ffixed;fts+=ftstep) {
for (fts=0.9; fts<=0.9; fts+=ftstep) {  // increment by fraction of transitions
  
  // Better to initialize all file names here

  sprintf(filename, "data/N%d_Sw_walk%d_fts%1.2f.out",N,maxwalklength,fts);
  fprintf(stdout, filename); // "Filename"
  fpout = fopen(filename,"w");
  if (fpout == NULL) {   
    fprintf(stdout, "Couldn't open file");
  }

  // DFE SHIFTS

  // Uncomment if want to check the effect of shifting on DFE (be sure to uncomment area where writing to file as well)
  //sprintf(filename, "data/dfe_shift/N%d_Sw_walk%d_fts%1.2f_Step175_reversalIncl_shiftedDfeStatsv7.out",N,maxwalklength,fts);
  //fprintf(stdout, filename); // "Filename"
  //dfeShift = fopen(filename,"w");
  //if (dfeShift == NULL) {
  //  fprintf(stdout, "Couldn't open file");
  //}

  if (fbenFlag == 1){
    sprintf(filename, "data/tstv_dfe/fbenData/N%d_Sw_walk%d_fts%1.2f_uniformReversal_fben.out",N,maxwalklength,fts);
    fprintf(stdout, filename); // "Filename"
    fbenFile = fopen(filename,"w");
    if (fbenFile == NULL) {
      fprintf(stdout, "Couldn't open file UnifStat");
    }
  }

  if (uniformFlag == 1){
    sprintf(filename, "data/uniform/N%d_Sw_walk%d_fts%1.2f_uniformBiasv7.out",N,maxwalklength,fts);
    fprintf(stdout, filename); // "Filename"
    unifStat = fopen(filename,"w");
    if (unifStat == NULL) {
      fprintf(stdout, "Couldn't open file UnifStat");
    }
  }

  if (autocorFlag == 1){
    sprintf(filename, "data/uniform/autocor/N%d_Sw_walk%d_fts%1.2f_autocorUnifv2.out", N, maxwalklength, fts);
    fprintf(stdout, filename); // "Filename"
    unifAutoCor = fopen(filename,"w");
    if (unifAutoCor == NULL) {
      fprintf(stdout, "Couldn't open file unifAutoCor");
    }
  }

  // Open file to write out biases for this specific maxwalklength, ftsmu and N
 // sprintf(filenamets, "data/tstv_bias/N%d_Sw_walk%d_ftsmu%1.2f_v9.out",N,maxwalklength,fts); // change back to v6
 // fprintf(stdout, filenamets); // "Filename"
 // fpout2 = fopen(filenamets,"w");
 // if (fpout2 == NULL) {
 //   fprintf(stdout, "Couldn't open file"); 
 //}

  // file for DFE, changes for each n.
  // dfe file name, remove no_reversal, add reversalIncl
 // sprintf(dfefilenamebias, "data/tstv_dfe/no_replicates_included/dfe_N%d_Sw_walk%d_ftsmu%1.2f_noReversal_0_25_50_100_125_fben7.out",N,maxwalklength,fts);
 // dfeBias = fopen(dfefilenamebias,"w");
 // if (dfeBias == NULL) {
 //   fprintf(stdout, "Couldn't open file");
 // }

 // Kishore: For my purposes, I have been running code with nwalks = 1. 
 // I am not sure (have not tested) how my code will change if nwalks > 1
 // I think many arrays will be re-written though, so may want to keep that in mind.
 // To fix this, and if memory is not an issue, I would add another dimension to each array
 for (iwalk=0;iwalk<nwalks;iwalk++) {  // loop over walks

   // first, create a new landscape for this walk, by setting the w_i and neighbours
    for (i=0;i<=3;i++) for (j=0;j<=3;j++)
      for (k=0;k<N;k++) wis[i][j][k] = (rand()/(double)RAND_MAX);
   
    for (k=0;k<N;k++) neighs[k] = (int)(N*(rand()/(double)RAND_MAX)); // random epistatic neighbour
   
    // pick a new random ancestor sequence
    for (k=0;k<N;k++) anc[k] = (int)(4*(rand()/(double)RAND_MAX));
   
    // compute the fitness of the ancestor 
    w0 = 0;
    for (k=0;k<N;k++) 
      w0 += wis[anc[k]][anc[neighs[k]]][k];

   // step 0.  fitness at step 0 is w0, and nsteps = 0.
   wseq[0] = w0;
   nsteps = 0;

   // now try to make a new step in the walk, up to maxwalklength number of attempts
    for (itry = 1; itry<=maxwalklength; itry++) {

      if (nsteps >= MAXSTEPS) fprintf(stdout, "nsteps must be less than MAXSTEPS\n");

    // use either "mutation" or "transition/transversion" depending on which case.
    // Since the gc case is more complicated, we move it all to the subroutine.

    // sequence a real random number (current ts bias) 
    // every step, change transition matrix 
    // current bias stored and if rand < bias mutation rate then pick from normal at 0 with SD sigma 
    // current bias = current + r (sampled from normal)

     if (gc_flag) mutation(anc,tmpseq,ffixed,currentBias);
     else {
       pos = (int)(N*(rand()/(double)RAND_MAX));  // find position to mutate
       if ((rand()/(double)RAND_MAX)<fts)  // transition with probability fts
         transition(anc,pos,tmpseq);
       else
        transversion(anc,pos,tmpseq);
     }
     // compute fitness of the mutated sequence, tmpseq (if not fit, will be discarded so only temp)
     wtmp = 0;
     for (k=0;k<N;k++) 
       wtmp += wis[tmpseq[k]][tmpseq[neighs[k]]][k];
     // and then compute s, relative to the current wt fitness
     s = wtmp/wseq[nsteps]-1;
    if (s>0) {  // if beneficial mutation
      if ((rand()/(double)RAND_MAX)<2*s) {  // if it survives drift (prob 2s)
        nsteps = nsteps+1;  // we accept this new sequence as a step in the walk

        for (k=0;k<N;k++) anc[k] = tmpseq[k]; // new sequence replaces the ancestor
        wseq[nsteps] = wtmp;  // new fitness gets added to the list of fitnesses on this walk
      }
    }  else {  // if the mutation is deleterious or neutral, we might still accept it
    // with probability exp(2s)/Ne
      if ((rand()/(double)RAND_MAX)<oneoverNe*exp(2*s)) {
       nsteps = nsteps+1;
       for (k=0;k<N;k++) anc[k] = tmpseq[k];
       wseq[nsteps] = wtmp;
      }
    }
   }  // finish loop on itry, we are now finished this adaptive walk

   fwt[iwalk] = wseq[nsteps];  // save the final fitness of the wt at end of walk
   // now make the DFE for the wt
   npos=0; nneg = 0; sumpos=0; sumneg=0;
   for (imute=0;imute<nmutes;imute++) {  // loop over new mutations to add to wt DFE
      if (imute >= MAXMUTES) fprintf(stdout, "imute must never be larger or equal to MAXMUTES\n");
   // first, make a single substitution in the wt sequence, as described above
     if (gc_flag) mutation(anc,tmpseq,ffixed,fts);
     else {
       pos = (int)(N*(rand()/(double)RAND_MAX));
       if ((rand()/(double)RAND_MAX)<fts)
        transition(anc,pos,tmpseq);
       else
         transversion(anc,pos,tmpseq);
     }

     // then compute the fitness and then s for the mutated sequence, tmpseq
     wtmp = 0;
     for (k=0;k<N;k++) 
       wtmp += wis[tmpseq[k]][tmpseq[neighs[k]]][k];
     // add this s value to the dfe for the wt for this walk
     dfe[iwalk][imute] = wtmp/fwt[iwalk] - 1;
     // keep track of how many beneficial and deleterious mutations in this dfe
     if (dfe[iwalk][imute]>0) {npos++; sumpos+= dfe[iwalk][imute];}
     if (dfe[iwalk][imute]<0) {nneg++; sumneg+= dfe[iwalk][imute];}
   }

   nstepsall[iwalk] = nsteps;  // keep track of how many accepted steps in this walk
   fben[iwalk] = (double)npos/nmutes;  // fraction beneficial in wt DFE
   fdel[iwalk] = (double)nneg/nmutes;  // fraction deleterious in wt DFE
   if (npos>0) sbarben[iwalk] = sumpos/npos; // mean s of beneficial
   else sbarben[iwalk]=0;
   if (nneg>0) sbardel[iwalk] = sumneg/nneg;  // mean s of deleterious
   else sbardel[iwalk]=0;

    fprintf(stdout, "\nThe number of steps in this long walk is %d, at a ft of %1.2f. \n", nsteps, fts);

  //  srand(7);
  //  fprintf(stdout, "7");

    // Here on each iteration we will do numrepl trials. For confidence intervals and 
    // to minimize statistical errors
    for (int n = 0; n < numrepl; n++){

      // Reset Seed, not actually necessary but for experiments I added this
      srand(n + 1);

      // Start off with our first DFE, for each of the n replicates I make a new DFE at steps
      // specified in dfenumbers array
      // Counter for how many DFEs I have made
      dfeNum = 0;

      // Counts how far we have moved in the short walk
      // this is the number of accepted steps
      istep = 0;
      ishort = 0; // start of short walk, the is number of tries in short walk

      // The fitness at the beggining of the short walk is the fitness at the END of the long walk. 
      // Also, this initial fitness should be the same at the beggining of each replicate
      wseqShort[0] = wseq[nsteps];

      // Initialize our fraction of transitions as the currentBias
      // The goal is for this bias to randomly change
      currentBias = fts; 

      if (numreverse >= NREVERSE) fprintf(stdout, "Do not let istep >= NREVERSE\n");
             
      // Record bias into matrix as first row entry
      biasEvolution[istep][n] = currentBias;

      // create copy of ancestral sequence
      // ancEvolve is the sequence at this point in the short walk
      for (k=0;k<N;k++) ancEvolve[k] = anc[k];

      // if we are at a position in the short walk where we can make 
      // a DFE, then do that
      // Note this is not in the while loop later because I want
      // this DFE to be made before currentBias has evolved.
      if (istep == dfenumbers[dfeNum]){

        // We want to make a DFE at this bias.... and then depending 
        // Make more DFEs of this sequence at a SHIFTED transition bias. 

        // Make a ancestral sequence to make the DFE from
        // ancDFE is the ancestral sequence that we mutate to get the DFE
        for (int g = 0; g < N; g++) ancDFE[g] = ancEvolve[g];

        // What is the initial fitness of the sequence (at short walk = istep)
        winit = 0;
        for (k=0;k<N;k++) 
          winit += wis[ancDFE[k]][ancDFE[neighs[k]]][k];

        // now make the DFE for the istep time of the short walk
        // Note that this is actually different than the shift loop, because
        // I want more than JUST descriptive statistics
        npos = 0; nneg = 0; sumpos = 0; sumneg = 0;
        for (imute=0;imute<nmutes;imute++) {  // loop over new mutations to add to wt DFE
          //first, make a single substitution in the wt sequence, as described above
          if (gc_flag) mutation(ancDFE,tmpseq,ffixed,currentBias);
          else {
            pos = (int)(N*(rand()/(double)RAND_MAX));
            if ((rand()/(double)RAND_MAX)<currentBias)
              transition(ancDFE,pos,tmpseq);
            else
              transversion(ancDFE,pos,tmpseq);
          }

          // then compute the fitness and then s for the mutated sequence, tmpseq
          wtmp = 0;
          for (k=0;k<N;k++) 
            wtmp += wis[tmpseq[k]][tmpseq[neighs[k]]][k];

          if (dfeNum >= MAXDFES) fprintf(stdout, "Do not let dfenum >= MAXDFES\n");
          if (imute >= MAXMUTES) fprintf(stdout, "Do not let imute >= MAXMUTES\n");
          // add this s value to the dfe for the wt for this walk
          dfe_ftsmu[dfeNum][imute][n] = wtmp/winit - 1;

          // Added 2/25/2022
          s_temp = wtmp/winit - 1;

          if(s_temp > 0) {npos++; sumpos+=wtmp/winit - 1;}
          if(s_temp < 0) {nneg++; sumneg+=wtmp/winit - 1;}


          
        }

        //fprintf(stdout, "The current dfeNum is %d \n", dfeNum);
        fbenDFEs[dfeNum][n] = (double) npos/nmutes;
        // End Add 2/25/2022
        
      // Now we want to shift the bias and calculate the dfe again, shift = 0 already accounted for above
      // Update the shift array so that the first element is the bias at this point 
      shftBias[0] = currentBias;
      // Iterate across all shifts
      for(int ishift = 0; ishift < nshifts; ishift++){

          // Make a ancestral sequence to make the DFE from
          for (int g = 0; g < N; g++) ancDFE[g] = ancEvolve[g];

          // The bias to shift do and run the DFE from
          shiftedbias = shftBias[ishift];

          // What is the initial fitness of the sequence (at short walk = istep)
          winit = 0;
          for (k=0;k<N;k++) 
            winit += wis[ancDFE[k]][ancDFE[neighs[k]]][k];

          // now make the DFE for the istep time of the short walk
          // and collect descriptive stats.
           npos = 0; nneg = 0; sumpos = 0; sumneg = 0;
           for (imute=0;imute<nmutes;imute++) {  // loop over new mutations to add to wt DFE
            //first, make a single substitution in the wt sequence, as described above
            if (gc_flag) mutation(ancDFE,tmpseq,ffixed,shiftedbias);
            else {
              pos = (int)(N*(rand()/(double)RAND_MAX));
              if ((rand()/(double)RAND_MAX)<shiftedbias)
                transition(ancDFE,pos,tmpseq);
              else
                transversion(ancDFE,pos,tmpseq);
            }

            // then compute the fitness and then s for the mutated sequence, tmpseq
            wtmp = 0;
            for (k=0;k<N;k++) 
              wtmp += wis[tmpseq[k]][tmpseq[neighs[k]]][k];
              // add this s value to the dfe for the wt for this walk
            s_temp = wtmp/winit - 1;

          // Keep track of data
          if(s_temp > 0) {npos++; sumpos+=wtmp/winit - 1;}
          if(s_temp < 0) {nneg++; sumneg+=wtmp/winit - 1;}

          }

          // dfe_shift_data contains:
          // currentBias, shiftedBias, sbarben, sbardel, fben, fdel
          dfe_shift_data[ishift][0][n] = currentBias; dfe_shift_data[ishift][1][n] = shiftedbias; 
          dfe_shift_data[ishift][4][n] = (double) npos/nmutes; dfe_shift_data[ishift][5][n] = (double) nneg/nmutes;

          // Must handle seperate cases where npos = 0 or nneg = 0
          if (npos > 0) dfe_shift_data[ishift][2][n] = sumpos/npos;
          else dfe_shift_data[ishift][2][n] = 0;
          if (nneg > 0) dfe_shift_data[ishift][3][n] = sumneg/nneg; 
          else dfe_shift_data[ishift][3][n] = 0;

       } // end of loop on dfe shifts

        // We have just done the first DFEs so increment
        dfeNum++;
      }  // end of if (istep = dfenumbers[dfeNum])

      // Now we have made DFEs initially, time to mutate the bias
      if (uniformFlag == 1) {
        // Want to start before evolving bias
        uniformFit[istep][0][n] = currentBias;
        uniformFit[istep][1][n] = wseq[nsteps]; // fitness at this point
      }

      // loop over mutations from the evolved ancestral sequence
      while (istep < numreverse){

        // placehold the bias before we update it, so that we can revert
        // back later
        oldBias = currentBias;
        
        // Evolve the bias with a gaussian deviate
        if (uniformFlag == 0) currentBias = currentBias + sigma * gasdev(&idum);

        // Evolve the bias with a uniform deviation (extreme shifts more likely)
        else {
          // Shift bias with probability 0.1
          if ((rand() / (double) RAND_MAX) < 0.1) currentBias = ((rand()/(double)RAND_MAX));
        }

        // biases above 1 and below 0 make no sense, so adjust for this
        if (currentBias > 1) currentBias = 1;
        if (currentBias < 0) currentBias = 0;

        // This will essentially make it so that we do not evolve the bias
        if (noReverseFlag == 1) currentBias = fts;

        //fprintf(stdout, "The current bias is %1.5f\n", currentBias);

        // mutate with this new bias
        // on first iterate ancEvolve was made before the while loop and if we accept this bias
        // (leads to favorable s) then we will use this new sequence as ancEvolve
        if (gc_flag) mutation(ancEvolve,tmpseq,ffixed,currentBias);
        else {
          pos = (int)(N*(rand()/(double)RAND_MAX));  // find position to mutate
          if ((rand()/(double)RAND_MAX)<currentBias)  // transition with probability currentBias
            transition(ancEvolve,pos,tmpseq);
          else
            transversion(ancEvolve,pos,tmpseq);
        }

      // compute fitness of the mutated sequence, tmpseq (if not fit, will be discarded so only temp)
      wtmp = 0;
      for (k=0;k<N;k++) 
        wtmp += wis[tmpseq[k]][tmpseq[neighs[k]]][k];
      // and then compute s, relative to the current wt fitness
      s = wtmp/wseqShort[istep]-1;

      if (s>0) {  // if beneficial mutation
        // fitness-weighted adaptive walk 
        // previously this was 2*s
        if ((rand()/(double)RAND_MAX)< 2*s) {  // if it survives drift (prob 2s)
          wseqShort[istep + 1] = wtmp;  // new fitness gets added to the list of fitnesses on this walk
          // Record wtmp if uniform deviates and record pos in short walk
          if (uniformFlag == 1) {
            uniformFit[istep  + 1][0][n] = currentBias; // Bias at this point
            uniformFit[istep + 1][1][n] = wtmp;  // Enter fitness
            uniformFit[istep + 1][2][n] = ishort;
          }          

          for (k=0;k<N;k++) ancEvolve[k] = tmpseq[k]; // new sequence replaces the ancestor

          // make dfe at this step if it is a 'dfe' step
          // IF ISTEP is in dfenumbers array, MAKE DFE
          if (dfeNum < ndfes && istep == dfenumbers[dfeNum]){

            // We want to make a DFE at this bias.... and then depending on what the flag is
            // Make more DFEs of this sequence at a SHIFTED transition bias. 

            // Make an ancestral sequence to make the DFE from
            for (int g = 0; g < N; g++) ancDFE[g] = ancEvolve[g];

            // What is the initial fitness of the sequence (at short walk = istep)
            winit = 0;
            for (k=0;k<N;k++) 
              winit += wis[ancDFE[k]][ancDFE[neighs[k]]][k];  // changed from tmpseq

            // now make the DFE for the istep time of the short walk
            npos = 0; nneg = 0; sumpos = 0; sumneg = 0;
            for (imute=0;imute<nmutes;imute++) {  // loop over new mutations to add to wt DFE
              //first, make a single substitution in the wt sequence, as described above
              if (gc_flag) mutation(ancDFE,tmpseq,ffixed,currentBias);
              else {
                pos = (int)(N*(rand()/(double)RAND_MAX));
                if ((rand()/(double)RAND_MAX)<currentBias)
                  transition(ancDFE,pos,tmpseq);
                else
                  transversion(ancDFE,pos,tmpseq);
              }

              // then compute the fitness and then s for the mutated sequence, tmpseq
              wtmp = 0;
              for (k=0;k<N;k++) 
                wtmp += wis[tmpseq[k]][tmpseq[neighs[k]]][k];

              if (dfeNum >= MAXDFES) fprintf(stdout, "Do not let dfenum >= MAXDFES\n");
              if (imute >= MAXMUTES) fprintf(stdout, "Do not let imute >= MAXMUTES\n");

              // add this s value to the dfe for the wt for this walk
              dfe_ftsmu[dfeNum][imute][n] = wtmp/winit - 1;
              // Added 2/25/2022
              s_temp = wtmp/winit - 1;

              if(s_temp > 0) {npos++; sumpos+=wtmp/winit - 1;}
              if(s_temp < 0) {nneg++; sumneg+=wtmp/winit - 1;}
          }

         // fprintf(stdout, "The current dfeNum is %d", dfeNum);
          fbenDFEs[dfeNum][n] = (double) npos/nmutes;
          // End Add 2/25/2022

           // Now we want to shift the bias and calculate the dfe again, shift = 0 already accounted for above
           // First, find DFE with the current bias
           shftBias[0] = currentBias;
           for(int ishift = 0; ishift < nshifts; ishift++){
              // Make a ancestral sequence to make the DFE from
              for (int g = 0; g < N; g++) ancDFE[g] = ancEvolve[g];

              shiftedbias = shftBias[ishift];

              // What is the initial fitness of the sequence (at short walk = istep)
              winit = 0;
              for (k=0;k<N;k++) 
                winit += wis[ancDFE[k]][ancDFE[neighs[k]]][k];

              // now make the DFE for the istep time of the short walk
              npos = 0; nneg = 0; sumpos = 0; sumneg = 0;
              for (imute=0;imute<nmutes;imute++) {  // loop over new mutations to add to wt DFE
                //first, make a single substitution in the wt sequence, as described above
                if (gc_flag) mutation(ancDFE,tmpseq,ffixed,shiftedbias);
                else {
                  pos = (int)(N*(rand()/(double)RAND_MAX));
                  if ((rand()/(double)RAND_MAX)<shiftedbias)
                    transition(ancDFE,pos,tmpseq);
                  else
                    transversion(ancDFE,pos,tmpseq);
                }

                // then compute the fitness and then s for the mutated sequence, tmpseq
                wtmp = 0;
                for (k=0;k<N;k++) 
                  wtmp += wis[tmpseq[k]][tmpseq[neighs[k]]][k];

                // add this s value to the dfe for the wt for this walk
                s_temp = wtmp/winit - 1;

                // Keep track of data
                if(s_temp > 0) {npos++; sumpos+=wtmp/winit - 1;}
                if(s_temp < 0) {nneg++; sumneg+=wtmp/winit - 1;}

              }

              // The first row in our data matrix for shifts is
              // currentBias, shiftedBias sbarben, sbardel, fben, fdel
              dfe_shift_data[ishift][0][n] = currentBias; dfe_shift_data[ishift][1][n] = shiftedbias;
              dfe_shift_data[ishift][4][n] = (double) npos/nmutes; dfe_shift_data[ishift][5][n] = (double) nneg/nmutes;

              // Must handle seperate cases where npos = 0 or nneg = 0
              if (npos > 0) dfe_shift_data[ishift][2][n] = sumpos/npos;
              else dfe_shift_data[ishift][2][n] = 0;
              if (nneg > 0) dfe_shift_data[ishift][3][n] = sumneg/nneg; 
              else dfe_shift_data[ishift][3][n] = 0;

           }

            dfeNum++;
          }

          istep = istep+1;  // we accept this new sequence as a step in the walk
          // save bias
          if (istep > MAXSTEPS) {
            fprintf(stdout, "ERROR, MAXSTEPS is too small");
          }
          biasEvolution[istep][n] = currentBias;

        } 
        else {
          // Added recently
          currentBias = oldBias;
       }
      }  else {  // if the mutation is deleterious or neutral, we might still accept it
      // with probability exp(2s)/Ne
        if ((rand()/(double)RAND_MAX)<oneoverNe*exp(2*s)) {
          wseqShort[istep + 1] = wtmp;  // new fitness gets added to the list of fitnesses on this walk

          // Input this as an accepted step if looking at uniform dist
          if (uniformFlag == 1) {
            uniformFit[istep  + 1][0][n] = currentBias; // Bias at this point
            uniformFit[istep + 1][1][n] = wtmp;  // Enter fitness
            uniformFit[istep + 1][2][n] = ishort;
          }  
          for (k=0;k<N;k++) ancEvolve[k] = tmpseq[k];

          // Make dfe if needed
          // IF ISTEP is in dfenumbers array, MAKE DFE
          // also must make sure that the dfe number we are currently making is always less than ndfess
          if (dfeNum < ndfes && istep == dfenumbers[dfeNum]){
            // MAKE DFE

            // We want to make a DFE at this bias.... and then depending on what the flag is
            // Make more DFEs of this sequence at a SHIFTED transition bias. 

            // Make a ancestral sequence to make the DFE from
            for (int g = 0; g < N; g++) ancDFE[g] = ancEvolve[g];

            // What is the initial fitness of the sequence (at short walk = istep)
            winit = 0;
            for (k=0;k<N;k++) 
              winit += wis[ancDFE[k]][ancDFE[neighs[k]]][k];

            // now make the DFE for the istep time of the short walk
            npos = 0; nneg = 0; sumpos = 0; sumneg = 0;
            for (imute=0;imute<nmutes;imute++) {  // loop over new mutations to add to wt DFE
              //first, make a single substitution in the wt sequence, as described above
              if (gc_flag) mutation(ancDFE,tmpseq,ffixed,currentBias);
              else {
                pos = (int)(N*(rand()/(double)RAND_MAX));
                if ((rand()/(double)RAND_MAX)<currentBias)
                  transition(ancDFE,pos,tmpseq);
                else
                  transversion(ancDFE,pos,tmpseq);
              }

              // then compute the fitness and then s for the mutated sequence, tmpseq
              wtmp = 0;
              for (k=0;k<N;k++) 
                wtmp += wis[tmpseq[k]][tmpseq[neighs[k]]][k];

              if (dfeNum >= MAXDFES) fprintf(stdout, "Do not let dfenum >= MAXDFES");
              if (imute >= MAXMUTES) fprintf(stdout, "Do not let imute >= MAXMUTES");

              // add this s value to the dfe for the wt for this walk
              dfe_ftsmu[dfeNum][imute][n] = wtmp/winit - 1;

              // Added 2/25/2022
              s_temp = wtmp/winit - 1;

              if(s_temp > 0) {npos++; sumpos+=wtmp/winit - 1;}
              if(s_temp < 0) {nneg++; sumneg+=wtmp/winit - 1;}
            }
            //fprintf(stdout, "The current dfeNum is %d", dfeNum);

            fbenDFEs[dfeNum][n] = (double) npos/nmutes;
            // End 2/25/2022

            // Now we want to shift the bias and calculate the dfe again, shift = 0 already accounted for above
            shftBias[0] = currentBias;
            for(int ishift = 0; ishift < nshifts; ishift++){
              // Make a ancestral sequence to make the DFE from
              for (int g = 0; g < N; g++) ancDFE[g] = ancEvolve[g];

              shiftedbias = shftBias[ishift];

              // What is the initial fitness of the sequence (at short walk = istep)
              winit = 0;
              for (k=0;k<N;k++) 
                winit += wis[ancDFE[k]][ancDFE[neighs[k]]][k];

              // now make the DFE for the istep time of the short walk
              npos = 0; nneg = 0; sumpos = 0; sumneg = 0;
              for (imute=0;imute<nmutes;imute++) {  // loop over new mutations to add to wt DFE
                //first, make a single substitution in the wt sequence, as described above
                if (gc_flag) mutation(ancDFE,tmpseq,ffixed,shiftedbias);
                else {
                  pos = (int)(N*(rand()/(double)RAND_MAX));
                  if ((rand()/(double)RAND_MAX)<shiftedbias)
                    transition(ancDFE,pos,tmpseq);
                  else
                    transversion(ancDFE,pos,tmpseq);
                }

                // then compute the fitness and then s for the mutated sequence, tmpseq
                wtmp = 0;
                for (k=0;k<N;k++) 
                  wtmp += wis[tmpseq[k]][tmpseq[neighs[k]]][k];
                // add this s value to the dfe for the wt for this walk
                s_temp = wtmp/winit - 1;

                // Keep track of data
                if(s_temp > 0) {npos++; sumpos+=wtmp/winit - 1;}
                if(s_temp < 0) {nneg++; sumneg+=wtmp/winit - 1;}

              }

              // The first row in our data matrix for shifts is
              // currentBias, shiftedBias, sbarben, sbardel, fben, fdel
              dfe_shift_data[ishift][0][n] = currentBias; dfe_shift_data[ishift][1][n] = shiftedbias; 
              dfe_shift_data[ishift][4][n] = (double) npos/nmutes; dfe_shift_data[ishift][5][n] = (double) nneg/nmutes;

              // Must handle seperate cases where npos = 0 or nneg = 0
              if (npos > 0) dfe_shift_data[ishift][2][n] = sumpos/npos;
              else dfe_shift_data[ishift][2][n] = 0;
              if (nneg > 0) dfe_shift_data[ishift][3][n] = sumneg/nneg; 
              else dfe_shift_data[ishift][3][n] = 0;

            }
            
            dfeNum++;
          }
          // istep only increases if the mutation is favorable!
          istep = istep + 1;   
          // save bias
          biasEvolution[istep][n] = currentBias; 
       //   wseqShort[istep] = wtmp;
        }

        else {
          // If we have a unfit sequence, then we go back to the old bias
          currentBias = oldBias;
          // Do NOT save this fitness (wtmp) value. This is because it was not an accepted step.
          }
        } // end of else (s <= 0) clause

      ishort++; // have moved another step in short walk
      } // End of iteration over short walk



    }

    // KISHORE OUTPUT 

    if (uniformFlag == 1){

      // Output Results
      for (int irev = 0; irev < numreverse; irev++){
        for(int istat = 0; istat < 3; istat++){

          int irep = 0;
          double avgShift = 0;  // Record average shift across replicates (first few amount)

          // comment out this line if want the three-pronged outputs (bias, fitness, step in short walk)
          if (istat == 0){

          while (irep < numrepl){
            avgShift += uniformFit[irev][istat][irep];
            var[irep] = uniformFit[irev][istat][irep]; // to compute variance
            irep++;
            // uncomment below if want three-prong output
            // fprintf(unifStat, "%f", uniformFit[irev][istat][irep]) to write out all data
          }

          avgShift = avgShift / irep;
          fprintf(stdout, "The average shift: %f\n", avgShift);

        fprintf(unifStat, "%f,%f,", avgShift, computeVariance(var, numrepl));
        }
        }

      fprintf(unifStat, "\n");
      }

      fclose(unifStat);
    }

    if (autocorFlag == 1){
      int lag = 1;

      for(int irep = 0; irep < numrepl; irep++){
        for (int istat = 0; istat < 3; istat++){
          for(int irev = 0; irev < numreverse; irev++){
              if (istat == 0){
                // Look at the autocorrelations of just the bias
                autocor[irev] = uniformFit[irev][istat][irep];
              }
          }

        }
        fprintf(unifAutoCor, "%f\n", compute_autoc(autocor,numreverse,lag));

      }

      fclose(unifAutoCor);
    }



  if (fbenFlag == 1){
    for (int i = 0; i < numrepl; i++){
      for (int idfe = 0; idfe < ndfes; idfe++){
        fprintf(fbenFile, "%f,", fbenDFEs[idfe][i]);
      }
      fprintf(fbenFile, "\n");
    }

    fclose(fbenFile);
  }

    // output dfe results during short walk 
  //  for (int irep = 0; irep < numrepl; irep++){
  //    for (int imute = 0; imute < nmutes; imute++){
  //        for (int idfe = 0; idfe < ndfes; idfe++){
  //          fprintf(dfeBias, "%f,", dfe_ftsmu[idfe][imute][irep]);
  //        }
  //        fprintf(dfeBias, "\n");
  //    }
  //  }

   //fclose(dfeBias);

    // output results
 //   for (int i = 0; i < numreverse; i++){
 //     for (int j = 0; j < numrepl; j++){
 //       fprintf(fpout2, "%f,", biasEvolution[i][j]);
 //     }
 //     fprintf(fpout2, "\n");
 //   }

 //  fclose(fpout2);

  // Output summarized DFE statistics for shifts
  // Skip the first shift (this is just the current Bias and NOT actually shifted)
  for (int ishift = 1; ishift < nshifts; ishift++) {
    // Ignore the 'current bias' statistic, since this will shift on each replicate
    for (int istat = 1; istat <= 5; istat++) {
      double statbar = 0; int repl = 0;
      for(int n = 0; n < numrepl; n++) {
          // Compute average statistics
          // dfe_shift_data has: currentBias, shiftedBias, sbarben, sbardel, fben, fdel 
          repl++; statbar+=dfe_shift_data[ishift][istat][n];
          if (istat == 4) {     // The fourth statistic (really the 5th since 0 inclusive) is fben
            var[n] = dfe_shift_data[ishift][istat][n];
          }
      }
      // Put average statistics into array
      // mean(shiftedBias), mean(sbarben), mean(sbardel), mean(fben), mean(fdel)
      dfe_shift_data_summary[ishift - 1][istat - 1] = (double) statbar / repl;
    }
    // mean(shiftedBias), mean(sbarben), mean(sbardel), mean(fben), mean(fdel), var(fben)
    dfe_shift_data_summary[ishift - 1][5] = computeVariance(var, numrepl);
  }
// output results as the form:
// Each row represents a different shifted bias with data points: shiftedBias, mean(sbarben), mean(sbardel), mean(fben), mean(fdel), var(fben) across
// all n replicates
  //for (int i = 0; i < nshifts - 1; i++){
  //  for (int j = 0; j < 6; j++){
  //    fprintf(dfeShift, "%f,", dfe_shift_data_summary[i][j]);
  //  }
  //  fprintf(dfeShift, "\n");
 // }

  //fclose(dfeShift);
  
// END KISHORE OUTPUT

  // now compute dfe for a mutator with a different ts fraction from the wt
  // The mutator has the same sequence as the wt, and differs only in ts:tv bias
  // (or gc:at bias).
  ispec = 0;   

  // loop over ts fraction of mutator, edit the line below for gc:at case as needed
  //for (ftsmu = 0.1; ftsmu <= 0.9; ftsmu += 0.8) {   // Kishore: Just do 0.1, 0.9
  for (ftsmu=ftstep;ftsmu<=1-ftstep;ftsmu+=ftstep) { // Original Case in code, 0.1, 0.2, 0.3...0.9
  // for (ftsmu=0.05; ftsmu<=1-ffixed; ftsmu+=ftstep) {  // example for gc:at case

    // make a dfe for this mutator, exactly as described above for the wt
    // but only for ts = 0.9
    npos=0; nneg = 0; sumpos=0; sumneg=0;
    for (imute=0;imute<nmutes;imute++) { 

      if (gc_flag) mutation(anc,tmpseq,ffixed,ftsmu);
      else {
        pos = (int)(N*(rand()/(double)RAND_MAX));
        if ((rand()/(double)RAND_MAX)<currentBias)
          transition(anc,pos,tmpseq);
        else {       
          transversion(anc,pos,tmpseq);
        }
      }
     wtmp = 0;
     for (k=0;k<N;k++) 
       wtmp += wis[tmpseq[k]][tmpseq[neighs[k]]][k];
     dfemu = wtmp/fwt[iwalk] - 1;
     if (dfemu>0) {npos++; sumpos+= dfemu;}
     if (dfemu<0) {nneg++; sumneg+= dfemu;}

    


    }  // end of loop on mutations


   fbenmu[iwalk][ispec] = (double)npos/nmutes;
   fdelmu[iwalk][ispec] = (double)nneg/nmutes;
   if (npos>0) sbarbenmu[iwalk][ispec] = sumpos/npos;
   else sbarbenmu[iwalk][ispec]=0;
   if (nneg>0) sbardelmu[iwalk][ispec] = sumneg/nneg;
   else sbardelmu[iwalk][ispec] = 0;
   ispec++;
   }  // end of loop on ispec (stepping through mutator ts values)

   
   // at this point, we have completed an adaptive walk across this landscape
   // We have created a dfe for the wt
   // We have created a dfe for a set of mutators that have the same sequence
   // as the wt, but differ in their ts fraction.
   // Now we repeat this be creating a new landscape, new ancestor, and walking again.
}  // end of loop on iwalks

// output the shifted data for DFEs

// keep track of the average fraction beneficial, etc, across all walks
sums[1] = 0; sums[2] = 0; sums[3] = 0; sums[4] = 0; sums[5]=0; sums[6]=0;
int nzb=0; int nzd=0;
for (iwalk=0;iwalk<nwalks;iwalk++) {
  sums[1]+= fben[iwalk]; sums[2]+= sbarben[iwalk];
  sums[3]+= fdel[iwalk]; sums[4]+= sbardel[iwalk];
  sums[5]+= fwt[iwalk];  sums[6]+= nstepsall[iwalk];
  if (sbarben[iwalk]>0) nzb++;  // number of beneficial mean s values that are non-zero
  if (sbardel[iwalk]<0) nzd++; // ditto, deleterious
}
// print results to output file, first for wt and then for mutants
// note first row of output file contains: fts, ftsstep, fwt, nsteps
// following rows contain fben, sbarben, fdel, sbardel for wt and then all mutators
// put your output filename here, in this example "Sw" means s-weighted walk

// PUT FILE OPENING HERE
fprintf(fpout, "%f %f %f %f\n", fts, ftstep, sums[5]/nwalks, sums[6]/nwalks); 
fprintf(fpout,"%f %f %f %f\n",sums[1]/nwalks,sums[2]/nzb,sums[3]/nwalks,sums[4]/nzd);

for (i=0;i<ispec;i++) {
sums[1] =0; sums[2]=0; sums[3]= 0; sums[4] = 0; nzb=0; nzd=0; 
for (iwalk=0;iwalk<nwalks;iwalk++) {
  sums[1]+= fbenmu[iwalk][i]; sums[2]+= sbarbenmu[iwalk][i];
  sums[3]+= fdelmu[iwalk][i]; sums[4]+= sbardelmu[iwalk][i];
  if (sbarbenmu[iwalk][i]>0) nzb++;
  if (sbardelmu[iwalk][i]<0) nzd++;
}
fprintf(fpout,"%f %f %f %f\n",sums[1]/nwalks,sums[2]/nzb,sums[3]/nwalks,sums[4]/nzd);
}

fclose(fpout);


// uncomment below for more detailed output as required

 //fprintf(stdout,"Note: printing all data to output files\n");
 //fpout = fopen("wt_data.txt","w");
 //for (iwalk=0;iwalk<nwalks;iwalk++)
 //  fprintf(fpout,"%f %f %f %f %f %d\n",fben[iwalk],sbarben[iwalk],fdel[iwalk],sbardel[iwalk],fwt[iwalk],nstepsall[iwalk]);
 //fclose(fpout);
 //fpout = fopen("mt_data.txt","w");
 //for  (iwalk=0;iwalk<nwalks;iwalk++)
 //  for (i=0;i<ispec;i++)
 //    fprintf(fpout," %f %f %f %f\n",fbenmu[iwalk][i],sbarbenmu[iwalk][i],fdelmu[iwalk][i],sbardelmu[iwalk][i]);
 //fclose(fpout);



}  // loop on fts
 }  // loop on walklength



}


double computeVariance(double arr[MAXREPL], int n)
{
  int count = 0; double sum = 0; double sumsq = 0;
  for (int i = 0; i < n; i++){
    count++; sum+=arr[i]; sumsq+=arr[i]*arr[i];
  }
  double var = (sumsq - sum*sum / (double) count) / (double) (count - 1);
  return var;

}

// Function to compute autocorrelations
double compute_autoc(double arr[MAXREPL], int N, int lag)
{
  double   autocv;      // Autocovariance value
  double   ac_value;    // Computed autocorrelation value to be returned
  int      i, j;        // Loop counter
  double mean = 0; 

  // Loop to compute mean
  for (j = 0; j < N; j++){
    mean += arr[j];
  }

  // Average by division
  mean = mean / N;

  // Loop to compute autocovariance
  autocv = 0.0;
  for (i=0; i<(N - lag); i++)
    autocv = autocv + ((arr[i] - mean) * (arr[i+lag] - mean));
  autocv = (1.0 / (N - lag)) * autocv;

  // Autocorrelation is autocovariance divided by variance
  ac_value = autocv / computeVariance(arr, N);

  return ac_value;
}


void transition(int seq[MAXN],int pos,int newseq[MAXN])
// give sequence seq, make a random transition at position pos to form newseq
{
for (int i=0;i<N;i++) newseq[i]=seq[i]; // first make a copy of seq into newseq
// then make a transition at position pos 
switch(seq[pos]) {
  case 0  : newseq[pos] = 1; break;
  case 1  : newseq[pos] = 0; break;
  case 2  : newseq[pos] = 3; break;
  case 3  : newseq[pos] = 2; break;
}
}

void transversion(int seq[MAXN],int pos,int newseq[MAXN])
// give sequence seq, make a random transversion at position pos to form newseq
{
for (int i=0;i<N;i++) newseq[i]=seq[i];  // first copy the seq to newseq
// then make a transversion at position pos (depends on base at position pos)
switch(seq[pos]) {
  case 0  :
  if ((rand()/(double)RAND_MAX)<0.5) newseq[pos] = 2; 
  else newseq[pos] = 3;
   break;
  case 1  :
  if ((rand()/(double)RAND_MAX)<0.5) newseq[pos] = 2;
  else newseq[pos] = 3;
  break;
  case 2  :
  if ((rand()/(double)RAND_MAX)<0.5) newseq[pos] = 0;
  else newseq[pos] = 1;
  break;
  case 3  :
  if ((rand()/(double)RAND_MAX)<0.5) newseq[pos] = 0;
  else newseq[pos] = 1;
  break;
}
}


void mutation(int seq[MAXN],int newseq[MAXN],double ffixed, double fts)
// given some fixed AT->AT and GC->GC fraction "ffixed",
// and GC->AT fraction "fts"
// make a mutation to seq to form newseq
{
  int i, pos;
  double muprob;
  
for (i=0;i<N;i++) newseq[i]=seq[i];  // first copy seq to newseq
 muprob = rand()/(double)RAND_MAX;
 if (muprob<ffixed) {  // we want AT to AT or GC to GC
  pos = (int)(N*(rand()/(double)RAND_MAX));
  switch(seq[pos]) {
    case 0: newseq[pos] = 1; break;
    case 1: newseq[pos] = 0; break;
    case 2: newseq[pos] = 3; break;
    case 3: newseq[pos] = 2; break;
   }
  }
 else {
  if (muprob<(ffixed+fts)) {   // we want GC to AT
   pos = (int)(N*(rand()/(double)RAND_MAX));
   while (seq[pos]<2)  pos = (int)(N*(rand()/(double)RAND_MAX));
   if ((rand()/(double)RAND_MAX)<0.5) newseq[pos] = 0;
   else newseq[pos] = 1;
  }
  else {// we want AT to GC
   pos = (int)(N*(rand()/(double)RAND_MAX));
   while (seq[pos]>1)  pos = (int)(N*(rand()/(double)RAND_MAX));
   if ((rand()/(double)RAND_MAX)<0.5) newseq[pos] = 2;
   else newseq[pos] = 3;
  }
 }
}


float gasdev(long *idum)
{
	float ran1(long *idum);
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;

	if (*idum < 0) iset=0;
	if  (iset == 0) {
		do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}
void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}


float gamdev(int ia, long *idum)
{
	float ran1(long *idum);
	void nrerror(char error_text[]);
	int j;
	float am,e,s,v1,v2,x,y;

	if (ia < 1) nrerror("Error in routine gamdev");
	if (ia < 6) {
		x=1.0;
		for (j=1;j<=ia;j++) x *= ran1(idum);
		x = -log(x);
	} else {
		do {
			do {
				do {
					v1=ran1(idum);
					v2=2.0*ran1(idum)-1.0;
				} while (v1*v1+v2*v2 > 1.0);
				y=v2/v1;
				am=ia-1;
				s=sqrt(2.0*am+1.0);
				x=s*y+am;
			} while (x <= 0.0);
			e=(1.0+y*y)*exp(am*log(x/am)-s*y);
		} while (ran1(idum) > e);
	}
	return x;
}
float ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

