/* 
 * simulation of sequence length of sequence-by-synthesis (SBS) 
 * with flow cycle number f fixed, 
 * and comparison between the simulation results and theoretical results
 * 
 * Yong Kong
 * Yale University
 *
 * $Id: sbs_ffcm_simul.c,v 1.4 2011/08/18 14:03:20 yk336 Exp $
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <sys/time.h>
#include <ctype.h>
#include <time.h>
#include <dirent.h>

#include <assert.h>


#define MAX_RJ      20            /* max of rj */
#define MAX_C       256           /* upper limit of incorporation cycle */
#define MAX_N       1024          /* max of sequence length */


long long baseCnt[4];             /* base cnt */
double nSum = 0.0, nSum2 = 0.0;   /* seq length */

double ndist[MAX_N];              /* length counts */
double ndist4[MAX_N][4];          /* length counts of the 4 nt flows */


void usage(const char *name) {
  fprintf(stderr, "%s: -f <flow cycles> -r <repeats> -a <prob of a> -b <prob of b> -c <prob of c> -d <prob of g> -A a0,a1,... -B b0,b1,...  -C c0,c1,... -D d0,d1,...\n", name);
  exit(-1);
}

/* simulation for one sequence */
void seq_gen(int cycle, double *cumProb, double (*cumPr)[MAX_C], int *pr_cnt) {
  int f = 1;   /* cycle count */
  int n = 0;   /* sequence length */
  int i = 0, j = 0, ii, ok, lastj = 0;
  double rand;

  i = 0;  /* flow  */
  while (1) { /* cycle */
    /* create a base of the template */
    ok = 0;
    rand = drand48();
    //    lastj = j;  /* the last base that ends up in cycle f */
    for (j = 0; j < 4; j++) { /* j: base of the template */
      if (rand < cumProb[j]) {
	ok = 1;
	baseCnt[j]++;
	break;
      }
    }
    /* now j is the base to be incorporated */

    if (ok == 0) {
      printf("The base is not incorporated! (rand = %f)\n", rand);
      exit(-1);
    } 

    /* 
     * flush the flows, until the base is incorporated, or the cycle
     *  limit is reached 
     */
    while (f <= cycle) {
      // printf("<%1d:%d>", i, f);
      if (i != j) { /* flow and base do not match */
	i++;        /* i: flow */
	
	if (i > 3) {
	  i = 0;
	  f++;
	  if (f > cycle) {
	    break;
	  }
	}
      } else { /* incorporate the base based on probs */
	ok = 0;
	rand = drand48();
	/* ii: index of gi  */
	for (ii = 0; ii < pr_cnt[i]; ii++) {
	  if (rand < cumPr[i][ii]) { /* the base will be incorporated */ 
	    ok = 1;
	    break;
	  } else {
	    f++;
	    if (f > cycle) {
	      break;
	    }
	  }
	}

	if (ok == 0 && f <= cycle) {
	  printf("The base is not incorporated! (f: %d n: %d base: %d flow: %d g_cnt: %d ii: %d cumPr: %f rand: %f)\n", 
		 f, n, j, i, pr_cnt[i], ii, cumPr[i][ii], rand);
	  exit(-1);
	}      
	 
	if (f <= cycle) {
	  /* base got incorporated */
	  n++;    /* seq length */
	  lastj = j;
	  break;
	}
      }
    } /* while (f <= cycle) */
    if (f > cycle) {
      break;
    }
  }

  nSum += n;
  nSum2 += (n * n);

  if (n < MAX_N) {
    ndist[n]++;
    ndist4[n][lastj]++;
  }
}


/*
 * functions needed for theoretical results
 */

/* s2 */
double get_s2(double *prob) {
  double s23, s2;

  s23 = prob[2] + prob[3];
  s2 = prob[0] * (prob[1] + s23)
    +  prob[1] * s23
    +  prob[2] * prob[3];

  return s2;
}

/* s3 */
double get_s3(double *prob) {
  double s3;

  s3 = prob[0] * prob[1] * prob[2]
    +  prob[0] * prob[1] * prob[3]
    +  prob[0] * prob[2] * prob[3]
    +  prob[1] * prob[2] * prob[3];

  return s3;
}

/* s4 */
double get_s4(double *prob) {
  double s4;

  s4 = prob[0] * prob[1] * prob[2] * prob[3];

  return s4;
}



/* t1'(1) and "var" of t1 */
double get_t1p(double *prob, double (*pr)[MAX_C], int *pr_cnt, 
	       double *var_t1) {
  int i, j;
  double g, g2, g2sum, t1p;
  
  t1p = 0.0;
  g2sum = 0.0;
  for (i = 0; i < 4; i++) {
    g = 0.0;
    g2 = 0.0;
    for (j = 1; j < pr_cnt[i]; j++) {
      g += j * pr[i][j];
      g2 += j*(j-1) * pr[i][j];
    }
    g *= prob[i];
    g2 *= prob[i];

    t1p += g;
    g2sum += g2;
  }

  *var_t1 = g2sum + t1p - t1p*t1p;

  return t1p;
}

/* t1''(1)  */
double get_t1p2(double *prob, double (*pr)[MAX_C], int *pr_cnt) {
  int i, j;
  double g, t1p2;
  
  t1p2 = 0.0;
  for (i = 0; i < 4; i++) {
    g = 0.0;
    for (j = 2; j < pr_cnt[i]; j++) {
      g += j * (j-1) * pr[i][j];
    }
    g *= prob[i];

    t1p2 += g;
  }

  return t1p2;
}


/* t1'''(1)  */
double get_t1p3(double *prob, double (*pr)[MAX_C], int *pr_cnt) {
  int i, j;
  double g, t1p3;
  
  t1p3 = 0.0;
  for (i = 0; i < 4; i++) {
    g = 0.0;
    for (j = 3; j < pr_cnt[i]; j++) {
      g += j * (j-1) * (j-2) * pr[i][j];
    }
    g *= prob[i];

    t1p3 += g;
  }

  return t1p3;
}


/* t2'(1) -- see yale notebook 3 p12 */
double get_t2p(double *prob, double (*pr)[MAX_C], int *pr_cnt, double t1p) {
  int i, j;
  double g[4], t, t2p;
  
  for (i = 0; i < 4; i++) {
    g[i] = 0.0;
    for (j = 1; j < pr_cnt[i]; j++) {
      g[i] += j * pr[i][j];
    }
    g[i] *= prob[i]; /* g_i'(1) */
  }

  t = 0.0;
  for (i = 0; i < 4; i++) 
    t += prob[i] * g[i];

  t2p = t1p - t;

  return t2p;
}

/* t2''(1) */
double get_t2p2(double *prob, double (*pr)[MAX_C], int *pr_cnt) {
  int i, j;
  double g1[4], g2[4], t2p2;
  
  for (i = 0; i < 4; i++) {
    g1[i] = 0.0;
    g2[i] = 0.0;
    for (j = 1; j < pr_cnt[i]; j++) {
      g1[i] += j * pr[i][j];
      g2[i] += j * (j-1) * pr[i][j];
    }
    g1[i] *= prob[i]; /* g_i'(1) */
    g2[i] *= prob[i]; /* g_i''(1) */
  }

  t2p2 = 0.0;
  for (i = 0; i < 4; i++) {
    t2p2 += (g2[i] * (1.0 - prob[i]));
  }

  t2p2 += 2.0 * (g1[0] * g1[1] + g1[0] * g1[2] + g1[0] * g1[3]
		 + g1[1] * g1[2] + g1[1] * g1[3]
		 + g1[2] * g1[3]);

  return t2p2;
}


/* t3'(1) */
double get_t3p(double *prob, double (*pr)[MAX_C], int *pr_cnt) {
  int i, j;
  double g[4], t3p;
  
  for (i = 0; i < 4; i++) {
    g[i] = 0.0;
    for (j = 1; j < pr_cnt[i]; j++) {
      g[i] += j * pr[i][j];
    }
    g[i] *= prob[i]; /* g_i'(1) */
  }

  t3p = g[0] * (prob[1] * prob[2] + prob[1] * prob[3] + prob[2] * prob[3]) +
        g[1] * (prob[0] * prob[2] + prob[0] * prob[3] + prob[2] * prob[3]) +
        g[2] * (prob[0] * prob[1] + prob[0] * prob[3] + prob[1] * prob[3]) +
        g[3] * (prob[0] * prob[1] + prob[0] * prob[2] + prob[1] * prob[2]);

  return t3p;
}






int main(int argc, char *argv[]) { 
  int i, j, k;
  int cycle = 10, repeat = 1000, ch, rsum = 0, peak, zerocnt;
  int fixedSeed = 1, pad = 0;
  double prob[] = {0.25, 0.25, 0.25, 0.25};  /* prob of 4 nucleotides */
  double cumProb[4];                         /* cumulative probs */
  double pr[4][MAX_C];                       /* incorporation probs */
  double cumPr[4][MAX_C];

  double sum = 0.0;
  char *t;
  int pr_cnt[4];                  /* upper index of incorporation probs */

  double nAvg, nVar;
  double nAvgth, nVarth;
  double nfac, nfac2, nfac3, nfac4, vfac, tmp2, u, v;
  double s2, s3, s4, t1p, var_t1, t2p, t1p2, t1p3, t2p2, t3p;

  char fn[1024];
  FILE *fp;

  while ((ch = 
	  getopt(argc, argv, 
		 "a:b:c:d:sf:r:p:A:B:C:D:"
		 )) != -1)
    switch (ch) {
    case 'a':
      prob[0] = atof(optarg);             break;
    case 'b':
      prob[1] = atof(optarg);             break;      
    case 'c':
      prob[2] = atof(optarg);             break;
    case 'd':
      prob[3] = atof(optarg);             break;      
    case 's':
      fixedSeed = 0;                      break;   
    case 'p':
      pad       = 1;                      break;   
    case 'f':
      cycle = atoi(optarg);              break;   
    case 'r':
      repeat = atoi(optarg);              break;   
      

    case 'A':
      for (i = 0, t = strtok(optarg, " ,"); t; t = strtok(NULL, " ,"), i++)
	pr[0][i] = atof(t);
      
      pr_cnt[0] = i;
      break;   

    case 'B':
      for (i = 0, t = strtok(optarg, " ,"); t; t = strtok(NULL, " ,"), i++)
	pr[1][i] = atof(t);
 
      pr_cnt[1] = i;
      break;   

    case 'C':
      for (i = 0, t = strtok(optarg, " ,"); t; t = strtok(NULL, " ,"), i++)
	pr[2][i] = atof(t);
 
      pr_cnt[2] = i;
      break;   

    case 'D':
      for (i = 0, t = strtok(optarg, " ,"); t; t = strtok(NULL, " ,"), i++)
	pr[3][i] = atof(t);

      pr_cnt[3] = i;
      break;   


    default:
      usage(argv[0]);
    }
  
  /* make sure the probs add up to 1 */
  for (i = 0; i < 4; i++) {
    if (prob[i] < 0.0 || prob[i] > 1.0) {
      usage(argv[0]);
    }
    sum += prob[i];
  }

  if (fabs(sum - 1.0) > 1.0e-10)
    usage(argv[0]);


  for (i = 0; i < 4; i++) {
    sum = 0.0;
    for (j = 0; j < pr_cnt[i]; j++) {

      if (pr[i][j] < 0.0 || pr[i][j] > 1.0) {
	usage(argv[0]);
      }
      
      sum += pr[i][j];
    }
    if (fabs(sum - 1.0) > 1.0e-10) {
      printf("%c\t%f\n", i + 65, sum);
      usage(argv[0]);
    }
  }

  /* set up a unique file name */
  sprintf(fn, "sbs_ffcm_f_%d_repeat_%d", cycle, repeat);
  for (i=0; i<4; i++) {
    sprintf(fn, "%s_%c%4.3f", fn, i+65, prob[i]);
    for (j=0; j<pr_cnt[i]; j++)
      sprintf(fn, "%s_%3.2f", fn, pr[i][j]);
  }
  printf("%s\n", fn);
  fp = fopen(fn, "w");
  if (fp == NULL) {
    fprintf(stderr, "cannot open file %s to write!\n", fn);
    exit(-1);
  }


  fprintf(fp, "######## Input parameters ###########################################\n");
  if (fixedSeed) {
    srand48(12345);
    fprintf(fp, "#\tseed: fixed as 12345\n");
  } else {
    srand48(time(NULL));
    fprintf(fp, "#\tseed: use time\n");
  }

  fprintf(fp, "#\tflow cycle: %d\n", cycle);
  fprintf(fp, "#\trepeat: %d\n", repeat);


  /* set up cumulative probs */
  cumProb[0] = prob[0];
  fprintf(fp, "\n#\tnucleotide composition probabilities\tcumulative\n");  
  fprintf(fp, "#\t%c\t%f\t%f\n", 'a', prob[0], cumProb[0]);
  for (i = 1; i < 4; i++) {
    cumProb[i] = cumProb[i-1] + prob[i];
    fprintf(fp, "#\t%c\t%f\t%f\n", i+97, prob[i], cumProb[i]);
  }

  fprintf(fp, "\n#\tnucleotide incorporation probabilities\tcumulative\n"); 
  for (i = 0; i < 4; i++) {
    cumPr[i][0] = pr[i][0];
    fprintf(fp, "#\t%c\t%d\t%f\t%f\n", i+97, 0, pr[i][0], cumPr[i][0]);
    for (j = 1; j < pr_cnt[i]; j++) {
      cumPr[i][j] = cumPr[i][j-1] + pr[i][j];
      fprintf(fp, "#\t%c\t%d\t%f\t%f\n", i+97, j, pr[i][j], cumPr[i][j]);
    }
    fprintf(fp, "\n");
  }

  
  /* init */
  for (i = 0; i < 4; i++) {
    baseCnt[i] = 0;
  }

  
  for (j = 0; j < MAX_N; j++) {
    ndist[j] = 0.0;
    for (i = 0; i < 4; i++)
      ndist4[j][i] = 0.0; 
  }


  /*** run the simulation ***/
  for (k = 0; k < repeat; k++) {
    seq_gen(cycle, cumProb, cumPr, pr_cnt);
  }

  
  fprintf(fp, "#####################################################################\n");
  fprintf(fp, "######## Theoretical results of sequence length distribution ########\n");

  /* avg and var from theory */
  s2 = get_s2(prob);
  s3 = get_s3(prob);
  s4 = prob[0] * prob[1] * prob[2] * prob[3]; 
  t1p  = get_t1p(prob, pr, pr_cnt, &var_t1);
  t2p  = get_t2p(prob, pr, pr_cnt, t1p);
  t1p2 = get_t1p2(prob, pr, pr_cnt);
  t1p3 = get_t1p3(prob, pr, pr_cnt);
  t2p2 = get_t2p2(prob, pr, pr_cnt);
  t3p  = get_t3p(prob, pr, pr_cnt);

  nfac = s2 + t1p;
  vfac = 2.0*s3 + s2 - 3.0*s2*s2 + t1p2 + t1p - t1p*t1p + 2.0*t2p 
    - 4.0*t1p*s2;  /* w in jcb p824 */

  nfac2 = nfac  * nfac;
  nfac3 = nfac2 * nfac;
  nfac4 = nfac3 * nfac;

  // 4/15/2011
  u = nfac;
  v = 2*t2p + t1p2 + 2*s3;

  tmp2 = v*v;


  /* 4/25/2010 */
  nAvgth = 1.0/nfac * cycle 
    + 0.5/(nfac2) * (-2.0*s2*s2 + 2.0*s3 
		     - 4.0*t1p*s2 - 2.0*t1p*t1p + t1p2 + 2.0*t2p); 
  
  // 4/15/2011 a newer version
  nVarth = 1.0/(nfac3) * (v - (3*s2 + t1p -1)*u) * cycle 
    + 1.0/12.0/(nfac4) * (
			  6.0*u*v*(3*u-4*s2) + 15.0*v*v
			  -8.0*u*(
				  6.0*t3p + t1p3 + 3.0*t2p2 + 6.0*s4
				  + 3.0*u*(t1p2 + t2p)
				  )
			  );
 


  fprintf(fp, "\tAvg\t%f\tVar\t%f\n", 
	  nAvgth, 
	  nVarth);

  if (nAvgth >= 0.5*MAX_N) {
    fprintf(stderr, "MAX_N should be increased!\n");
    exit(-1);
  }


  fprintf(fp, "######## Simulation results of sequence length distribution #########\n");
  /* avg and var from simulation */
  nAvg = nSum / (double)repeat;
  nVar = nSum2 / (double)repeat -  nAvg * nAvg;
  fprintf(fp, "\tAvg\t%f\tVar\t%f\n", 
	  nAvg, 
	  nVar);

  fprintf(fp, "#####################################################################\n");


  fprintf(fp, "\n");
  fprintf(fp, "######### base counts \n");
  for (i = 0; i < 4; i++) {
    fprintf(fp, "#\t%c\t%lld\n", i+97, baseCnt[i]);
  }



  /* use a simple peak finder to truncate the unneccessary printout */
  fprintf(fp, "\n");
  fprintf(fp, "######### length distribution########################################\n");
  fprintf(fp, "n\tabcd_count\tabcd_frac\ta_count\ta_frac\tb_count\tb_frac\tc_count\tc_frac\td_count\td_frac\n");
  peak = 0;
  zerocnt = 0;
  for (j=0; j<MAX_N; j++) {
    fprintf(fp, "%d\t%f\t%f\t", j, ndist[j], ndist[j]/repeat);
    for (i = 0; i < 4; i++)
      fprintf(fp, "%f\t%f\t", ndist4[j][i], ndist4[j][i]/repeat);
    fprintf(fp, "\n");
    rsum += ndist[j];

    if (ndist[j] > 0) 
      peak = 1;

    if (peak && ndist[j] == 0.0)
      zerocnt++;

    if (zerocnt > 100)
      break;
  }
  //  fprintf(fp, "Total sequences: %d\n", rsum);




  fclose(fp);
  exit (0);
}
