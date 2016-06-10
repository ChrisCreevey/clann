
/* version 3.6. (c) Copyright 1993-2004 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, Andrew Keeffe,
   and Dan Fineman.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#ifdef OSX_CARBON
#include <Carbon/Carbon.h>
#endif

#include <stdio.h>
#include <signal.h>
#ifdef WIN32
#include <windows.h>
/* for console code (clear screen, text color settings) */
CONSOLE_SCREEN_BUFFER_INFO savecsbi;
HANDLE hConsoleOutput;

void phyClearScreen();
void phySaveConsoleAttributes();
void phySetConsoleAttributes();
void phyRestoreConsoleAttributes();
void phyFillScreenColor();
#endif

#include "phylip.h"

#ifndef OLDC
static void crash_handler(int signum);

#endif
#if defined(OSX_CARBON) && defined(__MWERKS__)
boolean fixedpath = false;
#endif
FILE *infile, *outfile, *intree, *intree2, *outtree, *weightfile, *catfile, *ancfile, *mixfile, *factfile;
long spp, words, bits;
boolean ibmpc, ansi, tranvsp;
naym *nayme;                     /* names of species */

static void crash_handler(int sig_num)
{ /* when we crash, lets print out something usefull */
  printf("ERROR:  ");
  switch(sig_num) {
#ifdef SIGSEGV
    case SIGSEGV:
      puts("This program has caused a Segmentation fault.");
      break;
#endif /* SIGSEGV */
#ifdef SIGFPE
    case SIGFPE:
      puts("This program has caused a Floating Point Exception");
      break;
#endif  /* SIGFPE */
#ifdef SIGILL
    case SIGILL:
      puts("This program has attempted an illegal instruction");
      break;
#endif  /* SIGILL */
#ifdef SIGPIPE 
    case SIGPIPE:
      puts("This program tried to write to a broken pipe");
      break;
#endif  /* SIGPIPE */
#ifdef SIGBUS
    case SIGBUS:
      puts("This program had a bus error");
      break;
#endif /* SIGBUS */
  }   
  if (sig_num == SIGSEGV) {
    puts(
    "       This may have been caused by an incorrectly formatted input file");
    puts(
        "       or input tree file.  You should check those files carefully.");
    puts("       If this seems to be a bug, please mail joe@gs.washington.edu");
  }
  else {
    puts("       Most likely, you have encountered a bug in the program.");
    puts("       Since this seems to be a bug, please mail joe@gs.washington.edu");
  }
  puts("       with the name of the program, your computer system type,");
  puts("       a full description of the problem, and with the input data file.");
  puts("       (which should be in the body of the message, not as an Attachment).");

#ifdef WIN32
  puts ("Press Enter or Return to close program.");
  puts("  You may have to press Enter or Return twice.");
  getchar ();
  getchar ();
  phyRestoreConsoleAttributes();
#endif
  abort();
}


void init(int argc, char** argv) 
{ /* initialization routine for all programs 
   * anything done at the beginig for every program should be done here */ 
 
  /* set up signal handler for 
   * segfault,floating point exception, illeagal instruction, bad pipe, bus error
   * there are more signals that can cause a crash, but these are the most common
   * even these aren't found on all machines.  */
#ifdef SIGSEGV
  signal(SIGSEGV, crash_handler);
#endif /* SIGSEGV */
#ifdef SIGFPE
  signal(SIGFPE, crash_handler);
#endif /* SIGFPE */
#ifdef SIGILL
  signal(SIGILL, crash_handler);
#endif /* SIGILL */
#ifdef SIGPIPE
  signal(SIGPIPE, crash_handler);
#endif /* SIGPIPE */
#ifdef SIGBUS
  signal(SIGBUS, crash_handler);
#endif /* SIGBUS */

#ifdef WIN32
  phySetConsoleAttributes();
  phyClearScreen();
#endif

}

void scan_eoln(FILE *f) 
{ /* eat everything to the end of line or eof*/
  char ch;

  while (!eoff(f) && !eoln(f)) 
    gettc(f);
  if (!eoff(f)) 
    ch = gettc(f);
}


boolean eoff(FILE *f)
{ /* check for end of file */
    int ch;

    if (feof(f)) 
      return true;
    ch = getc(f);
    if (ch == EOF) {
      ungetc(ch, f);
      return true;
    }
    ungetc(ch, f);
    return false;
}  /*eoff*/


boolean eoln(FILE *f)
{ /* check for end of line or eof*/
    register int ch;

    ch = getc(f);
    if (ch == EOF)
      return true;
    ungetc(ch, f);
    return ((ch == '\n') || (ch == '\r'));
}  /*eoln*/


int filexists(char *filename)
{ /* check whether file already exists */
  FILE *fp;
  fp =fopen(filename,"r");
  if (fp) {
    fclose(fp);
    return 1;
  } else
    return 0;
}  /*filexists*/


const char* get_command_name (const char *vektor)
{ /* returns the name of the program from vektor without the whole path */
  char *last_slash;

  /* Point to the last slash... */
  last_slash = strrchr (vektor, DELIMITER);

  if (last_slash)
    /* If there was a last slash, return the character after it */
    return last_slash + 1;
  else
    /* If not, return the vector */
    return vektor;

}  /*get_command_name*/


void getstryng(char *fname)
{ /* read in a file name from stdin and take off newline if any */

  fname = fgets(fname, 100, stdin);
  if (strchr(fname, '\n') != NULL)
    *strchr(fname, '\n') = '\0';
} /* getstryng */


void countup(long *loopcount, long maxcount)
{ /* count how many times this loop has tried to read data, bail out
     if exceeds maxcount */

  (*loopcount)++;
  if ((*loopcount) >= maxcount) {
    printf("\nERROR: Made %ld attempts to read input in loop. Aborting run.\n",
            *loopcount);
    exxit(-1);
  }
} /* countup */


void openfile(FILE **fp,const char *filename,const char *filedesc,
              const char *mode,const char *application, char *perm)
{ /* open a file, testing whether it exists etc. */
  FILE *of;
  char file[FNMLNGTH];
  char filemode[2];
  char input[FNMLNGTH];
  char ch;
  const char *progname_without_path;
  long loopcount, loopcount2;
#if defined(OSX_CARBON) && defined(__MWERKS__)
  ProcessSerialNumber myProcess;
  FSRef bundleLocation;
  unsigned char bundlePath[FNMLNGTH];

  if(!fixedpath){
    /* change path to the bundle location instead of root directory */
    GetCurrentProcess(&myProcess);
    GetProcessBundleLocation(&myProcess, &bundleLocation);
    FSRefMakePath(&bundleLocation, bundlePath, FNMLNGTH);
    chdir((const char*)bundlePath);
    chdir(".."); /* get out of the .app directory */
    
    fixedpath = true;
  }
#endif

  progname_without_path = get_command_name(application);

  strcpy(file,filename);
  strcpy(filemode,mode);
  loopcount = 0;
  while (1){
    if (filemode[0] == 'w' && filexists(file)){
      printf("\n%s: the file \"%s\" that you wanted to\n",
          progname_without_path, file);
      printf("     use as %s already exists.\n", filedesc);
      printf("     Do you want to Replace it, Append to it,\n");
      printf("     write to a new File, or Quit?\n");
      loopcount2 = 0;
      do {
        printf("     (please type R, A, F, or Q) \n");
#ifdef WIN32
        phyFillScreenColor();
#endif
        fgets(input, sizeof(input), stdin);
        ch  = input[0];
        uppercase(&ch);
        countup(&loopcount2, 10);
      } while (ch != 'A' && ch != 'R' && ch != 'F' && ch != 'Q');
      if (ch == 'Q')
        exxit(-1);
      if (ch == 'A') {
        strcpy(filemode,"a");
        continue;
      }
      else if (ch == 'F') {
        file[0] = '\0';
        loopcount2 = 0;
        while (file[0] =='\0') {
          printf("Please enter a new file name> ");
          getstryng(file);
          countup(&loopcount2, 10);
        }
        strcpy(filemode,"w");
        continue;
      }
    }
    of = fopen(file,filemode);
    if (of)
      break;
    else {
      switch (filemode[0]){

      case 'r':
        printf("%s: can't find %s \"%s\"\n", progname_without_path,
            filedesc, file);
        file[0] = '\0';
        loopcount2 = 0;
        while (file[0] =='\0'){
          printf("Please enter a new file name> ");
          countup(&loopcount2, 10);
          getstryng(file);}
        break;

      case 'w':
      case 'a':
        printf("%s: can't write %s \"%s\"\n", progname_without_path,
            filedesc, file);
        file[0] = '\0';
        loopcount2 = 0;
        while (file[0] =='\0'){
          printf("Please enter a new file name> ");
          countup(&loopcount2, 10);
          getstryng(file);}
        continue;
      default:
     printf("There is some error in the call of openfile. Unknown mode.\n");
        exxit(-1);
      }
    }
    countup(&loopcount, 20);
  }
  *fp = of;
  if (perm != NULL)
    strcpy(perm,file);
} /* openfile */


void cleerhome()
{ /* home cursor and clear screen, if possible */
#ifdef WIN32
  if(ibmpc || ansi){
    phyClearScreen();
  } else {
    printf("\n\n");
  }
#else
  printf("%s", ((ibmpc || ansi) ? ("\033[2J\033[H") : "\n\n"));
#endif
} /* cleerhome */


double randum(longer seed)
{ /* random number generator -- slow but machine independent
  This is a multiplicative congruential 32-bit generator
  x(t+1) = 1664525 * x(t) mod 2^32, one that passes the
  Coveyou-Macpherson and Lehmer tests, see Knuth ACP vol. 2
  We here implement it representing each integer in base-64
  notation -- i.e. as an array of 6 six-bit chunks   */
  long i, j, k, sum;
  longer mult, newseed;
  double x;

  mult[0] = 13;   /* these four statements set the multiplier */
  mult[1] = 24;   /* -- they are its "digits" in a base-64    */
  mult[2] = 22;   /*    notation: 1664525 = 13*64^3+24*64^2   */
  mult[3] = 6;    /*                         +22*64+6         */
  for (i = 0; i <= 5; i++)
    newseed[i] = 0;
  for (i = 0; i <= 5; i++) {
    sum = newseed[i];
    k = i;
    if (i > 3)
      k = 3;
    for (j = 0; j <= k; j++)
      sum += mult[j] * seed[i - j];
    newseed[i] = sum;
    for (j = i; j <= 4; j++) {
      newseed[j + 1] += newseed[j] / 64;
      newseed[j] &= 63;
    }
  }
  memcpy(seed, newseed, sizeof(longer));
  seed[5] &= 3;
  x = 0.0;
  for (i = 0; i <= 5; i++)
    x = x / 64.0 + seed[i];
  x /= 4.0;
  return x;
}  /* randum */


void randumize(longer seed, long *enterorder)
{ /* randomize input order of species */
  long i, j, k;

  for (i = 0; i < spp; i++) {
    j = (long)(randum(seed) * (i+1));
    k = enterorder[j];
    enterorder[j] = enterorder[i];
    enterorder[i] = k;
  }
} /* randumize */


double normrand(longer seed)
{/* standardized Normal random variate */
  double x;

  x = randum(seed)+randum(seed)+randum(seed)+randum(seed)
       + randum(seed)+randum(seed)+randum(seed)+randum(seed)
       + randum(seed)+randum(seed)+randum(seed)+randum(seed)-6.0;
  return(x);
} /* normrand */ 


long readlong(const char *prompt)
{ /* read a long */
  long res, loopcount;
  char string[100];

  loopcount = 0;
  do {
    printf("%s",prompt);
    getstryng(string);
    if (sscanf(string,"%ld",&res) == 1)
      break;
    countup(&loopcount, 10);
   } while (1);
  return res;
}  /* readlong */


void uppercase(Char *ch)
{ /* convert ch to upper case */
  *ch = (islower (*ch) ? toupper(*ch) : (*ch));
}  /* uppercase */


void initseed(long *inseed, long *inseed0, longer seed)
{ /* input random number seed */
  long i, loopcount;

  loopcount = 0;
  do {
    printf("Random number seed (must be odd)?\n");
    scanf("%ld%*[^\n]", inseed);
    getchar();
    countup(&loopcount, 10);
  } while (((*inseed) < 0) || ((*inseed) & 1) == 0);
  *inseed0 = *inseed;
  for (i = 0; i <= 5; i++)
    seed[i] = 0;
  i = 0;
  do {
    seed[i] = *inseed & 63;
    *inseed /= 64;
    i++;
  } while (*inseed != 0);
}  /*initseed*/


void initjumble(long *inseed, long *inseed0, longer seed, long *njumble)
{ /* input number of jumblings for jumble option */
  long loopcount;

  initseed(inseed, inseed0, seed);
  loopcount = 0;
  do {
    printf("Number of times to jumble?\n");
    scanf("%ld%*[^\n]", njumble);
    getchar();
    countup(&loopcount, 10);
  } while ((*njumble) < 1);
}  /*initjumble*/


void initoutgroup(long *outgrno, long spp)
{ /* input outgroup number */
  long loopcount;
  boolean done;

  loopcount = 0;
  do {
    printf("Type number of the outgroup:\n");
    scanf("%ld%*[^\n]", outgrno);
    getchar();
    done = (*outgrno >= 1 && *outgrno <= spp);
    if (!done) {
      printf("BAD OUTGROUP NUMBER: %ld\n", *outgrno);
      printf("  Must be in range 1 - %ld\n", spp);
    }
    countup(&loopcount, 10);
  } while (done != true);
}  /*initoutgroup*/


void initthreshold(double *threshold)
{ /* input threshold for threshold parsimony option */
  long loopcount;
  boolean done;

  loopcount = 0;
  do {
    printf("What will be the threshold value?\n");
    scanf("%lf%*[^\n]", threshold);
    getchar();
    done = (*threshold >= 1.0);
    if (!done)
      printf("BAD THRESHOLD VALUE:  it must be greater than 1\n");
    else
      *threshold = (long)(*threshold * 10.0 + 0.5) / 10.0;
    countup(&loopcount, 10);
  } while (done != true);
}  /*initthreshold*/


void initcatn(long *categs)
{ /* initialize category number for rate categories */
  long loopcount;

  loopcount = 0;
  *categs = 0;
  do {
    printf("Number of categories (1-%d)?\n", maxcategs);
    scanf("%ld%*[^\n]", categs);
    getchar();
    countup(&loopcount, 10);
  } while (*categs > maxcategs || *categs < 1);
}  /*initcatn*/


void initcategs(long categs, double *rate)
{ /* initialize category rates for HMM rates */
  long i, loopcount, scanned;
  char line[100], rest[100];
  boolean done;

  loopcount = 0;
  for (;;){
    printf("Rate for each category? (use a space to separate)\n");
    getstryng(line);
    done = true;
    for (i = 0; i < categs; i++){
      scanned = sscanf(line,"%lf %[^\n]", &rate[i],rest);
      if ((scanned < 2 && i < (categs - 1)) ||
       (scanned < 1 && i == (categs - 1))){
     printf("Please enter exactly %ld values.\n",categs);
     done = false;
     break;
      }
      strcpy(line,rest);
    }
    if (done)
      break;
    countup(&loopcount, 100);
  }
}  /*initcategs*/


void initprobcat(long categs, double *probsum, double *probcat)
{ /* input probabilities of rate categores for HMM rates */
  long i, loopcount, scanned;
  boolean done;
  char line[100], rest[100];

  loopcount = 0;
  do {
    printf("Probability for each category?");
    printf(" (use a space to separate)\n");
    getstryng(line);
    done = true;
    for (i = 0; i < categs; i++){
      scanned = sscanf(line,"%lf %[^\n]",&probcat[i],rest);
      if ((scanned < 2 && i < (categs - 1)) ||
       (scanned < 1 && i == (categs - 1))){
     done = false;
     printf("Please enter exactly %ld values.\n",categs);
     break;}
      strcpy(line,rest);
    }
    if (!done)
      continue;
    *probsum = 0.0;
    for (i = 0; i < categs; i++)
      *probsum += probcat[i];
    if (fabs(1.0 - (*probsum)) > 0.001) {
      done = false;
      printf("Probabilities must add up to");
      printf(" 1.0, plus or minus 0.001.\n");
    }
    countup(&loopcount, 100);
  } while (!done);
}  /*initprobcat*/


void lgr(long m, double b, raterootarray lgroot)
{ /* For use by initgammacat.  Get roots of m-th Generalized Laguerre
     polynomial, given roots of (m-1)-th, these are to be
     stored in lgroot[m][] */
  long i;
  double upper, lower, x, y;
  boolean dwn;   /* is function declining in this interval? */

  if (m == 1) {
    lgroot[1][1] = 1.0+b;
  } else {
    dwn = true;
    for (i=1; i<=m; i++) {
      if (i < m) {
        if (i == 1)
          lower = 0.0;
        else
          lower = lgroot[m-1][i-1];
        upper = lgroot[m-1][i];
      } else {                 /* i == m, must search above */
        lower = lgroot[m-1][i-1];
        x = lgroot[m-1][m-1];
        do {
          x = 2.0*x;
          y = glaguerre(m, b,x);
        } while ((dwn && (y > 0.0)) || ((!dwn) && (y < 0.0)));
        upper = x;
      }
      while (upper-lower > 0.000000001) {
        x = (upper+lower)/2.0;
        if (glaguerre(m, b, x) > 0.0) {
          if (dwn)
            lower = x;
          else
            upper = x;
        } else {
          if (dwn)
            upper = x;
          else
            lower = x;
        }        
      }
      lgroot[m][i] = (lower+upper)/2.0;
      dwn = !dwn;                /* switch for next one */
    }
  }
} /* lgr */


double logfac (long n)
{ /* log(n!) values were calculated with Mathematica
     with a precision of 30 digits */
  long i;
  double x;

  switch (n)
    {
    case 0:
      return 0.;
    case 1:
      return 0.;
    case 2:
      return 0.693147180559945309417232121458;
    case 3:
      return 1.791759469228055000812477358381;
    case 4:
      return 3.1780538303479456196469416013;
    case 5:
      return 4.78749174278204599424770093452;
    case 6:
      return 6.5792512120101009950601782929;
    case 7:
      return 8.52516136106541430016553103635;
    case 8:
      return 10.60460290274525022841722740072;
    case 9:
      return 12.80182748008146961120771787457;
    case 10:
      return 15.10441257307551529522570932925;
    case 11:
      return 17.50230784587388583928765290722;
    case 12:
      return 19.98721449566188614951736238706;
    default:
      x = 19.98721449566188614951736238706;
      for (i = 13; i <= n; i++)
        x += log(i);
      return x;
    }
}

                        
double glaguerre(long m, double b, double x)
{ /* Generalized Laguerre polynomial computed recursively.
     For use by initgammacat */
  long i;
  double gln, glnm1, glnp1; /* L_n, L_(n-1), L_(n+1) */

  if (m == 0)
    return 1.0;
  else {
    if (m == 1)
      return 1.0 + b - x;
    else {
      gln = 1.0+b-x;
      glnm1 = 1.0;
      for (i=2; i <= m; i++) {
        glnp1 = ((2*(i-1)+b+1.0-x)*gln - (i-1+b)*glnm1)/i;
        glnm1 = gln;
        gln = glnp1;
      }
      return gln;
    }
  }
} /* glaguerre */


void initlaguerrecat(long categs, double alpha, double *rate, double *probcat)
{ /* calculate rates and probabilities to approximate Gamma distribution
     of rates with "categs" categories and shape parameter "alpha" using
     rates and weights from Generalized Laguerre quadrature */
  long i;
  raterootarray lgroot; /* roots of GLaguerre polynomials */
  double f, x, xi, y;

  alpha = alpha - 1.0;
  lgroot[1][1] = 1.0+alpha;
  for (i = 2; i <= categs; i++)
    lgr(i, alpha, lgroot);                   /* get roots for L^(a)_n */
  /* here get weights */
  /* Gamma weights are (1+a)(1+a/2) ... (1+a/n)*x_i/((n+1)^2 [L_{n+1}^a(x_i)]^2)  */
  f = 1;
  for (i = 1; i <= categs; i++)
    f *= (1.0+alpha/i);
  for (i = 1; i <= categs; i++) {
    xi = lgroot[categs][i];
    y = glaguerre(categs+1, alpha, xi);
    x = f*xi/((categs+1)*(categs+1)*y*y);
    rate[i-1] = xi/(1.0+alpha);
    probcat[i-1] = x;
  }
} /* initlaguerrecat */


double hermite(long n, double x)
{ /* calculates hermite polynomial with degree n and parameter x */
  /* seems to be unprecise for n>13 -> root finder does not converge*/
  double h1 = 1.;
  double h2 = 2. * x;
  double xx = 2. * x;
  long i;

  for (i = 1; i < n; i++) {
    xx = 2. * x * h2 - 2. * (i) * h1;
    h1 = h2;
    h2 = xx;
  }
  return xx;
} /* hermite */


void root_hermite(long n, double *hroot)
{ /* find roots of Hermite polynmials */
  long z;
  long ii;
  long start;

  if (n % 2 == 0) {
    start = n/2;
    z = 1;
  } else {
    start = n/2 + 1;
    z=2;
    hroot[start-1] = 0.0;
  }
  for (ii = start; ii < n; ii++) {         /* search only upwards*/
    hroot[ii] = halfroot(hermite,n,hroot[ii-1]+EPSILON, 1./n);
    hroot[start - z] = -hroot[ii];
    z++;
  }
} /* root_hermite */


double halfroot(double (*func)(long m, double x), long n, double startx,
                double delta)
{ /* searches from the bound (startx) only in one direction
     (by positive or negative delta, which results in
     other-bound=startx+delta)
     delta should be small.
     (*func) is a function with two arguments  */
  double xl;
  double xu;
  double xm;
  double fu;
  double fl;
  double fm = 100000.;
  double gradient;
  boolean dwn;

  /* decide if we search above or below startx and escapes to trace back
     to the starting point that most often will be
     the root from the previous calculation */
  if (delta < 0) {
    xu = startx;
    xl = xu + delta;
  } else {
    xl = startx;
    xu = xl + delta;
  }
  delta = fabs(delta);
  fu = (*func)(n, xu);
  fl = (*func)(n, xl);
  gradient = (fl-fu)/(xl-xu);
  while(fabs(fm) > EPSILON) {        /* is root outside of our bracket?*/
    if ((fu<0.0 && fl<0.0) || (fu>0.0 && fl > 0.0)) {
      xu += delta;
      fu = (*func)(n, xu);
      fl = (*func)(n, xl);
      gradient = (fl-fu)/(xl-xu);
      dwn = (gradient < 0.0) ? true : false;
    } else {
      xm = xl - fl / gradient;
      fm = (*func)(n, xm);
      if (dwn) {
        if (fm > 0.) {
          xl = xm;
          fl = fm;
        } else {
          xu = xm;
          fu = fm;
        }
      } else {
        if (fm > 0.) {
          xu = xm;
          fu = fm;
        } else {
          xl = xm;
          fl = fm;
        }
      }
      gradient = (fl-fu)/(xl-xu);
    }
  }
  return xm;
} /* halfroot */


void hermite_weight(long n, double * hroot, double * weights)
{
  /* calculate the weights for the hermite polynomial at the roots
     using formula from Abramowitz and Stegun chapter 25.4.46 p.890 */
  long i;
  double hr2;
  double numerator;

  numerator = exp(0.6931471805599 * ( n-1.) + logfac(n)) / (n*n);
  for (i = 0; i < n; i++) {
    hr2 = hermite(n-1, hroot[i]);
    weights[i] = numerator / (hr2*hr2);
  }
} /* hermiteweight */


void inithermitcat(long categs, double alpha, double *rate, double *probcat)
{ /* calculates rates and probabilities */
  long i;
  double *hroot;
  double std;

  std = SQRT2 /sqrt(alpha);
  hroot = (double *) Malloc((categs+1) * sizeof(double));
  root_hermite(categs, hroot);         /* calculate roots */
  hermite_weight(categs, hroot, probcat);  /* set weights */
  for (i=0; i<categs; i++) {           /* set rates */
    rate[i] = 1.0 + std * hroot[i];
    probcat[i] = probcat[i];
  }
  free(hroot);
} /* inithermitcat */


void initgammacat (long categs, double alpha, double *rate, double *probcat)
{ /* calculate rates and probabilities to approximate Gamma distribution
   of rates with "categs" categories and shape parameter "alpha" using
   rates and weights from Generalized Laguerre quadrature or from
   Hermite quadrature */

  if (alpha >= 100.0)
    inithermitcat(categs, alpha, rate, probcat);
  else
    initlaguerrecat(categs, alpha, rate, probcat);
} /* initgammacat */


void inithowmany(long *howmany, long howoften)
{/* input how many cycles */
  long loopcount;

  loopcount = 0;
  do { 
    printf("How many cycles of %4ld trees?\n", howoften);
    scanf("%ld%*[^\n]", howmany);
    getchar();
    countup(&loopcount, 10);
  } while (*howmany <= 0);
}  /*inithowmany*/



void inithowoften(long *howoften)
{ /* input how many trees per cycle */
  long loopcount;

  loopcount = 0;
  do {
    printf("How many trees per cycle?\n");
    scanf("%ld%*[^\n]", howoften);
    getchar();
    countup(&loopcount, 10);
  } while (*howoften <= 0);
}  /*inithowoften*/


void initlambda(double *lambda)
{ /* input patch length parameter for autocorrelated HMM rates */
  long loopcount;

  loopcount = 0;
  do {
    printf("Mean block length of sites having the same rate (greater than 1)?\n");
    scanf("%lf%*[^\n]", lambda);
    getchar();
    countup(&loopcount, 10);
  } while (*lambda <= 1.0);
  *lambda = 1.0 / *lambda;
}  /*initlambda*/


void initfreqs(double *freqa, double *freqc, double *freqg, double *freqt)
{ /* input frequencies of the four bases */
  char input[100];
  long scanned, loopcount;

  printf("Base frequencies for A, C, G, T/U (use blanks to separate)?\n");
  loopcount = 0;
  do {
    getstryng(input);
    scanned = sscanf(input,"%lf%lf%lf%lf%*[^\n]", freqa, freqc, freqg, freqt);
    if (scanned == 4)
      break;
    else
      printf("Please enter exactly 4 values.\n");
    countup(&loopcount, 100);
  } while (1);
}  /* initfreqs */


void initratio(double *ttratio)
{ /* input transition/transversion ratio */
  long loopcount;

  loopcount = 0;
  do {
    printf("Transition/transversion ratio?\n");
    scanf("%lf%*[^\n]", ttratio);
    getchar();
    countup(&loopcount, 10);
  } while (*ttratio < 0.0);
}  /* initratio */


void initpower(double *power)
{
  printf("New power?\n");
  scanf("%lf%*[^\n]", power);
  getchar();
}  /*initpower*/


void initdatasets(long *datasets)
{
  /* handle multi-data set option */
  long loopcount;
  boolean done;

  loopcount = 0;
  do {
    printf("How many data sets?\n");
    scanf("%ld%*[^\n]", datasets);
    getchar();
    done = (*datasets > 1);
      if (!done)
     printf("Bad data sets number:  it must be greater than 1\n");
    countup(&loopcount, 10);
  } while (!done);
} /* initdatasets */


void justweights(long *datasets)
{
  /* handle multi-data set option by weights */
  long loopcount;
  boolean done;

  loopcount = 0;
  do {
    printf("How many sets of weights?\n");
    scanf("%ld%*[^\n]", datasets);
    getchar();
    done = (*datasets >= 1);
      if (!done)
     printf("BAD NUMBER:  it must be greater than 1\n");
    countup(&loopcount, 10);
  } while (!done);
} /* justweights */


void initterminal(boolean *ibmpc, boolean *ansi)
{
  /* handle terminal option */
  if (*ibmpc) {
    *ibmpc = false;
    *ansi = true;
  } else if (*ansi)
      *ansi = false;
    else
      *ibmpc = true;
}  /*initterminal*/


void initnumlines(long *screenlines)
{
  long loopcount;

  loopcount = 0;
  do {
    *screenlines = readlong("Number of lines on screen?\n");
    countup(&loopcount, 10);
  } while (*screenlines <= 12);
}  /*initnumlines*/


void initbestrees(bestelm *bestrees, long maxtrees, boolean glob)
{
  /* initializes either global or local field of each array in bestrees */
  long i;

  if (glob)
    for (i = 0; i < maxtrees; i++)
      bestrees[i].gloreange = false;
  else
    for (i = 0; i < maxtrees; i++)
      bestrees[i].locreange = false;
} /* initbestrees */


void newline(FILE *filename, long i, long j, long k)
{
  /* go to new line if i is a multiple of j, indent k spaces */
  long m;

  if ((i - 1) % j != 0 || i <= 1)
    return;
  putc('\n', filename);
  for (m = 1; m <= k; m++)
    putc(' ', filename);
}  /* newline */


void inputnumbersold(long *spp, long *chars, long *nonodes, long n)
{
  /* input the numbers of species and of characters */

  if (fscanf(infile, "%ld%ld", spp, chars) != 2 || *spp <= 0 || *chars <= 0) {
    printf(
    "ERROR: Unable to read the number of species or characters in data set\n");
    printf(
      "The input file is incorrect (perhaps it was not saved text only).\n");
  }
  *nonodes = *spp * 2 - n;
}  /* inputnumbersold */


void inputnumbers(long *spp, long *chars, long *nonodes, long n)
{
  /* input the numbers of species and of characters */

  if (fscanf(infile, "%ld%ld", spp, chars) != 2 || *spp <= 0 || *chars <= 0) {
    printf(
    "ERROR: Unable to read the number of species or characters in data set\n");
    printf(
      "The input file is incorrect (perhaps it was not saved text only).\n");
  }
  *nonodes = *spp * 2 - n;
}  /* inputnumbers */


void inputnumbers2(long *spp, long *nonodes, long n)
{
  /* read species number */

  if (fscanf(infile, "%ld", spp) != 1 || *spp <= 0) {
    printf("ERROR: Unable to read the number of species in data set\n");
    printf(
      "The input file is incorrect (perhaps it was not saved text only).\n");
  }
  fprintf(outfile, "\n%4ld Populations\n", *spp);
  *nonodes = *spp * 2 - n;
}  /* inputnumbers2 */


void inputnumbers3(long *spp, long *chars)
{
  /* input the numbers of species and of characters */

  if (fscanf(infile, "%ld%ld", spp, chars) != 2 || *spp <= 0 || *chars <= 0) {
    printf(
    "ERROR: Unable to read the number of species or characters in data set\n");
    printf(
       "The input file is incorrect (perhaps it was not saved text only).\n");
    exxit(-1);
  }
}  /* inputnumbers3 */


void samenumsp(long *chars, long ith)
{
  /* check if spp is same as the first set in other data sets */
  long cursp, curchs;

  if (eoln(infile)) 
    scan_eoln(infile);
  fscanf(infile, "%ld%ld", &cursp, &curchs);
  if (cursp != spp) {
    printf(
         "\n\nERROR: Inconsistent number of species in data set %ld\n\n", ith);
    exxit(-1);
  }
  *chars = curchs;
} /* samenumsp */


void samenumsp2(long ith)
{
  /* check if spp is same as the first set in other data sets */
  long cursp;

  if (eoln(infile)) 
    scan_eoln(infile);
  if (fscanf(infile, "%ld", &cursp) != 1) {
    printf("\n\nERROR: Unable to read number of species in data set %ld\n",
      ith);
    printf(
      "The input file is incorrect (perhaps it was not saved text only).\n");
    exxit(-1);
  }
  if (cursp != spp) {
    printf(
      "\n\nERROR: Inconsistent number of species in data set %ld\n\n", ith);
    exxit(-1);
  }
} /* samenumsp2 */


void readoptions(long *extranum, const char *options)
{ /* read option characters from input file */
  Char ch;

  while (!(eoln(infile))) {
    ch = gettc(infile);
    uppercase(&ch);
    if (strchr(options, ch) != NULL)
     (* extranum)++;
    else if (!(ch == ' ' || ch == '\t')) {
      printf("BAD OPTION CHARACTER: %c\n", ch);
      exxit(-1);
    }
  }
  scan_eoln(infile);
}  /* readoptions */


void matchoptions(Char *ch, const char *options)
{  /* match option characters to those in auxiliary options line */

  *ch = gettc(infile);
  uppercase(ch);
  if (strchr(options, *ch) == NULL) {
    printf("ERROR: Incorrect auxiliary options line");
    printf(" which starts with %c\n", *ch);
    exxit(-1);
  }
}  /* matchoptions */


void inputweightsold(long chars, steptr weight, boolean *weights)
{
  Char ch;
  int i;

  for (i = 1; i < nmlngth ; i++)
    getc(infile);
 
  for (i = 0; i < chars; i++) {
    do {
      if (eoln(infile)) 
        scan_eoln(infile);
      ch = gettc(infile);
      if (ch == '\n')
        ch = ' ';
    } while (ch == ' ');
    weight[i] = 1;
    if (isdigit(ch))
      weight[i] = ch - '0';
    else if (isalpha(ch)) {
      uppercase(&ch);
      weight[i] = ch - 'A' + 10;
    } else {
      printf("\n\nERROR: Bad weight character: %c\n\n", ch);
      exxit(-1);
    }
  }
  scan_eoln(infile);
  *weights = true;
} /*inputweightsold*/


void inputweights(long chars, steptr weight, boolean *weights)
{
  /* input the character weights, 0-9 and A-Z for weights 0 - 35 */
  Char ch;
  long i;

  for (i = 0; i < chars; i++) {
    do {
      if (eoln(weightfile)) 
        scan_eoln(weightfile);
      ch = gettc(weightfile);
      if (ch == '\n')
        ch = ' ';
    } while (ch == ' ');
    weight[i] = 1;
    if (isdigit(ch))
      weight[i] = ch - '0';
    else if (isalpha(ch)) {
      uppercase(&ch);
      weight[i] = ch - 'A' + 10;
    } else {
      printf("\n\nERROR: Bad weight character: %c\n\n", ch);
      exxit(-1);
    }
  }
  scan_eoln(weightfile);
  *weights = true;
}  /* inputweights */


void inputweights2(long a, long b, long *weightsum,
        steptr weight, boolean *weights, const char *prog)
{
  /* input the character weights, 0 or 1 */
  Char ch;
  long i;

  *weightsum = 0;
  for (i = a; i < b; i++) {
    do {
      if (eoln(weightfile))
        scan_eoln(weightfile);
      ch = gettc(weightfile);
    } while (ch == ' ');
    weight[i] = 1;
    if (ch == '0' || ch == '1')
      weight[i] = ch - '0';
    else {
      printf("\n\nERROR: Bad weight character: %c -- ", ch);
      printf("weights in %s must be 0 or 1\n", prog);
      exxit(-1);
    }
    *weightsum += weight[i];
  }
  *weights = true;
  scan_eoln(weightfile);
}  /* inputweights2 */


void printweights(FILE *filename, long inc, long chars,
        steptr weight, const char *letters)
{
  /* print out the weights of sites */
  long i, j;
  boolean letterweights;

  letterweights = false;
  for (i = 0; i < chars; i++)
    if (weight[i] > 9)
      letterweights = true;
  fprintf(filename, "\n    %s are weighted as follows:",letters);
  if (letterweights)
    fprintf(filename, " (A = 10, B = 11, etc.)\n");
  else
    putc('\n', filename);
  for (i = 0; i < chars; i++) {
    if (i % 60 == 0) {
      putc('\n', filename);
      for (j = 1; j <= nmlngth + 3; j++)
        putc(' ', filename);
    }
    if (weight[i+inc] < 10)
      fprintf(filename, "%ld", weight[i + inc]);
    else
      fprintf(filename, "%c", 'A'-10+(int)weight[i + inc]);
    if ((i+1) % 5 == 0 && (i+1) % 60 != 0)
      putc(' ', filename);
  }
  fprintf(filename, "\n\n");
}  /* printweights */


void inputcategs(long a, long b, steptr category, long categs,const char *prog)
{
  /* input the categories, 1-9 */
  Char ch;
  long i;

  for (i = a; i < b; i++) {
    do {
      if (eoln(catfile)) 
        scan_eoln(catfile);
      ch = gettc(catfile);
    } while (ch == ' ');
    if ((ch >= '1') && (ch <= ('0'+categs)))
      category[i] = ch - '0';
    else {
      printf("\n\nERROR: Bad category character: %c", ch);
      printf(" -- categories in %s are currently 1-%ld\n", prog, categs);
      exxit(-1);
    }
  }
  scan_eoln(catfile);
}  /* inputcategs */


void printcategs(FILE *filename, long chars, steptr category,
                 const char *letters)
{
  /* print out the sitewise categories */
  long i, j;

  fprintf(filename, "\n    %s are:\n",letters);
  for (i = 0; i < chars; i++) {
    if (i % 60 == 0) {
      putc('\n', filename);
      for (j = 1; j <= nmlngth + 3; j++)
        putc(' ', filename);
    }
    fprintf(filename, "%ld", category[i]);
    if ((i+1) % 10 == 0 && (i+1) % 60 != 0)
      putc(' ', filename);
  }
  fprintf(filename, "\n\n");
}  /* printcategs */


void inputfactors(long chars, Char *factor, boolean *factors)
{
  /* reads the factor symbols */
  long i;

  for (i = 0; i < (chars); i++) {
    if (eoln(factfile)) 
      scan_eoln(factfile);
    factor[i] = gettc(factfile);
    if (factor[i] == '\n')
      factor[i] = ' ';
  }
  scan_eoln(factfile);
  *factors = true;
}  /* inputfactors */


void printfactors(FILE *filename, long chars, Char *factor, const char *letters)
{
  /* print out list of factor symbols */
  long i;

  fprintf(filename, "Factors%s:\n\n", letters);
  for (i = 1; i <= nmlngth - 5; i++)
    putc(' ', filename);
  for (i = 1; i <= (chars); i++) {
    newline(filename, i, 55, nmlngth + 3);
    putc(factor[i - 1], filename);
    if (i % 5 == 0)
      putc(' ', filename);
  }
  putc('\n', filename);
}  /* printfactors */


void headings(long chars, const char *letters1, const char *letters2)
{
  long i, j;

  putc('\n', outfile);
  j = nmlngth + (chars + (chars - 1) / 10) / 2 - 5;
  if (j < nmlngth - 1)
    j = nmlngth - 1;
  if (j > 37)
    j = 37;
  fprintf(outfile, "Name");
  for (i = 1; i <= j; i++)
    putc(' ', outfile);
  fprintf(outfile, "%s\n", letters1);
  fprintf(outfile, "----");
  for (i = 1; i <= j; i++)
    putc(' ', outfile);
  fprintf(outfile, "%s\n\n", letters2);
}  /* headings */


void initname(long i)
{
  /* read in species name */
  long j;

  for (j = 0; j < nmlngth; j++) {
    if (eoff(infile) | eoln(infile)){
      printf("\n\nERROR: end-of-line or end-of-file");
      printf(" in the middle of species name for species %ld\n\n", i+1);
      exxit(-1);
    }
    nayme[i][j] = gettc(infile);
    if ((nayme[i][j] == '(') || (nayme[i][j] == ')') || (nayme[i][j] == ':')
        || (nayme[i][j] == ',') || (nayme[i][j] == ';') || (nayme[i][j] == '[')
        || (nayme[i][j] == ']')) {
      printf("\nERROR: Species name may not contain characters ( ) : ; , [ ] \n");
      printf("       In name of species number %ld there is character %c\n\n",
              i+1, nayme[i][j]);
      exxit(-1);
    }
  }
} /* initname */


void findtree(boolean *found,long *pos,long nextree,long *place,bestelm *bestrees)
{
  /* finds tree given by array place in array bestrees by binary search */
  /* used by dnacomp, dnapars, dollop, mix, & protpars */
  long i, lower, upper;
  boolean below, done;

  below = false;
  lower = 1;
  upper = nextree - 1;
  (*found) = false;
  while (!(*found) && lower <= upper) {
    (*pos) = (lower + upper) / 2;
    i = 3;
    done = false;
    while (!done) {
      done = (i > spp);
      if (!done)
        done = (place[i - 1] != bestrees[(*pos) - 1].btree[i - 1]);
      if (!done)
        i++;
    }
    (*found) = (i > spp);
    if (*found)
      break;
    below = (place[i - 1] <  bestrees[(*pos )- 1].btree[i - 1]);
    if (below)
      upper = (*pos) - 1;
    else
      lower = (*pos) + 1;
  }
  if (!(*found) && !below)
    (*pos)++;
}  /* findtree */


void addtree(long pos,long *nextree,boolean collapse,long *place,bestelm *bestrees)
{
  /* puts tree from array place in its proper position in array bestrees */
  /* used by dnacomp, dnapars, dollop, mix, & protpars */
  long i;

  for (i = *nextree - 1; i >= pos; i--){
    memcpy(bestrees[i].btree, bestrees[i - 1].btree, spp * sizeof(long));
    bestrees[i].gloreange = bestrees[i - 1].gloreange;
    bestrees[i - 1].gloreange = false;
    bestrees[i].locreange = bestrees[i - 1].locreange;
    bestrees[i - 1].locreange = false;
    bestrees[i].collapse = bestrees[i - 1].collapse;
  }
  for (i = 0; i < spp; i++)
    bestrees[pos - 1].btree[i] = place[i];
    bestrees[pos - 1].collapse = collapse;
  (*nextree)++;
}  /* addtree */


long findunrearranged(bestelm *bestrees, long nextree, boolean glob)
{
  /* finds bestree with either global or local field false */
  long i;

  if (glob) {
    for (i = 0; i < nextree - 1; i++)
      if (!bestrees[i].gloreange)
        return i;
  } else {
    for (i = 0; i < nextree - 1; i++)
      if (!bestrees[i].locreange)
        return i;
  }
  return -1;
} /* findunrearranged */


boolean torearrange(bestelm *bestrees, long nextree)
{ /* sees if any best tree is yet to be rearranged */

  if (findunrearranged(bestrees, nextree, true) >= 0)
    return true;
  else if (findunrearranged(bestrees, nextree, false) >= 0)
    return true;
  else
    return false;
} /* torearrange */


void reducebestrees(bestelm *bestrees, long *nextree)
{
  /* finds best trees with collapsible branches and deletes them */
  long i, j;

  i = 0;
  j = *nextree - 2;
  do {
    while (!bestrees[i].collapse && i < *nextree - 1) i++;
    while (bestrees[j].collapse && j >= 0) j--;
    if (i < j) {
      memcpy(bestrees[i].btree, bestrees[j].btree, spp * sizeof(long));
      bestrees[i].gloreange = bestrees[j].gloreange;
      bestrees[i].locreange = bestrees[j].locreange;
      bestrees[i].collapse = false;
      bestrees[j].collapse = true;
    }
  } while (i < j);
  *nextree = i + 1;
} /* reducebestrees */


void shellsort(double *a, long *b, long n)
{ /* Shell sort keeping a, b in same order */
  /* used by dnapenny, dolpenny, & penny */
  long gap, i, j, itemp;
  double rtemp;

  gap = n / 2;
  while (gap > 0) {
    for (i = gap + 1; i <= n; i++) {
      j = i - gap;
      while (j > 0) {
     if (a[j - 1] > a[j + gap - 1]) {
       rtemp = a[j - 1];
       a[j - 1] = a[j + gap - 1];
       a[j + gap - 1] = rtemp;
       itemp = b[j - 1];
       b[j - 1] = b[j + gap - 1];
       b[j + gap - 1] = itemp;
     }
     j -= gap;
      }
    }
    gap /= 2;
  }
}  /* shellsort */


void getch(Char *c, long *parens, FILE *treefile)
{ /* get next nonblank character */

  do {
    if (eoln(treefile)) 
      scan_eoln(treefile);
    (*c) = gettc(treefile);

    if ((*c) == '\n' || (*c) == '\t')
      (*c) = ' ';
  } while ( *c == ' ' && !eoff(treefile) );
  if ((*c) == '(')
    (*parens)++;
  if ((*c) == ')')
    (*parens)--;
}  /* getch */


void getch2(Char *c, long *parens)
{ /* get next nonblank character */
  do {
    if (eoln(intree)) 
      scan_eoln(intree);
    *c = gettc(intree);
    if (*c == '\n' || *c == '\t')
      *c = ' ';
  } while (!(*c != ' ' || eoff(intree)));
  if (*c == '(')
   (*parens)++;
  if (*c == ')')
    (*parens)--;
}  /* getch2 */


void findch(Char c, Char *ch, long which)
{ /* scan forward until find character c */
  boolean done;
  long dummy_parens;
  done = false;
  while (!done) {
    if (c == ',') {
      if (*ch == '(' || *ch == ')' || *ch == ';') {
        printf(
      "\n\nERROR in user tree %ld: unmatched parenthesis or missing comma\n\n",
          which);
        exxit(-1);
      } else if (*ch == ',')
        done = true;
    } else if (c == ')') {
      if (*ch == '(' || *ch == ',' || *ch == ';') {
        printf("\n\nERROR in user tree %ld: ", which);
        printf("unmatched parenthesis or non-bifurcated node\n\n");
        exxit(-1);
      } else {
        if (*ch == ')')
          done = true;
      }
    } else if (c == ';') {
      if (*ch != ';') {
        printf("\n\nERROR in user tree %ld: ", which);
        printf("unmatched parenthesis or missing semicolon\n\n");
        exxit(-1);
      } else
        done = true;
    }
    if (*ch != ')' && done)
      continue;
   getch(ch, &dummy_parens, intree);
  }
}  /* findch */


void findch2(Char c, long *lparens, long *rparens, Char *ch)
{ /* skip forward in user tree until find character c */
  boolean done;
  long dummy_parens;
  done = false;
  while (!done) {
    if (c == ',') {
      if (*ch == '(' || *ch == ')' || *ch == ':' || *ch == ';') {
        printf("\n\nERROR in user tree: ");
        printf("unmatched parenthesis, missing comma");
        printf(" or non-trifurcated base\n\n");
     exxit(-1);
      } else if (*ch == ',')
        done = true;
    } else if (c == ')') {
      if (*ch == '(' || *ch == ',' || *ch == ':' || *ch == ';') {
        printf(
   "\n\nERROR in user tree: unmatched parenthesis or non-bifurcated node\n\n");
        exxit(-1);
      } else if (*ch == ')') {
        (*rparens)++;
        if ((*lparens) > 0 && (*lparens) == (*rparens)) {
          if ((*lparens) == spp - 2) {
           getch(ch, &dummy_parens, intree);
            if (*ch != ';') {
              printf( "\n\nERROR in user tree: ");
              printf("unmatched parenthesis or missing semicolon\n\n");
              exxit(-1);
            }
          }
        }
     done = true;
      }
    }
    if (*ch != ')' && done)
      continue;
    if (*ch == ')')
     getch(ch, &dummy_parens, intree);
  }
}  /* findch2 */


void processlength(double *valyew, double *divisor, Char *ch, 
        boolean *minusread, FILE *treefile, long *parens)
{ /* read a branch length from a treefile */
  long digit, ordzero;
  boolean pointread;

  ordzero = '0';
  *minusread = false;
  pointread = false;
  *valyew = 0.0;
  *divisor = 1.0;
  getch(ch, parens, treefile);
  digit = (long)(*ch - ordzero);
  while ( ((digit <= 9) && (digit >= 0)) || *ch == '.' || *ch == '-') {
    if (*ch == '.' )
      pointread = true;
    else if (*ch == '-' )
      *minusread = true;
    else {
      *valyew = *valyew * 10.0 + digit;
      if (pointread)
        *divisor *= 10.0;
    }
    getch(ch, parens, treefile);
    digit = (long)(*ch - ordzero);
  }
  if (*minusread)
    *valyew = -(*valyew);
}  /* processlength */


void writename(long start, long n, long *enterorder)
{ /* write species name and number in entry order */
  long i, j;

  for (i = start; i < start+n; i++) {
    printf(" %3ld. ", i+1);
    for (j = 0; j < nmlngth; j++)
      putchar(nayme[enterorder[i] - 1][j]);
    putchar('\n');
    fflush(stdout);
  }
}  /* writename */


void memerror()
{
  printf("Error allocating memory\n");
  exxit(-1);
}  /* memerror */


void odd_malloc(long x)
{ /* error message if attempt to malloc too little or too much memory */
  printf ("ERROR: a function asked for an inappropriate amount of memory:");
  printf ("  %ld bytes\n", x);
  printf ("       This can mean one of two things:\n");
  printf ("       1.  The input file is incorrect");
  printf (" (perhaps it was not saved as Text Only),\n");
  printf ("       2.  There is a bug in the program.\n");
  printf ("       Please check your input file carefully.\n");
  printf ("       If it seems to be a bug, please mail joe@gs.washington.edu\n");
  printf ("       with the name of the program, your computer system type,\n");
  printf ("       a full description of the problem, and with the input data file.\n");
  printf ("       (which should be in the body of the message, not as an Attachment).\n");

  /* abort() can be used to crash */
  
  exxit(-1);
}


MALLOCRETURN *mymalloc(long x)
{ /* wrapper for malloc, allowing error message if too little, too much */
  MALLOCRETURN *new_block;

  if ((x <= 0) ||
      (x > TOO_MUCH_MEMORY))
    odd_malloc(x);

  new_block = (MALLOCRETURN *)calloc(1,x);

  if (!new_block) {
    memerror();
    return (MALLOCRETURN *) new_block;
  } else
    return (MALLOCRETURN *) new_block;
} /* mymalloc */


void gnu(node **grbg, node **p)
{ /* this and the following are do-it-yourself garbage collectors.
     Make a new node or pull one off the garbage list */

  if (*grbg != NULL) {
    *p = *grbg;
    *grbg = (*grbg)->next;
  } else
    *p = (node *)Malloc(sizeof(node));

  (*p)->back       = NULL;
  (*p)->next       = NULL;
  (*p)->tip        = false;
  (*p)->times_in_tree = 0.0;
  (*p)->r          = 0.0;
  (*p)->theta      = 0.0;
  (*p)->x          = NULL;
  (*p)->protx           = NULL;        /* for the sake of proml     */
}  /* gnu */


void chuck(node **grbg, node *p)
{ /* collect garbage on p -- put it on front of garbage list */
  p->back = NULL;
  p->next = *grbg;
  *grbg = p;
}  /* chuck */


void zeronumnuc(node *p, long endsite)
{
  long i,j;

  for (i = 0; i < endsite; i++)
    for (j = (long)A; j <= (long)O; j++)
      p->numnuc[i][j] = 0;
} /* zeronumnuc */


void zerodiscnumnuc(node *p, long endsite)
{
  long i,j;

  for (i = 0; i < endsite; i++)
    for (j = (long)zero; j <= (long)seven; j++)
      p->discnumnuc[i][j] = 0;
} /* zerodiscnumnuc */


void allocnontip(node *p, long *zeros, long endsite)
{ /* allocate an interior node */
  /* used by dnacomp, dnapars, & dnapenny */

  p->numsteps = (steptr)Malloc(endsite*sizeof(long));
  p->oldnumsteps = (steptr)Malloc(endsite*sizeof(long));
  p->base = (baseptr)Malloc(endsite*sizeof(long));
  p->oldbase = (baseptr)Malloc(endsite*sizeof(long));
  p->numnuc = (nucarray *)Malloc(endsite*sizeof(nucarray));
  memcpy(p->base, zeros, endsite*sizeof(long));
  memcpy(p->numsteps, zeros, endsite*sizeof(long));
  memcpy(p->oldbase, zeros, endsite*sizeof(long));
  memcpy(p->oldnumsteps, zeros, endsite*sizeof(long));
  zeronumnuc(p, endsite);
}  /* allocnontip */


void allocdiscnontip(node *p, long *zeros, unsigned char *zeros2, long endsite)
{ /* allocate an interior node */
  /* used by pars */

  p->numsteps = (steptr)Malloc(endsite*sizeof(long));
  p->oldnumsteps = (steptr)Malloc(endsite*sizeof(long));
  p->discbase = (discbaseptr)Malloc(endsite*sizeof(unsigned char));
  p->olddiscbase = (discbaseptr)Malloc(endsite*sizeof(unsigned char));
  p->discnumnuc = (discnucarray *)Malloc(endsite*sizeof(discnucarray));
  memcpy(p->discbase, zeros2, endsite*sizeof(unsigned char));
  memcpy(p->numsteps, zeros, endsite*sizeof(long));
  memcpy(p->olddiscbase, zeros2, endsite*sizeof(unsigned char));
  memcpy(p->oldnumsteps, zeros, endsite*sizeof(long));
  zerodiscnumnuc(p, endsite);
}  /* allocdiscnontip */


void allocnode(node **anode, long *zeros, long endsite)
{ /* allocate a node */
  /* used by dnacomp, dnapars, & dnapenny */

  *anode = (node *)Malloc(sizeof(node));
  allocnontip(*anode, zeros, endsite);
}  /* allocnode */


void allocdiscnode(node **anode, long *zeros, unsigned char *zeros2, 
        long endsite)
{ /* allocate a node */
  /* used by pars */

  *anode = (node *)Malloc(sizeof(node));
  allocdiscnontip(*anode, zeros, zeros2, endsite);
}  /* allocdiscnontip */


void gnutreenode(node **grbg, node **p, long i, long endsite, long *zeros)
{ /* this and the following are do-it-yourself garbage collectors.
     Make a new node or pull one off the garbage list */

  if (*grbg != NULL) {
    *p = *grbg;
    *grbg = (*grbg)->next;
    memcpy((*p)->numsteps, zeros, endsite*sizeof(long));
    memcpy((*p)->oldnumsteps, zeros, endsite*sizeof(long));
    memcpy((*p)->base, zeros, endsite*sizeof(long));
    memcpy((*p)->oldbase, zeros, endsite*sizeof(long));
    zeronumnuc(*p, endsite);
  } else
    allocnode(p, zeros, endsite);
  (*p)->back = NULL;
  (*p)->next = NULL;
  (*p)->tip = false;
  (*p)->visited = false;
  (*p)->index = i;
  (*p)->numdesc = 0;
  (*p)->sumsteps = 0.0;
}  /* gnutreenode */


void gnudisctreenode(node **grbg, node **p, long i,
        long endsite, long *zeros, unsigned char *zeros2)
{ /* this and the following are do-it-yourself garbage collectors.
     Make a new node or pull one off the garbage list */

  if (*grbg != NULL) {
    *p = *grbg;
    *grbg = (*grbg)->next;
    memcpy((*p)->numsteps, zeros, endsite*sizeof(long));
    memcpy((*p)->oldnumsteps, zeros, endsite*sizeof(long));
    memcpy((*p)->discbase, zeros2, endsite*sizeof(unsigned char));
    memcpy((*p)->olddiscbase, zeros2, endsite*sizeof(unsigned char));
    zerodiscnumnuc(*p, endsite);
  } else
    allocdiscnode(p, zeros, zeros2, endsite);
  (*p)->back = NULL;
  (*p)->next = NULL;
  (*p)->tip = false;
  (*p)->visited = false;
  (*p)->index = i;
  (*p)->numdesc = 0;
  (*p)->sumsteps = 0.0;
}  /* gnudisctreenode */


void chucktreenode(node **grbg, node *p)
{ /* collect garbage on p -- put it on front of garbage list */

  p->back = NULL;
  p->next = *grbg;
  *grbg = p;
}  /* chucktreenode */


void setupnode(node *p, long i)
{ /* initialization of node pointers, variables */

  p->next = NULL;
  p->back = NULL;
  p->times_in_tree = (double) i * 1.0;
  p->index = i;
  p->tip = false;
}  /* setupnode */


long count_sibs (node *p)
{ /* Count the number of nodes in a ring, return the total number of */
  /* nodes excluding the one passed into the function (siblings)     */
  node *q;
  long return_int = 0;

  if (p->tip) {
   printf ("Error: the function count_sibs called on a tip.  This is a bug.\n");
    exxit (-1);
  }

  q = p->next;
  while (q != p) {
    if (q == NULL) {
      printf ("Error: a loop of nodes was not closed.\n");
      exxit (-1);
    } else {
      return_int++;
      q = q->next;
    }
  }

  return return_int;
}  /* count_sibs */


void inittrav (node *p)
{ /* traverse to set pointers uninitialized on inserting */
  long i, num_sibs;
  node *sib_ptr;
  
  if (p == NULL)
    return;
  if (p->tip)
    return;
  num_sibs = count_sibs (p);
  sib_ptr  = p;
  for (i=0; i < num_sibs; i++) {
    sib_ptr              = sib_ptr->next;
    sib_ptr->initialized = false;
    inittrav(sib_ptr->back);
  }
} /* inittrav */


void commentskipper(FILE ***intree, long *bracket)
{ /* skip over comment bracket contents in reading tree */
  char c;
  
  c = gettc(**intree);
  
  while (c != ']') {
    
    if(feof(**intree)) {
      printf("\n\nERROR: Unmatched comment brackets\n\n");
      exxit(-1);
    }

    if(c == '[') {
      (*bracket)++;
      commentskipper(intree, bracket);
    }
    c = gettc(**intree);
  }
  (*bracket)--;
}  /* commentskipper */


long countcomma(FILE **treefile, long *comma)
{
  /* Modified by Dan F. 11/10/96 */ 

  /* The next line inserted so this function leaves the file pointing
     to where it found it, not just re-winding it. */
  long orig_position = ftell (*treefile);

  Char c;
  long  lparen = 0;
  long bracket = 0;
  (*comma) = 0;


  for (;;){
    c = getc(*treefile);
    if (feof(*treefile))
      break;
    if (c == ';')
      break;
    if (c == ',')
      (*comma)++;
    if (c == '(')
         lparen++;
    if (c == '[') {
      bracket++;
      commentskipper(&treefile, &bracket);
    }
  }

  /* Don't just rewind, */
  /* rewind (*treefile); */
  /* Re-set to where it pointed when the function was called */

  fseek (*treefile, orig_position, SEEK_SET);

  return lparen + (*comma);
}  /*countcomma*/
/* countcomma rewritten so it passes back both lparen+comma to allocate nodep
   and a pointer to the comma variable.  This allows the tree to know how many
   species exist, and the tips to be placed in the front of the nodep array */


long countsemic(FILE **treefile)
{ /* Used to determine the number of user trees.  Return
     either a: the number of semicolons in the file outside comments
     or b: the first integer in the file */
  Char c;
  long return_val, semic = 0;
  long bracket = 0;
  
  /* Eat all whitespace */
  c = gettc(*treefile);
  while ((c == ' ')  ||
      (c == '\t') ||
      (c == '\n')) {
    c = gettc(*treefile);
  }

  /* Then figure out if the first non-white character is a digit; if
     so, return it */
  if (isdigit (c)) {
    ungetc(c, *treefile);
    fscanf((*treefile), "%ld", &return_val);
  } else {

    /* Loop past all characters, count the number of semicolons
       outside of comments */
    for (;;){
      c = fgetc(*treefile);
      if (feof(*treefile))
     break;
      if (c == ';')
     semic++;
      if (c == '[') {
     bracket++;
     commentskipper(&treefile, &bracket);
      }
    }
    return_val = semic;
  }

  rewind (*treefile);
  return return_val;
}  /* countsemic */


void hookup(node *p, node *q)
{ /* hook together two nodes */

  p->back = q;
  q->back = p;
}  /* hookup */


void link_trees(long local_nextnum, long nodenum, long local_nodenum,
        pointarray nodep)
{
  if(local_nextnum == 0)
    hookup(nodep[nodenum],nodep[local_nodenum]);
  else if(local_nextnum == 1)
      hookup(nodep[nodenum], nodep[local_nodenum]->next);
    else if(local_nextnum == 2)
        hookup(nodep[nodenum],nodep[local_nodenum]->next->next);
      else
        printf("Error in Link_Trees()");
} /* link_trees() */


void allocate_nodep(pointarray *nodep, FILE **treefile, long  *precalc_tips)  
{ /* pre-compute space and allocate memory for nodep */

  long numnodes;      /* returns number commas & (    */
  long numcom = 0;        /* returns number commas */
  
  numnodes = countcomma(treefile, &numcom) + 1;
  *nodep      = (pointarray)Malloc(2*numnodes*sizeof(node *));

  (*precalc_tips) = numcom + 1;        /* this will be used in placing the
                                          tip nodes in the front region of
                                          nodep.  Used for species check?  */
} /* allocate_nodep -plc */


void malloc_pheno (node *p, long endsite, long rcategs)
{ /* Allocate the phenotype arrays; used by dnaml */
  long i;

  p->x  = (phenotype)Malloc(endsite*sizeof(ratelike));
  for (i = 0; i < endsite; i++)
    p->x[i]  = (ratelike)Malloc(rcategs*sizeof(sitelike));
} /* malloc_pheno */


void malloc_ppheno (node *p,long endsite, long rcategs)
{
  /* Allocate the phenotype arrays; used by proml */
  long i;

  p->protx  = (pphenotype)Malloc(endsite*sizeof(pratelike));
  for (i = 0; i < endsite; i++)
    p->protx[i]  = (pratelike)Malloc(rcategs*sizeof(psitelike));
} /* malloc_ppheno */


long take_name_from_tree (Char *ch, Char *str, FILE *treefile)
{
  /* This loop takes in the name from the tree.
     Return the length of the name string.  */

  long name_length = 0;

  do {
    if ((*ch) == '_')
      (*ch) = ' ';
    str[name_length++] = (*ch);
    if (eoln(treefile)) 
      scan_eoln(treefile);
    (*ch) = gettc(treefile);
    if (*ch == '\n')
      *ch = ' ';
  } while ((*ch) != ':' && (*ch) != ',' && (*ch) != ')' &&
        (*ch) != '[' && (*ch) != ';' && name_length <= MAXNCH);
  return name_length;
}  /* take_name_from_tree */


void match_names_to_data (Char *str, pointarray treenode, node **p, long spp)
{
  /* This loop matches names taken from treefile to indexed names in
     the data file */

  boolean found;
  long i, n;

  n = 1;  
  do {
    found = true;
    for (i = 0; i < nmlngth; i++) {
      found = (found && ((str[i] == nayme[n - 1][i]) ||
        (((nayme[n - 1][i] == '_') && (str[i] == ' ')) ||
        ((nayme[n - 1][i] == ' ') && (str[i] == '\0')))));
    }
    
    if (found)
      *p = treenode[n - 1];
    else
      n++;

  } while (!(n > spp || found));
  
  if (n > spp) {
    printf("\n\nERROR: Cannot find species: ");
    for (i = 0; (str[i] != '\0') && (i < MAXNCH); i++)
      putchar(str[i]);
    printf(" in data file\n\n");
    exxit(-1);
  }
}  /* match_names_to_data */


void addelement(node **p, node *q, Char *ch, long *parens, FILE *treefile,
        pointarray treenode, boolean *goteof, boolean *first, pointarray nodep,
        long *nextnode, long *ntips, boolean *haslengths, node **grbg,
        initptr initnode)
{
  /* Recursive procedure adds nodes to user-defined tree
     This is the main (new) tree-reading procedure */

  node *pfirst;
  long i, len = 0, nodei = 0;
  boolean notlast;
  Char str[MAXNCH];
  node *r;

  if ((*ch) == '(') {
    (*nextnode)++;          /* get ready to use new interior node */
    nodei = *nextnode;      /* do what needs to be done at bottom */
    (*initnode)(p, grbg, q, len, nodei, ntips,
                  parens, bottom, treenode, nodep, str, ch, treefile);
    pfirst      = (*p);
    notlast = true;
    while (notlast) {          /* loop through immediate descendants */
      (*initnode)(&(*p)->next, grbg, q,
                   len, nodei, ntips, parens, nonbottom, treenode,
                   nodep, str, ch, treefile);
                       /* ... doing what is done before each */
      r = (*p)->next;
      getch(ch, parens, treefile);      /* look for next character */
      
       /* handle blank names */
      if((*ch) == ',' || (*ch) == ':'){
        ungetc((*ch), treefile);
        *ch = 0;
      } else if((*ch)==')'){
        ungetc((*ch), treefile);
        (*parens)++;
        *ch = 0;
      }
 
      addelement(&(*p)->next->back, (*p)->next, ch, parens, treefile,
        treenode, goteof, first, nodep, nextnode, ntips,
        haslengths, grbg, initnode);        /* call self recursively */

      (*initnode)(&r, grbg, q, len, nodei, ntips,
                    parens, hslength, treenode, nodep, str, ch, treefile);
                           /* do what is done after each about length */
      pfirst->numdesc++;               /* increment number of descendants */
      *p = r;                         /* make r point back to p */

      if ((*ch) == ')') {
        notlast = false;
        do {
          getch(ch, parens, treefile);
        } while ((*ch) != ',' && (*ch) != ')' &&
           (*ch) != '[' && (*ch) != ';' && (*ch) != ':');
      }
    }
    
    (*p)->next = pfirst;
    (*p)       = pfirst;

  } else if ((*ch) != ')') {       /* if it's a species name */
    for (i = 0; i < MAXNCH; i++)   /* fill string with nulls */
      str[i] = '\0';

    len = take_name_from_tree (ch, str, treefile);  /* get the name */

    if ((*ch) == ')')
      (*parens)--;         /* decrement count of open parentheses */
    (*initnode)(p, grbg, q, len, nodei, ntips,
                  parens, tip, treenode, nodep, str, ch, treefile);
                              /* do what needs to be done at a tip */
  } else
    getch(ch, parens, treefile);
  if (q != NULL)
    hookup(q, (*p));                    /* now hook up */
  (*initnode)(p, grbg, q, len, nodei, ntips, 
                parens, iter, treenode, nodep, str, ch, treefile);
                              /* do what needs to be done to variable iter */
  if ((*ch) == ':')
    (*initnode)(p, grbg, q, len, nodei, ntips, 
                  parens, length, treenode, nodep, str, ch, treefile);
                                   /* do what needs to be done with length */
  else if ((*ch) != ';' && (*ch) != '[')
    (*initnode)(p, grbg, q, len, nodei, ntips, 
                  parens, hsnolength, treenode, nodep, str, ch, treefile);
                             /* ... or what needs to be done when no length */
  if ((*ch) == '[')
    (*initnode)(p, grbg, q, len, nodei, ntips,
                  parens, treewt, treenode, nodep, str, ch, treefile);
                             /* ... for processing a tree weight */
  else if ((*ch) == ';')     /* ... and at end of tree */
    (*initnode)(p, grbg, q, len, nodei, ntips,
                  parens, unittrwt, treenode, nodep, str, ch, treefile);
}  /* addelement */


void treeread (FILE *treefile, node **root, pointarray treenode,
        boolean *goteof, boolean *first, pointarray nodep, 
        long *nextnode, boolean *haslengths, node **grbg, initptr initnode)
{
  /* read in user-defined tree and set it up */
  char  ch;
  long parens = 0;
  long ntips = 0;
  
  (*goteof) = false;
  (*nextnode) = spp;

  /* eat blank lines */
  while (eoln(treefile) && !eoff(treefile)) 
    scan_eoln(treefile);

  if (eoff(treefile)) {
    (*goteof) = true;
    return;
  } 

  getch(&ch, &parens, treefile);

  while (ch != '(') {
    /* Eat everything in the file (i.e. digits, tabs) until you
       encounter an open-paren */
    getch(&ch, &parens, treefile);
  }
  (*haslengths) = true; 
  addelement(root, NULL, &ch, &parens, treefile,
         treenode, goteof, first, nodep, nextnode, &ntips,
         haslengths, grbg, initnode);
  
  /* Eat blank lines and end of current line*/ 
  do {
    scan_eoln(treefile);
  }
  while (eoln(treefile) && !eoff(treefile));
  
  (*first) = false;
  if (parens != 0) {
    printf("\n\nERROR in tree file: unmatched parentheses\n\n");
    exxit(-1);
  }
}  /* treeread */


void addelement2(node *q, Char *ch, long *parens, FILE *treefile,
        pointarray treenode, boolean lngths, double *trweight, boolean *goteof,
        long *nextnode, long *ntips, long no_species, boolean *haslengths)
{
  /* recursive procedure adds nodes to user-defined tree
     -- old-style bifurcating-only version */
  node *pfirst = NULL, *p;
  long i, len, current_loop_index;
  boolean notlast, minusread;
  Char str[MAXNCH];
  double valyew, divisor;

  if ((*ch) == '(') {

    current_loop_index = (*nextnode) + spp;
    (*nextnode)++;

    /* This is an assignment of an interior node */
    p = treenode[current_loop_index];
    pfirst = p;
    notlast = true;
    while (notlast) {
      /* This while loop goes through a circle (triad for
      bifurcations) of nodes */
      p = p->next;
      /* added to ensure that non base nodes in loops have indices */
      p->index = current_loop_index + 1;
      
      getch(ch, parens, treefile);
      
      addelement2(p, ch, parens, treefile, treenode, lngths, trweight,
        goteof, nextnode, ntips, no_species, haslengths);

      if ((*ch) == ')') {
        notlast = false;
        do {
          getch(ch, parens, treefile);
        } while ((*ch) != ',' && (*ch) != ')' &&
           (*ch) != '[' && (*ch) != ';' && (*ch) != ':');
      }
    }

  } else if ((*ch) != ')') {
    for (i = 0; i < MAXNCH; i++) 
      str[i] = '\0';
    len = take_name_from_tree (ch, str, treefile);
    match_names_to_data (str, treenode, &p, spp);
    pfirst = p;
    if ((*ch) == ')')
      (*parens)--;
    (*ntips)++;
    strncpy (p->nayme, str, len);
  } else
    getch(ch, parens, treefile);
  
  if ((*ch) == '[') {    /* getting tree weight from last comment field */
    if (!eoln(treefile)) {
      fscanf(treefile, "%lf", trweight);
      getch(ch, parens, treefile);
      if (*ch != ']') {
        printf("\n\nERROR: Missing right square bracket\n\n");
        exxit(-1);
      }
      else {
        getch(ch, parens, treefile);
        if (*ch != ';') {
          printf("\n\nERROR: Missing semicolon after square brackets\n\n");
          exxit(-1);
        }
      }
    }
  }
  else if ((*ch) == ';') {
    (*trweight) = 1.0 ;
    if (!eoln(treefile))
      printf("WARNING: tree weight set to 1.0\n");
  }
  else
    (*haslengths) = ((*haslengths) && q == NULL);
  
  if (q != NULL)
    hookup(q, pfirst);
  if (q != NULL) {
    if (q->branchnum < pfirst->branchnum)
      pfirst->branchnum = q->branchnum;
    else
      q->branchnum = pfirst->branchnum;
  }

  if ((*ch) == ':') {
    processlength(&valyew, &divisor, ch,
       &minusread, treefile, parens);
    if (q != NULL) {
      if (!minusread)
        q->oldlen = valyew / divisor;
      else
        q->oldlen = 0.0;
      if (lngths) {
        q->v = valyew / divisor;
        q->back->v = q->v;
        q->iter = false;
        q->back->iter = false;
        q->back->iter = false;
    }
    }
  }
  
}  /* addelement2 */


void treeread2 (FILE *treefile, node **root, pointarray treenode,
        boolean lngths, double *trweight, boolean *goteof,
        boolean *haslengths, long *no_species)
{
  /* read in user-defined tree and set it up
     -- old-style bifurcating-only version */
  char  ch;
  long parens = 0;
  long ntips = 0;
  long nextnode;
  
  (*goteof) = false;
  nextnode = 0;

  /* Eats all blank lines at start of file */
  while (eoln(treefile) && !eoff(treefile)) 
    scan_eoln(treefile);

  if (eoff(treefile)) {
    (*goteof) = true;
    return;
  } 

  getch(&ch, &parens, treefile);

  while (ch != '(') {
    /* Eat everything in the file (i.e. digits, tabs) until you
       encounter an open-paren */
    getch(&ch, &parens, treefile);
  }

  addelement2(NULL, &ch, &parens, treefile, treenode, lngths, trweight,
          goteof, &nextnode, &ntips, (*no_species), haslengths);
  (*root) = treenode[*no_species];

  /*eat blank lines */
  while (eoln(treefile) && !eoff(treefile)) 
    scan_eoln(treefile);

  (*root)->oldlen = 0.0;

  if (parens != 0) {
    printf("\n\nERROR in tree file:  unmatched parentheses\n\n");
    exxit(-1);
  }
}  /* treeread2 */


void exxit(int exitcode)
{
#ifdef WIN32
  if (exitcode == 0)
#endif
    exit (exitcode);
#ifdef WIN32
  else {
    puts ("Hit Enter or Return to close program.");
    puts("  You may have to hit Enter or Return twice.");
    getchar ();
    getchar ();
    phyRestoreConsoleAttributes();
    exit (exitcode);
  }
#endif
} /* exxit */


char gettc(FILE* file) 
{ /* catch eof's so that other functions not expecting an eof
   * won't have to worry about it */
  int ch;

  ch=getc(file);

  if (ch == EOF )   {
    puts("Unexpected End of File");
    exxit(-1);
  }

  if ( ch == '\r' ) {
    ch = getc(file);
    if ( ch != '\n' )
      ungetc(ch,file);
    ch = '\n';
  }
  return ch;
} /* gettc */

void unroot(tree *t, long nonodes) 
{
  /* used by fitch, restml and contml */
  if (t->start->back == NULL) { 
    if (t->start->next->back->tip)
      t->start = t->start->next->next->back;
    else  t->start = t->start->next->back;
  }
  if (t->start->next->back == NULL) {
    if (t->start->back->tip)
      t->start = t->start->next->next->back;
    else t->start = t->start->back;
  }
  if (t->start->next->next->back == NULL)  {
    if (t->start->back->tip)
      t->start = t->start->next->back;
    else t->start = t->start->back;
  }
    

  unroot_r(t->start,t->nodep,nonodes);
  unroot_r(t->start->back, t->nodep, nonodes);
}

void unroot_here(node* root, node** nodep, long nonodes)
{
  node* tmpnode;
  double newl;
  /* used by unroot */
  /* assumes bifurcation this is ok in the programs that use it */
 
  
  newl = root->next->oldlen + root->next->next->oldlen;
  root->next->back->oldlen = newl;
  root->next->next->back->oldlen = newl;

  newl = root->next->v + root->next->next->v;
  root->next->back->v = newl;
  root->next->next->back->v = newl;

  root->next->back->back=root->next->next->back;
  root->next->next->back->back = root->next->back;
  
  while ( root->index != nonodes ) {
    tmpnode = nodep[ root->index ];
    nodep[root->index] = root;
    root->index++;
    root->next->index++;
    root->next->next->index++;
    nodep[root->index - 2] = tmpnode;
    tmpnode->index--;
    tmpnode->next->index--;
    tmpnode->next->next->index--;
  }
}

void unroot_r(node* p, node** nodep, long nonodes) 
{
  /* used by unroot */
  node *q;
  
  if ( p->tip) return;

  q = p->next;
  while ( q != p ) {
    if (q->back == NULL)
      unroot_here(q,nodep,nonodes);
    else unroot_r(q->back,nodep,nonodes);
    q = q->next;
  }
}

void clear_connections(tree *t, long nonodes) 
{
  long i;
  for ( i = 0 ; i < nonodes ; i++) {
    if ( i > spp) {
      t->nodep[i]->next->back = NULL;
      t->nodep[i]->next->v = 0;
      t->nodep[i]->next->next->back = NULL;
      t->nodep[i]->next->next->v = 0;
    }
    t->nodep[i]->back = NULL;
    t->nodep[i]->v = 0;
  }
}

#ifdef WIN32
void phySaveConsoleAttributes()
{
  GetConsoleScreenBufferInfo( hConsoleOutput, &savecsbi );
} /* PhySaveConsoleAttributes */


void phySetConsoleAttributes()
{
  hConsoleOutput = GetStdHandle(STD_OUTPUT_HANDLE);

  phySaveConsoleAttributes();

  SetConsoleTextAttribute(hConsoleOutput, 
    BACKGROUND_GREEN | BACKGROUND_BLUE | BACKGROUND_INTENSITY);
} /* phySetConsoleAttributes */


void phyRestoreConsoleAttributes()
{
  COORD coordScreen = { 0, 0 };
  DWORD cCharsWritten;
  DWORD dwConSize; 

  dwConSize = savecsbi.dwSize.X * savecsbi.dwSize.Y;

  SetConsoleTextAttribute(hConsoleOutput, savecsbi.wAttributes);

  FillConsoleOutputAttribute( hConsoleOutput, savecsbi.wAttributes,
         dwConSize, coordScreen, &cCharsWritten );
} /* phyRestoreConsoleAttributes */


void phyFillScreenColor()
{
  COORD coordScreen = { 0, 0 };
  DWORD cCharsWritten;
  CONSOLE_SCREEN_BUFFER_INFO csbi; /* to get buffer info */ 
  DWORD dwConSize; 

  GetConsoleScreenBufferInfo( hConsoleOutput, &csbi );
  dwConSize = csbi.dwSize.X * csbi.dwSize.Y;

  FillConsoleOutputAttribute( hConsoleOutput, csbi.wAttributes,
         dwConSize, coordScreen, &cCharsWritten );
} /* PhyFillScreenColor */


void phyClearScreen()
{
   COORD coordScreen = { 0, 0 };    /* here's where we'll home the
                                       cursor */ 
   DWORD cCharsWritten;
   CONSOLE_SCREEN_BUFFER_INFO csbi; /* to get buffer info */ 
   DWORD dwConSize;                 /* number of character cells in
                                       the current buffer */ 

   /* get the number of character cells in the current buffer */ 

   GetConsoleScreenBufferInfo( hConsoleOutput, &csbi );
   dwConSize = csbi.dwSize.X * csbi.dwSize.Y;

   /* fill the entire screen with blanks */ 

   FillConsoleOutputCharacter( hConsoleOutput, (TCHAR) ' ',
      dwConSize, coordScreen, &cCharsWritten );

   /* get the current text attribute */ 

   GetConsoleScreenBufferInfo( hConsoleOutput, &csbi );

   /* now set the buffer's attributes accordingly */ 

   FillConsoleOutputAttribute( hConsoleOutput, csbi.wAttributes,
      dwConSize, coordScreen, &cCharsWritten );

   /* put the cursor at (0, 0) */ 

   SetConsoleCursorPosition( hConsoleOutput, coordScreen );
   return;
} /* PhyClearScreen */
#endif

