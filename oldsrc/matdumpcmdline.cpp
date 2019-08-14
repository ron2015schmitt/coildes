/************************************************************************* 
* 
*   File Name    :  cmdline.cpp
*   Platform     :  Red Hat LINUX 
*   Author       :  Ron Schmitt
*   Date         :  16 May 2003
* 
*
*   SYNOPSIS     
*     This code parses the input from the unix command line
*
**************************************************************************/


#include <cstdio>
#include <cstring>
#include <cstdlib>


#include "matdumpcmdline.hpp"

bool printdir;
char *vnames[255];
int vn;
char *mnames[255];
int mn;

const int NUM_FILES = 1;

const int NUM_OPTIONS = 3;

const int MAX_PARMS_PER_OPTION = 1;


char *myname = "";
// Actual text for each command line option
char *opts[NUM_OPTIONS] = 
  {
    "-d", "-v", "-m"
  };


// enumerated types for option
typedef enum 
  {
    opt_dir, opt_vec, opt_mat
  } t_opts;


// number of parameters for each option
int opts_parms[NUM_OPTIONS] = 
  {
    0,1,1
  };


// whether or not this is a required option (yes this is an oxymoron)
int opts_required[NUM_OPTIONS] = 
  {
    0,0,0
  };


// text that describes each option
char *opts_help[NUM_OPTIONS] = 
  {
    "print out the names of all variables stored in given file",
    "print out the contents of the vector <vname>",
    "print out the contents of the matrix <mname>"
  };


// text that for labeling the parameters of each option
char *opts_parm_text[NUM_OPTIONS][MAX_PARMS_PER_OPTION] = 
  {
    {""},
    {"<vname>"},
    {"<mname>"},
  };



// Synopsis for the command
char synopsis[] = {
  "\n\tThis program prints out information from .mat matlab files:\n"
  "\t\t example:  matdump -v x -m A results.mat\n"
  "\t\t \n"
};



/*************************************************************************
 *
 *   FUNCTION     : opt_func
 *   INPUTS       :  
 *     opt - option type
 *     parms - array of parameter strings that were input with the given option
 *   OUTPUTS      :  
 *     none
 *   I/O          :
 *    none
 *
 *   SYNOPSIS
 *
 *   This function is called by the command line parsing engine.
 *   It is called once for each option that is parsed. This function
 *   then performs the appropriate action for the option.
 *
 **************************************************************************/

void opts_func(t_opts opt, char *parms[])
{

#ifdef DEBUG
  {
    int debug;
    printf("found: %s",opts[opt]);
    for(debug=0; debug < opts_parms[opt]; debug++)
      printf(" '%s'",parms[debug]);
    printf("\n");
  }
#endif

  switch(opt) {
  case opt_dir:
    printdir = true;
    break;
  case opt_vec:
    vnames[vn++] = parms[0];
    break;
  case opt_mat:
    mnames[mn++] = parms[0];
    break;


  } // switch(opt)

  
}




/*-------------------------COMMAND PARSING ENGINE-----------------------*/
/************************************************************************
 This code is generic, i.e. it contains no specific code to this project
 and thus allows easier portability

 This code can parse out all your options and filenames from the command line.
 Options must start with either a '+' or a '-'.
************************************************************************/


int opts_found[NUM_OPTIONS];


void help()
{
  int i,j;

  fprintf(stderr,"------------------------------------------------------------------------------\n"); 
  fprintf(stderr,"Help for %s \n",myname);
  fprintf(stderr,"\nUsage:\n");
  fprintf(stderr," %s [-options]",myname);
  if (NUM_FILES == -1)
    fprintf(stderr," file1 file2 ... fileN");
  else
    for(i=0; i<NUM_FILES; i++)
      fprintf(stderr," file%i.mat",i+1);

  fprintf(stderr,"\n");

  fprintf(stderr,"\nSynopsis:\n%s\n", synopsis);
    
  fprintf(stderr,"Options:\n");
  fprintf(stderr," -h,--help\n\tDisplays this help screen.\n");
  for (i=0; i< NUM_OPTIONS; i++)
  {
    fprintf(stderr," %s ",opts[i]);
      
    for(j=0; j<opts_parms[i]; j++)
    {
      if (j>0)
	fprintf(stderr,",");
      fprintf(stderr,opts_parm_text[i][j],j);
    }

    fprintf(stderr,"\n\t");
    if (opts_required[i])
      fprintf(stderr,"*");
    fprintf(stderr,"%s\n",opts_help[i]);
  }
  fprintf(stderr,
	  "\nNotes:\n"
 "Options may be placed in any order.\n"
 "Asterisks (*'s) denote required options.\n"
);
  fprintf(stderr,"\n------------------------------------------------------------------------------\n"); 
  exit(1); 

}




int nfiles;
char *filenames[255];

void parse_cmd (int argc, char *argv[])
{
  int i,j;
  int not_option;
  
  myname = argv[0];
  nfiles = 0;

  for (j=0; j< NUM_OPTIONS; j++)
    opts_found[j] = 0;

  for(i=1; i<argc; i++)
  {
    not_option = 1;
    for(j=0; j<NUM_OPTIONS; j++)
    {  
      if (opts_parms[j] == 0)
      {
	/* check for normal option */
	if ( (strcmp(argv[i],opts[j])) == 0) 
	{
	  not_option = 0;
	  opts_found[j] = 1;
	  opts_func((t_opts)j,NULL);
	  break;
	}
      } else
      {
	/* check for option that has parmeters*/
	if ( strncmp(argv[i],opts[j],strlen(opts[j])) == 0) 
	{
	  char stemp[80];
	  char *sptr = argv[i] + strlen(opts[j]);
	  char *sparms[5];
	  int k;

	  strncpy(stemp,sptr,78);
	  strcat(stemp,",");
	  sptr = stemp;

	  if ((sparms[0] = strtok(sptr,",")) == NULL)
	  {
	    sptr = argv[++i];
	    if ((sparms[0] = strtok(sptr,",")) == NULL)
	    {
	      fprintf(stderr,"%s: Bad input (not enough parameters) for option %s\n",myname,opts[j]);
	      exit(1);
	    }
	  }
	  for(k=1;k<opts_parms[j]; k++)
	  {
	    if ((sparms[k] = strtok(NULL,",")) == NULL)
	    {
	      fprintf(stderr,"%s: Bad input (not enough parameters) for option %s\n",myname,opts[j]);
	      exit(1);
	    }
	  } 
	  if ((sparms[k] = strtok(NULL,",")) != NULL)
	  {
	    fprintf(stderr,"%s: Bad input (too many parameters) for option %s\n",myname,opts[j]);
	    exit(1);
	  }

	  opts_found[j] = 1;
	  not_option = 0;
	  opts_func((t_opts)j,sparms);
	  break;
	} 
      }  /* else (check for opts with parms) */
    
    } /*for(j=o) */
    
    if (not_option)
    {
      if ( ( strcmp(argv[i],"-h") == 0 ) || ( strcmp(argv[i],"-help") == 0 ) || ( strcmp(argv[i],"--help") == 0 ) )
	help();

      if ( ( argv[i][0] == '-' ) || ( argv[i][0] == '+' ))
      {
	fprintf(stderr,"\n%s: Illegal option: %s\n\n", myname, argv[i]);
	exit(1);
      }
      filenames[nfiles++] = argv[i];
    }
    
  } /* for(i=0) */

  
  for (j=0; j< NUM_OPTIONS; j++)
    if (opts_required[j] && !opts_found[j])
    {
      fprintf(stderr,"\n%s: Required option %s not found on command line.\n\n", myname, opts[j]);
      exit(1);
    }
  

#ifdef DEBUG
  {
    int debug;
    fprintf(stderr,"file list:");
    for(debug=0; debug < nfiles; debug++)
      fprintf(stderr," '%s'",filenames[debug]);
    fprintf(stderr,"\n");
  }
#endif

  if ((NUM_FILES != -1) && (nfiles != NUM_FILES))
  {
    fprintf(stderr,"\n%s: Expecting (%d) filenames on the command line\n",myname,NUM_FILES);
    exit(1);
  }

} // parse_cmd()
