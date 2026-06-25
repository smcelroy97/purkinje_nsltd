TITLE getenv

COMMENT
get environment vars for parallel batch simulation
ENDCOMMENT

NEURON {
   SUFFIX nothing

   GLOBAL k


 }

ASSIGNED {

   k
 

  
  
}

VERBATIM
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
ENDVERBATIM 

FUNCTION loadenv() {

VERBATIM

  char *j;


  j = getenv("PBS_ARRAYID");

  
  k =atof(j);

   
 
ENDVERBATIM 
 

 }

