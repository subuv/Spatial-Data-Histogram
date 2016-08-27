///
//	Filename: spatialDataHistogram.c
//	returns: integer.
//			 success code of 0
//			 failure code of 1
// Output: written onto a file /tmp/SDH.log
///

//***************************************************************************
// Fall 2015 Advanced Databases Final project
// SPATIAL DATA HISTOGRAM generation using complex aggregates in postgresql.
// Authors: Subramanian Viswanathan, Sneha Rajan, Raviteja Yerrapati
//***************************************************************************

/* Header files
// Postgres.h should always be included first in any source file, because it declares a number of things that you will need from now.
// Definitions for the Postgres function manager and function-call interface. 
// The file fmgr.h must be included by all Postgres modules that either define or call fmgr-callable functions.
// We calculate distance between particles and we need some mathematical operations to be done. Hence math.h is added 
*/
#include "postgres.h"
#include "fmgr.h"
#include <string.h>
#include <math.h>

/* To ensure that a dynamically loaded object file is not loaded into an incompatible server, we add a magic block.
// This allows the server to detect obvious incompatibilities, such as code compiled for a different major version of PostgreSQL 
*/
PG_MODULE_MAGIC;

/* Creating a structure histogram which should hold the counts of all the buckets */
typedef struct hist_entry {
	unsigned long long int d_cnt;   
} bucket;

/* Class level variables declared
// histogram of type bucket is created in which all the distance are finally stored 
*/
bucket *histogram;

int i, j;
double dist = 0.0;

/* The user defined arguments are stored in these variables below respectively. */
int32 sdhNoOfSamples = 0;
int32 sdhBucketWidth = 0;

/* This is the of the 3D space considered */
int32 boxSize = 23000;

/* variables to calculate the inter particle distances */
int32 *xPos;
int32 *yPos;
int32 *zPos;
int h_pos = 0;

/* total number of buckets in the histogram */
int noOfBuckets = 0;

/* File pointer to write histogram output */
FILE *fp;

/* Success Code */
int32 successCode = 0;

/* A macro is created. The calling convention relies on macros to suppress most of the complexity of passing arguments and results. */
Datum spatialDataHistogram(PG_FUNCTION_ARGS);

/* The macro call is done in the name of the source file. It is generally written just before the function itself. */
PG_FUNCTION_INFO_V1(spatialDataHistogram);

/* The actual function starts here. The function name is spatialDataHistogram.
// PG_FUNCTION_ARGS accepts the runtime arguments given by the user. 
// In this program we take the arguments from the table in which the arguments values are specified 
*/
Datum spatialDataHistogram(PG_FUNCTION_ARGS) {
	double xa, xb, ya, yb, za, zb;
	
	/* An initial check to verify if the user has provided the two required inputs:
	// 1. First argument should be the number of samples
	// 2. Second argument should be user defined bucket width 
	*/
   	if( PG_ARGISNULL(0) || PG_ARGISNULL(1)) {
   		ereport(ERROR,
            (
             	errcode(ERRCODE_NULL_VALUE_NOT_ALLOWED),
             	errmsg("Null values are not allowed"),
             	errdetail("The number of samples and bucket width inputs are not specified."),
             	errhint("Provide the number of samples and bucket width")
            )
        );
        successCode = 1;
      	PG_RETURN_NULL();
   	}

   	/* PG_GETARG_INT32() macro is used to fetch the arguments. */
   	sdhNoOfSamples = PG_GETARG_INT32(0);
   	sdhBucketWidth = PG_GETARG_INT32(1);
   	
	if (sdhNoOfSamples > 0 && sdhBucketWidth > 0) {
		noOfBuckets = (int)(boxSize * 1.732 / sdhBucketWidth) + 1;

		/* Palloc is used to allocate memory to the following variables */
		xPos = (int32 *)palloc((sizeof(int32)*sdhNoOfSamples));
		yPos = (int32 *)palloc((sizeof(int32)*sdhNoOfSamples));
		zPos = (int32 *)palloc((sizeof(int32)*sdhNoOfSamples));
		histogram = (bucket *)palloc(sizeof(bucket)*noOfBuckets);

		/* memset initialises the variables with their default values */
		memset((int32 *) xPos, -1, sizeof(int32));
		memset((int32 *) yPos, -1, sizeof(int32));
		memset((int32 *) zPos, -1, sizeof(int32));
		memset((bucket *) histogram, 0, (sizeof(bucket)*noOfBuckets));

		/* assign the coordinates of the particles in the box */
		srand(1);

		/* all the X, Y, Z coordinates of the 3D particles are assigned and stored in the corresponding arrays. */
		for(i = 0;  i < sdhNoOfSamples; i++) {
			xPos[i] = ((double)(rand()) / RAND_MAX) * boxSize;
			yPos[i] = ((double)(rand()) / RAND_MAX) * boxSize;
			zPos[i] = ((double)(rand()) / RAND_MAX) * boxSize;
		}

		/* Let us build the histogram now */
		for(i = 0; i < sdhNoOfSamples; i++) {
			for(j = i+1; j < sdhNoOfSamples; j++) {
				xa = xPos[i];
				xb = xPos[j];
				
				ya = yPos[i];
				yb = yPos[j];
				
				za = zPos[i];
				zb = zPos[j];	
				
				/* Inter-particle Distances are calculated using the below function. SQRT function is used here. */
				dist = sqrt((xa - xb)*(xa - xb) + (ya - yb)*(ya - yb) + (za - zb)*(za - zb));
				
				/* The position of the distance in the histogram is calculated */
				h_pos = (int) (dist/sdhBucketWidth);
				
				histogram[h_pos].d_cnt++;
			}
		}

		/* free the dynamic variable pointers */
		pfree(xPos);
		pfree(yPos);
		pfree(zPos);
		
		/* Variables for building the output array. */
		char outputText[15000];
		char tempOutput[15000] = "=";
		int32 tempOutputLen = 0;
		text *outputFinal = NULL;
	
		/* memset initialises the variables with their default values */
		memset((char *) outputText, ' ', sizeof(char));
		memset((char *) tempOutput, ' ', sizeof(char));

		sprintf(tempOutput, "%s %d %s %d %s","============================================================================================== \n\t THE SPATIAL DATA HISTOGRAM FOR", sdhNoOfSamples, "SAMPLES AND", sdhBucketWidth, "BUCKET WIDTH IS:\t");
		strcat(tempOutput, "\n==============================================================================================");

		/* Write the output histogram from the structure bucket to the text variable */
   		for(i=0; i<noOfBuckets; i++) {
			if(i%5==0){ 
				/* prints 5 buckets in a row */
				sprintf(outputText, "\n%02d: ", i);
				strcat(tempOutput, outputText);
			}
			sprintf(outputText, "%15lld ", histogram[i].d_cnt);
			strcat(tempOutput, outputText);
	  		
	  		if(i != noOfBuckets-1) {
				sprintf(outputText, "%s", "| ");
				strcat(tempOutput, outputText);
			}
		}
		
		/* free the dynamic variable pointers */
		pfree(histogram);

		strcat(tempOutput, "\n==============================================================================================");

		tempOutputLen = strlen(tempOutput);

		/* Allocate memory and set data structure size for the final output */
   		outputFinal = (text *)palloc(tempOutputLen);
   		memset((text *) outputFinal, ' ', sizeof(text));
		
		SET_VARSIZE(outputFinal, tempOutputLen + VARHDRSZ);

   		/* Construct the final output to be written onto the file */
		strcpy(VARDATA(outputFinal), tempOutput);

		/* Write the output to a file */
		fp = fopen("/tmp/SDH.log", "a+");
		fwrite(VARDATA(outputFinal), 1, VARSIZE(outputFinal)-VARHDRSZ, fp);
		fprintf(fp, "\n");
		fclose(fp);

		/* free the dynamic variable pointers */
		pfree(outputFinal);
	}
	/* Return success code */
	PG_RETURN_INT32(successCode);
}
