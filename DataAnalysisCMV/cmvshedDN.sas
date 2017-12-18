***********************************************************************
*
* ASCII Datafile Name:  shedcmv.oct_96_fin.dat
*
* Create Date:          October 12, 1996
*
* Protocol:             ACTG 181
*
* Purpose:              Last Negative and First Positive CMV Cultures,
*                       and CD4 Counts and Dates for Positive Cultures
*
* Programmer:		Eva Stamenovic
*
**********************************************************************;

/* Version adapted by Daniel Nevo, July 5, 2017 
Added proc export to have the data in .csv form
Original file was found on http://hedwig.mgh.harvard.edu/biostatistics/node/32
in interval
*/

options nocenter errors=3 linesize=90 ps=500 ;

data cmvshed ;
   length PID 4 STUDY $ 4 SEX $ 1 RACE CMVIND sheddind
	cmvdiag digcd4 DEATHDAT STARTRX OFFSTDAT DIGDATE DIGCD4DT
        RAND081  RAND181  SBAS181D SBASCD4D FIRCD4DT LASCD4DT
	SBASCD4 FIRSTCD4 LASTCD4
	blposind BLOODEDT BLDNEGDT BLDPOSDT BPOCD4DT BPOCD4 
	urposind URINEEDT URNNEGDT URNPOSDT UPOCD4DT UPOCD4 4 ;
   format DEATHDAT STARTRX OFFSTDAT DIGDATE DIGCD4DT  
	RAND081  RAND181  SBAS181D SBASCD4D FIRCD4DT LASCD4DT
	BLOODEDT BLDNEGDT BLDPOSDT BPOCD4DT
	URINEEDT URNNEGDT URNPOSDT UPOCD4DT date7. ;

   infile	'/folders/myfolders/cmvshed.dat' 
		lrecl=80 missover firstobs=60 ;

   input
	/* line #1 */
	  @7 PID @14 STUDY @18 SEX 
	 @20 RACE @22 CMVIND @24 SHEDDIND @26 DEATHCEN 
	 @28 cmvdiag @30 digcd4 4.
         @35 (DEATHDAT STARTRX OFFSTDAT DIGDATE DIGCD4DT )(date7. +1) /
	/* line #2 */
	 @14 (RAND081 RAND181  SBAS181D SBASCD4D FIRCD4DT
	     LASCD4DT ) (date7. +1) SBASCD4 FIRSTCD4 LASTCD4 /
	/* line #3 */
         @14 blposind 
	 @16 (BLOODEDT BLDNEGDT BLDPOSDT BPOCD4DT )(date7. +1) BPOCD4 /
	/* line #4 */
         @14 urposind 
         @16 (URINEEDT URNNEGDT URNPOSDT UPOCD4DT )(date7. +1) UPOCD4 ;

   label /* line #1 */ 
	PID     = 'Patient ID'
	STUDY    = 'Protocol Number'
	SEX      = 'Gender'
	RACE     = 'Race/Ethnicity'
	CMVIND   = 'CMV Baseline-Diag. Indicator'
	SHEDDIND = '181 CMV Shedding Indicator'
	DEATHCEN = '081 Death Censor'
	DEATHDAT = '081 Death Date'
	CMVDIAG  = 'CMV Diagnosis Indicator'
	DIGDATE  = 'CMV Diagnosis Date'
	DIGCD4   = 'CD4 Count for CMV Diagnosis'
	DIGCD4DT = 'CD4 Date for CMV Diagnosis'
	STARTRX  = '181 Protocol Rx Started'
	OFFSTDAT = '181 Offstudy Date'

	/* line #2 */
	RAND081  = '081 Randomization Date'
	RAND181  = '181 Randomization Date'
	SBAS181D = '181 Baseline Date'
	SBASCD4D = '181 Baseline CD4 Date'
	SBASCD4  = '181 Baseline CD4 Count'
	FIRCD4DT = 'Earliest 181 CD4 Date'
	FIRSTCD4 = 'Earliest 181 CD4 Count'
	LASCD4DT = 'Last 181 CD4 Date'
	LASTCD4  = 'Last 181 CD4 Count'

	/* line #3 */
	BLPOSIND = '181 Blood Shedding Indicator'
	BLOODEDT = 'Blood Earliest Date'
	BLDNEGDT = 'Date Last Negative Blood'
	BLDPOSDT = 'Date 1st Positive Blood'
	BPOCD4DT = 'CD4 Date of Positive Blood Shed.'
	BPOCD4   = 'CD4 Count of Positive Blood Shed.'

	/* line #4 */
	URPOSIND = '181 Urine Shedding Indicator'
	URINEEDT = 'Urine Earliest Date'
	URNNEGDT = 'Date Last Negative Urine'
	URNPOSDT = 'Date 1st Positive Urine'
	UPOCD4DT = 'CD4 Date of Positive Urine Shed.'
	UPOCD4   = 'CD4 Count of Positive Urine Shed.' ;
run ;

proc contents ;
	title 'ACTG 181' ;
	title2 'Last Negative and First Positive CMV Cultures';
run ;
proc freq ;
        tables STUDY SEX RACE CMVIND CMVDIAG DEATHCEN 
		SHEDDIND BLPOSIND URPOSIND BLPOSIND*URPOSIND 
                / missing norow nocol nopercent missing ;
run ;

proc export data=cmvshed outfile="/folders/myfolders/IC_data/cmvshed.txt"
	dbms=tab
   	replace;
   run;

proc print uniform ; * label; *by pid ;
	var PID SEX RACE CMVIND DEATHCEN DEATHDAT
            RAND081 RAND181 STARTRX OFFSTDAT 
            BLOODEDT BLDNEGDT BLDPOSDT BPOCD4DT BPOCD4
            URINEEDT URNNEGDT URNPOSDT UPOCD4DT UPOCD4 ;
run ;



