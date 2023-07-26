
* Replication material for 'Stratified Pathways to Italy’s “Latest-Late” Transition to Adulthood'

* Date of the last update: 2023-07-26

* Title: Stratified Pathways to Italy’s “Latest-Late” Transition to Adulthood

* Author: Badolato Luca (badolato.3@osu.edu), Department of Sociology and Institue for Population Research, The Ohio State University

* Journal: Advances in Life Course Research


* Description: Script to clean the data, select the variable of interest, and save the cleaned data as data.dta

* Clear the enviroment and set the working directory
clear

cd ""

* Load the dataset as downloaded from https://www.europeansocialsurvey.org/ (Round 9, version 03)
use "ESS9e03.dta", replace 

* Keep only Italy
keep if cntry == "IT" 

* Rename variables of interest: 
* YW (Year Work), YL (Year Leaving), YC ( Year Cohabitation or Marriage), and YB (Year child was Born)
rename gndr gender

rename pdempyr YW 
rename lvpntyr YL
rename lvptnyr YC
rename fcldbrn YB

tab YW, miss
tab YL, miss
tab YC, miss
tab YB, miss

* Assign the value 0 for not applicable answers (for example "never left the parental home")
* and drop the other observations with missing values.
foreach var of varlist YW YC YB {

  drop if `var' == .b | `var' == .c | `var'== .d
  replace `var' = 0 if `var' == .a

}

drop if YL == .a | YL == .b | YL == .c | YL== .d | YL == 1111
drop if yrbrn == .a | yrbrn == .b | yrbrn == .c | yrbrn== .d 

drop if agea < 30

* Generate YWage, YCage, YLage, and YBage, the age when the specific event ocours.
* If it has never occured, 0 is the default value.
foreach var of varlist YW YL YC YB {
  gen `var'age = `var' - yrbrn
  replace `var'age = 0 if `var'age < 0
}

* For each possible age (13:30), create the variables "age"Y"event" equal to 1 if
* the event occurs and equal to 0 if not. 
foreach i of num 13/30 {

       gen age`i'YW = ""
	   gen age`i'YL = ""
	   gen age`i'YC = ""
	   gen age`i'YB = ""
	   
	   replace age`i'YW = "1" if YWage <= `i' & YWage != 0
       replace age`i'YW = "0" if YWage > `i' | YWage == 0
	   
	   replace age`i'YL = "1" if YLage <= `i' & YLage != 0
       replace age`i'YL = "0" if YLage > `i' | YLage == 0
	   
	   replace age`i'YC = "1" if YCage <= `i' & YCage != 0
       replace age`i'YC = "0" if YCage > `i' | YCage == 0
	   
	   replace age`i'YB = "1" if YBage <= `i' & YBage != 0
       replace age`i'YB = "0" if YBage > `i' | YBage == 0
		}
 
* Finally, generate the timing variables and put together the four 
* events to obtain the unique sequence with the states of interest, of the form "WLCB". 
foreach i of num 13/30 {
     
       gen age`i' = age`i'YW+age`i'YL+age`i'YC+age`i'YB
	  
        }
	
* Define cohorts
gen cohort = ""
replace cohort = "1928 - 1942" if yrbrn >= 1928 &  yrbrn <= 1942
replace cohort = "1943 - 1957" if yrbrn >= 1943 &  yrbrn <= 1957
replace cohort = "1958 - 1967" if yrbrn >= 1958 &  yrbrn <= 1967
replace cohort = "1968 - 1977" if yrbrn >= 1968 &  yrbrn <= 1977
replace cohort = "1978 - 1988" if yrbrn >= 1978 &  yrbrn <= 1988

foreach i of num 13/30 {
 tab age`i'YW
 tab age`i'YL
 tab age`i'YC
 tab age`i'YB
}

foreach i of num 13/30 {
 tab age`i'YB
 }
 
* Save data
save data, replace

* End

