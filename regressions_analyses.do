
* Replication material for 'Stratified Pathways to Italy’s “Latest-Late” Transition to Adulthood'

* Date of the last update**: 2023-07-26

* Title: Stratified Pathways to Italy’s “Latest-Late” Transition to Adulthood

* Author: Badolato Luca (badolato.3@osu.edu), Department of Sociology and Institue for Population Research, The Ohio State University

* Journal: Advances in Life Course Research


* Description: script to compute the regression analyses and replicate Figure 4, Table 5, Figure A4, Table A2, Table A3, Table A5, Table A6. Takes as input data.dta, data_cluster_optimal_matching.dta, and data_cluster_optimal_matching_gender_specific.dta. 

clear all

cd "Insert working directory"

* Load data
use "data_cluster_optimal_matching.dta", clear
merge 1:1 _n using "data.dta"

gen cluster = ""
replace cluster = "EE" if clust == 1
replace cluster = "NT" if clust  == 2
replace cluster = "FT" if clust  == 3
replace cluster = "LE" if clust  == 4
replace cluster = "EI" if clust == 5
replace cluster = "TI" if clust  == 6
replace cluster = "EST" if clust  == 7

encode cohort, gen(cohort_c)
encode cluster, gen(cluster_c)

* Baseline model: cohort only

* Women
mlogit cluster_c i.cohort_c [aweight=anweight] if gender==2, base(5) 

estimates store m0w

margins, dydx(cohort) post
estimates restore m0w

margins cohort_c, atmeans predict(outcome(1)) post
estimates store m01w 
estimates restore m0w

margins cohort_c, atmeans predict(outcome(2)) post
estimates store m02w
estimates restore m0w

margins cohort_c, atmeans predict(outcome(3)) post
estimates store m03w
estimates restore m0w

margins cohort_c, atmeans predict(outcome(4)) post
estimates store m04w
estimates restore m0w

margins cohort_c, atmeans predict(outcome(5)) post
estimates store m05w
estimates restore m0w

margins cohort_c, atmeans predict(outcome(6)) post
estimates store m06w
estimates restore m0w

margins cohort_c, atmeans predict(outcome(7)) post
estimates store m07w
estimates restore m0w

* Men
mlogit cluster_c i.cohort_c [aweight=anweight] if gender==1, base(5) 

estimates store m0m
margins, dydx(cohort) post

estimates restore m0m
margins cohort_c, atmeans predict(outcome(1)) post
estimates store m01m
estimates restore m0m

margins cohort_c, atmeans predict(outcome(2)) post
estimates store m02m
estimates restore m0m

margins cohort_c, atmeans predict(outcome(3)) post
estimates store m03m
estimates restore m0m

margins cohort_c, atmeans predict(outcome(4)) post
estimates store m04m
estimates restore m0m

margins cohort_c, atmeans predict(outcome(5)) post
estimates store m05m
estimates restore m0m

margins cohort_c, atmeans predict(outcome(6)) post
estimates store m06m
estimates restore m0m

margins cohort_c, atmeans predict(outcome(7)) post
estimates store m07m
estimates restore m0m

coefplot m01w m01m, vertical recast(connected) ciopts(recast(rcap)) offset(0) ///
ysize(6.5) xsize(11) graphregion(margin(large)) ///
ytitle("Predicted probability", size(medium)) ///
xlabel(, labsize(small)) ///
legend(rows(1) size(small) label(2 "Women") label(4 "Men")) ///
xtitle("Cohort", size(medium)) ///
graphregion(color(white)) 

coefplot m02w m02m, vertical recast(connected) ciopts(recast(rcap)) offset(0) ///
ysize(6.5) xsize(11) graphregion(margin(large)) ///
ytitle("Predicted probability", size(medium)) ///
xlabel(, labsize(small)) ///
legend(rows(1) size(small) label(2 "Women") label(4 "Men")) ///
xtitle("Cohort", size(medium)) ///
graphregion(color(white)) 

coefplot m03w m03m, vertical recast(connected) ciopts(recast(rcap)) offset(0) ///
ysize(6.5) xsize(11) graphregion(margin(large)) ///
ytitle("Predicted probability", size(medium)) ///
xlabel(, labsize(small)) ///
legend(rows(1) size(small) label(2 "Women") label(4 "Men")) ///
xtitle("Cohort", size(medium)) ///
graphregion(color(white)) 

coefplot m04w m04m, vertical recast(connected) ciopts(recast(rcap)) offset(0) ///
ysize(6.5) xsize(11) graphregion(margin(large)) ///
ytitle("Predicted probability", size(medium)) ///
xlabel(, labsize(small)) ///
legend(rows(1) size(small) label(2 "Women") label(4 "Men")) ///
xtitle("Cohort", size(medium)) ///
graphregion(color(white)) 

coefplot m05w m05m, vertical recast(connected) ciopts(recast(rcap)) offset(0) ///
ysize(6.5) xsize(11) graphregion(margin(large)) ///
ytitle("Predicted probability", size(medium)) ///
xlabel(, labsize(small)) ///
legend(rows(1) size(small) label(2 "Women") label(4 "Men")) ///
xtitle("Cohort", size(medium)) ///
graphregion(color(white)) 

coefplot m06w m06m, vertical recast(connected) ciopts(recast(rcap)) offset(0) ///
ysize(6.5) xsize(11) graphregion(margin(large)) ///
ytitle("Predicted probability", size(medium)) ///
xlabel(, labsize(small)) ///
legend(rows(1) size(small) label(2 "Women") label(4 "Men")) ///
xtitle("Cohort", size(medium)) ///
graphregion(color(white)) 

coefplot m07w m07m, vertical recast(connected) ciopts(recast(rcap)) offset(0) ///
ysize(6.5) xsize(11) graphregion(margin(large)) ///
ytitle("Predicted probability", size(medium)) ///
xlabel(, labsize(small)) ///
legend(rows(1) size(small) label(2 "Women") label(4 "Men")) ///
xtitle("Cohort", size(medium)) ///
graphregion(color(white)) 

* Model: baseline + parental education

* Education of parents as ISCED scale
gen educ_father = eiscedf if eiscedf>=1 & eiscedf<=8 
replace educ_father = 1 if eiscedf==2
replace educ_father = 2 if eiscedf==3
replace educ_father = 3 if eiscedf==4 | eiscedf==5
replace educ_father = 4 if eiscedf==6 | eiscedf==7 | eiscedf==8

tab educ_father

gen educ_mother = eiscedm if eiscedm>=1 & eiscedm<=7 
replace educ_mother = 1 if eiscedm==2
replace educ_mother = 2 if eiscedm==3
replace educ_mother = 3 if eiscedm==4 | eiscedm==5
replace educ_mother = 4 if eiscedm==6 | eiscedm==7 | eiscedm==8
tab educ_mother

gen educ_parents=max(educ_father, educ_mother)
tab educ_parents

* Women
mlogit cluster_c i.cohort_c i.educ_parents  [aweight=anweight]  if gender==2, base(5) 
estimates store m1w

margins, dydx(cohort) post
eststo

estimates restore m1w

margins, dydx(educ_parents) post
eststo

* Men
mlogit cluster_c i.cohort_c i.educ_parents  [aweight=anweight]  if gender==1, base(5) 
estimates store m1m

margins, dydx(cohort) post
eststo

estimates restore m1m

margins, dydx(educ_parents) post
eststo

* Robustness check: baseline + parental occupation 
*

generate occupation_father = occf14b

*farm workers and unskilled workers
replace occupation_father = 1 if occf14b==8 | occf14b==9
*semi-skilled workers
replace occupation_father = 2 if occf14b==7 
* sales occupations, skilled worker, clerical occupations and service occupations
replace occupation_father = 3 if occf14b==3 | occf14b==4 | occf14b==5 | occf14b==6
*Higher administrator occupations, Professional and technical occupations
replace occupation_father = 4 if  occf14b==2 | occf14b==1 

generate occupation_mother = occm14b

*farm workers and unskilled workers
replace occupation_mother = 1 if occm14b==8 | occm14b==9
*semi-skilled workers
replace occupation_mother = 2 if occm14b==7 
* sales occupations, skilled worker, clerical occupations and service occupations
replace occupation_mother = 3 if occm14b==3 | occm14b==4 | occm14b==5 | occm14b==6
*Higher administrator occupations, Professional and technical occupations
replace occupation_mother = 4 if  occm14b==2 | occm14b==1 

spearman occupation_father occupation_mother

generate occupation_parents=max(occupation_mother, occupation_father)

* Women
mlogit cluster_c i.cohort_c i.occupation_parents [aweight=anweight] if gender==2, base(5) 
estimates store m5w

margins, dydx(cohort) post
eststo

estimates restore m5w

margins, dydx(social_class) post
eststo

* Men
mlogit cluster_c i.cohort_c i.occupation_parents [aweight=anweight] if gender==1, base(5) 
estimates store m5m

margins, dydx(cohort) post
eststo

estimates restore m5m

margins, dydx(social_class) post
eststo

* Additional analyses: gender-specific typologies

* Men
use "data_cluster_optimal_matching_gender_specific.dta", replace

* Keep only men 
keep if gender == 1

gen cluster = ""
replace cluster = "NT" if clust == 1
replace cluster = "LE" if clust  == 2
replace cluster = "FT" if clust  == 3
replace cluster = "EE" if clust  == 4
replace cluster = "TI" if clust == 5

encode cohort, gen(cohort_c)
encode cluster, gen(cluster_c)

gen educ_father = eiscedf if eiscedf>=1 & eiscedf<=8 
replace educ_father = 1 if eiscedf==2
replace educ_father = 2 if eiscedf==3
replace educ_father = 3 if eiscedf==4 | eiscedf==5
replace educ_father = 4 if eiscedf==6 | eiscedf==7 | eiscedf==8

tab educ_father

gen educ_mother = eiscedm if eiscedm>=1 & eiscedm<=7 
replace educ_mother = 1 if eiscedm==2
replace educ_mother = 2 if eiscedm==3
replace educ_mother = 3 if eiscedm==4 | eiscedm==5
replace educ_mother = 4 if eiscedm==6 | eiscedm==7 | eiscedm==8
tab educ_mother

gen educ_parents=max(educ_father, educ_mother)
tab educ_parents

mlogit cluster_c i.cohort_c [aweight=anweight], base(5) 

estimates store m0m
margins, dydx(cohort) post

estimates restore m0m
margins cohort_c, atmeans predict(outcome(1)) post
estimates store m01m
estimates restore m0m

margins cohort_c, atmeans predict(outcome(2)) post
estimates store m02m
estimates restore m0m

margins cohort_c, atmeans predict(outcome(3)) post
estimates store m03m
estimates restore m0m

margins cohort_c, atmeans predict(outcome(4)) post
estimates store m04m
estimates restore m0m

margins cohort_c, atmeans predict(outcome(5)) post
estimates store m05m
estimates restore m0m

coefplot m01m, vertical recast(connected) ciopts(recast(rcap)) offset(0) ///
ysize(6.5) xsize(11) graphregion(margin(large)) ///
ytitle("Predicted probability", size(medium)) ///
xlabel(, labsize(small)) ///
legend(rows(1) size(small) label(2 "Men")) ///
xtitle("Cohort", size(medium)) ///
graphregion(color(white)) 

coefplot m02m, vertical recast(connected) ciopts(recast(rcap)) offset(0) ///
ysize(6.5) xsize(11) graphregion(margin(large)) ///
ytitle("Predicted probability", size(medium)) ///
xlabel(, labsize(small)) ///
legend(rows(1) size(small) label(2 "Men")) ///
xtitle("Cohort", size(medium)) ///
graphregion(color(white)) 

coefplot m03m, vertical recast(connected) ciopts(recast(rcap)) offset(0) ///
ysize(6.5) xsize(11) graphregion(margin(large)) ///
ytitle("Predicted probability", size(medium)) ///
xlabel(, labsize(small)) ///
legend(rows(1) size(small) label(2 "Men")) ///
xtitle("Cohort", size(medium)) ///
graphregion(color(white)) 

coefplot m04m, vertical recast(connected) ciopts(recast(rcap)) offset(0) ///
ysize(6.5) xsize(11) graphregion(margin(large)) ///
ytitle("Predicted probability", size(medium)) ///
xlabel(, labsize(small)) ///
legend(rows(1) size(small) label(2 "Men")) ///
xtitle("Cohort", size(medium)) ///
graphregion(color(white)) 

coefplot m05m, vertical recast(connected) ciopts(recast(rcap)) offset(0) ///
ysize(6.5) xsize(11) graphregion(margin(large)) ///
ytitle("Predicted probability", size(medium)) ///
xlabel(, labsize(small)) ///
legend(rows(1) size(small) label(2 "Men")) ///
xtitle("Cohort", size(medium)) ///
graphregion(color(white)) 

* Adding parental education
mlogit cluster_c i.cohort_c i.educ_parents  [aweight=anweight], base(5) 
estimates store m1m

margins, dydx(cohort) post

estimates restore m1m

margins, dydx(educ_parents) post

* Women
use "data_cluster_optimal_matching_gender_specific.dta", replace

* Keep only women 
keep if gender == 2

gen cluster = ""
replace cluster = "LE" if clust == 1
replace cluster = "EST" if clust  == 2
replace cluster = "NT" if clust  == 3
replace cluster = "EI" if clust  == 4
replace cluster = "FT" if clust == 5

encode cohort, gen(cohort_c)
encode cluster, gen(cluster_c)

gen educ_father = eiscedf if eiscedf>=1 & eiscedf<=8 
replace educ_father = 1 if eiscedf==2
replace educ_father = 2 if eiscedf==3
replace educ_father = 3 if eiscedf==4 | eiscedf==5
replace educ_father = 4 if eiscedf==6 | eiscedf==7 | eiscedf==8

tab educ_father

gen educ_mother = eiscedm if eiscedm>=1 & eiscedm<=7 
replace educ_mother = 1 if eiscedm==2
replace educ_mother = 2 if eiscedm==3
replace educ_mother = 3 if eiscedm==4 | eiscedm==5
replace educ_mother = 4 if eiscedm==6 | eiscedm==7 | eiscedm==8
tab educ_mother

gen educ_parents=max(educ_father, educ_mother)
tab educ_parents

mlogit cluster_c i.cohort_c [aweight=anweight], base(5) 

estimates store m0m
margins, dydx(cohort) post

estimates restore m0m
margins cohort_c, atmeans predict(outcome(1)) post
estimates store m01m
estimates restore m0m

margins cohort_c, atmeans predict(outcome(2)) post
estimates store m02m
estimates restore m0m

margins cohort_c, atmeans predict(outcome(3)) post
estimates store m03m
estimates restore m0m

margins cohort_c, atmeans predict(outcome(4)) post
estimates store m04m
estimates restore m0m

margins cohort_c, atmeans predict(outcome(5)) post
estimates store m05m
estimates restore m0m

coefplot m01m, vertical recast(connected) ciopts(recast(rcap)) offset(0) ///
ysize(6.5) xsize(11) graphregion(margin(large)) ///
ytitle("Predicted probability", size(medium)) ///
xlabel(, labsize(small)) ///
legend(rows(1) size(small) label(2 "Women")) ///
xtitle("Cohort", size(medium)) ///
graphregion(color(white)) 

coefplot m02m, vertical recast(connected) ciopts(recast(rcap)) offset(0) ///
ysize(6.5) xsize(11) graphregion(margin(large)) ///
ytitle("Predicted probability", size(medium)) ///
xlabel(, labsize(small)) ///
legend(rows(1) size(small) label(2 "Women")) ///
xtitle("Cohort", size(medium)) ///
graphregion(color(white)) 

coefplot m03m, vertical recast(connected) ciopts(recast(rcap)) offset(0) ///
ysize(6.5) xsize(11) graphregion(margin(large)) ///
ytitle("Predicted probability", size(medium)) ///
xlabel(, labsize(small)) ///
legend(rows(1) size(small) label(2 "Women")) ///
xtitle("Cohort", size(medium)) ///
graphregion(color(white)) 

coefplot m04m, vertical recast(connected) ciopts(recast(rcap)) offset(0) ///
ysize(6.5) xsize(11) graphregion(margin(large)) ///
ytitle("Predicted probability", size(medium)) ///
xlabel(, labsize(small)) ///
legend(rows(1) size(small) label(2 "Women")) ///
xtitle("Cohort", size(medium)) ///
graphregion(color(white)) 

coefplot m05m, vertical recast(connected) ciopts(recast(rcap)) offset(0) ///
ysize(6.5) xsize(11) graphregion(margin(large)) ///
ytitle("Predicted probability", size(medium)) ///
xlabel(, labsize(small)) ///
legend(rows(1) size(small) label(2 "Women")) ///
xtitle("Cohort", size(medium)) ///
graphregion(color(white)) 

* Adding parental education
mlogit cluster_c i.cohort_c i.educ_parents [aweight=anweight], base(5) 
estimates store m1m

margins, dydx(cohort) post

estimates restore m1m

margins, dydx(educ_parents) post

* End
