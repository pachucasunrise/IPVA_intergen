*09.07.24
*Analyses v.10

****************************
*SET UP

timer on 1

*Set drives

ssc install mim

*******************************
timer on 2
*TABLES 1-4
*Table 1: Association of parental IPVA with offspring IPVA, overall and by IPVA type, crude, adjusted A, adjusted B (A + ACEs): victimisation as outcome
*Table 2: As for Table 1, but perpetration as outcome
*Table 3: Prevalence, single PAFs (joint PAFs handled in a separate loop), and RDs: victimisation as outcome
*Table 4: As for Table 3, but perpetration as outcome

*The following exposures are mother-reported, '_acts' are those that only use acts-based measures (sensitivity analysis)
global exposure "dv_mrep psych_mrep cc_mrep phys_mrep dv_mrep_acts psych_mrep_acts phys_mrep_acts"
global outcome "vic_1821 per_1821"

putexcel set "tables", sheet("T1_T4") modify
putexcel A1="Exposure" B1="Outcome" C1="Adjusted" D1="beta" E1="SE" F1="P Value" G1="Sex" H1="no_perc" I1="yes_perc" J1="RR" K1="LCI" L1="UCI" M1="PAF" N1="LPAF" O1="UPAF" P1="pc" Q1="prev" R1="RD" S1="LRD" T1="URD" U1="Risk ref"

local x=1
*Looping over sex, types of DV exposure, and vic/perp as an outcome
*kz021 is sex variable in ALSPAC		
forvalues s=1/2{
	foreach exp of global exposure{
		foreach out of global outcome{
			use "imputed_sample_bysex_`exp'_`out'.dta", replace
			keep if kz021==`s'
			
			*Generate variables to be filled in later
			gen imp_risk_ref = .
			gen imp_rd = .
			gen imp_rd_se = .
			gen imp_rd_diff_sq = .
			gen imp_rd_se_sq = .
			
			*Unadjusted
			mim: glm `out' `exp', family(poisson) link(log) vce(robust)
			mim, storebv
			
			local beta = _b[`exp']
			local se = _se[`exp']
			
			test `exp'		
			local pvalue = `r(p)'
			
			*Percentages
			bysort _mj `exp': egen denom = count(aln) if _mj != 0
			bysort _mj `exp': egen num2 = count(aln) if _mj != 0 & `out'==1
			bysort _mj `exp': egen num = min(num2)
			gen prop = num/denom if _mj != 0
			bysort _mj `exp': egen seq = seq() if _mj != 0 
			summ prop if seq == 1 & `exp'==0
			local prop_no = r(mean)
			summ prop if seq == 1 & `exp'==1
			local prop_yes = r(mean)
			drop denom num num2 prop seq
			
			*Now calculating proportion of exposure among cases
			bysort _mj `out': egen denom = count(aln) if _mj != 0
			bysort _mj `out': egen num2 = count(aln) if _mj != 0 & `exp'==1
			bysort _mj `out': egen num = min(num2)
			gen prop = num/denom if _mj != 0
			bysort _mj `out': egen seq = seq() if _mj != 0 
			summ prop if seq == 1 & `out'==1
			local pc = r(mean)
			drop denom num num2 prop seq
			
			local rr=2.718^(`beta')
			local lci=2.718^(`beta'-1.96*`se')
			local uci=2.718^(`beta'+1.96*`se')
			local paf=`pc'*(1-(1/`rr'))*100
			local lpaf=`pc'*(1-(1/`lci'))*100
			local upaf=`pc'*(1-(1/`uci'))*100
			
			bysort _mj: egen denom = count(aln) if _mj != 0
			bysort _mj: egen num2 = count(aln) if _mj != 0 & `exp'==1
			bysort _mj: egen num = min(num2)
			gen prop = num/denom if _mj != 0
			bysort _mj: egen seq = seq() if _mj != 0 
			summ prop if seq == 1
			local prev = r(mean)
			drop denom num num2 prop seq
			
			local x=`x'+1
			putexcel A`x'="`exp'" B`x'="`out'" C`x'="No" ///
			D`x'=`beta' E`x'=`se' F`x'=`pvalue' G`x'="`s'" H`x'=`prop_no'*100 I`x'=`prop_yes'*100 J`x'=`rr' K`x'=`lci' L`x'=`uci' M`x'=`paf' N`x'=`lpaf' O`x'=`upaf' P`x'=`pc'*100 Q`x'=`prev'*100
			
			*Adjusted model A (without prenatal IPVA)
			xi: mim: glm `out' `exp' marstat matage birthwgt parity smokpreg i.highested i.hhsocclas, family(poisson) link(log) vce(robust)
			mim, storebv
			
			local beta = _b[`exp']
			local se = _se[`exp']
			
			test `exp'		
			local pvalue = `r(p)'
			
			local rr=2.718^(`beta')
			local lci=2.718^(`beta'-1.96*`se')
			local uci=2.718^(`beta'+1.96*`se')
			
			*PAFs
			local paf=`pc'*(1-(1/`rr'))*100
			local lpaf=`pc'*(1-(1/`lci'))*100
			local upaf=`pc'*(1-(1/`uci'))*100
		
			*Predicted probs
			*To average later
			bysort _mj `exp': egen seq = seq() if _mj != 0 

			*Estimate risk difference within each imputation
			forvalues i=1/35{
				logistic `out' `exp' marstat matage birthwgt parity smokpreg i.highested i.hhsocclas preg_mrep ///
				if _mj==`i'
				*Risks of outcome overall within current imputation i (different estimated risks across the individuals)
				predict pred_prob if _mj==`i', pr
				
				*Risks of outcome amongst the unexposed
				summ pred_prob if `exp'==0 & _mj==`i', det
				*Average risk (and variance) of outcome amongst the unexposed (save this for later)
				local pred_prob_0 = r(mean)
				*local pred_se_0 = r(sd)/(sqrt(N))
				local pred_var_0 = r(Var)
				
				*Risks of outcome amongst the exposed
				summ pred_prob if `exp'==1 & _mj==`i', det
				local pred_prob_1 = r(mean)
				local pred_var_1 = r(Var)
				
				*Only one row's value of imp_risk_ref, imp_rd, imp_rd_se is replaced for that particular imputation i
				*(will have 35 values replaced across whole spreadsheet when loop has finished)
				replace imp_risk_ref = `pred_prob_0' if `exp'==1 & seq==1 & _mj==`i'
				replace imp_rd = `pred_prob_1'-`pred_prob_0' if `exp'==1 & seq==1 & _mj==`i'
				*prob1 and prob0 have no covariance (as probability of same outcome, conditional on exposure status. Noone can be both)
				replace imp_rd_se = sqrt(`pred_var_1'+`pred_var_0') if `exp'==1 & seq==1 & _mj==`i'
				drop pred_prob
				}
			
			summ imp_risk_ref, det
			local imp_risk_ref=r(mean)
			
			summ imp_rd, det
			local rd_pool=r(mean)
			
			*Now replace 35 values of imp_rd_diff_sq (as there will be 35 non-missing values of imp_rd)
			replace imp_rd_diff_sq=(imp_rd-`rd_pool')^2
			egen totsq = total(imp_rd_diff_sq)
			*Divide by n-1=35-1
			local vb = totsq/34
			replace imp_rd_se_sq=(imp_rd_se)^2
			egen totse = total(imp_rd_se_sq)
			*Divide by n=35
			local vw = totse/35
			local rd_se_pool=`vw'+`vb'+(`vb'/35)
			
			local lrd_pool=`rd_pool'-1.96*`rd_se_pool'
			local urd_pool=`rd_pool'+1.96*`rd_se_pool'
			drop seq totsq totse

			local x=`x'+1
			putexcel A`x'="`exp'" B`x'="`out'" C`x'="Yes minus prenatal" D`x'=`beta' E`x'=`se' F`x'=`pvalue' G`x'="`s'" H`x'=`prop_no'*100 I`x'=`prop_yes'*100 J`x'=`rr' K`x'=`lci' L`x'=`uci' M`x'=`paf' N`x'=`lpaf' O`x'=`upaf' P`x'=`pc'*100 Q`x'=`prev'*100 R`x'=`rd_pool' S`x'=`lrd_pool' T`x'=`urd_pool' U`x'=`imp_risk_ref'			
			
			*Adjusted model A
			xi: mim: glm `out' `exp' marstat matage birthwgt parity smokpreg i.highested i.hhsocclas preg_mrep, family(poisson) link(log) vce(robust)
			mim, storebv
			
			local beta = _b[`exp']
			local se = _se[`exp']
			
			test `exp'		
			local pvalue = `r(p)'
			
			local rr=2.718^(`beta')
			local lci=2.718^(`beta'-1.96*`se')
			local uci=2.718^(`beta'+1.96*`se')
			
			*PAFs
			local paf=`pc'*(1-(1/`rr'))*100
			local lpaf=`pc'*(1-(1/`lci'))*100
			local upaf=`pc'*(1-(1/`uci'))*100
		
			*Predicted probs
			*To average later
			bysort _mj `exp': egen seq = seq() if _mj != 0 

			forvalues i=1/35{
				logistic `out' `exp' marstat matage birthwgt parity smokpreg i.highested i.hhsocclas preg_mrep ///
				if _mj==`i'
				*Risks of outcome overall within current imputation i (different estimated risks across the individuals)
				predict pred_prob if _mj==`i', pr
				
				summ pred_prob if `exp'==0 & _mj==`i', det
				*Average risk (and variance) of outcome amongst the unexposed (save this for later)
				local pred_prob_0 = r(mean)
				*local pred_se_0 = r(sd)/(sqrt(N))
				local pred_var_0 = r(Var)
				
				summ pred_prob if `exp'==1 & _mj==`i', det
				local pred_prob_1 = r(mean)
				local pred_var_1 = r(Var)
				
				replace imp_risk_ref = `pred_prob_0' if `exp'==1 & seq==1 & _mj==`i'
				replace imp_rd = `pred_prob_1'-`pred_prob_0' if `exp'==1 & seq==1 & _mj==`i'
				*prob1 and prob0 have no covariance (as probability of same outcome, conditional on exposure status. Noone can be both)
				replace imp_rd_se = sqrt(`pred_var_1'+`pred_var_0') if `exp'==1 & seq==1 & _mj==`i'
				drop pred_prob
				}
			
			summ imp_risk_ref, det
			local imp_risk_ref=r(mean)
			
			summ imp_rd, det
			local rd_pool=r(mean)
			
			replace imp_rd_diff_sq=(imp_rd-`rd_pool')^2
			egen totsq = total(imp_rd_diff_sq)
			*Divide by n-1=35-1
			local vb = totsq/34
			replace imp_rd_se_sq=(imp_rd_se)^2
			egen totse = total(imp_rd_se_sq)
			*Divide by n=35
			local vw = totse/35
			local rd_se_pool=`vw'+`vb'+(`vb'/35)
			
			local lrd_pool=`rd_pool'-1.96*`rd_se_pool'
			local urd_pool=`rd_pool'+1.96*`rd_se_pool'
			drop seq
			
			local x=`x'+1
			putexcel A`x'="`exp'" B`x'="`out'" C`x'="Yes" D`x'=`beta' E`x'=`se' F`x'=`pvalue' G`x'="`s'" H`x'=`prop_no'*100 I`x'=`prop_yes'*100 J`x'=`rr' K`x'=`lci' L`x'=`uci' M`x'=`paf' N`x'=`lpaf' O`x'=`upaf' P`x'=`pc'*100 Q`x'=`prev'*100 R`x'=`rd_pool' S`x'=`lrd_pool' T`x'=`urd_pool' U`x'=`imp_risk_ref'												  	
		
			*Adjusted for ACEs
			xi: mim: glm `out' `exp' marstat matage birthwgt parity smokpreg i.highested i.hhsocclas preg_mrep ///
			i.bullying_0_1 i.emotional_ab i.emotional_ne i.physical_abu i.sexual_abuse i.mentl_hlth_p i.substance_ho i.parent_convi, family(poisson) link(log) vce(robust)
			mim, storebv
			
			local beta = _b[`exp']
			local se = _se[`exp']
			
			test `exp'		
			local pvalue = `r(p)'
			
			local rr=2.718^(`beta')
			local lci=2.718^(`beta'-1.96*`se')
			local uci=2.718^(`beta'+1.96*`se')
			local paf=`pc'*(1-(1/`rr'))*100
			local lpaf=`pc'*(1-(1/`lci'))*100
			local upaf=`pc'*(1-(1/`uci'))*100

			local x=`x'+1
			putexcel A`x'="`exp'" B`x'="`out'" C`x'="Yes + ACEs" ///
			D`x'=`beta' E`x'=`se' F`x'=`pvalue' G`x'="`s'" H`x'=`prop_no'*100 I`x'=`prop_yes'*100 J`x'=`rr' K`x'=`lci' L`x'=`uci' M`x'=`paf' N`x'=`lpaf' O`x'=`upaf' P`x'=`pc'*100 Q`x'=`prev'*100
			
			*Specific adjustment with same covariates as we have in the joint PAFs analysis (next section)
			*Binary household social class:
			cap recode hhsocclas (1 2 3 = 0) (4 5 = 1), gen(hhsocclas_bin)
			
			xi: mim: glm `out' `exp' marstat matage smokpreg hhsocclas_bin preg_mrep, family(poisson) link(log) vce(robust)
			mim, storebv
			
			local beta = _b[`exp']
			local se = _se[`exp']
			
			test `exp'		
			local pvalue = `r(p)'
			
			local rr=2.718^(`beta')
			local lci=2.718^(`beta'-1.96*`se')
			local uci=2.718^(`beta'+1.96*`se')
			local paf=`pc'*(1-(1/`rr'))*100
			local lpaf=`pc'*(1-(1/`lci'))*100
			local upaf=`pc'*(1-(1/`uci'))*100
			
			local x=`x'+1
			putexcel A`x'="`exp'" B`x'="`out'" C`x'="Yes (reduced set)" ///
			D`x'=`beta' E`x'=`se' F`x'=`pvalue' G`x'="`s'" H`x'=`prop_no'*100 I`x'=`prop_yes'*100 J`x'=`rr' K`x'=`lci' L`x'=`uci' M`x'=`paf' N`x'=`lpaf' O`x'=`upaf' P`x'=`pc'*100 Q`x'=`prev'*100
			
			*Restricted exposure
			use "imputed_sample_bysex_`exp'_`out'.dta", replace
			keep if kz021==`s'

			xi: mim: glm `out' `exp' marstat birthwgt parity i.highested matage smokpreg hhsocclas_bin preg_mrep, family(poisson) link(log) vce(robust)
			mim, storebv
			
			local beta = _b[`exp']
			local se = _se[`exp']
			
			test `exp'		
			local pvalue = `r(p)'
			
			local rr=2.718^(`beta')
			local lci=2.718^(`beta'-1.96*`se')
			local uci=2.718^(`beta'+1.96*`se')

			local x=`x'+1
			putexcel A`x'="`exp'" B`x'="`out'" C`x'="Yes (restricted exp)" ///
			D`x'=`beta' E`x'=`se' F`x'=`pvalue' G`x'="`s'" J`x'=`rr' K`x'=`lci' L`x'=`uci'
			}
		}
	}

timer off 2
timer list


********************************************************************************
*Now PAFs for other factors (only adjustment A needed)
*Noting that now household social class binarized and no longer adjusting for birthweight, parity, highest education
*Joint PAFs, based on Rockhill et al, 1998
*When several risk factors are being considered simultaneously, the exposure categories arise from a complete cross classification of the risk factors under consideration.
timer on 3

putexcel set "tables", sheet("Joint PAFs") modify
putexcel A1="Joint exposure level" B1="Outcome" C1="Adjusted" D1="beta" E1="SE" F1="P Value" G1="Sex" H1="no_perc" I1="yes_perc" J1="RR" K1="LCI" L1="UCI" M1="PAF" N1="LPAF" O1="UPAF" P1="pc" Q1="Exposure 1" R1="Exposure 2" S1="prev" T1="RD" U1="RD_SE" V1="LRD" W1="URD" X1="Risk ref"

global exposure1 "dv_mrep psych_mrep cc_mrep phys_mrep"
global exposure2 "maltreat aces3"
global outcome "vic_1821 per_1821" 

local x=1
forvalues s=1/2{
	foreach exp1 of global exposure1{
		foreach exp2 of global exposure2{
			foreach out of global outcome{	
				use "imputed_sample_bysex_dv_mrep_`out'.dta", clear
				keep if kz021==`s'		
				
				*Note that each 'interaction' term is gp_`v'_`w' where `v' is dv_mrep, etc., and `w' is maltreat, restoftrio, and aces3
				*0 = "Neither or just maltreat/rest/aces3", 1 = "dv_mrep without maltreat/rest/aces3", 2 = "Both"
				*gen group = gp_`exp1'_`exp2'
				
				egen group2 = group(`exp1' `exp2')
				*Labels: group2 1 "Neither" 2 "0 `exp1', 1 `exp2'" 3 "1 `exp1' 0 `exp2'" 4 "Both"
				recode group2 (1 2 = 0) (3 = 1) (4 = 2), gen(group)
				label define group 0 "No DVA" 1 "DVA only" 2 "Both"
				label values group group

				*Binary household social class:
				cap recode hhsocclas (1 2 3 = 0) (4 5 = 1), gen(hhsocclas_bin)
		
				xi: mim: glm `out' i.group marstat matage smokpreg hhsocclas_bin preg_mrep, family(poisson) link(log) vce(robust)
								
				mim, storebv
				
				*Coefs, SEs, ps
				forvalues i=1/2{
					local beta`i' = _b[_Igroup_`i']
					local se`i' = _se[_Igroup_`i']
					test _Igroup_`i'
					local p`i' = `r(p)'
					}
				
				*Proportion each level with outcome
				bysort _mj group: egen denom = count(aln) if _mj != 0
				bysort _mj group: egen num2 = count(aln) if _mj != 0 & `out'==1
				bysort _mj group: egen num = min(num2)
				gen prop = num/denom if _mj != 0
				bysort _mj group: egen seq = seq() if _mj != 0 
				
				forvalues i=0/2{
					summ prop if seq == 1 & group==`i'
					local prop_`i' = r(mean)
					}
				drop denom num num2 prop* seq
				
				*Proportion of exposure among cases
				forvalues i=1/2{
					bysort _mj `out': egen denom = count(aln) if _mj != 0
					bysort _mj `out': egen num2 = count(aln) if _mj != 0 & group==`i'
					bysort _mj `out': egen num = min(num2)
					gen prop = num/denom if _mj != 0
					bysort _mj `out': egen seq = seq() if _mj != 0 
					summ prop if seq == 1 & `out'==1
					local pc`i' = r(mean)
					drop denom num num2 prop seq

					local rr`i'=2.718^(`beta`i'')
					local lci`i'=2.718^(`beta`i''-1.96*`se`i'')
					local uci`i'=2.718^(`beta`i''+1.96*`se`i'')
					local paf`i'=`pc`i''*(1-(1/`rr`i''))*100
					local lpaf`i'=`pc`i''*(1-(1/`lci`i''))*100
					local upaf`i'=`pc`i''*(1-(1/`uci`i''))*100
					
					bysort _mj: egen denom = count(aln) if _mj != 0
					bysort _mj: egen num2 = count(aln) if _mj != 0 & group==`i'
					bysort _mj: egen num = min(num2)
					gen prop = num/denom if _mj != 0
					bysort _mj: egen seq = seq() if _mj != 0 
					summ prop if seq == 1
					local prev`i' = r(mean)*100
					drop denom num num2 prop seq
					}
							
				*Risk differences, based on: 
				xi: mim: glm `out' i.group marstat birthwgt parity i.highested matage smokpreg hhsocclas_bin preg_mrep, family(binomial) link(identity)
				mim, storebv

				forvalues i=1/2{
					local rd`i'=_b[_Igroup_`i']
					local rd_se`i' = _se[_Igroup_`i']
					local rd_lrd`i'=`rd`i''-1.96*`rd_se`i''
					local rd_urd`i'=`rd`i''+1.96*`rd_se`i''
					}
				
				local x=`x'+1
				putexcel A`x'="`exp1'=1 `exp2'=0" B`x'="`out'" C`x'="Yes" D`x'=`beta1' E`x'=`se1' F`x'=`p1' G`x'="`s'" H`x'=1-`prop_1' I`x'=`prop_1' J`x'=`rr1' K`x'=`lci1' L`x'=`uci1' M`x'=`paf1' N`x'=`lpaf1' O`x'=`upaf1' P`x'=`pc1' Q`x'="`exp1'" R`x'="`exp2'" S`x'=`prev1' T`x'=`rd1' U`x'=`rd_se1' V`x'=`rd_lrd1' W`x'=`rd_urd1' X`x'=`prop_0'
				 
				local x=`x'+1
				putexcel A`x'="`exp1'=1 `exp2'=1" B`x'="`out'" C`x'="Yes" D`x'=`beta2' E`x'=`se2' F`x'=`p2' G`x'="`s'" H`x'=1-`prop_2' I`x'=`prop_2' J`x'=`rr2' K`x'=`lci2' L`x'=`uci2' M`x'=`paf2' N`x'=`lpaf2' O`x'=`upaf2' P`x'=`pc2' Q`x'="`exp1'" R`x'="`exp2'" S`x'=`prev2' T`x'=`rd2' U`x'=`rd_se2' V`x'=`rd_lrd2' W`x'=`rd_urd2' X`x'=`prop_0'
				}		
			}
		}
	}

*Have to do rest of trio separately because in per_1821 as we cannot estimate for psych+restoftrio and cc+restoftrio
global exposure1 "dv_mrep psych_mrep cc_mrep phys_mrep"
global exposure2 "restoftrio"
global outcome "vic_1821" 
forvalues s=1/2{
	foreach exp1 of global exposure1{
		foreach exp2 of global exposure2{
			foreach out of global outcome{	
				use "imputed_sample_bysex_dv_mrep_`out'.dta", clear
				keep if kz021==`s'		
				
				egen group2 = group(`exp1' `exp2')
				*Labels: group2 1 "Neither" 2 "0 `exp1', 1 `exp2'" 3 "1 `exp1' 0 `exp2'" 4 "Both"
				recode group2 (1 2 = 0) (3 = 1) (4 = 2), gen(group)
				label define group 0 "No DVA" 1 "DVA only" 2 "Both"
				label values group group

				cap recode hhsocclas (1 2 3 = 0) (4 5 = 1), gen(hhsocclas_bin)

				xi: mim: glm `out' i.group marstat matage smokpreg hhsocclas_bin preg_mrep, family(poisson) link(log)
				mim, storebv
				
				*Coefs, SEs, ps
				forvalues i=1/2{
					local beta`i' = _b[_Igroup_`i']
					local se`i' = _se[_Igroup_`i']
					test _Igroup_`i'
					local p`i' = `r(p)'
					}
				
				*Proportion each level with outcome
				bysort _mj group: egen denom = count(aln) if _mj != 0
				bysort _mj group: egen num2 = count(aln) if _mj != 0 & `out'==1
				bysort _mj group: egen num = min(num2)
				gen prop = num/denom if _mj != 0
				bysort _mj group: egen seq = seq() if _mj != 0 
				
				forvalues i=0/2{
					summ prop if seq == 1 & group==`i'
					local prop_`i' = r(mean)
					}
				drop denom num num2 prop* seq
				
				*Proportion of exposure among cases
				forvalues i=1/2{
					bysort _mj `out': egen denom = count(aln) if _mj != 0
					bysort _mj `out': egen num2 = count(aln) if _mj != 0 & group==`i'
					bysort _mj `out': egen num = min(num2)
					gen prop = num/denom if _mj != 0
					bysort _mj `out': egen seq = seq() if _mj != 0 
					summ prop if seq == 1 & `out'==1
					local pc`i' = r(mean)
					drop denom num num2 prop seq

					local rr`i'=2.718^(`beta`i'')
					local lci`i'=2.718^(`beta`i''-1.96*`se`i'')
					local uci`i'=2.718^(`beta`i''+1.96*`se`i'')
					local paf`i'=`pc`i''*(1-(1/`rr`i''))*100
					local lpaf`i'=`pc`i''*(1-(1/`lci`i''))*100
					local upaf`i'=`pc`i''*(1-(1/`uci`i''))*100
					
					bysort _mj: egen denom = count(aln) if _mj != 0
					bysort _mj: egen num2 = count(aln) if _mj != 0 & group==`i'
					bysort _mj: egen num = min(num2)
					gen prop = num/denom if _mj != 0
					bysort _mj: egen seq = seq() if _mj != 0 
					summ prop if seq == 1
					local prev`i' = r(mean)*100
					drop denom num num2 prop seq
					}
							
				*Risk differences, based on: 
				xi: mim: glm `out' i.group marstat birthwgt parity i.highested matage smokpreg hhsocclas_bin preg_mrep, family(binomial) link(identity)
				mim, storebv

				forvalues i=1/2{
					local rd`i'=_b[_Igroup_`i']
					local rd_se`i' = _se[_Igroup_`i']
					local rd_lrd`i'=`rd`i''-1.96*`rd_se`i''
					local rd_urd`i'=`rd`i''+1.96*`rd_se`i''
					}
				
				local x=`x'+1
				putexcel A`x'="`exp1'=1 `exp2'=0" B`x'="`out'" C`x'="Yes" D`x'=`beta1' E`x'=`se1' F`x'=`p1' G`x'="`s'" H`x'=1-`prop_1' I`x'=`prop_1' J`x'=`rr1' K`x'=`lci1' L`x'=`uci1' M`x'=`paf1' N`x'=`lpaf1' O`x'=`upaf1' P`x'=`pc1' Q`x'="`exp1'" R`x'="`exp2'" S`x'=`prev1' T`x'=`rd1' U`x'=`rd_se1' V`x'=`rd_lrd1' W`x'=`rd_urd1' X`x'=`prop_0'
				 
				local x=`x'+1
				putexcel A`x'="`exp1'=1 `exp2'=1" B`x'="`out'" C`x'="Yes" D`x'=`beta2' E`x'=`se2' F`x'=`p2' G`x'="`s'" H`x'=1-`prop_2' I`x'=`prop_2' J`x'=`rr2' K`x'=`lci2' L`x'=`uci2' M`x'=`paf2' N`x'=`lpaf2' O`x'=`upaf2' P`x'=`pc2' Q`x'="`exp1'" R`x'="`exp2'" S`x'=`prev2' T`x'=`rd2' U`x'=`rd_se2' V`x'=`rd_lrd2' W`x'=`rd_urd2' X`x'=`prop_0'
				}		
			}
		}
	}
	
global exposure1 "dv_mrep phys_mrep"
global exposure2 "restoftrio"
global outcome "per_1821" 
forvalues s=1/2{
	foreach exp1 of global exposure1{
		foreach exp2 of global exposure2{
			foreach out of global outcome{	
				use "imputed_sample_bysex_dv_mrep_`out'.dta", clear
				keep if kz021==`s'		
				
				egen group2 = group(`exp1' `exp2')
				*Labels: group2 1 "Neither" 2 "0 `exp1', 1 `exp2'" 3 "1 `exp1' 0 `exp2'" 4 "Both"
				recode group2 (1 2 = 0) (3 = 1) (4 = 2), gen(group)
				label define group 0 "No DVA" 1 "DVA only" 2 "Both"
				label values group group

				cap recode hhsocclas (1 2 3 = 0) (4 5 = 1), gen(hhsocclas_bin)
		
				xi: mim: glm `out' i.group marstat matage smokpreg hhsocclas_bin preg_mrep, family(poisson) link(log) vce(robust)
				mim, storebv
				
				*Coefs, SEs, ps
				forvalues i=1/2{
					local beta`i' = _b[_Igroup_`i']
					local se`i' = _se[_Igroup_`i']
					test _Igroup_`i'
					local p`i' = `r(p)'
					}
				
				*Proportion each level with outcome
				bysort _mj group: egen denom = count(aln) if _mj != 0
				bysort _mj group: egen num2 = count(aln) if _mj != 0 & `out'==1
				bysort _mj group: egen num = min(num2)
				gen prop = num/denom if _mj != 0
				bysort _mj group: egen seq = seq() if _mj != 0 
				
				forvalues i=0/2{
					summ prop if seq == 1 & group==`i'
					local prop_`i' = r(mean)
					}
				drop denom num num2 prop* seq
				
				*Proportion of exposure among cases
				forvalues i=1/2{
					bysort _mj `out': egen denom = count(aln) if _mj != 0
					bysort _mj `out': egen num2 = count(aln) if _mj != 0 & group==`i'
					bysort _mj `out': egen num = min(num2)
					gen prop = num/denom if _mj != 0
					bysort _mj `out': egen seq = seq() if _mj != 0 
					summ prop if seq == 1 & `out'==1
					local pc`i' = r(mean)
					drop denom num num2 prop seq

					local rr`i'=2.718^(`beta`i'')
					local lci`i'=2.718^(`beta`i''-1.96*`se`i'')
					local uci`i'=2.718^(`beta`i''+1.96*`se`i'')
					local paf`i'=`pc`i''*(1-(1/`rr`i''))*100
					local lpaf`i'=`pc`i''*(1-(1/`lci`i''))*100
					local upaf`i'=`pc`i''*(1-(1/`uci`i''))*100
					
					bysort _mj: egen denom = count(aln) if _mj != 0
					bysort _mj: egen num2 = count(aln) if _mj != 0 & group==`i'
					bysort _mj: egen num = min(num2)
					gen prop = num/denom if _mj != 0
					bysort _mj: egen seq = seq() if _mj != 0 
					summ prop if seq == 1
					local prev`i' = r(mean)*100
					drop denom num num2 prop seq
					}
							
				*Risk differences, based on: 
				xi: mim: glm `out' i.group marstat birthwgt parity i.highested matage smokpreg hhsocclas_bin preg_mrep, family(binomial) link(identity)
				mim, storebv

				forvalues i=1/2{
					local rd`i'=_b[_Igroup_`i']
					local rd_se`i' = _se[_Igroup_`i']
					local rd_lrd`i'=`rd`i''-1.96*`rd_se`i''
					local rd_urd`i'=`rd`i''+1.96*`rd_se`i''
					}
				
				local x=`x'+1
				putexcel A`x'="`exp1'=1 `exp2'=0" B`x'="`out'" C`x'="Yes" D`x'=`beta1' E`x'=`se1' F`x'=`p1' G`x'="`s'" H`x'=1-`prop_1' I`x'=`prop_1' J`x'=`rr1' K`x'=`lci1' L`x'=`uci1' M`x'=`paf1' N`x'=`lpaf1' O`x'=`upaf1' P`x'=`pc1' Q`x'="`exp1'" R`x'="`exp2'" S`x'=`prev1' T`x'=`rd1' U`x'=`rd_se1' V`x'=`rd_lrd1' W`x'=`rd_urd1' X`x'=`prop_0'
				 
				local x=`x'+1
				putexcel A`x'="`exp1'=1 `exp2'=1" B`x'="`out'" C`x'="Yes" D`x'=`beta2' E`x'=`se2' F`x'=`p2' G`x'="`s'" H`x'=1-`prop_2' I`x'=`prop_2' J`x'=`rr2' K`x'=`lci2' L`x'=`uci2' M`x'=`paf2' N`x'=`lpaf2' O`x'=`upaf2' P`x'=`pc2' Q`x'="`exp1'" R`x'="`exp2'" S`x'=`prev2' T`x'=`rd2' U`x'=`rd_se2' V`x'=`rd_lrd2' W`x'=`rd_urd2' X`x'=`prop_0'
				}		
			}
		}
	}
	
timer off 3
	

*******************************************************************************
*SUPP TABLE S3: descriptive stats
timer on 4

forvalues s=1/2{
	putexcel set "tables", sheet("ST3_characs_`s'") modify
	putexcel A1="COHORT" A2="variable" B2="label" C2="n" D2="p" E2="extra" F1="ALSPAC" F2="variable" G2="label" H2="n" I2="p" J2="extra"

	local x=2
	use "imputed_sample_bysex_dv_mrep_vic_1821.dta", replace

	keep if kz021==`s'
	summ _mj
	local nimp=r(max)

	*Continuous vars
	foreach v of varlist matage birthwgt{ 
		local med_tot=0
		local lq_tot=0
		local uq_tot=0
		forvalues j=1/`nimp'{
			summ `v' if _mj==`j' & `v'>=0 & `v'!=., det
			local median = `r(p50)'
			local lq = `r(p25)'
			local uq = `r(p75)'
			local med_tot = `med_tot'+`median'
			local lq_tot = `lq_tot'+`lq'
			local uq_tot = `uq_tot'+`uq'
			}
		local x=`x'+1
		putexcel A`x'="`v'" B`x'="." C`x'=(`med_tot'/`nimp') D`x'=(`lq_tot'/`nimp') E`x'=(`uq_tot'/`nimp')
		}

	local var_parent "marstat smokpreg parity mated pated highested hhsocclas mentl_hlth_p substance_ho parent_convi parental_sep"
	local var_offspring "kz021 ethnicity maltreat emotional_ab emotional_ne physical_abu sexual_abu bullying_0_1 restoftrio aces3 aces4 dv_3plus aces4_dv_mal aces4_trio dv_mal dv_trio"
	
	*Sum
	egen sum=rowtotal(emotional_ne emotional_ab physical_abu bullying_0_1 substance_ho mentl_hlth_p parent_convi parental_sep sexual_abuse)
	*Already have 3+ ACEs variable
	
	*Generate 4+ ACEs (any)
	cap gen aces4=sum
	replace aces4=0 if sum<4
	replace aces4=1 if sum>=4 & sum!=.
	*Parental IPVA + 3+ other ACEs (i.e. 4+ including IPVA)
	cap gen dv_3plus=0 if aces4==0 | dv_mrep==0
	replace dv_3plus=1 if aces4==1 & dv_mrep==1
	*4+ ACEs where two are DV and maltreat
	cap gen aces4_dv_mal=0 if aces4==0 | dv_mrep==0 | maltreat==0
	replace aces4_dv_mal=1 if aces4==1 & dv_mrep==1 & maltreat==1
	*4+ ACEs where three are toxic trio
	cap gen aces4_trio=0 if aces4==0 | dv_mrep==0 | substance_ho==0 | mentl_hlth_p==0
	replace aces4_trio=1 if aces4==1 & dv_mrep==1 & substance_ho==1 & mentl_hlth_p==1
		
	*Parental IPVA + maltreatment
	cap gen dv_mal=0 if dv_mrep==0 | maltreat==0
	replace dv_mal=1 if dv_mrep==1 & maltreat==1

	*Parental IPVA + rest of toxic trio (i.e. entire trio)
	cap gen dv_trio=0 if dv_mrep==0 | restoftrio==0
	replace dv_trio=1 if dv_mrep==1 & restoftrio==1
	
	local var "`var_parent' `var_offspring'" 

	foreach v of local var{ 
		*replace `v'=-99 if `v'==.
		tab `v' if `v'>=0 & `v'!=. & _mj==1, m
		local total=r(N)
		levelsof `v', local(levels)
		foreach i of local levels{
			local num_tot=0
			forvalues j=1/`nimp'{
				tab `v' if `v'==`i' & `v'>=0 & `v'!=. & _mj==`j', matcell(matrix)
				local num_tot = `num_tot'+matrix[1,1]
				}
			local x=`x'+1
			putexcel A`x'="`v'" B`x'=`i' C`x'=`num_tot'/`nimp' D`x'=(`num_tot'/(`nimp'*`total')*100) E`x'=`total'
			}
		}

	*Have to open new dataset each time for reporting DV %s in Table 1
	local var_dv "dv_mrep phys_mrep psych_mrep cc_mrep preg_mrep" 

	foreach v of local var_dv{ 
		use "imputed_sample_bysex_`v'_vic_1821.dta", replace
		cap rename sex kz021
		keep if kz021==`s'
		*replace `v'=-99 if `v'==.
		tab `v' if `v'>=0 & `v'!=. & _mj==1, m
		local total=r(N)
		levelsof `v', local(levels)
		foreach i of local levels{
			local num_tot=0
			forvalues j=1/`nimp'{
				tab `v' if `v'==`i' & `v'>=0 & `v'!=. & _mj==`j', matcell(matrix)
				local num_tot = `num_tot'+matrix[1,1]
				}
			local x=`x'+1
			putexcel A`x'="`v'" B`x'=`i' C`x'=`num_tot'/`nimp' D`x'=(`num_tot'/(`nimp'*`total')*100) E`x'=`total'
			}
		}

	local var_ipva "vic per"
	foreach v of local var_ipva{ 
		use "imputed_sample_bysex_dv_mrep_`v'_1821.dta", replace
		keep if kz021==`s'
		*replace `v'=-99 if `v'==.
		tab `v'_1821 if `v'_1821>=0 & `v'_1821!=. & _mj==1, m
		local total=r(N)
		levelsof `v'_1821, local(levels)
		
		foreach i of local levels{
			local num_tot=0
			forvalues j=1/`nimp'{
				tab `v'_1821 if `v'_1821==`i' & `v'_1821>=0 & `v'_1821!=. & _mj==`j', matcell(matrix)
				local num_tot = `num_tot'+matrix[1,1]
				}
			local x=`x'+1
			putexcel A`x'="`v'_1821" B`x'=`i' C`x'=`num_tot'/`nimp' D`x'=(`num_tot'/(`nimp'*`total')*100) E`x'=`total'
		}
		
		*replace `v'_emo=-99 if `v'_emo==.
		tab `v'_emo_1821 if `v'_emo_1821>=0 & `v'_emo_1821!=. & _mj==1, m
		local total=r(N)
		levelsof `v'_emo_1821, local(levels)
		foreach i of local levels{
			local num_tot=0
			forvalues j=1/`nimp'{
				tab `v'_emo_1821 if `v'_emo_1821==`i' & `v'_emo_1821>=0 & `v'_emo_1821!=. & _mj==`j', matcell(matrix)
				local num_tot = `num_tot'+matrix[1,1]
				}
			local x=`x'+1
			putexcel A`x'="`v'_emo_1821" B`x'=`i' C`x'=`num_tot'/`nimp' D`x'=(`num_tot'/(`nimp'*`total')*100) E`x'=`total'
		}
		
		*replace `v'_emo_co=-99 if `v'_emo_co==.
		cap tab `v'_emo_co_1821 if `v'_emo_co_1821>=0 & `v'_emo_co_1821!=. & _mj==1, m
		local total=r(N)
		levelsof `v'_emo_co_1821, local(levels)
		foreach i of local levels{
			local num_tot=0
			forvalues j=1/`nimp'{
				cap tab `v'_emo_co_1821 if `v'_emo_co_1821==`i' & `v'_emo_co_1821>=0 & `v'_emo_co_1821!=. & _mj==`j', matcell(matrix)
				local num_tot = `num_tot'+matrix[1,1]
				}
			local x=`x'+1
			putexcel A`x'="`v'_emo_co_1821" B`x'=`i' C`x'=`num_tot'/`nimp' D`x'=(`num_tot'/(`nimp'*`total')*100) E`x'=`total'
		}
		
		*replace `v'_phys=-99 if `v'_phys==.
		tab `v'_phys_1821 if `v'_phys_1821>=0 & `v'_phys_1821!=. & _mj==1, m
		local total=r(N)
		levelsof `v'_phys_1821, local(levels)
		foreach i of local levels{
			local num_tot=0
			forvalues j=1/`nimp'{
				tab `v'_phys_1821 if `v'_phys_1821==`i' & `v'_phys_1821>=0 & `v'_phys_1821!=. & _mj==`j', matcell(matrix)
				local num_tot = `num_tot'+matrix[1,1]
				}
			local x=`x'+1
			putexcel A`x'="`v'_phys_1821" B`x'=`i' C`x'=`num_tot'/`nimp' D`x'=(`num_tot'/(`nimp'*`total')*100) E`x'=`total'
		}
		
		*replace `v'_sex=-99 if `v'_sex==.
		tab `v'_sex_1821 if `v'_sex_1821>=0 & `v'_sex_1821!=. & _mj==1, m
		local total=r(N)
		levelsof `v'_sex_1821, local(levels)
		foreach i of local levels{
			local num_tot=0
			forvalues j=1/`nimp'{
				tab `v'_sex_1821 if `v'_sex_1821==`i' & `v'_sex_1821>=0 & `v'_sex_1821!=. & _mj==`j', matcell(matrix)
				local num_tot = `num_tot'+matrix[1,1]
				}
			local x=`x'+1
			putexcel A`x'="`v'_sex_1821" B`x'=`i' C`x'=`num_tot'/`nimp' D`x'=(`num_tot'/(`nimp'*`total')*100) E`x'=`total'
			}
		}

	*ALSPAC
	local var_parent "marstat smokpreg parity mated pated highested hhsocclas mentl_hlth_p substance_ho parent_convi parental_sep"
	local var_offspring "kz021 ethnicity maltreat emotional_ab emotional_ne physical_abu sexual_abu bullying_0_1 restoftrio aces3 aces4 dv_3plus aces4_dv_mal aces4_trio dv_mal dv_trio"
	
	putexcel set "tables", sheet("ST3_characs_`s'") modify
	local x=2
	use  "\\rdsfcifs.acrc.bris.ac.uk\MRC-IEU-research\projects\ieu2\p6\046\working\data\cohort\p3 - intergen\imputed_alspac_bysex_dv_mrep_vic_1821.dta", replace
	keep if kz021==`s'
	summ _mj
	local nimp=r(max)
	foreach v of varlist matage birthwgt{ 
		local med_tot=0
		local lq_tot=0
		local uq_tot=0
		forvalues j=1/`nimp'{
			summ `v' if _mj==`j' & `v'>=0 & `v'!=., det
			local median = `r(p50)'
			local lq = `r(p25)'
			local uq = `r(p75)'
			local med_tot = `med_tot'+`median'
			local lq_tot = `lq_tot'+`lq'
			local uq_tot = `uq_tot'+`uq'
			}
		local x=`x'+1
		putexcel F`x'="`v'" G`x'="." H`x'=(`med_tot'/`nimp') I`x'=(`lq_tot'/`nimp') J`x'=(`uq_tot'/`nimp')
		}

	*Sum
	egen sum=rowtotal(emotional_ne emotional_ab physical_abu bullying_0_1 substance_ho mentl_hlth_p parent_convi parental_sep sexual_abuse)
	*Already have 3+ ACEs variable
	
	*Generate 4+ ACEs (any)
	cap gen aces4=sum
	replace aces4=0 if sum<4
	replace aces4=1 if sum>=4 & sum!=.
	*Parental IPVA + 3+ other ACEs (i.e. 4+ including IPVA)
	cap gen dv_3plus=0 if aces4==0 | dv_mrep==0
	replace dv_3plus=1 if aces4==1 & dv_mrep==1
	*4+ ACEs where two are DV and maltreat
	cap gen aces4_dv_mal=0 if aces4==0 | dv_mrep==0 | maltreat==0
	replace aces4_dv_mal=1 if aces4==1 & dv_mrep==1 & maltreat==1
	*4+ ACEs where three are toxic trio
	cap gen aces4_trio=0 if aces4==0 | dv_mrep==0 | substance_ho==0 | mentl_hlth_p==0
	replace aces4_trio=1 if aces4==1 & dv_mrep==1 & substance_ho==1 & mentl_hlth_p==1
		
	*Parental IPVA + maltreatment
	cap gen dv_mal=0 if dv_mrep==0 | maltreat==0
	replace dv_mal=1 if dv_mrep==1 & maltreat==1

	*Parental IPVA + rest of toxic trio (i.e. entire trio)
	cap gen dv_trio=0 if dv_mrep==0 | restoftrio==0
	replace dv_trio=1 if dv_mrep==1 & restoftrio==1
	
	local var "`var_parent' `var_offspring'" 

	foreach v of local var{ 
		*replace `v'=-99 if `v'==.
		tab `v' if `v'>=0 & `v'!=. & _mj==1, m
		local total=r(N)
		levelsof `v', local(levels)
		foreach i of local levels{
			local num_tot=0
			forvalues j=1/`nimp'{
				tab `v' if `v'==`i' & `v'>=0 & `v'!=. & _mj==`j', matcell(matrix)
				local num_tot = `num_tot'+matrix[1,1]
				}
			local x=`x'+1
			putexcel F`x'="`v'" G`x'=`i' H`x'=`num_tot'/`nimp' I`x'=(`num_tot'/(`nimp'*`total')*100) J`x'=`total'
			}
		}

	*Have to open new dataset each time for reporting DV %s in Table 1
	local var_dv "dv_mrep phys_mrep psych_mrep cc_mrep preg_mrep" 

	foreach v of local var_dv{ 
		use "imputed_alspac_bysex_`v'_vic_1821.dta", replace
		cap rename sex kz021 
		keep if kz021==`s'
		*replace `v'=-99 if `v'==.
		tab `v' if `v'>=0 & `v'!=. & _mj==1, m
		local total=r(N)
		levelsof `v', local(levels)
		foreach i of local levels{
			local num_tot=0
			forvalues j=1/`nimp'{
				tab `v' if `v'==`i' & `v'>=0 & `v'!=. & _mj==`j', matcell(matrix)
				local num_tot = `num_tot'+matrix[1,1]
				}
			local x=`x'+1
			putexcel F`x'="`v'" G`x'=`i' H`x'=`num_tot'/`nimp' I`x'=(`num_tot'/(`nimp'*`total')*100) J`x'=`total'
			}
		}

	local var_ipva "vic per"
	foreach v of local var_ipva{ 
		use "imputed_sample_bysex_dv_mrep_`v'_1821.dta", replace
		cap rename sex kz021 
		keep if kz021==`s'
		*replace `v'=-99 if `v'==.
		tab `v'_1821 if `v'_1821>=0 & `v'_1821!=. & _mj==1, m
		local total=r(N)
		levelsof `v'_1821, local(levels)
		
		foreach i of local levels{
			local num_tot=0
			forvalues j=1/`nimp'{
				tab `v'_1821 if `v'_1821==`i' & `v'_1821>=0 & `v'_1821!=. & _mj==`j', matcell(matrix)
				local num_tot = `num_tot'+matrix[1,1]
				}
			local x=`x'+1
			putexcel F`x'="`v'" G`x'=`i' H`x'=`num_tot'/`nimp' I`x'=(`num_tot'/(`nimp'*`total')*100) J`x'=`total'
		}
		
		*replace `v'_emo=-99 if `v'_emo==.
		tab `v'_emo_1821 if `v'_emo_1821>=0 & `v'_emo_1821!=. & _mj==1, m
		local total=r(N)
		levelsof `v'_emo_1821, local(levels)
		foreach i of local levels{
			local num_tot=0
			forvalues j=1/`nimp'{
				tab `v'_emo_1821 if `v'_emo_1821==`i' & `v'_emo_1821>=0 & `v'_emo_1821!=. & _mj==`j', matcell(matrix)
				local num_tot = `num_tot'+matrix[1,1]
				}
			local x=`x'+1
			putexcel F`x'="`v'_emo" G`x'=`i' H`x'=`num_tot'/`nimp' I`x'=(`num_tot'/(`nimp'*`total')*100) J`x'=`total'
		}
		
		*replace `v'_emo_co=-99 if `v'_emo_co==.
		cap tab `v'_emo_co_1821 if `v'_emo_co_1821>=0 & `v'_emo_co_1821!=. & _mj==1, m
		local total=r(N)
		levelsof `v'_emo_co_1821, local(levels)
		foreach i of local levels{
			local num_tot=0
			forvalues j=1/`nimp'{
				cap tab `v'_emo_co_1821 if `v'_emo_co_1821==`i' & `v'_emo_co_1821>=0 & `v'_emo_co_1821!=. & _mj==`j', matcell(matrix)
				local num_tot = `num_tot'+matrix[1,1]
				}
			local x=`x'+1
			putexcel F`x'="`v'_emo_co" G`x'=`i' H`x'=`num_tot'/`nimp' I`x'=(`num_tot'/(`nimp'*`total')*100) J`x'=`total'
		}
		
		*replace `v'_phys=-99 if `v'_phys==.
		tab `v'_phys_1821 if `v'_phys_1821>=0 & `v'_phys_1821!=. & _mj==1, m
		local total=r(N)
		levelsof `v'_phys_1821, local(levels)
		foreach i of local levels{
			local num_tot=0
			forvalues j=1/`nimp'{
				tab `v'_phys_1821 if `v'_phys_1821==`i' & `v'_phys_1821>=0 & `v'_phys_1821!=. & _mj==`j', matcell(matrix)
				local num_tot = `num_tot'+matrix[1,1]
				}
			local x=`x'+1
			putexcel F`x'="`v'_phys" G`x'=`i' H`x'=`num_tot'/`nimp' I`x'=(`num_tot'/(`nimp'*`total')*100) J`x'=`total'
		}
		
		*replace `v'_sex=-99 if `v'_sex==.
		tab `v'_sex_1821 if `v'_sex_1821>=0 & `v'_sex_1821!=. & _mj==1, m
		local total=r(N)
		levelsof `v'_sex_1821, local(levels)
		foreach i of local levels{
			local num_tot=0
			forvalues j=1/`nimp'{
				tab `v'_sex_1821 if `v'_sex_1821==`i' & `v'_sex_1821>=0 & `v'_sex_1821!=. & _mj==`j', matcell(matrix)
				local num_tot = `num_tot'+matrix[1,1]
				}
			local x=`x'+1
			putexcel F`x'="`v'_sex" G`x'=`i' H`x'=`num_tot'/`nimp' I`x'=(`num_tot'/(`nimp'*`total')*100) J`x'=`total'
			}
		}
	}

timer off 4
timer list


*SUPP TABLE S4
*Complete cases
local outcome "vic_1821 per_1821"
local exposure "dv_mrep psych_mrep cc_mrep phys_mrep"

putexcel set "tables", sheet("complete_case") modify
putexcel A1="Exposure" B1="Outcome" C1="RR" D1="LCI" E1="UCI" F1="N" G1="Sex"

local x=1
foreach out of local outcome{
	foreach exp of local exposure{
	use "preimputation.dta"

		forvalues i=1/2{
			glm `out' `exp' marstat matage birthwgt parity smokpreg i.highested i.hhsocclas preg_mrep if _mj==1 & kz021==`i', family(poisson) link(log) vce(robust)
					
			local beta = _b[`exp']
			local se = _se[`exp']
			
			local rr=2.718^(`beta')
			local lci=2.718^(`beta'-1.96*`se')
			local uci=2.718^(`beta'+1.96*`se')
			
			local total=e(N)
			
			local x=`x'+1
			putexcel A`x'="`exp'" B`x'="`out'" C`x'=`rr' D`x'=`lci' E`x'=`uci' F`x'=`total' G`x'=`i'
			}
		}
	}

*Proportions of missing data in exposures, outcomes, and covariates, prior to multiple imputation  
use "preimputation_update.dta", replace
keep if kz021==1
mdesc marstat mum_sexmin smokpreg matage parity highested hhsocclas dv_mrep phys_mrep psych_mrep cc_mrep ethnicity sexmin_ever birthwgt vic_1821 vic_phys_1821 vic_emo_1821 vic_emo_co_1821 vic_sex_1821 per_1821 per_phys_1821 per_emo_1821 per_emo_co_1821 per_sex_1821 maltreat emotional_ab emotional_ne physical_abu sexual_abu bullying_0_1 mentl_hlth_p substance_ho parent_convi parental_sep
*paste into "tables", sheet("complete_case_characs") using text import delimited
use "preimputation_update.dta", replace
keep if kz021==2
mdesc marstat mum_sexmin smokpreg matage parity highested hhsocclas dv_mrep phys_mrep psych_mrep cc_mrep ethnicity sexmin_ever birthwgt vic_1821 vic_phys_1821 vic_emo_1821 vic_emo_co_1821 vic_sex_1821 per_1821 per_phys_1821 per_emo_1821 per_emo_co_1821 per_sex_1821 maltreat emotional_ab emotional_ne physical_abu sexual_abu bullying_0_1 mentl_hlth_p substance_ho parent_convi parental_sep
*paste into "tables", sheet("complete_case_characs") using text import delimited


********************************************************************************	
*27.02.23
*ANALYSES FOR VALUES WITHIN MAIN TEXT
timer on 5

*Testing interaction of sex
use "imputed_sample_bysex_dv_mrep_vic_1821.dta", clear
xi: mim: glm vic_1821 dv_mrep##sex marstat matage birthwgt parity smokpreg i.highested i.hhsocclas preg_mrep, family(poisson) link(log) vce(robust)

use "imputed_sample_bysex_dv_mrep_per_1821.dta", clear
xi: mim: glm per_1821 dv_mrep##sex marstat matage birthwgt parity smokpreg i.highested i.hhsocclas preg_mrep, family(poisson) link(log) vce(robust)
*p= vic: 0.184 perp: 0.213
						
*Overlap between psych_mrep and phys_mrep
use "imputed_sample_bysex_general.dta", clear
tab psych_mrep phys_mrep if _mj==0, m 
bysort _mj: egen denom = count(aln) if _mj != 0
bysort _mj: egen num2 = count(aln) if _mj != 0 & psych_mrep==1 & phys_mrep==1
bysort _mj: egen num = min(num2)
local di `num'
bysort _mj: egen seq = seq() if _mj != 0 
summ num if seq == 1
local mean_num = r(mean)
drop denom num num2 seq
di `mean_num'
*Nov 2024: n=370 in non-imputed, 641 in imputed

*How many cases of DVA also have CC?
use "preimputation.dta", clear
tab dv_mrep cc_mrep, row m
tab dv_mrep psych_mrep, row m
tab dv_mrep phys_mrep, row m
*16%, 82% 51%, respectively.

*How many answering caregivers are the mother on their own (whether bio or not)?
tab mumonly, m
*n = 3143 (i.e. 97%)

*What proportion of 3+ ACEs are maltreatment?
tab aces3 maltreat, row m
*91%
				 
*91% of dv_mrep + 3+ other ACEs cases included maltreatment. See what the PAFs are when separating my maltreatment and not
use "imputed_sample_bysex_dv_mrep_per_1821.dta", clear
keep if sex == 1
egen group2 = group(dv_mrep aces3 maltreat)
*Labels: group2 1 "None at all" 2 "dv_mrep == 1 aces3 == 0 maltreat == 0" 3 "" 4 "" 5 "" 6 ""
recode group2 (1 2 = 0) (3 = 1) (4 = 2), gen(group)
label define group 0 "No DVA" 1 "DVA only" 2 "Both"
label values group group
cap recode hhsocclas (1 2 3 = 0) (4 5 = 1), gen(hhsocclas_bin)

xi: mim: glm `out' i.group marstat matage smokpreg hhsocclas_bin preg_mrep, family(poisson) link(log) vce(robust)
mim, storebv

putexcel set "tables", sheet("maltreat_aces_sens") modify
putexcel A1="Exposure" B1="Outcome" C1="Adjusted" D1="beta" E1="SE" F1="P Value" G1="Sex" H1="no_perc" I1="yes_perc" J1="RR" K1="LCI" L1="UCI" M1="PAF" N1="LPAF" O1="UPAF" P1="pc" Q1="prev" R1="RD" S1="LRD" T1="URD" U1="Risk ref"
local exp1 "dv_mrep"
local exp2 "aces3"
forvalues s = 1/2{
		use "imputed_sample_bysex_dv_mrep_`out'.dta", clear
		keep if kz021==`s'		
		keep if maltreat==`j'
		*Note that each 'interaction' term is gp_`v'_`w' where `v' is dv_mrep, etc., and `w' is maltreat, restoftrio, and aces3
		*0 = "Neither or just maltreat/rest/aces3", 1 = "dv_mrep without maltreat/rest/aces3", 2 = "Both"
		*gen group = gp_`exp1'_`exp2'
		
		egen group2 = group(`exp1' `exp2')
		*Labels: group2 1 "Neither" 2 "0 `exp1', 1 `exp2'" 3 "1 `exp1' 0 `exp2'" 4 "Both"
		recode group2 (1 2 = 0) (3 = 1) (4 = 2), gen(group)
		label define group 0 "No DVA" 1 "DVA only" 2 "Both"
		label values group group
		cap recode hhsocclas (1 2 3 = 0) (4 5 = 1), gen(hhsocclas_bin)
			
		xi: mim: glm `out' i.group marstat matage smokpreg hhsocclas_bin preg_mrep, family(binomial) link(log)
		mim, storebv
		
		*Coefs, SEs, ps
		forvalues i=1/2{
			local beta`i' = _b[_Igroup_`i']
			local se`i' = _se[_Igroup_`i']
			test _Igroup_`i'
			local p`i' = `r(p)'
			}
		
		*Proportion each level with outcome
		bysort _mj group: egen denom = count(aln) if _mj != 0
		bysort _mj group: egen num2 = count(aln) if _mj != 0 & `out'==1
		bysort _mj group: egen num = min(num2)
		gen prop = num/denom if _mj != 0
		bysort _mj group: egen seq = seq() if _mj != 0 
		
		forvalues i=0/2{
			summ prop if seq == 1 & group==`i'
			local prop_`i' = r(mean)
			}
		drop denom num num2 prop* seq
		
		*Proportion of exposure among cases
		forvalues i=1/2{
			bysort _mj `out': egen denom = count(aln) if _mj != 0
			bysort _mj `out': egen num2 = count(aln) if _mj != 0 & group==`i'
			bysort _mj `out': egen num = min(num2)
			gen prop = num/denom if _mj != 0
			bysort _mj `out': egen seq = seq() if _mj != 0 
			summ prop if seq == 1 & `out'==1
			local pc`i' = r(mean)
			drop denom num num2 prop seq

			local rr`i'=2.718^(`beta`i'')
			local lci`i'=2.718^(`beta`i''-1.96*`se`i'')
			local uci`i'=2.718^(`beta`i''+1.96*`se`i'')
			local paf`i'=`pc`i''*(1-(1/`rr`i''))*100
			local lpaf`i'=`pc`i''*(1-(1/`lci`i''))*100
			local upaf`i'=`pc`i''*(1-(1/`uci`i''))*100
			
			bysort _mj: egen denom = count(aln) if _mj != 0
			bysort _mj: egen num2 = count(aln) if _mj != 0 & group==`i'
			bysort _mj: egen num = min(num2)
			gen prop = num/denom if _mj != 0
			bysort _mj: egen seq = seq() if _mj != 0 
			summ prop if seq == 1
			local prev`i' = r(mean)*100
			drop denom num num2 prop seq
			}
					
		*Risk differences, based on: 
		xi: mim: glm `out' i.group marstat birthwgt parity i.highested matage smokpreg hhsocclas_bin preg_mrep, family(binomial) link(identity)
		mim, storebv

		forvalues i=1/2{
			local rd`i'=_b[_Igroup_`i']
			local rd_se`i' = _se[_Igroup_`i']
			local rd_lrd`i'=`rd`i''-1.96*`rd_se`i''
			local rd_urd`i'=`rd`i''+1.96*`rd_se`i''
			}	
		}

timer off 5

********************************************************************************	
timer off 1
timer list