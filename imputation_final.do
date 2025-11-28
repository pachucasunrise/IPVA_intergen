*09.10.22

local parental_dv "dv_mrep phys_mrep psych_mrep dv_mrep_acts phys_mrep_acts psych_mrep_acts"
local child_dv "vic_1821 per_1821"
local ace_cats "maltreat restoftrio aces3"
local ace_extra "bullying_0_1 parental_sep emotional_ab emotional_ne physical_abu sexual_abuse mentl_hlth_p substance_ho parent_convi"

foreach v of local parental_dv{	
	foreach p of local child_dv {
		use full_alspac_update.dta, clear
		ice `v' `p' `ace_extra' `ace_cats' preg_mrep marstat matage parity highested hhsocclas mum_depress aggress birthwgt smokpreg mum_sexmin, ///
		by(kz021) m(35) ///
		cmd(parity highested hhsocclas: ologit, `v' `p' `ace_extra' preg_mrep marstat aggress mum_sexmin: logit, matage mum_depress: regress) ///
		passive(maltreat: (emotional_ab==1 | emotional_ne==1 | physical_abu==1 | sexual_abuse==1) \ restoftrio: (mentl_hlth_p==1 & substance_ho==1) \ aces3: ((emotional_ne + emotional_ab + physical_abu + bullying_0_1 + substance_ho + mentl_hlth_p + parent_convi + parental_sep + sexual_abuse)>=3 & (emotional_ne + emotional_ab + physical_abu + bullying_0_1 + substance_ho + mentl_hlth_p + parent_convi + parental_sep + sexual_abuse) != .)) ///
		eqdrop(`v' `p' `ace_extra' restoftrio aces3 preg_mrep marstat matage birthwgt parity smokpreg highested hhsocclas aggress: maltreat, `v' `p' `ace_extra' maltreat aces3 preg_mrep marstat matage birthwgt parity smokpreg highested hhsocclas aggress: restoftrio, `v' `p' `ace_extra' maltreat restoftrio preg_mrep marstat matage birthwgt parity smokpreg highested hhsocclas aggress: aces3, `v' `p' `ace_extra' `ace_cats' preg_mrep marstat matage parity highested hhsocclas aggress mum_depress smokpreg: birthwgt, `v' `p' `ace_extra' `ace_cats' preg_mrep marstat matage parity highested hhsocclas aggress mum_depress birthwgt: smokpreg) ///
		seed(1234) saving(imputed_alspac_bysex_`v'_`p'.dta, replace) 
		
		*No need to impute vic or perp for study sample as everything complete
		use preimputation_update, clear
		ice `v' `p' `ace_extra' `ace_cats' preg_mrep marstat matage parity highested hhsocclas mum_depress aggress birthwgt smokpreg mum_sexmin, ///
		by(kz021) m(35) ///
		cmd(parity highested hhsocclas: ologit, `v' `p' `ace_extra' preg_mrep marstat aggress mum_sexmin: logit, matage mum_depress: regress) ///
		passive(maltreat: (emotional_ab==1 | emotional_ne==1 | physical_abu==1 | sexual_abuse==1) \ restoftrio: (mentl_hlth_p==1 & substance_ho==1) \ aces3: ((emotional_ne + emotional_ab + physical_abu + bullying_0_1 + substance_ho + mentl_hlth_p + parent_convi + parental_sep + sexual_abuse)>=3 & (emotional_ne + emotional_ab + physical_abu + bullying_0_1 + substance_ho + mentl_hlth_p + parent_convi + parental_sep + sexual_abuse) != .)) ///
		eqdrop(`v' `p' `ace_extra' restoftrio aces3 preg_mrep marstat matage birthwgt parity smokpreg highested hhsocclas aggress: maltreat, `v' `p' `ace_extra' maltreat aces3 preg_mrep marstat matage birthwgt parity smokpreg highested hhsocclas aggress: restoftrio, `v' `p' `ace_extra' maltreat restoftrio preg_mrep marstat matage birthwgt parity smokpreg highested hhsocclas aggress: aces3, `v' `p' `ace_extra' `ace_cats' preg_mrep marstat matage parity highested hhsocclas aggress mum_depress smokpreg: birthwgt, `v' `p' `ace_extra' `ace_cats' preg_mrep marstat matage parity highested hhsocclas aggress mum_depress birthwgt: smokpreg) ///
		seed(1234) saving(imputed_sample_bysex_`v'_`p'.dta, replace) 
		}
	}
