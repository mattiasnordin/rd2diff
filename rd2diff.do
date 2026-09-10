* non-parametric p-values & confidence intervals for discrete RD setting 
* double-differences approach 
* performs inference on micro data (default) or data aggregated to the running variable 
* micro regression allows for clustering (except at level of running variable) & controls 

cap program drop rd2diff
program define rd2diff, rclass 

	syntax varlist( min=2 max=2 numeric) [if]  [aweight fweight pweight /], cutoff(real) ///
	[ siglevel(real 0.05) cluster(varlist min=1 max=1) controls(varlist) sem(varlist min=1 max=1) saveto(string) ] 
	
	preserve
	
	local yvar : word 1 of `varlist'
	local runvar : word 2 of `varlist' 
	
	qui {			
			// select analysis sample 
		if "`if'"!="" {
			keep `if'
		}
		
		drop if `runvar'==. 
		
		levelsof `runvar', hexadecimal 
		local suppts1 = r(r)
		
		foreach var in `yvar' `controls' `exp' `cluster' {			
			drop if `var'==. 
		}		
		if "`sem'"=="" {
			bys `runvar': drop if _N<2 
		}
		
		levelsof `runvar', hexadecimal
		local suppts2 = r(r)
					
			// check for contradictions 
		if `suppts2'<`suppts1' {
			noisily di as err "Insufficient number of observations at one or more support points, possibly due to dropping missing values; cannot compute variance."
			exit 498
		}		
		if "`sem'"!="" {
			levelsof `runvar' 
			cap assert r(r)==r(N)
			if _rc!=0 {
				noisily di as err "Data apparently not aggregated at level of running variable."
				exit 498
			}
			else {
				if "`weight'"!="" | "`controls'"!="" | "`cluster'"!="" {
					noisily di as err "Weight, list of controls, or clustering variable has been supplied, but data are aggregated at level of running variable."
					exit 498
				}
			}
		}
						
		sum `runvar'
		local Nttl1 = r(N)  
		if `cutoff'>r(max) | `cutoff'<=r(min) {
			noisily di as err "Cut-off lies outside permissible range."
			exit 498
		}	
		if "`yvar'"=="`runvar'" {
			noisily di as err "Running variable same as outcome variable."
			exit 498
		}				
		
			// index running variable, construct selector matrix 
		levelsof `runvar', hexadecimal local(lvls)
		
		tempvar g 
		gen `g' = .
		tempname X Mu V R Xsel stats Tmin Tmax Tstat pval BTElo BTEup CIlo CIup
			
		local i = 1
		foreach l of local lvls {
			
			replace `g' = `i' if `runvar'==`l' // running variable as index 
			if "`sem'"=="" {
				tempvar rv_`i' // dummies to be used in regression 
				gen `rv_`i'' = ( `runvar'==`l' )
			}
			if `i'==1 { // running variable to matrix 
				matrix `X' = `l'
			}
			else if `i'>1 {
				matrix `X' = `X' \ `l' 
			}		
			local ++i 		
		}	
		cap assert `g'!=. 
		if _rc>0 {
			noisily di as err "Failed to index running variable."
			exit 498
		}	
		local G `=rowsof(`X')'
		
		sum `runvar' if `runvar'>=`cutoff' 
		local Ntreat = r(N)
		sum `g' if `runvar'==r(min) 
		local c = r(mean)
		
		mata: Rmatrix("`X'",`c')
		
		matrix `R' = r(R)
		matrix `Xsel' = r(Xsel)
				
			// get covariance matrix 
		if "`sem'"=="" {
			
			if "`exp'"!="" {
				local wt "[`weight' = `exp']"
			}			
			if "`cluster'"=="" {
				local se "robust"
			}
			else {
				local se "cluster(`cluster')"
			}				
			
			reg `yvar' `rv_1'-`rv_`G'' `controls' `wt', nocons `se'
			local Nttl2 = e(N)
			
			if `Nttl2'<`Nttl1' {
				noisily di as err "Regression drops observations unexpectedly."
				exit 498
			}
			
			matrix `Mu' = ( e(b)[1,1..`G'] )'
			matrix `V' = e(V)[1..`G',1..`G']		
		}
		else {		
			levelsof `runvar' 
			cap assert r(r)==r(N)
			if _rc!=0 {
				noisily di "Data apparently not aggregated at level of running variable."
				exit 198
			}
			
			tempvar Mu V
			
			gen `Mu' = `yvar'
			gen `V' = `sem'^2
			
			mkmat `Mu', mat(`Mu')
			mkmat `V', mat(`V')	
			matrix `V' = diag(`V')
		}
		
			// get estimates
		mata: allstats("`Mu'","`V'","`R'",`siglevel')
		
		matrix `stats' = `Xsel', r(stats) 		
		matrix colnames `stats' = runvr DD S T CIlo CIup
		if "`saveto'"!="" {			
			clear 			
			svmat `stats', n(col)	
			save "`saveto'", replace 
		}	
		return matrix stats = `stats' 
		
		foreach sc in Tmin Tmax Tstat pval BTElo BTEup CIlo CIup {
			return scalar `sc' = r(`sc')
		}
	}	
	
	local Gtreat = `G' - `c' + 1
	local Lvl : di ceil(100*(1-`siglevel'))
	
	di ""
	
	di as text "Dependent variable: " _col(22) abbrev("`yvar'",15) _col(39) "Cut-off" _col(59) " = " _col(62) as result %9.2g `cutoff'
	di as text "Running variable: " _col(22) abbrev("`runvar'",15) _col(39) "Support points" _col(59) " = " _col(62) as result %9.2gc `G'
	di as text  _col(39) "Suppt pts, treated" _col(59) " = " _col(62) as result %9.2gc `Gtreat'
	if "`sem'"=="" {
		di as text  _col(39) "Observations" _col(59) " = " _col(62) as result %9.2gc `Nttl1'
		di as text  _col(39) "Obs, treated" _col(59) " = " _col(62) as result %9.2gc `Ntreat'
	}
	di as text "{hline 27}{c TT}{hline 42}"
	di as text _col(12) "Minimum T-ratio {c |} " _col(25) as result %9.5fc return(Tmin)	
	di as text _col(12) "Maximum T-ratio {c |} " _col(25) as result %9.5fc return(Tmax)	
	di as text _col(13) "Test statistic {c |} " _col(25) as result %9.5fc return(Tstat)	
	di as text _col(20) "p-value {c |} " _col(25) as result %9.5fc return(pval)	
	di as text _col(0) "Bounds on treatment effect {c |} " 	_col(50) "[" as result %9.5fc return(BTElo) "," as result %9.5fc return(BTEup) "]"
	di as text _col(4) "`Lvl'% confidence interval {c |} " 	_col(50) "[" as result %9.5fc return(CIlo) "," as result %9.5fc return(CIup) "]"
	di as text "{hline 27}{c BT}{hline 42}"		
end 	
	
		// Mata functions 
mata:
	mata clear 
	mata set matastrict on 
	
	void Rmatrix( string scalar M, real scalar c ) // exractor matrices 
	{
		
		real colvector X
		real scalar G
		real scalar cgap 
		real colvector phi 
		real matrix R
		real colvector Xsel 
		real scalar i
		
		X = st_matrix(M)
		
		G = rows(X)
		
		cgap = X[c] - X[c-1]
		
		phi = J(G,1,0)
		
		for (i=2; i<=G; i++) {
			phi[i-1] = cgap/( X[i] - X[i-1] )
		}
		
		R = diag( phi ) -  ( J(G-1,1,0),diag(phi[1..G-1]) \ J(1,G,0) ) + ( J(G,c-2,0) , -J(G,1,1), J(G,1,1), J(G,G-c,0) )
		
		if (c==2) {
			R = R[c..G-1,.] 
			Xsel = X[3..G] 
		}
		else if (c==G) {
			R = R[1..c-2,.]
			Xsel = X[2..c-1] 
		}
		else {
			R = R[1..c-2,.] \ R[c..G-1,.]
			Xsel = X[2..c-1] \ X[c+1..G] 
		}
		
		st_matrix("r(R)",R)
		st_matrix("r(Xsel)",Xsel)
	}	

	void allstats(string scalar Mu, string scalar V, string scalar R, real scalar alph) // estimation 
	{
		
		real matrix DD
		real colvector S
		real colvector T
		real scalar z
		real colvector CIlo
		real colvector CIup
		real scalar Tstat		
		
		DD = st_matrix(R)*st_matrix(Mu) 
		S = sqrt(diagonal(st_matrix(R)*st_matrix(V)*st_matrix(R)'))
		T = DD:/S 
		
		z = invnormal(1-alph/2)
		
		CIlo = DD - z*S 
		CIup = DD + z*S 
				
		Tstat = ( min(T)>0 )*min(T) + ( max(T)<0 )*abs(max(T))		
				
		st_matrix("r(stats)",( DD, S, T, CIlo, CIup ))
		
		st_numscalar("r(Tmin)",min(T))
		st_numscalar("r(Tmax)",max(T))
		st_numscalar("r(Tstat)",Tstat )		
		st_numscalar("r(pval)",2*( 1  - normal(Tstat)))
		st_numscalar("r(BTElo)",min(DD))
		st_numscalar("r(BTEup)",max(DD))
		st_numscalar("r(CIlo)",min(CIlo))
		st_numscalar("r(CIup)",max(CIup))	
	}	
end