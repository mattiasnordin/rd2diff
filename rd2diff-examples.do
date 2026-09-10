* how to use rd2diff 

cd "" // set working directory, which must contain the file rd2diff.do 
		
do rd2diff // load the command 

clear 

	// create basic features of dataset 
set obs 20

gen x = rnormal() 
gen g = _n 
gen n = 50 + rbinomial(400,0.5)

expand n

	// an example with ordinal running variable 
local cut = 12

gen d1 = g>=`cut' 
gen y1 = g + 0.5*d1 + rnormal() 

rd2diff y1 g, cutoff(`cut') 
rd2diff y1 g if inrange(g,`cut'-3,`cut'+2), cutoff(`cut')  // restrict the bandwidth 

	// an example with discrete running variable whose support points are drawn from the reals 
local cut = 0

gen d2 = x>=`cut' 
gen y2 = x + 0.5*d2 + 0.1*rnormal() 

rd2diff y2 x, cutoff(`cut') 

	// an example with a control variable 
local cut = 12
	
gen d3 = g>=`cut'	
	
gen w = d3 + rnormal()
gen y3 = g + w + 0.5*rnormal()

rd2diff y3 g, cutoff(`cut') 
rd2diff y3 g, cutoff(`cut') controls(w)

	// estimation on aggregate data
local cut = 12

rd2diff y1 g, cutoff(`cut') 
	
collapse y1 (sem) y1se = y1, by(g)

rd2diff y1 g, cutoff(`cut') sem(y1se)