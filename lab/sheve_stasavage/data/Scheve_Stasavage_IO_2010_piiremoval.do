/*
	Script to deal with ages under 18 and listed union membership. 
	Produces an edited version of the dta and csv files. 

*/

// Directory
cd "D:\Users\jw2282\studies\D062\"

use Scheve_Stasavage_IO_2010_AIPO0263INCOME3W.dta, clear

// The variable to bottom code is AGE
replace AGE = 18 if AGE < 18

save Scheve_Stasavage_IO_2010_AIPO0263INCOME3W_EDITED.dta, replace
outsheet using Scheve_Stasavage_IO_2010_AIPO0263INCOME3W_EDITED.csv, comma replace
 

use Scheve_Stasavage_IO_2010__AIPO0242FW.dta, clear

replace AGE = 18 if AGE < 18

save Scheve_Stasavage_IO_2010_AIPO0242FW_EDITED.dta, replace
outsheet using Scheve_Stasavage_IO_2010_AIPO0242FW_EDITED.csv, comma replace


use Scheve_Stasavage_IO_2010__AIPO0242FW.dta, clear


replace AGE = 18 if AGE < 18

save Scheve_Stasavage_IO_2010__AIPO0242FW_EDITED.dta, replace
outsheet using Scheve_Stasavage_IO_2010__AIPO0242FW_EDITED.csv, comma replace
