#!/bin/csh
#---------------------------------------------------------------------
# Setting |
#---------------------------------------------------------------------

set sdate=(1979)
set edate=(2023)

#---------------------------------------------------------------------

@ iyr = ${sdate[1]}

while(${iyr} <= ${edate[1]})

    @ jyr = $iyr + 1

    csh make_atm.csh ${iyr} 1 1 ${jyr} 1 1
    csh make_fflux.csh ${iyr} 1 1 ${jyr} 1 1
    
    @ iyr++
    
end
