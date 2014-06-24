c
c  icslap.h
c
c  Definitions used by icslap.f
c

	integer lenw, leniw, nelt, nelt_lt, isym
	parameter(nelt = 100000,ndvar=2, nelt_lt=(nelt+ndvar)/2+1,
     &            leniw=nelt_lt + ndvar+11,
     &            lenw=nelt_lt+5*ndvar,
     &            isym=0)
        double precision rwork(lenw)
        integer iwork(leniw) 	
