C Return p given quant and df

      double precision function pchis(quant,df)
      double precision quant,df
      real pgamm
      real rquant, rdf

      rquant = quant / 2.0
      rdf = df / 2.0
      pchis = pgamm (rquant, rdf)
      return
      end


      real function pgamm(quant,eta)
      real quant,eta
      real e,term,sum,gamln,eps,r,test
      integer last,j
      precis=10e-7
      if(.not.(quant .le. 0.))goto 23000
      call chi_error (1)
      call exit (1)
23000 continue
      if(.not.(eta .le. 0.))goto 23002
      call chi_error (2)
      call exit (1)
23002 continue
      e = eta
      term = 1.
      sum = 1.
      if(.not.(quant .ge. eta .and. quant .ge. 7.0))goto 23004
      last = ifix(eta)+10
      if(.not.(quant .lt. 11.))goto 23006
      last = ifix(quant+eta)-1
23006 continue
      do 23008 j = 1,last 
      e = e-1.
      term = term*e/quant
      sum = sum+term
23008 continue
      term = term*(e-1.)
      sum = sum+term/(quant-e+2.)
      pgamm = 1.-exp(-quant+(eta-1.)*alog(quant)+alog(sum)-gamln(eta))
      goto 23005
23004 continue
      eps = precis
23010 continue
      e = e+1.
      r = quant/e
      term = term*r
      sum = sum+term
      if(.not.(r .le. 0.5))goto 23013
      test = term/sum
      if(.not.(test .le. eps))goto 23015
      goto 23012
23015 continue
23013 continue
23011 goto 23010
23012 continue
      pgamm = exp(eta*alog(quant)-quant+alog(sum)-gamln(eta))/eta
23005 continue
      return
      end
      real function gam(eta)
      real eta
      real a(8),feta,g1feta,test,prod,term
      data a/-0.577191652,0.988205891,-.897056937,0.918206857,
     & -0.756704078,0.482199394,-0.193527818,0.035868343/
      if(.not.(eta .gt. 33.))goto 23017
      call chi_error (3)
      call exit (1)
23017 continue
      if(.not.(eta .le. 0.0))goto 23019
      call chi_error (4)
      call exit (1)
23019 continue
      feta = amod(eta,1.)
      if(.not.(feta .gt. 0.))goto 23021
      g1feta = 1.+feta*(a(1)+feta*(a(2)+feta*(a(3)+feta*(a(4)+feta*(a(5)
     & +feta*(a(6)+feta*(a(7)+feta*a(8))))))))
      goto 23022
23021 continue
      g1feta = 1.
23022 continue
      if(.not.(eta .lt. 1.))goto 23023
      gam = g1feta/feta
      goto 23024
23023 continue
      if(.not.(eta .le. 2.))goto 23025
      gam = g1feta
      goto 23026
23025 continue
      test = eta-.5
      prod = 1.
      term = 1.+feta
23027 continue
      prod = prod*term
      term = term+1.
23028 if(.not.(term .ge. test))goto 23027
      gam = prod*g1feta
23026 continue
23024 continue
      return
      end
      real function gamln(eta)
      real eta
      real gamma,xlsq2p,etam,gam
      data xlsq2p/.91893854/
      if(.not.(eta .le. 0.0))goto 23030
      call chi_error (5)
      call exit (1)
23030 continue
      if(.not.(eta .le. 30.0))goto 23032
      gamma=alog(gam(eta))
      goto 23033
23032 continue
      etam = eta-1.0
      gamma = xlsq2p+(etam+.5)*alog(etam)-etam+1./(12.*etam)+1./(288.*
     & etam*etam)
23033 continue
      gamln = (gamma)
      return
      end
      subroutine chi_error(err)
      integer err
      print *, "chi_error: err =", err
      return
      end
