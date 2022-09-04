function E=HillEffect(c,E0,Emax,EC50,Hcoeff)
    E=E0+Emax.*c.^Hcoeff/(EC50.^Hcoeff+c.^Hcoeff);
end