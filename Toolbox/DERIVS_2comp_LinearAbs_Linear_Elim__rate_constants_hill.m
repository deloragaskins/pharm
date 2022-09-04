function dydt = DERIVS_2comp_LinearAbs_Linear_Elim__rate_constants_hill(t,y,deriv_params)
p=deriv_params;
%%%%PD params
%p.n_factor
%p.kn_killmax
%p.kn_kill50
%kn_a
%p.kn_e

%PK params
%p.ka
%p.k_elim
%p.k12
%p.k21
%p.V1
%p.V2

a=y(1);
c1=y(2);
c2=y(3);
n=y(4);

kn_kill=Hill_Effect(p.n_factor*y(2),0,p.kn_killmax,p.kn_kill50,1);

dydt = [ - 1*p.ka*a, 
          p.F*p.ka*(a/p.V1) - p.k_elim*c1 - p.k12*c1*p.V1 + p.k21*c2*p.V2/p.V1,        
                                             +p.k12*c1*p.V1/p.V2 -  p.k21*c2*p.V2,
                                             p.kn_a*log(p.n_carry/n)-(p.kn_e+kn_kill)*n];    
end