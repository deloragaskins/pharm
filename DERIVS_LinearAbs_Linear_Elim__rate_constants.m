FUNCTION dcdt = DERIVS_LinearAbs_Linear_Elim__rate_constants(a,c1,c2,params)
    %required parameters:
    %p.ka
    %p.k_elim
    %p.k12
    %p.k21
    %p.V1
    %p.V2
    p=params;
    dcdt = [ - 1*p.ka*a, 
             p.F*p.ka*(a/p.V1) - p.k_elim*c1 - p.k12*c1*p.V1 + p.k21*c2*p.V2/p.V1,        
                                             +p.k12*c1*p.V1/p.V2 -  p.k21*c2*p.V2 ];

end