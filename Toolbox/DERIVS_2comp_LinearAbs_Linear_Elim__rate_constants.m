function  dcdt = DERIVS_2comp_LinearAbs_Linear_Elim__rate_constants(t,c,deriv_params)
    %required parameters:
    %p.ka
    %p.k_elim
    %p.k12
    %p.k21
    %p.V1
    %p.V2
    
    a=c(1);
    c1=c(2);
    c2=c(3);
    
    p=deriv_params;
    dcdt = [ - 1*p.ka*a, 
             p.F*p.ka*(a/p.V1) - p.k_elim*c1 - p.k12*c1*p.V1 + p.k21*c2*p.V2/p.V1,        
                                             +p.k12*c1*p.V1/p.V2 -  p.k21*c2*p.V2 ];
    end