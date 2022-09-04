%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = DERIVS_1comp_LinearAbs_Linear_Elim__rate_constants(t, y, params)
dydt = [-params.ka*y(1), 
         params.ka*y(1) - params.k_elim*y(2) ];
end
