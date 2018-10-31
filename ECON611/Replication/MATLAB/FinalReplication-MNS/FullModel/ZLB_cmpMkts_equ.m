function [ resid , J_L, J, J_P] = ZLB_cmpMkts_equ( X, eps )

[ resid , J_L, J, J_P] = ZLB_equ( X, [eps; 1;1] );



end

