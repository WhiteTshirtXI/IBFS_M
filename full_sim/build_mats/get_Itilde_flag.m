function Itilde_flag = get_Itilde_flag( nb )

Itilde_flag = sparse( 2*nb, 3*nb );

indr = 1 : nb;
indcx = 1 : 3 : 3*nb-2;

Itilde_flag = Itilde_flag + sparse( indr, indcx, ...
    ones(size(indr)), 2*nb, 3*nb);

indcy = 2 : 3 : 3*nb -1;
Itilde_flag = Itilde_flag + sparse( nb + indr, indcy, ...
    ones(size(indr)), 2*nb, 3*nb);
