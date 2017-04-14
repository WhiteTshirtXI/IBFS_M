function Q = get_Q_flag( nb, h_el_0, parms )

%Take a vector of nodal x and y forces and turn them into the global force
%vector for the FEM formulation...

nel = nb - 1;

Q = zeros( 3*nb, 2*nb );

for jj = 1 : 2*nb
   
    F = zeros( 3*nb, 1 );
    fv = zeros( 2*nb, 1 );
    fv(jj) = 1;
    
    fx = fv( 1 : nb); fy = fv( nb + (1 : nb ) );

    for j = 1 : nel

        %Build element force vector:
        fxsmall = fx(j : j + 1); %Get x forces at that element
        fysmall = fy(j : j + 1); %Get y forces at that element
        F_e = get_Fe( fxsmall, fysmall, h_el_0 );

        %Assemble force vector:
        index = get_ind_flag( j );
        F = F_assemble( F, F_e, index );


    end
    
    %BCs

    %clamped?
    if parms.clamped == 'T'

        bc_type = [1 1 1];

    else

        bc_type = [1 1 0];

    end

    count = 0;

    if parms.inverted == 'T'

        for j = 3*nb - 2 : 3*nb

            count = count + 1;

            if bc_type(count) == 1
                F(j) = 0;
            end

        end

    else

        for j = 1 : 3

            count = count + 1;

            if bc_type(count) == 1
                F(j) = 0;
            end

        end


    end
    
    Q(:, jj) = F;
    

end

