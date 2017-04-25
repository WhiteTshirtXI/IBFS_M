clear all, close all, clc

%script that runs the IB code

%If it_start = 0, change variables in "if" statement to desired quantities.
%Otherwise set the timestep at which you want to restart from file.
it_start = 0; %Time step to start at (0 = new simulation)
it_stop = 30000; %save every it_save time steps
it_save = 1000; %Time step to end at

if ( it_start == 0 )
    parms.it_start = it_start; 
    parms.it_stop = it_stop; 
    parms.it_save = it_save; 
    parms.Re = 200; %Reynolds #
    parms.m = 200; %# of x points on finest domain
    parms.len = 2; %length of domain in x-direction
    parms.n = 64; %# of y points on finest domain (length in y-dirn is len/m * n )
    parms.dt = 0.001; %time step size in physical time
    parms.mg = 5; %# of domains
    parms.U_body = -1; %velocity of body in lab frame
    parms.offx = 0.2; %offset in x dirn (on fine domain, x-grid runs from
                      %-parms.offx to len-parms.offx).
    parms.offy = 0.32; %offset in y dirn (same as offx but in y-dirn)
    
    parms.body = 'flg'; %body ('cyl' for cylinder with r=.5, 'flg' for 
                        %flag with length 1)
        %If body is not one of those two flags, specify body in file
        %body.mat that will be read by code.
    parms.stationary = 'F'; %Is body moving?
    parms.deform = 'T'; %Is body deformable?

    
    %Parameters for flag (unused if parms.deform = 'F')
    parms.R_rho = 20;
    parms.R_E = 1e3;
    parms.R_sh = 1e-7;
    parms.R_th = 0.01;
    parms.inverted = 'F';
    parms.clamped = 'F';

else
    file_start = ['outputs/runvars_it_',num2str(it_start),'.mat'];
    if exist( file_start , 'file') == 2
        load( ['outputs/runvars_it_',num2str(it_start),'.mat'] )
        parms.it_start = it_start;
        parms.it_stop = it_stop;
        parms.it_save = it_save;

    else
        error( 'File does not exist at desired start time!')
    end

end
%---

if it_start == 0
    soln = [];
end

mats = [];

%run main function
tic
main_fun( soln, mats, parms )

simtime = toc;
display(['Spent ',num2str(simtime),' secs on run.'])
