clear all, close all, clc

%script that calls the IB code

%If it_start = 0, change variables in "if" statement to desired quantities.
%Otherwise set the timestep at which you want to restart from file.
it_start = 0; %Time step to start at (0 = new simulation)
it_stop = 1000; %save every it_save time steps
it_save = 1000; %Time step to end at

%set body
parms.body = 'cyl';
parms.L = 1; %Main length scale in problem (e.g. for cylinder, diameter).

%Is body moving?
parms.stationary = 'T';

if ( it_start == 0 )
    parms.it_start = it_start; 
    parms.it_stop = it_stop; 
    parms.it_save = it_save; 
    parms.Re = 50; %Reynolds #
    parms.m = 100; %# of x points on finest domain
    parms.len = 3; %length of domain in x-direction
    parms.n = 100; %# of y points on finest domain (length in y-dirn is len/m * n )
    parms.dt = 0.001; %time step size in physical time
    parms.mg = 5; %# of domains
    parms.U_body = -1; %velocity of body in lab frame
    parms.offx = 0.8; %offset in x dirn (on fine domain, x-grid runs from
                      %-parms.offx to len-parms.offx).
    parms.offy = 1.5; %offset in y dirn (same as offx but in y-dirn)

else
    file_start = ['outputs/runvars_it_',num2str(it_start),'.mat'];
    if exist( file_start , 'file') == 2
        load( ['outputs/runvars_it_',num2str(it_start),'.mat'] )
        parms.it_start = it_start;
        parms.it_stop = it_stop;
        parms.it_save = it_save;

    else
        error( 'File does not exist at desired start point!')
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
