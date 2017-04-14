clear all, close all, clc

%script that calls the IB code

%--check to make sure we have a base flow

    file_base = '../base_flow/base.mat';

    %don't have a base flow
    if exist( file_base , 'file') ~= 2
        error( ['Base flow required! Run get_base.m in the '...
            ' base_flow directory to compute it'])
    else
        load(file_base)
    end

%--


it_start = 1700; %Time step to start at (0 = new simulation)
it_stop = 40000; %save every it_save time steps
it_save = 100; %Time step to end at
dt = 0.001; %time step size in physical time



if it_start == 0
    soln_lin = [];
    
else
    
    load(['outputs/runvars_it_',num2str(it_start),'.mat'])
    
    parms.it_start = it_start;
    parms.it_stop = it_stop;
    parms.it_save = it_save;
    parms.dt = dt;
end

%run main function
main_fun( soln_lin, soln, mats, parms )
