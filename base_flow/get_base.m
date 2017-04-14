clear all, close all, clc

%compute a base flow satisfying the steady equations of motion using a
%Newton-Raphson method.

file_start = 'base.mat';
%--Restarting from earlier guess?

    %don't already have a guess
    if exist( file_start , 'file') ~= 2
    
        parms.Re = 100; %Reynolds #
        parms.m = 200; %# of x points on finest domain
        parms.len = 3; %length of domain in x-direction
        parms.n = 80; %# of y points on finest domain (length in y-dirn is len/m * n )
        parms.mg = 6; %# of domains
        parms.U_body = -1; %velocity of body in lab frame
        parms.offx = 0.2; %offset in x dirn (on fine domain, x-grid runs from
                          %-parms.offx to len-parms.offx).
        parms.offy = 0.6; %offset in y dirn (same as offx but in y-dirn)
        parms.tol = 1e-8; %tolerance for terminating N-R
        
        %Parameters for flag
        parms.R_rho = 30;
        parms.R_E = 1e3;
        parms.R_sh = 1e-7;
        parms.R_th = 0.01;
        parms.inverted = 'T';
        parms.clamped = 'T';
        
        mats = []; soln = [];
        
    %otherwise read from file
    else
        load( file_start )
        %Parameters for flag
        parms.R_rho = 50;
        parms.R_E = 1e5;
        parms.R_sh = 4e-6;
        parms.R_th = 0.01;
        parms.inverted = 'T';
        parms.clamped = 'T';
    end
        
%--

%run main function
main_fun( parms , mats, soln )
    
    
    


