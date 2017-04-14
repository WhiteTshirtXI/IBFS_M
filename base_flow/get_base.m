clear all, close all, clc

%compute a base flow satisfying the steady equations of motion using a
%Newton-Raphson method.

file_start = 'base.mat';
%--Restarting from earlier guess?

    %don't already have a guess
    if exist( file_start , 'file') ~= 2

        %set body
        parms.body = 'cyl';

        %Is body moving?
        parms.stationary = 'T';
        parms.Re = 47; %Reynolds #
        parms.m = 160; %# of x points on finest domain
        parms.len = 5; %length of domain in x-direction
        parms.n = 120; %# of y points on finest domain (length in y-dirn is len/m * n )
        parms.mg = 6; %# of domains
        parms.U_body = -1; %velocity of body in lab frame
        parms.offx = 0.8; %offset in x dirn (on fine domain, x-grid runs from
                          %-parms.offx to len-parms.offx).
        parms.offy = 1.5; %offset in y dirn (same as offx but in y-dirn)
        parms.tol = 1e-8; %tolerance for terminating N-R
        
        mats = []; soln = [];
        
    %otherwise read from file
    else
        load( file_start )
        parms.Re = 55; %Reynolds #
    end
        
%--

%run main function
main_fun( parms , mats, soln )
    
    
    


