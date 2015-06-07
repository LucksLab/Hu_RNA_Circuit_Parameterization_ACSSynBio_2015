% ------------------------------------------------------------------------------------- %
% Copyright (c) 2015 Lucks Lab,
% School of Chemical and Biomolecular Engineering,
% Cornell University, Ithaca NY 14853 USA.
%
% Details about this software can be found in 
% Hu et al. “Generating effective models and parameters for RNA genetic circuits” 
% ACS SynBio, 2015, http://dx.doi.org/10.1021/acssynbio.5b00077 .
%
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions: 
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.
%
% Double_Repression_Analytical_SSM.m
% A function that solves the sensitivity matrix analytically with given ODE and parameters 
% ------------------------------------------------------------------------------------- %

function [ SSM ] = Double_Repression_Analytical_SSM( DF,x,nstep,t_inc,Prev_ID)

P=DF.Initial_Parameters;   %Obtain the initial guess vector from struct
P_size=DF.Num_Parameters;  
C=DF.Constants;           
CSTR_LV=DF.Construct;

[ J_mtrx ] = Double_Repression_Jacobian( P,C,x,nstep ); % get the jacobian matrix
[ P_mtrx ] = Double_Repression_Pmatrix (P,C,x,nstep );  % get the pmatrix
z0=[0;0;0;0;0;0;0;0];  % Initial value for z thats used in 0de15s 

SSM=zeros((nstep),(P_size)); 


for i=2:(nstep+1)
    print=['solving for SSM at point ',num2str(i),'/',num2str(nstep+1)];
    disp(print);
    
    for j=1:(P_size)

        
        tspan=0:(t_inc):((i-1)*t_inc);  % to save computational time, only calculate z to the time point of interest
        [t,z]=ode15s(@senfunction,tspan,z0,[],J_mtrx,P_mtrx,t_inc,j); 
 
        SSM(i,j)=z(i,8);     % store z into SSM, here only take z(8), which is Gm
   
    end

end

 [SSM]=Normalize_SSM (DF,SSM,P,x,nstep);  % Identifiablity was estimated using an normalized SSM

if CSTR_LV==1
    
     SSM(:,1)=0; %P1 and P2 are in fixed ratio with P3, also P1=P2=0 when there is only one level present
     SSM(:,2)=0;
     
     for i=1:size(Prev_ID)
        SSM(:,Prev_ID(i))=0; %skips previous identified parameter and mark according columns 0
     end
    
elseif CSTR_LV==2
    
     SSM(:,1)=0;
     SSM(:,2)=0;
     
     for i=1:size(Prev_ID)
        SSM(:,Prev_ID(i))=0;
     end
     
elseif CSTR_LV==3

     SSM(:,1)=0;
     SSM(:,2)=0;
     
     Prev_ID(Prev_ID==5|Prev_ID==8|Prev_ID==13)=1; % if 5 8 13 are previously identified by CSTR_LV2, they are not counted as identified in this level

     for i=1:size(Prev_ID)
        SSM(:,Prev_ID(i))=0;
     end
     
 elseif CSTR_LV==4  
     SSM(:,1)=0;
     SSM(:,2)=0;
     
     for i=1:size(Prev_ID)
        SSM(:,Prev_ID(i))=0;
     end
     

end

