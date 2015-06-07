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
% Identifiability.m
% This program uses SSM to to determine identifiability of parameters.
% ------------------------------------------------------------------------------------- % 

function [List,pset] = Identifiability (SSM)
%this function uses SSM to to determine identifiability of parameters

P_i_name=char('beta_2','beta_1','beta_m','k_2','k_1','k_E','d_2','d_1','d_m','alpha','k_I','rm2','rm1','p_c','P_t'); % The list of all parameters in order

[r,c]=size(SSM); % r=time steps, c=15
X=[];
pset=[];


% McAuley procedure doi:10.1081/PRE-120024426
R=SSM; % the first time, use SSM to do column sumComputeIdentifiability

for j=1:c
     for i=1:c
         M(i)=R(:,i)'*R(:,i); % the square sum of each column
     end
     
     [a,pos]=max(M);          % finds the indices of the maximum values of M, and returns them in output vector pos. 
     
     if a>1.0e-08                 % this is the tolerance
         X=[X SSM(:,pos)];      % the colomn that has the largest SS magnitude
         pset=[pset; pos];      % give the index of the parameter
         Shat=X*inv(X'*X)*X'*SSM; % Find the prediction SSM
         R=SSM-Shat;            % residual mtrx, now the residual matrix is the new mtrx that we find the next identifiable parameter, return to the top of the J loop  
     end
 end

pset=unique(pset);                % identifiable parameters

c=size(pset);
for i=1:c
    k=pset(i);
    List(i,:)=P_i_name(k,:);
end

end

