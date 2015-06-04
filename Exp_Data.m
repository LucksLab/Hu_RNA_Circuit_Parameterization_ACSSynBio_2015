% ------------------------------------------------------------------------------------- %
% Copyright (c) 2015 Lucks Lab,
% School of Chemical and Biomolecular Engineering,
% Cornell University, Ithaca NY 14853 USA.
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions: % The above copyright notice and this permission notice shall be included in
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
% Exp_Data.m
% A function loads experimental data from text file and passes to main
% script
% ------------------------------------------------------------------------------------- %


function [ Exp_Data_DF] = Exp_Data( CSTR_LV )

if CSTR_LV==1
    Data= dlmread('Exp_Data_1LV.txt', '\t'); 
    [r,c]=size(Data);
    time=Data(:,1);
    GFP=Data(:,2:c);  

    stepsz=size(time,1);

    for i=1:(stepsz)

        Upbound(i)=max(GFP(i,:));
        lowbound(i)=min(GFP(i,:));
        avg(i)=sum(GFP(i,:))/(c-1); 

    end


elseif CSTR_LV==3
   Data= dlmread('Exp_Data_2LV+3LV.txt', '\t'); 
    [r,c]=size(Data);
    time=Data(:,1);
    GFP=Data(:,2:c);  

    stepsz=size(time,1);

    for i=1:(stepsz)

        Upbound(i)=max(GFP(i,:));
        lowbound(i)=min(GFP(i,:));
        avg(i)=sum(GFP(i,:))/(c-1); 

    end


elseif CSTR_LV==2
   Data= dlmread('Exp_Data_2LV+1LV.txt', '\t'); 
    [r,c]=size(Data);
    time=Data(:,1);
    GFP=Data(:,2:c);  

    stepsz=size(time,1);

    for i=1:(stepsz)

        Upbound(i)=max(GFP(i,:));
        lowbound(i)=min(GFP(i,:));
        avg(i)=sum(GFP(i,:))/(c-1); 

    end


elseif CSTR_LV==4
    Data= dlmread('Exp_Data_3LV+2LV+1LV.txt', '\t'); 
    [r,c]=size(Data);
    time=Data(:,1);
    GFP=Data(:,2:c);  
    stepsz=size(time,1);

    for i=1:(stepsz)

        Upbound(i)=max(GFP(i,:));
        lowbound(i)=min(GFP(i,:));
        avg(i)=sum(GFP(i,:))/(c-1); 

    end
end

Exp_Data_DF.Data=GFP;
Exp_Data_DF.avg=avg';
Exp_Data_DF.Upbound=Upbound';
Exp_Data_DF.Lowbound=lowbound';
Exp_Data_DF.timestep=stepsz;
end

