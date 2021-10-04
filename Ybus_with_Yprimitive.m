%------------------------------------------------------------------------%
% START %
%------------------------------------------------------------------------%
% This is program for Creating Ybus Matrix when Yprimitive is known Method %
%------------------------------------------------------------------------%

clc;
clear all;

%------------------------------------------------------------------------%
% User Inputs %
%------------------------------------------------------------------------%

n=input('Enter number of Buses in the Power System : ');
l=input('Enter the number of Links in the Power System Graph : ');
A=input('Enter the Incidence Matrix of the Power System : ');
%Enter nodes starting with Slack Bus followed by all PQ Buses and at last
%all PV Buses

s1=size(A); % Checking Size of User entered Incidence Matrix
if (l~=s1(1,1) || n~=s1(1,2))
    fprintf('Size of Incidence Matrix is incorrect \n');
    goto(16)
    return;
end

Yprim=input('Enter the Primitive Ybus Matrix of the Power System : ');

s2=size(Yprim); % Checking Size of User entered Incidence Matrix
if (l~=s1(1,1) || l~=s1(1,2))
    fprintf('Size of Primitive Ybus Matrix is incorrect \n');
    goto(27)
    return;
end

%------------------------------------------------------------------------%
% Calculation of Ybus Matrix %
%------------------------------------------------------------------------%

Ybus=A'*Yprim*A; % Calculating Ybus Matrix

%------------------------------------------------------------------------%
% END %
%------------------------------------------------------------------------%


