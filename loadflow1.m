%------------------------------------------------------------------------%
%START%
%------------------------------------------------------------------------%
% This is program for Load Flow studies using GAUSS-SEIDEL Method %
%------------------------------------------------------------------------%

clc;
clear all;

%------------------------------------------------------------------------%
% User Inputs %
%------------------------------------------------------------------------%

n=input('Enter number of Buses in the Power System : ');
l=input('Enter the number of Links in the Power System Graph : ');
PVnum=input('Enter number of PV Buses : ');
A=input('Enter the Incidence Matrix of the Power System : ');
%Enter nodes starting with Slack Bus followed by all PV Buses and at last
%all PQ Buses

s1=size(A); % Checking Size of User entered Incidence Matrix
if (l~=s1(1,1) || n~=s1(1,2))
    fprintf('Size of Incidence Matrix is incorrect \n');
    goto(17)
    return;
end

LP=input('Enter the Matrix of Link Parameters : ');
% Column 1- Link Number, Column 2- Series Impedance 
% Column 3- Shunt Admittance (of lines should be the 1/2 the total i.e
% equivalen t to one leg of the Pi-Circuit model of transmission line
% Any Unspecified value should be initialised to zero

s2=size(LP);% Checking Size of User entered Link Parameter Matrix
if (l~=s2(1,1) || 3~=s2(1,2))
    fprintf('Size of Link Parameter Matrix is incorrect \n');
    goto(28)
    return;
end

SP=input('Enter System Information Matrix : ');
% Column 1- Bus Number(n) , Column 2- Real Power p.u.(Pi)
% Column 3- Reactive Power p.u.(Qi), Column 4- Max Qi p.u. 
% Column 5- Min Qi p.u. Column 6-Bus Voltage p.u. (Vi), Column 7- Bus Type
% Bus Types : 1- Slack Bus, 2- PQ Bus, 3- PV Bus
% Any value which is not specified should be initialised to zero, provide
% flat start for unspecified voltages
%Enter nodes starting with Slack Bus followed by all PV Buses and at last
%all PQ Buses

s3=size(SP);% Checking Size of User entered System Parameter Matrix
if (n~=s3(1,1) || 7~=s3(1,2))
    fprintf('Size of Line Parameter Matrix is incorrect \n');
    goto(41)
    return;
end   

%------------------------------------------------------------------------%

% Creating Ybus Matrix  %
%------------------------------------------------------------------------%

Ybus=zeros(n);% Initialising Ybus Matrix
for i=1:n
    for j=1:n
        if (i==j)
            for i1=1:l
                if (A(i1,i)==1 || A(i1,i)==-1)
                    Ybus(i,j)=Ybus(i,j)+(1/LP(i1,2))+LP(i1,3);
                end
            end
        else
            for j1=1:l
                if ((A(j1,i)==1 || A(j1,i)==-1)&&(A(j1,j)==1 || A(j1,j)==-1))
                    Ybus(i,j)=-1/(LP(j1,2));
                end
            end
        end
    end
end

%------------------------------------------------------------------------%

%Creating Matrices of size [n 2] of the Bus Parameters i.e intialising 
%Parameter vectors 
%------------------------------------------------------------------------%

PQnum=n-PVnum-1; % Claculating number of PQ Buses

% Column 1- For calculation of present iteration
% Column 2- For storing values of last iteration

V=zeros(n,2);
Del=zeros(n,2);
P=zeros(n,2);
Q=zeros(n,2);
for i2=1:n
    V(i2,1)=SP(i2,6);
    Del(i2,1)=angle(V(i2,1));
    P(i2,1)=real(SP(i2,2));
    Q(i2,1)=real(SP(i2,3));
    
    V(i2,2)=SP(i2,6);
    Del(i2,2)=angle(V(i2,2));
    P(i2,2)=abs(SP(i2,2));
    Q(i2,2)=abs(SP(i2,3));
end
while_loop_counter=0;% Setting iteration count
ee=1; % Initialising while condition
    
%------------------------------------------------------------------------%

% Gauss-Seidel Iterations %
%------------------------------------------------------------------------%

while (ee==1)
    while_loop_counter=while_loop_counter+1;% Incrasing iteration count
    
   
    
    % Calculating Voltage vector for PV Buses %
 %------------------------------------------------------------------------%
    for i4=2:PVnum+1
        YV2=0+0i; % Initialising summation part of the equation to zero
        for k2=1:n
            YV2=YV2+(Ybus(i4,k2)*V(k2,1));
        end
        Q(i4,1)=-imag(conj(V(i4,1))*YV2);% Qi calculation
        
        VQ=1*(cos(Del(i4,1))+sin(Del(i4,1))*i);% Flat Start for 
        
        %  Checking of violation of Qi limits%
%------------------------------------------------------------------------%        
        if (Q(i4,1)>SP(i4,4)) % Checking Qi>Qimax
            Q(i4,1)=SP(i4,4);
            
            YV3=0+0i; %initialising Summation part of the equation
           
            for k3=1:n
                YV3=YV3+(Ybus(i4,k3)*V(k3,1));
            end
            V(i4,1)=(1/Ybus(i4,i4))*(((P(i4,1)-Q(i4,1)*i)/(conj(VQ)))+...
            (Ybus(i4,i4)*V(i4,1))-YV3);% Vi calculation
            Del(i4,1)=angle(V(i4,1));% Updating Del Value
            
        elseif (Q(i4,1)<SP(i4,5)) % Checking Qi<Qimin
                Q(i4,1)=SP(i4,5);
            
                YV4=0+0i; %initialising Summation part of the equation
           
                for k4=1:n
                    YV4=YV4+(Ybus(i4,k4)*V(k4,1));
                end
                V(i4,1)=(1/Ybus(i4,i4))*(((P(i4,1)-Q(i4,1)*i)/(conj(VQ)))+...
                (Ybus(i4,i4)*V(i4,1))-YV4);% Vi calculation
                Del(i4,1)=angle(V(i4,1));% Updating Del Value
                
        else % Qi limits not violated
            
            YV5=0+0i; %initialising Summation part of the equation
           
            for k5=1:n
                YV5=YV5+(Ybus(i4,k5)*V(k5,1));
            end
            Del(i4,1)=angle((1/Ybus(i4,i4))*(((P(i4,1)-Q(i4,1)*i)/(conj(V(i4,1))))+...
            (Ybus(i4,i4)*V(i4,1))-YV5));% Calculating Del Value
            V(i4,1)=abs(SP(i4,6))*(cos(Del(i4,1)+sin(Del(i4,1))*i));% Updating Vi Value
        end
    end
    
     % Calculating Voltage vector for PQ Buses %
%------------------------------------------------------------------------%    
  
    for i3=PVnum+2:n
        YV1=0+0i; %initialising Summation part of the equation
        for k1=1:n
            YV1=YV1+(Ybus(i3,k1)*V(k1,1));
        end
        V(i3,1)=(1/Ybus(i3,i3))*(((P(i3,1)-Q(i3,1)*i)/(conj(V(i3,1))))+...
            (Ybus(i3,i3)*V(i3,1))-YV1);% Vi calculation
        Del(i3,1)=angle(V(i3,1));% Updating Del Value 
    end
%------------------------------------------------------------------------%
    % Printing Bus Parameters after each iteration%
%------------------------------------------------------------------------%
    fprintf('Iteration number = %d\n\n',while_loop_counter);
    for m1=2:n
        fprintf('The Bus %d : V%d = %.4f p.u.; Del%d = %.4f deg; P%d = %.4f p.u.; Q%d = %.4f p.u.\n',m1,m1,abs(V(m1,1)),m1,radtodeg(Del(m1,1)),m1,P(m1,1),m1,Q(m1,1));
    end
    fprintf('\n\n'); % Keeping space between Different Sections
%------------------------------------------------------------------------%
    % Updating PV Bus Voltages to original magnitudes %
%------------------------------------------------------------------------%
    for i6=PQnum+2:n
        V(i6,1)=abs(SP(i6,6))*(cos(Del(i6,1))+sin(Del(i6,1))*i);
    end        
%--- ---------------------------------------------------------------------%
    % Calculating Volatage errors  %
%------------------------------------------------------------------------%
    ee1=0; % Initialising Count for Error Free Vi
    for i5=1:n
        e=abs(V(i5,1)-V(i5,2));
        if (e<=0.001)
            ee1=ee1+1;
        end
    end
    if (ee1==n) % Updating While Condition
        ee=0;
    else
        ee=1;
    end

%------------------------------------------------------------------------%    
    % Copying updated values to second column of individual matrices to %
    % have last iteration values  %
%------------------------------------------------------------------------%
    for i7=1:n
        V(i7,2)=V(i7,1);
        Del(i7,2)=Del(i7,1);
        P(i7,2)=P(i7,1);
        Q(i7,2)=Q(i7,1);
    end
end
  
%------------------------------------------------------------------------%

% Calculating Slack Bus Power %
%------------------------------------------------------------------------% 
YV6=0+0i;

for k6=1:n
    YV6=YV6+(Ybus(1,k6)*V(k6,1));
end

SlackbusPower=V(1,1)*conj(YV6);
RSlackbusPower=real(SlackbusPower);
ImSlackbusPower=imag(SlackbusPower);
P(1,1)=RSlackbusPower;
Q(1,1)=ImSlackbusPower;
%------------------------------------------------------------------------%

% Calculating Line Currents %
%------------------------------------------------------------------------%
I=zeros(n);% Initialising Ybus Matrix

for i8=1:n
    for j2=1:n
        if (i8==j2)
           continue;
        else
            for j3=1:l
                if ((A(j3,i8)==1 || A(j3,i8)==-1)&&(A(j3,j2)==1 || A(j3,j2)==-1))
                    I(i8,j2)=((V(i8,1))-V(j2,1)*(-Ybus(i8,j2)))+(V(i8,1)*LP(j3,3)); % Calculation of Line Current
                end
            end
        end
    end
end

%------------------------------------------------------------------------%

% Calculating Line Power Flow Matrix %
%------------------------------------------------------------------------%
Si=zeros(n); % Initialising Line Power Flow Matrix

for i9=1:n
    for j4=1:n
        if (i9==j4)
            continue;
        else
            Si(i9,j4)=V(i9,1)*conj(I(i9,j4)); % Calculating Line Power Flow Matrix
        end
    end
end

%------------------------------------------------------------------------%

% Calculating Line Losses %
%------------------------------------------------------------------------%

S=Si+Si';

%------------------------------------------------------------------------%
% Printing Bus Parameters At last Iteration %
%------------------------------------------------------------------------%
fprintf('Number of Iterations required to complete Load Flow Studies = %d\n\n',while_loop_counter);
    for m2=1:n
        if (m2==1)
           fprintf('The Slack Bus parameters  : V%d = %.4f p.u.; Del%d = %.4f deg; P%d = %.4f p.u.; Q%d = %.4f p.u.\n',m2,abs(V(m2,1)),m2,radtodeg(Del(m2,1)),m2,P(m2,1),m2,Q(m2,1));
        else
            fprintf('The Bus %d : V%d = %.4f p.u.; Del%d = %.4f deg; P%d = %.4f p.u.; Q%d = %.4f p.u.\n',m1,m1,abs(V(m1,1)),m1,radtodeg(Del(m1,1)),m1,P(m1,1),m1,Q(m1,1));
        end            
    end
    fprintf('\n\n'); % Keeping space between Different Sections
    
%------------------------------------------------------------------------%
% Printing Line Currents %
%------------------------------------------------------------------------%    
 for m3=1:n
     for n3=1:n
         if (n3<=m3)
             continue;
         else
             fprintf('The Line Current from Node%d to Node%d = (%.4f) + (%.4f)j\n',m3,n3,real(I(m3,n3)),imag(I(m3,n3)));
         end
     end
 end
 fprintf('\n\n'); %  Keeping space between Different Sections
 
%------------------------------------------------------------------------%
% Printing Line Power Flows %
%------------------------------------------------------------------------%
for m5=1:n
    for n5=1:n
        if (n5<=m5)
            continue;
        else
            fprintf('The Power from Node%d to Node%d = (%.4f) + (%.4f)i\n',m5,n5,real(Si(m5,n5)),imag(Si(m5,n5)));
        end
    end
end
fprintf('\n\n'); %  Keeping space between Different Sections
%------------------------------------------------------------------------%
% Printing Line Losses %
%------------------------------------------------------------------------%
for m4=1:n
    for n4=1:n
        if (n4<=m4)
            continue;
        else
            fprintf('The Power Loss In Line %d-%d = (%.4f) + (%.4f)i\n',m4,n4,real(S(m4,n4)),imag(S(m4,n4)));
        end
    end
end
%------------------------------------------------------------------------%
% END %
%------------------------------------------------------------------------%
 
 


