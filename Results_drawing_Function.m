function [Drawing_comparison,MSE_SINDy,MSE_SINDy_Constraint,MSE_UN_SINDy] = Results_drawing_Function(out,Pre_library, Reduced_library, SINDY_UN_coefficients, SINDY_coefficients, SINDY_constrained_coefficients, Symbolic_Lagrange_equations,Removed_Range,DataSet_Library,DataSet_Torque,UN_DataSet_Library,UN_DataSet_Torque,joints)

% This function is used for generatimg the figures that compare the the R-SIEL equations with SINDy, SINDy with constrained library (SINDy C_L), Euler-Lagrange equation, simulation.

% (Inputs)
% out: The output simulation data from SIMULINK (SCARA_simulink_Data, KUKA_simulink_Data, SCARA_simulink_Data)

% Pre_library: The prelimnary library (KUKA Libraries, SCARA Libraries, UR5 Libraries)

% Reduced_library: The constrained library (KUKA Libraries, SCARA Libraries, UR5 Libraries)

% SINDY_UN_coefficients: The coefficents of original SINDy (SINDY_results_KUKA, SINDY_results_SCARA, SINDY_results_UR5)

% SINDY_coefficients: The coefficents of SINDy C-L (SINDY_results_KUKA, SINDY_results_SCARA, SINDY_results_UR5)

% SINDY_constrained_coefficients: The coefficents of SINDy R-SIEL (SINDY_results_KUKA, SINDY_results_SCARA, SINDY_results_UR5)

% Symbolic_Lagrange_equations: The Euler-Lagrange equations in symbolic
% form (Symbolic_Lagrange_equations_KUKA, Symbolic_Lagrange_equations_SCARA, Symbolic_Lagrange_equations_UR5)

% Removed_Range: The removed range from the simulation data sets (The
% number of observations in the simulation data are 9951), this parameter is
% introduced in order to reduce the time of drawing.

% DataSet_Library: Data sets for the constrained library (Google drive link)

% DataSet_Torque: Data sets for the constrained torques (Google drive link)

% UN_DataSet_Library: Data sets for the constrained library (Google drive link)

% UN_DataSet_Torque: Data sets for the constrained torques (Google drive link)

% joints: "R" for revolute and "P" for prismatic

% (Outputs)
% Drawing_comparison: The comparison figures
% MSE_SINDy: Mean squared error for SINDy with constrained library (SINDy C-L)
% MSE_SINDy_Constraint: Mean squared error for R-SIEL
% MSE_UN_SINDy: Mean squared error for original SINDy



n=length(Reduced_library);
% Mean squared error calculation for each algorithm
MSE_SINDy=zeros(n,1);
MSE_SINDy_Constraint=zeros(n,1);
MSE_UN_SINDy=zeros(n,1);

for i=1:n

MSE_SINDy(i)=sum((DataSet_Torque(:,i)-DataSet_Library{i}*SINDY_coefficients{i}).^2)/(size(DataSet_Library{i},1));  % SINDy C-L

MSE_SINDy_Constraint(i)=sum((DataSet_Torque(:,i)-DataSet_Library{i}*SINDY_constrained_coefficients{i}).^2)/(size(DataSet_Library{i},1)); % R-SIEL

MSE_UN_SINDy(i)=sum((UN_DataSet_Torque(:,i)-UN_DataSet_Library{i}*SINDY_UN_coefficients{i}).^2)/(size(UN_DataSet_Library{i},1)); % SINDy

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equations (SINDy)
SINDy_UN_Equations=cell(1,n);

for i=1:n
    SINDy_UN_Equations{i}=[SINDY_UN_coefficients{i} Pre_library{i}(:,2)];
end
%%%
Equa_UN_SINDY=SINDy_UN_Equations;
for i=1:n
    Equa_UN_SINDY{i}=SINDy_UN_Equations{i}(double(SINDy_UN_Equations{i}(:,1))~=0,:);
end
%%%
SINDy_UN_Equations_sym=cell(1,n);
for i=1:n
    SINDy_UN_Equations_sym{i}=str2sym(Equa_UN_SINDY{i});
end
%%
Symbolic_SINDy_UN_equations=cell(1,n);
for i=1:n
    Symbolic_SINDy_UN_equations{i}=sum(prod(SINDy_UN_Equations_sym{i},2));
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Equations (SINDy C-L)
SINDy_Equations=cell(1,n);

for i=1:n
    SINDy_Equations{i}=[SINDY_coefficients{i} Reduced_library{i}(:,2)];
end
%%%
Equa_SINDY=SINDy_Equations;
for i=1:n
    Equa_SINDY{i}=SINDy_Equations{i}(double(SINDy_Equations{i}(:,1))~=0,:);
end
%%%
SINDy_Equations_sym=cell(1,n);
for i=1:n
    SINDy_Equations_sym{i}=str2sym(Equa_SINDY{i});
end
%%
Symbolic_SINDy_equations=cell(1,n);
for i=1:n
    Symbolic_SINDy_equations{i}=sum(prod(SINDy_Equations_sym{i},2));
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equations (R-SIEL)
Constrained_SINDy_Equations=cell(1,n);

for i=1:n
    Constrained_SINDy_Equations{i}=[SINDY_constrained_coefficients{i} Reduced_library{i}(:,2)];
end
%%%
Equa_COS_SINDY=Constrained_SINDy_Equations;
for i=1:n
    Equa_COS_SINDY{i}=Constrained_SINDy_Equations{i}(double(Constrained_SINDy_Equations{i}(:,1))~=0,:);
end
%%%
Constrained_SINDy_Equations_sym=cell(1,n);
for i=1:n
    Constrained_SINDy_Equations_sym{i}=str2sym(Equa_COS_SINDY{i});
end
%%
Symbolic_Constrained_SINDy_equations=cell(1,n);
for i=1:n
    Symbolic_Constrained_SINDy_equations{i}=sum(prod(Constrained_SINDy_Equations_sym{i},2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate the substitution functions
syms theta_1(t)

suffix=zeros(n,1);
for i=1:n
    suffix(i)=string(i);
end

for i=1:n
    if joints(i)=="R"
        VARIABLES(i)=str2sym("theta_"+suffix(i)+"(t)");
    else
        VARIABLES(i)=str2sym("d"+suffix(i)+"(t)"); 
    end
end

syms s
substitute_functions=s;
for i=1:n
substitute_functions(i)=diff(VARIABLES(i), t, t);
substitute_functions(i+n)=diff(VARIABLES(i), t);
substitute_functions(i+(2*n))=VARIABLES(i);
end

% Collect the simulation data from simulink
simulation_Data=cell(1,3);
ju=3;
for k=1:3
simulation_Data{k}=zeros(length(out.tout),n);
jump=0;
for i=1:n
simulation_Data{k}(:,i)=out.yout{ju+jump}.Values.Data;
jump=jump+3;
end
ju=ju-1;
end
simulation_Data_total=zeros(length(out.tout),3*n);
jump=n;
jump1=1;
for i=1:3
simulation_Data_total(:,jump1:jump)=simulation_Data{i};
jump=jump+n;
jump1=jump1+n;
end

% remove data
simulation_Data_total(Removed_Range,:)=[];
Time=out.tout;
Time(Removed_Range)=[];

% storage cells
Torque_simulation=cell(1,n);   
TorqueLagrange=cell(1,n);
TorqueSindy=cell(1,n);
Torque_UN_Sindy=cell(1,n);
TorqueCons_Sindy=cell(1,n);

for i=1:n
    bb=out.yout{3*n+i}.Values.Data;
    bb(Removed_Range)=[];
    Torque_simulation{i}=bb;
    TorqueLagrange{i}=zeros(length(Time),1);
    Torque_UN_Sindy{i}=zeros(length(Time),1);
    TorqueSindy{i}=zeros(length(Time),1);
    TorqueCons_Sindy{i}=zeros(length(Time),1);
end

% Generate values for drawing
parfor k=1:n
for i=1:length(Time)
    
       TorqueLagrange{k}(i)=subs(Symbolic_Lagrange_equations{k},substitute_functions,[simulation_Data_total(i,:)]);
       Torque_UN_Sindy{k}(i)=subs(Symbolic_SINDy_UN_equations{k},substitute_functions,[simulation_Data_total(i,:)]);
       TorqueSindy{k}(i)=subs(Symbolic_SINDy_equations {k},substitute_functions,[simulation_Data_total(i,:)]);
       TorqueCons_Sindy{k}(i)=subs(Symbolic_Constrained_SINDy_equations{k},substitute_functions,[simulation_Data_total(i,:)]);
       
end
end

% generate comparison figures
Drawing_comparison=tiledlayout(n,1);

for i=1:n

    if joints(i)=="R"
    nexttile;
    plot(Time,TorqueLagrange{i},'g',Time,Torque_simulation{i},'r:',Time,TorqueSindy{i},'m--',Time,Torque_UN_Sindy{i},'c.',Time,TorqueCons_Sindy{i},'b-.','LineWidth',2);
    
    if i==1
    legend('Euler-lagrange','Simulation','SINDy C-L','SINDy','R-SIEL')
    end

    xlabel('Time (s)','FontSize',12,'FontWeight','bold')
    ylabel(char("Torque of Joint-" +i +"(N.m)"),'FontSize',12,'FontWeight','bold')
    else
    nexttile;
    plot(Time,TorqueLagrange{i},'g',Time,Torque_simulation{i},'r:',Time,TorqueSindy{i},'m--',Time,Torque_UN_Sindy{i},'c.',Time,TorqueCons_Sindy{i},'b-.','LineWidth',2);
    
    if i==1
    legend('Euler-lagrange','Simulation','SINDy C-L','SINDy','R-SIEL')
    end
    
    xlabel('Time (s)','FontSize',12,'FontWeight','bold')
    ylabel(char("Force of Joint-" +i +"(N)"),'FontSize',12,'FontWeight','bold')        
    end

end
end
