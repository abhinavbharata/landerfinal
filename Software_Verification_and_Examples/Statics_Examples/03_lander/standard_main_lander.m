 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%Static calculation of a tensegrity tower %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;close all;
% Specify material properties
[consti_data,Eb,Es,sigmab,sigmas,rho_b,rho_s]=material_lib('PLA','Ninja-Flex');
material{1}='linear_elastic'; % index for material properties:'linear_elastic' multielastic plastic
material{2}=0; % index for considering slack of string (1) for yes,(0) for no (for compare with ANSYS)

% cross section design cofficient
thick=6e-3;        % thickness of hollow bar
hollow_solid=0;          % use hollow bar or solid bar in minimal mass design (1)hollow (0)solid
c_b=0.1;           % coefficient of safty of bars 0.5
c_s=0.1;           % coefficient of safty of strings 0.3

substep=100;                                     %substep
lumped=0;               % use lumped matrix 1-yes,0-no
saveimg=0;              % save image or not (1) yes (0)no
savedata=0;             % save data or not (1) yes (0)no
savevideo=0;            % make video(1) or not(0)
gravity=0;              % consider gravity 1 for yes, 0 for no
savePath=fullfile(fileparts(mfilename('fullpath')),'data_temp'); %Save files in same folder as this code

%% N C of the structure
% Manually specify node positions of a lander.
L = 1.45;
th1 = -atan(1/2);
N = [
%% Stage-one
0.3625 0 0;             0.3625 0 1.45;           -0.3625 0 0;            -0.3625 0 1.45;
0 -0.7250 1.0875;       0 0.7250 1.0875;          0 -0.7250 0.3625;       0 0.7250 0.3625;
0.7250 -0.3625 0.7250; -0.7250 -0.3625 0.7250;    0.7250 0.3625 0.7250;  -0.7250 0.3625 0.7250;

%% Stage-two
1.1065 -0.7309 0.7459;  1.1065 -0.7309 2.1959;   0.3815 -0.7309 0.7459;   0.3815 -0.7309 2.1959;
0.7440 -1.4559 1.8334;  0.7440 -0.0059 1.8334;   0.7440 -1.4559 1.1084;   0.7440 -0.0059 1.1084;
1.4690 -1.0934 1.4709;  0.0190 -1.0934 1.4709;   1.4690 -0.3684 1.4709;   0.0190 -0.3684 1.4709;

%% Stage-three
1.8553 -1.4682 1.4952; 1.8553 -1.4682 2.9452; 1.1303 -1.4682 1.4952; 1.1303 -1.4682 2.9425;
1.4925 -2.1932 2.5827; 1.4925 -0.7432 2.5827; 1.4925 -2.1932 1.8577; 1.4925 -0.7432 1.8577;
2.2178 -1.8307 2.2202; 0.7678 -1.8307 2.2202; 2.2178 -1.1057 2.2202; 0.7678 -1.1057 2.2202
]';
N = [1 0 0;0 cos(th1) -sin(th1);0 sin(th1) cos(th1)]*N;
N(3,:)=N(3,:)+2;
% N = [cos(th2) -sin(th2) 0;sin(th2) cos(th2) 0;0 0 1]*N;
C_b_in = [
          1 2;3 4;5 6;7 8;9 10;11 12;
          13 14;15 16;17 18;19 20;21 22;23 24;
          25 26;27 28;29 30;31 32;33 34;35 36
          ];

C_s_in = [
 2 6;2 11;4 5;4 6;4 10;4 12;
 1 7;1 8;1 9;1 11;3 7;3 8;3 10;3 12;
 5 10;7 9;7 10;6 11;6 12;8 11;8 12;
 2 20;2 24;9 20;9 15;5 15;5 24;
 %%
 14 18;14 23;16 17;16 18;16 22;
 16 24;13 19;13 20;13 21;13 23;15 19;15 22;
 17 22;19 21;19 22;18 23;18 24;20 23;
 14 32;14 36;21 32;21 27;17 27;17 36;
 %% 
 26 30;26 35;26 29;26 33;28 29;28 30;28 34;28 36;25 31;
 25 32;25 33;25 35;27 31;27 34;29 33;29 34;
 31 33;31 34;30 35;30 36;32 35
];
% Convert the above matrices into full connectivity matrices.
C_b = tenseg_ind2C(C_b_in,N);%%
C_s = tenseg_ind2C(C_s_in,N);
C=[C_b;C_s];
[ne,nn]=size(C);        % ne:No.of element;nn:No.of node
% Plot the structure to make sure it looks right
tenseg_plot(N,C_b,C_s);

%% Boundary constraints
pinned_X=([])'; pinned_Y=([])'; pinned_Z=([])';
[Ia,Ib,a,b]=tenseg_boundary(pinned_X,pinned_Y,pinned_Z,nn);

%% Group information
%generate group index
% gr={[1:6];[7:ne]};       % number of elements in one group
gr=[];                     % if no group is used
Gp=tenseg_str_gp(gr,C);    % generate group matrix

%% self-stress design
%Calculate equilibrium matrix and member length
[A_1a,A_1ag,A_2a,A_2ag,l,l_gp]=tenseg_equilibrium_matrix1(N,C,Gp,Ia);

%SVD of equilibrium matrix
[U1,U2,V1,V2,S]=tensegL_svd(A_1ag);

%external force in equilibrium design
w0=zeros(numel(N),1); w0a=Ia'*w0;

%prestress design
index_gp=1;                   % number of groups with designed force
fd=-1e6;                        % force in bar is given as -1000
[q_gp,t_gp,q,t]=tenseg_prestress_design(Gp,l,l_gp,A_1ag,V2,w0a,index_gp,fd);    %prestress design

%% cross sectional design
index_b=find(t<0);              % index of bar in compression
index_s=setdiff(1:ne,index_b);	% index of strings
[A_b,A_s,A_gp,A,r_b,r_s,r_gp,radius,E,l0,rho,mass]=tenseg_minimass(t,l,Gp,sigmas,sigmab,Eb,Es,index_b,index_s,c_b,c_s,rho_b,rho_s,thick,hollow_solid);

% Plot the structure with radius
R3Ddata.Bradius=0.01*interp1([min(radius),max(radius)],[0.2,0.8],r_b);
R3Ddata.Sradius=0.01*interp1([min(radius),max(radius)],[0.2,0.8],r_s);
R3Ddata.Nradius=0.01*ones(nn,1);
tenseg_plot(N,C_b,C_s,[],[],[],'Tensegrity Lander',R3Ddata);

%% mass matrix and damping matrix
M=tenseg_mass_matrix(mass,C,lumped); % generate mass matrix

% damping matrix
d=0; % damping coefficient
D=d*2*max(sqrt(mass.*E.*A./l0))*eye(3*nn); % critical damping

%% mode analysis
num_plt=7:9; % mode index to plot
[V_mode,omega]=tenseg_mode(Ia,C,C_b,C_s,N(:),E,A,l0,M,num_plt,saveimg);

%% external force, forced motion of nodes, shrink of strings
% calculate external force and
ind_w=[];w=[];
ind_dnb=[]; dnb0=[];
% ind_dl0=[1;2]; dl0=[1;1];   %extent rest length of bar
ind_dl0=[7;8;9;10]; dl0=[-0.5;-0.5;-0.5;-0.5];   %extent rest length of bar
[w_t,dnb_t,l0_t,Ia_new,Ib_new]=tenseg_load_prestress(substep,ind_w,w,ind_dnb,dnb0,ind_dl0,dl0,l0,b,gravity,[0;9.8;0],C,mass);

%% input file of ANSYS
ansys_input_gp_00(N,C,A_gp,t_gp,b,Eb,Es,rho_b,rho_s,Gp,index_s,find(t_gp>0),ind_w,w,ind_dnb,dnb0,fullfile(savePath,'lander_ansys'));
% ansys_input_gp(N,C,A_gp,t_gp,b,Eb,Es,rho_b,rho_s,Gp,index_s,find(t_gp>0),'lander');

%% equilibrium calculation
% input data
data.N=N; data.C=C; data.ne=ne; data.nn=nn; data.Ia=Ia_new; data.Ib=Ib_new;
data.E=E; data.A=A; data.l0=l0; data.index_b=index_b; data.index_s=index_s;
data.consti_data=consti_data;   data.material=material; %constitue info
data.w_t=w_t;  % external force
data.dnb_t=dnb_t;% forced movement of pinned nodes
data.l0_t=l0_t;% forced movement of pinned nodes
data.substep=substep;    % substep

% nonlinear analysis
data_out=static_solver(data);        %solve equilibrium using mNewton method
% data_out=static_solver2(data);        %solve equilibrium using mNewton method
% data_out{i}=equilibrium_solver_pinv(data);        %solve equilibrium using mNewton method

t_t=data_out.t_out;          %member force in every step
n_t=data_out.n_out;          %nodal coordinate in every step
N_out=data_out.N_out;

%% plot member force
tenseg_plot_result(1:substep,t_t(1:3,:),{'element 1','element 2','element 3'},{'Load step','Force (N)'},fullfile(savePath,'plot_member_force.png'),saveimg);

%% Plot nodal coordinate curve X Y
tenseg_plot_result(1:substep,n_t([3*4-2,3*4],:),{'4X','4Z'},{'Time (s)','Coordinate (m)'},fullfile(savePath,'plot_coordinate.png'),saveimg);

%% Plot final configuration
% tenseg_plot_catenary(reshape(n_t(:,end),3,[]),C_b,C_s,[],[],[0,0],[],R3Ddata,l0_t(index_s,end));

%% save output data
if savedata==1
    save (fullfile(savePath,['tensegrity_ball',material{1},'.mat']));
end

%% make video of the dynamic
name=fullfile(savePath,['lander',material{1},num2str(material{2})]);
% tenseg_video(n_t,C_b,C_s,[],substep,name,savevideo);
tenseg_video_slack(n_t,C_b,C_s,l0_t,index_s,[],[],[],min(substep,50),name,savevideo,material{2});

%% linearized dynaimcs
[A_lin,B_lin]=tenseg_lin_mtrx(C,N(:),E,A,l0,M,D,Ia,A_1a);
