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
savevideo=1;            % make video(1) or not(0)
gravity=0;              % consider gravity 1 for yes, 0 for no
savePath=fullfile(fileparts(mfilename('fullpath')),'data_temp'); %Save files in same folder as this code

%% N C of the structure
% Manually specify node positions of a lander.
L = 1.45;
th1 = -atan(1/2);

N=[
%% FIRST stage one and stage two 
-0.7754	-0.1008	-0.023;-0.6855	0.0155	-0.0172;-0.7952	-0.0728	-0.0887;-0.7059	0.0437	-0.0793;
-0.786	0.0214	-0.1057;-0.7449 -0.0329	  0.035; -0.7414 -0.0203 -0.1339;-0.6956 -0.0724 -0.0116;
-0.8039	0.0125	-0.0177;-0.7142	-0.1022	-0.093; -0.7647	0.0409	-0.0086;-0.6704	-0.0441	-0.0849;
%% SECOND
-0.6632	-0.1026	0.0178;-0.5765	0.0167	0.0153;-0.6819	-0.0817	-0.0511;-0.5972	0.0399	-0.0514;
-0.6733	0.0076	-0.0673;-0.6345	-0.0338	0.0718;-0.626	-0.0281	-0.1014;-0.5851	-0.0706	0.0301;
-0.6985	-0.0234	0.0199;-0.5964	-0.1014	-0.0563;-0.6569	0.0316	0.0232;-0.5531	-0.0389	-0.0573;
 
  ]';
N = [1 0 0;0 cos(th1) -sin(th1);0 sin(th1) cos(th1)]*N;
N(3,:)=N(3,:)+2;
C_b_in = [
    1 2;3 4;5 6;7 8;9 10;11 12;13 14;15 16;17 18;
    19 20;21 22;23 24;
    %25 26;27 28;29 30;31 32;33 34;35 36 
          ];
C_s_in = [ 
    1 6;1 8;1 9;1 10;2 6;2 11;2 17;2 21;3 5;3 7;3 9;3 10;4 5;4 7;4 11;4 12;
    5 9;5 11;6 9;6 11;7 10;7 12;8 10;8 15;8 21;12 15;12 17;13 18;13 20;13 21;
    13 22;14 18;14 23;15 19;15 22;16 17;16 19;16 23;16 24;17 23;14 24;14 20;20 22;20 24;18 21;18 23;19 22;19 24
%%
%     1 6;1 8;1 9;1 10;2 6;2 11;2 17;2 21;3 5;3 7;3 9;3 10;4 5;4 7;4 11;4 12;
%     5 9;5 11;6 9;6 11;7 10;7 12;8 10;8 15;8 21;12 15;12 17;13 18;13 20;13 21;
%     13 22;14 18;14 23;15 19;15 22;16 17;16 19;16 23;16 24;17 23;20 22;18 21;
%     18 23;19 22;19 24;
%     14 33;14 29;20 33;20 27;24 29;29 35;29 28;24 27;30 25;30 26;30 33;30 35;
%     31 34;31 36;32 34;32 36;31 27;31 28;32 25;32 26;28 35;28 36;26 35;26 36;
%     25 33;25 34;27 28;

%%    
%     1 6;1 8;1 9;1 10;
%     2 6;2 11;2 17;2 21;
%     3 5;3 7;3 9;3 10;
%     4 5;4 7;4 11;4 12;
%     5 9;5 11;
%     6 9;6 11;
%     7 10;7 12;
%     8 10;8 15;8 21;
%     12 15;12 17;
%     13 18;13 20;13 21;13 22;
%     14 18;14 29;14 23;14 33;
%     15 19;15 22;
%     16 17;16 19;16 23;16 24;
%     17 23;18 21;18 23;
%     19 22;19 24;
%     20 22;20 27;20 33;
%     24 27;24 29;
%     25 30;25 32;25 33;25 34;
%     26 32;26 35;26 36;
%     27 31;27 34;
%     28 29;28 31;28 35;28 36;
%     29 35;
%     30 26;30 33;30 35;
%     31 34;31 36;
%     32 34;32 36
];
% Convert the above matrices into full connectivity matrices.
C_b = tenseg_ind2C(C_b_in,N);%%
C_s = tenseg_ind2C(C_s_in,N);
C=[C_b;C_s];
[ne,nn]=size(C);        % ne:No.of element;nn:No.of node
% Plot the structure to make sure it looks right
tenseg_plot(N,C_b,C_s);
grid on
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
fd=-1e2;                        % force in bar is given as -1000
[q_gp,t_gp,q,t]=tenseg_prestress_design(Gp,l,l_gp,A_1ag,V2,w0a,index_gp,fd);    %prestress design

%% cross sectional design
index_b=find(t<0);             % index of bar in compression
index_s=setdiff(1:ne,index_b);	% index of strings
[A_b,A_s,A_gp,A,r_b,r_s,r_gp,radius,E,l0,rho,mass]=tenseg_minimass(t,l,Gp,sigmas,sigmab,Eb,Es,index_b,index_s,c_b,c_s,rho_b,rho_s,thick,hollow_solid);


% Plot the structure with radius
R3Ddata.Bradius=interp1([min(radius),max(radius)],[0.0005,0.001],r_b);
R3Ddata.Sradius=interp1([min(radius),max(radius)],[0.0001,0.0002],r_s);
R3Ddata.Nradius=0.01*ones(nn,1);
tenseg_plot(N,C_b,C_s,[],[],[],'Tensegrity Lander',R3Ddata);
grid on
% %% mass matrix and damping matrix
% M=tenseg_mass_matrix(mass,C,lumped); % generate mass matrix
% 
% % damping matrix
% d=0; % damping coefficient
% D=d*2*max(sqrt(mass.*E.*A./l0))*eye(3*nn); % critical damping
% 
% %% mode analysis
% num_plt=7:9; % mode index to plot
% [V_mode,omega]=tenseg_mode(Ia,C,C_b,C_s,N(:),E,A,l0,M,num_plt,saveimg);
% 
% %% external force, forced motion of nodes, shrink of strings
% % calculate external force and
% ind_w=[];w=[];
% ind_dnb=[]; dnb0=[];
% % ind_dl0=[1;2]; dl0=[1;1];   %extent rest length of bar
% ind_dl0=[7;8;9;10]; dl0=[-0.5;-0.5;-0.5;-0.5];   %extent rest length of bar
% [w_t,dnb_t,l0_t,Ia_new,Ib_new]=tenseg_load_prestress(substep,ind_w,w,ind_dnb,dnb0,ind_dl0,dl0,l0,b,gravity,[0;9.8;0],C,mass);
% 
% %% input file of ANSYS
% ansys_input_gp_00(N,C,A_gp,t_gp,b,Eb,Es,rho_b,rho_s,Gp,index_s,find(t_gp>0),ind_w,w,ind_dnb,dnb0,fullfile(savePath,'lander_ansys'));
% % ansys_input_gp(N,C,A_gp,t_gp,b,Eb,Es,rho_b,rho_s,Gp,index_s,find(t_gp>0),'lander');
% 
% %% equilibrium calculation
% % input data
% data.N=N; data.C=C; data.ne=ne; data.nn=nn; data.Ia=Ia_new; data.Ib=Ib_new;
% data.E=E; data.A=A; data.l0=l0; data.index_b=index_b; data.index_s=index_s;
% data.consti_data=consti_data;   data.material=material; %constitue info
% data.w_t=w_t;  % external force
% data.dnb_t=dnb_t;% forced movement of pinned nodes
% data.l0_t=l0_t;% forced movement of pinned nodes
% data.substep=substep;    % substep
% 
% % nonlinear analysis
% data_out=static_solver(data);        %solve equilibrium using mNewton method
% % data_out=static_solver2(data);        %solve equilibrium using mNewton method
% % data_out{i}=equilibrium_solver_pinv(data);        %solve equilibrium using mNewton method
% 
% t_t=data_out.t_out;          %member force in every step
% n_t=data_out.n_out;          %nodal coordinate in every step
% N_out=data_out.N_out;
% 
% %% plot member force
% tenseg_plot_result(1:substep,t_t(1:3,:),{'element 1','element 2','element 3'},{'Load step','Force (N)'},fullfile(savePath,'plot_member_force.png'),saveimg);
% 
% %% Plot nodal coordinate curve X Y
% tenseg_plot_result(1:substep,n_t([3*4-2,3*4],:),{'4X','4Z'},{'Time (s)','Coordinate (m)'},fullfile(savePath,'plot_coordinate.png'),saveimg);
% 
% %% Plot final configuration
% tenseg_plot_catenary(reshape(n_t(:,end),3,[]),C_b,C_s,[],[],[0,0],[],R3Ddata,l0_t(index_s,end));
% 
% %% save output data
% if savedata==1
%     save (fullfile(savePath,['tensegrity_ball',material{1},'.mat']));
% end
% 
% %% make video of the dynamic
% name=fullfile(savePath,['lander',material{1},num2str(material{2})]);
% % tenseg_video(n_t,C_b,C_s,[],substep,name,savevideo);
% tenseg_video_slack(n_t,C_b,C_s,l0_t,index_s,[],[],[],min(substep,50),name,savevideo,material{2});
% 
% %% linearized dynaimcs
%  [A_lin,B_lin]=tenseg_lin_mtrx(C,N(:),E,A,l0,M,D,Ia,A_1a);
% 
