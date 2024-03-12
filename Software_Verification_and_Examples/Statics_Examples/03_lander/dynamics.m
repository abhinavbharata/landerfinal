clc;clear;close all;
% Global variable
[consti_data,Eb,Es,sigmab,sigmas,rho_b,rho_s]=material_lib('PLA','Ninja-Flex');
material{1}='linear_elastic'; % index for material properties: multielastic, plastic.
material{2}=0; % index for considering slack of string (1) for yes,(0) for no (for comparision with ANSYS)

% cross section design cofficient
thick=6e-3;        % thickness of hollow bar
hollow_solid=0;    % use hollow bar or solid bar in minimal mass design (1) hollow (0)solid
c_b=0.1;           % coefficient of safty of bars 0.5
c_s=0.1;           % coefficient of safty of strings 0.3

% dynamic analysis set
amplitude=50;           % amplitude of external force of ground motion
period=0.5;             % period of seismic

dt=1e-5;                % time step in dynamic simulation
auto_dt=1;              % use(1 or 0) auto time step, converengency is guaranteed if used
tf=2;                   % final time of dynamic simulation
out_dt=0.001;           % output data interval(approximately, not exatly)
lumped=0;               % use lumped matrix 1-yes,0-no
saveimg=0;              % save image or not (1) yes (0)no
savedata=0;             % save data or not (1) yes (0)no
savevideo=1;            % make video(1) or not(0)
gravity=0;              % consider gravity 1 for yes, 0 for no
% move_ground=0;        % for earthquake, use pinned nodes motion(1) or add inertia force in free node(0)
savePath=fullfile(fileparts(mfilename('fullpath')),'data_temp'); %Save files in same folder as this code

%% N C of the structure
% Manually specify node positions of double layer prism.
L = 0.145;
th1 = -atan(1/2);
N = [ %% First Stage
0.03625	0	0;0.03625	0	0.145;-0.03625	0	0;-0.03625	0	0.145;
0	-0.0725	0.10875;0	0.0725	0.10875;0	-0.0725	0.03625;0	0.0725	0.03625;
0.0725	-3.63E-02	0.0725;-0.0725	-3.63E-02	0.0725;0.0725	0.03625	0.0725;-0.0725	0.03625	0.0725;
%% Second Stage
0.11065	-0.07309	0.07459;0.11065	-0.07309	0.21959;0.03815	-0.07309	0.07459;0.03815	-0.07309	0.21959;
0.0744	-0.14559	0.18334;0.0744	-0.00059	0.18334;0.0744	-0.14559	0.11084;0.0744	-0.00059	0.11084;
0.1469	-0.10934	0.14709;0.0019	-0.10934	0.14709;0.1469	-0.03684	0.14709;0.0019	-0.03684	0.14709;
%% Third Stage
0.18553	-0.14682	0.14952;0.18553	-0.14682	0.29452;0.11303	-0.14682	0.14952;0.11303	-0.14682	0.29425;
0.14925	-0.21932	0.25827;0.14925	-0.07432	0.25827;0.14925	-0.21932	0.18577;0.14925	-0.07432	0.18577;
0.22178	-0.18307	0.22202;0.07678	-0.18307	0.22202;0.22178	-0.11057	0.22202;0.07678	-0.11057	0.22202
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
title('A tensegrity lander');

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
R3Ddata.Bradius=0.01*interp1([min(radius),max(radius)],[0.2,0.5],r_b);
R3Ddata.Sradius=0.01*interp1([min(radius),max(radius)],[0.05,0.1],r_s);
R3Ddata.Nradius=0.005*ones(nn,1);
tenseg_plot(N,C_b,C_s,[],[],[],'Tensegrity Lander',R3Ddata);

%% input file of ANSYS
ansys_input_gp(N,C,A_gp,t_gp,b,Eb,Es,rho_b,rho_s,Gp,index_s,find(t_gp>0),fullfile(savePath,'lander'));

%% mass matrix and damping matrix
M=tenseg_mass_matrix(mass,C,lumped); % generate mass matrix
% damping matrix
d=0; % damping coefficient
D=d*2*max(sqrt(mass.*E.*A./l0))*eye(3*nn); % critical damping

%% mode analysis
num_plt=7:9; % mode index to plot
[V_mode,omega]=tenseg_mode(Ia,C,C_b,C_s,N(:),E,A,l0,M,num_plt,saveimg);

%% input data for dynamic analysis
% time step
if auto_dt
    dt=pi/(8*max(omega)); % time step dt is 1/8 of the smallest period, guarantee convergence in solving ODE
end
tspan=0:dt:tf;

out_tspan=interp1(tspan,tspan,0:out_dt:tf, 'nearest','extrap');  % output data time span

% calculate external force and forced motion of nodes
[~,~,~,~,dz_a_t]=tenseg_ex_force(tspan,b,'impulse',gravity,[0;0;9.8],C,mass,3*(4:9)-1,0,period);
% [w_t,dnb_t,dnb_d_t,dnb_dd_t,dz_a_t]=tenseg_ex_force(tspan,a,b,'step',gravity,C,mass,3*(4:9)-2,1e5,period);
% [w_t,dnb_t,dnb_d_t,dnb_dd_t,dz_a_t]=tenseg_ex_force(tspan,a,b,'ramp',gravity,C,mass,3*(4:9),-1e6,period);
% [w_t,dnb_t,dnb_d_t,dnb_dd_t,dz_a_t]=tenseg_ex_force(tspan,a,b,'vib_force',gravity,C,mass,3*(4:9)-2,50,period);
% [w_t,dnb_t,dnb_d_t,dnb_dd_t,dz_a_t]=tenseg_ex_force(tspan,a,b,'vib_nodes',gravity,C,mass,3*(1:3)-2,50,period);
w_t='w_t_ground_force';         % external force given in a script
dnb_t=zeros(b,1);dnb_d_t=zeros(b,1);dnb_dd_t=zeros(b,1); % zero forced motion

% give initial speed of free coordinates
n0a_d=4*kron(ones(numel(a)/3,1),[0;0;-1]); %initial speed in X direction
% n0a_d(1:end/2)=-n0a_d(1:end/2);

% input data
data.N=N; data.C=C; data.ne=ne; data.nn=nn; data.Ia=Ia; data.Ib=Ib;
data.E=E; data.A=A; data.l0=l0; data.index_b=index_b; data.index_s=index_s;
data.consti_data=consti_data;   data.material=material; %constitue info
data.w_t=w_t; % external force
data.dnb_t=dnb_t; data.dnb_d_t=dnb_d_t;  data.dnb_dd_t=dnb_dd_t; % forced movement of pinned nodes
data.n0a_d=n0a_d; %initial speed of free coordinates
data.M=M;data.D=D;
data.tf=tf;data.dt=dt;data.tspan=tspan;data.out_tspan=out_tspan;

%% dynamic analysis
% solve dynamic equation
data_out=dynamic_solver(data);        %solve ODE of dynamic equation
% time history of structure
t_t=data_out.t_t;   %time history of members' force
n_t=data_out.n_t;   %time history of nodal coordinate
l_t=data_out.l_t;   %time history of members' length
nd_t=data_out.nd_t;   %time history of nodal coordinate

%% save output data
if savedata==1
    save (['FEM_lander_',material{1},'.mat']);
end

%% plot member force
tenseg_plot_result(out_tspan,t_t(7:8,:),{'element 7','element 8'},{'Time (s)','Force (N)'},fullfile(savePath,'member_force.png'),saveimg);

%% Plot nodal coordinate curve X Y
tenseg_plot_result(out_tspan,n_t(3*8,:),{'8Z'},{'Time (s)','Coordinate (m)'},fullfile(savePath,'Z_coordinate.png'),saveimg);
tenseg_plot_result(out_tspan,n_t(3*8-1,:),{'8Y'},{'Time (s)','Coordinate (m)'},fullfile(savePath,'Y_coordinate.png'),saveimg);
tenseg_plot_result(out_tspan,n_t(3*8-2,:),{'8X'},{'Time (s)','Coordinate (m)'},fullfile(savePath,'X_coordinate.png'),saveimg);
%% plot external force information
% tenseg_plot_exforce(Ib,tspan,w_t,[4:6],dnb_t,dnb_d_t,dnb_dd_t,[1],saveimg);

%% make video of the dynamic
name=fullfile(savePath,['lander','tf_',num2str(tf),material{1}]);
% tenseg_video(n_t,C_b,C_s,[],100,name,savevideo);

nx_t=n_t(1:3:end); min_x=min(min(nx_t));max_x=max(max(nx_t));d_x=max_x-min_x+1e-3;
ny_t=n_t(2:3:end); min_y=min(min(ny_t));max_y=max(max(ny_t));d_y=max_y-min_y+1e-3;
nz_t=n_t(3:3:end); min_z=min(min(nz_t));max_z=max(max(nz_t));d_z=max_z-min_z+1e-3;
d=max(max(n_t))-min(min(n_t));
axislim=[min_x-0.2*d_x,max_x+0.2*d_x,min_y-0.2*d_y,max_y+0.2*d_y,min_z,max_z+0.2*d_z];

if savevideo==1
    figure(99);
    set(gcf,'Position',get(0,'ScreenSize'));
    for p = 1:floor(size(n_t,2)/100):size(n_t,2)
        tenseg_plot(reshape(n_t(:,p),3,[]),C_b,C_s,99);hold on
        axis(axislim)
        [X,Y] = meshgrid(-0.8:0.05:0.8,-0.5:0.05:1);
        Z = -0.1*ones(size(X));
        surf(X,Y,Z,'MeshStyle','none','FaceColor','c');
        surf(X,Y,Z);
        tenseg_savegif_forever(name);
        hold off;
    end
    close
end

%% linearized dynamics
[A_lin,B_lin]=tenseg_lin_mtrx(C,N(:),E,A,l0,M,D,Ia,A_1a);