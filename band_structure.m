%% band_structure
%
% script to caulculate band structure from simple tight-binding
% Hamiltonians. Hamiltonians must be defined for each lattice. Variables of
% the matrix must be kx and ky. 
%
% contains sections for infinite 2D lattices, semi-infinite lattices 
% (ribbons) and Floquet topological insulators. The sections to be executed
% can be selected in input parameters via booleans. Currently only for the
% Lieb-kagome transition lattice ('trans') all functions are available.

tic

%   clear all
% clear imin zmin diracp autoDP cones EV ev
clearvars -except cvar vs8 conf8 vs conf vs1 vs2 conf1 conf2 conemovement coneenergy icm theta nexp startcol scanvar kp
close all %-except Figure 999

global bound1 bound2 stepsize steps ew lpos D N point_names points plotcount; % variables needed for function Bandstruktur1D

folder='I:\Desktop\Simulation - phi\save-load\'; % folder used to save and load files (like positions of DP)

%% ---------------------------input parameters-----------------------------

lattice='trans'; %choose lattice.   trans, transition, transition2, transition3, transition4, trans_st, trans5, dislocated, rect, graphene, check, disloc7, trans_qr, trans_q, trans_i, sq1, sq2 
calculate_lattice=true; %manually decide if sections depending on 2D Hamiltonian are executed.

t=1; % coupling constant, usually 1
a=1; % lattice constant, usually 1

theta=15; % shearing angle for transition lattices
nexp=4; % exponent of nnn interaction: larger nexp corresponds to waeker nnn interaction

delta=3e-3; % additional dislocation of some elements in trans_st lattice. 0<delta<0.5. should be small
gamma=0.5; % used in trans_qr lattices

steps=2000; % ~ 500-3000 depending on RAM, graphics card and desired accuracy
bound1=-2*pi; % define range of k for bandstructure calculation
bound2=2*pi; % 2*pi/a 1.8*pi 0-1.5*pi
stepsize=abs(bound1-bound2)/steps;

autosaveDP=false; % decides if section autosaveDP is executed. Will overwrite files if they already exist
calculateEV=true; % decides if section calculateEV is executed.

cutoff=15; %15 % used for calculating DP: Local minima *cutoff* times larger than the global minimum are ignored. Might have to be adjusted sometimes

if autosaveDP
    warning('autosave activated')
end

%% -----------------------ribbon input parameters--------------------------

ribbon_size=15; %width of the ribbon in unit cells in finite direction. Also used for Floquet ribbon.

ribbon=lattice; %choose lattice for ribbon calculation.   trans, lieb_y_bearded, trans2_flat, kagome3_mixed, trans4_st_mixed, trans3_mixed, kagome2_flat, trans2_mixed, trans2_st_bearded, trans_st_bearded, graphene_zigzag2, graphene_zigzag3, lieb_st_y_bearded, graphene_y_zigzag, trans_y_bearded, kagome_x, kagome2_x, kagome3_x, trans2
calculate_ribbon=true; %manually decide if sections depending on ribbon Hamiltonian are executed.

steps_r=1000;
bound1_r=0*pi/a;
bound2_r=2*pi/a;
stepsize_r=abs(bound1_r-bound2_r)/steps_r;

gamma0=0.2886;

%% -------------------Floquet ribbon input parameters----------------------

floquet=lattice; %choose lattice for Floquet ribbon calculation.   trans, graphene_y_zigzag, trans4_mixed
calculate_floquet = true; %manually decide if sections depending on Floquet ribbon Hamiltonian are executed.

a_f=20e-6;

steps_f=2000;
bound1_f=0*pi/a_f;
bound2_f=2*pi/a_f;
stepsize_f=abs(bound1_f-bound2_f)/steps_f;

% Helix parameters
Z=1e-2; % Z=1e-2; % helix period [m]
omega=2*pi/Z; % angular frequency of the helix
omega2=omega; % set to 1 or omega (omega is probably correct)
R=0.4*a_f; % 0.32*a % helix radius
z1=0; % time evolution operator is calculated over one helix period
z2=Z;
t_resolution=100; % number of steps for calculation of time evolution operator
zstep=Z/t_resolution;

% parameters needed for peierl's substitution
n0=1.45; % refractive index of the medium
lambda=633e-9; %633e-9 % wavelength of excitation light [m]
k0=2*pi*n0/lambda; % Influences the Peierl's substitution

%% ---------------- check if all Hamiltonians are defined -----------------
% checks for the selected lattice if 2D, ribbon and Floquet ribbon
% Hamiltonians are defined. If any are undefined the respective sections
% are skipped.

[valid_lattice, valid_ribbon, valid_floquet] = Valid(lattice, ribbon, floquet);

%% ----------------------lattice-specific parameters-----------------------
% input parameters which depend on the lattice:
% - N: Dimension of the Hamiltonian, mandatory
% - point_names, points: names and coordinates of high-symmetry points.
%   required for Bandstruktur1D. If not set, 1D bandstructure won't be
%   calculated.

if calculate_lattice && valid_lattice % check if following sections shall/can be executed
    
    point_names=[];
    points=[];
    
    if strcmp(lattice,'transition') || strcmp(lattice,'transition2') || strcmp(lattice,'transition3') || strcmp(lattice,'transition4') || strcmp(lattice,'trans') || strcmp(lattice,'trans_st') || strcmp(lattice,'trans5') || strcmp(lattice,'trans_qr') || strcmp(lattice,'trans_q') || strcmp(lattice,'trans_i') || strcmp(lattice,'sq1') || strcmp(lattice,'sq2')
        
        if strcmp(lattice,'trans5')
            N=5;
        elseif strcmp(lattice,'sq1')
            N=1;
        elseif strcmp(lattice,'sq2')
            N=2;
        else
            N=3;
        end
        
        load([folder, 'elsepoints.mat']);
        load([folder, 'elsepoints05.mat']);
        point_names=["\Gamma","M","K","X"];
        
        if theta==0
            point_names=["\Gamma","M","X"];
            points=[[0, 0];[pi, pi];[pi, 0]];
        elseif theta==5
            points=[[0,0];[3.14159265358979,2.87873928454845];[3.39344970623451,2.60388554194302];[3.14159265358979,-0.274853742605432]];
        elseif theta==10
            points=[[0,0];[3.14159265358979,2.63610923693645];[3.60640983615649,2.08216168971307];[3.14159265358979,-0.553947547223382]];
        elseif theta==15
            points=[[0,0];[3.14159265358979,2.41062882833589];[3.78751870139358,1.56884161385895];[3.14159265358979,-0.841787214476933]];
        elseif theta==20
            points=[[0,0];[3.14159265358979,2.19976685802782];[3.94224231223721,1.05632064393176];[3.14159265358979,-1.14344621409606]];
        elseif theta==25
            points=[[0,0];[3.14159265358979,2.00141525117335];[4.07486791236014,0.536466538000949];[3.14159265358979,-1.46494871317240]];
        elseif theta==30
            points=[[0,0];[3.14159265358980,1.81379936423422];[4.18879020478639,0];[3.14159265358979,-1.81379936423422]];
        elseif mod(theta,1)==0 && ~(theta<0 || theta>30)
            points=elsepoints(:,:,theta);
        elseif mod(theta,1)==0.5 && ~(theta<0 || theta>30)
            points=elsepoints05(:,:,theta+0.5);
        end
        
    elseif strcmp(lattice,'dislocated')
        
        N=4;
        
        steps=1000;
        bound1=-0.7*pi/a;
        bound2=0.7*pi/a;
        stepsize=abs(bound1-bound2)/steps;
        
    elseif strcmp(lattice,'rect')
        N=4;
        point_names=["\Gamma","X","M","Y"];
        if theta==0
            points=[[0,0];[2.09439510239320,0];[2.09439510239320,3.14159265358979];[0,3.14159265358979]];
        end
        
    elseif strcmp(lattice,'graphene') || strcmp(lattice,'check')
        
        N=2;
        points=((2*pi)/(3*a))*[[1,1/sqrt(3)]; [1,-1/sqrt(3)]; [-1,1/sqrt(3)]; [-1,-1/sqrt(3)]; [0,1]; [0,-1]];
        
    elseif strcmp(lattice,'disloc7')
        N=7;
        
    else
        error('invalid lattice string')
    end
    
    D=length(point_names);
    
    %% -----------------------calculate band structure-------------------------
    % calculate and diagonalize the Hamiltonian on a meshgrid of k-values. its
    % eigenvalues form bands which can be plotted over all kx, ky to obtain the
    % band structure. from the Hamiltonian's eigenvectors the system's
    % eigenstates can be calulculated.
    
    [x,y]=meshgrid(bound1:stepsize:bound2-stepsize); % create meshgrid
    ew=zeros(N,steps,steps); % initialize matrix for eigenvalues
    
    for n=1:steps % iterate over meshgrid
        for m=1:steps
            kx=x(n,m); % set variables of Hamiltonian to current value
            ky=y(n,m);
            
            M = Hamiltonian(lattice, kx, ky, t, a, theta, nexp, delta, gamma);
            % function to calculate Hamiltonian at the end of this file.
            % Hamiltonian must be defined for each lattice.
            
            v=eig(M); % calculate eigenvalues
            for d=1:N
                ew(d,n,m)=v(d); % write eigenvalues into prepared matrix
            end
        end
    end
    
    %% -----------------------plot 1D band structure---------------------------
    % plots 1D bandstructure along high symmetry lines. Uses function
    % Bandstruktur1D from the end of this file and external function
    % subplot_tight. Written specifically for trans lattice. Must be adapted
    % for use with other lattices.
    
    
    if ~isempty(points)
        fig1=figure('Name','1D-Bandstruktur');
        fig1.Units='normalized';
        fig1.OuterPosition=[0 0 1 1]; % fullscreen
        ma=[0,0];
        lpos=0;
        plotcount=1;
        for d=1:D
            subplot_tight(1,D,d,ma); % function that allows seamlessly connected subplots
            if d==1
                lin1=Bandstruktur1D(points(d,:), points(mod(d,D)+1,:));
            elseif d==2
                lin2=Bandstruktur1D(points(d,:), points(mod(d,D)+1,:));
            elseif d==3
                lin3=Bandstruktur1D(points(d,:), points(mod(d,D)+1,:));
            elseif d==4
                lin4=Bandstruktur1D(points(d,:), points(mod(d,D)+1,:));
            else
                Bandstruktur1D(points(d,:), points(mod(d,D)+1,:));
            end
            
        end
        title("\theta = " + num2str(theta) + '°');
        %     saveas(fig1, [num2str(nexp) 'lin' num2str(theta) '.png']);
    end
    
    %% -----------------------plot 2D band structure---------------------------
    % visualize band structure as a surface plot.
    
    fig2=figure('Name','2D-Bandstruktur');
    fig2.Units='normalized';
    fig2.OuterPosition=[0 0 1 1];
    shading interp
    for d=1:N
        s=mesh(x,y,squeeze(ew(d,:,:)));
        hold on
        
        s.EdgeAlpha=1.0;
        s.FaceLighting='gouraud';
        s.FaceColor='interp';
        colormap jet;
    end
    set(gca,'xtick',[-2*pi:pi:2*pi]) % where to set the tick marks
    set(gca,'xticklabels',{'-2\pi','-\pi','0','\pi','2\pi'},'FontSize',14) % give them user-defined labels
    set(gca,'ytick',[-2*pi:pi:2*pi]) % where to set the tick marks
    set(gca,'yticklabels',{'-2\pi','-\pi','0','\pi','2\pi'},'FontSize',14) % give them user-defined labels
    view(83.862060314988838,2.092363018197072); %view(55,6); % viewing angle
    title("\theta = " + num2str(theta) + '°');
    xlabel("k_x/a")
    ylabel("k_y/a")
    zlabel("\beta_n [a.u.]")
    % saveas(fig2,[num2str(nexp),'s',num2str(theta),'t.png']);
    
    hold off
    clear CO
    
    %% --------------------------find Dirac points-----------------------------
    % finds DP by finding minima of the difference between bands. Uses
    % external functions extrema.m and extrema2.m from the matlab file exchange.
    % Also calculates EV at DP to use for lattice excitation.
    
    
    
    ic=1;
    for d=1:N-1
        diff(d,:,:)=abs(ew(d,:,:)-ew(d+1,:,:)); % calculate difference between bands
        [~,~,zmin,imin]=extrema2(squeeze(diff(d,:,:))); % minima of difference -> DP
        for i=1:length(imin) % save relevant values of DP in variable
            if zmin(i)/min(zmin) <= cutoff
                diracp(ic,1)=x(imin(i));
                diracp(ic,2)=y(imin(i));
                diracp(ic,3)=zmin(i);
                diracp(ic,4)=imin(i);
                diracp(ic,5)=d;
                %             trialvalues(ic)=ew(exband,round((y(imin(i))+trial(2)-bound1)/stepsize)+1,round((x(imin(i))+trial(1)-bound1)/stepsize)+1);
                
                kx=x(imin(i));
                ky=y(imin(i));
                
                if calculateEV % calculate EV at DP.
                    
                    M = Hamiltonian(lattice, kx, ky, t, a, theta, nexp, delta, gamma);
                    
                    [EV(:,:,ic),~]=eig(M);
                end
                
                ic=ic+1;
            end
        end
    end
    
    %% -------------------------plot Dirac points------------------------------
    % Highlights DP in plots of individual bands and contour plots. The top and
    % bottom band are omitted. DP on top of the bands are marked in blue, on
    % the bottom in red. Cones are numbered to help decide which to save into
    % files.
    
    for d=2:N-1
        
        figure('Name',['diracp-points, band ',num2str(d)]);
        fig=gcf;
        fig.Units='normalized';
        fig.OuterPosition=[0 0 1 1];
        
        s2=mesh(x,y,squeeze(ew(d,:,:)));
        s2.EdgeAlpha=1.0;
        s2.FaceLighting='gouraud';
        s2.FaceColor='interp';
        colormap jet;
        
        set(gca,'xtick',[-2*pi:pi:2*pi]) % where to set the tick marks
        set(gca,'xticklabels',{'-2\pi','-\pi','0','\pi','2\pi'},'FontSize',14) % give them user-defined labels
        set(gca,'ytick',[-2*pi:pi:2*pi]) % where to set the tick marks
        set(gca,'yticklabels',{'-2\pi','-\pi','0','\pi','2\pi'},'FontSize',14) % give them user-defined labels
        
        hold on
        for i=1:size(diracp,1) % i=1:size(diracp,1)
            if (diracp(i,5) == d || diracp(i,5) == d-1)
                str={['X: ',num2str(diracp(i,1))],['Y: ',num2str(diracp(i,2))],['cone ',num2str(i)]};
                zc=ew(d,round((y(diracp(i,4))-bound1)/stepsize)+1,round((x(diracp(i,4))-bound1)/stepsize)+1);
                %             if i>size(diracp,1)/2
                if diracp(i,5) == d
                    plot3(diracp(i,1),diracp(i,2),zc,'.b','markersize',20);
                else
                    plot3(diracp(i,1),diracp(i,2),zc,'.r','markersize',20);
                end
                text(diracp(i,1)+0.02*abs(bound1-bound2),diracp(i,2)+0.02*abs(bound1-bound2),zc,str);
            end
        end
        hold off
        
        figure('Name',['contour-plot, band ',num2str(d)]);
        fig=gcf;
        fig.Units='normalized';
        fig.OuterPosition=[0 0 1 1];
        contour(x,y,squeeze(ew(d,:,:)));
        pbaspect([1 1 1]);
        
        set(gca,'xtick',[-2*pi:pi:2*pi]) % where to set the tick marks
        set(gca,'xticklabels',{'-2\pi','-\pi','0','\pi','2\pi'},'FontSize',14) % give them user-defined labels
        set(gca,'ytick',[-2*pi:pi:2*pi]) % where to set the tick marks
        set(gca,'yticklabels',{'-2\pi','-\pi','0','\pi','2\pi'},'FontSize',14) % give them user-defined labels
        
        hold on
        for i=1:size(diracp,1)
            if (diracp(i,5) == d || diracp(i,5) == d-1)
                str={['X: ',num2str(diracp(i,1))],['Y: ',num2str(diracp(i,2))],['cone ',num2str(i)]};
                %             if i>size(diracp,1)/2
                if diracp(i,5) == d
                    plot3(diracp(i,1),diracp(i,2),diracp(i,3),'.b','markersize',20);
                else
                    plot3(diracp(i,1),diracp(i,2),diracp(i,3),'.r','markersize',20);
                end
                text(diracp(i,1)+0.1,diracp(i,2)+0.1,0,str);
            end
        end
        hold off
        
    end
    
    fn=3000+nexp; % number used to automatically differentiate files with Dirac points. Set initial filenumber for following sections. Can be changed from the command line if necessary.
    clear cn; % ensures the next 2 sectons aren't executing when running the full script.
    
    %% -----------------------plot dirac points with BZ------------------------
    % plots Dirac points along with the first Brilouin zone (outline and
    % high-symmetry points). Must be adapted for each lattice. Currently only
    % defined for lattice 'trans' and theta 0, 10, 15, 20, 30.
    
    if strcmp(lattice,'trans') && (theta == 0 || theta == 10 || theta == 20 || theta == 30 || theta == 15)
        
        %for d=1:N-1
        
        d = 2;
        
        load(['kp-' num2str(theta) '.mat']);
        load(['xp-' num2str(theta) '.mat']);
        load(['yp-' num2str(theta) '.mat']);
        if theta ~= 0
            load(['mp-' num2str(theta) '.mat']);
        end
        if strcmp(lattice, 'graphene')
            kp=diracp;
        end
        
        figure('Name', '1st Brillouin zone');
        fig=gcf;
        fig.Units='normalized';
        fig.OuterPosition=[0 0 1 1];
        contour(x,y,squeeze(ew(d,:,:)));
        pbaspect([1 1 1]);
        
        title("\theta = " + num2str(theta) + '°');
        xlabel("k_x/a")
        ylabel("k_y/a")
        set(get(gca,'YLabel'),'Rotation',0)
        
        set(gca,'xtick',[-2*pi:pi:2*pi]) % where to set the tick marks
        set(gca,'xticklabels',{'-2\pi','-\pi','0','\pi','2\pi'},'FontSize',20) % give them user-defined labels
        set(gca,'ytick',[-2*pi:pi:2*pi]) % where to set the tick marks
        set(gca,'yticklabels',{'-2\pi','-\pi','0','\pi','2\pi'},'FontSize',20) % give them user-defined labels
        
        hold on
        plot3(0,0,0,'.k','markersize',30); %20
        text(0.1,0.1,"\Gamma",'FontSize',20)
        for i=1:size(diracp,1)
            if (diracp(i,5) == d || diracp(i,5) == d-1)
                str={['X: ',num2str(diracp(i,1))],['Y: ',num2str(diracp(i,2))],['cone ',num2str(i)]};
                %             if i>size(diracp,1)/2
                if diracp(i,5) == d
                    plot3(diracp(i,1),diracp(i,2),diracp(i,3),'.b','markersize',40); %25
                else
                    if theta ~= 30 && theta ~= 0
                        plot3(diracp(i,1),diracp(i,2),diracp(i,3),'.r','markersize',40);
                    end
                end
                
                %             text(diracp(i,1)+0.1,diracp(i,2)+0.1,0,str);
            end
        end
        
        for i=1:size(xp,1)
            plot3(xp(i,1),xp(i,2),0,'.k','markersize',30);
            text(xp(i,1)+0.1,xp(i,2)+0.1,0,'X','FontSize',20);
        end
        
        for i=1:size(yp,1)
            plot3(yp(i,1),yp(i,2),0,'.k','markersize',30);
            text(yp(i,1)+0.1,yp(i,2)+0.1,0,'Y','FontSize',20);
        end
        
        if theta ~= 0
            for i=1:size(mp,1)
                plot3(mp(i,1),mp(i,2),0,'.k','markersize',30);
                text(mp(i,1)+0.1,mp(i,2)+0.1,0,'M','FontSize',20);
            end
        end
        
        for i=1:size(kp,1)
            
            if theta == 0
                text(kp(i,1)+0.1,kp(i,2)+0.1,0,'M','FontSize',20)
            else
                if mod(i,2)==1
                    text(kp(i,1)+0.1,kp(i,2)+0.1,0,"K",'FontSize',20)
                else
                    text(kp(i,1)+0.1,kp(i,2)+0.1,0,"K'",'FontSize',20)
                end
            end
            
            px(i)=kp(i,1);
            py(i)=kp(i,2);
            pz(i)=0;
            
            if theta ~= 0 %&& theta ~= 30
                plot3(kp(i,1),kp(i,2),0,'.k','markersize',30); %20
            end
        end
        line(px,py,pz,'Color','black','LineStyle','--');
        line([kp(size(kp,1),1) kp(1,1)], [kp(size(kp,1),2) kp(1,2)], [0 0],'Color','black','LineStyle','--')
        
        hold off
        
    end
    
    fn=2000+nexp; % set initial filenumber for following sections. Can be changed from the command line if necessary.
    clear cn; % ensures the next 2 sectons aren't executing when running the full script.
    
    %end
    
    %% ------------------------save DP and EV to file--------------------------
    % Execute this section manually to save DP and EV to a file. This and the
    % next section are not run when running the full script, but have to be run
    % individually. Cone numbers can be seen in the previous section's plots and
    % have to be entered into the variable "cn". Files are saved under the name
    % *theta*cones*fn*, e.g. 25cones2003. If the file already exists the
    % section will not be executed.
    
    
    if exist('cn','var')
        
        cn=[12 15 9 13 16 10 14 11]; % enter cones to be saved here
        cones=zeros(length(cn),5);
        
        if (exist([folder,num2str(theta),'cones',num2str(fn),'.mat'],'file') || exist([folder,num2str(theta),'cones',num2str(fn),'.png'],'file') || exist([folder,num2str(theta),'cones',num2str(fn),'.fig'],'file') || exist([folder,num2str(theta),'ev',num2str(fn),'.mat'],'file'))
            warning('file already exists, please delete it if you want to overwrite it, or change fn otherwise')
        else
            if ~isempty(cn)
                for i=1:length(cn)
                    cones(i,:)=diracp(cn(i),:);
                    ev(:,:,i)=EV(:,:,cn(i));
                end
                save([folder,num2str(theta),'cones',num2str(fn),'.mat'],'cones'); % save k-space coordinates of DP
                if calculateEV
                    save([folder,num2str(theta),'ev',num2str(fn),'.mat'],'ev'); % save EV at DP
                end
                
                if cones(i,5)==1
                    d=2;
                else
                    d=cones(i,5);
                end
                
                figure('Name',[num2str(theta),'cones',num2str(fn)]); % contour plot with DP to show which cones are saved
                fig=gcf;
                fig.Units='normalized';
                fig.OuterPosition=[0 0 1 1];
                contour(x,y,squeeze(ew(d,:,:)));
                pbaspect([1 1 1]);
                
                hold on
                for i=1:size(cones,1)
                    str={['X: ',num2str(cones(i,1))],['Y: ',num2str(cones(i,2))],['cone ',num2str(i)]};
                    if i>size(cones,1)/2 % first half of DP are marked in blue, rest in red, to facilitate exciting specific DP in simulation
                        plot3(cones(i,1),cones(i,2),0,'.b','markersize',20);
                    else
                        plot3(cones(i,1),cones(i,2),0,'.r','markersize',20);
                    end
                    text(cones(i,1)+0.1,cones(i,2)+0.1,0,str);
                    
                end
                hold off
                saveas(gcf,[folder,num2str(theta),'cones',num2str(fn),'.png']) % save as .png as a reminder which cones were saved
                %             saveas(gcf,[folder,num2str(theta),'cones',num2str(fn),'.fig'])
                
                fn=fn+1000 % automatical increment of filenumber after running the section. Can be changed from the command line if necessary.
            end
        end
    end
    
    %% ------------------------auto-save DP and EV-----------------------------
    % auto-saves DP to file when executing full script. Only executed if
    % "autosaveDP" is set to true in input prameters. Does not overwrite files.
    
    if autosaveDP
        if exist([folder,num2str(theta),'cones',num2str(fn),'.mat'],'file')
            warning([num2str(theta),'ev',num2str(fn),'.mat',' already exists, please delete it if you want to overwrite it, or change fn otherwise'])
        else
            
            ic=1;
            
            for i=1:size(diracp,1)
                diracp(i,6)=sqrt(diracp(i,1)^2 + diracp(i,2)^2);
                diracp(i,7)=i;
                if diracp(i,5)==2
                    autoDP(ic,1)=diracp(i,1);
                    autoDP(ic,2)=diracp(i,2);
                    autoDP(ic,3)=diracp(i,6);
                    autoDP(ic,4)=diracp(i,7);
                    ic=ic+1;
                end
            end
            
            autoDP=sortrows(autoDP,3);
            
            for i=1:6
                cones(i,:)=autoDP(i,:);
                ev(:,:,i)=EV(:,:,autoDP(i,4));
            end
            
            figure('Name',[num2str(theta),'cones',num2str(fn)]);
            fig=gcf;
            fig.Units='normalized';
            fig.OuterPosition=[0 0 1 1];
            contour(x,y,squeeze(ew(d,:,:)));
            pbaspect([1 1 1]);
            
            hold on
            for i=1:size(cones,1)
                str={['X: ',num2str(cones(i,1))],['Y: ',num2str(cones(i,2))],['cone ',num2str(i)]};
                plot3(cones(i,1),cones(i,2),0,'.b','markersize',20);
                text(cones(i,1)+0.1,cones(i,2)+0.1,0,str);
            end
            hold off
            
            saveas(gcf,[folder,num2str(theta),'cones',num2str(fn),'.png'])
            %         saveas(gcf,[folder,num2str(theta),'cones',num2str(fn),'.fig'])
            save([folder,num2str(theta),'cones',num2str(fn),'.mat'],'cones');
            if calculateEV
                save([folder,num2str(theta),'ev',num2str(fn),'.mat'],'ev');
            end            
        end
    end
    
end

%% -------------------calculate ribbon band structure----------------------
% calculates ribbon band structure for the selected lattice. Further
% information in the corresponding function.

if calculate_ribbon && valid_ribbon % check if section shall/can be executed
    
    [matrix_size, ~] = Ribbon(ribbon, ribbon_size, t, a, theta, nexp, delta, gamma0, 0);
    
    kxr = zeros(1, steps_r+1);
    ewr = zeros(steps_r+1, matrix_size);
    evr = zeros(steps_r+1, matrix_size, matrix_size);
    
    m=1;
    for k=bound1_r:stepsize_r:bound2_r
        
        [~, M] = Ribbon(lattice, ribbon_size, t, a, theta, nexp, delta, gamma0, k);
        % function to calculate ribbon Hamiltonian at the end of this file.
        % Hamiltonian must be defined for each lattice.
        
        kxr(m)=k;
        [EV,EW]=eig(M);
        ewr(m,:)=diag(EW);
        evr(m,:,:)=EV;
        m=m+1;
        
    end
    
    figure('Name', 'ribbon band structure');
    title(append("\theta = ", num2str(theta), "°"))
    % fullscreen();
    hold on
    box on
    
    fig=gcf;
    fig.Units='normalized';
    fig.OuterPosition=[0 0 1 1];
    pbaspect([1 1 1]);
    
    xlabel("k_xa")
    ylabel("\beta [a.u.]")
    set(gca,'xtick',[0:pi/2:2*pi])
    % set(gca,'xticklabels',{'0','\pi/2','\pi','3\pi/2','2\pi'})
    set(gca,'xticklabels',{'0','\pi/4','\pi/2','3\pi/4','\pi'})
    set(gca,'xLim',[0,2*pi])
    
    for n=1:size(ewr,2)
        plot(kxr,ewr(:,n),'k');
    end
    
end

%% ---------------calculate Floquet ribbon band structure------------------
% calculates Floquet ribbon band structure for the selected lattice. 
% Further information in the corresponding function.

if calculate_floquet && valid_floquet % check if section shall/can be executed
    
    [matrix_size, ~] = Floquet(floquet, ribbon_size, a_f, theta, nexp, k0, R, omega, omega2, 0);
    
    xf=linspace(bound1_f,bound2_f,steps_f); % k-values in chosen range
    l=zeros(matrix_size,steps_f); % matrix for eigenvalues
    
    for m=1:steps_f
        
        k=xf(m);
        
        [~, M] = Floquet(lattice, ribbon_size, a_f, theta, nexp, k0, R, omega, omega2, k);
        % function to calculate Floquet ribbon Hamiltonian at the end of this file.
        % Hamiltonian must be defined for each lattice.
        
        U=eye(matrix_size); % initialize time-evolution operator
        
        for z=zstep:zstep:Z % calculate time-evolution operator. feval calculates M at time z
            U=expm(-1i*feval(M,z)*zstep) * U;
        end
        %         U
        [~,Ueig]=eig(U); % calculate eigenvectors (C) and eigenvalues (Ueig)
        %         ExpectM_0 = C' * M0 * C; % sorting of eigenvalues. Doesn't work with nnn-int.act.
        %         [ExpSort,Index]=sort(diag(ExpectM_0));
        ewf=angle(diag(Ueig))/Z; % calculate quasienergies from eigenvalues
        [ewf, sorti]=sort(ewf);
        %         C=C(:,sorti);
        l(:,m)=ewf; % put quasienergies in prepared matrix
        %         l(:,m)=l(Index(:),m); % sorting of eigenvalues. Doesn't work with nnn-int.act.
        %         Dk1(m,:,:)=C(:,:); % eigenvectors
        
        if mod(m,100)==0 % indicates progress
            disp(['Floquet calculation: progress = ', num2str((m/steps_f)*100), '%'])
        end
        
    end
    
    fig=figure('Name', 'Floquet ribbon band structure');
    pbaspect([1 1 1]);
    fig.Units='normalized';
    fig.OuterPosition=[0 0 1 1];
    
    xf=xf*2e-5;
    l=l/100;
    
%     for li=1:size(l,1)
%         for lj=1:size(l,2)
%             if l(li,lj) < - 280
%                 l(li,lj)=l(li,lj) + 100*2*pi;
%             end
%         end
%     end
    
    for i=1:matrix_size
        
        plot(xf,l(i,:),'.','Color','black','markersize',2); % bands should be plotted as dots, because sorting of ew does not work with nnn
        
        xlabel("k_xa")
        ylabel("\beta [a.u.]")
        %     set(gca,'xtick',[0:pi/2:2*pi])
        %     set(gca,'xticklabels',{'0','\pi/2','\pi','3\pi/2','2\pi'})
        % %     set(gca,'xticklabels',{'0','\pi/4','\pi/2','3\pi/4','\pi'})
        %     set(gca,'xLim',[0,2*pi])
        
        hold on
    end
    pbaspect([1 1 1]);
    
    title(['R=' num2str(R) ' m   ' 'Z=' num2str(Z) ' m   ' 'theta=' num2str(theta) '°   ' 'n=' num2str(nexp) ''])
%     title(append("\theta = ", num2str(theta), "°"))
%     title(append('R = ', num2str(R/a_f), " \cdot a"))
    % title(append('R = ', num2str(R/a), 'a'))
    hold off
    
    % saveas(gcf,['0F' num2str(cvar) '.png']);
    
end

    %% ------------------------------------------------------------------------
  
    cn=[]; % enables running section "save DP to file"
    disp(['total elapsed time t=' num2str(toc) 'sec .']) % displays total runtime

%% function Valid
% checks if the chosen lattice string is valid for a 2D-Hamiltonian, ribbon
% and/or Floquet-ribbon. Each new lattice added to any of the other
% functions has to be added here too.

function[valid_lattice, valid_ribbon, valid_floquet] = Valid(lattice, ribbon, floquet)

%check 2D Hamiltonian
if strcmp(lattice, 'transition') || strcmp(lattice, 'transition2') || strcmp(lattice, 'transition3') || strcmp(lattice, 'transition4') || strcmp(lattice, 'trans') || strcmp(lattice, 'trans_st') || strcmp(lattice, 'trans5') || strcmp(lattice, 'dislocated') || strcmp(lattice, 'rect') || strcmp(lattice, 'graphene') || strcmp(lattice, 'check') || strcmp(lattice, 'disloc7') || strcmp(lattice, 'trans_qr') || strcmp(lattice, 'trans_q') || strcmp(lattice, 'trans_i') || strcmp(lattice, 'sq1') || strcmp(lattice, 'sq2')
    valid_lattice = true;
else
    valid_lattice = false;
    warning('2D Hamiltonian undefined')
end

%check ribbon Hamiltonian
if strcmp(ribbon, 'lieb_y_bearded') || strcmp(ribbon, 'trans2_flat') || strcmp(ribbon, 'kagome3_mixed') || strcmp(ribbon, 'trans') || strcmp(ribbon, 'trans4_st_mixed') || strcmp(ribbon, 'trans3_mixed') || strcmp(ribbon, 'kagome2_flat') || strcmp(ribbon, 'trans2_mixed') || strcmp(ribbon, 'trans2_st_bearded') || strcmp(ribbon, 'trans_st_bearded') || strcmp(ribbon, 'graphene_zigzag2') || strcmp(ribbon, 'graphene_zigzag3') || strcmp(ribbon, 'lieb_st_y_bearded') || strcmp(ribbon, 'graphene_y_zigzag') || strcmp(ribbon, 'trans_y_bearded') || strcmp(ribbon, 'kagome_x') || strcmp(ribbon, 'kagome2_x') || strcmp(ribbon, 'kagome3_x') || strcmp(ribbon, 'trans2')
    valid_ribbon = true;
else
    valid_ribbon = false;
    warning('ribbon Hamiltonian undefined')
end

%check Floquet ribbon Hamiltonian
if strcmp(floquet, 'trans') || strcmp(floquet, 'graphene_y_zigzag') || strcmp(floquet, 'trans4_mixed')
    valid_floquet = true;
else
    valid_floquet = false;
    warning('Floquet ribbon Hamiltonian undefined')
end

end

%% function Bandstruktur1D
%
% Calculates the Bandstructure along lines between points in k-space.
% Calculates line between the points and takes corresponding k-values from
% meshgrids and eigenvalues from 2D-bandstructure.
%
% algorithm by Jan Wichmann.

function[ewl]=Bandstruktur1D(startp, endp)

global bound1 bound2 stepsize steps ew lpos D N point_names points plotcount;

pointdiffs=zeros(1,D); % calculate the distance between Dirac points
for d = (1:D)
    pointdiffs(d)=norm(points(mod(d,D)+1,:)-points(d,:));
end
for d = 1:D
    xasp(d)=(0.4*pointdiffs(d))/max(pointdiffs); % set aspect ratio of subplots depending on distance between Dirac points
end

x1=((startp(1)-bound1)/stepsize);
y1=((startp(2)-bound1)/stepsize);
x2=((endp(1)-bound1)/stepsize);
y2=((endp(2)-bound1)/stepsize);

% calculate line through matrix ew between the points given as arguments

slope=(y2-y1)/(x2-x1); % slope of the line  thorugh the matrix, defined by the two points

if abs(slope)<1 && x2-x1>0 || x2-x1<0 % case slope < 1: x is independent variable, y = x*slope. case slope > 1: y is independent variable, x = y*slope^-1.
    slope=(y2-y1)/(x2-x1);
    
    l=zeros(1,abs(round(x2-x1))+1);
    k=zeros(1,abs(round(x2-x1))+1);
    ewl=zeros(N,abs(round(x2-x1))+1);
    
    for j=1:abs(round(x2-x1))+1
        
        l(j)=sign(x2-x1)*(j-1)+round(x1); % increase/decrease the independent variable in steps of 1
        k(j)=round(((l(j)-l(1)))*slope+y1); % linear equation
        
        for d=1:N
            ewl(d,j)=ew(d,k(j),l(j));
        end
    end
else
    slope=slope^-1; % case slope > 1. use inverse of slope and y as independent variable
    
    l=zeros(1,abs(round(y2-y1))+1);
    k=zeros(1,abs(round(y2-y1))+1);
    ewl=zeros(N,abs(round(y2-y1))+1);
    
    for j=1:abs(round(y2-y1))+1
        
        k(j)=sign(y2-y1)*(j-1)+round(y1); % linear equation
        l(j)=round(((k(j)-k(1)))*slope+x1); % increase/decrease the independent variable in steps of 1
        
        for d=1:N
            ewl(d,j)=ew(d,k(j),l(j));
        end
    end
end

% ---------------------------------plot------------------------------------
%
% plot resulting energy bands in 1D

% xasp=(0.4*pointdiffs(plotcount))/max(pointdiffs);
yasp=1; % =sum(xasp);
xxasp=xasp(plotcount)/sum(xasp);
hold on
for d=1:N
    testplot=plot(squeeze(ewl(d,:)));
    uistack(testplot,'top')
end
hold off
pbaspect([xxasp yasp 1]);
ax=gca;
ax.TickDir='out';
% axis([1 j -3 4]);
axis([1 j -3 4.5]);
ax.FontSize=32;
set(ax,'box','on');
ylabel("\beta [a.u.]")

if plotcount==1
    ax.XTick=[1,j];
    ax.XTickLabel={point_names(1),point_names(2)};
elseif plotcount==D
    ax.XTick=[j];
    ax.XTickLabel={point_names(1)};
else
    ax.XTick=[j];
    ax.XTickLabel={point_names(plotcount+1)};
end

if plotcount~=1
    ax.YTick=[];
    ax.YAxisLocation = 'right';
    ylabel([])
end

% if plotcount==D
%     ax.YAxisLocation = 'right';
% elseif plotcount~=1
%     ax.YTick=[];
% end

ax.Position=[lpos+0.1,0,xxasp*0.3125,1];
lpos=lpos+xxasp*0.3125;
plotcount=plotcount+1;

end

%% function Hamiltonian
% calculates tight-binding Hamiltonian. Hamiltonian must be defined for
% each individual lattice

function[M] = Hamiltonian(lattice, kx, ky, t, a, theta, nexp, delta, gamma)

gamma1=(exp(1-sqrt(2-2*sind(theta))))^nexp;  % t_nnn/t_nn for transition lattices
gamma2=(exp(1-sqrt(2+2*sind(theta))))^nexp;

gamma_AB_p=exp(-2*nexp*delta); % additional gammas for trans_st
gamma_AB_m=exp(2*nexp*delta);
gamma_AC_1_p=exp(nexp*(1-(1+2*delta)*sqrt(2-2*sind(theta))));
gamma_AC_1_m=exp(nexp*(1-(1-2*delta)*sqrt(2-2*sind(theta))));
a_AC_2_p=norm([0.5*(sind(theta)+1)+delta*(sind(theta)-1) cosd(theta)*(0.5+delta)]);
a_AC_2_m=norm((0.5+delta)*[1 0]+(0.5-delta)*[sind(theta) cosd(theta)]);
gamma_AC_2_p=exp(nexp*(1-(a_AC_2_p*2)));
gamma_AC_2_m=exp(nexp*(1-(a_AC_2_m*2)));

if strcmp(lattice,'transition') % insert Hamiltonian here, variables: kx, ky
    M=2*t*[0,cos(kx/1),cos((kx*sind(theta)+ky*cosd(theta))/1);cos(kx/1),0,gamma1*cos(kx*(1-sind(theta))/1 - ky*cosd(theta)/1) + gamma2*cos(kx*(1-sind(theta))/1 + ky*cosd(theta)/1);cos((kx*sind(theta)+ky*cosd(theta))/1),gamma1*cos(kx*(1-sind(theta))/1 - ky*cosd(theta)/1) + gamma2*cos(kx*(1-sind(theta))/1 + ky*cosd(theta)/1),0];
elseif strcmp(lattice,'transition2')
    M=2*t*[0,cos(kx/2),gamma1*cos(kx*(1-sind(theta))/2 - ky*cosd(theta)/2) + gamma2*cos(kx*(1+sind(theta))/2 + ky*cosd(theta)/2);cos(kx/2),0,cos((kx*sind(theta)+ky*cosd(theta))/2);gamma1*cos(kx*(1-sind(theta))/2 - ky*cosd(theta)/2) + gamma2*cos(kx*(1+sind(theta))/2 + ky*cosd(theta)/2),cos((kx*sind(theta)+ky*cosd(theta))/2),0];
elseif strcmp(lattice,'transition3')
    M=2*t*[0,(1+exp(-1i*(kx*sind(theta)+ky*cosd(theta)))),gamma1*(1+exp(-1i*(kx*(1-sind(theta)) - ky*cosd(theta)))) + gamma2*(1+exp(-1i*(kx*(1-sind(theta))/2 + ky*cosd(theta)/2)));(1+exp(+1i*(kx*sind(theta)+ky*cosd(theta)))),0,(1+exp(-1i*kx));gamma1*(1+exp(+1i*(kx*(1-sind(theta)) - ky*cosd(theta)))) + gamma2*(1+exp(+1i*(kx*(1-sind(theta))/2 + ky*cosd(theta)/2))),(1+exp(+1i*kx)),0];
elseif strcmp(lattice,'transition4')
    AB=cos((kx*sind(theta)+ky*cosd(theta))/2);
    BC=cos(kx/2);
    AC=gamma1*cos(kx*(1-sind(theta))/2 - ky*cosd(theta)/2) + gamma2*cos(kx*(1-sind(theta))/2 + ky*cosd(theta)/2);
    M=2*t*[0,BC,AB;conj(BC),0,AC;conj(AB),conj(AC),0];
elseif strcmp(lattice,'trans')
   
    AB=cos((kx*sind(theta)+ky*cosd(theta))/2);
    BC=cos(kx/2);
    AC=gamma1*cos(kx*(1-sind(theta))/2 - ky*cosd(theta)/2) + gamma2*cos(kx*(1+sind(theta))/2 + ky*cosd(theta)/2);
    
    M=2*t*[0,AB,AC;conj(AB),0,BC;conj(AC),conj(BC),0];
    
elseif strcmp(lattice,'trans_st')
    
    AB=gamma_AB_m*exp(-1i*(0.5-delta)*(kx*sind(theta)+ky*cosd(theta)))+gamma_AB_p*exp(1i*(0.5+delta)*(kx*sind(theta)+ky*cosd(theta)));
    BC=gamma_AB_m*exp(1i*(0.5-delta)*kx)+gamma_AB_p*exp(-1i*(0.5+delta)*kx);
    AC1=gamma_AC_1_m*exp(1i*(0.5-delta)*(kx*(1-sind(theta)) - ky*cosd(theta)))+gamma_AC_1_p*exp(-1i*(0.5+delta)*(kx*(1-sind(theta)) - ky*cosd(theta)));
    AC2=gamma_AC_2_p*exp(1i*(kx*(0.5*(sind(theta)+1)+delta*(sind(theta)-1))+ky*(cosd(theta)*(0.5+delta))))+gamma_AC_2_m*exp(-1i*(kx*((0.5+delta)+(0.5-delta)*sind(theta))+ky*((0.5-delta)*cosd(theta))));
    AC=AC1+AC2;
  
    M=2*t*[0,BC,AB;conj(BC),0,AC;conj(AB),conj(AC),0];
    
elseif strcmp(lattice,'trans5')
    
    Hy=exp((1i/3)*(kx*sind(theta)+ky*cosd(theta)));
    Hx=exp(1i*kx/3);
    Hv=gamma1*exp((1i/3)*(kx*(sind(theta)-1)+ky*cosd(theta)));
    Hw=gamma2*exp((1i/3)*(kx*(sind(theta)+1)+ky*cosd(theta)));
    
    M=t*[0 conj(Hy) Hy Hw Hv;
        Hy 0 conj(Hy) conj(Hv) conj(Hw);
        conj(Hy) Hy 0 Hx conj(Hx);
        conj(Hw) Hv conj(Hx) 0 Hx;
        conj(Hv) Hw Hx conj(Hx) 0];
    
elseif strcmp(lattice,'dislocated')
    M=[0,exp(a*1i*kx),0,2*cos(a*ky);exp(-a*1i*kx),0,exp(a*1i*kx),0;0,exp(-a*1i*kx),0,exp(a*1i*kx);2*cos(a*ky),0,exp(-a*1i*kx),0];
    
elseif strcmp(lattice,'rect')
    AB=2*cos((kx*sind(theta)+ky*cosd(theta))/2);
    AC=gamma1*exp((1i/2)*((1-sind(theta))*kx-cosd(theta)*ky)) + gamma2*exp((1i/2)*((1+sind(theta))*kx+cosd(theta)*ky));
    AD=gamma1*exp((-1i/2)*((1-sind(theta))*kx-cosd(theta)*ky)) + gamma2*exp((-1i/2)*((1+sind(theta))*kx+cosd(theta)*ky));
    BC=exp((1i/2)*kx);
    BD=exp((-1i/2)*kx);
    CD=exp((1i/2)*kx);
    M=t*[0 AB AC AD;
        conj(AB) 0 BC BD;
        conj(AC) conj(BC) 0 CD;
        conj(AD) conj(BD) conj(CD) 0];
    
elseif strcmp(lattice,'graphene')
    
    AB=t*exp(1i*ky*a/(2*sqrt(3)))*(exp(1i*sqrt(3)*ky*a/2)+2*cos(kx*a/2));
    M=[0,AB;conj(AB),0];
    
elseif strcmp(lattice,'check')
    
    AB=2*t*(cos(sind(45)*(kx+ky))+cos(sind(45)*(kx-ky)));
    M=[0,AB;conj(AB),0];
    
elseif strcmp(lattice,'disloc7')
    
    AB=exp(-1i*a*(kx*sind(theta)+ky*cosd(theta))/4);
    AC=0;
    AD=exp(1i*a*(kx*sind(theta)+ky*cosd(theta))/4);
    AE=gamma2*exp((1i/4)*(kx*(1+sind(theta)) + ky*cosd(theta)));
    AF=0;
    AG=gamma2*exp((-1i/4)*(kx*(1+sind(theta)) + ky*cosd(theta)));
    BC=AB;
    BD=0;
    BE=0;
    BF=0;
    BG=exp(-1i*a*kx/4);
    CD=AB;
    CE=gamma1*exp((1i/4)*(kx*(1-sind(theta)) - ky*cosd(theta)));
    CF=0;
    CG=gamma1*exp((-1i/4)*(kx*(1-sind(theta)) - ky*cosd(theta)));
    DE=exp(1i*a*kx/4);
    DF=0;
    DG=0;
    EF=DE;
    EG=0;
    FG=DE;
    
    M=[0 AB AC AD AE AF AG;
        conj(AB) 0 BC BD BE BF BG;
        conj(AC) conj(BC) 0 CD CE CF CG;
        conj(AD) conj(BD) conj(CD) 0 DE DF DG;
        conj(AE) conj(BE) conj(CE) conj(DE) 0 EF EG;
        conj(AF) conj(BF) conj(CF) conj(DF) conj(EF) 0 FG;
        conj(AG) conj(BG) conj(CG) conj(DG) conj(EG) conj(FG) 0;];
    
elseif strcmp(lattice,'trans_qr')
    
    AB=cos(kx/2);
    BC=cos(ky/2);
    AC=gamma*cos((kx+ky)/2);
    M=2*t*[0,AB,AC;conj(AB),0,BC;conj(AC),conj(BC),0];
    
elseif strcmp(lattice,'trans_q')
    
    AB=1+exp(1i*kx);
    BC=1+exp(1i*ky);
    AC=gamma*(1+exp(1i*(kx+ky)));
    M=t*[0,AB,AC;conj(AB),0,BC;conj(AC),conj(BC),0];
    
elseif strcmp(lattice,'trans_i')
    
    AB=1+exp(0.5*1i*(sind(theta)*kx+cosd(theta)*ky));
    BC=1+exp(0.5*1i*kx);
    AC=gamma1*(1+exp(1i*kx*(1-sind(theta))/2 - 1i*ky*cosd(theta)/2)) + gamma2*(1+exp(1i*kx*(1+sind(theta))/2 + 1i*ky*cosd(theta)/2));
    M=t*[0,AB,AC;conj(AB),0,BC;conj(AC),conj(BC),0];
    
elseif strcmp(lattice,'sq1')
    
    M=2*t*(cos(kx*a)+cos(ky*a));
    
elseif strcmp(lattice,'sq2')
    
    M=2*t*[0,(cos(kx*a)+cos(ky*a));(cos(kx*a)+cos(ky*a)),0];
    
else
    
    error('invalid lattice');
    
end

end

%% function Ribbon
% a unit cell is defined in the lattice which is repeated ribbon_size times
% in the finite direction. The definition of the unit cell decides the
% edges of the lattice. The elements in the unit cell are numbered. The
% matrix element M(n,m), for every n in the unit cell, is then:
% -   1 if there is a connection between n and m in the same unit cell,
% -   exp(1i*a*k) if there is a connection between n and m in the next unit cell,
% -   exp(-1i*a*k) if there is a connection between n and m in the previous unit cell,
% where "same, next, previous" unit cell  refers to the infinite direction.
% This is repeated ribbon_size times to fill the whole matrix of dimension
% unitcell_size*ribbon_size.

function[matrix_size, M] = Ribbon(ribbon, ribbon_size, t, a, theta, nexp, delta, gamma0, k)

gamma1=(exp(1-sqrt(2-2*sind(theta))))^nexp;  % t_nnn/t_nn for transition lattices
gamma2=(exp(1-sqrt(2+2*sind(theta))))^nexp;

gamma_p=exp(-2*nexp*delta); % for staggered trans lattices
gamma_m=exp(2*nexp*delta);
gamma_AC1_p=exp(nexp*(1-(1+2*delta)*sqrt(2-2*sind(theta))));
gamma_AC1_m=exp(nexp*(1-(1-2*delta)*sqrt(2-2*sind(theta))));
a_AC2=norm([0.5*(sind(theta)+1)+delta*(sind(theta)-1) cosd(theta)*(0.5+delta)]);
gamma_AC2=exp(nexp*(1-(a_AC2*2)));
a_AC3=norm(1.5*[sind(theta)+(1/3)-(2/3)*delta, cosd(theta)]);
gamma_AC3=0*exp(nexp*(1-(a_AC3*2)));
% a_AC_2_p=norm([0.5*(sind(theta)+1)+delta*(sind(theta)-1) cosd(theta)*(0.5+delta)]);
% a_AC_2_m=norm((0.5+delta)*[1 0]+(0.5-delta)*[sind(theta) cosd(theta)]);
% gamma_AC2=exp(nexp*(1-(a_AC_2_p*2)));
% gamma_AC2=exp(nexp*(1-(a_AC_2_m*2)));

if strcmp(ribbon, 'lieb_y_bearded')
    
    matrix_size=ribbon_size*3+1;
    M=zeros(matrix_size);
    
    for n=1:3:matrix_size
        
        M(n,n+1)=1;
        M(n+1,n+2)=exp(1i*k/2);
        M(n,n+2)=gamma0*(1+exp(1i*k/2));
        
        M(n+1,n)=1;
        M(n+2,n+1)=exp(-1i*k/2);
        M(n+2,n)=gamma0*(1+exp(-1i*k/2));
        
        M(n+1,n+3)=1;
        M(n+2,n+3)=gamma0*(1+exp(1i*k/2));
        
        M(n+3,n+1)=1;
        M(n+3,n+2)=gamma0*(1+exp(-1i*k/2));
        
    end
    
    matrix_size=size(M,1);
    
elseif strcmp(ribbon, 'trans2_flat')
    
    unitcell_size=3;
    matrix_size=ribbon_size*unitcell_size;
    M=zeros(matrix_size);
    
    for n=1:unitcell_size:matrix_size
        
        M(n,n+1)=1+exp(1i*a*k);
        if n<=matrix_size-unitcell_size
            M(n,n+2)=gamma1+gamma2*exp(1i*a*k);
            M(n+1,n+2)=1;
            M(n+2,n+3)=1;
            M(n+2,n+4)=gamma2+gamma1*exp(-1i*a*k);
        end
        
        M(n+1,n)=1+exp(-1i*a*k);
        if n<=matrix_size-unitcell_size
            M(n+2,n)=gamma1+gamma2*exp(-1i*a*k);
            M(n+2,n+1)=1;
            M(n+3,n+2)=1;
            M(n+4,n+2)=gamma2+gamma1*exp(1i*a*k);
        end
        
    end
    
    matrix_size=size(M,1);
    
elseif strcmp(ribbon, 'kagome3_mixed')
    
    unitcell_size=6;
    matrix_size=ribbon_size*unitcell_size;
    M=zeros(matrix_size);
    
    for n=1:unitcell_size:matrix_size
        
        M(n,n+1)=1+exp(-1i*a*k);
        M(n,n+2)=exp(-1i*a*k);
        M(n+1,n+2)=1;
        M(n+2,n+3)=1;
        M(n+2,n+4)=exp(1i*a*k);
        M(n+3,n+4)=1+exp(1i*a*k);
        M(n+3,n+5)=1;
        M(n+4,n+5)=1;
        if n<matrix_size-unitcell_size
            M(n+5,n+6)=1;
            M(n+5,n+7)=1;
        end
        
        M(n+1,n)=1+exp(1i*a*k);
        M(n+2,n)=exp(1i*a*k);
        M(n+2,n+1)=1;
        M(n+3,n+2)=1;
        M(n+4,n+2)=exp(-1i*a*k);
        M(n+4,n+3)=1+exp(-1i*a*k);
        M(n+5,n+3)=1;
        M(n+5,n+4)=1;
        if n<matrix_size-unitcell_size
            M(n+6,n+5)=1;
            M(n+7,n+5)=1;
        end
        
    end
    
    matrix_size=size(M,1);
    
elseif strcmp(ribbon, 'trans') % works best. flat edge at the bottom, zig-zag edge at the top.
    
    unitcell_size=6;
    matrix_size=ribbon_size*unitcell_size;
    M=zeros(matrix_size);
    
    for n=1:unitcell_size:matrix_size
        
        M(n,n+1)=1+exp(-1i*a*k);
        M(n,n+2)=exp(-1i*a*k);
        M(n+1,n+2)=gamma1+gamma2*exp(-1i*a*k);
        M(n+2,n+3)=1;
        M(n+2,n+4)=gamma2+gamma1*exp(1i*a*k);
        M(n+3,n+4)=1+exp(1i*a*k);
        M(n+3,n+5)=1;
        M(n+4,n+5)=gamma1+gamma2*exp(-1i*a*k);
        if n<matrix_size-unitcell_size
            M(n+5,n+6)=1;
            M(n+5,n+7)=gamma1+gamma2*exp(-1i*a*k);
        end
        
        M(n+1,n)=1+exp(1i*a*k);
        M(n+2,n)=exp(1i*a*k);
        M(n+2,n+1)=gamma1+gamma2*exp(1i*a*k);
        M(n+3,n+2)=1;
        M(n+4,n+2)=gamma2+gamma1*exp(-1i*a*k);
        M(n+4,n+3)=1+exp(-1i*a*k);
        M(n+5,n+3)=1;
        M(n+5,n+4)=gamma1+gamma2*exp(1i*a*k);
        if n<matrix_size-unitcell_size
            M(n+6,n+5)=1;
            M(n+7,n+5)=gamma1+gamma2*exp(1i*a*k);
        end
        
    end
    
    matrix_size=size(M,1);
    
elseif strcmp(ribbon, 'trans4_st_mixed')
    
    unitcell_size=6;
    matrix_size=ribbon_size*unitcell_size;
    M=zeros(matrix_size);
    
    for n=1:unitcell_size:matrix_size
        
        M(n,n+1)=gamma_p+gamma_m*exp(-1i*a*k);
        M(n,n+2)=gamma_m*exp(-1i*a*k)  +  gamma_AC3;
        M(n+1,n+2)=gamma_AC1_m+gamma_AC2*exp(-1i*a*k);
        M(n+2,n+3)=gamma_p;
        M(n+2,n+4)=gamma_AC2+gamma_AC1_p*exp(1i*a*k);
        M(n+3,n+4)=gamma_m+gamma_p*exp(1i*a*k);
        M(n+3,n+5)=gamma_m  +  gamma_AC3*exp(1i*a*k);
        M(n+4,n+5)=gamma_AC1_m+gamma_AC2*exp(-1i*a*k);
        if n<matrix_size-unitcell_size
            M(n+5,n+6)=gamma_p;
            M(n+5,n+7)=gamma_AC1_p+gamma_AC2*exp(-1i*a*k);
        end
        
        M(n+1,n)=gamma_p+gamma_m*exp(1i*a*k);
        M(n+2,n)=gamma_m*exp(1i*a*k)  +  gamma_AC3;
        M(n+2,n+1)=gamma_AC1_m+gamma_AC2*exp(1i*a*k);
        M(n+3,n+2)=gamma_p;
        M(n+4,n+2)=gamma_AC2+gamma_AC1_p*exp(-1i*a*k);
        M(n+4,n+3)=gamma_m+gamma_p*exp(-1i*a*k);
        M(n+5,n+3)=gamma_m  +  gamma_AC3*exp(-1i*a*k);
        M(n+5,n+4)=gamma_AC1_m+gamma_AC2*exp(1i*a*k);
        if n<matrix_size-unitcell_size
            M(n+6,n+5)=gamma_p;
            M(n+7,n+5)=gamma_AC1_p+gamma_AC2*exp(1i*a*k);
        end
        
    end
    
    matrix_size=size(M,1);
    
elseif strcmp(ribbon, 'trans3_mixed')
    
    unitcell_size=12;
    matrix_size=ribbon_size*unitcell_size;
    M=zeros(matrix_size);
    
    for n=1:unitcell_size:matrix_size
        
        M(n,n+2)=gamma_m;
        M(n,n+3)=gamma_AC1_m;
        M(n,n+8)=gamma_p*exp(1i*a*k);
        M(n,n+9)=gamma_AC2*exp(1i*a*k);
        M(n+1,n+3)=gamma_AC2;
        M(n+1,n+4)=gamma_m;
        M(n+1,n+5)=gamma_AC1_m;
        M(n+1,n+9)=gamma_AC1_p*exp(1i*a*k);
        M(n+1,n+10)=gamma_p*exp(1i*a*k);
        M(n+1,n+11)=gamma_AC2*exp(1i*a*k);
        M(n+2,n+3)=gamma_m;
        M(n+2,n+6)=gamma_p;
        M(n+3,n+4)=gamma_p;
        M(n+3,n+6)=gamma_AC2;
        M(n+3,n+7)=gamma_AC1_p;
        M(n+4,n+5)=gamma_m;
        M(n+4,n+7)=gamma_p;
        M(n+5,n+7)=gamma_AC2;
        M(n+6,n+8)=gamma_m;
        M(n+6,n+9)=gamma_AC1_m;
        M(n+7,n+9)=gamma_AC2;
        M(n+7,n+10)=gamma_m;
        M(n+7,n+11)=gamma_AC1_m;
        M(n+8,n+9)=gamma_m;
        M(n+9,n+10)=gamma_p;
        M(n+10,n+11)=gamma_m;
        if n<matrix_size-unitcell_size
            M(n+5,n+12)=gamma_AC2;
            M(n+5,n+14)=gamma_p;
            M(n+5,n+18)=gamma_AC1_p;
            M(n+11,n+12)=gamma_AC1_p*exp(-1i*a*k);
            M(n+11,n+18)=gamma_AC2;
            M(n+11,n+20)=gamma_p;
        end
        
        M(n+2,n)=gamma_m;
        M(n+3,n)=gamma_AC1_m;
        M(n+8,n)=gamma_p*exp(-1i*a*k);
        M(n+9,n)=gamma_AC2*exp(-1i*a*k);
        M(n+3,n+1)=gamma_AC2;
        M(n+4,n+1)=gamma_m;
        M(n+5,n+1)=gamma_AC1_m;
        M(n+9,n+1)=gamma_AC1_p*exp(-1i*a*k);
        M(n+10,n+1)=gamma_p*exp(-1i*a*k);
        M(n+11,n+1)=gamma_AC2*exp(-1i*a*k);
        M(n+3,n+2)=gamma_m;
        M(n+6,n+2)=gamma_p;
        M(n+4,n+3)=gamma_p;
        M(n+6,n+3)=gamma_AC2;
        M(n+7,n+3)=gamma_AC1_p;
        M(n+5,n+4)=gamma_m;
        M(n+7,n+4)=gamma_p;
        M(n+7,n+5)=gamma_AC2;
        M(n+8,n+6)=gamma_m;
        M(n+9,n+6)=gamma_AC1_m;
        M(n+9,n+7)=gamma_AC2;
        M(n+10,n+7)=gamma_m;
        M(n+11,n+7)=gamma_AC1_m;
        M(n+9,n+8)=gamma_m;
        M(n+10,n+9)=gamma_p;
        M(n+11,n+10)=gamma_m;
        if n<matrix_size-unitcell_size
            M(n+12,n+5)=gamma_AC2;
            M(n+14,n+5)=gamma_p;
            M(n+18,n+5)=gamma_AC1_p;
            M(n+12,n+11)=gamma_AC1_p*exp(1i*a*k);
            M(n+18,n+11)=gamma_AC2;
            M(n+20,n+11)=gamma_p;
        end
        
    end
    
    matrix_size=size(M,1);
    
elseif strcmp(ribbon, 'kagome2_flat')
    
    unitcell_size=3;
    matrix_size=ribbon_size*unitcell_size;
    M=zeros(matrix_size);
    
    for n=1:unitcell_size:matrix_size
        
        M(n,n+1)=1+exp(-1i*a*k);
        if n<=matrix_size-unitcell_size
            M(n,n+2)=1;
            M(n+1,n+2)=1;
            M(n+2,n+3)=1+exp(-1i*a*k);
            M(n+2,n+4)=1+exp(1i*a*k);
        end
        
        M(n+1,n)=1+exp(1i*a*k);
        if n<=matrix_size-unitcell_size
            M(n+2,n)=1;
            M(n+2,n+1)=1;
            M(n+3,n+2)=1+exp(1i*a*k);
            M(n+4,n+2)=1+exp(-1i*a*k);
        end
        
    end
    
    matrix_size=size(M,1);
    
elseif strcmp(ribbon, 'trans2_mixed')
    
    unitcell_size=3;
    matrix_size=ribbon_size*unitcell_size;
    M=zeros(matrix_size);
    
    for n=1:unitcell_size:matrix_size
        
        M(n,n+1)=1+exp(1i*a*k);
        M(n,n+2)=gamma1+gamma2*exp(1i*a*k);
        M(n+1,n+2)=1;
        if n<=matrix_size-unitcell_size
            M(n+2,n+3)=1;
            M(n+2,n+4)=gamma2+gamma1*exp(-1i*a*k);
        end
        
        M(n+1,n)=1+exp(-1i*a*k);
        M(n+2,n)=gamma1+gamma2*exp(-1i*a*k);
        M(n+2,n+1)=1;
        if n<=matrix_size-unitcell_size
            M(n+3,n+2)=1;
            M(n+4,n+2)=gamma2+gamma1*exp(1i*a*k);
        end
        
    end
    
    matrix_size=size(M,1);
    
elseif strcmp(ribbon, 'trans2_st_bearded')
    
    unitcell_size=3;
    matrix_size=ribbon_size*unitcell_size+1;
    M=zeros(matrix_size);
    
    for n=1:unitcell_size:matrix_size
        
        M(n,n+1)=gamma_p;
        M(n+1,n+2)=gamma_m+gamma_p*exp(1i*a*k);
        M(n,n+2)=gamma_AC2+gamma_AC1_p*exp(-1i*a*k);
        M(n+1,n+3)=gamma_m;
        M(n+2,n+3)=gamma_AC1_p+gamma_AC2*exp(1i*a*k);
        
        M(n+1,n)=gamma_p;
        M(n+2,n+1)=gamma_m+gamma_p*exp(-1i*a*k);
        M(n+2,n)=gamma_AC2+gamma_AC1_p*exp(1i*a*k);
        M(n+3,n+1)=gamma_m;
        M(n+3,n+2)=gamma_AC1_p+gamma_AC2*exp(-1i*a*k);
        
    end
    
    matrix_size=size(M,1);
    
elseif strcmp(ribbon, 'trans_st_bearded')
    
    unitcell_size=3;
    matrix_size=ribbon_size*unitcell_size+1;
    M=zeros(matrix_size);
    
    for n=1:unitcell_size:matrix_size
        
        M(n,n+1)=1;
        M(n+1,n+2)=1+exp(1i*a*k);
        M(n,n+2)=gamma2+gamma1*exp(-1i*a*k);
        M(n+1,n+3)=1;
        M(n+2,n+3)=gamma1+gamma2*exp(1i*a*k);
        
        M(n+1,n)=1;
        M(n+2,n+1)=1+exp(-1i*a*k);
        M(n+2,n)=gamma2+gamma1*exp(1i*a*k);
        M(n+3,n+1)=1;
        M(n+3,n+2)=gamma1+gamma2*exp(-1i*a*k);
        
    end
    
    matrix_size=size(M,1);
    
elseif strcmp(ribbon, 'graphene_zigzag2')
    
    unitcell_size=4;
    matrix_size=ribbon_size*unitcell_size;
    M=zeros(matrix_size);
    
    for n=1:unitcell_size:matrix_size
        
        M(n,n+1)=1+exp(1i*k);
        M(n+1,n+2)=1;
        M(n+2,n+3)=1+exp(-1i*k);
        if n<=matrix_size-unitcell_size
            M(n+3,n+4)=1;
        end
        
        M(n+1,n)=conj(M(n,n+1));
        M(n+2,n+1)=conj(M(n+1,n+2));
        M(n+3,n+2)=conj(M(n+2,n+3));
        if n<=matrix_size-unitcell_size
            M(n+4,n+3)=1;
        end
        
    end
    
    matrix_size=size(M,1);
    
elseif strcmp(ribbon, 'graphene_zigzag3')
    
    unitcell_size=4;
    matrix_size=ribbon_size*unitcell_size;
    M=zeros(matrix_size);
    
    for n=1:unitcell_size:matrix_size
        
        M(n,n+2)=1+exp(1i*k);
        M(n+1,n+3)=1+exp(1i*k);
        M(n+2,n+3)=1;
        if n<=matrix_size-unitcell_size
            M(n+1,n+4)=1;
        end
        
        M(n+2,n)=conj(M(n,n+2));
        M(n+3,n+1)=conj(M(n+1,n+3));
        M(n+3,n+2)=conj(M(n+2,n+3));
        if n<=matrix_size-unitcell_size
            M(n+4,n+1)=1;
        end
        
    end
    
    matrix_size=size(M,1);
    
elseif strcmp(ribbon, 'lieb_st_y_bearded')
    
    matrix_size=ribbon_size*3+1;
    M=zeros(matrix_size);
    
    for n=1:3:matrix_size
        
        M(n,n+1)=1;
        M(n+1,n+2)=exp(1i*k*(0.5-delta));
        M(n,n+2)=gamma0*(1+exp(1i*k*(0.5-delta)));
        
        M(n+1,n)=1;
        M(n+2,n+1)=exp(-1i*k*(0.5-delta));
        M(n+2,n)=gamma0*(1+exp(-1i*k*(0.5-delta)));
        
        M(n+1,n+3)=1;
        M(n+2,n+3)=gamma0*(1+exp(1i*k*(0.5+delta)));
        
        M(n+3,n+1)=1;
        M(n+3,n+2)=gamma0*(1+exp(-1i*k*(0.5+delta)));
        
    end
    
    matrix_size=size(M,1);
    
elseif strcmp(ribbon, 'graphene_y_zigzag')
    
    matrix_size=ribbon_size*2;
    M=zeros(matrix_size);
    
    for n=1:2:matrix_size
        
        M(n,n+1)=t*(1+exp(-1i*k));
        M(n+1,n)=t*(1+exp(1i*k));
        if (n+2<=ribbon_size*2)
            M(n+1,n+2)=t;
            M(n+2,n+1)=t;
        end
    end
    
    matrix_size=size(M,1);
    
elseif strcmp(ribbon, 'trans_y_bearded')
    
    matrix_size=ribbon_size*3+1;
    M=zeros(matrix_size);
    
    for n=1:3:matrix_size
        
        M(n,n+1)=1+exp(-1i*sind(theta)*k/2);
        M(n+1,n+2)=exp(1i*k/2);
        M(n,n+2)=gamma1*(1+exp(1i*(1-sind(theta))*k/2));
        
        M(n+1,n)=1+exp(1i*sind(theta)*k/2);
        M(n+2,n+1)=exp(-1i*k/2);
        M(n+2,n)=gamma1*(1+exp(-1i*(1-sind(theta))*k/2));
        
        M(n+1,n+3)=1+exp(-1i*sind(theta)*k/2);
        M(n+2,n+3)=gamma2*(1+exp(1i*(1+sind(theta))*k/2));
        
        M(n+3,n+1)=1+exp(1i*sind(theta)*k/2);
        M(n+3,n+2)=gamma2*(1+exp(-1i*(1+sind(theta))*k/2));
        
    end
    
    matrix_size=size(M,1);
    
elseif strcmp(ribbon, 'kagome_x')
    
    matrix_size=ribbon_size*6;
    M=zeros(matrix_size);
    
    for n=1:6:matrix_size
        
        M(n,n+1)=exp(-1i*k);
        M(n+1,n+2)=1+exp(-1i*k);
        M(n+2,n+3)=1+exp(1i*k);
        M(n+3,n+4)=exp(1i*k);
        M(n+3,n+5)=1+exp(1i*k);
        M(n+4,n+5)=1+exp(-1i*k);
        
        M(n+1,n)=exp(1i*k);
        M(n+2,n+1)=1+exp(1i*k);
        M(n+3,n+2)=1+exp(-1i*k);
        M(n+4,n+3)=exp(-1i*k);
        M(n+5,n+3)=1+exp(-1i*k);
        M(n+5,n+4)=1+exp(+1i*k);
        
    end
    
    matrix_size=size(M,1);
    
elseif strcmp(ribbon, 'kagome2_x')
    
    matrix_size=ribbon_size*4;
    M=zeros(matrix_size);
    
    for n=1:4:matrix_size
        
        M(n,n+1)=exp(1i*k);
        M(n+1,n+2)=exp(1i*k);
        M(n+2,n+3)=1+exp(1i*k);
        if (n<=matrix_size-4)
            M(n+3,n+4)=1+exp(1i*k);
        end
        M(n+1,n)=exp(-1i*k);
        M(n+2,n+1)=exp(-1i*k);
        M(n+3,n+2)=1+exp(-1i*k);
        if (n<=matrix_size-4)
            M(n+4,n+3)=1+exp(-1i*k);
        end
    end
    
    matrix_size=size(M,1);
    
elseif strcmp(ribbon, 'kagome3_x')
    
    matrix_size=ribbon_size*3;
    M=zeros(matrix_size);
    
    for n=1:3:matrix_size
        
        M(n,n+1)=exp(1i*k); %12
        M(n+1,n+2)=1+exp(-1i*k/2); %23
        M(n,n+2)=1+exp(1i*k/2); %13
        if (n<=matrix_size-3)
            M(n+2,n+3)=1+exp(1i*k/2); %34
        end
        
        M(n+1,n)=exp(-1i*k); %21
        M(n+2,n+1)=1+exp(1i*k/2); %32
        M(n+2,n)=1+exp(-1i*k/2); %31
        if (n<=matrix_size-3)
            M(n+3,n+2)=1+exp(-1i*k/2); %43
        end
        
    end
    
    matrix_size=size(M,1);
    
elseif strcmp(ribbon,'trans2')
    
    AB=exp((1i*k*sind(theta))/2);
    BC=exp(1i*k/2);
    AC=gamma1*exp(1i*k*(1-sind(theta))/2) + gamma2*exp(1i*k*(1+sind(theta))/2);
    
    matrix_size=ribbon_size*3;
    M=zeros(matrix_size);
    
    for n=1:3:matrix_size
        
        M(n,n+1)=AB;
        M(n,n+2)=AC;
        M(n+1,n+2)=BC;
        
        M(n+1,n)=conj(AB);
        M(n+2,n)=conj(AC);
        M(n+2,n+1)=conj(BC);
        
    end
    
    matrix_size=size(M,1);
    
end

end

%% function Floquet
% Ribbon bandstructure calculation for floquet topological insulators
% (helical waveguides).
%
% The matrix has to be separated into matrices for each possible
% "time dependance" (z dependance), because in matlab a function has to be
% defined in a single line of code. There is a different time dependance
% for each possible vector connecting two waveguides. The time-dependant
% Hamiltonian is then used to calculate the time evolution operator. From
% the eigenvalues of that quasienergies are calculated, which correspond to
% the energy bands. These are however unsorted, so they should be plotted as
% dots instead of lines.

function[matrix_size, M] = Floquet(floquet, ribbon_size, a, theta, nexp, k0, R, omega, omega2, k)

gamma1=(exp(1-sqrt(2-2*sind(theta))))^nexp;  % t_nnn/t_nn for transition lattices
gamma2=(exp(1-sqrt(2+2*sind(theta))))^nexp;

if strcmp(floquet, 'trans4_mixed') || strcmp(floquet, 'trans') % set unitcell_size
    unitcell_size=6;
    matrix_size=ribbon_size*unitcell_size;
elseif strcmp(floquet, 'graphene_y_zigzag')
    unitcell_size=2;
    matrix_size=ribbon_size*unitcell_size;
end

if strcmp(floquet, 'trans4_mixed') % obsolete
    
    Ms=zeros(matrix_size);
    Mtm=zeros(matrix_size);
    Mtp=zeros(matrix_size);
    
    for n=1:unitcell_size:matrix_size
        
        Ms(n,n+1)=1;
        Mtm(n,n+1)=1;
        Mtm(n,n+2)=1;
        Ms(n+1,n+2)=gamma1;
        Mtm(n+1,n+2)=gamma2;
        Ms(n+2,n+3)=1;
        Ms(n+2,n+4)=gamma2;
        Mtp(n+2,n+4)=gamma1;
        Ms(n+3,n+4)=1;
        Mtp(n+3,n+4)=1;
        Ms(n+3,n+5)=1;
        Ms(n+4,n+5)=gamma1;
        Mtm(n+4,n+5)=gamma2;
        if n<matrix_size-unitcell_size
            Ms(n+5,n+6)=1;
            Ms(n+5,n+7)=gamma1;
            Mtm(n+5,n+7)=gamma2;
        end
        
        Ms(n+1,n)=1;
        Mtp(n+1,n)=1;
        Mtp(n+2,n)=1;
        Ms(n+2,n+1)=gamma1;
        Mtp(n+2,n+1)=gamma2;
        Ms(n+3,n+2)=1;
        Ms(n+4,n+2)=gamma2;
        Mtm(n+4,n+2)=gamma1;
        Ms(n+4,n+3)=1;
        Mtm(n+4,n+3)=1;
        Ms(n+5,n+3)=1;
        Ms(n+5,n+4)=gamma1;
        Mtp(n+5,n+4)=gamma2;
        if n<matrix_size-unitcell_size
            Ms(n+6,n+5)=1;
            Ms(n+7,n+5)=gamma1;
            Mtp(n+7,n+5)=gamma2;
        end
        
    end
    
    M=@(t) (Ms + (Mtp.*exp(1i.*a.*(k+R.*sin(omega.*t)))) + (Mtm.*exp(-1i.*a.*(k+R.*sin(omega.*t)))));
    %         M1=@(t) Mtp.*exp(1i*a*(k+R*omega*sin(omega*t)));
    %         M2=@(t) Mtm.*exp(-1i*a*(k+R*omega*sin(omega*t)));
    %         M=@(t) (Ms+M1(t)+M2(t));
    
elseif strcmp(floquet, 'trans') % latest
    
    % initialize matrices. One matrix for each possible vector
    % connecting waveguides.
    Mxp=zeros(matrix_size);
    Mxm=zeros(matrix_size);
    Myp=zeros(matrix_size);
    Mym=zeros(matrix_size);
    Mac1p=zeros(matrix_size);
    Mac1m=zeros(matrix_size);
    Mac2p=zeros(matrix_size);
    Mac2m=zeros(matrix_size);
    
    for n=1:unitcell_size:matrix_size
        
        Mxp(n,n+1)=exp(-1i*a*k);
        Mxp(n+3,n+4)=1;
        %             Mxp(n+6,n+7)=exp(-1i*a*k);
        
        Mxp(n+1,n)=1;
        Mxp(n+4,n+3)=exp(-1i*a*k);
        %             Mxp(n+7,n+6)=1;
        
        Mxm(n,n+1)=1;
        Mxm(n+3,n+4)=exp(1i*a*k);
        %             Mxm(n+6,n+7)=1;
        
        Mxm(n+1,n)=exp(1i*a*k);
        Mxm(n+4,n+3)=1;
        %             Mxm(n+7,n+6)=exp(1i*a*k);
        
        Myp(n,n+2)=exp(-1i*a*k);
        Myp(n+2,n+3)=1;
        Myp(n+3,n+5)=1;
        if n<matrix_size-unitcell_size
            Myp(n+5,n+6)=1;
        end
        
        Mym(n+2,n)=exp(1i*a*k);
        Mym(n+3,n+2)=1;
        Mym(n+5,n+3)=1;
        if n<matrix_size-unitcell_size
            Mym(n+6,n+5)=1;
        end
        
        Mac1m(n+1,n+2)=gamma1;
        Mac1m(n+4,n+5)=gamma1;
        if n<matrix_size-unitcell_size
            Mac1m(n+5,n+7)=gamma1;
        end
        Mac1m(n+2,n+4)=gamma1*exp(1i*a*k);
        
        Mac1p(n+2,n+1)=gamma1;
        Mac1p(n+5,n+4)=gamma1;
        if n<matrix_size-unitcell_size
            Mac1p(n+7,n+5)=gamma1;
        end
        Mac1p(n+4,n+2)=gamma1*exp(-1i*a*k);
        
        Mac2p(n+1,n+2)=gamma2*exp(-1i*a*k);
        Mac2p(n+2,n+4)=gamma2;
        Mac2p(n+4,n+5)=gamma2*exp(-1i*a*k);
        if n<matrix_size-unitcell_size
            Mac2p(n+5,n+7)=gamma2*exp(-1i*a*k);
        end
        
        Mac2m(n+2,n+1)=gamma2*exp(1i*a*k);
        Mac2m(n+4,n+2)=gamma2;
        Mac2m(n+5,n+4)=gamma2*exp(1i*a*k);
        if n<matrix_size-unitcell_size
            Mac2m(n+7,n+5)=gamma2*exp(1i*a*k);
        end
        
    end
    
    % time-dependent ribbon-matrix.
    %         M = @(t) (Mxp*exp(1i*((a*e)/(hq*c))*k0*R*omega2*sin(omega*t)) + Mxm*exp(-1i*((a*e)/(hq*c))*k0*R*omega2*sin(omega*t)) + Myp*exp(1i*((a*e)/(hq*c))*k0*R*omega2*(sind(theta)*sin(omega*t)-cosd(theta)*cos(omega*t))) + Mym*exp(-1i*((a*e)/(hq*c))*k0*R*omega2*(sind(theta)*sin(omega*t)-cosd(theta)*cos(omega*t))) + Mac1p*exp(1i*((a*e)/(hq*c))*k0*R*omega2*((1-sind(theta))*sin(omega*t)+cosd(theta)*cos(omega*t))) + Mac1m*exp(-1i*((a*e)/(hq*c))*k0*R*omega2*((1-sind(theta))*sin(omega*t)+cosd(theta)*cos(omega*t))) + Mac2p*exp(1i*((a*e)/(hq*c))*k0*R*omega2*((1+sind(theta))*sin(omega*t)-cosd(theta)*cos(omega*t))) + Mac2m*exp(-1i*((a*e)/(hq*c))*k0*R*omega2*((1+sind(theta))*sin(omega*t)-cosd(theta)*cos(omega*t))))*(-1)*145;
    M = @(t) (Mxp*exp(1i*a*k0*R*omega2*sin(omega*t)) + Mxm*exp(-1i*a*k0*R*omega2*sin(omega*t)) + Myp*exp(1i*a*k0*R*omega2*(sind(theta)*sin(omega*t)-cosd(theta)*cos(omega*t))) + Mym*exp(-1i*a*k0*R*omega2*(sind(theta)*sin(omega*t)-cosd(theta)*cos(omega*t))) + Mac1p*exp(1i*a*k0*R*omega2*((1-sind(theta))*sin(omega*t)+cosd(theta)*cos(omega*t))) + Mac1m*exp(-1i*a*k0*R*omega2*((1-sind(theta))*sin(omega*t)+cosd(theta)*cos(omega*t))) + Mac2p*exp(1i*a*k0*R*omega2*((1+sind(theta))*sin(omega*t)-cosd(theta)*cos(omega*t))) + Mac2m*exp(-1i*a*k0*R*omega2*((1+sind(theta))*sin(omega*t)-cosd(theta)*cos(omega*t))))*(-1)*145;
    %         M = @(t) (Mxp*exp(-1i*a*R*omega2*cos(omega*t)) + Mxm*exp(1i*a*R*omega2*cos(omega*t)) + Myp*exp(1i*a*R*omega2*(-sind(theta)*cos(omega*t)+cosd(theta)*sin(omega*t))) + Mym*exp(-1i*a*R*omega2*(-sind(theta)*cos(omega*t)+cosd(theta)*sin(omega*t))) + Mac1p*exp(1i*a*R*omega2*(-(1-sind(theta))*cos(omega*t)-cosd(theta)*sin(omega*t))) + Mac1m*exp(-1i*a*R*omega2*(-(1-sind(theta))*cos(omega*t)-cosd(theta)*sin(omega*t))) + Mac2p*exp(1i*a*R*omega2*(-(1+sind(theta))*cos(omega*t)+cosd(theta)*sin(omega*t))) + Mac2m*exp(-1i*a*R*omega2*(-(1+sind(theta))*cos(omega*t)+cosd(theta)*sin(omega*t)))); % alternative definition of R_ij=r_j-r_i or R_ij=r_i-r_j
    
elseif strcmp(floquet, 'graphene_y_zigzag')
    
    Mxp=zeros(matrix_size);
    Mxm=zeros(matrix_size);
    Myp=zeros(matrix_size);
    Mym=zeros(matrix_size);
    
    for n=1:unitcell_size:matrix_size
        
        Myp(n,n+1)=exp(-1i*a*k);
        Mym(n+1,n)=exp(1i*a*k);
        if (n+2<=matrix_size)
            Mxp(n+1,n+2)=1;
            Mxm(n+2,n+1)=1;
        end
        
    end
    
    M = @(t) (Mxp*exp(1i*a*k0*R*omega*sin(omega*t)) + Mxm*exp(-1i*a*k0*R*omega*sin(omega*t)) + Myp*exp(1i*0.5*a*(k0*R*omega*sin(omega*t)+sqrt(3)*k0*R*omega*cos(omega*t))) + Mym*exp(-1i*0.5*a*(k0*R*omega*sin(omega*t)+sqrt(3)*k0*R*omega*cos(omega*t))))*145;
    
end

end
