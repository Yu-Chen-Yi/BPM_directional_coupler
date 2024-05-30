classdef ChenYi_BPM
    properties
        lambda double = 1.55;
        xmin double = -10;
        xmax double = 10;
        xmesh double = 0.0025;
        zmesh double = 0.005;
        zsample double = 20000;
        Waveguide struct = struct(...
            'core_index', 3.4777, ...
            'cladding_index', 1.444, ...
            'gap', 0.15, ...
            'Width', 0.45, ...
            'profile',[]);
        
        Source struct = struct(...
            'amplitude', 1, ...
            'postion', 0.15/2+0.45/2, ...
            'bandwidth', 0.1, ...
            'Efield', 0);
        Result struct = struct(...
            'x',[],...
            'z',[],...
            'I',[],...
            'Lc',[],...
            'Delta_n',[]);
    end
    methods
        function [lambda, xa, xb, dx, deltaz, zmax, nSi, nSiO2, gap, wg_width, A, x0, W0] = getParameters(obj)
            lambda = obj.lambda;
            xa = obj.xmin;
            xb = obj.xmax;
            dx = obj.xmesh;
            deltaz = obj.zmesh;
            zmax = obj.zsample;
            nSi = obj.Waveguide.core_index;
            nSiO2 = obj.Waveguide.cladding_index;
            gap = obj.Waveguide.gap;
            wg_width = obj.Waveguide.Width;
            A = obj.Source.amplitude;
            x0 = obj.Source.postion;
            W0 = obj.Source.bandwidth;
        end
        function obj = StructureView(obj)
            [~, xa, xb, dx, deltaz, zmax, nSi, nSiO2, gap, wg_width, A, x0, W0] = getParameters(obj);
            x = xa:dx:xb;
            obj.Result.x = x;
            z = (0:10:zmax)*deltaz;
            obj.Result.z = z;
            %%-- Set 1st waveguide position coordinate -------%%
            x1 = -wg_width-gap/2;       % left side of 1st waveguide step  (um)
            x2 = -gap/2;                % right side of 1st waveguide step (um)
            %%-- Set 2st waveguide position coordinate -------%%
            x3 = gap/2;                 % left side of 2st waveguide step  (um)
            x4 = wg_width+gap/2;        % right side of 2st waveguide step (um)
            % ---------- Gaussian source Generation ---------------------------
            Efield = A*exp (-((x+x0)/W0).^2);  % Generate gaussian source
            obj.Source.Efield = Efield;
            figure(1);
            subplot(3,1,1)
            % Plots the gaussian source
            colororder({'r','b'})
            yyaxis left
            plot(x,Efield,'r-','linewidth',3);grid on; xlim([xa xb]);xlim([x1*3 x4*3]);
            goodfigure('Gaussian source & Profile of index','w',[0 0.1 0.7 0.7]);
            goodplot2('\itx (\mum)','\itAmplitude (V/m)',sprintf('%d*e^{-(x+%g/%g)^2}',A,x0,W0),24)
            set(gca,'xDir','reverse');
            plot_darkmode;
            % --------------- Calculation the average index -----------------------
            n = ((x >= x1) & (x <= x2)) | ((x >= x3) & (x <= x4));
            n = +n*nSi + +(~n)*nSiO2;
            obj.Waveguide.profile = repmat(n,[zmax/10+1,1]);
            yyaxis right
            plot(x,n,'b-','linewidth',3);xlim([xa xb]);ylim([0 nSi]);xlim([x1*3 x4*3]);
            goodplot2('\itx (\mum)','\itRefractive index',sprintf('%d*e^{-(x+%g/%g)^2}',A,x0,W0),24)
            set(gca,'xDir','reverse');
            legend('Source','Refractive index');
            plot_darkmode;
            % --------------- Calculation propagate angular frequency with attenuation -----------------
            PML = (x <= x1*1.5) | (x >= x4*1.5);
            PML = +PML*1000;
            figure(1);
            subplot(3,1,2)
            plot(x,PML,'g','linewidth',3);
            goodplot2('\itx (\mum)','\itAttenuation constant','attenuation region',24); xlim([xa xb]);
            xlim([x1*3 x4*3]);
            set(gca,'xDir','reverse');
            legend('Attenuation constant')
            plot_darkmode;
            subplot(3,1,3)
            imagesc(x,z,obj.Waveguide.profile);colormap jet;colorbar;
            goodplot2('\itx (\mum)','\itz (\mum)','',24)
            set(gca,'xDir','reverse');
            plot_darkmode;
        end
        function obj = BPM_Run(obj)
            [lambda0, xa, xb, dx, deltaz, zmax, nSi, nSiO2, gap, wg_width, A, x0, W0] = getParameters(obj);
            x = xa:dx:xb;
            z = (0:10:zmax)*deltaz;
            x1 = -wg_width-gap/2;       % left side of 1st waveguide step  (um)
            x2 = -gap/2;                % right side of 1st waveguide step (um)
            x3 = gap/2;                 % left side of 2st waveguide step  (um)
            x4 = wg_width+gap/2;        % right side of 2st waveguide step (um)
            k0=(2*pi/lambda0);           % Wave number (1/um)
            % ---------- Gaussian source Generation ---------------------------
            Efield = A*exp (-((x+x0)/W0).^2);  % Generate gaussian source
            % --------------- Calculation the average index -----------------------
            navg = (nSi + nSiO2)/2 * 1.18;
            n = ((x >= x1) & (x <= x2)) | ((x >= x3) & (x <= x4));
            n = +n*nSi + +(~n)*nSiO2;
            num_samples = length(x);
            Lx = xb - xa;               % Width : Space domain (um)
            % --------------- Calculation of spactial frequency and propagate angular frequency -------------------
            kx = linspace(-num_samples/2,num_samples/2-1,num_samples);
            kx = (2*pi/Lx)*fftshift(kx);
            phase1 = exp((1i*deltaz*kx.^2)./(navg*k0 + sqrt(max(0,navg^2*k0*2 - kx.^2))));
            PML = (x <= x1*1.5) | (x >= x4*1.5);
            PML = +PML*1000;
            phase2 = exp(-(PML + 1i*(n - navg)*k0)*deltaz);
            %% --------------- Propagation forloop cycle------------------------
            E = zeros(length(z),length(x));
            for k = 0:1:zmax/10
                for j = 1:1:10 % Propagate 10 steps to save the next element of E(z,x)
                    % ASBPM
                    Efield = ifft((fft(Efield).*phase1)).*phase2;
                end
                E(k+1,:) = Efield;
            end
            I = abs(E).^2;
            Imax = max(max(I(zmax/10/2:end,:)));
            I = I/Imax;
            [~,idx_WG1] = min(abs(x-x1));
            [~,idx_WG2] = min(abs(x-x2));
            P1 = sum(I(:,idx_WG1:idx_WG2),2);
            norm_val = max(P1(zmax/10/2:end));
            P1 = P1/norm_val;
            %% fit Lc & attenuation decay
            fitnessFunction = @(x) sum(abs(P1' - (x(1)*cos(z/x(2)*pi).*exp(-z*x(3))).^2).^2);
            numVariables = 3;
            lowerBound = [1, 3  , 0.001];
            upperBound = [5, 200, 0.01];
            maxGenerations = 200;
            populationSize = 700;
            [p,~] = greyWolfOptimization(fitnessFunction, numVariables, lowerBound, upperBound, maxGenerations, populationSize);
            %% normalize
            I = 1/p(1)^2*I.*exp(z'*p(3)*2);
            Lc = p(2);
            obj.Result.I = I;
            obj.Result.Lc = Lc;
            obj.Result.Delta_n = 2*lambda0/Lc;
        end
        function obj = Field_SideView(obj)
            x = obj.Result.x;
            z = obj.Result.z;
            I = obj.Result.I;
            Lc = obj.Result.Lc;
            [~, ~, ~, ~, ~, ~, ~, ~, gap, wg_width, ~, ~, ~] = getParameters(obj);
            x1 = -wg_width-gap/2;       % left side of 1st waveguide step  (um)
            x4 = wg_width+gap/2;        % right side of 2st waveguide step (um)
            figure(3);
            %% draw figure 311
            % ---------- Meshgrid x & z coordinate for plotting propagated field -------------------
            [X,Z] = meshgrid (x,z);
            subplot(3,1,1)
            mesh(X,Z,I,'FaceColor','interp','edgecolor','none');
            color = fliplr(jet);
            color = [color(1,:);color(51,:);color(102,:);color(153,:);color(204,:);color(256,:)];
            mycolormap = customcolormap([linspace(0,0.9,6) 1], [color;0 0 0]);
            colormap(mycolormap);
            camlight right;grid on;set(gca,'xDir','reverse');set(gca,'yDir','normal');
            view(70,60);xlim([x1*1.5 x4*1.5]);
            ylim([0 Lc*5]);
            zlim([0 1]);
            if exist('clim', 'file') == 2
                % ㄏノclimㄧ计
                clim([0 1]);
            else
                % ㄏノcaxisㄧ计
                caxis([0 1]);
            end
            colorbar;
            goodplot2('\itx axis (\mum)','\itz axis (\mum)','Magnitude of the Propagated Gaussian Source',16);
            goodfigure('Top view Propagation mode','w',[0.6 0 0.4 1]);
            ylim('tight')
            plot_darkmode;
        end
        function obj = Field_topView(obj)
            x = obj.Result.x;
            z = obj.Result.z;
            I = obj.Result.I;
            Lc = obj.Result.Lc;
            [~, ~, ~, ~, ~, ~, ~, ~, gap, wg_width, ~, ~, ~] = getParameters(obj);
            x1 = -wg_width-gap/2;       % left side of 1st waveguide step  (um)
            x3 = gap/2;                 % left side of 2st waveguide step  (um)
            x4 = wg_width+gap/2;        % right side of 2st waveguide step (um)
            %% draw figure 312
            figure(3);
            subplot(3,1,2)
            imagesc(z,x,I');hold on;
            set(gca,'YDir','normal');
            rectangle('Position',[min(z) x1 max(z)-min(z) wg_width],'EdgeColor','w','linewidth',1);
            rectangle('Position',[min(z) x3 max(z)-min(z) wg_width],'EdgeColor','w','linewidth',1);
            hold off
            ylim([x1*3 x4*3]);
            xlim([0 Lc*5])
            set(gca,'Layer','top');
            color = fliplr(jet);
            color = [color(1,:);color(51,:);color(102,:);color(153,:);color(204,:);color(256,:)];
            mycolormap = customcolormap([linspace(0,0.9,6) 1], [color;0 0 0]);
            colormap(mycolormap);colorbar;
            if exist('clim', 'file') == 2
                % ㄏノclimㄧ计
                clim([0 1]);
            else
                % ㄏノcaxisㄧ计
                caxis([0 1]);
            end
            goodplot2('\itz axis (\mum)','\itx axis (\mum)','Magnitude of the Propagated Gaussian Source',16);
            xlim('tight')
            plot_darkmode;
        end
        function obj = Coupling_efficeincy(obj)
            x = obj.Result.x;
            z = obj.Result.z;
            I = obj.Result.I;
            Lc = obj.Result.Lc;
            [~, ~, ~, ~, ~, zmax, ~, ~, gap, wg_width, ~, ~, ~] = getParameters(obj);
            %%-- Set 1st waveguide position coordinate -------%%
            x1 = -wg_width-gap/2;       % left side of 1st waveguide step  (um)
            x2 = -gap/2;                % right side of 1st waveguide step (um)
            %%-- Set 2st waveguide position coordinate -------%%
            x3 = gap/2;                 % left side of 2st waveguide step  (um)
            x4 = wg_width+gap/2;        % right side of 2st waveguide step (um)
            %% draw figure 313
            [~,idx_WG1] = min(abs(x-x1));
            [~,idx_WG2] = min(abs(x-x2));
            [~,idx_WG3] = min(abs(x-x3));
            [~,idx_WG4] = min(abs(x-x4));
            P1 = sum(I(:,idx_WG1:idx_WG2),2);
            norm_val = max(P1(zmax/10/2:end));
            P1 = P1/norm_val;
            P2 = sum(I(:,idx_WG3:idx_WG4),2);
            P2 = P2/norm_val;
            figure(3);
            subplot(3,1,3)
            set(gcf,'name','BPM')
            plot(z,P1,'color','r','linewidth',2);hold on;
            plot(z,P2,'color','b','linewidth',2);hold off;
            goodplot2('Lc (\mum)','Coupling ratio','Beam propagation method',16);
            xlim([0 Lc*5])
            ylim([-0.05 1.05])
            legend('BPM-P_1','BPM-P_2','Location','best');
            xlim('tight')
            plot_darkmode;
            colorbar;
        end
    end
end