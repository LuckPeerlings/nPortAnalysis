classdef PronyMethod < handle
    %This class is used to determine the impedance using the so called
    %Prony method.
    %The details of this method are taken from the paper
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    %TODO List
    % - Normalization or any other method to be able to define epsilon in a
    %   logical matter
    % - Fix the the mylittle prony
    % - Fix the code so that it works for inputs of P which are matrices
    % - Include code to be able to get the data on regular intervals
    % - Include a comparison of the different methods with the inclusion of
    %   noise on the input data.
    % - Comment the code with references to the paper where necessary
    
    properties
        %Vector of the measured pressures
        Frequency
        P
        MicSpacing %The microphone spacing
        NrModes
        Amplitude
        kappa
        Parameter
        Epsilon
    end
    
    methods
        %Empty Constructor 
        %Constructor
        function obj = PronyMethod(Frequency,P,MicSpacing,NrModes,Parameter)
            if nargin > 1 
                obj.Frequency = Frequency; 
                obj.P = P;
                obj.MicSpacing = MicSpacing;
                obj.NrModes = NrModes;
                obj.Parameter = Parameter;
            elseif nargin == 1
                obj.TestFunctionNicolas
            else
                
                obj.TestClass;
            end
        end
        
        function obj = DetermineImpedance(obj)
            
        end
        
        function obj = DetermineWaveNumber(obj)
            %Calculating the wavenumber and the amplitude of the modes for
            %each frequency
            for ff = 1:size(obj.P,2)
                    [Amplitude_temp,kappa_temp] = obj.ESPRIT(obj.P(:,ff).',obj.dX,obj.NrModes, 0.1); % 1e-10
                    obj.Amplitude(:,ff) = [Amplitude_temp ; zeros(obj.NrModes - length(Amplitude_temp), 1)];
                    obj.kappa(:,ff) = [kappa_temp ; zeros(obj.NrModes - length(kappa_temp), 1)];
%                     [Amplitude_temp,kappa_temp] = obj.MyLittleProny(obj.P(:,ff).',obj.dX,obj.NrModes, 1e-10);
%                     obj.Amplitude(:,ff) = [Amplitude_temp ; zeros(obj.NrModes - length(Amplitude_temp), 1)];
%                     obj.kappa(:,ff) = [transpose(kappa_temp) ; zeros(obj.NrModes - length(kappa_temp), 1)];
%                     [obj.Amplitude(:,ff),obj.kappa(:,ff)] = obj.PronyCoefficients(obj.P(:,ff).',obj.dX,obj.NrModes);                end
            end
        end
        
        function obj = PlotWaveNumberComplexDomain(obj,ModeNr)
            nn = ModeNr;
            
            [~,I] = sort(obj.Frequency);
            figure;
            AX  = axes();
            hold on;
            PronyMethod.ArrowPlotComplexDomain(AX,obj.kappa(nn,I))
            axis square;
            
            ylabel('Imag. part of wavenumber')  
            xlabel('Real part of wavenumber')
           
        end
        function obj = PlotAmplitude(obj)
            %Function to plot the obtained real and imaginary part of each
            %of the modes and the relative amplitude of the modes compared to the largest mode present. 
            
            figure
            for nn = 1:obj.NrModes
                subplot(obj.NrModes+1,1,nn)
                plot(obj.Frequency,real(obj.kappa(nn,:)),'x');
                hold all
                plot(obj.Frequency,imag(obj.kappa(nn,:)),'o');
            end
            legend('Real part of wavenumber','Imag. part of wavenumber')  
            
           
            AX1 = subplot(obj.NrModes+1,1,obj.NrModes+1);            
            for ff = 1:size(obj.Frequency,2)
                plot(obj.Frequency(ff),abs(obj.Amplitude(:,ff))/norm(obj.Amplitude(:,ff)),'x')
                hold on         
                AX1.ColorOrderIndex = 1;
                ylabel('Rel. Amplitude')
                xlabel('Frequency')
            end
        end
        
        function TestFunctionNicolas(obj)
        %function to determine whether the equispacing function works well
        
        %different signals to be tested 
        X_test=0.01+transpose(0:0.055:5*0.055);
        t= transpose(0:0.01:1);
        func = @(x) exp(1i*(-11.510242292516846 + 6.292327191264097i)*x);
        Testsignal=func(X_test);
        Testsignal2=func(t);
        %function call: 
        % 1. Parameter: KxN Vector of space locations
        % 2. Parameter: KxN Vector of testfunction
        % 3. Parameter: Value of L: K>=L>=N with N= nonharmonic bandwidth
        % 4. Parameter: Vaue of b:  Variance of the gaussian
        
        [equiPressure, equiX] =PronyMethod.equispacing(X_test,Testsignal,6,12.3);
        [Testdata, equigrid] =PronyMethod.equispacing(t,Testsignal2,99,12.3);
        
%        MSE of the approximated Data and the Testdata
         error=equiPressure-func(equiX);

        figure
        subplot(3,2,1)
        plot(equigrid,real(Testdata),'x')
        ylabel('Reconstructed Amplitude')
        xlabel('Equispaced Steps')
        hold on
        plot(t,real(Testsignal2))
        plot(equiX,real(equiPressure),'o')
        title('Real part')
      
        
        %imaginary Part
        
        subplot(3,2,2)
        plot(equigrid,imag(Testdata),'x')
        ylabel('Reconstructed')
        xlabel('Equispaced Steps')
        hold on
        plot(t,imag(Testsignal2))
        plot(equiX,imag(equiPressure),'o')
        title('Imaginary part')
        legend('generated Data','Testsignal','Measurement points')
        
        
        % Plot errrors of real Part  
        subplot(3,2,3)
        plot(equiX,real(error))
        ylabel('Error of reconstruction')
        xlabel('Equispaced Steps')
        hold on
        plot(t,zeros(length(t)));
        
        % Plot errrors of imag Part  
        subplot(3,2,4)
        plot(equiX,imag(error))
        ylabel('Error of reconstruction')
        xlabel('Equispaced Steps')        
        hold on
        plot(t,zeros(length(t)));
        
        % plot of the input with non equispaced nodes 
        %real
        
        subplot(3,2,5)
        plot(X_test,real(Testsignal),'o')
        ylabel('Testfunction')
        xlabel('Non-Equispaced Steps')
        hold on
        plot(t,real(Testsignal2))
        
        %imaginary
        subplot(3,2,6)
        plot(X_test,imag(Testsignal),'o')
        ylabel('Testfunction')
        xlabel('Non-Equispaced Steps')
        hold on
        plot(t,imag(Testsignal2))
        
        %best parameter estimate for b
        kx=[4.48291066436601 - 0.0232247424295838i;5.41782308166276 - 0.0252505740248546i;6.34711835706620 - 0.0271515381332884i;7.27292078959503 - 0.0289422972352744i;8.19639320338978 - 0.0306373246438771i;9.11822807268444 - 0.0322490790633878i;10.0388646289507 - 0.0337879280000863i;10.9585955258497 - 0.0352624496801973i;11.8776236209035 - 0.0366797698925568i;12.7960941782502 - 0.0380458519566265i;13.7141140953456 - 0.0393657289147092i;14.6317638826141 - 0.0406436852143972i;15.5491054014463 - 0.0418833981827764i;16.4661870190230 - 0.0430880485163045i;17.3830471353964 - 0.0442604071291219i;18.2997166540144 - 0.0454029039623454i;19.2162207483640 - 0.0465176829647735i;20.1325801487590 - 0.0476066463992921i;21.0488120952079 - 0.0486714908444128i;21.9649310535897 - 0.0497137366808734i;22.8809492612396 - 0.0507347524250974i;23.7968771477129 - 0.0517357749536686i;24.7127236629481 - 0.0527179264258585i;25.6284965358633 - 0.0536822285329635i;26.5442024800763 - 0.0546296145681435i;27.4598473590062 - 0.0555609397073330i;28.3754363194589 - 0.0564769898124638i;29.2909739005355 - 0.0573784890067363i;30.2064641230508 - 0.0582661062236577i;31.1219105634366 - 0.0591404608938080i;32.0373164152010 - 0.0600021279033990i;32.9526845403354 - 0.0608516419348758i;33.8680175125510 - 0.0616895012807165i;34.7833176538282 - 0.0625161712061832i;35.6985870654663 - 0.0633320869242929i;36.6138276545808 - 0.0641376562360900i];
        b=linspace(1,100,100);
        nn=length(kx);
        predicted=zeros(6,length(b),nn);
        real1=zeros(6,length(b),nn);
        fonc = @(x,y) exp(1i*(y)*x);
        
        
        for bb=1:length(b)
            for kk=1:nn
            [predicted(:,bb,kk),~]=PronyMethod.equispacing(X_test,fonc(X_test,kx(kk)),6,b(bb));
            real1(:,bb,kk)=fonc(X_test,kx(kk));
            end
        end
        argument=predicted-real1;
        handover=sqrt(dot(argument,argument));
        mse=(1/nn)*dot(sqrt(handover),sqrt(handover),3);
        
        figure;
        plot(b,mse)
       end
        
        function TestClass(obj)
            % Set the properties of the class
            obj.NrModes = 12;
            obj.MicSpacing = 0.1;
            obj.Epsilon = 1e-10;
                                   
            %Function to determine whether the prony method is implemented
            %correctly
             %1-MicSpacing; %0.15
            %The wavenumbers have to be large enough that there is a
            %variation accross the microphones
            
            %The wavenumbers are be sorted using the size of the absolute
            %value.
            k = [2+1i,1+2i,5+4i,3,6i, 1.5+1.5i]*5;
            A = [10-8i,5+9i,-5-6i,-5+4i,2+2i,1-1i];
            
            %Determine pressure with the above exponentials at the
            %positions X.
            X = 0:obj.MicSpacing:2.4;
            Pressure = zeros(1,length(X));
            for ii = 1:length(k) %NrModes
                Pressure = Pressure + A(ii) * exp(1i * k(ii) * X);
            end
            
            obj.P = Pressure;
            
            [C, kappa] = PronyMethod.PronyCoefficients(obj.P,obj.MicSpacing,obj.NrModes);
            
%             [C, kappa] = PronyMethod.ESPRIT(obj.P,obj.MicSpacing,obj.NrModes,obj.Epsilon);
             
%             [C, kappa] = PronyMethod.MatrixPencil(obj.P,obj.MicSpacing,obj.NrModes,obj.Epsilon);
%             
%             [C, kappa] = PronyMethod.MyLittleProny(obj.P,obj.MicSpacing,obj.NrModes,obj.Epsilon);

              
               
               
            
            
            %Sort the modes by amplitude
            figure
            subplot(3,2,1)
            plot(X,real(Pressure));
            xlabel('Distance')
            ylabel('Real Pressure')

            subplot(3,2,2)
            plot(X,imag(Pressure));
            xlabel('Distance')
            ylabel('Imag Pressure')
            
            subplot(3,2,3)
            plot(real(kappa),'x')
            hold all
            plot(real(k),'o')
            ylabel('Real Wavenumber')
            xlabel('Mode nr')
            legend('Educed','Input')
            
            subplot(3,2,4)
            plot(imag(kappa),'x')
            hold all
            plot(imag(k),'o')
            ylabel('Imag Wavenumber')
            xlabel('Mode nr')   
            legend('Educed','Input')
                        
            subplot(3,2,5)
            plot(real(C),'x')
            hold all
            plot(real(A),'o')
            ylabel('Real Amplitude')
            xlabel('Mode nr')
            legend('Educed','Input')
            
            subplot(3,2,6)
            plot(imag(C),'x')
            hold all
            plot(imag(A),'o')
            ylabel('Imag Amplitude')
            xlabel('Mode nr')
            legend('Educed','Input')
        end
    end
    
    methods(Static)
        function [Amplitude, kappa] = PronyCoefficients(P,dX,NrModes)
            %Input
            %P is the complex microphone pressure which should be a row
            %vectors
            %dX is the microphone seperation
            %NrModes is the amount of modes that have to be decomposed
            M = length(P);
            if NrModes>M/2
                error('The number of modes can not exceed the half of the measurement points')
            end
            %The hankel matrix is constructed, in order to solve for the
            %exponential functions
            H = toeplitz(P(NrModes:(M-1)),fliplr(P(1:(NrModes))));
            h_M = P(NrModes+1:end).';
            %Obtaining the values of the complex exponential functions
            Lambda = roots([1; -pinv(H)*h_M]);
%             Solution = pinv(transpose(H) * H) * transpose(H) * h_M;
%             Lambda = roots([1; -Solution]);
            
            %Construct the vandermonde matrix to obtain the coefficients
            V = PronyMethod.Vandermonde(M,Lambda);
            %Obtain the coefficients for each of the exponential
            C = V\P.';
            Amplitude = C;
            %The wavenumbers associated with the coefficients
            kappa = -1i*log(Lambda)/dX;
            
            [C,I] = sort(abs(C),'descend');
            Amplitude=Amplitude(I);
            kappa=kappa(I);
        end       
        
        function [Amplitudes, kappa] = MatrixPencil(P, dX, L, epsilon)
            % L = number of modes that are calculated
            % M = number of modes we look at
            
            N = length(P)/2;
            if L > N
                error('The number of modes can not exceed the half of the measurement points')
            elseif L < 3
                error('The number of modes must be superior or equal to 3')
            end
            
            H = PronyMethod.Hankel(2*N-L, L+1, 0, P);
            
            [Q, R, MPi] = qr(H);
            
            M = 1;
            Rabs = abs(R);
            [size_1, size_2] = size(Rabs);
            for ii = 1:min(size_1, size_2)-1
                if Rabs(M+1, M+1) > epsilon * Rabs(1, 1)
                    if M < L
                        M = M+1;
                    end
                end
            end
            
            S = R * transpose(MPi);
                        
            T_s0 = PronyMethod.Matrix_T(S, M, L, 0);
            T_s1 = PronyMethod.Matrix_T(S, M, L, 1);
            
            D = diag(diag(R(1:M, 1:M)));
            
            T_s0 = pinv(D) * T_s0;
            T_s1 = pinv(D) * T_s1;
            
            F = pinv(transpose(T_s0)) * transpose(T_s1);
            
            Z = eig(F);
            
            if length(Z) < M
                disp('There are less eigenvalues than modes')
            elseif length(Z) > M
                disp('There are more eigenvalues than modes')
            end
            
            f_i = log(Z);
            
            V = PronyMethod.Vandermonde(2*N, Z);
            
            C = V \ transpose(P);
            
            Amplitudes = C;
            [C,I] = sort(abs(C),'descend');
            Amplitudes = Amplitudes(I);
            
            
            for ii = 1:length(f_i)
                kappa(ii) = -1i * f_i(ii) / dX;
            end
            
            kappa = kappa(I);
            
        end
        
        function [Amplitudes, kappa] = MyLittleProny(P, dX, L, epsilon)
            % L = number of modes that are calculated
            
            N = length(P)/2;
            if L > N
                error('The number of modes can not exceed the half of the measurement points')
            end
            
            H = PronyMethod.Hankel(2*N-L, L, 0, P);
            h = zeros(2*N-L, 1);
            for ii = 1:length(h)
                h(ii) = -P(ii+L);
            end
            
            H_t = toeplitz(P(L:(2*N-1)),fliplr(P(1:L)));
            h_M = P(L+1:end).';
            
            q = H\h;
            
            Z_0 = roots([1 ; q]);
            
            V_0 = PronyMethod.Vandermonde(2*N, Z_0);
            
            C_0 = V_0\transpose(P);
            
            index = [];
            for ii = 1:length(C_0)
                if abs(C_0(ii)) > epsilon
                    index = [index, ii];
                end
            end
            
            Z_1 = Z_0; %(index);
            
            V_1 = PronyMethod.Vandermonde(2*N, Z_1);
            
            C_1 = C_0; %V_1\transpose(P);
            
            f_i = log(Z_1);
            
            Amplitudes = C_1;
            [C_1,I] = sort(abs(C_1),'descend');
            Amplitudes = Amplitudes(I);
            
            for ii = 1:length(f_i)
                kappa(ii) = -1i * f_i(ii) / dX;
            end
            
            kappa = kappa(I);
            
        end
        
        function [Amplitudes, kappa] = ESPRIT(P, dX, L, epsilon)
            % L = number of modes that are calculated
            
            N = length(P)/2;
            if L > N
                error('The number of modes can not exceed the half of the measurement points')
            end
            
            H = PronyMethod.Hankel(2*N-L, L+1, 0, P);
            
            [~, D, W] = svd(H);
            
            W = ctranspose(W);
            
            sigma = diag(D);
            a = [];
            
            for ii = 2:length(sigma)
                if sigma(ii) < epsilon * sigma(1)
                    a = [a, ii];
                end
            end
            
            if isempty(a) == 1
                M = L;
            else
                M = a(1)-1;
            end
            
            W_rect = W(1:M, 1:L+1);
            
            W_0 = W_rect(1:M, 1:L);
            W_1 = W_rect(1:M, 2:L+1);
            
            F = pinv(transpose(W_0)) * transpose(W_1);
            
            Z = eig(F);
            
            V = PronyMethod.Vandermonde(2*N, Z);
            
            C = V\transpose(P);
            
            f_i = log(Z);
            
            kappa = -1i * f_i / dX;
            
            Amplitudes = C;
            [~,I] = sort(abs(C),'descend');
            Amplitudes = Amplitudes(I);
            
            kappa = kappa(I);
            
        end
        
        function [AX] = ArrowPlotComplexDomain(AX,Data)
            
            for ii = 1:length(Data)-1
                if ii == 1
                quiver(real(Data(ii)),imag(Data(ii)),...
                       real(Data(ii+1)-Data(ii)), ...
                       imag(Data(ii+1)-Data(ii)),0,'o','filled','MarkerFaceColor','r','MaxHeadSize',0.5);
                else
                                    quiver(real(Data(ii)),imag(Data(ii)),...
                       real(Data(ii+1)-Data(ii)), ...
                       imag(Data(ii+1)-Data(ii)),0,'ko','filled','MaxHeadSize',0.5);
                end                
            end
            plot(real(Data(ii+1)),imag(Data(ii+1)),'rx');
                       
        end
        
        function [H] = Hankel(M, N, s, P)
            H = zeros(M, N);
            for ii = 1:M
                for jj = 1:N
                    H(ii, jj) = P((ii - 1) + (jj - 1) + s + 1);
                end
            end
        end
        
        function [V] = Vandermonde(M, Z)
            for ii = 1:M 
                   V(ii, :) = Z(:) .^ (ii - 1);
            end
        end
        
        function [T] = Matrix_T(S, M, L, s)
           T = S(1:M, (1+s):(L+s));
        end
        
        function [S_1] = Matrix_S(S_0, L, s)
            S_1 = S_0(1:length(S_0(:,1)), (1+s):(L+s));
        end
        
       function [equiPressure, equiX] = equispacing(X,P,n,b)
         %This function approximates the complex pressure function with
         %nonequispaced steps and transforms it to a equispaced stepsized
         
         % Input 
         % 1. Parameter: KxN Vector of space locations
         % 2. Parameter: KxN Vector of testfunction
         % 3. Parameter: Value of L: K>=L>=N with N= nonharmonic bandwidth
         % 4. Parameter: b= Variance of the gaussian
        
         % Output
         % 1. Parameter: Kx1 vector equiPressure= approximated Pressure with homogenious
         %               step size 
         % 2. Parameter: Kx1 vector equiX= vector of new homogenious step size
         
           
         % 1) Initialization of parameters for the algorithm
         N=length(X);
         equiX=zeros(N,1);
                 
         % 2) Preconditioning of necessary functions
         %Adapting a new step size to gain a uniform step size
         
         step_size=(max(X)-min(X))/(N-1); 
         equiX(1)=X(1);
                  
         for ii=2:N
              equiX(ii)=equiX(1)+(ii-1)*step_size;             
         end
             
         %Preparing the window function and the stepsize     
         argumentX=zeros(N,n);
         augmentedX=zeros(N,n);
         
         l=-n/2:(n/2)-1;
         
         %argumentX is going to be the new argument used for the
         %windowfunction with eiquispaced input
         for ll=1:n
            argumentX(:,ll)=((X-max(X)/2)/max(X))-((1/n)*l(ll)); 
         end
         
         %augmented X is going to be the new argument used for the
         %windowfunction with eiquispaced input
         for ll=1:n
            augmentedX(:,ll)=((equiX-max(X)/2)/max(X))-((1/n)*l(ll)); 
         end
         
         
         %truncated window function 
         truncatedphi= PronyMethod.twindow(argumentX,b,n);
         truncatedphi2= PronyMethod.twindow(augmentedX,b,n);
         
%        Computation of the support coefficients h
         h=truncatedphi\P; 
         
%        Computation of the pressure with equispaced nodes
         equiPressure = truncatedphi2*h;

        end
        
        function truncatedphi= twindow(argument,b,n)
            % this functions sets up the truncated window function for a
            % gaussian bell window function 
            % inputs: argument: Matrix of different locations in space per node
            %         b: Variance of the gaussian
            %         n: number of parameters h
           
        
            % Initialization of parameters
            truncatedphi=zeros(size(argument,1),n);  
            %truncated phi is a matrix Nxn  containing the different values
            %of the gaussian bell depending on the parameters b,n
            for ii=1:size(argument,1)
                for jj=1:n
                    truncatedphi(ii,jj)=exp(-((n*argument(ii,jj))^2)/b);
                end
            end
         end

    end
    
end
