classdef ZeroNumerical < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    %Todo list: Fix when a zero is present on the boundary. 
    %Option is to find the four neighbouring points, remove them from the
    %search list and replace the Domain with a different subdivision.
    properties
        function_handle
        diff_delta = 1e-6;
        initialZeros = [];
        initialDomain = []
        iterationNumber = 0;
        domainList = [];
        randSeedNewtonSearch = 5;
        found_zeros = [];
        uniqueZeroTol = 1e-3;
        expand_ratio = 1;
        domainSubDivision = 2;
        maxIteration = 500;
        PlotDomains = true;
        NewtonSearchOptions
        FminUncSearchOptions
    end
    
    methods        
        function obj = ZeroNumerical(function_handle,initialDomain)
            if nargin == 0
                %Construct an example object
                z_scale = 0.5;
                obj.function_handle = @(z) (z./z_scale).^50 + (z./z_scale).^12-5*sin(20*(z./z_scale)).*cos(12*(z./z_scale)) - 1 ;
                obj.function_handle = @(z) (z).^50 + 0.9;
                D.Origin = -3.1-3.1i;
                D.Width = 5;
                D.Height = 5;
                obj.initialDomain = D;
                obj.NewtonSearchOptions = optimoptions('fsolve','Display','none','Algorithm','levenberg-marquardt','ScaleProblem','Jacobian','TolX',1e-20,'TolFun',1e-20);
            else
                obj.function_handle = function_handle;
                obj.initialDomain = initialDomain;
            end
        end
        
        function plotZeros(obj)
            if isempty(obj.found_zeros)
                warning('No zeros calculated')
                return
            end
            figure;
            axes(); hold all;
            for ii = 1:length(obj.found_zeros)
                plot(real(obj.found_zeros(ii)),imag(obj.found_zeros(ii)), 'o')
            end
            obj.PlotDomain(obj.initialDomain, gca);
        end
        function [] = PlotFunValuesDomain(obj,Domain,nPoints)
            ReAxis = linspace(real(Domain.Origin), real(Domain.Origin)+Domain.Width,nPoints);
            ImAxis = linspace(imag(Domain.Origin), imag(Domain.Origin)+Domain.Height,nPoints);
            
            [Re,Im] = meshgrid(ReAxis,ImAxis);
            FunVal = obj.function_handle(Re+1i*Im);
            A = log10(min(abs( FunVal(:) )));
            B = log10(max(abs( FunVal(:) )));            
            lev1 = logspace(A,B,10);            
            figure
            contour(Re,Im,abs(FunVal),lev1,'k'); hold all
            contour(Re,Im,abs(FunVal),lev2,'r'); 
            figure;
            contour(Re,Im,abs(FunVal),lev2,'r');
        end
        function [NewDomainList,IntVal] = NewDomain(obj,Domain,Iteration)  
            NrIterations = 4;
            ZerosOnVertices = true;
            qq = 1;
            while ZerosOnVertices && qq < NrIterations
                IntVal = [];
                ZeroOnVertex = [];
                NewDomainList = [];
                ZerosOnVertices = false;                
                if mod(Iteration,2) == 0 
                    NewDomainList = obj.SubDivideSearchDomain(Domain,1+qq,1);
                else
                    NewDomainList = obj.SubDivideSearchDomain(Domain,1,1+qq);
                end
                for ii = 1:length(NewDomainList)
                    [IntVal(ii),ZeroOnVertex(ii)] = obj.IntegrationSearchDomain(NewDomainList(ii));
                    if ZeroOnVertex(ii) ~= 0
                        ZerosOnVertices = true;
                        qq = qq +1 ;
                        break
                    end
                end
            end
            if qq == NrIterations
                error('No subdivision possible which does not contain a zero on one of the vertices')
            end
        end
        function ZerosList = NewtonSearch(obj,Domain)                        
            %Perform a Newton's search in the domain
            nn = 1;
            foundzero = [];
            for ii = 1:obj.randSeedNewtonSearch
                z0 = obj.StartValueSearchDomain(Domain);
                if imag(z0) == 0
                    %If the initial value has no imaginary
                    %component, add a small component to search
                    %for complex zeros.
                    z0 = z0 + 1e3*1i;
                end
                fun = obj.function_handle;
                
                %Scale the function to have a better convergence of the
                %solver.
                funscaled = @(z) 1/fun(z0)*fun(z);
                [z,fval,exitflag,output] = fsolve(funscaled,z0,obj.NewtonSearchOptions);
                if exitflag > 0
                    % Check if the found zero is within the domain
                    if obj.InDomain(Domain, z)
                        foundzero(nn) = z;
                        nn = nn+1;
                    end
                end
            end 
            [UniqueZeros,Ia,Ib] = uniquetol(abs(foundzero),obj.uniqueZeroTol);
            ZerosList = foundzero(Ia);
        end
        
        function searchZeros(obj)
            %First determine the winding number in the initial domain
            [obj.initialZeros,ZeroOnVertex] = obj.IntegrationSearchDomain(obj.initialDomain);
            if ZeroOnVertex ~= 0
                error('The boundary of the initial domain contains a zero')
            end
             if round(obj.initialZeros) == 0
                 error('There are no zeros in the initial domain')
             end
             
            %Subdivide the initial domain, add the subdomains to the domain list and calculate the winding number.
            %When the winding number is 1, find the zero using Newton's
            %method and remove the subdomain from the list if a zero is
            %found within the search domain.
            %Iterate until there is no subdomain on the list.
            
            Iteration = 0;
            obj.found_zeros = [];
            ax = [];
            DomainList = obj.initialDomain;
            while ~isempty(DomainList) || Iteration < obj.maxIteration
                IntVal = [];
                NewDomainList = [];
                ZerosList = [];
                
                %Subdivide contour and Integrate around each sub contour
                disp('Subdividing domain and integrating around domain')   
                for dd = 1:length(DomainList)                                     
                    [NewDomains,NewIntVal] = obj.NewDomain(DomainList(dd),Iteration);
                    NewDomainList = [NewDomainList, NewDomains];
                    IntVal = [IntVal,NewIntVal];
                end 
                
                if obj.PlotDomains
                    cla(ax)
                    for dd = 1:length(NewDomainList)
                    ax = obj.PlotDomain(NewDomainList(dd), ax);  
                    ax = obj.PlotDomainVal(NewDomainList(dd), round(IntVal(dd)), ax);  
                    end
                    set(ax,'Xlim',[real(obj.initialDomain.Origin), real(obj.initialDomain.Origin) + obj.initialDomain.Width],...
                           'Ylim',[imag(obj.initialDomain.Origin), imag(obj.initialDomain.Origin) + obj.initialDomain.Height]);
                    drawnow;
                end
                ZerosList = [];
                IntVal = round(IntVal);
                RemoveDomain = [];
                for dd = 1:length(NewDomainList)                      
                    if IntVal(dd) == 0    
                            RemoveDomain = [RemoveDomain,dd];
                    end
                    if IntVal(dd) == 1                        
                        NewZeroes = obj.NewtonSearch(NewDomainList(dd));
                        if length(NewZeroes) == 1
                            RemoveDomain = [RemoveDomain,dd];
                            ZerosList = [ZerosList,NewZeroes];
                        end
                        if length(NewZeroes) > 1
                            warning('No unique zero found')
                        end
                    end
                end    
                NewDomainList(RemoveDomain) = [];
                Iteration = Iteration + 1;                
                DomainList = NewDomainList;
                
                %Update the class properties
                obj.found_zeros = [obj.found_zeros,ZerosList];                
                obj.iterationNumber = Iteration;
                obj.domainList = NewDomainList;
            end        
            
        end
        
       
        function [ax] = PlotDomain(obj,Domain, ax)
            if isempty(ax)
                figure;
                ax = axes(gcf); hold on;
            end
            
            IntPoints = obj.DomainToIntegrationPoints(Domain);
            for ii = 1:length(IntPoints)-1
                line([real(IntPoints(ii)) real(IntPoints(ii+1))],[imag(IntPoints(ii)) imag(IntPoints(ii+1))])
            end
        end
        
        function [ax] = PlotDomainVal(obj,Domain, Val, ax)
            if isempty(ax)
                figure;
                ax = axes(gcf); hold on;
            end
            text(real(Domain.Origin) + 0.5*Domain.Width, imag(Domain.Origin) + 0.5*Domain.Height,num2str(Val));         
        end
        
        function val = NumIntegrand(obj, z)
            val = 1/(2*pi*1i)*(obj.function_handle(z+obj.diff_delta)-obj.function_handle(z-obj.diff_delta))./(2*obj.diff_delta)./...
                obj.function_handle(z);
        end
        
        function z0 = StartValueSearchDomain(obj,Domain)
            z0 = Domain.Origin + (0.3*(rand(1)-0.5)+0.5)*Domain.Width + 1i*(0.3*(rand(1)-0.5)+0.5)*Domain.Height;
        end
        
        function bool = InDomain(obj,Domain,z)
            bool = false;
            %Check wether the complex number is the domain bounded by the
            %Domain.
            if (real(z) > real(Domain.Origin)) && (real(z) < real(Domain.Origin) +Domain.Width) ...
                    && (imag(z) > imag(Domain.Origin)) && (imag(z) < imag(Domain.Origin) +Domain.Height)
                bool = true;
            end
        end        

        function SubDividedDomain = SubDivideSearchDomain(obj,Domain,nRe,nIm)
            %Function to subdivide the search Domain in nRe*nIm sections
            for ii = 1:nRe
                for jj = 1:nIm
                    SubDividedDomain(nIm*(ii-1) + jj).Origin = Domain.Origin + (ii-1)*Domain.Width/nRe + 1i*(jj-1)*Domain.Height/nIm;
                    SubDividedDomain(nIm*(ii-1) + jj).Width = Domain.Width/nRe;
                    SubDividedDomain(nIm*(ii-1) + jj).Height = Domain.Height/nIm;
                end
            end
        end
        
        function IntPoints = DomainToIntegrationPoints(obj,Domain)
            %This function returns the points in the complex domain
            %that build up the Domain, based on an origin on the lower
            %left corner, the width (real part) and height (imaginary
            %part) of the rectangle.
            %The points are ordered counter clock wise from the origin;
            %Domain(1) = Origin
            %Domain(2) = Width
            %Domain(3) = Height
            IntPoints(1) = Domain.Origin;
            IntPoints(2) = Domain.Origin + Domain.Width;
            IntPoints(3) = Domain.Origin + Domain.Width + 1i*Domain.Height;
            IntPoints(4) = Domain.Origin + 1i*Domain.Height;
            IntPoints(5) = Domain.Origin;
        end
        
        function [IntVal,ZeroOnVertex] = IntegrationSearchDomain(obj,Domain)
            IntVal = 0;
            ZeroOnVertex = 0;
            IntPoints = DomainToIntegrationPoints(obj,Domain);
            warning(''); %Clear the warning
            for ii = 1:length(IntPoints)-1
                    IntVal = IntVal + obj.NumIntegration(IntPoints(ii),IntPoints(ii+1));
                    %Check for the last warning, if it is related to the
                    %integral, the chance is large that the function f'/f
                    %is not defined
                    %Return the side of the rectangle where the warning was
                    %issued.
                    [warnmsg, msgid] = lastwarn;
                    if strcmp(msgid,'MATLAB:integral:MaxIntervalCountReached')
                        ZeroOnVertex = ii;
                        IntVal = Inf;
                        break
                    end
            end
        end
        
        function Int = NumIntegration(obj,a,b)
            fun = @(z) obj.NumIntegrand(z);
            Int = integral(fun,a,b);            
        end
        
        
    end
end


