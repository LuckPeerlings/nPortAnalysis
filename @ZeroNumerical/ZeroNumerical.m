classdef ZeroNumerical < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties 
        %The handle to the function that has to be evaluated
        function_handle        
        %The winding number after the integration of the initial domain
        initialZeros = [];
        %The initial domain, containing three fields, origin (complex), width (real) and height (real) of the domain to be searched. Can be an array to create domains excluding for example poles.
        initialDomain = [];
        %How many complete iterations have been performed
        iterationNumber;
        %The domains that still have zeros in them
        domainList = [];
        %Number of random points within each domain when searching for the zeros;
        randSeedNewtonSearch = 5;
        %List of the zeros that have been found, not ordered
        found_zeros = [];
        %Tolerance when defining an unique zero when comparing the obtained values from the zero search
        uniqueZeroTol = 1e-3;
        %How often a domain will redivided in order to have no zeros on the
        %boundaries.
        subDivisionIterations = 5;
        
        %Show a figure of the domains at each iteration step;
        PlotDomains = true; 
        %The options to use for the newton search method (using fsolve)
        NewtonSearchOptions = optimoptions('fsolve','Display','none','Algorithm','levenberg-marquardt','ScaleProblem','Jacobian','TolX',1e-20,'TolFun',1e-20);
        %Maximum number of complete iterations, subdividing the boundaries and searching for zeros 
        maxIteration = 100; 
        %Step taken to determine the derivative of the functions, based on the central difference method
        diff_delta = 1e-6;  
    end
    
    methods        
        function obj = ZeroNumerical(function_handle,initialDomain)
            if nargin == 0
                %Construct an example object
                obj.function_handle = @(z) z.^50 + z.^12-5*sin(20*z).*cos(12*z) - 1 ;
                D.Origin = -6.1-6.1i;
                D.Width = 10;
                D.Height = 10;
                obj.initialDomain = D;                
            else
                obj.function_handle = function_handle;
                obj.initialDomain = initialDomain;
            end
        end
        
        function searchZeros(obj)
            %Main method that finds the zeros in the domain.
            
            %First determine the winding number in the initial domain
            %obj.initialDomain = obj.SubDivideDomain(obj.initialDomain,100,100);
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
            while ~isempty(DomainList) && Iteration < obj.maxIteration
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
                            disp(NewZeroes)
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
        
        function [] = PlotFunValuesDomain(obj,Domain,nPoints)
            %Method to plot the absolute value of the evaluating function
            %in the desired domain.
            %For debugging purposes.
            if isempty(Domain)
                Domain = obj.initialDomain;
            end
            if isempty(nPoints)
                nPoints = 1000;
            end
            ReAxis = linspace(real(Domain.Origin), real(Domain.Origin)+Domain.Width,nPoints);
            ImAxis = linspace(imag(Domain.Origin), imag(Domain.Origin)+Domain.Height,nPoints);
            
            [Re,Im] = meshgrid(ReAxis,ImAxis);
            FunVal = obj.function_handle(Re+1i*Im);
            A = log10(min(abs( FunVal(:) )));
            B = log10(max(abs( FunVal(:) )));            
            lev1 = logspace(A,B,30);            
            figure
            contourf(Re,Im,abs(FunVal),lev1);
            
            figure
            PhasePlot(Re+1i*Im,FunVal,'j');
        end
        
        function [NewDomainList,IntVal] = NewDomain(obj,Domain,Iteration) 
            %Subdivide the input Domain, such that the new boundaries do
            %not have any zeros on them.
            %If there is a zero, the domain is redivided with one extra
            %division
            
            ZerosOnVertices = true;
            qq = 1;
            while ZerosOnVertices && qq < obj.subDivisionIterations
                IntVal = [];
                ZeroOnVertex = [];
                ZerosOnVertices = false;                
                if mod(Iteration,2) == 0 
                    NewDomainList = obj.SubDivideDomain(Domain,1+qq,1);
                else
                    NewDomainList = obj.SubDivideDomain(Domain,1,1+qq);
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
            if qq == obj.subDivisionIterations
                error('No subdivision possible for which a zero is not present on one of the vertices')
            end
        end
        
        function ZerosList = NewtonSearch(obj,Domain)                        
            %Perform a Newton's search in the domain to find the zero's
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
                [z,~,exitflag,~] = fsolve(funscaled,z0,obj.NewtonSearchOptions);
                if exitflag > 0
                    % Check if the found zero is within the domain
                    if obj.InDomain(Domain, z)
                        foundzero(nn) = z;
                        nn = nn+1;
                    end
                end
            end 
            %Determine if the found zeros are unique and only return the
            %unique zeros
            [~,Ia,~] = uniquetol(abs(foundzero),obj.uniqueZeroTol);
            ZerosList = foundzero(Ia);
        end
        
        
        
        function plotZeros(obj)
            %Plot the found zeros in the complex domain, together with the
            %bounding box of the search domain
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
            xlabel('Real axis')
            ylabel('Imaginary axis')
            title('Plot of the zeros in the domain')
        end
        
        function [ax] = PlotDomain(obj,Domain, ax)
            %Method to function the plot the domain in the axes given by
            %the axes handle "ax"
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
            %Method to plot the winding number in the axes given by "ax".
            %The position is in the middle of the domain.
            if isempty(ax)
                figure;
                ax = axes(gcf); hold on;
            end
            text(real(Domain.Origin) + 0.5*Domain.Width, imag(Domain.Origin) + 0.5*Domain.Height,num2str(Val));         
        end
        
        function val = NumIntegrand(obj, z)
            %The integrand f'(z)/f(z) that is evaluated to obtain the
            %number of zeros in the domain
            %The derivative is obtained using the central difference method
            val = 1/(2*pi*1i)*(obj.function_handle(z+obj.diff_delta)-obj.function_handle(z-obj.diff_delta))./(2*obj.diff_delta)./...
                obj.function_handle(z);
            
            %val = 1/(2*pi*1i).*obj.function_handle(z);
        end
        
        function z0 = StartValueSearchDomain(obj,Domain)
            %Obtain a random starting value for the Newton's search that is
            %within the domain.
            z0 = Domain.Origin + (0.3*(rand(1)-0.5)+0.5)*Domain.Width + 1i*(0.3*(rand(1)-0.5)+0.5)*Domain.Height;
        end
        
        function bool = InDomain(obj,Domain,z)
            %Method to check wether the complex number is the domain bounded by the
            %Domain.
            bool = false;
            if (real(z) > real(Domain.Origin)) && (real(z) < real(Domain.Origin) +Domain.Width) ...
                    && (imag(z) > imag(Domain.Origin)) && (imag(z) < imag(Domain.Origin) +Domain.Height)
                bool = true;
            end
        end        

        function SubDividedDomain = SubDivideDomain(obj,Domain,nRe,nIm)
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
            %that builds up the Domain, based on an origin on the lower
            %left corner, the width (real part) and height (imaginary
            %part) of the rectangle.
            %The points are ordered counter clock wise from the origin;
            %The fifth point is added to easily determine the integrals
            %when using a loop.
            IntPoints(1) = Domain.Origin;
            IntPoints(2) = Domain.Origin + Domain.Width;
            IntPoints(3) = Domain.Origin + Domain.Width + 1i*Domain.Height;
            IntPoints(4) = Domain.Origin + 1i*Domain.Height;
            IntPoints(5) = Domain.Origin;
        end
        
        function [IntVal,ZeroOnVertex] = IntegrationSearchDomain(obj,Domain)
            %Perform the integration over the boundaries of the domain.
            IntVal = 0;
            ZeroOnVertex = 0;
            IntPoints = DomainToIntegrationPoints(obj,Domain);
            lastwarn(''); %Clear the warning
            for ii = 1:length(IntPoints)-1
                    IntVal = IntVal + obj.NumIntegration(IntPoints(ii),IntPoints(ii+1));
                    %Check for the last warning, if it is related to the
                    %integral, the chance is large that the function f'/f
                    %is not defined
                    %Return the side of the rectangle where the warning was
                    %issued and stop the loop as the obtained value of the
                    %winding number will not be correct.
                    [~, msgid] = lastwarn;
                    if strcmp(msgid,'MATLAB:integral:MaxIntervalCountReached')
                        ZeroOnVertex = ii;
                        IntVal = Inf;
                        break
                    end
            end
        end
        
        function Int = NumIntegration(obj,a,b)
            %Method that evaluates the numerical integral, defined by the
            %integrand.
            fun = @(z) obj.NumIntegrand(z);
            Int = integral(fun,a,b);            
        end
        
        
    end
end


