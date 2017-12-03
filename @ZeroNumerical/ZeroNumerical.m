classdef ZeroNumerical < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    %Todo list: Fix when a zero is present on the boundary. 
    %Option is to find the four neighbouring points, remove them from the
    %search list and replace the square with a different subdivision.
    properties
        function_handle
        diff_delta = 1e-6;
        initialZeros = [];
        initialDomain = []
        iterationNumber = 0;
        domainList = [];
        randSeedNewtonSearch = 5;
        zeros = [];
        uniqueZeroTol = 1e-6;
        expand_ratio = 1;
        domainSubDivision = 2;
        maxIteration = 500;
        PlotDomains = true;
    end
    
    methods        
        function obj = ZeroNumerical(function_handle,initialDomain)
            if nargin == 0
                %Construct an example object
                obj.function_handle = @(z) z.^50 + z.^12-5*sin(20*z).*cos(12*z) - 1 ;
                D.Origin = -5.0-4.0i;
                D.Width = 10;
                D.Height = 8;
                obj.initialDomain = D;
            else
                obj.function_handle = function_handle;
                obj.initialDomain = initialDomain;
            end
        end
        
        function plotZeros(obj)
            if isempty(obj.zeros)
                warning('No zeros calculated')
                return
            end
            figure;
            axes(); hold all;
            for ii = 1:length(obj.zeros)
                plot(real(obj.zeros(ii)),imag(obj.zeros(ii)), 'o')
            end
        end
                
            
            
        function searchZeros(obj)
            %First determine the winding number in the initial domain
            [obj.initialZeros,ZeroOnVertex] = obj.IntegrationSearchArea(obj.initialDomain);
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
            
            SubDividedDomain = obj.SubDivideSearchArea(obj.initialDomain,2,1);
            obj.domainList = SubDividedDomain;
            Iteration = 0;
            obj.zeros = [];
            ax = [];

            
           
            while ~isempty(obj.domainList) || Iteration < obj.maxIteration
                DomainList = obj.domainList;
                NewDomainList = [];
                ZerosList = [];
                
                %Subdivide contour and Integrate around each sub contour                            
                for dd = 1:length(DomainList)
                    disp('Integrating around domain')
                    [IntVal(dd),ZeroOnVertex(dd)] = obj.IntegrationSearchArea(DomainList(dd));
                    
                    if ZeroOnVertex(dd) ~= 0
                        
                    end
                end 
                
                if obj.PlotDomains
                    for dd = 1:length(DomainList)
                    ax = obj.PlotDomain(DomainList(dd), ax);  
                    %ax = obj.PlotDomainVal(obj.domainList(dd), IntVal(dd), ax);  
                    end
                    set(ax,'Xlim',[real(obj.initialDomain.Origin), real(obj.initialDomain.Origin) + obj.initialDomain.Width],...
                           'Ylim',[imag(obj.initialDomain.Origin), imag(obj.initialDomain.Origin) + obj.initialDomain.Height]);
                    drawnow;
                end

                IntVal = round(IntVal);
                for dd = 1:length(DomainList)    
                    if IntVal(dd) == 1
                        %Perform a Newton's search in the domain
                        nn = 1;
                        foundzero = [];
                        for ii = 1:obj.randSeedNewtonSearch
                            z0 = obj.StartValueSearchArea(DomainList(dd));
                            if imag(z0) == 0
                                %If the initial value has no imaginary
                                %component, add a small component to search
                                %for complex zeros.
                                z0 = z0 + 1e3*1i;
                            end
                            disp('Searching for zero')
                            z = fsolve(obj.function_handle,z0);                            
                            % Check if the found zero is within the domain
                            if obj.InDomain(DomainList(dd), z)
                                foundzero(nn) = z;
                                nn = nn+1;
                            end
                        end 
                        if length(uniquetol(abs(foundzero),obj.uniqueZeroTol)) >= 1
                            if length(uniquetol(abs(foundzero),obj.uniqueZeroTol)) > 1
                                warning('No unique zero found')
                            end
                            [UniqueZeros,Ia,Ib] = uniquetol(abs(foundzero),obj.uniqueZeroTol);
                            ZerosList = [ZerosList,foundzero(Ia) ];
                        elseif isempty(uniquetol(abs(foundzero),obj.uniqueZeroTol))
                            %If no zeros have been found within the
                            %domain, subdivide the domain.
                            warning('No zeros found within the domain, subdividing domain')
                            NewDomainList = [NewDomainList, obj.SubDivideSearchArea(DomainList(dd),mod(Iteration,2)+1,mod(Iteration+1,2)+1)];
                        end
                    elseif IntVal(dd) == Inf
                        %If IntVal returns infinity, a zero is present on
                        %the border of the integration. Expand the
                        %Integration area and add it to the list for the
                        %next iteration                        
                        disp('Expanding integration area')
                        ExpandedArea = obj.ExpandArea(DomainList(dd),obj.expand_ratio);
                        NewDomainList = [NewDomainList,ExpandedArea];
                    elseif IntVal(dd) > 1
                        %Subdivide the domain
                        disp('Subdividing integration area')
                        NewDomainList = [NewDomainList, obj.SubDivideSearchArea(DomainList(dd),mod(Iteration,2)+1,mod(Iteration+1,2)+1)];
                        
                    elseif IntVal(dd) == -1
                        disp('Expanding integration area')
                        ExpandedArea = obj.ExpandArea(DomainList,obj.expand_ratio);
                        NewDomainList = [NewDomainList,ExpandedArea];
                        continue
                    end
                end
                
                Iteration = Iteration + 1;
                %Update the class properties
                obj.zeros = [obj.zeros,ZerosList];                
                obj.iterationNumber = Iteration;
                obj.domainList = NewDomainList;
            end
            
            
        end
        function [NewDomainList, NewIntVal] = ReDivisionDomain(obj, DomainList, IntVal, VertexError)
            %Function to redivision the integration area
            
            %Find the domains that have an error
            [I] = find(VertexError ~= 0);
            %Obtain the origins of their parent rectangle. Before the
            %division there was no error
            for ii=1:length(I)
                Parent(ii) = DomainList(I(ii)).ParentOrigin;
            end
            %Remove the rectangles with the same parent rectangle from the
            %domain list.
            DomainList(I) = [];
            %Obtain the unique parent rectangles
            [ParentOrigin] = unique(ParentOrigin);
            %Redivide the domain, calculate the integral and return the data if it is correct 
            
            for ii = 1:length(ParentOrigin)
                
                [NewDomainList, obj.SubDivideSearchArea(DomainList(dd),mod(Iteration,2)+1,mod(Iteration+1,2)+1)];
            end
                
            for dd = 1:length(obj.domainList)
                if VertexError(dd) ~= 0
                    %Determine on which vertice of the rectangle the
                    %zero is placed and find the origins of the side
                    %rectangle.
                    %There always has to be one
                    if VertexError(dd) == 1
                        NeighBouringOrigin = DomainList(dd).Origin - 1i*DomainList(dd).Height;
                        NewDomain.Origin = NeighBouringOrigin;
                    end
                    if VertexError(dd) == 2
                        NeighBouringOrigin = DomainList(dd).Origin + DomainList(dd).Width;
                        NewDomain.Origin = DomainList(dd).Origin;
                    end
                    if VertexError(dd) == 3
                        NeighBouringOrigin = DomainList(dd).Origin + 1i*DomainList(dd).Height;
                        NewDomain.Origin = DomainList(dd).Origin;
                    end
                    if VertexError(dd) == 4
                        NeighBouringOrigin = DomainList(dd).Origin - DomainList(dd).Width;
                        NewDomain.Origin = NeighBouringOrigin;
                    end
                    
                    for nn = 1:length(DomainList)
                        if DomainList(nn).Origin == NeighBouringOrigin
                            
                            break
                        end
                    end
                    
                    NewDomain.Origin = NeighBouringOrigin;
                    %Add the two domains together, and subdivide them
                    if VertexError(dd) == 1 || VertexError(dd) == 3
                        NewDomain.Height = 2*DomainList(dd).Height;
                        NewDomain.Width = 1*DomainList(dd).Width;
                    end
                    if VertexError(dd) == 2 || VertexError(dd) == 4
                        NewDomain.Height = 1*DomainList(dd).Height;
                        NewDomain.Width = 2*DomainList(dd).Width;
                    end
                    qq = 1;
                    NrMaxSubDivisions = 3;
                    Success = false;
                    
                    while ~Success && qq < NrMaxSubDivisions
                        IntVal_New = [];
                        VertexError_New = [];
                        Success = true;
                        SubDividedNewDomain = obj.SubDivideSearchArea(obj.initialDomain,obj.domainSubDivision + qq);
                        %                             for jj = 1:length(SubDividedNewDomain)
                        %                                 ax = obj.PlotDomain(SubDividedNewDomain(jj), ax);
                        %                             end
                        for jj = 1:length(SubDividedNewDomain)
                            [IntVal_New(jj),VertexError_New(jj)] = obj.IntegrationSearchArea(SubDividedNewDomain(jj));
                            if VertexError_New(jj) ~= 0
                                %Restart the loop and add one more to the
                                %subdivision of the domain
                                qq = qq +1;
                                Success = false;
                                break
                            end
                        end
                    end
                    if ~Success
                        error('Integration along the edges does not work')
                    else
                        %Remove the old domains and replace them with
                        %the new data
                        obj.domainList([dd,nn]) = [];
                        IntVal([dd,nn]) = [];
                        IntVal = [IntVal, IntVal_New];
                        obj.domainList = [obj.domainList, SubDividedNewDomain];
                    end
                end
            end
        
        function [ax] = PlotDomain(obj,Square, ax)
            if isempty(ax)
                figure;
                ax = axes(gcf); hold on;
            end
            
            IntPoints = obj.SquareToIntegrationPoints(Square);
            for ii = 1:length(IntPoints)-1
                line([real(IntPoints(ii)) real(IntPoints(ii+1))],[imag(IntPoints(ii)) imag(IntPoints(ii+1))])
            end
        end
        
        function [ax] = PlotDomainVal(obj,Square, Val, ax)
            if isempty(ax)
                figure;
                ax = axes(gcf); hold on;
            end
            text(real(Square.Origin) + 0.5*Square.Width, imag(Square.Origin) + 0.5*Square.Height,num2str(Val));         
        end
        
        function val = NumIntegrand(obj, z)
            val = 1/(2*pi*1i)*(obj.function_handle(z+obj.diff_delta)-obj.function_handle(z-obj.diff_delta))./(2*obj.diff_delta)./...
                obj.function_handle(z);
        end
    
        function NewSquare = ExpandArea(obj, Square,ExpansionRatio)
            NewSquare.Origin = Square.Origin-0.5*Square.Width*ExpansionRatio - 1i*0.5*Square.Height*ExpansionRatio;
            NewSquare.Width = Square.Width*ExpansionRatio;
            NewSquare.Height = Square.Height*ExpansionRatio;
        end
        
        function z0 = StartValueSearchArea(obj,Square)
            z0 = Square.Origin + rand(1)*Square.Width + 1i*rand(1)*Square.Height;
        end
        
        function bool = InDomain(obj,Square,z)
            bool = false;
            %Check wether the complex number is the domain bounded by the
            %Square.
            if (real(z) > real(Square.Origin)) && (real(z) < real(Square.Origin) +Square.Width) ...
                    && (imag(z) > imag(Square.Origin)) && (imag(z) < imag(Square.Origin) +Square.Height)
                bool = true;
            end
        end
        

        function SubDividedSquare = SubDivideSearchArea(obj,Square,nRe,nIm)
            %Function to subdivide the search area in nRe*nIm sections
            for ii = 1:nRe
                for jj = 1:nIm
                    SubDividedSquare(nIm*(ii-1) + jj).Origin = Square.Origin + (ii-1)*Square.Width/nRe + 1i*(jj-1)*Square.Height/nIm;
                    SubDividedSquare(nIm*(ii-1) + jj).Width = Square.Width/nRe;
                    SubDividedSquare(nIm*(ii-1) + jj).Height = Square.Height/nIm;
                    SubDividedSquare.ParentOrigin = Square.Origin;
                    SubDividedSquare.ParentWidth = Square.Width;
                    SubDividedSquare.ParentHeight = Square.Height;
                end
            end
        end
        
        function IntPoints = SquareToIntegrationPoints(obj,Square)
            %This function returns the points in the complex domain
            %that build up the square, based on an origin on the lower
            %left corner, the width (real part) and height (imaginary
            %part) of the rectangle.
            %The points are ordered counter clock wise from the origin;
            %Square(1) = Origin
            %Square(2) = Width
            %Square(3) = Height
            IntPoints(1) = Square.Origin;
            IntPoints(2) = Square.Origin + Square.Width;
            IntPoints(3) = Square.Origin + Square.Width + 1i*Square.Height;
            IntPoints(4) = Square.Origin + 1i*Square.Height;
            IntPoints(5) = Square.Origin;
        end
        
        function [IntVal,ZeroOnVertex] = IntegrationSearchArea(obj,Square)
            IntVal = 0;
            ZeroOnVertex = 0;
            IntPoints = SquareToIntegrationPoints(obj,Square);
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


