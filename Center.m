% Center class.
classdef Center < handle
    properties(SetAccess = private)
        x
        y
        pres
    end
    
    methods
        %% Initialize variables stored at cell center.
        function obj = Center(domain)
            % set the grid
            obj.x = linspace(-0.5, domain.nx+2-1.5, domain.nx+2)*domain.dx;
            obj.y = linspace(-0.5, domain.ny+2-1.5, domain.ny+2)*domain.dy;
            % obj.pressure
            obj.pres = zeros(domain.nx+2, domain.ny+2);
        end
        
        
        %% Calculate the pressure field.
        function solvePressure(obj, domain, param, fluid, face)
            % initialize variables
            [temp1, temp2] = deal(zeros(domain.nx+2, domain.ny+2));
            % calculate source term and the coefficient for obj.pressure
            rhoTemp = fluid.rho;
            largeNum = 1000;
            rhoTemp(1:domain.nx+2,1) = largeNum;
            rhoTemp(1:domain.nx+2,domain.ny+2) = largeNum;
            rhoTemp(1,1:domain.ny+2) = largeNum;
            rhoTemp(domain.nx+2,1:domain.ny+2) = largeNum;
            for i=2:domain.nx+1
                for j=2:domain.ny+1
                    temp1(i,j) = (0.5/param.dt)* ...
                        ((face.uTemp(i,j)-face.uTemp(i-1,j  ))/domain.dx+ ...
                         (face.vTemp(i,j)-face.vTemp(i  ,j-1))/domain.dy);
                    temp2(i,j) = 1.0/((1./domain.dx)* ...
                        (1./(domain.dx*(rhoTemp(i+1,j  )+rhoTemp(i,j)))+...
                         1./(domain.dx*(rhoTemp(i-1,j  )+rhoTemp(i,j))))+ ...
                        (1./domain.dy)* ...
                        (1./(domain.dy*(rhoTemp(i  ,j+1)+rhoTemp(i,j)))+...
                         1./(domain.dy*(rhoTemp(i  ,j-1)+rhoTemp(i,j)))));
                end
            end

            % construct the obj.pressure field using SOR
            for it=1:param.maxIter
                oldPres = obj.pres;
                for i=2:domain.nx+1
                    for j=2:domain.ny+1
                        obj.pres(i,j) = (1.0-param.beta)*obj.pres(i,j)+ ...
                            param.beta*temp2(i,j)* ...
                            ((1./domain.dx)* ...
                            (obj.pres(i+1,j)/(domain.dx* ...
                            (rhoTemp(i+1,j)+ rhoTemp(i,j)))+ ...
                             obj.pres(i-1,j)/(domain.dx* ...
                            (rhoTemp(i-1,j)+ rhoTemp(i,j))))+ ...
                            (1./domain.dy)* ...
                            (obj.pres(i,j+1)/(domain.dy* ...
                            (rhoTemp(i,j+1)+ rhoTemp(i,j)))+...
                             obj.pres(i,j-1)/(domain.dy* ...
                            (rhoTemp(i,j-1)+ rhoTemp(i,j))))-temp1(i,j));
                    end
                end
                if max(max(abs(oldPres-obj.pres))) < param.maxErr
                    break
                end
            end
        end
        
    end
end
