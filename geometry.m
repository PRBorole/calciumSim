classdef geometry
    methods
        function stencil = stencilMakerCalcium(pt,nx)
        dx = getdx(pt,nx);
        stencil = zeros(nx);
        stencil(1,1) = -2/dx(1)^2;
        stencil(1,2) = 2/dx(1)^2;
        stencil(nx,nx-1) = 2/dx(nx-1)^2;
        stencil(nx,nx) = -2/dx(nx-1)^2;
        for i=2:nx-1
            c1 = 2/(dx(i-1)*(dx(i-1)+dx(i)));
            c2 = -2/(dx(i)*dx(i-1));
            c3 = 2/(dx(i)*(dx(i-1)+dx(i)));
           stencil(i,i-1:i+1) = [c1 c2 c3];
        end
        end

        function dx = getdx(pt,nx)
        dx(nx-1,nx-1) = 0;
        for i = 1:nx-1
            dx();
        end
        end
        
        function createGeometry(nx,type)
            if type == "1D"
                geom = zeros(nx);
                geom(1,1) = 1;
                geom(1,2) = 1;
                geom(nx,nx-1) = 1;
                geom(nx,nx) = 1;
                for i = 2:nx-1
                    geom(i,i-1:i+1) = [1 1 1];
                end
            elseif type == "Branch_1"
                %Here branch is created at center node going outside
                %only 4 points for now
                geom = zeros(nx);
                geom(1,1) = 1;
                geom(1,2) = 1;
                geom(nx,nx-1) = 1;
                geom(nx,nx) = 1;
                for i = 2:nx-1
                    geom(i,i-1:i+1) = [1 1 1];
                end
                
                if mod(nx,2) == 1
                geom
            elseif type == "Branch_2"
                
            elseif type == "Branch_In"
                
            end
            
        end
    end
end

