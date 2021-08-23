function [matrixx,matrixy] = q4w6newscheme(numberofpointsinrho,numberofpointsintime,tau)

        % We look to implement Curve shortening flow using finite elements where we
        % have parametrised some intial curve as vectors of (x,y). This will be for
        % the postive semi unit circle attached orthonally to the x axis.

                              %   Initial conditions
                                
        % For \rho in [0,1] we can parameterise the positive semi unit circle as
        % [cos\pi\rho, sin\pi\rho]^{T} depending on how many angle segments we choose
        % to implement. We will include a forcing term (f) so that it
        % can be used such that forcing is included.

        drho = (1-0)/numberofpointsinrho;
        dt = (tau-0)/numberofpointsintime;
 
        pointrho = zeros(1,numberofpointsinrho+1);
        pointtime = zeros(1,numberofpointsintime+1);
        
            for i = 1:numberofpointsinrho+1
                pointrho(i) = (i-1)*drho;
            end
        
            for j=1:numberofpointsintime+1
                pointtime(j) = (j-1)*dt;
            end

        % First create matrices for the x and y component of each vector

        matrixx = zeros(numberofpointsintime+1,numberofpointsinrho+1); % x component 
        matrixy = zeros(numberofpointsintime+1,numberofpointsinrho+1); % y component
        
        % Create column vectors for RHS linear coefficients of prev
        % iteration
      
        

        % Set up initial conditions into the vector matrices above with
        % matrixx_0 = cos(pi*rho), matrixy_0 = sin(pi*rho).
        
        for i = 1:numberofpointsinrho+1
            matrixx(1,i) = 4*cos(pi*pointrho(i)); %*4
            matrixy(1,i) = 4*sin(pi*pointrho(i)); %*4
        end
         
        % We now look to build all of our matrix systems. Since these
        % change on each loop, we now start our iterated loop sequence over
        % time.
        
        j=1;
        
           %F(x,y,t) = cos(2pi*t)*y - sin(2pi*t)*x
           
          
           %FX = @(x,y,t) 0;
           %FY = @(x,y,t) 1;
           %FT = @(x,y,t) 0;
           
           FX = @(x,y,t)       -sin(2*pi*t);
           FY = @(x,y,t)       cos(2*pi*t);
           FT = @(x,y,t)       -2*pi*( (y*(sin(2*pi*t))) + (x*(cos(2*pi*t))));
           
           
           %FX = @(x,y,t)       0;
           %FY = @(x,y,t)       1;
           %FT = @(x,y,t)       -1;
           
           
           for j=1:numberofpointsintime
                
            %-------------------------------------------------------------%
            
        
                    % Forcing term to start we will say
                   
                     p = @(j) 0;
                    
                    % The below was for week 5 - no forcing in week 6.
        
                    % If we want to use the condition p = 2/ r(t) where r(t) denotes
                    % the radius at time t , then we need to insert this into the loop so this
                    % will change every iteration.
                    
                    % How to calculate the radius, we will just use the first point (rho = 0) since the forcing is constant and
                    % we are starting with the unit circle where the curvature is constant.
                   
            
                    %Start with building the left hand side matrix to be
                    %inverted
                    
                    matrixtobeinverted = zeros((2*numberofpointsinrho)+2,(2*numberofpointsinrho)+2);
                    RHSIvec = zeros(2*(numberofpointsinrho+1),1);
                    RHSforcingvec = zeros(2*(numberofpointsinrho+1),1);
                    
                    % Help simplify the below vectors and re write them as
                    % written
                    
                    x1 = matrixx(j,1);
                    x2 = matrixx(j,2);
                    y1 = matrixy(j,1);
                    y2 = matrixy(j,2);
                    
                    xend = matrixx(j,numberofpointsinrho+1);
                    xendmin1 = matrixx(j,numberofpointsinrho);
                    yend = matrixy(j,numberofpointsinrho+1);
                    yendmin1 = matrixy(j,numberofpointsinrho);
                    
                    t = pointtime(j);
                    
                    FX1 = FX(x1,y1,t); 
                    FY1 = FY(x1,y1,t); 
                    FT1 = FT(x1,y1,t);
                    
                    FXend = FX(xend,yend,t); 
                    FYend = FY(xend,yend,t); 
                    FTend = FT(xend,yend,t);

                    
                    for i=2:numberofpointsinrho 
                        
                    % Define the q vectors
                    
                    q(i-1) = sqrt((matrixx(j,i) - matrixx(j,i-1))^2 + (matrixy(j,i) - matrixy(j,i-1))^2);
                    q(i) = sqrt((matrixx(j,i+1) - matrixx(j,i))^2 + (matrixy(j,i+1) - matrixy(j,i))^2);
                    q(1) = sqrt((matrixx(j,2) - matrixx(j,1))^2 + (matrixy(j,2) - matrixy(j,1))^2);
                    q(numberofpointsinrho) = sqrt((matrixx(j,numberofpointsinrho+1) - matrixx(j,numberofpointsinrho))^2 + (matrixy(j,numberofpointsinrho+1) - matrixy(j,numberofpointsinrho))^2);
            
                    % Define the coefficients of each term first
                    
                    % ------   LHS First half of matrix  ------
                    
                    LHcIminus1 = -1/drho ;
                    LHCI = 2/drho + (1/(2*dt*drho))*((q(i-1))^2 +(q(i))^2);
                    LHCIplus1 = -1/drho ;
            
                    
                    % ------ RHS forcing term Coefficients term ----- %
                    
                    % First perpindicular term
                    
                    xnonperp1 = matrixx(j,i) - matrixx(j,i-1);
                    ynonperp1 = matrixy(j,i) - matrixy(j,i-1);
                    xnonperp2 = matrixx(j,i+1) - matrixx(j,i);
                    ynonperp2 = matrixy(j,i+1) - matrixy(j,i);
                    
                    
                    % Using the perp operator in the clock wise direction
                    
                    perp1 = [ynonperp1;-xnonperp1];
                    perp2 = [ynonperp2;-xnonperp2];
                    
                    
                    RHforcingterm = (perp1 + perp2); % NOTE THIS IS A VECTOR
                    
                     % ------- RHS Forcing vector -------------%
                   
                    RHSforcingvec(i) = dot(perp1,[1;0])*p(j)*(q(i-1)/(2*drho)) + dot(perp2,[1;0])*p(j)*(q(i)/(2*drho)) ;
                    RHSforcingvec((numberofpointsinrho+1)+i) = dot(perp1,[0;1])*p(j)*(q(i-1)/(2*drho)) + dot(perp2,[0;1])*p(j)*(q(i)/(2*drho));
                    
                    % ----- LHS to be inverted matrix ------ %
                    
                    matrixtobeinverted(i,i-1) = LHcIminus1 ; matrixtobeinverted((numberofpointsinrho+1)+i,(numberofpointsinrho+1)+i-1) = LHcIminus1 ;
                    matrixtobeinverted(i,i) = LHCI ;         matrixtobeinverted((numberofpointsinrho+1)+i,(numberofpointsinrho+1)+i) = LHCI ;
                    matrixtobeinverted(i,i+1)= LHCIplus1 ;   matrixtobeinverted((numberofpointsinrho+1)+i,(numberofpointsinrho+1)+i+1)= LHCIplus1 ;

                    % ------------------------ RHS vectors -----------------------------------%
                    
                    % -------- FT vector --------------- %

                    FTvec = zeros(2*(numberofpointsinrho+1),1);
                    
                    FTvec(1) = -FT1*FX1*(((q(1))^2)/(2*drho))*dt;
                    
                    FTvec(numberofpointsinrho+1) = -FTend*FXend*((q(numberofpointsinrho))^2)/(2*drho)*dt ;
                    
                    FTvec(numberofpointsinrho+2) = FT1*FY1*((q(1))^2)/(2*drho)*dt ;
                    
                    FTvec(2*(numberofpointsinrho+1)) = FTend*FYend*((q(numberofpointsinrho))^2)/(2*drho)*dt;
                    
                    
                    % ----- RHS explicit terms ---- %
                    
                    RHSIvec(i) = (1/(2*dt*drho))*((q(i-1))^2 + (q(i))^2)*matrixx(j,i);
                    RHSIvec((numberofpointsinrho+1)+i) = (1/(2*dt*drho))*((q(i-1))^2 + (q(i))^2)*matrixy(j,i);
                    
                    end
                    
                    % ---- Boundary conditions for LHS to be inverted matrix --- %
                    
                    % First boundary condition for LHS side %
                    
                    matrixtobeinverted(1,1)= (((q(1))^2)/(2*drho))*((FX1)^2 + (FY1)^2) + ((FY1)^2)/(drho);  %(x_1^n)
                    matrixtobeinverted(1,2)= -((FY1)^2)/(drho);                                                %(x_2^n)
                    matrixtobeinverted(1,numberofpointsinrho+2)=  -(FY1*FX1)/(drho);                           %(y_1^n)
                    matrixtobeinverted(1,numberofpointsinrho+3)= (FY1*FX1)/(drho);                             %(y_2^n)

                    % J+1 boundary condition for LHS side % 

                    matrixtobeinverted(numberofpointsinrho+1,numberofpointsinrho) = -((FYend)^2)/(drho);  %(x_{J}^n)
                    matrixtobeinverted(numberofpointsinrho+1,numberofpointsinrho+1) = (((q(numberofpointsinrho))^2)/(2*drho))*((FXend)^2 + (FYend)^2) + ((FYend)^2)/(drho);  %(x_{J+1}^n)
                    matrixtobeinverted(numberofpointsinrho+1,(2*numberofpointsinrho)+1) = (FXend*FYend)/(drho); %(y_{J}^n)
                    matrixtobeinverted(numberofpointsinrho+1,2*(numberofpointsinrho+1))= -(FXend*FYend)/(drho);  %(y_{J+1}^n)
                    
                    % J+2 boundary condition for LHS side %
                    
                    matrixtobeinverted(numberofpointsinrho+2,1)= (FX1*FY1)/(drho);    %(x_1^n)
                    matrixtobeinverted(numberofpointsinrho+2,2)= -(FX1*FY1)/(drho);  %(x_2^n)
                    matrixtobeinverted(numberofpointsinrho+2,numberofpointsinrho+2)= -(((q(1))^2)/(2*drho))*((FX1)^2 + (FY1)^2) - ((FX1)^2)/(drho) ;   %(y_1^n)
                    matrixtobeinverted(numberofpointsinrho+2,numberofpointsinrho+3)= ((FX1)^2)/(drho); %(y_2^n)

                    % 2(J+1) boundary condition for LHS side % 
                    
                    matrixtobeinverted(2*(numberofpointsinrho+1),numberofpointsinrho)= -(FXend*FYend)/(drho);   %(x_{J}^n)
                    matrixtobeinverted(2*(numberofpointsinrho+1),numberofpointsinrho+1) = (FXend*FYend)/(drho);  %(x_{J+1}^n)
                    matrixtobeinverted(2*(numberofpointsinrho+1),(2*numberofpointsinrho)+1) = ((FXend)^2)/(drho);   %(y_{J}^n)
                    matrixtobeinverted(2*(numberofpointsinrho+1),(2*(numberofpointsinrho+1))) = -(((q(numberofpointsinrho))^2)/(2*drho))*((FXend)^2 + (FYend)^2) - ((FXend)^2)/(drho);   %(y_{J+1}^n)
                    
                    % ---- Boundary explicit terms (without F_t terms) and including forcing terms ----
                    
                    RHSIvec(1) = (((q(1))^2)/(2*drho))*((FX1)^2 + (FY1)^2)*x1 + ((FY1*q(1)*p(j))/(2*drho))*(FX1*(x2 -x1) + FY1*(y2-y1)) ;
                    
                    RHSIvec(numberofpointsinrho+1) =(((q(numberofpointsinrho))^2)/(2*drho))*((FXend)^2 + (FYend)^2)*xend + ((FYend*q(numberofpointsinrho)*p(j))/(2*drho))*(FXend*(xend-xendmin1) + FYend*(yend - yendmin1)); 
                    
                    RHSIvec(numberofpointsinrho+2) = -(((q(1))^2)/(2*drho))*((FX1)^2 + (FY1)^2)*y1 + ((FX1*q(1)*p(j))/(2*drho))*(FX1*(x2 -x1) + FY1*(y2-y1)) ;
                    
                    RHSIvec(2*(numberofpointsinrho+1)) = -(((q(numberofpointsinrho))^2)/(2*drho))*((FXend)^2 + (FYend)^2)*yend + ((FYend*q(numberofpointsinrho)*p(j))/(2*drho))*(FXend*(xend-xendmin1) + FYend*(yend - yendmin1));
                    
                    % --- Combine RHS known vector values together --- %
                    
                    RHSvec = RHSIvec + FTvec + RHSforcingvec;
                    Vecofnewvalues = inv(matrixtobeinverted)*RHSvec;
                    
                    for i=1:numberofpointsinrho+1
                    matrixx(j+1,i) = Vecofnewvalues(i);
                    matrixy(j+1,i) = Vecofnewvalues((numberofpointsinrho+1)+i);    
                    end
           
            
            end
            
            plot(matrixx(1,:),matrixy(1,:))
            hold on
            scatter(matrixx(numberofpointsintime/2,:),matrixy(numberofpointsintime/2,:),'r')
            hold on
            scatter(matrixx(numberofpointsintime+1,:),matrixy(numberofpointsintime+1,:),'g')
          
            
end