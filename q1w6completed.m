function [matrixx,matrixy] = q1w6completed(numberofpointsinrho,numberofpointsintime,tau)

        % We look to implement Curve shortening flow using finite elements where we
        % have parametrised some intial curve as vectors of (x,y). This will be for
        % the postive semi unit circle attached orthonally to the x axis.

                              %   Initial conditions
                                
        % For \rho in [0,1] we can parameterise the positive semi unit circle as
        % [cos\pi\rho, sin\pi\rho]^{T} depending on how many angle segments we choose
        % to implement. We will include a forcing term (f) so that it
        % can be used such that forcing is included.


        % Create the mesh for rho in [0,1] and time in [0,tau] where tau is the end
        % point (note as tau -> 0.5 we enter problems for p = 0 (no
        % forcing))
        
        % FOR Q1, FX FY and FT are all constants. Since the level set is
        % F(x,y) = y.
        
        FX = 0;
        FY = 1;
        FT = 0;

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
      
        RHSIvec = zeros(2*(numberofpointsinrho+1),1);
        RHSforcingvec = zeros(2*(numberofpointsinrho+1),1);

        % Set up initial conditions into the vector matrices above with
        % matrixx_0 = cos(pi*rho), matrixy_0 = sin(pi*rho).
        
        for i = 1:numberofpointsinrho+1
            matrixx(1,i) = cos(pi*pointrho(i));
            matrixy(1,i) = sin(pi*pointrho(i));
        end
        
        % We now look to build all of our matrix systems. Since these
        % change on each loop, we now start our iterated loop sequence over
        % time.
        
        j=1;
        
        
        
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
                    
                    for i=2:numberofpointsinrho 
                        
                    % Define the q vectors
                    
                    q(i-1) = sqrt((matrixx(j,i) - matrixx(j,i-1))^2 + (matrixy(j,i) - matrixy(j,i-1))^2);
                    q(i) = sqrt((matrixx(j,i+1) - matrixx(j,i))^2 + (matrixy(j,i+1) - matrixy(j,i))^2);
            
                    % Define the coefficients of each term first
                    
                    % ------   LHS First half of matrix  ------
                    
                    LHcIminus1 = -1/(q(i-1));
                    LHCI =(1/(2*dt))*(q(i-1) + q(i)) + 1/(q(i-1)) + 1/(q(i));
                    LHCIplus1 = -1/(q(i));
            
                    
                    % ------ RHS forcing term Coefficients term ----- %
                    
                    % First perpindicular term
                    
                    xnonperp1 = matrixx(j,i) - matrixx(j,i-1);
                    ynonperp1 = matrixy(j,i) - matrixy(j,i-1);
                    xnonperp2 = matrixx(j,i+1) - matrixx(j,i);
                    ynonperp2 = matrixy(j,i+1) - matrixy(j,i);
                    
                    % Using the perp operator in the clock wise direction
                    
                    perp1 = [ynonperp1;-xnonperp1];
                    perp2 = [ynonperp2;-xnonperp2];
                    
                    RHforcingterm = 0.5*(perp1 + perp2)*p(j); % NOTE THIS IS A VECTOR
                    
                    % ----- LHS to be inverted matrix ------ %
                    
                    matrixtobeinverted(i,i-1) = LHcIminus1 ; matrixtobeinverted((numberofpointsinrho+1)+i,(numberofpointsinrho+1)+i-1) = LHcIminus1 ;
                    matrixtobeinverted(i,i) = LHCI ;         matrixtobeinverted((numberofpointsinrho+1)+i,(numberofpointsinrho+1)+i) = LHCI ;
                    matrixtobeinverted(i,i+1)= LHCIplus1 ;   matrixtobeinverted((numberofpointsinrho+1)+i,(numberofpointsinrho+1)+i+1)= LHCIplus1 ;

                    % ------------------------ RHS vectors -----------------------------------%
                    
                    % ------- Forcing vector -------------%
                   
                    RHSforcingvec(i) = dot(RHforcingterm,[1;0]);
                    RHSforcingvec((numberofpointsinrho+1)+i) = dot(RHforcingterm,[0;1]);
                    
                    % -------- FT vector --------------- %
                    % Since FT = 0, this vector is just 0's %
                    
                    FTvec = zeros(2*(numberofpointsinrho+1),1);
                    
                    % ----- RHS explicit terms ---- %
                    
                    RHSIvec(i) = (1/(2*dt))*(q(i-1) + q(i)) * matrixx(j,i);
                    RHSIvec((numberofpointsinrho+1)+i) = (1/(2*dt))*(q(i-1) + q(i)) * matrixy(j,i);
                    
                    end
                    
                    % ---- Boundary conditions for LHS to be inverted matrix --- %
                    
                    % First boundary condition for LHS side %
                    
                    q(1) = sqrt((matrixx(j,2) - matrixx(j,1))^2 + (matrixy(j,2) - matrixy(j,1))^2);
                    
                    matrixtobeinverted(1,1)= (q(1)/(2*dt))*((FX)^2 + (FY)^2) + (FY/(q(1)))*(FY);  %(x_1^n)
                    matrixtobeinverted(1,2)=  -(FY/(q(1)))*(FY);                        %(x_2^n)
                    matrixtobeinverted(1,numberofpointsinrho+2)= -(FY/(q(1)))*(FX);   %(y_1^n)
                    matrixtobeinverted(1,numberofpointsinrho+3)= (FY/(q(1)))*(FX);  %(y_2^n)

                    % J+1 boundary condition for LHS side % ------ PROBLEM?
                    
                    q(numberofpointsinrho) = sqrt((matrixx(j,numberofpointsinrho+1) - matrixx(j,numberofpointsinrho))^2 + (matrixy(j,numberofpointsinrho+1) - matrixy(j,numberofpointsinrho))^2);
                    
                    matrixtobeinverted(numberofpointsinrho+1,numberofpointsinrho)= -(FY/(q(numberofpointsinrho)))*(FY);  %(x_{J}^n)
                    matrixtobeinverted(numberofpointsinrho+1,numberofpointsinrho+1)= (q(numberofpointsinrho)/(2*dt))*((FX)^2 + (FY)^2) + (FY/(q(numberofpointsinrho)))*FY; %(x_{J+1}^n)
                    matrixtobeinverted(numberofpointsinrho+1,(2*numberofpointsinrho)+1)= (FY/(q(numberofpointsinrho)))*(FX);   %(y_{J}^n)
                    matrixtobeinverted(numberofpointsinrho+1,(2*numberofpointsinrho+1))= (FY/(q(numberofpointsinrho)))*(-FX);  %(y_{J+1}^n)
                    
                    % J+2 boundary condition for LHS side %
                    
                    matrixtobeinverted(numberofpointsinrho+2,1)= (FX/(q(1)))*(FY);  %(x_1^n)
                    matrixtobeinverted(numberofpointsinrho+2,2)=  -(FX/(q(1)))*(FY); %(x_2^n)
                    matrixtobeinverted(numberofpointsinrho+2,numberofpointsinrho+2)= - (q(1)/(2*(dt)))*((FY)^2 +(FX)^2) - (FX/(q(1)))*(FX);   %(y_1^n)
                    matrixtobeinverted(numberofpointsinrho+2,numberofpointsinrho+3)= (FX/(q(1)))*(FX);  %(y_2^n)

                    % 2(J+1) boundary condition for LHS side % ------% PROBLEM?
                    
                    matrixtobeinverted(2*(numberofpointsinrho+1),numberofpointsinrho)= -(FX/q(numberofpointsinrho))*FY;  %(x_{J}^n)
                    matrixtobeinverted(2*(numberofpointsinrho+1),numberofpointsinrho+1)= (FX/q(numberofpointsinrho))*FY; %(x_{J+1}^n)
                    matrixtobeinverted(2*(numberofpointsinrho+1),(2*numberofpointsinrho)+1)= (FX/q(numberofpointsinrho))*FX;   %(y_{J}^n)
                    matrixtobeinverted(2*(numberofpointsinrho+1),(2*(numberofpointsinrho+1))) = - (q(numberofpointsinrho)/(2*dt))*((FX)^2 + (FY)^2) - (FX/q(numberofpointsinrho))*FX;  %(y_{J+1}^n)
                    
                    % ---- Boundary explicit terms (without F_t terms) and including forcing terms ----
                    
                    RHSIvec(1)= ((q(1)/(2*(dt)))*((FX)^2 + (FY)^2)*matrixx(j,1)) + ((p(j)*FY)/2)*((FY*(matrixy(j,2) - matrixy(j,1))) + FX*(matrixx(j,2) - matrixx(j,1)));
                    
                    RHSIvec(numberofpointsinrho+1) = (FY*p(j)/2)*(FY*(matrixy(j,numberofpointsinrho+1) - matrixy(j,numberofpointsinrho)) + FX*(matrixx(j,numberofpointsinrho+1) - matrixx(j,numberofpointsinrho))) + (q(numberofpointsinrho)/(2*dt))*((FX)^2 + (FY)^2)*matrixx(j,numberofpointsinrho+1);
                    
                    RHSIvec(numberofpointsinrho+2) = -(q(1)/(2*dt))*((FX)^2 + (FY)^2)*matrixy(j,1) + (p(j)*FX/2)*(FY*(matrixy(j,2) - matrixy(j,1)) + FX*(matrixx(j,2)-matrixx(j,1)));
                    
                    RHSIvec(2*(numberofpointsinrho+1)) = (p(j)*FX/2)*(FY*(matrixy(j,numberofpointsinrho+1) - matrixy(j,numberofpointsinrho)) + FX*(matrixx(j,numberofpointsinrho+1) - matrixx(j,numberofpointsinrho))) - (q(numberofpointsinrho)/(2*dt))*((FY)^2 + (FX)^2)*matrixy(j,numberofpointsinrho+1) ;
                    
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
            plot(matrixx(numberofpointsintime/2,:),matrixy(numberofpointsintime/2,:))
            hold on
            plot(matrixx(numberofpointsintime+1,:),matrixy(numberofpointsintime+1,:))
end