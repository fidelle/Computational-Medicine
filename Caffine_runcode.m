% Initialize structure
p = struct();

% Constants
p.q=0;
%p.V = 45;           % L (volume of distribution)
p.ka = 5.94;        % hr-1 (rate constant for absorption for oral dosing)
%p.ke = 0.58;        % hr-1 (rate constant for elimination)
tspan = [0:0.1: 5];     % Time span for simulation, hrs
tspan2=[5:0.1:10];
p.kc=log(2)/7;
p.i = 1;            % Number of doses given, used for mass balance
weights = [112, 133, 174, 210, 255] * 0.453592; % weights of each subject 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Switch-Case
simulation = 'optim';
switch simulation 
    
case 'optim'
        
        % ----- Optimization -----
        % Experimental Data
        t_exp = [0.25, 1.5, 3, 4.75, 6, 8];
        
        % concentration (mg/L)
        y_exp = [3.5, 3.5, 3.0, 2.2, 2.5; 
            4.0, 3.5, 3.5, 2.5, 2.5,;
            3.1, 2.5, 3.0, 2.0, 2.0; 
            2.5, 2.0, 3.0, 1.5, 1.5;
            4.5, 3.5, 4.8, 3.0, 3.0;
            3.5, 2.5, 4.0, 2.0, 2.5];
        
        
        tspan = [0:0.1: 5];     % Time span for simulation, hrs
        tspan2=[5:0.1:10];
        % Set parameters
        p.dose = 152;         % mg (drug input)
        %p.V = 0.7           % L (volume of distribution)
        %tspan = 0:0.1:8;   % Time span for simulation, hrs
        %p.D0 = dose/p.V;    % Initial drug concentration, mg/L
        
        % Set initial conditions
        y0 = zeros(3,5);
        y0(:,1) = [0 0 p.dose];   % Subject 1
        y0(:,2) = [0 0 p.dose];   % Subject 2
        y0(:,3) = [0 0 p.dose];   % Subject 3
        y0(:,4) = [0 0 p.dose];   % Subject 4
        y0(:,5) = [0 0 p.dose];   % Subject 5
        
        
        
        k0(1,1:5)=4; %/hr absorbtion 
        k0(2,1:5)=0.5; %rate constant for clearance
        k0(3,1:5)=40; %l/kg vol of distribution

% %       

        
        % The cost function file requires multiple inputs to allow it to
        % run inside the loop. Because lsqnonlin requires the cost function
        % to only take the parameters to be optimized as input, this
        % anonymous function gives the cost function the other parameters
        cost = @(k0) caffeine_cost(k0, t_exp, y_exp, tspan,tspan2, y0, p, weights);
        
        % Optimization parameters
        lower = [0 0 0];      % Lower bound
        upper = [Inf Inf Inf];      % Upper bound
        options = optimoptions('lsqnonlin', 'Display', 'iter', ...
            'StepTolerance', 1e-4, 'MaxFunctionEvaluations', 6e2);
        
        % Run lsqnonlin for the anonymous cost function
        optimal = lsqnonlin(cost, k0, lower, upper, options);
        
        % Set the k values to the optimization results
        p.vr1 = optimal(1);
        p.vr2 = optimal(2);
        p.vr3 = optimal(3);
        p.vr4 = optimal(4);
        p.vr5 = optimal(5);
        p.kc = optimal(6); 

        
%         %%getting vals for 2c
%         Y=[];
%         for i = 1:5
% %             p.ka=optimal(1,i);
% %             p.kc=optimal(2,i);
% %             p.v=optimal(3,i);
%         %    
%         %    % p.ka=k0(1,i);
%             p.ka=5.94
%             p.kc=optimal(1,i);
%             p.v=optimal(2,i);
% 
%             % Run the ODE solver
%             [T1i, Y1i, ~, ~] = caffeine_main(tspan,y0(:,i),p);
% 
%             y0new=Y1i(end,:);
%             y0new(end,3)=y0new(end,3)+92;
% 
% 
%             [T2i, Y2i, ~, ~] = caffeine_main(tspan2,y0new,p);
% 
%             Ti=[T1i;T2i];
%             Yi=[Y1i;Y2i];
% 
%             % Select the output from the central compartment
%             Yi = Yi(:,1);
% 
%             % Store output
%             Y = cat(2,Y,Yi);
% 
%         end

        
%         % ----- Save Output -----
%         t_exp = t_exp';     % Make a column vector to match y_exp
%         
%         %save exp.mat t_exp y_exp
%         save optim2b.mat Ti Y
%         


    %case 'localSens'
%% local sensitivity 
      %% run base case
      vd= [35.5616   42.2294   55.2475   66.6780   80.9662];

      dose1=152;
      dose2=92;
      
      p0.q=0;
      
      %p0.dose1=dose1;
      p0.dose=dose1;
      
      p0.dose2=dose2;
      
      p0.i=1;
      
      %tspan = [0:0.1: 5];     % Time span for simulation, hrs
      tspan1=[0.5:0.1:5];
      tspan2=[5:0.1:10];
      
      tspan30=[0:0.1:0.5];
      
      
      
      y0=[0 0 p0.dose]; %initial conditions
      
      baseAUC=zeros(1,5); %output conditions matrix initialized
      baseAUC30=zeros(1,5); %auc after 30mins
      baseCost=zeros(1,5); %cost
      
      
     Y=[]; %store y vals
      
      
      for i=1:5 
          
          p0.ka=optimal(1,i);
          p0.kc=optimal(2,i);
          p0.v=optimal(3,i);
          %p.v=vd(m);
        [Td0i, Yd0i, ~, AUCd0] = caffeine_main(tspan30, y0, p0); %for 30 mins
%         
          y30i = Yd0i(end,:);%save output to var 
%                     
          y30_update=Yd0i(end,:); %update output
%          
          AUC30base= AUCd0;
          
          baseAUC30(i)=baseAUC30(i)+ AUC30base; %gets auc after 30mins
%         

        
        [Td1i, Yd1i, ~, AUCd1] = caffeine_main(tspan1, y30_update, p0);
        
                
         % Select initial conditions for next time period
         y0i = Yd1i(end,:);
                    
         y0_update=Yd1i(end,:);
         y0_update(end,3)=y0_update(end,3)+p0.dose2;
         
         T1toti=[Td0i;Td1i]; %gets values after 30mins up to 5 hrs
         Y1toti=[Yd0i;Yd1i];

    
        [Td2i, Yd2i, ~, AUCd2] = caffeine_main(tspan2,y0_update,p0);%time after 2nd dose

         Ti=[T1toti;Td2i];
         Yi=[Y1toti;Yd2i];
         
          Y = cat(2,Y,Yi); %concatenate values for each patient

         
         AUCbase=AUCd1+AUCd2+AUCd0;
         
         baseAUC(i) = baseAUC(i) + AUCbase;
         
         
        y_diff = zeros(length(t_exp),1);

           % For each experimental time point
        for k = 1:length(t_exp)

            % Difference between the output and experimental time
            t_diff = abs(Ti - t_exp(k));

            % Finding the value of T that minimizes the difference
            [~, t_index] = min(t_diff);

            % Find the difference between the experimental concentrations and model
            % output at that time point
           y_diff(k,:) = Yi(t_index, 1) - y_exp(k,i); %ydiff per patient
           %disp(size(Yi(t_index,1)))

        end
         
          SSE = norm(y_diff,2)^2;
          baseCost(i)=baseCost(i)+SSE;

         
      end 
         
         
%      
         
         
         
         
         
     
      % ----- Sensitivity -----
      
         
%       %simulation parameters
        delta = 0.1;
        params = {'p.dose','p.dose2',  'p.ka', 'p.kc', 'p.v'};
        dose1=152;
        dose2=92;
        
      
       
         AUC = zeros(length(params),5);
         AUC30=zeros(length(params),5);
         parCost=zeros(length(params),5);

         p.q=0;
         p.dose=dose1;
         p.dose2=dose2;
         p.i=1;


         Y=[];
     for m=1:5 %loop over patients
%            
            
           
           p.ka=optimal(1,m);
           p.kc=optimal(2,m);
           p.v=optimal(3,m);

%            

        tspan1=[0.5:0.1:5];
        tspan2=[5:0.1:10];

        tspan30=[0:0.1:0.5];
      
         % y02=[0,0,p.dose];
%            
%          % Initialize output matrix (rows = parameters, columns = dosing)
%           
%           
             for z=1:length(params) %loop over param
                eval(sprintf('%s = %s * (delta+1);', params{z}, params{z})); %turns to interger

                eval(sprintf('disp(%s)', params{z}));
                y02=[0,0,p.dose];
                
                
                
                  [t30m, y30m, ~, auc30m] = caffeine_main(tspan30, y02, p); %for 30 mins
    %         
                  Y30parm = y30m(end,:);%save output to var 
        %                     
                  y30_updatem=y30m(end,:); %update output
        %          
                  auc30param= auc30m;

                  AUC30(z,m)=AUC30(z,m)+ auc30param; %gets auc after 30mins
                
%                   [Td0, Yd0, ~, AUCd0] = caffeine_main(tspan30, y02, p); %for 30 mins
%         %         
%                   y30i = Yd0i(end,:);%save output to var 
%         %                     
%                   y30_update=Yd0i(end,:); %update output
%         %          
%                   AUC30base= AUCd0;
% 
%                  AUC30(z,m)=AUC30(z,m)+ AUC30base; %gets auc after 30mins
             
                
%                 
%                 
%                  
%                 
                [t1m, y1m, ~, auc1] = caffeine_main(tspan1, y30_updatem, p);
                
                 % Select initial conditions for next time period
                 y02m = y1m(end,:);

                 y02_update=y1m(end,:);
                 y02_update(end,3)=y02_update(end,3)+p.dose2;
                 
                 t1totm=[t30m;t1m]; %gets values after 30mins up to 5 hrs
                 y1totm=[y30m;y1m];



                [t2m, y2m, ~, auc2] = caffeine_main(tspan2,y02_update,p);

                 T2m=[t1totm;t2m];
                 Y2m=[y1totm;y2m];
                 
                 AUCpar=auc1+auc2+auc30param;
                 
                AUC(z,m) = AUC(z,m) + AUCpar;
                
                
                
                
                %%% cost -----------
                
                 y_diffpar = zeros(length(t_exp),1);

                   % For each experimental time point
                for k = 1:length(t_exp)

                    % Difference between the output and experimental time
                    t_diff = abs(T2m - t_exp(k));

                    % Finding the value of T that minimizes the difference
                    [~, t_index] = min(t_diff);

                    % Find the difference between the experimental concentrations and model
                    % output at that time point
                   y_diffpar(k,:) = Y2m(t_index, 1) - y_exp(k,m); %ydiff per patient
                   %disp(size(Yi(t_index,1)))

                end

                  SSEpar= norm(y_diffpar,2)^2;
                  
                  parCost(z,m)=parCost(z,m)+SSEpar;
              
                eval(sprintf('%s = %s / (delta+1);', params{z}, params{z}));
                
             end
              % Change parameter back to original value
             

         end 
%---------------Calculate change in AUC----------
        deltaAUC = AUC - baseAUC;
        
        % Convert to percentage
        percentAUC = deltaAUC./baseAUC;
        
        % Percentage change in each parameter
        percentP = [0.1 0.1 0.1 0.1 0.1]';
        
        % Percent change in AUC divided by percent change in parameter
        AUCchange = percentAUC./percentP;
 
%----------------Calculate change in AUC 30mins----------------
        deltaAUC30 = AUC30 - baseAUC30;
        
        % Convert to percentage
        percentAUC30 = deltaAUC30./baseAUC30;
        
        % Percentage change in each parameter
        
        % Percent change in AUC divided by percent change in parameter
        AUC30change = percentAUC30./percentP;
        
        
%----------------Calculate change in cost----------------
        deltaCost = parCost - baseCost;
        
        % Convert to percentage
        percentCost = deltaCost./baseCost;
        
        % Percentage change in each parameter
        
        % Percent change in AUC divided by percent change in parameter
        Costchange = percentCost./percentP;
        
        
        save AUCchange.mat AUCchange
        save AUC30change.mat AUC30change
        save costchange.mat Costchange

 
 end 
      
      
  %% run base case
      
  
      
          