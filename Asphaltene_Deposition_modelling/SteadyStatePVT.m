function m = SteadyStatePVT(m,u,d,maxiter)

if isempty(maxiter),
    maxiter = 10; 
end

nst = m.th.nst;

if size(d,2)~=nst,
   error('wrong number of well sections')
end
m.u=u; m.d=d;

   
    dx  = m.th.L;  g = m.th.g; rhoL = m.th.rhoL; R=m.th.R; M=m.th.Ma;

%%
    p            = zeros(nst+1,1); 
    p(nst+1)     = u; %WHP in [bar]
    
% read pvt data file
    [pp, TT, rhogEq, rhoLEq, wgEq]=pvtData();  
    
% interpolate the equilibrium values
    %rhogE = interp2(pp,TT,rhogEq,(p(nst+1)*1e+5),(m.th.Ti(nst)-273));
    rhoLE = interp2(pp,TT,rhoLEq,(p(nst+1)*1e+5),(m.th.Ti(nst)-273));
    wgE = interp2(pp,TT,wgEq,(p(nst+1)*1e+5),(m.th.Ti(nst)-273));
    rhogE = interp2(pp,TT,rhogEq,(p(nst+1)*1e+5),(m.th.Ti(nst)-273)); 
    
%%    
    rhoL          = rhoLE; m.th.rhoL=rhoL; % const in the entire pipe

    %qL = zeros(nst,1);
    %qL(1) = d(1,1)*dx/rhoL;
    %for i=2:nst, qL(i)=qL(i-1)+d(1,i)*dx/rhoL; end   
   
    mdot=zeros(nst+1,1); % total mass flow rate [kg/s]
    mdot(1) = d(1,1)*dx + d(2,1)*dx;
    for i=2:nst, 
    mdot(i) = mdot(i-1) + d(1,i)*dx + d(2,i)*dx; end
    
    Fg = zeros(nst+1,1); % gas mass flow rate [kg/s]

    mdot(nst+1)   = mdot(nst);
    Fg(nst+1)     = mdot(nst+1)*wgE;

    rhoL(nst+1)  = rhoL;
    rhog(nst+1)  = p(nst+1)*1e5*M/(R*m.th.Z(nst)*m.th.Ti(nst)); % independent to equilibrium
    rhog(nst+1)  = rhogE;
    qL           = zeros(nst+1,1);
    qg           = zeros(nst+1,1);
    alfag        = zeros(nst+1,1);
    alfaL        = zeros(nst+1,1);
    
    qg(nst+1)    = Fg(nst+1)/rhog(nst+1);   
    qL(nst+1)    = (mdot(nst+1) - Fg(nst+1))/rhoL(nst+1);    
    
    ubsf=0.35;
    ubs          = ubsf*sqrt(g*m.th.IDt(nst))*sqrt(sin(m.th.incl(nst)))*(1+cos(m.th.incl(nst)))^1.2;
    alfag(nst+1) = qg(nst+1)/(m.th.DriftFlux.C0 * (qL(nst+1)+ qg(nst+1)) + m.th.A(nst) * ubs);
    alfaL(nst+1) = 1-alfag(nst+1);
    
    if alfag(nst+1) < 1e-6,
      momip1       = rhoL(nst+1)*qL(nst+1)*qL(nst+1)/alfaL(nst+1);
    else
      momip1       = rhog(nst+1)*qg(nst+1)*qg(nst+1)/alfag(nst+1) + rhoL(nst+1)*qL(nst+1)*qL(nst+1)/alfaL(nst+1); 
    end

    momip2       = momip1;
    
 %%   
    for i=nst:-1:1,

      Fg(i)  = Fg(i+1);
      qg(i)    = Fg(i)/rhog(i+1);
      qL(i)    = (mdot(i) - Fg(i))/rhoL(i+1); 
      
      ubs      = ubsf*sqrt(g*m.th.IDt(i))*sqrt(sin(m.th.incl(i)))*(1+cos(m.th.incl(i)))^1.2;
      alfag(i) = qg(i)/(m.th.DriftFlux.C0 * (qL(i)+ qg(i)) + m.th.A(i) * ubs);
      alfaL(i) = 1-alfag(i); 
      %
      
      if alfag(i) < 1e-6,
        momi       = rhoL(i+1)*qL(i)*qL(i)/alfaL(i);
      else
        momi       = rhog(i)*qg(i)*qg(i)/alfag(i) + rhoL(i+1)*qL(i)*qL(i)/alfaL(i); 
      end      
      
      rhomix   = alfaL(i)*rhoL(i+1) + alfag(i)*rhog(i);
      umix     = (1/m.th.A(i))*(qL(i) + qg(i));
          
      dpgrav   = 1e-5*dx*rhomix*g*sin(m.th.incl(i));
      dpfric   = 1e-5*dx*(0.5*(m.th.Friction.fTP(i)/m.th.IDt(i))*rhomix*abs(umix)*umix);
         
      dpmom    = 1e-5*(1/m.th.A(i)^2)*(momip1-momi);  %nb; no multiplication by dx
      p(i)     = p(i+1)+ dpmom+ dpgrav + dpfric;      
      rhog(i)  = p(i)*1e5*M/(R*m.th.Z(i)*m.th.Ti(i));
      
      % Equilibrium calc. 
      wgE = interp2(pp,TT,wgEq,(p(i)*1e+5),(m.th.Ti(i)-273));
      rhoL(i) = interp2(pp,TT,rhoLEq,(p(i)*1e+5),(m.th.Ti(i)-273));      
      rhog(i) = interp2(pp,TT,rhogEq,(p(i)*1e+5),(m.th.Ti(i)-273));
            
      Fg(i)     = mdot(i)*wgE;     
      qg(i)    = Fg(i)/rhog(i);     
      qL(i)    = (mdot(i) - Fg(i))/rhoL(i);      
      alfag(i) = qg(i)/(m.th.DriftFlux.C0 * (qL(i)+ qg(i)) + m.th.A(i) * ubs);    
      alfaL(i) = 1-alfag(i);
 
      if alfag(i) < 1e-6,
        momi       = rhoL(i)*qL(i)*qL(i)/alfaL(i);
      else
        momi       = rhog(i)*qg(i)*qg(i)/alfag(i) + rhoL(i)*qL(i)*qL(i)/alfaL(i); 
      end
           
      momip1   = momi;      
      %
    end

    for i=1:nst, m.x(3*(i-1)+3) = qL(i); end
    for i=1:nst, m.x(3*(i-1)+2) = p(i); end
    for i=1:nst, m.x(3*(i-1)+1) = alfaL(i); end    
            
%%        
    for iter=1:maxiter,
      iter
      momip1     = momip2;
      p_err      = zeros(nst,1); Fg_err = zeros(nst,1); rhog_err =zeros(nst,1); % 
      
      for i=nst:-1:1,
        % 
        p_prev   = p(i); Fg_prev = Fg(i); rhog_prev = rhog(i); % 
        qg(i)    = Fg(i)/rhog(i);     
        qL(i)    = (mdot(i) - Fg(i))/rhoL(i);      
        ubs      = ubsf*sqrt(g*m.th.IDt(i))*sqrt(sin(m.th.incl(i)))*(1+cos(m.th.incl(i)))^1.2;
        alfag(i) = qg(i)/(m.th.DriftFlux.C0 * (qL(i)+ qg(i)) + m.th.A(i) * ubs);    
        alfaL(i) = 1-alfag(i);
                
        rhomix   = alfaL(i)*rhoL(i) + alfag(i)*rhog(i);
        umix     = (1/m.th.A(i))*(qL(i) + qg(i));
                
        dpgrav   = 1e-5*dx*rhomix*g*sin(m.th.incl(i));
        dpfric   = 1e-5*dx*(0.5*(m.th.Friction.fTP(i)/m.th.IDt(i))*rhomix*abs(umix)*umix);

        if alfag(i) < 1e-6,
          momi       = rhoL(i)*qL(i)*qL(i)/alfaL(i);
        else
          momi       = rhog(i)*qg(i)*qg(i)/alfag(i) + rhoL(i)*qL(i)*qL(i)/alfaL(i); 
        end        
        
        dpmom    = 1e-5*(1/m.th.A(i)^2)*(momip1-momi);
        p(i)     = p(i+1)+ dpmom + dpgrav + dpfric;      
        rhog(i)  = p(i)*1e5*M/(R*m.th.Z(i)*m.th.Ti(i));

       % Equilibrium calc. 
        wgE = interp2(pp,TT,wgEq,(p(i)*1e+5),(m.th.Ti(i)-273)); 
        rhoL(i)= interp2(pp,TT,rhoLEq,(p(i)*1e+5),(m.th.Ti(i)-273));  
        rhog(i) = interp2(pp,TT,rhogEq,(p(i)*1e+5),(m.th.Ti(i)-273));
                
        Fg(i)     = mdot(i)*wgE;     
        qg(i)    = Fg(i)/rhog(i);     
        qL(i)    = (mdot(i) - Fg(i))/rhoL(i);      
        alfag(i) = qg(i)/(m.th.DriftFlux.C0 * (qL(i)+ qg(i)) + m.th.A(i) * ubs);    
        alfaL(i) = 1-alfag(i);
 
        if alfag(i) < 1e-6,
          momi       = rhoL(i)*qL(i)*qL(i)/alfaL(i);
        else
          momi       = rhog(i)*qg(i)*qg(i)/alfag(i) + rhoL(i)*qL(i)*qL(i)/alfaL(i); 
        end
           
        momip1   = momi;
        p_err(i)      = p(i) - p_prev;
        Fg_err(i)      = Fg(i) - Fg_prev;        
        rhog_err(i)      = rhog(i) - rhog_prev;
        %

        uGs(i)  = qg(i)/m.th.A(i);
        uLs(i)  = qL(i)/m.th.A(i);
        
      end
      
      for i=1:nst, m.x(3*(i-1)+3) = qL(i); end
      for i=1:nst, m.x(3*(i-1)+2) = p(i); end
      for i=1:nst, m.x(3*(i-1)+1) = alfaL(i); end      
            
      p_norm(iter)=norm(p_err); Fg_norm(iter)=norm(Fg_err); rhog_norm(iter)=norm(rhog_err); 
      %check derivatives
      %[dxdt] = TwoPhaseRLwell(0,m.x,u,d,m.th,nst);  normdxdt(iter) = norm(dxdt);
      if iter==10,
      piter10=p; alfagiter10=alfag; rhogiter10=rhog; end 
  
    end
    
    m.th.rhoL=rhoL(1);
    
    %figure; semilogy(p_norm); grid; title(['Norm:dp=p(new) - p(prev) vs. niter']); 
    %figure; semilogy(Fg_norm); grid; title(['Norm:dFg=Fg(new) - Fg(prev) vs. niter']); 
    %figure; semilogy(rhog_norm); grid; title(['Norm:drhog=rhog(new) - rhog(prev) vs. niter']);

%%    
    % geometry information
    
    %figure; plot((1:nst+1)',[m.th.A'; m.th.A(nst)],'-b','MarkerEdgeColor','k','MarkerSize',2); grid;
    %hold on
    %plot((1:nst+1)',[m.th.IDt'; m.th.IDt(nst)],'-r','MarkerEdgeColor','k','MarkerSize',2);
    %plot((1:nst+1)',[m.th.incl'; m.th.incl(nst)],'-g','MarkerEdgeColor','k','MarkerSize',2);
    %title(['Geometry variables: ALLI01']);
    %ylabel('geometry variables');
    %xlabel('cell number');
    %legend('Cross-sectional area, m2','Inner diameter, m','inclination','Location','NorthEast');
    %hold off
    
    % pressure   
    
    %figure; plot(((0:nst)*dx)',piter10,'-or','MarkerSize',4); grid;
    %hold on        
    %plot(((0:nst)*dx)',p,'-oB','MarkerSize',4);   
    %title (sprintf('Pressure drop: ALLI01; WHP=%0.3f; mdotL=%2.3f; mdotG=%2.3f',u,sum(d(1,:))*dx, sum(d(2,:))*dx));
    %ylabel('Pressure, bar');
    %xlabel('Pipe length, m');   
    %legend('iter=10',['iter=' num2str(maxiter)],'Location','NorthEast'); 
    %hold off
    
    % gas vf
    
    %figure; plot(((0:nst)*dx)',alfag,'-ob','MarkerSize',4); grid;    
    %hold on
    %title (sprintf('Gas volume fraction: ALLI01; WHP=%0.3f; mdotL=%2.3f; mdotG=%2.3f',u,sum(d(1,:))*dx, sum(d(2,:))*dx));
    %ylabel('VF, -');
    %xlabel('Pipe length, m');   
    %hold off    
    % gas density
    
    %figure; plot(((0:nst)*dx)',rhog,'-oc','MarkerSize',4); grid;    
    %hold on;
    %title (sprintf('Gas density: ALLI01; WHP=%0.3f; mdotL=%2.3f; mdotG=%2.3f',u,sum(d(1,:))*dx, sum(d(2,:))*dx));    
    %ylabel('rho, -kg/m3');
    %xlabel('Pipe length, m'); 
    %hold off
 
 %%
     uG = zeros(nst, 1); uL = zeros(nst, 1);

    for i=1:nst,
    
        if (alfag(i)>0)
         uG(i)=uGs(i)/alfag(i);
        end
  
        if (alfaL(i)>0)
         uL(i)=uLs(i)/alfaL(i);
        end
           
    end
    % output for matfile
    matname  = '.\driftFluxRes.mat';    
    
    rhog = rhog';
    rhoL = rhoL';
    uGs = uGs';
    uLs = uLs';
    PipeL= ((0:nst)*dx)';
     
    savename=['save(matname,''PipeL'',''p'',''alfag'''];
    savename=[savename ',''rhog'',''rhoL''']; 
    savename=[savename ',''mdot'',''Fg'''];     
    savename=[savename ',''uGs'',''uLs'''];     
    savename=[savename ',''uG'',''uL'''];
    savename=[savename ');'];
    eval(savename) 
    
