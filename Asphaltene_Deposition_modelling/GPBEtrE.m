function [ dndt ] = GPBEtrE(t,n,da,dc,betaA,Sbr,npr,nst)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%load('driftFluxRes_sp1.mat');
load('driftFluxResdxVC500.mat');
%load('driftFluxRExC10s2.mat');
%dx=4.;
%dx=4049.47335514323/500;
%dx=8.09894671028646;
dx=PipeL(2);


dpar=da(1);
imax=da(2);
ca=da(3);
rhoa=da(4);
mpr=pi/6.*(dpar^3)*rhoa;
Navo=6.023e+23;
phia=5e-6;
%phia=1.;

n0 =[];
n0=[n0 (npr)*1e-3/Navo];
for i=2:imax
    n0=[n0 0];
end
    

%dndt= zeros(imax:nst,1);
dndt = [];

%%
    %alfaL(i),p(i),qL(i),IDt(i)



for i=1:nst
    A(i)=pi/4*((IDt(i))^2);
    uL(i)=qL(i)/(alfaL(i)*A(i));
end

    p(nst+1)=p(nst);

for i=1:nst
    dp(i)=p(i+1)-p(i);
    %Bpr(i)=-(0.001125)*(dp(i)/dx)/mpr*(1e-3/Navo);   % [kmol/m3 m]; production of primary particles
    Bpr(i)=0.; 
    % the minimum pressure at which precipitation occurs
    %if p(i) > 200.
    Bpr(i)=-(4.4025e-5)*(dp(i)/dx)/mpr*(1e-3/Navo);   % [kmol/m3-liquid m]; production of primary particles; Eskin(2011); Oil-C
    %end
    dc_Bpr(i)=dc(1,i)/rhoL(i)*(ca*rhoa)/(A(i))/(alfaL(i))/mpr*(1e-3/Navo); % [kmol/m3-liquid s]; inflow of primary particles; Oil-C      
end
%%

    BA(1:2,1:imax)=0.;
    DA(1:2,1:imax)=0.; 
    BB(1:imax)=0.;
    DB(1:imax)=0.;
    S(1:imax)=0.;
    
    betaA0=betaA(1:imax,1:imax)*phia*(Navo/1e-3); % m3/kmol s
    Sbr0=Sbr(1:imax,1);
    

for j = 1:imax
            if  j > 2
                for k = 1:j-2
                BA(1,j)= BA(1,j)+(2^(k-j+1))*betaA0(j-1,k)*n(j-1)*n(k);
                end
            end
    
            if j > 1
                BA(2,j)=1./2.*betaA0(j-1,j-1)*n(j-1)*n(j-1);
                %DB(j)=Sbr0(j)*n(j);             
            end
 
            %if j < imax
            %BB(j)=2.*Sbr0(j+1)*n(j+1);
            %end                    
                        
            %S(j) = BA(1,j) + BA(2,j) + BB(j) - DB(j) -uL(1)/dx*(n(j)-n0(j));             

            S(j)= -uL(1)/dx*(n(j)-n0(j));
end

    



for j=1:(imax-1)
            
            if j > 1
        
                for k = 1:j-1
                    DA(1,j)= DA(1,j)+(2^(k-j))*betaA0(j,k)*n(k);
                end    
                DA(1,j) = DA(1,j)* n(j);
            end          
        
            for k = j:(imax-1)
                DA(2,j) = DA(2,j)+betaA0(j,k)*n(k);
            end
    
            DA(2,j) = DA(2,j) *n(j);       
        
        
            %S(j) = S(j) - DA(1,j) - DA(2,j);               

end

%S(1)=S(1)+Bpr(1)*uL(1);

for j=1:imax

dndt = [dndt; S(j)];

end

%%
for i=2:nst
    
    % i:i+1
    betaA0=betaA(1:imax,1+(i-1)*imax:i*imax)*phia*(Navo/1e-3);
    Sbr0=Sbr(1:imax,i);

    BA(1:2,1:imax)=0.;
    DA(1:2,1:imax)=0.; 
    BB(1:imax)=0.;
    DB(1:imax)=0.;
    S(1:imax)=0.;



        for j = 1:imax
            if  j > 2
                for k = 1:j-2
                    BA(1,j)= BA(1,j)+(2^(k-j+1))*betaA0(j-1,k)*n(j-1+(i-1)*imax)*n(k+(i-1)*imax);
                end
            end
    
            if j > 1
                BA(2,j)=1./2.*betaA0(j-1,j-1)*n(j-1+(i-1)*imax)*n(j-1+(i-1)*imax);
                
                %DB(j)=Sbr0(j)*n(j+(i-1)*imax); 
            end
            
            %if j < imax
            %BB(j)=2.*Sbr0(j+1)*n(j+1+(i-1)*imax);
            %end
                
            S(j) = BA(1,j) + BA(2,j) + BB(j) - DB(j);
            
            %drhoqLn=rhoL(i)*alfaL(i)*uL(i)*A(i)*n(j+(i-1)*imax)-rhoL(i-1)*alfaL(i-1)*uL(i-1)*A(i-1)*n(j+(i-2)*imax);
            
            % correct, 06/06/13            
            dqLn=alfaL(i)*uL(i)*A(i)*n(j+(i-1)*imax)-alfaL(i-1)*uL(i-1)*A(i-1)*n(j+(i-2)*imax);
            
            
            %S(j) = S(j) - drhoqLn/dx/(rhoL(i)*alfaL(i)*A(i)); 
            
            %S(j) = - drhoqLn/dx/(rhoL(i)*alfaL(i)*A(i));
            
            % correct, 06/06/13
            %S(j) = S(j) - dqLn/dx/(alfaL(i)*A(i));
            
            S(j) = - dqLn/dx/(alfaL(i)*A(i));
            
        end

        for j=1:(imax-1)
            
            if j > 1        
                for k = 1:j-1
                    DA(1,j)= DA(1,j)+(2^(k-j))*betaA0(j,k)*n(k+(i-1)*imax);
                end    
                DA(1,j) = DA(1,j)*n(j+(i-1)*imax);
            end          
        
            for k = j:(imax-1)
                DA(2,j) = DA(2,j)+betaA0(j,k)*n(k+(i-1)*imax);
            end
    
            DA(2,j) = DA(2,j) *n(j+(i-1)*imax);                      

            %dndt(j+(i-1)*imax) = dndt(j+(i-1)*imax)- DA(1,j) - DA(2,j);             
            %S(j) = S(j)- DA(1,j) - DA(2,j); 
            
        end
        
    %S(1)=S(1)+Bpr(i)*uL(i)+dc_Bpr(i);
    %S(1)=S(1)+Bpr(i)*uL(i);
        
        
    for j=1:imax

    dndt = [dndt; S(j)];

    end        
        
        
        
        
end


end








    
    


