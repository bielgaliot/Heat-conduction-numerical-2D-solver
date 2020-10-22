%clc
clear 
close all
tic
%% Introducció de dades
%Físiques
rho1 = 1500; %[kg/m^3]
rho2 = 1600; %[kg/m^3]
rho3 = 1900; %[kg/m^3]
rho4 = 2500; %[kg/m^3]

cp1 = 750; %[J/(kgK)]
cp2 = 770; %[J/(kgK)]
cp3 = 810; %[J/(kgK)]
cp4 = 930; %[J/(kgK)]

k1 = 170; %[W/(mK)]
k2 = 140; %[W/(mK)]
k3 = 200; %[W/(mK)]
k4 = 140; %[W/(mK)]

%Geomètriques
p1x = 0.5; %[m]
p1y = 0.4; %[m]
p2x = 0.5; %[m]
p2y = 0.7; %[m]
p3x = 1.1; %[m]
p3y = 0.8; %[m]

W = 1; %[m]

%De contron
T_isotherm = 23+273.15;
Q_flow = 60;
Tg = 33+273.15;
alpha = 9;
To = 8+273.15;

%Numèriques
N = 40;
delta = 1e-6; %Criteri de convergència
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta = 1;        %SELECTOR ESQUEMA DE RESOLUCIÓ: 0 = EXPLÍCIT | 0.5 = C-N | 1 = IMPLÍCIT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fr = 1;
n=1;
timestep = 1; %[s/timestep]
tmax=10000;   %temps final de la simulació en segons
max_timesteps=round(tmax/timestep);
t = linspace(0,max_timesteps*timestep,max_timesteps+1); %vector temps
p=1; %comptador de timesteps avaluats

%% Definició grid i càlculs previs
rho = zeros(N+2,N+2);
cp = zeros(N+2,N+2);
k = zeros(N+2,N+2);

dx = p3x/N;
dy = p3y/N;

Vp = zeros(N+2,N+2);
dPE = zeros(N+2,N+2);
dPW = zeros(N+2,N+2);
dPN = zeros(N+2,N+2);
dPS = zeros(N+2,N+2);
dPe = zeros(N+2,N+2);
dPw = zeros(N+2,N+2);
dPn = zeros(N+2,N+2);
dPs = zeros(N+2,N+2);
Se = zeros(N+2,N+2);
Sw = zeros(N+2,N+2);
Sn = zeros(N+2,N+2);
Ss = zeros(N+2,N+2);

aE = zeros(N+2,N+2); 
aW = zeros(N+2,N+2);
aN = zeros(N+2,N+2);
aS = zeros(N+2,N+2);
aP = zeros(N+2,N+2);
bP = zeros(N+2,N+2);

sumQPn = zeros(N+2,N+2);

%Distribució de densitats, conductivitats i calors i creació de matrius de V, dPX, Sx.
for i=2:N+1
     for j=2:N+1
       
         if i<round((N*(p1x/p3x))) && j<round((N*(p1y/p3y)))
         rho(i,j)=rho1;
         cp(i,j)=cp1;
         k(i,j)=k1;
         end
         if i>=round((N*(p1x/p3x))) && j<round((N*(p2y/p3y)))
         rho(i,j)=rho2;
         cp(i,j)=cp2;
         k(i,j)=k2;
         end
         if i<round((N*(p1x/p3x))) && j>=round((N*(p1y/p3y)))
         rho(i,j)=rho3;
         cp(i,j)=cp3;
         k(i,j)=k3;
         end
         if i>=round((N*(p1x/p3x))) && j>=round((N*(p2y/p3y)))
         rho(i,j)=rho4;
         cp(i,j)=cp4;
         k(i,j)=k4;
         end
         
         Vp(i,j) = dx*dy*W; 
         dPE(i,j) = dx;
         dPW(i,j) = dx;
         dPN(i,j) = dy;
         dPS(i,j) = dy;
         dPe(i,j) = dx/2;
         dPw(i,j) = dx/2;
         dPn(i,j) = dy/2;
         dPs(i,j) = dy/2;
         
         Se(i,j) = dy*W;
         Sw=Se;
         Sn(i,j) = dx*W;
         Ss=Sn; 
     end
end

%Càlcul k harmònica
    k_har = zeros(N+2,N+2);
    k_harx = zeros(N+2,N+2);
    k_hary = zeros(N+2,N+2);
    
    for i=2:N+1
        for j=2:N+1
           k_harx(i,j) = dPE(i,j)./((dPe(i,j)./k(i,j))+(dPe(i+1,j)./k(i+1,j)));
           k_hary(i,j) = dPN(i,j)./((dPn(i,j)./k(i,j))+(dPn(i,j+1)./k(i,j+1)));
           k_har(i,j) = (k_harx(i,j) + k_hary(i,j))/2;
           
           %k per als contorns:
           if i==(N+1)
               k_har(i,j) = k2;
           end
           if j==(N+1) && i<round(N*(p2x/p3x))
               k_har(i,j) = k3;
           end
           if j==(N+1) && i==(round(N*(p2x/p3x))-1)
               k_har(i,j) = 182.35;
           end
           if j==(N+1) && i>=round(N*(p2x/p3x))
               k_har(i,j) = k4;
           end
       
        end
    end
 

%% Mapa inicial de temperatures

T = zeros(N+2,N+2,max_timesteps);
T(:,:,:) = To;
T_asterisc = zeros(N+2,N+2,max_timesteps);
T_data = zeros(N+2,N+2,max_timesteps);
T_anterior = zeros(N+2,N+2,max_timesteps);


%% Resolució i iteració GS

while p<=max_timesteps %p com a comptador de timesteps. Equivalent a: Següent increment?

    T(:,:,n) = T(:,:,n+1); %Update Tn
    
    while max(max(abs(T_asterisc(:,:,n+1)-T(:,:,n+1))))>delta %Comprovació d'error per a fer convergir.
         if fr==1
             T_asterisc(:,:,n+1) = T(:,:,n) + fr*(T(:,:,n+1)-T(:,:,n));
         else
        delta=1e-4;
        T(:,:,n)=T_asterisc(:,:,n+1);
        T_asterisc(:,:,n+1) = T(:,:,n) + fr*(T(:,:,n+1)-T(:,:,n));
         end
             for i=2:N+1
                for j=2:N+1
                    
                   %Suma de calors a l'instant anterior
                   sumQPn(i,j) = -k_har(i-1,j)*(T(i,j,n)-T(i-1,j,n))*Sw(i,j)/(dPW(i,j)) ...
                                 -k_har(i+1,j)*(T(i,j,n)-T(i+1,j,n))*Se(i,j)/(dPE(i,j)) ...  
                                 -k_har(i,j-1)*(T(i,j,n)-T(i,j-1,n))*Ss(i,j)/(dPS(i,j)) ... 
                                 -k_har(i,j+1)*(T(i,j,n)-T(i,j+1,n))*Sn(i,j)/(dPN(i,j));     

                   %Nodes interns
                   aE(i,j) = beta*(k_har(i+1,j)*Se(i,j))/dPE(i,j);
                   aW(i,j) = beta*(k_har(i-1,j)*Sw(i,j))/dPW(i,j);
                   aN(i,j) = beta*(k_har(i,j+1)*Sn(i,j))/dPN(i,j);
                   aS(i,j) = beta*(k_har(i,j-1)*Ss(i,j))/dPS(i,j);
                   aP(i,j) = aE(i,j) + aW(i,j) + aN(i,j) +  aS(i,j) + rho(i,j)*Vp(i,j)*cp(i,j)/timestep;
                   bP(i,j) = (1-beta)*sumQPn(i,j) + rho(i,j)*Vp(i,j)*cp(i,j)*T(i,j,n)/timestep;

                   %Condicions de contorn
                   %Top
                   aE(i,N+1) = 0;
                   aW(i,N+1) = 0;
                   aN(i,N+1) = 0;
                   aS(i,N+1) = beta*(k_har(i,N+1-1)*Ss(i,N+1))/dPs(i,N+1);
                   aP(i,N+1) = beta*(k_har(i,N+1-1)*Ss(i,N+1))/dPs(i,N+1) + rho(i,N+1)*Vp(i,N+1)*cp(i,N+1)/timestep;
                   bP(i,n+1) = Q_flow*Sn(i,N+1) + rho(i,N+1)*Vp(i,N+1)*cp(i,N+1)*T(i,N+1,n)/timestep; 


                   %Left qconv=qcond
                   aE(2,j) = beta*(k_har(2+1,j))/dPe(2,j);
                   aW(2,j) = 0;
                   aS(2,j) = 0;
                   aN(2,j) = 0;
                   aP(2,j) = beta*(k_har(2+1,j))/dPe(2,j) + alpha;
                   bP(2,j) = alpha*Tg;
                   
                   %Resolució
                   T(i,j,n+1) = (aE(i,j)*T(i+1,j,n+1) + aW(i,j)*T(i-1,j,n+1) + aS(i,j)*T(i,j-1,n+1) + aN(i,j)*T(i,j+1,n+1)...
                                + bP(i,j))/aP(i,j);
                   
                   %Condicions de contorn de T uniforme
                   %Bottom         
                   T(i,2,n+1) = 23+273.15; 
                   
                   %Right
                   T(N+1,j,n+1) = 273.15 + 8 + 0.005*t(p);
                   
                   %Cantonades
                   T(N+1,N+1,n+1)=T(N,N,n+1);
                   T(2,2,n+1)=T(3,3,n+1);
                   T(2,N+1,n+1)=T(3,N,n+1);
                   T(N+1,2,n+1)=T(N,3,n+1);
        
                end
             end 
   
    end
  
    T_data(:,:,p)=T(:,:,n+1); %%Guarda la temperatura a cada instant excepte per t=0 on T=To.

    T_asterisc(:,:,n+1) = T(:,:,n); %Actualització T estimada

    p=p+1; %seguent timestep
end

%% Tractament de dades
%Treure els bordes de la matriu de dades i passar a ºC
T_data_noborder = T_data-273.15;
T_data_noborder(1,:,:)=[];
T_data_noborder(end,:,:)=[];
T_data_noborder(:,1,:)=[];
T_data_noborder(:,end,:)=[];

T_trans = flip(permute(T_data_noborder, [2 1 3])); %Transposar la matriu de dades per al canvi d'eixos.
T_trans_5000s = transpose(T_data_noborder(:,:,round(5000/timestep))); %Transposar la matriu de dades per al canvi d'eixos. Imatge a 5000s (input en nº de timesteps)

%% Verificació als 5000s
h = pcolor(T_trans_5000s);
shading interp;
colorbar;
set(h, 'EdgeColor', 'none');
colormap jet
xticklabels(0.1375:0.1375:1.1);
yticklabels({'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8'})
set(gca,'TickLabelInterpreter','latex');
xlabel('$X(m)$','Interpreter','latex')
ylabel('$Y(m)$','Interpreter','latex')

%% Arxiu .dat
results = zeros(length(t)-1,3);
t(1)=[];
%Guardar valors d'interès
results(:,1) = t;
results(:,2) = squeeze(T_trans(round((0.65/1.1)*N),round((0.56/0.8)*N),:)); %omplir columnes de la martiu amb les dades
results(:,3) = squeeze(T_trans(round((0.74/1.1)*N),round((0.72/0.8)*N),:));

results = [0 To-273.15 To-273.15; results]; %Afegir condicions inicials 
results = round(results.*10000)./10000; %Arrodonir a 4 decimals

%Creació arxiu
delete GALIOT.dat %Eliminar-lo per sobreescriure dades a cada execució.
writematrix(results,'GALIOT.dat','Delimiter','tab');
toc