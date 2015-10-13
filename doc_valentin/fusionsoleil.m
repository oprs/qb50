%Simulation de l'observateur non linéaire.
function [Xfinal]=fusionsoleil(duree)
	% clear all

	clc


	%% ===Tables de B===
	global B_tab;
	B_tab = textread('dataMag.csv','','delimiter',';');

	% disp(B_tab);

	%% ===Définition des variables et initialisation===

	%On doit redéclarer ces variables pour qu'elles soient globales
	global T
	T=duree;
	global Torbite
	Torbite=5400;
	global Ttot
	Ttot=duree;%La durée que je rajoute à la fin pour voir si ça reconverge
	%Attention ici, changement de définition
	global TIME_CONTROLE ;%Inutile de déclarer cette valeur mais certain scripts
	TIME_CONTROLE=duree;%de François utilisent cette variable globale
	global moment
	moment=zeros(1,3);

	fe=1;%Fréquence d'échantillonnage 
	global freq
	freq=fe;

	%% ===Plot des courbes===
	plotomega=1;
	plotbiais=1;
	plotquat=1;%plot q
	plottheta=1;%pointing knowlege error
	% plotthetaref=1;%pointing accuracy error
	plotaccuracyquat=0;
	plotmoment=0;
	% plotpuissance=1;
	plotforces=0;
	plotenergie=0;
	plotlyapu=0;

	%% Dynamiques perturbatrices
	global gradient 
	global frot 
	global circuit
	gradient=1*0;
	frot=1*0;
	circuit=1*0;

	%% ====Gains du filtre===
	global kB kS kI kP
	kB=0.3; kS=0.1; kI=3e-4; kP=1;%Attention, kI sera à choisir en fonction de la dérive de biais de gyro (à estimer)


	%% =Variance des bruits=

	%Bruit sur les gyros
	global sigmagyro;
	sigmagyro=0.011*(pi/180)*sqrt(fe)*4;
	% fprintf('Ecart-type du bruit de gyro : %f\n', sigmagyro);

	%Bruit sur le magnéto
	global sigmamagneto
	sigmamagneto=5e-3*sqrt(fe)*1e-4*4;%En Tesla

	%Bruit sur le capteur solaire
	global sigmasolaire
	sigmasolaire=5e-2*sqrt(fe)*4;

	%Coefficient aérodynamique global
	global Cd
	Cd=2.2;


	%===Etat Initial===
	% Torbite=5400;
	fprintf('==============================================\n');

	fprintf('Initialisation partiellement aléatoire des conditions initiales\n');


	%Quaternion définissantla matrice de rotation
% 	q = [2*(rand(1,4)-.5)];
	q=[3 6 -4 -20];
	%q=[0 0 0 1];
	q=q/norm(q);
	disp('Quaternion initial :');
	disp(q);
	theta0=norm(2*asin(norm(q(1:3)))*180/pi);
	fprintf('Erreur initiale de %f ° \n', theta0);
	qref=[0 0 0 1];
	omegaref=[0 2*pi/Torbite 0];

	%Vecteur rotation du satellite
	omega =[0 2*pi/Torbite 0]+1e-6*randn(1,3);
% 	omega=[30 -30 30]*pi/180*0+5e-4*randn(1,3);
	fprintf('\nOmega inital x1000:\n');
	disp(1000*omega);
	fprintf('\nOmega de référence x1000:\n');
	disp(1000*omegaref);
	fprintf('Norme de l erreur initiale sur omega : %e \n', norm(omega-omegaref));
	%Biais de mesure sur omega
	b = 5e-6*randn(1,3);
	fprintf('\nBiais :\n');
	disp(b);
	%Observateur pour le quaternion
	qo =[0 0 1 1];
    qo=qo/norm(qo);
	%Observateur pour le biais :
	bo = zeros(1,3);
	%Vecteur global stockant le moment magnétique généré pour le tracé des
	%courbes
	moment=zeros(1,3);


	%Définition de l'état X initial

	X = [q omega b qo bo qref omegaref]';

	% X = [q omega b qo bo qref omegaref]';


	%=Constantes pour l'intégration=

	pas=1/fe;%Période d'échantillonnage
	n=int64(Ttot*fe);%nombre d'itérations
	vectTemps= 0:pas:Ttot;%Vecteur sotckant les points de temps
	global vecteurTemps
	vecteurTemps=vectTemps; % c'est sale, c'est le vecteur utilisé dans 
	%f pour trouver index
	vectX(:,1)=X;%Vecteur stockant l'état X à chaque itération

	vecjour=ones(size(vectTemps));
	veccontrol=zeros(size(vectTemps));
	vectmoment=zeros(3, length(vectTemps));
	%===Résolution de l'équation différentielle avec RK4===



	t=0;%initialisation du temps
	for i=2:(n+1)%Subtilité avec le nombre d'itérations,faut être sûr qu'on a le même...
		            %nombre de points de temps et d'itérations.
		
		touslesxpourcents=5;
		if(mod(100*i/touslesxpourcents,n)==0)
		    fprintf('\nLoading...%d%%\n',100*i/(n+1));
		end

	  
		k1 = f(t,X);
		k2 = f(t + pas/2 , X + (pas/2)*k1);
		k3 = f(t + pas/2 , X + (pas/2)*k2);
		k4 = f(t + pas , X + pas*k3);
		X = X + (pas/6) * (k1 + 2*k2 + 2*k3 + k4);
		vecjour(i)=decroissancejour(t);
		veccontrol(i)=controle(t);
		t=t+pas;
		vectX(:,i)=X;
		vectmoment(:,i)=moment;
	end;
	vectTemps=vectTemps/Torbite;
	%===Output de la fonction===
	Xfinal=vectX(:,n+1);


	%===Estimation des erreurs===

	%matrice des données à tracer, une colonne par point
	%Convergence pour le biais sur omega
	vectb=[vectX(8,:)-vectX(15,:);...
		vectX(9,:)-vectX(16,:);...
		vectX(10,:)-vectX(17,:)];
	



	%Convergence pour le quaternion (pas sûr de l'ordre, quel inverse faut-il prendre?)
	vectqknow=zeros(4,n);
	vectqacc=zeros(4,n);
	for i=1:n+1
	%     vectqknow(:,i)=quatprod(vectX(1:4,i)',invQ(vectX(11:14,i)'));
	%     vectqacc(:,i)=quatprod(vectX(1:4,i)',invQ(vectX(18:21,i)'));
		vectqknow(:,i)=quatprod(invQ(vectX(11:14,i)'),vectX(1:4,i)');
		vectqacc(:,i)=quatprod(invQ(vectX(18:21,i)'),vectX(1:4,i)');
	end;

	%Erreur angulaire
	theta=zeros(1,n+1);
	thetaref=zeros(1,n+1);
	for i=1:n+1
	%     theta(i)=norm(2*asin(norm(vectq(1:3,i)))*180/pi);
		theta(i)=norm(acos( QrotInv(vectqknow(:,i),[1;0;0])'*[1;0;0] )*180/pi);
		thetaref(i)=norm(acos( QrotInv(vectqacc(:,i),[1;0;0])'*[1;0;0] )*180/pi);
	%     thetaref(i)=norm(2*asin(norm(vectqref(1:3,i)))*180/pi);
	end


	%===Tracé des courbes===

	if(plotbiais==1)
		%Tracé de vectb
		figure('color', 'white');
	%     const=max(vectb(:));
		plot(vectTemps, vectX(8:10,:), vectTemps, vectX(15:17,:), 'linewidth', 2);
		title('$Convergence de l''estimation du biais$', 'interpreter', 'latex', 'fontsize',13);
		xlabel('Time $[T_0]$', 'interpreter', 'latex', 'fontsize',13);
		ylabel('Différence entre observation et biais réel', 'fontsize',13);
		legend('bx', 'by', 'bz', 'box', 'boy', 'boz');
		grid
	end


	if(plotquat==1)
		%Tracé de vectq
		figure('color', 'white');
		
		subplot(211);
		plot(vectTemps,vectqknow, 'linewidth', 2);
		title('Quaternion knowledge error', 'interpreter', 'latex', 'fontsize',13);
		xlabel('Time $[T_0]$', 'interpreter', 'latex', 'fontsize',13);
		ylabel('$\widehat{q}^{-1} q$', 'interpreter', 'latex', 'fontsize',13);
		legend('err1', 'err2', 'err3', 'err4');
		grid
		
		subplot(212);
		plot(vectTemps,vectqacc, 'linewidth', 2);
		title('Quaternion accuracy error', 'interpreter', 'latex', 'fontsize',13);
		xlabel('Time $[T_0]$', 'interpreter', 'latex', 'fontsize',13);
		ylabel('$q_{ref}^{-1} q$', 'interpreter', 'latex', 'fontsize',13);
		legend('err1', 'err2', 'err3', 'err4');
		grid
	end

	vectomega=(vectX(5:7,:));
	%vectomegao=(vectX(31:33,:));
	if(plotomega==1)
		figure('color', 'white');
		subplot(211);
		plot(vectTemps,vectomega, '-', 'linewidth', 2);
		title('Evolution of $\omega$', 'interpreter', 'latex', 'fontsize',13);
		xlabel('Time $[T_0]$', 'interpreter', 'latex', 'fontsize',13);
		ylabel('$\omega$', 'interpreter', 'latex', 'fontsize',13);
		legend('X', 'Y', 'Z');
		grid
		
		subplot(212);
		plot(vectTemps,.5*log10(  sum( ( vectomega-omegaref'*ones(1,length(vectTemps)) ).^2 )  ), vectTemps, .5*log10(sum(vectb.^2)), 'linewidth', 2);
		title('$\log(||\omega-\omega_{ref}||) et \log(||\widehat{b}-b||)$', 'interpreter', 'latex', 'fontsize',13);
		xlabel('Time $[T_0]$', 'interpreter', 'latex', 'fontsize',13);
		ylabel('$\log(\epsilon)$', 'interpreter', 'latex', 'fontsize',13);
		legend('log(error(omega))', 'log(error(b))');
		grid
	end

	if(plottheta==1)
	%     const=max(theta);
		const=10;
		%Pointing knowledge error [°]
		figure('color', 'white');
		%vectTemps, 2.1*veccontrol, 'k',
		
		subplot(211);
		plot( vectTemps, const*vecjour, 'r--',  vectTemps, theta, vectTemps, thetaref, 'linewidth', 2);
		xlabel('Time $[T_0]$', 'interpreter', 'latex', 'fontsize',13);
		ylabel('Pointing  error [�]');
		legend('Light Intensity', 'Pointing knowledge error', 'Pointing accuracy error');
		grid
		title('Pointing error ', 'fontsize', 14);
		
		
		subplot(212);
		plot(vectTemps, log10(theta/180), vectTemps, log10(thetaref/180), 'linewidth', 2);
		xlabel('Time $[T_0]$', 'interpreter', 'latex', 'fontsize',13);
		ylabel('log(Pointing error [°]/180)');
		legend('log(Pointing knowledge error/180)', 'log(Pointing accuracy error/180)');
		grid
		title('log(Pointing error/180)');
		


	end
	% 
	if(plotaccuracyquat==1)
		figure('color', 'white');
		plot(vectTemps,vectqacc, 'linewidth', 2);
		title('Accuracy quaternion error');
		xlabel('Time $[T_0]$', 'interpreter', 'latex', 'fontsize',13);
		ylabel('q*qref^-1');
		grid
	end
	% figure;
	% plot(vectTemps,vectX(5:7,:)-vectX(22:24,:), 'linewidth', 2);
	% title('Evolution de omega-omegaref');
	% xlabel('Temps(s)');
	% ylabel('omega-omegaref');
	% legend('x', 'y', 'z');
	% grid

	% if(plotthetaref==1)
	%     %Pointing accuracy error
	%     figure('color', 'white');
	%     plot(vectTemps, thetaref, 'linewidth', 2);
	%     xlabel('Time [T0]');
	%     ylabel('Pointing accuracy error [°]');
	%     legend('Pointing accuracy error');
	%     grid
	%     title('Pointing accuracy error ');
	% end

	% vectmoment=vectX(25:27,:);
	if(plotmoment==1)
		% Torque
		figure('color', 'white');
		plot(vectTemps, vectmoment, 'linewidth', 2);
		xlabel('Time $[T_0]$', 'interpreter', 'latex', 'fontsize',13);
		ylabel('Torque $[A.m^2]$', 'interpreter', 'latex', 'fontsize',13);
		legend('X', 'Y', 'Z');
		grid
		title('Torque generated by the coils', 'interpreter', 'latex', 'fontsize',13);
	end

	%Calcul de l'intensité/puissance parcourant les bobines
	%il faudrait faire une fonction spécifique quand les données seront sues
	surfacex=(75e-3)^2; surfacey=surfacex; surfacez=surfacex;
	Nx=1500; Ny=Nx; Nz=Nx;
	vectintensite=[vectmoment(1,:)/(surfacex*Nx); vectmoment(2,:)/(surfacey*Ny); vectmoment(3,:)/(surfacez*Nz)];

	%calcul de la consommation 
	resistance=60;
	vectpuissance=zeros(1,length(vectintensite));
	vectenergie=vectpuissance;
	for compt=1:length(vectintensite)
		vectpuissance(compt)=resistance*(vectintensite(1,compt)^2+vectintensite(2,compt)^2+vectintensite(3,compt)^2);
		if(compt>1)
		    vectenergie(compt)=vectpuissance(compt)*pas+vectenergie(compt-1);
		end
		
	end


	if(plotenergie==1)
		figure('color', 'white');
		
		subplot(211)
		plot(vectTemps, vectpuissance, 'linewidth', 2);
		xlabel('Time $[T_0]$', 'interpreter', 'latex', 'fontsize',13);
		ylabel('Power [Watt]');
	%     legend('X', 'Y', 'Z');
		grid
		title('Power consumption');
		
		subplot(212)
		plot(vectTemps, vectenergie, 'linewidth', 2);%vectenergie/(3600*7.2) pour les A.h
		xlabel('Time $[T_0]$', 'interpreter', 'latex', 'fontsize',13);
		ylabel('Energy [J]');
		grid
		title('Power consumption');
	end



	%Recalcul du moment aero
	%Ce calcul est fait une seconde fois pour ne pas avoir une autre variable
	%globale...
	n1 = [1;0;0];
	n2 = [0;1;0];
	n3 = [0;0;1];

	surf1=.1*.1; surf2=.1*.2; surf3=surf2;
	central=[-1e-2, 0, 0];
	vectmomaero=zeros(3,length(vectTemps));
	vectmomgrad=vectmomaero;
	% disp(size(vectmomaero));
	% disp(size(vectTemps));
	om0 = sqrt(6.67e-11*5.97e24/((6370e3+250e3)^(3)));
	I=getInertie();
	v=7.8e3;
	for j=1:length(vectTemps)
		
		uv = QrotInv(vectqacc(:,j)',[1 0 0]); %expression de la vitesse dans le repère satellite
		uv = uv/norm(uv);
		Kfrot1=0.5*Cd*airdensity(vectTemps(j)*Torbite)*v^2*surf1*dot(uv,n1)*cross3(uv, central);
		Kfrot2=0.5*Cd*airdensity(vectTemps(j)*Torbite)*v^2*surf2*dot(uv,n2)*cross3(uv, central);
		Kfrot3=0.5*Cd*airdensity(vectTemps(j)*Torbite)*v^2*surf3*dot(uv,n3)*cross3(uv, central);
		Kfrot = Kfrot1+Kfrot2+Kfrot3;
		vectmomaero(:,j)=Kfrot';
		qorbsat=vectqacc(:,j);
		zrot = [(2*(-qorbsat(4)*qorbsat(2)+qorbsat(1)*qorbsat(3)));(2*(qorbsat(2)*qorbsat(3)+qorbsat(4)*qorbsat(1)));(1-2*qorbsat(1)^2-2*qorbsat(2)^2)];

		vectmomgrad(:,j) = 3*om0^2*monCross(zrot,I*zrot)';
	end


	if (plotforces==1)
		figure('color', 'white');
		
		subplot(221);
		plot(vectTemps, vectmomaero, 'linewidth', 2);
		xlabel('Time $[T_0]$', 'interpreter', 'latex', 'fontsize',13);
		ylabel('Torque $[N.m]$', 'interpreter', 'latex', 'fontsize',13);
		legend('aeroX', 'aeroY', 'aeroZ');
		grid
		title('Aerodynamical disturbance torque', 'interpreter', 'latex', 'fontsize',13);
		
		subplot(222);
		plot(vectTemps, log10(airdensity(vectTemps*Torbite)), 'linewidth', 2);
		xlabel('Time $[T_0]$', 'interpreter', 'latex', 'fontsize',13);
		ylabel('Air density (log) $[kg.m^{-3}]$', 'interpreter', 'latex', 'fontsize',13);
		grid
		title('Air density', 'interpreter', 'latex', 'fontsize',13);
		
		
		subplot(223);
		plot(vectTemps, vectmomgrad, 'linewidth', 2);
		xlabel('Time $[T_0]$', 'interpreter', 'latex', 'fontsize',13);
		ylabel(' Torque $[N.m]$', 'interpreter', 'latex', 'fontsize',13);
		legend('gradX', 'gradY', 'gradZ');
		grid
		title('Gradient gravity disturbance torque', 'interpreter', 'latex', 'fontsize',13);
		
		subplot(224);
		plot(vectTemps, vectmoment, 'linewidth', 2);
		xlabel('Time $[T_0]$', 'interpreter', 'latex', 'fontsize',13);
		ylabel('Torque $[A.m^2]$', 'interpreter', 'latex', 'fontsize',13);
		legend('X', 'Y', 'Z');
		grid
		title('Torque generated by the coils', 'interpreter', 'latex', 'fontsize',13);
	end


	if(plotlyapu==1)
		tracerlyapu(vectqacc, vectomega);
	end

	
		fprintf('==============================================\n');
end

%% ===Equation différentielle vérifiée par X====

function dX=f(t,X)

	global sigmamagneto
	global sigmasolaire
	global kB kP kI kS
	global Cd
	global moment

	%On redécompse X en ses différentes valeurs
	q=X(1:4)';
	omega=X(5:7)';
	b=X(8:10)';
	qo=X(11:14)';
	bo=X(15:17)';
	qref=X(18:21)';
	omegaref=X(22:24)';


	%Test au cas où la norme du quaternion explose (problème de schéma
	%numérique)
	if(norm(qo)>1.3)
		disp('-------------------------------------------');
		fprintf('t=%.1f\n', t);
		disp('q=');
		disp(q);
		disp('qo=');
		disp(qo);
		return
	end


	%% Lecture des mesures des capteurs
	%Mesure de omega
	omegam=mesure(omega, b);



	%Magnetometer
	Bgc=champmag(t);
	Bs = QrotInv(q,Bgc)+sigmamagneto*randn(1,3);
	Bspasnorm=Bs;
	if(norm(Bs)~=0)
		Bs=Bs/norm(Bs);
	end
	% 
	% if(isnan(Bs(1)))
	%     disp('Bs is NaN:');
	%     disp(Bs);
	%     disp('Bgc');
	%     disp(Bgc);
	%     disp('Quaternion');% C'est lui !
	%     disp(q);
	%     return
	%     
	% end


	Bso=QrotInv(qo,Bgc); %"Estimée" de Bs cf notation de l'article
	Bsopasnorm=Bso; %Notatioin un peu sale : le vecteur que l'on ne renorme pas
	if(norm(Bso)~=0)
		Bso=Bso/norm(Bso);
	end

	%Sun sensors
	Sgc=Sgctime(t);
	Ss = QrotInv(q,Sgc)+sigmasolaire*randn(1,3);
% 	Ss = estimerS9(Ss, sigmasolaire)';
	Ss = Ss/norm(Ss);
	Sso= QrotInv(qo,Sgc);%"Estimée" de Ss idem
	Sso = Sso/norm(Sso);%on norme ces vecteurs

	%On récupère les constantes kB, kS, kI, kP définies en haut, et on modifie
	%kB et kS en fonction du soleil et du contrôle
	%On cherche ici des gains optimaux pour minimiser l'erreur. Lorsqu'on a
	%convergé, il est bon d'avoir des gains petits, et il est souhaitable d'en
	%avoir des gros auparavant pour converger rapidement, mais comment les
	%choisir ? kP, kI, kB et kS sont des constantes globales qui "marchent
	% %assez bien" pour T=5400s

	%% Choix des gains de l'observateyur
	kPi=kP;

	% kIi=0.003/10;
	kIi=kI;
	kBi=kB*(1-controle(t));%+kSi;
	% kBi=kB;
	kSi=kS*decroissancejour(t)/8;
	global Torbite

	if t>Torbite/3
		kPi=kP/3;
		kIi=kIi/10;
		kSi=kS*decroissancejour(t)/12;
		kBi=kB/3*(1-controle(t));
	end

	%%%%%%%Control%%%%%%


	%% CALCUL DU CONTROLE U
	%calcul de l'ecart par rapport a l'etat de reference
	if detumbling(t)==1
		    Bgc=Bgctime(t);
		    k=10000;
		    m = -k*cross3(QrotInv(qo,Bgc),omegam-bo-omegaref);
	else
		if(controle(t)==1)      
		    dq = quatprod(invQ(qo), qref);
		    m=(-4.8e-5*cross3(Bsopasnorm, omegam-bo-omegaref)-3e-7*cross3(Bsopasnorm, dq(1:3)))/norm(Bsopasnorm)^2;
	%             m=(-1e-4*cross3(Bsopasnorm, omegam-bo-omegaref)-2.5e-10*cross3(Bsopasnorm, dq(1:3)))/norm(Bsopasnorm)^2;
		else
		    m=zeros(1,3);
		end
	end

	%%%%%%%%%%%Erreur sur le moment=Erreur de contrôle (Skype X-Mines)
	%Ici c'est une erreur de contrôle statique, a essayer avec un bruit blanc
	% donc la variance est adéquate
	%m(1)=1.1*m(1);
	%m(2)=0.95*m(2);
	%m(3)=1.03*m(3);
	%%%%%%%%%

	if(max(abs(m))>0.3)
		m = 0.3/max(abs(m))*m;
	end
	moment=m';
	Kt = cross3(m',Bspasnorm(1:3)');
	Ks= Kt';
	%%%%%%%%%%%%%%%%%%%%
	%moment généré par les bobines
	% Ks=momentsat(t,q);
	%Evaluation du moment dû au gradient de gravité
	%qorbsat quaternion passage de orbital à body
	%% Calcul du gradient
	global gradient
	qorbsat = quatprod(invQ(qref),q);
	if(gradient)
	   
		%définition du quaternion de rotation entre le repere orbital et le
		%repere du satellite
		
		%definition du vecteur associe (cf rapport)
		zrot = [(2*(-qorbsat(4)*qorbsat(2)+qorbsat(1)*qorbsat(3)));(2*(qorbsat(2)*qorbsat(3)+qorbsat(4)*qorbsat(1)));(1-2*qorbsat(1)^2-2*qorbsat(2)^2)];
		om0 = sqrt(6.67e-11*5.97e24/((6370e3+250e3)^(3)));
		I=getInertie();
		Kgrad = 3*om0^2*monCross(zrot,I*zrot);
		Ks = Kgrad'+Ks;
	end
	%% Calcul des frottements
	global frot
	if(frot)
		v= (7.8e3);
		coeffdrag2=.5*Cd*airdensity(t)*(v^2+0*(5e-2*norm(omega))^2);
		n1 = [1;0;0];
		n2 = [0;1;0];
		n3 = [0;0;1];
		uv = QrotInv(qorbsat,[1 0 0]); %expression de la vitesse dans le repère satellite
		uv = uv/norm(uv); 
		surf1=.1*.1; surf2=.1*.2; surf3=surf2;
		central=[-1e-2, 0, 0];
		Kfrot1=coeffdrag2*surf1*dot(uv,n1)*cross3(uv, central);
		Kfrot2=coeffdrag2*surf2*dot(uv,n2)*cross3(uv, central);
		Kfrot3=coeffdrag2*surf3*dot(uv,n3)*cross3(uv, central);
		Kfrot = Kfrot1+Kfrot2+Kfrot3;
		Ks = Kfrot+Ks;
		
	end

	%% On ajoute un couple constant : celui des circuits
	global circuit
	Ks=Ks+cross3([0,0,4e-4],Bspasnorm)*circuit;
	%% Observateur
	%Notations du papier de recherche
	wmes=kBi*cross3(Bs,Bso)+kSi*cross3(Ss,Sso);


	[dX(22:24), dX(18:21)]=physic(t, omegaref, qref, 0);%l'état de référence
	[dX(5:7), dX(1:4)]=physic(t, omega, q, Ks');%cf fonction physic
	dX(8:10)=[0 0 0] + 5e-4*pi/180*randn(1,3)/10*0;%le biais constant (lentement variable)
	dX(11:14)=1/2*quatprod(qo,[omegam-bo+kPi*wmes 0])+(1-norm(qo)^2)*qo;%observateur
	%de l'article sur le quaternion + un terme de normalisation
	dX(15:17)=-kIi*wmes;%observateur de l'article


	% dX(25:33)=zeros(size(X(25:33)));
	dX=dX';

end


function Bgc=Bgctime(t)%en gauss
	global Torbite;
	global B_tab;%Le champ est en nT
	% compteur_mag=1+mod(floor(360*t/T)+90,360);
	compteur_mag=1+floor(360/Torbite*mod(t,Torbite));
	Bgc = 1e-9*[B_tab(compteur_mag,2) B_tab(compteur_mag,3) B_tab(compteur_mag,4)];
end


%Fonction qui donne le vecteur solaire dans le ref géo. Pour le moment je
%l'ai inventée
function Sgc=Sgctime(t)
	global T
	% Sgc=(ones(1,3))/sqrt(3);
	Sgc=[-1 0 0];
end


%La fonction de mesure du gyro, soumise à un bruit et à un biais
function omegam=mesure(omega, b)
	global sigmagyro;
    omegam=omega+b+sigmagyro*randn(1,3);
end

%moment généré par les bobines
% function Ks=momentsat(t,q)
% Bgc=Bgctime(t);%Le moment est basé sur le champ réel extérieur
% Bs=quatprod(quatprod(invQ(q),[Bgc 0]),q);
% Bs=Bs(1:3);
% ms=1e-3*[1,1,1]*controle(t); % ici penser à changer le moment en fonction du contrôle
% Ks=-cross3(ms,Bs);
% end


% function res=decroissancejour(t)
% global T
% res=1/(1+exp(-(T2/2-t)*20/T2));
% res=0;
% res=1/(1+exp(-(T/3.5-t)))+1/(1+exp((3*T/4-t)));
% end

function res=decroissancejour(t)
	global Torbite
	tt = mod(t,Torbite);
	if(tt<Torbite/2)
		res = luminosite(tt);
	else
		res = luminosite(Torbite-tt);
	end
%res=1/(1+exp(-(Torbite/2-t)*20/Torbite));
end

function res=luminosite(tt)
	global Torbite
	temp = 80;
	if(tt<Torbite/4-temp)
		res = 1;
	else
		if(tt>Torbite/4+temp)
		    res = 0;
		else
		    res = -tt/(2*temp)+0.5+Torbite/(8*temp);
		end
	end
end



function cont=controle(t)
	cont=0;
	global Torbite
	% if( (t>=Torbite && t<5*Torbite) || (t>7*Torbite && t<8*Torbite) || (t>10*Torbite && t<11*Torbite) )
	% if t>1*Torbite && (mod(t, 100)/100)>0.3
	% if (mod(t, 100)/100)>0.3
	%     cont=1;
	% end
	if (mod(t, 200)/200)>0.3 && t>Torbite/2
		cont=1;
	end
end

function res=detumbling(t)%Fonction qui déclenche le détumbling classique
	res=0;
	global Torbite
	% if( (t>=Torbite && t<5*Torbite) || (t>7*Torbite && t<8*Torbite) || (t>10*Torbite && t<11*Torbite) )
	if t<2*Torbite
		res=0;
	end
end

function res=airdensity(t)%Evolution de la densité de l'air
	global T
	% res=ones(size(t)).*(10.^(2*t/T-12));
	res=ones(size(t))*1e-12;
end

function res=cross3(x,y)%Une version plus rapide du produit vectoriel en dim 3
	res=zeros(size(x));
	res(1)=x(2)*y(3)-x(3)*y(2);
	res(2)=y(1)*x(3)-y(3)*x(1);
	res(3)=x(1)*y(2)-x(2)*y(1);
end


