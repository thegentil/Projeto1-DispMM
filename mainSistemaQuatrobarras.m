function [] = mainSistemaQuatrobarras()

% Definindo os inputs

theta0 = 30*pi/180;
vin = 0.25; % m/s
l = 0.3; % m
g = 9.81; % m/s^2
m = 8; % Massa da barra 1
m2 = 5; % Massa das barras 2,3,4 e 5
mr = 3; % Massa do rolete
Gy = 300; %For�a
lvar = linspace(0,1.5,100);
m2var = linspace(0.5,10,100);

tmax = l*cos(theta0)/vin;

% Criando a s�rie de tempo (linspace)

tempo = linspace(0,tmax,120);

% Equacionando theta

x = cos(theta0)-vin.*tempo/l;

theta = acos(x);
theta1 = (90*pi)/180 - theta;

% CINEM�TICA

% Equacionando as velocidades 

Vbx = vin;
Vby = vin.*cot(theta);

Vg = 2*vin.*cot(theta); % Vgy

wbf = (-vin)./(l.*sin(theta));
wgb = (vin)./(l.*sin(theta));

% Equacionando as acelera��es 

alphabf = ((vin./l)^2).*(cos(theta)./(sin(theta).^3));
alphagb = -(wgb.^2).*cot(theta);


Ab = (-alphabf).*l.*cos(theta) - (wbf.^2).*l.*sin(theta); % Aby

Ag = Ab + alphagb.*l.*cos(theta) - (wgb.^2).*l.*sin(theta);

Acm4x = (-0.5).*alphagb.*l.*sin(theta) - 0.5.*(wgb.^2).*l.*cos(theta);

Acm4y = Ab + 0.5.*alphagb.*l.*cos(theta) - 0.5.*(wgb.^2).*l.*sin(theta);

Acm5x = (-0.5).*alphabf.*l.*sin(theta) + 0.5.*(wbf.^2).*l.*cos(theta);

Acm5y = (-0.5).*alphabf.*l.*cos(theta) - 0.5.*(wbf.^2).*l.*sin(theta);

Ary = Ab;

Arx = 0;

% DIN�MICA VARIANDO O PAR�METRO l (comprimento das barras 4 e 5) considerando a
% angula��o da barra de 30 graus constante (theta0) e a massa da
% barra perfeita.

wbflvar = (-vin)./(lvar.*sin(theta0));
alphabflvar = ((vin./lvar).^2).*(cos(theta0)./(sin(theta0).^3));
wgblvar = (vin)./(lvar.*sin(theta0));
alphagblvar = -(wgblvar.^2).*cot(theta0);
Ablvar = (-alphabflvar).*lvar.*cos(theta0) - (wbflvar.^2).*lvar.*sin(theta0);
Acm4xlvar = (-0.5).*alphagblvar.*lvar.*sin(theta0) - 0.5.*(wgblvar.^2).*lvar.*cos(theta0);
Acm4ylvar = Ablvar + 0.5.*alphagblvar.*lvar.*cos(theta0) - 0.5.*(wgblvar.^2).*lvar.*sin(theta0);
Arylvar = Ablvar;
Acm5xlvar = (-0.5).*alphabflvar.*lvar.*sin(theta0) + 0.5.*(wbflvar.^2).*lvar.*cos(theta0);
Acm5ylvar = (-0.5).*alphabflvar.*lvar.*cos(theta0) - 0.5.*(wbflvar.^2).*lvar.*sin(theta0);
Arxlvar = 0;
Gylvar = 300;

Fy4lvar = m2.*Acm4ylvar + Gylvar;
Fy5lvar = mr.*(Arylvar - g) + Fy4lvar;
Fx5lvar = ((m2.*lvar.*alphagblvar)/(3.*sin(theta0))) - ((Fy5lvar.*cos(theta0))/sin(theta0));
Oxlvar = m2.*Acm5xlvar + Fx5lvar;
Oylvar = m2.*Acm5ylvar + Fy5lvar;
Gxlvar = ((m2.*lvar.*alphabflvar)./(12.*sin(theta0))) - ((Fy4lvar.*cos(theta0))./(2.*sin(theta0))) - ((Gylvar.*cos(theta0))./(2.*sin(theta0))) + m2.*Acm4xlvar;
Fx4lvar = Gxlvar - m2.*Acm4xlvar;
Fmlvar = -(Fx5lvar + Fx4lvar + mr.*Arxlvar).*2;

% DIN�MICA VARIANDO O PAR�METRO m (massa das barras 4 e 5) considerando a
% angula��o da barra de 30 graus constante (theta0) e o comprimento de
% barra perfeito.
m2

wbfmvar = (-vin)./(l.*sin(theta0));
alphabfmvar = ((vin./l).^2).*(cos(theta0)./(sin(theta0).^3));
wgbmvar = (vin)./(l.*sin(theta0));
alphagbmvar = -(wgbmvar.^2).*cot(theta0);
Abmvar = (-alphabfmvar).*l.*cos(theta0) - (wbfmvar.^2).*l.*sin(theta0);
Acm4xmvar = (-0.5).*alphagbmvar.*l.*sin(theta0) - 0.5.*(wgbmvar.^2).*l.*cos(theta0);
Acm4ymvar = Abmvar + 0.5.*alphagbmvar.*l.*cos(theta0) - 0.5.*(wgbmvar.^2).*l.*sin(theta0);
Arymvar = Abmvar;
Acm5xmvar = (-0.5).*alphabfmvar.*l.*sin(theta0) + 0.5.*(wbfmvar.^2).*l.*cos(theta0);
Acm5ymvar = (-0.5).*alphabfmvar.*l.*cos(theta0) - 0.5.*(wbfmvar.^2).*l.*sin(theta0);
Arxmvar = 0;
Gymvar = 300;

Fy4mvar = m2var.*Acm4ymvar + Gymvar;
Fy5mvar = mr.*(Arymvar - g) + Fy4mvar;
Fx5mvar = ((m2var.*l.*alphagbmvar)/(3.*sin(theta0))) - ((Fy5mvar.*cos(theta0))/sin(theta0));
Oxmvar = m2var.*Acm5xmvar + Fx5mvar;
Oymvar = m2var.*Acm5ymvar + Fy5mvar;
Gxmvar = ((m2var.*l.*alphabfmvar)./(12.*sin(theta0))) - ((Fy4mvar.*cos(theta0))./(2.*sin(theta0))) - ((Gymvar.*cos(theta0))./(2.*sin(theta0))) + m2.*Acm4xmvar;
Fx4mvar = Gxmvar - m2var.*Acm4xmvar;
Fmmvar = -(Fx5mvar + Fx4mvar + mr.*Arxmvar).*2;

% DIN�MICA COM OS PAR�METROS PERFEITOS

Fy4 = m2.*Acm4y + Gy;
Fy5 = mr.*(Ary - g) + Fy4;
Fx5 = ((m2.*l.*alphagb)/(3.*sin(theta))) - ((Fy5.*cos(theta))/sin(theta));
Ox = m2.*Acm5x + Fx5;
Oy = m2.*Acm5y + Fy5;
Gx = ((m2.*l.*alphabf)/(12*sin(theta))) - ((Fy4.*cos(theta))/(2.*sin(theta))) - ((Gy.*cos(theta))/(2.*sin(theta))) + m2.*Acm4x;
Fx4 = Gx - m2.*Acm4x;
Fm = -(Fx5 + Fx4 + mr.*Arx).*2;

F4 = (Fy4.^2 + Fx4.^2).^(1/2);
F5 = (Fx5.^2 + Fy5.^2).^(1/2);
O = (Ox.^2 + Oy.^2).^(1/2);
G = (Gx.^2 + Gy.^2).^(1/2);
Gplot = 2.*G;


    

% Janela 2: Vby x t | Vby x Theta

figure(2)

whitebg('k')

subplot (1,2,1)
plot(tempo,Vby,'r','LineWidth',2)
title('Velocidade by x Tempo')
xlabel('Tempo (s)')
ylabel('Vby (m/s)')
grid
axis tight

subplot(1,2,2)
plot(theta*180/pi,Vby,'w','LineWidth',2)
title('Velocidade by x �ngulo')
xlabel('Theta (graus)')
ylabel('Vby (m/s)')
grid
axis tight

% Janela 3: wbf x t | wbf x Theta

figure(3)

whitebg('k')

subplot (1,2,1)
plot(tempo,wbf, 'r','LineWidth',2)
title('Veloc. �ng. barra 5 x Tempo')
xlabel('Tempo (s)')
ylabel('W barra 5 (m/s)')
grid
axis tight

subplot(1,2,2)
plot(theta*180/pi,wbf,'w','LineWidth',2)
title('Veloc. �ng. barra 5 x �ngulo')
xlabel('Theta (graus)')
ylabel('W barra 5 (m/s)')
grid
axis tight

% Janela 4: wgb x t | wgb x Theta

figure(4)

whitebg('k')

subplot (1,2,1)
plot(tempo,wgb,'r','LineWidth',2)
title('Veloc. �ng. barra 4 x Tempo')
xlabel('Tempo (s)')
ylabel('W barra 4 (m/s)')
grid
axis tight

subplot(1,2,2)
plot(theta*180/pi,wgb, 'w','LineWidth',2)
title('Veloc. �ng. barra 4 x �ngulo')
xlabel('Theta (graus)')
ylabel('W barra 4 (m/s)')
grid
axis tight

% Janela 5: Vg x t | Vg x Theta

figure(5)

whitebg('k')

subplot (1,2,1)
plot(tempo,Vg,'r','LineWidth',2)
title('Velocidade gy x Tempo')
xlabel('Tempo (s)')
ylabel('Vg (m/s)')
grid
axis tight

subplot(1,2,2)
plot(theta*180/pi,Vg, 'w','LineWidth',2)
title('Velocidade gy x �ngulo')
xlabel('Theta (graus)')
ylabel('Vg (m/s)')
grid
axis tight


% Janela 6: Ag x t | Ag x Theta

figure (6)

whitebg('k')

subplot (1,2,1)
plot(tempo,Ag,'r','LineWidth',2)
title('Acelera��o gy x Tempo')
xlabel('Tempo (s)')
ylabel('Ag (m/s^2)')
grid
axis tight

subplot(1,2,2)
plot(theta*180/pi,Ag,'w','LineWidth',2)
title('Acelera��o gy x �ngulo')
xlabel('Theta (graus)')
ylabel('Ag (m/s^2)')
grid
axis tight

% Janela 7: Ab x t | Ab x Theta

figure (7)

whitebg('k') 

subplot (1,2,1)
plot(tempo,Ab,'r','LineWidth',2)
title('Acelera��o by x Tempo')
xlabel('Tempo (s)')
ylabel('Ab (m/s^2)')
grid
axis tight

subplot (1,2,2)
plot(theta*180/pi,Ab,'w','LineWidth',2)
title('Acelera��o by x �ngulo')
xlabel('Theta (graus)')
ylabel('Ab (m/s^2)')
grid
axis tight

figure (8)

whitebg('k') 

subplot (1,2,1)
plot(tempo,alphabf,'r','LineWidth',2)
title('Alpha barra 5 x Tempo')
xlabel('Tempo (s)')
ylabel('Alpha barra 5 (rad/s^2)')
grid
axis tight

subplot (1,2,2)
plot(theta*180/pi,alphabf,'w','LineWidth',2)
title('Alpha barra 5 x �ngulo')
xlabel('Theta (graus)')
ylabel('Alpha barra 5 (rad/s^2)')
grid
axis tight

figure (9)

whitebg('k') 

subplot (1,2,1)
plot(tempo,alphagb,'r','LineWidth',2)
title('Alpha barra 4 x Tempo')
xlabel('Tempo (s)')
ylabel('Alpha barra 4 (rad/s^2)')
grid
axis tight

subplot (1,2,2)
plot(theta*180/pi,alphagb,'w','LineWidth',2)
title('Acelera��o barra 4 x �ngulo')
xlabel('Theta (graus)')
ylabel('Alpha barra 4 (rad/s^2)')
grid
axis tight

figure (10)

whitebg('k') 

subplot (1,2,1)
plot(theta*180/pi,F4,'r','LineWidth',2)
title('For�a Barra 4 x Theta')
xlabel('Theta (graus)')
ylabel('For�a Barra 4 (N)')
grid
axis tight

subplot (1,2,2)
plot(theta*180/pi,F5,'w','LineWidth',2)
title('For�a Barra 5 x Theta')
xlabel('Theta (graus)')
ylabel('For�a Barra 5 (N)')
grid
axis tight

figure (11)

whitebg('k') 

subplot (1,2,1)
plot(theta*180/pi,O,'r','LineWidth',2)
title('For�a Ponto O x Theta')
xlabel('Theta (graus)')
ylabel('For�a Ponto O (N)')
grid
axis tight

subplot (1,2,2)
plot(theta*180/pi,G,'w','LineWidth',2)
title('For�a Ponto G x Theta')
xlabel('Theta (graus)')
ylabel('For�a Ponto G (N)')
grid
axis tight

figure (12)

whitebg('k') 

subplot (1,2,1)
plot(tempo,Fm,'r','LineWidth',2)
title('For�a de entrada x Tempo')
xlabel('Tempo (s)')
ylabel('For�a de entrada (N)')
grid
axis tight

subplot (1,2,2)
plot(theta*180/pi,Fm,'w','LineWidth',2)
title('For�a de entrada x Theta')
xlabel('Theta (graus)')
ylabel('For�a de entrada (N)')
grid
axis tight


figure (13)

whitebg('k') 

subplot (1,2,1)
plot(tempo,Gplot,'r','LineWidth',2)
title('For�a de sa�da x Tempo')
xlabel('Tempo (s)')
ylabel('For�a de sa�da (N)')
grid
axis tight

subplot (1,2,2)
plot(tempo,Fm,'w','LineWidth',2)
title('For�a de entrada x Tempo')
xlabel('Tempo (s)')
ylabel('For�a de entrada (N)')
grid
axis tight

figure (14)

whitebg('k') 

plot(lvar,Fmlvar,'r','LineWidth',2)
title('For�a de entrada x Comprimento das barras 4 e 5')
xlabel('Comprimento da barra (m)')
ylabel('For�a de entrada (N)')
grid
axis tight

figure (15)

whitebg('k') 

plot(m2var,Fmmvar,'r','LineWidth',2)
title('For�a de entrada x Massa das barras 4 e 5')
xlabel('Massa das barras (kg)')
ylabel('For�a de entrada (N)')
grid
axis tight

end