function [t_theo t_shell ft f1 f2 fa fr X N_bolts t_head A_actual pitch Bolt_circle_dia] = PED2(OpPressure,fmax,Di,J,Length,CA,head,Sf,rho,rhoflange,apex,h,m,Ya,N,fpermissible)
Wp = (1.1)*OpPressure; %Design pressure is 1.1 times operating pressure%
t = (Wp*Di)/((2*fmax*J)-Wp);
t_theo = t;
if(t<=4) %Minimum thickness of shell = 4mm %
    t_shell=4;
else
    t_shell = t;
end
Do = Di + 2*t_shell;
 
ft = (Wp*(Di+t_shell))/(2*t_shell); %hoop stress due to internal pressure%
 
f1 = (Wp*Di)/(4*t_shell); %axial stress due to internal pressure%
V_inside = (3.14*Di*Di*0.25*Length)/1000000;
Davg = (Do+Di)/2;
M_shell = 3.14*Davg*Length*(t_shell+CA)*rho*(10^-6); %Mass of shell%
V_sf = 3.14*0.25*Di*Di*Sf*(10^-9)*2; %Volume of straight flange%
 
 
 bo = N/2;
   if(bo<=6.3)
       b = bo;
   else
       b = 2.5*sqrt(bo);
   end
   r = sqrt((Ya-(m*Wp))/(Ya-((m+1)*Wp)));
   Num = 2*N*r;
   Denom = r-1;
   Go = Num/Denom; 
   Gi = (2*N)/(r-1); %Inner diameter of gasket%%
   Gavg = (Go + Gi)/2
   Wm1 = 3.14*b*Gavg*Ya; %Gasket reaction under atmospheric pressure%
   Wm2 = 3.14*2*b*Gavg*m*Wp + 3.14*Gavg*Gavg*Wp*0.25; %Gasket reaction under operating pressure%
   Nb = ((Gavg/10))/2.5;
   q = Nb/4; %Number of bolts rounded off to nearest multiple of 4%
   q1 = ceil(q);
   q2 = floor(q);
   q3 = (q1+q2)/2;
   q4 = q3*4;
   if(Nb<q4)
       N_bolts = 4*q2;
   else
       N_bolts = 4*q1;
   end
   if(Wm1>Wm2)
       Amt = Wm1/fpermissible;       Wm = Wm1;
   else
       Amt = Wm2/fpermissible; 
       Wm = Wm2;
   end
   A1b = Amt/N_bolts; %Area of one bolt%
  X1 = exp(log(A1b/0.51)/2.09);
  X2 = ceil(X1);
  u1 = rem(X2,2);
  if(u1==0) %standard bolt diameter%
      X = X2;
  else
      X = X2+1;
  end
  A_actual = N_bolts*0.51*(X^2.09); %total area of bolt%
  Abmax = (2*3.14*Ya*Gavg*N)/fpermissible;
  
  if(A_actual<Abmax)
      disp('Bolt size is acceptable');
  else
       disp('Bolt size is not acceptable');
  end
  pitch = 4.5*X; %Pitch of bolts is 4.5 times bolt diameter%
  Bolt_circle_dia = (pitch*N_bolts)/3.14; %Bolt Circle diamter%
  Flange_dia = (Bolt_circle_dia + 2*X)/1000;
  h_G = (Bolt_circle_dia-Gavg)/2;
  H = 3.14*0.25*Gavg*Gavg*Wp;
  y1 = (1.5*Wm*h_G*Gavg)/(H*Gavg);
  k = 1/(0.3+y1);
  tf = Gavg*sqrt(Wp/(k*fpermissible)); %Thickness of flange%
  Flange_Mass = 2*3.14*0.25*((Flange_dia)^2-(Di/1000)^2)*(tf/1000)*rhoflange;
  
 
if(strcmp(head,'A'))
     t_h = (Wp*Di*1.77)/(2*fmax*J); %Minimum thickness of head = 4mm%
      if(t_h<4)
          t_head =4;
      else
          t_head = t_h;
      end
          
          V_head = (0.081*Di*Di*Di)*(10^-9)*2; %Volume of head%
          B = (1.024*Do) + 0.67*0.06*Di + 2*Sf;
          M_head = 3.14*B*B*0.25*t_head*rho*(10^-9)*2;
          
          V_total = ceil(V_sf + V_head + V_inside); 
          M_content = 1000*V_total; 
          M_total = Flange_Mass + M_shell + M_head; %Total mass of vessel%
          M_final = 1.15*M_total + M_content; %Mass taken for calculation is 1.15 times actual mass of vessel%
      f2 = (M_final*9.8)/(3.14*t_shell*(Di+t_shell)); %Compressive stress due to weight of vessel and content%
  fa = f1 - f2; %total stress in axial direction%
  fr = sqrt(ft*ft - ft*fa + fa*fa); %Equivalent stress combining all other stress on basis of shear strain energy theory%
  
  if( fr<fmax && fa<fmax && ft<fmax)
      disp('thickness is ok')
  else
      disp('one of the stress values is over the limit,increase thickness or change material')
  end
%conical head%
 else     
      h_apex = (apex*pi)/360; %Apex angle in radians%
 
    t_c = (Wp*Di)/(2*fmax*J*cos(h_apex));
     if(t_c<4)
         t_head = 4;
     else
         t_head = t_c;
     end
    R1 = (Di/apex)*180;
    d1 = Di - ((Di*Di*h*h)/(R1*R1 - 0.25*Di*Di)); %inner diameter of cone%
    E1 = R1 - sqrt(h*h +0.25*((Di-d1)^2));
    A1 = 2*R1*sin(h_apex)+ 25;
    B1 = R1 - E1*cos(h_apex) + 25;
    
    M_cone = A1*B1*t_head*rho*2*(10^-9); %mass of conical head%
    V_cone = 0.262*h*(Di*Di +d1*d1 + Di*d1)*(10^-9)*2; %volume of conical head%
    
     V_total = ceil(V_sf + V_cone + V_inside);
          M_content = 1000*V_total;
          M_total = Flange_Mass + M_shell + M_cone;
          M_final = 1.15*M_total + M_content;
      f2 = (M_final*9.8)/(3.14*t_shell*(Di+t));
  fa = f1 - f2;
  fr = sqrt(ft*ft - ft*fa + fa*fa);
  
  if( fr<fmax && fa<fmax && ft<fmax)
      disp('thickness is ok')
  else
      disp('one of the stress values is over the limit,increase thickness or change material')
  end
          
end
 
end

