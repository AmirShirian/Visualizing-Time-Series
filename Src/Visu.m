clc
close all
clear
dt=0.1;			%Sampling Period

%-------------------------Choosing Case
Case=0;

if Case==0
	 t_final=100;
	 T=0:dt:t_final;
	 
	 f1=30;
	 f2=28;
	 f3=20;
	 
	 rb=60;
	 Win=100;
	 
	 a0=[zeros(1,rb+Win),ones(1,length(T)-(rb+Win))];
	 a=a0.*T;
	 b0=[zeros(1,rb+2*Win),ones(1,length(T)-(rb+2*Win))];
	 b=b0.*T;
	 c0=[zeros(1,rb+3*Win),ones(1,length(T)-(rb+3*Win))];
	 c=c0.*T;
	 d0=[zeros(1,rb+4*Win),ones(1,length(T)-(rb+4*Win))];
	 d=d0.*T;
	 f4=0.1;
	 A1=0*sin(f4*T)+2;
	 %---------------------------------------------------
	 in_sig=(40*cos(f1*a)+50*cos(f2*b)-10*cos(f3*c)+10+80*sin(26*d)+1).*A1;					%Input Signal
end

if Case==1
	Data
end 


n_max=100;		%Maximmun Number of Modes
n_start=1;		%Starting Number of Modes
%---------------------Initializing RLS Parameters
Q=length(in_sig)-n_max;
Y=(in_sig(n_max+1:end))';
U=zeros(Q,n_max);
for j=1:n_max
    U(:,j)=in_sig(n_max+1-j:end-j);
end
%---
PBig0=ones(n_max+1)*0.01;
Lan=0.1:0.005:1;			%Range of Landas
L_Lan=length(Lan);

for j=1:L_Lan
    Pr(j).PBig=PBig0;
    Pr(j).Lan=Lan(j);
end
Pr(1).PWin=PBig0;

len_Pic = (2*pi*0.5)/dt; 		%Length of Picture
ddt=0.05;						%Step beetwen each freq
PIC = [];

for id=1:Q
    yd=Y(id);
    ud=(U(id,:))';
    Ther_max=-1;
    index_max=0;
    %Forgetting factor
    for j=1:L_Lan			%Update Cov Matrix with RLS Algorithm
        PBig=Pr(j).PBig;
        lan=Lan(j);
        PBig=PBig*lan+[ud;yd]*[ud;yd]';
        Pr(j).PBig=PBig;
        His=[];
        for n=1:n_max;
            Pd=PBig(1:n,1:n);
            His=[His,log(cond(Pd))];
            Fd=PBig(1:n,end);
            PB=PBig([1:n,n_max+1],[1:n,n_max+1]);
            index=(log(cond(PB))-max(His));
            if (index>index_max)
                nb=n;
                Pb=Pd;
                Fb=Fd;
                index_max=index;
                lanb=lan;
            end
        end
    end
    %moving Window
    PWin=Pr(1).PWin;
    if id>n_max
        yb=Y(id-n_max);
        ub=(U(id-n_max,:))';
    else
        yb=yd*0;
        ub=ud*0;
    end
    
    PWin=PWin+[ud;yd]*[ud;yd]'-[ub;yb]*[ub;yb]';
    Pr(1).PWin=PWin;
    for n=1:n_max;
        Pd=PWin(1:n,1:n);
        Fd=PWin(1:n,end);
        PB=PWin([1:n,n_max+1],[1:n,n_max+1]);
        index=(log(cond(PB))-log(cond(Pd)));
        if (index>index_max)
            nb=n;
            Pb=Pd;
            Fb=Fd;
            index_max=index;
            lanb=1;
        end
    end
    IND(id)=index_max;
    TET=(Pb+0.1*eye(length(Pb)))\Fb;		%Poles of signal
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    ut = ud(1:nb);
    N_ = length(TET);
    z = [];
    tet = [1,-TET'];
    for i=1:N_
        l = 0;
        c_tet = 1;
        for j=i:-1:1
            l = l + ut(j)*tet(c_tet);
            c_tet = c_tet + 1;
        end
        z = [ z ; l ];
    end
    Zeros = roots([z' , 0]);				%Zeros of Signal
    [ Res , FerZ , k ] = residuez([z',0],[1,-TET']);		
    for i=1:length(FerZ)		%Normalize of discrete Poles
        FerZ(i) = FerZ(i)/abs(FerZ(i));
    end
    FerS_ = (1/dt)*log(FerZ);		%Convert to Continouse Poles 
    C = abs(Res);
    %%%%%%%%%%%%%%%%%%%%%
    Fer=roots([1,-TET']);
    FerS=imag(1/dt*log(Fer))
    %%%%%%%%%%%%%%%%%%%%%
    Fer = FerS_;
    M = fi(Fer);
    nx = numerictype(1,20,13);
    F_q = double(quantize(M,nx,'Nearest'));        %%% quantize frequencies
    [ l , ~ ] = size(F_q);
    pic = zeros(length(0:ddt:len_Pic),1);
    for r=1:l
        pic(length(0:ddt:abs(Fer(r)))) = pic(length(0:ddt:abs(Fer(r)))) + C(r);
    end
    PIC = [ PIC , pic ];
    %%%%%%%%%%%%%%%%%
    plot(ones(length(FerS),1)*id,FerS,'*')		%Drowing Picture each Step
    hold on
    drawnow
    N(n_max+1+id)=nb;
    L(n_max+1+id)=lanb;
end


figure
% plot((a0+b0+c0+d0)*2+1,'g')
% hold on
plot(N,'r')
% hold on
figure()
plot(L,'k')

figure
plot(LAN)
title('Landa')
figure
a=image(uint8(PIC(:,1:end)),'CDataMapping','scaled');
Ytick = linspace(0,pi/dt,10);
set(gca,'Ydir','normal' );
r=gca;
u = num2cell(uint8(linspace(0,pi/dt,length(r.YTickLabel)+1)));
u(1)=[];
r.YTickLabel = u;
% u = num2cell(uint8(linspace(0,t_final,length(r.XTickLabel)+1)));
% u(1)=[];
% r.XTickLabel = u;
% xlabel Time(s)
ylabel Freq
colorbar


