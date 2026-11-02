% Dynamics Coursework

close all;
clear;
clc;

% Define constants (all units in m, N, kg)
Tnlim = [0.35,0.7];
ulim = 0.03;
tensilelim = [15,24.1]*(10^6);
shearlim = tensilelim/sqrt(3);
E = [1.07,2.03]*(10^9);
rho = [1.04,1.23]*(10^3);
psi = 0;
n = 2;

% Define constants and calculate beta, f0 and w
x = [0.01,0.01,0.02];
f = [1.5,2.5,3.5];
w = 2*pi*f;
f0 = x.*(w.^2);

%% Initiate & Iteration to determine Tn, umax and maxTensileStress

% Create vectors to iterate through
db = 0.001;
dd = 0.001;
dh = 0.001;
b = 0.015:db:0.042;
d = 0.002:dd:0.003;
h = 0.133:dh:0.133;

% Initiate matrices to store Tn, umax, tensileStress, mtotal
Tn = NaN(length(b),length(b),length(d),length(d),length(h),length(E));
umax = Tn;
maxTensileStress = Tn;
maxShearStress = Tn;
buckling = Tn; 
sumd = 0; 
sumb = 0; 
for i = length(d):-1:1
    sumd = sumd+i; 
end
for i = length(b):-1:length(b)-10
    sumb = sumb+i; 
end
iterationNumber = sumb*sumd*length(h)*length(E);

l = 0;
for i1 = 1:11
    for i2 = i1:length(b)
        for i3 = 1:length(d)
            for i4 = i3:length(d)
                for i5 = 1:length(h)
                    for i6 = 1:length(E)
                        % Show progress
                        clc;
                        l = l+1;
                        fprintf('PART 1 --- Loading: %.2f/100\n',l/iterationNumber*100);
                        
                        % Calculate the natural period
                        fA = @(x)((b(i1)+(b(i2)-b(i1))/h(i5)*x).*(d(i3)+(d(i4)-d(i3))/h(i5)*x));
                        m = (n*0.25)+(pi*(0.025^2)*0.01+0.05*0.025*0.017-0.01*2*pi*(0.007^2)/4)*rho(i6)/2+(integral(fA,0,h(i5)/2))*rho(i6);
                        L = h(i5);
                        fK = @(x)(((L-x).^2)./((b(i2)-(b(i2)-b(i1))/h(i5)*x).*((d(i4)-(d(i4)-d(i3))/h(i5)*x).^3)/12));
                        K = E(i6)/integral(fK,0,L);
                        wn = sqrt(K/m);
                        Tn(i1,i2,i3,i4,i5,i6) = 2*pi/wn;
                        
                        % Continue if Tn does not meet criteria
                        if Tn(i1,i2,i3,i4,i5,i6)>=Tnlim(2) && Tn(i1,i2,i3,i4,i5,i6)<=Tnlim(1)
                            continue;
                        end
                        
                        % Calculate umax for all 3 cases
                        beta = w./wn;
                        ust = f0./(wn.^2);
                        R = 1./sqrt(((1-beta.^2).^2)+((2*psi*beta).^2));
                        umax1 = ust.*R;
                        umax(i1,i2,i3,i4,i5,i6) = max(umax1);
                        
                        % Continue if umax is larger than ulim
                        if umax(i1,i2,i3,i4,i5,i6)>=ulim
                            continue;
                        end
                        
                        % Calculate the maximum tensile stress
                        Vmax = K*max([ust,umax(i1,i2,i3,i4,i5,i6)]);
                        fy = @(x)((d(i4)-(d(i4)-d(i3))/h(i5)*x)/2);
                        fI = @(x)((b(i2)-(b(i2)-b(i1))/h(i5)*x).*((d(i4)-(d(i4)-d(i3))/h(i5)*x).^3)/12);
                        fM = @(x)(Vmax*(L-x));
                        x1 = 0:0.001:L;
                        tensileStress = fM(x1).*fy(x1)./fI(x1);
                        maxTensileStress(i1,i2,i3,i4,i5,i6) = max(tensileStress);
                        
                        % Calculate the maximum shear stress
                        maxShearStress(i1,i2,i3,i4,i5,i6) = Vmax/fA(0);
                        
                        % Calculate the weight
                        f1 = @(x)(((b(i1)+(b(i2)-b(i1))/h(i5)*x).*((d(i3)+(d(i4)-d(i3))/h(i5)*x).^3)/12).*((sin(pi*x/h(i5)/2)).^2));
                        f2 = @(x)((cos(pi*x/h(i5)/2)).^2);
                        P = E(i6)*(pi^2)/(h(i5)^2)*integral(f1,0,h(i5))/integral(f2,0,h(i5))/4;
                        if (m*9.81)>=P
                            buckling(i1,i2,i3,i4,i5,i6) = 0;
                        else
                            buckling(i1,i2,i3,i4,i5,i6) = 1;
                        end
                    end
                end
            end
        end
    end
end

%%

err = 0.0001;
merr = 0.001; 
tensileFOS1 = 1.2;
shearFOS1 = 1;
tensileFOS2 = 1.2;
shearFOS2 = 1;
mtotal = NaN(length(b),length(b),length(d),length(d),length(h),length(E));

l = 0;
for i1 = 1:16
    for i2 = i1:length(b)
        for i3 = 1:length(d)
            for i4 = i3:length(d)
                for i5 = 1:length(h)
                    for i6 = 1:length(E)
                        % Show progress
                        clc;
                        l = l+1;
                        fprintf('PART 2 --- Loading: %.2f/100\n',l/iterationNumber*100);
                        
                        % Check if the chosen dimensions without error meet
                        % the criteria
                        if Tn(i1,i2,i3,i4,i5,i6)<Tnlim(2) && Tn(i1,i2,i3,i4,i5,i6)>Tnlim(1) && umax(i1,i2,i3,i4,i5,i6)<ulim && tensileFOS1*maxTensileStress(i1,i2,i3,i4,i5,i6)<tensilelim(i6) && shearFOS1*maxShearStress(i1,i2,i3,i4,i5,i6)<shearlim(i6) && buckling(i1,i2,i3,i4,i5,i6)==1
                            if err==0
                                % Calculate the total mass
                                fA = @(x)((b(i1)+(b(i2)-b(i1))/h(i5)*x).*(d(i3)+(d(i4)-d(i3))/h(i5)*x));
                                mtotal(i1,i2,i3,i4,i5,i6) = (pi*(0.025^2)*0.01+0.05*0.025*0.017-0.01*2*pi*(0.007^2)/4+0.01*pi*(0.055^2)-4*0.01*pi*(0.01^2)/4)*rho(i6)/2+integral(fA,0,h(i5))*rho(i6);
                            else
                                % Initiate empty vectors
                                TnC = zeros(1,1+2^5);
                                umaxC = TnC;
                                maxTensileStressC = TnC;
                                maxShearStressC = TnC;
                                bucklingC = TnC; 
                                TnC(1) = Tn(i1,i2,i3,i4,i5,i6);
                                umaxC(1) = umax(i1,i2,i3,i4,i5,i6);
                                maxTensileStressC(1) = maxTensileStress(i1,i2,i3,i4,i5,i6);
                                maxShearStressC(1) = maxShearStress(i1,i2,i3,i4,i5,i6);
                                bucklingC(1) = buckling(i1,i2,i3,i4,i5,i6); 
                                
                                % Calculate Tn, umax and maxTensileStress
                                % considering error
                                k = 1;
                                for j1 = -1:2:1
                                    for j2 = -1:2:1
                                        for j3 = -1:2:1
                                            for j4 = -1:2:1
                                                for j5 = -1:2:1
                                                    k = k+1;
                                                    
                                                    % Calculate the natural period
                                                    fA = @(x)(((b(i1)+j1*err)+((b(i2)+j2*err)-(b(i1)+j1*err))/(h(i5)+j5*err)*x).*((d(i3)+j3*err)+((d(i4)+j4*err)-(d(i3)+j3*err))/(h(i5)+j5*err)*x));
                                                    m = (n*0.25)+(pi*(0.025^2)*0.01+0.05*0.025*0.017-0.01*2*pi*(0.007^2)/4)*rho(i6)/2+(integral(fA,0,(h(i5)+j5*err)/2))*rho(i6);
                                                    L = h(i5)+j5*err;
                                                    fK = @(x)(((L-x).^2)./(((b(i2)+j2*err)-((b(i2)+j2*err)-(b(i1)+j1*err))/(h(i5)+j5*err)*x).*(((d(i4)+j4*err)-((d(i4)+j4*err)-(d(i3)+j3*err))/(h(i5)+j5*err)*x).^3)/12));
                                                    K = E(i6)/integral(fK,0,L);
                                                    wn = sqrt(K/m);
                                                    TnC(k) = 2*pi/wn;
                                                    
                                                    % Calculate umax for all 3 cases
                                                    beta = w./wn;
                                                    ust = f0./(wn.^2);
                                                    R = 1./sqrt(((1-beta.^2).^2)+((2*psi*beta).^2));
                                                    umax1 = ust.*R;
                                                    umaxC(k) = max(umax1);
                                                    
                                                    % Calculate the maximum tensile stress
                                                    Vmax = K*max([ust,umaxC(k)]);
                                                    fy = @(x)(((d(i4)+j4*err)-((d(i4)+j4*err)-(d(i3)+j3*err))/(h(i5)+j5*err)*x)/2);
                                                    fI = @(x)(((b(i2)+j2*err)-((b(i2)+j2*err)-(b(i1)+j1*err))/(h(i5)+j5*err)*x).*(((d(i4)+j4*err)-((d(i4)+j4*err)-(d(i3)+j3*err))/(h(i5)+j5*err)*x).^3)/12);
                                                    fM = @(x)(Vmax*(L-x));
                                                    x1 = 0:0.001:L;
                                                    tensileStress = fM(x1).*fy(x1)./fI(x1);
                                                    maxTensileStressC(k) = max(tensileStress);
                                                    
                                                    % Calculate the maximum shear stress
                                                    maxShearStressC(k) = Vmax/fA(0);
                                                    
                                                    % Calculate the weight
                                                    f1 = @(x)((((b(i1)+j1*err)+((b(i2)+j2*err)-(b(i1)+j1*err))/(h(i5)+j5*err)*x).*(((d(i3)+j3*err)+((d(i4)+j4*err)-(d(i3)+j3*err))/(h(i5)+j5*err)*x).^3)/12).*((sin(pi*x/(h(i5)+j5*err)/2)).^2));
                                                    f2 = @(x)((cos(pi*x/(h(i5)+j5*err)/2)).^2);
                                                    P = E(i6)*(pi^2)/((h(i5)+j5*err)^2)*integral(f1,0,(h(i5)+j5*err))/integral(f2,0,(h(i5)+j5*err))/4;
                                                    if (m*9.81)>=P
                                                        bucklingC(k) = 0;
                                                    else
                                                        bucklingC(k) = 1;
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                                
                                if all(TnC<Tnlim(2)) && all(TnC>Tnlim(1)) && all(umaxC<ulim) && all(tensileFOS2*maxTensileStressC<tensilelim(i6)) && all(shearFOS2*maxShearStressC<shearlim(i6)) && all(bucklingC)
                                    % Calculate the total mass
                                    fA = @(x)((b(i1)+(b(i2)-b(i1))/h(i5)*x).*(d(i3)+(d(i4)-d(i3))/h(i5)*x));
                                    mtotal(i1,i2,i3,i4,i5,i6) = (pi*(0.025^2)*0.01+0.05*0.025*0.017-0.01*2*pi*(0.007^2)/4+0.01*pi*(0.055^2)-4*0.01*pi*(0.01^2)/4)*rho(i6)/2+integral(fA,0,h(i5))*rho(i6);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
fprintf('Number of solutions = %i\n',nnz(~isnan(mtotal)));

%% Find the dimensions of the structure with the minimum mass

% Determine the index to the minimum mass
[A1,i6] = min(mtotal,[],6);
[A2,i5] = min(A1,[],5);
[A3,i4] = min(A2,[],4);
[A4,i3] = min(A3,[],3);
[A5,i2] = min(A4,[],2);
[mmin,i1] = min(A5,[],1);
i2 = i2(i1);
i3 = i3(i1,i2);
i4 = i4(i1,i2,i3);
i5 = i5(i1,i2,i3,i4);
i6 = i6(i1,i2,i3,i4,i5);
fprintf('The minimum mass is %fkg: b1=%.3fm b2=%.3fm d1=%.3fm d2=%.3fm h=%.3fm material=%i n=%i\n',mmin,b(i1),b(i2),d(i3),d(i4),h(i5),i6,n);

% Calculate the natural period
fA = @(x)((b(i1)+(b(i2)-b(i1))/h(i5)*x).*(d(i3)+(d(i4)-d(i3))/h(i5)*x));
m = (n*0.25)+(pi*(0.025^2)*0.01+0.05*0.025*0.017-0.01*2*pi*(0.007^2)/4)*rho(i6)/2+(integral(fA,0,h(i5)/2))*rho(i6);
L = h(i5);
fK = @(x)(((L-x).^2)./((b(i2)-(b(i2)-b(i1))/h(i5)*x).*((d(i4)-(d(i4)-d(i3))/h(i5)*x).^3)/12));
K = E(i6)/integral(fK,0,L);
wn = sqrt(K/m);
TnS = 2*pi/wn;
fprintf('The natural period is %.4fs. \n',TnS);

% Calculate umax for all 3 cases
beta = w./wn;
ust = f0./(wn.^2);
R = 1./sqrt(((1-beta.^2).^2)+((2*psi*beta).^2));
umax1 = ust.*R;
umaxS = max(umax1);
fprintf('The maximum dynamic displacement for load case %i is %.6f. \n',[[1,2,3];umax1]);
fprintf('The maximum displacement is %.6f. \n',max([ust,umaxS]));

% Calculate the maximum tensile stress
Vmax = K*max([ust,umaxS]);
fy = @(x)((d(i4)-(d(i4)-d(i3))/h(i5)*x)/2);
fI = @(x)((b(i2)-(b(i2)-b(i1))/h(i5)*x).*((d(i4)-(d(i4)-d(i3))/h(i5)*x).^3)/12);
fM = @(x)(Vmax*(L-x));
x1 = 0:0.001:L;
tensileStress = fM(x1).*fy(x1)./fI(x1);
maxTensileStressS = max(tensileStress);
fprintf('The maximum tensile stress is %.2f. \n',maxTensileStressS);
tensileFOS = tensilelim(i6)/maxTensileStressS;
fprintf('The factor of safety for tensile strength is %.3f. \n',tensileFOS);

% Calculate the maximum shear stress
maxShearStressS = Vmax/fA(0);
fprintf('The maximum shear stress is %.2f. \n',maxShearStressS);
shearFOS = shearlim(i6)/maxShearStressS;
fprintf('The factor of safety for shear strength is %.3f. \n',shearFOS);

% Calculate buckling load
f1 = @(x)(((b(i1)+(b(i2)-b(i1))/h(i5)*x).*((d(i3)+(d(i4)-d(i3))/h(i5)*x).^3)/12).*((sin(pi*x/h(i5)/2)).^2));
f2 = @(x)((cos(pi*x/h(i5)/2)).^2);
P = E(i6)*(pi^2)/(h(i5)^2)*integral(f1,0,h(i5))/integral(f2,0,h(i5))/4;
fprintf('The weight is %.2f and the buckling load is %.2f. \n',m*9.81,P); 