%Liel Giat
%ID - 313612095

%Part 1

n = -1000:1000;
a_n = rectangularPulse(-100,100,n);% Function Of Rectangular Pulse
a_n (1101) = 0;% we want that in the edges will be 0, not 0.5 (the function fading slowly)
a_n (901) = 0;
f1 = figure;
stem(n,a_n,'g');% plot of discrete signal
title('Rectangular Pulse Function Of a_n');
ylabel('a_n');
xlabel('n (Time)');

%Part 2

[a_k] = Fourier_coefficients(n,a_n);% call the function that calculate the coefficients of fourier
f2 = figure;
plot(n,a_k); %plot of continuous signal
title('Fourier Coefficients (a_k)');
ylabel('a_k');
xlabel('n');

%Part 3
% If we divide The imaginary of a_k by the real of a_k and the result will
% be less then 1 in a million - a_k is real
if var(imag(a_k))/var(real(a_k)) < 1/100000000
    fprintf ('a_k is real')
end
% I also will show that a_k is real with a plot 
% (compare to the complex conjugate)
k = -1000:1000;
nk = n'*k;% nk = n * k
N = length(n);
ex = exp(1i*2*pi*nk/N);% this is the complex conjugate (+i ,not -i)
% ex - we will also use it for the inverse transform
C_c = (a_n./N)*ex;%Complex conjugate, so we will can to check if a_k real
R_n = flip(a_k);%Reverse of a_k, so we will can to check symmetry
f3 = figure;

% the complex conjugate of a_k
subplot(2,2,1);
plot (n,C_c,'b');
title('Complex conjugate');
ylabel('C_c');
xlabel('n');

% the complex conjugate of a_k - a_k => need to be zero in the plot
subplot(2,2,2);
Real = C_c - a_k;
plot (n,Real,'c');
title('if all 0 - a_k is real');
ylabel('C_c - a_k');
xlabel('n');

% the reverse of a_k
subplot(2,2,3);
plot (n,R_n,'r');
title('Reverse Of a_k');
ylabel('R_n');
xlabel('n');

% the reverse of a_k - a_k => need to be zero in the plot
subplot(2,2,4);
symmetri = R_n - a_k;
plot (n,symmetri,'g');
title('if all 0 - a_k is symmetrical');
ylabel('R_n - a_k');
xlabel('n');

%Part 4

%From the analytical formula
a_k2 = (exp(-1i*2*pi*k*100/N) - exp(1i*2*pi*k*99/N)) ./ (N*(exp(-1i*2*pi*k/N)-1));
f4 = figure;
plot (n,a_k2,'--r*');
hold on %function that help us to show 2 different Y on the same plot
plot(n,a_k,'.-g');% what we calaculate in part 2
legend('a_k2','a_k')% show us which line of who
hold off
title('Fourier Coefficients (a_k/a_k2)');
ylabel('a_k and a_k2');
xlabel('n');

%Part 5

% from the formula on the page
b_k = a_k.*exp(-1i*2*pi*150*k/N);
b_n = b_k*ex; % the inverse transform
f5 = figure;
subplot(2,1,1); % function that help us to show how many plots that we want in the same figure
stem(n,b_n,'c');
title('Time Shifting of b_n');
ylabel('b_n');
xlabel('n');

% I shift a_n 150 places forward (like b_n) to check if they the same
a_n2 = rectangularPulse(50,250,n);
a_n2 (1251) = 0;
a_n2 (1051) = 0;
subplot(2,1,2);
stem(n,a_n2,'g');
title('Time Shifting of a_n');
ylabel('a_n');
xlabel('n');

%Part 6

% from the formula on the page
c_k = a_k.*(1-exp(-1i*k*2*pi/N));
c_n = c_k*ex;% the inverse transform
f6 = figure;
stem(n,c_n,'m');
title('Time Derivative');
ylabel('c_n');
xlabel('n');

%Part 7

% from the formula on the page
d_k = (a_k.^2).*N;
d_n = d_k*ex;% the inverse transform
f7 = figure;
stem(n,d_n,'r');
title('Time Convolution');
ylabel('d_n');
xlabel('n');

%Part 8

% from the formula on the page, the left side of Parseval
A1 = abs(d_k);
A2 = (A1.^2);
A = sum(A2);
% from the formula on the page, the right side of Parseval
B1 = abs(d_n);
B2 = (B1.^2);
B = sum(B2)*(1/N);
% checking if thy equal each other
if (round (A/B) == 1)
    fprintf('Parseval Works');
else
    fprintf('Parseval does not work');
end

%Part 9

% from the formula on the page
e_n = a_n .* b_n;
[e_k] = Fourier_coefficients(n,e_n);% the function that we use in part 2 (fourier transform)
f8 = figure;
subplot(2,2,1);
plot(n,abs(e_k),'k');
title('Frequency Convolution');
ylabel('e_k - Amplitude');
xlabel('n');

subplot(2,2,2);
plot(n,angle(e_k),'k');
title('Frequency Convolution');
ylabel('e_k - Phase');
xlabel('n');

f_k = zeros(1,N);% we put 0 in all f_k
% 3 loops that run all over k and l
for k = 1:N % in each k we need to calaulate the formula on the page
    for L = 1:N % in each k we need to calculate with all the values of L
        shft = k - L + 1; % the shift k-L in the formula
        if(shft <= 0)% the shft must be positive
            shft = N + shft;
        end
        f_k(k) = f_k(k) + a_k(L)*b_k(shft); % the sigma of f_k
    end
end
f_k = circshift(f_k,1001);%move each value in f_k 1001 places
subplot(2,2,3);
plot(n,abs(f_k),'m');
title('Cyclic Convolution');
ylabel('f_k - Amplitude');
xlabel('n');

subplot(2,2,4);
plot(n,angle(f_k),'m');
title('Cyclic Convolution');
ylabel('f_k - Phase');
xlabel('n');
k = -1000:1000;% the value of k is changing during the 'for' (this is for the rest of the parts)

%Part 10

% from the formula on the page
g_n = a_n.*cos(2*pi*500*(n/N));
[g_k] = Fourier_coefficients(n,g_n);% the function that we use in part 2 (fourier transform)
f9 = figure;
plot (n,g_k);
title('Doubling In Cos In Time');
ylabel('g_k');
xlabel('n');

%Part 11

% from the formula on the page
N2 = 10005;
f_k  = zeros(1,N2);% we change all the values in f_k to 0 and also extend f_k
%we start in the third element to compare f_k to a_k and then we will
%compare every fifth elemnt in f_k to a_k untili the end
f_k (3:5:end) = a_k; 
f_k = (1/5)*f_k;% from the formula on the page
n2 = -5002:5002;
k2 = -5002:5002;
nk2 = n2'*k2;
ex2 = exp(1i*2*pi*nk2/N2);% nk2 = n * k
f_n = f_k*ex2;
f10 = figure;
stem(n2,f_n,'r');
title('Shrink Time');
ylabel('f_n');
xlabel('n2');


%Part 12

% 6 options of M between 700 - 1000
M1 = -1000:1000;
M2 = -900:900;
M3 = -850:850;
M4 = -800:800;
M5 = -750:750;
M6 = -700:700;
% M act here like k
nM1 = M1'*n;
nM2 = M2'*n;
nM3 = M3'*n;
nM4 = M4'*n;
nM5 = M5'*n;
nM6 = M6'*n;
% The inverse transform
a_m1 = a_k*exp(1i*2*pi*nM1/N);
a_m2 = a_k(100:1900)*exp(1i*2*pi*nM2/N);% a_k(100:1900) = a_k(-900:900)
a_m3 = a_k(150:1850)*exp(1i*2*pi*nM3/N);% a_k(150:1850) = a_k(-850:850)
a_m4 = a_k(200:1800)*exp(1i*2*pi*nM4/N);% a_k(200:1800) = a_k(-800:800)
a_m5 = a_k(250:1750)*exp(1i*2*pi*nM5/N);% a_k(250:1750) = a_k(-750:750)
a_m6 = a_k(300:1700)*exp(1i*2*pi*nM6/N);% a_k(300:1700) = a_k(-700:700)

f11 = figure;

subplot(3,2,1);
stem(n,a_m1);
title('-1000:1000');
ylabel('a_m1');
xlabel('n');

subplot(3,2,2);
stem(n,a_m2);
title('-900:900');
ylabel('a_m2');
xlabel('n');

subplot(3,2,3);
stem(n,a_m3);
title('-850:850');
ylabel('a_m3');
xlabel('n');

subplot(3,2,4);
stem(n,a_m4);
title('-800:800');
ylabel('a_m4');
xlabel('n');

subplot(3,2,5);
stem(n,a_m5);
title('-750:750');
ylabel('a_m5');
xlabel('n');

subplot(3,2,6);
stem(n,a_m6);
title('-700:700');
ylabel('a_m6');
xlabel('n');

%Part 13

%A
h_n = a_n.*sin(2*pi*500*(n/N));% from the formula on the page
 f12 = figure;
subplot(2,2,[1 2]);
stem(n,h_n,'g');
title('Doubling In sin In Time');
ylabel('h_n');
xlabel('n');

%B
[h_k] = Fourier_coefficients(n,h_n);% the function that we use in part 2 (fourier transform)
subplot(2,2,3);
plot(n,real(h_k),'c');
title('Real numbers');% plot of the real numbers
ylabel('h_k');
xlabel('n');

subplot(2,2,4);
plot(n,imag(h_k),'r');% plot of the imaginary numbers
title('Imaginary numbers');
ylabel('h_k');
xlabel('n');

%C
k = -999:999;% (N-2)/2 = 999
n = -999:999;
H_k = -1i*sign(k);% from the formula on the page
g_k = g_k(2:2000); % we need the same  dimensions for matrix like k
hTilda_k = H_k.*g_k;% from the formula on the page
f13 = figure;
subplot(2,2,1);
plot(n,real(hTilda_k));% plot of the real numbers
title('Real numbers');
ylabel('hTilda_k');
xlabel('n');

subplot(2,2,2);
plot(n,imag(hTilda_k),'m');
title('Imaginary numbers');% plot of the imaginary numbers
ylabel('hTilda_k');
xlabel('n');

%D
N = 1999;
hTilda_n = hTilda_k*exp(1i*2*pi*(n'*k)/N);% from the formula on the page
subplot(2,2,[3 4]);
stem(n,hTilda_n,'k');
title('Hilbert transform');
ylabel('hTilda_n ');
xlabel('n');


% fourier transform of a_n
% input - n and a_n , output - a_k
function [a_k] = Fourier_coefficients (n,a_n)
k = -1000:1000;
nk = n'*k;   
N = length(n); % = 2001
ex1 = exp(-1i*2*pi*nk/N);% nk = n * k
a_k = (a_n./N)*ex1;
end
