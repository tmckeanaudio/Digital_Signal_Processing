%% *EECE5666 (DSP) : Homework-7*
% *Due on April 19, 2022 by 11:59 pm via submission portal.* 
% 
% *NAME*: McKean, Tyler
% 
% 
%% *Instructions*
%% 
% # You are required to complete this assignment using Live Editor.
% # Enter your MATLAB script in the spaces provided. If it contains a plot, 
% the plot will be displayed after the script.
% # All your plots must be properly labeled and should have appropriate titles 
% to get full credit.
% # Use the equation editor to typeset mathematical material such as variables, 
% equations, etc. 
% # After completeing this assignment, export this Live script to PDF and submit 
% the PDF file through the provided submission portal.
% # You will have only one attempt to submit your assignment. Make every effort 
% to submit the correct and completed PDF file the first time.
% # Please submit your homework before the due date/time. A late submission 
% after midnight of the due date will result in loss of points at a rate of 10% 
% per hour until 8 am the following day, at which time the solutions will be published.
%% *Default Plot Parameters*

set(0,'defaultfigurepaperunits','points','defaultfigureunits','points');
set(0,'defaultaxesfontsize',10); set(0,'defaultaxeslinewidth',1.5);
set(0,'defaultaxestitlefontsize',1.4,'defaultaxeslabelfontsize',1.2);
%% 
% 
%% Problem 7.1
% We want to design a lowpass analog Chebyshev-I filter that has a 0.5 dB or 
% better ripple at 10 Hz and at least 45 dB of attenuation at 20 Hz.

clc; close all; clear;
%% 
% *(a)* Using the design procedure on Page 640 of the textbook (or that in Example 
% 11.3) obtain the system function in a rational function form.
% 
% *Solution:* Follow the following steps. Perform numerical calculations using 
% MATLAB below each step.
% 
% *Step-0*: Determine the analog passband ripple parameter $\epsilon$ and stopband 
% attenuation parameter $A.$
% 
% We can determine the values of parameters $\epsilon$ and $A$ by
% 
% $$\epsilon = \sqrt{10^{Ap(0.1)} -1} = \sqrt{10^{(0.5)(0.1)} -1} = 0.3493$$
% 
% $$A = 10^{As(0.05)} = 10^{(45)(0.05)} = 177.83$$

Ap = 0.5; As = 45; Omegap = 2*pi*10; Omegas = 2*pi*20;
epsilon = sqrt(10^(Ap/10)-1), A = 10^(As/20)
%% 
% *Step-1*: Compute the parameters $\alpha$ and $\beta$ using (11.50):
% 
% Using (11.50), we calculate the values of $\alpha$ and $\beta$ to be
% 
% $$\alpha = \frac{\Omega_s}{\Omega_p} = \frac{2\pi(20)}{2\pi(10)} = 2$$
% 
% $$\beta = \frac{1}{\epsilon}\sqrt{A^2 -1} = \frac{1}{0.3493}\sqrt{(177.83)^2 
% - 1} = 509.0734$$

alpha = Omegas/Omegap, beta = sqrt((A^2 -1)/epsilon^2)
%% 
% *Step-2:* Compute order $N$ using (11.49) and round upwards to the nearest 
% integer:
% 
% $$N \geq \frac{\ln(\beta + \sqrt{\beta^2 -1})}{\ln(\alpha + \sqrt{\alpha^2 
% -1})}  = \frac{\ln((509) + \sqrt{(509)^2 -1})}{\ln(2 + \sqrt{2^2 -1})} $$

N = ceil(log(beta + sqrt(beta^2 - 1))/(log(alpha + sqrt(alpha^2 - 1))))
%% 
% *Step-3:* Set $\Omega_\mathrm{c}=\Omega_\mathrm{p}$ and compute $a$ and $b$ 
% using (11.44) and (11.45):
% 
% We have a passband frequency requirement of 10Hz, thus 
% 
% $$\Omega_c = \Omega_p \rightarrow 20\pi \frac{rad}{sec}$$
% 
% We must first calculate the value of $\gamma$ in order to compute values $a$ 
% and $b$
% 
% $$\gamma ={\left(\frac{1}{\epsilon }+\sqrt{1+\frac{1}{\epsilon^2 }}\right)}^{\frac{1}{N}} 
% ={\left(\frac{1}{0\ldotp 3493}+\sqrt{1+\frac{1}{{\left(0\ldotp 3493\right)}^2 
% }}\right)}^{\frac{1}{6}} =1\ldotp 3441\;$$ 
% 
% $$\begin{array}{l}a=\frac{1}{2}\left(\gamma -\gamma^{-1} \right)=\frac{1}{2}\left(1\ldotp 
% 3441-1\ldotp {3441}^{-1} \right)=0\ldotp 3\\b=\frac{1}{2}\left(\gamma +\gamma^{-1} 
% \right)=\frac{1}{2}\left(1\ldotp 3441+1\ldotp {3441}^{-1} \right)=1\ldotp 044\end{array}$$

gamma = (1/epsilon + sqrt(1+1/(0.3493)^2))^(1/N)
a = 0.5*(gamma - 1/gamma), b = 0.5*(gamma + 1/gamma)
%% 
% *Step-4:* Compute the pole locations using (11.41) and (11.42):
% 
% We compute the pole locations using the formulas
% 
% $\sigma_k =\left\lbrack \Omega_c \sinh \left(\phi \right)\right\rbrack \cos 
% \left(\theta_k \right)$ with $a\;$$\triangleq$ $\sinh \left(\phi \right)$
% 
% $\Omega_k =\left\lbrack \Omega_c \cosh \left(\phi \right)\right\rbrack \sin 
% \left(\theta_k \right)$ with $b\;$$\triangleq$ $\cosh \left(\phi \right)$ with 
% the angle being
% 
% $$\theta_k =\frac{\pi }{2}+\frac{2k-1}{2N}\pi ,\;\;\;\textrm{with}\;k=1,2,\ldotp 
% \ldotp \ldotp ,2N$$
% 
% Substituting $\Omega_{c\;} ,\;a,\;\textrm{and}\;b\;\textrm{into}\;\textrm{the}\;\textrm{first}\;\textrm{formulas}\;\textrm{we}\;\textrm{get}$:
% 
% $$\sigma_k =\;$$$$20\pi \left(0\ldotp 3\right)\cos \left(\theta_k \right)$$
% 
% $$\Omega_k =20\pi \left(1\ldotp 044\right)\sin \left(\theta_k \right)$$

Omegac = Omegap;
% Step-4: Calculations of Poles
k = 1:N; thetak = pi/2+(2*k-1)*pi/(2*N);
sigmak = (a*Omegac)*cos(thetak); Omegak = (b*Omegac)*sin(thetak);
sk = cplxpair(sigmak + 1j*Omegak)
%% 
% Taking only the first 6 poles calculated are the left-half elliptic shaped 
% poles of the Chebyshev-I filter, which are
% 
% $$\begin{array}{l}s_1 =-4\ldotp 8789+\textrm{j63}\ldotp 3635\\s_2 =-13\ldotp 
% 3295+\textrm{j46}\ldotp 3853\\s_3 =-18\ldotp 2085+\textrm{j16}\ldotp 9782\\s_4 
% =-18\ldotp 2085-\textrm{j16}\ldotp 9782\\s_5 =-13\ldotp 3295-\textrm{j46}\ldotp 
% 3853\\s_6 =-4\ldotp 8789-\textrm{j63}\ldotp 3635\end{array}$$
% 
% *Step-5:* Compute the filter gain $G$ and the system function $H(\j\Omega)$ 
% from (11.43): Since $N$ is even, 
% 
% We can compute the $\textrm{filter}\;\textrm{gain}\;G\;\textrm{by}\;\textrm{the}\;\textrm{product}\;\textrm{of}\;\textrm{the}\;\textrm{first}\;N\;\textrm{poles}\;\textrm{by}\;\left(\frac{1}{\sqrt{1+\epsilon^2 
% }}\right),\textrm{which}\;\textrm{yields}:$

Rp = 1/sqrt(1+epsilon^2);
D = real(poly(sk));
G = D(end)*Rp 
%% 
% Thus, $\textrm{the}\;\textrm{filter}\;\textrm{gain}\;G=5\ldotp 5046$$\mathrm{e}{+09}$
% 
% Finally, determine the system function:
% 
% $$H_C \left(s\right)=\frac{5\ldotp 5046\textrm{e09}}{\left(s+4\ldotp 8789-\textrm{j63}\ldotp 
% 3635\right)\left(s+13\ldotp 3295-\textrm{j46}\ldotp 3853\right)\left(s+18\ldotp 
% 2085-\textrm{j16}\ldotp 9782\right)\left(s+18\ldotp 2085+\textrm{j16}\ldotp 
% 9782\right)\left(s+13\ldotp 3295+\textrm{j46}\ldotp 3853\right)\left(s+4\ldotp 
% 8789+\textrm{j63}\ldotp 3635\right)}$$

s = tf('s');
Hc = G/((s + abs(sk(1)))*(s + abs(sk(2)))*(s + abs(sk(3)))*(s + abs(sk(4)))*(s + abs(sk(5)))*(s + abs(sk(6))))
%% 
% Thus, the system function $H\left(j\Omega \right)\;\textrm{is}\;\textrm{described}\;\textrm{by}\;\textrm{the}\;\textrm{rational}\;\textrm{function}\;\textrm{above}$
% 
% 
% 
% *(b)* Verify your design using the |*cheb1ord*| and |*cheby1*| functions.
% 
% *Solution:* 

[N,Omegac] = cheb1ord(Omegap,Omegas,Ap,As,'s')
Fc = Omegac/(2*pi)
%% 
% Thus, using the cheb1ord function yields the same result for $N = 6$ for the 
% order of the filter with the cutoff frequency of 
% 
% $F_c = 10$Hz

[C,D] = cheby1(N,Ap,Omegac,'s');
%% 
% 
% 
% *(c)* Provide plots of impulse response, amplitude response, log-magnitude 
% response in dB, and phase response in one figure using two rows and two columns.
% 
% *MATLAB script:* 

% Impulse Response
trsys = tf(C,D);
[h t] = impulse(trsys);
% Frequency Response
Fmax = 40; F = linspace(0,Fmax,1001);
omega = 2*pi*F;
H = freqs(C,D,omega);
Hmag = abs(H);
Hdb = mag2db(Hmag);
% Phase Response
Hpha = angle(H);
Hgd = -diff(unwrap(Hpha))./diff(omega); 
Hgd = [Hgd Hgd(end)];
% Amplitude Response
[Hr,wr] = zerophase(C,D,512,'whole');
% Plot Results
figure('Units','inches','Position',[0,0,12,4]);
% Impulse Response Plot
subplot(2,2,1),plot(t,h,'LineWidth',1.5),title('Impulse Response'),xlabel('Time Index t/s'),ylabel('\it{h(t)}')
xlim([0 t(end)])
% Magnitude Response Plot
subplot(2,2,2),plot(F,Hdb,'LineWidth',1.5),title('Log-Magnitude Response (dB)'), grid on 
ylim([-80 10]), yticks([-45 -0.5]), xticks([10 20])
xlabel('Frequency (Hz)'),ylabel('Decibels')
% Amplitude Response Plot
subplot(2,2,3), plot(wr/pi,Hr,'LineWidth',1.5), title('Amplitude Response')
xlabel('\omega/\pi'), ylabel('Amplitude'), ylim([0.85 1.025])
% Phase Response Plot
subplot(2,2,4), plot(F,Hpha,'LineWidth',1.5), title('Phase Delay Response'), xlabel('Frequency (Hz)')
ylabel('Radians')
%% 
% 
%% Problem 7.2
% *Text Problem 11.21 (Page 694)*
% 
% Consider a $9^{\mathrm{th}$-order analog Butterworth lowpass filter $H_\mathrm{c}(s)$ 
% with $3$-dB cutoff frequency of 10 Hz.
% 
% 
% 
% *(a)* Determine and graph pole locations of $H_\mathrm{c}(s).$
% 
% *Solution:* .

clc; close all; clear;
N = 9; Fc = 10; Omegac = 2*pi*Fc;
[C,D] = butter(N,Omegac,'s');
sk = roots(D)
%% 
% Thus, the 9 poles of the Butterworth lowpass filter are:
% 
% $$s_{1,9} =-10\ldotp 9101\pm \textrm{j61}\ldotp 8773$$
% 
% $$s_{2,8\;} =-31\ldotp 4159\pm \textrm{j54}\ldotp 414$$
% 
% $$s_{3,7} =-48\ldotp 132\pm \textrm{j40}\ldotp 3875$$
% 
% $$s_{4,6} =-59\ldotp 0426\pm \textrm{j21}\ldotp 489$$
% 
% $$s_5 =-62\ldotp 8319$$

figure, zplane(C,D), title('N = 9 Pole Locations of \it{H_c(s)}')
%% 
% 
% 
% *(b)* Plot the magnitude and log-magnitude responses over $[0,100]$ Hz range.
% 
% *MATLAB script:* 

% Frequency Response
Fmax = 100; F = linspace(0,Fmax,2001);
omega = 2*pi*F;
H = freqs(C,D,omega);
Hmag = abs(H);
Hdb = mag2db(Hmag);
% Plot Results
figure('Units','inches','Position',[0,0,12,4]);
% Magnitude Response Plot
subplot(1,2,1),plot(F,Hmag,'LineWidth',1.5),title('Magnitude Response'), grid on
xlim([0 100]), ylim([0 1.1]), yticks([0 0.707 1])
xlabel('Frequency (Hz)'),ylabel('Magnitude')
% Log-Magnitude Response Plot
subplot(1,2,2),plot(F,Hdb,'LineWidth',1.5),title('Log-Magnitude Response (dB)'), grid on
xlabel('Frequency (Hz)'),ylabel('Decibels'), xlim([0 100])
yticks([-40 -3])
%% 
% 
% 
% *(c)* Determine frequencies at which the attenuation is $30$ db, $40$ db, 
% and $50$ db.
% 
% *Solution:* 

figure
% Log-Magnitude Response Plot
plot(F,Hdb),title('Log-Magnitude Response (dB)'), grid on
xlabel('Frequency (Hz)'),ylabel('Decibels')
ylim([-55 -25]), yticks([-50 -40 -30]), xlim([11 22])
ax = gca;
chart = ax.Children(1);
datatip(chart,14.7,-30.1213,'Location','southwest');
datatip(chart,16.7,-40.0894,'Location','northeast');
datatip(chart,19,-50.1757,'Location','northeast');
%% 
% Observing the data points of the Log-Magnitude Response, we can see that the 
% corresponding frequencies are:
% 
% $$F_1 =14\ldotp 7\;\textrm{Hz}\;\textrm{with}\;\textrm{attenuation}\;\textrm{of}\;-30\textrm{dB}$$
% 
% $$F_2 =16\ldotp 7\;\textrm{Hz}\;\textrm{with}\;\textrm{attentuation}\;\textrm{of}-40\textrm{dB}\;$$
% 
% $$F_3 =19\;\textrm{Hz}\;\textrm{with}\;\textrm{attenuation}\;\textrm{of}-50\textrm{dB}$$
% 
% 
%% Problem 7.3 
% *Text Problem 11.31 (Page 695)*
% 
% A lowpass digital filter's specifications are given by:
% 
% $$\omega_\mathrm{p}=0.4\mathrm{\pi}, \quad A_\mathrm{p}=0.5\text{ dB},\quad 
% \omega_\mathrm{s}=0.55\mathrm{\pi}, \quad A_\mathrm{s}=50\text{ dB}.$$
%% 
% *(a)* Using bilinear transformation and Chebyshev-I approximation approach 
% obtain a system function $H(z)$ in the rational function form that satisfies 
% the above specifications.
% 
% *Solution:* 

clc; close all; clear;
Omegap = 0.4*pi; Omegas = 0.55*pi;
Ap = 0.5; As = 50; Td = 2; 
Omegap = (2/Td)*tan(Omegap/2); Omegas = (2/Td)*tan(Omegas/2);
[N,Omegac] = cheb1ord(Omegap,Omegas,Ap,As,'s')
[C,D] = cheby1(N,Ap,Omegac,'s');
%% 
% After using the Chebyshev-I approximation approach, we find that the CT system 
% function $H_c \left(s\right)$ is:
% 
% $$H_c \left(s\right)=\frac{0\ldotp 0874}{s^9 +1\ldotp 4358s^8 +4\ldotp 5838s^7 
% +4\ldotp 8208s^6 +6\ldotp 9361s^5 +5\ldotp 0495s^4 +3\ldotp 8733s^3 +1\ldotp 
% 6865s^2 +0\ldotp 5853s+0\ldotp 0874}$$
% 
% Then using the bilinear transformation function along with an arbitrary factor 
% of $T_d =2$, we derive the system function $H\left(z\right)$ to be:

[B,A] = bilinear(C,D,1/Td)
%% 
% $$H\left(z\right)=\frac{0\ldotp 0003+0\ldotp 0028z^{-1} +0\ldotp 0097z^{-2} 
% +0\ldotp 0194z^{-3} +0\ldotp 0243z^{-4} +0\ldotp 0194z^{-5} +0\ldotp 0097z^{-6} 
% +0\ldotp 0028z^{-7} +0\ldotp 0003z^{-8} }{1-3\ldotp 8656z^{-1} +8\ldotp 2625z^{-2} 
% -11\ldotp 6939z^{-3} +11\ldotp 7756z^{-4} -8\ldotp 5442z^{-5} +4\ldotp 3559z^{-6} 
% -1\ldotp 4343z^{-7} +0\ldotp 2381z^{-8} }$$
% 
% 
% 
% *(b)* Provide design plots in the form of log-magnitude, phase delay, group-delay, 
% and impulse responses.
% 
% *MATLAB script:* 

% Frequency Response
omega = linspace(0,pi,1001);
H = freqz(B,A,omega);
Hmag = abs(H);
Hdb = mag2db(Hmag);
% Impulse Response
N = 50; n = 0:N; x = (n==0); h = filter(B,A,x);
% Phase Response
Hpha = angle(H);
Hgd = -diff(unwrap(Hpha))./diff(omega); 
Hgd = [Hgd Hgd(end)];
%[sos G] = tf2sos(B,A);
%Hgd = grpdelay(sos,1001)';
Hgd = medfilt1(Hgd,3);
% Plot Results
figure('Units','inches','Position',[0,0,12,4]);
% Magnitude Response Plot
subplot(2,2,1),plot(omega/pi,Hdb,'LineWidth',1.5),title('Log-Magnitude Response (dB)'), grid on
xticks([0 0.4 0.55 1]), xticklabels({'0','0.4','0.55','1'}), yticks([-50 -0.5]), ylim([-80 10]),
xlabel('\omega/\pi'),ylabel('Magnitude')
% Impulse Response Plot
subplot(2,2,2),stem(n,h,'filled'),title('Impulse Response'),xlabel('Time Index t/s'),ylabel('\it{h(t)}')
% Group-Delay Plot
subplot(2,2,3), plot(omega/pi,Hgd,'LineWidth',1.5), title('Group-Delay Response'), xlabel('\omega/\pi')
ylabel('Samples'), ylim([0 30])
% Phase Delay Plot
subplot(2,2,4), plot(omega/pi,Hpha,'LineWidth',1.5), title('Phase Delay Response'), xlabel('\omega/\pi')
ylabel('Radians')
%% 
% 
% 
% *(c)* Determine the exact band-edge frequencies for the given specifications.
% 
% *Solution:* 

% Exact Band-Edge Frequencies
ind = find(Hdb > -Ap,1,'last'); w1 = omega(ind)/pi % Exact Passband Edge
ind = find(Hdb < -As,1,'first'); w2 = omega(ind)/pi % Exact Stopband Edge
%% 
% Thus, the exact band-edge frequencies are 
% 
% $\omega_p =0\ldotp 399\pi$ and $\omega_s =0\ldotp 548\pi$, which satisfy the 
% design requirements.
% 
% 
%% Problem 7.4
% *Text Problem 11.38*
% 
% A highpass filter specifications are given by:
% 
% $$\omega_\mathrm{s}=0.6\mathrm{\pi}, \quad A_\mathrm{s}=40\text{ dB},\quad 
% \omega_\mathrm{p}=0.8\mathrm{\pi}, \quad A_\mathrm{p}=1\text{ dB}.$$
% 
% 
% 
% *(a)* Using the Butterworth approximation obtain a system function $H(z)$ 
% in the cascade function form that satisfies the above specifications.
% 
% *Solution:* 

clc; close all; clear;
omegas = 0.6*pi; omegap = 0.8*pi;
As = 40; Ap = 1;
[N,Omegac] = buttord(omegap/pi,omegas/pi,Ap,As)
[bhp,ahp] = butter(N,Omegac,'high')
[b0,B,A] = dir2cas(bhp,ahp)
%% 
% Thus, the rational system function $H\left(z\right)$ in its cascade form is:
% 
% $$H\left(z\right)\;=\;\left(2\ldotp 0235\times {10}^{-4} \right)\times \left(\frac{1-1\ldotp 
% 9907+0\ldotp 9908z^{-2} }{1+1\ldotp 3114z^{-1} +0\ldotp 7441z^{-2} }\right)\times 
% \left(\frac{1-2\ldotp 0031z^{-1} +1\ldotp 0031z^{-2} }{1+1\ldotp 0657z^{-1} 
% +0\ldotp 4174z^{-2} }\right)\times \left(\frac{1-2\ldotp 0135z^{-1} +1\ldotp 
% 0135z^{-2} }{1+0\ldotp 9434z^{-1} +0\ldotp 2547z^{-2} }\right)\times \left(\frac{1-0\ldotp 
% 9927z^{-1} }{1+0\ldotp 4532z^{-1} }\right)$$

[sos G] = tf2sos(bhp,ahp);
%% 
% 
% 
% *(b)* Provide design plots in the form of log-magnitude, phase delay, group-delay, 
% and impulse responses.
% 
% *Solution:* 

% Frequency Response
omega = linspace(0,pi,2001);
H = freqz(bhp,ahp,omega);
Hmag = abs(H);
Hdb = mag2db(Hmag);
% Impulse Response
N = 50; n = 0:N; x = (n==0); h = filter(bhp,ahp,x);
% Phase Response
Hpha = angle(H);
Hgd = grpdelay(sos,2001)';
% Plot Results
figure('Units','inches','Position',[0,0,12,4]);
% Log-Magnitude Response Plot
subplot(2,2,1), plot(omega/pi,Hdb,'LineWidth',1.5),title('Log-Magnitude Response (dB)'), grid on
xticks([0 0.6 0.8 1]),xticklabels({'0','0.6','0.8','1'}), yticks([-40 1]), ylim([-80 10]),
xlabel('\omega/\pi'),ylabel('Magnitude')
% Impulse Response Plot
subplot(2,2,2),stem(n,h,'filled'),title('Impulse Response'),xlabel('Time Index t/s'),ylabel('\it{h(t)}')
ylim([-0.25 0.25])
% Group-Delay Plot
subplot(2,2,3), plot(omega/pi,Hgd,'LineWidth',1.5), title('Group-Delay Response'), xlabel('\omega/\pi')
ylabel('Samples'), ylim([0 15])
% Phase Delay Plot
subplot(2,2,4), plot(omega/pi,Hpha,'LineWidth',1.5), title('Phase Delay Response'), xlabel('\omega/\pi')
ylabel('Radians')
%% 
% 
% 
% *(c)* Determine the exact band-edge frequencies for the given specifications.
% 
% *Solution:* 

% Exact Band-Edge Frequencies
ind = find(Hdb > -Ap,1,'first'); w1 = omega(ind)/pi % Exact Passband Edge
ind = find(Hdb < -As,1,'last'); w2 = omega(ind)/pi % Exact Stopband Edge
%% 
% Thus, the exact band-edge frequencies are 
% 
% $\omega_p =0\ldotp 79\pi$ and $\omega_s =0\ldotp 5995\pi$, which meets the 
% design requirements for the highpass filter specifications.
% 
% 
%% Problem 7.5
% *Text Problem 11.43 (Page 698)*
% 
% A digital filter is specified by the following band parameters:
% 
% $$\begin{array}{ll}\text{Band-1: } [0,0.3\mathrm{\pi}], & \text{Attn. } = 
% 50 \text{ dB},\\\text{Band-2: } [0.4\mathrm{\pi},0.5\mathrm{\pi}], & \text{Attn. 
% } = 1 \text{ dB},\\\text{Band-3: } [0.6\mathrm{\pi},\mathrm{\pi}], & \text{Attn. 
% } = 50 \text{ dB}.\end{array}$$
% 
% 
% 
% *(a)* Using Chebyshev II approximation, obtain a system function $H(z)$ in 
% the rational function form that satisfies the above specifications.
% 
% *Solution:* 

clc; close all; clear;
omegas1 = 0.3*pi; omegap1 = 0.4*pi; 
omegas2 = 0.6*pi; omegap2 = 0.5*pi;
As1 = 50; Ap = 1; As2 = 50; As = max(As1,As2);
omegas = [omegas1 omegas2];
omegap = [omegap1 omegap2];
[N,Omegac] = cheb2ord(omegap/pi,omegas/pi,Ap,As)
[b,a] = cheby2(N,As,Omegac,'bandpass') 
%% 
% With the numerator and denominator coefficients now calculated, the system 
% function can be described as
% 
% $$H\left(z\right)=\frac{0\ldotp 0068-0\ldotp 0054z^{-1} +0\ldotp 0053z^{-2} 
% -0\ldotp 0055z^{-3} +0\ldotp 0114z^{-4} -0\ldotp 0055z^{-5} +0\ldotp 0053z^{-6} 
% -0\ldotp 0054z^{-7} +0\ldotp 0068z^{-8} }{1-1\ldotp 2238z^{-1} +3\ldotp 5354z^{-2} 
% -2\ldotp 8788z^{-3} +4\ldotp 2798z^{-4} -2\ldotp 2237z^{-5} +2\ldotp 1128z^{-6} 
% -0\ldotp 5604z^{-7} +0\ldotp 3536z^{-8} }$$
% 
% 
% 
% *(b)* Provide design plots in the form of magnitude, log-magnitude, group-delay, 
% and impulse responses.
% 
% *MATLAB script:* 

% Frequency Response
omega = linspace(0,pi,1001);
H = freqz(b,a,omega);
Hmag = abs(H);
Hdb = mag2db(Hmag);
% Impulse Response
N = 50; n = 0:N; x = (n==0); h = filter(b,a,x);
% Phase Response
Hpha = angle(H);
Hgd = -diff(unwrap(Hpha))./diff(omega); 
Hgd = [Hgd Hgd(end)];
Hgd = medfilt1(Hgd,3);
% Plot Results
figure('Units','inches','Position',[0,0,12,4]);
% Magnitude Response Plot
subplot(2,2,1),plot(omega/pi,Hmag,'LineWidth',1.5),title('Magnitude Response'), grid on
xticks([0 0.3 0.4 0.5 0.6 1]), yticks([0 0.707 1]), ylim([0 1.1]),
xlabel('\omega/\pi'),ylabel('Magnitude')
% Log-Magnitude Response Plot
subplot(2,2,2),plot(omega/pi,Hdb,'LineWidth',1.5),title('Magnitude Response (dB)'), grid on
xticks([0 0.3 0.4 0.5 0.6 1]), yticks([-50 -1]), ylim([-100 10]),
xlabel('\omega/\pi'),ylabel('Magnitude')
% Group-Delay Plot
subplot(2,2,3), plot(omega/pi,Hgd,'LineWidth',1.5), title('Group-Delay Response'), xlabel('\omega/\pi')
ylabel('Samples'), ylim([0 24])
% Impulse Response Plot
subplot(2,2,4),stem(n,h,'filled'),title('Impulse Response'),xlabel('Time Index t/s'),ylabel('\it{h(t)}')
ylim([-0.15 0.15])
%% 
% 
% 
% *(c)* Determine the exact band-edge frequencies for the given attenuation.
% 
% *Solution:* 

% Exact Band-Edge Frequencies
ind = find(Hdb > -Ap); wp1 = omega(ind)/pi; % Exact Passband Edges
LowerPassEdge = wp1(1), UpperPassEdge = wp1(end)
inds1 = find(Hdb(1:401) < -As,1,'last'); % Exact Lower Stopband Edge
LowerStopEdge = omega(inds1)/pi
inds2 = find(Hdb < -As); % Exact Upper Stopband Edge
LowerStopEdge = omega(602)/pi
%% 
% Thus the exact band-edge frequencies are:
% 
% $\omega_{\textrm{s1}} =0\ldotp 3\pi$, $\omega_{\textrm{s2}} =0\ldotp 601\pi$
% 
% $\omega_{\textrm{p1}} =0\ldotp 392\pi$, $\omega_{\textrm{p2}} =0\ldotp 497\pi$
% 
% 
%% Problem 7.6
% *Text Problem 11.66 (Page 701)*
% 
% A digital filter is specified by the following band parameters:
% 
% $$\begin{array}{ll}\text{Band-1: } [0,0.2\mathrm{\pi}], & \text{Attn. } = 
% 1 \text{ dB},\\\text{Band-2: } [0.35\mathrm{\pi},0.5\mathrm{\pi}], & \text{Attn. 
% } = 50 \text{ dB},\\\text{Band-3: } [0.65\mathrm{\pi},\mathrm{\pi}], & \text{Attn. 
% } = 1 \text{ dB}.\end{array}$$
% 
% *(a)* Using Butterworth approximation, obtain a system function $H(z)$ in 
% the cascade form that satisfies the above specifications.
% 
% *Solution:* 

clc; close all; clear;
wp1 = 0.2; ws1 = 0.35; Ap1 = 1;
wp2 = 0.65; ws2 = 0.5; As = 50;
omegap = [wp1 wp2]; 
omegas = [ws1 ws2];
[N,Omegac] = buttord(omegap,omegas,Ap1,As)
[b,a] = butter(N,Omegac,'stop');
[bo,B,A] = dir2cas(b,a)
%% 
% Here the coefficents for the cascade form of $H\left(z\right)$ can be written 
% as:
% 
% $$H\left(z\right)=\left(0\ldotp 0982\right)\times \left(\frac{1-0\ldotp 4752z^{-1} 
% +1\ldotp 0001z^{-2} }{1+0\ldotp 5928z^{-1} +0\ldotp 7650z^{-2} }\right)\times 
% \left(\frac{1-0\ldotp 4766z^{-1} +0\ldotp 9956z^{-2} }{1-0\ldotp 3322z^{-1} 
% +0\ldotp 4315z^{-2} }\right)\times \left(\frac{1-0\ldotp 4788z^{-1} +1\ldotp 
% 0045z^{-2} }{1-0\ldotp 0441z^{-1} +0\ldotp 2424z^{-2} }\right)\times \left(\frac{1-0\ldotp 
% 4815z^{-1} +0\ldotp 9956z^{-2} }{1-0\ldotp 5481z^{-1} +0\ldotp 2885z^{-2} }\right)\times 
% \left(\frac{1-0\ldotp 4838z^{-1} +1\ldotp 0044z^{-2} }{1-0\ldotp 9404z^{-1} 
% +0\ldotp 5150z^{-2} }\right)\times \left(\frac{1-0\ldotp 4851z^{-1} +0\ldotp 
% 999z^{-2} }{1-1\ldotp 2428z^{-1} +0\ldotp 8133z^{-2} }\right)$$
% 
% 
% 
% *(b)* Provide design plots in the form of magnitude, log-magnitude, group-delay, 
% and impulse responses.
% 
% *Solution:* 

% Frequency Response
omega = linspace(0,pi,1001);
H = freqz(b,a,omega);
Hmag = abs(H);
Hdb = mag2db(Hmag);
% Impulse Response
L = 50; n = 0:L; x = (n==0); h = filter(b,a,x);
% Phase Response
Hpha = angle(H);
Hgd = -diff(unwrap(Hpha))./diff(omega); 
Hgd = [Hgd Hgd(end)];
Hgd = medfilt1(Hgd,3);
% Plot Results
figure('Units','inches','Position',[0,0,12,4]);
% Magnitude Response Plot
subplot(2,2,1),plot(omega/pi,Hmag,'LineWidth',1.5),title('Magnitude Response'), grid on
xticks([0 0.2 0.35 0.5 0.65 1]), yticks([0 0.707 1]), ylim([0 1.1]),
xlabel('\omega/\pi'),ylabel('Magnitude')
% Log-Magnitude Response Plot
subplot(2,2,2),plot(omega/pi,Hdb,'LineWidth',1.5),title('Magnitude Response (dB)'), grid on
xticks([0 0.2 0.35 0.5 0.65 1]), yticks([-50 -1]), ylim([-100 10]),
xlabel('\omega/\pi'),ylabel('Magnitude')
% Group-Delay Plot
subplot(2,2,3), plot(omega/pi,Hgd,'LineWidth',1.5), title('Group-Delay Response'), xlabel('\omega/\pi')
ylabel('Samples'),xticks([0 0.2 0.35 0.5 0.65 1]), ylim([0 24])
% Impulse Response Plot
subplot(2,2,4),stem(n,h,'filled'),title('Impulse Response'),xlabel('Time Index t/s'),ylabel('\it{h(t)}')
ylim([-0.3 0.6])
%% 
% 
% 
% *(c)* Determine the exact band-edge frequencies for the given attenuation.
% 
% *Solution:* 

ind = find(Hdb(1:401) < -Ap1, 1,'first'); wplow = omega(ind)/pi;
ind = find(Hdb(1:401) > -As, 1,'last'); wslow = omega(ind)/pi;
ind = Hdb(500:502); wsupper = omega(502)/pi;
ind = find(Hdb(1:701) < -Ap1, 1,'last'); wpupper = omega(ind)/pi;
wplow,wslow,wsupper,wpupper
%% 
% Thus, the exact band-edge frequencies for this Butterworth BandStop IIR Filter 
% are:
% 
% $\omega_{\textrm{p1}} =0\ldotp 243\pi$, $\omega_{\textrm{s1}} =0\ldotp 349\pi$, 
% $\omega_{\textrm{s2}} =0\ldotp 501\pi$, and $\omega_{\textrm{p2}} =0\ldotp 631\pi$ 
% 
% 
%% Problem 7.7
% *Text Problem 11.70 (Page 702)*
% 
% An analog signal $x_\mathrm{c}(t)=5\sin(2\mathrm{\pi}250t)+10\sin(2\mathrm{\pi}300t)$is 
% to be processed using the effective continuous-time system of Figure 6.18 in 
% which the sampling frequency is 1 kHz.
% 
% 
% 
% *(a)* Design a minimum-order IIR digital filter that will suppress the $300$ 
% Hz component down to $50$ dB while pass the $250$ Hz component with attenuation 
% of less than $1$ dB. The digital filter should have an equiripple passband and 
% stopband. Determine the system function of the filter and plot its log-magnitude 
% response in dB.
% 
% *Solution:* 
% 
% Since the requirements call for a filter with an equiripple passband and stopband, 
% the only filter that can achieve the design specs are an analog Elliptic filter

clc; close all; clear;
% Design Requirements
Fs = 1000; Fs2 = Fs/2;
fs = 300/Fs; As = 50; 
fp = 250/Fs; Ap = 1;
omegas = 2*fs;
omegap = 2*fp;
[N,omegac] = ellipord(omegap,omegas,Ap,As)
[b1,a1] = ellip(N,Ap,As,omegac)
%% 
% Deriving the coefficients of the Digital Elliptic Lowpass IIR Filter shows 
% that the system function is:
% 
% $$H\left(z\right)=\frac{0\ldotp 0481+0\ldotp 1381z^{-1} +0\ldotp 2542z^{-2} 
% +0\ldotp 301z^{-3} +0\ldotp 2542z^{-4} +0\ldotp 1381z^{-5} +0\ldotp 0481z^{-6} 
% }{1-1\ldotp 166z^{-1} +2\ldotp 2689z^{-2} -1\ldotp 8296z^{-3} +1\ldotp 5033z^{-4} 
% -0\ldotp 6909z^{-5} +0\ldotp 2401z^{-6} }$$

% Frequency Response
omega = linspace(0,pi,2001);
H = freqz(b1,a1,omega);
Hmag = abs(H);
Hdb = mag2db(Hmag);
% Log-Magnitude Response Plot
figure('Units','inches','Position',[0,0,6,3]);
plot(omega/pi,Hdb,'LineWidth',1.5),title('Log-Magnitude Response (dB)'), grid on
xticks([0 0.5 0.6 1]), yticks([-50 -1]), ylim([-100 10]),
xlabel('\omega/\pi'),ylabel('Magnitude')
%% 
% 
% 
% *(b)* Process the signal $x_\mathrm{c}(t)$ through the effective analog system. 
% Generate sufficient samples so that the output response $y_\mathrm{c}(t)$ goes 
% into steady-state. Plot the steady-state $y_{\mathrm{ss}}(t)$ and comment on 
% the filtering result.
% 
% *Solution:* 

Ts = 1/Fs; t = 0:0.0001:.5;
% Continuous-Time Signal
xc = 5*sin(2*pi*250.*t)+10*sin(2*pi*300.*t);
% Discrete-Time Signal
n = 0:500; nTs = n*Ts;
xn = 5*sin(2*pi*250.*nTs)+10*sin(2*pi*300.*nTs);
% Plot the CT signal and Sampled DT signal
figure('Units','inches','Position',[0,0,12,4]);
plot(t(1:1000)*1000,xc(1:1000),'b','LineWidth',1.5), hold on, stem(nTs(1:101)*1000,xn(1:101),'--r','LineWidth',1.5), xlabel('Time (ms)'),ylabel('Amplitude')
title('CT Signal and Sampled DT Signal'), legend('x_c(t)','x[n]')
% Filter the Discrete-Signal x[n] through Elliptic Lowpass IIR Filter
yn = filter(b1,a1,xn);
% Reconstruct the analog signal through interpolation and analyze results 
yt = yn * sinc(Fs*(ones(length(n),1)*t-nTs'*ones(1,length(t))));
figure('Units','inches','Position',[0,0,12,3]);
plot(t(4001:5000)*1000,yt(4001:5000),'b','LineWidth',1.5), hold on, stem(n(401:500),yn(401:500),'--r','LineWidth',1.5)
xlabel('Time (msec)')
ylabel('Amplitude')
title('IIR Filtered DT and Reconstructed Steady-State CT signal'), legend('y_c(t)','y[n]')
%% 
% *Comment:* 
% 
% After analyzing the the reconstructed signal, $y_c \left(t\right)$, we see 
% that the resulting signal is a pure sinusoid with only one fundamental frequency. 
% The higher frequency of 300 Hz has been filtered out by the IIR filter and only 
% the 250 Hz component remains. Once the signal reaches a steady-state, there 
% is a small difference in magnitude as the reconstructed signal has a peak amplitude 
% of about 4.5, compared to the original sinusoidal component $\;5\sin \left(2\pi 
% 250\right)$.
% 
% 
% 
% *(c)* Repear parts (a) and (b) by designing an equiripple FIR filter. Compare 
% the orders of the two filters and their filtering results.
% 
% *Solution:* 

% Compute Absolute passband and stopband ripple values
[dp,ds] = spec_convert(Ap,As,'rel','abs')
% Estimate Filter using FIRPMORD function
f = [0.5 0.6]; % Band-edge array
a = [1 0]; % Band-edge desired gain
dev = [dp,ds]; % Band tolerance
[M,fo,ao,W] = firpmord(f,a,dev); M
% Filter Design using FIRPM function
[h,delta] = firpm(M,fo,ao,W);
% Frequency Response
omega = linspace(0,pi,2001);
H = freqz(h,1,omega);
Hmag = abs(H);
Hdb = mag2db(Hmag);
% Log-Magnitude Response Plot
figure('Units','inches','Position',[0,0,6,3]);
plot(omega/pi,Hdb,'LineWidth',1.5),title('Log-Magnitude Response (dB)'), grid on
xticks([0 0.5 0.6 1]), yticks([-50 0]), ylim([-100 10]),
xlabel('\omega/\pi'),ylabel('Magnitude')
% Filter the Discrete-Signal x[n] through Equiripple Lowpass FIR Filter
yn = filter(h,1,xn);
% Reconstruct the analog signal through interpolation and analyze results 
yt = yn * sinc(Fs*(ones(length(n),1)*t-nTs'*ones(1,length(t))));
figure('Units','inches','Position',[0,0,12,3]);
plot(t(4001:5000)*1000,yt(4001:5000),'b','LineWidth',1.5), hold on 
stem(n(401:500),yn(401:500),'--r','LineWidth',1.5)
xlabel('Time (msec)')
ylabel('Amplitude')
title('FIR Filtered DT and Reconstructed Steady-State CT signal'), legend('y_c(t)','y[n]')
%% 
% *Comment:* 
% 
% After drawing comparisons between the two filters, the Elliptic IIR filter 
% had an order of $N=6$, while the FIR filter had order $N=34$, which is almost 
% 6 times the order of the IIR Filter. When comparing the number of multiplications 
% per output sample:
% 
% The IIR Filter had order $N=6$, thus this filter had $4\left(\frac{N}{2}\right)\to 
% 4\left(\frac{6}{2}\right)\to 4\left(3\right)$= 12 multiplications per output 
% sample
% 
% The FIR Filter had order $M=34$, thus the number of multiplications per output 
% sample was $\frac{M}{2}\to \frac{34}{2}=17$. 
% 
% When comparing the exact band-edge frequencies, the IIR filter had a passband 
% edge at $\omega_p =0\ldotp 4995\pi$ and a sharper transition band with a stopband 
% edge at $\omega_s =0\ldotp 5575\pi$
% 
% Whereas the FIR filter had a narrower stopband width with a passband edge 
% at $\omega_p =0\ldotp 505\pi$ and stopband edge at $\omega_s =0\ldotp 5995\pi$ 
% 
% The resulting reconstructed signals had a difference in amplitude, as the 
% FIR filter produced a signal that was closer to the sinusoidal component of 
% $5\sin \left(2\pi 250\right)$, whereas the IIR had a smaller peak amplitude 
% due to its sharper transition band cutting off frequencies before the desired 
% stopband-edge frequency.
% 
% Overall, the IIR filter achieved the design requirements using about 6 times 
% less coefficients and 5 multiplcations/output sample less than the FIR filter 
% while also providing a shorter transition band that provided the necessary attenuations 
% in the corresponding frequency bands. 
% 
% 
%% Problem 7.8
% *Text Problem 11.71 (Page 703)* 
% 
% Consider the following bandpass digital filter specifications:
% 
% $$\begin{array}{rcl}\text{Stopband-1}&:& [0,0.4\pi], \quad&\text{Attn.}&=& 
% 40\ \text{dB} \\\text{Passband}&:& [0.45\pi,0.55\pi], \quad&\text{Attn.}&=& 
% 0.5\ \text{dB} \\\text{Stopband-2}&:& [0.65\pi,\pi], \quad&\text{Attn.}&=& 50\ 
% \text{dB} \end{array}$$
% 
% 
% 
% *(a)* Design a minimum order FIR filter to satisfy the  above specifications. 
% Plot its magnitude, log-magnitude (dB), and group-delay responses in one figure 
% using 3 rows and 1 column.
% 
% *MATLAB script*: 

clc; close all; clear;
omegas1 = 0.4*pi; omegap1 = 0.45*pi;
omegas2 = 0.65*pi; omegap2 = 0.55*pi; 
As1 = 40; As2 = 50; As = max(As1,As2);
Ap = 0.5; 
[dp,ds1] = spec_convert(Ap,As1,'rel','abs')
[dp,ds2] = spec_convert(Ap,As2,'rel','abs')
% Estimate Filter using FIRPMORD function
f = [omegas1,omegap1,omegap2,omegas2]/pi; % Band-edge array
a = [0 1 0]; % Band-edge desired gain
dev = [ds1 dp ds2]; % Band tolerance
[M,fo,ao,W] = firpmord(f,a,dev); M
% Filter Design using FIRPM function
[h,delta] = firpm(M,fo,ao,W);
n = 0:M; omega = linspace(0,pi,1001);
% Magnitude Response
H = freqz(h,1,omega); Hmag = abs(H);
Hdb = mag2db(Hmag);
% Phase and Group-Delay
Hpha = angle(H);
Hgd = -diff(unwrap(Hpha))./diff(omega); 
Hgd = [Hgd Hgd(end)];
Hgd = medfilt1(Hgd,3);
% Plot Results
figure
% Magnitude Response Plot
subplot(3,1,1),plot(omega/pi,Hmag,'LineWidth',1.5),title('Magnitude Response'), grid on
xticks([0 0.4 0.45 0.55 0.65 1]), yticks([0 0.707 1]), ylim([0 1.1]),
xlabel('\omega/\pi'),ylabel('Magnitude')
% Log-Magnitude Response Plot
subplot(3,1,2),plot(omega/pi,Hdb,'LineWidth',1.5),title('Magnitude Response (dB)'), grid on
xticks([0 0.4 0.45 0.55 0.65 1]), yticks([-50 -40 -0.5]), ylim([-80 10]),
xlabel('\omega/\pi'),ylabel('Magnitude')
% Group-Delay Plot
subplot(3,1,3), plot(omega/pi,Hgd,'LineWidth',1.5), title('Group-Delay Response'), xlabel('\omega/\pi')
ylabel('Samples'),xticks([0 0.4 0.45 0.55 0.65 1]), ylim([32 33])
%% 
% 
% 
% *(b)* Design a minimum order IIR filter to satisfy the above specifications. 
% Plot its magnitude, log-magnitude (dB), and group-delay responses in one figure 
% using 3 rows and 1 column. From your plots determine the exact band-edge frequencies.
% 
% *MATLAB script*: 

Omegas = [omegas1 omegas2];
Omegap = [omegap1 omegap2];
[N,Omegac] = ellipord(Omegap/pi,Omegas/pi,Ap,As); N
[b,a] = ellip(N,Ap,As,Omegac,'bandpass');
omega = linspace(0,pi,1001);
% Magnitude Response
H = freqz(b,a,omega); Hmag = abs(H);
Hdb = mag2db(Hmag);
% Phase and Group-Delay
Hpha = angle(H);
Hgd = -diff(unwrap(Hpha))./diff(omega); 
Hgd = [Hgd Hgd(end)];
Hgd = medfilt1(Hgd,3);
% Plot Results
figure
% Magnitude Response Plot
subplot(3,1,1),plot(omega/pi,Hmag,'LineWidth',1.5),title('Magnitude Response'), grid on
xticks([0 0.4 0.45 0.55 0.65 1]), yticks([0 0.707 1]), ylim([0 1.1]),
xlabel('\omega/\pi'),ylabel('Magnitude')
% Log-Magnitude Response Plot
subplot(3,1,2),plot(omega/pi,Hdb,'LineWidth',1.5),title('Magnitude Response (dB)'), grid on
xticks([0 0.4 0.45 0.55 0.65 1]), yticks([-50 -0.5]), ylim([-80 10]),
xlabel('\omega/\pi'),ylabel('Magnitude')
% Group-Delay Plot
subplot(3,1,3), plot(omega/pi,Hgd,'LineWidth',1.5), title('Group-Delay Response'), xlabel('\omega/\pi')
ylabel('Samples'),xticks([0 0.4 0.45 0.55 0.65 1]),ylim([0 100])
%% 
% 
% 
% *(c)* Compare the two filter designs in terms of their responses and orders.
% 
% *Comparison*: 
% 
% The FIR filter design was designed using the Parks-McClellan algorithm which 
% produced an equiripple FIR filter with 
% 
% minimum order $M=65$. The design produced a linear-phase which resulted in 
% a constant group-delay. The magnitude response was unsatisfactory towards the 
% design specifications, as the first stopband did not meet the attenuation requirement 
% of 40dB, due to differencebetween the lower and upper transition bands. The 
% number of multiplications per output sample for the FIR filter was $\frac{M}{2}\to 
% \frac{65}{2}=32\ldotp 5$
% 
% The IIR filter was designed using a digital Elliptic filter and achieved a 
% minimum order of $N=5\ldotp$ This filter was able to achieve attenuation requirements 
% in both stopbands and the passband, but resulted in a nonlinear phase and group-delay 
% response. The number of multiplications per output sample for the IIR filter 
% was $4\left(\frac{N}{2}\right)\to 4\left(\frac{5}{2}\right)\to 4\left(2\ldotp 
% 5\right)=10$
% 
% Thus, the IIR filter was able to satisfy the design requirements with order 
% $N=5$, which resulted in only 10 multiplications per output sample. This is 
% about 3 times less the number of multiplications, so the IIR filter not only 
% satisfied the magnitude requirements, but is it also computationally quicker. 
% 
% 

function [A,B] = spec_convert(C,D,typein,typeout)
%  typein: 'abs' or 'rel' or 'ana'
% typeout: 'abs' or 'rel' or 'ana'
%     C,D: input specifications
%     A,B: output specifications

% Enter your function code below
if nargin > 4
    error('too many input arguments')
end

switch typein
    case 'abs'
        switch typeout
            case 'rel'
                Ap = 20*log10((1+C)/(1-C)); % Relative Output: Passband ripple
                As = floor(20*log10((1+C)/D)); % Relative Output: Stopband Attenuation
                A = Ap; B = As;
            case 'ana'
                Ap = 20*log10((1+C)/(1-C)); % Relative Output: Passband ripple
                As = floor(20*log10((1+C)/D)); % Relative Output: Stopband Attenuation
                A = sqrt(10^(Ap/10)-1); % Analog Output: Passband
                B = floor(10^(As/20)); % Analog Output: Stopband                
        end
    case 'rel'
        switch typeout
            case 'abs'
                dp = (10^(C/20)-1)/(1+10^(C/20)); % Absolute Output: Passband Error
                ds = (1+dp)/(10^(D/20)); % Absolute Output: Stopband Error
                A = dp; B = ds;
            case 'ana'
                epsilon = sqrt(10^(C/10)-1); % Analog Output: Passband
                B = floor(10^(D/20)); % Analog Output: Stopband
                A = epsilon;
        end
   case 'ana'
        switch typeout
            case 'rel'
                Ap = round(10*log10(C^(2)+1),2); % Relative Output: Passband ripple (in dB)
                As = 20*log10(D); % Relative Output: Stopband Attenuation (in dB)
                A = Ap; B = As;
            case 'abs'
                Ap = round(10*log10(C^(2)+1),2); % Relative Output: Passband ripple (in dB)
                As = 20*log10(D); % Relative Output: Stopband Attenuation (in dB)
                dp = (10^(Ap/20)-1)/(1+10^(Ap/20)); % Absolute Output: Passband Error
                ds = (1+dp)/(10^(Ap/20)); % Absolute Output: Stopband Error
                A = dp; B = ds;
        end     
end
end