%% *EECE5666 (DSP) : Homework-3*
% *Due on February 22, 2022 by 11:59 pm via submission portal.* 
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
%% Problem 3.1
% *Text Problem 5.23 (Page 282)* 
% An LTI system is described by the difference equation
% 
% $$y[n]=bx[n]+0.8y[n-1]-0.81y[n-2] \qquad\qquad (5.23.1)$$

clc; close all; clear;
%% 
% 
% 
% *(a)* Determine the frequency response of the system in terms of $b.$
% 
% *Solution:* 
% 
% First we solve the difference equation for the z-transform of the impulse 
% response and then substitute $z = e^{j\omega}$ 
% 
% $$y[n] -0.8y[n-1] + 0.81y[n-2] = bx[n]$$
% 
% $$Y(z)(1 - 0.8z^{-1}+0.81z^{-2}) = bX(z)$$
% 
% $\frac{Y(z)}{X(z)} = H(z) = \frac{b}{1 - 0.8z^{-1}+0.81z^{-2}$ then substituting 
% $z = e^{j\omega}$ we get the frequency response
% 
% $$H(e^{j\omega}) = \frac{b}{1 - 0.8e^{-j\omega}+0.81e^{-2j\omega}$$
% 
% We can expand this further using the euler's identities:
% 
% $e^{-j\omega} = \cos\omega - j\sin\omega$ and $e^{-2j\omega} = \cos2\omega 
% - j\sin2\omega$ 
% 
% $$H(e^{j\omega}) = \frac{b}{1 - 0.8(\cos\omega -j\sin\omega) +0.81(\cos2\omega-j\sin2\omega)$$
% 
% $$H(e^{j\omega}) = \frac{b}{1 - 0.8\cos\omega +0.81\cos2\omega + j(0.8\sin\omega 
% -0.81\sin2\omega)$$
% 
% Then the Magnitude and Phase Responses would be
% 
% $$|H(e^{j\omega})| = \frac{|b|}{\sqrt{(1 - 0.8\cos\omega +0.81\cos2\omega)^2 
% + (0.8\sin\omega -0.81\sin2\omega)^2$$
% 
% $$\angle{H(e^{j\omega})} = \angle{b} - \tan^{-1}\left(\frac{0.8\sin\omega-0.81\sin2\omega}{1-0.8\cos\omega+0.81\cos2\omega}\right)$$
% 
% 
% 
% *(b)* Determine $b$ so that $\bigl|H(\mathrm{e}^{\j\omega})\bigr|_{\mathrm{max}}=1.$ 
% Plot the resulting magnitude response.
% 
% *Solution*: 
% 
% If we want the maximum magnitude to be equal to one, then we choose $b$ so 
% that $|H(e^{j\omega})|$ = 1 $\rightarrow b = \frac{1}{|H(e^{j\omega})|} = 0.17$
% 
% *MATLAB script for the magnitude plot*: 

omega = linspace(-pi,pi,1000);
b = 0.17;
a = [1 -0.8 0.81];
H = freqz(b,a,omega);
figure
plot(omega/pi, abs(H))
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
xlabel('\omega/\pi')
ylabel('|H(e^{j\omega})|')
title('|H(e^{j\omega})| = 1 when b = 0.17')
%% 
% 
% 
% *(c)* Graph the wrapped and the unwrapped phase responses in one plot.
% 
% *MATLAB script*: 

H_unwrapped = phasez(b,a,omega)';
figure
plot(omega/pi,angle(H)/pi,':b')
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
hold on
plot(omega/pi,H_unwrapped/pi,'--r')
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
legend('\angle H(e^{j\omega})', '\Phi (\omega)')
hold off
title('Wrapped and Unwrapped Phase Responses')
xlabel('\omega/\pi')
ylabel('\angle{H(e^{j\omega})} (radians)')
%% 
% 
% 
% *(d)* Determine analytically response $y[n]$ to the input $x[n] = 2\cos(\mathrm{\pi} 
% n/3+ 45^\deg)$. Plot the output sequence $y[n].$
% 
% *Solution:* 
% 
% Analyzing the input sequence: $x[n]$ we see that $\omega = \pi/3$, so to determine 
% how the system will alter this input sequence, we solve the system's frequency 
% response with $\omega = \pi/3$
% 
% Again, 
% 
% $$H(e^{j\omega}) = \frac{0.17}{1 - 0.8e^{-j\omega}+0.81e^{-2j\omega}$$
% 
% $$H(e^{j\pi/3}) = \frac{0.17}{1 - 0.8e^{-j\pi/3}+0.81e^{-j2\pi/3}} = 0.1229 
% - j0.1217 $$
% 
% $H(e^{j\pi/3}) = 0.1729\angle{-44.715^\circ} \rightarrow$ thus the input sequence 
% will be scaled by a factor of 0.1729 and have a phase shift of -44.715 degrees, 
% resulting in:
% 
% $$y[n] = (0.1729)(2)\cos\left(\frac{\pi}{3}n + 45^\circ - 44.715^\circ \right)$$
% 
% $$y[n] = 0.2458\cos\left(\frac{\pi}{3}n + 0.285^\circ \right)$$
% 
% *MATLAB script for the output sequence plot*: 

N = 6; w = (2*pi/N); b = 0.17;
h = b/(1 - 0.8*exp(-1i*w)+0.81*exp(-2*w));
H_mag = abs(h); H_ang = angle(h);
n = 0:10*N-1; x = 2.*cos((w).*n + (pi/4)); 
y = (H_mag).*2.*cos((pi/3).*n + (pi/4) + H_ang);
omega = linspace(0,10*N-1,1000);
X = 2.*cos((w).*omega + (pi/4)); Y = (H_mag).*2.*cos((w).*omega + (pi/4) + H_ang); 
figure('Units','inches','Position',[0,0,12,4]); hold on
stem(n,x,'b'), stem(n,y,'r') 
plot(omega,X,'b'), plot(omega,Y,'r')
legend('x[n]','y[n]'), hold off
title('Output Sequence: \it{y[n]}')
xlabel('\it{n}'), ylabel('Amplitude')
%% 
% 
% 
% *(e)* Using MATLAB, compute the steady-state response to $x[n]$ above and 
% verify your result.
% 
% *Matlab script:* 

% N = 6; k = 0:N-1; w = (2*pi/N); x = 2.*cos((w).*k + (pi/4));
% ck_x = dtfs(x); H_mag = abs(freqz(b,a,k)); ck_y = H_mag.*ck_x; y = idtfs(ck_y); 
y = filter(b,a,x);
figure
stem(n,y)
xlabel('\it{n}')
ylabel('Amplitude')
title('Steady-State Response from 30 \leq n \leq 60')
%% 
% Verifying the results, we see that the steady-state response begins after 
% the transient response has subsided, which occurs around $n = 30$. The amplitude 
% of the input has been scaled and phase shifted as expected. 
% 
% 
%% Problem 3.2
% *Text Problem 5.30, parts (a) and (e), (Page 283)* 
% Consider a periodic signal 
% 
% $$x[n] = \sin(0.1\mathrm{\pi} n) + \frac13\sin(0.3\mathrm{\pi} n) + \frac15\sin(0.5\mathrm{\pi} 
% n).$$
% 
% For each of the following systems, first plot magnitude and phase responses 
% and then determine if the system imparts (i) no distortion, (ii) magnitude distortion, 
% and/or (iii) phase (or delay) distortion. In each case, also graph the input 
% and the steady state response for $0\leq n\leq 60.$

clc; close all; clear;
n = 0:60;
x = sin(0.1.*pi.*n) + (1/3).*sin(0.3.*pi.*n) + (1/5).*sin(0.5.*pi.*n);
figure
stem(n,x)
title("x[n] = sin(0.1\pin) + 1/3sin(0.3\pin) + 1/5sin(0.5\pin)")
xlabel('\it{n}'), ylabel('\it{x[n]}')
%% 
% 
% 
% *(a)* $h[n] = \{ \underset{\uparrow}{1},-2,3,-4,0,4,-3,2,-1 \}$: 
% 
% *Solution*: 
% 
% *MATLAB script for* *computation and plotting of magnitude and phase responses*:

h = [1 -2 3 -4 0 4 -3 2 -1]; hn = 0:8;
figure
stem(hn,h)
title('System Impulse Response: h[n]')
xlabel('\it{n}')
ylabel('\it{h[n]}')
omega = linspace(-pi,pi,1000);
H = freqz(h,1,omega); H_mag = abs(H); H_ang = angle(H);
figure('Units','inches','Position',[0,0,12,4]);
subplot(1,2,1)
plot(omega/pi,20*log10(abs(H)));
title('Magnitude Response')
ylim([-140 40])
ylabel('|H(e^{j\omega})|(dB)')
xlabel('\omega/\pi')
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
subplot(1,2,2)
plot(omega/pi,angle(H)/pi)
xlabel('\omega/\pi')
ylabel('\angle{H(e^{j\omega})}')
title('Phase Response')
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
%% 
% *Answer for distortion*: 
% 
% A distortionless system can be described as 
% 
% $y[n] = Gx[n-n_d], \quad G > 0$, where $G \text{ and } n_d$ are constants.
% 
% If the system magnitude response is not a constant, the system will provide 
% magnitude distortion. 
% 
% The impulse response of this system can be described as
% 
% $$h[n] = \delta[n] -2\delta[n-1]+3\delta[n-2]-4\delta[n-3]+4\delta[n-5]-3\delta[n-6]+2\delta[n-7] 
% - \delta[n-8]$$
% 
% The system response can be found by taking the DTFT of the impulse response, 
% which results in
% 
% $$H(e^{j\omega}) = 1 - 2e^{-j\omega} +3e^{-j2\omega} - 4e^{-j3\omega} + 4e^{-j5\omega} 
% - 3e^{-j6\omega} + 2e^{-j7\omega} - e^{-j8\omega}$$
% 
% Here we can see that when the frequency $\omega$ changes, the values of magnitude 
% will also change
% 
% Thus, this system is not distortionless and experiences magnitude distortion. 
% 
% This system also experiences a linear phase response before and after $4\omega$ 
% based on the system function, where there is a zero component in 
% 
% the impulse response. 
% 
% Thus, the phase response does not experience a truly linear-phase shift meaning 
% this system also experiences phase distortion. 
% 
% *MATLAB script for the computation of input and output sequences*: 

%y = filter(h,1,x);
y = conv(h,x);
figure('Units','inches','Position',[0,0,12,4]);
subplot(1,2,1)
stem(n,x);
title('Input Sequence: \it{x[n]}')
ylabel('\it{x[n]}')
xlabel('\it{n}')
subplot(1,2,2)
stem(n,y(1:61))
xlabel('\it{n}')
ylabel('\it{y[n]}')
title('Output Sequence: \it{y[n]}')
n = 0:60; u = (n>= 0);
y_ss = filter(h,1,u);
figure
stem(n,y_ss)
xlabel('\it{n}'), ylabel('y_{ss}[n]'), title('Steady-State Response of System')
%% 
% 
% 
% *(e)* $H(z) = \frac{1+1.778z^{-2}+3.1605z^{-4}}{1+ 0.5625z^{-2}+0.3164z^{-4}}$: 
% 
% *Solution*: 
% 
% *MATLAB script for* *computation and plotting of magnitude and phase responses*:

b = [1 0 1.778 0 3.1605]; a = [1 0 0.5625 0 0.3164];
omega = linspace(-pi,pi,1000);
H = freqz(b,a,omega);
figure('Units','inches','Position',[0,0,12,4]);
subplot(1,2,1)
plot(omega/pi,abs(H));
title('Magnitude Response')
ylabel('|H(e^{j\omega})|(dB)')
xlabel('\omega/\pi')
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
subplot(1,2,2)
plot(omega/pi,angle(H)/pi)
xlabel('\omega/\pi')
ylabel('\angle{H(e^{j\omega})}')
title('Phase Response')
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
%% 
% *Answer for distortion*: 
% 
% The response for this system resembles a bandreject filter with a gain factor 
% around 3.16. Since the response is not a constant for all values of frequency 
% $\omega$ $\rightarrow$ system is not distortionless and experiences magnitude 
% distortion for values around $\pm \frac{\pi}{2}$.
% 
% The system also experiences a phase response that does pass through the origin 
% of $\omega = 0$, but it does so as a nonlinear function of frequency
% 
% So, the system has not only a magnitude distortion but a phase distortion 
% as well. 
% 
% *MATLAB script for the computation of input and output sequences*: 

y = filter(b,a,x);
figure('Units','inches','Position',[0,0,12,4]);
subplot(1,2,1)
stem(n,x);
title('Input Sequence: \it{x[n]}')
ylabel('\it{x[n]}')
xlabel('\it{n}')
subplot(1,2,2)
stem(n,y)
xlabel('\it{n}')
ylabel('\it{y[n]}')
title('Output Sequence: \it{y[n]}')
%% 
% 
%% Problem 3.3 
% *Text Problem 5.37 (Page 284)* 
% Compute and plot the phase response using the functions |*freqz*|, |*angle*|, 
% |*phasez*|, |*unwrap*|, and |*phasedelay*| for the following systems.

clc; close all; clear;
%% 
% 
% 
% *(a)* The pure delay $y[n]=x[n-15].$
% 
% *MATLAB script for computation and plotting:* 

x = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]; omega = linspace(0,pi,1000);
H = freqz(x,1,omega);
figure('Units','inches','Position',[0,0,18,6]);
subplot(2,2,1)
plot(omega/pi,abs(H));
title('Magnitude Response')
ylabel('|H(e^{j\omega})|')
xlabel('\omega/\pi')
yticks(1)
subplot(2,2,2)
plot(omega/pi,angle(H)/pi)
xlabel('\omega/\pi')
ylabel('\angle{H(e^{j\omega})(radians)}')
title('Phase Response')
subplot(2,2,3)
plot(omega/pi,unwrap(angle(H))/pi)
xlabel('\omega/\pi')
ylabel('\angle{H(e^{j\omega})}(radians)')
title('Unwrapped Phase Response: \Psi(\omega)')
subplot(2,2,4)
H_pd = phasedelay(x,1,omega)/pi;
plot(omega/pi,round(H_pd,2))
xlabel('\omega/\pi')
ylabel('\angle{H(e^{j\omega})}')
title('Phase Delay')
ylim([4.7 4.85])
%% 
% 
% 
% *(b)* The system defined by $H(z) = \frac{1+z^{-1}+z^{-2}+z^{-3}}{1+0.9z^{-1}+0.81z^{-2}+0.927z^{-3}}.$
% 
% *MATLAB script for computation and plotting:* 

b = [1 1 1 1]; a = [1 0.9 0.81 0.927]; omega = linspace(-pi,pi,1000);
H = freqz(b,a,omega);
figure('Units','inches','Position',[0,0,18,6]);
subplot(2,2,1)
plot(omega/pi,abs(H));
title('Magnitude Response')
ylabel('|H(e^{j\omega})|')
xlabel('\omega/\pi')
subplot(2,2,2)
plot(omega/pi,angle(H)/pi)
xlabel('\omega/\pi')
ylabel('\angle{H(e^{j\omega})}')
title('Phase Response')
subplot(2,2,3)
plot(omega/pi,unwrap(angle(H))/pi)
xlabel('\omega/\pi')
ylabel('\angle{H(e^{j\omega})}(radians)')
title('Unwrapped Phase Response: \Psi(\omega)')
H_pd = phasedelay(b,a,omega);
subplot(2,2,4)
plot(omega/pi,H_pd/pi)
xlabel('\omega/\pi')
ylabel('\angle{H(e^{j\omega})}')
title('Phase Delay Response')
%% 
% 
%% Problem 3.4
% *Text Problem 5.40 (Page 285)* 
% Consider a second-order IIR notch filter specification that satisfies the 
% following requirements: (1) the magnitude response has notches at $\omega_{1,2}=\pm2\mathrm{\pi}/3$; 
% (2) The maximum magnitude response is 1; (3) the magnitude response is approximately 
% $1/\sqrt{2}$ at frequencies $\omega_{1,2}\pm0.01$.
% 
% 
% 
% *(a)* Using the pole-zero placement approach determine locations of two poles 
% and two zeros of the required filter and then compute its system function $H(z).$
% 
% *Solution*: 
% 
% The formula for a Notch Filter is 
% 
% $$H(z) = b_0\frac{(1-e^{j\omega_0}z^{-1})(1-e^{-j\omega_0}z^{-1})}{(1-re^{j\omega_0}z^{-1})(1-re^{-j\omega_0}z^{-1})$$
% 
% Because we want to block the frequencies at $\pm\frac{2\pi}{3}$, we set are 
% zeros to $z_{1} = e^{j2\pi/3}$ and $z_2 = e^{-j2\pi/3$
% 
% Then if we want the maximum magnitude response to be equal to 1 then $b_0 
% $ must be equal to $\frac{1}{|H(e^{j\omega})|$
% 
% Finally, for the magnitude response to be approximately $\frac{1}{\sqrt{2}$ 
% at frequencies $\omega_{1,2}\pm0.01$, we set our poles equal to $p_1 = 0.969e^{j\frac{2\pi}{3}} 
% \text{ and }  0.969e^{-j\frac{2\pi}{3}} $
% 
% Which resolves into the system function:
% 
% $$H(z) = b_0\frac{1-2\cos\left(\frac{2\pi}{3}\right)z^{-1} + z^{-2}}{1-\left(1.938\cos\left(\frac{2\pi}{3}\right)\right)z^{-1} 
% +0.939z^{-2}$$
% 
% 
% 
% *(b)* Graph the magnitude response of the filter and verify the given requirements.
% 
% *MATLAB script*: 

clc; close all; clear;
w0 = 2*pi/3; r1 = 0.969; %1/sqrt(2); %r2 = 1/sqrt(2);
 % a1 = [1 -r*exp(1i*0.01)]; a2 = [1 -r*exp(-1i*0.01)]; a = a.*conv(a1,a2);
%b1 = [1 -exp(1i*w0)]; b2 = [1 -exp(-1i*w0)]; a1 = [1 -r*exp(1i*0.01)]; a2 = [1 -r*exp(-1i*0.01)]; b = conv(b1,b2); a = conv(a1,a2);
omega = linspace(-pi,pi,1000);
b = [1 -2*cos(w0) 1]; %b2 = [1 -2*cos(0) 1]; b = conv(b1,b2);
a = [1 -(2*r1)*cos(w0) r1^2]; %a2 = [1 -(2*r2)*cos(0) r2^2]; a = conv(a1,a2); 
H = freqz(b,a,omega);
b0 = 1/max(abs(H)); H = b0.*H;
figure
zplane(b,a), title('Pole-Zero Pattern of System')
figure
plot(omega/pi, abs(H))
xticks(-1:1/3:1), xlabel('\omega/\pi'), ylabel('|H(e^{j\omega})|'), title('Magnitude Response')
xticklabels({'-\pi','-2\pi/3','\pi/3','0','\pi/3','2\pi/3','\pi'});
%% 
% *Verification:* We will probe the magnitude responses at $\omega_{1,2}\pm0.01$ 
% and verify that it is approximates $1/\sqrt{2}=0.7071.$

figure
x_pos = (2/3) + 0.01; x_min = (2/3) - 0.01; y_pos = 1/sqrt(2);
plot(omega/pi, abs(H)), hold on
plot(x_pos,y_pos,'r*'), plot(x_min,y_pos,'r*')
plot(-x_pos,y_pos,'r*'), plot(-x_min,y_pos,'r*')
xticks(-1:1/3:1), xlabel('\omega/\pi'), ylabel('|H(e^{j\omega})|'), title('Magnitude Response at \omega_{1,2} \pm 0.01')
xticklabels({'-\pi','-2\pi/3','\pi/3','0','\pi/3','2\pi/3','\pi'});
hold off;
xlim([-0.734 0.734])
ylim([0.6724 0.7378])
%% 
% In the above plot, we can see the points at $\omega_{1,2} \pm 0.01$ where 
% the values of $\frac{1}{\sqrt{2}}$ are being highlighted by red asterisks. This 
% confirms that the magnitude response is equal to $\frac{1}{\sqrt{2}}$ at $\omega_{1,2} 
% \pm 0.01$.
% 
% 
% 
% *(c)* Graph phase and group-delay responses in one plot.
% 
% *MATLAB script*: 

H_phase = phasez(b,a,omega);
figure('Units','inches','Position',[0,0,10,5]);
plot(omega/pi,H_phase/pi, '--b')
xlabel('\omega/\pi')
ylabel('\angle{H(e^{j\omega})}(radians)')
title('Phase & Group Delay Response')
%xticklabels({'-\pi','-4\pi/5','-3\pi/5','-2\pi/5','-1\pi/5','0','1\pi/5','2\pi/5','3\pi/5','4\pi/5','\pi'});
hold on
H_gd = grpdelay(b,a,omega)';
plot(omega/pi,H_gd/pi,'-b'), legend('Principle Value', '\Psi (\omega)'), hold off
%% 
% 
%% Problem 3.5
% *Text Problem 5.55, parts (c) and (d) (Page 288)* 
% Determine the system function, magnitude response, and phase response of the 
% following systems and use the pole-zero pattern to explain the shape of their 
% magnitude response.

clc; close all; clear;
%% 
% 
% 
% *(c)* $y[n] = x[n] - x[n-4] - 0.6561y[n-4]$
% 
% *Solution*: 
% 
% Using the difference equation to derive the system function, we find that
% 
% $$y[n] + 0.6561y[n-4] = x[n] - x[n-4] \rightarrow \frac{Y(z)}{X(z)} = H(z) 
% = \frac{1 - z^{-4}}{1 + 0.6561z^{-4}}$$
% 
% The system function can then be described by substituting $z = e^{j\omega}$ 
% from the above equation yielding:
% 
% $$H(e^{j\omega}) = \frac{1 - e^{-j4\omega}}{1 + 0.6561e^{-j4\omega}}$$
% 
% *MATLAB script for various system responses and plots*: 

b = [1 0 0 0 -1]; a = [1 0 0 0 0.6561]; omega = linspace(-pi,pi,1000);
H = freqz(b,a,omega); figure; zplane(b,a);
figure('Units','inches','Position',[0,0,12,4]);
subplot(1,2,1)
plot(omega/pi,20*log10(abs(H)));
title('Magnitude Response')
ylabel('|H(e^{j\omega})|(dB)')
xlabel('\omega/\pi')
ylim([-52 20])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
subplot(1,2,2)
plot(omega/pi,angle(H)/pi)
xlabel('\omega/\pi')
ylabel('\angle{H(e^{j\omega})} (Radians)')
title('Phase Response')
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
%% 
% *Explanation of plots:* 
% 
% From the pole-zero plots, we see that there are 4 poles at $\pm\frac{\pi}{4} 
% \text{ and}\pm\frac{3\pi}{4}$ with a magnitude of 0.6561 for each pole location. 
% This correlates to these frequencies being boosted in the magnitude response 
% by a gain of about 15dB at each pole location. There were also 4 zero locations 
% in 
% 
% the pole-zero plot located at $\left[0,\frac{\pi}{2},\pi, \frac{3\pi}{2}\right]$, 
% which results in a notch in the magnitude response at these corresponding frequencies, 
% as seen in the magnitude response plot above. 
% 
% 
% 
% *(d)* $y[n]=x[n]-x[n-1]+0.99y[n-1]-0.9801y[n-2]$
% 
% *Solution*: 
% 
% Using the difference equation to derive the system function, we find that
% 
% $$y[n] -0.99y[n-1] +0.9801y[n-2] = x[n] -x[n-1] \rightarrow \frac{Y(e^{j\omega})}{X(e^{j\omega})}  
% = H(e^{j\omega})$$
% 
% $$H(e^{j\omega}) = \frac{1 - e^{-j\omega}}{1 - 0.99e^{-j\omega}+0.9801e^{-j2\omega}$$
% 
% *MATLAB script for various system responses and plots*: 

b = [1 -1]; a = [1 -0.99 0.9801]; omega = linspace(-pi,pi,1000);
H = freqz(b,a,omega); figure; zplane(b,a);
figure('Units','inches','Position',[0,0,12,4]);
subplot(1,2,1)
plot(omega/pi,abs(H));
title('Magnitude Response')
ylabel('|H(e^{j\omega})|')
xlabel('\omega/\pi')
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
subplot(1,2,2)
plot(omega/pi,angle(H)/pi)
xlabel('\omega/\pi')
ylabel('\angle{H(e^{j\omega})} (Radians)')
title('Phase Response')
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
%% 
% *Explanation of plots:* 
% 
% The pole-zero plots show that there are two pole locations that have a magnitude 
% of about 0.99 at frequencies $\omega_{1,2} = \pm\frac{\pi}{3}$, which results 
% in the magnitude response have a large and sharp boost at these particular frequencies 
% in the magnitude response plot. There is also a zero that is 
% 
% present on the unit circle at an angle of $\omega = 0$, which corresponds 
% to a notch at this frequency in the magnitude response. This zero placement
% 
% explains the large notch in magnitude when $\omega = 0$, respectively.
% 
% 
%% Problem 3.6
% An ideal highpass filter is described in the frequency-domain by
% 
% $$H_\mathrm{d}( \mathrm{e}^{\j\omega }) =\left\{\begin{array}{c@{\quad}c}1\cdot 
% \mathrm{e}^{-\j\alpha \omega }, & \omega _\mathrm{c}<\left| \omega \right| \le 
% \mathrm{\pi}  \\0, & \left| \omega \right| \le \omega _\mathrm{c} \end{array}\right.$$
% 
% where $$\omega _\mathrm{c}$$ is the cutoff frequency and $\alpha$ is called 
% the phase delay.
% 
% 
% 
% *(a)* Determine the ideal impulse response $h_\mathrm{d}[n]$ using the DTFT 
% synthesis equation.
% 
% *Solution:* 
% 
% The ideal impulse response $h_d[n]$ can be extracted using the IDTFT equation:
% 
% $$h_d[n] = \frac{1}{2\pi}\int_{-\pi}^{\pi}H_d(e^{j\omega})e^{j\omega n}d\omega 
% $$ 
% 
% $$h_d[n] = \frac{1}{2\pi}\left[\int_{-\pi}^{-\omega_c}e^{-j\alpha\omega}e^{j\omega 
% n}d\omega + \int_{\omega_c}^{\pi}e^{-j\alpha\omega}e^{j\omega n}d\omega\right]$$
% 
% $$h_d[n] = \frac{1}{2\pi}\left[\int_{-\pi}^{-\omega_c}e^{j\omega(n-\alpha)}d\omega 
% + \int_{\omega_c}^{\pi}e^{j\omega(n-\alpha)}d\omega\right]$$
% 
% $$h_d[n] = \frac{1}{2\pi}\left(\left[\frac{e^{j\omega(n-\alpha)}}{j(n-\alpha)}\right]_{\omega 
% = -\pi}^{\omega = -\omega_c} + \left[\frac{e^{j\omega(n-\alpha)}}{j(n-\alpha)}\right]_{\omega 
% = \omega_c}^{\omega = \pi}\right)$$
% 
% $$h_d[n] = \frac{1}{\pi(n-\alpha)}\left[\frac{e^{-j\omega_c(n-\alpha)} - e^{-j\pi(n-\alpha)} 
% + e^{j\pi(n-\alpha)} - e^{j\omega_c(n-\alpha)}}{2j}\right]$$
% 
% $$h_d[n] = \frac{1}{\pi(n-\alpha)}\left[\frac{e^{j\pi(n-\alpha)} - e^{-j\pi(n-\alpha)}}{2j} 
% - \frac{e^{j\omega_c(n-\alpha)} - e^{-j\omega_c(n-\alpha)}}{2j}\right]$$
% 
% $$h_d[n] = \frac{\sin\pi(n-\alpha)-\sin\omega_c(n-\alpha)}{\pi(n-\alpha)$$
% 
% 
% 
% *(b)* Determine and plot the truncated impulse response $h[n] =\left\{\begin{array}{c@{\quad}c}h_\mathrm{d}[n] 
% , & 0\le n\le N-1 \\0, & \hbox{otherwise}\end{array}\right.,$ for $$N=31$$, 
% $$\alpha =15$$, and $\omega _\mathrm{c}=0.5\mathrm{\pi}.$
% 
% *MATLAB script for computation and plot*: 

clc; close all; clear;
N = 31; alpha = 15; wc = 0.5*pi; n = 0:N-1;
a = sin((pi).*(n-alpha)); b = sin(wc.*(n-alpha)); c = pi.*(n-alpha);
h = (a-b)./(c); h(n==alpha)=(pi-wc)/pi;
figure
stem(n,h,'b')
xlabel('\it{n}')
ylabel('\it{h[n]}')
title('Truncated Impulse Response: \it{h[n]}')
%% 
% *(c)* Determine and plot the frequency response function $H( \mathrm{e}^{\j\omega 
% })$, and compare it with the ideal highpass filter response $H_\mathrm{d}( \mathrm{e}^{\j\omega 
% }).$ 
% 
% *MATLAB script*: 

omega = linspace(-pi,pi,1000); H = freqz(h,1,omega); 
figure
plot(omega/pi,abs(H),'b'); 
title('Magnitude Response of H(e^{j\omega})')
ylabel('|H(e^{j\omega})|')
xlabel('\omega/\pi')
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
%% 
% *Comment on your observations*: 
% 
% We can observe that the frequency response function $H(e^{j\omega})$is experiencing 
% a phenomenon called the Gibb's phenomenon in which the Fourier sums overshoot 
% at a jump discontinuity which results in the rippling effect in the magnitude 
% response. This realized highpass filter cannot have the clear and distinct cutoff 
% boundaries such as those in the ideal highpass filter $H_d(e^{j\omega})$. Since 
% the cutoff frequency is $\omega_c = \pm\frac{\pi}{2}$ our realized highpass 
% filter has a transition band that slowly rises/falls at the cutoff frequency, 
% instead of having a clear and sharp frequency cut that occurs in the ideal highpass 
% filter. The ideal highpass filter would require an infinite amount of harmonics 
% to eliminate the Gibb's phenomenon, which makes it impractical for real-world 
% applications. 
% 
% 
%% Problem 3.7
% *Text Problem 6.22 part (b) (Page 346)* 
% Signal $$x_\mathrm{c}(t)= 3 + 2\sin(16\mathrm{\pi} t) + 10\cos(24\mathrm{\pi} 
% t)$$ is sampled at a rate of $$F_\mathrm{s}=20$ Hz to obtain the discrete-time 
% signal $x[n].$ 

clc; close all; clear;
%% 
% 
% 
% *(i)* Determine the spectra $X(\mathrm{e}^{\j\omega})$of $x[n].$
% 
% *Solution:* 
% 
% The CTFT of our original signal, $x_c(t)$is described as
% 
% $$X_c(j2\pi F) = 3\delta - j\delta(F+8)+j\delta(F-8)+5\delta(F+12)+5\delta(F-12)$$
% 
% Our new sampled signal, can be obtained by 
% 
% $$x[n] = x_c(nT) =x_c(n/F) = 3 + 2\sin(16\pi n/20)+10\cos(24\pi n/20)$$
% 
% $$= 3 + 2\sin(0.8\pi n)+10\cos(1.2\pi n) = 3 + 2\sin(0.8\pi n)+10\cos(2\pi 
% n - 0.8\pi n) $$
% 
% $$= 3 + 2\sin(0.8\pi n)+10\cos(-0.8\pi n)  \rightarrow 3 + 2\sin(0.8\pi n)+10\cos(0.8\pi 
% n) $$
% 
% Thus, $x[n] = 3 + 2\sin(0.8\pi n)+10\cos(0.8\pi n) $
% 
% The spectra $X(e^{j\omega})$ of $x[n]$ can be obtained from the equation:
% 
% $$X_s(F) = F_s \sum X_c(F - kF_s) = 20\sum X_c(F-20k)$$
% 
% $$= 20 \sum\left[3\delta - j\delta(F+8-20k)+j\delta(F-8-20k) + 5\delta(F+8-20k)+5\delta(F-8-20k)\right]$$
% 
% $$= 20 \sum\left[3\delta +(5- j)\delta(F+8-20k)+(5+j)\delta(F-8-20k)\right]$$
% 
% 
% 
% *(ii)* Plot magnitude of $X(\mathrm{e}^{\j\omega})$ as a function of $$\omega$$ 
% in $$\frac{\text{rad}}{\text{sam}}$$ and as a function of $$F$$ in Hz.
% 
% *Solution:* You can hand plot or use a plotting app or use MATLAB for the 
% magnitude plot.
% 
% 
% 
% *MATLAB script for the magnitude plot*: You can use MATLAB to numerically 
% compute and plot the magnitude $\vert X(\mathrm{e}^{\j\omega})\vert$over $-\mathrm{\pi}<\omega\leq\mathrm{\pi}$. 
% The amplitudes will be approximate since we cannot synthesize a true impulse 
% but should be proportional to the areas under the impulses. Therefore, plot 
% normalized magnitudes. The important issues here are the location of spikes 
% and their proportional heights.

F = 20; Ts = 1/20; n = -10:10; w1 = 16*pi; w2 = 24*pi;
% Discrete Time Signal 
x = 3 + 2*sin(w1*(n*Ts)) + 10*cos(w2*(n*Ts));
% Discrete-time Fourier Transform
K = 500; k = 0:1:K; w = pi*k/K;
X = x * exp(-1i*n'*w); X = real(X);
w = [-fliplr(w), w(2:K+1)]; X = [fliplr(X), X(2:K+1)];
figure
plot(w/pi,abs(X/F)-0.65)
title('Magnitude Response of X(e^{j\omega})')
ylabel('|X(e^{j\omega})|')
xlabel('\omega/\pi')
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
%% 
% 
% 
% *(iii)* Explain if $x_\mathrm{c}(t)$ can be recovered from $x[n].$
% 
% *Explanation:* 
% 
% If we were to derive a reconstructed signal $y_c(t)$ from the above figure 
% after lowpass filtering, we would obtain
% 
% $$y_c(t) = 3 + 2\sin(16\pi t) + 10\cos(16\pi t) \ne x_c(t)$$
% 
% Because our sampling frequency, $F_s = 20 Hz$ did not meet the Nyquist rate 
% of twice the highest frequency component, the reconstructed signal experienced 
% aliasing in the $10\cos(24\pi t)$component and altered the original signal after 
% sampling. If the sampling frequency was $F_s > 24Hz$ we would have been able 
% to recover $x_c(t)$ from our sample signal $x[n].$ 
% 
% 
%% Problem 3.8
% *Text Problem 6.26 (Page 347)* 
% Consider a continuous-time signal $x_\mathrm{c}(t)= 3\cos(2\mathrm{\pi} F_{1} 
% t + 45^{\deg})+3\sin(2\mathrm{\pi} F_{2} t).$It is sampled at $$t=0.001n$$ to 
% obtain $$x[n]$$ which is then applied to an ideal DAC to obtain another continuous-time 
% signal $y_{\mathrm{r}}(t).$

clc; close all; clear;
%% 
% 
% 
% *(a)* For $$F_{1}=150$$ Hz and $$F_{2}=400$$ Hz, determine $$x[n]$$ and graph 
% its samples along with the signal $x_\mathrm{c}(t)$ in one plot (choose few 
% cycles of the $x_\mathrm{c}(t)$ signal).
% 
% *Solution:* 
% 
% $$x_c(t) = 3\cos(300\pi t + \pi/4) + 3\sin(800\pi t)$$
% 
% $$x[n] = x_c(nT) = 3\cos\left(\frac{300\pi n}{1000} + \frac{\pi}{4}\right) 
% + 3\sin\left(\frac{800\pi n}{1000}\right)$$
% 
% $$x[n] = 3\cos\left(0.3\pi n+ \frac{\pi}{4}\right) + 3\sin\left(0.8\pi n\right)$$

F1 = 150; F2 = 400; w1 = 2*pi*F1; w2 = 2*pi*F2;
Fs = 1000; Ts = 1/Fs; n = -20:20;
% Continuous-time Signal 
Dt = 0.00001; t = -0.02:Dt:0.02; xc = 3.*cos(w1*t + pi/4) + 3*sin(w2*t);
% Discrete Time Signal 
xn = 3.*cos(w1*(n*Ts) + pi/4) + 3*sin(w2*(n*Ts));
figure('Units','inches','Position',[0,0,12,3]);
plot(t*1000,xc,'--b'), hold on, stem(n*Ts*1000,xn,'r')
xlabel('Time(msec) and Samples(n)')
ylabel('Amplitude')
title('CT signal and Sampled DT signal - Ts = 1 msec'), legend('x_c(t)','x[n]')
%% 
% 
% 
% *(b)* Determine $y_{\mathrm{r}}(t)$ for the above $$x[n]$$ as a sinusoidal 
% signal. Graph and compare it with $x_\mathrm{c}(t).$
% 
% *Solution:* 
% 
% $$y_r(t) = x[n]|_{n = tFs} = x[n]|_{n=1000t}$$
% 
% $$y_r(t) = 3\cos(0.3\pi(1000t) + \pi/4) + 3\sin(0.8\pi(1000t))$$
% 
% $$y_r(t) = 3\cos(300\pi t + \pi/4) + 3\sin(800\pi t) = x_c(t)$$
% 
% Thus, no aliasing occurred during the reconstruction process.
% 
% *MATLAB script*: 

Fs = 1000; Ts = 1/Fs; n = -20:20; t= -0.02:Dt:0.02; nTs = n*Ts;
% Reconstructed Signal 
yr = xn * sinc(Fs*(ones(length(n),1)*t-nTs'*ones(1,length(t))));
figure('Units','inches','Position',[0,0,12,3]);
plot(t*1000,yr,'b'), hold on, stem(n,xn,'--r')
xlabel('Time (msec)')
ylabel('Amplitude')
title('Sampled and Reconstructed Signals - Ts = 1 msec'), legend('y_r(t)','x[n]')
figure('Units','inches','Position',[0,0,12,3]);
plot(t*1000,yr,'b'), hold on, plot(t*1000,xc,'--b')
xlabel('Time (msec)')
ylabel('Amplitude')
title('Original CT signal and Reconstructed Signal'), legend('y_r(t)','x_c(t)')
%% 
% 
% 
% *(c)* Repeat (a) and (b) for $$F_{1}=300$$ Hz and $$F_{2}=700$$ Hz. 
% 
% *Solution*: 
% 
% $$x_c(t) = 3\cos(600\pi t + \pi/4) + 3\sin(1400\pi t)$$
% 
% $$x[n] = x_c(nT) = 3\cos\left(\frac{600\pi n}{1000} + \frac{\pi}{4}\right) 
% + 3\sin\left(\frac{1400\pi n}{1000}\right)$$
% 
% $$x[n] = 3\cos\left(0.6\pi n+ \frac{\pi}{4}\right) + 3\sin\left(1.4\pi n\right)$$
% 
% $x[n] = 3\cos\left(0.6\pi n+ \frac{\pi}{4}\right) + 3\sin\left(2\pi n -1.4\pi 
% n\right) \rightarrow 3\sin\left(2\pi n -1.4\pi n\right) = 3\sin\left(-0.6\pi 
% n\right) \rightarrow -3\sin\left(0.6\pi n\right)$ (odd symmetry)
% 
% $$x[n] = 3\cos\left(0.6\pi n+ \frac{\pi}{4}\right) - 3\sin\left(0.6\pi n\right)$$ 
% 
% The high frequency component $$F_{2}=700$$ Hz experiences aliasing because 
% $F_s \ne F_s >2F_2$
% 
% $$y_r(t) = x[n]|_{n = tFs} = x[n]|_{n=1000t}$$
% 
% $$y_r(t) = 3\cos(0.6\pi(1000t) + \pi/4) - 3\sin(0.6\pi(1000t))$$
% 
% $$y_r(t) = 3\cos(600\pi t + \pi/4) - 3\sin(600\pi t) \ne x_c(t)$$
% 
% Thus, aliasing occurred during the reconstruction process because of the sampling 
% rate not meeting the Nyquist rate. 

F1 = 300; F2 = 700; w1 = 2*pi*F1; w2 = 2*pi*F2;
Fs = 1000; Ts = 1/Fs; n = -20:20; nTs = n*Ts;
% Continuous-time Signal 
Dt = 0.00001; t = -0.02:Dt:0.02; xc = 3.*cos(w1*t + pi/4) + 3*sin(w2*t);
% Discrete Time Signal 
xn = 3.*cos(w1*(n*Ts) + pi/4) + 3*sin(w2*(n*Ts));
figure('Units','inches','Position',[0,0,12,3]);
plot(t*1000,xc,'--b'), hold on, stem(n*Ts*1000,xn,'r')
xlabel('t in msec')
ylabel('Amplitude')
title('CT signal and Sampled DT signal - Ts = 1 msec'), legend('x_c(t)','x[n]')
%% 
% *MATLAB script*: 

% Reconstructed Signal 
yr = xn * sinc(Fs*(ones(length(n),1)*t-nTs'*ones(1,length(t))));
figure('Units','inches','Position',[0,0,12,3]);
plot(t*1000,yr,'b'), hold on, stem(n,xn,'--r')
xlabel('t in msec')
ylabel('Amplitude')
title('Sampled and Reconstructed Signals - Ts = 1 msec'), legend('y_r(t)','x[n]')
figure('Units','inches','Position',[0,0,12,3]);
plot(t*1000,yr,'b'), hold on, plot(t*1000,xc,'--r')
xlabel('t in msec')
ylabel('Amplitude')
title('Original CT signal and Reconstructed Signal'), legend('y_r(t)','x_c(t)')
%% 
% *Comment on your results*:
% 
% Analyzing both the reconstructed plots, for parts (a) and (b), the sampling 
% frequency $F_s = 1000$ satisfied the Nyquist rate, which was able to recover 
% the original continuous signal after reconstruction. However, for parts (c) 
% and (d), the sampling rate was not twice the value of the highest frequency 
% component $F_2 = 700Hz$, which resulted in aliasing after reconstructing the 
% signal from the sampled signal. The higher frequencies being aliased resulted 
% in a lower frequency sinusoid, as seen in the plot above.  
% 
% 
%% Problem 3.9
% *Text Problem 6.23 part (c) (Page 348)* 
% Signal $$x\mathrm{c}(t)= 5\mathrm{e}^{\j40t} + 3\mathrm{e}^{-\j70t}$$ is sampled 
% periodically with $T=0.1$ sec to obtain the discrete-time signal $$x[n].$ Determine 
% the spectra $$X(\mathrm{e}^{\j\omega})$$ of $$x[n]$$ and plot its magnitude 
% as a function of $$\omega$$ in $$\frac{\text{rad}}{\text{sam}}.$ 

clc; close all; clear;
%% 
% *Solution:* You can hand plot or use a plotting app or use MATLAB for the 
% magnitude plot.
% 
% 
% 
% *MATLAB script for the magnitude plot*: You can use MATLAB to numerically 
% compute and plot the magnitude $\vert X(\mathrm{e}^{\j\omega})\vert$over $-\mathrm{\pi}<\omega\leq\mathrm{\pi}$. 
% The amplitudes will be approximate since we cannot synthesize a true impulse 
% but should be proportional to the areas under the impulses. Therefore, plot 
% normalized magnitudes. The important issues here are the location of spikes 
% and their proportional heights.

Ts = 0.1; n = 0:100; nTs = n*Ts;
x = 5*exp(40i*nTs) + 3*exp(-70i*nTs);
w = -pi:pi/1000:pi;
X = x*exp(-1i*n'*w);
X_mag = abs(X);
X_phase = angle(X);

subplot(2,1,1)
plot(w/pi, X_mag/10);
grid on
xlabel('\omega/\pi (rad/sample)');
ylabel('|X(e^{j\omega})|');
title('Magnitude Spectrum');
subplot(2,1,2)
plot(w/pi, X_phase);
grid on
xlabel('\omega/\pi (rad/sample)');
ylabel('\angle X(e^{j\omega}) (radians)');
title('Phase Spectrum');
%% 
% *Explain if* $$x_\mathrm{c}(t)$$ *can be recovered from* $x[n]$: 
% 
% The frequency components of the signal $x_c(t)$ are $F_1 = \frac{20}{\pi} 
% \approx 6.4Hz \text{ and } F_2 = \frac{35}{\pi} \approx 11.1Hz$. With our sampling 
% interval being $T_S = 0.1\text{s}$ this means our sampling frequency is $F_s 
% = 10Hz$. So, aliasing will occur during the sampling process and will ultimately 
% lose the $F_2$ component because the sampling interval is too low. The lower 
% frequency component $F_1$ will be sampled, but aliasing will occur that will 
% make $x_c(t)$ be impossible to reconstruct from our sample signal $x[n]$.
% 
% 
%% Problem 3.10
% *Text Problem 6.37 (Page 348)* 
% Signal $x_\mathrm{c}(t)$ with spectra $X(\j2\mathrm{\pi}F)$ shown below is 
% sampled at a rate of $$F_\mathrm{s}=6$$ Hz to obtain the discrete-time signal 
% $$x[n].$
% 
% 

clc; close all; clear;
%% 
% 
% 
% *(i)* Determine the spectra $$X(\mathrm{e}^{\j\omega})$$ of $$x[n].$ 
% 
% *Solution:* 
% 
% First, we see that $X_c(j2\pi F) = 1-\frac{|F|^2}{25}, \quad \text{for }|F| 
% \leq5$ and 0, everywhere else. 
% 
% We can relate the spectra $X_c(j2\pi F)$ to $X_c(j\Omega)$ since $\Omega = 
% 2\pi F$, then using the equation
% 
% $X(e^{j\omega}) |_{\omega = \Omega T} = \frac1T \sum_{k = -\infty}^{\infty} 
% X_c\left(j\Omega - j\frac{2\pi}{T}k\right)$ along with our sample interval $F_s 
% = 6Hz \rightarrow T_s = 1/F_s = \frac16$
% 
% we can approximate our new sampled spectra $X(e^{j\omega})$ using the substitution 
% $\omega = \Omega T$ to get:
% 
% $$X(e^{j\omega}) = 6 \sum_{k = -\infty}^{\infty}X_c[j(6\omega - \frac{\pi}{3}k)]$$
% 
% where 
% 
% $$X(e^{j\omega}) = 6-\frac{|\frac{\omega}{2\pi}|^2}{25}, \quad \text{for  
% } |6\omega| \leq 10\pi$$
% 
% 
% 
% *(ii)* Plot it as a function of $$\omega$$ in rad/sam. 
% 
% *Solution:* You can hand plot or use a plotting app or use MATLAB for the 
% magnitude plot.
% 
% 
% 
% The plot of $X(e^{j\omega})$results in an overlapped version of the original 
% spectra. The highlighted rippling top-edges of the plot are the resulting spectrum 
% shape. 
% 
% 
% 
% *(iii)* Explain if $x_\mathrm{c}(t)$ can be recovered from $$x[n].$
% 
% *Solution:* 
% 
% We would not be able to recover the original signal $x_c(t)$ from the sampled 
% signal $x[n]$, because the sampling frequency $F_s = 6Hz$ did not meet the Nyquist 
% rate, which should have been $\geq 10Hz$. The original spectra starts to overlap 
% because of the sampling process around $3Hz$, where edges of the semi-circles 
% were summed together and produced a new aliased spectra shown in $X(e^{j\omega})$ 
% above. 
% 
%