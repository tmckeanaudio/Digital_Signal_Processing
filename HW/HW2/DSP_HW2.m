%% *EECE5666 (DSP) : Homework-2*
% *Due on February 8, 2022 by 11:59 pm via submission portal.* 
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
%% Problem 2.1
% *Text Problem 3.45, parts (b) and (d) only (Page 130)* 
% Determine the $z$-transform and sketch the pole-zero plot with the ROC for 
% each of the following sequences.

clc; close all; clear;
%% 
% 
% 
% *(b)* $x[n] = (1/2)^n u[n+1]+3^n u[-n-1]$:   
% 
% *Solution*: 
% 
% $$X(z) = \sum_{n = -\infty}^{\infty} \left((1/2)^n u[n+1]+3^n u[-n-1]\right)z^{-n}$$
% 
% $$X(z) = \sum_{n = -\infty}^{\infty} (1/2)^n u[n+1]z^{-n}+ \sum_{n = -\infty}^{\infty} 
% 3^n u[-n-1]z^{-n}$$
% 
% $$X(z) = \sum_{n = -1}^{\infty} \left(\frac{1/2}{z}\right)^n+ \sum_{n = -\infty}^{-1} 
% \left(\frac{z}{3}\right)^{n}$$
% 
% $$X(z) = \frac{1}{1 - 0.5z^{-1}} - \frac{1}{1 - 3z^{-1}} \text{ , } \quad 
% \text{ROC: } \frac{1}{2} < |z| < 3$$
% 
% *MATLAB script for pole-zero and ROC plot*:

p = [0.5;3]; z = [0;0];
figure
zplane(z,p)
%% 
% 
% 
% *(d)* $x[n]=|n|(1/2)^{|n|}$
% 
% *Solution*: 
% 
% $$X(z) = \sum_{n = - \infty}^{\infty}\left(|n|(1/2)^{|n|}\right)z^{-n}$$
% 
% $$X(z) = \sum_{n = -\infty}^{-1} n\Big(\frac12\Big)^{n}z^{-n} + 0 + \sum_{n 
% = 1}^{\infty} n\Big(\frac12\Big)^{n}z^{-n}$$
% 
% $$X(z) = \frac{\frac12z^{-1}}{\Big(1 - \frac12z^{-1}\Big)^{2}} + \frac{\frac12z^{-1}}{\Big(1 
% - \frac12z^{-1}\Big)^{2}} $$
% 
% $$X(z) = \frac{z^{-1}}{1 - z^{-1} + \frac14z^{-2}}$$
% 
% *MATLAB script for pole-zero and ROC plot*:

b = [0 1]; a = [1 -1 1/4];
[A,p,C] = residuez(b,a)
zplane(b,a)
%% 
% The ROC for this sequence would then be $\text{ROC: } |z|> \frac12  $
% 
% 
%% Problem 2.2
% Consider the $z$-transform expression:
% 
% $$X(z) = \frac{(z-0.91)\bigl(z^2+0.3z+0.4\bigr)} {\bigl(z+1.5 \bigr) \bigl(z^2-0.6z+0.6\bigr)}.$$

clc; close all; clear;
%% 
% 
% 
% *(a)* Provide a zero-pole plot of $$X(z)$$.
% 
% *MATLAB script*: 

b = [1 -0.61 0.127 -0.364]; a = [1 0.9 -0.3 0.9];
z = roots(b), p = roots(a)
figure
zplane(z,p)
%% 
% 
% 
% *(b)* List all possible regions of convergence (ROCs) for this $z$-transform.
% 
% *Solution:* 
% 
% Given the calculated poles at $p_1 = -1.5, p_2  = 0.3 + j0.7141, \text{and 
% } p_3 = 0.3 - j0.7141$
% 
% There are 3 possible regions of convergence being: 
% 
% $$ \text{ROC}_1 = |z| > 1.5$$
% 
% $$ \text{ROC}_2 = 0.3 < |z| < 1.5$$
% 
% $$ \text{ROC}_3 = |z| < 0.3$$
% 
% 
% 
% *(c)* Determine the inverse $z$-transform so that the resulting sequence is 
% absolutely summable. This sequence $x[n]$ should be a real-valued sequence. 
% Provide a |*stem*| plot of $x[n]$.
% 
% *Solution* 
% 
% *MATLAB script for sequence plot*:

[A,p,C] = residuez(b,a)
Ma = (abs(A))', Mp = (abs(p))'
Aa = (angle(A)), Ap = (angle(p))
%% 
% $X(z) = -0.4044 + \frac{0.9426}{\Big(1 + \frac32z^{-1}\Big)} + \frac{\Big(0.2309 
% + j0.1643\Big)}{1 - (0.3+j0.7141)z^{-1}} + \frac{\Big(0.2309 - j0.1643\Big)}{1 
% - (0.3-j0.7141)z^{-1}} $  or 
% 
% $$X(z) = -0.4044 + \frac{0.9426}{\Big(1 + \frac32z^{-1}\Big)} + \frac{0.2834e^{j35.4^\circ}}{1 
% - 0.78e^{j67.2^\circ}z^{-1}} + \frac{0.2834e^{-j35.4^\circ}}{1 - 0.78e^{-j67.2^\circ}z^{-1}} 
% $$
% 
% After computing the magnitude and phases of the complex poles and coefficients, 
% we get the input sequence:
% 
% $$x[n] = -0.4044\delta[n] + 0.9426\Big(-\frac32\Big)^{n}u[-n-1] + 0.56\Big(0.78\Big)^{n}\cos\Big(67.2^\circ 
% n + 35.4^\circ\Big)u[n]$$

n = 0:50;
x1 = (n==0); x2 = (n>=0); x3 = (n<=-1);
x = -0.4044*x1 + 0.9426.*((-3/2).^(n)).*x3 + 0.56.*((0.78).^(n)).*cos(Ap(2).*n + Aa(2)).*x2;
figure
stem(n,x)
xlabel("\it{n}")
ylabel("\it{x[n]}")
title("Input sequence: \it{x[n]}")
%% 
% 
%% Problem 2.3 
% *Text Problem 3.47, parts (b) and (e), (Page 131)* 
% Given the $z$-transform pair $x[n]\leftrightarrow X(z)=z^{-1}/(1+0.8 z^{-1})$with 
% ROC: $|z|>0.8$, use the $z$-transform properties to determine the $z$-transform 
% of the following sequences:
% 
% 
% 
% *(b)* $y[n]=x[3-n] = x[-(n-3)]$: 
% 
% *Solution*: 
% 
% The output signal $y[n]$ exhibits properties of the Time-Shift and Folding 
% operations performed on the input signal $x[n]$ where,
% 
% $x[-n] \leftrightarrow X(1/z) = \frac{(z^{-1})^{-1}}{1 + 0.8(z^{-1})^{-1}} 
% = \frac{z}{1+0.8z}$ would be the Folding operation performed on $x[n]$ and then 
% 
% $x[n-k] \leftrightarrow z^{-k}X(z) = z^{3}\frac{z}{1+0.8z}  = \frac{z^4}{1+0.8z$ 
% with $\text{ROC} = |z| < 0.8$
% 
% 
% 
% *(e)* $y[n]=x[n]\ast x[2-n]$:
% 
% *Solution:* 
% 
% $x[2-n]$ is a folded and shifted sequence of the input signal $x[n]$ and can 
% be derived from the following operations:
% 
% $x[2 -n] \leftrightarrow z^2X(1/z) = z^2\frac{(z^{-1})^{-1}}{1 + 0.8(z^{-1})^{-1}} 
% = \frac{z^3}{1+0.8z}$   then the operation 
% 
% $y[n]=x[n]\ast x[2-n] \leftrightarrow X(z)X(1/z)z^2$ which means the z-transform 
% is the product of the z-transforms for both input sequences giving:
% 
% $$\left(\frac{z^{-1}}{1 + 0.8z^{-1}}\right)\left(\frac{z^{3}}{1 + 0.8z^}\right) 
% = \frac{z^{2}}{(1 + 0.8z^{-1})(1+0.8z)}$$
% 
% 
%% Problem 2.4
% An LTI system described by the following impulse response
% 
% $$h[n] = n\Bigl(\frac13\Bigr)^{n}u[n]+\Bigl(-\frac14\Bigr)^{n}u[n].$$

clc; close all; clear;
%% 
% 
% 
% *(a)* Determine the system function representation.
% 
% *Solution*: 
% 
% The system function is defined as $H(z) = \sum_{n = -\infty}^{\infty} h[n]z^{-n}$ 
% so,
% 
% $$H(z) = Z\{h[n]\} = Z\{n\Bigl(\frac13\Bigr)^{n}u[n]\}+Z\{\Bigl(-\frac14\Bigr)^{n}u[n]\}$$
% 
% $H(z) = \frac{\frac13z^{-1}}{\Big(1-\frac13z^{-1}\Big)^{2}} + \frac{1}{\Big(1 
% + \frac14z^{-1}\Big)$     which can be expanded to 
% 
% $$H(z) = \frac{\frac13z^{-1}}{\Big(1-\frac23z^{-1}+\frac19z^{-2}\Big)} + \frac{1}{\Big(1 
% + \frac14z^{-1}\Big)$$ 
% 
% then in order to make it a single rational fraction we multiple the denominators 
% together and the numerators by the other's denominator to get
% 
% $H(z) = \frac{\frac13z^{-1}\Big(1-\frac23z^{-1}+\frac19z^{-2}\Big) + 1\Big(1-\frac23z^{-1}+\frac19z^{-2}\Big)}{\Big(1-\frac23z^{-1}+\frac19z^{-2}\Big)\Big(1 
% + \frac14z^{-1}\Big)} = \frac{1 - \frac13z^{-1} + \frac{7}{36}z^{-2}}{1- \frac{5}{12}z^{-1} 
% - \frac{1}{18}z^{-2} + \frac{1}{36}z^{-3}$ hence,
% 
% $$H(z) = \frac{1 - \frac13z^{-1} + \frac{7}{36}z^{-2}}{1- \frac{5}{12}z^{-1} 
% - \frac{1}{18}z^{-2} + \frac{1}{36}z^{-3}$$
% 
% 
% 
% *(b)* Determine the difference equation representation.
% 
% *Solution*: 
% 
% $$H(z) = \frac{B(z)}{A(z)} = \frac{1 - \frac13z^{-1} + \frac{7}{36}z^{-2}}{1- 
% \frac{5}{12}z^{-1} - \frac{1}{18}z^{-2} + \frac{1}{36}z^{-3}$$
% 
% $$Y(z) = - \sum_{k = 1}^{N} a_kz^{-k}Y(z) + \sum_{k = 0}^{M} b_kz^{-k}X(z)$$
% 
% $Y(z) = \frac{5}{12}z^{-1}Y(z) + \frac{1}{18}z^{-2}Y(z) - \frac{1}{36}z^{-3}Y(z) 
% + X(z) - \frac13z^{-1}X(z) + \frac{7}{36}z^{-2}X(z)$ then taking the inverse 
% z-transform to solve for $y[n] $ we get:
% 
% $$y[n] =  \frac{5}{12}y[n-1] + \frac{1}{18}y[n-2]- \frac{1}{36}y[n-3] + x[n] 
% - \frac13x[n-1] + \frac{7}{36}x[n-2]$$
% 
% 
% 
% *(c)* Determine the pole-zero plot.
% 
% *MATLAB script*:  

b = [1, -1/3, 7/36]; a = [1, -5/12, -1/18, 1/36];
[A,p,C] = residuez(b,a)
zeros = roots(b), poles = roots(a)
figure
zplane(b,a)
%% 
% Because this an LTI system the ROC will then be $\text{ROC: } |z|>\frac13$
% 
% 
% 
% *(d)* Determine the output sequence $y[n]$ when the input is $$x[n]=\Bigl(\frac14\Bigr)^{n}u[n]$$.
% 
% *Solution*: 
% 
% $X(z) = \frac{1}{1 - \frac14z^{1}$ and $H(z) = \frac{1 - \frac13z^{-1} + \frac{7}{36}z^{-2}}{1- 
% \frac{5}{12}z^{-1} - \frac{1}{18}z^{-2} + \frac{1}{36}z^{-3}$ then 
% 
% $$Y(z) = X(z)H(z) = \Big(\frac{1}{1-\frac14z^{-1}}\Big) \Big(\frac{1 - \frac13z^{-1} 
% + \frac{7}{36}z^{-2}}{1- \frac{5}{12}z^{-1} - \frac{1}{18}z^{-2} + \frac{1}{36}z^{-3}}\Big) 
% = \frac{1 - \frac13z^{-1}+ \frac{7}{36}z^{-2}}{1 - \frac23z^{-1} + \frac{7}{144}z^{-2} 
% + \frac{1}{24}z^{-3} - \frac{1}{144}z^{-4}$$
% 
% We can solve for $y[n] $ by taking the partial fraction expansion of $Y(z) 
% $followed by the inverse z-transform

b = [1, -1/3, 7/36]; a = [1, -2/3, 7/144, 1/24, -1/144];
[A,p,C] = residuez(b,a);
A = A', p = p'
%% 
% Using the partial fraction coefficients along with the calculates poles we 
% get an output signal of 
% 
% $$y[n] = 12.5\Big(\frac14\Big)^{n}u[n] + \frac12\Big(-\frac{1}{4}\Big)^{n}u[n] 
% + 4\Big(\frac13\Big)^{n}u[n] - 16\Big(\frac13\Big)^nu[n]$$
% 
% 
%% Problem 2.5
% *Text Problem 3.57 (Page 132)* 
% Determine the impulse response of the system described by
% 
% $$y[n]+\frac{11}{6}y[n-1]+\frac12y[n-2]=2x[n].$$
% 
% for all possible regions of convergence.
% 
% *Solution:* 
% 
% Taking the z-transform of the expression above we get:
% 
% $$Y(z) + \frac{11}{6}z^{-1}Y(z) + \frac12z^{-2}Y(z) = 2X(z)$$
% 
% $$(1 + \frac{11}{6}z^{-1} + \frac12z^{-2})Y(z) = 2X(z)$$
% 
% $$H(z) = \frac{Y(z)}{X(z)} = \frac{2}{1 + \frac{11}{6}z^{-1} + \frac12z^{-2}$$
% 
% Solving for the partial fraction coefficients and poles we get: 

b = 2; a = [1, 11/6, 1/2];
[A,p,C] = residuez(b,a); A = A', p = p'
zplane(b,a)
%% 
% Which gives us coefficients 2.5714 and -0.5714 with poles $p_1 = \frac13 \text{ 
% and } p_2 = \frac32$
% 
% $H(z) = -\frac{0.5714}{1 + \frac13z^{-1}} + \frac{2.5714}{1 + \frac32z^{-1}$ 
% then by taking the inverse z-transform we get 3 possible regions of convergence
% 
% $$h[n] = 2.5714\Big(-\frac32\Big)^{n}u[n] - 0.5714\Big(-\frac13\Big)^{n}u[n], 
% \quad |z|>\frac32 \text{ (causal)$$
% 
% $$h[n] = -2.5714\Big(-\frac32\Big)^{n}u[-n-1] + 0.5714\Big(-\frac13\Big)^{n}u[-n-1], 
% \quad |z|<\frac13 \text{ (anticausal)$$
% 
% $$h[n] = -2.5714\Big(-\frac32\Big)^{n}u[-n-1] - 0.5714\Big(-\frac13\Big)^{n}u[n], 
% \quad \frac13 <|z|<\frac32 \text{ (two-sided)$$
% 
% 
%% Problem 2.6
% *Text Problem 3.63 (Page 133)* 
% Consider the following LCCDE
% 
% $$y[n] = 2\cos(\omega_0)y[n-1]-y[n-2]$$
% 
% with no input but with initial conditions $y[-1]=0$ and $y[-2] = -A\sin(\omega_0).$
% 
% 
% 
% *(a)* Show that the solution of the above LCCDE is given by $y[n] = A\sin[(n+1)\omega_0] 
% u[n].$ This system is known as a _digital oscillator_.
% 
% *Solution:* 
% 
% $$y[n] = 2\cos(\omega_0)y[n-1] - y[n-2]$$ 
% 
% We can apply the z-transform to the LCCDE to obtain $Y(z) :$
% 
% $$Y(z) - 2\cos(\omega_0)z^{-1}\Big(Y(z) + y[-1]z\Big) + z^{-2}\Big( Y(z) + 
% y[-1]z + y[-2]z^2\Big) = 0$$
% 
% $$(1 - 2\cos(\omega_0)z^{-1} + z^{-2})Y(z) -A\sin(\omega_0) = 0$$ 
% 
% $Y(z) =  \frac{A\sin(\omega_0)}{1 - 2\cos(\omega_0)z^{-1} + z^{-2}$ noticing 
% the numerator doesn't have a $z$ term we can write the z-transform as if it 
% was being time-shifted by +1 so, 
% 
% $Y(z) =  \frac{z^1\big(A\sin(\omega_0)z^{-1}\big)}{1 - 2\cos(\omega_0)z^{-1} 
% + z^{-2}$ which after taking the inverse transform to obtain $y[n]$ we observe
% 
% $y[n] = A\sin[(n+1)\omega_0]u[n]$ which verifies the solution to the LCCDE.
% 
% 
% 
% *(b)* For $A=2$ and $\omega_0=0.1\mathrm{\pi}$, verify the operation of the 
% above digital oscillator using MATLAB.
% 
% *MATLAB script*: 

clc; close all; clear;
A = 2; w0 = pi/10;
a = [1, -A*cos(w0),1]; n = 0:100; x = (n==0);
y = filter(1,a,x);
figure
plot(n,y)
xlabel("\it{n}")
ylabel("\it{y[n]}")
title("Impulse Response of Digital Oscillator for \it{0 \leq n \leq 100}")
%% 
% 
%% Problem 2.7
% *Text Problem 4.38, parts (a) and (d) only, (Page 197)* 
% Determine whether or not each of the following signals is periodic. If a signal 
% is periodic, determine its fundamental period.
% 
% 
% 
% *(a)* $x_1(t)=|\sin(7\mathrm{\pi}t)|\cos(11\mathrm{\pi}t)$
% 
% *Solution:* 
% 
% Since $|\sin(7\pi t)| $ has only positive values, it becomes periodic from 
% $0 <t<\pi$ instead of $2\pi$
% 
% Thus, $x_1(t) = \sin\Big(\frac72\pi t\Big)\cos\Big(11\pi t \Big) $ then we 
% can use a trigonometric identity for the product of these sinusoids giving
% 
% $$x_1(t) = \sin\Big(\frac72\pi t\Big)\cos\Big(11\pi t \Big)  = \frac12\Big[\sin\Big(\frac{29\pi}{2}t\Big) 
% - \sin\Big(\frac{15\pi}{2}t\Big)\Big]$$
% 
% The first sinusoid is periodic with period $T_1 = \frac{2\pi}{\omega_0} = 
% \frac{2\pi}{29\pi/2} = \frac{4}{29}$ and then
% 
% The second sinusoid is periodic with period $T_2 = \frac{2\pi}{\omega_0} = 
% \frac{2\pi}{15\pi/2} = \frac{4}{15}$
% 
% Using the fact that the sum of two periodic signals result in a periodic signal, 
% hence $x_1(t)$ is thus periodic.
% 
% The fundamental period of $x_1(t) $would then be the least common multiple 
% of both $T_1 \text{ and } T_2$
% 
% $$T = LCM(T_1,T_2) = LCM\Big(\frac{4}{29},\frac{4}{15}\Big) = 4$$
% 
% Thus, the fundamental period of $x_1(t) \text{ is } T  = 4$
% 
% 
% 
% *(d)* $x_4[n]=\mathrm{e}^{\j\mathrm{\pi}n/7}+\mathrm{e}^{\j\mathrm{\pi}n/11}:$
% 
% *Solution:* 
% 
% Observing the fundamental frequencies of each of the exponentials of $x_4[n]$ 
% we see
% 
% $\mathrm{e}^{\j\mathrm{\pi}n/7} \rightarrow \omega_0 = \frac\pi7 \rightarrow 
% f = \frac{\omega_0}{2\pi} = \frac{1}{14} \text{Hz}$  and $\mathrm{e}^{\j\mathrm{\pi}n/11} 
% \rightarrow \omega_0 = \frac\pi{11} \rightarrow f = \frac{\omega_0}{2\pi} = 
% \frac{1}{22} \text{Hz}$
% 
% both frequencies are rational values of $\omega_0$ , thus the sequence $x_4[n]$ 
% is periodic. 
% 
% It's fundamental period can be calculated from the Least Common Multiple of 
% each frequency, so 
% 
% $T_1 = 1/F = 14 Hz$ and $T_2 = 1/F = 22 Hz$ then fundamental period would 
% be 
% 
% $$N = LCM(14,22) = 154$$
% 
% Thus, the sequence $x_4[n] $ is periodic with fundamental period $= 154$
% 
% 
%% Problem 2.8
% *Text Problem 4.45, parts (c) and (d) only, (Page 198)* 
% Given that $x[n]$ is a periodic sequence with fundamental period $N$ and Fourier 
% coefficients $a_k$, determine the Fourier coefficients of the following sequences 
% in terms of $a_k$.
% 
% Note that
% 
% $$a_k = \frac1N \sum_{n=0}^{N-1} x[n] \mathrm{e}^{-\j2\mathrm{\pi}kn/N}. \qquad\qquad 
% (2.8.1)$$
% 
% 
% 
% *(c)* $x_3[n]=3\cos(2\mathrm{\pi}5n/N)x[-n],\; N>5:$ 
% 
% *Solution*: 
% 
% Using the analysis equation of the Fourier Series, we find that
% 
% $$x_3[n] = 3F\{\cos(2\pi5n/N)x[-n]\} = 3\left[\frac12e^{j\left(\frac{2\pi}{N}5n\right)} 
% + \frac12e^{-j\left(\frac{2\pi}{N}5n\right)}\right]x[-n] = \left[\frac32e^{j\left(\frac{2\pi}{N}5n\right)} 
% + \frac32e^{-j\left(\frac{2\pi}{N}5n\right)}\right]x[-n]$$
% 
% To obtain the Fourier coefficients of the sequence we can then use the following 
% properties of the Fourier Series:
% 
% $$x[-n] \leftrightarrow a_{-k}$$
% 
% $$e^{j\left(\frac{2\pi}{N}Mn\right)} \leftrightarrow a_{k-M}$$
% 
% $$e^{-j\left(\frac{2\pi}{N}Mn\right)} \leftrightarrow a_{k+M}$$
% 
% then let $c_k $ be the Fourier coefficients in terms of $a_k$ we get
% 
% $$c_k = \frac32\left(a_{-k-5} + a_{-k+5}\right)$$
% 
% 
% 
% *(d)* $x_4[n] = x[n]+x^\ast[-n]:$
% 
% *Solution*: 
% 
% Using the conjugate and folding properties of the Fourier series, we obtain 
% the coefficients represented by $c_k$ in terms of $a_k $ 
% 
% $$c_k = a_k + a_k^*$$
% 
% Since our sequence is an evenly symmetry sequence our final result of the 
% Fourier coefficients would be 
% 
% $$c_k = 2Re\{a_k\}$$
% 
% 
%% Problem 2.9
% *Text Problem 4.49, partsa (c) and (d) only, (Page 198)* 
% Determine sequences corresponding to each of the following Fourier transforms.
% 
% 
% 
% *(c)* $X_3(\mathrm{e}^{\j\omega}) = \j\mathrm{e}^{-\j4\omega}\bigl[2+3\cos(\omega)+\cos(2\omega)\bigr]:$
% 
% *Solution*: 
% 
% Distributing the $je^{-j4\omega}$ term to the components inside the brackets 
% we get: 
% 
% $$X_3(e^{jw}) = 2je^{-j4\omega} + 3je^{-j4\omega}\Big[\frac{e^{-j\omega}+e^{j\omega}}{2}\Big] 
% + je^{-j4\omega}\Big[\frac{e^{-j2\omega}+e^{j2\omega}}{2}\Big]$$
% 
% $$X_3(e^{jw}) = 2je^{-j4\omega} + \frac32je^{-j5\omega} + \frac32je^{-j3\omega} 
% + \frac12je^{-j6\omega} + \frac12je^{-j2\omega} $$
% 
% Applying the inverse Discrete Fourer Transform to $X_3(e^{jw})$ results in:
% 
% $$x_3[n] = \frac12j\delta(n - 2) + \frac32j\delta(n - 3) + 2j\delta(n - 4) 
% + \frac32j\delta(n - 5) + \frac12j\delta(n - 6)$$
% 
% 
% 
% *(d)* $X_4(\mathrm{e}^{\j\omega})=\left\{\begin{array}{l}2,\quad 0\le |\omega|\le 
% \mathrm{\pi}/8 \\1,\quad \mathrm{\pi}/8 \leq |\omega|\leq 3\mathrm{\pi}/4 \\0,\quad 
% 3\mathrm{\pi}/4\le |\omega|\le\mathrm{\pi}\end{array}\right.$
% 
% *Solution*: 
% 
% $$x_4[n] = \frac{1}{2\pi} \int_{<2\pi>}X(e^{j\omega})e^{j\omega n} d\omega$$
% 
% $$x_4[n] = \frac{1}{2\pi}\Big[\int_{\frac{-\pi}{8}}^0 2e^{j\omega n}d\omega 
% + \int_{0}^{\frac{\pi}{8}} 2e^{j\omega n}d\omega + \int_{-\frac{3\pi}{4}}^{\frac{\pi}{8}} 
% e^{j\omega n}d\omega + \int_{\frac{\pi}{8}}^{\frac{3\pi}{4}} e^{j\omega n}d\omega\Big]$$
% 
% $$x_4[n] = \frac{1}{2\pi}\Big[ 2 \Big(\frac{1 - e^{-j\frac{\pi}{8}n}}{jn}\Big) 
% + 2\Big(\frac{e^{j\frac{\pi}{8}n}-1}{jn}\Big) + \Big(\frac{e^{j\frac{\pi}{8}}-e^{-j\frac{3\pi}{4}n}}{jn}\Big) 
% + \Big(\frac{e^{j\frac{3\pi}{4}}-e^{j\frac{\pi}{8}n}}{jn}\Big)\Big]$$
% 
% $$x_4[n] = \frac{1}{2\pi}\left[\frac{2}{jn}(1-1) + \frac{2}{jn}\Big(e^{j\frac{\pi}{8}n} 
% - e^{-j\frac{\pi}{8}n}\Big) + \frac{1}{jn}\Big(e^{j\frac{\pi}{8}n} - e^{j\frac{\pi}{8}n}\Big) 
% + \frac{1}{jn}\Big(e^{j\frac{3\pi}{4}n} - e^{-j\frac{3\pi}{4}n}\Big) \right]$$
% 
% Using the euler's identity for $\sin$ and canceling some terms we obtain
% 
% $x_4[n] = \frac{1}{\pi n} \Big[2\sin\Big(\frac{\pi}{8}n\Big) + \sin\Big(\frac{3\pi}{4}n\Big)\Big]$ 
% is the corresponding sequence.
% 
% 
%% Problem 2.10
% *Text Problem 4.53 (Page 199)* 
% *Note*: There are two corrections in the errata sheet. Please follow them.
% 
% Let a sinusoidal pulse be given by $x(n) =\left( \cos \omega _0n\right) \bigl(u[n]-u[n-M]\bigr).$
% 
% 
% 
% *(a)* Using the frequency-shifting property of the DTFT, show that the real-part 
% of DTFT of $x(n)$ is given by
% 
% $$X_{\mathrm{R}}(\mathrm{e}^{\j\omega}) = \frac12 \cos\left[\frac{(\omega-\omega_0)(M-1)}{2}\right] 
% \left[\frac{\sin\left(\frac{(\omega-\omega_0)M}{2}\right)} {\sin\left(\frac{\omega-\omega_0}{2}\right)} 
% \right]+ \frac12 \cos\left[\frac{(\omega+\omega_0)(M-1)}{2}\right] \left[\frac{\sin\left(\frac{(\omega+\omega_0)M}{2}\right)} 
% {\sin\left(\frac{\omega+\omega_0}{2}\right)} \right] \qquad\qquad (2.10.1)$$
% 
% *Solution*: 
% 
% The frequency-shifting property of the DTFT says that $\rightarrow x[n]e^{j\omega_0n} 
% \leftrightarrow X\Big(e^{j(\omega - \omega_0)}\Big)$
% 
% Using euler's formula, we can represent the sinusoidal pulse $x(n)$ as 
% 
% $$\cos(\omega_0n) = \frac{e^{j\omega_0n} + e^{-j\omega_0n}}{2}$$ $$\rightarrow 
% x(n) = \frac12e^{j\omega_0n}(u[n]-u[n-M]) + \frac12e^{-j\omega_0n}(u[n]-u[n-M])$$
% 
% Here, the square pulse is being frequency-shifted by a factor of $\omega_0$ 
% so we can apply this when taking the DTFT of the pulse train, such that
% 
% $u[n] - u[n-M] \leftarrow\text{DTFT}\rightarrow \sum_{n=0}^{M-1}e^{-j\omega 
% n} = \frac{1-e^{-j\omega M}}{1-e^{-j\omega}}$ then apply the frequency shifts 
% to get
% 
% $$X(e^{j\omega}) = \frac12\left[\frac{1-e^{-jM(\omega-\omega_0) }}{1-e^{-j(\omega-\omega_0)}} 
% + \frac{1-e^{-jM(\omega+\omega_0) }}{1-e^{-j(\omega+\omega_0)}}\right]$$
% 
% We can reduce this term even further by factoring out exponential terms
% 
% $$\frac{1-e^{-jM(\omega-\omega_0) }}{1-e^{-j(\omega-\omega_0)}} \rightarrow 
% \frac{e^{-jM\Big(\frac{\omega - \omega_0}{2}\Big)}}{e^{-j\Big(\frac{\omega - 
% \omega_0}{2}\Big)}}\left[\frac{e^{jM\Big(\frac{\omega - \omega_0}{2}\Big)} - 
% e^{-jM\Big(\frac{\omega - \omega_0}{2}\Big)}}{e^{j\Big(\frac{\omega - \omega_0}{2}\Big)} 
% - e^{-j\Big(\frac{\omega - \omega_0}{2}\Big)}}\right] \rightarrow e^{-j\left(\frac{(\omega-\omega_0)(M-1)}{2}\right)}\left[\frac{2j\sin\left(\frac{(\omega-\omega_0)M}{2}\right)}{2j\sin\left(\frac{(\omega-\omega_0)}{2}\right)}\right]$$  
% 
% Thus, $X(e^{j\omega}) = \frac12e^{-j\left(\frac{(\omega-\omega_0)(M-1)}{2}\right)}\left[\frac{\sin\left(\frac{(\omega-\omega_0)M}{2}\right)}{\sin\left(\frac{(\omega-\omega_0)}{2}\right)}\right] 
% + \frac12e^{-j\left(\frac{(\omega+\omega_0)(M-1)}{2}\right)}\left[\frac{\sin\left(\frac{(\omega+\omega_0)M}{2}\right)}{\sin\left(\frac{(\omega+\omega_0)}{2}\right)}\right]$ 
% 
% Then if we only want to take the $Re\{X(e^{j\omega})\}$ we can expand out 
% the exponential terms to cosine terms resulting in
% 
% $$Re\{X(e^{j\omega})\} = X_R(e^{j\omega}) = \frac12\cos\left(\frac{(\omega-\omega_0)(M-1)}{2}\right)\left[\frac{\sin\left(\frac{(\omega-\omega_0)M}{2}\right)}{\sin\left(\frac{(\omega-\omega_0)}{2}\right)}\right] 
% + \frac12\cos\left(\frac{(\omega+\omega_0)(M-1)}{2}\right)\left[\frac{\sin\left(\frac{(\omega+\omega_0)M}{2}\right)}{\sin\left(\frac{(\omega+\omega_0)}{2}\right)}\right]$$ 
% 
% 
% 
% *(b)* Compute and plot $X_\mathrm{R}(\mathrm{e}^{\j\omega})$ for $\omega_0=\mathrm{\pi}/2$ 
% and $M=5,\;15,\;25,\;100.$ Use the plotting interval of $[-\mathrm{\pi},\mathrm{\pi}].$ 
% Comment on your results.
% 
% *MATLAB script*: 

clc; close all; clear;
Xr = zeros(4,629);
w0 = pi/2; omega = -pi:1/100:pi; M = 5;
A = sin(((omega-w0).*M)/2)./sin((omega-w0)./2); B = sin(((omega+w0).*M)./2)./sin((omega+w0)./2);
Xr(1,:) = 0.5.*((cos(((omega-w0).*(M-1))./2)).*A + (cos(((omega+w0).*(M-1))/2)).*B);
M = 15;
Xr(2,:) = 0.5.*((cos(((omega-w0).*(M-1))./2)).*A + (cos(((omega+w0).*(M-1))./2)).*B);
M = 25;
Xr(3,:) = 0.5.*((cos(((omega-w0).*(M-1))./2)).*A + (cos(((omega+w0).*(M-1))./2)).*B);
M = 100;
Xr(4,:) = 0.5.*((cos(((omega-w0).*(M-1))./2)).*A + (cos(((omega+w0).*(M-1))./2)).*B);
figure
subplot(4,1,1)
plot(omega,Xr(1,:))
xlim([-pi pi])
xticks([-pi -pi/2 0 pi/2 pi]);
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
ylabel("X_R(e^{j\omega})"), title("X_R(e^{j\omega}) for \it{M} = 5");
subplot(4,1,2)
plot(omega,Xr(2,:))
xlim([-pi pi])
xticks([-pi -pi/2 0 pi/2 pi]);
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
ylabel("X_R(e^{j\omega})"), title("X_R(e^{j\omega}) for \it{M} = 15");
subplot(4,1,3)
plot(omega,Xr(3,:))
xlim([-pi pi])
xticks([-pi -pi/2 0 pi/2 pi]);
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
ylabel("X_R(e^{j\omega})"), title("X_R(e^{j\omega}) for \it{M} = 25");
subplot(4,1,4)
plot(omega,Xr(4,:))
xlim([-pi pi])
xticks([-pi -pi/2 0 pi/2 pi]);
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
ylabel("X_R(e^{j\omega})"), title("X_R(e^{j\omega}) for \it{M} = 100");
%% 
% Comparing the results for $M = 5, 15, 25, 100$, we observe the real signals 
% are symmetrical across zero and are oscillating at higher frequencies as M increases. 
% Then since the cosine term is multiplying with a square-pulse train, we end 
% up with symmetrical, frequency-shifted sinc functions.
%% 
%