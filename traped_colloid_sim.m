
clear all ; close all ; clc ;

%% Load constants       ==================================================================
constants

%% Initialise variables ==================================================================
particle = zeros ( frames , 1 ) ;                                           % Particle positions in 1 D
k = 1e-6 ;                                                                  % Trap stiffness, assuming Hook
tic                                                                         % runtime

%% Start simulation ======================================================================
sigma = sqrt( 2 * k_B * T * dt / mobility ) ;                               % [ m kg / s ] 

noise  = sigma * randn( frames - 1 , 1 ) ;                                  % correlated excitation

for frame_index  = 1 : frames - 1                                           % Calculate Langevin

    V = - k * particle( frame_index ) ;                                     % Potential acting on particle
    dparticle = mobility * V * dt + mobility * noise( frame_index , 1 ) ;
    particle( frame_index + 1 ) = particle( frame_index ) + dparticle ;
end

t = dt * ones( frames - 1 , 1) ; t( 1 ) = 0 ; t = cumsum( t ) ;             % Build timestamps

particle( 1 , : ) = [ ] ;                                                   % Trim starting zeros

trap_stiff = k_B * T / mean( particle.^2 ) ;

%% Calculate power spectrum ===========================================================
% Calculate power spectrum for entire data set, samp_f and smooth using blocking
fNyq    =   frame_index / max( t ) / 2 ;                                    % Nyquist frequency
x_fft = fft( particle ) ;                                                   % Fourier transform data
x_0 = fftshift( x_fft ) ;                                                   % Shift values to centre on zero freq
x_ps = x_0 .* conj( x_0 ) ;                                                 % Power spectra
freq = ( [ 1 : frame_index ]  / max( t ) ) ;                                % 0-centered frequency range

ind  = find( freq <= fNyq ) ;                                               % only to the Nyquist f
freq = rot90( freq( ind ) ) ;                                               % rotate to keap everthing n x 1
x_ps = x_ps( ind ) ;

% Blocking data to smooth out noise
count = 0 ;
x_ps_smooth = zeros( floor( length( x_ps ) / block ) , 1 ) ;
freq_smooth = zeros( floor( length( x_ps ) / block ) , 1 ) ;

for b = 1 : block : length( x_ps ) - block
    count = count + 1 ;
    x_block = sum( x_ps( b : b + block ) ) / block ;
    freq_block = sum( freq( b : b + block ) ) / block ;
    x_ps_smooth( count ) = x_block ;
    freq_smooth(count)=freq_block;
end
    
% Scalling power spectra according to comp. phys. comm. 159 225 (2004), note that var(x) = area under spectra for this scalling
PS = mean( dt )^2 * x_ps_smooth / max( t ) ;

toc
	

%% Plotting ==========================================================================
range = 1 : 1e4 ;
ptsize = 2 ;
linewdth = 0 ;

scrsz = get( 0 , 'ScreenSize' ) ;
scr_wdth = scrsz( 3 ) ;
scr_hght = scrsz( 4 ) ;

%Particle position-------------------------------------------------------------------------------------------------------------------------
fig1 = figure( 'OuterPosition' , [ 0 0 scr_wdth / 2 scr_hght / 2 ] ) ;

subplot ( 2 , 2 , 1 ) ;
plot1 = plot( t( range ) , particle( range ) , 'k' ) ;

subplot ( 2 , 2 , 2 ) ;
[ dist , dist_x ] = hist( particle , sqrt( frames - 1 ) ) ;
U_x = - log( dist / sum( dist) ) ;
plot( dist_x , U_x, 'o',...
	'MarkerEdgeColor','g',...
	'MarkerFaceColor','g',...
	'MarkerSize', ptsize ) ;

subplot ( 2 , 2 , 3 ) ;
loglog( freq_smooth(  freq_smooth <= freq_trim  , 1 ) , PS(  freq_smooth <= freq_trim  , 1 ) , '.r' ) ;
hold on;
xlabel('{Frequency (Hz)}');
ylabel('{P(f) (m s)}');

y_top = 1e-14; n_tex = 8; y_tex = 0 : n_tex ;                               % Position of top text box; number of text boxes; Array of y positions of text boxes
y_pos = y_top * exp( -y_tex ) ;                                             % Array of text box positions linearly spaed on a log scale
x_pos = 2 ;

n_guess = viscosity ;                                                       % Initial guess for viscosity in mPa
PS_par_guess = [ a , viscosity, k ] ;                                       % [ particle radius, viscosity, trap stiffness]
fit_opt = statset('MaxIter',1e6);
PS_par_fit = nlinfit( freq_smooth(  freq_smooth <= freq_trim  , 1 ) , PS(  freq_smooth <= freq_trim  , 1 ) ,  @Lorentzian , PS_par_guess , fit_opt );

PS_fit = Lorentzian( PS_par_guess , freq_smooth(  freq_smooth <= freq_trim  , 1 ) ) ;

loglog( freq_smooth(  freq_smooth <= freq_trim  , 1 ) , PS_fit , '-b' ) ;

hold off