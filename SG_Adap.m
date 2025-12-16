
%%%% This code provides an efficient implementation of a Lagrange multiplier-based approach to solve systems of nonlinear wave equations with an adaptive time-stepping strategy.

    clear
    clc
    close
    delta = 0.5;
    x_span = [ - 10, 10 ]; y_span = x_span; t_span = [ 0, 10 ]; 
    hx = 1; hy = hx;
dt_min = 1 / 100; dt_max = 6e-2; alpha = 2;
    t_start = cputime;

    rho1 = 1 / delta ^ 2;
    rho2 = 1;

    x_axis = x_span( 1 ) : hx : x_span( 2 ) - hx;
    y_axis = y_span( 1 ) : hx : y_span( 2 ) - hy;

    beta1 = 1; beta2 = 1;

%% 
    hxy = hx * hy;

    Mx = length( x_axis );
    My = length( y_axis );
    N = length( t_span( 1 ) : dt_min : t_span( 2 ) );

    a = 320 / 393; b = 310 / 393;
    ex = ones( Mx, 1 ); ey = ones( My, 1 );
    Ix = eye( Mx ); Iy = eye( My ); Ixy = eye( Mx * My );
    
    Ax = spdiags( [ 23 / 2358 * ex, 344 / 1179 * ex, ex, 344 / 1179 * ex, 23 / 2358 * ex ], - 2 : 2, Mx, Mx );
    Ax( 1, end - 1 ) = 23 / 2358; Ax( 1, end ) = 344 / 1179;
    Ax( 2, end ) = 23 /2358;
    Ax( end - 1, 1 ) = 23 / 2358;
    Ax( end, 1 ) = 344 / 1179; Ax( end, 2 ) = 23 / 2358;
    invAx = Ax \ Ix;
    
    Ay = spdiags( [ 23 / 2358 * ey, 344 / 1179 * ey, ey, 344 / 1179 * ey, 23 / 2358 * ey ], - 2 : 2, My, My );
    Ay( 1, end - 1 ) = 23 / 2358; Ay( 1, end ) = 344 / 1179;
    Ay( 2, end ) = 23 /2358;
    Ay( end - 1, 1 ) = 23 / 2358;
    Ay( end, 1 ) = 344 / 1179; Ay( end, 2 ) = 23 / 2358;
    invAy = Ay \ Iy;

    DIx = spdiags( [ ex, - 2 * ex, ex ], - 1 : 1, Mx, Mx ) / hx ^ 2;
    DIx( 1, end ) = 1 / hx ^ 2; DIx( end, 1 ) = 1 / hx ^ 2;

    DIy = spdiags( [ ey, - 2 * ey, ey ], - 1 : 1, My, My ) / hy ^ 2;
    DIy( 1, end ) = 1 / hy ^ 2; DIy( end, 1 ) = 1 / hy ^ 2;

    DIx_h = spdiags( [ ex, 0 * ex, - 2 * ex, 0 * ex, ex ], - 2 : 2, Mx, Mx ) / ( 4 * hx ^ 2 );
    DIx_h( 1, end - 1 ) = 1 / ( 4 * hx ^ 2 );
    DIx_h( 2, end ) = 1 / ( 4 * hx ^ 2 );
    DIx_h( end - 1, 1 ) = 1 / ( 4 * hx ^ 2 );
    DIx_h( end, 2 ) = 1 / ( 4 * hx ^ 2 );

    DIy_h = spdiags( [ ey, 0 * ey, - 2 * ey, 0 * ey, ey ], - 2 : 2, My, My ) / ( 4 * hy ^ 2 );
    DIy_h( 1, end - 1 ) = 1 / ( 4 * hy ^ 2 );
    DIy_h( 2, end ) = 1 / ( 4 * hy ^ 2 );
    DIy_h( end - 1, 1 ) = 1 / ( 4 * hy ^ 2 );
    DIy_h( end, 2 ) = 1 / ( 4 * hy ^ 2 );

    coef = kron( Iy, invAx ) * kron( Iy, ( a * DIx + b * DIx_h ) ) + kron( invAy, Ix ) * kron( ( a * DIy + b * DIy_h ), Ix );
%%%%%%%%%%%%%%%%
    [ u, v, phi, psi ] = deal( zeros( Mx, My, N ) );
    [ eta, t_axis ] = deal( zeros( 1, N ) );
    eta( 1 ) = 1; t_axis( 1 ) = t_span( 1 );
    [ iter_count, dt, rho ] = deal( zeros( 1, N - 1 ) );
    dt( 1 ) = dt_min; rho( 1 ) = 0;

    M = zeros( Mx, My );
    for i = 1 : Mx
        for j = 1 : My
            M( i, j ) = x_axis( i ) ^ 2 + y_axis( j ) ^ 2;
        end
    end

    u( :, :, 1 ) = 4 * atan( exp( 3 - 5 * sqrt( M ) ) );
    v( :, :, 1 ) = 0.25 * u( :, :, 1 );
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
    beta1_coef = beta1 * coef;
    beta2_coef = beta2 * coef;
    rho1_beta1_coef = beta1_coef / rho1;
    rho2_beta2_coef = beta2_coef / rho2;

    for k = 1 : N - 1

        if k ~= 1
            dt( k ) = max( dt_min, dt_max * min( 1 / sqrt( 1 + alpha * ( err_u_k_kp1_2 / dt( k - 1 ) ) ^ 2 ), 1 / sqrt( 1 + alpha * ( err_v_k_kp1_2 / dt( k - 1 ) ) ^ 2 ) ) );
            rho( k ) = dt( k ) / dt( k - 1 );
        end
        t_axis( k + 1 ) = t_axis( k ) + dt( k );

        if t_axis( k + 1 ) >= t_span( 2 )
            t_axis( k + 1 : end ) = [];
            dt( k : end ) = []; rho( k : end ) = [];
            eta( k + 1 : end ) = [];
            break
        end

        rhok = rho( k );
        dtk = dt( k );

        dt2_4_rho_1 = 0.25 * dtk ^ 2 / rho1;
        dt2_4_rho_2 = 0.25 * dtk ^ 2 / rho2;
        dt2_4 = 0.25 * dtk ^ 2;
    
        uA1 = ( Ixy - dt2_4 * rho1_beta1_coef );
        uA2 = 2 * Ixy - uA1;
        uA1 = uA1 \ Ixy;
    
        vA1 = ( Ixy - dt2_4 * rho2_beta2_coef );
        vA2 = 2 * Ixy - vA1;
        vA1 = vA1 \ Ixy;
        
        uk = u( :, :, k ); phi_k = phi( :, :, k );
        vk = v( :, :, k ); psi_k = psi( :, :, k );

        ukr = reshape( uk, [], 1 ); phi_k_r = reshape( phi_k, [], 1 );
        vkr = reshape( vk, [], 1 ); psi_k_r = reshape( psi_k, [], 1 );

        if k == 1
            Fu_k_r = sin( ukr - vkr );
            Fv_k_r = - Fu_k_r;
        else
            Fu_k_r = reshape( sin( ( 1 + 0.5 * rhok ) * uk - 0.5 * rhok * u( :, :, k - 1 ) - ( 1 + 0.5 * rhok ) * vk + 0.5 * rhok * v( :, :, k - 1 ) ), [], 1 );
            Fv_k_r = - Fu_k_r;
        end

        qk = - uA1 * ( dt2_4_rho_1 * Fu_k_r );
        pk = uA1 * ( uA2 * ukr + dtk * phi_k_r ) + eta( k ) * qk;

        gk = - vA1 * ( dt2_4_rho_2 * Fv_k_r );
        fk = vA1 * ( vA2 * vkr + dtk * psi_k_r ) + eta( k ) * gk;

        Xk = 1; count = 0;
        tol_k = 1; F_k = 1;

        while ( abs( tol_k ) > 1e-15 ) && ( abs( F_k ) > 1e-15)

            u_kp1_r = pk + Xk * qk;
            v_kp1_r = fk + Xk * gk;

            R_k = hxy * sum( Fu_k_r .* ( u_kp1_r - ukr ) + Fv_k_r .* ( v_kp1_r - vkr ) );
            e_k = ( eta( k ) + Xk ) * 0.5;

            F_k = hxy * sum( cos( ukr - vkr ) - cos( u_kp1_r - v_kp1_r ) ) - e_k * R_k;
            J_k = hxy * sum( sin( u_kp1_r - v_kp1_r ) .* ( qk - gk ) ) - 0.5 * R_k - e_k * hxy * sum( Fu_k_r .* qk + Fv_k_r .* gk );

            tol_k = - J_k \ F_k;
            Xk = Xk + tol_k;

            count = count + 1;

        end

        eta( k + 1 ) = Xk;
        u_kp1_r = pk + Xk * qk;
        v_kp1_r = fk + Xk * gk;

        u( :, :, k + 1 ) = reshape( u_kp1_r, Mx, My );
        v( :, :, k + 1 ) = reshape( v_kp1_r, Mx, My );

        phi( :, :, k + 1 ) = 2 * ( u( :, :, k + 1 ) - uk ) / dtk - phi_k;
        psi( :, :, k + 1 ) = 2 * ( v( :, :, k + 1 ) - vk ) / dtk - psi_k;

        iter_count( k ) = count;

        diff_u_k_kp1 = u( :, :, k + 1 ) - uk;
        diff_v_k_kp1 = v( :, :, k + 1 ) - vk;

        err_u_k_kp1_2 = sqrt( hxy * sum( sum( diff_u_k_kp1 .^ 2 ) ) );
        err_v_k_kp1_2 = sqrt( hxy * sum( sum( diff_v_k_kp1 .^ 2 ) ) );

    end

%%%%%%%%%%%%%%%%%%%%%
    max_iter = max( iter_count );

    eta_err = abs( eta - 1 );
    max_eta_err = max( eta_err );

    avg_dt = ( dt_max + dt_min ) / 2;

    [ X, Y ] = meshgrid( x_axis, x_axis );
    for kk = [ 0 1, 3, 5, 7, 9]

        [ ~, idx ] = min( abs( t_axis - kk ) );

        figure
        [ ~, h1 ] = contourf( X', Y', sin(u( :, :, idx )/2), 50 );
        hold on
        set( h1, 'linestyle', 'none' )
        title( sprintf( 'T = %f', t_axis( idx ) ) )
        hold off
      xlabel('x')
    ylabel('y')
% % %         figure
% % %        [ ~, h2 ] = contourf( X', Y', v( :, :, idx ), 50 );
% % % %  mesh( X', Y', v( :, :, idx ) );
% % %         hold on
% % %          set( h2, 'linestyle', 'none' )
% % %         title( sprintf( 'T = %f', t_axis( idx ) ) )
% % %         hold off
% % % xlabel('x')
% % %     ylabel('y')
    end

    t_end = cputime;
    t_elapsed = t_end - t_start;