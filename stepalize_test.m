function stepalize_test

    log = @(x) fprintf(x);

    tried = 0;
    failed = 0;
    
    [t, f] = test_arg_parsing(log);
    tried = tried + t; failed = failed + f;
    
    [t, f] = test_unconstrained(log);
    tried = tried + t; failed = failed + f;
    
    [t, f] = test_eig_constrained(log);
    tried = tried + t; failed = failed + f;
    
    [t, f] = test_response_constrained_delay(log);
    tried = tried + t; failed = failed + f;

    [t, f] =  test_response_constrained_no_delay(log);
    tried = tried + t; failed = failed + f;
    
    log(['\n\nTotal Tests: ' num2str(tried) '\nTotal failures: ' num2str(failed) '\n\n'])
end


function [tried, failed] = test_unconstrained(log)
    % Test a few basic systems, both single-output and multi-output.
    tried = 0;
    failed = 0;
    
    log(['\nTesting unconstrained realization\n', ...
           '---------------------------------------------------------']);
    
    function assertAlmostEqualEig(A0, A1)
        tried = tried + 1;
        maxDif = max(abs(sort(eig(A0)) - sort(eig(A1))));
        if maxDif > 1e-12
            log(['fail\n  - Max difference in eigenvalues: ' num2str(maxDif) '.\n']);
            failed = failed + 1;
        else
            log('pass\n');
        end
    end

    function assertAlmostEqualResponse(y0, y1)
        tried = tried + 1;
        maxDif = max(max(abs(y0 - y1)));
        if maxDif > 1e-12
            log(['fail\n  - Max difference in response: ' num2str(maxDif) '.\n']);
            failed = failed + 1;
        else
            log('pass\n');
        end
    end

    function assertEstimateSuccess(G0, y0, Ge)
        log('...eigenvalues...');
        assertAlmostEqualEig(Ge.a, G0.a);
    
        ye = step(Ge);
        log('...response...');
        assertAlmostEqualResponse(y0, ye);
    end

    
    % Stable single-output system.
    A0 = [0.3, 0.1; -0.1, 0.3];
    B0 = [1; 0];
    C0 = [1, 0];
    D0 = 0.2;
    
    log('\nTest unconstrained SISO system w/ 1 time delay...\n');
    G0 = ss(A0, B0, C0, 0, []);
    y0 = step(G0);
    [Ae, Be, Ce, De] = stepalize(y0, 'Order', 2);
    Ge = ss(Ae, Be, Ce, De, []);
    assertEstimateSuccess(G0, y0, Ge);
    
    log('\nTest unconstrained SISO system w/ 0 time delay...\n');
    G0 = ss(A0, B0, C0, D0, []);
    y0 = step(G0);
    [Ae, Be, Ce, De] = stepalize(y0, 'Order', 2, 'TimeDelay', 0);
    Ge = ss(Ae, Be, Ce, De, []);
    assertEstimateSuccess(G0, y0, Ge);
    
    log('\nTest unconstrained SISO system w/ 2 time delays...\n');
    G0 = ss(A0, B0, C0, 0, []);
    G0.InputDelay = 1;
    y0 = step(G0);
    [Ae, Be, Ce, De] = stepalize(y0, 'Order', 2, 'TimeDelay', 2);
    Ge = ss(Ae, Be, Ce, De, []);
    Ge.InputDelay = 1;
    assertEstimateSuccess(G0, y0, Ge);
    
    % Stable multi-output system.
    A0 = [0.3, 0.1; -0.1, 0.3];
    B0 = [1; 0];
    C0 = [1, 0; 0, 1];
    D0 = [0.2; 0.1];
    
    log('\nTest unconstrained SIMO system w/ 1 time delay\n');
    G0 = ss(A0, B0, C0, 0, []);
    y0 = step(G0);
    [Ae, Be, Ce, De] = stepalize(y0, 'Order', 2);
    Ge = ss(Ae, Be, Ce, De, []);
    assertEstimateSuccess(G0, y0, Ge);
    
    log('\nTest unconstrained SIMO system w/ 0 time delay\n');
    G0 = ss(A0, B0, C0, D0, []);
    y0 = step(G0);
    [Ae, Be, Ce, De] = stepalize(y0, 'Order', 2, 'TimeDelay', 0);
    Ge = ss(Ae, Be, Ce, De, []);
    assertEstimateSuccess(G0, y0, Ge);

    log('\nTest unconstrained SIMO system w/ 2 time delays...\n');
    G0 = ss(A0, B0, C0, 0, []);
    G0.InputDelay = 1;
    y0 = step(G0);
    [Ae, Be, Ce, De] = stepalize(y0, 'Order', 2, 'TimeDelay', 2);
    Ge = ss(Ae, Be, Ce, De, []);
    Ge.InputDelay = 1;
    assertEstimateSuccess(G0, y0, Ge);
    
    log(['\n' num2str(failed) ' out of ' num2str(tried) ' tests failed.\n']);
end


function [tried, failed] = test_eig_constrained(log)
    tried = 0;
    failed = 0;
    
    log(['\nTesting eigenvalue-constrained realization\n', ...
           '---------------------------------------------------------']);
       
    function assertAlmostEqualEig(A0, A1)
        tried = tried + 1;
        maxDif = max(abs(sort(eig(A0)) - sort(eig(A1))));
        if maxDif > 1e-10
            log(['fail\n  - Max difference in eigenvalues: ' num2str(maxDif) '.\n']);
            failed = failed + 1;
        else
            log('pass\n');
        end
    end
    
    function assertAlmostEqualResponse(y0, y1)
        tried = tried + 1;
        maxDif = max(max(abs(y0 - y1)));
        if maxDif > 1e-10
            log(['fail\n  - Max difference in response: ' num2str(maxDif) '.\n']);
            failed = failed + 1;
        else
            log('pass\n');
        end
    end

    function assertEstimateEqual(G0, y0, Ge)
        log('...eigenvalues...');
        assertAlmostEqualEig(Ge.a, G0.a);
    
        ye = step(Ge);
        log('...response...');
        assertAlmostEqualResponse(y0, ye);
    end

    function assertEstimateInSet(Ge, f)
        tried = tried + 1;
        log('...constraints...');
        lambda_e = eig(Ge);
        
        result = false;
        for k = 1:length(lambda_e)
            if any(eig(f(lambda_e(k))) < 0)
                result = true;
                break;
            end
        end
        
        if result
            log('fail\n');
        else
            log('pass\n');
        end
    end

    
    % LMI region for strictly positive eigenvalues. 
    % Default delta_p is 1e-4.
    delta_p = 1e-4;
    alpha_p = delta_p*[2, 0; 0, -2];
    beta_p = [0, 0; 0, 1];
    f_p = @(z) alpha_p + beta_p*z + beta_p'*z';
    

    log('\nTest positive real eigenvalue constraint - system in set...\n');
    A0 = [0.3, 0.1; -0.1, 0.3];
    B0 = [1; 0];
    C0 = [1, 0; 0, 1];
    D0 = [0; 0];
    G0 = ss(A0, B0, C0, D0, []);
    y0 = step(G0);
    [Ae, Be, Ce, De] = stepalize(y0, 'Order', 2, 'Eigenvalues', 'Positive');
    Ge = ss(Ae, Be, Ce, De, []);
    assertEstimateEqual(G0, y0, Ge);
    assertEstimateInSet(Ge, f_p);
    
    
    log('\nTest positive real eigenvalue constraint - system not in set...\n');
    A0 = [-0.3, 0.1; -0.1, -0.3];
    B0 = [1; 0];
    C0 = [1, 0; 0, 1];
    D0 = [0; 0];
    G0 = ss(A0, B0, C0, D0, []);
    y0 = step(G0);
    [Ae, Be, Ce, De] = stepalize(y0, 'Order', 2, 'Eigenvalues', 'Positive');
    Ge = ss(Ae, Be, Ce, De, []);
    assertEstimateInSet(Ge, f_p);
    
    % LMI region for positive eigenvalues w/ custom delta. 
    delta_p = 0.25;
    alpha_p = delta_p*[2, 0; 0, -2];
    beta_p = [0, 0; 0, 1];
    f_p = @(z) alpha_p + beta_p*z + beta_p'*z';
    

    log('\nTest positive real eigenvalue constraint & custom delta - system in set...\n');
    A0 = [0.3, 0.1; -0.1, 0.3];
    B0 = [1; 0];
    C0 = [1, 0; 0, 1];
    D0 = [0; 0];
    G0 = ss(A0, B0, C0, D0, []);
    y0 = step(G0);
    [Ae, Be, Ce, De] = stepalize(y0, 'Order', 2, 'Eigenvalues', 'Positive', 'DeltaP', delta_p);
    Ge = ss(Ae, Be, Ce, De, []);
    assertEstimateEqual(G0, y0, Ge);
    assertEstimateInSet(Ge, f_p);
    
    
    log('\nTest positive real eigenvalue constraint & custom delta - system not in set...\n');
    A0 = [0.2, 0.1; -0.1, 0.2];
    B0 = [1; 0];
    C0 = [1, 0; 0, 1];
    D0 = [0; 0];
    G0 = ss(A0, B0, C0, D0, []);
    y0 = step(G0);
    [Ae, Be, Ce, De] = stepalize(y0, 'Order', 2, 'Eigenvalues', 'Positive', 'DeltaP', delta_p);
    Ge = ss(Ae, Be, Ce, De, []);
    assertEstimateInSet(Ge, f_p);
    
    
    % LMI region for strictly stable eigenvalues. 
    % Default delta_s is 1e-4.
    delta_s = 1e-4;
    alpha_s = (1 - delta_s)*eye(2);
    beta_s = [0, 1; 0, 0];
    f_s = @(z) alpha_s + beta_s*z + beta_s'*z';
    

    log('\nTest stable eigenvalue constraint - system in set...\n');
    A0 = [0.3, 0.1; -0.1, 0.3];
    B0 = [1; 0];
    C0 = [1, 0; 0, 1];
    D0 = [0; 0];
    G0 = ss(A0, B0, C0, D0, []);
    y0 = step(G0);
    [Ae, Be, Ce, De] = stepalize(y0, 'Order', 2, 'Eigenvalues', 'Stable');
    Ge = ss(Ae, Be, Ce, De, []);
    assertEstimateEqual(G0, y0, Ge);
    assertEstimateInSet(Ge, f_s);
    
    
    log('\nTest stable eigenvalue constraint - system not in set...\n');
    A0 = [1, 0.1; -0.1, 1];
    B0 = [1; 0];
    C0 = [1, 0; 0, 1];
    D0 = [0; 0];
    G0 = ss(A0, B0, C0, D0, []);
    y0 = step(G0, 70);
    [Ae, Be, Ce, De] = stepalize(y0, 'Order', 2, 'Eigenvalues', 'Stable');
    Ge = ss(Ae, Be, Ce, De, []);
    assertEstimateInSet(Ge, f_s);

    
    % LMI region for strictly stable eigenvalues w/ custom delta. 
    delta_s = 0.4;
    alpha_s = (1 - delta_s)*eye(2);
    beta_s = [0, 1; 0, 0];
    f_s = @(z) alpha_s + beta_s*z + beta_s'*z';
    

    log('\nTest stable eigenvalue constraint & custom delta - system in set...\n');
    A0 = [0.3, 0.1; -0.1, 0.3];
    B0 = [1; 0];
    C0 = [1, 0; 0, 1];
    D0 = [0; 0];
    G0 = ss(A0, B0, C0, D0, []);
    y0 = step(G0);
    [Ae, Be, Ce, De] = stepalize(y0, 'Order', 2, 'Eigenvalues', 'Stable', 'DeltaS', delta_s);
    Ge = ss(Ae, Be, Ce, De, []);
    assertEstimateEqual(G0, y0, Ge);
    assertEstimateInSet(Ge, f_s);
    
    
    log('\nTest stable eigenvalue constraint & custom delta - system not in set...\n');
    A0 = [0.5, 0.1; -0.1, 0.5];
    B0 = [1; 0];
    C0 = [1, 0; 0, 1];
    D0 = [0; 0];
    G0 = ss(A0, B0, C0, D0, []);
    y0 = step(G0, 70);
    [Ae, Be, Ce, De] = stepalize(y0, 'Order', 2, 'Eigenvalues', 'Stable', 'DeltaS', delta_s);
    Ge = ss(Ae, Be, Ce, De, []);
    assertEstimateInSet(Ge, f_s);
    

    % LMI region for real eigenvalues. 
    % Default delta_r is 1e-8.
    delta_r = 1e-8;
    alpha_r = delta_r*eye(2);
    beta_r = [0, 0.5; -0.5, 0];
    f_r = @(z) alpha_r + beta_r*z + beta_r'*z';
    

    log('\nTest real eigenvalue constraint - system in set...\n');
    A0 = [-0.3, 0; 0, 0.8];
    B0 = [1; 1];
    C0 = [1, 1; 0, 1];
    D0 = [0; 0];
    G0 = ss(A0, B0, C0, D0, []);
    y0 = step(G0);
    [Ae, Be, Ce, De] = stepalize(y0, 'Order', 2, 'Eigenvalues', 'Real');
    Ge = ss(Ae, Be, Ce, De, []);
    assertEstimateEqual(G0, y0, Ge);
    assertEstimateInSet(Ge, f_r);
    
    
    log('\nTest real eigenvalue constraint - system not in set...\n');
    A0 = [0.3, 0.1; -0.1, 0.3];
    B0 = [1; 0];
    C0 = [1, 0; 0, 1];
    D0 = [0; 0];
    G0 = ss(A0, B0, C0, D0, []);
    y0 = step(G0, 70);
    [Ae, Be, Ce, De] = stepalize(y0, 'Order', 2, 'Eigenvalues', 'Real');
    Ge = ss(Ae, Be, Ce, De, []);
    assertEstimateInSet(Ge, f_s);
    
    
    % LMI region for real eigenvalues w/ custom delta.
    delta_r = 0.2;
    alpha_r = delta_r*eye(2);
    beta_r = [0, 0.5; -0.5, 0];
    f_r = @(z) alpha_r + beta_r*z + beta_r'*z';
    

    log('\nTest real eigenvalue constraint - system in set...\n');
    A0 = [0.3, 0.1; -0.1, 0.3];
    B0 = [1; 1];
    C0 = [1, 1; 0, 1];
    D0 = [0; 0];
    G0 = ss(A0, B0, C0, D0, []);
    y0 = step(G0);
    [Ae, Be, Ce, De] = stepalize(y0, 'Order', 2, 'Eigenvalues', 'Real', 'DeltaR', delta_r);
    Ge = ss(Ae, Be, Ce, De, []);
    assertEstimateEqual(G0, y0, Ge);
    assertEstimateInSet(Ge, f_r);
    
    
    log('\nTest real eigenvalue constraint - system not in set...\n');
    A0 = [0.3, 0.4; -0.4, 0.3];
    B0 = [1; 0];
    C0 = [1, 0; 0, 1];
    D0 = [0; 0];
    G0 = ss(A0, B0, C0, D0, []);
    y0 = step(G0);
    [Ae, Be, Ce, De] = stepalize(y0, 'Order', 2, 'Eigenvalues', 'Real', 'DeltaR', delta_r);
    Ge = ss(Ae, Be, Ce, De, []);
    assertEstimateInSet(Ge, f_s);
    
    
    % All three:
    f = @(z) [f_r(z), zeros(2, 4); zeros(2), f_p(z), zeros(2); zeros(2, 4), f_s(z)];
    log('\nTest real, positive, and stable eigenvalue constraint - system in set...\n');
    A0 = [0.3, 0; 0, 0.8];
    B0 = [1; 1];
    C0 = [1, 1; 0, 1];
    D0 = [0; 0];
    G0 = ss(A0, B0, C0, D0, []);
    y0 = step(G0);
    [Ae, Be, Ce, De] = stepalize(y0, 'Order', 2, 'Eigenvalues', {'Real', 'Positive', 'Stable'});
    Ge = ss(Ae, Be, Ce, De, []);
    assertEstimateEqual(G0, y0, Ge);
    assertEstimateInSet(Ge, f);
    
    log('\nTest real, positive, and stable eigenvalue constraint - system not in set...\n');
    A0 = [-1, 0.5; -0.5, -1];
    B0 = [1; 1];
    C0 = [1, 1; 0, 1];
    D0 = [0; 0];
    G0 = ss(A0, B0, C0, D0, []);
    y0 = step(G0);
    [Ae, Be, Ce, De] = stepalize(y0, 'Order', 2, 'Eigenvalues', {'Real', 'Positive', 'Stable'});
    Ge = ss(Ae, Be, Ce, De, []);
    assertEstimateInSet(Ge, f);
    
    log(['\n' num2str(failed) ' out of ' num2str(tried) ' tests failed.\n']);
end


function [tried, failed] = test_response_constrained_delay(log)
    % Test a few basic systems, both single-output and multi-output.
    tried = 0;
    failed = 0;
    
    log(['\nTesting steady-state-constrained realization with time delay\n', ...
           '------------------------------------------------------------']);
    
    function assertAlmostEqualSteadyState(yss0, Ge)
        tried = tried + 1;
        ysse = (Ge.c/(eye(size(Ge.a)) - Ge.a)*Ge.b + Ge.d)';
        maxDif = max(max(abs(yss0 - ysse)));
        if maxDif > 1e-8
            log(['fail\n  - Max difference in steady-state: ' num2str(maxDif) '.\n']);
            failed = failed + 1;
        else
            log('pass\n');
        end
    end
    
    function assertAlmostEqualResponse(y0, y1)
        tried = tried + 1;
        maxDif = max(max(abs(y0 - y1)));
        if maxDif > 1e-8
            log(['fail\n  - Max difference in response: ' num2str(maxDif) '.\n']);
            failed = failed + 1;
        else
            log('pass\n');
        end
    end

    function assertEstimateEqual(y0, yss, Ge)
        log('...steady-state...');
        assertAlmostEqualSteadyState(yss, Ge);
    
        ye = step(Ge);
        log('...response...');
        assertAlmostEqualResponse(y0, ye);
    end

    function assertNoUndershoot(Ge)
        log('...undershoot...');
        tried = tried + 1;
        ye = step(Ge);
        yss = Ge.c/(eye(order(Ge)) - Ge.a)*Ge.b + Ge.d;
        if any(any(ye.*repmat(sign(yss'), size(ye, 1), 1) < 0))
            log('fail\n');
            failed = failed + 1;
        else
            log('passed\n');
        end
    end

    log('\nTest steady-state constraint - use true value...\n');
    A0 = [0.3, 0.4; -0.4, 0.3];
    B0 = [1; 0];
    C0 = [1, 0; 0, 1];
    D0 = [0; 0];
    G0 = ss(A0, B0, C0, D0, []);
    y0 = step(G0);
    yss = (C0/(eye(2) - A0)*B0 + D0)';
    [Ae, Be, Ce, De] = stepalize(y0, 'Order', 2, 'SteadyState', yss);
    Ge = ss(Ae, Be, Ce, De, []);
    assertEstimateEqual(y0, yss, Ge);
    
    log('\nTest steady-state constraint - use wrong value...\n');
    yss = [3, -1];
    [Ae, Be, Ce, De] = stepalize(y0, 'Order', 2, 'SteadyState', yss);
    Ge = ss(Ae, Be, Ce, De, []);
    log('...steady-state...');
    assertAlmostEqualSteadyState(yss, Ge);
    
    log(['\nTesting undershoot-constrained realization\n', ...
           '---------------------------------------------------------']);
       
    % Test when the true system has no undershoot.
    log('\nTest undershoot constraint - true has no undershoot...\n');
    [Ae, Be, Ce, De] = stepalize(y0, 'Order', 2, 'NoUndershoot', true);
    Ge = ss(Ae, Be, Ce, De, []);
    yss = (C0/(eye(2) - A0)*B0 + D0)';
    assertEstimateEqual(y0, yss, Ge);
    
    % This has some serious undershoot. Make sure we can estimate it
    % unconstrained first.
    log('\nTest undershoot constraint - true has undershoot, no constraint...\n');
    G0u = ss(zpk(1.2, [0.3+0.4i, 0.3-0.4i], 1, []));
    G0u = ss(G0u.a, G0u.b, [G0u.c; 1, 0], [0; 0], []);
    y0u = step(G0u);
    yssu = (G0u.c/(eye(2) - G0u.a)*G0u.b + G0u.d)';
    [Ae, Be, Ce, De] = stepalize(y0u, 'Order', 2, 'NoUndershoot', false);
    Ge = ss(Ae, Be, Ce, De, []);
    assertEstimateEqual(y0u, yssu, Ge);
    
    % Now see if the constraints work.
    log('\nTest undershoot constraint - true has undershoot with constraint...\n');
    [Ae, Be, Ce, De] = stepalize(y0u, 'Order', 2, 'NoUndershoot', true);
    Ge = ss(Ae, Be, Ce, De, []);
    assertNoUndershoot(Ge);
    
    % It turns out this problem is infeasible, so make sure we can figure
    % that out.
    log('\nTest infeasible constraint...');
    tried = tried + 1;
    try
        stepalize(y0u, 'Order', 2, 'NoUndershoot', true, 'SteadyState', yssu);
    catch ex
        if strcmpi(ex.identifier, 'STEPALIZE:YALMIPError')
            log('pass');
        else
            log('failed');
            failed = failed + 1;
        end
    end
end


function [tried, failed] = test_response_constrained_no_delay(log)
    % Test a few basic systems, both single-output and multi-output.
    tried = 0;
    failed = 0;
    
    log(['\nTesting steady-state-constrained realization with no time delay\n', ...
           '---------------------------------------------------------------']);
    
    function assertAlmostEqualSteadyState(yss0, Ge)
        tried = tried + 1;
        ysse = (Ge.c/(eye(size(Ge.a)) - Ge.a)*Ge.b + Ge.d)';
        maxDif = max(max(abs(yss0 - ysse)));
        if maxDif > 1e-8
            log(['fail\n  - Max difference in steady-state: ' num2str(maxDif) '.\n']);
            failed = failed + 1;
        else
            log('pass\n');
        end
    end
    
    function assertAlmostEqualResponse(y0, y1)
        tried = tried + 1;
        maxDif = max(max(abs(y0 - y1)));
        if maxDif > 1e-8
            log(['fail\n  - Max difference in response: ' num2str(maxDif) '.\n']);
            failed = failed + 1;
        else
            log('pass\n');
        end
    end

    function assertEstimateEqual(y0, yss, Ge)
        log('...steady-state...');
        assertAlmostEqualSteadyState(yss, Ge);
    
        ye = step(Ge);
        log('...response...');
        assertAlmostEqualResponse(y0, ye);
    end

    function assertNoUndershoot(Ge)
        log('...undershoot...');
        tried = tried + 1;
        ye = step(Ge);
        yss = Ge.c/(eye(order(Ge)) - Ge.a)*Ge.b + Ge.d;
        if any(any(ye.*repmat(sign(yss'), size(ye, 1), 1) < 0))
            log('fail\n');
            failed = failed + 1;
        else
            log('passed\n');
        end
    end

    log('\nTest steady-state constraint - use true value...\n');
    A0 = [0.3, 0.4; -0.4, 0.3];
    B0 = [1; 0];
    C0 = [1, 0; 0, 1];
    D0 = [0; -0.3];
    G0 = ss(A0, B0, C0, D0, []);
    y0 = step(G0);
    yss = (C0/(eye(2) - A0)*B0 + D0)';
    [Ae, Be, Ce, De] = stepalize(y0, 'Order', 2, 'TimeDelay', 0, 'SteadyState', yss);
    Ge = ss(Ae, Be, Ce, De, []);
    assertEstimateEqual(y0, yss, Ge);

    
    log('\nTest steady-state constraint - use wrong value...\n');
    yss = [3, -1];
    [Ae, Be, Ce, De] = stepalize(y0, 'Order', 2, 'SteadyState', yss);
    Ge = ss(Ae, Be, Ce, De, []);
    log('...steady-state...');
    assertAlmostEqualSteadyState(yss, Ge);
    
    
    log(['\nTesting undershoot-constrained realization\n', ...
           '---------------------------------------------------------']);
       
    % Test when the true system has no undershoot.
    log('\nTest undershoot constraint - true has no undershoot...\n');
    [Ae, Be, Ce, De] = stepalize(y0, 'Order', 2, 'TimeDelay', 0, 'NoUndershoot', true);
    Ge = ss(Ae, Be, Ce, De, []);
    yss = (C0/(eye(2) - A0)*B0 + D0)';
    assertEstimateEqual(y0, yss, Ge);

    
    % This has some serious undershoot. Make sure we can estimate it
    % unconstrained first.
    log('\nTest undershoot constraint - true has undershoot, no constraint...\n');
    G0u = ss(zpk(1.2, [0.3+0.4i, 0.3-0.4i], 1, []));
    G0u = ss(G0u.a, G0u.b, [G0u.c; 1, 0], [0.1; 0.2], []);
    y0u = step(G0u);
    yssu = (G0u.c/(eye(2) - G0u.a)*G0u.b + G0u.d)';
    [Ae, Be, Ce, De] = stepalize(y0u, 'Order', 2, 'TimeDelay', 0, 'NoUndershoot', false);
    Ge = ss(Ae, Be, Ce, De, []);
    assertEstimateEqual(y0u, yssu, Ge);
    
    % Now see if the constraints work.
    log('\nTest undershoot constraint - true has undershoot with constraint...\n');
    [Ae, Be, Ce, De] = stepalize(y0u, 'Order', 2, 'TimeDelay', 0, 'NoUndershoot', true);
    Ge = ss(Ae, Be, Ce, De, []);
    assertNoUndershoot(Ge);
    
    
    % It turns out this problem is infeasible, so make sure we can figure
    % that out.
    log('\nTest undershoot constraint with steady-state...\n');
    tried = tried + 1;
	[Ae, Be, Ce, De] = stepalize(y0u, 'Order', 2, 'TimeDelay', 0, 'NoUndershoot', true, 'SteadyState', yssu);
    Ge = ss(Ae, Be, Ce, De, []);
	assertNoUndershoot(Ge);
end


function [tried, failed] = test_arg_parsing(log)
% TODO: finish
    tried = 0;
    failed = 0;
    
    % This asserts that the expected exception is thrown.
    function assertException(fnc, id)
        tried = tried + 1;
        try
            fnc();
        catch ex
            if strcmp(ex.identifier, id)
                log('pass\n');
            else
                log(['fail\n  - Expected: ' id '\n  - Raised:   ' ...
                        ex.identifier '\n']);
                failed = failed + 1;
            end
        end
    end

    % This checks to make sure non-numeric types are not allowed for a
    % given argument.
    function assertNonNumericException(fnc)
        tried = tried + 1;
        
        % Test for cell types.
        log('......cell...');
        x = {};
        assertException(@() fnc(x), 'MATLAB:invalidType');
        
        log('......struct...');
        x = struct();
        assertException(@() fnc(x), 'MATLAB:invalidType');
        
        log('......logical...');
        x = true;
        assertException(@() fnc(x), 'MATLAB:invalidType');
        
        log('......char...');
        x = 'blah';
        assertException(@() fnc(x), 'MATLAB:invalidType');
    end
        

    log('\nTest output signal "y" argument:\n');
    log('...non-numeric argument...\n');
    assertNonNumericException(@(x) stepalize(x));
    
    log('...empty argument...');
    y = [];
    assertException(@() stepalize(y), 'MATLAB:expectedNonempty');
       
    log('...complex argument...');
    y = rand(1) + rand(1)*1i;
    assertException(@() stepalize(y), 'MATLAB:expectedReal');
    
    log('...NaN argument...');
    y = NaN;
    assertException(@() stepalize(y), 'MATLAB:expectedNonNaN');

    log('...Inf argument...');
    y = Inf;
    assertException(@() stepalize(y), 'MATLAB:expectedFinite');
    
    
    % We have multiple optional arguments that are matrix dimensions. This
    % one function tests them all.
    function assertArgIsDimension(arg)
        log(['\nTest "' arg '" argument:\n']);
        
        log('...non-numeric argument...\n');
        assertNonNumericException(@(x) stepalize(0, arg, x));
        
        log('...zero argument...');
        x = 0;
        assertException(@() stepalize(0, arg, x), 'MATLAB:expectedPositive');
        
        log('...negative argument...');
        x = -1;
        assertException(@() stepalize(0, arg, x), 'MATLAB:expectedPositive');
        
        log('...non-itegeter argument...');
        x = 1.1;
        assertException(@() stepalize(0, arg, x), 'MATLAB:expectedInteger');
        
        log('...vector argument...');
        x = [1, 1];
        assertException(@() stepalize(0, arg, x), 'MATLAB:expectedScalar');
        
        log('...matrix argument...');
        x = ones(2);
        assertException(@() stepalize(0, arg, x), 'MATLAB:expectedScalar');
    end

    assertArgIsDimension('Order');
    assertArgIsDimension('BlockRows');
    
    log('\nTest SteadyState argument...\n');
    assertNonNumericException(@(x) stepalize(0, 'SteadyState', x));
    
    log('...dimensions...');
    y = zeros(3, 2);
    yss = zeros(3, 1);
    assertException(@() stepalize(y, 'SteadyState', yss), 'STEPALIZE:incorrectSteadyStateSize');
    
    log('...empty argument...');
    yss = [];
    assertException(@() stepalize(y, 'SteadyState', yss), 'MATLAB:expectedVector');
       
    log('...complex argument...');
    yss = rand(1) + rand(1)*1i;
    assertException(@() stepalize(y, 'SteadyState', yss), 'MATLAB:expectedReal');
    
    log('...NaN argument...');
    yss = NaN;
    assertException(@() stepalize(y, 'SteadyState', yss), 'MATLAB:expectedNonNaN');

    log('...Inf argument...');
    yss = Inf;
    assertException(@() stepalize(y, 'SteadyState', yss), 'MATLAB:expectedFinite');
    
    
    log(['\n' num2str(failed) ' out of ' num2str(tried) ' tests failed.\n']);
end

