function rmask = twoSpeaker(sig, sid, type, nGau, bW, snr_criterion, nStep, workFolder)

% parameters
nChan = 128;

% temporay file for store mixture
prefix = [workFolder,'/',datestr(now,30),num2str(round(rand()*100))];
while (size(dir([prefix,'*']))>0)
    pause(1);
    prefix = [workFolder,'/',datestr(now,30),num2str(round(rand()*100))];
end
sigFN = [prefix,'.in'];
dlmwrite(sigFN,sig);

% separation
switch type
    case 'ReddyRaj07'
        ver = 0;
        cmd = sprintf('!../c/twoSpk/twoSpk %f %d %.1f %d %d %d %s.in %s.out',...
            0, nChan, ver, sid(1), sid(2), nGau, prefix, prefix);
        eval(cmd);
        rmask{1} = load([prefix,'.out'])';
        rmask{2} = 1-rmask{1};
    
    case {'MMSE','MAP','acoustDym'}
        switch type
            case 'MMSE'
                ver = 1;
            case 'MAP'
                ver = 2;
            case 'acoustDym'
                ver = 4;
        end
        switch ver
            case {1,2}
                cmd = sprintf('!../c/twoSpk/twoSpk %f %d %.1f %d %d %d %s.in %s.out',...
                    0, nChan, ver, sid(1), sid(2), nGau, prefix, prefix);                
            case 4
                cmd = sprintf('!../c/twoSpk/twoSpk %f %d %.1f %d %d %d %s.in %d %s.out',...
                    0, nChan, ver, sid(1), sid(2), nGau, prefix, bW, prefix);                
        end
        eval(cmd);
        rmask{1} = load([prefix,'.out'])';
        rmask{2} = 1-rmask{1};       
        
    case {'MAP_iter','MMSE_iter'}        
        addpath('snr'); % note this includes other libraries        
        addpath('cochleagram'); % note this includes other libraries        
        switch type
            case 'MMSE_iter'
                ver = 1;
            case 'MAP_iter'
                ver = 2;
        end        
        % initial mask estimation
        step = 1;
        cmd = sprintf('!../c/twoSpk/twoSpk %f %d %.1f %d %d %d %s.in %s.out',...
             0, nChan, ver, sid(1), sid(2), nGau, prefix, prefix);
        eval(cmd);
        % store results
        est_isnr{step} = 0;
        tmasks{step} = load([prefix,'.out'])';        
        % iterative estimation
        snr_diff = inf;
        while step<nStep && abs(snr_diff)>snr_criterion            
            step = step +1;            
            % input SNR estimation
            params.tmpDir = workFolder;
            params.broadband = '';
            params.ebm = tmasks{step-1};
            [dummy, tmp] = snrEst(sig, params);
            est_isnr{step} = tmp.broadband;
            snr_diff = est_isnr{step} - est_isnr{step-1};
            for i=1:step
                fprintf('* Step %d, estimated input SNR=%f dB *\n', i, est_isnr{i});
            end
            fprintf('\n');            
            % mask estimation
            cmd = sprintf('!../c/twoSpk/twoSpk %f %d %.1f %d %d %d %s.in %s.out',...
                est_isnr{step}, nChan, ver, sid(1), sid(2), nGau, prefix, prefix);
            eval(cmd);
            tmasks{step} = load([prefix,'.out'])';
        end
        rmask{1} = tmasks{end};
        rmask{2} = 1-rmask{1};
        
    case {'acoustDym_iter'}
        ver = 4;                
        step = 1;               
        % initial mask estimation
        cmd = sprintf('!../c/twoSpk/twoSpk %f %d %.1f %d %d %d %s.in %d %s.out',...
            0, nChan, ver, sid(1), sid(2), nGau, prefix, bW, prefix);
        eval(cmd);
        % store results
        est_isnr{step} = 0;
        tmasks{step} = load([prefix,'.out'])';                    
        % iterative estimation
        snr_diff = inf;
        while step<nStep && abs(snr_diff)>snr_criterion
            step = step +1;            
            % input SNR estimation
            params.tmpDir = workFolder;
            params.broadband = '';
            params.ebm = tmasks{step-1};
            [dummy, tmp] = snrEst(sig, params);
            est_isnr{step} = tmp.broadband;
            snr_diff = est_isnr{step} - est_isnr{step-1};            
            for i=1:step
                fprintf('* Step %d, estimated input SNR=%f dB *\n', i, est_isnr{i});
            end
            fprintf('\n');
            % mask estimation            
            cmd = sprintf('!../c/twoSpk/twoSpk %f %d %.1f %d %d %d %s.in %d %s.out',...
                est_isnr{step}, nChan, ver, sid(1), sid(2), nGau, prefix, bW, prefix);            
            eval(cmd);
            tmasks{step} = load([prefix,'.out'])';
        end
      
        rmask{1} = tmasks{end};
        rmask{2} = 1-rmask{1};
end
dlmwrite('ratio_mask_1.ascii', rmask{1});
dlmwrite('ratio_mask_2.ascii', rmask{2});
delete([prefix,'.*']);
