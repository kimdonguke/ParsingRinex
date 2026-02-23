%%%% SV Clock Correction %%%%
function [ttx_epoch, ttx_table] = clock_Correction(Obs,Nav)  % ttx_cell{i,1} = time, ttx_cell{i,2}=satellite number
    %% get Obs data %%
    %sObs=table2array(Obs);
    c=299792458;                    % speed of light(meter/sec)
    tolerance = 1e-12;              %convergence threshold
    mu=3.986005e+14;                % earth's gravitational constant meter^3/sec^2
    F=-2*sqrt(mu)/c^2;              % 상대성 효과 보정 힘
    L1 = 1575.42;                   % L1CA 주파수
    L2 = 1227.60;                   % L2 주파수
    mu2 = (L1/L2)^2;                % TGD 계수
    
    Obs = groupingTable(Obs, 'Time');
    
    ttx_list = [];
    ttx_epoch=cell(height(Obs),1);
    
    
    for i=1:height(Obs) % iteration by whole grouped data
        
        obs_i=Obs{i};
        trx = obs_i.Time(1);
        PRN=obs_i.PRN;
        t_k=0;
        index_num=0;
        
        TTX=nan(length(PRN),1);
        P=Obs{i}.P1;     %Pseudorange
        SNR = Obs{i}.SNR1; % 
        % code = Obs{i}.code;
        
        %%% algorithm problem on here
        for j=1:length(PRN) % iteration by each satellite, calculate t_tx
             sNav=Nav{PRN(j)};
             index_num=get_NavIndexbyTime(sNav,trx);
    
             M0=sNav.M0(index_num);
             A=sNav.A(index_num);
             toe=sNav.toe(index_num);
             toc=sNav.toc(index_num);
             deln=sNav.deln(index_num);
             e=sNav.e(index_num); 
             af0=sNav.af0(index_num);
             af1=sNav.af1(index_num);
             af2=sNav.af2(index_num);
             TGD=sNav.tgd(index_num);
             % if obs_i.code == 23
             %     TGD = mu2*TGD;
             % end
             
             ttx_new=obs_i.Time(j)-obs_i.P1(j)/c;   % define ttx_new to first iteration set ttx_new
             t_SV=ttx_new;
             ttx=0;
             % fprintf("ttx_new = ")
             % disp(fff(ttx_new))
             max_iter=1;
            
             while abs(ttx_new-ttx)>tolerance && max_iter<200
                max_iter=max_iter+1;

                ttx=ttx_new;
                t_k=ttx_new-toe; %  time from reference epoch
                M_k=M0+(sqrt(mu/(A^3))+deln)*t_k; % M_k = M_0 + n*t_k /// n=n_0+del n
                E_k=M_k;
                for k=1:10
                    E_kOld=E_k;
                    E_k=M_k+e*sin(E_k);
                    if norm(E_k-E_kOld)<tolerance
                        break
                    end
                end
                delt_r=F*sqrt(A)*e*sin(E_k);
                delt_SV=af0+af1*(ttx-toc)+af2*(ttx-toc)^2+delt_r-TGD;
                ttx_new=t_SV-delt_SV;
             end
             P(j)=P(j)+c*delt_SV;
             TTX(j)=ttx_new;
        end

        ttx_epoch{i}=table(PRN,TTX,P,SNR);
        ttx_list=[ttx_list; [PRN,TTX,P,SNR]];
    end
    PRN = ttx_list(:,1);
    TTX = ttx_list(:,2);
    P = ttx_list(:,3);
    SNR = ttx_list(:,4);
    ttx_table = table(PRN,TTX,P,SNR);
end