clear;clc;
folder_list = dir('./');
cur_dir = pwd();
num_of_dir = length(folder_list);
for k = 1:num_of_dir
    if (folder_list(k).isdir)
        num_of_dirname = length(folder_list(k).name);
        if num_of_dirname < 5
            continue;
        end
        if all(folder_list(k).name(1:4)=='1938')
            Flag = work_in_single_folder(folder_list(k).name, cur_dir);
        end
    end
end

function Flag = work_in_single_folder(folder_name, root_dir)
chdir(folder_name);
delete('记录.txt');
diary('记录.txt');
diary on;
namelist = dir('./');
num_of_file = length(namelist);
for k = 1:num_of_file
    num_of_name = length(namelist(k).name);
    if num_of_name < 10
        continue;
    end
    if all(namelist(k).name(1:5)=='RW200') && all(namelist(k).name(num_of_name-3:num_of_name)=='.TXT')
        for t=6:num_of_name-3
            if all(namelist(k).name(t:t+1)=='常值')
                cal_step(namelist(k).name);
            end
            if all(namelist(k).name(t:t+1)=='正弦')
                cal_sine(namelist(k).name);
            end
            if all(namelist(k).name(t:t+1)=='传函')
                cal_tf(namelist(k).name);
            end
        end
    end
end
diary off;
chdir(root_dir);
Flag = strcat(folder_name, ' Success!');
end

function [t_up,t_adjust,t_peak,overshoot,max_err,mean_err,var_err,sigma_3] = cal_step(name_of_datafile)
    disp(strcat('processing:',name_of_datafile,'...'));
    [Counter,Cmd,Ans,~,~] = textread(name_of_datafile,'%f%f%f%f%f');
    timestep = 0.1;
    time = timestep * (Counter - Counter(1));
    Num = length(Ans);
    for k = 2:Num-20
        if abs(Ans(k)-Ans(k-1)) > 10
            for t = 1:20
                if Ans(k+t) - Ans(k-1) < 1
                    break;
                end
            end
            Ans(k:k+t-1) = (Ans(k+t)-Ans(k-1))/(time(k+t)-time(k-1))*(time(k:k+t-1)-time(k-1))+Ans(k-1);
        end
    end
    final = mean(Ans(Num-249:Num-49));
    peak = max(Ans);
    t_peak = time(Ans == peak);
    if length(t_peak) > 1
        t_peak = mean(t_peak);
    end
    overshoot = (peak - final) / final;
    for k = Num-1:-1:1
        if Ans(k)>1.05*final || Ans(k)<0.95*final
            t_adjust = time(k+1);
            break;
        end
    end
    if k == 1
        t_adjust = inf;
    end
    for k=1:Num
        if Ans(k)>0.1*final 
            left_bound = time(k);
            break;
        end
    end
    for k=k+1:Num
        if Ans(k)>0.9*final 
            right_bound = time(k);
            break;
        end
    end
    t_up = right_bound - left_bound;
    err = Ans - Cmd;
    max_err = max(abs(err(Num-249:Num-49)));
    mean_err = mean(err(Num-249:Num-49));
    var_err = var(err(Num-249:Num-49));
    sigma_3 = 3 * sqrt(var_err);
    avg_tor = cal_torq(time,Ans);
    disp(name_of_datafile);
    disp(strcat('上升时间 = ',num2str(t_up),' 秒'));
    disp(strcat('调节时间 = ',num2str(t_adjust),' 秒'));
    disp(strcat('峰值时间 = ',num2str(t_peak),' 秒'));
    disp(strcat('超调量 = ',num2str(overshoot*100),'%'));
    disp(strcat('最大误差 = ',num2str(max_err),' rps'));
    disp(strcat('误差均值 = ',num2str(mean_err),' rps'));
    disp(strcat('误差方差 = ',num2str(var_err),' rps'));
    disp(strcat('3σ = ',num2str(sigma_3),' rps'));
    if length(avg_tor)==3
        disp(strcat('0rps～50rps 平均力矩 = ',num2str(avg_tor(1)),' Nm'));
        disp(strcat('50rps～100rps 平均力矩 = ',num2str(avg_tor(2)),' Nm'));
        disp(strcat('100rps～125rps 平均力矩 = ',num2str(avg_tor(3)),' Nm'));
    elseif length(avg_tor)==2
        disp(strcat('0rps～50rps 平均力矩 = ',num2str(avg_tor(1)),' Nm'));
        disp(strcat('50rps～100rps 平均力矩 = ',num2str(avg_tor(2)),' Nm'));
    else
        disp(strcat('0rps～50rps 平均力矩 = ',num2str(avg_tor(1)),' Nm'));
    end
    disp(' ');
    epsilon = sqrt(1 / ((-pi / log(overshoot))^2 + 1));
    omega_n = 3.5 / epsilon / t_adjust;
    s = tf('s');
    s1 = -epsilon * omega_n + omega_n * sqrt(epsilon^2 - 1);
    s2 = -epsilon * omega_n - omega_n * sqrt(epsilon^2 - 1);
    func = 1 / (s - s1) / (s - s2)
    figure(1)
    subplot(2,1,1);
    hold on;grid on;
    plot(time,Cmd,'-.blue');
    plot(time,Ans,'-red');
    legend('指令','响应','Location','SouthEast');
    set(get(gca, 'XLabel'), 'String', '时间(s)');
    set(get(gca, 'YLabel'), 'String', '转速(rps)');
    subplot(2,1,2);
    hold on;grid on;
    plot(time,err,'-red');
    set(get(gca, 'XLabel'), 'String', '时间(s)');
    set(get(gca, 'YLabel'), 'String', '转速(rps)');
    k = length(name_of_datafile);
    name_of_datafile(k-2:k) = 'png';
    saveas(gcf,name_of_datafile);
    % figure(2)
    % subplot(2,1,1);
    % hold on;grid on;
    % plot(time,Cur,'-red');
    % set(get(gca, 'XLabel'), 'String', '时间(s)');
    % set(get(gca, 'YLabel'), 'String', '电流(mA)');
    % subplot(2,1,2);
    % hold on;grid on;
    % plot(time,Tem,'-red');
    % set(get(gca, 'XLabel'), 'String', '时间(s)');
    % set(get(gca, 'YLabel'), 'String', '温度({}^{\circ}C)');
    % saveas(gcf,strcat(name_of_datafile,'电流与温度.png'));
    close all;
end
function [max_err,mean_err,var_err,sigma_3] = cal_sine(name_of_datafile)
    disp(strcat('processing:',name_of_datafile,'...'));
    [Counter,Cmd,Ans,~,~] = textread(name_of_datafile,'%f%f%f%f%f');
    Num = length(Ans);
    timestep = 0.1;
    time = timestep * (Counter - Counter(1));
    for k = 2:Num-20
        if abs(Ans(k)-Ans(k-1)) > 10
            for t = 1:20
                if Ans(k+t) - Ans(k-1) < 1
                    break;
                end
            end
            Ans(k:k+t-1) = (Ans(k+t)-Ans(k-1))/(time(k+t)-time(k-1))*(time(k:k+t-1)-time(k-1))+Ans(k-1);
        end
    end
    err = Ans - Cmd;
    max_err = max(abs(err(Num-249:Num-49)));
    mean_err = mean(err(Num-249:Num-49));
    var_err = var(err(Num-249:Num-49));
    sigma_3 = 3 * sqrt(var_err);
    disp(strcat('最大误差 = ',num2str(max_err),' rps'));
    disp(strcat('误差均值 = ',num2str(mean_err),' rps'));
    disp(strcat('误差方差 = ',num2str(var_err),' rps'));
    disp(strcat('3σ = ',num2str(sigma_3),' rps'));
    disp(' ');
    figure(1)
    subplot(2,1,1);
    hold on;grid on;
    plot(time,Cmd,'-.blue');
    plot(time,Ans,'-red');
    legend('指令','响应');
    set(get(gca, 'XLabel'), 'String', '时间(s)');
    set(get(gca, 'YLabel'), 'String', '转速(rps)');
    subplot(2,1,2);
    hold on;grid on;
    plot(time,err,'-red');
    set(get(gca, 'XLabel'), 'String', '时间(s)');
    set(get(gca, 'YLabel'), 'String', '转速(rps)');
    k = length(name_of_datafile);
    name_of_datafile(k-2:k) = 'png';
    saveas(gcf,name_of_datafile);
    % figure(2)
    % subplot(2,1,1);
    % hold on;grid on;
    % plot(time,Cur,'-red');
    % set(get(gca, 'XLabel'), 'String', '时间(s)');
    % set(get(gca, 'YLabel'), 'String', '电流(mA)');
    % subplot(2,1,2);
    % hold on;grid on;
    % plot(time,Tem,'-red');
    % set(get(gca, 'XLabel'), 'String', '时间(s)');
    % set(get(gca, 'YLabel'), 'String', '温度({}^{\circ}C)');
    % saveas(gcf,strcat(name_of_datafile,'电流与温度.png'));
    close all;
end

function avg_tor = cal_torq(time,Ans)
    len = length(time);
    t1=0;t2=0;t3=0;
    for k=1:len
        if Ans(k) > 50
            t1 = k;
            break;
        end
    end
    for k=k+1:len
        if Ans(k) > 100
            t2 = k;
            break;
        end
    end
    for k=k+1:len
        if Ans(k) > 125
            t3 = k;
            break;
        end
    end
    torq = diff(Ans) / 0.1 * 2 * pi;
    if t2==0
        avg_tor = mean(torq(1:t1)) * 6.095e-6;
    elseif t3==0
        avg_tor = [mean(torq(1:t1)),mean(torq(t1+1:t2))] * 6.095e-6;
    else
        avg_tor = [mean(torq(1:t1)),mean(torq(t1+1:t2)),mean(torq(t2+1:t3))] * 6.095e-6;
    end
end

function [] = cal_tf(name_of_datafile)
    disp(strcat('processing:',name_of_datafile,'...'));
    [Time, Cmd, Ans, ~, ~]=textread(name_of_datafile,'%f%f%f%f%f');
    timestep = 0.1;
    t = timestep * (Time-Time(1));
    err = Ans - Cmd;
    figure(1);
    subplot(2,1,1);
    hold on;grid on;
    plot(t,Cmd,'-blue','linewidth',1.5);
    plot(t,Ans,'-.red' ,'linewidth',1.5);
    xlabel('时间(s)');
    ylabel('转/秒');
    legend('指令','响应');
    subplot(2,1,2);
    hold on;grid on;
    plot(t,err,'-red','linewidth',1);
    xlabel('时间(s)');
    ylabel('转/秒');
    saveas(gcf,strcat('响应与指令.png'));
    time_stamp = t(diff([Cmd(1); Cmd])~=0);
    for ind = 1:length(time_stamp)-1
        if ind == 1
            left = time_stamp(ind)/2;
        else
            left = (time_stamp(ind - 1) + time_stamp(ind)) / 2;
        end
        right = (time_stamp(ind + 1) + time_stamp(ind)) / 2;
        IOcharactor(50, left, right, t, Cmd, Ans)
    end
end

function func = IOcharactor(picnum ,left_bound, right_bound, t, Cmd, Ans)
    %IOcharactor - Description
    %
    % Syntax: output = IOcharactor(left_bound, right_bound, t, Cmd, Ans)
    %
    % Long description
        head = find(t == left_bound);
        tail = find(t == right_bound);
        tim = t(head:tail);
        cmd = Cmd(head:tail);
        res = Ans(head:tail);
        figure(picnum);
        hold on; grid on;
        plot(tim, cmd, '-blue', 'linewidth', 1.5);
        plot(tim, res, '-.red', 'linewidth', 1.5);
        xlabel('时间(s)');
        ylabel('转/秒');
        legend('指令', '响应');
        stage = cmd(end) - cmd(1);
        overshoot = (max(res) - cmd(end)) / stage;
        shot = find(diff([cmd(1); cmd]) == max(diff([cmd(1); cmd])));
        findadjtime = true;
        findrisetime = true;
    
        for po = tail - head + 1:-1:1
    
            if abs(res(po) - cmd(end)) > stage * 0.05 && findadjtime
                findadjtime = false;
                adjtime = tim(po) - tim(shot);
            end
    
            if abs(res(po) - cmd(end)) > stage * 0.2 && findrisetime
                findrisetime = false;
                risetime = tim(po) - tim(shot);
            end
    
        end
    
        epsilon = sqrt(1 / ((-pi / log(overshoot))^2 + 1));
        omega_n = 3.5 / epsilon / adjtime;
        s = tf('s');
        s1 = -epsilon * omega_n + omega_n * sqrt(epsilon^2 - 1);
        s2 = -epsilon * omega_n - omega_n * sqrt(epsilon^2 - 1);
        disp(strcat('时间区间:', num2str(left_bound), "s~", num2str(right_bound), 's'));
        disp(strcat('阶跃高度:', num2str(stage * 2 * pi), "rad/s"));
        disp(strcat('超调量:', num2str(overshoot * 100), '%'));
        disp(strcat('80%上升时间:', num2str(risetime), 's'));
        disp(strcat('5%调节时间:', num2str(adjtime), 's'));
        disp('传递函数:');
        func = 1 / (s - s1) / (s - s2);
    end
