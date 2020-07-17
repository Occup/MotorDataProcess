close all;
clear;
clc;
delete('记录.txt');
diary('记录.txt');
diary on;
[Time, Cmd, Ans, Cur, Tem]=textread('RW200_2019-9-30_10-29-29_-速度-频率.TXT','%f%f%f%f%f');
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

% figure(2);
% subplot(2,1,1);
% hold on;grid on;
% plot(t,Cur,'-red','linewidth',1.5);
% xlabel('时间(s)');
% ylabel('电流(mA)');
% subplot(2,1,2);
% hold on;grid on;
% plot(t,Tem,'-red' ,'linewidth',1.5);
% xlabel('时间(s)');
% ylabel('温度({}^{\circ}C)');
% saveas(gcf,strcat('电流与温度.png'));

IOcharactor(2, 30, 50, t, Cmd, Ans)
IOcharactor(3, 50, 70, t, Cmd, Ans)
IOcharactor(4, 70, 100, t, Cmd, Ans)
diary off;
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
    saveas(gcf,strcat(num2str(left_bound),'s~',num2str(right_bound),'s.png'));
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