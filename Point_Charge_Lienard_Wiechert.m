function Point_Charge_Lienard_Wiechert()
close all
format long
% Initialization and variable declarations
fig = figure('Position', [300 100 400 400], 'WindowButtonDownFcn', @myBDCallback, 'WindowButtonUpFcn', @myBUCallback);
fig.KeyPressFcn = @simulatior_break;
fig.Color = 'w';
q = 1;
k = 8.9875e9;
xq = 1;
yq = 1;
h = 0.01;
[x, y] = meshgrid(0:h:2, 0:h:2);
state = 0;

c = 1; % speed of light

position_str = zeros(5000,3,1);
v_str = zeros(5000,3,1);

position_str(1,1:2,1:length(xq)) = [xq; yq];


% Calculation of electric field
EX = zeros(size(x));
EY = zeros(size(x));
for i = 1:length(q)
    rx = x - xq(i);
    ry = y - yq(i);
    r = sqrt(rx.^2 + ry.^2);
    EX = EX + k * q(i) * rx ./ r.^3;
    EY = EY + k * q(i) * ry ./ r.^3;
end
E = sqrt(EX.^2 + EY.^2);

% Video Writer Settings
%video = VideoWriter('Point_Charge.avi', 'Uncompressed AVI');
%open(video)

disp("Press escape key to exit.")

% Drawing a point charge
ss{1} = scatter(2.1, 1.5, 200, 'red', 'filled');
hold on
ss{2} = scatter(2.1, 0.5, 200, 'blue', 'filled');
rectangle('Position', [0.1 0.1 1.9 1.8], 'EdgeColor', 'k')
for k = 1:length(q)
    s{k} = scatter(xq(k), yq(k), 200, 'filled');
    if q(k) > 0
        s{k}.MarkerFaceColor = [1 0 0];
    else
        s{k}.MarkerFaceColor = [0 0 1];
    end
end

pulse = '+';
minus = '-';
for k = 1:length(q)
    if q(k) > 0
        T{k} = text(xq(k) - 0.038, yq(k), pulse, 'FontSize', 20, 'Color', 'w');
    else
        T{k} = text(xq(k) - 0.038, yq(k), minus, 'FontSize', 20, 'Color', 'w');
    end
end
text(ss{1}.XData - 0.038, ss{1}.YData, pulse, 'FontSize', 20, 'Color', 'w');
text(ss{2}.XData - 0.038, ss{2}.YData, minus, 'FontSize', 20, 'Color', 'w');

% Calculate electric lines of force
n = 15;
dl = h * 10;
theta = linspace(0, 2 * pi, n).';
N = 200;
X = zeros(n, N, length(q));
Y = zeros(n, N, length(q));
for k = 1:length(q)
    X(:, 1, k) = xq(k) .* ones(size(theta)) + dl * cos(theta);
    Y(:, 1, k) = yq(k) .* ones(size(theta)) + dl * sin(theta);
end
switch_state = ones(n, k);
draw_state = ones(n,k);
i = 2;
while any(switch_state == 1, 'all')
    EX_tmp = ones(n, 1);
    EY_tmp = ones(n, 1);
    E_tmp = ones(n, 1);
    for k = 1:length(q)
        switch_state(X(:, i - 1, k) > 2.2, k) = 0;
        switch_state(Y(:, i - 1, k) > 2.2, k) = 0;
        switch_state(X(:, i - 1, k) < -0.2, k) = 0;
        switch_state(Y(:, i - 1, k) < -0.2, k) = 0;
        for ii = 1:n
            if switch_state(ii, k) == 0
                continue
            end
            if any(sqrt((xq - X(ii, i - 1, k)).^2 + (yq - Y(ii, i - 1, k)).^2) < dl * 0.99)
                switch_state(ii, k) = 0;
                draw_state(ii,k) = 0;
                continue
            end            
            x_idx = round(X(ii, i - 1, k) / h);
            y_idx = round(Y(ii, i - 1, k) / h);

            if x_idx <= 0 || y_idx <= 0 || x_idx >= length(x(:, 1)) || y_idx >= length(y(:, 1))
                switch_state(ii, k) = -1;
                continue
            end
            EX_tmp(ii) = EX(y_idx,x_idx);
            EY_tmp(ii) = EY(y_idx,x_idx);
            E_tmp(ii) = E(y_idx,x_idx);
        end
        X(:, i, k) = X(:, i - 1, k) + sign(q(k)) * (EX_tmp ./ E_tmp) .* switch_state(:, k) * dl * 0.1;
        Y(:, i, k) = Y(:, i - 1, k) + sign(q(k)) * (EY_tmp ./ E_tmp) .* switch_state(:, k) * dl * 0.1;
    end
    i = i + 1;
end

% Drawing lines of electric force
for ii = 1:n
    for k = 1:length(q)
        handle{ii, k} = plot(X(ii, 1:i - 1, k), Y(ii, 1:i - 1, k), 'k-', 'LineWidth', 0.8, 'LineStyle', '-.');
        if q(k) < 0 && draw_state(ii,k) == 1 || q(k) > 0
            handle{ii, k}.Visible = 'on';
        else
            handle{ii, k}.Visible = 'off';
        end
    end
end

percent = 0;
percent_msg = "$||v||/c$: " + sprintf('%.4f',percent);
percent_text = text(0.6,0,percent_msg,'Interpreter','latex','FontSize',15);


% Plot Settings
axis([0.1 2.3 0.1 1.9])
daspect([1 1 1])
axis off
set(gca, 'LooseInset', get(gca, 'TightInset'))

% Main loop
step = 1;
fin = 0;
while fin ~= 1
    step = step + 1;
    EX = zeros(size(x));
    EY = zeros(size(x));
    for i = 1:length(q)
        if i == state    
            % Electric field with Liénard–Wiechert potentials
            Cp = get(gca,'CurrentPoint');
            position_str(step,:,state) = [Cp(1,1) Cp(1,2) 0];
            v = (position_str(step,:,state) - position_str(step-1,:,state));
            v_mag = norm(v);
            if v_mag > c
                v = v/v_mag;
            end

            v_str(step,:,state) = v;
            a = (v_str(step,:,state) - v_str(step-1,:,state));

            
            Ex = zeros(size(x));
            Ey = zeros(size(x));
            v_mag = norm(v);
            percent = v_mag/c;
            percent_text.String = "$||v||/c$: " + sprintf('%.4f',percent);
            for ii = 1:numel(x)
                R_vec = [x(ii) - position_str(step,1,state), y(ii) - position_str(step,2,state), 0];
                R_mag = norm(R_vec);

                term1 = (1 - (v_mag/c)^2) * (R_vec - dot(R_vec, v)/c * v) / (R_mag - dot(R_vec, v)/c)^3;
                term2 = cross(R_vec, cross(R_vec - dot(R_vec, v)/c, a)) / (R_mag - dot(R_vec, v)/c)^3;
    
                E_LW = q(state) * (term1 + term2/c);
    
                Ex(ii) = E_LW(1);
                Ey(ii) = E_LW(2);
            end
            EX = EX + Ex;
            EY = EY + Ey;
        else
            rx = x - xq(i);
            ry = y - yq(i);
            r = sqrt(rx.^2 + ry.^2);
            EX = EX + k * q(i) * rx ./ r.^3;
            EY = EY + k * q(i) * ry ./ r.^3;
        end        
    end
   
    E = sqrt(EX.^2 + EY.^2);

    X = zeros(n, N, length(q));
    Y = zeros(n, N, length(q));
    for k = 1:length(q)
        X(:, 1, k) = xq(k) .* ones(size(theta)) + dl * cos(theta);
        Y(:, 1, k) = yq(k) .* ones(size(theta)) + dl * sin(theta);
    end
    switch_state = ones(n, k);
    draw_state = ones(n, k);
    i = 2;
    while any(switch_state == 1, 'all')   
        EX_tmp = ones(n, 1);
        EY_tmp = ones(n, 1);
        E_tmp = ones(n, 1);
        for k = 1:length(q)
            switch_state(X(:, i - 1, k) > 2.2, k) = 0;
            switch_state(Y(:, i - 1, k) > 2.2, k) = 0;
            switch_state(X(:, i - 1, k) < -0.2, k) = 0;
            switch_state(Y(:, i - 1, k) < -0.2, k) = 0;
            for ii = 1:n
                if switch_state(ii, k) == 0
                    continue
                end
                if any(sqrt((xq - X(ii, i - 1, k)).^2 + (yq - Y(ii, i - 1, k)).^2) < dl * 0.99)
                    switch_state(ii, k) = 0;
                    draw_state(ii,k) = 0;
                    continue
                end                
                x_idx = round(X(ii, i - 1, k) / h);
                y_idx = round(Y(ii, i - 1, k) / h);

                if x_idx <= 0 || y_idx <= 0 || x_idx >= length(x(:, 1)) || y_idx >= length(y(:, 1))
                    switch_state(ii, k) = 0;
                    continue
                end
                EX_tmp(ii) = EX(y_idx,x_idx);
                EY_tmp(ii) = EY(y_idx,x_idx);
                E_tmp(ii) = E(y_idx,x_idx);
            end

            X(:, i, k) = X(:, i - 1, k) + sign(q(k)) * (EX_tmp ./ E_tmp) .* switch_state(:, k) * dl * 0.1;
            Y(:, i, k) = Y(:, i - 1, k) + sign(q(k)) * (EY_tmp ./ E_tmp) .* switch_state(:, k) * dl * 0.1;
        end
        i = i + 1;
        if i > 10000
            break
        end
    end
    if length(q) ~= length(s)
        for ii = 1:n
            k = length(q);
            if q(k) < 0 && draw_state(ii,k) == 1 || q(k) > 0
                handle{ii, k} = plot(X(ii, 1:i - 1, k), Y(ii, 1:i - 1, k), 'k-', 'LineWidth', 0.8, 'LineStyle', '-.');
            else
                handle{ii, k} = plot(X(ii, 1:i - 1, k), Y(ii, 1:i - 1, k), 'w-', 'LineWidth', 0.8, 'LineStyle', '-.');
            end
        end
        if q(k) > 0
            s{end + 1} = scatter(xq(end), yq(end), 200, 'red', 'filled');
            T{k} = text(xq(k) - 0.038, yq(k), pulse, 'FontSize', 20, 'Color', 'w');
        else
            s{end + 1} = scatter(xq(end), yq(end), 200, 'blue', 'filled');
            T{k} = text(xq(k) - 0.038, yq(k), minus, 'FontSize', 20, 'Color', 'w');
        end
    end
    for k = 1:length(q)
        xq(k) = s{k}.XData;
        yq(k) = s{k}.YData;
        T{k}.Position(1) = s{k}.XData - 0.038;
        T{k}.Position(2) = s{k}.YData;
    end
    for ii = 1:n
        for k = 1:length(q)
            handle{ii, k}.XData = X(ii, 1:i - 1, k);
            handle{ii, k}.YData = Y(ii, 1:i - 1, k);
            if q(k) < 0 && draw_state(ii,k) == 1 || q(k) > 0
                handle{ii, k}.Visible = 'on';
            else
                handle{ii, k}.Visible = 'off';
            end
        end
    end

    position_str(step,1:2,:) = [xq; yq];

    drawnow
    %frame = getframe(gcf);
    %writeVideo(video, frame)
end
%close(video)

% Set various callback functions
    function myBDCallback(src, ~)
        set(src, 'WindowButtonMotionFcn', @myBMCallback);
        function myBMCallback(~, ~)
            Cp = get(gca, 'CurrentPoint');
            for kk = 1:length(q)
                if Cp(1, 1) < xq(kk) + 0.2 && Cp(1, 1) > xq(kk) - 0.2 && Cp(1, 2) < yq(kk) + 0.2 && Cp(1, 2) > yq(kk) - 0.2 && state == 0 || state == kk
                    state = kk;
                    set(s{kk}, 'XData', Cp(1, 1), 'YData', Cp(1, 2))
                end
            end
            if Cp(1, 1) < ss{1}.XData + 0.4 && Cp(1, 1) > ss{1}.XData - 0.4 && Cp(1, 2) < ss{1}.YData + 0.4 && Cp(1, 2) > ss{1}.YData - 0.4 && state == 0
                q = [q 1];
                xq = [xq Cp(1, 1)];
                yq = [yq Cp(1, 2)];
                state = length(q);
                position_str(:,:,state) = 0;
                v_str(:,:,state) = 0;
            elseif Cp(1, 1) < ss{2}.XData + 0.4 && Cp(1, 1) > ss{2}.XData - 0.4 && Cp(1, 2) < ss{2}.YData + 0.4 && Cp(1, 2) > ss{2}.YData - 0.4 && state == 0
                q = [q -1];
                xq = [xq Cp(1, 1)];
                yq = [yq Cp(1, 2)];
                state = length(q);
                position_str(:,:,state) = 0;
                v_str(:,:,state) = 0;
            end
            for kk = 1:length(q) - (length(q) ~= length(s))
                xq(kk) = s{kk}.XData;
                yq(kk) = s{kk}.YData;
                T{kk}.Position(1) = s{kk}.XData - 0.038;
                T{kk}.Position(2) = s{kk}.YData;
            end
        end
    end

    function myBUCallback(src, ~)
        state = 0;
        set(src, 'WindowButtonMotionFcn', '');
    end

    function simulatior_break(~, evantdata)
        if strcmp(evantdata.Key, 'escape')
            fin = 1;
        end
    end
end
